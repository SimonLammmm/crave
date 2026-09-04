#### Shared analysis primitives ####
# The pieces used by more than one analysis: screen filtering, hit selection, the
# preludes shared by the Correlate tabs, dendrogram construction, and the cutoff and
# hypergeometric machinery behind the set analyses.

#### Screen filtering ####

# The Correlate sidebar filter fields, mapped to the comparisons column each one
# constrains. Adding a filter now means adding one row here plus one input.
CORRELATE_FILTER_MAP <- c(
  library          = "Library",
  source           = "Source",
  organism         = "Organism",
  kind             = "Kind",
  timepoint        = "Timepoint",
  days_grown_diff  = "Days grown (diff)",
  treatment_diff   = "Treatment (diff)",
  dose_diff        = "Dose (diff)",
  knockout_diff    = "Knockout (diff)",
  cellline_diff    = "Cell line (diff)",
  days_grown_ref   = "Days grown (ref)",
  treatment_ref    = "Treatment (ref)",
  dose_ref         = "Dose (ref)",
  knockout_ref     = "Knockout (ref)",
  cellline_ref     = "Cell line (ref)",
  contrast         = "FriendlyID"
)

#' Apply the Correlate screen filters and return matching comparison IDs.
#'
#' @param data App data object.
#' @param filters Named list; names are the keys of CORRELATE_FILTER_MAP plus
#'   `citation` and `custom`. NULL or empty entries are ignored.
applyScreenFilter <- function(data, filters) {
  cmp <- data$comparisons
  if (nrow(cmp) == 0) return(character(0))

  keep <- rep(TRUE, nrow(cmp))

  # Citation is matched against the tag-stripped label, as fixed substring.
  if (length(filters$citation) && any(nzchar(filters$citation))) {
    plain <- stripHtml(as.character(cmp$Citation))
    hit <- rep(FALSE, nrow(cmp))
    for (needle in filters$citation) {
      if (nzchar(needle)) hit <- hit | grepl(needle, plain, fixed = TRUE)
    }
    keep <- keep & hit
  }

  for (key in names(CORRELATE_FILTER_MAP)) {
    want <- filters[[key]]
    if (length(want) == 0) next
    col <- CORRELATE_FILTER_MAP[[key]]
    keep <- keep & (as.character(cmp[[col]]) %in% as.character(want))
  }

  if (length(filters$custom) == 1 && nzchar(filters$custom)) {
    keep <- keep & matchTextFilter(cmp$FriendlyID, filters$custom, filters$customRegex)
  }

  unique(as.character(cmp$`Comparison ID`[keep]))
}

#' Match a free-text filter against a character vector.
#'
#' @param regex TRUE to treat `pattern` as a Perl-compatible regular expression,
#'   FALSE to match it as a literal substring.
#' @return A logical vector. An invalid regular expression matches nothing rather
#'   than raising, so a half-typed pattern narrows the result instead of breaking
#'   the tab; validateTextFilter() is what tells the user about it.
matchTextFilter <- function(x, pattern, regex = FALSE) {
  if (!isTRUE(regex)) return(grepl(pattern, x, fixed = TRUE))
  tryCatch(
    grepl(pattern, x, perl = TRUE),
    error = function(e) rep(FALSE, length(x))
  )
}

#' Check a Perl regular expression, returning NULL if valid or the error message.
validateTextFilter <- function(pattern, regex = FALSE) {
  if (!isTRUE(regex) || length(pattern) != 1 || !nzchar(pattern)) return(NULL)
  tryCatch({
    grepl(pattern, "", perl = TRUE)
    NULL
  }, error = function(e) conditionMessage(e))
}

#### Quantile normalisation ####

#' Quantile normalise scores so that every screen shares one score distribution.
#'
#' Screens seldom measure identical gene sets, so the textbook sort-the-columns-and
#' -average approach does not apply: the columns have different lengths and
#' different missing entries. Instead each screen's values are placed on the
#' (rank - 0.5) / n probability scale, a reference quantile function is formed by
#' averaging the screens' own quantile functions over a shared probability grid, and
#' every value is mapped back through that reference by interpolation.
#'
#' Normalisation is rank-preserving within a screen, so FDRs and rank-based cutoffs
#' are unaffected; only score magnitudes, and therefore score-based cutoffs, change.
#'
#' @param df Long-format statistics, one row per gene per comparison.
#' @param valueCol Column holding the score.
#' @param byCol Column identifying the screen.
quantileNormaliseScores <- function(df, valueCol = "score", byCol = "comparison_id",
                                    gridSize = 1024L) {
  if (nrow(df) == 0) return(df)
  if (!all(c(valueCol, byCol) %in% names(df))) return(df)

  dt <- data.table::as.data.table(df)
  screens <- unique(dt[[byCol]])
  if (length(screens) < 2) {
    log_info("Quantile normalisation needs at least two screens; leaving scores as they are.")
    return(df)
  }

  probs <- (seq_len(gridSize) - 0.5) / gridSize

  # Each screen's empirical quantile function, sampled on the shared grid.
  perScreen <- vapply(screens, function(s) {
    v <- dt[[valueCol]][dt[[byCol]] == s]
    v <- v[!is.na(v)]
    if (length(v) == 0) return(rep(NA_real_, gridSize))
    stats::quantile(v, probs = probs, names = FALSE, type = 7)
  }, numeric(gridSize))

  reference <- rowMeans(perScreen, na.rm = TRUE)
  if (all(is.na(reference))) {
    log_warn("Quantile normalisation found no usable scores; leaving them as they are.")
    return(df)
  }

  # Plotting position of every value within its own screen.
  dt[, .qnP := {
    v  <- .SD[[1L]]
    ok <- !is.na(v)
    out <- rep(NA_real_, .N)
    if (any(ok)) out[ok] <- (rank(v[ok], ties.method = "average") - 0.5) / sum(ok)
    out
  }, by = c(byCol), .SDcols = valueCol]

  mapped <- stats::approx(x = probs, y = reference,
                          xout = dt$.qnP[!is.na(dt$.qnP)], rule = 2)$y
  dt[!is.na(.qnP), (valueCol) := mapped]
  dt[, .qnP := NULL]

  log_info(paste0("Quantile normalised ", nrow(dt), " scores across ",
                  length(screens), " screens."))
  as_tibble(dt)
}

#### Hit selection ####

#' Add within-comparison score ranks.
#'
#' rank1 ranks ascending (strongest hypersensitivity first), rank2 descending.
addRanks <- function(df, by = "comparison") {
  if (nrow(df) == 0) {
    df$rank1 <- numeric(0)
    df$rank2 <- numeric(0)
    return(df)
  }
  dt <- data.table::as.data.table(df)
  dt[, c("rank1", "rank2") := list(rank(score), rank(-score)), by = c(by)]
  as_tibble(dt)
}

#' Logical mask of rows that count as hits under the given statistic and cutoff.
#'
#' NA statistics are never hits.
hitMask <- function(df, stat, cutoff) {
  n <- nrow(df)
  if (n == 0) return(logical(0))
  cutoff <- asNumericOr(cutoff, NA_real_)
  if (is.na(cutoff)) return(rep(FALSE, n))
  stat <- stat %||% "Score"
  out <- switch(
    stat,
    Score = !is.na(df$score) & abs(df$score) >= abs(cutoff),
    FDR   = !is.na(df$fdr)   & df$fdr <= cutoff,
    Rank  = (!is.na(df$rank1) & df$rank1 <= cutoff) |
            (!is.na(df$rank2) & df$rank2 <= cutoff),
    rep(FALSE, n)
  )
  out[is.na(out)] <- FALSE
  out
}

#' Pick genes of interest automatically when the user supplied none.
#'
#' @param df Must contain gene, score and (for Rank) rank1/rank2.
#' @param cap Maximum number of genes to keep; the strongest are kept.
autoGoi <- function(df, stat, cutoff, cap = CRAVE_MAX_AUTO_GOI) {
  hits <- df[hitMask(df, stat, cutoff), , drop = FALSE]
  if (nrow(hits) == 0) return(character(0))
  best <- hits %>%
    group_by(gene) %>%
    summarise(score = max(abs(score), na.rm = TRUE), .groups = "drop")
  if (nrow(best) > cap) {
    log_warn(paste0("Too many genes pass the cutoff of ", cutoff,
                    "; keeping the strongest ", cap, "."))
    best <- best %>% slice_max(order_by = score, n = cap, with_ties = FALSE)
  }
  unique(best$gene)
}

#' Every gene measured, capped by strength but not filtered by the cutoff.
#'
#' The counterpart to autoGoi() for analyses that want the whole library rather than
#' the hits. With an infinite cap this is just the distinct gene list, and it skips
#' the grouped maximum entirely — worth doing when the library runs to tens of
#' thousands of genes.
allGenesByStrength <- function(df, cap = Inf) {
  if (nrow(df) == 0) return(character(0))
  if (!is.finite(cap)) return(unique(df$gene))
  best <- df %>%
    group_by(gene) %>%
    summarise(score = max(abs(score), na.rm = TRUE), .groups = "drop")
  if (nrow(best) > cap) {
    log_warn(paste0("Library has ", nrow(best), " genes; keeping the strongest ",
                    cap, "."))
    best <- best %>% slice_max(order_by = score, n = cap, with_ties = FALSE)
  }
  unique(best$gene)
}

#' The prelude shared by the Gene query and Clustergram tabs.
#'
#' Finds the screens that contain at least one hit, then fetches every gene for
#' those screens and attaches the comparison metadata.
#'
#' @return list(status, screens, data, goi) where status is one of
#'   "ok", "no_screens", "no_hits", "too_many_screens".
correlateHitData <- function(data, screens, genes, stat, cutoff, method,
                             maxScreens = CRAVE_MAX_SCREENS_PLOT,
                             cap = CRAVE_MAX_AUTO_GOI,
                             normalise = FALSE) {
  if (length(screens) == 0) return(list(status = "no_screens"))

  probe <- fetchStat(data, screens, genes, method, normalise = normalise) %>%
    transmute(gene = gene_id, comparison = comparison_id, fdr = fdr, score = score) %>%
    addRanks()
  include <- unique(probe$comparison[hitMask(probe, stat, cutoff)])

  if (length(include) == 0) return(list(status = "no_hits"))
  if (length(include) > maxScreens) {
    return(list(status = "too_many_screens", n = length(include)))
  }

  fetched <- fetchStat(data, include, NULL, method, normalise = normalise) %>%
    transmute(gene = gene_id, `Comparison ID` = comparison_id, fdr = fdr, score = score) %>%
    addRanks(by = "Comparison ID") %>%
    filter(!is.na(fdr), !is.na(score), gene != "X") %>%
    left_join(data$comparisons, by = "Comparison ID")

  goi <- if (length(genes) > 0) genes else autoGoi(fetched, stat, cutoff, cap)

  list(status = "ok", screens = include, data = fetched, goi = goi)
}

#' The prelude shared by the Heatmap, Network and Reduce tabs.
#'
#' @param applyCutoff TRUE to keep only genes passing the sidebar cutoff, which is
#'   what the matrix analyses want. FALSE to keep every gene measured, subject only
#'   to `cap`; Reduce uses this so an embedding is drawn against the whole library.
#'   Ignored when the user has named genes of interest explicitly.
#' @return A wide gene x comparison score matrix as a tibble, or NULL.
correlateWideScores <- function(data, screens, genes, stat, cutoff, method,
                               cap = CRAVE_MAX_AUTO_GOI, normalise = FALSE,
                               applyCutoff = TRUE) {
  raw <- fetchStat(data, screens, genes, method, normalise = normalise)
  if (nrow(raw) == 0) return(NULL)

  wide <- raw %>%
    transmute(gene = gene_id, comparison = comparison_id, score = score) %>%
    tidyr::pivot_wider(id_cols = gene, names_from = comparison, values_from = score)

  if (length(genes) == 0) {
    ranked <- raw %>%
      transmute(gene = gene_id, comparison = comparison_id, fdr = fdr, score = score) %>%
      filter(gene != "X")
    gois <- if (isTRUE(applyCutoff)) {
      # addRanks() is only needed for a Rank cutoff, and it is a grouped rank over
      # every row, so it is skipped when the cutoff is not being applied at all.
      autoGoi(addRanks(ranked), stat, cutoff, cap)
    } else {
      allGenesByStrength(ranked, cap)
    }
    wide <- wide %>% filter(gene %in% gois)
  }
  wide <- wide %>% filter(gene != "X")
  if (nrow(wide) == 0) return(NULL)
  wide
}

#### Gene-gene similarity ####

#' Pairwise cosine similarity between the columns of a matrix.
#'
#' Missing values are handled pairwise, as cor(use = "pairwise.complete.obs") does.
#' Complete matrices take a vectorised path; only the ragged case falls back to the
#' pairwise loop.
cosineSimilarity <- function(m) {
  labels <- colnames(m)
  p <- ncol(m)

  if (!anyNA(m)) {
    nrm <- sqrt(colSums(m * m))
    out <- crossprod(m) / outer(nrm, nrm)
    out[!is.finite(out)] <- NA_real_
    dimnames(out) <- list(labels, labels)
    return(out)
  }

  out <- matrix(NA_real_, p, p, dimnames = list(labels, labels))
  for (i in seq_len(p)) {
    for (j in i:p) {
      ok <- !is.na(m[, i]) & !is.na(m[, j])
      if (!any(ok)) next
      a <- m[ok, i]
      b <- m[ok, j]
      den <- sqrt(sum(a * a)) * sqrt(sum(b * b))
      v <- if (den > 0) sum(a * b) / den else NA_real_
      out[i, j] <- v
      out[j, i] <- v
    }
  }
  out
}

#' Gene-gene similarity matrix for a wide gene x comparison table.
#'
#' @param method One of CRAVE_SIMILARITY_METHODS.
similarityMatrix <- function(wide, method = CRAVE_SIMILARITY_DEFAULT) {
  genes <- wide$gene
  m <- t(as.matrix(wide[, setdiff(names(wide), "gene"), drop = FALSE]))
  colnames(m) <- genes

  if (identical(method, "Cosine similarity")) return(cosineSimilarity(m))

  if (identical(method, "Kendall") && length(genes) > CRAVE_KENDALL_WARN_GENES) {
    log_warn(paste0("Kendall's tau over ", length(genes),
                    " genes is quadratic in the number of pairs and will be slow."))
  }
  suppressWarnings(cor(m, use = "pairwise.complete.obs", method = tolower(method)))
}

#' Gene-gene similarity in long form.
similarityLong <- function(wide, method = CRAVE_SIMILARITY_DEFAULT,
                           nameA = "Gene A", nameB = "Gene B",
                           valueName = "Similarity") {
  genes <- wide$gene
  cm <- similarityMatrix(wide, method)
  out <- as_tibble(cm)
  names(out) <- genes
  out[[nameA]] <- genes
  tidyr::pivot_longer(out, -all_of(nameA), names_to = nameB, values_to = valueName)
}

#### Dendrograms ####

#' Build a dendrogram layer positioned to sit alongside a heatmap.
#'
#' Consolidates three copies of the same twenty-line construction: cluster with
#' agnes, convert to ggdendro segments, extend the top branch into a root stub,
#' rescale into the heatmap's coordinate space, and optionally transpose.
#'
#' @param wide Wide table whose first column holds the labels.
#' @param labelCol Name of the label column.
#' @param offset Extent of the heatmap along the axis the dendrogram sits beside.
#' @param flip TRUE to draw the dendrogram rotated and transpose its segments,
#'   which is what a dendrogram on the y (gene) axis needs.
#' @param method One of CRAVE_DENDRO_METHODS.
#' @param metric "euclidean", "manhattan", or (Clara only) "jaccard".
#' @param k Number of clusters, used only by Clara. Must satisfy 2 <= k < nrow.
#' @return list(layer, order, clusters) or NULL if clustering was not possible.
#'   `clusters` is a named integer vector of cluster assignments under Clara and
#'   NULL otherwise.
dendrogramLayer <- function(wide, labelCol, offset, flip = FALSE,
                            method = CRAVE_DENDRO_DEFAULT,
                            metric = CRAVE_DISTANCE_DEFAULT,
                            k = NULL,
                            rootLength = 2, heightScale = 0.1) {
  if (!needPkg("cluster", "ggdendro")) return(NULL)
  tryCatch({
    m <- as.data.frame(wide)
    labels <- as.character(m[[labelCol]])
    m[[labelCol]] <- NULL
    m[is.na(m)] <- 0
    rownames(m) <- labels
    if (nrow(m) < 3 || ncol(m) < 1) return(NULL)

    built <- if (identical(method, "Clara")) {
      claraTree(m, labels, k, metric)
    } else {
      agnesTree(m, method, metric)
    }
    if (is.null(built)) return(NULL)

    # Clara with fewer than three clusters yields an ordering but no tree.
    if (is.null(built$segments)) {
      return(list(layer = NULL, order = built$order, clusters = built$clusters))
    }

    dendro <- built$segments
    maxY <- max(dendro$y, dendro$yend)
    top  <- dendro[dendro$y == dendro$yend & dendro$y == maxY, , drop = FALSE]
    mid  <- mean(c(top$x, top$xend))
    root <- data.frame(x = mid, xend = mid, y = maxY, yend = maxY + rootLength)

    dd <- rbind(dendro[, c("x", "y", "xend", "yend")],
                root[,     c("x", "y", "xend", "yend")])
    scaleTo <- max(dd$y, dd$yend)
    dd$y    <- dd$y    / scaleTo * offset * heightScale + offset + 0.5
    dd$yend <- dd$yend / scaleTo * offset * heightScale + offset + 0.5
    if (flip) {
      dd <- data.frame(x = dd$y, y = dd$x, xend = dd$yend, yend = dd$xend)
    }
    layer <- built$proto
    layer$data <- dd
    list(layer = layer, order = built$order, clusters = built$clusters)
  }, error = function(e) {
    log_warn(paste0("Dendrogram could not be built: ", conditionMessage(e)))
    NULL
  })
}

#' Hierarchical clustering with agnes.
#'
#' @return list(segments, proto, order, clusters)
agnesTree <- function(m, method, metric) {
  agnesName <- agnesMethod(method)
  if (is.null(agnesName)) {
    log_warn(paste0("Unrecognised clustering method '", method, "'."))
    return(NULL)
  }
  if (!metric %in% CRAVE_DISTANCE_AGNES) {
    log_warn(paste0("agnes() does not support the '", metric,
                    "' metric; using euclidean instead."))
    metric <- "euclidean"
  }
  args <- list(x = m, method = agnesName, metric = metric)
  # Only the flexible (Lance-Williams) linkage takes a parameter, and it is
  # mandatory when that linkage is chosen.
  if (identical(agnesName, "flexible")) args$par.method <- CRAVE_AGNES_FLEXIBLE_PAR
  ag <- do.call(cluster::agnes, args)

  gd <- ggdendro::ggdendrogram(stats::as.dendrogram(ag))
  list(segments = gd$layers[[2]]$data,
       proto    = gd$layers[[2]],
       order    = ag$order.lab,
       clusters = NULL)
}

#' Partitioning with Clara, plus a dendrogram over the resulting medoids.
#'
#' Clara assigns clusters but builds no tree, so there is nothing to draw directly.
#' Each cluster is given a contiguous block of the axis (alphabetical within the
#' block), the k medoids are clustered hierarchically to order and relate the blocks,
#' and the resulting k-leaf tree is stretched across the axis so each leaf sits above
#' the centre of its block.
claraTree <- function(m, labels, k, metric) {
  kk <- suppressWarnings(as.integer(k %||% NA))
  if (is.na(kk) || kk < 2 || kk >= nrow(m)) {
    log_warn(paste0(
      "Clara needs 2 <= k < ", nrow(m), " for this axis but k = ",
      k %||% "unset", "; drawing no dendrogram on it."))
    return(NULL)
  }
  if (!metric %in% CRAVE_DISTANCE_CLARA) metric <- "euclidean"

  cl <- cluster::clara(m, k = kk, metric = metric, stand = FALSE,
                       samples = 50, pamLike = TRUE)
  assignment <- setNames(as.integer(cl$clustering), labels)

  medoids <- as.data.frame(cl$medoids)
  rownames(medoids) <- paste0("cluster", seq_len(nrow(medoids)))

  # A tree needs at least three leaves; below that, order clusters by number.
  if (nrow(medoids) < 3) {
    ord <- order(assignment, labels)
    return(list(segments = NULL, proto = NULL,
                order = labels[ord], clusters = assignment))
  }

  medoidMetric <- if (metric %in% CRAVE_DISTANCE_AGNES) metric else "euclidean"
  ag <- cluster::agnes(medoids, method = "ward", metric = medoidMetric)
  clusterSeq <- ag$order  # cluster ids, in the order the medoid tree puts them

  blocks <- lapply(clusterSeq, function(cid) sort(labels[assignment == cid]))
  itemOrder <- unlist(blocks, use.names = FALSE)
  sizes <- vapply(blocks, length, integer(1))
  # Centre of each cluster's block on the item axis.
  centres <- cumsum(sizes) - (sizes - 1) / 2

  gd <- ggdendro::ggdendrogram(stats::as.dendrogram(ag))
  seg <- gd$layers[[2]]$data
  # Leaf i of the medoid tree sits at x = i; stretch that onto the block centres.
  # Internal nodes are means of leaf positions, so a monotone interpolation over the
  # same mapping keeps them above their children.
  remap <- stats::approxfun(seq_along(centres), centres, rule = 2)
  seg$x    <- remap(seg$x)
  seg$xend <- remap(seg$xend)

  list(segments = seg, proto = gd$layers[[2]],
       order = itemOrder, clusters = assignment)
}

#### Divergent fill scales ####

#' Choose a fill scale that keeps zero at the midpoint colour.
divergentFill <- function(values, low, mid, high) {
  values <- values[!is.na(values)]
  if (length(values) == 0) return(scale_fill_gradient2(low = low, mid = mid, high = high))
  lo <- min(values); hi <- max(values)
  if (hi < 0 || lo > 0) {
    scale_fill_gradient2(low = low, mid = mid, high = high)
  } else {
    scale_fill_gradientn(colours = c(low, mid, high),
                         values = scales::rescale(c(lo, 0, hi)))
  }
}

#### Cutoff sets for Overlap and Enrichment ####

#' Reshape fetched statistics into the per-screen frame the set analyses expect.
#'
#' `rank` is computed once and reused, so the `rank` column and the Rank statistic
#' always describe the same ordering despite ties.method = "random".
setAnalysisFrame <- function(data, comparisonIds, method, stat, normalise = FALSE) {
  stat <- stat %||% "Score"
  base <- fetchStat(data, comparisonIds, NULL, method, normalise = normalise) %>%
    left_join(data$comparisons, by = c("comparison_id" = "Comparison ID")) %>%
    mutate(screen = FriendlyID) %>%
    group_by(screen) %>%
    mutate(
      sgn = sign(score),
      rnk = rank(score, ties.method = "random"),
      nn  = n()
    ) %>%
    ungroup()

  base$val <- switch(
    stat,
    Score = base$score,
    FDR   = base$fdr,
    Rank  = abs(base$rnk) - base$nn / 2,
    base$score
  )

  base %>%
    transmute(screen, gene = gene_id, sign = sgn, rank = rnk, value = val) %>%
    filter(gene != "X") %>%
    group_by(screen)
}

#' Apply a direction and a cutoff to a per-screen frame.
#'
#' `alldata` must still be grouped by screen so that max(rank) and slice_max()
#' operate within each screen.
applyDirectionCutoff <- function(alldata, stat, cutoff, direction) {
  stat   <- stat %||% "Score"
  cutoff <- asNumericOr(cutoff, NA_real_)
  if (is.na(cutoff)) return(alldata %>% filter(FALSE))

  out <- if (identical(stat, "Rank")) {
    switch(direction %||% "",
      `Hypersensitivity hits` = alldata %>% filter(rank < max(rank) / 2),
      `Suppressing hits`      = alldata %>% filter(rank > max(rank) / 2),
      alldata)
  } else {
    switch(direction %||% "",
      `Hypersensitivity hits` = alldata %>% filter(sign == -1),
      `Suppressing hits`      = alldata %>% filter(sign == 1),
      alldata)
  }

  topN <- max(1L, as.integer(round(cutoff, 0)))
  switch(
    stat,
    Score = out %>% filter(!is.na(value) & abs(value) >= abs(cutoff)),
    FDR   = out %>% filter(!is.na(value) & value <= cutoff),
    Rank  = out %>% slice_max(abs(value), n = topN, with_ties = FALSE),
    out
  )
}

#' Split a cutoff frame into one gene set per screen.
hitSetsByScreen <- function(cutoffdata) {
  if (nrow(cutoffdata) == 0) return(list())
  d <- data.table::as.data.table(cutoffdata)
  sets <- split(d$gene, d$screen)
  lapply(sets, unique)
}

#' Pairwise hypergeometric tests between hit sets.
#'
#' The universe for each pair is the set of genes measured in both screens, and both
#' the overlap count and the population sizes are taken from within it. The
#' upper-tail probability is P(X >= q), matching the Enrichment tab.
pairwiseHyper <- function(sets, measuredByScreen) {
  nms <- names(sets)
  k <- length(nms)
  if (k < 2) return(tibble())

  rows <- vector("list", k * (k - 1) / 2)
  idx <- 0L
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i <= j) next
      idx <- idx + 1L
      universe <- intersect(measuredByScreen[[nms[i]]], measuredByScreen[[nms[j]]])
      hits_i <- intersect(sets[[i]], universe)
      hits_j <- intersect(sets[[j]], universe)
      nn <- length(universe)
      hk <- length(hits_i)
      hm <- length(hits_j)
      hn <- nn - hm
      hq <- length(intersect(hits_i, hits_j))
      hp <- if (nn == 0 || hk == 0 || hm == 0) NA_real_ else
        stats::phyper(hq - 1, hm, hn, hk, lower.tail = FALSE)
      rows[[idx]] <- tibble(
        `Screen 1`               = c(nms[i], nms[j]),
        `Screen 2`               = c(nms[j], nms[i]),
        `Hits 1`                 = c(hk, hm),
        `Hits 2`                 = c(hm, hk),
        `Hits 1 & 2`             = hq,
        `Hypergeometric p-value` = hp,
        `Total genes`            = nn
      )
    }
  }
  out <- bind_rows(rows)
  if (nrow(out) == 0) return(out)
  out %>%
    arrange(`Screen 1`, `Screen 2`) %>%
    mutate(
      `Hypergeometric FDR`     = signif(p.adjust(`Hypergeometric p-value`, method = "BH"), 4),
      `Hypergeometric p-value` = signif(`Hypergeometric p-value`, 4)
    ) %>%
    relocate(`Screen 1`, `Screen 2`, `Hits 1`, `Hits 2`, `Hits 1 & 2`,
             `Hypergeometric p-value`, `Hypergeometric FDR`, `Total genes`)
}

#' Genes measured in each screen, for use as hypergeometric universes.
measuredGenesByScreen <- function(alldata) {
  if (nrow(alldata) == 0) return(list())
  d <- data.table::as.data.table(alldata)
  lapply(split(d$gene, d$screen), unique)
}

#### Venn compartments ####

#' Partition hit genes into exclusive Venn compartments.
#'
#' A compartment is the set of genes that are hits in exactly one combination of
#' screens: with screens A, B and C the compartments are A only, B only, C only,
#' A and B, A and C, B and C, and A, B and C. Genes hit in no screen belong to no
#' compartment.
#'
#' Membership is only meaningful over a common gene set, so compartments are
#' restricted to `universe` — the genes measured in every screen under comparison.
#'
#' @param sets Named list of hit gene sets, one per screen.
#' @param universe Character vector of genes eligible for compartment membership.
#' @param sep Separator joining screen names into a compartment label.
#' @return list(compartments = named list of gene vectors, labels = ordered labels)
vennCompartments <- function(sets, universe, sep = CRAVE_COMPARTMENT_SEP) {
  screens <- names(sets)
  if (length(screens) == 0) return(list(compartments = list(), labels = character(0)))

  eligible <- intersect(unique(unlist(sets, use.names = FALSE)), universe)
  if (length(eligible) == 0) return(list(compartments = list(), labels = character(0)))

  # membership[i, j] is TRUE when gene i is a hit in screen j.
  membership <- vapply(sets, function(s) eligible %in% s, logical(length(eligible)))
  if (length(eligible) == 1) membership <- matrix(membership, nrow = 1,
                                                 dimnames = list(NULL, screens))

  label <- apply(membership, 1, function(row) paste(screens[row], collapse = sep))
  keep <- nzchar(label)
  if (!any(keep)) return(list(compartments = list(), labels = character(0)))

  compartments <- split(eligible[keep], label[keep])

  # Order by how many screens the compartment spans, then by size: the exclusive
  # single-screen compartments come first, the deepest intersection last.
  degree <- vapply(names(compartments),
                   function(l) length(strsplit(l, sep, fixed = TRUE)[[1]]),
                   integer(1))
  ord <- order(degree, -vapply(compartments, length, integer(1)), names(compartments))

  list(compartments = compartments[ord], labels = names(compartments)[ord])
}

#' Hypergeometric enrichment of ontology classes within each group of genes.
#'
#' Shared by the per-screen and the per-compartment forms of the Enrichment tab.
#'
#' @param groups Named list of gene sets to test.
#' @param universe Genes eligible for the test; the population.
#' @param ontology Two-column table of symbol and class.
#' @param upperTail TRUE for enrichment, FALSE for depletion.
#' @return A tibble with one row per group x class.
classEnrichment <- function(groups, universe, ontology, upperTail = TRUE) {
  if (length(groups) == 0 || length(universe) == 0) return(tibble())

  ont <- data.table::as.data.table(ontology)[, .(symbol = as.character(symbol),
                                                 class  = as.character(class))]
  ont <- unique(ont[symbol %in% universe & !is.na(class)])
  if (nrow(ont) == 0) return(tibble())

  classes <- sort(unique(ont$class))
  N <- length(unique(ont$symbol))
  classMembers <- split(ont$symbol, ont$class)

  rows <- lapply(names(groups), function(g) {
    hits <- intersect(groups[[g]], ont$symbol)
    k <- length(hits)
    tibble(
      group = g,
      class = classes,
      x = vapply(classes, function(cl) length(intersect(hits, classMembers[[cl]])), integer(1)),
      m = vapply(classes, function(cl) length(classMembers[[cl]]), integer(1)),
      k = k
    )
  })

  out <- bind_rows(rows) %>%
    mutate(
      N = N,
      n = N - m,
      phyper = stats::phyper(x - as.integer(upperTail), m, n, k,
                             lower.tail = !upperTail)
    )
  out %>%
    mutate(
      adj.phyper = signif(stats::p.adjust(phyper, method = "BH"), 4),
      phyper     = signif(phyper, 4),
      recall     = if_else(m > 0, x / m, NA_real_),
      sig        = if_else(!is.na(adj.phyper) & adj.phyper < 0.05, "*", "")
    )
}
