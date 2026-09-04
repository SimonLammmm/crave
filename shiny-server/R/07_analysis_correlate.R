#### Correlate analyses ####

#' Resolve Correlate plot dimensions and heatmap colours from the options.
correlateStyle <- function(opts) {
  dims <- plotDims(opts, CORRELATE_DEFAULTS)
  list(
    height = dims$height,
    width  = dims$width,
    low    = opts$heatmap_low  %||% CORRELATE_DEFAULTS$heatmap_low,
    mid    = opts$heatmap_mid  %||% CORRELATE_DEFAULTS$heatmap_mid,
    high   = opts$heatmap_high %||% CORRELATE_DEFAULTS$heatmap_high,
    fill   = opts$genequery_fillscheme %||% CORRELATE_DEFAULTS$genequery_fillscheme
  )
}

#' Translate a correlateHitData() failure status into a rendered message.
hitDataMessage <- function(res) {
  switch(
    res$status,
    no_screens = msgResult("No screens found with the specified filters.\nConsider relaxing filters."),
    no_hits    = msgResult("No hits found among screens with the specified filters and cutoffs.\nConsider relaxing filters and/or cutoffs."),
    too_many_screens = msgResult(paste0(
      "Too many screens to display (", res$n, ", maximum ", CRAVE_MAX_SCREENS_PLOT,
      ").\nConsider making filters and/or cutoffs more stringent.")),
    msgResult("Nothing to plot.")
  )
}

#### Gene query ####

runGeneQuery <- function(data, screens, genes, stat, cutoff, method,
                        plotType = "Violin plot", normalise = FALSE, opts = list()) {
  res <- correlateHitData(data, screens, genes, stat, cutoff, method,
                          normalise = normalise)
  if (!identical(res$status, "ok")) {
    return(c(hitDataMessage(res), list(brush = NULL)))
  }

  style <- correlateStyle(opts)
  plotdata <- res$data
  plotdata$goi <- plotdata$gene %in% res$goi
  plotdata$sig <- hitMask(plotdata, stat, cutoff)

  p <- ggplot(plotdata, aes(
    y = score, x = FriendlyID,
    fdr = fdr, score = score, rank_synth = rank1, rank_supp = rank2,
    DaysDiff = `Days grown (diff)`, TreatmentDiff = `Treatment (diff)`,
    DoseDiff = `Dose (diff)`, KnockoutDiff = `Knockout (diff)`,
    CellLineDiff = `Cell line (diff)`, DaysRef = `Days grown (ref)`,
    TreatmentRef = `Treatment (ref)`, DoseRef = `Dose (ref)`,
    KnockoutRef = `Knockout (ref)`, CellLineRef = `Cell line (ref)`,
    Library = Library
  ))

  myGeom <- if (identical(plotType, "Boxplot")) geom_boxplot else geom_violin
  # NB the comparisons table has no plain `Cell line`; it is `Cell line (diff)`.
  fillCol <- switch(
    style$fill,
    Treatment  = "Treatment (diff)",
    Knockout   = "Knockout (diff)",
    `Cell line` = "Cell line (diff)",
    Library    = "Library",
    "Treatment (diff)"
  )
  p <- p + myGeom(na.rm = TRUE, mapping = aes(colour = .data[[fillCol]]))
  p <- p +
    geom_point(data = plotdata %>% filter(goi), aes(colour = gene, shape = sig)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") + ylab("Score")

  displaydata <- plotdata %>%
    filter(goi) %>%
    select(-`Comparison ID`, -`Experiment ID`, -Contrast) %>%
    relocate(Gene = gene, Contrast = FriendlyID, FDR = fdr, Score = score)

  list(
    p = ggplotly(p, height = style$height, width = style$width) %>%
      layout(dragmode = "select"),
    plotdata = displaydata,
    height = style$height,
    brush = displaydata %>%
      transmute(x = as.numeric(factor(Contrast)), y = Score, Gene = Gene)
  )
}

#### Clustergram ####

runClustergram <- function(data, screens, genes, stat, cutoff, method,
                           geneDendrogram = TRUE, screenDendrogram = TRUE,
                           dendroMethod = CRAVE_DENDRO_DEFAULT,
                           distMethod = CRAVE_DISTANCE_DEFAULT,
                           k = CRAVE_CLARA_K_DEFAULT,
                           normalise = FALSE, opts = list()) {
  res <- correlateHitData(data, screens, genes, stat, cutoff, method,
                          normalise = normalise)
  if (!identical(res$status, "ok")) return(hitDataMessage(res))

  style <- correlateStyle(opts)
  plotdata <- res$data
  plotdata$goi <- plotdata$gene %in% res$goi
  plotdata$sig <- if_else(hitMask(plotdata, stat, cutoff), "*", "")
  plotdata <- plotdata[plotdata$goi, , drop = FALSE]

  if (nrow(plotdata) == 0) {
    return(msgResult("None of the genes of interest were measured in these screens."))
  }

  nGenes   <- length(unique(plotdata$gene))
  nScreens <- length(unique(plotdata$FriendlyID))

  geneDendro <- NULL
  if (nGenes >= 3 && isTRUE(geneDendrogram)) {
    wide <- plotdata %>%
      tidyr::pivot_wider(id_cols = gene, names_from = FriendlyID, values_from = score)
    geneDendro <- dendrogramLayer(wide, "gene", offset = nScreens, flip = TRUE,
                                  method = dendroMethod, metric = distMethod, k = k)
  }
  plotdata$gene <- if (!is.null(geneDendro)) {
    factor(plotdata$gene, ordered = TRUE, levels = geneDendro$order)
  } else {
    factor(plotdata$gene, levels = sort(unique(plotdata$gene), decreasing = TRUE))
  }

  screenDendro <- NULL
  if (nScreens >= 3 && isTRUE(screenDendrogram)) {
    wide2 <- plotdata %>%
      tidyr::pivot_wider(id_cols = FriendlyID, names_from = gene, values_from = score)
    screenDendro <- dendrogramLayer(wide2, "FriendlyID", offset = nGenes, flip = FALSE,
                                    method = dendroMethod, metric = distMethod, k = k)
  }
  plotdata$FriendlyID <- if (!is.null(screenDendro)) {
    factor(plotdata$FriendlyID, ordered = TRUE, levels = screenDendro$order)
  } else {
    # Sort the screens actually present, rather than relying on a filter input.
    factor(plotdata$FriendlyID, ordered = TRUE, levels = sort(unique(plotdata$FriendlyID)))
  }

  # Clara assigns clusters; surface them in the results table so the blocks on the
  # axes can be related back to gene and screen identities.
  if (!is.null(geneDendro$clusters)) {
    plotdata$`Gene cluster` <- unname(geneDendro$clusters[as.character(plotdata$gene)])
  }
  if (!is.null(screenDendro$clusters)) {
    plotdata$`Screen cluster` <- unname(screenDendro$clusters[as.character(plotdata$FriendlyID)])
  }

  tileAes <- aes(
    y = gene, x = FriendlyID, fdr = fdr, fill = score,
    rank_synth = rank1, rank_supp = rank2, label = sig,
    DaysDiff = `Days grown (diff)`, TreatmentDiff = `Treatment (diff)`,
    DoseDiff = `Dose (diff)`, KnockoutDiff = `Knockout (diff)`,
    CellLineDiff = `Cell line (diff)`, DaysRef = `Days grown (ref)`,
    TreatmentRef = `Treatment (ref)`, DoseRef = `Dose (ref)`,
    KnockoutRef = `Knockout (ref)`, CellLineRef = `Cell line (ref)`,
    Library = Library
  )

  p <- ggplot() +
    geom_tile(data = plotdata, mapping = tileAes) +
    geom_text(data = plotdata, mapping = tileAes, nudge_y = -0.15) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") + ylab("") +
    divergentFill(plotdata$score, style$low, style$mid, style$high)

  # Clara returns an ordering but no tree, so the layer can legitimately be absent.
  if (!is.null(geneDendro$layer))   p <- p + geneDendro$layer
  if (!is.null(screenDendro$layer)) p <- p + screenDendro$layer

  displaydata <- plotdata %>%
    select(-`Comparison ID`, -`Experiment ID`, -Contrast) %>%
    relocate(Gene = gene, Contrast = FriendlyID, FDR = fdr, Score = score)

  list(
    p = ggplotly(p, height = style$height, width = style$width),
    plotdata = displaydata,
    height = style$height
  )
}

#### Correlation heatmap ####

runHeatmap <- function(data, screens, genes, stat, cutoff, method,
                       similarity = CRAVE_SIMILARITY_DEFAULT,
                       dendroMethod = CRAVE_DENDRO_DEFAULT,
                       distMethod = CRAVE_DISTANCE_DEFAULT,
                       k = CRAVE_CLARA_K_DEFAULT,
                       normalise = FALSE, opts = list()) {
  if (length(screens) < 2) {
    return(msgResult("At least two screens are needed for a correlation heatmap."))
  }
  wide <- correlateWideScores(data, screens, genes, stat, cutoff, method,
                              normalise = normalise)
  if (is.null(wide) || nrow(wide) < 2) {
    return(msgResult("Not enough genes to plot a heatmap.\nAdjust the genes of interest or the cutoff."))
  }

  style <- correlateStyle(opts)
  valueName <- similarity
  plotdata <- similarityLong(wide, method = similarity, valueName = valueName)
  nB <- length(unique(plotdata$`Gene B`))

  dendro <- if (length(unique(plotdata$`Gene A`)) >= 3) {
    dendrogramLayer(
      tidyr::pivot_wider(plotdata, names_from = `Gene B`, values_from = all_of(valueName)),
      "Gene A", offset = nB, flip = FALSE,
      method = dendroMethod, metric = distMethod, k = k
    )
  } else NULL

  levs <- if (!is.null(dendro)) dendro$order else sort(unique(plotdata$`Gene A`))
  plotdata$`Gene A` <- factor(plotdata$`Gene A`, ordered = TRUE, levels = levs)
  plotdata$`Gene B` <- factor(plotdata$`Gene B`, ordered = TRUE, levels = levs)

  p <- ggplot(plotdata) +
    geom_tile(aes(x = `Gene A`, y = `Gene B`, fill = .data[[valueName]])) +
    scale_y_discrete(limits = rev(levels(plotdata$`Gene B`))) +
    theme_classic() +
    xlab("") + ylab("") +
    labs(fill = valueName) +
    divergentFill(plotdata[[valueName]], style$low, style$mid, style$high)

  if (!is.null(dendro$layer)) {
    p <- p + dendro$layer +
      theme(axis.line.y = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      annotate(x = 0, xend = 0, y = 0, yend = nB + 0.5, geom = "segment")
  }

  if (!is.null(dendro$clusters)) {
    plotdata$`Cluster (Gene A)` <- unname(dendro$clusters[as.character(plotdata$`Gene A`)])
  }

  list(
    p = ggplotly(suppressWarnings(p), height = style$height, width = style$width),
    plotdata = plotdata,
    height = style$height
  )
}

#### Correlation network ####

runNetwork <- function(data, screens, genes, stat, cutoff, method,
                       similarity = CRAVE_SIMILARITY_DEFAULT,
                       corrCutoff = defaults_corr_cutoff_value,
                       normalise = FALSE) {
  if (length(screens) < 2) {
    return(list(p = msgNetwork("At least two screens are needed for a network."),
                plotdata = tibble(), renderer = "visNetwork"))
  }
  if (!needPkg("igraph", "visNetwork")) {
    return(list(p = msgNetwork("The igraph package is not installed."),
                plotdata = tibble(), renderer = "visNetwork"))
  }
  wide <- correlateWideScores(data, screens, genes, stat, cutoff, method,
                              normalise = normalise)
  if (is.null(wide) || nrow(wide) < 2) {
    return(list(p = msgNetwork("Not enough genes to build a network.\nAdjust the genes of interest or the cutoff."),
                plotdata = tibble(), renderer = "visNetwork"))
  }

  edges <- similarityLong(wide, method = similarity, nameA = "from", nameB = "to",
                          valueName = "corr") %>%
    filter(to != from, !is.na(corr), abs(corr) >= asNumericOr(corrCutoff, 0.7))

  if (nrow(edges) == 0) {
    return(list(p = msgNetwork("No interactions found.\nAdjust the screen filters, genes of interest, or correlation cutoff."),
                plotdata = tibble(), renderer = "visNetwork"))
  }

  # No set.seed() here: visNetwork lays the graph out in the browser, and seeding
  # the global RNG would make every later random draw in the process deterministic,
  # including Exorcise task ids and the "random" tie-breaking in rank().
  nodes <- unique(c(edges$from, edges$to))
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  igraph::V(g)$color.background <- "#F2E95C"
  igraph::V(g)$color.border     <- "#F0D337"
  igraph::V(g)$size             <- 15
  igraph::V(g)$borderWidth      <- 3

  p <- visNetwork::visIgraph(g) %>%
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE))

  out <- edges %>% transmute(`Gene A` = from, `Gene B` = to, value = corr)
  names(out)[3] <- similarity

  list(p = p, plotdata = out, renderer = "visNetwork")
}

#### Reduce: dimensionality reduction ####

# Minimum number of screens (that is, matrix columns) each method needs.
CRAVE_REDUCE_MIN_SCREENS <- c("UMAP" = 4L, "PCA" = 2L, "t-SNE" = 3L)

#' Project genes into two dimensions.
#'
#' @param reduceMethod One of CRAVE_REDUCE_METHODS.
#' @return list(p, plotdata, height, brush) or a message result.
runReduce <- function(data, screens, genes, stat, cutoff, method,
                      reduceMethod = CRAVE_REDUCE_DEFAULT,
                      normalise = FALSE, opts = list()) {
  reduceMethod <- reduceMethod %||% CRAVE_REDUCE_DEFAULT
  if (!reduceMethod %in% CRAVE_REDUCE_METHODS) {
    return(msgResult(paste0("Unrecognised reduction method '", reduceMethod, "'.")))
  }
  # `[` not `[[`: a missing name gives NA rather than raising.
  minScreens <- unname(CRAVE_REDUCE_MIN_SCREENS[reduceMethod])
  if (is.na(minScreens)) minScreens <- 2L
  if (length(screens) < minScreens) {
    return(msgResult(paste0("At least ", minScreens, " screens are needed for ",
                            reduceMethod, ".")))
  }

  needed <- switch(reduceMethod,
                   "UMAP"  = c("umap", "mice"),
                   "t-SNE" = c("Rtsne", "mice"),
                   "PCA"   = "mice")
  if (!needPkg(needed)) {
    return(msgResult(paste0("The ", paste(needed, collapse = " and "),
                            " package(s) are needed for ", reduceMethod, ".")))
  }

  wide <- correlateWideScores(data, screens, genes, stat, cutoff, method,
                              cap = CRAVE_MAX_AUTO_GOI_REDUCE, normalise = normalise,
                              applyCutoff = CRAVE_REDUCE_USES_CUTOFF)
  if (is.null(wide) || nrow(wide) < 4) {
    return(msgResult(paste0("Not enough genes to compute a ", reduceMethod,
                            ".\nAdjust the genes of interest or the cutoff.")))
  }
  log_info(paste0("Reduce: projecting ", nrow(wide), " genes across ",
                  ncol(wide) - 1L, " screens with ", reduceMethod, "."))

  style <- correlateStyle(opts)
  # Genes are the observations and screens the variables, for both the imputation
  # and the projection that follows it: a gene's missing score in one screen is
  # predicted from its scores in the others, fitted across all the genes. Screens
  # are renamed to syntactic identifiers because mice() builds model formulae from
  # its column names.
  geneNames <- wide$gene
  m <- as.data.frame(wide[, setdiff(names(wide), "gene"), drop = FALSE])
  names(m) <- paste0("s", seq_along(m))
  rownames(m) <- NULL

  projected <- tryCatch({
    # All three methods need a complete matrix, so imputation comes first — but only
    # if anything is actually missing, which matters now that the whole library can
    # be in play. m = 1 because complete() takes the first imputation and the other
    # four were being computed and thrown away.
    imputed <- if (anyNA(m)) {
      as.matrix(mice::complete(
        mice::mice(m, m = 1, method = "lasso.norm", printFlag = FALSE)))
    } else {
      as.matrix(m)
    }
    if (anyNA(imputed)) stop("imputation left missing values")
    switch(
      reduceMethod,
      "UMAP" = {
        cfg <- umap::umap.defaults
        cfg$n_neighbors  <- min(15L, max(2L, nrow(imputed) - 1L))
        cfg$min_dist     <- 0.3
        cfg$random_state <- 101079
        list(coords = umap::umap(imputed, config = cfg)$layout,
             labels = c("UMAP1", "UMAP2"))
      },
      "PCA" = {
        pca <- stats::prcomp(imputed, center = TRUE, scale. = FALSE)
        if (ncol(pca$x) < 2) stop("PCA produced fewer than two components")
        pct <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
        list(coords = pca$x[, 1:2, drop = FALSE],
             labels = sprintf("PC%d (%.1f%% of variance)", 1:2, pct[1:2]))
      },
      "t-SNE" = {
        n <- nrow(imputed)
        # Rtsne requires 3 * perplexity < n - 1.
        perplexity <- min(CRAVE_TSNE_MAX_PERPLEXITY, floor((n - 2) / 3))
        if (perplexity < 2) {
          stop(paste0("t-SNE needs more genes: ", n,
                      " gives a perplexity below 2. Relax the cutoff."))
        }
        ts <- withSeed(101079, Rtsne::Rtsne(
          as.matrix(imputed), dims = 2, perplexity = perplexity,
          check_duplicates = FALSE, pca = TRUE, verbose = FALSE))
        list(coords = ts$Y, labels = c("t-SNE 1", "t-SNE 2"))
      }
    )
  }, error = function(e) {
    log_warn(paste0(reduceMethod, " failed: ", conditionMessage(e)))
    NULL
  })

  if (is.null(projected)) {
    return(msgResult(paste0(
      "The ", reduceMethod, " could not be computed for this selection.\n",
      "Try more screens, or relax the cutoff to include more genes.")))
  }

  result <- tibble(
    Gene = geneNames,
    Dim1 = as.numeric(projected$coords[, 1]),
    Dim2 = as.numeric(projected$coords[, 2])
  )

  # One point per gene, even though a gene may belong to several ontology classes.
  # The ontology is collapsed to one row per gene *before* the join, so the join is
  # one-to-one by construction and `result` cannot gain rows: collapsing afterwards
  # would leave the row count depending on a group-order assumption, and the brush
  # table is built from these same rows.
  #
  # The lexicographically first class is the one plotted and coloured; the rest are
  # reported in Others so nothing is silently hidden. arrange() rather than sort()
  # because dplyr orders in the C locale, so the choice does not vary with the
  # deployment's locale settings.
  #
  # NB the input column is `cls`, not `Class`. summarise() evaluates its arguments in
  # order and each new column shadows any existing one of the same name, so
  # `summarise(Class = Class[1], Others = paste(Class[-1], ...))` reads the new
  # length-1 Class in the second expression and makes Others empty for every gene.
  # Naming the input differently means neither output can shadow it and the two
  # expressions can be written in either order.
  result$Class  <- NA_character_
  result$Others <- ""
  if (nrow(data$ontology) > 0) {
    byGene <- data$ontology %>%
      transmute(Gene = as.character(symbol), cls = as.character(class)) %>%
      filter(!is.na(cls), nzchar(cls)) %>%
      distinct() %>%
      arrange(Gene, cls) %>%
      summarise(Class  = cls[1],
                Others = paste(cls[-1], collapse = ", "),
                .by = Gene)
    result <- result %>%
      select(-Class, -Others) %>%
      left_join(byGene, by = "Gene") %>%
      mutate(Others = if_else(is.na(Others), "", Others))
  }

  p <- ggplot(result, aes(x = Dim1, y = Dim2, Gene = Gene, colour = Class,
                          Others = Others)) +
    geom_point() +
    theme_classic() +
    xlab(projected$labels[1]) +
    ylab(projected$labels[2]) +
    ggtitle(reduceMethod)

  list(
    p = ggplotly(p, height = style$height, width = style$width) %>%
      layout(dragmode = "select"),
    plotdata = result,
    height = style$height,
    brush = result %>% transmute(x = Dim1, y = Dim2, Gene = Gene)
  )
}

#### Enrichment ####

runEnrichment <- function(data, screens, method, stat, cutoff, direction, tail,
                          compartments = FALSE, normalise = FALSE) {
  if (nrow(data$ontology) == 0) {
    return(msgResult("No ontology is loaded for this dataset."))
  }
  if (!isNumericString(as.character(cutoff))) {
    return(msgResult("The cutoff must be a number."))
  }
  if (length(screens) == 0) {
    return(msgResult("No screens found with the specified filters.\nConsider relaxing filters."))
  }

  if (isTRUE(compartments)) {
    return(runEnrichmentCompartments(data, screens, method, stat, cutoff,
                                     direction, tail, normalise = normalise))
  }

  alldata    <- setAnalysisFrame(data, screens, method, stat, normalise = normalise)
  cutoffdata <- applyDirectionCutoff(alldata, stat, cutoff, direction)

  allC <- alldata %>%
    ungroup() %>%
    left_join(data$ontology, by = c("gene" = "symbol"), relationship = "many-to-many") %>%
    filter(!is.na(class))
  cutC <- cutoffdata %>%
    ungroup() %>%
    left_join(data$ontology, by = c("gene" = "symbol"), relationship = "many-to-many") %>%
    filter(!is.na(class))

  if (nrow(cutC) == 0) {
    return(msgResult("No significant hits with the given cutoff.\nConsider relaxing the cutoff."))
  }

  # x, m, N and k for every screen x class cell, as four grouped counts.
  allU <- unique(data.table::as.data.table(allC)[, .(screen, gene, class)])
  cutU <- unique(data.table::as.data.table(cutC)[, .(screen, gene, class)])

  m_tab <- allU[, .(m = .N), by = .(screen, class)]
  x_tab <- cutU[, .(x = .N), by = .(screen, class)]
  N_tab <- unique(allU[, .(screen, gene)])[, .(N = .N), by = screen]
  k_tab <- unique(cutU[, .(screen, gene)])[, .(k = .N), by = screen]

  grid <- data.table::CJ(screen = unique(allU$screen),
                         class  = unique(as.character(data$ontology$class)),
                         unique = TRUE)
  grid <- merge(grid,  m_tab, by = c("screen", "class"), all.x = TRUE)
  grid <- merge(grid,  x_tab, by = c("screen", "class"), all.x = TRUE)
  grid <- merge(grid,  N_tab, by = "screen", all.x = TRUE)
  grid <- merge(grid,  k_tab, by = "screen", all.x = TRUE)
  for (cl in c("m", "x", "N", "k")) grid[is.na(get(cl)), (cl) := 0L]
  grid[, n := N - m]

  upper <- !identical(tail, "Lower tail (depleted classes)")
  grid[, phyper := stats::phyper(x - as.integer(upper), m, n, k, lower.tail = !upper)]
  grid[, adj.phyper := signif(stats::p.adjust(phyper, method = "BH"), 4)]
  grid[, phyper := signif(phyper, 4)]

  hypers <- as_tibble(grid)

  plotdata <- cutC %>%
    group_by(screen, class) %>%
    summarise(
      genes = gsub("(.+?;.+?;.+?;.+?;.+?;.+?);(.+?)", "\\1\n\\2",
                   paste0(sort(gene), collapse = "; ")),
      .groups = "drop"
    ) %>%
    left_join(hypers, by = c("screen", "class")) %>%
    mutate(recall = if_else(m > 0, x / m, NA_real_),
           sig = if_else(!is.na(adj.phyper) & adj.phyper < 0.05, "*", ""))

  p <- ggplot(plotdata, aes(
    x = screen, y = class, fill = recall, label = sig, symbols = genes,
    phyperx = x, phyperk = k, phyperm = m, phypern = n,
    phyper = phyper, adj.phyper = adj.phyper
  )) +
    geom_tile() +
    geom_text(nudge_y = -0.3) +
    scale_fill_gradient(low = "#FFFFFF", high = "#FF0087", breaks = c(0, 1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") + ylab("")

  list(
    p = ggplotly(p),
    plotdata = plotdata %>%
      transmute(Contrast = screen, Class = class, Genes = genes,
                `Hypergeometric p-value` = phyper,
                `Hypergeometric FDR` = adj.phyper),
    height = 900
  )
}

#' Enrichment within the exclusive Venn compartments of the selected screens.
#'
#' Instead of one column per screen, there is one column per combination of screens:
#' genes hit in A alone, in A and B but not C, and so on. The x axis is drawn as an
#' upset matrix so each column states which screens it represents.
#'
#' Rendered with renderPlot rather than plotly, because ggupset's combination-matrix
#' axis has no plotly equivalent.
runEnrichmentCompartments <- function(data, screens, method, stat, cutoff,
                                      direction, tail, normalise = FALSE) {
  if (length(screens) > CRAVE_MAX_VENN_SETS) {
    return(c(msgResult(paste0(
      "Compartment enrichment is limited to ", CRAVE_MAX_VENN_SETS,
      " screens; ", length(screens), " are selected.\n",
      "Tighten the screen filters, or switch compartments off to get one column per screen.")),
      list(renderer = "plotly")))
  }
  if (length(screens) < 2) {
    return(c(msgResult("Compartment enrichment needs at least two screens."),
             list(renderer = "plotly")))
  }
  if (!needPkg("ggupset")) {
    return(c(msgResult("The ggupset package is needed to draw the compartment axis."),
             list(renderer = "plotly")))
  }

  alldata    <- setAnalysisFrame(data, screens, method, stat, normalise = normalise)
  cutoffdata <- applyDirectionCutoff(alldata, stat, cutoff, direction)
  sets       <- hitSetsByScreen(cutoffdata)

  if (length(sets) < 2) {
    return(c(msgResult("Fewer than two screens have hits at this cutoff.\nConsider relaxing the cutoff."),
             list(renderer = "plotly")))
  }

  # Compartment membership is only defined over genes every screen measured.
  measured <- measuredGenesByScreen(alldata)
  universe <- Reduce(intersect, measured)
  if (length(universe) == 0) {
    return(c(msgResult("The selected screens share no measured genes, so compartments cannot be formed."),
             list(renderer = "plotly")))
  }

  vc <- vennCompartments(sets, universe)
  if (length(vc$compartments) == 0) {
    return(c(msgResult("No hits fall in any compartment.\nConsider relaxing the cutoff."),
             list(renderer = "plotly")))
  }

  upper <- !identical(tail, "Lower tail (depleted classes)")
  enr <- classEnrichment(vc$compartments, universe, data$ontology, upperTail = upper)
  if (nrow(enr) == 0) {
    return(c(msgResult("No ontology classes overlap the compartment genes."),
             list(renderer = "plotly")))
  }

  geneList <- vapply(vc$compartments, function(g) paste(sort(g), collapse = "; "),
                     character(1))
  enr <- enr %>%
    mutate(
      group = factor(group, levels = vc$labels),
      Genes = unname(geneList[as.character(group)])
    )

  p <- ggplot(enr, aes(x = group, y = class, fill = recall, label = sig)) +
    geom_tile() +
    geom_text(nudge_y = -0.3) +
    scale_fill_gradient(low = "#FFFFFF", high = "#FF0087",
                        breaks = c(0, 1), na.value = "#F5F5F5") +
    ggupset::axis_combmatrix(sep = CRAVE_COMPARTMENT_SEP) +
    theme_classic() +
    xlab("") + ylab("") +
    labs(fill = "Recall")

  list(
    p = p,
    renderer = "plot",
    height = 900,
    plotdata = enr %>%
      transmute(
        Compartment = gsub(CRAVE_COMPARTMENT_SEP, " & ", as.character(group), fixed = TRUE),
        Class = class, Genes,
        `Hits in compartment` = k, `Hits in class` = x,
        `Class size` = m, `Universe` = N,
        `Hypergeometric p-value` = phyper,
        `Hypergeometric FDR` = adj.phyper
      )
  )
}

#### Pendragonator ####

runPendragonator <- function(data, screens, genes, method, statName, n,
                             goiOnly = FALSE, normalise = FALSE) {
  if (length(screens) == 0) {
    return(tibble(Result = "No screens found with the specified filters."))
  }
  fetched <- fetchStat(data, screens, if (isTRUE(goiOnly)) genes else NULL, method,
                       normalise = normalise)
  if (nrow(fetched) == 0) {
    return(tibble(Result = "No data found for the specified filters."))
  }

  dt <- data.table::as.data.table(fetched)
  dt <- dt[gene_id != "X" & !grepl("Non-targeting", gene_id)]
  if (nrow(dt) == 0) {
    return(tibble(Result = "No gene-level data found for the specified filters."))
  }

  # Smaller is stronger: FDR directly, or the negated absolute score.
  dt[, .stat := if (identical(statName, "FDR")) fdr else -abs(score)]
  dt <- dt[!is.na(.stat)]
  if (nrow(dt) == 0) {
    return(tibble(Result = paste0("No usable ", statName, " values for the specified filters.")))
  }

  # One best row per gene, deterministically: order, then take the first per gene.
  data.table::setorder(dt, .stat)
  dt <- dt[, .SD[1], by = gene_id]
  data.table::setorder(dt, .stat)
  dt <- utils::head(dt, max(1L, as.integer(asNumericOr(n, 2000))))

  as_tibble(dt) %>%
    left_join(data$comparisons, by = c("comparison_id" = "Comparison ID")) %>%
    transmute(
      Gene = gene_id, Score = score, FDR = fdr,
      Method = methodName(analysis_type_id),
      Contrast,
      `Days grown (diff)`, `Treatment (diff)`, `Dose (diff)`,
      `Knockout (diff)`, `Cell line (diff)`,
      `Days grown (ref)`, `Treatment (ref)`, `Dose (ref)`,
      `Knockout (ref)`, `Cell line (ref)`,
      Library, Organism, Citation, Source
    )
}

#### Bulk download ####

runBulkDownload <- function(data, screens, filename = "export.tsv.gz",
                            normalise = FALSE) {
  delim <- if (grepl("\\.(tsv|txt)(\\.|$)", filename)) "\t" else ","

  if (length(screens) == 0) {
    return(list(console = "No screens found with the specified filters.",
                table = tibble("That didn't work." = ""),
                filename = filename, delim = delim))
  }
  if (length(screens) > CRAVE_MAX_SCREENS_DOWNLOAD) {
    return(list(
      console = paste0("Query too large (", length(screens), " screens, maximum ",
                       CRAVE_MAX_SCREENS_DOWNLOAD,
                       "). Limit your search, or contact the authors if such a large query is necessary."),
      table = tibble("That didn't work." = ""),
      filename = filename, delim = delim))
  }

  table <- fetchStat(data, screens, NULL, NULL, normalise = normalise) %>%
    left_join(data$comparisons, by = c("comparison_id" = "Comparison ID")) %>%
    transmute(
      Citation, Contrast = FriendlyID, `Gene symbol` = gene_id,
      # methodName() decodes all four analysis types.
      Method = methodName(analysis_type_id),
      Score = score, FDR = fdr,
      `Upper-tail p-value` = pos_p, `Lower-tail p-value` = neg_p,
      `Days grown (diff)`, `Treatment (diff)`, `Dose (diff)`,
      `Knockout (diff)`, `Cell line (diff)`,
      `Days grown (ref)`, `Treatment (ref)`, `Dose (ref)`,
      `Knockout (ref)`, `Cell line (ref)`,
      Library, Source
    )

  found <- sort(unique(stripHtml(as.character(table$Citation))))
  console <- paste0(
    "Found ", nrow(table), " entries among ", length(found), " screens. ",
    "Screens found: ", paste(found, collapse = ", "), ". ",
    "Showing the first six entries here. Click the download button to download them all."
  )
  list(console = console, table = table, filename = filename, delim = delim)
}

#### Guide library ####

#' Filter the guide library database.
#'
#' Every user-supplied value is a bound parameter, never interpolated into the SQL.
runGuideLibrary <- function(con, filters) {
  spec <- list(
    list(values = filters$geneclasses,  sql = "[Gene Type] IN (%s)"),
    list(values = filters$targets,      sql = "[Target] IN (%s)"),
    list(values = filters$chemistry,    sql = "[Chemistry] IN (%s)"),
    list(values = filters$assembly,     sql = "[Assembly] IN (%s)"),
    list(values = filters$pam,          sql = "[PAM] IN (%s)"),
    list(values = filters$library,      sql = "[Library] IN (%s)"),
    list(values = filters$organism,     sql = "[Organism] IN (%s)")
  )

  clauses <- character(0)
  params  <- list()
  for (s in spec) {
    v <- s$values
    if (length(v) == 0) next
    clauses <- c(clauses, sprintf(s$sql, paste(rep("?", length(v)), collapse = ",")))
    params  <- c(params, as.list(as.character(v)))
  }
  if (length(filters$chromosomes) > 0) {
    clauses <- c(clauses, paste0("(",
      paste(rep("[Alignment] LIKE ?", length(filters$chromosomes)), collapse = " OR "), ")"))
    params <- c(params, as.list(paste0("chr", filters$chromosomes, ":%")))
  }

  if (length(clauses) == 0) {
    return(tibble(` ` = "Please select at least one option in the sidebar."))
  }

  sql <- paste("SELECT * FROM guideLibrary WHERE", paste(clauses, collapse = " AND "))
  out <- tryCatch(
    DBI::dbGetQuery(con, sql, params = params),
    error = function(e) {
      log_error(paste0("Guide library query failed: ", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(out)) return(tibble(` ` = "The guide library query failed."))
  if (nrow(out) == 0) return(tibble(` ` = "No guides match those filters."))

  for (cl in intersect(c("PAM", "Target", "Assembly", "Library", "Chemistry", "Organism"),
                       names(out))) {
    out[[cl]] <- factor(out[[cl]])
  }
  as_tibble(out)
}
