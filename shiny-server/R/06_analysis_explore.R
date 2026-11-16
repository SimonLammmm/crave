#### Explore analyses ####
# Volcano, Rank, Biplot, Overlap and ROC. Each is a pure function of the app data
# and an explicit argument list, rather than of Shiny's `input` object, so they
# can be called and tested outside a session and are unaffected by module
# namespacing.

#' Resolve Explore plot styling from the customisation options.
exploreStyle <- function(opts) {
  colours <- c(
    opts$point_colour %||% EXPLORE_DEFAULTS$point_colour,
    opts$goi_colour   %||% EXPLORE_DEFAULTS$goi_colour
  )
  width <- asNumericOr(opts$width, EXPLORE_DEFAULTS$width)
  if (width == 0) width <- NULL
  height <- asNumericOr(opts$height, EXPLORE_DEFAULTS$height)
  if (height == 0) height <- EXPLORE_DEFAULTS$height
  list(colours = colours, width = width, height = height)
}

#' A scattergl of one screen: shared by the Volcano and Rank plots.
#'
#' @param plotdata Must contain `goi` plus the columns named by xvar/yvar.
explorePointPlot <- function(plotdata, xvar, yvar, style,
                             title, xtitle, ytitle, xtype, ytype,
                             yreverse = FALSE) {
  yaxis <- list(title = ytitle, type = ytype)
  if (yreverse) yaxis$autorange <- "reversed"
  plot_ly(
    plotdata,
    x = plotdata[[xvar]], y = plotdata[[yvar]],
    hovertext = plotdata$hovertext,
    color = plotdata$goi,
    colors = style$colours[seq_len(length(unique(plotdata$goi)))],
    type = "scattergl", mode = "markers",
    width = style$width, height = style$height
  ) %>%
    layout(
      title  = title,
      xaxis  = list(title = xtitle, type = xtype),
      yaxis  = yaxis,
      legend = list(title = list(text = "<b>GOI?</b>"), x = 0.1, y = 0.9),
      dragmode = "select"
    )
}

#### Volcano ####

plotVolcano <- function(data, contrast, method, goiGenes = character(0), opts = list()) {
  if (length(contrast) != 1 || !nzchar(contrast) ||
      length(method) != 1 || !nzchar(method)) {
    return(c(msgResult("Select a screen."), list(brush = NULL)))
  }
  ids <- comparisonIdsFor(data, contrast)
  raw <- fetchStat(data, ids, NULL, method)
  if (nrow(raw) == 0) {
    return(c(msgResult("No data for this screen with the selected method."), list(brush = NULL)))
  }

  style <- exploreStyle(opts)
  plotdata <- raw %>%
    transmute(gene = gene_id, fdr = fdr, score = score) %>%
    filter(!is.na(fdr), !is.na(score), gene != "X") %>%
    mutate(
      goi = if_else(gene %in% goiGenes, "Yes", "No"),
      hovertext = paste0(gene, "\nScore: ", score, "\nFDR: ", fdr, "\nGOI?: ", goi)
    ) %>%
    arrange(desc(goi), score)

  p <- explorePointPlot(
    plotdata, "score", "fdr", style,
    title  = contrastFormatter(contrast),
    xtitle = paste0(method, " score"),
    ytitle = paste0(method, " FDR"),
    xtype  = axisType(opts$x_axis_log, allow_log = FALSE, auto = "linear"),
    ytype  = axisType(opts$y_axis_log, allow_log = TRUE,  auto = "log"),
    yreverse = TRUE
  )

  list(
    p = p,
    plotdata = plotdata %>% transmute(Gene = gene, Score = score, FDR = fdr, `GOI?` = goi),
    height = style$height,
    brush = plotdata %>% transmute(x = score, y = fdr, Gene = gene)
  )
}

#### Rank ####

plotRank <- function(data, contrast, method, goiGenes = character(0), opts = list()) {
  if (length(contrast) != 1 || !nzchar(contrast) ||
      length(method) != 1 || !nzchar(method)) {
    return(c(msgResult("Select a screen."), list(brush = NULL)))
  }
  ids <- comparisonIdsFor(data, contrast)
  raw <- fetchStat(data, ids, NULL, method)
  if (nrow(raw) == 0) {
    return(c(msgResult("No data for this screen with the selected method."), list(brush = NULL)))
  }

  style <- exploreStyle(opts)
  plotdata <- raw %>%
    transmute(gene = gene_id, fdr = fdr, score = score) %>%
    filter(!is.na(fdr), !is.na(score), gene != "X") %>%
    mutate(
      rank = rank(score, ties.method = "random"),
      goi  = if_else(gene %in% goiGenes, "Yes", "No"),
      hovertext = paste0(gene, "\nRank: ", rank, "\nScore: ", score,
                         "\nFDR: ", fdr, "\nGOI?: ", goi)
    ) %>%
    arrange(desc(goi), score)

  p <- explorePointPlot(
    plotdata, "rank", "score", style,
    title  = contrastFormatter(contrast),
    xtitle = "Rank",
    ytitle = paste0(method, " score"),
    xtype  = axisType(opts$x_axis_log, allow_log = TRUE,  auto = "linear"),
    ytype  = axisType(opts$y_axis_log, allow_log = FALSE, auto = "linear")
  )

  list(
    p = p,
    plotdata = plotdata %>% transmute(Gene = gene, Rank = rank, Score = score,
                                      FDR = fdr, `GOI?` = goi),
    height = style$height,
    brush = plotdata %>% transmute(x = rank, y = score, Gene = gene)
  )
}

#### Biplot ####

plotBiplot <- function(data, contrastX, contrastY, methodX, methodY,
                       goiGenes = character(0), opts = list()) {
  valid <- length(contrastX) == 1 && nzchar(contrastX) &&
           length(contrastY) == 1 && nzchar(contrastY) &&
           length(methodX)   == 1 && nzchar(methodX) &&
           length(methodY)   == 1 && nzchar(methodY) &&
           !identical(c(contrastX, methodX), c(contrastY, methodY))
  if (!valid) {
    return(c(msgResult("Select two screens."), list(brush = NULL)))
  }

  xs <- fetchStat(data, comparisonIdsFor(data, contrastX), NULL, methodX) %>%
    transmute(gene = gene_id, fdrX = fdr, scoreX = score)
  ys <- fetchStat(data, comparisonIdsFor(data, contrastY), NULL, methodY) %>%
    transmute(gene = gene_id, fdrY = fdr, scoreY = score)

  plotdata <- inner_join(xs, ys, by = "gene") %>%
    filter(!is.na(fdrX), !is.na(scoreX), !is.na(fdrY), !is.na(scoreY),
           gene != "X", !grepl("Non-targeting", gene)) %>%
    mutate(
      goi = if_else(gene %in% goiGenes, "Yes", "No"),
      hovertext = paste0(
        gene,
        "\n", methodX, " score (x-axis): ", scoreX,
        "\n", methodX, " FDR (x-axis): ",   fdrX,
        "\n", methodY, " score (y-axis): ", scoreY,
        "\n", methodY, " FDR (y-axis): ",   fdrY,
        "\nGOI?: ", goi
      )
    )

  if (nrow(plotdata) == 0) {
    return(c(msgResult("No genes are shared between these two screens."), list(brush = NULL)))
  }

  style <- exploreStyle(opts)
  square <- opts$y_equals_x %||% EXPLORE_DEFAULTS$y_equals_x
  limits <- if (isTRUE(square)) {
    c(min(plotdata$scoreX, plotdata$scoreY), max(plotdata$scoreX, plotdata$scoreY))
  } else {
    c(max(min(plotdata$scoreX), min(plotdata$scoreY)),
      min(max(plotdata$scoreX), max(plotdata$scoreY)))
  }

  p <- plot_ly(
    plotdata, x = ~scoreX, y = ~scoreY, hovertext = ~hovertext, color = ~goi,
    colors = style$colours[seq_len(length(unique(plotdata$goi)))],
    type = "scattergl", mode = "markers",
    width = style$width, height = style$height
  ) %>%
    layout(
      shapes = list(list(type = "line", xref = "x", yref = "y",
                         x0 = limits[1], x1 = limits[2],
                         y0 = limits[1], y1 = limits[2],
                         line = list(color = "black"))),
      title = paste0("<i>x</i>: ", contrastFormatter(contrastX),
                     "\n<i>y</i>: ", contrastFormatter(contrastY)),
      xaxis = list(title = paste0(methodX, " score: ", contrastX),
                   type = axisType(opts$x_axis_log, allow_log = FALSE)),
      yaxis = list(title = paste0(methodY, " score: ", contrastY),
                   type = axisType(opts$y_axis_log, allow_log = FALSE)),
      legend = list(title = list(text = "<b>GOI?</b>"), x = 0.1, y = 0.9),
      dragmode = "select"
    )

  brush <- plotdata %>% transmute(x = scoreX, y = scoreY, Gene = gene)

  table <- plotdata %>%
    transmute(gene, scoreX, fdrX, scoreY, fdrY,
              delta = round(abs(scoreX - scoreY), 3),
              rmsd  = round(sqrt(abs(scoreX^2 - scoreY^2)), 3),
              goi) %>%
    arrange(desc(delta))
  names(table) <- c(
    "Gene",
    paste0(methodX, " score: ", contrastX), paste0(methodX, " FDR: ", contrastX),
    paste0(methodY, " score: ", contrastY), paste0(methodY, " FDR: ", contrastY),
    "Delta(score)", "RMSD(score)", "GOI?"
  )

  list(p = p, plotdata = table, height = style$height, brush = brush)
}

#### Overlap: Venn and Upset ####

#' Overlap analysis between two or more screens.
#'
#' @param style "Venn diagram" or "Upset plot".
#' @return list(p, plotdata, height, renderer) where renderer is "plotly" or "plot".
plotOverlap <- function(data, contrasts, method, stat, cutoff, direction,
                        style = "Venn diagram") {
  ids <- comparisonIdsFor(data, contrasts)
  if (length(ids) < 2) {
    return(c(msgResult("Select two or more screens."), list(renderer = "plotly")))
  }

  alldata    <- setAnalysisFrame(data, ids, method, stat)
  cutoffdata <- applyDirectionCutoff(alldata, stat, cutoff, direction)
  sets       <- hitSetsByScreen(cutoffdata)

  if (length(sets) < 2) {
    return(c(msgResult("Fewer than two screens have hits at this cutoff.\nConsider relaxing the cutoff."),
             list(renderer = "plotly")))
  }

  hypers <- pairwiseHyper(sets, measuredGenesByScreen(alldata))

  if (identical(style, "Upset plot")) {
    if (!needPkg("ggupset")) {
      return(c(msgResult("The ggupset package is not installed."), list(renderer = "plotly")))
    }
    wide <- cutoffdata %>%
      ungroup() %>%
      select(screen, gene) %>%
      mutate(value = TRUE) %>%
      tidyr::pivot_wider(names_from = "screen", values_from = "value")
    cols <- intersect(names(sets), names(wide))
    mat <- as.matrix(wide[, cols, drop = FALSE])
    mat[is.na(mat)] <- FALSE
    # lapply rather than apply(): scale_x_upset() needs a list column, and apply()
    # collapses to a plain vector whenever every row has the same number of members.
    wide$combination <- lapply(seq_len(nrow(mat)), function(i) cols[as.logical(mat[i, ])])
    p <- ggplot(wide, aes(x = combination)) +
      geom_bar() +
      ggupset::scale_x_upset() +
      theme_classic() +
      xlab("") + ylab("Hits")
    return(list(p = p, plotdata = hypers, height = 900, renderer = "plot"))
  }

  if (length(sets) > CRAVE_MAX_VENN_SETS) {
    return(c(msgResult(paste0(
      "Too many screens for a Venn diagram (", length(sets), " selected, maximum ",
      CRAVE_MAX_VENN_SETS, ").\nSwitch the plot type to Upset.")),
      list(renderer = "plotly")))
  }
  if (!needPkg("ggVennDiagram")) {
    return(c(msgResult("The ggVennDiagram package is not installed."), list(renderer = "plotly")))
  }

  gg <- ggVennDiagram::ggVennDiagram(sets, label = "both", show_intersect = TRUE,
                                     force_upset = FALSE)
  # ggVennDiagram with show_intersect returns a plotly object; recolour it and
  # break the set labels onto two lines.
  for (i in seq_along(gg$x$data)) {
    gg$x$data[[i]]$fillcolor <- "rgba(255,255,255,1)"
    if (all(gg$x$data[[i]]$text %in% names(sets))) {
      gg$x$data[[i]]$textfont$size <- 16
      gg$x$data[[i]]$text <- sub(": ", ":\n", gg$x$data[[i]]$text)
    }
  }
  list(p = gg, plotdata = hypers, height = 900, renderer = "plotly")
}

#### ROC ####

#' Sweep order for the query screens: the most base-hit-like gene first.
#'
#' The base screen's direction orients the sweep, because it has to -- a screen that
#' enriched a gene the base screen depleted has not recovered it, and an undirected
#' order would say it had. It never removes a gene from the sweep.
rocSweepKey <- function(df, stat, direction) {
  down <- identical(direction, "Hypersensitivity hits")
  up   <- identical(direction, "Suppressing hits")
  if (identical(stat, "FDR")) {
    # An FDR carries no direction of its own, so on a directed question a gene scored
    # the wrong way is ordered by how confidently it went the wrong way: 1 is the
    # hinge where the two directions meet and neither is evidence of anything.
    wrong <- if (down) df$score > 0 else if (up) df$score < 0 else rep(FALSE, nrow(df))
    ifelse(!is.na(wrong) & wrong, 2 - df$fdr, df$fdr)
  } else if (identical(stat, "Rank")) {
    if (down) df$rank1 else if (up) df$rank2 else pmin(df$rank1, df$rank2)
  } else {
    if (down) df$score else if (up) -df$score else -abs(df$score)
  }
}

#' Sweep one query screen gene by gene and return the curve it traces.
#'
#' @param key Sweep order from rocSweepKey(), smallest first.
#' @param inBase TRUE where the gene is one of the base screen's hits.
#' @return NULL if the screen holds no positives or no negatives, otherwise the sweep
#'   order, the per-gene rates in that order, the counts behind them, and the
#'   de-duplicated points the AUC is taken over.
rocSweep <- function(key, inBase) {
  # A gene the query screen could not score is swept last rather than dropped: it was
  # measured and not recovered, and dropping it would shrink whichever denominator it
  # belongs to and flatter the screen.
  key[!is.finite(key)] <- Inf
  ord    <- order(key)
  key    <- key[ord]
  inBase <- inBase[ord]

  P <- sum(inBase)
  N <- length(inBase) - P
  if (P == 0 || N == 0) return(NULL)

  # Tied statistics move as a block: every gene sharing a value takes the counts from
  # the end of its block, so the curve crosses it diagonally instead of treating an
  # arbitrary sort order as evidence. This is the Mann-Whitney AUC with average ranks
  # -- a tie counts as half -- and it is what lets the gene-level table below be
  # integrated straight back into the AUC reported on the plot.
  runs <- rle(key)$lengths
  ends <- cumsum(runs)
  tpr  <- cumsum(inBase)[rep(ends, runs)] / P
  fpr  <- cumsum(!inBase)[rep(ends, runs)] / N

  list(
    order = ord, tpr = tpr, fpr = fpr, P = P, N = N,
    # (0, 0) belongs at the start of the curve. (1, 1) needs no anchor: the sweep
    # covers every gene in the screen, so it arrives there on its own.
    curve = tibble(fpr = c(0, fpr[ends]), tpr = c(0, tpr[ends]))
  )
}

#' Trapezoid over points already in sweep order.
rocAuc <- function(curve) {
  sum(diff(curve$fpr) * (head(curve$tpr, -1) + tail(curve$tpr, -1)) / 2)
}

plotRoc <- function(data, baseContrast, headContrasts, method, stat, cutoff, direction) {
  if (length(baseContrast) < 1 || !nzchar(baseContrast[1])) {
    return(msgResult("Please select a base screen."))
  }
  if (length(headContrasts) < 1) {
    return(msgResult("Please select at least one query screen."))
  }
  if (any(baseContrast %in% headContrasts)) {
    return(msgResult("The base screen must not also appear in the query screens."))
  }
  if (is.na(asNumericOr(cutoff, NA_real_))) return(msgResult("The cutoff must be a number."))

  # The statistic, the cutoff and the direction define the positives, and they do that
  # in the base screen alone. Applied to a query screen as well they would delete
  # positives and negatives alike from the sweep, and recall would reach 1 having
  # never looked at the genes that screen happened to score the other way.
  baseIds <- comparisonIdsFor(data, baseContrast)
  if (length(baseIds) == 0) return(msgResult("The base screen is not in the loaded data."))
  baseGenes <- unique(applyDirectionCutoff(
    setAnalysisFrame(data, baseIds, method, stat), stat, cutoff, direction
  )$gene)
  if (length(baseGenes) == 0) {
    return(msgResult("No hits in the base screen at this cutoff.\nPlease relax the cutoff."))
  }

  headIds <- comparisonIdsFor(data, headContrasts)
  if (length(headIds) == 0) {
    return(msgResult("None of the query screens are in the loaded data."))
  }
  headFrame <- fetchStat(data, headIds, NULL, method) %>%
    filter(gene_id != "X") %>%
    # One row per gene: a screen reporting a gene twice would count twice in both
    # denominators and once in the sweep order, for no reason the reader could see.
    distinct(comparison_id, gene_id, .keep_all = TRUE)
  if (nrow(headFrame) == 0) {
    return(msgResult("No data in any of the query screens for this method."))
  }
  headFrame <- headFrame %>%
    group_by(comparison_id) %>%
    # Ranked within the screen. Pooling the query screens would make one position in
    # the sweep mean a different thing in each of them.
    mutate(rank1 = rank(score, ties.method = "min"),
           rank2 = rank(-score, ties.method = "min")) %>%
    ungroup()
  headFrame$key <- rocSweepKey(headFrame, stat, direction)

  nBase  <- length(baseGenes)
  sweeps <- Filter(Negate(is.null), lapply(split(headFrame, headFrame$comparison_id), function(this) {
    sw <- rocSweep(this$key, this$gene_id %in% baseGenes)
    if (is.null(sw)) return(NULL)
    this   <- this[sw$order, ]
    screen <- friendlyIdsFor(data, this$comparison_id[1])[1]
    if (sw$P < nBase) {
      log_warn(paste0("ROC: ", nBase - sw$P, " of ", nBase, " base hits are absent from ",
                      screen, "; recall is measured against the ", sw$P, " it can reach."))
    }
    list(
      screen = screen, auc = rocAuc(sw$curve), P = sw$P,
      curve = sw$curve,
      genes = tibble(
        Contrast       = screen,
        `Gene symbol`  = this$gene_id,
        Score          = this$score,
        FDR            = this$fdr,
        `Hit in base`  = this$gene_id %in% baseGenes,
        # Rounded where the rest of the app rounds, at the point of reporting. Six
        # figures still rebuild the AUC to within a part in a million.
        Recall         = signif(sw$tpr, 6),
        FPR            = signif(sw$fpr, 6)
      )
    )
  }))
  if (length(sweeps) == 0) {
    return(msgResult(paste0("None of the query screens contain any of the ", nBase,
                            " base hits, so there is nothing to recover.")))
  }

  # How many of the base hits each screen could reach belongs on the plot: a curve
  # drawn over 40 of 60 positives is not the curve the legend would otherwise claim.
  legend <- vapply(sweeps, function(s) paste0(
    s$screen, "\nAUC: ", signif(s$auc, 4),
    if (s$P < nBase) paste0("\n", s$P, " of ", nBase, " hits present") else ""
  ), character(1))

  curve <- bind_rows(lapply(seq_along(sweeps), function(i) {
    sweeps[[i]]$curve %>% mutate(legend = legend[i])
  }))

  p <- ggplot(curve, aes(x = fpr, y = tpr, colour = legend)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_colour_discrete(name = "Contrast") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_classic() +
    xlab("False positive rate\n(1 - specificity)") +
    ylab("True positive rate\n(recall)") +
    ggtitle(paste0(
      "Recovery of hits in <b>", baseContrast[1], "</b>",
      "<br><sup>", nBase, " hit", if (nBase == 1) "" else "s", ": ", stat %||% "Score",
      " better than ", asNumericOr(cutoff, NA_real_), ", ",
      switch(direction %||% "",
             `Hypersensitivity hits` = "hypersensitivity",
             `Suppressing hits`      = "suppressing",
             "either direction"),
      ", ", method, "</sup>"
    ))

  list(
    p = ggplotly(p),
    # The gene-level sweep, not the curve's coordinates: this says which genes were
    # recovered and in what order, so the curve and the AUC can both be rebuilt from
    # it -- the trapezoid over (FPR, Recall) down these rows is the AUC on the plot --
    # and the hit set can be read off the flagged rows.
    plotdata = bind_rows(lapply(sweeps, function(s) s$genes)),
    height = 900
  )
}
