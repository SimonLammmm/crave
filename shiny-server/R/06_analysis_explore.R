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

#' Precision-recall sweep for one query screen against a reference hit set.
#'
#' Sorts once and uses cumulative counts, so the whole threshold range costs a
#' single pass over the query screen.
#'
#' @param ascending TRUE when smaller statistics are stronger (FDR, Rank), so the
#'   selection is `stat <= s`; FALSE for Score, where it is `stat >= s`.
rocSweep <- function(headStat, inBase, steps, ascending) {
  ok <- !is.na(headStat)
  headStat <- headStat[ok]
  inBase   <- inBase[ok]
  if (length(headStat) == 0) return(list(q = rep(0, length(steps)), tot = rep(0, length(steps))))

  ord <- order(headStat, decreasing = !ascending)
  sorted <- headStat[ord]
  cumIn  <- cumsum(inBase[ord])

  idx <- if (ascending) {
    findInterval(steps, sorted)
  } else {
    findInterval(-steps, -sorted)
  }
  q <- ifelse(idx > 0, cumIn[pmax(idx, 1L)], 0)
  list(q = q, tot = idx)
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
  if (!needPkg("DescTools")) {
    return(msgResult("The DescTools package is not installed, so AUC cannot be computed."))
  }

  cutoffValue <- asNumericOr(cutoff, NA_real_)
  if (is.na(cutoffValue)) return(msgResult("The cutoff must be a number."))

  #' Apply the direction filter and derive the sweep statistic.
  prepare <- function(df, applyCutoff) {
    if (identical(direction, "Hypersensitivity hits")) df <- df %>% filter(score < 0)
    else if (identical(direction, "Suppressing hits")) df <- df %>% filter(score > 0)
    if (identical(stat, "Score")) {
      df <- df %>% mutate(stat = abs(score))
      if (applyCutoff) df <- df %>% filter(abs(score) >= abs(cutoffValue))
    } else if (identical(stat, "FDR")) {
      df <- df %>% mutate(stat = fdr)
      if (applyCutoff) df <- df %>% filter(!is.na(fdr) & fdr <= cutoffValue)
    } else {
      df <- df %>% mutate(rank = rank(abs(score)), stat = rank)
      if (applyCutoff) df <- df %>% filter(rank <= cutoffValue)
    }
    df
  }

  baseHits <- prepare(
    fetchStat(data, comparisonIdsFor(data, baseContrast), NULL, method), TRUE
  )
  if (nrow(baseHits) < 1) {
    return(msgResult("No hits in the base screen at this cutoff.\nPlease relax the cutoff."))
  }

  headHits <- prepare(
    fetchStat(data, comparisonIdsFor(data, headContrasts), NULL, method), FALSE
  )
  if (nrow(headHits) < 1) {
    return(msgResult("No data in any of the query screens with these settings."))
  }

  ascending <- !identical(stat, "Score")
  baseGenes <- unique(baseHits$gene_id)
  nBase <- nrow(baseHits)

  headSplit <- split(headHits, headHits$comparison_id)
  roc <- bind_rows(lapply(names(headSplit), function(h) {
    this <- headSplit[[h]]
    maxStat <- suppressWarnings(max(this$stat, na.rm = TRUE))
    if (!is.finite(maxStat)) return(NULL)
    steps <- seq(from = 0, to = maxStat, length.out = CRAVE_ROC_STEPS)
    sweep <- rocSweep(this$stat, this$gene_id %in% baseGenes, steps, ascending)
    k <- nrow(this)
    precision <- 1 - ((sweep$tot - sweep$q) / k)
    recall    <- sweep$q / nBase
    # Normalise to rates.
    if (max(precision, na.rm = TRUE) != 0) precision <- precision / max(precision, na.rm = TRUE)
    if (max(recall,    na.rm = TRUE) != 0) recall    <- recall    / max(recall,    na.rm = TRUE)
    tibble(step = steps, precision = precision, recall = recall) %>%
      bind_rows(tibble(step = NA_real_, precision = c(0, 1), recall = c(1, 0))) %>%
      mutate(head = h, screen = friendlyIdsFor(data, h)[1])
  }))

  if (is.null(roc) || nrow(roc) == 0) {
    return(msgResult("The precision-recall curve could not be computed."))
  }

  auc <- roc %>%
    group_by(head) %>%
    summarise(auc = suppressWarnings(as.numeric(DescTools::AUC(1 - precision, recall))),
              .groups = "drop")

  roc <- roc %>%
    left_join(auc, by = "head") %>%
    mutate(legend = paste0(screen, "\nAUC: ", signif(auc, 4)))

  p <- ggplot(roc, aes(x = 1 - precision, y = recall, colour = legend, auc = auc)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0) +
    scale_colour_discrete(name = "Contrast") +
    theme_classic() +
    xlab("1 - Precision\n(False positive rate)") +
    ylab("Recall\n(True positive rate)") +
    ggtitle(paste0("Precision-recall for <b>", baseContrast[1], "</b>"))

  list(
    p = ggplotly(p),
    plotdata = roc %>%
      transmute(Contrast = screen, `ROC-AUC` = auc, Precision = precision, Recall = recall) %>%
      distinct(),
    height = 900
  )
}
