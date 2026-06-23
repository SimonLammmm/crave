#### Shiny server: bi-plot function ####
plotBiplot <- function(input) {
  if (input$comparison_X != "" & input$comparison_Y != "" & input$analysis_method_se_b_x != "" & input$analysis_method_se_b_y != "" & !identical(c(input$comparison_X, input$analysis_method_se_b_x), c(input$comparison_Y, input$analysis_method_se_b_y))) {
    # Plot a biplot
    comparison_x <- comparisons$`Comparison ID`[input$comparison_X == comparisons$FriendlyID]
    comparison_y <- comparisons$`Comparison ID`[input$comparison_Y == comparisons$FriendlyID]
    plotdataX <- fetchDs(comparison_x, NULL, input$analysis_method_se_b_x) %>%
      transmute(gene = gene_id, fdrX = fdr, scoreX = score)
    plotdataY <- fetchDs(comparison_y, NULL, input$analysis_method_se_b_y) %>%
      transmute(gene = gene_id, fdrY = fdr, scoreY = score)
    plotdata <- inner_join(plotdataX, plotdataY, by = "gene")
    plotdata <- plotdata %>%
      filter(!is.na(fdrX)) %>%
      filter(!is.na(scoreX)) %>%
      filter(!is.na(fdrY)) %>%
      filter(!is.na(scoreY)) %>%
      filter(gene != "X") %>%
      filter(!grepl("Non-targeting", gene)) %>%
      mutate(goi = case_when(gene %in% input$genes_selected ~ "Yes", T ~ "No"),
             hovertext = paste0(gene, "\n", input$analysis_method_se_b_x, " score (x-axis): ", scoreX, "\n", input$analysis_method_se_b_x, " FDR (x-axis): ", fdrX,
                                "\n", input$analysis_method_se_b_y, " score (y-axis): ", scoreY, "\n", input$analysis_method_se_b_y, " FDR (y-axis): ", fdrY, "\nGOI?: ", goi))
    # apply colours
    colors = c(defaults_explore_customise_point_colour, defaults_explore_customise_goi_colour)
    if (!is.null(input$explore_customise_point_colour)) colors[1] = input$explore_customise_point_colour
    if (!is.null(input$explore_customise_goi_colour)) colors[2] = input$explore_customise_goi_colour
    # apply x axis
    scale_x = defaults_explore_customise_x_axis_log
    if (!is.null(input$explore_customise_x_axis_log)) scale_x = input$explore_customise_x_axis_log
    scale_x = case_when(scale_x == "Automatic" ~ "linear",
                        scale_x == "Logarithmic" ~ "linear", # not allow log x-axis
                        scale_x == "Linear" ~ "linear",
                        T ~ "linear")
    # apply y axis
    scale_y = defaults_explore_customise_y_axis_log
    if (!is.null(input$explore_customise_y_axis_log)) scale_y = input$explore_customise_y_axis_log
    scale_y = case_when(scale_y == "Automatic" ~ "linear",
                        scale_y == "Logarithmic" ~ "linear", # not allow log x-axis
                        scale_y == "Linear" ~ "linear",
                        T ~ "linear")
    # apply dimensions
    width = 0
    if (!is.null(input$explore_customise_width)) width = input$explore_customise_width
    if (width == 0) width = NULL # zero width means default
    height = 0
    if (!is.null(input$explore_customise_height)) height = input$explore_customise_height
    if (height == 0) height = defaults_explore_customise_height # zero height means default
    # apply square plot
    square = defaults_explore_customise_y_equals_x
    if (!is.null(input$explore_customise_y_equals_x)) square = input$explore_customise_y_equals_x
    if (square) {
      # square plot
      limits = c(min(plotdata$scoreX, plotdata$scoreY), max(plotdata$scoreX, plotdata$scoreY))
    } else {
      # scaled to each axis
      limits = c(max(min(plotdata$scoreX), min(plotdata$scoreY)), min(max(plotdata$scoreX), max(plotdata$scoreY)))
    }
    # plotly
    p <- plot_ly(plotdata, x = ~scoreX, y = ~scoreY, hovertext = ~hovertext, color = ~goi, colors = colors[1:length(unique(plotdata$goi))], type = "scattergl", mode = "markers", width = width, height = height) %>%
      layout(shapes = list(list(type = "line", x0 = ~limits[1], x1 = ~limits[2], xref = "x", y0 = ~limits[1], y1 = ~limits[2], yref = "y", line = list(color = "black"))),
             title = paste0("<i>x</i>: ",
                            contrastFormatter(input$comparison_X),
                            "\n<i>y</i>: ",
                            contrastFormatter(input$comparison_Y)),
             xaxis = list(title = paste0(input$analysis_method_se_b_x, " score: ", input$comparison_X),
                          type = scale_x),
             yaxis = list(title = paste0(input$analysis_method_se_b_y, " score: ", input$comparison_Y),
                          type = scale_y),
             legend = list(title = list(text = '<b>GOI?</b>'),
                           x = 0.1,
                           y = 0.9),
             dragmode = "select")

    # Export data
    plotdata <- plotdata %>% transmute(gene, scoreX, fdrX, scoreY, fdrY,
                                       delta = round(abs(scoreX - scoreY),3),
                                       rmsd = round(sqrt(abs(scoreX ^ 2 - scoreY ^ 2)),3),
                                       goi) %>%
      arrange(desc(delta))
    names(plotdata) <- c("Gene", paste0(input$analysis_method_se_b_x, " score: ", input$comparison_X), paste0(input$analysis_method_se_b_x, " FDR: ", input$comparison_X),
                         paste0(input$analysis_method_se_b_y, " score: ", input$comparison_Y), paste0(input$analysis_method_se_b_y, " FDR: ", input$comparison_Y),
                         "Delta(score)", "RMSD(score)",
                         "GOI?")
  } else {
    plotdata = tibble(x = 0, y = 0, label = "Select two screens.", Gene = NA)
    p <- plot_ly(data = plotdata, x = ~x, y = ~y, text = ~label, type = "scattergl", mode = "text") %>%
      layout(plot_bgcolor='#FFFFFF',
             xaxis = list(title='', zerolinecolor = '#FFFFFF', zerolinewidth = 0, gridcolor = '#FFFFFF', tickvals = list()),
             yaxis = list(title='', zerolinecolor = '#FFFFFF', zerolinewidth = 0, gridcolor = '#FFFFFF', tickvals = list()))
    height = 100
  }
  biplot <- list(plotdata = plotdata, p = p, height = height)
  return(biplot)
}