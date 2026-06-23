
#### Shiny server: volcano plot function ####
plotVolcano <- function(input) {
  if (input$comparison_volcano != "" & input$analysis_method_se_v != "") {
    # Plot a volcano plot
    comparison_volcano <- comparisons$`Comparison ID`[input$comparison_volcano == comparisons$FriendlyID]
    plotdata <- fetchDs(comparison_volcano, NULL, input$analysis_method_se_v) %>%
      transmute(gene = gene_id, fdr = fdr, score = score)
    plotdata <- plotdata %>%
      filter(!is.na(fdr)) %>%
      filter(!is.na(score)) %>%
      filter(gene != "X") %>%
      mutate(goi = case_when(gene %in% input$genes_selected ~ "Yes", T ~ "No"),
             hovertext = paste0(gene, "\nScore: ", score, "\nFDR: ", fdr, "\nGOI?: ", goi)) %>%
      arrange(score) %>%
      arrange(desc(goi))
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
    scale_y = case_when(scale_y == "Automatic" ~ "log",
                        scale_y == "Logarithmic" ~ "log",
                        scale_y == "Linear" ~ "linear",
                        T ~ "log")
    # apply dimensions
    width = 0
    if (!is.null(input$explore_customise_width)) width = input$explore_customise_width
    if (width == 0) width = NULL # zero width means default
    height = 0
    if (!is.null(input$explore_customise_height)) height = input$explore_customise_height
    if (height == 0) height = defaults_explore_customise_height # zero height means default
    # plotly
    p <- plot_ly(plotdata, x = ~score, y = ~fdr, hovertext = ~hovertext, color = ~goi, colors = colors[1:length(unique(plotdata$goi))], type = "scattergl", mode = "markers", width = width, height = height) %>%
      layout(title = contrastFormatter(input$comparison_volcano),
             xaxis = list(title = paste0(input$analysis_method_se_v, " score"),
                          type = scale_x),
             yaxis = list(title = paste0(input$analysis_method_se_v, " FDR"),
                          type = scale_y, autorange = "reversed"),
             legend = list(title = list(text = '<b>GOI?</b>'),
                           x = 0.1,
                           y = 0.9),
             dragmode = "select")
    # Export data
    plotdata <- plotdata %>% transmute(Gene = gene, Score = score, FDR = fdr, `GOI?` = goi)
  } else {
    plotdata = tibble(x = 0, y = 0, label = "Select a screen.", Gene = NA, Score = NA, FDR = NA)
    p <- plot_ly(data = plotdata, x = ~x, y = ~y, text = ~label, type = "scattergl", mode = "text") %>%
      layout(plot_bgcolor='#FFFFFF',
             xaxis = list(title='', zerolinecolor = '#FFFFFF', zerolinewidth = 0, gridcolor = '#FFFFFF', tickvals = list()),
             yaxis = list(title='', zerolinecolor = '#FFFFFF', zerolinewidth = 0, gridcolor = '#FFFFFF', tickvals = list()))
    height = 100
  }
  volcano <- list(plotdata = plotdata, p = p, height = height)
  return(volcano)
}