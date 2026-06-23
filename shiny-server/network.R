
#### Shiny server: network function ####
plot_network <- function(input) {
  finalScreens <- gq_apply_filter(input)
  if (length(finalScreens) >= 2 & length(finalScreens) <= Inf) {
    # Fetch scores
    plotdata0 <- fetchDs(finalScreens, input$genes_queried, input$analysis_method_gq)
    plotdata <- plotdata0 %>%
      transmute(gene = gene_id, comparison = comparison_id, score = score) %>%
      pivot_wider(id_cols = gene, names_from = comparison, values_from = score)
    if(length(input$genes_queried) == 0) {
      fdrs <- plotdata0 %>%
        transmute(gene = gene_id, comparison = comparison_id, fdr = fdr, score = score) %>%
        group_by(comparison) %>%
        mutate(rank1 = rank(score), rank2 = rank(-score)) %>%
        ungroup() %>%
        filter(gene != "X") %>%
        filter(comparison %in% finalScreens)
      gois <- fdrs %>%
        filter(case_when(input$gq_stat == "Score" ~ abs(score) >= abs(input$gq_cutoff),
                         input$gq_stat == "FDR" ~ fdr <= input$gq_cutoff,
                         input$gq_stat == "Rank" ~ rank1 <= input$gq_cutoff | rank2 <= input$gq_cutoff)) %>%
        group_by(gene) %>% summarise(score = max(abs(score), na.rm = T))
      # Limit auto-detected GOIs to 100
      if(nrow(gois) > 250) {
        log_warn(paste0("Too many genes with cutoff better than ", input$gq_cutoff, ". Accepting the best 250 scores only."))
        gois <- gois %>% slice_max(order_by = score, n = 250)
      }
      gois <- gois %>% select(gene) %>% unlist() %>% unique()
      plotdata <- plotdata %>% filter(gene %in% gois)
    }
    my_genes <- plotdata$gene
    plotdata <- plotdata %>% select(-gene) %>% t()
    colnames(plotdata) <- my_genes
    # Correlate
    plotdata <- cor(plotdata, use = "pairwise.complete.obs", method = "spearman")
    plotdata <- plotdata %>%
      as_tibble() %>%
      setNames(my_genes) %>%
      mutate(from = my_genes) %>%
      pivot_longer(-from, names_to = "to", values_to = "corr")
    plotdata <- plotdata %>%
      filter(to != from) %>%
      filter(abs(corr) >= input$corr_cutoff)
    # Network
    if(nrow(plotdata > 0)) {
      set.seed(3141)
      nodes <- unique(c(plotdata$from, plotdata$to))
      g <- graph_from_data_frame(plotdata, directed = F, vertices = nodes)
      V(g)$color.background <- "#F2E95C"
        V(g)$color.border <- "#F0D337"
          V(g)$size <- as.double(15)
          V(g)$borderWidth <- as.double(3)
          p <- ggraph(g, layout = "nicely") + # Plot network
            geom_edge_link(aes(width = abs(corr)), colour = "#00007E") +
            geom_node_circle(aes(r = size / 120, colour = color.border, fill = color.background), size = 200 / length(nodes)) +
            geom_node_text(aes(label = name), size = 4) +
            coord_equal() +
            scale_colour_manual(breaks = c("#5B9D73", "#00007E", "#F0D337"), values = c("#5B9D73", "#00007E", "#F0D337")) +
            scale_fill_manual(breaks = c("#C64095", "#FFFFFF", "#F2E95C"), values = c("#C64095", "#FFFFFF", "#F2E95C")) +
            theme(panel.background = element_rect(fill = "transparent"), legend.position = "none")
          pp <- visIgraph(g) %>% # Plot interactive
            visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T))
          plotdata <- plotdata %>% transmute(`Gene A` = from, `Gene B` = to, `Correlation coefficient` = corr)
    } else {
      pp <- visNetwork(nodes = tibble(label = "No interactions found.\nAdjust screen filters or GOI."))
      plotdata <- tibble()
    }
  } else if (length(finalScreens) > Inf) {
    pp <- visNetwork(nodes = tibble(label = "Too many screens found with the specified filters.\nConsider making filters more stringent."))
    plotdata <- tibble()
  } else {
    pp <- visNetwork(nodes = tibble(label = "Not enough screens found with the specified filters."))
    plotdata <- tibble()
  }
  network <- list(plot = pp, table = plotdata)
  return(network)
}