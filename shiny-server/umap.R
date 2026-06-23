
#### Shiny server: UMAP function ####
plot_umap <- function(input) {
  finalScreens <- gq_apply_filter(input)
  if (length(finalScreens) >= 4 & length(finalScreens) <= Inf) {
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
      if(nrow(gois) > 100) {
        log_warn(paste0("Too many genes with cutoff better than ", input$gq_cutoff, ". Accepting the best 250 scores only."))
        gois <- gois %>% slice_max(order_by = score, n = 100)
      }
      gois <- gois %>% select(gene) %>% unlist() %>% unique()
      plotdata <- plotdata %>% filter(gene %in% gois)
    }
    plotdata <- plotdata %>% filter(gene != "X") 
    my_genes <- gsub("-", "temporarydash", plotdata$gene)
    plotdata <- plotdata %>% select(-gene)
    rownames(plotdata) <- my_genes
    # Impute
    my_screens <- names(plotdata)
    names(plotdata) <- 1:length(plotdata)
    umap_ok <- F
    try({
      plotdata <- t(complete(mice(t(plotdata), method = "lasso.norm")))
      names(plotdata) <- my_screens
      custom.config <- umap.defaults
      custom.config$n_neighbors=3
      custom.config$min_dist=0.3
      custom.config$random_state = 101079
      plotdata <- umap(plotdata, config = custom.config)$layout
      plotdata <- plotdata %>%
        as_tibble() %>%
        transmute(Gene = gsub("temporarydash", "-", my_genes), UMAP1 = V1, UMAP2 = V2)
      if(nrow(ontology) == 0) {
        plotdata$Class <- NA
      } else {
        plotdata <- plotdata %>%
          left_join(ontology, by = c("Gene" = "symbol")) %>%
          transmute(Gene, UMAP1, UMAP2, Class = class)
      }
      umap_ok <- T
    })
    if(umap_ok) {
      # Plot
      p <- ggplot(plotdata, aes(x = UMAP1, y = UMAP2, Gene = Gene, colour = Class)) +
        geom_point()
      pp <- ggplotly(p) %>% layout(dragmode = "select")
      # Handle when umap failed
    } else {
      plotdata <- tibble(x = 0, y = 0, label = "Couldn't do it - server error.\nPlease try again.")
      pp <- ggplotly(ggplot(plotdata, aes(x = x, y = y, label = label)) + geom_text())
    }
  } else if (length(finalScreens) > Inf) {
    plotdata <- tibble(x = 0, y = 0, label = "Too many screens found with the specified filters.\nConsider making filters more stringent.")
    pp <- ggplotly(ggplot(plotdata, aes(x = x, y = y, label = label)) + geom_text())
  } else {
    plotdata <- tibble(x = 0, y = 0, label = "Not enough screens found with the specified filters.")
    pp <- ggplotly(ggplot(plotdata, aes(x = x, y = y, label = label)) + geom_text())
  }
  umap <- list(plot = pp, table = plotdata)
  return(umap)
}