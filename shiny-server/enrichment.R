

#### Shiny server: enrichment analysis function ####
plotEnrichment <- function(input) {
  if (input$analysis_method_gq != "" & grepl("^-*\\d+\\.*\\d*$", input$enrichment_cutoff_value) ) {
    # Apply screen filters
    finalScreens <- gq_apply_filter(input)
    # Enrichment analysis
    enrichment_cutoff_value <- as.numeric(input$enrichment_cutoff_value)
    alldata <- fetchDs(finalScreens, NULL, input$analysis_method_gq) %>%
      left_join(comparisons, by = c("comparison_id" = "Comparison ID")) %>%
      mutate(screen = FriendlyID) %>%
      group_by(screen) %>%
      transmute(gene = gene_id,
                sign = sign(score),
                rank = rank(score, ties.method = "random"),
                value = case_when(input$enrichment_cutoff_stat == "Score" ~ score,
                                  input$enrichment_cutoff_stat == "FDR" ~ fdr,
                                  input$enrichment_cutoff_stat == "Rank" ~ abs(rank)-n()/2)) %>%
      filter(gene != "X")
    
    if(input$enrichment_cutoff_direction == "Hypersensitivity hits" & input$enrichment_cutoff_stat == "Rank") cutoffdata <- alldata %>% filter(rank < max(rank)/2)
    if(input$enrichment_cutoff_direction == "Suppressing hits" & input$enrichment_cutoff_stat == "Rank") cutoffdata <- alldata %>% filter(rank > max(rank)/2)
    if(input$enrichment_cutoff_direction == "Hypersensitivity hits" & input$enrichment_cutoff_stat != "Rank") cutoffdata <- alldata %>% filter(sign == -1)
    if(input$enrichment_cutoff_direction == "Suppressing hits" & input$enrichment_cutoff_stat != "Rank") cutoffdata <- alldata %>% filter(sign == 1)
    
    if(input$enrichment_cutoff_stat == "Score") cutoffdata <- cutoffdata %>% filter(abs(value) >= abs(enrichment_cutoff_value))
    if(input$enrichment_cutoff_stat == "FDR") cutoffdata <- cutoffdata %>% filter(value <= enrichment_cutoff_value)
    if(input$enrichment_cutoff_stat == "Rank") cutoffdata <- cutoffdata %>% slice_max(abs(value), n = round(enrichment_cutoff_value, digits = 0))
    
    alldata <- alldata %>%
      left_join(ontology, by = c("gene" = "symbol"), relationship = "many-to-many") %>%
      filter(!is.na(class))
    
    cutoffdata <- cutoffdata %>%
      left_join(ontology, by = c("gene" = "symbol"), relationship = "many-to-many") %>%
      filter(!is.na(class))
    
    if(nrow(cutoffdata) > 0) {
      # Hypergeometric tests
      hypers <- foreach(screen = unique(alldata$screen), .combine = "bind_rows") %do% {
        foreach(class = unique(ontology$class), .combine = "bind_rows") %do% {
        tibble(screen = screen,
               class = class,
               x = length(unique(cutoffdata$gene[cutoffdata$class == class & cutoffdata$screen == screen])),
               m = length(unique(alldata$gene[alldata$class == class & alldata$screen == screen])),
               n = length(unique(alldata$gene[alldata$screen == screen])) - m,
               k = length(unique(cutoffdata$gene[cutoffdata$screen == screen])),
               phyper = phyper(x, m, n, k, input$enrichment_cutoff_tail == "Lower tail (depleted classes)"))
        }
      }
      hypers = hypers %>% 
        mutate(adj.phyper = signif(p.adjust(hypers$phyper, method = "BH")), 4,
               phyper = signif(phyper, 4))
      
      # Plot
      plotdata <- cutoffdata %>%
        group_by(screen, class) %>%
        summarise(genes = gsub("(.+?;.+?;.+?;.+?;.+?;.+?);(.+?)", "\\1\n\\2", paste0(sort(gene), collapse = "; "))) %>%
        left_join(hypers, by = c("screen", "class")) %>%
        mutate(recall = x / m,
               sig = case_when(adj.phyper < 0.05 ~ "*",
                               T ~ "")) %>%
        ungroup()
      
      p <- ggplot(plotdata, aes(x = screen, y = class, fill = recall, label = sig, symbols = genes,
                                phyperx = x, phyperk = k, phyperm = m, phypern = n, phyper = phyper, adj.phyper = adj.phyper)) +
        geom_tile() +
        geom_text(nudge_y = -0.3) +
        scale_fill_gradient(low = "#FFFFFF", high = "#FF0087", breaks = c(0, 1)) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      p <- ggplotly(p)
      
      plotdata <- plotdata %>% transmute(Contrast = screen, Class = class, Genes = genes, `Hypergeometric p-value` = phyper, `Hypergeometric FDR` = adj.phyper)
    } else {
      # Fail if no significant hits
      plotdata <- tibble(x = 0, y = 0, label = "No significant hits with given cutoff.\nConsider relaxing the cutoff.")
      p <- ggplotly(ggplot(plotdata, aes(x = x, y = y, label = label)) + geom_text())
    }
  } else {
    # Fail if invalid inputs
    plotdata <- tibble(x = 0, y = 0, label = "Error with input form.")
    p <- ggplotly(ggplot(plotdata, aes(x = x, y = y, label = label)) + geom_text())
  }
  enrichment <- list(plotdata = plotdata, p = p)
  return(enrichment)
}