#### Shiny server: gene query function ####
genequery <- function(input) {
  # Apply screen filters
  finalScreens <- gq_apply_filter(input)
  if (length(finalScreens) > 0 & length(finalScreens) <= Inf) {
    # Obtain screens of interest from GOI below cutoff
    include <- fetchDs(finalScreens, input$genes_queried, input$analysis_method_gq) %>%
      transmute(gene = gene_id, comparison = comparison_id, fdr = fdr, score = score) %>%
      group_by(comparison) %>%
      mutate(rank1 = rank(score), rank2 = rank(-score)) %>%
      ungroup() %>%
      filter(case_when(input$gq_stat == "Score" ~ abs(score) >= abs(input$gq_cutoff),
                       input$gq_stat == "FDR" ~ fdr <= input$gq_cutoff,
                       input$gq_stat == "Rank" ~ rank1 <= input$gq_cutoff | rank2 <= input$gq_cutoff)) %>%
      select(comparison) %>% unique() %>% unlist() %>% setNames(NULL)
    if (length(include) > 0 & length(include) <= 120) {
      # For screens of interest, obtain score and fdr
      fetchdata <- fetchDs(include, NULL, input$analysis_method_gq) %>%
        transmute(gene = gene_id, `Comparison ID` = comparison_id, fdr = fdr, score = score) %>%
        group_by(`Comparison ID`) %>%
        mutate(rank1 = rank(score), rank2 = rank(-score)) %>%
        ungroup() %>%
        filter(!is.na(fdr)) %>%
        filter(!is.na(score)) %>%
        filter(!is.na(rank1)) %>%
        filter(!is.na(rank2)) %>%
        filter(gene != "X") %>%
        left_join(comparisons, by = "Comparison ID")
      if(length(input$genes_queried) > 0) {
        myGoi <- input$genes_queried
      } else {
        myGoi <- fetchdata %>%
          filter(case_when(input$gq_stat == "Score" ~ abs(score) >= abs(input$gq_cutoff),
                           input$gq_stat == "FDR" ~ fdr <= input$gq_cutoff,
                           input$gq_stat == "Rank" ~ rank1 <= input$gq_cutoff | rank2 <= input$gq_cutoff)) %>%
          group_by(gene) %>% summarise(score = max(abs(score), na.rm = T))
        # Limit auto-detected GOIs
        if(nrow(myGoi) > 250) {
          log_warn(paste0("Too many genes with cutoff better than ", input$gq_cutoff, ". Accepting the top 250 scores only."))
          myGoi <- myGoi %>% slice_max(order_by = score, n = 250)
        }
        myGoi <- myGoi %>% select(gene) %>% unlist() %>% unique()
      }
      plotdata <- fetchdata %>%
        mutate(goi = gene %in% myGoi,
               sig = case_when(input$gq_stat == "Score" & abs(score) >= input$gq_cutoff ~ T,
                               input$gq_stat == "FDR" & fdr <= input$gq_cutoff ~ T,
                               input$gq_stat == "Rank" & (rank1 <= input$gq_cutoff | rank2 <= input$gq_cutoff) ~ T,
                               T ~ F))
      p <- ggplot(plotdata, aes(y = score, x = `FriendlyID`, fdr = fdr, score = score, rank_synth = rank1, rank_supp = rank2, DaysDiff = `Days grown (diff)`, TreatmentDiff = `Treatment (diff)`, DoseDiff = `Dose (diff)`, KnockoutDiff = `Knockout (diff)`, CellLineDiff = `Cell line (diff)`, DaysRef = `Days grown (ref)`, TreatmentRef = `Treatment (ref)`, DoseRef = `Dose (ref)`, KnockoutRef = `Knockout (ref)`, CellLineRef = `Cell line (ref)`, Library = `Library`))
      if(is.null(input$correlate_customise_genequery_fillscheme)) {
        fillscheme <- defaults_correlate_customise_genequery_fillscheme
      } else fillscheme <- input$correlate_customise_genequery_fillscheme
      # plot type settings
      plottype <- input$genequery_violin
      if(plottype == "Violin plot") my_geom <- geom_violin
      if(plottype == "Boxplot") my_geom <- geom_boxplot
      # fill settings
      if(fillscheme == "Treatment") p <- p + my_geom(na.rm = T, mapping = aes(colour = `Treatment (diff)`))
      if(fillscheme == "Knockout") p <- p + my_geom(na.rm = T, mapping = aes(colour = `Knockout (diff)`))
      if(fillscheme == "Cell line") p <- p + my_geom(na.rm = T, mapping = aes(colour = `Cell line`))
      if(fillscheme == "Library") p <- p + my_geom(na.rm = T, mapping = aes(colour = `Library`))
      p <- p + geom_point(data = plotdata %>% filter(goi), aes(colour = gene, shape = sig)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("") +
        ylab("Score")
      displaydata <- plotdata %>%
        filter(goi) %>%
        select(-`Comparison ID`, -`Experiment ID`, -Contrast) %>%
        relocate(Gene = gene, Contrast = FriendlyID, FDR = fdr, Score = score)
      if(is.null(input$correlate_customise_height)) { height <- defaults_correlate_customise_height
      } else height <- input$correlate_customise_height
      if(is.null(input$correlate_customise_width_auto)) { width <- NULL
      } else if (input$correlate_customise_width_auto) { width <- NULL
      } else width <- input$correlate_customise_width
      p <- ggplotly(p, height = height, width = width) %>% layout(dragmode = "select")
    } else if (length(include) > 120) {
      displaydata <- tibble(x = 0, y = 0, label = "Too many screens to display.\nConsider making filters and/or cutoffs more stringent.")
      p <- ggplotly(ggplot(displaydata, aes(x = x, y = y, label = label)) + geom_text())
    } else {
      displaydata <- tibble(x = 0, y = 0, label = "No hits found among screens with the specified filters and cutoffs.\nConsider relaxing filters and/or cutoffs.")
      p <- ggplotly(ggplot(displaydata, aes(x = x, y = y, label = label)) + geom_text())
    }
  } else if (length(finalScreens) > Inf) {
    displaydata <- tibble(x = 0, y = 0, label = "Too many screens found with the specified filters.\nConsider making filters more stringent.")
    p <- ggplotly(ggplot(displaydata, aes(x = x, y = y, label = label)) + geom_text())
  } else {
    displaydata <- tibble(x = 0, y = 0, label = "No screens found with the specified filters.\nConsider relaxing filters.")
    p <- ggplotly(ggplot(displaydata, aes(x = x, y = y, label = label)) + geom_text())
  }
  genequery <- list(displaydata = displaydata, p = p)
  return(genequery)
}