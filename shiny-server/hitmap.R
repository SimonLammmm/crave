#### Shiny server: gene query function ####
hitmap <- function(input) {
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
      # plotdata
      plotdata <- fetchdata %>%
        mutate(goi = gene %in% myGoi,
               #goi = case_when(grepl(paste0("^", myGoi, "$", collapse = "|"), gene, ignore.case = T) ~ T, T ~ F),
               sig = case_when(input$gq_stat == "Score" & abs(score) >= input$gq_cutoff ~ "*",
                               input$gq_stat == "FDR" & fdr <= input$gq_cutoff ~ "*",
                               input$gq_stat == "Rank" & (rank1 <= input$gq_cutoff | rank2 <= input$gq_cutoff) ~ "*",
                               T ~ "")) %>%
        filter(goi)
      # Dendrogram
      clustergram_ok <- length(unique(plotdata$gene)) >= 3 & input$hitmap_gene_dendrogram # Check if there is enough plotdata for a clustergram and the user wants one
      if(clustergram_ok) { # Create dendrogram on x if ok
        d <- plotdata %>%
          pivot_wider(id_cols = gene, names_from = FriendlyID, values_from = score)
        d[is.na(d)] <- 0
        rownames <- d$gene
        d$gene <- NULL
        suppressWarnings({ rownames(d) <- rownames })
        d <- d %>% agnes(method = "ward")
        order <- d$order.lab
        d <- d %>% as.dendrogram() %>% ggdendrogram(rotate = T)
        dendro <- d$layers[[2]]$data # Manipulate dendrogram for overplotting heatmap
        dendroMaxX <- max(dendro$x, dendro$xend)
        dendroMinX <- min(dendro$x, dendro$xend)
        dendroMaxY <- max(dendro$y, dendro$yend)
        heatmapOffset <- length(unique(plotdata$FriendlyID))
        rootLength <- 2
        dendroHeightScale <- 0.1
        topbranch <- dendro %>% filter(y == yend & y == max(y))
        topbranchmidpoint <- mean(c(topbranch$x, topbranch$xend))
        root <- data.frame(x = topbranchmidpoint,
                           xend = topbranchmidpoint,
                           y = dendroMaxY,
                           yend = dendroMaxY + rootLength)
        dendroWithRoot <- rbind(dendro, root)
        dendroMaxY <- max(dendroWithRoot$y, dendroWithRoot$yend)
        dendroWithRoot$y <- dendroWithRoot$y / dendroMaxY * heatmapOffset * dendroHeightScale + heatmapOffset + 0.5
        dendroWithRoot$yend <- dendroWithRoot$yend / dendroMaxY * heatmapOffset * dendroHeightScale + heatmapOffset + 0.5
        dendroWithRoot <- data.frame(x = dendroWithRoot$y, y = dendroWithRoot$x, xend = dendroWithRoot$yend, yend = dendroWithRoot$xend)
        d$layers[[2]]$data <- dendroWithRoot
        plotdata$gene <- factor(plotdata$gene, ordered = T, levels = order) # Sort x axis by dendrogram order
      } else {
        plotdata$gene <- factor(plotdata$gene, levels = sort(unique(plotdata$gene), decreasing = T)) # If no dendrogram, order the genes axis alphabetically
      }
      # Dendrogram
      clustergram2_ok <- length(unique(plotdata$FriendlyID)) >= 3 & input$hitmap_screen_dendrogram # Check if there is enough plotdata for a clustergram and the user wants one
      if(clustergram2_ok) { # Create dendrogram on x if ok
        d2 <- plotdata %>%
          pivot_wider(id_cols = FriendlyID, names_from = gene, values_from = score)
        d2[is.na(d2)] <- 0
        rownames <- d2$FriendlyID
        d2$FriendlyID <- NULL
        suppressWarnings({ rownames(d2) <- rownames })
        d2 <- d2 %>% agnes(method = "ward")
        order2 <- d2$order.lab
        d2 <- d2 %>% as.dendrogram() %>% ggdendrogram()
        dendro <- d2$layers[[2]]$data # Manipulate dendrogram for overplotting heatmap
        dendroMaxX <- max(dendro$x, dendro$xend)
        dendroMinX <- min(dendro$x, dendro$xend)
        dendroMaxY <- max(dendro$y, dendro$yend)
        heatmapOffset <- length(unique(plotdata$gene))
        rootLength <- 2
        dendroHeightScale <- 0.1
        topbranch <- dendro %>% filter(y == yend & y == max(y))
        topbranchmidpoint <- mean(c(topbranch$x, topbranch$xend))
        root <- data.frame(x = topbranchmidpoint,
                           xend = topbranchmidpoint,
                           y = dendroMaxY,
                           yend = dendroMaxY + rootLength)
        dendroWithRoot <- rbind(dendro, root)
        dendroMaxY <- max(dendroWithRoot$y, dendroWithRoot$yend)
        dendroWithRoot$y <- dendroWithRoot$y / dendroMaxY * heatmapOffset * dendroHeightScale + heatmapOffset + 0.5
        dendroWithRoot$yend <- dendroWithRoot$yend / dendroMaxY * heatmapOffset * dendroHeightScale + heatmapOffset + 0.5
        d2$layers[[2]]$data <- dendroWithRoot
        plotdata$FriendlyID <- factor(plotdata$FriendlyID, ordered = T, levels = order2) # Sort y axis by dendrogram order
      } else {
        plotdata$FriendlyID <- factor(plotdata$FriendlyID, ordered = T, levels = unique(input$gq_filter_contrast))
      }
      # plot
      p <- ggplot() +
        geom_tile(data = plotdata, mapping = aes(y = gene, x = `FriendlyID`, fdr = fdr, fill = score, rank_synth = rank1, rank_supp = rank2, label = sig, DaysDiff = `Days grown (diff)`, TreatmentDiff = `Treatment (diff)`, DoseDiff = `Dose (diff)`, KnockoutDiff = `Knockout (diff)`, CellLineDiff = `Cell line (diff)`, DaysRef = `Days grown (ref)`, TreatmentRef = `Treatment (ref)`, DoseRef = `Dose (ref)`, KnockoutRef = `Knockout (ref)`, CellLineRef = `Cell line (ref)`, Library = `Library`)) +
        geom_text(data = plotdata, mapping = aes(y = gene, x = `FriendlyID`, fdr = fdr, fill = score, rank_synth = rank1, rank_supp = rank2, label = sig, DaysDiff = `Days grown (diff)`, TreatmentDiff = `Treatment (diff)`, DoseDiff = `Dose (diff)`, KnockoutDiff = `Knockout (diff)`, CellLineDiff = `Cell line (diff)`, DaysRef = `Days grown (ref)`, TreatmentRef = `Treatment (ref)`, DoseRef = `Dose (ref)`, KnockoutRef = `Knockout (ref)`, CellLineRef = `Cell line (ref)`, Library = `Library`), nudge_y = -0.15) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("") +
        ylab("")
      if(is.null(input$correlate_customise_heatmap_low)) { mylow = defaults_correlate_customise_heatmap_low
      } else mylow <- input$correlate_customise_heatmap_low
      if(is.null(input$correlate_customise_heatmap_mid)) { mymid = defaults_correlate_customise_heatmap_mid
      } else mymid <- input$correlate_customise_heatmap_mid
      if(is.null(input$correlate_customise_heatmap_high)) { myhigh = defaults_correlate_customise_heatmap_high
      } else myhigh <- input$correlate_customise_heatmap_high
      if(max(plotdata$score, na.rm = T) < 0 | min(plotdata$score, na.rm = T) > 0) {
        p <- p + scale_fill_gradient2(low = mylow, mid = mymid, high = myhigh)
      } else {
        p <- p + scale_fill_gradientn(colours = c(mylow, mymid, myhigh), values = rescale(c(min(plotdata$score, na.rm = T), 0, max(plotdata$score, na.rm = T))))
      }
      if (clustergram_ok) { # Overplot dendrogram if exist
        p <- p + d$layers[[2]]
      }
      if (clustergram2_ok) { # Overplot dendrogram if exist
        p <- p + d2$layers[[2]]
      }
      displaydata <- plotdata %>%
        filter(goi) %>%
        select(-`Comparison ID`, -`Experiment ID`, -Contrast) %>%
        relocate(Gene = gene, Contrast = FriendlyID, FDR = fdr, Score = score)
      if(is.null(input$correlate_customise_height)) { height <- defaults_correlate_customise_height
      } else height <- input$correlate_customise_height
      if(is.null(input$correlate_customise_width_auto)) { width <- NULL
      } else if (input$correlate_customise_width_auto) { width <- NULL
      } else width <- input$correlate_customise_width
      p <- ggplotly(p, height = height, width = width)
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
  hitmap <- list(displaydata = displaydata, p = p)
  return(hitmap)
}