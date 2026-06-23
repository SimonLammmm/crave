#### Shiny server: heatmap function ####
plot_heatmap <- function(input) {
  # Apply screen filters
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
      mutate(`Gene A` = my_genes) %>%
      pivot_longer(-`Gene A`, names_to = "Gene B", values_to = "Correlation coefficient")
    # Dendrogram
    clustergram_ok <- length(unique(plotdata$`Gene A`)) >= 3 # Check if there is enough plotdata for a clustergram
    if(clustergram_ok) { # Create dendrogram on x if ok
      d <- plotdata %>%
        pivot_wider(names_from = `Gene B`, values_from = `Correlation coefficient`)
      d[is.na(d)] <- 0
      rownames <- d$`Gene A`
      d$`Gene A` <- NULL
      suppressWarnings({ rownames(d) <- rownames })
      d <- d %>% agnes(method = "ward")
      order <- d$order.lab
      d <- d %>% as.dendrogram() %>% ggdendrogram()
      dendro <- d$layers[[2]]$data # Manipulate dendrogram for overplotting heatmap
      dendroMaxX <- max(dendro$x, dendro$xend)
      dendroMinX <- min(dendro$x, dendro$xend)
      dendroMaxY <- max(dendro$y, dendro$yend)
      heatmapOffset <- length(unique(plotdata$`Gene B`))
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
      d$layers[[2]]$data <- dendroWithRoot
    }
    if (clustergram_ok) { # Factorise plotdata with clustergram order if exist
      plotdata$`Gene A` <- factor(plotdata$`Gene A`, ordered = T, levels = order) # Sort x axis by dendrogram order
      plotdata$`Gene B` <- factor(plotdata$`Gene B`, ordered = T, levels = order) # Sort y axis by dendrogram order
    } else { # If dendrogram not exist, sort both axes by alphabet order
      plotdata$`Gene A` <- factor(plotdata$`Gene A`, ordered = T, levels = sort(unique(plotdata$`Gene A`)))
      plotdata$`Gene B` <- factor(plotdata$`Gene B`, ordered = T, levels = sort(unique(plotdata$`Gene B`)))
    }
    # Plot
    p <- suppressWarnings( { # Plot heatmap
      ggplot(plotdata) +
        geom_tile(aes(x = `Gene A`, y = `Gene B`, fill = `Correlation coefficient`)) +
        scale_y_discrete(limits = rev(levels(plotdata$`Gene B`))) +
        theme_classic() +
        xlab("") +
        ylab("")
    } )
    if(is.null(input$correlate_customise_heatmap_low)) { mylow = defaults_correlate_customise_heatmap_low
    } else mylow <- input$correlate_customise_heatmap_low
    if(is.null(input$correlate_customise_heatmap_mid)) { mymid = defaults_correlate_customise_heatmap_mid
    } else mymid <- input$correlate_customise_heatmap_mid
    if(is.null(input$correlate_customise_heatmap_high)) { myhigh = defaults_correlate_customise_heatmap_high
    } else myhigh <- input$correlate_customise_heatmap_high
    if(max(plotdata$`Correlation coefficient`, na.rm = T) < 0 | min(plotdata$`Correlation coefficient`, na.rm = T) > 0) {
      p <- p + scale_fill_gradient2(low = mylow, mid = mymid, high = myhigh)
    } else {
      p <- p + scale_fill_gradientn(colours = c(mylow, mymid, myhigh), values = rescale(c(min(plotdata$`Correlation coefficient`, na.rm = T), 0, max(plotdata$`Correlation coefficient`, na.rm = T))))
    }
    if (clustergram_ok) { # Overplot dendrogram if exist; fix y-axis
      p <- p + d$layers[[2]] +
        theme(axis.line.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        annotate(x=0, xend=0, y=0, yend = heatmapOffset + 0.5, geom = "segment")
    }
    if (nrow(plotdata) == 0) { # Deal when plotdata is empty
      plotdata <- tibble(x = 0, y = 0, label = "Not enough genes to plot a heatmap.\nAdjust GOI.")
      p <- ggplot(plotdata, aes(x = x, y = y, label = label)) +
        geom_text(size = 6) +
        theme_classic() +
        theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) +
        xlab("") +
        ylab("")
      plotdata <- tibble()
    }
    pp <- ggplotly(p)  # Plotly
  } else if (length(finalScreens) > Inf) {
    plotdata <- tibble(x = 0, y = 0, label = "Too many screens found with the specified filters.\nConsider making filters more stringent.")
    pp <- ggplotly(ggplot(plotdata, aes(x = x, y = y, label = label)) + geom_text())
  } else { # Deal when no screens found
    plotdata <- tibble(x = 0, y = 0, label = "Not enough screens found with the specified filters.")
    pp <- ggplotly(ggplot(plotdata, aes(x = x, y = y, label = label)) + geom_text())
  }
  heatmap <- list(plot = pp, table = plotdata)
  return (heatmap)
}