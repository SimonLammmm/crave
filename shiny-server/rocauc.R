#### Shiny server: ROC-AUC ####

plotRoc <- function(input) {
  # check for valid inputs: at least one base screen, at least one head screen, base screen not equal to head screen, at least one hit in base screen at the cutoff specified
  if(length(input$rocauc_base) < 1) {
    # fail if less than 1 base screen
    table = tibble(x = 0, y = 0, label = "Please select a base screen.")
    p = ggplotly(ggplot(table, aes(x = x, y = y, label = label)) + geom_label() + theme_classic())
    return(list(plotdata = table, p = p))
  } else if(length(input$rocauc_head) < 1) {
    # fail if less than 1 head screen
    table = tibble(x = 0, y = 0, label = "Please select a head screen.")
    p = ggplotly(ggplot(table, aes(x = x, y = y, label = label)) + geom_label() + theme_classic())
    return(list(plotdata = table, p = p))
  } else if(input$rocauc_base %in% input$rocauc_head) {
    # fail if base screen exists in head screens
    table = tibble(x = 0, y = 0, label = "Please make sure the base screen does not appear in the head screens.")
    p = ggplotly(ggplot(table, aes(x = x, y = y, label = label)) + geom_label() + theme_classic())
    return(list(plotdata = table, p = p))
  } else {
    # obtain hits in base screen at the cutoff specified
    # obtain hits
    baseHits = fetchDs(comparison = comparisons$`Comparison ID`[comparisons$FriendlyID %in% input$rocauc_base], gene = NULL, method = input$rocauc_method)
    # apply direction specified
    if (input$rocauc_cutoff_direction == "Hypersensitivity hits") {
      baseHits = baseHits %>%
        filter(score < 0)
    } else if (input$rocauc_cutoff_direction == "Suppressing hits") {
      baseHits = baseHits %>%
        filter(score > 0)
    }
    # filter by stat and cutoff specified
    # define stat to step along - strongest hits have smaller values
    if (input$rocauc_cutoff_stat == "Score") {
      baseHits = baseHits %>%
        filter(abs(score) >= abs(as.numeric(input$rocauc_cutoff_value))) %>%
        mutate(stat = abs(score))
    } else if (input$rocauc_cutoff_stat == "FDR") {
      baseHits = baseHits %>%
        filter(fdr <= as.numeric(input$rocauc_cutoff_value)) %>%
        mutate(stat = fdr)
    } else if (input$rocauc_cutoff_stat == "Rank") {
      baseHits = baseHits %>%
        mutate(rank = rank(abs(score))) %>%
        filter(rank <= as.numeric(input$rocauc_cutoff_value)) %>%
        mutate(stat = rank)
    }
    if (nrow(baseHits) < 1) {
      # fail if no remaining hits in base screen
      table = tibble(x = 0, y = 0, label = "No hits in the base screen with the specified cutoff.\nPlease relax the cutoff.")
      p = ggplotly(ggplot(table, aes(x = x, y = y, label = label)) + geom_label() + theme_classic())
      return(list(plotdata = table, p = p))
    }
    # obtain hits in head screens without applying cutoff
    headHits = fetchDs(comparison = comparisons$`Comparison ID`[comparisons$FriendlyID %in% input$rocauc_head], gene = NULL, method = input$rocauc_method)
    # head screens with no hits are automatically excluded
    # apply direction specified
    if (input$rocauc_cutoff_direction == "Hypersensitivity hits") {
      headHits = headHits %>%
        filter(score < 0)
    } else if (input$rocauc_cutoff_direction == "Suppressing hits") {
      headHits = headHits %>%
        filter(score > 0)
    }
    # calculate stat in the same way as for the base screen
    # without applying the cutoff
    if (input$rocauc_cutoff_stat == "Score") {
      headHits = headHits %>%
        #filter(abs(score) >= abs(as.numeric(input$rocauc_cutoff_value))) %>%
        mutate(stat = abs(score))
    } else if (input$rocauc_cutoff_stat == "FDR") {
      headHits = headHits %>%
        #filter(fdr <= as.numeric(input$rocauc_cutoff_value)) %>%
        mutate(stat = fdr)
    } else if (input$rocauc_cutoff_stat == "Rank") {
      headHits = headHits %>%
        mutate(rank = rank(abs(score))) %>%
        #filter(rank <= as.numeric(input$rocauc_cutoff_value)) %>%
        mutate(stat = rank)
    }
    if (nrow(headHits) < 1) {
      # fail if no remaining hits in any head screens
      table = tibble(x = 0, y = 0, label = "No hits in any of the head screens.\nPlease contact the app developer.")
      p = ggplotly(ggplot(table, aes(x = x, y = y, label = label)) + geom_label() + theme_classic())
      return(list(plotdata = table, p = p))
    }
    # calculate precision and recall for each head screen along the entire range of scores/FDRs as specified
    # loop through head screens
    roc = foreach(h = unique(headHits$comparison_id), .combine = "bind_rows") %do% {
      thisHeadHits = headHits %>% filter(comparison_id == h)
      steps = seq(from = 0, to = max(thisHeadHits$stat, na.rm = T), length.out = 2000)
      n = nrow(baseHits) # n: number of base hits
      k = nrow(thisHeadHits) # k: number of head hits
      # step through
      t = foreach(s = steps, .combine = "bind_rows") %do% {
        if (input$rocauc_cutoff_stat == "Score") {
          # If the stat is "Score", find hits with Score s or better
          stepHeadHits = thisHeadHits %>% filter(stat >= s)
        } else {
          # If the stat is "FDR" or "Rank", find hits with FDR or Rank s or better
          stepHeadHits = thisHeadHits %>% filter(stat <= s)
        }
        q = sum(stepHeadHits$gene_id %in% baseHits$gene_id) # q: number of hits at this s that are also in the base hits
        r = nrow(stepHeadHits) - q # r: number of hits at this s that aren't in the base hits
        precision = 1 - (r / k) # precision: freedom from non-base hits at this s
        recall = q / n # recall: recovered base hits at this s as a proportion of total base hits
        tibble(step = s, precision = precision, recall = recall)
      }
      # convert precision and recall to rate
      t$precision = t$precision / max(t$precision)
      t$recall = t$recall / max(t$recall)
      # add the start and end points
      t = bind_rows(t,
                tibble(precision = c(0, 1), recall = c(1, 0)))
      t$head = h
      t$screen = comparisons$FriendlyID[comparisons$`Comparison ID` == h]
      t
    }
    # calculate area under the curve
    # loop through head screens
    auc = foreach(h = unique(headHits$comparison_id), .combine = "bind_rows") %do% {
      thisRoc = roc %>% filter(head == h)
      suppressWarnings({
        thisAuc = AUC(1-thisRoc$precision, thisRoc$recall)
      })
      tibble(head = h, auc = as.numeric(thisAuc))
    }
    roc = roc %>%
      left_join(auc, by = "head") %>%
      mutate(legend = paste0(screen, "\nAUC: ", signif(auc, 4)))
    # plot ROC curve
    p = ggplot(roc, aes(x = 1-precision, y = recall, colour = legend, auc = auc)) +
      geom_line() +
      geom_abline(slope = 1, intercept = 0) +
      scale_colour_discrete(name = "Contrast") +
      theme_classic() +
      xlab("1 - Precision\n(False positive rate)") +
      ylab("Recall\n(True positive rate)") +
      ggtitle(paste0("Precision-recall for <b>", input$rocauc_base, "</b>"))
    # return plot, table
    p = ggplotly(p)
    p
    table = roc %>%
      transmute(Contrast = screen, `ROC-AUC` = auc, Precision = precision, Recall = recall) %>%
      unique()
    return(list(plotdata = table, p = p, height = 900))
  }
}