
plotUpset <- function(input) {
  # Get comparison IDs from contrasts
  comparisonIds = comparisons$`Comparison ID`[comparisons$FriendlyID %in% input$comparison_overlap]
  # Cast cutoff value to numeric
  cutoff_value <- as.numeric(input$overlap_cutoff_value)
  # Fetch scores from the database
  alldata = fetchDs(comparisonIds, NULL, input$overlap_cutoff_method) %>%
    left_join(comparisons, by = c("comparison_id" = "Comparison ID")) %>%
    mutate(screen = FriendlyID) %>%
    group_by(screen) %>%
    transmute(gene = gene_id,
              sign = sign(score),
              rank = rank(score, ties.method = "random"),
              # Fetch the specified score
              value = case_when(input$overlap_cutoff_stat == "Score" ~ score,
                                input$overlap_cutoff_stat == "FDR" ~ fdr,
                                input$overlap_cutoff_stat == "Rank" ~ abs(rank)-n()/2)) %>%
    filter(gene != "X")
  # Apply cutoff
  if(input$overlap_cutoff_direction == "Hypersensitivity hits" & input$overlap_cutoff_stat == "Rank") cutoffdata <- alldata %>% filter(rank < max(rank)/2)
  if(input$overlap_cutoff_direction == "Suppressing hits"      & input$overlap_cutoff_stat == "Rank") cutoffdata <- alldata %>% filter(rank > max(rank)/2)
  if(input$overlap_cutoff_direction == "Hypersensitivity hits" & input$overlap_cutoff_stat != "Rank") cutoffdata <- alldata %>% filter(sign == -1)
  if(input$overlap_cutoff_direction == "Suppressing hits"      & input$overlap_cutoff_stat != "Rank") cutoffdata <- alldata %>% filter(sign == 1)
  if(input$overlap_cutoff_stat == "Score") cutoffdata <- cutoffdata %>% filter(abs(value) >= abs(cutoff_value))
  if(input$overlap_cutoff_stat == "FDR")   cutoffdata <- cutoffdata %>% filter(value <= cutoff_value)
  if(input$overlap_cutoff_stat == "Rank")  cutoffdata <- cutoffdata %>% slice_max(abs(value), n = round(cutoff_value, digits = 0))
  # Determine Venn groups per screen
  venns = foreach(contrast = unique(cutoffdata$screen), .final = function(x) setNames(x, unique(cutoffdata$screen))) %do% {
    cutoffdata %>%
      filter(screen == contrast) %>%
      select(gene) %>%
      unlist() %>%
      unique()
  }
  # Fail if less than 2 groups
  if(length(venns) < 2) {
    hypers = tibble(x = 0, y = 0, label = "Select two or more screens.")
    gg = ggplotly(ggplot(hypers, aes(x = x, y = y, label = label)) +
                    geom_text())
  } else {
    # Do hypergeometric tests
    hypers = foreach(i = 1:length(venns), .combine = "bind_rows") %do% {
      foreach(j = 1:length(venns), .combine = "bind_rows") %do% {
        if (i <= j) NULL
        else {
          # Universe is the set of genes detected in both experiments
          universe = alldata %>%
            filter(screen %in% names(venns)[c(i, j)]) %>%
            group_by(gene) %>%
            summarise(n = length(unique(screen))) %>%
            filter(n == 2) %>%
            select(gene) %>%
            unlist() %>%
            unique()
          # Determine hits per experiment within the universe
          hits_1 = venns[[i]][venns[[i]] %in% universe]
          hits_2 = venns[[j]][venns[[j]] %in% universe]
          # Perform hypergeometric test
          hyper_nn = length(universe)
          hyper_k = length(hits_1)
          hyper_m = length(hits_2)
          hyper_n = hyper_nn - hyper_m
          hyper_q = sum(venns[[i]] %in% venns[[j]])
          hyper_p = phyper(hyper_q, hyper_m, hyper_n, hyper_k, lower.tail = F)
          tibble(`Screen 1` = c(names(venns)[i], names(venns)[j]),
                 `Screen 2` = c(names(venns)[j], names(venns)[i]),
                 `Hits 1` = c(hyper_k, hyper_m),
                 `Hits 2` = c(hyper_m, hyper_k),
                 `Hits 1 & 2` = hyper_q,
                 `Hypergeometric p-value` = hyper_p,
                 `Total genes` = hyper_nn)
        }
      }
    } %>% arrange(`Screen 2`) %>% arrange(`Screen 1`)
    # Round and adjust p-values
    hypers = hypers %>%
      mutate(`Hypergeometric FDR` = signif(p.adjust(`Hypergeometric p-value`, method = "BH"), 4),
             `Hypergeometric p-value` = signif(`Hypergeometric p-value`, 4)) %>%
      relocate(`Screen 1`, `Screen 2`, `Hits 1`, `Hits 2`, `Hits 1 & 2`, `Hypergeometric p-value`, `Hypergeometric FDR`, `Total genes`)
    
    # Combinator helper for upset plot
    combinator <- function(plotdata, categories) {
      combination = list()
      for (i in 1:length(categories)){
        combination[[i]] = case_when(plotdata[[categories[i]]] ~ categories[i], T ~ NA)
      }
      combination2 = list()
      for (j in 1:nrow(plotdata)) {
        combination2[[j]] = NA
        for (i in 1:length(categories)) {
          combination2[[j]] = c(combination2[[j]], combination[[i]][j])
          combination2[[j]] = combination2[[j]][!is.na(combination2[[j]])] 
        }
      }
      combination2
    }
    
    # Plot
    plotdata <- cutoffdata %>%
      select(screen, gene) %>%
      mutate(value = T) %>%
      pivot_wider(names_from = "screen", values_from = "value")
    
    plotdata <- plotdata %>% mutate(combination = combinator(plotdata, names(venns)))
    
    p =
      ggplot(plotdata, aes(x = combination)) +
      geom_bar() +
      scale_x_upset() +
      theme_classic()
    
    gg = ggplotly(p)
  }
  # Return
  overlap <- list(plotdata = hypers, p = p, height = 900)
}