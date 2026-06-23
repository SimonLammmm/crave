pendragonate <- function(input) {
  # Fetch
  if(input$pendragonator_goi) {
    pendragonated = fetchDs(comparison = gq_apply_filter(input),
                               gene = input$genes_queried,
                               method = input$analysis_method_gq)
  } else {
    pendragonated = fetchDs(comparison = gq_apply_filter(input),
                               gene = NULL,
                               method = input$analysis_method_gq)
  }
  pendragonated = pendragonated %>%
    # Filter trash
    filter(gene_id != "X") %>%
    filter(!grepl("Non-targeting", gene_id)) %>%
    # Find the single best p-value for each gene
    group_by(gene_id) %>%
    mutate(stat = case_when(input$pendragonator_stat == "Score" ~ -abs(score),
                            input$pendragonator_stat == "FDR" ~ fdr),
           order = rank(stat, ties.method = "random")) %>%
    arrange(order) %>%
    filter(order == min(order)) %>%
    select(-order) %>%
    ungroup() %>%
    arrange(stat) %>%
    head(input$pendragonator_n) %>%
    left_join(comparisons, by = c("comparison_id" = "Comparison ID")) %>%
    # Return the experiment that produced the hit
    transmute(Gene = gene_id, Score = score, FDR = fdr,
              Method = case_when(analysis_type_id == 1 ~ "MAGeCK", analysis_type_id == 2 ~ "DrugZ"),
              Contrast,
              `Days grown (diff)`, `Treatment (diff)`, `Dose (diff)`, `Knockout (diff)`, `Cell line (diff)`,
              `Days grown (ref)`, `Treatment (ref)`, `Dose (ref)`, `Knockout (ref)`, `Cell line (ref)`,
              Library, Organism, Citation, Source)
  return(pendragonated)
}