#### Shiny server: bulk download function ####
bulkDownload <- function(input) {
  # Apply screen filters
  finalScreens <- gq_apply_filter(input)
  # Obtain data
  if (length(finalScreens) > 0 & length(finalScreens) <= 1000) {
    # If screens found, then proceed
    table <- fetchDs(finalScreens, NULL, NULL) %>%
      left_join(comparisons, by = c("comparison_id" = "Comparison ID")) %>%
      # Make stats readable
      mutate(Method = case_when(analysis_type_id == 1 ~ "MAGeCK score",
                                analysis_type_id == 2 ~ "DrugZ normZ")) %>%
      # Rename fields
      transmute(Citation, Contrast = FriendlyID, `Gene symbol` = gene_id, Method, Score = score, FDR = fdr, `Upper-tail p-value` = pos_p,
                `Lower-tail p-value` = neg_p,
                `Days grown (diff)`, `Treatment (diff)`, `Dose (diff)`, `Knockout (diff)`, `Cell line (diff)`,
                `Days grown (ref)`, `Treatment (ref)`, `Dose (ref)`, `Knockout (ref)`, `Cell line (ref)`,
                Library, Source)
    screens <- sort(unique(gsub("<.+?>", "", table$Citation)))
    console <- paste0("Found ", nrow(table), " entries among ", length(screens), " screens.
                      Screens found: ", paste0(screens, collapse = ", "), ".
                      Showing the first six entries here. Click the download button to download them all.")
  } else if (length(finalScreens) > 1000) {
    console <- "Query too large. Limit your search, or contact the authors if such a large query is necessary."
    table <- tibble("That didn't work." = "")
  }
    else {
    # If no screens found, then fail
    console <- "No screens found with the specified filters."
    table <- tibble("That didn't work." = "")
  }
  # Interpret file format
  filename <- input$download_filename
  delim <- case_when(grepl("\\.tsv.*$", filename) ~ "\t",
                     grepl("\\.txt.*$", filename) ~ "\t",
                     T ~ ",")
  # Return
  return(list(console = console, table = table, filename = filename, delim = delim))
}