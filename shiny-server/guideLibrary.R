# Filter a guide library according to inputs
# Todo: implement this as an SQL database
guideLibraryinate <- function(input) {
  
  # Build SQL query
  q <- "SELECT * FROM guideLibrary WHERE "
  cond <- NULL
  # Pick gene classes if specified
  if (length(input$guideLibrary_geneclasses) > 0) {
    cond <- c(cond, paste0("( ", paste0("[Gene Type] = ", paste0("'", input$guideLibrary_geneclasses, "'", collapse = " OR [Gene Type] = ")), " )"))
  }
  # Pick targets if specified
  if (length(input$guideLibrary_targets) > 0) {
    cond <- c(cond, paste0("( ", paste0("[Target] = ", paste0("'", input$guideLibrary_targets, "'", collapse = " OR [Target] = ")), " )"))
  }
  # Pick chromosomes if specified
  if (length(input$guideLibrary_chromosomes) > 0) {
    cond <- c(cond, paste0("( ", paste0("[Alignment] LIKE ", paste0("'chr", input$guideLibrary_chromosomes, ":%'", collapse = " OR [Alignment] LIKE ")), " )"))
  }
  # Pick chemistry if specified
  if (length(input$guideLibrary_chemistry) > 0) {
    cond <- c(cond, paste0("( ", paste0("[Chemistry] = ", paste0("'", input$guideLibrary_chemistry, "'", collapse = " OR [Chemistry] = ")), " )"))
  }
  # Pick assembly if specified
  if (length(input$guideLibrary_assembly) > 0) {
    cond <- c(cond, paste0("( ", paste0("[Assembly] = ", paste0("'", input$guideLibrary_assembly, "'", collapse = " OR [Assembly] = ")), " )"))
  }
  # Pick PAM if specified
  if (length(input$guideLibrary_pam) > 0) {
    cond <- c(cond, paste0("( ", paste0("[PAM] = ", paste0("'", input$guideLibrary_pam, "'", collapse = " OR [PAM] = ")), " )"))
  }
  # Pick library if specified
  if (length(input$guideLibrary_library) > 0) {
    cond <- c(cond, paste0("( ", paste0("[Library] = ", paste0("'", input$guideLibrary_library, "'", collapse = " OR [Library] = ")), " )"))
  }
  # Pick organism if specified
  if (length(input$guideLibrary_organism) > 0) {
    cond <- c(cond, paste0("( ", paste0("[Organism] = ", paste0("'", input$guideLibrary_organism, "'", collapse = " OR [Organism] = ")), " )"))
  }
  conds <- paste0(cond, collapse = " AND ")
  query <- paste0(q, conds)
  # Fail if nothing selected
  if(query == q) {
    guideLibraryFiltered <- tibble(V1 = "Please select at least one option in the sidebar.")
  } else {
    guideLibraryFiltered <- sqldf(connection = guideLibrary, x = query)
    # Refactor
    guideLibraryFiltered$PAM = factor(guideLibraryFiltered$PAM)
    guideLibraryFiltered$Target = factor(guideLibraryFiltered$Target)
    guideLibraryFiltered$Assembly = factor(guideLibraryFiltered$Assembly)
    guideLibraryFiltered$Library = factor(guideLibraryFiltered$Library)
    guideLibraryFiltered$Chemistry = factor(guideLibraryFiltered$Chemistry)
    guideLibraryFiltered$Organism = factor(guideLibraryFiltered$Organism)
  } 
  # Return
  return(guideLibraryFiltered)
}