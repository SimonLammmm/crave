#### Common functions ####
# Negative log scale
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

# fread with glob
fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)

# Generate SQL queries by string building
queryDs0 <- function(comparison, gene, method) {
  q <- "SELECT * FROM stat WHERE "
  cond <- NULL
  # Pick comparisons if specified
  if (!is.null(comparison))    cond <- c(cond, paste0("( ", paste0("[comparison_id] = ", paste0("'", comparison, "'", collapse = " OR [comparison_id] = ")), " )"))
  # Select from a table according to the method and stats specified
  if(!is.null(method)) {
    if (method == "DrugZ")       cond <- c(cond, paste0("[analysis_type_id] = 2"))
    else if (method == "MAGeCK") cond <- c(cond, paste0("[analysis_type_id] = 1"))
  }
  # Return genes (rows) if specified
  if (!is.null(gene))          cond <- c(cond, paste0("( ", paste0("[gene_id] = ", paste0("'", gene, "'", collapse = " OR [gene_id] = ")), " )"))
  conds <- paste0(cond, collapse = " AND ")
  q <- paste0(q, conds)
  return (q)
}

# Generate SQL queries by parameterisation
queryDs <- function(comparison, gene, method) {
  # SELECT * FROM stat WHERE [analysis_type_id] IN (?) AND ( [comparison_id] IN (?, ..., ?) ) AND ( [gene_id] IN (?, ..., ?) )
  # parameters: (?, ..., ?)
  q <- "SELECT * FROM stat WHERE "
  cond <- NULL
  p <- NULL
  # Pick comparisons if specified
  if (!is.null(comparison)) {
    cond <- c(cond,
              paste0("[comparison_id] IN (",
                     paste(rep("?", length(comparison)), collapse = ","),
                     ")")
    )
    p <- c(p, comparison)
  }
  # Select from a table according to the method and stats specified
  if (!is.null(method)) {
    method = unique(method)
    method = method[method %in% c("MAGeCK", "DrugZ", "Chronos", "Manual")]
    cond <- c(cond,
              paste0("[analysis_type_id] IN (",
                     paste(rep("?", length(method)), collapse = ","),
                     ")")
    )
    isMAGeCK = "MAGeCK" %in% method
    isDrugZ = "DrugZ" %in% method
    isChronos = "Chronos" %in% method
    isManual = "Manual" %in% method
    if (isMAGeCK) p <- c(p, "1")
    if (isDrugZ) p <- c(p, "2")
    if (isChronos) p <- c(p, "3")
    if (isManual) p <- c(p, "4")
  }
  # Return genes (rows) if specified
  if (!is.null(gene)) {
    cond <- c(cond,
              paste0("[gene_id] IN (",
                     paste(rep("?", length(gene)), collapse = ","),
                     ")")
    )
    p <- c(p, gene)
  }
  # Finish building
  conds <- paste0(cond, collapse = " AND ")
  q <- paste0(q, conds)
  p <- setNames(p, NULL)
  return (list(q = q, p = p))
}

# Fetch data from databases
fetchDs <- function(comparison, gene, method) {
  # Build an empty frame with the right columns from any loaded database
  fetched <- bind_rows(lapply(cons, function(con) sqldf("SELECT * FROM stat LIMIT 0", connection = con)))
  
  log_info(paste0("Database query detected with comparisons: ", paste(comparison, collapse = ","), "; genes: ",
                  paste(gene, collapse = ","), "; method: ", paste(method, collapse = ","), "."))
  # Query each dataset's database in turn
  for (source in names(cons)) {
    # choose the appropriate database to query each comparison, if specified
    if(!is.null(comparison)) {
      comparison_this <- comparisons %>%
        filter(Source == source) %>%
        filter(`Comparison ID` %in% comparison) %>%
        select(`Comparison ID`) %>%
        unlist()
      # and do so if yes
      if(length(comparison_this) > 0) {
        log_info(paste0("Querying ", source, " database with comparisons: ", paste(comparison_this, collapse = ","), "."))
        query = queryDs(comparison_this, gene, method)
        resq <- dbSendQuery(statement = query$q,
                            params = query$p,
                            conn = cons[[source]])
        fetched <- bind_rows(fetched,
                             dbFetch(resq)) %>%
          mutate(fdr = as.numeric(as.character(signif(fdr, 4))), score = round(score, 4)) # rounding
        dbClearResult(resq)
      }
    } else {
      # if not specified, query all
      comparison_this <- NULL
      log_info(paste0("Querying ", source, " database."))
      query = queryDs(comparison_this, gene, method)
      resq <- dbSendQuery(statement = query$q,
                          params = query$p,
                          conn = cons[[source]])
      fetched <- bind_rows(fetched,
                           dbFetch(resq)) %>%
        mutate(fdr = as.numeric(as.character(signif(fdr, 4))), score = round(score, 4)) # rounding
      dbClearResult(resq)
    }
  }
  log_info("Query done.")
  return(fetched)
}

#### Shiny server: correlate common apply filters function ####
gq_apply_filter <- function(input) {
  theseScreens <- comparisons
  # If no overrides, then apply filters
  # Experiment
  if (!is.null(input$gq_filter_citation)) theseScreens <- theseScreens %>% filter(grepl(input$gq_filter_citation, `Citation`, fixed = T))
  if (!is.null(input$gq_filter_library)) theseScreens <- theseScreens %>% filter(`Library` %in% input$gq_filter_library)
  if (!is.null(input$gq_filter_source)) theseScreens <- theseScreens %>% filter(`Source` %in% input$gq_filter_source)
  if (!is.null(input$gq_filter_organism)) theseScreens <- theseScreens %>% filter(`Organism` %in% input$gq_filter_organism)
  if (!is.null(input$gq_filter_kind)) theseScreens <- theseScreens %>% filter(`Kind` %in% input$gq_filter_kind)
  if (!is.null(input$gq_filter_timepoint)) theseScreens <- theseScreens %>% filter(`Timepoint` %in% input$gq_filter_timepoint)
  # Numerator levels
  if (!is.null(input$gq_filter_days_grown_diff)) theseScreens <- theseScreens %>% filter(`Days grown (diff)` %in% input$gq_filter_days_grown_diff)
  if (!is.null(input$gq_filter_treatment_diff)) theseScreens <- theseScreens %>% filter(`Treatment (diff)` %in% input$gq_filter_treatment_diff)
  if (!is.null(input$gq_filter_dose_diff)) theseScreens <- theseScreens %>% filter(`Dose (diff)` %in% input$gq_filter_dose_diff)
  if (!is.null(input$gq_filter_knockout_diff)) theseScreens <- theseScreens %>% filter(`Knockout (diff)` %in% input$gq_filter_knockout_diff)
  if (!is.null(input$gq_filter_cellline_diff)) theseScreens <- theseScreens %>% filter(`Cell line (diff)` %in% input$gq_filter_cellline_diff)
  # Denominator levels
  if (!is.null(input$gq_filter_days_grown_ref)) theseScreens <- theseScreens %>% filter(`Days grown (ref)` %in% input$gq_filter_days_grown_ref)
  if (!is.null(input$gq_filter_treatment_ref)) theseScreens <- theseScreens %>% filter(`Treatment (ref)` %in% input$gq_filter_treatment_ref)
  if (!is.null(input$gq_filter_dose_ref)) theseScreens <- theseScreens %>% filter(`Dose (ref)` %in% input$gq_filter_dose_ref)
  if (!is.null(input$gq_filter_knockout_ref)) theseScreens <- theseScreens %>% filter(`Knockout (ref)` %in% input$gq_filter_knockout_ref)
  if (!is.null(input$gq_filter_cellline_ref)) theseScreens <- theseScreens %>% filter(`Cell line (ref)` %in% input$gq_filter_cellline_ref)
  # Custom
  if (!is.null(input$gq_filter_contrast)) theseScreens <- theseScreens %>% filter(`FriendlyID` %in% input$gq_filter_contrast)
  if (!is.null(input$gq_filter_custom)) if(input$gq_filter_custom != "") theseScreens <- theseScreens %>% filter(grepl(input$gq_filter_custom, FriendlyID, fixed = T))
  
  # Return
  finalScreens <- theseScreens %>% select(`Comparison ID`) %>% unique() %>% unlist() %>% setNames(NULL)
  return(finalScreens)
}

#### Shiny server: correlate apply gene filters function
gq_gene_filter <- function(input) {
  theseGenes <- NULL
  # If genes are specified, then simply accept those genes
  if (!is.null(input$genes_queried)) {
    theseGenes <- c(theseGenes, input$genes_queried)
    # If not specified, then filter from the picker
  }
  theseGenes <- unique(theseGenes)
  return(theseGenes)
}

# Format a contrast
contrastFormatter <- function(comparison) {
  stem <- sub("^(.+?): .+$", "\\1", comparison)
  comparison <- sub("^.+?: ", "", comparison)
  if(grepl("^Day", comparison)) {
    formatted <- sub("^D", "Essentialome from d", comparison) %>% contrastFormatter()
  } else if(grepl("^\\(", comparison)) {
    formatted <- sub("^\\(.+?\\)\\s*", "", comparison) %>% contrastFormatter()
  } else if(grepl("^Essentialome", comparison)) {
    formatted <- sub("^(Essentialome.*?)((, | at | in | under )(.+))?$", "<b>\\1</b><i>\\2</i>", comparison)
  } else if(grepl("^Sorted", comparison)) {
    formatted <- sub("^(Sorted.*?)((, | at | in | under )(.+))?$", "<b>\\1</b><i>\\2</i>", comparison)
  } else {
    formatted <- sub("^(.+?)((, | at | in | under )(.+))?$", "<b>Effect of \\1</b><i>\\2</i>", comparison)
  }
  formatted <- paste0(formatted, ": ", stem)
  return(formatted) 
}

# Toupper first initial
toupper_first_initial <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}