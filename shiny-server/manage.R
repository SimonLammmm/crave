#### Metadata management ####
## In-app editing of experiments_metadata.csv.gz and comparisons_metadata.csv.gz.
## Gated by `enable_editing` (config.R). Writes back to the raw CSVs (with a .bak
## backup) and re-derives everything via load(). Cell editing only; key columns
## (IDs, Source) are locked because the database links to them.

## --- Helpers -----------------------------------------------------------------

# Resolve a dataset's directory from the config registry
manage_dataset_path <- function(nm) {
  entry <- Find(function(d) d$name == nm, datasets)
  if (is.null(entry)) return(NULL)
  entry$path
}

# Read the raw (untransformed) metadata CSVs for one dataset
manage_read <- function(nm) {
  path <- manage_dataset_path(nm)
  exp_path  <- file.path(path, "experiments_metadata.csv.gz")
  comp_path <- file.path(path, "comparisons_metadata.csv.gz")
  exp  <- if (file.exists(exp_path))  as.data.frame(data.table::fread(exp_path,  colClasses = "character")) else data.frame()
  comp <- if (file.exists(comp_path)) as.data.frame(data.table::fread(comp_path, colClasses = "character")) else data.frame()
  list(exp = exp, comp = comp, exp_path = exp_path, comp_path = comp_path)
}

# Write metadata back, taking a .bak backup of each file first
manage_write <- function(exp_df, comp_df, exp_path, comp_path) {
  if (file.exists(exp_path))  file.copy(exp_path,  paste0(exp_path,  ".bak"), overwrite = TRUE)
  if (file.exists(comp_path)) file.copy(comp_path, paste0(comp_path, ".bak"), overwrite = TRUE)
  data.table::fwrite(exp_df,  exp_path)   # .gz extension => gzip-compressed
  data.table::fwrite(comp_df, comp_path)
}

# Validate one dataset's metadata against itself and against database.db
# Returns an HTML tagList describing any problems found.
manage_validate <- function(nm, exp_df, comp_df) {
  issues <- list()

  expIDs  <- unique(exp_df[["Experiment ID"]])
  compExp <- unique(comp_df[["Experiment ID"]])
  metaCmp <- unique(comp_df[["Comparison ID"]])

  # 1. Comparisons whose Experiment ID has no experiment metadata row
  orphan_comp <- setdiff(compExp, expIDs)
  if (length(orphan_comp) > 0)
    issues <- c(issues, list(paste0("Comparisons reference Experiment IDs with no experiment metadata: ",
                                     paste(orphan_comp, collapse = ", "))))

  # 2. Experiments with no comparisons
  orphan_exp <- setdiff(expIDs, compExp)
  if (length(orphan_exp) > 0)
    issues <- c(issues, list(paste0("Experiments have no comparisons: ",
                                     paste(orphan_exp, collapse = ", "))))

  # 3. Metadata <-> database.db comparison_id consistency
  if (!is.null(cons[[nm]])) {
    dbCmp <- tryCatch(
      DBI::dbGetQuery(cons[[nm]], "SELECT DISTINCT comparison_id FROM stat")$comparison_id,
      error = function(e) NULL)
    if (!is.null(dbCmp)) {
      # Only check metadata -> database (catches metadata describing screens that
      # don't exist). The reverse is intentionally NOT checked: the stat table
      # legitimately contains self-comparisons, reversed-direction contrasts, and
      # per-sample reference points that are not meant to have metadata rows.
      meta_not_db <- setdiff(metaCmp, dbCmp)
      if (length(meta_not_db) > 0)
        issues <- c(issues, list(paste0("Comparison IDs in metadata but not in database.db: ",
                                         paste(meta_not_db, collapse = ", "))))
    } else {
      issues <- c(issues, list("Could not read comparison_id from database.db (no 'stat' table?)."))
    }
  }

  # 4. Source column should match the dataset name
  if ("Source" %in% names(comp_df)) {
    bad_src <- setdiff(unique(comp_df[["Source"]]), nm)
    bad_src <- bad_src[!is.na(bad_src) & bad_src != ""]
    if (length(bad_src) > 0)
      issues <- c(issues, list(paste0("Comparisons have a Source other than \"", nm, "\": ",
                                       paste(bad_src, collapse = ", "))))
  }

  if (length(issues) == 0) {
    tags$div(style = "color: #1a7f37;",
             tags$strong("\u2713 No problems found."),
             tags$p(paste0(length(expIDs), " experiments, ", nrow(comp_df), " comparisons checked.")))
  } else {
    tags$div(style = "color: #b35900;",
             tags$strong(paste0("\u26a0 ", length(issues), " issue(s) found:")),
             tags$ul(lapply(issues, tags$li)))
  }
}

## --- UI -----------------------------------------------------------------------
## Only built when editing is enabled; otherwise NULL (navbarPage skips NULL tabs).

manageTab <- if (isTRUE(get0("enable_editing", ifnotfound = FALSE))) {
  tabPanel(
    "Manage", icon = icon("table-list"),
    fluidPage(
      tags$h3("Metadata management"),
      tags$p("Edit experiment and comparison metadata, then save to write it back to the dataset and reload. Key columns (IDs, Source) are locked. A .bak backup is written on every save."),
      fluidRow(
        column(4, selectInput("manage_dataset", "Dataset", choices = NULL)),
        column(4, tags$br(), actionButton("manage_validate", "Validate", icon = icon("circle-check"))),
        column(4, tags$br(), actionButton("manage_save", "Save changes", icon = icon("floppy-disk"), class = "btn-warning"))
      ),
      uiOutput("manage_validation"),
      textOutput("manage_status"),
      tags$hr(),
      tags$h4("Experiments metadata"),
      DTOutput("manage_experiments"),
      tags$hr(),
      tags$h4("Comparisons metadata"),
      DTOutput("manage_comparisons")
    )
  )
} else NULL

## --- Server -------------------------------------------------------------------
## Call manageServer(input, output, session) once inside server().

manageServer <- function(input, output, session) {
  if (!isTRUE(get0("enable_editing", ifnotfound = FALSE))) return(invisible())
  if (length(cons) == 0) return(invisible())

  # Prime with the first dataset so the tables render with correct columns
  init_nm <- names(cons)[1]
  .init   <- manage_read(init_nm)
  rv <- reactiveValues(exp = .init$exp, comp = .init$comp,
                       exp_path = .init$exp_path, comp_path = .init$comp_path)

  # Lock key columns (0-based indices, rownames = FALSE)
  exp_lock  <- which(names(.init$exp)  %in% c("Experiment ID")) - 1
  comp_lock <- which(names(.init$comp) %in% c("Experiment ID", "Comparison ID", "Source")) - 1

  updateSelectInput(session, "manage_dataset", choices = names(cons), selected = init_nm)

  output$manage_experiments <- renderDT(
    isolate(rv$exp),
    editable = list(target = "cell", disable = list(columns = exp_lock)),
    rownames = FALSE, filter = "top", server = TRUE,
    options = list(iDisplayLength = 25, scrollX = TRUE)
  )
  output$manage_comparisons <- renderDT(
    isolate(rv$comp),
    editable = list(target = "cell", disable = list(columns = comp_lock)),
    rownames = FALSE, filter = "top", server = TRUE,
    options = list(iDisplayLength = 25, scrollX = TRUE)
  )
  proxy_exp  <- dataTableProxy("manage_experiments")
  proxy_comp <- dataTableProxy("manage_comparisons")

  # Switch dataset -> reload raw CSVs into the tables
  observeEvent(input$manage_dataset, {
    req(input$manage_dataset)
    d <- manage_read(input$manage_dataset)
    rv$exp <- d$exp; rv$comp <- d$comp
    rv$exp_path <- d$exp_path; rv$comp_path <- d$comp_path
    replaceData(proxy_exp,  rv$exp,  rownames = FALSE, resetPaging = FALSE, clearSelection = "all")
    replaceData(proxy_comp, rv$comp, rownames = FALSE, resetPaging = FALSE, clearSelection = "all")
    # Clear any leftover column filters / global search from the previous dataset
    updateSearch(proxy_exp,  keywords = list(global = "", columns = rep("", ncol(rv$exp))))
    updateSearch(proxy_comp, keywords = list(global = "", columns = rep("", ncol(rv$comp))))
    output$manage_validation <- renderUI(NULL)
    output$manage_status <- renderText("")
  }, ignoreInit = TRUE)

  # Persist cell edits into the in-memory copies
  observeEvent(input$manage_experiments_cell_edit, {
    rv$exp <- editData(rv$exp, input$manage_experiments_cell_edit, proxy_exp, rownames = FALSE)
  })
  observeEvent(input$manage_comparisons_cell_edit, {
    rv$comp <- editData(rv$comp, input$manage_comparisons_cell_edit, proxy_comp, rownames = FALSE)
  })

  # Validate the current (edited, unsaved) data
  observeEvent(input$manage_validate, {
    output$manage_validation <- renderUI(manage_validate(input$manage_dataset, rv$exp, rv$comp))
  })

  # Save -> backup + write + reload
  observeEvent(input$manage_save, {
    res <- tryCatch({
      manage_write(rv$exp, rv$comp, rv$exp_path, rv$comp_path)
      load()  # re-derive experiments/comparisons/etc. from the new CSVs
      "ok"
    }, error = function(e) conditionMessage(e))
    if (identical(res, "ok")) {
      output$manage_status <- renderText(paste0(
        "Saved ", input$manage_dataset,
        " and reloaded. Use \"Refresh data\" to update the other tabs. Backups written as *.csv.gz.bak."))
    } else {
      output$manage_status <- renderText(paste0("Save failed: ", res))
    }
  })
}
