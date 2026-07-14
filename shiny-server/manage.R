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

# Remove one or more experiments from a dataset's database.db via the Exorcise
# image's `crispr-screen-viewer remove` subcommand. Backs up database.db first.
# Returns list(ok, console).
manage_remove_experiments <- function(dataset_dir, exp_ids) {
  dataset_dir <- normalizePath(dataset_dir, mustWork = FALSE)
  db <- file.path(dataset_dir, "database.db")
  if (file.exists(db)) file.copy(db, paste0(db, ".bak"), overwrite = TRUE)
  
  if (!nzchar(get0("exorcise_docker", ifnotfound = ""))) {
    return(list(ok = FALSE, console = "exorcise_docker is not set in config.R."))
  }
  platform <- get0("exorcise_platform", ifnotfound = "")
  platflag <- if (nzchar(platform)) paste0(" --platform ", platform) else ""
  
  ids <- paste(shQuote(exp_ids), collapse = " ")
  command <- paste0(
    "docker run --rm", platflag,
    " -v ", dataset_dir, ":/d ",
    exorcise_docker,
    " crispr-screen-viewer remove /d ", ids,
    " 2>&1"
  )
  
  console <- tryCatch(
    paste(system(command, intern = TRUE), collapse = "\n"),
    error = function(e) paste("Failed to run docker:", conditionMessage(e))
  )
  ok <- !grepl("Traceback \\(most recent call last\\)", console)
  list(ok = ok, console = console)
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
      DTOutput("manage_comparisons"),
      tags$hr(),
      tags$h4("Remove experiment"),
      tags$p("Deletes an experiment's rows from the database (stat/comparison tables). A .bak backup of database.db is written first. This cannot be undone except by restoring the backup."),
      fluidRow(
        column(6, selectizeInput("remove_exp_ids", "Experiment ID(s) to remove", choices = NULL, multiple = TRUE)),
        column(6, tags$br(), actionButton("remove_exp_go", "Remove experiment(s)", icon = icon("trash"), class = "btn-danger"))
      ),
      textOutput("remove_exp_status"),
      tags$pre(style = "max-height: 200px; overflow-y: auto; font-size: 11px;",
               textOutput("remove_exp_console"))
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
  
  # Client-side tables that re-render from rv (no server-side search state to stick)
  output$manage_experiments <- renderDT({
    datatable(rv$exp, editable = list(target = "cell", disable = list(columns = exp_lock)),
              rownames = FALSE, options = list(iDisplayLength = 25, scrollX = TRUE))
  }, server = FALSE)
  
  output$manage_comparisons <- renderDT({
    datatable(rv$comp, editable = list(target = "cell", disable = list(columns = comp_lock)),
              rownames = FALSE, options = list(iDisplayLength = 25, scrollX = TRUE))
  }, server = FALSE)
  
  # Switch dataset -> reload raw CSVs; tables re-render reactively from rv
  observeEvent(input$manage_dataset, {
    req(input$manage_dataset)
    d <- manage_read(input$manage_dataset)
    rv$exp <- d$exp; rv$comp <- d$comp
    rv$exp_path <- d$exp_path; rv$comp_path <- d$comp_path
    output$manage_validation <- renderUI(NULL)
    output$manage_status <- renderText("")
    updateSelectizeInput(session, "remove_exp_ids", choices = sort(unique(d$exp[["Experiment ID"]])), server = TRUE)
  }, ignoreInit = FALSE)
  
  # Persist cell edits into the in-memory copies (coerce the single edited cell)
  observeEvent(input$manage_experiments_cell_edit, {
    info <- input$manage_experiments_cell_edit
    tmp <- rv$exp
    tmp[info$row, info$col + 1] <- info$value
    rv$exp <- tmp
  })
  observeEvent(input$manage_comparisons_cell_edit, {
    info <- input$manage_comparisons_cell_edit
    tmp <- rv$comp
    tmp[info$row, info$col + 1] <- info$value
    rv$comp <- tmp
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
  
  ## --- Remove experiment (task 5 stage 3) ---
  # Confirm before removing -- this modifies database.db permanently (backup aside).
  observeEvent(input$remove_exp_go, {
    req(input$remove_exp_ids)
    showModal(modalDialog(
      title = "Confirm removal",
      paste0("Remove ", length(input$remove_exp_ids), " experiment(s) (",
             paste(input$remove_exp_ids, collapse = ", "), ") from ",
             input$manage_dataset, "? A backup of database.db will be written first, ",
             "but this action modifies the live dataset."),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("remove_exp_confirm", "Remove", class = "btn-danger")
      )
    ))
  })
  
  observeEvent(input$remove_exp_confirm, {
    removeModal()
    req(input$remove_exp_ids, input$manage_dataset)
    dataset_dir <- manage_dataset_path(input$manage_dataset)
    output$remove_exp_status <- renderText(paste0("Removing ", length(input$remove_exp_ids), " experiment(s)..."))
    res <- manage_remove_experiments(dataset_dir, input$remove_exp_ids)
    output$remove_exp_console <- renderText(res$console)
    if (res$ok) {
      output$remove_exp_status <- renderText(paste0(
        "Removed ", paste(input$remove_exp_ids, collapse = ", "), " from ", input$manage_dataset,
        ". Backup at database.db.bak. Use \"Refresh data\" to update the other tabs."))
      # Refresh the experiment ID choices (removed ones should drop out on next dataset read)
      d <- manage_read(input$manage_dataset)
      updateSelectizeInput(session, "remove_exp_ids", choices = sort(unique(d$exp[["Experiment ID"]])), server = TRUE)
    } else {
      output$remove_exp_status <- renderText("Remove failed — see log below. database.db.bak has the pre-removal backup.")
    }
  })
}