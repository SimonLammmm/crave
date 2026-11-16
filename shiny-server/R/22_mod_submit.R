#### Submit module ####
# Upload a counts matrix and an analysis workbook, run the analyses the workbook
# asks for, and put the results into a dataset.
#
# The work happens in 21_pipeline.R; this module is the form in front of it. Two
# things shape the layout. The first is that a submission is refused before a
# container starts rather than after: the workbook is read as soon as it is
# uploaded, and what CRAVE understood is shown back, so a wrong file is obvious
# while it is still cheap. The second is that the job outlives the click -- a
# Chronos run is minutes -- so the panel polls a log rather than blocking, and
# reconnects to whatever job the dataset is running if the browser is reopened.

#### UI ####

submitUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Submit", icon = icon("upload"),
    sidebarLayout(
      sidebarPanel(
        width = 4,
        tags$h4(tags$strong("Dataset")),
        selectizeInput(ns("dataset"), "Add to", choices = NULL),
        tags$hr(),
        tags$h4(tags$strong("Files")),
        fileInput(ns("workbook"), "Experimental design (.xlsx)",
                  multiple = FALSE, accept = c(".xlsx")),
        fileInput(ns("counts"), "Counts matrix (.tsv, .csv, optionally .gz)",
                  multiple = FALSE,
                  accept = c(".tsv", ".csv", ".txt", ".gz", "text/plain")),
        shinyWidgets::materialSwitch(
          ns("replace"), "Replace this experiment if it is already in the dataset",
          value = FALSE, status = "primary", right = TRUE
        ),
        tags$br(),
        actionButton(ns("submit"), "Submit", icon = icon("play"), class = "btn-success"),
        actionButton(ns("reset"), "Clear", icon = icon("ban"), class = "btn-danger"),
        tags$hr(),
        tags$p(tags$strong("Help")),
        tags$p("The", tags$strong("counts matrix"), "is a table whose first column",
               "holds unique guide ids, whose second holds the gene symbols CRAVE",
               "will show, and whose remaining columns are one sample each."),
        tags$p("The", tags$strong("experimental design"), "workbook describes the",
               "screen in four sheets: Experiment details, Sample details,",
               "Control groups and Analyses. The Analyses sheet is what chooses",
               "the methods; MAGeCK, DrugZ and Chronos read the counts matrix, and",
               "a manual analysis supplies its own table instead."),
        tags$p("Everything submitted is kept inside the dataset, so an experiment",
               "can later be rebuilt from its workbook without running the",
               "analyses again."),
        tags$p("The", tags$strong("raw analysis output"), "download is everything",
               "MAGeCK, DrugZ and Chronos themselves wrote for an experiment,",
               "as a .tar.gz. The tabulated tables are left out: those are what",
               "the dataset is built from, so every plot and every data download",
               "elsewhere in CRAVE already comes from them."),
        tags$p("Analyses run in the background. You can leave this tab, or close",
               "the browser, and come back to it.")
      ),
      mainPanel(
        tags$h4("Submit"),
        tags$p("Run ", extLinkTight(DOI_MAGECK, "MAGeCK"), ", ",
               extLinkTight(DOI_DRUGZ, "DrugZ"), " and ",
               extLinkTight(DOI_CHRONOS, "Chronos"),
               " over a screen and add the results to a dataset, using the ",
               extLinkTight(URL_EXORCISE_REPO, "Exorcise"), " pipeline."),
        tags$br(),
        tags$h5(tags$strong("What was uploaded")),
        uiOutput(ns("summary")),
        tags$br(),
        tags$h5(tags$strong("Raw analysis output")),
        # Hidden rather than disabled while there is nothing to download: a
        # picker with no choices beside a button that errors is worse than the
        # line underneath saying why neither is there yet.
        shinyjs::hidden(tags$div(
          id = ns("rawBlock"),
          selectizeInput(ns("rawExperiment"), "Experiment", choices = NULL,
                         width = "100%"),
          downloadButton(ns("rawDownload"), "Download raw output",
                         icon = icon("download"))
        )),
        uiOutput(ns("rawSummary")),
        tags$br(),
        tags$h5(tags$strong("Job")),
        uiOutput(ns("status")),
        tags$pre(style = "max-height: 28em; overflow: auto; font-size: 85%",
                 textOutput(ns("log"), container = tags$code))
      )
    )
  )
}

#### Server ####

#' @param data Reactive holding the app data object.
#' @param config Deployment configuration.
#' @param onChanged Called after a job finishes successfully, to reload the data.
submitServer <- function(id, data, config, onChanged) {
  moduleServer(id, function(input, output, session) {

    jobId   <- reactiveVal(NULL)
    lastRun <- reactiveVal(NULL)   # the state a finished job was last seen in
    lastExp <- reactiveVal(NULL)   # the experiment this tab last submitted
    rawSeen <- reactiveVal(0)      # bumped when the archives could have changed

    #' The datasets this tab can write to.
    writableDatasets <- reactive({
      ds <- data()$datasets
      names(Filter(function(d) isTRUE(d$writable), ds))
    })

    observeEvent(data(), {
      choices <- writableDatasets()
      updateSelectizeInput(session, "dataset", choices = choices,
                           selected = keepOrFirst(isolate(input$dataset), choices))
    })

    dataset <- reactive({
      nm <- input$dataset
      if (is.null(nm) || !nzchar(nm)) return(NULL)
      data()$datasets[[nm]]
    })

    #### Reading the uploads ####
    # Read once, when the file arrives, rather than on every render: a counts
    # matrix is tens of megabytes and the summary below is redrawn often.
    workbook <- reactive({
      up <- input$workbook
      if (is.null(up)) return(NULL)
      readWorkbook(up$datapath[1])
    })

    counts <- reactive({
      up <- input$counts
      if (is.null(up)) return(NULL)
      readCountsHead(up$datapath[1])
    })

    #' Everything wrong with the submission as it stands, in the order it matters.
    problems <- reactive({
      out <- character(0)
      ds <- dataset()
      if (is.null(ds)) out <- c(out, "Choose a dataset to add to.")
      wb <- workbook()
      if (is.null(wb)) {
        out <- c(out, "Upload an experimental design workbook.")
      } else if (!isTRUE(wb$ok)) {
        out <- c(out, wb$errors)
      }
      cn <- counts()
      needsCounts <- isTRUE(wb$ok) && isTRUE(wb$needsCounts)
      if (needsCounts && is.null(cn)) {
        out <- c(out, paste0(
          "Upload a counts matrix: this workbook asks for ",
          paste(intersect(wb$methods, CRAVE_COUNTS_METHODS), collapse = ", "),
          ", which read one."))
      } else if (!is.null(cn) && !isTRUE(cn$ok)) {
        out <- c(out, cn$errors)
      }
      if (isTRUE(wb$ok) && !is.null(ds) && wb$experimentId %in% existingExperiments() &&
          !isTRUE(input$replace)) {
        out <- c(out, paste0(
          "'", wb$experimentId, "' is already in this dataset. Tick \"Replace\" to ",
          "rebuild it from these files, or change the Experiment name in the ",
          "workbook to add it alongside."))
      }
      held <- if (is.null(ds)) NULL else datasetLock(ds$path)
      if (!is.null(held) && !isTRUE(held$stale)) {
        out <- c(out, paste0("A job is already running against this dataset (",
                             held$jobId, ")."))
      }
      out
    })

    existingExperiments <- reactive({
      ds <- dataset()
      if (is.null(ds)) return(character(0))
      as.character(datasetExperiments(data(), ds$name)$`Experiment ID`)
    })

    #### What was uploaded ####
    output$summary <- renderUI({
      wb <- workbook()
      cn <- counts()
      if (is.null(wb) && is.null(cn)) {
        return(tags$p(tags$em("Nothing uploaded yet.")))
      }
      rows <- list()
      if (!is.null(wb) && isTRUE(wb$ok)) {
        known <- wb$experimentId %in% existingExperiments()
        rows <- c(rows, list(
          tags$li(tags$strong("Experiment: "), wb$experimentId,
                  if (known) tags$span(style = "color: #b36b00",
                                       " (already in this dataset)") else NULL),
          tags$li(tags$strong("Analyses: "), paste(wb$methods, collapse = ", ")),
          if (length(wb$groups)) tags$li(tags$strong("Control groups: "),
                                         paste(wb$groups, collapse = ", ")) else NULL,
          if (!is.null(wb$countsFile)) {
            # The saved name follows the upload's compression, so say which name
            # it will be rather than leaving the workbook's to be assumed.
            saved <- if (is.null(cn)) wb$countsFile else
              countsTargetName(wb$countsFile, isTRUE(cn$gzipped))$name
            tags$li(
              tags$strong("Counts file named in the workbook: "), wb$countsFile,
              tags$span(style = "color: #666",
                        if (identical(saved, wb$countsFile))
                          " — the upload is saved under this name"
                        else paste0(" — the upload is gzipped, so it is saved as ",
                                    saved, ", which Exorcise resolves from that name")))
          } else NULL
        ))
      }
      if (!is.null(cn) && isTRUE(cn$ok)) {
        rows <- c(rows, list(tags$li(
          tags$strong("Counts columns: "), length(cn$columns),
          " (", length(cn$samples), " sample(s): ",
          summariseItems(cn$samples), ")")))
      }
      bad <- problems()
      tagList(
        if (length(rows)) tags$ul(rows) else NULL,
        if (length(bad)) tags$div(
          style = paste0("background-color: #fff9e6; padding: 10px; ",
                         "border-left: 6px solid #fec000"),
          tags$ul(lapply(bad, tags$li))
        ) else tags$p(style = "color: #2e7d32",
                      icon("circle-check"), " Ready to submit.")
      )
    })

    #### Submitting ####
    observeEvent(input$reset, {
      jobId(NULL)
      lastRun(NULL)
      shinyjs::reset("workbook")
      shinyjs::reset("counts")
      shinyWidgets::updateMaterialSwitch(session, "replace", value = FALSE)
    })

    observeEvent(input$submit, {
      bad <- problems()
      if (length(bad)) {
        showNotification(bad[1], type = "warning", duration = 8)
        return()
      }
      ds <- dataset()
      wb <- workbook()

      withProgress(message = "Staging the submission...", value = NULL, {
        staged <- stageSubmission(ds, wb, input$workbook$datapath[1],
                                  if (is.null(input$counts)) NULL else input$counts$datapath[1])
        if (!isTRUE(staged$ok)) {
          showNotification(staged$message, type = "error", duration = 10)
          return()
        }
        started <- submitExperimentJob(ds, config, staged$workbook,
                                       replace = isTRUE(input$replace))
        if (!isTRUE(started$ok)) {
          showNotification(started$message, type = "error", duration = 10)
          return()
        }
        jobId(started$jobId)
        lastRun(NULL)
        lastExp(wb$experimentId)
        showNotification(paste0("Submitted. ", wb$experimentId,
                                " is being analysed; this panel will follow it."),
                         type = "message", duration = 6)
      })
    })

    #### Raw analysis output ####
    # Recomputed only when the set of archives could have changed -- the dataset
    # changed, or a job stopped -- and deliberately not on the job poll, which
    # would move the picker under the hand of anyone using it twice a second.
    rawExperiments <- reactive({
      rawSeen()
      ds <- dataset()
      if (is.null(ds)) return(character(0))
      rawResultExperiments(ds$path)
    })

    observe({
      ids <- rawExperiments()
      shinyjs::toggle("rawBlock", condition = length(ids) > 0)
      # The experiment this tab just submitted is the one to offer; failing that,
      # whatever was already chosen. Isolated, so that choosing one by hand does
      # not re-enter this observer and overwrite the choice.
      want <- isolate(lastExp())
      if (is.null(want) || !want %in% ids) want <- isolate(input$rawExperiment)
      updateSelectizeInput(session, "rawExperiment", choices = ids,
                           selected = keepOrFirst(want, ids))
    })

    output$rawSummary <- renderUI({
      ids <- rawExperiments()
      if (length(ids) == 0) {
        return(tags$p(tags$em(paste0(
          "Nothing to download yet. MAGeCK, DrugZ and Chronos write their own ",
          "output into the dataset as they run, and it is offered here once an ",
          "experiment has some."))))
      }
      ds <- dataset()
      id <- input$rawExperiment
      if (is.null(ds) || is.null(id) || !nzchar(id) || !id %in% ids) return(NULL)
      s <- rawResultSummary(ds$path, id)
      tags$p(style = "color: #666",
             paste0(paste(s$methods, collapse = ", "), ": ", s$files,
                    " file(s), ", prettyBytes(s$bytes), " before compression."))
    })

    output$rawDownload <- downloadHandler(
      filename = function() {
        paste0(safeFilename(input$rawExperiment) %||% "experiment",
               "-raw-output.tar.gz")
      },
      content = function(f) {
        ds <- dataset()
        ok <- !is.null(ds) && rawResultsArchive(ds$path, input$rawExperiment, f)
        if (!isTRUE(ok)) {
          stop("There is no raw analysis output in this dataset for '",
               input$rawExperiment, "'.", call. = FALSE)
        }
      }
    )

    #### Following the job ####
    # Polls only while something is running, so an idle Submit tab costs nothing.
    status <- reactive({
      ds <- dataset()
      if (is.null(ds)) return(NULL)
      id <- jobId() %||% {
        held <- datasetLock(ds$path)
        if (!is.null(held) && !isTRUE(held$stale)) held$jobId else NULL
      }
      if (is.null(id)) return(NULL)
      st <- jobStatus(ds$path, id)
      if (identical(st$state, "running")) invalidateLater(CRAVE_JOB_POLL_MS, session)
      st$jobId <- id
      st
    })

    # One reload when a job finishes, not one per poll.
    observeEvent(status(), {
      st <- status()
      if (is.null(st) || identical(st$state, "running")) return()
      if (identical(lastRun(), st$jobId)) return()
      lastRun(st$jobId)
      # Either way: a failed run can still have left some of its methods' output
      # behind, and that is often the thing worth looking at.
      rawSeen(rawSeen() + 1)
      if (identical(st$state, "ok")) {
        showNotification("Analysis finished and the dataset was updated. Reloading.",
                         type = "message", duration = 6)
        onChanged()
      } else {
        showNotification("The job did not finish. The log below says why.",
                         type = "error", duration = 10)
      }
    })

    output$status <- renderUI({
      st <- status()
      if (is.null(st)) return(tags$p(tags$em("No job running.")))
      switch(
        st$state,
        running = tags$p(icon("spinner"), tags$strong(" Running"),
                         " — ", st$jobId,
                         tags$br(), tags$em("Started ", format(st$started), ".")),
        ok      = tags$p(style = "color: #2e7d32", icon("circle-check"),
                         tags$strong(" Finished"), " — ", st$jobId),
        failed  = tags$p(style = "color: #c62828", icon("circle-exclamation"),
                         tags$strong(" Failed"), " — ", st$jobId,
                         if (!is.na(st$exit)) paste0(" (exit ", st$exit, ")") else ""),
        tags$p(tags$em("No job running."))
      )
    })

    output$log <- renderText({
      st <- status()
      if (is.null(st) || length(st$log) == 0) return("")
      paste(st$log, collapse = "\n")
    })
  })
}

#### Staging ####

#' Put the uploads where the Exorcise container will look for them.
#'
#' The counts table is saved under the name the workbook declares, not the name it
#' was uploaded with: exorcise-analyse resolves the counts file by that name inside
#' --counts-dir, so a file called anything else is simply not found. The one
#' adjustment is compression -- stageCountsFile() has the reasoning. The workbook
#' is saved under its experiment id, which is what makes it findable later by the
#' Manage tab.
stageSubmission <- function(ds, wb, workbookPath, countsPath) {
  if (!ensureSubmitDirs(ds$path)) {
    return(list(ok = FALSE, message = paste0(
      "Could not write to '", file.path(ds$path, CRAVE_SUBMIT_DIR),
      "'. The dataset directory has to be writable to submit into it.")))
  }
  p <- submitPaths(ds$path)
  target <- file.path(p$workbooks, paste0(safeFilename(wb$experimentId), ".xlsx"))

  ok <- tryCatch(file.copy(workbookPath, target, overwrite = TRUE),
                 error = function(e) FALSE)
  if (!isTRUE(ok)) {
    return(list(ok = FALSE, message = "The workbook could not be saved into the dataset."))
  }

  if (!is.null(countsPath) && !is.null(wb$countsFile)) {
    # Under a name that matches the bytes, which is not always the name the
    # workbook declares -- see stageCountsFile().
    staged <- stageCountsFile(countsPath, p$counts, wb$countsFile)
    if (!isTRUE(staged$ok)) return(list(ok = FALSE, message = staged$message))
  }
  log_info(paste0("Staged submission of '", wb$experimentId, "' into dataset '",
                  ds$name, "'."))
  list(ok = TRUE, workbook = target)
}
