#### Manage module ####
# What a dataset holds, and the three ways it may be changed.
#
#   Relabel   replace the metadata tables with edited copies. CRAVE writes these
#             itself, so replaceMetadataTables() refuses anything that would alter
#             which experiments and comparisons exist -- it can rename, and nothing
#             else.
#   Rebuild   republish an experiment from a corrected workbook, reusing the result
#             tables kept from its submission. No analysis is re-run.
#   Remove    drop experiments from the database and the metadata together.
#
# Rebuild and Remove are exorcise-database invocations, which rewrite the database
# and both metadata files in one pass; Relabel is checked against the database
# before it is allowed. That is the whole of the consistency argument, and the
# integrity panel at the top is there to show it holding.

#### UI ####

manageUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Manage", icon = icon("database"),
    fluidPage(
      tags$br(),
      tags$h4("Manage"),
      tags$p("Inspect and change what is in a dataset. Every change rewrites the",
             "statistics and the metadata that describes them together, so the two",
             "cannot drift apart."),
      inputGrid(
        selectizeInput(ns("dataset"), "Dataset", choices = NULL),
        tagList(
          tags$label("Consistency"),
          uiOutput(ns("integrity")),
          actionLink(ns("recheck"), "Re-check")
        )
      ),
      uiOutput(ns("lockNotice")),
      tags$hr(),

      tags$h4("Experiments"),
      tags$p("Select rows to act on them."),
      DTOutput(ns("experiments")),
      tags$br(),
      actionButton(ns("remove"), "Remove selected", icon = icon("trash"),
                   class = "btn-danger"),
      tags$hr(),

      tags$h4("Rebuild an experiment from a corrected workbook"),
      tags$p("Upload a replacement experimental design for an experiment that is",
             "already here. Its analyses are not run again: the result tables kept",
             "from the original submission are republished under the new",
             "description. Use this to correct a contrast, a treatment or a",
             "citation at source."),
      inputGrid(
        fileInput(ns("workbook"), "Experimental design (.xlsx)",
                  multiple = FALSE, accept = c(".xlsx")),
        tagList(
          uiOutput(ns("rebuildSummary")),
          actionButton(ns("rebuild"), "Rebuild", icon = icon("rotate"),
                       class = "btn-success")
        )
      ),
      tags$hr(),

      tags$h4("Metadata tables"),
      tags$p("Download the tables, edit them, and upload them back to relabel what",
             "is already here. Rows may be rewritten but not added or removed, and",
             "no comparison may change experiment: those are what the statistics",
             "are keyed on. To add an experiment use", tags$strong("Submit", .noWS = noWS),
             "; to drop one use", tags$strong("Remove", .noWS = noWS), "."),
      inputGrid(
        tagList(
          downloadButton(ns("downloadExperiments"), "Download experiments_metadata"),
          tags$br(), tags$br(),
          fileInput(ns("uploadExperiments"), "Replacement experiments_metadata",
                    multiple = FALSE, accept = c(".csv", ".gz", ".tsv"))
        ),
        tagList(
          downloadButton(ns("downloadComparisons"), "Download comparisons_metadata"),
          tags$br(), tags$br(),
          fileInput(ns("uploadComparisons"), "Replacement comparisons_metadata",
                    multiple = FALSE, accept = c(".csv", ".gz", ".tsv"))
        )
      ),
      actionButton(ns("applyMetadata"), "Apply metadata", icon = icon("floppy-disk"),
                   class = "btn-success"),
      uiOutput(ns("metadataNotice")),
      tags$hr(),

      tags$h4("Job"),
      uiOutput(ns("status")),
      tags$pre(style = "max-height: 24em; overflow: auto; font-size: 85%",
               textOutput(ns("log"), container = tags$code)),
      tags$br()
    )
  )
}

#### Server ####

#' @param data Reactive holding the app data object.
#' @param config Deployment configuration.
#' @param onChanged Called after the dataset has been changed, to reload it.
manageServer <- function(id, data, config, onChanged) {
  moduleServer(id, function(input, output, session) {

    jobId   <- reactiveVal(NULL)
    lastRun <- reactiveVal(NULL)
    notice  <- reactiveVal(NULL)
    recheck <- reactiveVal(0)

    observeEvent(data(), {
      choices <- names(Filter(function(d) isTRUE(d$writable), data()$datasets))
      updateSelectizeInput(session, "dataset", choices = choices,
                           selected = keepOrFirst(isolate(input$dataset), choices))
    })

    dataset <- reactive({
      nm <- input$dataset
      if (is.null(nm) || !nzchar(nm)) return(NULL)
      data()$datasets[[nm]]
    })

    #### Integrity ####
    integrity <- reactive({
      recheck()
      ds <- dataset()
      if (is.null(ds)) return(NULL)
      datasetIntegrity(data(), ds$name)
    })

    observeEvent(input$recheck, recheck(recheck() + 1))

    output$integrity <- renderUI({
      res <- integrity()
      if (is.null(res)) return(tags$p(tags$em("No dataset selected.")))
      if (isTRUE(res$ok)) {
        tags$p(style = "color: #2e7d32", icon("circle-check"), " ", res$message)
      } else {
        tags$p(style = "color: #c62828", icon("triangle-exclamation"), " ", res$message)
      }
    })

    #### Lock ####
    lock <- reactive({
      ds <- dataset()
      if (is.null(ds)) return(NULL)
      status()  # so that a finishing job clears this notice too
      datasetLock(ds$path)
    })

    output$lockNotice <- renderUI({
      held <- lock()
      if (is.null(held)) return(NULL)
      if (!isTRUE(held$stale)) {
        return(noticeBanner(tagList(
          "A job is running against this dataset (", held$jobId,
          "). Changes are refused until it finishes."), "info"))
      }
      noticeBanner(tagList(
        "A previous job (", held$jobId, ") left a lock behind without finishing. ",
        actionLink(session$ns("clearLock"), "Clear it"),
        " to allow changes again."), "warning")
    })

    observeEvent(input$clearLock, {
      ds <- dataset()
      if (is.null(ds)) return()
      clearDatasetLock(ds$path)
      showNotification("Lock cleared.", type = "message", duration = 4)
    })

    #### Experiments ####
    experiments <- reactive({
      ds <- dataset()
      if (is.null(ds)) return(tibble(` ` = "No dataset selected."))
      out <- datasetExperiments(data(), ds$name)
      if (nrow(out) == 0) tibble(` ` = "This dataset holds no experiments.") else out
    })

    output$experiments <- renderDT(
      experiments(), filter = "top", server = TRUE, escape = FALSE,
      options = list(pageLength = 25, scrollX = TRUE)
    )

    selectedExperiments <- reactive({
      rows <- input$experiments_rows_selected
      tbl <- experiments()
      if (length(rows) == 0 || !"Experiment ID" %in% names(tbl)) return(character(0))
      as.character(tbl$`Experiment ID`[rows])
    })

    #### Remove ####
    observeEvent(input$remove, {
      ds <- dataset()
      picked <- selectedExperiments()
      if (is.null(ds) || length(picked) == 0) {
        showNotification("Select the experiments to remove first.", type = "warning")
        return()
      }
      tbl <- experiments()
      affected <- sum(tbl$Comparisons[tbl$`Experiment ID` %in% picked], na.rm = TRUE)
      showModal(modalDialog(
        title = "Remove experiments",
        easyClose = TRUE,
        tags$p("These ", length(picked), " experiment(s) and the ", affected,
               " comparison(s) belonging to them will be removed from ",
               tags$strong(ds$name, .noWS = noWS), ":"),
        tags$ul(lapply(picked, tags$li)),
        tags$p("Their statistics are deleted from the database and their rows from",
               "the metadata, in one step. The files they were submitted with are",
               "kept, so the same experiment can be submitted again."),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(session$ns("removeConfirm"), "Remove", icon = icon("trash"),
                       class = "btn-danger")
        )
      ))
    })

    observeEvent(input$removeConfirm, {
      removeModal()
      ds <- dataset()
      picked <- selectedExperiments()
      res <- removeExperimentsJob(
        ds, config, picked,
        present = as.character(datasetExperiments(data(), ds$name)$`Experiment ID`))
      if (!isTRUE(res$ok)) {
        showNotification(res$message, type = "error", duration = 10)
        return()
      }
      jobId(res$jobId)
      lastRun(NULL)
      showNotification("Removing. This panel will follow the job.",
                       type = "message", duration = 5)
    })

    #### Rebuild ####
    rebuildWorkbook <- reactive({
      up <- input$workbook
      if (is.null(up)) return(NULL)
      readWorkbook(up$datapath[1])
    })

    output$rebuildSummary <- renderUI({
      wb <- rebuildWorkbook()
      ds <- dataset()
      if (is.null(wb)) return(tags$p(tags$em("No workbook uploaded.")))
      if (!isTRUE(wb$ok)) {
        return(tags$div(style = "color: #c62828", tags$ul(lapply(wb$errors, tags$li))))
      }
      present <- if (is.null(ds)) character(0) else
        as.character(datasetExperiments(data(), ds$name)$`Experiment ID`)
      if (!wb$experimentId %in% present) {
        return(tags$p(style = "color: #c62828",
                      "'", wb$experimentId, "' is not in this dataset. Rebuilding ",
                      "republishes an experiment that is already here; use Submit ",
                      "to add a new one."))
      }
      haveResults <- dir.exists(file.path(submitPaths(ds$path)$results,
                                          wb$experimentId, "tables"))
      tagList(
        tags$p(tags$strong("Will rebuild: "), wb$experimentId),
        if (!haveResults) tags$p(style = "color: #c62828",
          "No result tables are kept for this experiment, so there is nothing to ",
          "republish. It was added outside CRAVE, or its submission area has been ",
          "cleaned. Submit it again with its counts matrix instead.") else NULL
      )
    })

    observeEvent(input$rebuild, {
      ds <- dataset()
      wb <- rebuildWorkbook()
      if (is.null(ds) || is.null(wb) || !isTRUE(wb$ok)) {
        showNotification("Upload a valid workbook first.", type = "warning")
        return()
      }
      present <- as.character(datasetExperiments(data(), ds$name)$`Experiment ID`)
      if (!wb$experimentId %in% present) {
        showNotification(paste0("'", wb$experimentId, "' is not in this dataset."),
                         type = "warning", duration = 8)
        return()
      }
      staged <- stageSubmission(ds, wb, input$workbook$datapath[1], NULL)
      if (!isTRUE(staged$ok)) {
        showNotification(staged$message, type = "error", duration = 10)
        return()
      }
      res <- rebuildExperimentJob(ds, config, staged$workbook)
      if (!isTRUE(res$ok)) {
        showNotification(res$message, type = "error", duration = 10)
        return()
      }
      jobId(res$jobId)
      lastRun(NULL)
      showNotification("Rebuilding from the retained results.", type = "message",
                       duration = 5)
    })

    #### Metadata ####
    metadataPath <- function(which) {
      ds <- dataset()
      if (is.null(ds)) return(NULL)
      file.path(ds$path, paste0(which, "_metadata.csv.gz"))
    }

    output$downloadExperiments <- downloadHandler(
      filename = function() paste0(dataset()$name %||% "dataset", "-experiments_metadata.csv.gz"),
      content  = function(f) file.copy(metadataPath("experiments"), f, overwrite = TRUE)
    )
    output$downloadComparisons <- downloadHandler(
      filename = function() paste0(dataset()$name %||% "dataset", "-comparisons_metadata.csv.gz"),
      content  = function(f) file.copy(metadataPath("comparisons"), f, overwrite = TRUE)
    )

    observeEvent(input$applyMetadata, {
      ds <- dataset()
      if (is.null(ds)) return()
      res <- replaceMetadataTables(
        ds,
        experimentsFile = if (is.null(input$uploadExperiments)) NULL
                          else input$uploadExperiments$datapath[1],
        comparisonsFile = if (is.null(input$uploadComparisons)) NULL
                          else input$uploadComparisons$datapath[1]
      )
      notice(res)
      if (isTRUE(res$ok)) {
        showNotification(res$message, type = "message", duration = 8)
        recheck(recheck() + 1)
        onChanged()
      } else {
        showNotification("The metadata was not replaced; see the panel.",
                         type = "error", duration = 8)
      }
    })

    output$metadataNotice <- renderUI({
      res <- notice()
      if (is.null(res)) return(NULL)
      noticeBanner(res$message, if (isTRUE(res$ok)) "info" else "warning")
    })

    #### Following a job ####
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

    observeEvent(status(), {
      st <- status()
      if (is.null(st) || identical(st$state, "running")) return()
      if (identical(lastRun(), st$jobId)) return()
      lastRun(st$jobId)
      if (identical(st$state, "ok")) {
        showNotification("The dataset was changed. Reloading.", type = "message",
                         duration = 6)
        recheck(recheck() + 1)
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
        running = tags$p(icon("spinner"), tags$strong(" Running"), " — ", st$jobId),
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
