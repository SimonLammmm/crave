#### Guides module ####
# The guide library lives in its own SQLite file and can be large, so neither the
# connection nor the five DISTINCT queries that populate the pickers run at
# start-up. They are opened on the first visit to the tab, by the first user to get
# there, and then shared by every session for the life of the process.

.guideCache <- new.env(parent = emptyenv())

#' Open (once) and return the guide library connection and its picker choices.
guideLibraryHandle <- function(path) {
  if (is.null(path)) return(NULL)
  key <- paste0("gl:", path)
  if (!is.null(.guideCache[[key]])) return(.guideCache[[key]])

  handle <- tryCatch({
    con <- DBI::dbConnect(RSQLite::SQLite(), path, flags = RSQLite::SQLITE_RO)
    distinct1 <- function(col) {
      sort(DBI::dbGetQuery(con, paste0("SELECT DISTINCT [", col, "] AS v FROM guideLibrary"))$v)
    }
    list(
      con = con,
      choices = list(
        targets  = distinct1("Target"),
        assembly = distinct1("Assembly"),
        pam      = distinct1("PAM"),
        library  = distinct1("Library"),
        organism = distinct1("Organism")
      )
    )
  }, error = function(e) {
    log_error(paste0("Could not open the guide library at ", path, ": ", conditionMessage(e)))
    NULL
  })

  .guideCache[[key]] <- handle
  handle
}

GUIDE_GENE_CLASSES <- c("Protein-coding", "Non-coding RNA", "Pseudogene",
                        "Base edit", "Array", "Other", "Non-targeting")
GUIDE_CHROMOSOMES  <- c(as.character(1:22), "X", "Y")
GUIDE_CHEMISTRY    <- c("Knockout", "Interference", "Activation",
                        "Cytosine base editor", "Adenine base editor")

guidesUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Guides", icon = icon("worm"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        tags$h4(tags$strong("Target controls")),
        selectizeInput(ns("targets"), "Targets", choices = NULL, multiple = TRUE),
        checkboxGroupInput(ns("geneclasses"), "Gene classes",
                           choices = GUIDE_GENE_CLASSES, inline = TRUE),
        checkboxGroupInput(ns("chromosomes"), "Chromosomes",
                           choices = GUIDE_CHROMOSOMES, inline = TRUE),
        tags$hr(),
        tags$h4(tags$strong("Library controls")),
        checkboxGroupInput(ns("chemistry"), "CRISPR chemistry", choices = GUIDE_CHEMISTRY),
        selectizeInput(ns("assembly"), "Assembly", choices = NULL, multiple = TRUE),
        selectizeInput(ns("pam"),      "PAM",      choices = NULL, multiple = TRUE),
        selectizeInput(ns("library"),  "Library",  choices = NULL, multiple = TRUE),
        selectizeInput(ns("organism"), "Organism", choices = NULL, multiple = TRUE),
        tags$br(), tags$br(),
        actionButton(ns("example"), "Example", icon = icon("fire")),
        actionButton(ns("reset"), "Reset form", icon = icon("arrows-rotate"), class = "btn-danger"),
        tags$hr(),
        tags$p(tags$strong("Help")),
        tags$p("View CRISPR sgRNA sequences used in libraries featured on this portal."),
        tags$p("Use the controls to filter for targets and library designs.")
      ),
      mainPanel(
        tags$h4("Guides"),
        tags$p("View CRISPR sgRNA sequences used in libraries featured on this portal."),
        tags$br(),
        actionButton(ns("submit"), "Guides", icon = icon("worm"), class = "btn-success"),
        tags$br(), tags$br(),
        DTOutput(ns("table")),
        tags$br(),
        downloadButton(ns("download"), "Download", icon = icon("download"))
      )
    )
  )
}

guidesServer <- function(id, data, activeTab) {
  moduleServer(id, function(input, output, session) {

    pickerIds <- c("targets", "assembly", "pam", "library", "organism")
    ready <- reactiveVal(FALSE)

    # Populate the pickers the first time the tab is opened.
    observeEvent(activeTab(), {
      if (!identical(activeTab(), "Guides") || isTRUE(ready())) return()
      handle <- guideLibraryHandle(data()$files$guideLibrary)
      if (is.null(handle)) {
        showNotification("The guide library could not be opened.", type = "error")
        return()
      }
      withProgress(message = "Loading the guide library...", value = NULL, {
        for (p in pickerIds) {
          updateSelectizeInput(session, p, choices = handle$choices[[p]],
                               selected = character(0), server = TRUE,
                               options = list(maxOptions = CRAVE_MAX_SELECTIZE_SERVER))
        }
      })
      ready(TRUE)
    })

    resetForm <- function() {
      for (p in pickerIds) updateSelectizeInput(session, p, selected = character(0))
      for (g in c("geneclasses", "chromosomes", "chemistry")) {
        updateCheckboxGroupInput(session, g, selected = character(0))
      }
    }

    observeEvent(input$reset, resetForm())

    observeEvent(input$example, {
      resetForm()
      updateCheckboxGroupInput(session, "geneclasses", selected = "Protein-coding")
      updateCheckboxGroupInput(session, "chromosomes", selected = "10")
      updateCheckboxGroupInput(session, "chemistry",   selected = "Knockout")
      updateSelectizeInput(session, "organism", selected = "Human")
    })

    result <- eventReactive(input$submit, {
      handle <- guideLibraryHandle(data()$files$guideLibrary)
      if (is.null(handle)) return(tibble(` ` = "The guide library is not available."))
      withProgress(message = "Searching guides...", value = NULL, {
        runGuideLibrary(handle$con, list(
          targets     = input$targets,
          geneclasses = input$geneclasses,
          chromosomes = input$chromosomes,
          chemistry   = input$chemistry,
          assembly    = input$assembly,
          pam         = input$pam,
          library     = input$library,
          organism    = input$organism
        ))
      })
    }, ignoreInit = TRUE)

    output$table <- renderDT(
      {
        if (input$submit == 0) tibble(` ` = "Fill the form in the sidebar and then click the Guides button.")
        else result()
      },
      filter = "top", server = TRUE, escape = FALSE,
      options = list(pageLength = 100, scrollX = TRUE)
    )

    output$download <- downloadHandler(
      filename = "guides.csv.gz",
      content  = function(f) data.table::fwrite(result(), f, sep = ",")
    )

    stateIds <- c(pickerIds, "geneclasses", "chromosomes", "chemistry")

    list(
      state = function() captureInputs(input, stateIds),
      restore = function(st) {
        if (!is.list(st)) return(invisible(NULL))
        # These pickers are server = TRUE, so the choices have to be resent with the
        # selection or selectize discards a value it has not paged in.
        handle <- guideLibraryHandle(data()$files$guideLibrary)
        for (p in intersect(pickerIds, names(st))) {
          if (!is.null(handle)) {
            updateSelectizeInput(session, p, choices = handle$choices[[p]],
                                 selected = st[[p]], server = TRUE)
          } else {
            updateSelectizeInput(session, p, selected = st[[p]])
          }
        }
        for (g in intersect(c("geneclasses", "chromosomes", "chemistry"), names(st))) {
          updateCheckboxGroupInput(session, g, selected = st[[g]])
        }
        invisible(NULL)
      }
    )
  })
}
