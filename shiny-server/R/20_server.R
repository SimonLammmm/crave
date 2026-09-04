#### Top-level server ####
# The wiring layer. It owns four things and nothing else: the shared data handle, the
# message bus between tabs, the two pseudo-tabs (Refresh data and Save/Load), and the
# plotly brush, which has to be read at the root session because plotly does not
# namespace its event inputs.

#' Build the server function.
#'
#' @param config Deployment configuration.
#' @param dataStore An app-level reactiveVal holding the app data object. Shared
#'   across sessions, so one user refreshing the data refreshes it for everyone.
craveServer <- function(config, dataStore) {
  function(input, output, session) {

    if (!is.null(config$privacy_toast)) {
      showNotification(ui = config$privacy_toast, duration = 8, type = "message")
    }

    data <- reactive(dataStore())
    caps <- isolate(dataStore()$caps)

    # Cross-tab message bus. Explore and Correlate each publish their current
    # selections here so the other can copy them, and the brush handler publishes
    # into it too.
    bus <- reactiveValues(
      exploreGenes = character(0), exploreContrasts = character(0),
      correlateGenes = character(0), correlateContrasts = character(0),
      brushTable = NULL, brush = NULL
    )

    navigate <- function(tab) updateTabsetPanel(session, "main", selected = tab)

    #### Tabs ####
    homeServer("home", caps, navigate)
    correlateHandle <- correlateServer("correlate", data, bus)
    exploreHandle   <- exploreServer("explore", data, bus)

    activeTab <- reactive(input$main)

    guidesHandle <- if (isTRUE(caps$guides)) guidesServer("guides", data, activeTab) else NULL
    if (isTRUE(caps$libraries)) librariesServer("libraries", data)
    if (isTRUE(caps$ontology))  ontologyServer("ontology", data)
    exorciseHandle <- if (isTRUE(caps$exorcise)) exorciseServer("exorcise", config) else NULL

    #### Plot brushing ####
    # plotly reports selections through a fixed, un-namespaced input id, so it can
    # only be read here. The module that drew the plot has told us, via the bus,
    # which coordinates map to which genes.
    observeEvent(input$`plotly_brushed-A`, {
      published <- bus$brushTable
      if (is.null(published) || is.null(published$table)) return()
      raw <- input$`plotly_brushed-A`
      xs <- suppressWarnings(as.numeric(unlist(strsplit(
        sub("^.+x\\\":\\[(.+?)\\],\\\"y.+$", "\\1", raw), ","))))
      ys <- suppressWarnings(as.numeric(unlist(strsplit(
        sub("^.+y\\\":\\[(.+?)\\].+$", "\\1", raw), ","))))
      if (length(xs) < 2 || length(ys) < 2 || anyNA(c(xs[1:2], ys[1:2]))) return()
      tbl <- published$table
      inside <- !is.na(tbl$x) & !is.na(tbl$y) &
        tbl$x >= min(xs) & tbl$x <= max(xs) &
        tbl$y >= min(ys) & tbl$y <= max(ys)
      genes <- unique(tbl$Gene[inside])
      if (length(genes) == 0) return()
      bus$brush <- list(target = published$target, genes = genes, at = Sys.time())
    })

    #### Refresh data and Save/Load pseudo-tabs ####
    # Both are navbar entries that act as buttons: they bounce the user back to
    # where they were and then do their job.
    lastTab <- reactiveVal("Home")
    pseudoTabs <- c("Refresh data", "Save/Load")

    observeEvent(input$main, {
      if (!(input$main %in% pseudoTabs)) lastTab(input$main)
    })

    saveLoadOpened <- reactiveVal(0)

    observeEvent(input$main, {
      if (!(input$main %in% pseudoTabs)) return()
      navigate(lastTab())
      if (identical(input$main, "Save/Load")) {
        saveLoadOpened(saveLoadOpened() + 1)
      } else {
        withProgress(message = "Reloading datasets...", value = NULL, {
          old <- dataStore()
          fresh <- tryCatch(
            loadCraveData(config$datasets, config$exorcise_root, config$exorcise_docker,
                          create_indexes = isTRUE(config$create_indexes)),
            error = function(e) {
              log_error(paste0("Data refresh failed: ", conditionMessage(e)))
              NULL
            }
          )
          if (is.null(fresh) || isFALSE(fresh$caps$explore)) {
            showNotification("The data could not be reloaded. The previous data are still in use.",
                             type = "error", duration = 8)
          } else {
            dataStore(fresh)
            closeCraveData(old)
            showNotification("Datasets reloaded.", type = "message", duration = 4)
          }
        })
      }
    })

    modules <- Filter(Negate(is.null), list(
      explore   = exploreHandle,
      correlate = correlateHandle,
      guides    = guidesHandle,
      exorcise  = exorciseHandle
    ))
    saveLoadServer("saveload", modules, saveLoadOpened)

    #### Logging ####
    # Logging every input change installs an observer per input and writes a line
    # on every keystroke, so it is off unless the deployer asks for it.
    if (isTRUE(config$enable_input_logging)) log_shiny_input_changes(input)
  }
}
