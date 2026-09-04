#### Reusable "Customise plot" modal module ####
# Explore and Correlate share this. The fields are declared as data (see
# EXPLORE_CUSTOMISE_FIELDS and CORRELATE_CUSTOMISE_FIELDS) and the module builds the
# modal, mirrors edits into a store that survives the modal closing, and exposes a
# fully populated list of values — so an analysis never has to cope with a field
# the user has not visited yet.

#' Describe one customisation field.
#'
#' @param type One of "numeric", "checkbox", "colour", "select".
cfField <- function(id, type, label, choices = NULL) {
  list(id = id, type = type, label = label, choices = choices)
}

EXPLORE_CUSTOMISE_FIELDS <- list(
  cfField("height",       "numeric",  "Plot height (pixels)"),
  cfField("width",        "numeric",  "Plot width (pixels)"),
  cfField("point_colour", "colour",   "Point colour"),
  cfField("goi_colour",   "colour",   "Genes of interest colour"),
  cfField("x_axis_log",   "select",   "x-axis scale", c("Linear", "Logarithmic", "Automatic")),
  cfField("y_axis_log",   "select",   "y-axis scale", c("Linear", "Logarithmic", "Automatic")),
  cfField("y_equals_x",   "checkbox", "Equal x and y axes (biplot only)")
)

CORRELATE_CUSTOMISE_FIELDS <- list(
  cfField("height",               "numeric",  "Plot height (pixels)"),
  cfField("width",                "numeric",  "Plot width (pixels)"),
  cfField("width_auto",           "checkbox", "Fit plot to window width"),
  cfField("point_colour",         "colour",   "Point colour"),
  cfField("genequery_fillscheme", "select",   "Colour violins by (gene query only)",
          c("Treatment", "Knockout", "Cell line", "Library")),
  cfField("heatmap_high",         "colour",   "Positive correlation colour (heatmap only)"),
  cfField("heatmap_mid",          "colour",   "Zero correlation colour (heatmap only)"),
  cfField("heatmap_low",          "colour",   "Negative correlation colour (heatmap only)")
)

#' Button that opens the customisation modal.
customiseUI <- function(id, label = "Customise plot") {
  ns <- NS(id)
  actionButton(ns("open"), label, icon = icon("brush"))
}

#' Server for the customisation modal.
#'
#' @param fields List of cfField() descriptions.
#' @param defaults Named list of default values, keyed by field id.
#' @return list(values = reactive(named list), restore = function(list))
customiseServer <- function(id, fields, defaults) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    store <- do.call(reactiveValues, defaults)

    buildInput <- function(f, value) {
      switch(
        f$type,
        numeric  = numericInput(ns(f$id), f$label, value = value, min = 0, max = Inf),
        checkbox = checkboxInput(ns(f$id), f$label, value = isTRUE(value)),
        colour   = colourpicker::colourInput(ns(f$id), f$label, value = value),
        select   = selectizeInput(ns(f$id), f$label, choices = f$choices, selected = value),
        NULL
      )
    }

    showCustomise <- function() {
      showModal(modalDialog(
        title = "Customise plot",
        easyClose = TRUE,
        footer = NULL,
        lapply(fields, function(f) buildInput(f, store[[f$id]])),
        tags$hr(),
        actionButton(ns("reset"), "Restore defaults", icon = icon("arrows-rotate")),
        modalButton("Save options", icon = icon("floppy-disk"))
      ))
      # Width is meaningless when the plot is told to fill the window.
      # shinyjs applies the module namespace itself, so the id must be bare.
      if ("width_auto" %in% vapply(fields, `[[`, character(1), "id")) {
        if (isTRUE(store$width_auto)) shinyjs::disable("width")
      }
    }

    observeEvent(input$open, showCustomise())

    # Mirror every edit into the store so values survive the modal being closed
    # and re-opened.
    lapply(fields, function(f) {
      observeEvent(input[[f$id]], {
        store[[f$id]] <- input[[f$id]]
      }, ignoreInit = TRUE)
    })

    observeEvent(input$reset, {
      for (nm in names(defaults)) store[[nm]] <- defaults[[nm]]
      removeModal()
      showCustomise()
    })

    observeEvent(input$width_auto, {
      if (isTRUE(input$width_auto)) shinyjs::disable("width")
      else shinyjs::enable("width")
    }, ignoreInit = TRUE)

    list(
      values = reactive(reactiveValuesToList(store)),
      restore = function(state) {
        if (!is.list(state)) return(invisible(NULL))
        for (nm in intersect(names(state), names(defaults))) {
          if (!is.null(state[[nm]])) store[[nm]] <- state[[nm]]
        }
        invisible(NULL)
      }
    )
  })
}
