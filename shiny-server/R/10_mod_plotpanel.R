#### Reusable plot panel module ####
# The structure every analysis tab shares: a submit button, a plot, a data table and
# download handlers for both. An analysis tab supplies its controls, its prose, and a
# function returning list(p, plotdata, height, renderer, brush).

#' UI for a plot panel.
#'
#' @param id Module id.
#' @param title Heading shown above the panel.
#' @param description Character vector; one <p> per element.
#' @param controls Optional UI for the tab's own inputs, shown above the button.
#' @param submitLabel,submitIcon Label and icon for the run button.
#' @param tableFirst TRUE to place the table above the plot.
plotPanelUI <- function(id, title, description = character(0), controls = NULL,
                        submitLabel = "Run", submitIcon = "play",
                        tableFirst = FALSE) {
  ns <- NS(id)
  plotBlock  <- uiOutput(ns("plot"))
  tableBlock <- DTOutput(ns("table"))
  tagList(
    tags$br(),
    tags$h4(title),
    lapply(description, function(d) tags$p(HTML(d))),
    tags$br(),
    controls,
    actionButton(ns("submit"), submitLabel, icon = icon(submitIcon), class = "btn-success"),
    tags$br(),
    if (tableFirst) tagList(tableBlock, tags$br(), plotBlock)
    else            tagList(plotBlock,  tags$br(), tableBlock),
    tags$br(),
    downloadButton(ns("downloadPlot"), "Download plot", icon = icon("download")),
    downloadButton(ns("downloadData"), "Download data", icon = icon("download"))
  )
}

#' Server for a plot panel.
#'
#' @param id Module id.
#' @param compute Zero-argument function returning the analysis result. Called in
#'   a reactive context, so it may read inputs from the enclosing module.
#' @param busyMessage Progress message shown while `compute` runs.
#' @param onResult Optional function(result) called after each successful run;
#'   used to publish brushable coordinates back to the parent.
#' @return The result reactive.
plotPanelServer <- function(id, compute, busyMessage = "Working...", onResult = NULL) {
  moduleServer(id, function(input, output, session) {

    result <- eventReactive(input$submit, {
      withProgress(message = busyMessage, value = NULL, {
        tryCatch(compute(), error = function(e) {
          log_error(paste0("Analysis '", id, "' failed: ", conditionMessage(e)))
          c(msgResult(paste0("Something went wrong running this analysis.\n",
                             conditionMessage(e))),
            list(renderer = "plotly"))
        })
      })
    }, ignoreInit = TRUE)

    output$plot <- renderUI({
      res <- result()
      renderer <- res$renderer %||% "plotly"
      if (identical(renderer, "none")) return(NULL)
      p <- switch(
        renderer,
        plotly     = renderPlotly({ res$p }),
        plot       = renderPlot({ res$p }),
        visNetwork = renderVisNetwork({ res$p }),
        renderPlotly({ res$p })
      )
      attr(p, "outputArgs") <- list(height = res$height %||% 900)
      p
    })

    output$table <- renderDT(
      {
        res <- result()
        tbl <- res$plotdata %||% tibble()
        if (nrow(tbl) == 0) tibble(` ` = "No data.") else tbl
      },
      server = TRUE, filter = "top", escape = FALSE,
      options = list(scrollX = TRUE, pageLength = 25)
    )

    output$downloadData <- downloadHandler(
      filename = function() paste0(id, ".csv"),
      content  = function(f) data.table::fwrite(result()$plotdata %||% tibble(), f, sep = ",")
    )

    output$downloadPlot <- downloadHandler(
      filename = function() {
        switch(result()$renderer %||% "plotly",
               plot = paste0(id, ".png"),
               none = paste0(id, ".txt"),
               paste0(id, ".html"))
      },
      content = function(f) {
        res <- result()
        switch(
          res$renderer %||% "plotly",
          plot = ggplot2::ggsave(f, plot = res$p, width = 12, height = 8, dpi = 150),
          none = writeLines("This analysis produces a table only.", f),
          htmlwidgets::saveWidget(res$p, f, selfcontained = TRUE)
        )
      }
    )

    if (!is.null(onResult)) {
      observeEvent(result(), onResult(result()))
    }

    result
  })
}

#' Lay a list of inputs out in a two-column table.
inputGrid <- function(..., width = "50%", ncol = 2) {
  items <- list(...)
  items <- items[!vapply(items, is.null, logical(1))]
  if (length(items) == 0) return(NULL)
  rows <- split(items, ceiling(seq_along(items) / ncol))
  tags$table(
    style = "width: 100%",
    lapply(rows, function(cells) {
      tags$tr(
        style = "vertical-align: top",
        lapply(cells, function(cell) {
          tags$td(style = paste0("padding: 6px; width: ", width), cell)
        })
      )
    })
  )
}
