#### Simple table tabs: Libraries, Ontology ####
# Declarative outputs inside a tabPanel, so Shiny suspends them until the tab is
# first viewed and neither table costs anything at start-up.

librariesUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Libraries", icon = icon("book-open"),
    fluidPage(
      tags$h4("Libraries"),
      tags$p("View details of CRISPR libraries featured on this portal."),
      tags$br(),
      DTOutput(ns("table"))
    )
  )
}

librariesServer <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    libraries <- reactive({
      f <- data()$files$libraries
      if (is.null(f)) return(tibble(` ` = "No library metadata is loaded for this dataset."))
      raw <- readTableGlob(f)
      if (is.null(raw)) return(tibble(` ` = "The library metadata file could not be read."))

      needed <- c("Library", "Guide length", "Number of guides", "Number of targets",
                  "Species", "CRISPR chemistry", "PAM", "Citation", "Lab",
                  "Source link", "Reference")
      missing <- setdiff(needed, names(raw))
      if (length(missing)) {
        return(tibble(` ` = paste0(
          "The library metadata file is missing the column(s): ",
          paste(missing, collapse = ", "), ".")))
      }

      # Wrap a label in a link when the URL looks like one. Written as plain
      # vector operations rather than a helper called inside transmute(), which
      # relied on dplyr's data mask reaching into the helper's promises.
      linkify <- function(url, label) {
        ok <- grepl("^http", url)
        out <- as.character(label)
        out[ok] <- paste0("<a href=", url[ok], " target=\"_blank\">", label[ok], "</a>")
        out
      }
      libraryLabel <- linkify(raw$`Source link`, raw$Library)
      noSource <- !grepl("^http", raw$`Source link`)
      libraryLabel[noSource] <- linkify(raw$Reference, raw$Library)[noSource]

      tibble(
        Library             = libraryLabel,
        `Guide length`      = raw$`Guide length`,
        `Number of guides`  = raw$`Number of guides`,
        `Number of targets` = raw$`Number of targets`,
        Species             = raw$Species,
        `CRISPR chemistry`  = raw$`CRISPR chemistry`,
        PAM                 = raw$PAM,
        Citation            = linkify(raw$Reference, raw$Citation),
        Lab                 = raw$Lab
      )
    })

    output$table <- renderDT(
      libraries(), filter = "top", server = TRUE, escape = FALSE,
      options = list(pageLength = 100, scrollX = TRUE)
    )
  })
}

ontologyUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Ontology", icon = icon("sitemap"),
    fluidPage(
      tags$h4("Ontology"),
      tags$p("View members of gene classes."),
      tags$br(),
      DTOutput(ns("table"))
    )
  )
}

ontologyServer <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    output$table <- renderDT({
      ont <- data()$ontology
      if (is.null(ont) || nrow(ont) == 0) {
        tibble(` ` = "No ontology is loaded for this dataset.")
      } else {
        ont %>% mutate(class = factor(class))
      }
    }, filter = "top", server = TRUE, escape = FALSE,
       options = list(pageLength = 100, scrollX = TRUE))
  })
}

#### Legal ####

legalUI <- function(text) {
  tabPanel("Legal notices", icon = icon("gavel"), fluidPage(text))
}
