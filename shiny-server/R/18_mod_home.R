#### Home page ####

#' A coloured notice banner, or NULL if there is nothing to say.
noticeBanner <- function(content, kind = c("info", "warning", "emergency")) {
  if (is.null(content)) return(NULL)
  kind <- match.arg(kind)
  style <- switch(
    kind,
    info      = list(bg = "#e6f0fc", border = "#0469e3", fg = "black", icon = "circle-info"),
    warning   = list(bg = "#fff9e6", border = "#fec000", fg = "black", icon = "triangle-exclamation"),
    emergency = list(bg = "#FF3127", border = "#FF3127", fg = "white", icon = "circle-exclamation")
  )
  tags$div(
    style = paste0(
      "background-color: ", style$bg, "; padding: 12px; ",
      "border-style: none none none solid; border-color: ", style$border,
      "; border-width: 6px"
    ),
    tags$p(style = paste0("color: ", style$fg), icon(style$icon), content)
  )
}

#' The Home tab.
#'
#' @param caps Capability flags from the app data object.
#' @param config The deployment configuration.
homeUI <- function(id, caps, config) {
  ns <- NS(id)

  # Jump buttons and feature descriptions for the tabs that are actually enabled.
  features <- list(
    list(cap = TRUE,           tab = "Correlate", label = "Correlate", icon = "chart-line",
         desc = ", where you can inspect genes across many screens in the database."),
    list(cap = TRUE,           tab = "Explore",   label = "Explore",   icon = "microscope",
         desc = ", where you can visualise particular screens of interest."),
    list(cap = caps$guides,    tab = "Guides",    label = "Guides",    icon = "worm",
         desc = paste0(", where you can search the contents of sgRNA libraries featured on ",
                       config$portal_name, ".")),
    list(cap = caps$libraries, tab = "Libraries", label = "Libraries", icon = "book-open",
         desc = ", where you can see information on those libraries."),
    list(cap = caps$ontology,  tab = "Ontology",  label = "Ontology",  icon = "sitemap",
         desc = paste0(", where you can see the members of the gene classes used on ",
                       config$portal_name, ".")),
    list(cap = caps$exorcise,  tab = "Exorcise",  label = "Exorcise",  icon = "ghost",
         desc = ", where you can reannotate spCas9 guide sequences with targets in any CRISPR chemistry.")
  )
  enabled <- Filter(function(f) isTRUE(f$cap), features)

  tabPanel(
    "Home", icon = icon("house"),
    noticeBanner(config$notices_info, "info"),
    noticeBanner(config$notices_warning, "warning"),
    noticeBanner(config$notices_emergency, "emergency"),
    tags$div(style = "text-align: right", tags$p(textOutput(ns("motd")))),
    tags$div(
      style = "text-align: center",
      tags$p("Welcome to"),
      h1(tags$strong(config$portal_name), style = "font-size: 96px"),
      h4(config$portal_subtitle),
      lapply(enabled, function(f) {
        actionButton(ns(paste0("jump_", f$tab)), f$label, icon = icon(f$icon))
      }),
      tags$br(), tags$br(), tags$hr()
    ),
    tags$div(
      tags$h4("Information"),
      HTML(config$front_page_html),
      tags$p("The currently available features are:"),
      tags$ul(lapply(enabled, function(f) {
        tags$li(tags$strong(f$label, .noWS = noWS), f$desc)
      })),
      tags$br(),
      tags$h4("Application"),
      tags$p(paste0(config$portal_subtitle, ", this version ", CRAVE_VERSION,
                    ", last updated ", CRAVE_UPDATED, ".")),
      config$portal_footer,
      tags$br()
    )
  )
}

#' Server for the Home tab.
#'
#' @param navigate Function(tabName) that switches the top-level tabset.
homeServer <- function(id, caps, navigate) {
  moduleServer(id, function(input, output, session) {
    for (tab in c("Correlate", "Explore", "Guides", "Libraries", "Ontology", "Exorcise")) {
      # `local` captures the loop variable; without it every button would jump to
      # the last tab in the list.
      local({
        thisTab <- tab
        observeEvent(input[[paste0("jump_", thisTab)]], navigate(thisTab))
      })
    }

    # One message of the day per hour, shared by everyone connected.
    #
    # Chosen arithmetically, not with set.seed() + sample(). This runs on the tab
    # every visitor lands on, and seeding the global RNG here would make every later
    # random draw in the process deterministic per hour — including the Exorcise
    # task ids that name files in a shared directory.
    output$motd <- renderText({
      slot <- as.integer(as.numeric(Sys.time()) %/% 3600)
      CRAVE_MOTDS[[1L + slot %% length(CRAVE_MOTDS)]]
    })
  })
}
