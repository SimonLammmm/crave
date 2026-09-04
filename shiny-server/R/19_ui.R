#### Top-level UI ####

#' Build the navbar brand: a supplied logo if one exists, otherwise the name.
portalBrand <- function(config) {
  if (!is.null(config$portal_logo) &&
      file.exists(file.path("www", config$portal_logo))) {
    return(tags$img(src = config$portal_logo, height = "24px"))
  }
  HTML(paste0(
    '<a style="text-decoration:none;cursor:default;" class="active" href="#">',
    htmltools::htmlEscape(config$portal_name), '</a>'
  ))
}

#' Assemble the whole UI.
#'
#' Tabs whose data are not present are omitted rather than rendered empty.
craveUI <- function(caps, config) {
  navbarPage(
    theme = shinytheme(config$portal_theme),
    portalBrand(config),
    windowTitle = config$portal_name,
    id = "main",
    shinyjs::useShinyjs(),

    homeUI("home", caps, config),
    correlateUI("correlate", caps),
    exploreUI("explore"),
    if (isTRUE(caps$guides))    guidesUI("guides")       else NULL,
    if (isTRUE(caps$libraries)) librariesUI("libraries") else NULL,
    if (isTRUE(caps$ontology))  ontologyUI("ontology")   else NULL,
    if (isTRUE(caps$exorcise))  exorciseUI("exorcise")   else NULL,

    tabPanel("Refresh data", icon = icon("arrows-rotate"),
             fluidPage(tags$p("Refreshing the data connection. Please wait..."))),
    tabPanel("Save/Load", icon = icon("folder-open")),
    legalUI(config$legal_text),

    header = fluidPage(
      busy_start_up(loader = spin_epic("self-building-square", color = "#ff0087"),
                    text = "Almost there...", mode = "auto"),
      add_busy_bar(color = "#FF0000")
    )
  )
}

#' A minimal UI used when CRAVE cannot start.
craveErrorUI <- function(message) {
  fluidPage(
    tags$head(tags$title("CRAVE: startup error")),
    tags$div(
      style = "max-width: 40em; margin: 4em auto; font-family: sans-serif",
      tags$h3("CRAVE could not start"),
      tags$p(HTML(message))
    )
  )
}
