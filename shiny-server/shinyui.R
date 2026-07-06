#### Shiny UI: preamble ####

# Decide which features are available given the datasets specified
if(canGuideLibrary) {
  jumpGuideLibrary <- actionButton(inputId = "jumpGuideLibrary", label = "Guides", icon = icon("worm"))
  descGuideLibrary <- tags$li(tags$strong("Guides", .noWS = noWS), paste0(", where you can search the contents of sgRNA libraries featured on ", portal_name, "."))
} else {
  jumpGuideLibrary <- ""
  descGuideLibrary <- ""
}
if(canLibraries) {
  jumpLibraries <- actionButton(inputId = "jumpLibraries", label = "Libraries", icon = icon("book-open"))
  descLibraries <- tags$li(tags$strong("Libraries", .noWS = noWS), ", where you can see information on those libraries.")
} else {
  jumpLibraries <- ""
  descLibraries <- ""
}
if(canOntology) {
  jumpOntology <- actionButton(inputId = "jumpOntology", label = "Ontology", icon = icon("sitemap"))
  descOntology <- tags$li(tags$strong("Ontology", .noWS = noWS), paste0(", where you can see the members of the gene classes used on ", portal_name, "."))
} else {
  jumpOntology <- ""
  descOntology <- ""
}
if(canExorcise) {
  jumpExorcise <- actionButton(inputId = "jumpExorcise", label = "Exorcise", icon = icon("ghost"))
  descExorcise <- tags$li(tags$strong("Exorcise", .noWS = noWS), ", where you can reannotate spCas9 guide sequences with targets in any CRISPR chemistry.")
} else {
  jumpExorcise <- ""
  descExorcise <- ""
}

# Decide whether to render a custom logo
# Use the default logo
logo <- HTML(paste0('<a style="text-decoration:none;cursor:default;" class="active" href="#">', portal_name, '</a>'))

# Replace with a supplied logo if specified in config.R and file exists
if(!is.null(portal_logo)) {
  if(file.exists(paste0("www/", portal_logo))) {
    logo <- tags$img(src=portal_logo, height="24px")
  }
}

ui <-
  navbarPage(theme = shinytheme(portal_theme),
             logo,
             windowTitle = portal_name,
             id = "main",
             shinyjs::useShinyjs(),
             
             #### Shiny UI: Homepage ####
             tabPanel("Home", icon = icon("house"),
                      noticebanner_info,
                      noticebanner_warning,
                      noticebanner_emergency,
                      tags$div(style = "text-align: right",
                               tags$p(textOutput("motd"))
                      ),
                      tags$div(
                        style = "text-align: center",
                        tags$p("Welcome to"),
                        h1(tags$strong(portal_name), style = "font-size: 96px"),
                        h4(portal_subtitle),
                        actionButton(inputId = "jumpCorrelate", label = "Correlate", icon = icon("chart-line")),
                        actionButton(inputId = "jumpExplore", label = "Explore", icon = icon("microscope")),
                        jumpGuideLibrary,
                        jumpLibraries,
                        jumpOntology,
                        jumpExorcise,
                        tags$br(),
                        tags$br(),
                        tags$hr()
                      ),
                      tags$div(
                        tags$h4("Information"), 
                        HTML(front_page_html),
                        tags$p("The currently available features are:"),
                        tags$ul(
                          tags$li(tags$strong("Correlate", .noWS = noWS), ", where you can inspect genes across many screens in the database."),
                          tags$li(tags$strong("Explore", .noWS = noWS), ", where you can visualise particular screens of interest."),
                          descGuideLibrary,
                          descLibraries,
                          descOntology,
                          descExorcise
                        ),
                        tags$br(),
                        tags$h4("Application"),
                        tags$p(paste0(portal_subtitle, ", this version ", ver, ", last updated ", updated, ".")),
                        #### Deployer-specific content (credits, logos, contact) lives in config.R: portal_footer ####
                        portal_footer,
                        tags$br()
                      )
             ),
             
             #### Shiny UI: Correlate ####
             tabPanel("Correlate", icon = icon("chart-line"),
                      sidebarLayout(
                        correlateSidebar,
                        mainPanel(
                          tabsetPanel(
                            genequeryTab,
                            clustergramTab,
                            heatmapTab,
                            networkplotTab,
                            umapTab,
                            enrichmentanalysisTab,
                            ddromeTab,
                            datadownloadTab,
                          )
                        )
                      )
             ),
             
             #### Shiny UI: Explore ####
             tabPanel("Explore", icon = icon("microscope"),
                      sidebarLayout(
                        exploreSidebar,
                        mainPanel(
                          tabsetPanel(
                            experimentSelectTab,
                            comparisonSelectTab,
                            volcanoplotTab,
                            rankplotTab,
                            biplotTab,
                            overlapTab,
                            rocaucTab
                          )
                        )
                      )
             ),
             #### Shiny UI: Guides ####
             guideLibraryTab,
             #### Shiny UI: Libraries ####
             libraryTab,
             #### Shiny UI: Ontology ####
             ontologyTab,
             #### Shiny UI: Exorcise ####
             exorciseTab,
             manageTab,
             analysisTab,
             #### Shiny UI: Refresh data connection ####
             tabPanel("Refresh data", icon = icon("arrows-rotate"),
                      fluidPage(
                        tags$p("Refreshing data connection.\nPlease wait...")
                      )),
             #### Shiny UI: Save/Load ####
             tabPanel("Save/Load", icon = icon("folder-open"),
                      #actionButton(inputId = 'saveload', label = 'Save/Load', icon = icon("folder-open"), class = "btn-success"),
             ),
             #### Shiny UI: Legal notices ####
             tabPanel("Legal notices", icon = icon("gavel"),
                      fluidPage(
                        legalTab
                      )
             ),
             header = fluidPage(busy_start_up(loader = spin_epic("self-building-square", color = "#ff0087"), text = "Almost there...", mode = "auto"), # Wait screen
                                add_busy_bar(color = "#FF0000")) # Busy indicator
  )