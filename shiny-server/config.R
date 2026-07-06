datasets <- list(
  list(
    name             = "CRAVE",
    path             = "~/Desktop/crave-data/crave-20260619/",
    guides           = TRUE,
    citation_from_id = TRUE
  ),
  list(
    name             = "DDRcs",
    path             = "~/Desktop/crave-data/ddrcs-20260718/",
    guides           = FALSE,
    citation_from_id = FALSE
  )
)

exorcise_root     <- ""
exorcise_docker   <- "simonlammmm/exorcise:1.5.3.14_arm64"
exorcise_platform <- "linux/amd64"
enable_analysis   <- TRUE
portal_theme      <- "cerulean"

portal_name     <- "CRAVE"
portal_subtitle <- "CRISPR results app for visualisation and exploration"
front_page_html <- "<p>CRAVE enables the analysis of hits within and between CRISPR screens.</p>"

portal_footer <- tagList(
  tags$h4("Contact"),
  tags$p("sl681@cam.ac.uk")
)

privacy_toast <- tags$div(
  tags$strong("Privacy information"),
  tags$p("By continuing to use the site, you agree to the Privacy Policy in the Legal tab.")
)

legal_text <- tags$div(
  tags$h1("Privacy policy"),
  tags$p("Placeholder — replace before public deployment."),
  tags$hr(),
  tags$h1("Terms and conditions"),
  tags$p("Placeholder — replace before public deployment."),
  tags$hr(),
  tags$h1("Disclaimer"),
  tags$p("Placeholder — replace before public deployment.")
)

## Metadata editing — TRUE only on trusted/internal instances (writes to disk)
enable_editing <- TRUE
