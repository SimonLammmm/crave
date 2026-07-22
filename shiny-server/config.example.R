#### CRAVE deployment config ####
## This is your local working config (dev/Mac paths). Keep it out of git.
## For the public repo, copy to config.example.R with /data/... paths.

## Datasets — one entry per database directory
# `name` is the display name in CRAVE.
# `path` is the absolute path of the directory containing CRAVE format datasets
# containing a database.db file. If using docker compose, this is the absolute path
# in the Docker container, specified in docker-compose.yml
# Set `guides` = TRUE if you also have a guideLibrary.db file.
# Most users will leave this FALSE.
# Set `citation_from_id` = TRUE if you want to dynamically generate citations
# from experiment IDs.
# Most users will leave this FALSE.

datasets <- list(
  list(
    name             = "Dataset 1",
    path             = "/data/dataset-1",
    guides           = FALSE,
    citation_from_id = FALSE
  )
  # Add additional entries if you have additional database directories
  #,
  #list(
  #  name             = "Dataset 2",
  #  path             = "/data/dataset-2",
  #  guides           = FALSE,
  #  citation_from_id = FALSE
  #)
)

## Exorcise location - should point to the directory containing .2bit files
# This should be an absolute path. If using docker compose, this should
# be the absolute path in the Docker container, specified in docker-compose.yml
# For example, if your bind mount is "/opt/exorcise/:/data/exorcise/",
# then this variable should be "/data/exorcise/"
exorcise_root <- "/data/exorcise/"

## Exorcise host location (docker/docker compose only)
# If using docker/docker compose, this variable needs to be the absolute path
# of the Exorcise location on the host. For example, if your bind mount is
# "/opt/exorcise/:data/exorcise", then this variable should be "/opt/exorcise/"
# If not using docker/docker compose, then this is the same as exorcise_root
exorcise_hostroot <- exorcise_root

## Exorcise Docker image name - docker pull simonlammmm/exorcise
# If using, obtain Exorcise by running `docker pull simonlammmm/exorcise:latest`
# This variable should be simonlammmm/exorcise:latest unless you pulled a different tag
exorcise_docker <- "simonlammmm/exorcise:latest"

## Look and feel
# Choose a shinytheme. Possible values:
# "cerulean", "cosmo", "cyborg", "darkly", "flatly", "journal", "lumen", "paper",
# "readable", "sandstone", "simplex", "slate", "spacelab", "superhero", "united", "yeti"
portal_theme <- "cerulean"

## Branding
# Optional logo to be displayed in the top bar. Will be rendered at 24px height
portal_logo     <- NULL
# Branding text to be shown on the front page
portal_name     <- "CRAVE"
portal_subtitle <- "CRISPR results app for visualisation and exploration"
front_page_html <- "<p>CRAVE enables the analysis of hits within and between CRISPR screens.</p>"

## Front-page footer (credits, logos, contact) — rendered at the bottom of the Home page
portal_footer <- tagList(
  tags$h4("Abstract"),
  tags$p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")
)

## Privacy toast. This appears when the portal is accessed
privacy_toast <- {
  tags$strong("Privacy information")
tags$div(
  tags$p("Your privacy is important to us. We may collect usage information for:"),
  tags$ul(
    tags$li("necessary purposes,"),
    tags$li("monitoring purposes, and"),
    tags$li("analytics purposes.")
  ),
  tags$p("By continuing to use the site, you agree to the Privacy Policy which can be found in the Legal tab."))
}

## Temporary notices to be shown as banners on the front page
notices_info <- NULL
notices_warning <- NULL
notices_emergency <- NULL

## Legal text
legal_text <- 
  tags$div(
    tags$h1("Privacy policy"),
    tags$p("Lorem ipsum dolor sit amet consectetur adipiscing elit. Consectetur adipiscing elit quisque faucibus ex sapien vitae. Ex sapien vitae pellentesque sem placerat in id. Placerat in id cursus mi pretium tellus duis. Pretium tellus duis convallis tempus leo eu aenean."),
    tags$hr(),
    tags$h1("Terms and conditions"),
    tags$p("Lorem ipsum dolor sit amet consectetur adipiscing elit. Consectetur adipiscing elit quisque faucibus ex sapien vitae. Ex sapien vitae pellentesque sem placerat in id. Placerat in id cursus mi pretium tellus duis. Pretium tellus duis convallis tempus leo eu aenean."),
    tags$hr(),
    tags$h1("Disclaimer"),
    tags$p("Lorem ipsum dolor sit amet consectetur adipiscing elit. Consectetur adipiscing elit quisque faucibus ex sapien vitae. Ex sapien vitae pellentesque sem placerat in id. Placerat in id cursus mi pretium tellus duis. Pretium tellus duis convallis tempus leo eu aenean."),
  )
