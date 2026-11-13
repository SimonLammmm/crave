#### CRAVE deployment config ####
## Copy this file to config.R and edit it. Every setting below is optional except
## `datasets`; anything you leave out falls back to the default shown in the
## comment. CRAVE validates this file at start-up and, if something is wrong,
## serves a page explaining what rather than failing silently.

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

# Set `hostpath` if you are using Docker AND you want the Submit or Manage tabs.
# Those tabs run the Exorcise container against the dataset, and CRAVE asks the
# Docker daemon on the host to do it, so the daemon resolves the bind mount on the
# host rather than inside CRAVE's container. `hostpath` is the same directory as
# `path`, named the way the host names it: with a bind mount of
# "/opt/crave/dataset-1:/data/dataset-1", `path` is the second and `hostpath` the
# first. Running CRAVE natively, leave it out; the two are the same.

datasets <- list(
  list(
    name             = "Dataset 1",
    path             = "/data/dataset-1",
    #hostpath        = "/opt/crave/dataset-1",
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

## Features
# Three parts of CRAVE a deployment can withhold.
#
# Submit lets a visitor upload a counts matrix and an experimental design
# workbook, run MAGeCK, DrugZ and Chronos over them, and add the results to a
# dataset. Manage lets them relabel a dataset's metadata, rebuild an experiment
# from a corrected workbook, and remove experiments. Both write to the dataset
# directory, and neither asks who is asking -- CRAVE has no idea who its visitors
# are. Leave them FALSE on any deployment the public can reach, and put an
# authenticating proxy in front of one where you turn them on.
#
# Both also need `exorcise_docker` to name an image CRAVE can run, and a dataset
# directory that is writable. If either is missing the tab stays hidden and the
# reason is in the startup log.
enable_submit <- FALSE
enable_manage <- FALSE

# The Pendragonator sub-tab of Correlate, which ranks the strongest responding
# genes across the screens a filter selects. Switch it off to withhold it.
enable_pendragonator <- TRUE

## Server options
# Largest file a user may upload, in megabytes. Exorcise, Save/Load and -- where
# it is enabled -- Submit all use this, and a counts matrix is the big one: a
# genome-wide screen with a dozen samples runs to tens of megabytes.
max_upload_mb <- 100

# Host and port are deliberately NOT set here. Whoever launches the app owns them:
# RStudio's "Run App" chooses a port and expects the app to use it, and
# shiny-server and ShinyProxy assign one and pass it in. Docker passes them in its
# CMD. If you need to pin them for a bare Rscript deployment, do it at the call
# site instead:
#   Rscript -e "shiny::runApp('shiny-server', host = '0.0.0.0', port = 3838)"

## Performance
# On first start CRAVE adds two indexes to each dataset's `stat` table, which turn
# the gene-led queries behind Correlate and Pendragonator from full table scans
# into index lookups. Building them sorts the whole table, so on a very large
# dataset the first start takes several minutes; every start after that is fast.
# Set FALSE to skip it (queries will be slower), or leave TRUE and let it run once.
# Ignored anyway if the dataset directory is read-only.
create_indexes <- TRUE

## Diagnostics
# Log every input change. Installs an observer per input and writes a line on
# every keystroke, so keep it FALSE unless you are debugging.
enable_input_logging <- FALSE

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
