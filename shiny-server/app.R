#### CRAVE ####
# Contrasts-based Resource for Analysis, Visualisation, and Exploration
#
# This file is the entry point and nothing else. Everything it needs lives in R/,
# sourced in numeric order:
#
#   00_packages.R            eager and deferred package loading
#   01_constants.R           version, limits, defaults, static reference data
#   02_utils.R               small shared helpers
#   03_data.R                dataset loading, on-disk metadata cache, indexing
#   04_query.R               parameterised statistics queries
#   05_analysis_common.R     primitives shared between analyses
#   06_analysis_explore.R    volcano, rank, biplot, overlap, ROC
#   07_analysis_correlate.R  gene query, clustergram, heatmap, network, UMAP,
#                            enrichment, Pendragonator, bulk download, guides
#   08_analysis_exorcise.R   Exorcise
#   09_config.R              config.R loading and validation
#   10..18                   Shiny modules
#   19_ui.R, 20_server.R     assembly
#
# See CHANGELOG.md for what changed in 5.1.2.

#### 1. Load the codebase, exactly once ####
# Shiny (>= 1.5.0) sources every .R file in an app directory's R/ subdirectory
# before it evaluates app.R. Running app.R line by line at the console, or under an
# older Shiny, it does not. 00_packages.R sets the crave.sourced option, so this
# loop runs only when Shiny has not already done the work; otherwise every file is
# sourced twice, into two different environments.
if (!isTRUE(getOption("crave.sourced"))) {
  for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
}

log_threshold(INFO)

#### 2. Start up, announcing each stage before it begins ####
# Each stage logs before it runs, not after, so a stall is always attributable to a
# named step. The index build in particular can take minutes on a large dataset.
stage <- local({
  t0 <- Sys.time()
  function(msg) log_info(sprintf(
    "[%6.1fs] %s", as.numeric(difftime(Sys.time(), t0, units = "secs")), msg))
})

stage(paste0("Starting CRAVE ", CRAVE_VERSION, "."))

stage("Reading config.R ...")
loaded <- loadCraveConfig("config.R")
config <- loaded$config
for (w in loaded$warnings) log_warn(w)
for (e in loaded$errors)   log_error(e)

craveData <- NULL
if (length(loaded$errors) == 0) {
  stage(paste0("Loading ", length(config$datasets), " dataset(s) ..."))
  craveData <- loadCraveData(
    config$datasets, config$exorcise_root, config$exorcise_docker,
    create_indexes = isTRUE(config$create_indexes)
  )
}

#### 3. Go / no-go ####
startup <- if (length(loaded$errors) > 0) {
  list(ok = FALSE, message = paste0(
    "<ul>", paste0("<li>", loaded$errors, "</li>", collapse = ""), "</ul>",
    "<p>A template is available at <a href=\"",
    "https://github.com/SimonLammmm/crave/blob/main/shiny-server/config.example.R",
    "\">config.example.R</a>.</p>"))
} else if (!isTRUE(craveData$caps$explore)) {
  list(ok = FALSE, message = paste0(
    "<p>The configured datasets loaded no experiments. Check that each dataset ",
    "directory contains database.db, experiments_metadata.csv.gz and ",
    "comparisons_metadata.csv.gz, and that at least one experiment has a ",
    "non-empty Citation.</p>"))
} else {
  list(ok = TRUE)
}

#### 4. Build the UI and server, then call shinyApp() once, at top level ####
if (isTRUE(startup$ok)) {
  # App-level, so a refresh in one session is seen by all of them.
  dataStore <- reactiveVal(craveData)
  onStop(function() {
    log_info("Shutting down; closing dataset connections.")
    closeCraveData(isolate(dataStore()))
  })
  craveAppUI     <- craveUI(craveData$caps, config)
  craveAppServer <- craveServer(config, dataStore)
  stage("Ready.")
} else {
  craveAppUI     <- craveErrorUI(startup$message)
  craveAppServer <- function(input, output, session) {}
  stage("Serving the startup error page.")
}

# NB no shiny.host / shiny.port here, deliberately. Whoever launches the app owns
# those: RStudio's Run App picks a port and expects the app to use it, and
# shiny-server and ShinyProxy assign one and pass it in. Overriding them from
# inside app.R makes the app bind an address the launcher is not watching, which
# looks exactly like a hang. Docker passes them explicitly in its CMD.
options(shiny.maxRequestSize = (config$max_upload_mb %||% 100) * 1024 * 1024)

shinyApp(ui = craveAppUI, server = craveAppServer)
