#### Run analyses (crispr_pipeline) ####
## Stage 1: upload a counts matrix and an experiment-design workbook (.xlsx),
## run the Exorcise Docker image's crispr_pipeline, then show the result tables
## and offer the full output for download.
##
## Requires (config.R):
##   exorcise_docker    - image name, e.g. "simonlammmm/exorcise:1.5.3.14_arm64"
##   exorcise_platform  - optional docker --platform value (e.g. "linux/amd64" on
##                        Apple Silicon running an amd64 image); "" to omit
##   enable_analysis    - TRUE to expose the tab (deployer opt-in)
##
## The design workbook's Analyses sheet names the counts file; the uploaded
## counts file is saved under that expected name so the pipeline finds it.

## --- Helpers -----------------------------------------------------------------

# Read the counts filename that the design workbook expects (Analyses sheet)
analysis_expected_counts <- function(xlsx_path) {
  out <- tryCatch({
    sheets <- readxl::excel_sheets(xlsx_path)
    if (!("Analyses" %in% sheets)) return(NULL)
    an <- readxl::read_excel(xlsx_path, sheet = "Analyses")
    cf <- an[["Counts file"]]
    cf <- cf[!is.na(cf) & nzchar(cf)]
    if (length(cf) == 0) NULL else basename(cf[1])
  }, error = function(e) NULL)
  out
}

# Read the experiment name (used as the output subdirectory by the pipeline)
analysis_experiment_name <- function(xlsx_path) {
  out <- tryCatch({
    ed <- suppressMessages(readxl::read_excel(xlsx_path, sheet = "Experiment details",
                                              col_names = FALSE))
    # key/value sheet: find the row whose first column is "Experiment name"
    key <- as.character(ed[[1]])
    val <- as.character(ed[[2]])
    nm <- val[which(key == "Experiment name")[1]]
    if (is.na(nm) || !nzchar(nm)) NULL else nm
  }, error = function(e) NULL)
  out
}

# Build and run the crispr_pipeline docker command in a temp working dir.
# Returns list(ok, console, workdir, expname).
analysis_run <- function(counts_src, counts_name, design_src) {
  taskid  <- paste0(sample(c(LETTERS[1:6], 0:9), 8, TRUE), collapse = "")
  workdir <- file.path(tempdir(), paste0("analysis_", taskid))
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  
  design_dst <- file.path(workdir, "design.xlsx")
  counts_dst <- file.path(workdir, counts_name)
  file.copy(design_src, design_dst, overwrite = TRUE)
  file.copy(counts_src, counts_dst, overwrite = TRUE)
  
  platform <- get0("exorcise_platform", ifnotfound = "")
  platflag <- if (nzchar(platform)) paste0(" --platform ", platform) else ""
  
  command <- paste0(
    "docker run --rm", platflag,
    " -v ", workdir, ":/data ",
    exorcise_docker,
    " crispr_pipeline",
    " --counts /data/", counts_name,
    " --output-dir /data/output",
    " --overwrite",
    " /data/design.xlsx",
    " 2>&1"
  )
  
  console <- tryCatch(
    paste(system(command, intern = TRUE), collapse = "\n"),
    error = function(e) paste("Failed to run docker:", conditionMessage(e))
  )
  
  expname <- analysis_experiment_name(design_dst)
  tables_dir <- if (!is.null(expname)) file.path(workdir, "output", expname, "tables") else NULL
  ok <- !is.null(tables_dir) && dir.exists(tables_dir) &&
    length(list.files(tables_dir, pattern = "\\.csv$")) > 0
  
  list(ok = ok, console = console, workdir = workdir, expname = expname)
}

## --- UI -----------------------------------------------------------------------

analysisTab <- if (isTRUE(get0("enable_analysis", ifnotfound = FALSE))) {
  tabPanel(
    "Run analysis", icon = icon("gears"),
    fluidPage(
      tags$h3("Run analysis"),
      tags$p("Upload a counts matrix and an experiment-design workbook, then run the CRISPR pipeline (drugZ / MAGeCK) via the Exorcise Docker image."),
      fluidRow(
        column(5,
               fileInput("analysis_counts", "Counts matrix (.tsv / .tsv.gz)"),
               fileInput("analysis_design", "Experiment design (.xlsx)"),
               actionButton("analysis_run", "Run analysis", icon = icon("play"), class = "btn-success"),
               tags$br(), tags$br(),
               uiOutput("analysis_downloads")
        ),
        column(7,
               tags$h4("Log"),
               tags$pre(style = "max-height: 300px; overflow-y: auto; font-size: 11px;",
                        textOutput("analysis_console"))
        )
      ),
      tags$hr(),
      tags$h4("drugZ results"),
      DTOutput("analysis_drugz"),
      tags$hr(),
      tags$h4("MAGeCK results"),
      DTOutput("analysis_mageck")
    )
  )
} else NULL

## --- Server -------------------------------------------------------------------

analysisServer <- function(input, output, session) {
  if (!isTRUE(get0("enable_analysis", ifnotfound = FALSE))) return(invisible())
  
  result_dir <- reactiveVal(NULL)
  
  observeEvent(input$analysis_run, {
    req(input$analysis_counts, input$analysis_design)
    
    if (!nzchar(get0("exorcise_docker", ifnotfound = ""))) {
      output$analysis_console <- renderText("exorcise_docker is not set in config.R. Set it to the Exorcise image name.")
      return()
    }
    
    # The counts file must be named as the design's Analyses sheet expects
    expected <- analysis_expected_counts(input$analysis_design$datapath)
    counts_name <- if (!is.null(expected)) expected else input$analysis_counts$name
    
    output$analysis_console <- renderText("Running... this can take a while for large screens.")
    
    res <- analysis_run(input$analysis_counts$datapath, counts_name, input$analysis_design$datapath)
    output$analysis_console <- renderText(res$console)
    
    if (res$ok) {
      result_dir(file.path(res$workdir, "output", res$expname))
      tdir <- file.path(res$workdir, "output", res$expname, "tables")
      
      dz <- list.files(tdir, pattern = "drugz.*\\.csv$", full.names = TRUE)
      mg <- list.files(tdir, pattern = "mageck.*\\.csv$", full.names = TRUE)
      
      output$analysis_drugz <- renderDT({
        if (length(dz) == 0) return(NULL)
        datatable(data.table::fread(dz[1]), rownames = FALSE,
                  options = list(iDisplayLength = 15, scrollX = TRUE))
      }, server = TRUE)
      output$analysis_mageck <- renderDT({
        if (length(mg) == 0) return(NULL)
        datatable(data.table::fread(mg[1]), rownames = FALSE,
                  options = list(iDisplayLength = 15, scrollX = TRUE))
      }, server = TRUE)
      
      output$analysis_downloads <- renderUI({
        downloadButton("analysis_download_zip", "Download all results (.zip)", class = "btn-info")
      })
    } else {
      result_dir(NULL)
      output$analysis_drugz  <- renderDT(NULL)
      output$analysis_mageck <- renderDT(NULL)
      output$analysis_downloads <- renderUI(tags$em("No results produced. See the log above."))
    }
  })
  
  output$analysis_download_zip <- downloadHandler(
    filename = function() paste0("analysis-results-", Sys.Date(), ".zip"),
    content = function(f) {
      rd <- result_dir()
      req(rd)
      old <- setwd(dirname(rd))
      on.exit(setwd(old), add = TRUE)
      utils::zip(zipfile = f, files = basename(rd), flags = "-r9X")
    }
  )
}