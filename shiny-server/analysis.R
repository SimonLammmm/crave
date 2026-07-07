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

# Assemble the 5-sheet design workbook from in-app inputs and write it to path.
# details   : named list of Experiment-details key/value pairs
# samples   : data.frame with the Sample details columns
# controls  : data.frame with Group/Control sample/Test sample/Contrast
# analyses  : data.frame with the Analyses columns
analysis_build_workbook <- function(path, details, samples, controls, analyses) {
  wb <- openxlsx::createWorkbook()
  
  # Experiment details: key/value, two columns, no header
  openxlsx::addWorksheet(wb, "Experiment details")
  ed <- data.frame(key = names(details), value = unlist(details, use.names = FALSE),
                   stringsAsFactors = FALSE)
  openxlsx::writeData(wb, "Experiment details", ed, colNames = FALSE)
  
  openxlsx::addWorksheet(wb, "Sample details")
  openxlsx::writeData(wb, "Sample details", samples, colNames = TRUE)
  
  openxlsx::addWorksheet(wb, "Control groups")
  openxlsx::writeData(wb, "Control groups", controls, colNames = TRUE)
  
  openxlsx::addWorksheet(wb, "Analyses")
  openxlsx::writeData(wb, "Analyses", analyses, colNames = TRUE)
  
  # gubbins helper sheet: pipeline tolerates it present; write sample names
  openxlsx::addWorksheet(wb, "gubbins")
  openxlsx::writeData(wb, "gubbins", data.frame(x = samples$Sample), colNames = FALSE)
  
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}

# Read sample (column) names from a counts matrix header: guide, gene, <samples...>
analysis_counts_samples <- function(counts_path) {
  out <- tryCatch({
    hdr <- readLines(counts_path, n = 1)
    # detect tab vs whitespace
    parts <- strsplit(hdr, "\t")[[1]]
    if (length(parts) < 3) parts <- strsplit(trimws(hdr), "\\s+")[[1]]
    parts[-(1:2)]  # drop guide + gene
  }, error = function(e) character(0))
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
      tags$p("Run the CRISPR pipeline (drugZ / MAGeCK) via the Exorcise Docker image. Either upload an experiment-design workbook, or build one in the app."),
      tabsetPanel(
        id = "analysis_mode",
        
        ## Mode 1: upload an existing design workbook
        tabPanel(
          "Upload design",
          tags$br(),
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
          )
        ),
        
        ## Mode 2: build the design in-app
        tabPanel(
          "Build design",
          tags$br(),
          fileInput("build_counts", "Counts matrix (.tsv / .tsv.gz) — sample names are read from its header"),
          tags$h4("Experiment details"),
          fluidRow(
            column(4, textInput("build_expname", "Experiment name (required)", "")),
            column(4, textInput("build_library", "Library", "")),
            column(4, selectInput("build_organism", "Organism", c("Human", "Mouse", "Other"), "Human"))
          ),
          fluidRow(
            column(6, textInput("build_doi", "DOI", "")),
            column(6, textInput("build_citation", "Citation", ""))
          ),
          textAreaInput("build_description", "Description", "", width = "100%"),
          tags$hr(),
          tags$h4("Sample details"),
          tags$p("One row per sample (seeded from the counts header). Treatment / KO / Cell line accept known values or new ones."),
          DTOutput("build_samples"),
          tags$hr(),
          tags$h4("Control groups"),
          tags$p("Define which samples are control vs test for each comparison."),
          DTOutput("build_controls"),
          tags$hr(),
          tags$h4("Analysis"),
          fluidRow(
            column(4, checkboxGroupInput("build_methods", "Methods", c("drugz", "mageck"), c("drugz", "mageck"))),
            column(4, numericInput("build_pseudocount", "Add pseudocount", 5, min = 0)),
            column(4, checkboxInput("build_paired", "Paired", FALSE))
          ),
          actionButton("build_run", "Build design & run", icon = icon("play"), class = "btn-success"),
          tags$br(), tags$br(),
          uiOutput("build_downloads"),
          tags$h4("Log"),
          tags$pre(style = "max-height: 300px; overflow-y: auto; font-size: 11px;",
                   textOutput("build_console"))
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
  
  # Shared: given a completed analysis_run() result, populate tables + downloads.
  show_results <- function(res, console_out, downloads_out, download_id) {
    output[[console_out]] <- renderText(res$console)
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
      output[[downloads_out]] <- renderUI(
        downloadButton(download_id, "Download all results (.zip)", class = "btn-info"))
    } else {
      result_dir(NULL)
      output$analysis_drugz  <- renderDT(NULL)
      output$analysis_mageck <- renderDT(NULL)
      output[[downloads_out]] <- renderUI(tags$em("No results produced. See the log above."))
    }
  }
  
  ## --- Mode 1: upload design ---
  observeEvent(input$analysis_run, {
    req(input$analysis_counts, input$analysis_design)
    if (!nzchar(get0("exorcise_docker", ifnotfound = ""))) {
      output$analysis_console <- renderText("exorcise_docker is not set in config.R.")
      return()
    }
    expected <- analysis_expected_counts(input$analysis_design$datapath)
    counts_name <- if (!is.null(expected)) expected else input$analysis_counts$name
    output$analysis_console <- renderText("Running... this can take a while for large screens.")
    res <- analysis_run(input$analysis_counts$datapath, counts_name, input$analysis_design$datapath)
    show_results(res, "analysis_console", "analysis_downloads", "analysis_download_zip")
  })
  
  ## --- Mode 2: build design in-app ---
  # Dropdown suggestions from loaded comparisons metadata (fall back to blanks)
  sugg <- function(col) {
    v <- tryCatch(unique(comparisons[[col]]), error = function(e) NULL)
    v <- v[!is.na(v) & nzchar(as.character(v))]
    sort(as.character(v))
  }
  treat_opts <- sugg("Treatment"); ko_opts <- sugg("KO"); cell_opts <- sugg("Cell line")
  
  build_samples_rv <- reactiveVal(NULL)
  build_controls_rv <- reactiveVal(data.frame(
    Group = "main", `Control sample` = "", `Test sample` = "", Contrast = "",
    check.names = FALSE, stringsAsFactors = FALSE))
  
  # Seed the sample table from the uploaded counts header
  observeEvent(input$build_counts, {
    samps <- analysis_counts_samples(input$build_counts$datapath)
    if (length(samps) == 0) return()
    df <- data.frame(
      Replicate = samps, Sample = sub("[0-9]+$", "", samps), Label = "",
      Treatment = "No treatment", Dose = "N/A", `Growth inhibition %` = "Unspecified",
      `Days grown` = "", Doublings = "Unspecified", Trajectory = "Unspecified",
      Clone = "Unspecified", `Cell line` = "Unspecified", KO = "",
      Notes = "", check.names = FALSE, stringsAsFactors = FALSE)
    build_samples_rv(df)
  })
  
  output$build_samples <- renderDT({
    df <- build_samples_rv()
    if (is.null(df)) df <- data.frame(Message = "Upload a counts matrix to seed samples.")
    datatable(df, editable = "cell", rownames = FALSE,
              options = list(dom = "t", scrollX = TRUE))
  }, server = FALSE)
  observeEvent(input$build_samples_cell_edit, {
    df <- build_samples_rv(); info <- input$build_samples_cell_edit
    df[info$row, info$col + 1] <- info$value; build_samples_rv(df)
  })
  
  output$build_controls <- renderDT({
    datatable(build_controls_rv(), editable = "cell", rownames = FALSE,
              options = list(dom = "t", scrollX = TRUE))
  }, server = FALSE)
  observeEvent(input$build_controls_cell_edit, {
    df <- build_controls_rv(); info <- input$build_controls_cell_edit
    df[info$row, info$col + 1] <- info$value; build_controls_rv(df)
  })
  
  observeEvent(input$build_run, {
    req(input$build_counts)
    if (!nzchar(input$build_expname)) {
      output$build_console <- renderText("Experiment name is required.")
      return()
    }
    if (is.null(build_samples_rv())) {
      output$build_console <- renderText("No samples — upload a counts matrix first.")
      return()
    }
    counts_name <- input$build_counts$name
    
    details <- list(
      "Experiment name" = input$build_expname,
      "Date published"  = as.character(Sys.Date()),
      "Library"         = input$build_library,
      "Multiplicity of infection" = "0",
      "Representation average"    = "0",
      "Experiment description (a few sentences)" = input$build_description,
      "Notes"    = "",
      "DOI"      = input$build_doi,
      "Citation" = input$build_citation,
      "Organism" = input$build_organism
    )
    analyses <- data.frame(
      `Control group` = build_controls_rv()$Group,
      Paired = input$build_paired,
      Method = paste(input$build_methods, collapse = ","),
      `Counts file` = counts_name,
      `Add pseudocount` = input$build_pseudocount,
      Arguments = "", Name = "",
      check.names = FALSE, stringsAsFactors = FALSE)
    
    workdir <- file.path(tempdir(), paste0("build_", paste0(sample(c(LETTERS[1:6], 0:9), 8, TRUE), collapse = "")))
    dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
    design_path <- file.path(workdir, "design.xlsx")
    tryCatch(
      analysis_build_workbook(design_path, details, build_samples_rv(),
                              build_controls_rv(), analyses),
      error = function(e) { output$build_console <- renderText(paste("Failed to build workbook:", conditionMessage(e))); NULL }
    )
    if (!file.exists(design_path)) return()
    
    output$build_console <- renderText("Building design and running... this can take a while.")
    res <- analysis_run(input$build_counts$datapath, counts_name, design_path)
    show_results(res, "build_console", "build_downloads", "analysis_download_zip")
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