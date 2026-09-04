#### Save/Load module ####
# Each module reports and restores its own state, so this module only has to collect
# and distribute. The screen-dependent pickers preserve their values across a choices
# update, so restoring a key takes a single pass.

CRAVE_KEY_FORMAT <- 2L

# Old flat input ids, mapped to (module, new id). Lets keys saved by CRAVE 4.x
# still load.
LEGACY_KEY_MAP <- list(
  explore = c(
    genes_selected_text = "genes_text",
    genes_selected = "genes",
    comparisons_selected = "contrasts",
    experiments_table_rows_selected = "_experiment_rows",
    comparison_selector_citation = "sel_citation",
    comparison_selector_kind = "sel_kind",
    comparison_selector_library = "sel_library",
    comparison_selector_timepoint = "sel_timepoint",
    comparison_selector_endpoint = "sel_endpoint",
    comparison_selector_organism = "sel_organism",
    comparison_selector_source = "sel_source",
    comparison_selector_daysgrown_diff = "sel_daysgrown_diff",
    comparison_selector_daysgrown_ref = "sel_daysgrown_ref",
    comparison_selector_treatment_diff = "sel_treatment_diff",
    comparison_selector_treatment_ref = "sel_treatment_ref",
    comparison_selector_dose_diff = "sel_dose_diff",
    comparison_selector_dose_ref = "sel_dose_ref",
    comparison_selector_knockout_diff = "sel_knockout_diff",
    comparison_selector_knockout_ref = "sel_knockout_ref",
    comparison_selector_cell_line_diff = "sel_cellline_diff",
    comparison_selector_cell_line_ref = "sel_cellline_ref",
    comparison_volcano = "volcano_contrast",
    analysis_method_se_v = "volcano_method",
    comparison_rank = "rank_contrast",
    analysis_method_se_r = "rank_method",
    comparison_X = "biplot_x",
    comparison_Y = "biplot_y",
    analysis_method_se_b_x = "biplot_method_x",
    analysis_method_se_b_y = "biplot_method_y",
    comparison_overlap = "overlap_contrasts",
    overlap_cutoff_method = "overlap_method",
    overlap_cutoff_stat = "overlap_stat",
    overlap_cutoff_value = "overlap_value",
    overlap_cutoff_direction = "overlap_direction",
    rocauc_base = "roc_base",
    rocauc_method = "roc_method",
    rocauc_cutoff_stat = "roc_stat",
    rocauc_cutoff_value = "roc_value",
    rocauc_cutoff_direction = "roc_direction",
    rocauc_head = "roc_head"
  ),
  correlate = c(
    genes_queried_text = "genes_text",
    genes_queried = "genes",
    gq_stat = "stat",
    gq_cutoff = "cutoff",
    analysis_method_gq = "method",
    gq_filter_citation = "f_citation",
    gq_filter_kind = "f_kind",
    gq_filter_library = "f_library",
    gq_filter_timepoint = "f_timepoint",
    gq_filter_organism = "f_organism",
    gq_filter_source = "f_source",
    gq_filter_days_grown_diff = "f_days_grown_diff",
    gq_filter_days_grown_ref = "f_days_grown_ref",
    gq_filter_treatment_diff = "f_treatment_diff",
    gq_filter_treatment_ref = "f_treatment_ref",
    gq_filter_dose_diff = "f_dose_diff",
    gq_filter_dose_ref = "f_dose_ref",
    gq_filter_knockout_diff = "f_knockout_diff",
    gq_filter_knockout_ref = "f_knockout_ref",
    gq_filter_cellline_diff = "f_cellline_diff",
    gq_filter_cellline_ref = "f_cellline_ref",
    gq_filter_contrast = "f_contrast",
    gq_filter_custom = "f_custom",
    genequery_violin = "gq_plottype",
    corr_cutoff = "corr_cutoff",
    enrichment_cutoff_stat = "enr_stat",
    enrichment_cutoff_value = "enr_value",
    enrichment_cutoff_direction = "enr_direction",
    enrichment_cutoff_tail = "enr_tail",
    pendragonator_n = "pdg_n",
    pendragonator_stat = "pdg_stat",
    download_filename = "dl_filename"
  ),
  exorcise = c(
    exorcise_seq = "seq",
    exorcise_orig = "orig",
    exorcise_control_string = "control_string",
    exorcise_pam = "pam",
    exorcise_genome = "genome",
    exorcise_mode = "mode",
    exorcise_mode_advanced = "mode_advanced"
  ),
  guides = c(
    guideLibrary_targets = "targets",
    guideLibrary_geneclasses = "geneclasses",
    guideLibrary_chromosomes = "chromosomes",
    guideLibrary_chemistry = "chemistry",
    guideLibrary_assembly = "assembly",
    guideLibrary_pam = "pam",
    guideLibrary_library = "library",
    guideLibrary_organism = "organism"
  )
)

#' Translate a CRAVE 4.x flat key into the module-keyed format.
translateLegacyKey <- function(flat) {
  out <- list()
  for (module in names(LEGACY_KEY_MAP)) {
    map <- LEGACY_KEY_MAP[[module]]
    part <- list()
    for (old in names(map)) {
      if (!is.null(flat[[old]])) part[[map[[old]]]] <- flat[[old]]
    }
    if (length(part)) out[[module]] <- part
  }
  out
}

#' Read a key file, accepting both the current and the CRAVE 4.x formats.
readKeyFile <- function(path) {
  # Current format: a plain RDS.
  asRds <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.list(asRds) && identical(asRds$format, CRAVE_KEY_FORMAT)) {
    return(list(state = asRds$state, legacy = FALSE))
  }
  # CRAVE 4.x format: an .RData holding one object called `key`.
  env <- new.env(parent = emptyenv())
  ok <- tryCatch({ base::load(path, envir = env); TRUE }, error = function(e) FALSE)
  if (ok && !is.null(env$key)) {
    return(list(state = translateLegacyKey(env$key), legacy = TRUE))
  }
  NULL
}

#' Server for the Save/Load modal.
#'
#' @param modules Named list of module handles, each with $state() and $restore().
#' @param opened Reactive that becomes truthy when the user asks for the dialog.
saveLoadServer <- function(id, modules, opened) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    notice <- reactiveVal("No key uploaded.")

    observeEvent(opened(), {
      showModal(modalDialog(
        title = "Save/Load",
        easyClose = TRUE,
        footer = NULL,
        tags$h4("Save"),
        tags$p("Download the state of your current session as a key file."),
        downloadButton(ns("save"), "Download key", class = "btn-success"),
        tags$br(), tags$hr(),
        tags$h4("Load"),
        tags$p("Restore the state of a previous session by uploading a key file."),
        fileInput(ns("upload"), "Upload key", accept = c(".rds", ".txt", ".RData")),
        textOutput(ns("console")),
        tags$br(),
        actionButton(ns("load"), "Load key", icon = icon("folder-open"), class = "btn-success"),
        tags$hr(),
        modalButton("Close", icon = icon("xmark"))
      ))
    }, ignoreInit = TRUE)

    output$console <- renderText(notice())

    observeEvent(input$upload, ignoreNULL = FALSE, {
      notice(if (is.null(input$upload)) "No key uploaded."
             else paste0("Key uploaded: ", input$upload$name, "."))
    })

    output$save <- downloadHandler(
      filename = function() paste0("crave-key-", format(Sys.Date(), "%Y%m%d"), ".rds"),
      content = function(f) {
        state <- lapply(modules, function(m) m$state())
        saveRDS(list(format = CRAVE_KEY_FORMAT,
                     version = CRAVE_VERSION,
                     saved = Sys.time(),
                     state = state), f)
      }
    )

    observeEvent(input$load, {
      if (is.null(input$upload)) {
        notice("Please upload a key.")
        return()
      }
      parsed <- readKeyFile(input$upload$datapath[1])
      if (is.null(parsed)) {
        notice("That file is not a CRAVE key, or it is corrupt.")
        return()
      }
      for (nm in intersect(names(modules), names(parsed$state))) {
        tryCatch(modules[[nm]]$restore(parsed$state[[nm]]),
                 error = function(e) log_warn(paste0(
                   "Could not restore state for '", nm, "': ", conditionMessage(e))))
      }
      notice(paste0(
        "State loaded from ", input$upload$name,
        if (isTRUE(parsed$legacy)) " (converted from an older CRAVE key format)." else "."
      ))
    })
  })
}
