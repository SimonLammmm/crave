#### Deployment configuration ####
# config.R is sourced into its own environment and merged over a complete set of
# defaults, then validated, so a config file that omits a setting produces a message
# the deployer can act on rather than an "object not found" from inside the UI.

CRAVE_CONFIG_DEFAULTS <- list(
  datasets              = list(),
  exorcise_root         = NULL,
  exorcise_hostroot     = NULL,
  exorcise_docker       = NULL,
  portal_theme          = "cerulean",
  portal_logo           = NULL,
  portal_name           = "CRAVE",
  portal_subtitle       = "CRISPR results app for visualisation and exploration",
  front_page_html       = "",
  portal_footer         = NULL,
  privacy_toast         = NULL,
  notices_info          = NULL,
  notices_warning       = NULL,
  notices_emergency     = NULL,
  legal_text            = NULL,
  create_indexes        = TRUE,
  enable_input_logging  = FALSE,
  max_upload_mb         = 100,
  # Tabs and sub-tabs a deployment can withhold. The two that write to a dataset
  # are off unless asked for: a public portal has no business offering them, and a
  # deployer who does want them has to say so rather than discover them.
  enable_submit         = FALSE,
  enable_manage         = FALSE,
  enable_pendragonator  = TRUE
)

# Settings that belong to whatever launches the app, not to CRAVE. Recognised only
# so that a config.R still setting them gets told, rather than being ignored in
# silence.
CRAVE_CONFIG_IGNORED <- c("shiny_host", "shiny_port")

VALID_SHINY_THEMES <- c(
  "cerulean", "cosmo", "cyborg", "darkly", "flatly", "journal", "lumen", "paper",
  "readable", "sandstone", "simplex", "slate", "spacelab", "superhero", "united", "yeti"
)

#' Source config.R and merge it over the defaults.
#'
#' @return list(config, errors, warnings)
loadCraveConfig <- function(path = "config.R") {
  errors <- character(0)
  warnings <- character(0)

  if (!file.exists(path)) {
    return(list(config = NULL, errors = paste0(
      "config.R was not found. Copy config.example.R to config.R, edit it, and ",
      "put it next to app.R. If you are using Docker Compose, bind it at /app/config.R."
    ), warnings = warnings))
  }

  env <- new.env(parent = globalenv())
  ok <- tryCatch({ sys.source(path, envir = env); TRUE },
                 error = function(e) {
                   errors <<- c(errors, paste0("config.R could not be read: ",
                                               htmltools::htmlEscape(conditionMessage(e))))
                   FALSE
                 })
  if (!ok) return(list(config = NULL, errors = errors, warnings = warnings))

  config <- CRAVE_CONFIG_DEFAULTS
  for (nm in ls(env, all.names = FALSE)) {
    value <- get(nm, envir = env)
    if (is.function(value)) next
    config[[nm]] <- value
  }

  # -- Validate ---------------------------------------------------------------
  if (!is.list(config$datasets) || length(config$datasets) == 0) {
    errors <- c(errors, "config.R defines no datasets. Add at least one entry to `datasets`.")
  } else {
    for (i in seq_along(config$datasets)) {
      ds <- config$datasets[[i]]
      label <- ds$name %||% paste0("entry ", i)
      if (is.null(ds$path) || !nzchar(ds$path)) {
        errors <- c(errors, paste0("Dataset '", label, "' has no `path`."))
        next
      }
      if (!dir.exists(ds$path)) {
        errors <- c(errors, paste0(
          "Dataset '", label, "' points at '", htmltools::htmlEscape(ds$path),
          "', which does not exist. Paths must be absolute, and inside the container ",
          "if you are using Docker."))
      } else if (!file.exists(file.path(ds$path, "database.db"))) {
        errors <- c(errors, paste0(
          "Dataset '", label, "' has no database.db in '",
          htmltools::htmlEscape(ds$path), "'."))
      }
      if (is.null(ds$name) || !nzchar(ds$name)) {
        warnings <- c(warnings, paste0("Dataset ", i, " has no `name`; using its path."))
        config$datasets[[i]]$name <- basename(ds$path)
      }
      if (!is.null(ds$hostpath) && !startsWith(ds$hostpath, "/")) {
        warnings <- c(warnings, paste0(
          "Dataset '", label, "' has a relative `hostpath`. It names the directory ",
          "as the Docker daemon sees it on the host, so it has to be absolute. ",
          "Ignoring it."))
        config$datasets[[i]]$hostpath <- NULL
      }
    }
    names <- vapply(config$datasets, function(d) d$name %||% "", character(1))
    if (anyDuplicated(names)) {
      errors <- c(errors, "Two or more datasets share the same `name`. Names must be unique.")
    }
  }

  stale <- intersect(CRAVE_CONFIG_IGNORED, ls(env, all.names = FALSE))
  if (length(stale)) {
    warnings <- c(warnings, paste0(
      "config.R sets ", paste(stale, collapse = " and "),
      ", which CRAVE no longer applies. The host and port belong to whatever ",
      "launches the app: RStudio's Run App, shiny-server and ShinyProxy all choose ",
      "them, and Docker passes them in its CMD. For a bare Rscript deployment, use ",
      "shiny::runApp('shiny-server', host = ..., port = ...) instead."))
  }

  # -- Submit and Manage ------------------------------------------------------
  # Both write to a dataset through the Exorcise container, so both need an image
  # to run and a directory they can write to. Saying so here, once, beats a tab
  # that appears and then fails on the first click.
  wantsWriting <- isTRUE(config$enable_submit) || isTRUE(config$enable_manage)
  if (wantsWriting) {
    asked <- paste(c(if (isTRUE(config$enable_submit)) "Submit",
                     if (isTRUE(config$enable_manage)) "Manage"), collapse = " and ")
    if (is.null(config$exorcise_docker) || !nzchar(config$exorcise_docker)) {
      warnings <- c(warnings, paste0(
        asked, " is enabled but `exorcise_docker` names no image. Both tabs run ",
        "exorcise-analyse and exorcise-database in that image, so they will stay ",
        "hidden until it is set."))
    }
    writable <- vapply(config$datasets, function(d) {
      !is.null(d$path) && dir.exists(d$path) && file.access(d$path, 2) == 0
    }, logical(1))
    if (length(writable) && !any(writable)) {
      warnings <- c(warnings, paste0(
        asked, " is enabled but no configured dataset directory is writable. ",
        "Both tabs rewrite the dataset in place, so they will stay hidden."))
    } else if (any(!writable)) {
      readonly <- vapply(config$datasets[!writable], function(d) d$name %||% d$path,
                         character(1))
      warnings <- c(warnings, paste0(
        "Dataset(s) ", paste(readonly, collapse = ", "), " are not writable, so ",
        asked, " will not offer them."))
    }
  }

  if (!config$portal_theme %in% VALID_SHINY_THEMES) {
    warnings <- c(warnings, paste0(
      "portal_theme '", config$portal_theme, "' is not a shinytheme; using 'cerulean'."))
    config$portal_theme <- "cerulean"
  }

  if (is.null(config$exorcise_hostroot)) config$exorcise_hostroot <- config$exorcise_root

  if (is.null(config$legal_text)) {
    config$legal_text <- tags$p("No legal text has been configured for this deployment.")
  }

  list(config = config, errors = errors, warnings = warnings)
}
