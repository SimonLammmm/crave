#### Dataset loading ####
# The dominant cost of starting CRAVE, so three things are done to keep it down:
# the metadata wrangle is cached to an .rds keyed on the size and mtime of the
# source files; the gene symbol list is filtered in a single pass; and the Exorcise
# probe short-circuits unless the reference directory actually holds a .2bit.

#' Locate a writable cache directory for a dataset.
#'
#' Prefers a hidden directory inside the dataset so the cache survives container
#' restarts; falls back to the session tempdir for read-only datasets.
craveCacheDir <- function(path_dataset) {
  preferred <- file.path(path_dataset, ".crave-cache")
  if (!dir.exists(preferred)) {
    ok <- isTRUE(suppressWarnings(dir.create(preferred, showWarnings = FALSE, recursive = TRUE)))
    if (!ok && !dir.exists(preferred)) {
      fallback <- file.path(tempdir(), "crave-cache")
      dir.create(fallback, showWarnings = FALSE, recursive = TRUE)
      return(fallback)
    }
  }
  if (file.access(preferred, 2) == 0) return(preferred)
  fallback <- file.path(tempdir(), "crave-cache")
  dir.create(fallback, showWarnings = FALSE, recursive = TRUE)
  fallback
}

#' Fingerprint a set of files by size and mtime.
fileFingerprint <- function(paths) {
  info <- file.info(paths)
  paste0(
    "v", CRAVE_CACHE_VERSION, "-",
    paste(basename(paths), info$size, as.integer(info$mtime), sep = ":", collapse = "|")
  )
}

#' Ensure the indexes CRAVE's query patterns need exist on a dataset database.
#'
#' The `stat` table ships with only its primary key, (comparison_id, gene_id,
#' analysis_type_id). That serves Explore well, because Explore always filters by
#' comparison first. Correlate and Pendragonator filter by gene across many
#' comparisons, which the primary key cannot serve, so those queries fall back to
#' a full table scan. Adding a gene-leading index turns them into index lookups.
#'
#' Silently does nothing if the database file is read-only, or if the deployer has
#' set `create_indexes <- FALSE` in config.R.
ensureStatIndexes <- function(con, label = "", enabled = TRUE) {
  if (!isTRUE(enabled)) {
    log_info(paste0("Index creation disabled by config.R; skipping for ", label, "."))
    return(invisible(FALSE))
  }
  tryCatch({
    existing <- DBI::dbGetQuery(
      con, "SELECT name FROM sqlite_master WHERE type = 'index'"
    )$name
    wanted <- list(
      crave_idx_stat_gene = "CREATE INDEX crave_idx_stat_gene ON stat (gene_id, analysis_type_id)",
      crave_idx_stat_type = "CREATE INDEX crave_idx_stat_type ON stat (analysis_type_id, comparison_id)"
    )
    todo <- setdiff(names(wanted), existing)
    if (length(todo) == 0) {
      # Stated rather than left implicit: every other exit from this function logs,
      # so without this line "already indexed" would be indistinguishable from
      # "this code never ran" when reading a startup log.
      log_info(paste0("Indexes already present on ", label, "; skipping."))
      return(invisible(TRUE))
    }
    # Building an index means sorting the whole stat table: on a large dataset that
    # is minutes of work. It happens once per dataset.
    log_info(paste0(
      "Creating ", length(todo), " index(es) on ", label,
      ". This sorts the whole stat table and can take several minutes on a large ",
      "dataset; it happens once. Set create_indexes <- FALSE in config.R to skip it."
    ))
    for (nm in todo) {
      DBI::dbExecute(con, wanted[[nm]])
      log_info(paste0("Created index ", nm, " on ", label, "."))
    }
    log_info(paste0("Running ANALYZE on ", label, " ..."))
    DBI::dbExecute(con, "ANALYZE")
    log_info(paste0("Indexing of ", label, " done."))
    invisible(TRUE)
  }, error = function(e) {
    log_info(paste0(
      "Could not create indexes on ", label, " (", conditionMessage(e),
      "). CRAVE will still work, but gene-led queries will be slower. ",
      "Make the dataset directory writable to allow indexing."
    ))
    invisible(FALSE)
  })
}

# Columns a CRAVE dataset must provide. Checked up front so a malformed dataset
# produces one clear message rather than an "object not found" from the middle of a
# dplyr pipeline.
CRAVE_REQUIRED_EXPERIMENT_COLS <- c("Experiment ID", "Citation", "Organism", "DOI")
CRAVE_REQUIRED_COMPARISON_COLS <- c(
  "Experiment ID", "Comparison ID", "Contrast", "Library",
  "Treatment", "KO", "Dose", "Cell line", "Days grown",
  "ControlTreatment", "ControlKO", "ControlDose", "ControlCell line", "ControlDays grown"
)

#' Stop with a readable message if required columns are absent.
requireColumns <- function(df, needed, what, dataset) {
  missing <- setdiff(needed, names(df))
  if (length(missing) == 0) return(invisible(TRUE))
  stop(sprintf(
    "Dataset '%s': %s is missing the column(s) %s. CRAVE needs: %s.",
    dataset, what, paste(sQuote(missing), collapse = ", "),
    paste(needed, collapse = ", ")
  ), call. = FALSE)
}

#' Wrangle raw experiment and comparison metadata into CRAVE's internal shape.
#'
#' Kept separate from loadDataset() so that its result can be memoised on disk: this
#' is the expensive half of loading a dataset, and its inputs rarely change.
#'
#' @param ds One entry from the `datasets` list in config.R.
#' @return list(experiments, comparisons)
wrangleMetadata <- function(experiments, comparisons, ds) {
  # With citation_from_id, Citation is synthesised below rather than read, so it
  # is not required to be present in the file.
  neededExperiments <- if (isTRUE(ds$citation_from_id)) {
    setdiff(CRAVE_REQUIRED_EXPERIMENT_COLS, "Citation")
  } else {
    CRAVE_REQUIRED_EXPERIMENT_COLS
  }
  requireColumns(experiments, neededExperiments,
                 "experiments_metadata.csv.gz", ds$name)
  requireColumns(comparisons, CRAVE_REQUIRED_COMPARISON_COLS,
                 "comparisons_metadata.csv.gz", ds$name)

  if (isTRUE(ds$citation_from_id)) {
    experiments <- experiments %>%
      transmute(
        `Experiment ID`, Organism, DOI,
        Citation = sub("^(.+)_(NVS\\d+|HS\\d+|ST_J\\d+|FGC\\d+|SLX\\d+)(.*?)$",
                       "\\1 (\\2\\3)", `Experiment ID`)
      )
  }

  experiments <- experiments %>%
    filter(Citation != "") %>%
    select(`Experiment ID`, Citation, Organism, DOI) %>%
    mutate(Citation = if_else(
      grepl("https://doi.org", DOI),
      paste0("<a href=", sub("^.*?\\((.+)\\).*?$", "\\1", DOI),
             " target=\"_blank\">", Citation, "</a>"),
      Citation
    )) %>%
    left_join(comparisons, by = "Experiment ID") %>%
    summarise(
      Treatments   = paste(sort(unique(Treatment)),    collapse = ", "),
      Knockouts    = paste(sort(unique(KO)),           collapse = ", "),
      `Cell lines` = paste(sort(unique(`Cell line`)),  collapse = ", "),
      Libraries    = paste(sort(unique(Library)),      collapse = ", "),
      Source       = ds$name,
      .by = c(`Experiment ID`, Citation, Organism)
    )

  # "Unspecified" standardisation, applied column-wise.
  blankToUnspecified <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- "Unspecified"
    x
  }
  blankToNumeric <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- NA_character_
    suppressWarnings(as.numeric(x))
  }
  textCols <- c("Dose", "Treatment", "KO", "Cell line",
                "ControlDose", "ControlTreatment", "ControlKO", "ControlCell line",
                "Library")
  numCols  <- c("Days grown", "ControlDays grown")

  comparisons <- comparisons %>%
    left_join(experiments %>% select(`Experiment ID`, Citation, Organism),
              by = "Experiment ID") %>%
    filter(Citation != "") %>%
    mutate(across(all_of(textCols), blankToUnspecified)) %>%
    mutate(across(all_of(numCols),  blankToNumeric)) %>%
    # Standardise inert treatments
    mutate(
      Treatment        = if_else(Treatment        %in% c("DMSO", "N/A"), "No treatment", Treatment),
      ControlTreatment = if_else(ControlTreatment %in% c("DMSO", "N/A"), "No treatment", ControlTreatment)
    ) %>%
    # Standardise plasmid
    mutate(
      `ControlCell line` = if_else(ControlKO == "Plasmid", "Plasmid", `ControlCell line`),
      ControlKO          = if_else(`ControlCell line` == "Plasmid", "Wildtype", ControlKO)
    ) %>%
    # Timepoint identity, relative to the other timepoints of the same arm.
    #
    # NB no na.rm here. If any comparison in an arm has a blank "Days grown", the
    # arm's min and max are both NA and every comparison in it is labelled
    # "Midpoint". Dropping the NAs would classify the rest of the arm instead, but
    # that changes the Timepoint column, its filter choices and therefore any saved
    # session key, so it is left alone.
    mutate(
      .dayMax = max(c(`Days grown`, `ControlDays grown`)),
      .dayMin = min(c(`Days grown`, `ControlDays grown`)),
      .by = c(`Experiment ID`, KO, `Cell line`, ControlKO, `ControlCell line`)
    ) %>%
    mutate(Timepoint = case_when(
      `Days grown` == .dayMax ~ "Final timepoint",
      `Days grown` == .dayMin ~ "Initial timepoint",
      TRUE                    ~ "Midpoint"
    )) %>%
    transmute(
      Citation = factor(Citation),
      `Experiment ID`,
      Contrast,
      Kind = factor(case_when(
        grepl("[Ee]ssentialome", Contrast)                       ~ "Essentialome",
        `ControlCell line` == "Plasmid" & `Cell line` != "Plasmid" ~ "Essentialome",
        Treatment != ControlTreatment                            ~ "Treatment",
        KO != ControlKO                                          ~ "Knockout",
        `Cell line` != `ControlCell line`                        ~ "Cell line",
        Dose != ControlDose                                      ~ "Dose",
        `Days grown` != `ControlDays grown`                      ~ "Essentialome",
        grepl("[Ss]orted", Contrast)                             ~ "Cell sorting",
        grepl("^CRISPR", Contrast)                               ~ "CRISPR targeting"
      )),
      Endpoint = factor(case_when(
        grepl(" on gH2A.X", Contrast) ~ sub("^.+ on (.+?)(, .+?| at .+?| in .+?| under .+?| from .+?| to .+?)*$", "\\1", Contrast),
        grepl(" on ", Contrast)       ~ toupper_first_initial(sub("^.+ on (.+?)(, .+?| at .+?| in .+?| under .+?| from .+?| to .+?)*$", "\\1", Contrast)),
        TRUE                          ~ "Default"
      )),
      Timepoint             = factor(Timepoint),
      `Days grown (diff)`   = factor(`Days grown`),
      `Dose (diff)`         = factor(Dose),
      `Treatment (diff)`    = factor(Treatment),
      `Knockout (diff)`     = factor(KO),
      `Cell line (diff)`    = factor(`Cell line`),
      `Days grown (ref)`    = factor(`ControlDays grown`),
      `Dose (ref)`          = factor(ControlDose),
      `Treatment (ref)`     = factor(ControlTreatment),
      `Knockout (ref)`      = factor(ControlKO),
      `Cell line (ref)`     = factor(`ControlCell line`),
      Library               = factor(Library),
      Organism              = factor(Organism),
      `Comparison ID`,
      FriendlyID            = paste0(stripHtml(Citation), ": ", Contrast),
      Source                = factor(ds$name)
    )

  list(experiments = experiments, comparisons = comparisons)
}

#' Filter a raw symbol list down to real gene symbols.
filterGeneSymbols <- function(symbols) {
  symbols <- symbols[!grepl(CRAVE_GENE_EXCLUDE, symbols)]
  symbols <- symbols[symbols != "X"]
  sort(unique(symbols))
}

#' Load one dataset declared in config.R.
#'
#' @param ds list(name, path, guides = FALSE, citation_from_id = FALSE)
#' @param create_indexes Whether to add CRAVE's indexes to the stat table if absent.
#' @return list or NULL if the dataset is missing.
loadDataset <- function(ds, create_indexes = TRUE) {
  path_dataset <- ds$path

  file_experiments <- file.path(path_dataset, "experiments_metadata.csv.gz")
  file_comparisons <- file.path(path_dataset, "comparisons_metadata.csv.gz")
  file_data        <- file.path(path_dataset, "database.db")
  file_ontology    <- file.path(path_dataset, "ontology.tsv.gz")
  file_guide       <- file.path(path_dataset, "guideLibrary.db")
  file_libs        <- file.path(path_dataset, "libraries.txt.gz")

  if (!file.exists(file_data)) {
    warning("Dataset ", file_data, " doesn't exist.")
    return(NULL)
  }

  log_info(paste0("Opening dataset '", ds$name, "' at ", file_data, "."))
  con <- DBI::dbConnect(RSQLite::SQLite(), file_data)
  ensureStatIndexes(con, label = ds$name, enabled = create_indexes)

  # -- Cached metadata wrangle -------------------------------------------------
  sources <- c(file_experiments, file_comparisons)
  if (file.exists(file_ontology)) sources <- c(sources, file_ontology)
  # The dataset's name and its citation_from_id flag are baked into the wrangled
  # output (as the Source column and as synthesised citations), so they have to be
  # part of the cache key. Otherwise renaming a dataset in config.R without
  # touching its files would reuse a cache carrying the old Source, and the
  # comparison-to-dataset routing in fetchStat() would silently find nothing.
  fingerprint <- paste0(fileFingerprint(sources), "|name=", ds$name,
                        "|cfi=", isTRUE(ds$citation_from_id))
  cacheFile <- file.path(
    craveCacheDir(path_dataset),
    paste0("meta-", gsub("[^A-Za-z0-9]+", "_", ds$name), ".rds")
  )

  cached <- NULL
  if (file.exists(cacheFile)) {
    cached <- tryCatch(readRDS(cacheFile), error = function(e) NULL)
    if (!identical(cached$fingerprint, fingerprint)) cached <- NULL
  }

  if (is.null(cached)) {
    log_info(paste0("Wrangling metadata for dataset '", ds$name, "'."))
    experiments <- readTableGlob(file_experiments, colClasses = "character")
    comparisons <- readTableGlob(file_comparisons, colClasses = "character")
    if (is.null(experiments) || is.null(comparisons)) {
      warning("Dataset ", ds$name, " is missing its metadata files.")
      DBI::dbDisconnect(con)
      return(NULL)
    }
    ontology <- if (file.exists(file_ontology)) {
      readTableGlob(file_ontology, colClasses = "character")
    } else {
      data.table()
    }
    wrangled <- tryCatch(wrangleMetadata(experiments, comparisons, ds),
                         error = function(e) {
                           warning(conditionMessage(e), call. = FALSE)
                           log_error(conditionMessage(e))
                           NULL
                         })
    if (is.null(wrangled)) {
      DBI::dbDisconnect(con)
      return(NULL)
    }
    cached <- list(
      fingerprint = fingerprint,
      experiments = wrangled$experiments,
      comparisons = wrangled$comparisons,
      ontology    = ontology
    )
    tryCatch(saveRDS(cached, cacheFile, compress = FALSE),
             error = function(e) log_info("Metadata cache is not writable; will re-wrangle next start."))
  } else {
    log_info(paste0("Reusing cached metadata for dataset '", ds$name, "'."))
  }

  # -- Gene symbols ------------------------------------------------------------
  symbols <- DBI::dbGetQuery(con, "SELECT symbol FROM gene")$symbol
  genes <- filterGeneSymbols(symbols)

  list(
    name         = ds$name,
    con          = con,
    experiments  = cached$experiments,
    comparisons  = cached$comparisons,
    ontology     = cached$ontology,
    genes        = genes,
    guideLibrary = if (isTRUE(ds$guides) && file.exists(file_guide)) file_guide else NULL,
    libraries    = if (isTRUE(ds$guides) && file.exists(file_libs))  file_libs  else NULL
  )
}

#' Probe whether the Exorcise Docker image and reference data are both present.
#'
#' Cheap check first: if the reference directory has no .2bit files, Exorcise
#' cannot run whatever Docker says, so skip the subprocess entirely.
detectExorcise <- function(root, image) {
  if (length(root) == 0 || length(image) == 0) return(FALSE)
  if (!nzchar(root) || !nzchar(image)) return(FALSE)
  if (!dir.exists(root)) return(FALSE)
  if (length(Sys.glob(file.path(root, "*.2bit"))) == 0) {
    log_info(paste0("No .2bit reference found under ", root, "; Exorcise disabled."))
    return(FALSE)
  }
  if (nzchar(Sys.which("docker")) == FALSE) {
    log_info("Docker CLI not found; Exorcise disabled.")
    return(FALSE)
  }
  status <- tryCatch(
    suppressWarnings(system2("docker", c("image", "inspect", shQuote(image)),
                             stdout = FALSE, stderr = FALSE)),
    error = function(e) 1L
  )
  ok <- identical(as.integer(status), 0L)
  if (!ok) log_info(paste0("Exorcise image '", image, "' not present locally; Exorcise disabled."))
  ok
}

#' Build the picker choice lists shared by the Explore and Correlate filters.
#'
#' These depend only on the dataset, so they are computed once at load time rather
#' than per session.
buildChoices <- function(experiments, comparisons) {
  cols <- c(
    Library             = "Library",
    Organism            = "Organism",
    Source              = "Source",
    Kind                = "Kind",
    Timepoint           = "Timepoint",
    Endpoint            = "Endpoint",
    `Days grown (diff)` = "Days grown (diff)",
    `Days grown (ref)`  = "Days grown (ref)",
    `Treatment (diff)`  = "Treatment (diff)",
    `Treatment (ref)`   = "Treatment (ref)",
    `Dose (diff)`       = "Dose (diff)",
    `Dose (ref)`        = "Dose (ref)",
    `Knockout (diff)`   = "Knockout (diff)",
    `Knockout (ref)`    = "Knockout (ref)",
    `Cell line (diff)`  = "Cell line (diff)",
    `Cell line (ref)`   = "Cell line (ref)"
  )
  out <- lapply(cols, function(cl) sort(unique(as.character(comparisons[[cl]]))))
  out$Citation           <- sort(unique(stripHtml(as.character(comparisons$Citation))))
  out$FriendlyID         <- comparisons$FriendlyID
  out$citationSummary    <- summariseItems(stripHtml(experiments$Citation))
  out
}

#' Load every dataset declared in config.R and assemble the app data object.
loadCraveData <- function(datasets, exorcise_root = NULL, exorcise_docker = NULL,
                          create_indexes = TRUE) {
  t0 <- Sys.time()

  loaded <- lapply(datasets, loadDataset, create_indexes = create_indexes)
  loaded <- loaded[!vapply(loaded, is.null, logical(1))]

  if (length(loaded) == 0) {
    return(list(
      cons = list(), experiments = tibble(), comparisons = tibble(),
      ontology = data.table(), genes = character(0),
      comparisonSource = character(0), statTemplate = NULL,
      files = list(), choices = list(),
      caps = list(explore = FALSE, guides = FALSE, libraries = FALSE,
                  ontology = FALSE, exorcise = FALSE)
    ))
  }

  names(loaded) <- vapply(loaded, `[[`, character(1), "name")

  cons        <- lapply(loaded, `[[`, "con")
  experiments <- bind_rows(lapply(loaded, `[[`, "experiments"))
  comparisons <- bind_rows(lapply(loaded, `[[`, "comparisons"))
  ontology    <- unique(bind_rows(lapply(loaded, `[[`, "ontology")))
  genes       <- sort(unique(unlist(lapply(loaded, `[[`, "genes"), use.names = FALSE)))

  # comparison_id -> dataset name, so fetchStat() can route a comparison to the
  # database that holds it in one indexed lookup.
  comparisonSource <- setNames(
    as.character(comparisons$Source),
    comparisons$`Comparison ID`
  )

  # Column template for the stat table, fetched once instead of once per query.
  statTemplate <- tryCatch(
    DBI::dbGetQuery(cons[[1]], "SELECT * FROM stat LIMIT 0"),
    error = function(e) NULL
  )

  guideFile <- Filter(Negate(is.null), lapply(loaded, `[[`, "guideLibrary"))
  libsFile  <- Filter(Negate(is.null), lapply(loaded, `[[`, "libraries"))

  out <- list(
    cons             = cons,
    experiments      = experiments,
    comparisons      = comparisons,
    ontology         = ontology,
    genes            = genes,
    comparisonSource = comparisonSource,
    statTemplate     = statTemplate,
    files = list(
      guideLibrary = if (length(guideFile)) guideFile[[1]] else NULL,
      libraries    = if (length(libsFile))  libsFile[[1]]  else NULL
    ),
    caps = list(
      explore   = nrow(experiments) > 0,
      guides    = length(guideFile) > 0,
      libraries = length(libsFile) > 0,
      ontology  = nrow(ontology) > 0,
      exorcise  = detectExorcise(exorcise_root, exorcise_docker)
    )
  )
  out$choices <- buildChoices(experiments, comparisons)

  log_info(paste0(
    "Loaded ", length(cons), " dataset(s), ", nrow(comparisons), " comparisons, ",
    length(genes), " gene symbols in ",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2), "s."
  ))
  out
}

#' Close every open dataset connection.
closeCraveData <- function(data) {
  for (con in data$cons) tryCatch(DBI::dbDisconnect(con), error = function(e) NULL)
  invisible(NULL)
}
