#### Exorcise ####
# Runs the Exorcise container over a set of guide sequences and returns plain data;
# the module decides how to render it.
#
# The Docker invocation is an argument vector passed to system2(), never a string
# handed to a shell, and the PAM and chemistry fields are validated against a
# character allow-list. Keep it that way: both are free-text user input.

#' Validate a PAM string.
validPam <- function(pam) {
  pam <- toupper(trimws(pam %||% ""))
  if (!nzchar(pam)) return("NGG")
  if (!grepl("^[ACGTN]+$", pam)) return(NULL)
  pam
}

#' Validate an Exorcise chemistry mode string.
#'
#' Accepts either a plain distance in nucleotides, or a base-editor spec of the
#' form beNM0123.
validModeString <- function(x) {
  x <- trimws(x %||% "")
  if (!nzchar(x)) return(NULL)
  if (grepl("^[0-9]{1,6}$", x)) return(x)
  if (grepl("^be[ACGT]{2}[0-9]{4}$", x)) return(x)
  NA_character_  # present but invalid
}

#' Run Exorcise over a set of guide sequences.
#'
#' @return list(status = "ok"|"error", message, input, result, taskId)
runExorcise <- function(params, config) {
  root     <- config$exorcise_root
  hostroot <- config$exorcise_hostroot
  image    <- config$exorcise_docker

  fail <- function(msg) list(status = "error", message = msg)

  tempDir <- file.path(root, "temp")
  if (!dir.exists(tempDir) && !dir.create(tempDir, showWarnings = FALSE, recursive = TRUE)) {
    return(fail("The Exorcise working directory is not writable. Please contact the administrator."))
  }

  # -- Gather input ------------------------------------------------------------
  if (is.null(params$upload)) {
    seqs  <- splitTokens(params$seq)
    origs <- splitTokens(params$orig)
    control <- params$control_string
  } else {
    crispick <- tryCatch(data.table::fread(params$upload$datapath[1]),
                         error = function(e) NULL)
    if (is.null(crispick) ||
        !all(c("sgRNA Sequence", "Input") %in% names(crispick))) {
      return(fail("That does not look like a CRISPick sgRNA designs file. It needs 'sgRNA Sequence' and 'Input' columns."))
    }
    seqs  <- as.character(crispick$`sgRNA Sequence`)
    origs <- as.character(crispick$Input)
    control <- "\\(NEG_CONTROL\\)"
  }

  # -- Validate ----------------------------------------------------------------
  pam <- validPam(params$pam)
  if (is.null(pam)) {
    return(fail("The PAM must contain only the characters A, C, G, T and N."))
  }
  if (length(seqs) == 0) {
    return(fail("No sequences submitted. Please revise your input."))
  }
  if (length(origs) > 0 && length(origs) != length(seqs)) {
    return(fail(paste0(
      "Submitted ", length(seqs), " sequences but ", length(origs),
      " annotations. Please supply one annotation per sequence, or none at all.")))
  }
  if (any(!grepl("^[ACGTacgt]+$", seqs))) {
    return(fail("Sequences must contain only A, C, G and T. Please revise your input."))
  }
  if (any(grepl("A{17}|C{17}|G{17}|T{17}", toupper(seqs)))) {
    return(fail("Mononucleotide run of 17 or more detected, which is not allowed. Please revise your input."))
  }
  if (any(nchar(seqs) + nchar(gsub("[^ACGTN]", "", pam)) < 21)) {
    return(fail("Short sequence detected. Sequence plus PAM must be at least 21 nucleotides. Please revise your input or PAM."))
  }

  advanced <- validModeString(params$mode_advanced)
  if (identical(advanced, NA_character_)) {
    return(fail("The advanced chemistry string must be a number (e.g. 1000) or a base-editor spec (e.g. beCA0110)."))
  }
  mode <- advanced %||% unname(EXORCISE_MODES[params$mode %||% "Knockout"])
  if (is.null(mode) || is.na(mode)) {
    return(fail("Unrecognised CRISPR chemistry. Please choose one from the list."))
  }

  genome <- unname(EXORCISE_GENOMES[params$genome %||% ""])
  exome  <- unname(EXORCISE_EXOMES[params$genome %||% ""])
  if (is.na(genome) || is.na(exome)) {
    return(fail("Unrecognised reference genome. Please choose one from the list."))
  }

  # -- Write the input file ----------------------------------------------------
  # Task ids name files in a directory shared by every concurrent user, so a
  # collision means two runs overwrite and then delete each other's data. The id
  # therefore does not rely on the RNG alone: it also carries the process id and a
  # microsecond timestamp, which are distinct even if the global random stream has
  # been seeded identically somewhere else.
  taskId <- paste0(
    format(Sys.getpid()), "-",
    format(as.numeric(Sys.time()) * 1e6, scientific = FALSE, trim = TRUE), "-",
    paste0(sample(c(LETTERS[1:6], 0:9), 6, replace = TRUE), collapse = "")
  )
  harmonise <- length(origs) == length(seqs) && length(origs) > 0
  exorciseInput <- if (harmonise) {
    data.table(seq = seqs, orig = origs)
  } else {
    data.table(seq = seqs)
  }

  infile  <- paste0("exorcise-input_", taskId, ".tsv")
  outdir  <- paste0("exorcise-output_", taskId)
  inPath  <- file.path(tempDir, infile)
  outPath <- file.path(tempDir, outdir)
  on.exit({
    suppressWarnings(file.remove(inPath))
    unlink(outPath, recursive = TRUE)
  }, add = TRUE)

  data.table::fwrite(exorciseInput, inPath, sep = "\t")

  # -- Run ---------------------------------------------------------------------
  args <- c(
    "run", "--rm",
    "-v", paste0(hostroot, ":/data/"),
    "-v", paste0(file.path(hostroot, "temp"), ":/temp/"),
    image, "exorcise",
    "-i", paste0("/temp/", infile),
    "-o", paste0("/temp/", outdir),
    "-z", pam,
    "-q", mode,
    "-g", "1",
    "-v", paste0("/data/", genome),
    "-w", paste0("/data/", exome)
  )
  if (harmonise) args <- c(args, "-n", "2")
  if (!is.null(control) && nzchar(control)) args <- c(args, "-c", control)

  log_info(paste0("Running Exorcise task ", taskId, " (", length(seqs), " sequences)."))
  status <- tryCatch(
    suppressWarnings(system2("docker", args, stdout = TRUE, stderr = TRUE)),
    error = function(e) {
      log_error(paste0("Exorcise failed to start: ", conditionMessage(e)))
      NULL
    }
  )

  outfile <- file.path(outPath, "exorcise.tsv")
  if (!file.exists(outfile)) {
    logs <- Sys.glob(file.path(outPath, "logfile*.log"))
    msg <- if (length(logs)) {
      tail <- tryCatch(readLines(logs[1], warn = FALSE), error = function(e) character(0))
      if (length(tail)) tail[length(tail)] else "Exorcise produced no output."
    } else if (length(status)) {
      paste(utils::tail(status, 3), collapse = " ")
    } else {
      "Exorcise produced no output."
    }
    log_warn(paste0("Exorcise task ", taskId, " failed: ", msg))
    return(list(status = "error", message = msg, input = exorciseInput, taskId = taskId))
  }

  result <- data.table::fread(outfile)
  log_info(paste0("Exorcise task ", taskId, " done, ", nrow(result), " rows."))
  list(
    status  = "ok",
    message = "Exorcise done, results are ready for download.",
    input   = exorciseInput,
    result  = result,
    taskId  = taskId
  )
}
