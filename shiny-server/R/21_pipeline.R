#### Dataset mutation: driving the Exorcise pipeline ####
#
# Everything here exists to keep one promise: database.db and the two metadata
# CSVs beside it always describe the same experiments.
#
# exorcise-database is the only program that writes any of the three, and it
# writes all three in one pass -- add_statistics(), commit, metadata.to_files() --
# so the rule CRAVE follows is that nothing else ever touches them. Experiments are
# added by running exorcise-analyse and then exorcise-database; they are removed
# with exorcise-database --remove; and the one mutation CRAVE performs itself,
# replacing the metadata tables with edited copies, is refused unless the incoming
# tables name exactly the experiments and comparisons the database already holds.
# That is checked in replaceMetadataTables() and is why it cannot desynchronise
# anything: it can relabel rows, and nothing else.
#
# Both long-running commands are run detached, through a shell script, rather than
# with system2(wait = TRUE). A MAGeCK or Chronos run is minutes to hours, and R is
# single threaded: waiting on it inside an observer would freeze every other
# session connected to the app, not just the one that pressed the button.

#' Keep a picker's current value if it is still one of the choices.
#'
#' Not `%||%`: an un-selected selectizeInput reports "", which is neither NULL nor
#' length zero, so `%||%` keeps it and the picker comes up empty with every choice
#' available but none chosen.
keepOrFirst <- function(current, choices) {
  if (length(choices) == 0) return(character(0))
  if (length(current) == 1 && !is.na(current) && nzchar(current) &&
      current %in% choices) current else choices[1]
}

#### Where a dataset keeps what has been submitted into it ####

#' The submission area inside a dataset directory.
#'
#' Results and workbooks are kept, not discarded, because they are what lets the
#' Manage tab rebuild an experiment's metadata later without re-running its
#' analyses: exorcise-database reads result tables from <results>/<id>/tables/,
#' which is exactly the layout exorcise-analyse writes.
submitPaths <- function(datasetPath) {
  root <- file.path(datasetPath, CRAVE_SUBMIT_DIR)
  list(
    root      = root,
    counts    = file.path(root, "counts"),
    workbooks = file.path(root, "workbooks"),
    results   = file.path(root, "results"),
    jobs      = file.path(root, "jobs"),
    lock      = file.path(root, "lock")
  )
}

#' Create the submission area, returning FALSE if the dataset is read-only.
ensureSubmitDirs <- function(datasetPath) {
  p <- submitPaths(datasetPath)
  for (d in c(p$root, p$counts, p$workbooks, p$results, p$jobs)) {
    if (!dir.exists(d)) {
      suppressWarnings(dir.create(d, recursive = TRUE, showWarnings = FALSE))
    }
    if (!dir.exists(d)) return(FALSE)
  }
  file.access(p$root, 2) == 0
}

#' TRUE if a dataset directory can be written to.
#'
#' Submit and Manage both rewrite the dataset in place, so a read-only mount
#' disqualifies a dataset from either regardless of what config.R allows.
datasetWritable <- function(path) {
  dir.exists(path) && file.access(path, 2) == 0
}

#' The path a file inside the dataset has when the dataset is mounted at /dataset.
#'
#' Every Exorcise invocation mounts the dataset directory and nothing else, so each
#' argument is this translation of a path CRAVE already knows. Deriving them rather
#' than pasting strings is what keeps a user-supplied filename from becoming a path.
containerPath <- function(datasetPath, path) {
  base <- sub("/+$", "", normalizePath(datasetPath, mustWork = FALSE))
  full <- normalizePath(path, mustWork = FALSE)
  if (!startsWith(full, base)) {
    stop("Refusing to name a path outside the dataset: ", path, call. = FALSE)
  }
  paste0("/dataset", substring(full, nchar(base) + 1))
}

#' The host path of a dataset directory.
#'
#' CRAVE reaches the Docker daemon through its socket, so every -v argument is
#' resolved by the daemon on the host, not inside CRAVE's own container. A dataset
#' bind-mounted at /data/dataset-1 from /opt/crave/dataset-1 has to be named by the
#' latter when the Exorcise container is handed it, exactly as exorcise_hostroot
#' already does for the Exorcise tab. A deployment running CRAVE natively sets no
#' hostpath, and the two are the same path.
datasetHostPath <- function(ds) {
  sub("/+$", "", ds$hostpath %||% ds$path)
}

#### Validating what was uploaded ####

#' Reduce a name from a workbook to something safe to use as a filename.
#'
#' The counts filename comes out of a cell the submitter typed, and is then used to
#' name a file inside the dataset. Basename and an allow-list, so it cannot climb
#' out of the directory or mean anything to a shell.
safeFilename <- function(x) {
  x <- basename(as.character(x %||% ""))
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- sub("^[._]+", "", x)
  if (!nzchar(x)) NULL else x
}

#' Read what CRAVE needs to know out of an analysis workbook.
#'
#' Advisory, not authoritative: exorcise-analyse parses the workbook again and its
#' reading is the one that counts. This exists so that an unusable workbook is
#' refused in a second, with a message naming what is wrong, instead of after a
#' container has started. The rules mirror exorcise_screens/workbook.py -- the four
#' required sheets, "Analysis name" as the legacy spelling of "Experiment name",
#' and comma-separated method and group cells.
#'
#' @return list(ok, experimentId, methods, countsFile, groups, errors)
readWorkbook <- function(path) {
  fail <- function(...) list(ok = FALSE, errors = c(...))
  if (!needPkg("readxl")) {
    return(fail("The readxl package is not installed, so workbooks cannot be read."))
  }
  sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) NULL)
  if (is.null(sheets)) {
    return(fail("That file could not be opened as an .xlsx workbook."))
  }
  missing <- setdiff(CRAVE_WORKBOOK_SHEETS, sheets)
  if (length(missing)) {
    return(fail(paste0(
      "The workbook is missing the sheet(s): ", paste(missing, collapse = ", "),
      ". An analysis workbook needs: ", paste(CRAVE_WORKBOOK_SHEETS, collapse = ", "), ".")))
  }

  details <- tryCatch(
    readxl::read_excel(path, sheet = "Experiment details", col_names = FALSE,
                       .name_repair = "minimal"),
    error = function(e) NULL)
  if (is.null(details) || ncol(details) < 2) {
    return(fail("The \"Experiment details\" sheet needs two columns: field name and value."))
  }
  field <- trimws(as.character(details[[1]]))
  value <- trimws(as.character(details[[2]]))
  pick <- function(label) {
    hit <- which(field == label & !is.na(value) & nzchar(value) & value != "NA")
    if (length(hit)) value[hit[1]] else NA_character_
  }
  experimentId <- pick("Experiment name") %||% NA_character_
  if (is.na(experimentId)) experimentId <- pick("Analysis name")
  if (is.na(experimentId) || !nzchar(experimentId)) {
    return(fail(paste0(
      "The workbook has no \"Experiment name\" in its \"Experiment details\" sheet. ",
      "That field is the experiment's identity in the database.")))
  }

  analyses <- tryCatch(readxl::read_excel(path, sheet = "Analyses"),
                       error = function(e) NULL)
  if (is.null(analyses) || !"Method" %in% names(analyses)) {
    return(fail("The \"Analyses\" sheet needs at least a \"Method\" column."))
  }
  cellValues <- function(column) {
    if (!column %in% names(analyses)) return(character(0))
    raw <- as.character(analyses[[column]])
    raw <- raw[!is.na(raw) & nzchar(trimws(raw))]
    sort(unique(trimws(unlist(strsplit(raw, ",")))))
  }
  methods <- tolower(cellValues("Method"))
  unknown <- setdiff(methods, CRAVE_ANALYSIS_METHODS)
  if (length(unknown)) {
    return(fail(paste0(
      "The \"Analyses\" sheet asks for method(s) CRAVE cannot run: ",
      paste(unknown, collapse = ", "), ". Available: ",
      paste(CRAVE_ANALYSIS_METHODS, collapse = ", "), ".")))
  }
  if (length(methods) == 0) {
    return(fail("The \"Analyses\" sheet names no methods, so there is nothing to run."))
  }

  countsFiles <- cellValues("Counts file")
  needsCounts <- length(intersect(methods, CRAVE_COUNTS_METHODS)) > 0
  if (needsCounts && length(countsFiles) != 1) {
    return(fail(paste0(
      "The workbook names ", length(countsFiles), " counts files (",
      paste(countsFiles, collapse = ", ") , "). MAGeCK, DrugZ and Chronos read one ",
      "counts table, so exactly one is supported.")))
  }

  # "Analysis version" moves the output under results/<id>/<version>/tables, which
  # is not where exorcise-database looks. Refusing is better than building an
  # experiment with no statistics in it.
  version <- pick("Analysis version")
  if (!is.na(version) && nzchar(version)) {
    return(fail(paste0(
      "The workbook sets \"Analysis version\" to '", version, "'. exorcise-analyse ",
      "then files its tables under results/<experiment>/", version, "/, where ",
      "exorcise-database does not look for them. Clear that field and resubmit.")))
  }

  list(
    ok           = TRUE,
    experimentId = experimentId,
    methods      = methods,
    countsFile   = if (length(countsFiles)) safeFilename(countsFiles[1]) else NULL,
    needsCounts  = needsCounts,
    groups       = cellValues("Control group"),
    errors       = character(0)
  )
}

#' Check that an uploaded counts matrix looks like one.
#'
#' Two columns at least, the first unique guide ids and the second gene symbols,
#' which is what example-analysis/commands.txt documents. Only the header and
#' the first rows are read: a counts matrix is tens of megabytes and the point is a
#' fast refusal, not a full validation.
readCountsHead <- function(path, n = 20) {
  head <- tryCatch(
    data.table::fread(file = path, nrows = n, header = TRUE, encoding = "UTF-8"),
    error = function(e) NULL)
  if (is.null(head)) {
    return(list(ok = FALSE, errors = "That file could not be read as a table."))
  }
  if (ncol(head) < 3) {
    return(list(ok = FALSE, errors = paste0(
      "A counts matrix needs guide ids, gene symbols and at least one sample ",
      "column; this has ", ncol(head), " column(s).")))
  }
  list(ok = TRUE, columns = names(head), samples = names(head)[-(1:2)],
       gzipped = isGzipFile(path), errors = character(0))
}

#### Staging a counts matrix ####
# fread() sniffs a gzip whatever the file is called, so CRAVE reads an upload
# either way. Exorcise does not: it resolves the counts file by the name the
# workbook declares and decides from that name's extension whether to decompress,
# so the name it is saved under has to agree with the bytes in it.

#' TRUE if this file begins with the gzip magic number.
#'
#' By content, not by name: the upload's own name is only what the browser sent,
#' and a file called .tsv is routinely a gzip.
isGzipFile <- function(path) {
  # suppressWarnings as well as tryCatch: file() warns rather than raises on a
  # directory or a missing path, and this is a probe, not a demand.
  con <- suppressWarnings(tryCatch(file(path, "rb"), error = function(e) NULL))
  if (is.null(con)) return(FALSE)
  on.exit(close(con), add = TRUE)
  magic <- tryCatch(readBin(con, "raw", n = 2L), error = function(e) raw(0))
  length(magic) == 2L && magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b)
}

#' Copy a file, compressing it on the way.
#'
#' In chunks through a connection rather than readBin(n = file.size()): a counts
#' matrix is tens of megabytes uncompressed and there is no reason to hold all of
#' it at once.
gzipCopy <- function(from, to, chunk = 1024L * 1024L) {
  input <- file(from, "rb")
  on.exit(close(input), add = TRUE)
  output <- gzfile(to, "wb")
  on.exit(close(output), add = TRUE)
  repeat {
    bytes <- readBin(input, "raw", n = chunk)
    if (length(bytes) == 0) break
    writeBin(bytes, output)
  }
  TRUE
}

#' The name a counts matrix has to be saved under for Exorcise to read it.
#'
#' Exorcise looks for the workbook's declared name, and maybe_its_gz() will reach
#' for a .gz the workbook did not mention -- but it never strips one. So a gzipped
#' upload is saved with .gz appended to the declared name and is found by that
#' fallback, while a plain upload declared as .gz has to be compressed on the way
#' in, because nothing would ever look for it uncompressed.
#'
#' @return list(name, compress), compress TRUE when the copy must be gzipped.
countsTargetName <- function(declared, gzipped) {
  declaredGz <- grepl("\\.gz$", declared, ignore.case = TRUE)
  if (isTRUE(gzipped) && !declaredGz) return(list(name = paste0(declared, ".gz"),
                                                  compress = FALSE))
  if (!isTRUE(gzipped) && declaredGz)  return(list(name = declared, compress = TRUE))
  list(name = declared, compress = FALSE)
}

#' Save an uploaded counts matrix into the dataset under a name Exorcise resolves.
#'
#' The counterpart name is removed when it is there. maybe_its_gz() prefers the
#' plain name, so a file left behind by an earlier submission -- including one
#' saved wrongly by a CRAVE before 6.0.0 -- would be picked ahead of the one just
#' written, and the submission would keep failing on a file the user had replaced.
stageCountsFile <- function(countsPath, countsDir, declared) {
  target <- countsTargetName(declared, isGzipFile(countsPath))
  path <- file.path(countsDir, target$name)
  other <- file.path(countsDir, if (identical(target$name, declared))
    paste0(declared, ".gz") else declared)

  ok <- tryCatch({
    if (isTRUE(target$compress)) gzipCopy(countsPath, path)
    else file.copy(countsPath, path, overwrite = TRUE)
  }, error = function(e) {
    log_warn(paste0("Could not stage the counts matrix: ", conditionMessage(e)))
    FALSE
  })
  if (!isTRUE(ok)) {
    return(list(ok = FALSE, message = "The counts matrix could not be saved into the dataset."))
  }

  if (file.exists(other) && !identical(normalizePath(other, mustWork = FALSE),
                                       normalizePath(path, mustWork = FALSE))) {
    unlink(other)
    log_info(paste0("Removed '", basename(other), "', which Exorcise would have ",
                    "resolved ahead of the counts matrix just submitted."))
  }
  if (!identical(target$name, declared)) {
    log_info(paste0("The workbook names '", declared, "'; the upload is ",
                    if (isTRUE(target$compress)) "not " else "",
                    "gzipped, so it was saved as '", target$name, "'."))
  }
  list(ok = TRUE, name = target$name)
}

#### Raw analysis output ####
# exorcise-analyse leaves two different things under results/<experiment>/. Each
# method writes its own output into a directory named after itself -- MAGeCK's
# sgRNA and gene summaries, DrugZ's per-comparison files, Chronos's model output
# -- and the tabulator then reduces all of that to one CSV per method under
# tables/, which is what exorcise-database reads.
#
# Only the method directories are offered for download. The tables are already in
# the database, and therefore behind every plot and every "Download data" in
# CRAVE, so putting them in the archive would hand back the same numbers a second
# time under a different name. The raw output is the part that is otherwise
# unreachable once the container has exited.

#' The methods that left raw output for one experiment.
#'
#' CRAVE_COUNTS_METHODS doubles as the list of directory names because
#' exorcise-analyse names each one after the method that wrote it. "manual" is
#' not among them -- a supplied table is not computed, so it leaves nothing --
#' and an experiment analysed that way has no raw output to offer.
rawResultMethods <- function(datasetPath, experimentId) {
  dirs <- file.path(submitPaths(datasetPath)$results, experimentId,
                    CRAVE_COUNTS_METHODS)
  CRAVE_COUNTS_METHODS[dir.exists(dirs)]
}

#' The experiments in a dataset that have raw output kept for them.
#'
#' Listed from the directory rather than derived from the metadata, for two
#' reasons: the folder name is whatever exorcise-analyse wrote, which is the
#' experiment id exactly as the workbook spelled it; and an id read off the disk
#' is one that exists, so nothing a submitter typed is ever used to build a path.
#' That is what makes rawResultsArchive() safe -- it only accepts a name from
#' this list.
rawResultExperiments <- function(datasetPath) {
  dir <- submitPaths(datasetPath)$results
  if (!dir.exists(dir)) return(character(0))
  ids <- list.dirs(dir, full.names = FALSE, recursive = FALSE)
  keep <- vapply(ids, function(id) length(rawResultMethods(datasetPath, id)) > 0,
                 logical(1))
  sort(ids[keep])
}

#' What an experiment's archive would hold.
#'
#' The uncompressed size, because it is what can be known without building the
#' archive, and because the point of reporting it is to warn before a long
#' download rather than to be exact about the gzip.
#'
#' @return list(methods, files, bytes)
rawResultSummary <- function(datasetPath, experimentId) {
  methods <- rawResultMethods(datasetPath, experimentId)
  dirs <- file.path(submitPaths(datasetPath)$results, experimentId, methods)
  files <- unlist(lapply(dirs, list.files, recursive = TRUE, full.names = TRUE,
                         all.files = TRUE))
  list(
    methods = methods,
    files   = length(files),
    bytes   = if (length(files)) sum(file.size(files), na.rm = TRUE) else 0
  )
}

#' Write one experiment's raw output to a .tar.gz.
#'
#' R's internal tar rather than the system one, so this works wherever CRAVE does
#' and does not depend on a tar binary being in the image.
#'
#' The paths are stored relative, so the archive unpacks as
#' <experiment>/<method>/... instead of recreating the server's directory tree on
#' whoever downloaded it. tar() stores the names it is handed, so that means
#' moving the working directory for the length of the call. R is single threaded
#' and tar() never returns to the event loop, so no other session can observe the
#' move; the target is resolved first, because it is a relative path once we have.
#'
#' @return TRUE if an archive was written.
rawResultsArchive <- function(datasetPath, experimentId, target) {
  if (!isTRUE(experimentId %in% rawResultExperiments(datasetPath))) return(FALSE)
  methods <- rawResultMethods(datasetPath, experimentId)
  target <- normalizePath(target, mustWork = FALSE)

  old <- setwd(submitPaths(datasetPath)$results)
  on.exit(setwd(old), add = TRUE)
  ok <- tryCatch(
    withCallingHandlers(
      identical(utils::tar(target, file.path(experimentId, methods),
                           compression = "gzip", tar = "internal"), 0L),
      warning = function(w) {
        # One warning is muffled and the rest are left to be logged: a tar header
        # stores the owner's id in an octal field too narrow for the ids macOS
        # hands out, so R substitutes nobody's and says so. It changes nothing
        # about the files in the archive, and it would otherwise be logged on
        # every download a native deployment served.
        if (grepl("invalid [ug]id value", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }),
    error = function(e) {
      log_error(paste0("Could not archive the raw output of '", experimentId,
                       "': ", conditionMessage(e)))
      FALSE
    })
  isTRUE(ok)
}

#### Jobs ####

#' A job id that is unique even if two sessions submit in the same second.
#'
#' Same reasoning as the Exorcise task ids: process id and a microsecond timestamp,
#' so it does not depend on the state of the random stream.
newJobId <- function(prefix = "job") {
  paste0(prefix, "-", format(Sys.time(), "%Y%m%d-%H%M%S"), "-",
         format(Sys.getpid()), "-",
         paste0(sample(c(LETTERS[1:6], 0:9), 4, replace = TRUE), collapse = ""))
}

#' Assemble the `docker run` argument vector for one Exorcise command.
#'
#' An argument vector, never a command string: every element is passed to the
#' process as-is, so a filename cannot become an argument or a shell operator.
exorciseDockerArgs <- function(ds, image, command, args) {
  c("run", "--rm",
    "-v", paste0(datasetHostPath(ds), ":/dataset"),
    image, command, args)
}

#' Write the shell script a job runs.
#'
#' One script with several docker invocations in it, rather than several detached
#' R processes, so that the steps are ordered, a failure stops the rest, and the
#' lock is released exactly once however the script ends.
writeJobScript <- function(jobDir, lockFile, commands) {
  script <- file.path(jobDir, "job.sh")
  quoted <- vapply(commands, function(argv) {
    paste(c("docker", shQuote(argv)), collapse = " ")
  }, character(1))
  writeLines(c(
    "#!/bin/sh",
    "# Written by CRAVE. One job: the lock is released by the trap, whatever happens.",
    "set -eu",
    paste0("JOB_DIR=", shQuote(jobDir)),
    paste0("LOCK=", shQuote(lockFile)),
    # The script redirects its own output rather than letting R do it: system2()
    # with stdout and stderr naming one file opens that file twice, and the two
    # streams then overwrite each other's lines.
    "exec >> \"$JOB_DIR/job.log\" 2>&1",
    "echo $$ > \"$JOB_DIR/pid\"",
    "printf 'running' > \"$JOB_DIR/status\"",
    "cleanup() {",
    "  code=$?",
    "  if [ \"$code\" -eq 0 ]; then printf 'ok' > \"$JOB_DIR/status\";",
    "  else printf 'failed:%s' \"$code\" > \"$JOB_DIR/status\"; fi",
    "  rm -f \"$LOCK\"",
    "}",
    "trap cleanup EXIT",
    "",
    quoted,
    ""
  ), script)
  Sys.chmod(script, "0755")
  script
}

#' Is a mutation already running against this dataset?
#'
#' The lock is a file naming the job that holds it. A job that died without its
#' trap running would otherwise hold the lock forever, so a lock whose job has a
#' terminal status, or whose process is gone, is reported as stale rather than
#' busy -- the Manage tab offers to clear those.
datasetLock <- function(datasetPath) {
  p <- submitPaths(datasetPath)
  if (!file.exists(p$lock)) return(NULL)
  jobId <- tryCatch(trimws(readLines(p$lock, warn = FALSE)[1]),
                    error = function(e) NA_character_)
  if (is.na(jobId) || !nzchar(jobId)) return(list(jobId = NA_character_, stale = TRUE))
  st <- jobStatus(datasetPath, jobId)
  list(jobId = jobId, stale = !identical(st$state, "running"), status = st)
}

#' Clear a lock left behind by a job that is no longer running.
clearDatasetLock <- function(datasetPath) {
  suppressWarnings(file.remove(submitPaths(datasetPath)$lock))
  invisible(NULL)
}

#' What a job is doing, and what it has said so far.
#'
#' @return list(state, exit, log, started) with state one of "running", "ok",
#'   "failed", "missing".
jobStatus <- function(datasetPath, jobId, logLines = CRAVE_JOB_LOG_LINES) {
  if (is.null(jobId) || is.na(jobId) || !nzchar(jobId)) {
    return(list(state = "missing", log = character(0)))
  }
  jobDir <- file.path(submitPaths(datasetPath)$jobs, basename(jobId))
  if (!dir.exists(jobDir)) return(list(state = "missing", log = character(0)))

  raw <- tryCatch(readLines(file.path(jobDir, "status"), warn = FALSE),
                  error = function(e) character(0))
  raw <- if (length(raw)) trimws(raw[1]) else ""
  log <- tryCatch(readLines(file.path(jobDir, "job.log"), warn = FALSE),
                  error = function(e) character(0))
  if (length(log) > logLines) log <- c("[...]", utils::tail(log, logLines))

  state <- if (startsWith(raw, "failed")) "failed"
           else if (identical(raw, "ok")) "ok"
           else if (identical(raw, "running")) "running"
           else "missing"

  # A script killed outright leaves "running" behind with no process to match it.
  if (identical(state, "running") && !jobProcessAlive(jobDir)) state <- "failed"

  list(
    state   = state,
    exit    = if (startsWith(raw, "failed:")) sub("^failed:", "", raw) else NA_character_,
    log     = log,
    jobDir  = jobDir,
    started = tryCatch(file.info(jobDir)$mtime, error = function(e) NA)
  )
}

#' TRUE while the shell that runs a job still exists.
#'
#' The pid file is written by the script itself, so for the first moment of a job
#' there is nothing to check. That window is treated as alive, but only for
#' CRAVE_JOB_START_GRACE seconds: a script that never ran at all would otherwise
#' leave a job that is reported as running for as long as the app is up.
jobProcessAlive <- function(jobDir) {
  pidFile <- file.path(jobDir, "pid")
  if (!file.exists(pidFile)) {
    age <- difftime(Sys.time(), file.info(jobDir)$mtime, units = "secs")
    return(!is.na(age) && age < CRAVE_JOB_START_GRACE)
  }
  pid <- tryCatch(suppressWarnings(as.integer(readLines(pidFile, warn = FALSE)[1])),
                  error = function(e) NA_integer_)
  if (is.na(pid)) return(TRUE)  # being written as we read it
  identical(0L, suppressWarnings(system2("kill", c("-0", pid),
                                         stdout = FALSE, stderr = FALSE)))
}

#' Start a job: take the lock, write the script, run it detached.
#'
#' @param commands List of docker argument vectors, run in order.
#' @return list(ok, jobId, message)
startJob <- function(ds, commands, label = "job") {
  if (!ensureSubmitDirs(ds$path)) {
    return(list(ok = FALSE, message = paste0(
      "The dataset directory '", ds$path, "' is not writable, so nothing can be ",
      "submitted into it.")))
  }
  p <- submitPaths(ds$path)
  held <- datasetLock(ds$path)
  if (!is.null(held) && !isTRUE(held$stale)) {
    return(list(ok = FALSE, message = paste0(
      "Another job (", held$jobId, ") is already running against this dataset. ",
      "One at a time, so that two of them cannot rewrite the database at once.")))
  }
  if (!is.null(held)) clearDatasetLock(ds$path)

  jobId <- newJobId(label)
  jobDir <- file.path(p$jobs, jobId)
  dir.create(jobDir, recursive = TRUE, showWarnings = FALSE)
  writeLines(jobId, p$lock)
  # Written here, not left to the script: between system2() returning and the
  # shell reaching its first line there is a window in which the job directory
  # exists with no status in it, and a caller polling in that window would be told
  # the job had finished. The script writes the same value again when it starts.
  cat("running", file = file.path(jobDir, "status"))

  script <- writeJobScript(jobDir, p$lock, commands)
  logFile <- file.path(jobDir, "job.log")
  writeLines(c(paste0("CRAVE ", CRAVE_VERSION, " starting ", label, " ", jobId,
                      " at ", format(Sys.time())), ""), logFile)

  log_info(paste0("Starting ", label, " ", jobId, " on dataset '", ds$name, "'."))
  started <- tryCatch({
    system2("sh", shQuote(script), wait = FALSE)
    TRUE
  }, error = function(e) {
    log_error(paste0("Could not start job ", jobId, ": ", conditionMessage(e)))
    FALSE
  })
  if (!started) {
    clearDatasetLock(ds$path)
    return(list(ok = FALSE, message = "The job could not be started. See the log."))
  }
  list(ok = TRUE, jobId = jobId, message = paste0("Started ", jobId, "."))
}

#### The two pipeline jobs ####

#' Add an experiment: run the analyses, then put the results in the database.
#'
#' The two commands are one job because the second is meaningless without the
#' first, and because a half-finished submission must not leave the lock held.
#'
#' @param replace TRUE to pass --update-existing, which replaces an experiment
#'   already in the database rather than skipping it.
submitExperimentJob <- function(ds, config, workbookPath, replace = FALSE) {
  p <- submitPaths(ds$path)
  image <- config$exorcise_docker

  analyse <- exorciseDockerArgs(ds, image, "exorcise-analyse", c(
    "--counts",     containerPath(ds$path, p$counts),
    "--output-dir", containerPath(ds$path, p$results),
    containerPath(ds$path, workbookPath)
  ))
  database <- exorciseDockerArgs(ds, image, "exorcise-database", c(
    "--out-dir",     containerPath(ds$path, ds$path),
    "--results-dir", containerPath(ds$path, p$results),
    "--counts-dir",  containerPath(ds$path, p$counts),
    if (isTRUE(replace)) "--update-existing",
    # Non-interactive by construction: exorcise-database refuses a destructive
    # step rather than prompting when stdin is not a terminal, and here there is
    # nobody to prompt. The confirmation happened in the UI instead.
    "--force-overwrite",
    containerPath(ds$path, workbookPath)
  ))
  startJob(ds, list(analyse, database), label = "submit")
}

#' Rebuild one experiment from its retained workbook and results.
#'
#' No analysis is re-run: exorcise-database reads the result tables that are still
#' under <submit>/results/<id>/tables, so a corrected workbook republishes the
#' metadata and the statistics it already produced. --update-existing is what makes
#' it a replacement rather than a skip.
rebuildExperimentJob <- function(ds, config, workbookPath) {
  p <- submitPaths(ds$path)
  database <- exorciseDockerArgs(ds, config$exorcise_docker, "exorcise-database", c(
    "--out-dir",     containerPath(ds$path, ds$path),
    "--results-dir", containerPath(ds$path, p$results),
    "--counts-dir",  containerPath(ds$path, p$counts),
    "--update-existing", "--force-overwrite",
    containerPath(ds$path, workbookPath)
  ))
  startJob(ds, list(database), label = "rebuild")
}

#' Remove experiments from a dataset.
#'
#' The ids are not free text: the caller passes what the user ticked in a table
#' built from the dataset's own metadata, and removeExperimentsJob() refuses
#' anything that is not currently in the database.
removeExperimentsJob <- function(ds, config, experimentIds, present) {
  experimentIds <- unique(as.character(experimentIds))
  unknown <- setdiff(experimentIds, present)
  if (length(unknown)) {
    return(list(ok = FALSE, message = paste0(
      "Not in this dataset, so not removing: ", paste(unknown, collapse = ", "), ".")))
  }
  if (length(experimentIds) == 0) {
    return(list(ok = FALSE, message = "No experiments were selected."))
  }
  args <- exorciseDockerArgs(ds, config$exorcise_docker, "exorcise-database", c(
    "--out-dir", containerPath(ds$path, ds$path),
    "--remove", "--force-overwrite",
    experimentIds
  ))
  startJob(ds, list(args), label = "remove")
}

#### Metadata replacement ####

#' Replace a dataset's metadata tables with edited copies.
#'
#' The only write CRAVE makes to a dataset without going through
#' exorcise-database, and therefore the only one that could put the metadata and
#' the statistics out of step. It cannot, because the incoming tables have to name
#' exactly the experiments and comparisons already present: same Experiment IDs,
#' same Comparison IDs, same mapping between them. Labels may change -- which is
#' what makes this the way to correct a citation or a treatment name -- and
#' nothing else may.
#'
#' Adding or dropping rows is refused rather than fixed up, because both have a
#' correct route: submit the experiment, or remove it.
#'
#' @param experimentsFile,comparisonsFile Uploaded replacements, or NULL to keep
#'   the table that is there.
#' @return list(ok, message, changed)
replaceMetadataTables <- function(ds, experimentsFile = NULL, comparisonsFile = NULL) {
  if (is.null(experimentsFile) && is.null(comparisonsFile)) {
    return(list(ok = FALSE, message = "No replacement table was uploaded."))
  }
  if (!datasetWritable(ds$path)) {
    return(list(ok = FALSE, message = paste0(
      "The dataset directory '", ds$path, "' is not writable.")))
  }
  held <- datasetLock(ds$path)
  if (!is.null(held) && !isTRUE(held$stale)) {
    return(list(ok = FALSE, message = paste0(
      "A job (", held$jobId, ") is running against this dataset. Wait for it to ",
      "finish before replacing the metadata.")))
  }

  expPath <- file.path(ds$path, "experiments_metadata.csv.gz")
  cmpPath <- file.path(ds$path, "comparisons_metadata.csv.gz")

  current <- list(
    experiments = readTableGlob(file.path(ds$path, "experiments_metadata.csv"),
                                colClasses = "character"),
    comparisons = readTableGlob(file.path(ds$path, "comparisons_metadata.csv"),
                                colClasses = "character")
  )
  if (is.null(current$experiments) || is.null(current$comparisons)) {
    return(list(ok = FALSE, message = "This dataset's metadata files could not be read."))
  }

  incoming <- current
  problems <- character(0)
  readUpload <- function(path, what) {
    out <- tryCatch(
      data.table::fread(file = path, colClasses = "character", encoding = "UTF-8"),
      error = function(e) NULL)
    if (is.null(out)) {
      problems <<- c(problems, paste0("The ", what, " upload could not be read as a table."))
    }
    out
  }
  if (!is.null(experimentsFile)) incoming$experiments <- readUpload(experimentsFile, "experiments")
  if (!is.null(comparisonsFile)) incoming$comparisons <- readUpload(comparisonsFile, "comparisons")
  if (length(problems)) return(list(ok = FALSE, message = paste(problems, collapse = " ")))

  # -- The checks that make this safe ----------------------------------------
  sameSet <- function(a, b) setequal(as.character(a), as.character(b))

  for (spec in list(
    list(tbl = "experiments", id = "Experiment ID", needed = CRAVE_REQUIRED_EXPERIMENT_COLS),
    list(tbl = "comparisons", id = "Comparison ID", needed = CRAVE_REQUIRED_COMPARISON_COLS)
  )) {
    new <- incoming[[spec$tbl]]
    old <- current[[spec$tbl]]
    label <- paste0(spec$tbl, "_metadata")
    if (!spec$id %in% names(new)) {
      problems <- c(problems, paste0(label, " has no '", spec$id, "' column."))
      next
    }
    missingCols <- setdiff(spec$needed, names(new))
    if (length(missingCols)) {
      problems <- c(problems, paste0(
        label, " is missing the column(s) ", paste(missingCols, collapse = ", "), "."))
    }
    if (anyDuplicated(new[[spec$id]])) {
      problems <- c(problems, paste0(label, " repeats a ", spec$id, "."))
    }
    if (!sameSet(new[[spec$id]], old[[spec$id]])) {
      added   <- setdiff(as.character(new[[spec$id]]), as.character(old[[spec$id]]))
      dropped <- setdiff(as.character(old[[spec$id]]), as.character(new[[spec$id]]))
      problems <- c(problems, paste0(
        label, " does not list the same ", spec$id, "s as the database. ",
        if (length(added))   paste0("Added: ", summariseItems(added), ". ") else "",
        if (length(dropped)) paste0("Missing: ", summariseItems(dropped), ". ") else "",
        "Replacing the metadata can relabel rows but cannot add or remove them: ",
        "submit an experiment to add one, and use Remove to drop one."))
    }
  }

  # The comparison-to-experiment mapping is what joins the two tables, and the
  # statistics table keys on it, so it has to survive the edit as well.
  if (all(c("Experiment ID", "Comparison ID") %in% names(incoming$comparisons))) {
    key <- function(d) paste0(d$`Comparison ID`, "\u001f", d$`Experiment ID`)
    if (!sameSet(key(incoming$comparisons), key(current$comparisons))) {
      problems <- c(problems, paste0(
        "comparisons_metadata moves a comparison to a different experiment. That ",
        "mapping is what the statistics are keyed on, so it cannot be edited here."))
    }
  }
  if (length(problems)) {
    return(list(ok = FALSE, message = paste(problems, collapse = " ")))
  }

  # -- Write, keeping a copy of what was there -------------------------------
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  backup <- file.path(submitPaths(ds$path)$root, "metadata-backups", stamp)
  dir.create(backup, recursive = TRUE, showWarnings = FALSE)
  ok <- tryCatch({
    for (f in c(expPath, cmpPath)) if (file.exists(f)) file.copy(f, backup, overwrite = TRUE)
    data.table::fwrite(incoming$experiments, expPath, compress = "gzip")
    data.table::fwrite(incoming$comparisons, cmpPath, compress = "gzip")
    TRUE
  }, error = function(e) {
    log_error(paste0("Metadata replacement failed: ", conditionMessage(e)))
    FALSE
  })
  if (!ok) {
    return(list(ok = FALSE, message = paste0(
      "The metadata could not be written. The previous tables are in ", backup, ".")))
  }

  log_info(paste0("Replaced metadata tables on dataset '", ds$name,
                  "'; previous copies in ", backup, "."))
  list(ok = TRUE, changed = c(
    if (!is.null(experimentsFile)) "experiments_metadata.csv.gz",
    if (!is.null(comparisonsFile)) "comparisons_metadata.csv.gz"),
    message = paste0(
      "Metadata replaced. The previous tables were kept in ",
      file.path(CRAVE_SUBMIT_DIR, "metadata-backups", stamp), "."))
}

#### Integrity ####

#' Compare a dataset's metadata against the statistics actually in its database.
#'
#' The check the two tabs lean on: every comparison the metadata describes should
#' have statistics, and every comparison with statistics should be described. A
#' dataset only built and modified through this app cannot drift, so anything this
#' reports came from outside CRAVE -- a hand-edited file, an interrupted build, or
#' a database written by a different version of Exorcise.
#'
#' @return list(ok, nMetadata, nStats, missingStats, undescribed, message)
datasetIntegrity <- function(data, datasetName) {
  con <- data$cons[[datasetName]]
  if (is.null(con)) {
    return(list(ok = FALSE, message = "That dataset is not loaded."))
  }
  cmp <- data$comparisons[as.character(data$comparisons$Source) == datasetName, , drop = FALSE]
  described <- unique(as.character(cmp$`Comparison ID`))
  withStats <- tryCatch(
    DBI::dbGetQuery(con, "SELECT DISTINCT comparison_id FROM stat")$comparison_id,
    error = function(e) {
      log_error(paste0("Integrity query failed: ", conditionMessage(e)))
      NULL
    })
  if (is.null(withStats)) {
    return(list(ok = FALSE, message = "The statistics table could not be read."))
  }
  withStats <- unique(as.character(withStats))

  missingStats <- setdiff(described, withStats)
  undescribed  <- setdiff(withStats, described)
  ok <- length(missingStats) == 0 && length(undescribed) == 0

  list(
    ok           = ok,
    nMetadata    = length(described),
    nStats       = length(withStats),
    missingStats = missingStats,
    undescribed  = undescribed,
    message = if (ok) {
      paste0("In step: ", length(described), " comparison(s), all of them with statistics.")
    } else {
      paste0(
        if (length(missingStats)) paste0(
          length(missingStats), " comparison(s) are described in the metadata but ",
          "have no statistics: ", summariseItems(missingStats), ". ") else "",
        if (length(undescribed)) paste0(
          length(undescribed), " comparison(s) have statistics but are not in the ",
          "metadata, so nothing in CRAVE can reach them: ",
          summariseItems(undescribed), ". ") else "",
        "Rebuilding or removing the experiments involved will put them back in step.")
    }
  )
}

#' The experiments in one dataset, with what the Manage tab needs to show.
datasetExperiments <- function(data, datasetName) {
  exps <- data$experiments[as.character(data$experiments$Source) == datasetName, , drop = FALSE]
  if (nrow(exps) == 0) return(tibble())
  cmp <- data$comparisons[as.character(data$comparisons$Source) == datasetName, , drop = FALSE]
  counts <- table(as.character(cmp$`Experiment ID`))
  tibble(
    `Experiment ID` = as.character(exps$`Experiment ID`),
    Citation        = stripHtml(as.character(exps$Citation)),
    Organism        = as.character(exps$Organism),
    Comparisons     = as.integer(counts[as.character(exps$`Experiment ID`)] %||% 0L),
    Libraries       = as.character(exps$Libraries)
  ) %>%
    mutate(Comparisons = if_else(is.na(Comparisons), 0L, Comparisons)) %>%
    arrange(`Experiment ID`)
}

#' Workbooks CRAVE has retained for a dataset, keyed by experiment id.
retainedWorkbooks <- function(datasetPath) {
  dir <- submitPaths(datasetPath)$workbooks
  if (!dir.exists(dir)) return(character(0))
  files <- list.files(dir, pattern = "\\.xlsx$", full.names = TRUE)
  setNames(files, sub("\\.xlsx$", "", basename(files)))
}
