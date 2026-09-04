#### Small shared helpers ####

#' Negative log scale for ggplot.
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv   <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}

#' Read a table whose exact extension is not known.
#'
#' CRAVE datasets ship metadata as .csv.gz but tolerate uncompressed files, so the
#' path is treated as a prefix. Deliberately not called `fread`: that name would
#' shadow data.table::fread and route every ordinary call through a Sys.glob().
readTableGlob <- function(path, ...) {
  hits <- Sys.glob(paste0(path, "*"))
  if (length(hits) == 0) return(NULL)
  data.table::fread(file = hits[1], ...)
}

#' Strip HTML tags from a character vector.
stripHtml <- function(x) gsub("<.+?>", "", x)

#' Split a free-text list into tokens on commas, spaces, tabs and newlines.
splitTokens <- function(x) {
  if (length(x) == 0) return(character(0))
  out <- unlist(strsplit(as.character(x), "[ \t\n\r,;]+"))
  out[nzchar(out)]
}

#' Return `x` unless it is NULL or empty, in which case return `default`.
`%||%` <- function(x, default) {
  if (is.null(x) || length(x) == 0) default else x
}

#' Numeric coercion that falls back to a default instead of producing NA.
asNumericOr <- function(x, default) {
  if (is.null(x) || length(x) == 0) return(default)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0 || any(is.na(v))) default else v
}

#' TRUE if `x` is a string that parses as a number.
isNumericString <- function(x) {
  length(x) == 1 && !is.na(x) && grepl("^\\s*-?\\d+\\.?\\d*([eE][-+]?\\d+)?\\s*$", x)
}

#' Human-readable summary of a character vector, truncated to a budget.
#'
#' Used for the "screens found" and "selected experiments" readouts. Works on the
#' vector directly rather than on a joined string, so it stays fast on large
#' selections and the leftover count is exact.
summariseItems <- function(x, max_chars = CRAVE_SUMMARY_CHARS) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  n <- length(x)
  if (n == 0) return("")
  if (n == 1) return(x)
  if (nchar(paste(x, collapse = "; ")) <= max_chars) {
    return(paste0(paste(x[-n], collapse = "; "), "; and ", x[n]))
  }
  widths <- cumsum(nchar(x) + 2L)
  keep   <- max(1L, sum(widths <= max_chars))
  paste0(paste(x[seq_len(keep)], collapse = "; "), "... and ", n - keep, " more.")
}

#' Format a contrast label for a plot title.
#'
#' The label is "<citation>: <contrast>". The citation is split off once, the
#' contrast body is normalised (recursively, because "Day 21" becomes
#' "Essentialome from day 21" and a leading parenthetical is stripped before the
#' rest is read), and the citation is appended once at the end.
contrastFormatter <- function(comparison) {
  if (length(comparison) == 0 || is.na(comparison[1]) || !nzchar(comparison[1])) return("")
  comparison <- comparison[1]
  stem <- sub("^(.+?): .+$", "\\1", comparison)
  body <- sub("^.+?: ", "", comparison)

  formatBody <- function(x, depth = 0L) {
    if (depth > 4L || !nzchar(x)) return(x)
    if (grepl("^Day", x)) {
      return(formatBody(sub("^D", "Essentialome from d", x), depth + 1L))
    }
    if (grepl("^\\(", x)) {
      stripped <- sub("^\\(.+?\\)\\s*", "", x)
      # Guard against a malformed label that the substitution cannot shorten.
      if (identical(stripped, x)) return(x)
      return(formatBody(stripped, depth + 1L))
    }
    if (grepl("^Essentialome", x)) {
      return(sub("^(Essentialome.*?)((, | at | in | under )(.+))?$",
                 "<b>\\1</b><i>\\2</i>", x))
    }
    if (grepl("^Sorted", x)) {
      return(sub("^(Sorted.*?)((, | at | in | under )(.+))?$",
                 "<b>\\1</b><i>\\2</i>", x))
    }
    sub("^(.+?)((, | at | in | under )(.+))?$", "<b>Effect of \\1</b><i>\\2</i>", x)
  }

  paste0(formatBody(body), ": ", stem)
}

#' Capitalise the first character.
toupper_first_initial <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#' Evaluate an expression with a fixed random seed, then restore the RNG.
#'
#' Some algorithms need a seed to be reproducible, but a bare set.seed() makes the
#' whole process's random stream deterministic from that point on — which is how two
#' sessions can end up generating identical Exorcise task ids. This confines the seed
#' to one expression.
withSeed <- function(seed, expr) {
  hadSeed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  oldSeed <- if (hadSeed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (hadSeed) assign(".Random.seed", oldSeed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

#### External links in prose ####
# htmltools separates a tag's siblings with a newline, which HTML collapses to the
# single space a sentence needs. Applying .noWS to a link removes that space and
# welds the link to the words on either side, so these helpers exist to make the
# right choice obvious at the call site.

#' A link in the middle of a sentence, keeping one space on each side.
extLink <- function(href, label) {
  tags$a(href = href, label, target = "_blank")
}

#' A link immediately followed by punctuation, suppressing only the trailing space.
extLinkTight <- function(href, label) {
  tags$a(href = href, label, target = "_blank", .noWS = "after")
}

#' Prose naming the analysis methods available, as raw HTML.
methodLinks <- function(chronos = TRUE) {
  paste0(
    "Choose between <a href=\"", DOI_DRUGZ, "\" target=\"_blank\">DrugZ</a>, ",
    "<a href=\"", DOI_MAGECK, "\" target=\"_blank\">MAGeCK</a>",
    if (chronos) paste0(", and <a href=\"", DOI_CHRONOS, "\" target=\"_blank\">Chronos</a>") else "",
    " methods."
  )
}

#### Placeholder plots ####
# Every analysis needs to render "nothing to show, here's why". These were
# open-coded twenty-odd times with slightly different styling.

#' A ggplot carrying a single centred message.
msgGg <- function(text) {
  d <- tibble(x = 0, y = 0, label = text)
  ggplot(d, aes(x = x, y = y, label = label)) +
    geom_text(size = 5) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank()) +
    xlab("") + ylab("")
}

#' A plotly carrying a single centred message, plus the empty table to match.
#'
#' @return list(p, plotdata, height) so callers can splice it straight into the
#'   same shape their success path returns.
msgResult <- function(text, height = 240) {
  list(
    p        = ggplotly(msgGg(text)),
    plotdata = tibble(x = 0, y = 0, label = text),
    height   = height
  )
}

#' A visNetwork carrying a single centred message.
msgNetwork <- function(text) {
  visNetwork(nodes = tibble(id = 1, label = text))
}

#### Plot sizing ####

#' Resolve height/width from a customisation option list.
#'
#' @param opts Named list of customisation values (see mod_customise.R).
#' @param defaults Named list of fallbacks.
plotDims <- function(opts, defaults) {
  height <- asNumericOr(opts$height, defaults$height)
  if (height == 0) height <- defaults$height
  width <- if (isTRUE(opts$width_auto)) NULL else asNumericOr(opts$width, defaults$width)
  if (!is.null(width) && width == 0) width <- NULL
  list(height = height, width = width)
}

#' Map a "Linear"/"Logarithmic"/"Automatic" choice to a plotly axis type.
#'
#' @param allow_log Whether a log axis is meaningful for this plot.
#' @param auto What "Automatic" means here.
axisType <- function(choice, allow_log = TRUE, auto = "linear") {
  choice <- choice %||% "Automatic"
  if (identical(choice, "Linear")) return("linear")
  if (identical(choice, "Logarithmic")) return(if (allow_log) "log" else "linear")
  auto
}

#### Save/Load state helpers ####
# Each module declares its own input ids; these two helpers read and write them.

#' Read a set of inputs into a plain list.
captureInputs <- function(input, ids) {
  out <- lapply(ids, function(i) input[[i]])
  names(out) <- ids
  out
}

#' Push a saved state back into a set of inputs.
#'
#' Widget type is inferred from the stored value, which covers checkboxes,
#' numeric inputs and selectize inputs. Text and slider inputs are ambiguous
#' (both can hold a bare number or string), so modules name those explicitly.
#'
#' @param textIds Ids backed by textInput or textAreaInput.
#' @param sliderIds Ids backed by sliderInput.
#' @param switchIds Ids backed by shinyWidgets::materialSwitch, which has its own
#'   updater and does not respond to updateCheckboxInput.
#' @param numericIds Ids backed by numericInput. Named explicitly because an empty
#'   numeric box arrives as logical NA, which would otherwise look like a checkbox.
#' @param serverChoices Named list of choice vectors for selectize inputs that were
#'   created with server = TRUE. Such an input only holds the options the client has
#'   paged in, so selectize discards a `selected` value it has never seen; the
#'   choices must be resent with it.
restoreInputs <- function(session, state, ids,
                          textIds = character(0), sliderIds = character(0),
                          switchIds = character(0), numericIds = character(0),
                          serverChoices = list()) {
  if (!is.list(state)) return(invisible(NULL))
  for (i in intersect(ids, names(state))) {
    v <- state[[i]]
    if (i %in% textIds) {
      updateTextAreaInput(session, i, value = v %||% "")
    } else if (i %in% sliderIds) {
      if (!is.null(v)) updateSliderInput(session, i, value = v)
    } else if (i %in% switchIds) {
      shinyWidgets::updateMaterialSwitch(session, i, value = isTRUE(v))
    } else if (i %in% numericIds) {
      updateNumericInput(session, i, value = if (length(v) == 1) v else NA)
    } else if (i %in% names(serverChoices)) {
      updateSelectizeInput(session, i, choices = serverChoices[[i]],
                           selected = v, server = TRUE)
    } else if (is.logical(v) && length(v) == 1 && !is.na(v)) {
      updateCheckboxInput(session, i, value = v)
    } else if (is.numeric(v) && length(v) == 1) {
      updateNumericInput(session, i, value = v)
    } else {
      updateSelectizeInput(session, i, selected = v)
    }
  }
  invisible(NULL)
}

