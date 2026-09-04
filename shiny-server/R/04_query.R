#### Statistics queries ####
# fetchStat() is the single entry point for reading the `stat` table. It routes each
# comparison to the dataset that holds it, chunks bound parameters so that no query
# exceeds SQLite's limit, and rounds once at the end.

# Conservative chunk size: stays under SQLITE_MAX_VARIABLE_NUMBER on old builds.
CRAVE_SQL_PARAM_CHUNK <- 400L

#' Split a vector into chunks of at most `size`.
chunkVector <- function(x, size = CRAVE_SQL_PARAM_CHUNK) {
  n <- length(x)
  if (n == 0) return(list())
  if (n <= size) return(list(x))
  split(x, ceiling(seq_len(n) / size))
}

#' Translate method names to analysis_type_id codes, dropping unknown names.
methodCodes <- function(method) {
  if (is.null(method) || length(method) == 0) return(integer(0))
  m <- unique(method)
  m <- m[m %in% names(CRAVE_METHODS)]
  unname(CRAVE_METHODS[m])
}

#' Build one parameterised SELECT against the `stat` table.
#'
#' @return list(sql, params)
buildStatQuery <- function(comparison = NULL, gene = NULL, codes = integer(0)) {
  clauses <- character(0)
  params  <- list()

  placeholders <- function(n) paste(rep("?", n), collapse = ",")

  if (length(comparison)) {
    clauses <- c(clauses, paste0("comparison_id IN (", placeholders(length(comparison)), ")"))
    params  <- c(params, as.list(as.character(comparison)))
  }
  if (length(codes)) {
    clauses <- c(clauses, paste0("analysis_type_id IN (", placeholders(length(codes)), ")"))
    params  <- c(params, as.list(as.integer(codes)))
  }
  if (length(gene)) {
    clauses <- c(clauses, paste0("gene_id IN (", placeholders(length(gene)), ")"))
    params  <- c(params, as.list(as.character(gene)))
  }

  sql <- "SELECT * FROM stat"
  if (length(clauses)) sql <- paste(sql, "WHERE", paste(clauses, collapse = " AND "))
  list(sql = sql, params = params)
}

#' Fetch statistics, chunking gene and comparison lists so no single query
#' exceeds SQLite's bound-parameter limit.
fetchStatOne <- function(con, comparison, gene, codes) {
  cmpChunks  <- if (length(comparison)) chunkVector(comparison) else list(NULL)
  geneChunks <- if (length(gene))       chunkVector(gene)       else list(NULL)
  out <- vector("list", length(cmpChunks) * length(geneChunks))
  i <- 0L
  for (cc in cmpChunks) {
    for (gc in geneChunks) {
      i <- i + 1L
      q <- buildStatQuery(cc, gc, codes)
      out[[i]] <- tryCatch(
        DBI::dbGetQuery(con, q$sql, params = q$params),
        error = function(e) {
          log_error(paste0("Statistics query failed: ", conditionMessage(e)))
          NULL
        }
      )
    }
  }
  data.table::rbindlist(Filter(Negate(is.null), out), use.names = TRUE, fill = TRUE)
}

#' Fetch statistics for the given comparisons, genes and methods.
#'
#' @param data The app data object from loadCraveData().
#' @param comparison Character vector of comparison IDs, or NULL for all.
#' @param gene Character vector of gene symbols, or NULL for all.
#' @param method Character vector of method names, or NULL for all methods.
#' @param normalise TRUE to quantile normalise scores across the screens returned by
#'   this query, so that every screen shares one score distribution. The reference
#'   distribution is therefore built from the comparisons asked for here, and changes
#'   if the caller's screen selection changes.
#' @return A tibble with the columns of the `stat` table.
fetchStat <- function(data, comparison = NULL, gene = NULL, method = NULL,
                      normalise = FALSE) {
  codes <- methodCodes(method)
  if (length(method) && length(codes) == 0) {
    log_warn(paste0("No recognised analysis method in: ",
                    paste(method, collapse = ", "), ". Querying all methods."))
  }

  log_info(paste0(
    "Statistics query: ", length(comparison %||% character(0)), " comparison(s), ",
    length(gene %||% character(0)), " gene(s), method(s): ",
    paste(method %||% "all", collapse = ",") , "."
  ))

  results <- list()

  emptyResult <- function() {
    as_tibble(if (is.null(data$statTemplate)) data.frame() else data$statTemplate)
  }

  # NULL means "no comparison restriction". An empty vector means "the caller asked
  # for a specific set of comparisons and none of them resolved", which must return
  # nothing rather than everything.
  if (!is.null(comparison) && length(comparison) == 0) {
    log_info("Statistics query asked for zero comparisons; returning no rows.")
    return(emptyResult())
  }

  if (is.null(comparison)) {
    # No comparison restriction: every dataset is in scope.
    for (nm in names(data$cons)) {
      results[[nm]] <- fetchStatOne(data$cons[[nm]], NULL, gene, codes)
    }
  } else {
    # Route each comparison to the dataset that holds it, in one indexed lookup
    # rather than one dplyr filter per dataset.
    comparison <- as.character(comparison)
    owner <- data$comparisonSource[comparison]
    known <- !is.na(owner)
    if (any(!known)) {
      log_warn(paste0("Ignoring ", sum(!known), " comparison ID(s) not present in any dataset."))
    }
    for (nm in unique(owner[known])) {
      ids <- comparison[known][owner[known] == nm]
      if (length(ids) == 0) next
      results[[nm]] <- fetchStatOne(data$cons[[nm]], ids, gene, codes)
    }
  }

  fetched <- data.table::rbindlist(
    Filter(function(x) !is.null(x) && nrow(x) > 0, results),
    use.names = TRUE, fill = TRUE
  )

  if (nrow(fetched) == 0) {
    log_info("Statistics query returned no rows.")
    return(emptyResult())
  }

  # Normalise before rounding, so the rounding applies to the values that are
  # actually reported. Normalisation is rank-preserving within a screen, so FDRs and
  # rank-based cutoffs are untouched; only score magnitudes change.
  if (isTRUE(normalise)) {
    fetched <- data.table::as.data.table(
      quantileNormaliseScores(fetched, valueCol = "score", byCol = "comparison_id")
    )
  }

  # Round once, at the end.
  fetched[, `:=`(fdr = signif(fdr, 4), score = round(score, 4))]
  log_info(paste0("Statistics query returned ", nrow(fetched), " rows."))
  as_tibble(fetched)
}

#' Which genes appear in the given comparisons?
#'
#' Used to narrow the Explore gene picker to genes actually present in the selected
#' screens. Only the datasets holding those comparisons are queried.
fetchGenesForComparisons <- function(data, comparison) {
  if (length(comparison) == 0) return(character(0))
  comparison <- as.character(comparison)
  owner <- data$comparisonSource[comparison]
  known <- !is.na(owner)
  out <- character(0)
  for (nm in unique(owner[known])) {
    ids <- comparison[known][owner[known] == nm]
    for (chunk in chunkVector(ids)) {
      sql <- paste0("SELECT DISTINCT gene_id FROM stat WHERE comparison_id IN (",
                    paste(rep("?", length(chunk)), collapse = ","), ")")
      got <- tryCatch(
        DBI::dbGetQuery(data$cons[[nm]], sql, params = as.list(chunk))$gene_id,
        error = function(e) {
          log_error(paste0("Gene lookup failed: ", conditionMessage(e)))
          character(0)
        }
      )
      out <- c(out, got)
    }
  }
  unique(out)
}

#' Map user-facing contrast labels (FriendlyID) to comparison IDs.
comparisonIdsFor <- function(data, friendlyIds) {
  if (length(friendlyIds) == 0) return(character(0))
  data$comparisons$`Comparison ID`[data$comparisons$FriendlyID %in% friendlyIds]
}

#' Map comparison IDs back to contrast labels.
friendlyIdsFor <- function(data, comparisonIds) {
  if (length(comparisonIds) == 0) return(character(0))
  data$comparisons$FriendlyID[data$comparisons$`Comparison ID` %in% comparisonIds]
}
