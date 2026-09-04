#### Package loading ####
# Only the packages needed to build the UI and load the datasets are attached at
# start-up. Everything else is attached on first use via needPkg(), which costs
# nothing on subsequent calls because library() short-circuits when a package is
# already attached.
#
# There is a third category below, CRAVE_PKGS_INSTALLED_ONLY: packages that must be
# installed but that CRAVE never loads, because another package reaches for them on
# our behalf. Read the comment there before removing anything from it.

Sys.setenv(
  `_R_CHECK_LENGTH_1_CONDITION_` = FALSE,
  `_R_CHECK_LENGTH_1_LOGIC2_`    = FALSE
)

# Packages required to serve the first page and to read the datasets.
CRAVE_PKGS_EAGER <- c(
  "shiny", "shinythemes", "shinyWidgets", "shinybusy", "shinyjs",
  "DT", "plotly", "visNetwork", "colourpicker",
  "tibble", "dplyr", "tidyr", "data.table",
  "ggplot2", "scales", "DBI", "RSQLite", "logger"
)

# Packages deferred until the analysis that needs them is actually run.
# Kept as a manifest so the Dockerfile and this file cannot drift apart.
CRAVE_PKGS_LAZY <- c(
  "cluster",        # agnes(): clustergram + correlation heatmap dendrograms
  "ggdendro",       # ggdendrogram(): the same dendrograms
  "igraph",         # correlation network
  "umap",           # Reduce tab, UMAP method
  "Rtsne",          # Reduce tab, t-SNE method
  "mice",           # Reduce tab imputation, all methods
  "ggVennDiagram",  # Overlap tab, Venn style
  "ggupset",        # Overlap tab, Upset style
  "DescTools"       # AUC(): ROC tab
)

# Packages that must be INSTALLED but are never loaded or attached by CRAVE.
#
# Do not delete these because nothing appears to call them. Grepping the source
# for callers will find none: they are reached indirectly, from inside another
# package.
#
#   R.utils - data.table::fread() hands compressed input to
#             R.utils::decompressFile(). Every CRAVE metadata file is gzipped
#             (experiments_metadata.csv.gz, comparisons_metadata.csv.gz,
#             ontology.tsv.gz, libraries.txt.gz), so without R.utils on the library
#             path fread() fails on all of them and no dataset will load. Recent
#             data.table versions decompress .gz natively via zlib and only need
#             R.utils for .bz2, but which of those applies depends on the
#             data.table build, so it stays installed either way.
CRAVE_PKGS_INSTALLED_ONLY <- c("R.utils")

suppressPackageStartupMessages({
  for (p in CRAVE_PKGS_EAGER) library(p, character.only = TRUE)
})

# Presence check only: system.file() looks on the library path without loading the
# namespace, so this costs nothing and does not attach anything. Failing here with
# a clear message beats failing later with an opaque fread() error on every file.
for (p in CRAVE_PKGS_INSTALLED_ONLY) {
  if (!nzchar(system.file(package = p))) {
    log_warn(paste0(
      "Package '", p, "' is not installed. It is never called by CRAVE directly, ",
      "but data.table::fread() needs it to read the gzipped metadata files, so ",
      "dataset loading is likely to fail. Install it with install.packages('", p, "')."
    ))
  }
}

#' Attach a deferred package on first use.
#'
#' @param ... Package names.
#' @return TRUE if every package is available, FALSE otherwise.
needPkg <- function(...) {
  pkgs <- c(...)
  ok <- TRUE
  for (p in pkgs) {
    if (isTRUE(p %in% .packages())) next
    if (!requireNamespace(p, quietly = TRUE)) {
      log_error(paste0("Package '", p, "' is not installed; the analysis that needs it is unavailable."))
      ok <- FALSE
      next
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  ok
}

options(warn = 1)

# Lets app.R tell whether Shiny has already auto-sourced R/ (it does so for app
# directories from Shiny 1.5.0 onwards) so that it does not source everything a
# second time, into a second environment.
options(crave.sourced = TRUE)
