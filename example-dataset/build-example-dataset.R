#### Build the example dataset that ships with CRAVE ####
#
# data/ is a small subset of DDRcs (https://sjlab.cruk.cam.ac.uk/ddrcs/), the
# lab's portal of published DNA damage response screens: real screens rather
# than invented numbers, so the plots a new user meets show real biology, and
# small enough to ship inside the image.
#
# What it keeps, and why:
#
#   * Human screens of five DNA-damaging agents -- olaparib, talazoparib,
#     cisplatin, etoposide and camptothecin -- that compare a treated arm with
#     an untreated one in the same genotype. This is the slice Correlate's
#     Example button selects, and it has to be big enough for the tabs that
#     correlate screens against one another: Adjacency needs two screens,
#     Reduce four. Wild-type screens are taken first, because the example also
#     filters on genotype, then others up to MAX_PER_TREATMENT.
#   * a few essentialomes from those same experiments, so that Kind has more
#     than one value, Explore has something to show on its own terms, and
#     Chronos scores are not missing everywhere.
#   * only genome-wide, symbol-keyed screens. A base-editing or variant-tiling
#     library keys its statistics by variant rather than by gene, and a focused
#     panel of a few hundred genes leaves most of a gene list blank; neither
#     makes a useful demonstration. MIN_EXAMPLE_GENES is the test.
#   * a gene subset rather than whole libraries: the genes the Correlate
#     example asks for, the genes named in the ontology, and a random sample of
#     those measured across most of the screens, up to N_GENES. Whole libraries
#     over thirty screens would be a database of a few hundred MB, which is not
#     something to put in a Docker image.
#
# Usage, from this directory:
#
#   Rscript build-example-dataset.R /path/to/source-dataset [outdir]
#
# The source is not shipped -- DDRcs is 25 GB. Point the script at any CRAVE
# dataset: the selection is by metadata, so it keeps whatever of it matches.

suppressPackageStartupMessages({
  library(data.table); library(DBI); library(RSQLite)
})

N_GENES            <- 1500L
SEED               <- 101079L
# TREATMENTS comes from the app's constants, below: the dataset is built to
# hold the screens the Correlate example asks for.
MAX_PER_TREATMENT  <- 5L
MAX_ESSENTIALOMES  <- 6L
# A screen has to measure this many of the Correlate example's 200 genes to
# count as genome-wide and symbol-keyed. Variant-tiling and focused libraries
# fall far below it.
MIN_EXAMPLE_GENES  <- 150L
# A gene has to be measured in this fraction of the screens to be sampled as
# filler: one seen in two screens out of thirty correlates with nothing and only
# pads the download.
MIN_COVERAGE       <- 0.8

args   <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Give the path to the CRAVE dataset to subset.")
srcDir <- args[1]
outDir <- if (length(args) >= 2) args[2] else "data"

# The example gene list lives with the app, so the dataset is built to cover
# exactly what the Example button asks for.
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "shiny-server", "R", "01_constants.R"), local = TRUE)
TREATMENTS <- CRAVE_EXAMPLE_TREATMENTS

cat("Reading ", srcDir, "\n", sep = "")
exps <- fread(file.path(srcDir, "experiments_metadata.csv.gz"), colClasses = "character")
cmps <- fread(file.path(srcDir, "comparisons_metadata.csv.gz"), colClasses = "character")
# The ontology is optional in a CRAVE dataset, but the Enrichment tab needs one
# and the example is built to exercise it.
ontFile <- file.path(srcDir, "ontology.tsv.gz")
ont <- if (file.exists(ontFile)) fread(ontFile, colClasses = "character") else {
  warning("No ontology.tsv.gz in the source; the Enrichment tab will have nothing to test.")
  data.table(symbol = character(0), class = character(0))
}

organism <- setNames(exps$Organism, exps$`Experiment ID`)
isEssentialome <- function(x) grepl("[Ee]ssentialome", x)

con <- dbConnect(SQLite(), file.path(srcDir, "database.db"), flags = SQLITE_RO)
invisible(dbExecute(con, "CREATE TEMP TABLE example_gene (id TEXT PRIMARY KEY)"))
invisible(dbExecute(con, "INSERT INTO example_gene VALUES (:id)",
                    params = list(id = CRAVE_EXAMPLE_GENES)))
exampleGeneCount <- function(cid) dbGetQuery(con,
  "SELECT COUNT(DISTINCT gene_id) AS n FROM stat
   WHERE comparison_id = ? AND gene_id IN (SELECT id FROM example_gene)",
  params = list(cid))$n

# Same genotype on both arms: a treated-versus-untreated contrast, whatever the
# background, rather than one genotype against another.
human <- cmps[organism[`Experiment ID`] == "Human" & KO == ControlKO]
treated <- human[!isEssentialome(Contrast) & Treatment != ControlTreatment &
                 Treatment %in% TREATMENTS]
cat("Checking ", nrow(treated), " candidate screens for gene-level coverage\n", sep = "")
treated[, covered := vapply(`Comparison ID`, exampleGeneCount, integer(1))]
treated <- treated[covered >= MIN_EXAMPLE_GENES]
# Wild-type first, so the example's genotype filter keeps as many as the source
# allows; then a stable order, so the same screens come out every time.
treated[, wtFirst := KO != "Wildtype"]
setorder(treated, Treatment, wtFirst, `Experiment ID`, `Cell line`, `Comparison ID`)
treated <- treated[, head(.SD, MAX_PER_TREATMENT), by = Treatment]

essentialomes <- human[isEssentialome(Contrast) &
                       `Experiment ID` %in% treated$`Experiment ID`]
essentialomes[, covered := vapply(`Comparison ID`, exampleGeneCount, integer(1))]
essentialomes <- essentialomes[covered >= MIN_EXAMPLE_GENES]
setorder(essentialomes, `Experiment ID`, Contrast)
essentialomes <- head(essentialomes, MAX_ESSENTIALOMES)

keep <- rbind(treated, essentialomes, fill = TRUE)
keep[, c("covered", "wtFirst") := NULL]
# A blank "Days grown" anywhere in an arm makes CRAVE label the whole arm
# "Midpoint", which would hide these screens behind the example's timepoint
# filter. See the Timepoint block in R/03_data.R.
keep <- keep[!is.na(suppressWarnings(as.numeric(`Days grown`)))]
cat("Keeping ", nrow(keep), " comparisons from ", uniqueN(keep$`Experiment ID`),
    " experiments, ", uniqueN(keep$`Cell line`), " cell lines\n", sep = "")
print(keep[, .(screens = .N, wildtype = sum(KO == "Wildtype")), by = Treatment])

#### The gene subset ####
ids <- keep$`Comparison ID`
invisible(dbExecute(con, "CREATE TEMP TABLE want_comparison (id TEXT PRIMARY KEY)"))
invisible(dbExecute(con, "INSERT INTO want_comparison VALUES (:id)", params = list(id = ids)))

measured <- setDT(dbGetQuery(con, "
  SELECT gene_id, COUNT(DISTINCT comparison_id) AS n FROM stat
  WHERE comparison_id IN (SELECT id FROM want_comparison) GROUP BY gene_id"))

named  <- intersect(unique(c(CRAVE_EXAMPLE_GENES, ont$symbol)), measured$gene_id)
common <- measured[n >= MIN_COVERAGE * nrow(keep), gene_id]
set.seed(SEED)
filler <- setdiff(common, named)
filler <- sample(filler, max(0L, min(length(filler), N_GENES - length(named))))
genes  <- sort(unique(c(named, filler)))
cat("Keeping ", length(genes), " genes: ", length(named), " named by the example or the ontology, ",
    length(filler), " sampled from the ", length(common),
    " measured in at least ", round(100 * MIN_COVERAGE), "% of the screens\n", sep = "")

invisible(dbExecute(con, "CREATE TEMP TABLE want_gene (id TEXT PRIMARY KEY)"))
invisible(dbExecute(con, "INSERT INTO want_gene VALUES (:id)", params = list(id = genes)))

#### Write the dataset ####
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
dbFile <- file.path(outDir, "database.db")
if (file.exists(dbFile)) file.remove(dbFile)
out <- dbConnect(SQLite(), dbFile)
# The same four tables exorcise-database writes. CRAVE reads its metadata from
# the two CSVs beside this file, so experiment and comparison stay empty here,
# exactly as they are in a dataset Exorcise built.
invisible(dbExecute(out, "CREATE TABLE gene (id VARCHAR NOT NULL PRIMARY KEY, symbol VARCHAR NOT NULL,
  official_id VARCHAR, symbol_with_ids VARCHAR, organism VARCHAR NOT NULL)"))
invisible(dbExecute(out, "CREATE TABLE experiment (stringid VARCHAR NOT NULL PRIMARY KEY, date DATETIME,
  library VARCHAR NOT NULL, doi VARCHAR, representation VARCHAR, moi VARCHAR, description VARCHAR,
  notes VARCHAR, reference VARCHAR, source VARCHAR, citation VARCHAR)"))
invisible(dbExecute(out, "CREATE TABLE comparison (stringid VARCHAR NOT NULL PRIMARY KEY, contrast VARCHAR NOT NULL,
  treatment_label VARCHAR NOT NULL, timepoint VARCHAR NOT NULL, cell VARCHAR NOT NULL,
  control_sample VARCHAR NOT NULL, test_sample VARCHAR NOT NULL, ko VARCHAR NOT NULL,
  control_treatment VARCHAR, control_ko VARCHAR, dose VARCHAR, gi VARCHAR, days_grown INTEGER,
  library VARCHAR NOT NULL, notes VARCHAR)"))
invisible(dbExecute(out, "CREATE TABLE stat (comparison_id VARCHAR NOT NULL, experiment_id VARCHAR NOT NULL,
  gene_id VARCHAR NOT NULL, analysis_type_id INTEGER NOT NULL, score FLOAT, fdr FLOAT, fdr10 FLOAT,
  pos_p FLOAT, neg_p FLOAT, PRIMARY KEY (comparison_id, gene_id, analysis_type_id))"))

dbWriteTable(out, "gene", dbGetQuery(con,
  "SELECT * FROM gene WHERE id IN (SELECT id FROM want_gene)"), append = TRUE)

total <- 0L
for (cid in ids) {
  st <- dbGetQuery(con, "SELECT * FROM stat WHERE comparison_id = ?
                         AND gene_id IN (SELECT id FROM want_gene)", params = list(cid))
  if (nrow(st)) { dbWriteTable(out, "stat", st, append = TRUE); total <- total + nrow(st) }
}
dbDisconnect(con)
invisible(dbExecute(out, "VACUUM"))
cat("Wrote ", total, " statistics rows, ",
    dbGetQuery(out, "SELECT COUNT(DISTINCT analysis_type_id) FROM stat")[[1]],
    " analysis method(s)\n", sep = "")
dbDisconnect(out)

fwrite(exps[`Experiment ID` %in% keep$`Experiment ID`],
       file.path(outDir, "experiments_metadata.csv.gz"), compress = "gzip")
fwrite(keep, file.path(outDir, "comparisons_metadata.csv.gz"), compress = "gzip")
if (nrow(ont)) fwrite(ont, file.path(outDir, "ontology.tsv.gz"), sep = "\t", compress = "gzip")

cat("\n", normalizePath(outDir), "\n", sep = "")
for (f in sort(list.files(outDir))) cat(sprintf("   %-34s %8.1f KB\n", f,
  file.info(file.path(outDir, f))$size / 1024))
cat(sprintf("   %-34s %8.1f MB\n", "total",
  sum(file.info(list.files(outDir, full.names = TRUE))$size) / 1024^2))
