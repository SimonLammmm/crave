# Changelog

## 5.1.2 — 2026-08-18

- **Fixed: the Reduce `Others` column was empty for every gene.** `summarise()`
  evaluates its arguments in order and each new column shadows any existing one of
  the same name, so `summarise(Class = Class[1], Others = paste(Class[-1], ...))`
  read the new length-1 `Class` in the second expression and `Class[-1]` was always
  `character(0)`. The ontology's class column is now called `cls` on input, so
  neither output can shadow it and the two expressions are order-independent —
  writing them in the working order would have fixed the symptom but left the next
  edit free to reintroduce it.
- Audited the rest of the codebase for the same pattern. Twelve other places create
  a column and read it in a later argument of the same `mutate()`/`transmute()`, and
  all twelve are deliberate: the plasmid reference chain in `wrangleMetadata()`,
  `n = N - m` feeding `phyper` in `classEnrichment()`, the FDR being derived from
  unrounded p-values before those are rounded, and the hover-text builders in the
  Explore plots.
- **Added `tools/rcheck.py`**, the static checker used while developing 5.0-5.1,
  standing in for an R linter where R is not installed. Python 3, no packages, exits
  non-zero on a hard failure so it can be a pre-commit hook. Checks bracket balance,
  `plotPanelUI`/`plotPanelServer` pairing, agreement between the package manifests
  and the Dockerfile, unresolved function calls, module input id cross-referencing,
  and the dplyr shadowing pattern above. See the README for what it can and cannot
  see.

---

## 5.1.1 — 2026-08-18

- **Reduce now projects the whole library.** The 250-gene cap suits the analyses
  that put a gene on an axis (Gene query, Clustergram) or on both axes (Heatmap,
  Network), where a matrix of hits is the point. An embedding is different: cluster
  structure only appears against the full gene set. Two constants govern this:
  `CRAVE_MAX_AUTO_GOI_REDUCE` (default `Inf`) and `CRAVE_REDUCE_USES_CUTOFF`
  (default `FALSE`). The second is the one that matters — lifting the cap alone did
  nothing, because gene selection filtered on the sidebar cutoff *before* the cap
  applied, so the tab still showed hits only. Reduce now takes every gene measured
  in the selected screens; a deployer wanting the old behaviour sets
  `CRAVE_REDUCE_USES_CUTOFF <- TRUE`, and one wanting a ceiling sets a finite cap,
  which keeps the strongest genes by absolute score.
- **Reduce plots one point per gene when a gene has several ontology classes.**
  Previously such genes were drawn once per class, overplotted at identical
  coordinates and duplicated in the results table. The lexicographically first class
  is now the one plotted and coloured, with the remainder listed in a new `Others`
  column so nothing is hidden. The ontology is collapsed to one row per gene before
  the join, so the join is one-to-one by construction rather than relying on
  post-hoc grouping — which matters because the brush table is built from the same
  rows.
- Reduce skips imputation when nothing is missing, and asks `mice()` for one
  imputation rather than five. `complete()` only ever used the first, so the other
  four were computed and discarded — negligible at 100 genes, not at 20,000.
- Reduce logs the gene and screen count before projecting, so a slow run is
  attributable.

---

## 5.1.0 — 2026-08-18

### New features

- **Quantile normalisation.** A switch at the top of the Correlate screen filter
  controls gives every screen the same score distribution. The reference is built
  from the screens the sidebar filters currently select, so it moves with the
  filters; results from two different filter selections are not directly
  comparable. Because screens rarely measure identical gene sets, the reference is
  formed by averaging the screens' own quantile functions over a shared probability
  grid rather than by sorting equal-length columns. Normalisation is rank-preserving
  within a screen, so FDRs and rank cutoffs are unaffected — only score magnitudes,
  and therefore score cutoffs, change. Applies to every Correlate analysis and to
  the bulk download.
- **Regex mode for the Correlate text filter.** A switch toggles between literal
  substring matching and Perl-compatible regular expressions. An invalid pattern
  matches nothing and reports the parse error under the field, rather than raising.
- **Similarity measures for Heatmap and Network.** Spearman (still the default),
  Pearson, Kendall, and cosine similarity. All four run from -1 to 1, so the
  divergent colour scale and the network's cutoff behave identically across them.
  Cosine takes a vectorised path on complete matrices and falls back to a pairwise
  loop only when values are missing. Kendall logs a warning above 300 genes, being
  quadratic in the number of pairs.
- **Clustering methods for Heatmap and Clustergram.** Agnes with Ward (default),
  average, single, complete, flexible and weighted linkage, plus Clara with a
  user-supplied k. Distance measures are euclidean and manhattan, with jaccard
  offered only for Clara — the distance dropdown updates itself to match the chosen
  method. Clara partitions rather than building a tree, so each cluster takes a
  contiguous block of the axis and the dendrogram drawn is over the cluster medoids,
  stretched so each leaf sits above the centre of its block. Cluster assignments are
  added to the results table.
- **UMAP is now Reduce**, with UMAP (default), PCA and t-SNE. Axis labels are
  per-method, and PCA reports the variance explained by each component. Each method
  states its own minimum screen count (UMAP 4, t-SNE 3, PCA 2) and t-SNE's
  perplexity is derived from the gene count so it always satisfies Rtsne's
  constraint. Result columns are now `Dim1`/`Dim2` rather than `UMAP1`/`UMAP2`.
  Adds `Rtsne` as a deferred dependency.
- **Venn-compartment enrichment.** A switch on the Enrichment tab replaces the
  one-column-per-screen layout with one column per exclusive combination of screens:
  genes hit in A alone, in A and B but not C, and so on. The x axis is drawn as an
  upset matrix naming the screens in each compartment. Limited to
  `CRAVE_MAX_VENN_SETS` (7) screens; above that it explains the limit instead of
  drawing. Compartment membership is only defined over genes every selected screen
  measured, so that intersection is used as the hypergeometric population.
  Rendered with `renderPlot`, since ggupset's combination-matrix axis has no plotly
  equivalent.

### Fixes

- **Links in prose lost the spaces around them.** `.noWS` was applied wholesale to
  anchors sitting mid-sentence, and it includes `"outside"`, which strips the
  whitespace on both sides — so "uses the … method" rendered as one run-on word.
  Two helpers now make the choice explicit: `extLink()` for a link inside a
  sentence, `extLinkTight()` for one followed immediately by punctuation. `noWS`
  itself remains, for the case it was meant for.
- **Save/Load silently dropped every server-side selectize selection.** The gene
  lists, the sixteen Correlate screen filters, the contrast pickers and the guide
  library pickers are all created with `server = TRUE`, which means the client holds
  only the options it has paged in — so a bare `selected` was discarded by
  selectize. Restoring a key left them empty. `restoreInputs()` now resends the
  choices alongside the selection for those inputs.
- **Reduce imputed in the wrong orientation.** `mice()` was handed the transposed
  matrix, making screens the observations and genes the variables: with the
  hundred-gene cap and a two-to-four-screen minimum that is a hundred predictors
  fitted on a handful of rows, which is why the tab needed a `try()` and a "server
  error" fallback. Genes are now the observations and screens the variables, so a
  gene's missing score in one screen is predicted from its scores in the others,
  fitted across all genes. **This changes UMAP output** — worth comparing against a
  known result.
- An empty numeric field arrives as logical `NA`, so Save/Load restored `cutoff`,
  `k` and the Pendragonator list length through `updateCheckboxInput`. Numeric ids
  are now declared explicitly.
- `CRAVE_REDUCE_MIN_SCREENS[[method]]` raised `subscript out of bounds` rather than
  falling back, for a method name outside the known set.
- A `NULL` distance selection made the metric observer's `if` condition
  length-zero.
- t-SNE needs a seed to be reproducible; it uses a new `withSeed()` helper that
  snapshots and restores `.Random.seed`, rather than seeding the process globally.

---

## 5.0.0 — 2026-08-17

Reorganisation and performance pass. No new features; the aim was faster start-up,
faster analyses, and a layout that is possible to navigate. The version history
that used to sit as a comment block at the top of `app.R` is preserved below.

### Layout

The application was 24 files in one directory: a 1584-line `server()`, a 765-line
block of UI element definitions, and one file per analysis that reached into
global variables. It is now:

```
shiny-server/
  app.R                    entry point only (82 lines)
  config.example.R
  R/
    00_packages.R          eager vs. deferred package loading
    01_constants.R         version, limits, defaults, static reference data
    02_utils.R             shared helpers
    03_data.R              dataset loading, metadata cache, SQLite indexing
    04_query.R             parameterised statistics queries
    05_analysis_common.R   primitives shared between analyses
    06_analysis_explore.R  volcano, rank, biplot, overlap, ROC
    07_analysis_correlate.R gene query, clustergram, heatmap, network, UMAP,
                           enrichment, Pendragonator, bulk download, guides
    08_analysis_exorcise.R Exorcise
    09_config.R            config.R loading and validation
    10_mod_plotpanel.R     reusable submit/plot/table/download panel
    11_mod_customise.R     reusable "Customise plot" modal
    12_mod_explore.R       Explore tab
    13_mod_correlate.R     Correlate tab
    14_mod_tables.R        Libraries, Ontology, Legal
    15_mod_guides.R        Guides tab
    16_mod_exorcise.R      Exorcise tab
    17_mod_saveload.R      Save/Load
    18_mod_home.R          Home
    19_ui.R                UI assembly
    20_server.R            module wiring
```

Each tab is a Shiny module with namespaced input ids. Analysis functions take
explicit arguments instead of Shiny's `input` object, so they are callable and
testable outside a session. Every `<<-` global is gone: shared data is passed in,
and cross-tab communication goes through one explicit message bus.

### Start-up time

- **Metadata caching.** The per-dataset metadata wrangle (a ~60-line dplyr
  pipeline) ran on every start. It is now memoised to `.crave-cache/*.rds` inside
  the dataset directory, keyed on the size and mtime of the source files plus the
  dataset name and its `citation_from_id` flag. Unchanged data skips the pipeline
  entirely. Falls back to the session tempdir if the dataset is read-only.
- **Deferred packages.** ~30 packages were attached before the first page was
  served. Only 19 are now, and the nine that a single tab needs (`cluster`,
  `ggdendro`, `igraph`, `umap`, `mice`, `ggVennDiagram`, `ggupset`, `DescTools`)
  are attached on first use.
- **Packages removed entirely.** `sqldf` (replaced by plain `DBI`, which also
  removes the `gsubfn`/`proto`/`chron` chain), `tidyverse` (only four of its
  members were used), `ggraph` (the network tab built a ggraph plot and then
  returned the visNetwork one instead, discarding it), `foreach` (every `%do%`
  loop is now vectorised), and ten packages that no code referenced:
  `heatmaply`, `leaflet`, `ggpointdensity`, `ggrepel`, `Hmisc`, `rmarkdown`,
  `optparse`, `ggdensity`, `tidyselect`, `shinydashboard`. Also the five
  Bioconductor packages, which nothing called.
- **`R.utils` stays installed but is never attached.** `data.table::fread()` hands
  compressed input to `R.utils::decompressFile()`, and every CRAVE metadata file is
  gzipped, so removing it stops any dataset from loading — but grepping the source
  for callers finds none, which makes it an easy thing to delete by mistake. It is
  now declared in `CRAVE_PKGS_INSTALLED_ONLY` in `R/00_packages.R`, which warns at
  start-up if it is absent, and the Dockerfile installs it in its own block with
  the reason attached.
- **Per-session work hoisted to app scope.** Every connected user was running
  `sort(unique(as.character(comparisons$<col>)))` for seventeen columns, in three
  separate places, and rebuilding the 600-element core-essentials vector. These
  depend only on the dataset and are computed once.
- **Exorcise probe.** Availability was tested with
  `system("docker images | grep ...")` on every start, which shells out, pipes
  through grep, and prints to the console. It now short-circuits when the
  reference directory holds no `.2bit` files and otherwise uses
  `docker image inspect` with output suppressed.
- **Gene symbol filtering** was eight sequential `grepl()` passes over the whole
  symbol table; it is one alternation, one pass.

### Query performance

- **SQLite indexes.** The `stat` table ships with only its primary key,
  `(comparison_id, gene_id, analysis_type_id)`. That serves Explore, which always
  filters by comparison first, but Correlate and Pendragonator filter by gene
  across many comparisons and so fell back to a full table scan. CRAVE now creates
  `(gene_id, analysis_type_id)` and `(analysis_type_id, comparison_id)` indexes on
  first start, and runs `ANALYZE`. Skipped silently if the dataset is read-only.
- **`fetchStat` (was `fetchDs`)** no longer runs `SELECT * FROM stat LIMIT 0`
  against every dataset per call to learn the column names, no longer re-filters
  the whole comparisons table once per dataset to route comparison IDs (a lookup
  vector does it), and no longer rounds the accumulated result inside the loop —
  which with N datasets rounded the first dataset's rows N times.
- **Bound-parameter chunking.** Gene and comparison lists are chunked at 400.
  SQLite caps bound parameters at 999 in many builds, so a query for a thousand
  genes could previously fail outright.
- **Enrichment.** The hypergeometric grid was a nested `foreach` over screens ×
  classes with four `length(unique(...))` scans of the whole table inside each
  cell — 2400 full scans for 30 screens and 20 classes. It is four grouped counts.
- **ROC.** The precision-recall sweep stepped through 2000 thresholds and
  re-filtered the whole query table at each one, per screen. It is now one sort
  plus a `cumsum` and a `findInterval`.
- **Overlap.** The per-pair universe was recomputed with a `group_by`/`summarise`
  over the whole table for every pair; it is now an intersection of precomputed
  per-screen gene sets. The Upset plot's set-membership construction was a nested
  per-row loop, quadratic in gene count; it is vectorised.
- **Pendragonator.** Best-row-per-gene was a grouped `rank()` with random tie
  breaking; it is an ordered `.SD[1]`, which is faster and reproducible.
- **Debouncing.** Typing in Correlate's free-text screen filter re-scanned the
  comparison table on every keystroke, and each change to the Explore comparison
  search fired seventeen selectize updates. Both are debounced.

### Bugs fixed

- **Enrichment FDR was never rounded and a spurious column appeared.**
  `mutate(adj.phyper = signif(p.adjust(...)), 4, ...)` had the `digits` argument
  stranded outside the `signif()` call, so `signif()` ran with its default and `4`
  became a column named `4`.
- **Overlap p-values were not well defined.** The overlap count came from the
  unrestricted hit lists while the population sizes came from the
  universe-restricted ones, so the overlap could exceed the sample size. The
  upper-tail probability was also `P(X > q)` rather than `P(X >= q)`, which is the
  convention the Enrichment tab already used (and which 4.27.2 fixed there).
- **The Upset plot crashed with fewer than two screens selected**, returning a
  variable that was never assigned on that path.
- **Gene query's "Cell line" violin colour scheme errored**, mapping to a
  `Cell line` column that does not exist (it is `Cell line (diff)`).
- **The clustergram's screen axis went blank** when the screen dendrogram was
  disabled: it fell back to the order of the *contrast filter* input, which is
  empty unless the user filtered by contrast, producing all-NA factor levels.
- **Chronos and Manual results showed as `NA`** in the bulk download and
  Pendragonator tables, which only decoded MAGeCK and DrugZ.
- **Save/Load never restored the Endpoint filter** (`key$comparison_selector_end`
  instead of `..._endpoint`).
- **Selecting a screen wiped every analysis tab's screen choice.** The per-analysis
  pickers were rebuilt with no `selected`, so changing the comparison selection
  reset them. This was also why Save/Load had to re-apply the whole key eight
  times via `shinyjs::click()` on its own button, hoping the cascade settled. The
  pickers now preserve valid selections and the retry loop is gone.
- **Seeding the global RNG.** The Home tab called `set.seed()` on the clock hour to
  pick a message of the day, and the network tab seeded it for a discarded layout.
  Both made the whole process's random stream deterministic — so two sessions
  started in the same hour generated *identical* Exorcise task ids and then
  overwrote and deleted each other's working files in the shared temp directory.
  Neither seeds now, and task ids also carry the process id and a microsecond
  timestamp.
- **Exorcise ignored the "Interference" chemistry**: the mode dictionary keyed that
  entry "Inhibition" while the UI offered "Interference", so the lookup returned
  `NA` and Exorcise received an empty `-q`.
- **Exorcise always claimed annotations had been supplied.** It compared
  `length(input$exorcise_orig) == length(input$exorcise_seq)`, and both are single
  strings from a textarea, so the test was always `1 == 1`. Exorcise was then told
  to harmonise against a column that did not exist when the box was empty.
- **`fetchDs` with an unresolvable contrast** did a full, ungated scan of every
  dataset's entire `stat` table. An empty comparison set now returns no rows.
- **`load()` shadowed `base::load`**, which is why the Save/Load tab had to spell
  out `base::load`. The dataset loader is now `loadCraveData()`.
- **The `Formaldehyde <- NULL` shim is gone.** It existed because `sqldf` resolves
  table names against the calling environment and collided with base R's built-in
  `Formaldehyde` dataset. Plain `DBI` has no such behaviour.
- `length(exorcise_root > 0)` (always 1) is now `length(exorcise_root) > 0`.
- Dead code removed: the `symbol` global was assembled from a field `loadOne()`
  never returned; `queryDs0()` was unused; `> Inf` guard branches could never be
  taken; `rm(genes_filtered)` targeted a local that did not exist.

### Security

- **Exorcise no longer builds a shell command by string concatenation.** The PAM
  and advanced-chemistry fields went straight into a line handed to `system()`, so
  a value like `NGG; rm -rf /` would have executed. The command is now an argument
  vector passed to `system2()` — no shell — and both fields are validated against a
  character allow-list.
- **The guide library query is parameterised.** It previously pasted user-supplied
  filter values directly into SQL.
- The guide library is opened read-only.

### Robustness

- `config.R` is sourced into its own environment and merged over a complete set of
  defaults, then validated: missing dataset paths, absent `database.db`, duplicate
  dataset names and unknown themes are reported on a startup page that says what is
  wrong, rather than surfacing as "object not found" from inside the UI code.
- Datasets missing required metadata columns are named in the error, and one bad
  dataset no longer prevents the others from loading.
- Analyses run inside `tryCatch` and report failures in the panel instead of
  greying out the session.
- Every long operation shows a progress message.
- Typed gene symbols that are not recognised are reported instead of silently
  dropped.
- Dataset connections are closed on shutdown and after a data refresh.

### Configuration

New optional `config.R` settings, all with defaults: `shiny_host`, `shiny_port`,
`max_upload_mb`, `enable_input_logging`. Input-change logging installed an observer
per input and wrote a line on every keystroke; it is now off unless asked for.

### Save/Load format

Keys are now RDS with a format marker and are keyed by module. Keys written by
CRAVE 4.x are detected and converted on load.

---

## Earlier history

| Author | Version | Date | Description |
|---|---|---|---|
| SL | 2.0.0 | 2023-08-03 | Complete rework. |
| SL | 2.1.0 | 2023-08-08 | Enable lazy-loading using SQL data structure. |
| SL | 2.1.1 | 2023-08-10 | Enable reference level KO in view. |
| SL | 2.1.2 | 2023-08-11 | Add busy indicator. Fix edge cases with screen selection. |
| SL | 2.1.3 | 2023-11-03 | Enable comparison contrast column support. |
| SL | 2.2.0 | 2023-11-03 | Enable case-insensitive search. |
| SL | 2.2.1 | 2023-11-03 | Change column order in comparison picker. |
| SL | 2.2.2 | 2023-11-03 | Add comment metadata feature. |
| SL | 2.2.3 | 2024-01-19 | Bug fixes. |
| SL | 2.3.1 | 2024-02-28 | Improve table and plot display. |
| SL | 3.0.0 | 2024-03-07 | Beta release CRAVE Correlate. |
| SL | 3.1.0 | 2024-03-07 | Add source data switcher. |
| SL | 3.1.1 | 2024-03-08 | Fix race condition between the two comparison pickers. |
| SL | 4.0.0 | 2024-06-07 | Add enrichment plot. Enable long style (DDRcs) dataset support. |
| SL | 4.1.0 | 2024-06-20 | Pre-populate gene list upon screen selection rather than showing them all. |
| SL | 4.2.0 | 2024-06-20 | Move enrichment plot to CRAVE Correlate. Enable bulk data download. |
| SL | 4.2.1 | 2024-06-26 | Support external dataset. |
| SL | 4.2.2 | 2024-06-27 | Fix bug where contrasts involving "Formaldehyde" couldn't be queried. |
| SL | 4.3.0 | 2024-06-28 | Enable external and internal datasets to be viewed together. |
| SL | 4.3.1 | 2024-07-01 | Enable html plot downloading. |
| SL | 4.3.2 | 2024-07-06 | Add MOTD, logging, code optimisations. |
| SL | 4.3.3 | 2024-07-09 | Increase tolerance before CRAVE refuses to plot a gene query plot. |
| SL | 4.3.4 | 2024-07-09 | Enable ontology colouring on UMAP. |
| SL | 4.3.5 | 2024-07-15 | Add heatmap of gene query results. |
| SL | 4.3.6 | 2024-07-16 | Add dendrograms to hits heatmaps. |
| SL | 4.3.7 | 2024-07-16 | Fix colour gradients for asymmetric divergent scales. |
| SL | 4.3.9 | 2024-08-16 | Add Pendragonator. |
| SL | 4.3.10 | 2024-08-19 | Prepare public-suitable CRAVE. |
| SL | 4.4 | 2024-08-20 | Security improvements. |
| SL | 4.4.1 | 2024-08-23 | Add citation filter to screen filter controls. |
| SL | 4.4.3 | 2024-09-16 | Re-enable customisation plot for CRAVE. |
| SL | 4.4.4 | 2024-09-21 | Improve gene and screen filtering options, enable brushing. |
| SL | 4.4.5 | 2024-11-28 | Deal nicely when data do not exist. |
| SL | 4.5 | 2025-01-16 | Add Exorcise. |
| SL | 4.6 | 2025-01-26 | Add support for Chronos. |
| SL | 4.6.1 | 2025-02-07 | Improve stability of Exorcise input form. |
| SL | 4.7 | 2025-05-17 | Add overlap analysis, minor tweaks to comparison picker. |
| SL | 4.8 | 2025-05-20 | Add comparison search by numerator and denominator levels. |
| SL | 4.8.1 | 2025-05-29 | Improve comparison search, add core essentialome shortcut. |
| SL | 4.8.2 | 2025-05-30 | Fix Correlate screen search for cell line. |
| SL | 4.9 | 2025-06-03 | Enable Exorcise upload CRISPick output, add guide library, remove the need for a dummy dataset, make switching between DDRcs/CRAVE easier. |
| SL | 4.9.1 | 2025-06-05 | Fix when Exorcise advanced chemistry string passed even when empty, move to using Exorcise Docker image. |
| SL | 4.10 | 2025-06-06 | Move Guide Library to SQL. |
| SL | 4.10.1 | 2025-06-06 | Clean up files after Exorcise to avoid collisions with other users. |
| SL | 4.10.2 | 2025-06-09 | Deal properly when internal dataset absent (i.e. in the case of DDRcs). |
| SL | 4.10.4 | 2025-06-11 | Implement externally bindable message banners on the front page. |
| SL | 4.10.5 | 2025-06-11 | Enable HTML in message banners. |
| SL | 4.10.6 | 2025-06-11 | Enable logging. |
| SL | 4.11 | 2025-06-17 | Adopt native Plotly for improved speed. |
| SL | 4.11.2 | 2025-06-18 | Make web links open in a new tab, improve Kinds. |
| SL | 4.12 | 2025-07-01 | Add libraries table, fix when 1 external experiment is selected and "clear form" is chosen, fix FGC accessions. |
| SL | 4.12.2 | 2025-07-04 | Add hypergeometric adjusted p-value to Enrichment and Overlap; fix bugs in Biplot. |
| SL | 4.13 | 2025-07-07 | Improve contrast picker in Correlate. |
| SL | 4.14 | 2025-07-09 | Add Explore ROC-AUC. |
| SL | 4.14.1 | 2025-07-09 | UI/UX improvements. |
| SL | 4.14.4 | 2025-07-30 | UI/UX improvements. |
| SL | 4.15 | 2025-07-31 | Save/load feature. |
| SL | 4.16 | 2025-10-09 | Enable score and rank cutoff in Correlate. |
| SL | 4.17 | 2025-10-10 | Implement upset plot. |
| SL | 4.18 | 2025-10-24 | Support "Manual" externally-analysed screens. |
| SL | 4.19 | 2025-10-27 | Enable Pendragonator for GOIs. |
| SL | 4.19.1 | 2025-11-06 | Add support for Endpoint selection for "Manual" method comparisons. |
| SL | 4.19.2 | 2025-11-18 | Add boxplot to gene query. |
| SL | 4.19.3 | 2026-01-27 | Enable enrichment for DDRcs. |
| SL | 4.20 | 2026-03-05 | Add ontology view. |
| SL | 4.21 | 2026-03-05 | Allow disabling dendrograms on the clustrgram. |
| AE | 4.22 | 2026-06-23 | Arbitrary dataset support via config.R; remove guides/libraries/ontology requirement. Custom branding (portal name, front-page text, footer) via config.R. |
| SL | 4.23 | 2026-06-23 | Intelligently enable/disable features that depend on ontology, library, guide library, and Exorcise available. |
| SL | 4.24 | 2026-06-23 | Deprecate genetypes, never implemented. |
| SL | 4.25 | 2026-07-06 | Basic config and dataset error checking. |
| SL | 4.26 | 2026-07-08 | Fix Exorcise data and work directories with config.R. |
| SL | 4.27 | 2026-07-22 | Add support for Exorcise 2.0.1. |
| SL | 4.27.1 | 2026-07-23 | Speed up check for specified Exorcise installation. |
| SL | 4.27.2 | 2026-08-04 | Fix inequality in computing upper-tail hypergeometric p-value. |
