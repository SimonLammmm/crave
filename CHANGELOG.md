# Changelog

## 6.0.1 — 2026-11-13

### Fixed

- **A comparison search field could only ever hold one value.** On the Explore
  tab, under Select comparison, picking a treatment, a cell line or any other
  search field left that field offering nothing but the value just picked, so a
  second could not be added and the OR the form implies was unreachable. The
  filtering was right all along — several values in one field have always been
  an OR — but each field's list of choices was rebuilt from the comparisons the
  whole form had already narrowed to, which for the field doing the narrowing is
  exactly the values picked in it.
- Each field's choices are now drawn with **that field's own filter lifted**: it
  offers everything still reachable once the other fields have had their say, so
  Camptothecin and Carboplatin can sit in Treatment (diff) together and the
  table shows the comparisons matching either. The other fields still narrow as
  they did, which is what keeps the form from offering combinations that match
  nothing.
- Values already picked stay in their field's list even when the rest of the
  form leaves them unreachable. selectize drops a selection that is not among
  the choices, so without this an over-narrow search would quietly rewrite
  itself instead of being there to undo.

### The ROC tab now measures what it says it measures

- **The false positive rate was divided by every gene in the query screen
  instead of by the non-hits.** The x axis therefore never reached 1, and the
  anchor point the curve was given at (1, 1) turned the gap into area: the
  reported AUC was roughly `(N/k)·AUC + P/k`. On a 221-gene worked example the
  tab read 0.762 and 0.804 where the answer is 0.665 and 0.662.
- **The direction and the cutoff were applied to the query screens as well as
  to the reference.** They now define the hits and nothing else. Every gene a
  query screen measured is swept, whichever way it scored, because a gene the
  reference depleted and the query screen enriched is a hit that screen failed
  to recover, not a gene to leave out of the arithmetic. Under the old filter
  20 of 60 positives disappeared from one example screen and 21 from the other,
  and a rescaling of the recall axis then drew the curve all the way to 1 as
  though nothing were missing.
- **Reference hits absent from a query screen are now counted and said out
  loud.** The legend reads "40 of 60 hits present" when a screen cannot reach
  them all, and a screen sharing no hits at all with the reference is named
  rather than silently given an AUC of 0.001.
- **Rank selected the weakest genes.** The ROC tab ranked by `rank(abs(score))`
  — rank 1 being the smallest magnitude — and pooled the query screens into one
  ranking, so a threshold meant a different thing in each. Hits are now chosen
  the way the Overlap tab chooses them, and query screens are ranked within
  themselves.
- **Ties, and the sweep itself.** The sweep ran over a fixed grid of 2,000
  thresholds and was integrated by `DescTools::AUC()`, which averages y wherever
  two points share an x and warns that it is doing so; the warning was
  suppressed. The sweep now steps from one distinct statistic to the next, tied
  genes crossing as a block, and is integrated directly. The result is the
  Mann-Whitney AUC with average ranks, to every digit.
- **The table under the plot is the sweep, gene by gene**, in the order the
  genes were reached, flagging which ones were reference hits and carrying the
  recall and false positive rate at each step. The old table held the curve's
  coordinates, which say nothing a reader can act on; this one names the genes
  recovered, and the trapezoid down its two rate columns is the AUC on the
  plot.

### A broader Correlate example, and a dataset it can show off on

- The Example button fills the form with **200 genes** rather than a handful,
  and selects five DNA-damaging agents at once — olaparib, talazoparib,
  cisplatin, etoposide and camptothecin — rather than one, since several values
  in one field are an OR. It sets a Score cutoff of 3, and it turns on the
  free-text filter with a regular expression that leaves out screen variants
  and alternative-library runs, which an underscore in a screen's name marks.
- Every Correlate tab now draws something from it. Two of them needed a tab
  control the example had no say over, so it now sets those too and Reset puts
  them back: the Network tab's similarity, because how many edges a threshold
  admits depends on how many screens are being correlated and the 0.7 default
  draws a hairball over few screens; and the Enrichment tab's cutoff, because
  at |normZ| >= 3 a handful of screens give too few hits for homologous
  recombination under a PARP inhibitor to survive multiple testing.

### A bigger example dataset

- **The dataset that ships with CRAVE is now thirty comparisons rather than
  two**, and carries an ontology. It is a subset of DDRcs: fifteen published
  experiments over twelve cell lines and eight libraries, with MAGeCK, DrugZ
  and — for the six essentialomes — Chronos scores for 1,500 genes, and 220
  ontology rows over nine DNA repair classes. Two comparisons could not
  demonstrate Adjacency, Network, Reduce or Enrichment at all; these can.
- `example-dataset/build-example-dataset.R` rebuilds it from any CRAVE dataset
  and documents what it selects and why. README has the summary.
- **The Exorcise worked example has moved to `example-analysis/`**, where it is
  no longer mixed up with the dataset CRAVE ships. The counts matrix, workbook
  and `commands.txt` are joined there by what they produce — the analysis
  results and the CRAVE dataset built from them — so every stage can be read
  without running anything. The commands themselves are now Exorcise 3.1.9's
  `exorcise-analyse` and `exorcise-database`, in place of 1.6's
  `crispr_pipeline` and `crispr-screen-viewer database`. The whole folder is in
  `.dockerignore`: it is a worked example, not something the image needs.

## 6.0.0 — 2026-10-08

Two new tabs write to a dataset rather than reading it, three parts of CRAVE can
now be withheld by configuration, the Network tab can find communities, and two
Correlate tabs are renamed — labels, module ids and all, which is the one breaking
change here and is described below. Major, because a CRAVE that can change its own data is
a different thing to deploy than one that cannot: the new tabs are off by default,
and the notes below on who may reach them are the reason.

### Community detection, on the Network tab

**The Network tab can group the genes it draws and ring each group.** A
"Community detection" dropdown offers Agnes under each of its six linkages,
Clara, k-means, non-negative matrix factorisation, a Gaussian mixture model,
Louvain and Leiden. "Off" is the default and leaves the tab exactly as it was: a
visNetwork to drag about, and an adjacency list of pairs.

- **The similarity and the distance are separate choices**, because they answer
  separate questions. The similarity measure — Spearman, Pearson, Kendall or
  cosine, all four supported — decides which genes are joined by an edge. The
  distance measure — euclidean, manhattan or jaccard — decides which genes are
  near each other for the purpose of grouping, as it already does for the
  dendrograms. A network of Spearman edges can therefore be grouped by jaccard
  distance.
- **Every method but the two graph ones honours the distance**, though not all of
  them can take one: Agnes and Clara are given it directly, k-means and the
  mixture model are given a classical MDS embedding of it, and the factorisation
  is given the affinity `1 - d/max(d)`. Each is said in the method's own comment
  rather than left to be discovered.
- **Louvain and Leiden group the network itself.** They are given the edges the
  tab draws, weighted by the absolute similarity exactly as the layout already
  weights them, and maximise modularity over that graph — so the similarity
  measure and the cutoff decide the communities by themselves, and no distance
  reaches them. They also find their own number of communities, so they take a
  **resolution** (default 1, above it for more and smaller communities, below it
  for fewer and larger) where the others take a k, and the diagram reports how
  many they found. Leiden is run against modularity rather than its own default
  constant Potts model, which reads the resolution as an absolute edge density and
  would return every gene as its own community at the default; under modularity
  the resolution means for Leiden what it means for Louvain, so the two can be
  compared. Both are seeded, so the same selection gives the same communities.
- **The controls shown follow the method.** k and the distance appear for the
  methods that take them, the resolution for Louvain and Leiden, and the
  membership threshold for the two probabilistic fits — rather than all of them
  appearing as soon as a method is chosen, because an inert control beside a
  dropdown invites the reading that it does something.
- **k** is asked for by every method that needs one — Agnes needs it to cut its
  tree, the factorisation as its rank, the mixture as its number of components. It
  is refused outside `2 <= k < genes`, with the bound in the message. The
  resolution is refused outside `0.05` to `10` the same way.
- **A resolution that leaves every gene alone is refused rather than drawn.** It
  is a valid partition and a useless picture, and on a network whose edges are all
  of similar weight it arrives abruptly; the message says which resolution did it
  and what to change.
- **A distance method and a graph method can disagree on the same selection**, and
  neither is wrong. The first groups genes whose scores are alike across the
  screens; the second groups genes the drawn network ties together, so genes that
  share a profile but produce few surviving edges are one community to Agnes and
  several to Louvain. Only the second is describing the picture on screen.
- **The factorisation and the mixture model are probabilistic.** Each gives every
  gene a membership in every community, so a gene can belong to two and the
  ellipses overlap. **Least community membership to be ringed** (default 0.7)
  decides how much of a gene must sit in a community for its ellipse to enclose
  it; the control appears only for those two methods, because for the others it
  could do nothing.
- **The adjacency list** keeps both triangles, as it already did — one row per
  ordered pair — and gains the community of each end, plus the membership of each
  end as a probability under a probabilistic method.
- **The diagram** is the same graph, laid out by Fruchterman-Reingold, with each
  community ringed by a 95% ellipse and each gene coloured by the community it
  sits most in. It is drawn through ggplot and plotly rather than visNetwork,
  which has no way to draw a region around a set of nodes; the trade is that the
  picture is no longer dragged about with the mouse, and hovering still names the
  gene, its community and its membership. With community detection off, the
  visNetwork is untouched.
- Communities are detected over the genes the network actually draws, not every
  gene in the score matrix: a community containing genes the network does not show
  could not be ringed, and would describe a different picture to the one on
  screen.
- The mixture model is fitted with diagonal covariances in at most ten dimensions,
  and its variances are floored. A full covariance per community over every screen
  in a selection has more parameters than a gene set has genes, and an
  unconstrained EM run ends by collapsing a component onto two identical genes
  rather than by finding a cluster.
- Every randomised step — both soft fits, the k-means start, the two graph sweeps
  and the layout — is seeded through `withSeed()`, so the same selection draws the
  same communities and the global random stream is left as it was found.

### Jaccard distance

- **Jaccard is now offered on the Heatmap and Adjacency tabs too**, alongside
  euclidean and manhattan. `agnes()` implements only the other two, so it is now
  handed a dissimilarity matrix rather than the data and a metric name — for
  euclidean and manhattan that is the same calculation it was doing itself, which
  is checked against `stats::dist()`.
- What jaccard means here is its continuous generalisation, the Tanimoto
  coefficient, so that it is defined on signed CRISPR scores; on binary values it
  reduces to the familiar ratio of set sizes, which is tested. `clara()`
  implements its own jaccard in C and is given the metric by name, so Clara's
  jaccard is cluster's and every other method's is CRAVE's.
### Breaking: Clustergram is now Heatmap, and Heatmap is now Adjacency

Renamed all the way down, not just on the tab. Leaving the old ids in place would
have meant `heatmap` naming the Adjacency tab while a tab labelled Heatmap sat
next to it under a different name — a double meaning that costs every future
reader of this module a detour, in exchange for a saved-key format that is
largely unused.

| | was | is |
|---|---|---|
| tab label | Clustergram | Heatmap |
| tab label | Heatmap | Adjacency |
| module id | `clustergram` | `heatmap` |
| module id | `heatmap` | `adjacency` |
| input ids | `cg_*` | `hm_*` |
| input ids | `hm_*` | `adj_*` |
| function | `runClustergram()` | `runHeatmap()` |
| function | `runHeatmap()` | `runAdjacency()` |

**What this does to a session key saved by 5.1.5 or earlier.** It still loads, and
everything outside those two tabs is restored as before: `restoreInputs()` matches
on id, so ids it no longer knows are skipped. The two tabs' clustering controls —
nine inputs — come back at their defaults instead of their saved values.

One case needed code rather than a note. An old key's `hm_dendro`, `hm_dist` and
`hm_k` described the tab now called Adjacency, and those are exactly the ids the
Heatmap tab uses now. Restoring them would have been worse than skipping them: the
values are valid for the control they would land on, so a user's old Adjacency
settings would silently reappear on a different tab. A key from that era is
recognisable by its `cg_*` ids, so the Correlate module drops those three when it
sees them, and says so in the log. No key is rejected, and nothing is restored
into the wrong tab.

Deployers who set `defaults_correlate_customise_heatmap_*` in config.R are
unaffected: those names are unchanged, because they are the colours of both tiled
plots rather than of either tab. Their labels in the Customise dialog now say so.

### Submit

**Upload a counts matrix and an experimental design workbook, run the analyses the
workbook asks for, and add the results to a dataset.** The work is Exorcise's:
`exorcise-analyse` runs MAGeCK, DrugZ and Chronos — whichever the Analyses sheet
names, which is what "where appropriate" comes to — and `exorcise-database` puts
the results and the metadata describing them into the dataset directory. CRAVE
contributes the form, the checks in front of it, and the job that sequences the
two.

- The workbook is read as soon as it is uploaded and what CRAVE understood is
  shown back: the experiment id, the methods, the control groups, the counts file
  the workbook names, and whether that experiment is already in the dataset.
  Refusing a wrong file takes a second, rather than a container start and a
  confusing log.
- The counts matrix is saved under the name the workbook declares, not the name it
  was uploaded with. `exorcise-analyse` resolves counts by the declared name
  inside `--counts-dir`, so a file saved under any other name is simply not found.
- **The one adjustment to that name is compression**, because Exorcise decides
  whether to decompress from the extension and `maybe_its_gz()` will reach for a
  `.gz` the workbook did not mention but never strips one. So a gzipped upload
  declared as `.tsv` is saved as `.tsv.gz` and found by that fallback, and a plain
  upload declared as `.tsv.gz` is compressed on the way in, since nothing would
  look for it uncompressed. Compression is read from the file's first two bytes,
  not from what the browser called it. Without this a gzipped upload was stored
  under a `.tsv` name and Exorcise stopped at validation with
  `UnicodeDecodeError: 'utf-8' codec can't decode byte 0x8b`; uploading the
  uncompressed matrix was the only way through.
- The counterpart name is removed when the matrix is staged, so a file left by an
  earlier submission cannot be resolved ahead of the one just uploaded — which is
  what otherwise makes a resubmission fail the same way as the first attempt.
- A workbook that sets "Analysis version" is refused, with the reason:
  `exorcise-analyse` then files its tables under
  `results/<experiment>/<version>/`, which is not where `exorcise-database` looks,
  and the experiment would be added with no statistics in it.
- Analyses run detached, through a shell script, and the panel polls a log. A
  Chronos run is minutes to hours and R is single threaded: waiting on it inside
  an observer would freeze every session connected to the app, not just the one
  that pressed the button. Closing the browser does not stop a job, and reopening
  the tab picks it back up.
- **The raw analysis output can be downloaded**, as a `.tar.gz` of everything
  MAGeCK, DrugZ and Chronos themselves wrote for one experiment — the method
  directories `exorcise-analyse` files them under — which unpacks as
  `<experiment>/<method>/…`. The sgRNA summaries, the per-comparison files and
  Chronos's model output are otherwise unreachable once the container has exited.
- **The tabulated tables are not in it.** They are what the dataset is built
  from, and so behind every plot and every "Download data" elsewhere in CRAVE;
  putting them in the archive would hand back the same numbers a second time
  under a different name.
- The picker lists the experiments in the dataset that have raw output, not every
  experiment: one analysed only as a manual table computed nothing and has none.
  It is read from the directory rather than from the metadata, which is also what
  keeps an experiment id typed into a workbook from being used to build a path.
  The experiment just submitted is the one offered, and the list refreshes when a
  job stops — whether it succeeded or not, since a failed run can still have left
  behind the output worth looking at.
- The archive is built with R's own tar, so it needs no `tar` binary in the
  image, and the size it will be before compression is shown next to the button.

### Manage

**What a dataset holds, and the three ways it may be changed.** Relabel its
metadata, rebuild an experiment from a corrected workbook, or remove experiments.

- **Rebuild** republishes an experiment from a corrected workbook without running
  its analyses again, because the result tables from its submission are kept.
  This is the route for correcting a contrast, a treatment or a citation at
  source — the workbook stays the thing the dataset is derived from.
- **Remove** drops the statistics and the metadata rows together.
- **Relabel** replaces the metadata tables with edited copies, and is the only
  write CRAVE makes to a dataset itself. It is refused unless the incoming tables
  name exactly the experiments and comparisons already present, with each
  comparison still under its own experiment. Labels may change; membership may
  not. Adding or dropping a row has a correct route, and this is not it. The
  previous tables are kept under `.crave-submit/metadata-backups/`.
- A consistency line at the top of the tab compares the metadata against the
  statistics actually in the database, and is computed from the files each time it
  is drawn. A dataset built and changed only through CRAVE cannot drift, so
  anything it reports came from outside: a hand-edited file, an interrupted
  build, or a database written by something else.

### How the two stay in step

`exorcise-database` writes `database.db` and both metadata CSVs in one pass, and
it is the only thing CRAVE lets near them. Adding, rebuilding and removing are all
invocations of it; relabelling is checked against the database before it is
allowed. That is the whole of the argument, and it is why the consistency line
has nothing to do under normal use.

One job at a time per dataset, enforced by a lock file naming the job that holds
it. A job that dies without releasing the lock is reported as stale rather than
busy, and Manage offers to clear it: a lock that outlives its job would otherwise
make the dataset permanently unchangeable.

### Configuration

- **`enable_submit`** and **`enable_manage`** (both `FALSE`) control the new tabs.
  **Neither tab authenticates anybody** — CRAVE does not know who its visitors
  are, so what it offers, it offers to everyone who can reach it. Leave them off
  on a public deployment and put an authenticating proxy in front of one where
  they are on.
- **`enable_pendragonator`** (`TRUE`) withholds the Pendragonator sub-tab of
  Correlate, including its line of sidebar help.
- **`hostpath`**, per dataset, names the dataset directory as the Docker daemon
  sees it. CRAVE reaches the daemon through its socket, so bind mounts are
  resolved on the host rather than inside CRAVE's container — the same thing
  `exorcise_hostroot` already does for the Exorcise tab. Running natively, leave
  it out.
- Both tabs also require `exorcise_docker` to name a runnable image and the
  dataset directory to be writable. A tab whose requirements are not met stays
  hidden, and the startup log says which one was missing rather than leaving it to
  be guessed at.

### Other changes

- **Exorcise 3.1 or newer is required** for Submit and Manage, which call
  `exorcise-analyse` and `exorcise-database`. The Exorcise tab is unaffected.
- `readxl` is a new deferred dependency, loaded the first time a workbook is read.
- `19_ui.R` and `20_server.R` are now `29_ui.R` and `30_server.R`, so that
  assembly still comes last with the new modules in between.
- `loadCraveData()` takes an optional `config`, which is what the capability flags
  of the withheld tabs are computed from. Omitting it still works and reports them
  unavailable, which for a console session they are.
- `tools/rcheck.py` finds modules at any two-digit number rather than 10 to 18, so
  the new ones are covered by the id cross-reference.

---

## 5.1.5 — 2026-09-30

- **Fixed: accented characters in the metadata files depended on the locale the
  app was started in.** `readTableGlob()` called `fread()` without an `encoding`,
  so it took fread's default, `"unknown"`: the bytes come back intact but the
  strings are left unmarked, and R then reads them as whatever the process locale
  says. On a UTF-8 host that is correct and nothing looks wrong. Started without a
  locale — a bare `Rscript`, a systemd unit, a shiny-server install with no
  `LANG` — it is silently wrong at the last step rather than the first: the bytes
  survive the read, `wrangleMetadata()`, `stripHtml()` and the cache, and are
  then destroyed on the way out, because `jsonlite` escapes each unmarked
  non-ASCII byte. DT and selectize send the browser
  `A<cc><81>lvarez-Quilo<cc><81>n et al (2020)` for `Álvarez-Quilón et al
  (2020)`. Every file read through this helper — both metadata files, the
  ontology and the library table — is written by `exorcise-database`, which
  writes UTF-8, so UTF-8 is now declared, and an accented citation renders
  identically under `LC_ALL=C` and `en_GB.UTF-8`. The parameter is exposed, so a
  dataset that genuinely is not UTF-8 can still say so.

  Values out of the database needed nothing: RSQLite marks its strings UTF-8
  already, under any locale.

- **`CRAVE_CACHE_VERSION` is 4L.** The metadata cache stores wrangled strings, so
  a cache written before this change holds unmarked ones, and the fingerprint is
  otherwise a function of the source files alone — a deployment that upgraded
  without touching its metadata would have gone on serving the unmarked copies.

## 5.1.4 — 2026-09-21

- **Added `build-and-push.sh`.** Builds, tags, pushes and optionally saves the
  image. `linux/amd64` only: CRAVE is served over HTTP from a server and nothing
  it depends on is architecture-sensitive, unlike Exorcise, which needs a
  multi-architecture manifest because BLAT and twoBitToFa are compiled binaries
  and its users run it on Apple silicon. Building one architecture keeps this to
  minutes instead of the hours QEMU would take. The script reads the version
  from `01_constants.R` rather than holding its own copy, so the image tag
  cannot disagree with what the running app reports.

- **Added `.dockerignore`.** The Dockerfile copies `shiny-server/` wholesale, so
  anything sitting in that directory at build time goes into the image. No
  released version has shipped a `config.R` — there is none in 5.1.2 — but a
  working tree acquires one as soon as you configure a local instance to test
  against, and building from that tree would bake it in. An image carrying one
  is also worse than an image without: `09_config.R` explains what to do when
  the file is missing, whereas a stale baked copy makes a container started
  without the bind mount come up quietly wrong.

  The example dataset, the images directory and assorted development leftovers
  are excluded too, which keeps the build context to the source rather than the
  8 MB the repository currently holds.

- **Version references brought back into line.** `01_constants.R` said 5.1.3
  while `docker-compose.yml` still pulled 5.1.2, and `app.R` and `README.md`
  both pointed at the 5.1.2 notes. The compose file now tracks the release and
  the two prose references no longer name a version, so they cannot go stale
  again.

## 5.1.3 — 2026-09-09

Reconstructed from the diff against 5.1.2; this release originally shipped
without an entry.

- **Fixed: selecting comparisons silently did nothing for most screens.** The
  "Selected comparisons" picker is created with `server = TRUE`, so the browser
  only holds the options it has been sent — the first
  `CRAVE_MAX_SELECTIZE_SERVER` at start-up, plus whatever a later search
  fetched. selectize cannot resolve a value whose option it has never seen and
  clears the control instead, so `updateSelectizeInput(selected = ...)` worked
  for screens near the top of the list and failed silently for everything
  below. Every write now goes through `setContrasts()`, which resends the
  choices alongside the selection. This affected selecting rows in the
  comparisons table, Select all from a search, and both Clear buttons.

- The remaining hardcoded `maxOptions = 1000` was replaced with
  `CRAVE_MAX_SELECTIZE_SERVER`, which is what every other picker already used.

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
