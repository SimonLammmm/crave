# CRAVE
Contrasts-based Resource for Analysis, Visualisation, and Exploration

## Abstract
CRAVE enables exploration of single contrasts between a pair of conditions in a biological experiment as well as analysis of trends across contrasts.

## Deployment
CRAVE can be deployed as a Docker container or natively in Rshiny.

### Docker Compose (recommended)
Clone this repo and navigate into it.
```
git clone https://github.com/SimonLammmm/crave.git
cd crave
```
Run CRAVE with Docker Compose
```
docker compose up -d
```
A preview of CRAVE is now accessible at http://localhost:8080.

Customise CRAVE by copying the shiny-server/config.example.R file as config.R and editing it. Update the docker-compose.yml file with any additional bind mounts you need.

### Docker
Clone this repo and navigate into it.
```
git clone https://github.com/SimonLammmm/crave.git
cd crave
```
Build the Docker image.
```
cd docker
docker build -t crave:latest .
```
Run CRAVE with Docker
```
docker run --rm crave:latest -v "./shiny-server/config.example.R:/app/config.R" -v "./example-dataset/data/:/data/dataset-1/" -p 8080:3838
```
A preview of CRAVE is now accessible at http://localhost:8080.

Customise CRAVE by copying the shiny-server/config.example.R file as config.R and editing it. At runtime, bind mount that instead of config.example.R. Adjust the command with any additional bind mounts you need.

### Rshiny
Run in native Rshiny on your local machine.
Clone this repo and navigate into it.
```
git clone https://github.com/SimonLammmm/crave.git
cd crave
```
Copy the shiny-server/config.example.R file as shiny-server/config.R.
```
cp shiny-server/config.example.R shiny-server/config.R
```
Install the dependencies listed in `shiny-server/R/00_packages.R`:
`CRAVE_PKGS_EAGER`, plus `CRAVE_PKGS_LAZY` if you want every tab, plus
`CRAVE_PKGS_INSTALLED_ONLY` — which CRAVE never loads itself but which
`data.table::fread()` needs in order to read the gzipped metadata files. Then run
Rshiny.
```
Rscript -e "shiny::runApp('shiny-server')"
```
A preview of CRAVE is now accessible at http://localhost:3838.

Customise CRAVE by editing the shiny-server/config.R file.

## Layout

```
shiny-server/
  app.R              entry point; sources R/ in numeric order
  config.example.R   template configuration, copy to config.R
  R/
    00_packages.R    which packages load at start-up and which load on demand
    01_constants.R   version, limits, defaults, static reference data
    02_utils.R       shared helpers
    03_data.R        dataset loading, metadata cache, SQLite indexing
    04_query.R       parameterised statistics queries
    05_analysis_common.R    primitives shared between analyses
    06_analysis_*.R  the analyses, as pure functions of explicit arguments
    09_config.R      config.R loading and validation
    1*_mod_*.R       one Shiny module per tab, plus two reusable ones
    21_pipeline.R    driving Exorcise to change a dataset on disk
    22_mod_submit.R  Submit
    23_mod_manage.R  Manage
    24_communities.R community detection for the Network tab
    29_ui.R          UI assembly
    30_server.R      module wiring
```

Each tab is a Shiny module, so its input ids are namespaced and local to it. The
analysis functions take explicit arguments rather than Shiny's `input` object, so
they can be called from the console:

```r
for (f in list.files("shiny-server/R", full.names = TRUE)) source(f)
config <- loadCraveConfig("shiny-server/config.R")$config
data   <- loadCraveData(config$datasets)
res    <- plotVolcano(data, contrast = data$comparisons$FriendlyID[1], method = "DrugZ")
res$plotdata
```

## The example data

Two folders, for the two directions. `example-dataset/` holds a dataset to look
at; `example-analysis/` shows where such a dataset comes from. Neither is in the
image — `.dockerignore` keeps both out of the build context — and both are
mounted or read from the host.

### example-dataset/

`example-dataset/data/` is what `docker-compose.yml` mounts, and what CRAVE
shows anyone who runs it without configuring a dataset of their own. It is a
subset of [DDRcs](https://sjlab.cruk.cam.ac.uk/ddrcs/), the lab's portal of
published DNA damage response screens — real screens rather than invented
numbers, so the first plots a new user sees mean something.

```
30 comparisons   15 experiments   12 cell lines   8 libraries
   24 treated-versus-untreated contrasts: olaparib, talazoparib, cisplatin,
      etoposide and camptothecin, in wild-type and in knockout backgrounds
    6 essentialomes from those same experiments
1500 genes       MAGeCK, DrugZ and (for the essentialomes) Chronos
 220 ontology rows over nine DNA repair classes
```

Both axes are trimmed, because whole libraries across thirty screens would be a
database of a few hundred MB rather than 19 MB. The genes kept are the ones
Correlate's Example button asks for, the ones the ontology names, and a random
sample of those measured in most of the screens. Nine of the thirty comparisons
are what that Example button selects, which is enough for every Correlate tab:
Adjacency needs two screens, Reduce four.

`build-example-dataset.R` in that directory rebuilds it, and documents the
selection rules at the top. The source it subsets is not shipped — DDRcs is
25 GB — but the script will take any CRAVE dataset and keep whatever of it
matches:

```
cd example-dataset
Rscript build-example-dataset.R /path/to/a/crave/dataset
```

### example-analysis/

A self-contained worked example of the other direction: how a counts matrix and
an experimental design become a CRAVE dataset, which is what the Submit tab does
from inside the app. Every stage is there to read without running anything:

```
example.counts.tsv.gz   a published olaparib screen, downsampled
example.xlsx            the experimental design: samples, contrasts, analyses
commands.txt            the two Exorcise commands that turn one into the other
results/                what exorcise-analyse wrote: DrugZ scores per contrast,
                        and the tabulator's combined table
data/                   the CRAVE dataset exorcise-database built from those
```

Re-running the commands overwrites `results/` and `data/`; `commands.txt`
explains the flags that make that work unattended.

## Community detection

The Network tab can group the genes it draws and ring each group on the diagram.
Two choices are made separately, and the tab asks for both:

- the **similarity measure** (Spearman, Pearson, Kendall, cosine) decides which
  genes are joined by an edge, and is what the adjacency list reports;
- the **distance measure** (euclidean, manhattan, jaccard) decides which genes are
  near each other for the purpose of forming communities, exactly as it does for
  the dendrograms on the Heatmap and Adjacency tabs.

Keeping them apart is what lets a network of Spearman edges be grouped by jaccard
distance. Most methods honour the distance, though they reach it differently:
Agnes and Clara take it directly, k-means and the mixture model are given a
classical MDS embedding of it, and the factorisation is given the affinity
`1 - d/max(d)`. Each is noted in that method's comment in `24_communities.R`.

Louvain and Leiden are the exception. They are graph methods: they partition the
network on screen rather than a distance between gene profiles, so the similarity
measure and the cutoff decide the communities by themselves and no distance is
involved. They also find their own number of communities, so they are given a
**resolution** rather than a k — above 1 for more and smaller communities, below 1
for fewer and larger — and the diagram says how many they found.

| Method | Communities per gene | Grouped on | Number of communities |
|---|---|---|---|
| Agnes (Ward, average, single, complete, flexible, weighted) | one | the distance | k, as a cut height |
| Clara, k-means | one | the distance | k |
| Non-negative matrix factorisation, Gaussian mixture model | a membership in each | the distance | k |
| Louvain, Leiden | one | the network's edges | found, from a resolution |

The factorisation and the mixture model are probabilistic: each gene gets a
membership in every community, so a gene can be ringed by more than one and the
ellipses overlap. **Least community membership to be ringed** (default 0.7) is
what decides how much of a gene has to sit in a community for its ellipse to
enclose it. The results table carries those memberships as probabilities; under
the hard methods it carries the community number alone.

The two families answer different questions, and on the same selection they can
disagree for good reason. A distance method groups genes whose scores are alike
across the screens; Louvain and Leiden group genes the drawn network ties
together, so a set of genes that share a profile but produce few surviving edges
will be one community to Agnes and several to Louvain. Both are describing the
data; only the second is describing the picture.

A resolution that leaves every gene in a community of its own is refused rather
than drawn, with the resolution named in the message: it is a valid partition and
a useless picture, and on a network whose edges are all of similar weight — two
screens, say, where every Spearman correlation is exactly 1 — it arrives abruptly.

Edges are weighted by the absolute similarity, as the network layout already is,
so a strongly anti-correlated pair is as firm a tie as a strongly correlated one.
Leiden is run against modularity rather than its own default constant Potts model,
which reads the resolution as an absolute edge density and would return every gene
as its own community at 1; under modularity the resolution means for Leiden what
it means for Louvain, so the two can be compared.

The adjacency list is long-format and keeps both triangles — one row per ordered
pair, so `A, B` and `B, A` both appear — with the community of each end.

Jaccard is the continuous generalisation (Tanimoto) of the set form, so that it is
defined on signed CRISPR scores and reduces to the familiar ratio exactly when the
values are binary. `clara()` implements its own jaccard in C and is given the
metric by name, so Clara's jaccard is cluster's and every other method's is
CRAVE's; they measure the same thing and may differ in detail.

## Submitting and managing data

Two tabs write to a dataset rather than reading it, and both are off until
`config.R` turns them on:

- **Submit** takes a counts matrix and an experimental design workbook, runs the
  analyses the workbook asks for — MAGeCK, DrugZ, Chronos, or a manual table —
  and adds the results to a dataset. It also hands back the raw output of those
  analyses, as described below.
- **Manage** relabels a dataset's metadata, rebuilds an experiment from a
  corrected workbook, and removes experiments.

Both need `exorcise_docker` to name an image CRAVE can run (3.1 or newer, for
`exorcise-analyse` and `exorcise-database`) and a dataset directory that is
writable. A tab whose requirements are not met stays hidden, and the startup log
says which one was missing.

**Neither tab authenticates anybody.** CRAVE does not know who its visitors are,
so anything it offers, it offers to everyone who can reach it. Leave both off on
a public deployment, and put an authenticating proxy in front of one where you
turn them on.

### Why the database and the metadata cannot drift apart

A CRAVE dataset is three files that have to agree: `database.db`, and the two
metadata tables beside it. `exorcise-database` writes all three in one pass, and
it is the only thing CRAVE lets near them — adding an experiment, rebuilding one
and removing one are all invocations of it.

The single exception is replacing the metadata tables with edited copies, which
CRAVE does itself. That is checked against the database first and refused unless
the incoming tables name exactly the experiments and comparisons already there,
with each comparison still under its own experiment. It can rename, and it can
do nothing else; adding or dropping a row has to go through Submit or Remove.

The Manage tab shows the result of that argument as a consistency line: how many
comparisons the metadata describes, and whether every one of them has statistics.
It is derived from the files each time it is drawn, so it reports what is there
rather than what should be.

### What is kept, and where

Submissions are kept inside the dataset, under `.crave-submit/`:

```
.crave-submit/
  workbooks/   the experimental design each experiment was added with
  counts/      the counts matrices, under the names the workbooks give them,
               with .gz added or dropped so the name matches the bytes
  results/     one directory per experiment: <experiment>/tables/, where
               exorcise-database looks, and <experiment>/<method>/, the raw
               output Submit offers for download
  jobs/        one directory per run: the script, its log, its exit status
  lock         names the job currently allowed to change this dataset
```

Each experiment's results directory holds both halves of what exorcise-analyse
produced: one directory per method, named after it, holding what that method
itself wrote, and `tables/`, holding the one CSV per method that the tabulator
reduced all of it to. `tables/` is what the database is built from.

Submit offers the method directories as a `.tar.gz`, one experiment at a time,
unpacking as `<experiment>/<method>/…`. The tables are deliberately left out:
they are already in the database, and so behind every plot and every "Download
data" elsewhere in CRAVE, whereas MAGeCK's sgRNA summaries, DrugZ's
per-comparison files and Chronos's model output are otherwise unreachable once
the container has exited. An experiment analysed only as a manual table has no
such output and is not offered.

Keeping the workbooks and results is what lets Manage rebuild an experiment
without running its analyses again. They are also what a rebuild needs, so a
deployment that prunes `.crave-submit/` to save space gives up rebuilding for the
experiments it prunes; Submit still works, with the counts matrix supplied again.

Analyses run detached, so a long Chronos run does not block the app, and one job
at a time per dataset — the lock is what stops two submissions rewriting the same
database at once. Closing the browser does not stop a job; reopening the tab
picks it back up.

## Performance notes

On first start against a given dataset, CRAVE writes two things into the dataset
directory:

- `.crave-cache/meta-<dataset>.rds`, the wrangled metadata, keyed on the size and
  modification time of the source files. Delete it at any time; it will be rebuilt.
- two indexes on the `stat` table, which make gene-led queries (Correlate,
  Pendragonator) index lookups rather than full scans.

Both are skipped, with a note in the log, if the dataset directory is read-only —
CRAVE still works, just more slowly. Making the directory writable once is worth it.

## Static checks

`tools/rcheck.py` stands in for an R linter, for editing environments where R is not
installed. It needs only Python 3 and no packages.

```
python3 tools/rcheck.py            # everything
python3 tools/rcheck.py --quiet    # just the checks that can fail
```

It exits 1 if a check marked `[HARD]` fails, so it works as a pre-commit hook or a
CI step:

- **[HARD]** every file's brackets balance
- **[HARD]** every `plotPanelUI()` has a matching `plotPanelServer()` and vice versa
- **[HARD]** the three package manifests in `R/00_packages.R` match the three
  `install.packages()` blocks in the Dockerfile
- **[ADVISORY]** calls to names that are neither defined in the sources nor on the
  script's allowlist of base and package functions
- **[ADVISORY]** per module, input ids read by the server against `ns()` ids
  declared in the UI
- **[ADVISORY]** columns created and then re-read inside one `mutate()` or
  `summarise()`, where dplyr's sequential evaluation means the later argument sees
  the *new* value

`tools/check-icons.R` is the companion that does need R. It checks every `icon()`
and `submitIcon` name against the Font Awesome set shiny ships, and exits 1 on one
that is missing. Worth running after any icon change: Font Awesome Pro has names
the free set does not, so a plausible name can leave a tab with no icon and only a
single start-up warning to say so.

```
Rscript tools/check-icons.R
```

rcheck.py is not an R parser: it strips comments and string literals and works on
the text that remains. It cannot see type errors or anything that depends on runtime values,
and the advisory sections have known false positives, described in their output.

See CHANGELOG.md for what changed in this release.
