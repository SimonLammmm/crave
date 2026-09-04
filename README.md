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
    19_ui.R          UI assembly
    20_server.R      module wiring
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

It is not an R parser: it strips comments and string literals and works on the text
that remains. It cannot see type errors or anything that depends on runtime values,
and the advisory sections have known false positives, described in their output.

See CHANGELOG.md for what changed in 5.1.2.
