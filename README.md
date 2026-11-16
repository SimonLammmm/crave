# CRAVE
Contrasts-based Resource for Analysis, Visualisation, and Exploration

## Abstract
CRAVE enables exploration of single contrasts between a pair of conditions in a biological experiment as well as analysis of trends across contrasts.

## Evaluation
(DDRcs)[https://sjlab.cruk.cam.ac.uk/ddrcs/], the DNA damage response CRISPR screen explorer, runs on CRAVE. Visit it to evaluate features without
deploying. DDRcs may be kept a few versions behind this repo and has some features disabled.

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

## Configuration
CRAVE ships with an example config file at `shiny-server/config.example.R`. Use
this file to point to CRAVE-compatible datasets, the Exorcise image if using, and
to enable and disable modules. Instructions are in the `config.example.R` file.

## Example data
The `example-dataset/` shipped with CRAVE is loaded when you run CRAVE out of the
box in Docker compose. Use this to demonstrate the features on CRAVE for evaluation.

The `example-analysis/` folder contains Exorcise commands to create your own CRAVE-
compatible dataset. CRAVE ships with example input files to demonstrate running
Exorcise. To run them on your own files, you'll need a counts matrix and an
experimental design .xlsx workbook, in the format provided in the example files.

### example-dataset/
Contains a subset of the DDRcs dataset.
```
  30 comparisons   15 experiments   12 cell lines   8 libraries
  24 treated-versus-untreated contrasts: olaparib, talazoparib, cisplatin,
     etoposide and camptothecin, in wild-type and in knockout backgrounds
   6 essentialomes from those same experiments
1500 genes       MAGeCK, DrugZ and (for the essentialomes) Chronos
 220 ontology rows over nine DNA repair classes
```

### example-analysis/
A self-contained worked example of how to create your own CRAVE-compatible
database using the (Exorcise)[https://github.com/SimonLammmm/exorcise] Docker image.

```
docker pull simonlammmm/exorcise:latest
```
CRAVE ships with these example files.
```
example.counts.tsv.gz   a published olaparib screen, downsampled
example.xlsx            the experimental design: samples, contrasts, analyses
commands.txt            the two Exorcise commands that turn one into the other
results/                what exorcise-analyse wrote: DrugZ scores per contrast,
                        and the tabulator's combined table
data/                   the CRAVE dataset exorcise-database built from those
```

## Submitting and managing data

If enabled in `config.R`, **Submit** and **Manage** tabs are enabled and wrap
(Exorcise)[https://github.com/SimonLammmm/exorcise] analysis and database features,
starting from a counts matrix. These tabs write to datasets attached via
options in your `config.R`. Enable the two tabs and name a Docker image in
`exorcise_docker`.

- **Submit** takes a counts matrix and an experimental design workbook, runs the
  analyses the workbook asks for — MAGeCK, DrugZ, Chronos, or a manual table —
  and adds the results to a dataset. It also hands back the raw output of those
  analyses, as described below.
- **Manage** relabels a dataset's metadata, rebuilds an experiment from a
  corrected workbook, and removes experiments.

Don't turn this feature on in deployment environments where you don't trust the
users, such as over the internet. CRAVE doesn't authenticate users, so any user
can access these tabs if they're enabled.

Analyses run detached, so a long run does not block the app, and one job
at a time per dataset. Closing the tab or browser does not stop a job.

See CHANGELOG.md to see what changed in this release.
