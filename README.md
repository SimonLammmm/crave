# CRAVE
Contrasts-based Resource for Analysis, Visualisation, and Exploration

## Abstract
CRAVE enables exploration of single contrasts between a pair of conditions in a biological experiment as well as analysis of trends across contrasts.

[DDRcs](https://sjlab.cruk.cam.ac.uk/app/ddrcs/), the DNA damage response CRISPR screen viewer, runs on CRAVE and can be visited to preview the portal.

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

### Docker
Clone this repo and navigate into it.
```
git clone https://github.com/SimonLammmm/crave.git
cd crave
```
Pull the Docker image from [Docker Hub](https://hub.docker.com/repository/docker/simonlammmm/crave/tags).
```
docker pull simonlammmm/crave:latest
```
Build the Docker image.
```
cd docker
docker build -t simonlammmm/crave:latest .
```
Run CRAVE with Docker
```
docker run --rm -v "./shiny-server/config.example.R:/app/config.R" -v "./example-dataset/data/:/data/dataset-1/" -p 8080:3838 crave:latest
```
A preview of CRAVE is now accessible at http://localhost:8080.

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
Install dependencies found in the shiny-server/app.R script in the usual way for your system and then run Rshiny.
```
Rscript shiny-server/app.R
```
A preview of CRAVE is now accessible at http://localhost:8080.

## Configuration

CRAVE is configured using a config.R file that is sourced when the app runs. Use config.R to specify
 * CRAVE dataset paths
 * Exorcise data path and Docker image
 * Branding and logos

The file is written in R and therefore you need to supply R code.

When using Docker or Docker Compose, any paths need to be container paths. You also need to bind mount volumes to those container paths. When using native Rshiny, the paths can be paths on your local machine.

When using branding logos with img src, the location for those logos needs to be /app/www for Docker or Docker Compose, or shiny-server/www for native Rshiny.

## Datasets

CRAVE requires a CRAVE-formatted dataset. An example dataset is provided in this repo under example-dataset/. You can build your own dataset of CRISPR contrasts data analysed with [DrugZ](https://doi.org/10.1186/s13073-019-0665-3), [MAGeCK](https://doi.org/10.1186/s13059-014-0554-4), and [Chronos](https://doi.org/10.1186/s13059-021-02540-7) using Exorcise ([GitHub](https://github.com/SimonLammmm/exorcise)) ([Docker Hub](https://hub.docker.com/repository/docker/simonlammmm/exorcise/tags)) ([Article](https://doi.org/10.1186/s13073-024-01414-4])).

Example inputs and commands to regenerate the example dataset are in the example-dataset/ directory.
