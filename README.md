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

cp shiny-server/config.example.R shiny-server/config.r
```
Install dependencies found in the shiny-server/app.R script in the usual way for your system and then run Rshiny.
```
Rscript shiny-server/app.R
```
A preview of CRAVE is now accessible at http://localhost:8080.

Customise CRAVE by editing the shiny-server/config.R file.
