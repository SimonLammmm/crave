# CRAVE
Contrasts-based Resource for Analysis, Visualisation, and Exploration

## Abstract
CRAVE enables exploration of single contrasts between a pair of conditions in a biological experiment as well as analysis of trends across contrasts.

## Deployment
CRAVE can be deployed as a Docker container or natively in Rshiny.

### Docker (recommended)
Clone this repo and navigate into it.
```
git clone https://github.com/SimonLammmm/crave.git
cd crave
```
Copy the shiny-server/config.example.R file as shiny-server/config.R
```
cp shiny-server/config.example.R shiny-server/config.R
```
Make customisation changes to the config gile. Specify where the datasets are. The dataset locations need to be relative to the /data volume that will be mounted in Docker at runtime.
```
nano shiny-server/config.R
```
Build the Docker image.
```
cd docker
docker build -t crave:latest .
```
Run the Docker container. Bind the directory containing the datasets into the /data directory in the container. Bind the Docker socket if you want to enable Exorcise.
```
docker run --rm -v path/to/folder/containing/datasets:/data -v /var/run/docker.sock:/var/run/docker.sock -p 3838:3838 crave
```
Access CRAVE in your browser on local port 3838, i.e. http://localhost:3838.

### Rshiny
Run in native Rshiny on your local machine.
Clone this repo and navigate into it.
```
git clone https://github.com/SimonLammmm/crave.git
cd crave
```
Copy the shiny-server/config.example.R file as shiny-server/config.R.
```
cp shiny-server/config.example.R shiny-server/config.r```
Make customisation changes to the config file. Specify where the datasets are.
```
nano shiny-server/config.R
```
Install dependencies found in the shiny-server/app.R script in the usual way for your system.
Run Rshiny.
```
Rscript shiny-server/app.R
```
Access CRAVE in your browser on local port 3838, i.e. http://localhost:3838.
