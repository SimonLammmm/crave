# initial fresh
FROM openanalytics/r-ver:latest

# install system packages
RUN apt-get update && apt-get install --no-install-recommends -y libcurl4-openssl-dev libssl-dev libxml2-dev libgdal-dev libproj-dev sqlite3 cmake libgeos-dev libjpeg-turbo8-dev libfontconfig1-dev zlib1g-dev libharfbuzz-dev libfribidi-dev pandoc libbz2-dev apt-transport-https ca-certificates curl software-properties-common libudunits2-dev && \
# install docker
  curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg && \
  echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null && \
  apt update && \
  apt-cache policy docker-ce && \
  apt install --no-install-recommends -y docker-ce && \
# install shiny
  rm -rf /var/lib/apt.lists/* && \
  R -q -e "options(warn=2); install.packages('shiny', repos = 'https://cran.ma.imperial.ac.uk')" && \
  locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

# install R packages
RUN R -q -e "options(warn=2); install.packages(c('shinyWidgets', 'shinydashboard', 'shinythemes', 'tibble', 'tidyverse', 'ggplot2', 'visNetwork', 'heatmaply', 'leaflet', 'ggpointdensity', 'igraph', 'Hmisc', 'BiocManager', 'ggrepel', 'umap', 'mice', 'colourpicker', 'ggdensity', 'data.table', 'dplyr', 'DT', 'tidyr', 'tidyselect', 'scales', 'plotly', 'RSQLite', 'sqldf', 'shinybusy', 'ggraph', 'R.utils', 'ggdendro', 'cluster', 'foreach', 'logger', 'optparse', 'ggVennDiagram', 'DescTools', 'ggupset'), repos = 'https://cran.ma.imperial.ac.uk')" && \
  R -q -e "BiocManager::install(c('ComplexHeatmap', 'GenomicRanges', 'plyranges', 'rtracklayer', 'Biostrings'))"

# install source code
COPY crave-dev/ /src

# bind data
VOLUME /data
# At runtime, bind datasets into here, as appropriate
# At runtime, bind Exorcise dataset here, if using
# At runtime, bind "/var/run/docker.sock:/var/run/docker.sock", if using Exorcise

# run
WORKDIR /src
EXPOSE 3838
CMD ["R", "-q", "-e", "shiny::runApp('/src')"]
