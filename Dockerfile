# pull r
FROM openanalytics/r-ver:4.5.3

# install system packages
RUN apt-get update && apt-get install --no-install-recommends -y \
      libcurl4-openssl-dev libssl-dev libxml2-dev sqlite3 cmake \
      libjpeg-turbo8-dev libfontconfig1-dev zlib1g-dev libharfbuzz-dev \
      libfribidi-dev libbz2-dev libpng-dev libtiff5-dev libglpk-dev \
      pandoc apt-transport-https ca-certificates curl gnupg lsb-release \
      software-properties-common libuv1-dev && \
# install the docker CLI, needed only if you enable Exorcise
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg && \
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" > /etc/apt/sources.list.d/docker.list && \
    apt-get update && \
    apt-get install --no-install-recommends -y docker-ce-cli && \
    rm -rf /var/lib/apt/lists/* && \
    locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

ARG CRAN=https://cran.ma.imperial.ac.uk

# The three install blocks below mirror the three manifests in
# shiny-server/R/00_packages.R. Keep them in step:
#
#   CRAVE_PKGS_EAGER          attached at start-up
#   CRAVE_PKGS_LAZY           attached on first use of the tab that needs them
#   CRAVE_PKGS_INSTALLED_ONLY installed but never loaded by CRAVE itself
#
# CHANGELOG.md records what earlier versions installed and why it went.

# CRAVE_PKGS_EAGER
RUN R -q -e "options(warn = 2); install.packages(c(\
      'shiny', 'shinythemes', 'shinyWidgets', 'shinybusy', 'shinyjs', \
      'DT', 'plotly', 'visNetwork', 'colourpicker', \
      'tibble', 'dplyr', 'tidyr', 'data.table', \
      'ggplot2', 'scales', 'DBI', 'RSQLite', 'logger'), \
      repos = '${CRAN}')"

# CRAVE_PKGS_LAZY
RUN R -q -e "options(warn = 2); install.packages(c(\
      'cluster', 'ggdendro', 'igraph', 'umap', 'Rtsne', 'mice', \
      'ggVennDiagram', 'ggupset', 'DescTools'), \
      repos = '${CRAN}')"

# CRAVE_PKGS_INSTALLED_ONLY
#
# DO NOT REMOVE because nothing appears to call them. They are reached from inside
# another package, so searching the source for callers finds nothing.
#
#   R.utils - data.table::fread() passes compressed input to
#             R.utils::decompressFile(). Every CRAVE metadata file is gzipped, so
#             without R.utils on the library path fread() fails on all of them and
#             no dataset loads at all.
#
# R/00_packages.R also warns at start-up if any of these are missing.
RUN R -q -e "options(warn = 2); install.packages('R.utils', repos = '${CRAN}')"

# install source code
COPY shiny-server/ /app

# bind data
VOLUME /data
# At runtime, bind datasets into here, as appropriate
# At runtime, bind the Exorcise dataset here, if using
# At runtime, bind "/var/run/docker.sock:/var/run/docker.sock", if using Exorcise

# run
# Host and port are passed here rather than set from inside app.R. Whoever launches
# the app owns them: RStudio's Run App picks a port and expects the app to use it,
# and shiny-server and ShinyProxy assign one and pass it in. An app.R that forces
# its own binds an address the launcher is not watching, which looks like a hang.
WORKDIR /app
EXPOSE 3838
CMD ["R", "-q", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]
