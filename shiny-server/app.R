
#### Versioning ####
# Author    Version       Date          Description
# SL        2.0.0         2023-08-03    Complete rework. 
# SL        2.1.0         2023-08-08    Enable lazy-loading using SQL data structure.
# SL        2.1.1         2023-08-10    Enable reference level KO in view.
# SL        2.1.2         2023-08-11    Add busy indicator. Fix edge cases with screen selection.
# SL        2.1.3         2023-11-03    Enable comparison contrast column support.
# SL        2.2.0         2023-11-03    Enable case-insensitive search.
# SL        2.2.1         2023-11-03    Change column order in comparison picker.
# SL        2.2.2         2023-11-03    Add comment metadata feature
# SL        2.2.3         2024-01-19    Bug fixes
# SL        2.3.1         2024-02-28    Improve table and plot display
# SL        3.0.0         2024-03-07    Beta release CRAVE Correlate
# SL        3.1.0         2024-03-07    Add source data switcher
# SL        3.1.1         2024-03-08    Fix race condition between the two comparison pickers
# SL        4.0.0         2024-06-07    Add enrichment plot. Enable long style (DDRcs) dataset support
# SL        4.1.0         2024-06-20    Pre-populate gene list upon screen selection rather than showing them all
# SL        4.2.0         2024-06-20    Move enrichment plot to CRAVE Correlate. Enable bulk data download
# SL        4.2.1         2024-06-26    Support external dataset
# SL        4.2.2         2024-06-27    Fix bug where contrasts involving "Formaldehyde" couldn't be queried
# SL        4.3.0         2024-06-28    Enable external and internal datasets to be viewed together
# SL        4.3.1         2024-07-01    Enable html plot downloading
# SL        4.3.2         2024-07-06    Add MOTD, logging, code optimisations
# SL        4.3.3         2024-07-09    Increase tolerance before CRAVE refuses to plot a gene query plot
# SL        4.3.4         2024-07-09    Enable ontology colouring on UMAP
# SL        4.3.5         2024-07-15    Add heatmap of gene query results
# SL        4.3.6         2024-07-16    Add dendrograms to hits heatmaps
# SL        4.3.7         2024-07-16    Fix colour gradients for asymmetric divergent scales
# SL        4.3.9         2024-08-16    Add Pendragonator
# SL        4.3.10        2024-08-19    Prepare public-suitable CRAVE
# SL        4.4           2024-08-20    Security improvements
# SL        4.4.1         2024-08-23    Add citation filter to screen filter controls
# SL        4.4.3         2024-09-16    Re-enable customisation plot for CRAVE
# SL        4.4.4         2024-09-21    Improve gene and screen filtering options, enable brushing
# SL        4.4.5         2024-11-28    Deal nicely when data do not exist
# SL        4.5           2025-01-16    Add Exorcise
# SL        4.6           2025-01-26    Add support for Chronos
# SL        4.6.1         2025-02-07    Improve stability of Exorcise input form
# SL        4.7           2025-05-17    Add overlap analysis, minor tweaks to comparison picker
# SL        4.8           2025-05-20    Add comparison search by numerator and denominator levels
# SL        4.8.1         2025-05-29    Improve comparison search, add core essentialome shortcut
# SL        4.8.2         2025-05-30    Fix Correlate screen search for cell line
# SL        4.9           2025-06-03    Enable Exorcise upload CRISPick output, add guide library, remove the need for a dummy dataset, make switching between DDRcs/CRAVE easier
# SL        4.9.1         2025-06-05    Fix when Exorcise advanced chemistry string passed even when empty, move to using Exorcise Docker image
# SL        4.10          2025-06-06    Move Guide Library to SQL
# SL        4.10.1        2025-06-06    Clean up files after Exorcise to avoid collisions with other users
# SL        4.10.2        2025-06-09    Deal properly when internal dataset absent (i.e. in the case of DDRcs)
# SL        4.10.4        2025-06-11    Implement externally bindable message banners on the front page
# SL        4.10.5        2025-06-11    Enable HTML in message banners
# SL        4.10.6        2025-06-11    Enable logging
# SL        4.11          2025-06-17    Adopt native Plotly for improved speed
# SL        4.11.2        2025-06-18    Make web links open in a new tab, improve Kinds
# SL        4.12          2025-07-01    Add libraries table, fix when 1 external experiment is selected and "clear form" is chosen, fix FGC accessions
# SL        4.12.2        2025-07-04    Add hypergeometric adjusted p-value to Enrichment and Overlap; fix bugs in Biplot
# SL        4.13          2025-07-07    Improve contrast picker in Correlate
# SL        4.14          2025-07-09    Add Explore ROC-AUC
# SL        4.14.1        2025-07-09    UI/UX improvements
# SL        4.14.4        2025-07-30    UI/UX improvements
# SL        4.15          2025-07-31    Save/load feature
# SL        4.16          2025-10-09    Enable score and rank cutoff in Correlate
# SL        4.17          2025-10-10    Implement upset plot
# SL        4.18          2025-10-24    Support "Manual" externally-analysed screens
# SL        4.19          2025-10-27    Enable Pendragonator for GOIs
# SL        4.19.1        2025-11-06    Add support for Endpoint selection for "Manual" method comparisons
# SL        4.19.2        2025-11-18    Add boxplot to gene query
# SL        4.19.3        2026-01-27    Enable enrichment for DDRcs
# SL        4.20          2026-03-05    Add ontology view
# SL        4.21          2026-03-05    Allow disabling dendrograms on the clustrgram
# AE        4.22          2026-06-23    Arbitrary dataset support via config.R; remove guides/libraries/ontology requirement
# AE        4.22          2026-06-23    Custom branding (portal name, front-page text, footer) via config.R
# SL        4.23          2026-06-23    Intelligently enable/disable features that depend on ontology, library, guide library, and Exorcise available
# SL        4.24          2026-06-23    Deprecate genetypes, never implemented
# SL        4.25          2026-07-06    Basic config and dataset error checking
# SL        4.26          2026-07-08    Fix Exorcise data and work directories with config.R
# SL        4.27          2026-07-22    Add support for Exorcise 2.0.1
# SL        4.27.1        2026-07-23    Speed up check for specified Exorcise installation

ver <<- "4.27.1"
updated <<- "2026-07-23"

#### Preamble ####

Sys.setenv(
  `_R_CHECK_LENGTH_1_CONDITION_` = F,
  `_R_CHECK_LENGTH_1_LOGIC2_` = F
)
suppressPackageStartupMessages({
  library(tibble)
  library(tidyverse)
  library(shiny)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(DT)
  library(tidyr)
  library(tidyselect)
  library(scales)
  library(plotly)
  library(shinythemes)
  library(RSQLite)
  library(sqldf)
  library(shinybusy)
  library(shinyWidgets)
  library(igraph)
  library(visNetwork)
  library(ggraph)
  library(R.utils)
  library(ggdendro)
  library(cluster)
  library(umap)
  library(mice)
  library(colourpicker)
  library(ggdensity)
  library(foreach)
  library(logger)
  library(rmarkdown)
  library(ggVennDiagram)
  library(DescTools)
  library(ggupset)
})
options(warn = 1)
noWS <<- c("before", "after", "outside", "after-begin", "before-end")

# Check that config.R exists
if(!file.exists("config.R")) {
  # If not, then serve an error page
  ui <- fluidPage("Error starting CRAVE: config.R not found. Please obtain a template config.example.R file at <https://github.com/SimonLammmm/crave/blob/a95e5a6034baee7837d999cfbba46adb0a35f48c/shiny-server/config.example.R>, configure it, and put it next to app.R. If you're using Docker Compose, edit docker-compose.yml and bind it at /app/config.R.")
  server <- function(input, output, session) {}
  options(shiny.host = "0.0.0.0", shiny.port = 3838)
  shinyApp(ui, server)
} else {
  
  # If config.R exists, load modules
  source("config.R")
  source("biplot.R")
  source("common.R")
  source("dataset.R")
  source("download.R")
  source("enrichment.R")
  source("exorcise.R")
  source("genequery.R")
  source("guideLibrary.R")
  source("hitmap.R")
  source("heatmap.R")
  source("init.R")
  source("legal.R")
  source("network.R")
  source("overlap.R")
  source("pendragonator.R")
  source("rank.R")
  source("rocauc.R")
  source("umap.R")
  source("upset.R")
  source("volcano.R")
  
  # Load and validate datasets
  load()
  if(!canExplore) {
    # If invalid, then serve an error page
    ui <- fluidPage("Error starting CRAVE: invalid datasets. Did you forget to specify a CRAVE dataset in config.R? If using Docker Compose, you need to bind config.R at /src/config.R and make sure the paths in config.R point to absolute paths in the container.")
    server <- function(input, output, session) {}
    options(shiny.host = "0.0.0.0", shiny.port = 3838)
    shinyApp(ui, server)
    
  } else {
    #### Run ####
    source("uielements.R")
    source("shinyui.R")
    source("shinyserver.R")
    options(shiny.host = "0.0.0.0", shiny.port = 3838, shiny.maxRequestSize = 100 * 1024 * 1024)
    shinyApp(ui, server)
  }
}
