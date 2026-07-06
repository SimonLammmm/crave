#### Front page notice banners ####
# Info banner
if(!is.null(notices_info)) {
  # Draw info banner if there is info to show
  noticebanner_info <- tags$div(
    style = "background-color: #e6f0fc; padding-top: 12px; padding-right: 12px; padding-bottom: 12px; padding-left: 12px; border-style: none none none solid; border-color: #0469e3; border-width: 6px",
    do.call(tags$p,
            list(style = "color: black",
                 icon("circle-info"),
                 notices_info))
  )
} else {
  noticebanner_info <- NULL
}
# Warning banner
if(!is.null(notices_warning)) {
  # Draw warning banner if there is a warning to show
  noticebanner_warning <- tags$div(
    style = "background-color: #fff9e6; padding-top: 12px; padding-right: 12px; padding-bottom: 12px; padding-left: 12px; border-style: none none none solid; border-color: #fec000; border-width: 6px",
    do.call(tags$p,
            list(style = "color: black",
                 icon("triangle-exclamation"),
                 notices_warning))
  )
} else {
  noticebanner_warning <- NULL
}
# Emergency banner
if(!is.null(notices_emergency)) {
  # Draw emergency banner if there is an emergency to show
  noticebanner_emergency <- tags$div(
    style = "background-color: #FF3127; padding-top: 12px; padding-right: 12px; padding-bottom: 12px; padding-left: 12px; border-style: none none none solid; border-color: #FF3127; border-width: 6px",
    do.call(tags$p,
            list(style = "color: white",
                 icon("circle-exclamation"),
                 notices_emergency))
  )
} else {
  noticebanner_emergency <- NULL
}

#### Shiny UI: correlate sidebar ####
correlateSidebar <-
  sidebarPanel(
    width = 3,
    tags$h4(tags$strong("Gene filter controls")),
    textAreaInput('genes_queried_text', label = "Type genes of interest", placeholder = "e.g. TP53, BRCA1, UBE2K", rows = 3),
    actionLink('genes_queried_text_validate', label = "Validate"),
    " · ",
    actionLink('genes_queried_essentials', label = "Core essentials"),
    " · ",
    actionLink('genes_copy_from_explore', label = "Copy GOIs from Explore"),
    tags$br(),
    tags$br(),
    selectizeInput('genes_queried', label = "Search genes of interest", choices = NULL, multiple = TRUE, options = list(placeholder = "e.g. TP53")),
    actionButton('gq_clear_genes', label = "Clear genes", icon = icon("ban"), class = "btn-danger"),
    tags$br(),
    tags$br(),
    selectizeInput('gq_stat', label = "Statistic", choices = c("Score", "FDR", "Rank"), selected = "Score"),
    numericInput('gq_cutoff', label = "Cutoff", value = defaults_gq_cutoff_selected),
    selectizeInput('analysis_method_gq', label = "Method", choices = c("DrugZ", "MAGeCK", "Chronos", "Manual")),
    tags$hr(),
    tags$h4(tags$strong("Screen filter controls")),
    tags$p(tags$strong("Experiment")),
    tags$table(style = "width: 100%",
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_citation', label = "Citation", choices = NULL, multiple = TRUE)),
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_kind', label = "Kind", choices = NULL, multiple = TRUE))
               ),
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_library', label = "Library", choices = NULL, multiple = TRUE)),
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_timepoint', label = "Timepoint", choices = NULL, multiple = TRUE))
               ),
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_organism', label = "Species", choices = NULL, multiple = TRUE)),
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_source', label = "Source", choices = NULL, multiple = TRUE))
               )
    ),
    tags$table(style = "width: 100%",
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px; width: 50%",
                               tags$p(tags$strong("Numerator levels"))),
                       tags$td(style = "padding: 6px; width: 50%",
                               tags$p(tags$strong("Denominator levels")))
               ),
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_days_grown_diff', label = "Days grown (diff)", choices = NULL, multiple = TRUE)),
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_days_grown_ref', label = "Days grown (ref)", choices = NULL, multiple = TRUE))
               ),
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_treatment_diff', label = "Treatment (diff)", choices = NULL, multiple = TRUE)),
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_treatment_ref', label = "Treatment (ref)", choices = NULL, multiple = TRUE))
               ),
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_dose_diff', label = "Dose (diff)", choices = NULL, multiple = TRUE)),
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_dose_ref', label = "Dose (ref)", choices = NULL, multiple = TRUE))
               ),
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_knockout_diff', label = "Knockout (diff)", choices = NULL, multiple = TRUE)),
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_knockout_ref', label = "Knockout (ref)", choices = NULL, multiple = TRUE))
               ),
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_cellline_diff', label = "Cell line (diff)", choices = NULL, multiple = TRUE)),
                       tags$td(style = "padding: 6px; width: 50%",
                               selectizeInput('gq_filter_cellline_ref', label = "Cell line (ref)", choices = NULL, multiple = TRUE))
               )
    ),
    tags$p(tags$strong("Custom")),
    selectizeInput('gq_filter_contrast', label = "Filter by contrast", choices = NULL, multiple = TRUE),
    actionLink('gq_filter_contrast_copy_from_explore', label = "Copy contrasts from Explore"),
    tags$br(),
    tags$br(),
    textInput('gq_filter_custom', label = "Filter by text input"),
    textOutput('gq_filtered_screens'),
    tags$br(),
    tags$br(),
    actionButton('gq_example', "Example", icon = icon("fire")),
    actionButton('gq_reset', "Reset form", icon = icon("arrows-rotate"), class = "btn-danger"),
    tags$br(),
    tags$br(),
    actionButton('correlate_customise', "Customise plot", icon = icon("brush")),
    tags$hr(),
    tags$p(tags$strong("Help")),
    tags$p("Use", tags$strong("gene filters"), "to choose genes of interest. Leave it blank to consider all genes. Use the FDR cutoff slider to constrain for hit strength."),
    tags$p("Use", tags$strong("screen filters"), "to choose screens of interest. Screens satisfying at least one filter in each field are considered. Screens of interest are previewed below the fields."),
    tags$p("Once you've chosen your filters, click a tab to run an analysis:"),
    tags$p(tags$strong("Gene query"), "and", tags$strong("Clustergram"), " show CRISPR scores across screens."),
    tags$p(tags$strong("Heatmap"), "and", tags$strong("Network"), " show pairwise correlation between genes"),
    tags$p(tags$strong("UMAP"), "uses the ", tags$a(href="https://doi.org/10.48550/arXiv.1802.03426", "uniform manifold approximation and projection", target="_blank", .noWS = noWS), " method to reveal clusters of similarly behaving genes"),
    tags$p(tags$strong("Enrichment"), "uses the hypergeometric test to reveal enrichment and depletion of ontological terms among hits, using the ontology described by ", tags$a(href = "https://doi.org/10.1016/j.cell.2020.05.040", "Olivieri et al (2020)", target="_blank", .noWS = noWS), "."),
    tags$p(tags$strong("Pendragonator"), "returns the top responding genes among your screens of interest and the conditions that generated them."),
    tags$p(tags$strong("Data download"), "lets you bulk download CRISPR scores.")
  )

#### Shiny UI: gene query ####
genequeryTab <-
  tabPanel("Gene query", icon = icon("magnifying-glass-chart"),
           tags$br(),
           tags$h4("Gene query"),
           tags$p("Generate violin plots of CRISPR scores and shows where your GOIs were hits among your screens of interest."),
           tags$p("Leave the genes of interest field empty to query the top 250 genes."),
           tags$p("GOIs that were hits with FDR better than the cutoff are shown as triangles. All other GOIs are shown as circles. Only screens with at least one triangle are shown."),
           tags$p("Double-click in the legend to isolate a plot element. Single-click to toggle plot elements on and off."),
           tags$br(),
           selectizeInput("genequery_violin", "Plot type", choices = c("Violin plot", "Boxplot"), selected = "Violin plot", multiple = F),
           tags$br(),
           actionButton('genequery_submit', "Gene query", icon = icon("magnifying-glass-chart"), class = "btn-success"),
           tags$br(),
           plotlyOutput('genequery_plot', height = NULL),
           DTOutput('genequery_table'),
           tags$br(),
           downloadButton('genequery_download_plot', "Download plot", icon = icon("download")),
           downloadButton('genequery_download', "Download data", icon = icon("download"))
  )

#### Shiny UI: clustergram ####
clustergramTab <-
  tabPanel("Clustergram", icon = icon("diagram-project"),
           tags$br(),
           tags$h4("Clustergram"),
           tags$p("Generate a heatmap of CRISPR scores for your GOIs among your screens of interest."),
           tags$p("Leave the genes of interest field empty to query the top 250 genes."),
           tags$p("Genes and screens are arranged into clusters of similar behaviour. Dendrograms are shown when there are at least three GOIs or at least three screens of interest. Clustering uses the ", tags$a(href="https://doi.org/10.1002/9780470316801", "agglomerative nesting (Agnes)", target="_blank", .noWS = noWS), " method."),
           tags$p("Bigger scores have greater colour intensity. Hits exceeding the FDR cutoff are shown with an asterisk."),
           tags$br(),
           actionButton('hitmap_submit', "Clustergram", icon = icon("diagram-project"), class = "btn-success"),
           checkboxInput('hitmap_gene_dendrogram', "Draw dendrogram on genes", T),
           checkboxInput('hitmap_screen_dendrogram', "Draw dendrogram on screens", T),
           tags$br(),
           plotlyOutput('hitmap_plot', height = NULL),
           tags$br(),
           DTOutput('hitmap_table'),
           tags$br(),
           downloadButton('hitmap_download_plot', "Download plot", icon = icon("download")),
           downloadButton('hitmap_download', "Download data", icon = icon("download"))
  )

#### Shiny UI: heatmap ####
heatmapTab <-
  tabPanel("Heatmap", icon = icon("temperature-half"),
           tags$br(),
           tags$h4("Heatmap"),
           tags$p("Generate a gene-gene correlation heatmap using evidence from your screens of interest."),
           tags$p("Genes are arranged into clusters of similar correlation coefficient. A dendrogram is shown when there are at least three GOIs. Clustering uses the ", tags$a(href="https://doi.org/10.1002/9780470316801", "agglomerative nesting (Agnes)", target="_blank", .noWS = noWS), " method."),
           tags$p("Stronger pairwise correlation coefficients have greater colour intensity. Correlation coefficients are computed using the Spearman method."),
           tags$br(),
           actionButton('heatmap_submit', "Heatmap", icon = icon("temperature-half"), class = "btn-success"),
           tags$br(),
           plotlyOutput('heatmap_plot', height = 900),
           tags$br(),
           DTOutput('heatmap_table'),
           tags$br(),
           downloadButton('heatmap_download_plot', "Download plot", icon = icon("download")),
           downloadButton('heatmap_download', "Download data", icon = icon("download"))
  )

#### Shiny UI: UMAP ####
umapTab <-
  tabPanel("UMAP", icon = icon("map"),
           tags$br(),
           tags$h4("UMAP"),
           tags$p("Generate a scatter plot of GOIs using the ", tags$a(href="https://doi.org/10.48550/arXiv.1802.03426", "uniform manifold approximation and projection (UMAP)", target="_blank", .noWS = noWS), " method"),
           tags$p("Missing values are dealt with by imputation."),
           tags$br(),
           actionButton('umap_submit', "UMAP", icon = icon("map"), class = "btn-success"),
           tags$br(),
           plotlyOutput('umap_plot', height = 900),
           tags$br(),
           DTOutput('umap_table'),
           tags$br(),
           downloadButton('umap_download_plot', "Download plot", icon = icon("download")),
           downloadButton('umap_download', "Download data", icon = icon("download"))
  )

#### Shiny UI: correlation network ####
networkplotTab <-
  tabPanel("Network", icon = icon("circle-nodes"),
           tags$br(),
           tags$h4("Network"),
           tags$p("Draw a network representation of gene-gene correlations using evidence from your screens of interest."),
           tags$p("Only pairwise correlations at least as strong as the cutoff are shown. Correlation coefficients are computed using the Spearman method."),
           tags$br(),
           sliderInput('corr_cutoff', label = "Correlation coefficient cutoff", min = defaults_corr_cutoff_min, max = defaults_corr_cutoff_max, step = defaults_corr_cutoff_step, value = defaults_corr_cutoff_value),
           actionButton('network_submit', "Network", icon = icon("circle-nodes"), class = "btn-success"),
           tags$br(),
           visNetworkOutput('network_plot', height = 900),
           tags$br(),
           DTOutput('network_table'),
           tags$br(),
           downloadButton('network_download_plot', "Download plot", icon = icon("download")),
           downloadButton('network_download', "Download data", icon = icon("download"))
  )

#### Shiny UI: enrichment analysis ####
if(canOntology) {
  enrichmentanalysisTab <-
    tabPanel("Enrichment", icon = icon("blender"),
             tags$br(),
             tags$h4("Enrichment"),
             tags$p("Generate a heatmap showing enrichment or depletion of ontological terms  among the hits within your screens of interest."),
             tags$p("Use the controls to define the cutoff for a hit. Use the controls to choose a direction and test tail."),
             tags$p("Significantly enriched or depleted terms are shown with an asterisk."),
             tags$p("Cells are coloured by the number of hits as a proportion of all genes in that ontological term (recall). More intense colour means greater recall."),
             tags$br(),
             tags$table(
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px", selectizeInput('enrichment_cutoff_stat', label = "Cutoff statistic", choices = c("Score", "FDR", "Rank"))),
                       tags$td(style = "padding: 6px", textInput('enrichment_cutoff_value', label = "Better than", value = 3))
               ),
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px", selectizeInput('enrichment_cutoff_direction', label = "Direction", choices = c("Hypersensitivity hits", "Suppressing hits"))),
                       tags$td(style = "padding: 6px", selectizeInput('enrichment_cutoff_tail', label = "Test tail", choices = c("Upper tail (enriched classes)", "Lower tail (depleted classes)")))
               ),
               tags$tr(style = "vertical-align: top",
                       tags$td(style = "padding: 6px", actionButton('enrichment_submit', 'Enrichment analysis', icon = icon("blender"), class = "btn-success"))
               )
             ),
             tags$br(),
             plotlyOutput('enrichment_plot', height = 900),
             tags$br(),
             DTOutput('enrichment_table'),
             tags$br(),
             downloadButton('enrichment_download_plot', "Download plot", icon = icon("download")),
             downloadButton('enrichment_download', "Download data")
    )
} else {
  enrichmentanalysisTab <- NULL
}

#### Shiny UI: DDRome ####
ddromeTab <-
  tabPanel("Pendragonator", icon = icon("dragon"),
           tags$br(),
           tags$h4("Pendragonator"),
           tags$p("Generate a list of genes in your screens of interest using Pendragon library design principles."),
           tags$p("Returns the strongest responding genes on either tail and the conditions that produced the hit."),
           tags$p("Tick the box to find the screens that produced the best hit among only your genes of interest."),
           tags$br(),
           tags$table(
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", 
                             numericInput('pendragonator_n', label = "List length", value = 2000, min = 1, max = 20000, step = 1)
                     ),
                     tags$td(style = "padding: 6px", 
                             selectizeInput('pendragonator_stat', label = "Statistic", choices = c("Score", "FDR"))
                     ),
                     tags$td(style = "padding: 6px", 
                             checkboxInput('pendragonator_goi', label = "Pendragonate only genes of interest", value = F)
                     )
             )
           ),
           tags$br(),
           actionButton('pendragonator_submit', "Pendragonate", icon = icon("dragon"), class = "btn-success"),
           tags$br(),
           tags$br(),
           DTOutput('pendragonator_table'),
           tags$br(),
           downloadButton('pendragonator_download', "Download", icon = icon("download"))
  )

#### Shiny UI: data download ####
datadownloadTab <-
  tabPanel("Data download", icon = icon("download"),
           tags$br(),
           tags$h4("Data download"),
           tags$p("Bulk download CRISPR data."),
           tags$p("Specify the file format by appending \".csv\" (comma-separated values) or \".tsv\" (tab-separated values). Append \".gz\" for gzipped output."),
           tags$p("Click \"Fetch\" to preview the output and then \"Download data\" to start the download."),
           tags$br(),
           textInput('download_filename', label = "Filename", value = "export.tsv.gz"),
           tags$br(),
           actionButton('download_submit', "Fetch", icon = icon("download"), class = "btn-success"),
           tags$br(),
           textOutput('download_console'),
           tags$br(),
           DTOutput('download_table'),
           tags$br(),
           downloadButton('download_download', "Download data", icon = icon("download"))
  )

### Shiny UI: explore sidebar ####
exploreSidebar <-
  sidebarPanel(
    width = 3,
    textAreaInput('genes_selected_text', label = "Type genes of interest", rows = 3),
    actionLink('genes_selected_text_validate', label = "Validate"),
    " · ",
    actionLink('genes_select_essentialome', label = "Core essentials"),
    " · ",
    actionLink('genes_copy_from_correlate', label = "Copy GOIs from Correlate"),
    tags$br(),
    tags$br(),
    selectizeInput('genes_selected', label = "Search genes of interest", choices = NULL, multiple = TRUE),
    actionButton('genes_clear', label = "Clear genes", icon = icon("ban"), class = "btn-danger"),
    tags$br(),
    tags$br(),
    selectizeInput('comparisons_selected', label = "Selected comparisons", choices = NULL, multiple = TRUE),
    actionButton('comparisons_table_selectNone', "Clear comparisons", icon = icon("ban"), class = "btn-danger"),
    tags$br(),
    tags$br(),
    actionButton('explore_customise', label = "Customise plot", icon = icon("brush")),
    tags$br(),
    tags$br(),
    tags$hr(),
    tags$p(tags$strong("Help")),
    tags$p("Use", tags$strong("Select experiment"), "and", tags$strong("Select comparison"), "tabs to choose your screen(s) of interest."),
    tags$p(tags$strong("Volcano plot"), "and", tags$strong("Rank plot"), "show the distribution of CRISPR scores in a single screen."),
    tags$p(tags$strong("Biplot"), " compares the distributions of CRISPR scores of two screens."),
    tags$p(tags$strong("Overlap"), " shows common hits between two or more screens and calculates the statistical significance of the overlap between pairs of screens using a hypergeometric test."),
    tags$p(tags$strong("ROC"), " shows how well hits from a reference screen are captured by other screens using precision-recall curves.")
  )

#### Shiny UI: experiment selector ####
experimentSelectTab <-
  tabPanel("Select experiment", icon = icon("flask"),
           tags$br(),
           tags$h4("Select experiment"),
           tags$p("Use the table to choose papers or NGS submissions of interest. The", tags$strong("Select comparison"), "tab will be populated with screens from only the selected experiments."),
           tags$p("Leave blank to show screens from the whole database."),
           tags$br(),
           tags$hr(),
           tags$br(),
           actionButton('experiments_table_selectAll', "Select all", icon = icon("check-double"), class = "btn-success"),
           actionButton('experiments_table_selectNone', "Select none", icon = icon("ban"), class = "btn-danger"),
           tags$br(),
           tags$br(),
           DTOutput('experiments_table')
  )

#### Shiny UI: comparison selector ####
comparisonSelectTab <-
  tabPanel("Select comparison", icon = icon("vials"),
           tags$br(),
           tags$h4("Select comparison"),
           tags$p("Use the form to filter screens by experimental design or conditions on the denominator (reference) or numerator (differential)."),
           tags$p("Use the table to choose screens of interest. Click the \"Select comparisons\" button to choose all of them."),
           tags$p("As you fill the form, the contents of the table updates."),
           tags$br(),
           tags$hr(),
           tags$br(),
           tags$p(tags$strong("Comparison search")),
           tags$p("Selected experiments:"),
           textOutput('comparison_selector_experiments'),
           tags$br(),
           tags$p(tags$strong("Experiment")),
           tags$table(
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_citation', label = "Citation", choices = NULL, multiple = T),
                     ),
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_kind', label = "Kind", choices = NULL, multiple = T),
                     )
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_library', label = "Library", choices = NULL, multiple = T)
                     ),
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_timepoint', label = "Timepoint", choices = NULL, multiple = T)
                     )
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_organism', label = "Species", choices = NULL, multiple = T)
                     ),
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_source', label = "Source", choices = NULL, multiple = T)
                     )
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_endpoint', label = "Endpoint", choices = NULL, multiple = T)
                     )
             )
           ),
           hr(),
           tags$table(
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             tags$p(tags$strong("Numerator levels"))),
                     tags$td(style = "padding: 6px",
                             tags$p(tags$strong("Denominator levels")))
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_daysgrown_diff', label = "Days grown (diff)", choices = NULL, multiple = T)),
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_daysgrown_ref', label = "Days grown (ref)", choices = NULL, multiple = T))
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_treatment_diff', label = "Treatment (diff)", choices = NULL, multiple = T)),
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_treatment_ref', label = "Treatment (ref)", choices = NULL, multiple = T))
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_dose_diff', label = "Dose (diff)", choices = NULL, multiple = T)),
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_dose_ref', label = "Dose (ref)", choices = NULL, multiple = T))
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_knockout_diff', label = "Knockout (diff)", choices = NULL, multiple = T)),
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_knockout_ref', label = "Knockout (ref)", choices = NULL, multiple = T))
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_cell_line_diff', label = "Cell line (diff)", choices = NULL, multiple = T)),
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_selector_cell_line_ref', label = "Cell line (ref)", choices = NULL, multiple = T))
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             actionButton('comparisons_selector_select', "Select comparisons", icon = icon("circle-check"), class = "btn-success"),
                             " · ",
                             actionLink('comparisons_selector_clear', "Clear form")),
                     tags$td(style = "padding: 6px",
                             textOutput('comparisons_selector_list')))
           ),
           tags$br(),
           tags$hr(),
           tags$br(),
           tags$p(tags$strong("All comparisons")),
           tags$br(),
           tags$br(),
           DTOutput('comparisons_table')
  )

#### Shiny UI: Screen explorer: Volcano plot ####
volcanoplotTab <-
  tabPanel("Volcano plot", icon = icon("volcano"),
           tags$br(),
           tags$h4("Volcano plot"),
           tags$p("Generate a volcano plot of CRISPR score against FDR for one screen."),
           tags$p("Choose between ", tags$a(href="https://doi.org/10.1186/s13073-019-0665-3", "DrugZ", target="_blank", .noWS = noWS), " and ", tags$a(href="https://doi.org/10.1186/s13059-014-0554-4", "MAGeCK", target="_blank", .noWS = noWS), " methods."),
           tags$br(),
           tags$table(
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", selectizeInput('comparison_volcano', label = "Select a screen", choices = NULL, options = list(placeholder = "Choose a screen"))),
                     tags$td(style = "padding: 6px", selectizeInput('analysis_method_se_v', label = "Method", choices = c("DrugZ", "MAGeCK", "Manual")))
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", actionButton('volcano_submit', 'Volcano plot', icon = icon("volcano"), class = "btn-success"))
             )
           ),
           tags$br(),
           uiOutput('volcano_plot'),
           tags$br(),
           DTOutput('volcano_table'),
           downloadButton('volcano_download_plot', "Download plot", icon = icon("download")),
           downloadButton('volcano_download', "Download data", icon = icon("download"))
  )

#### Shiny UI: Screen explorer: Rank plot ####
rankplotTab <-
  tabPanel("Rank plot", icon = icon("ranking-star"),
           tags$br(),
           tags$h4("Rank plot"),
           tags$p("Generate a rank plot of CRISPR score against FDR for one screen."),
           tags$p("Choose between ", tags$a(href="https://doi.org/10.1186/s13073-019-0665-3", "DrugZ", target="_blank", .noWS = noWS), ", ", tags$a(href="https://doi.org/10.1186/s13059-014-0554-4", "MAGeCK", target="_blank", .noWS = noWS), ", and ", tags$a(href="https://doi.org/10.1186/s13059-021-02540-7", "Chronos", target="_blank", .noWS = noWS), " methods."),
           tags$br(),
           tags$table(
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", selectizeInput('comparison_rank', label = "Select a screen", choices = NULL, options = list(placeholder = "Choose a screen"))),
                     tags$td(style = "padding: 6px", selectizeInput('analysis_method_se_r', label = "Method", choices = c("DrugZ", "MAGeCK", "Chronos", "Manual")))
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", actionButton('rank_submit', 'Rank plot', icon = icon("ranking-star"), class = "btn-success"))
             )
           ),
           tags$br(),
           uiOutput('rank_plot'),
           tags$br(),
           DTOutput('rank_table'),
           downloadButton('rank_download_plot', "Download plot", icon = icon("download")),
           downloadButton('rank_download', "Download data", icon = icon("download"))
  )

#### Shiny UI: Screen explorer: Biplot ####
biplotTab <-
  tabPanel("Biplot", icon = icon("compass-drafting"),
           tags$br(),
           tags$h4("Biplot"),
           tags$p("Plot CRISPR scores for two screens against each other."),
           tags$p("Choose between ", tags$a(href="https://doi.org/10.1186/s13073-019-0665-3", "DrugZ", target="_blank", .noWS = noWS), ", ", tags$a(href="https://doi.org/10.1186/s13059-014-0554-4", "MAGeCK", target="_blank", .noWS = noWS), ", and ", tags$a(href="https://doi.org/10.1186/s13059-021-02540-7", "Chronos", target="_blank", .noWS = noWS), " methods. You can choose a different method for each screen."),
           tags$br(),
           tags$table(
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px",
                             selectizeInput('comparison_X', label = "Select two screens", choices = NULL, options = list(placeholder = "Choose a screen to go on the x-axis")),
                             selectizeInput('comparison_Y', label = NULL, choices = NULL, options = list(placeholder = "Choose a screen to go on the y-axis"))),
                     tags$td(style = "padding: 6px",
                             selectizeInput('analysis_method_se_b_x', label = "Methods", choices = c("DrugZ", "MAGeCK", "Chronos", "Manual")),
                             selectizeInput('analysis_method_se_b_y', label = NULL, choices = c("DrugZ", "MAGeCK", "Chronos", "Manual")))
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", actionButton('bi_submit', 'Biplot', icon = icon("compass-drafting"), class = "btn-success"))
             )
           ),
           tags$br(),
           uiOutput('bi_plot'),
           tags$br(),
           DTOutput('bi_table'),
           downloadButton('bi_download_plot', "Download plot", icon = icon("download")),
           downloadButton('bi_download', "Download data", icon = icon("download"))
  )

#### Shiny UI: overlap analysis ####
overlapTab <-
  tabPanel("Overlap", icon = icon("circle-half-stroke"),
           tags$br(),
           tags$h4("Overlap"),
           tags$p("Draw a Venn diagram showing common hits between two or more screens."),
           tags$p("Use the controls to set the threshold for a hit."),
           tags$p("Choose between ", tags$a(href="https://doi.org/10.1186/s13073-019-0665-3", "DrugZ", target="_blank", .noWS = noWS), ", ", tags$a(href="https://doi.org/10.1186/s13059-014-0554-4", "MAGeCK", target="_blank", .noWS = noWS), ", and ", tags$a(href="https://doi.org/10.1186/s13059-021-02540-7", "Chronos", target="_blank", .noWS = noWS), " methods."),
           tags$br(),
           tags$table(
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", selectizeInput('comparison_overlap', label = "Select screens", choices = NULL, multiple = T)),
                     tags$td(style = "padding: 6px", selectizeInput('overlap_cutoff_method', label = "Method", choices = c("DrugZ", "MAGeCK", "Chronos", "Manual"))))),
           tags$table(
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", selectizeInput('overlap_cutoff_stat', label = "Statistic", choices = c("Score", "FDR", "Rank"), width = "100%")),
                     tags$td(style = "padding: 6px", textInput('overlap_cutoff_value', label = "Better than", value = 3, width = "100%")),
                     tags$td(style = "padding: 6px", selectizeInput('overlap_cutoff_direction', label = "Direction", choices = c("Hypersensitivity hits", "Suppressing hits"), width = "100%"))),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", selectizeInput('overlap_style', label = "Type of plot", choices = c("Venn diagram", "Upset plot"), width = "100%"))
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", actionButton('overlap_submit', 'Overlap', icon = icon("circle-half-stroke"), class = "btn-success"))
             )
           ),
           tags$br(),
           uiOutput('overlap_plot'),
           tags$br(),
           DTOutput('overlap_table'),
           tags$br(),
           downloadButton('overlap_download_plot', "Download plot", icon = icon("download")),
           downloadButton('overlap_download', "Download data"))

#### Shiny UI: ROC-AUC analysis ####
rocaucTab <-
  tabPanel("ROC", icon = icon("chart-area"),
           tags$br(),
           tags$h4("ROC"),
           tags$p("Draw a ROC curve showing the precision-recall of hits in a reference screen by query screens."),
           tags$p("Use the controls to set the threshold for a hit in the reference screen."),
           tags$p("Choose between ", tags$a(href="https://doi.org/10.1186/s13073-019-0665-3", "DrugZ", target="_blank", .noWS = noWS), ", ", tags$a(href="https://doi.org/10.1186/s13059-014-0554-4", "MAGeCK", target="_blank", .noWS = noWS), ", and ", tags$a(href="https://doi.org/10.1186/s13059-021-02540-7", "Chronos", target="_blank", .noWS = noWS), " methods."),
           tags$br(),
           tags$table(
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", selectizeInput('rocauc_base', label = "Select reference screen", choices = NULL)),
                     tags$td(style = "padding: 6px", selectizeInput('rocauc_method', label = "Method", choices = c("DrugZ", "MAGeCK", "Chronos", "Manual"))))),
           tags$table(
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", selectizeInput('rocauc_cutoff_stat', label = "Statistic", choices = c("Score", "FDR", "Rank"), width = "100%")),
                     tags$td(style = "padding: 6px", textInput('rocauc_cutoff_value', label = "Better than", value = 3, width = "100%")),
                     tags$td(style = "padding: 6px", selectizeInput('rocauc_cutoff_direction', label = "Direction", choices = c("Hypersensitivity hits", "Suppressing hits", "Both"), width = "100%"))),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", selectizeInput('rocauc_head', label = "Select query screens", choices = NULL, multiple = T))
             ),
             tags$tr(style = "vertical-align: top",
                     tags$td(style = "padding: 6px", actionButton('rocauc_submit', 'ROC', icon = icon("chart-area"), class = "btn-success"))
             )
           ),
           tags$br(),
           uiOutput('rocauc_plot'),
           tags$br(),
           DTOutput('rocauc_table'),
           tags$br(),
           downloadButton('rocauc_download_plot', "Download plot", icon = icon("download")),
           downloadButton('rocauc_download', "Download data"))

#### Shiny UI: Exorcise ####
if(canExorcise) {
  exorciseTab <-
    tabPanel("Exorcise", icon = icon("ghost"),
             sidebarLayout(
               sidebarPanel(
                 width = 4,
                 tags$h4(tags$strong("Sequences")),
                 fileInput('exorcise_upload', "Upload a CRISPick sgRNA designs file", multiple = F, accept = "text/plain"),
                 tags$p("or"),
                 textAreaInput('exorcise_seq', label = "Enter sequence(s), separated by commas or spaces", placeholder = "e.g. ATCGATCGATCGATCGATCG, GTACGTACGTACGTACGTAC", rows = 4),
                 textAreaInput('exorcise_orig', label = "Optional: existing annotation(s), one per sequence, separated by commas or spaces", placeholder = "e.g. TP53, negative_control", rows = 4),
                 textInput('exorcise_control_string', label = "Optional: existing annotation(s) that indicate controls", placeholder = "e.g. negative_control"),
                 tags$hr(),
                 tags$h4(tags$strong("Reannotation policy")),
                 textInput('exorcise_pam', label = "PAM sequence to append", value = "NGG", placeholder = "e.g. NGG"),
                 selectizeInput('exorcise_genome', label = "Select reference genome", choices = c("GRCh38 Ensembl 111", "GRCm39 Ensembl 111")),
                 selectizeInput('exorcise_mode', label = "Select CRISPR chemistry", choices = c("Knockout", "Interference", "Activation", "Cytosine base editor", "Adenine base editor")),
                 textInput('exorcise_mode_advanced', label = "Advanced: custom chemistry mode string", placeholder = "e.g. beCA0213"),
                 actionButton('exorcise_example', "Example", icon = icon("fire")),
                 actionButton('exorcise_reset', "Reset form", icon = icon("arrows-rotate"), class = "btn-danger"),
                 tags$br(),
                 tags$br(),
                 downloadButton('exorcise_download_input', "Download input", icon = icon("download")),
                 tags$hr(),
                 tags$p(tags$strong("Help")),
                 tags$p("Choose ", tags$strong("sequences"), " to Exorcise by uploading a CRISPick sgRNA designs file or by typing sequences into the box."),
                 tags$p("Sequences must be in [ACGT] and must not contain long runs of the same nucleotide. The length of the sequence plus protospacer-adjacent motif (PAM) must be at least 21 excluding any Ns."),
                 tags$p("Then, choose the ", tags$strong("reannotation policy"), " according to the genome to align to, the PAM to append, and the CRISPR chemistry."),
                 tags$p("The PAM can be in [ACGTN] and will be appended to each sequence during sequence search."),
                 tags$p("Knockout chemistry annotates guides with features occurring at the Cas9 cut site. Interference and activation chemistry finds features up to 500 nucleotides away in linear distance. Cytosine base editor chemistry edits C to G on both strands and assumes a base editing window of [2, 8]. Adenine base editor chemistry edits A to T on both strands and assumes a base editing window of [4, 9]."),
                 tags$p("Custom CRISPR chemistry can be a number indicating the number of nucleotides upstream and downstream of the alignment of the guide + PAM to look for features; or a string in the format `beNM0123`: where `N` and `M` indicate the base edited from and to, respectively, and both are in [ACTG]; and where `01` and `23` indicate the start and end positions of the base editing window, a one-based fully-closed range counting from the PAM-distal end."),
                 tags$p("Examples: `1000`: CRISPRa/i affecting genes up to 1 kb upstream/downstream in linear distance. `beCA0110`: base editing genomic C to edited A within the base editing window [1, 10].")
               ),
               mainPanel(
                 tags$h4("Exorcise"),
                 tags$p("Reannotate CRISPR guides by aligning to other reference genomes using the ", tags$a(href="https://doi.org/10.1186/s13073-024-01414-4", "Exorcise", target="_blank", .noWS = noWS), " algorithm."),
                 tags$p("This tool uses the implementation maintained at the ", tags$a(href="https://github.com/SimonLammmm/exorcise", "Exorcise GitHub repo", target="_blank", .noWS = noWS), "."),
                 tags$br(),
                 actionButton('exorcise_submit', "Exorcise", icon = icon("ghost"), class = "btn-success"),
                 tags$br(),
                 tags$br(),
                 textOutput('exorcise_console'),
                 tags$br(),
                 tags$br(),
                 downloadButton('exorcise_download_result', "Download result", icon = icon("download"))
               )
             )
    )
} else {
  exorciseTab <- NULL
}

#### Shiny UI: Guides ####
if(canGuideLibrary) {
  guideLibraryTab <-
    tabPanel("Guides", icon = icon("worm"),
             sidebarLayout(
               sidebarPanel(
                 width = 3,
                 tags$h4(tags$strong("Target controls")),
                 selectizeInput('guideLibrary_targets', "Targets", choices = NULL, multiple = T),
                 checkboxGroupInput('guideLibrary_geneclasses', "Gene classes",
                                    choices  = c("Protein-coding", "Non-coding RNA", "Pseudogene", "Base edit", "Array", "Other", "Non-targeting"), inline = T),
                 checkboxGroupInput('guideLibrary_chromosomes', "Chromosomes",
                                    choices  = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"), inline = T),
                 tags$hr(),
                 tags$h4(tags$strong("Library controls")),
                 checkboxGroupInput('guideLibrary_chemistry', "CRISPR chemistry",
                                    choices  = c("Knockout", "Interference", "Activation", "Cytosine base editor", "Adenine base editor")),
                 selectizeInput('guideLibrary_assembly', "Assembly", choices = NULL, multiple = T),
                 selectizeInput('guideLibrary_pam', "PAM", choices  = NULL, multiple = T),
                 selectizeInput('guideLibrary_library', "Library", choices = NULL, multiple = T),
                 selectizeInput('guideLibrary_organism', "Organism", choices  = NULL, multiple = T),
                 tags$br(),
                 tags$br(),
                 actionButton('guideLibrary_example', "Example", icon = icon("fire")),
                 actionButton('guideLibrary_reset', "Reset form", icon = icon("arrows-rotate"), class = "btn-danger"),
                 tags$hr(),
                 tags$p(tags$strong("Help")),
                 tags$p("View CRISPR sgRNA sequences used in libraries featured on CRAVE."),
                 tags$p("Use the controls to filter for targets and library designs.")
               ),
               mainPanel(
                 tags$h4("Guides"),
                 tags$p("View CRISPR sgRNA sequences used in libraries featured on CRAVE."),
                 tags$br(),
                 actionButton('guideLibrary_submit', "Guides", icon = icon("worm"), class = "btn-success"),
                 tags$br(),
                 tags$br(),
                 DTOutput('guideLibrary_table'),
                 downloadButton('guideLibrary_download', "Download", icon = icon("download"))
               )
             )
    )
} else {
  guideLibraryTab <- NULL
}

#### Shiny UI: Libraries ####
if(canLibraries) {
  libraryTab <-
    tabPanel("Libraries", icon = icon("book-open"),
             fluidPage(
               tags$h4("Libraries"),
               tags$p("View details of CRISPR libraries featured on CRAVE."),
               tags$br(),
               DTOutput('libraries_table')
             ))
} else {
  libraryTab <- NULL
}

#### Shiny UI: Ontology ####
if(canOntology) {
  ontologyTab <-
    tabPanel("Ontology", icon = icon("sitemap"),
             fluidPage(
               tags$h4("Ontology"),
               tags$p("View members of gene classes."),
               tags$br(),
               DTOutput('ontology_table')
             ))
} else {
  ontologyTab <- NULL
}

#### Shiny UI: legal notices ####
legalTab <-
  legal_text




