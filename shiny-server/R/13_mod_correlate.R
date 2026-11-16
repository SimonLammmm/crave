#### Correlate module ####

# The Correlate sidebar's screen filters. Each entry maps a namespaced input to
# the key of the precomputed choice list and to the filter name understood by
# applyScreenFilter(). Adding a filter is one line here.
CORRELATE_FILTER_FIELDS <- list(
  list(id = "citation",        label = "Citation",           key = "Citation",           filter = "citation"),
  list(id = "kind",            label = "Kind",               key = "Kind",               filter = "kind"),
  list(id = "library",         label = "Library",            key = "Library",            filter = "library"),
  list(id = "timepoint",       label = "Timepoint",          key = "Timepoint",          filter = "timepoint"),
  list(id = "organism",        label = "Species",            key = "Organism",           filter = "organism"),
  list(id = "source",          label = "Source",             key = "Source",             filter = "source"),
  list(id = "days_grown_diff", label = "Days grown (diff)",  key = "Days grown (diff)",  filter = "days_grown_diff"),
  list(id = "days_grown_ref",  label = "Days grown (ref)",   key = "Days grown (ref)",   filter = "days_grown_ref"),
  list(id = "treatment_diff",  label = "Treatment (diff)",   key = "Treatment (diff)",   filter = "treatment_diff"),
  list(id = "treatment_ref",   label = "Treatment (ref)",    key = "Treatment (ref)",    filter = "treatment_ref"),
  list(id = "dose_diff",       label = "Dose (diff)",        key = "Dose (diff)",        filter = "dose_diff"),
  list(id = "dose_ref",        label = "Dose (ref)",         key = "Dose (ref)",         filter = "dose_ref"),
  list(id = "knockout_diff",   label = "Knockout (diff)",    key = "Knockout (diff)",    filter = "knockout_diff"),
  list(id = "knockout_ref",    label = "Knockout (ref)",     key = "Knockout (ref)",     filter = "knockout_ref"),
  list(id = "cellline_diff",   label = "Cell line (diff)",   key = "Cell line (diff)",   filter = "cellline_diff"),
  list(id = "cellline_ref",    label = "Cell line (ref)",    key = "Cell line (ref)",    filter = "cellline_ref")
)

correlateFilterId <- function(f) paste0("f_", f$id)

#### UI ####

correlateSidebarUI <- function(ns, caps = list()) {
  filterInput <- function(f) {
    selectizeInput(ns(correlateFilterId(f)), f$label, choices = NULL, multiple = TRUE)
  }
  byId <- function(ids) {
    Filter(function(f) f$id %in% ids, CORRELATE_FILTER_FIELDS)
  }
  experimentFields <- byId(c("citation", "kind", "library", "timepoint", "organism", "source"))
  # Numerator and denominator fields, interleaved so each grid row is a diff/ref
  # pair.
  levelIds <- c("days_grown_diff", "days_grown_ref", "treatment_diff", "treatment_ref",
                "dose_diff", "dose_ref", "knockout_diff", "knockout_ref",
                "cellline_diff", "cellline_ref")
  levelFields <- CORRELATE_FILTER_FIELDS[match(levelIds,
    vapply(CORRELATE_FILTER_FIELDS, `[[`, character(1), "id"))]

  sidebarPanel(
    width = 3,
    tags$h4(tags$strong("Gene filter controls")),
    textAreaInput(ns("genes_text"), "Type genes of interest",
                  placeholder = "e.g. TP53, BRCA1, UBE2K", rows = 3),
    actionLink(ns("genes_validate"), "Validate"), " · ",
    actionLink(ns("genes_essentials"), "Core essentials"), " · ",
    actionLink(ns("genes_copy"), "Copy GOIs from Explore"),
    tags$br(), tags$br(),
    selectizeInput(ns("genes"), "Search genes of interest", choices = NULL,
                   multiple = TRUE, options = list(placeholder = "e.g. TP53")),
    actionButton(ns("genes_clear"), "Clear genes", icon = icon("ban"), class = "btn-danger"),
    tags$br(), tags$br(),
    selectizeInput(ns("stat"), "Statistic", choices = c("Score", "FDR", "Rank"),
                   selected = "Score"),
    numericInput(ns("cutoff"), "Cutoff", value = defaults_gq_cutoff_selected),
    selectizeInput(ns("method"), "Method", choices = CRAVE_METHOD_NAMES,
                   selected = "DrugZ"),
    tags$hr(),
    tags$h4(tags$strong("Screen filter controls")),
    shinyWidgets::materialSwitch(
      ns("normalise"), "Quantile normalise scores across screens",
      value = FALSE, status = "primary", right = TRUE
    ),
    tags$p(tags$em(
      "Gives every screen the same score distribution, using the screens selected",
      "below as the reference. Gene rankings within a screen are unchanged, so FDR",
      "and rank cutoffs behave the same; score cutoffs apply to normalised scores."
    ), style = "font-size: 85%"),
    tags$p(tags$strong("Experiment")),
    do.call(inputGrid, lapply(experimentFields, filterInput)),
    tags$table(
      style = "width: 100%",
      tags$tr(
        style = "vertical-align: top",
        tags$td(style = "padding: 6px; width: 50%", tags$p(tags$strong("Numerator levels"))),
        tags$td(style = "padding: 6px; width: 50%", tags$p(tags$strong("Denominator levels")))
      )
    ),
    do.call(inputGrid, lapply(levelFields, filterInput)),
    tags$p(tags$strong("Custom")),
    selectizeInput(ns("f_contrast"), "Filter by contrast", choices = NULL, multiple = TRUE),
    actionLink(ns("contrast_copy"), "Copy contrasts from Explore"),
    tags$br(), tags$br(),
    textInput(ns("f_custom"), "Filter by text input"),
    shinyWidgets::materialSwitch(
      ns("f_custom_regex"), "Treat as a regular expression",
      value = FALSE, status = "primary", right = TRUE
    ),
    tags$p(tags$em("Perl-compatible syntax, e.g. ^Olaparib|Talazoparib"),
           style = "font-size: 85%"),
    textOutput(ns("custom_regex_error")),
    textOutput(ns("filtered_screens")),
    tags$br(), tags$br(),
    actionButton(ns("example"), "Example", icon = icon("fire")),
    actionButton(ns("reset"), "Reset form", icon = icon("arrows-rotate"), class = "btn-danger"),
    tags$br(), tags$br(),
    customiseUI(ns("customise")),
    tags$hr(),
    tags$p(tags$strong("Help")),
    tags$p("Use", tags$strong("gene filters"), "to choose genes of interest. Leave it",
           "blank to consider all genes. Use the cutoff to constrain for hit strength."),
    tags$p("Use", tags$strong("screen filters"), "to choose screens of interest.",
           "Screens satisfying at least one filter in each field are considered.",
           "Screens of interest are previewed below the fields."),
    tags$p("Once you've chosen your filters, click a tab to run an analysis:"),
    tags$p(tags$strong("Gene query"), "and", tags$strong("Heatmap"),
           "show CRISPR scores across screens."),
    tags$p(tags$strong("Adjacency"), "and", tags$strong("Network"),
           "show pairwise similarity between genes.", tags$strong("Network"),
           "can also group the genes into communities and ring them."),
    # Both links are followed immediately by a comma, so the trailing space is
    # suppressed but the leading one is kept.
    tags$p(tags$strong("Reduce"), "projects genes into two dimensions with",
           extLinkTight(DOI_UMAP, "UMAP"), ", PCA or",
           extLinkTight(DOI_TSNE, "t-SNE"),
           ", revealing clusters of similarly behaving genes."),
    tags$p(tags$strong("Enrichment"), "uses the hypergeometric test to reveal",
           "enrichment and depletion of ontological terms among hits."),
    if (isTRUE(caps$pendragonator)) tags$p(
      tags$strong("Pendragonator"), "returns the top responding genes among",
      "your screens of interest and the conditions that generated them.") else NULL,
    tags$p(tags$strong("Data download"), "lets you bulk download CRISPR scores.")
  )
}

#' Clustering controls shared by the Heatmap and Adjacency tabs.
#'
#' `prefix` namespaces the three inputs so each tab keeps its own choice.
clusteringControls <- function(ns, prefix) {
  id <- function(suffix) ns(paste0(prefix, "_", suffix))
  tagList(
    inputGrid(
      selectizeInput(id("dendro"), "Clustering method",
                     choices = CRAVE_DENDRO_METHODS, selected = CRAVE_DENDRO_DEFAULT),
      selectizeInput(id("dist"), "Distance measure",
                     choices = CRAVE_DISTANCE_AGNES, selected = CRAVE_DISTANCE_DEFAULT)
    ),
    conditionalPanel(
      condition = sprintf("input['%s'] == 'Clara'", id("dendro")),
      numericInput(id("k"), "Number of clusters (k)",
                   value = CRAVE_CLARA_K_DEFAULT, min = 2, max = NA, step = 1),
      tags$p(tags$em(
        "Clara partitions rather than building a tree, so each cluster takes a",
        "contiguous block of the axis and the dendrogram shown is over the cluster",
        "medoids. k must be at least 2 and fewer than the number of items on the axis."
      ), style = "font-size: 85%")
    )
  )
}

correlateUI <- function(id, caps) {
  ns <- NS(id)

  enrichmentTab <- if (isTRUE(caps$ontology)) {
    tabPanel(
      "Enrichment", icon = icon("blender"),
      plotPanelUI(
        ns("enrichment"), "Enrichment",
        c("Generate a heatmap showing enrichment or depletion of ontological terms among the hits within your screens of interest.",
          "Use the controls to define the cutoff for a hit and to choose a direction and test tail.",
          "Significantly enriched or depleted terms are shown with an asterisk.",
          "Cells are coloured by the number of hits as a proportion of all genes in that ontological term (recall).",
          paste0("Switch on <strong>Venn compartments</strong> to test each exclusive ",
                 "combination of screens separately &mdash; genes hit in A alone, in A ",
                 "and B but not C, and so on &mdash; instead of one column per screen. ",
                 "The x axis then becomes an upset matrix naming the screens in each ",
                 "compartment. Limited to ", CRAVE_MAX_VENN_SETS, " screens.")),
        controls = tagList(
          inputGrid(
            selectizeInput(ns("enr_stat"), "Cutoff statistic", choices = c("Score", "FDR", "Rank")),
            textInput(ns("enr_value"), "Better than", value = CRAVE_ENRICHMENT_CUTOFF),
            selectizeInput(ns("enr_direction"), "Direction",
                           choices = c("Hypersensitivity hits", "Suppressing hits")),
            selectizeInput(ns("enr_tail"), "Test tail",
                           choices = c("Upper tail (enriched classes)", "Lower tail (depleted classes)"))
          ),
          shinyWidgets::materialSwitch(
            ns("enr_compartments"), "Venn compartments",
            value = FALSE, status = "primary", right = TRUE
          )
        ),
        submitLabel = "Enrichment analysis", submitIcon = "blender"
      )
    )
  } else NULL

  # Pendragonator is the one analysis a deployment may withhold without
  # withholding the tab it lives in, so it is built here and spliced in below.
  pendragonatorTab <- if (isTRUE(caps$pendragonator)) {
      tabPanel(
        "Pendragonator", icon = icon("dragon"),
        plotPanelUI(
          ns("pendragonator"), "Pendragonator",
          c("Generate a list of genes in your screens of interest using Pendragon library design principles.",
            "Returns the strongest responding genes on either tail and the conditions that produced the hit.",
            "Tick the box to find the screens that produced the best hit among only your genes of interest."),
          controls = inputGrid(
            numericInput(ns("pdg_n"), "List length", value = 2000, min = 1,
                         max = 20000, step = 1),
            selectizeInput(ns("pdg_stat"), "Statistic", choices = c("Score", "FDR")),
            checkboxInput(ns("pdg_goi"), "Pendragonate only genes of interest", value = FALSE),
            ncol = 3
          ),
          submitLabel = "Pendragonate", submitIcon = "dragon",
          tableFirst = TRUE
        )
      )
  } else NULL

  tabPanel(
    "Correlate", icon = icon("chart-line"),
    sidebarLayout(
      correlateSidebarUI(ns, caps),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Gene query", icon = icon("magnifying-glass-chart"),
            plotPanelUI(
              ns("genequery"), "Gene query",
              c("Generate violin plots of CRISPR scores and show where your GOIs were hits among your screens of interest.",
                paste0("Leave the genes of interest field empty to query the top ", CRAVE_MAX_AUTO_GOI, " genes."),
                "GOIs that pass the cutoff are shown as triangles, all others as circles. Only screens with at least one triangle are shown.",
                "Double-click in the legend to isolate a plot element. Single-click to toggle elements on and off."),
              controls = selectizeInput(ns("gq_plottype"), "Plot type",
                                        choices = c("Violin plot", "Boxplot"),
                                        selected = "Violin plot"),
              submitLabel = "Gene query", submitIcon = "magnifying-glass-chart"
            )
          ),

          tabPanel(
            "Heatmap", icon = icon("fire-flame-curved"),
            plotPanelUI(
              ns("heatmap"), "Heatmap",
              c("Generate a heatmap of CRISPR scores for your GOIs among your screens of interest.",
                paste0("Leave the genes of interest field empty to query the top ", CRAVE_MAX_AUTO_GOI, " genes."),
                paste0("Genes and screens are arranged into clusters of similar behaviour. ",
                       "Dendrograms are shown when there are at least three GOIs or at least ",
                       "three screens of interest. Choose the clustering method and distance ",
                       "measure below; <a href=\"", DOI_AGNES, "\" target=\"_blank\">",
                       "agglomerative nesting (Agnes)</a> with Ward linkage is the default."),
                "Bigger scores have greater colour intensity. Hits exceeding the cutoff are shown with an asterisk."),
              controls = tagList(
                checkboxInput(ns("hm_gene_dendro"), "Draw dendrogram on genes", TRUE),
                checkboxInput(ns("hm_screen_dendro"), "Draw dendrogram on screens", TRUE),
                clusteringControls(ns, "hm")
              ),
              submitLabel = "Heatmap", submitIcon = "fire-flame-curved"
            )
          ),

          tabPanel(
            "Adjacency", icon = icon("table-cells"),
            plotPanelUI(
              ns("adjacency"), "Adjacency",
              c("Generate a gene-gene similarity matrix using evidence from your screens of interest.",
                "Genes are arranged into clusters of similar behaviour. A dendrogram is shown when there are at least three GOIs.",
                "Spearman, Pearson and Kendall are rank or product-moment correlations; cosine similarity compares the direction of each gene's score vector. All four run from -1 to 1."),
              controls = tagList(
                selectizeInput(ns("adj_similarity"), "Similarity measure",
                               choices = CRAVE_SIMILARITY_METHODS,
                               selected = CRAVE_SIMILARITY_DEFAULT),
                clusteringControls(ns, "adj")
              ),
              submitLabel = "Adjacency", submitIcon = "table-cells"
            )
          ),

          tabPanel(
            "Network", icon = icon("circle-nodes"),
            plotPanelUI(
              ns("network"), "Network",
              c("Draw a network representation of gene-gene similarity using evidence from your screens of interest.",
                "Only pairs at least as strong as the cutoff are shown, in either direction.",
                paste0("Switch on <strong>community detection</strong> to group the genes ",
                       "drawn and ring each group. The similarity measure decides which ",
                       "genes are joined; the distance measure decides which are near ",
                       "each other for grouping, as it does on the Heatmap and Adjacency ",
                       "tabs. Agnes, Clara and k-means put each gene in one community; ",
                       "the factorisation and the mixture model give each gene a ",
                       "membership in every community, so their ellipses overlap and the ",
                       "results table carries the membership as a probability. Louvain and ",
                       "Leiden group the network itself rather than a distance between ",
                       "gene profiles, so they take a resolution in place of a k and find ",
                       "their own number of communities."),
                "The results table is an adjacency list: one row per ordered pair, so each pair appears twice, with the community of each end."),
              controls = tagList(
                inputGrid(
                  selectizeInput(ns("nw_similarity"), "Similarity measure",
                                 choices = CRAVE_SIMILARITY_METHODS,
                                 selected = CRAVE_SIMILARITY_DEFAULT),
                  sliderInput(ns("corr_cutoff"), "Similarity cutoff",
                              min = defaults_corr_cutoff_min,
                              max = defaults_corr_cutoff_max,
                              step = defaults_corr_cutoff_step,
                              value = defaults_corr_cutoff_value)
                ),
                selectizeInput(ns("nw_community"), "Community detection",
                               choices = CRAVE_COMMUNITY_METHODS,
                               selected = CRAVE_COMMUNITY_DEFAULT),
                # Each of the three blocks below is shown only for the methods
                # it applies to, rather than all of them whenever a method is
                # chosen: an inert control beside a dropdown invites the reading
                # that it does something. k and the distance are for the methods
                # that take them, the resolution is for Louvain and Leiden, and
                # the membership is for the two soft fits.
                #
                # The method lists are written into the JS by hand rather than
                # with a JSON writer: this is the only place CRAVE would need
                # one, and the names are its own constants, so there is nothing
                # here to escape.
                conditionalPanel(
                  condition = paste0(
                    "input['", ns("nw_community"), "'] != '", CRAVE_COMMUNITY_OFF,
                    "' && ['", paste(CRAVE_COMMUNITY_GRAPH, collapse = "','"),
                    "'].indexOf(input['", ns("nw_community"), "']) < 0"),
                  inputGrid(
                    numericInput(ns("nw_k"), "Number of communities (k)",
                                 value = CRAVE_CLARA_K_DEFAULT, min = 2, step = 1),
                    selectizeInput(ns("nw_dist"), "Distance measure",
                                   choices = CRAVE_DISTANCES,
                                   selected = CRAVE_DISTANCE_DEFAULT)
                  )
                ),
                conditionalPanel(
                  condition = paste0(
                    "['", paste(CRAVE_COMMUNITY_GRAPH, collapse = "','"),
                    "'].indexOf(input['", ns("nw_community"), "']) >= 0"),
                  numericInput(ns("nw_resolution"), "Resolution",
                               value = CRAVE_COMMUNITY_RESOLUTION_DEFAULT,
                               min = CRAVE_COMMUNITY_RESOLUTION_MIN,
                               max = CRAVE_COMMUNITY_RESOLUTION_MAX, step = 0.1),
                  tags$p(tags$em(
                    "These methods group the network drawn above, so the",
                    "similarity measure and the cutoff decide the communities and",
                    "no distance is used. They also choose how many communities",
                    "there are: raise the resolution for more and smaller ones,",
                    "lower it for fewer and larger. The diagram says how many it",
                    "found."), style = "font-size: 85%")
                ),
                conditionalPanel(
                  condition = paste0(
                    "['", paste(CRAVE_COMMUNITY_PROBABILISTIC, collapse = "','"),
                    "'].indexOf(input['", ns("nw_community"), "']) >= 0"),
                  numericInput(ns("nw_membership"),
                               "Least community membership to be ringed",
                               value = CRAVE_COMMUNITY_MEMBERSHIP_DEFAULT,
                               min = 0, max = 1, step = 0.05),
                  tags$p(tags$em(
                    "These methods give every gene a membership in every",
                    "community. A gene is ringed with each community it holds",
                    "at least this much of, so a gene may be ringed twice and",
                    "the ellipses may overlap."), style = "font-size: 85%")
                )
              ),
              submitLabel = "Network", submitIcon = "circle-nodes"
            )
          ),

          tabPanel(
            "Reduce", icon = icon("minimize"),
            plotPanelUI(
              ns("reduce"), "Reduce",
              c(paste0("Project your GOIs into two dimensions using ",
                       "<a href=\"", DOI_UMAP, "\" target=\"_blank\">UMAP</a>, ",
                       "principal component analysis, or ",
                       "<a href=\"", DOI_TSNE, "\" target=\"_blank\">t-SNE</a>."),
                "Missing values are dealt with by imputation before projection.",
                "UMAP needs at least four screens, t-SNE at least three, PCA at least two. t-SNE also needs enough genes to support a perplexity of 2 or more."),
              controls = selectizeInput(ns("reduce_method"), "Method",
                                        choices = CRAVE_REDUCE_METHODS,
                                        selected = CRAVE_REDUCE_DEFAULT),
              submitLabel = "Reduce", submitIcon = "minimize"
            )
          ),

          enrichmentTab,

          pendragonatorTab,

          tabPanel(
            "Data download", icon = icon("download"),
            tags$br(),
            tags$h4("Data download"),
            tags$p("Bulk download CRISPR data."),
            tags$p("Specify the file format by appending \".csv\" (comma-separated values)",
                   "or \".tsv\" (tab-separated values). Append \".gz\" for gzipped output."),
            tags$p("Click \"Fetch\" to preview the output and then \"Download data\" to start the download."),
            tags$br(),
            textInput(ns("dl_filename"), "Filename", value = "export.tsv.gz"),
            tags$br(),
            actionButton(ns("dl_submit"), "Fetch", icon = icon("download"), class = "btn-success"),
            tags$br(),
            textOutput(ns("dl_console")),
            tags$br(),
            DTOutput(ns("dl_table")),
            tags$br(),
            downloadButton(ns("dl_download"), "Download data", icon = icon("download"))
          )
        )
      )
    )
  )
}

#### Server ####

correlateServer <- function(id, data, bus) {
  moduleServer(id, function(input, output, session) {

    customise <- customiseServer("customise", CORRELATE_CUSTOMISE_FIELDS, CORRELATE_DEFAULTS)
    opts <- reactive(customise$values())

    #### Initialisation ####
    observeEvent(data(), {
      d <- data()
      for (f in CORRELATE_FILTER_FIELDS) {
        updateSelectizeInput(session, correlateFilterId(f), server = TRUE,
                             choices = d$choices[[f$key]])
      }
      # Source defaults to everything, so a fresh session sees all datasets.
      updateSelectizeInput(session, "f_source", server = TRUE,
                           choices = d$choices$Source, selected = d$choices$Source)
      updateSelectizeInput(session, "f_contrast", server = TRUE,
                           choices = d$choices$FriendlyID,
                           options = list(maxOptions = 1000))
      updateSelectizeInput(session, "genes", server = TRUE, choices = d$genes,
                           options = list(maxOptions = CRAVE_MAX_SELECTIZE_SERVER))
    })

    #### Filters ####
    filters <- reactive({
      out <- list(
        citation    = input$f_citation,
        custom      = input$f_custom,
        customRegex = isTRUE(input$f_custom_regex),
        contrast    = input$f_contrast
      )
      for (f in CORRELATE_FILTER_FIELDS) {
        out[[f$filter]] <- input[[correlateFilterId(f)]]
      }
      out
    })

    normalise <- reactive(isTRUE(input$normalise))

    # Debounced: typing in the free-text filter should not re-scan the comparison
    # table on every keystroke.
    screens <- debounce(reactive(applyScreenFilter(data(), filters())), 300)

    # A part-typed regular expression is not an error worth interrupting the user
    # for, but it does need explaining, because it silently matches nothing.
    output$custom_regex_error <- renderText({
      msg <- validateTextFilter(input$f_custom, isTRUE(input$f_custom_regex))
      if (is.null(msg)) "" else paste0("Invalid regular expression: ", msg)
    })

    output$filtered_screens <- renderText({
      found <- sort(friendlyIdsFor(data(), screens()))
      n <- length(found)
      if (n == 0) return("No screens found with the selected filters.")
      paste0(if (n == 1) "Screen found: " else "Screens found: ", summariseItems(found))
    })

    # jaccard is a clara-only metric, so the distance choices follow the method.
    for (prefix in c("hm", "adj")) {
      local({
        p <- prefix
        methodId <- paste0(p, "_dendro")
        distId   <- paste0(p, "_dist")
        observeEvent(input[[methodId]], {
          allowed <- if (identical(input[[methodId]], "Clara")) CRAVE_DISTANCE_CLARA
                     else CRAVE_DISTANCE_AGNES
          # %||% first: a NULL here would make the %in% test length zero, which
          # `if` refuses.
          current <- isolate(input[[distId]]) %||% CRAVE_DISTANCE_DEFAULT
          keep <- if (current %in% allowed) current else CRAVE_DISTANCE_DEFAULT
          updateSelectizeInput(session, distId, choices = allowed, selected = keep)
        }, ignoreInit = TRUE)
      })
    }

    #### Gene selection ####
    observeEvent(input$genes_validate, {
      txt <- input$genes_text %||% ""
      if (nchar(txt) < 3) return()
      pool  <- data()$genes
      typed <- unique(splitTokens(txt))
      found <- unique(pool[toupper(pool) %in% toupper(typed)])
      missed <- setdiff(toupper(typed), toupper(found))
      if (length(missed)) {
        showNotification(
          paste0(length(missed), " symbol(s) not recognised: ",
                 paste(utils::head(missed, 10), collapse = ", "),
                 if (length(missed) > 10) ", ..." else ""),
          type = "warning", duration = 6
        )
      }
      updateSelectizeInput(session, "genes", selected = unique(c(input$genes, found)),
                           choices = pool,
                           server = length(found) <= CRAVE_MAX_SELECTIZE_SERVER)
    })

    observeEvent(input$genes_essentials, {
      updateTextAreaInput(session, "genes_text",
                          value = paste(CRAVE_ESSENTIALS, collapse = ", "))
    })

    observeEvent(input$genes_clear, {
      updateSelectizeInput(session, "genes", selected = character(0),
                           choices = data()$genes, server = TRUE)
    })

    observeEvent(input$genes_copy, {
      updateSelectizeInput(session, "genes", selected = bus$exploreGenes,
                           choices = data()$genes, server = TRUE)
    })

    observeEvent(input$contrast_copy, {
      updateSelectizeInput(session, "f_contrast", selected = bus$exploreContrasts,
                           choices = data()$choices$FriendlyID, server = TRUE,
                           options = list(maxOptions = 1000))
    })

    observe({
      bus$correlateGenes <- input$genes
      bus$correlateContrasts <- input$f_contrast
    })

    observeEvent(bus$brush, {
      b <- bus$brush
      if (is.null(b) || !identical(b$target, "Correlate") || length(b$genes) == 0) return()
      updateSelectizeInput(session, "genes", selected = unique(c(input$genes, b$genes)),
                           choices = data()$genes, server = TRUE)
    }, ignoreInit = TRUE)

    #### Reset and example ####
    resetForm <- function() {
      d <- data()
      for (f in CORRELATE_FILTER_FIELDS) {
        updateSelectizeInput(session, correlateFilterId(f), selected = character(0))
      }
      updateSelectizeInput(session, "f_source", selected = d$choices$Source)
      updateSelectizeInput(session, "f_contrast", selected = character(0))
      updateTextInput(session, "f_custom", value = "")
      updateNumericInput(session, "cutoff", value = defaults_gq_cutoff_selected)
      updateSelectizeInput(session, "stat", selected = "Score")
      updateSelectizeInput(session, "method", selected = "DrugZ")
      shinyWidgets::updateMaterialSwitch(session, "f_custom_regex", value = FALSE)
      shinyWidgets::updateMaterialSwitch(session, "normalise", value = FALSE)
      # The two tab controls the form owns, because the example sets them:
      # leaving them out would make Reset unable to undo what Example did.
      updateSliderInput(session, "corr_cutoff", value = defaults_corr_cutoff_value)
      updateTextInput(session, "enr_value", value = CRAVE_ENRICHMENT_CUTOFF)
    }

    observeEvent(input$reset, resetForm())

    observeEvent(input$example, {
      d <- data()
      resetForm()
      updateTextAreaInput(session, "genes_text",
                          value = paste(CRAVE_EXAMPLE_GENES, collapse = ", "))
      updateSelectizeInput(session, "genes", selected = CRAVE_EXAMPLE_GENES,
                           choices = d$genes, server = TRUE)
      # Only pre-select example filter values that this dataset actually has.
      if ("Human" %in% d$choices$Organism) {
        updateSelectizeInput(session, "f_organism", selected = "Human")
      }
      if ("Treatment" %in% d$choices$Kind) {
        updateSelectizeInput(session, "f_kind", selected = "Treatment")
      }
      drugs <- intersect(CRAVE_EXAMPLE_TREATMENTS, d$choices$`Treatment (diff)`)
      if (length(drugs) > 0) {
        updateSelectizeInput(session, "f_treatment_diff", selected = drugs)
      }
      if ("Final timepoint" %in% d$choices$Timepoint) {
        updateSelectizeInput(session, "f_timepoint", selected = "Final timepoint")
      }
      if ("Wildtype" %in% d$choices$`Knockout (diff)`) {
        updateSelectizeInput(session, "f_knockout_diff", selected = "Wildtype")
      }
      if ("Wildtype" %in% d$choices$`Knockout (ref)`) {
        updateSelectizeInput(session, "f_knockout_ref", selected = "Wildtype")
      }
      # Leave out the screen variants and alternative-library runs, which are
      # what an underscore in a screen's name almost always marks: the regular
      # ones read "Aldo (NVS057): Etoposide", the variants "Alex (NVS163_PeCR):
      # Etoposide" or "...: GeneKO_G1_A549_48h". Matching is inclusive, so the
      # rule is spelled as the labels to keep -- no underscore anywhere.
      shinyWidgets::updateMaterialSwitch(session, "f_custom_regex", value = TRUE)
      updateTextInput(session, "f_custom", value = "^[^_]*$")
      # Score, hits at |normZ| >= 3: the threshold the Enrichment tab already
      # defaults to, and loose enough that every tab has something to draw.
      updateSelectizeInput(session, "stat", selected = "Score")
      updateNumericInput(session, "cutoff", value = 3)
      updateSelectizeInput(session, "method", selected = "DrugZ")
      # Tighter than the Network tab's own 0.7 default. How many edges a
      # similarity threshold admits depends on how many screens are being
      # correlated -- over the example dataset's nine, 0.7 draws 2,846 edges
      # between all 200 genes, which is a hairball; over the forty-odd a large
      # dataset offers, 0.9 draws eight. 0.85 is the value that leaves a
      # readable graph either way.
      updateSliderInput(session, "corr_cutoff", value = 0.85)
      # Looser than the Enrichment tab's own cutoff, the other way: a handful of
      # screens give a few dozen hits each, and at |normZ| >= 3 the classes that
      # are plainly enriched -- homologous recombination under a PARP inhibitor
      # -- do not survive multiple testing. At 2 they do.
      updateTextInput(session, "enr_value", value = "2")
      # Plot height belongs to the customise modal, whose inputs exist only while
      # it is open. The store behind them is what the plots read, and it is also
      # what the modal is rebuilt from, so this is the one way in. Taller than the
      # 900px default because the heatmap draws a row per gene of interest.
      customise$restore(list(height = 1500))
    })

    #### Analyses ####
    plotPanelServer("genequery", busyMessage = "Querying genes...",
      compute = function() {
        runGeneQuery(data(), screens(), input$genes, input$stat, input$cutoff,
                     input$method, input$gq_plottype,
                     normalise = normalise(), opts = opts())
      },
      onResult = function(res) publishBrush(bus, "Correlate", res$brush))

    plotPanelServer("heatmap", busyMessage = "Clustering...",
      compute = function() {
        runHeatmap(data(), screens(), input$genes, input$stat, input$cutoff,
                       input$method, input$hm_gene_dendro, input$hm_screen_dendro,
                       dendroMethod = input$hm_dendro, distMethod = input$hm_dist,
                       k = input$hm_k, normalise = normalise(), opts = opts())
      })

    plotPanelServer("adjacency", busyMessage = "Correlating genes...",
      compute = function() {
        runAdjacency(data(), screens(), input$genes, input$stat, input$cutoff,
                   input$method, similarity = input$adj_similarity,
                   dendroMethod = input$adj_dendro, distMethod = input$adj_dist,
                   k = input$adj_k, normalise = normalise(), opts = opts())
      })

    plotPanelServer("network", busyMessage = "Building network...",
      compute = function() {
        runNetwork(data(), screens(), input$genes, input$stat, input$cutoff,
                   input$method, similarity = input$nw_similarity,
                   corrCutoff = input$corr_cutoff, normalise = normalise(),
                   community = input$nw_community, communityK = input$nw_k,
                   communityDist = input$nw_dist,
                   communityResolution = input$nw_resolution,
                   membership = input$nw_membership,
                   opts = opts())
      })

    plotPanelServer("reduce", busyMessage = "Imputing and projecting...",
      compute = function() {
        runReduce(data(), screens(), input$genes, input$stat, input$cutoff,
                  input$method, reduceMethod = input$reduce_method,
                  normalise = normalise(), opts = opts())
      },
      onResult = function(res) publishBrush(bus, "Correlate", res$brush))

    plotPanelServer("enrichment", busyMessage = "Testing enrichment...",
      compute = function() {
        runEnrichment(data(), screens(), input$method, input$enr_stat,
                      input$enr_value, input$enr_direction, input$enr_tail,
                      compartments = isTRUE(input$enr_compartments),
                      normalise = normalise())
      })

    plotPanelServer("pendragonator", busyMessage = "Pendragonating...",
      compute = function() {
        list(
          plotdata = runPendragonator(data(), screens(), input$genes, input$method,
                                      input$pdg_stat, input$pdg_n, input$pdg_goi,
                                      normalise = normalise()),
          renderer = "none"
        )
      })

    #### Data download ####
    downloadResult <- eventReactive(input$dl_submit, {
      withProgress(message = "Fetching data...", value = NULL, {
        runBulkDownload(data(), screens(), input$dl_filename %||% "export.tsv.gz",
                        normalise = normalise())
      })
    }, ignoreInit = TRUE)

    output$dl_console <- renderText(downloadResult()$console)
    output$dl_table <- renderDT(
      utils::head(downloadResult()$table, 6),
      server = TRUE, filter = "top", escape = FALSE,
      options = list(scrollX = TRUE)
    )
    output$dl_download <- downloadHandler(
      filename = function() downloadResult()$filename,
      content  = function(f) {
        res <- downloadResult()
        data.table::fwrite(res$table, f, sep = res$delim)
      }
    )

    #### Save/Load ####
    stateIds <- c(
      "genes_text", "genes", "stat", "cutoff", "method", "normalise",
      vapply(CORRELATE_FILTER_FIELDS, correlateFilterId, character(1)),
      "f_contrast", "f_custom", "f_custom_regex",
      "gq_plottype",
      "hm_gene_dendro", "hm_screen_dendro", "hm_dendro", "hm_dist", "hm_k",
      "adj_similarity", "adj_dendro", "adj_dist", "adj_k",
      "nw_similarity", "corr_cutoff",
      "nw_community", "nw_k", "nw_dist", "nw_resolution", "nw_membership",
      "reduce_method",
      "enr_stat", "enr_value", "enr_direction", "enr_tail", "enr_compartments",
      "pdg_n", "pdg_stat", "pdg_goi", "dl_filename"
    )
    textStateIds    <- c("genes_text", "f_custom", "enr_value", "dl_filename")
    sliderStateIds  <- c("corr_cutoff")
    switchStateIds  <- c("normalise", "f_custom_regex", "enr_compartments")
    numericStateIds <- c("cutoff", "hm_k", "adj_k", "pdg_n", "nw_k", "nw_membership",
                         "nw_resolution")

    list(
      state = function() {
        st <- captureInputs(input, stateIds)
        st$`_customise` <- customise$values()
        st
      },
      restore = function(st) {
        # A key saved before 6.0.0 calls the Adjacency tab's clustering controls
        # hm_dendro, hm_dist and hm_k -- which are the names the Heatmap tab uses
        # now. Restoring them is worse than skipping them: the values are valid
        # for the control they would land on, so the old Adjacency settings would
        # silently appear on the Heatmap tab. Keys from that era are recognisable
        # by their cg_* ids, which nothing uses any more.
        if (any(c("cg_dendro", "cg_dist", "cg_k") %in% names(st))) {
          st[c("hm_dendro", "hm_dist", "hm_k")] <- NULL
          log_info("Loaded a pre-6.0.0 key; its clustering settings were left at the defaults.")
        }
        d <- data()
        # Every one of these was created with server = TRUE, so its choices have to
        # be resent alongside the restored selection.
        serverChoices <- list(genes = d$genes, f_contrast = d$choices$FriendlyID)
        for (f in CORRELATE_FILTER_FIELDS) {
          serverChoices[[correlateFilterId(f)]] <- d$choices[[f$key]]
        }
        restoreInputs(session, st, stateIds, textStateIds, sliderStateIds,
                      switchStateIds, numericStateIds, serverChoices)
        customise$restore(st$`_customise`)
        invisible(NULL)
      }
    )
  })
}
