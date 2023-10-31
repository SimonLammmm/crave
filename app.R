
# Author    Version       Date          Description
# SL        2.0.0         2023-08-03    Complete rework. 
# SL        2.1.0         2023-08-08    Enable lazy-loading using SQL data structure.
# SL        2.1.1         2023-08-10    Enable reference level KO in view.
# SL        2.1.2         2023-08-11    Add busy indicator. Fix edge cases with screen selection.
# SL        2.2.0         2023-08-29    Add exorcise service.

ver = "2.2.0"

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
})

options(digits = 2)
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)

queryDdrcs <- function(experiment, gene, method, type) {
  q <- "SELECT "
  if (is.null(experiment)) {
    q <- paste0(q, "* FROM ")
  } else {
    q <- paste0(q, "[gene], ", paste0("[", experiment, "]", collapse = ","), " FROM ")
  }
  if (method == "DrugZ") {
    q <- paste0(q, "drz_", type)
  } else if (method == "MAGeCK") {
    q <- paste0(q, "mag_", type)
  }
  if (!is.null(gene)) {
    q <- paste0(q, " WHERE [gene] = ", paste0("'", gene, "'", collapse = " OR [gene] = "))
  }
  return (q)
}

ui <-
  navbarPage(theme = shinytheme("flatly"),
             HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">CRAVE</a>'),
             windowTitle = "CRAVE",
             
             #### Homepage ####
             tabPanel("Home",
                      tags$div(
                        h1("CRAVE: CRISPR-screen Rshiny app for visualisation and exploration"),tags$br(),
                        tags$br(),tags$h4("Welcome to CRAVE"), 
                        "CRAVE is being developed to unify all CRISPR screens done in the lab and to allow unconventional exploration of trends across genetic backgrounds and treatments.",tags$br(),
                        tags$br(),
                        "This app is under development. Some features will not work. Please report any bugs to the authors, contactable as below.",tags$br(),
                        tags$br(),
                        tags$h4("Authors"),
                        "Dr Simon Lam, Steve Jackson Lab, Cancer Research UK Cambridge Institute, University of Cambridge, UK",tags$br(),
                        tags$br(),tags$br(),
                        tags$h4("Contact"),
                        "sl681@cam.ac.uk",tags$br(),tags$br(),
                        tags$a(href='https://www.cam.ac.uk//', target="_blank", tags$img(src = "Colour logo RGB_DM.jpg", width = "361px", height = "75px")),
                        tags$a(href='https://www.stevejacksonlab.org/', target="_blank", tags$img(src = "new_lab_logo.png", width = "306px", height = "75px")),
                        tags$a(href='https://www.cruk.cam.ac.uk/', target="_blank", tags$img(src = "CRUK_CAMBRIDGE_I_Pos_CMYK.jpg", width = "341px", height = "75px"))
                      )
             ),
             
             
             #### Gene query ####
             tabPanel("Gene query",
                      sidebarLayout(
                        sidebarPanel(
                          width = 3,
                          selectizeInput('genes_queried', label = "Genes of interest", choices = NULL, multiple = TRUE),
                          sliderInput('p_cutoff', label = "p-value cutoff", min = 0, max = 1, step = 0.001, value = 0.05),
                          selectizeInput('analysis_method_gq', label = "Method", choices = c("DrugZ", "MAGeCK")),
                          actionButton('genequery_submit', "Gene query")
                        ),
                        mainPanel(
                          fluidRow(
                            plotlyOutput('genequery_plot', height = 900)
                          ),
                          fluidRow(
                            DTOutput('genequery_table'),
                            downloadButton('genequery_download', "Download")
                          )
                        )
                      )
            ),
             
             #### Screen explorer ####
             tabPanel("Screen explorer",
                      sidebarLayout(
                        sidebarPanel(
                          width = 3,
                          selectizeInput('comparisons_selected', label = "Selected comparisons", choices = NULL, multiple = TRUE),
                          selectizeInput('genes_selected', label = "Genes of interest", choices = NULL, multiple = TRUE)
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Select experiment",
                                     fluidRow(
                                       actionButton('experiments_table_selectAll', "Select all experiments in view"),
                                       actionButton('experiments_table_selectNone', "Deselect all experiments"),
                                       DTOutput('experiments_table')
                                     )
                            ),
                            tabPanel("Select comparison",
                                     fluidRow(
                                       actionButton('comparisons_table_selectAll', "Select all comparisons in view"),
                                       actionButton('comparisons_table_selectNone', "Deselect all comparisons"),
                                       DTOutput('comparisons_table')
                                     )
                            ),
                            tabPanel("Volcano plot",
                                     fluidRow(
                                       column(
                                         width = 4.5,
                                         selectizeInput('comparison_volcano', label = "Select one comparison", choices = NULL)
                                         ),
                                       column(
                                         width = 4.5,
                                         selectizeInput('analysis_method_se', label = "Method", choices = c("DrugZ", "MAGeCK")),
                                         actionButton('volcano_submit', 'Generate volcano plot')
                                       )
                                     ),
                                     fluidRow(
                                       plotlyOutput('volcano_plot', height = 900),
                                       DTOutput('volcano_table'),
                                       downloadButton('volcano_download', "Download")
                                     )
                            ),
                            tabPanel("Rank plot",
                                     fluidRow(
                                       column(
                                         width = 4.5,
                                         selectizeInput('comparison_rank', label = "Select one comparison", choices = NULL)
                                         ),
                                       column(
                                         width = 4.5,
                                         selectizeInput('analysis_method_se', label = "Method", choices = c("DrugZ", "MAGeCK")),
                                         actionButton('rank_submit', 'Generate rank plot')
                                       )
                                     ),
                                     fluidRow(
                                       plotlyOutput('rank_plot', height = 900),
                                       DTOutput('rank_table'),
                                       downloadButton('rank_download', "Download")
                                     )
                            ),
                            tabPanel("Biplot",
                                     fluidRow(
                                       column(
                                         width = 4.5,
                                         selectizeInput('comparison_X', label = "Select X comparison", choices = NULL),
                                         selectizeInput('comparison_Y', label = "Select X comparison", choices = NULL)
                                         ),
                                       column(
                                         width = 4.5,
                                         selectizeInput('analysis_method_se', label = "Method", choices = c("DrugZ", "MAGeCK")),
                                         actionButton('bi_submit', 'Generate biplot')
                                       )
                                     ),
                                     fluidRow(
                                       plotlyOutput('bi_plot', height = 900),
                                       DTOutput('bi_table'),
                                       downloadButton('bi_download', "Download")
                                     )
                            )
                          )
                        )
                      ),
             #### Exorcise service ####
             tabPanel("Exorcise online",
                      fluidPage()
                      )
             ),
             add_busy_bar(color = "#FF0000")
  )


server <- function(input, output, session) {
  
  #### Initialise ####
  
  # Experiments table
  output$experiments_table <- renderDT(experiments, filter = "top", server = TRUE)
  experiments_table_proxy <- dataTableProxy("experiments_table")
  
  # Comparisons table
  output$comparisons_table <- renderDT(comparisons, filter = "top", server = TRUE)
  comparisons_table_proxy <- dataTableProxy("comparisons_table")
  
  # Volcano plot table
  output$volcano_table <- renderDT(tibble(), filter = "top", server = TRUE)
  
  # Rank plot table
  output$rank_table <- renderDT(tibble(), filter = "top", server = TRUE)
  
  # Biplot table
  output$bi_table <- renderDT(tibble(), filter = "top", server = TRUE)
  
  # Comparisons text picker
  updateSelectizeInput(session, "comparisons_selected", server = TRUE, choices = comparisons$`Comparison ID`)
  
  # Genes of interest text pickers
  updateSelectizeInput(session, "genes_queried", server = TRUE, choices = genes, selected = c("RNF4", "ATM"))
  updateSelectizeInput(session, "genes_selected", server = TRUE, choices = genes, selected = c("RNF4", "ATM"))
  
  
  
  
  #### React ####
  
  # When experiments are picked from the table
  observeEvent(input$experiments_table_rows_selected, {
    experiments_selected <- experiments$`Experiment ID`[input$experiments_table_rows_selected]
    # Filter the comparisons on the picked experiments
    comparisons_filtered <- comparisons %>% filter(`Experiment ID` %in% experiments_selected)
    output$comparisons_table <- renderDT(comparisons_filtered, filter = "top", server = TRUE)
    # Update the comparisons text picker choices
    updateSelectizeInput(session, "comparisons_selected", choices = comparisons_filtered$`Comparison ID`)
  })
  
  # When the "select all experiments in view" button is pressed
  observeEvent(input$experiments_table_selectAll, {
    # Select all visible rows in the experiments table
    selectRows(experiments_table_proxy, input$experiments_table_rows_all)
  })
  
  # When the "deselect all experiments" button is pressed
  observeEvent(input$experiments_table_selectNone, {
    # Clear the selection in the experiments table
    selectRows(experiments_table_proxy, NULL)
    # Remove the filter on the comparisons table
    output$comparisons_table <- renderDT(comparisons, filter = "top", server = TRUE)
    # Clear the comparisons text picker
    updateSelectizeInput(session, "comparisons_selected", server = TRUE, selected = NA, choices = comparisons$`Comparison ID`)
  })
  
  # When comparisons are picked from the table
  observeEvent(input$comparisons_table_rows_selected, {
    experiments_selected <- experiments$`Experiment ID`[input$experiments_table_rows_selected]
    if(length(experiments_selected) > 0) comparisons_filtered <- comparisons %>% filter(`Experiment ID` %in% experiments_selected)
    else comparisons_filtered <- comparisons
    comparisons_selected <- comparisons_filtered$`Comparison ID`[input$comparisons_table_rows_selected]
    # Update the comparisons text picker selection
    updateSelectizeInput(session, "comparisons_selected", selected = comparisons_selected)
  })
  
  # When comparisons are picked from the text picker
  observeEvent(input$comparisons_selected, {
    experiments_selected <- experiments$`Experiment ID`[input$experiments_table_rows_selected]
    if(length(experiments_selected) > 0) comparisons_filtered <- comparisons %>% filter(`Experiment ID` %in% experiments_selected)
    else comparisons_filtered <- comparisons
    comparisons_table_rows_selected <- which(comparisons_filtered$`Comparison ID` %in% input$comparisons_selected)
    # Select those rows in the comparisons table
    selectRows(comparisons_table_proxy, comparisons_table_rows_selected)
    # Update the individual comparison pickers
    updateSelectizeInput(session, "comparison_volcano", server = TRUE, choices = input$comparisons_selected)
    updateSelectizeInput(session, "comparison_rank", server = TRUE, choices = input$comparisons_selected)
    updateSelectizeInput(session, "comparison_X", server = TRUE, choices = input$comparisons_selected)
    updateSelectizeInput(session, "comparison_Y", server = TRUE, choices = input$comparisons_selected)
  })
  
  # When the "select all comparisons in view" button is pressed
  observeEvent(input$comparisons_table_selectAll, {
    # Select all visible rows in the comparisons table
    selectRows(comparisons_table_proxy, input$comparisons_table_rows_all)
  })
  
  # When the "deselect all experiments" button is pressed
  observeEvent(input$comparisons_table_selectNone, {
    # Clear the selection in the experiments table
    selectRows(comparisons_table_proxy, NULL)
    # Clear the comparisons text picker
    updateSelectizeInput(session, "comparisons_selected", selected = NA)
  })
  
  # When the volcano plot button is pressed
  observeEvent(input$volcano_submit, {
    if (input$comparison_volcano != "" & input$analysis_method_se != "") {
      # Plot a volcano plot
      plotdata_fdr <- sqldf(queryDdrcs(input$comparison_volcano, NULL, input$analysis_method_se, "fdr"), connection = con)
      plotdata_score <- sqldf(queryDdrcs(input$comparison_volcano, NULL, input$analysis_method_se, "score"), connection = con)
      plotdata <- inner_join(plotdata_fdr, plotdata_score, by = "gene")
      names(plotdata) <- c("gene", "fdr", "score")
      plotdata <- plotdata %>%
        filter(!is.na(fdr)) %>%
        filter(!is.na(score)) %>%
        mutate(goi = case_when(gene %in% input$genes_selected ~ T, T ~ F),
               label = case_when(goi ~ gene, T ~ "")) %>%
        arrange(goi)
      p <- ggplot(plotdata, aes(x = score, y = fdr, label = label, colour = goi, gene = gene)) +
        geom_point(aes(alpha = fdr)) +
        geom_text() +
        scale_y_continuous(trans = reverselog_trans()) +
        scale_alpha_continuous(trans = reverselog_trans())
      p <- ggplotly(p)
      output$volcano_table <- renderDT(plotdata, server = TRUE, filter = "top")
      output$volcano_plot <- renderPlotly({ p })
      output$volcano_download <- downloadHandler(
        filename = "export.csv",
        content = function(f) fwrite(plotdata, f, sep = ",")
      )
    }
  })
  
  # When the rank plot button is pressed
  observeEvent(input$rank_submit, {
    if (input$comparison_rank != "" & input$analysis_method_se != "") {
      # Plot a rank plot
      plotdata_fdr <- sqldf(queryDdrcs(input$comparison_rank, NULL, input$analysis_method_se, "fdr"), connection = con)
      plotdata_score <- sqldf(queryDdrcs(input$comparison_rank, NULL, input$analysis_method_se, "score"), connection = con)
      plotdata <- inner_join(plotdata_fdr, plotdata_score, by = "gene")
      names(plotdata) <- c("gene", "fdr", "score")
      plotdata <- plotdata %>%
        filter(!is.na(fdr)) %>%
        filter(!is.na(score)) %>%
        mutate(rank = rank(score),
               goi = case_when(gene %in% input$genes_selected ~ T, T ~ F),
               label = case_when(goi ~ gene, T ~ "")) %>%
        filter(!is.na(fdr)) %>%
        arrange(goi)
      p <- ggplot(plotdata, aes(x = rank, y = score, label = label, colour = goi, gene = gene)) +
        geom_point(aes(alpha = fdr)) +
        geom_text() +
        scale_alpha_continuous(trans = reverselog_trans())
      p <- ggplotly(p)
      output$rank_table <- renderDT(plotdata, server = TRUE, filter = "top")
      output$rank_plot <- renderPlotly({ p })
      output$rank_download <- downloadHandler(
        filename = "export.csv",
        content = function(f) fwrite(plotdata, f, sep = ",")
      )
    }
  })
  
  # When the biplot button is pressed
  observeEvent(input$bi_submit, {
    if (input$comparison_X != "" & input$comparison_Y != "" & input$analysis_method_se != "") {
      # Plot a biplot
      plotdataX_fdr <- sqldf(queryDdrcs(input$comparison_X, NULL, input$analysis_method_se, "fdr"), connection = con)
      plotdataX_score <- sqldf(queryDdrcs(input$comparison_X, NULL, input$analysis_method_se, "score"), connection = con)
      plotdataX <- inner_join(plotdataX_fdr, plotdataX_score, by = "gene")
      names(plotdataX) <- c("gene", "fdrX", "scoreX")
      plotdataY_fdr <- sqldf(queryDdrcs(input$comparison_Y, NULL, input$analysis_method_se, "fdr"), connection = con)
      plotdataY_score <- sqldf(queryDdrcs(input$comparison_Y, NULL, input$analysis_method_se, "score"), connection = con)
      plotdataY <- inner_join(plotdataY_fdr, plotdataY_score, by = "gene")
      names(plotdataY) <- c("gene", "fdrY", "scoreY")
      plotdata <- inner_join(plotdataX, plotdataY, by = "gene")
      plotdata <- plotdata %>%
        filter(!is.na(fdrX)) %>%
        filter(!is.na(scoreX)) %>%
        filter(!is.na(fdrY)) %>%
        filter(!is.na(scoreY)) %>%
        mutate(goi = case_when(gene %in% input$genes_selected ~ T, T ~ F),
               label = case_when(goi ~ gene, T ~ "")) %>%
        arrange(goi)
      p <- ggplot(plotdata, aes(x = scoreX, y = scoreY, label = label, colour = goi, gene = gene)) +
        geom_point() +
        geom_text() +
        geom_hline(yintercept = 0, linetype = "dotted") +
        geom_vline(xintercept = 0, linetype = "dotted") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed")
      p <- ggplotly(p)
      output$bi_table <- renderDT(plotdata, server = TRUE, filter = "top")
      output$bi_plot <- renderPlotly({ p })
      output$bi_download <- downloadHandler(
        filename = "export.csv",
        content = function(f) fwrite(plotdata, f, sep = ",")
      )
    }
  })
  
  # When the gene query button is pressed
  observeEvent(input$genequery_submit, {
    if (!is.null(input$genes_queried)) {
      # Plot gene(s)' scores across screens
      # Display only endpoint timepoints
      endpointScreens <- comparisons %>% filter(`Days (diff)` == `Days (ref)`) %>% select(`Comparison ID`) %>% unique() %>% unlist() %>% setNames(NULL)
      # Obtain screens of interest from GOI below cutoff
      include <- sqldf(queryDdrcs(endpointScreens, input$genes_queried, input$analysis_method_gq, "fdr"), connection = con) %>%
        pivot_longer(-gene, names_to = "Comparison ID", values_to = "fdr") %>%
        filter(fdr <= input$p_cutoff) %>%
        select(`Comparison ID`) %>% unique() %>% unlist() %>% setNames(NULL)
      if (length(include) > 0 ) {
        # For screens of interest, obtain score and fdr
        fdr <- sqldf(queryDdrcs(include, NULL, input$analysis_method_gq, "fdr"), connection = con) %>%
          pivot_longer(-gene, names_to = "Comparison ID", values_to = "fdr")
        score <- sqldf(queryDdrcs(include, NULL, input$analysis_method_gq, "score"), connection = con) %>%
          pivot_longer(-gene, names_to = "Comparison ID", values_to = "score")
        plotdata <- full_join(fdr, score, by = c("gene", "Comparison ID")) %>%
          filter(!is.na(fdr)) %>%
          filter(!is.na(score)) %>%
          left_join(comparisons, by = "Comparison ID") %>%
          mutate(goi = case_when(gene %in% input$genes_queried ~ T, T ~ F))
        p <- ggplot(plotdata, aes(y = score, x = `Comparison ID`)) +
          geom_violin(na.rm = T) +
          geom_point(data = plotdata %>% filter(goi), aes(colour = gene, fdr = fdr, ExperimentID = `Experiment ID`, TreatmentDiff = `Treatment (diff)`, DoseDiff = `Dose (diff)`, DaysDiff = `Days (diff)`, KnockoutDiff = `Knockout (diff)`, TreatmentRef = `Treatment (ref)`, DoseRef = `Dose (ref)`, DaysRef = `Days (ref)`, KnockoutRef = `Knockout (ref)`, CellLine = `Cell line`, Library = `Library`, Exorcised = `Exorcised`)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        displaydata <- plotdata %>%
          filter(goi)
        p <- ggplotly(p)
        output$genequery_plot <- renderPlotly({ p })
        output$genequery_table <- renderDT(displaydata, server = TRUE, filter = "top")
      } else {
        plotdata <- tibble(x = 0, y = 0, label = "No hits found.")
        output$genequery_plot <- renderPlotly(ggplotly(ggplot(plotdata, aes(x = x, y = y, label = label)) + geom_text()))
        output$genequery_table <- renderDT(plotdata, server = TRUE)
        output$genequery_download <- downloadHandler(
          filename = "export.csv",
          content = function(f) fwrite(plotdata, f, sep = ",")
        )
      }
    }
  })
  
}

#path_dataset <- "/srv/shiny-server/crave/www/app_data/"
path_dataset <- "~/bio/Projects/operation/ddrcs/www/ddrcs-internal-dev"


# I/O
file_experiments <- paste0(path_dataset, "/experiments_metadata.csv.gz")
file_comparisons <- paste0(path_dataset, "/comparisons_metadata.csv.gz")
file_data        <- paste0(path_dataset, "/ddrcs.db"                   )
con <- dbConnect(SQLite(), file_data)

experiments <- fread(file_experiments)
comparisons <- fread(file_comparisons)

# Wrangle experiments metadata
experiments <- experiments %>%
  select(`Experiment ID`) %>%
  left_join(comparisons, by = "Experiment ID") %>%
  group_by(`Experiment ID`) %>%
  summarise(Treatments = paste(unique(`Treatment`), collapse = ", "),
            Knockouts = paste(unique(`KO`), collapse = ", "),
            `Cell lines` = paste(unique(`Cell line`), collapse = ", "),
            Libraries = paste(unique(Library), collapse = ", "),
            Exorcised = paste(unique(Exorcised), collapse = ", "))

# Wrangle comparisons metadata
comparisons <- comparisons %>%
  transmute(`Experiment ID`,
            `Comparison ID`,
            `Treatment (diff)` = Treatment,
            `Dose (diff)` = Dose,
            `Days (diff)` = `Days grown`,
            `Knockout (diff)` = KO,
            `Treatment (ref)` = Treatment.ref,
            `Dose (ref)` = Dose.ref,
            `Days (ref)` = `Days grown.ref`,
            `Knockout (ref)` = KO.ref,
            `Cell line`,
            Library,
            Exorcised)

# Obtain genes list
genes <- sqldf("SELECT [gene] from drz_fdr", connection = con) %>% unlist() %>% setNames(NULL)
genes <- genes[!grepl("Non-targeting", genes)] # exclude Non-targeting
genes <- genes[!grepl("[ATCG]{19}", genes)] # exclude sequences
genes <- genes[genes != "X"] # exclude "X"

# Start the server
shinyApp(ui, server)

