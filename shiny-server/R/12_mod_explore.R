#### Explore module ####
# Owns the Explore tab: gene selection, experiment and comparison selection, and
# the five single-screen analyses. All input ids are namespaced, so the ids
# themselves are short and descriptive instead of being globally unique strings.

# The comparison search form: one row per field.
COMPARISON_SEARCH_FIELDS <- list(
  list(id = "citation",       label = "Citation",           key = "Citation",           col = "Citation"),
  list(id = "kind",           label = "Kind",               key = "Kind",               col = "Kind"),
  list(id = "library",        label = "Library",            key = "Library",            col = "Library"),
  list(id = "timepoint",      label = "Timepoint",          key = "Timepoint",          col = "Timepoint"),
  list(id = "organism",       label = "Species",            key = "Organism",           col = "Organism"),
  list(id = "source",         label = "Source",             key = "Source",             col = "Source"),
  list(id = "endpoint",       label = "Endpoint",           key = "Endpoint",           col = "Endpoint"),
  list(id = "daysgrown_diff", label = "Days grown (diff)",  key = "Days grown (diff)",  col = "Days grown (diff)"),
  list(id = "daysgrown_ref",  label = "Days grown (ref)",   key = "Days grown (ref)",   col = "Days grown (ref)"),
  list(id = "treatment_diff", label = "Treatment (diff)",   key = "Treatment (diff)",   col = "Treatment (diff)"),
  list(id = "treatment_ref",  label = "Treatment (ref)",    key = "Treatment (ref)",    col = "Treatment (ref)"),
  list(id = "dose_diff",      label = "Dose (diff)",        key = "Dose (diff)",        col = "Dose (diff)"),
  list(id = "dose_ref",       label = "Dose (ref)",         key = "Dose (ref)",         col = "Dose (ref)"),
  list(id = "knockout_diff",  label = "Knockout (diff)",    key = "Knockout (diff)",    col = "Knockout (diff)"),
  list(id = "knockout_ref",   label = "Knockout (ref)",     key = "Knockout (ref)",     col = "Knockout (ref)"),
  list(id = "cellline_diff",  label = "Cell line (diff)",   key = "Cell line (diff)",   col = "Cell line (diff)"),
  list(id = "cellline_ref",   label = "Cell line (ref)",    key = "Cell line (ref)",    col = "Cell line (ref)")
)

ANY_CHOICE <- "(any)"

# Columns hidden from the comparisons table: internal identifiers.
COMPARISONS_TABLE_HIDDEN <- c("Comparison ID", "Experiment ID", "FriendlyID")

#' Values of a comparison-table column as plain character, tags stripped for
#' Citation.
comparisonColumn <- function(cmp, col) {
  if (identical(col, "Citation")) stripHtml(as.character(cmp$Citation))
  else as.character(cmp[[col]])
}

#### UI ####

exploreSidebarUI <- function(ns) {
  sidebarPanel(
    width = 3,
    textAreaInput(ns("genes_text"), "Type genes of interest", rows = 3,
                  placeholder = "e.g. TP53, BRCA1, UBE2K"),
    actionLink(ns("genes_validate"), "Validate"), " · ",
    actionLink(ns("genes_essentials"), "Core essentials"), " · ",
    actionLink(ns("genes_copy"), "Copy GOIs from Correlate"),
    tags$br(), tags$br(),
    selectizeInput(ns("genes"), "Search genes of interest", choices = NULL, multiple = TRUE),
    actionButton(ns("genes_clear"), "Clear genes", icon = icon("ban"), class = "btn-danger"),
    tags$br(), tags$br(),
    selectizeInput(ns("contrasts"), "Selected comparisons", choices = NULL, multiple = TRUE),
    actionButton(ns("contrasts_clear"), "Clear comparisons", icon = icon("ban"), class = "btn-danger"),
    tags$br(), tags$br(),
    customiseUI(ns("customise")),
    tags$br(), tags$br(), tags$hr(),
    tags$p(tags$strong("Help")),
    tags$p("Use", tags$strong("Select experiment"), "and", tags$strong("Select comparison"),
           "tabs to choose your screen(s) of interest."),
    tags$p(tags$strong("Volcano plot"), "and", tags$strong("Rank plot"),
           "show the distribution of CRISPR scores in a single screen."),
    tags$p(tags$strong("Biplot"), "compares the distributions of CRISPR scores of two screens."),
    tags$p(tags$strong("Overlap"), "shows common hits between two or more screens and",
           "calculates the statistical significance of the overlap between pairs of",
           "screens using a hypergeometric test."),
    tags$p(tags$strong("ROC"), "shows how well hits from a reference screen are",
           "captured by other screens using precision-recall curves.")
  )
}

exploreSelectUI <- function(ns) {
  searchInput <- function(f) {
    selectizeInput(ns(paste0("sel_", f$id)), f$label, choices = ANY_CHOICE,
                   selected = ANY_CHOICE, multiple = TRUE)
  }
  experimentFields <- Filter(function(f) f$id %in%
    c("citation", "kind", "library", "timepoint", "organism", "source", "endpoint"),
    COMPARISON_SEARCH_FIELDS)
  levelFields <- Filter(function(f) !(f$id %in%
    c("citation", "kind", "library", "timepoint", "organism", "source", "endpoint")),
    COMPARISON_SEARCH_FIELDS)

  list(
    experiment = tabPanel(
      "Select experiment", icon = icon("flask"),
      tags$br(),
      tags$h4("Select experiment"),
      tags$p("Use the table to choose papers or NGS submissions of interest. The",
             tags$strong("Select comparison"),
             "tab will be populated with screens from only the selected experiments."),
      tags$p("Leave blank to show screens from the whole database."),
      tags$br(), tags$hr(), tags$br(),
      actionButton(ns("experiments_all"), "Select all", icon = icon("check-double"), class = "btn-success"),
      actionButton(ns("experiments_none"), "Select none", icon = icon("ban"), class = "btn-danger"),
      tags$br(), tags$br(),
      DTOutput(ns("experiments_table"))
    ),
    comparison = tabPanel(
      "Select comparison", icon = icon("vials"),
      tags$br(),
      tags$h4("Select comparison"),
      tags$p("Use the form to filter screens by experimental design or conditions on",
             "the denominator (reference) or numerator (differential)."),
      tags$p("Use the table to choose screens of interest. Click the",
             "\"Select comparisons\" button to choose all of them."),
      tags$p("As you fill the form, the contents of the table updates."),
      tags$br(), tags$hr(), tags$br(),
      tags$p(tags$strong("Comparison search")),
      tags$p("Selected experiments:"),
      textOutput(ns("selected_experiments")),
      tags$br(),
      tags$p(tags$strong("Experiment")),
      do.call(inputGrid, lapply(experimentFields, searchInput)),
      tags$hr(),
      tags$table(
        style = "width: 100%",
        tags$tr(
          style = "vertical-align: top",
          tags$td(style = "padding: 6px; width: 50%", tags$p(tags$strong("Numerator levels"))),
          tags$td(style = "padding: 6px; width: 50%", tags$p(tags$strong("Denominator levels")))
        )
      ),
      do.call(inputGrid, lapply(levelFields, searchInput)),
      inputGrid(
        tagList(
          actionButton(ns("sel_select"), "Select comparisons",
                       icon = icon("circle-check"), class = "btn-success"),
          " · ",
          actionLink(ns("sel_clear"), "Clear form")
        ),
        textOutput(ns("selected_comparisons"))
      ),
      tags$br(), tags$hr(), tags$br(),
      tags$p(tags$strong("All comparisons")),
      actionButton(ns("comparisons_all"), "Select all", icon = icon("check-double"), class = "btn-success"),
      actionButton(ns("comparisons_none"), "Select none", icon = icon("ban"), class = "btn-danger"),
      tags$br(), tags$br(),
      DTOutput(ns("comparisons_table"))
    )
  )
}

exploreUI <- function(id) {
  ns <- NS(id)
  sel <- exploreSelectUI(ns)

  tabPanel(
    "Explore", icon = icon("microscope"),
    sidebarLayout(
      exploreSidebarUI(ns),
      mainPanel(
        tabsetPanel(
          sel$experiment,
          sel$comparison,

          tabPanel(
            "Volcano plot", icon = icon("volcano"),
            plotPanelUI(
              ns("volcano"), "Volcano plot",
              c("Generate a volcano plot of CRISPR score against FDR for one screen.",
                methodLinks(chronos = FALSE)),
              controls = inputGrid(
                selectizeInput(ns("volcano_contrast"), "Select a screen", choices = NULL,
                               options = list(placeholder = "Choose a screen")),
                selectizeInput(ns("volcano_method"), "Method",
                               choices = c("DrugZ", "MAGeCK", "Manual"))
              ),
              submitLabel = "Volcano plot", submitIcon = "volcano"
            )
          ),

          tabPanel(
            "Rank plot", icon = icon("ranking-star"),
            plotPanelUI(
              ns("rank"), "Rank plot",
              c("Generate a rank plot of CRISPR score against FDR for one screen.",
                methodLinks()),
              controls = inputGrid(
                selectizeInput(ns("rank_contrast"), "Select a screen", choices = NULL,
                               options = list(placeholder = "Choose a screen")),
                selectizeInput(ns("rank_method"), "Method", choices = CRAVE_METHOD_NAMES)
              ),
              submitLabel = "Rank plot", submitIcon = "ranking-star"
            )
          ),

          tabPanel(
            "Biplot", icon = icon("compass-drafting"),
            plotPanelUI(
              ns("biplot"), "Biplot",
              c("Plot CRISPR scores for two screens against each other.",
                paste(methodLinks(), "You can choose a different method for each screen.")),
              controls = inputGrid(
                tagList(
                  selectizeInput(ns("biplot_x"), "Select two screens", choices = NULL,
                                 options = list(placeholder = "Choose a screen for the x-axis")),
                  selectizeInput(ns("biplot_y"), NULL, choices = NULL,
                                 options = list(placeholder = "Choose a screen for the y-axis"))
                ),
                tagList(
                  selectizeInput(ns("biplot_method_x"), "Methods", choices = CRAVE_METHOD_NAMES),
                  selectizeInput(ns("biplot_method_y"), NULL, choices = CRAVE_METHOD_NAMES)
                )
              ),
              submitLabel = "Biplot", submitIcon = "compass-drafting"
            )
          ),

          tabPanel(
            "Overlap", icon = icon("circle-half-stroke"),
            plotPanelUI(
              ns("overlap"), "Overlap",
              c("Show common hits between two or more screens.",
                "Use the controls to set the threshold for a hit.",
                methodLinks()),
              controls = tagList(
                inputGrid(
                  selectizeInput(ns("overlap_contrasts"), "Select screens", choices = NULL, multiple = TRUE),
                  selectizeInput(ns("overlap_method"), "Method", choices = CRAVE_METHOD_NAMES)
                ),
                inputGrid(
                  selectizeInput(ns("overlap_stat"), "Statistic", choices = c("Score", "FDR", "Rank")),
                  textInput(ns("overlap_value"), "Better than", value = 3),
                  selectizeInput(ns("overlap_direction"), "Direction",
                                 choices = c("Hypersensitivity hits", "Suppressing hits")),
                  selectizeInput(ns("overlap_style"), "Type of plot",
                                 choices = c("Venn diagram", "Upset plot")),
                  ncol = 2
                )
              ),
              submitLabel = "Overlap", submitIcon = "circle-half-stroke"
            )
          ),

          tabPanel(
            "ROC", icon = icon("chart-area"),
            plotPanelUI(
              ns("roc"), "ROC",
              c("Draw a ROC curve showing the precision-recall of hits in a reference screen by query screens.",
                "Use the controls to set the threshold for a hit in the reference screen.",
                methodLinks()),
              controls = tagList(
                inputGrid(
                  selectizeInput(ns("roc_base"), "Select reference screen", choices = NULL),
                  selectizeInput(ns("roc_method"), "Method", choices = CRAVE_METHOD_NAMES)
                ),
                inputGrid(
                  selectizeInput(ns("roc_stat"), "Statistic", choices = c("Score", "FDR", "Rank")),
                  textInput(ns("roc_value"), "Better than", value = 3),
                  selectizeInput(ns("roc_direction"), "Direction",
                                 choices = c("Hypersensitivity hits", "Suppressing hits", "Both")),
                  selectizeInput(ns("roc_head"), "Select query screens", choices = NULL, multiple = TRUE),
                  ncol = 2
                )
              ),
              submitLabel = "ROC", submitIcon = "chart-area"
            )
          )
        )
      )
    )
  )
}

#### Server ####

exploreServer <- function(id, data, bus) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    customise <- customiseServer("customise", EXPLORE_CUSTOMISE_FIELDS, EXPLORE_DEFAULTS)

    # Genes available in the currently selected screens, or all genes.
    availableGenes <- reactiveVal(character(0))
    # Guard against the table and the selectize input updating each other forever.
    syncing <- reactiveVal(FALSE)

    experimentsProxy <- dataTableProxy("experiments_table")
    comparisonsProxy <- dataTableProxy("comparisons_table")

    # Overrides the "screens found" readout with a warning when a selection is
    # refused, cleared on the next search change.
    selectionNotice <- reactiveVal(NULL)

    #' Reset every search field to "(any)" with choices drawn from `cmp`.
    resetSearchForm <- function(cmp) {
      for (f in COMPARISON_SEARCH_FIELDS) {
        updateSelectizeInput(
          session, paste0("sel_", f$id), server = FALSE, selected = ANY_CHOICE,
          choices = c(ANY_CHOICE, sort(unique(comparisonColumn(cmp, f$col))))
        )
      }
    }

    #### Initialisation, also re-run when the data are refreshed ####
    observeEvent(data(), {
      d <- data()
      availableGenes(d$genes)
      updateSelectizeInput(session, "contrasts", server = TRUE,
                           choices = d$choices$FriendlyID,
                           options = list(maxOptions = 1000))
      updateSelectizeInput(session, "genes", server = TRUE, choices = d$genes,
                           options = list(maxOptions = CRAVE_MAX_SELECTIZE_SERVER))
      resetSearchForm(d$comparisons)
    })

    #### Gene selection ####
    observeEvent(input$genes_validate, {
      txt <- input$genes_text %||% ""
      if (nchar(txt) < 3) return()
      pool  <- availableGenes()
      typed <- unique(splitTokens(txt))
      found <- unique(pool[toupper(pool) %in% toupper(typed)])
      missed <- setdiff(toupper(typed), toupper(found))
      if (length(missed)) {
        showNotification(
          paste0(length(missed), " symbol(s) not found in the selected screens: ",
                 paste(utils::head(missed, 10), collapse = ", "),
                 if (length(missed) > 10) ", ..." else ""),
          type = "warning", duration = 6
        )
      }
      updateSelectizeInput(session, "genes",
                           selected = unique(c(input$genes, found)),
                           choices = pool,
                           server = length(found) <= CRAVE_MAX_SELECTIZE_SERVER)
    })

    observeEvent(input$genes_essentials, {
      updateTextAreaInput(session, "genes_text",
                          value = paste(CRAVE_ESSENTIALS, collapse = ", "))
    })

    observeEvent(input$genes_clear, {
      updateSelectizeInput(session, "genes", selected = character(0),
                           choices = availableGenes(), server = TRUE)
    })

    observeEvent(input$genes_copy, {
      updateSelectizeInput(session, "genes", selected = bus$correlateGenes,
                           choices = availableGenes(), server = TRUE)
    })

    # Publish this tab's selections for the Correlate tab to copy.
    observe({
      bus$exploreGenes <- input$genes
      bus$exploreContrasts <- input$contrasts
    })

    # Genes brushed on a plot are added to the selection.
    observeEvent(bus$brush, {
      b <- bus$brush
      if (is.null(b) || !identical(b$target, "Explore") || length(b$genes) == 0) return()
      updateSelectizeInput(session, "genes",
                           selected = unique(c(input$genes, b$genes)),
                           choices = availableGenes(), server = TRUE)
    }, ignoreInit = TRUE)

    #### Experiments table ####
    output$experiments_table <- renderDT(
      data()$experiments %>% select(-`Experiment ID`),
      filter = "top", server = TRUE, escape = FALSE,
      options = list(pageLength = 100, scrollX = TRUE)
    )

    observeEvent(input$experiments_all,  selectRows(experimentsProxy, input$experiments_table_rows_all))
    observeEvent(input$experiments_none, selectRows(experimentsProxy, NULL))

    #' Experiments currently ticked, or all of them if none are.
    selectedExperiments <- reactive({
      d <- data()
      rows <- input$experiments_table_rows_selected
      if (length(rows) == 0) return(d$experiments)
      d$experiments[rows, , drop = FALSE]
    })

    #' Comparisons belonging to the ticked experiments.
    scopedComparisons <- reactive({
      d <- data()
      exps <- selectedExperiments()$`Experiment ID`
      d$comparisons %>% filter(`Experiment ID` %in% exps)
    })

    #### Comparisons table ####
    output$comparisons_table <- renderDT(
      scopedComparisons() %>% select(-all_of(COMPARISONS_TABLE_HIDDEN)),
      filter = "top", server = TRUE, escape = FALSE,
      options = list(pageLength = 100, scrollX = TRUE)
    )

    observeEvent(input$comparisons_all,  selectRows(comparisonsProxy, input$comparisons_table_rows_all))
    observeEvent(input$comparisons_none, {
      selectRows(comparisonsProxy, NULL)
      updateSelectizeInput(session, "contrasts", selected = character(0))
    })
    observeEvent(input$contrasts_clear, {
      selectRows(comparisonsProxy, NULL)
      updateSelectizeInput(session, "contrasts", selected = character(0))
    })

    # Table selection -> selectize. The one-shot `syncing` flag stops the two
    # directions from updating each other indefinitely.
    observeEvent(input$comparisons_table_rows_selected, ignoreNULL = FALSE, {
      if (isTRUE(syncing())) { syncing(FALSE); return() }
      picked <- scopedComparisons()$FriendlyID[input$comparisons_table_rows_selected]
      if (setequal(picked, input$contrasts %||% character(0))) return()
      syncing(TRUE)
      updateSelectizeInput(session, "contrasts", selected = picked)
    })

    # selectize -> table, plus the per-analysis screen pickers.
    observeEvent(input$contrasts, ignoreNULL = FALSE, {
      chosen <- input$contrasts %||% character(0)

      if (!isTRUE(syncing())) {
        rows <- which(scopedComparisons()$FriendlyID %in% chosen)
        current <- input$comparisons_table_rows_selected %||% integer(0)
        if (!setequal(rows, current)) {
          syncing(TRUE)
          selectRows(comparisonsProxy, rows)
        }
      } else {
        syncing(FALSE)
      }

      # Narrow the gene picker to genes measured in the selected screens.
      d <- data()
      ids <- comparisonIdsFor(d, chosen)
      pool <- if (length(ids) > 0 && length(ids) < CRAVE_MAX_SCREENS_PLOT) {
        withProgress(message = "Finding genes in the selected screens...", value = NULL, {
          intersect(d$genes, fetchGenesForComparisons(d, ids))
        })
      } else {
        d$genes
      }
      availableGenes(pool)
      keep <- intersect(input$genes %||% character(0), pool)
      updateSelectizeInput(session, "genes", server = TRUE, choices = pool,
                           selected = keep,
                           options = list(maxOptions = CRAVE_MAX_SELECTIZE_SERVER))

      # Per-analysis screen pickers all draw from the current selection.
      #
      # Each picker keeps whatever it already had, if that is still among the
      # selected screens, and otherwise falls back to a sensible default. Passing no
      # `selected` here would reset every analysis tab's screen choice whenever the
      # comparison selection changed.
      keepOr <- function(pid, fallback) {
        current <- intersect(isolate(input[[pid]]) %||% character(0), chosen)
        if (length(current)) current else fallback
      }
      updateSelectizeInput(session, "volcano_contrast", choices = chosen, server = TRUE,
                           selected = keepOr("volcano_contrast", chosen[1]))
      updateSelectizeInput(session, "rank_contrast", choices = chosen, server = TRUE,
                           selected = keepOr("rank_contrast", chosen[1]))
      updateSelectizeInput(session, "biplot_x", choices = chosen, server = TRUE,
                           selected = keepOr("biplot_x", chosen[1]))
      updateSelectizeInput(session, "biplot_y", choices = chosen, server = TRUE,
                           selected = keepOr("biplot_y",
                                             if (length(chosen) > 1) chosen[2] else chosen[1]))
      updateSelectizeInput(session, "overlap_contrasts", choices = chosen, server = TRUE,
                           selected = keepOr("overlap_contrasts",
                                             utils::head(chosen, CRAVE_MAX_VENN_SETS)))
      updateSelectizeInput(session, "roc_base", choices = chosen, server = TRUE,
                           selected = keepOr("roc_base", chosen[1]))
      updateSelectizeInput(session, "roc_head", choices = chosen, server = TRUE,
                           selected = keepOr("roc_head",
                                             if (length(chosen) > 1) chosen[2:min(length(chosen), 6)]
                                             else character(0)))
    })

    #### Comparison search form ####

    #' Current value of a search field, with the "(any)" sentinel removed.
    searchValue <- function(f) {
      want <- input[[paste0("sel_", f$id)]]
      want[want != ANY_CHOICE]
    }

    #' The comparisons that satisfy the search form, within the scoped set.
    searchedComparisons <- reactive({
      cmp <- scopedComparisons()
      for (f in COMPARISON_SEARCH_FIELDS) {
        want <- searchValue(f)
        if (length(want) == 0) next
        cmp <- cmp[comparisonColumn(cmp, f$col) %in% want, , drop = FALSE]
      }
      cmp
    })

    # Keep the form's choices, its "(any)" sentinel and the table's column filters
    # in step. Debounced, because one click would otherwise fan out into a
    # selectize update for every field in the form.
    searchTrigger <- debounce(
      reactive({
        list(
          fields = lapply(COMPARISON_SEARCH_FIELDS, function(f) input[[paste0("sel_", f$id)]]),
          rows   = input$experiments_table_rows_selected
        )
      }),
      250
    )

    observeEvent(searchTrigger(), ignoreNULL = FALSE, {
      selectionNotice(NULL)
      cmp <- searchedComparisons()
      tableCols <- setdiff(names(data()$comparisons), COMPARISONS_TABLE_HIDDEN)
      columnFilters <- rep("", length(tableCols) + 1L)

      for (f in COMPARISON_SEARCH_FIELDS) {
        inputId <- paste0("sel_", f$id)
        current <- input[[inputId]]

        # Drop the sentinel once something real is picked, and restore it when
        # the field is emptied.
        newSelected <- if (is.null(current) || length(current) == 0) {
          ANY_CHOICE
        } else if (length(current) > 1 && current[1] == ANY_CHOICE) {
          current[-1]
        } else if (current[length(current)] == ANY_CHOICE) {
          ANY_CHOICE
        } else {
          current
        }

        updateSelectizeInput(
          session, inputId, server = FALSE, selected = newSelected,
          choices = c(ANY_CHOICE, sort(unique(comparisonColumn(cmp, f$col))))
        )

        # Mirror the selection into the DT column filter. +1 for the row-name
        # column that DT prepends.
        colIdx <- match(f$col, tableCols)
        if (!is.na(colIdx)) {
          picked <- newSelected[newSelected != ANY_CHOICE]
          if (length(picked)) {
            columnFilters[colIdx + 1L] <-
              paste0("[", paste0("\"", picked, "\"", collapse = ","), "]")
          }
        }
      }

      updateSearch(comparisonsProxy,
                   keywords = list(global = NULL, columns = columnFilters))
    })

    output$selected_experiments <- renderText(
      summariseItems(stripHtml(as.character(selectedExperiments()$Citation)))
    )
    output$selected_comparisons <- renderText(
      selectionNotice() %||% summariseItems(sort(searchedComparisons()$FriendlyID))
    )

    observeEvent(input$sel_clear, resetSearchForm(scopedComparisons()))

    observeEvent(input$sel_select, {
      found <- searchedComparisons()$FriendlyID
      if (length(found) >= 100) {
        selectionNotice(paste0(
          "Won't select ", length(found), " comparisons. Please refine your search."))
        return()
      }
      updateSelectizeInput(session, "contrasts",
                           selected = unique(c(input$contrasts, found)))
    })

    #### Analyses ####
    opts <- reactive(customise$values())

    plotPanelServer("volcano", busyMessage = "Drawing volcano plot...",
      compute = function() {
        plotVolcano(data(), input$volcano_contrast, input$volcano_method,
                    input$genes, opts())
      },
      onResult = function(res) publishBrush(bus, "Explore", res$brush))

    plotPanelServer("rank", busyMessage = "Drawing rank plot...",
      compute = function() {
        plotRank(data(), input$rank_contrast, input$rank_method, input$genes, opts())
      },
      onResult = function(res) publishBrush(bus, "Explore", res$brush))

    plotPanelServer("biplot", busyMessage = "Drawing biplot...",
      compute = function() {
        plotBiplot(data(), input$biplot_x, input$biplot_y,
                   input$biplot_method_x, input$biplot_method_y, input$genes, opts())
      },
      onResult = function(res) publishBrush(bus, "Explore", res$brush))

    plotPanelServer("overlap", busyMessage = "Computing overlaps...",
      compute = function() {
        plotOverlap(data(), input$overlap_contrasts, input$overlap_method,
                    input$overlap_stat, input$overlap_value, input$overlap_direction,
                    input$overlap_style)
      })

    plotPanelServer("roc", busyMessage = "Computing precision-recall...",
      compute = function() {
        plotRoc(data(), input$roc_base, input$roc_head, input$roc_method,
                input$roc_stat, input$roc_value, input$roc_direction)
      })

    #### Save/Load ####
    stateIds <- c(
      "genes_text", "genes", "contrasts",
      paste0("sel_", vapply(COMPARISON_SEARCH_FIELDS, `[[`, character(1), "id")),
      "volcano_contrast", "volcano_method", "rank_contrast", "rank_method",
      "biplot_x", "biplot_y", "biplot_method_x", "biplot_method_y",
      "overlap_contrasts", "overlap_method", "overlap_stat", "overlap_value",
      "overlap_direction", "overlap_style",
      "roc_base", "roc_method", "roc_stat", "roc_value", "roc_direction", "roc_head"
    )
    textStateIds <- c("genes_text", "overlap_value", "roc_value")

    list(
      state = function() {
        st <- captureInputs(input, stateIds)
        st$`_customise` <- customise$values()
        st$`_experiment_rows` <- input$experiments_table_rows_selected
        st
      },
      restore = function(st) {
        d <- data()
        # The screen pickers were created with server = TRUE, so their choices have
        # to accompany the restored selection. The per-analysis pickers are keyed on
        # the restored contrast list rather than on the current one, so the restore
        # does not depend on which observer happens to run first.
        chosen <- st$contrasts %||% character(0)
        serverChoices <- list(
          genes             = d$genes,
          contrasts         = d$choices$FriendlyID,
          volcano_contrast  = chosen,
          rank_contrast     = chosen,
          biplot_x          = chosen,
          biplot_y          = chosen,
          overlap_contrasts = chosen,
          roc_base          = chosen,
          roc_head          = chosen
        )
        restoreInputs(session, st, stateIds, textStateIds,
                      serverChoices = serverChoices)
        customise$restore(st$`_customise`)
        selectRows(experimentsProxy, st$`_experiment_rows`)
        invisible(NULL)
      }
    )
  })
}

#' Publish brushable coordinates for the root session's brush handler.
publishBrush <- function(bus, target, brush) {
  if (is.null(brush) || nrow(brush) == 0) return(invisible(NULL))
  bus$brushTable <- list(target = target, table = brush)
  invisible(NULL)
}
