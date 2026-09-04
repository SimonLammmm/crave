#### Exorcise module ####

exorciseUI <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Exorcise", icon = icon("ghost"),
    sidebarLayout(
      sidebarPanel(
        width = 4,
        tags$h4(tags$strong("Sequences")),
        fileInput(ns("upload"), "Upload a CRISPick sgRNA designs file",
                  multiple = FALSE, accept = c(".txt", ".tsv", ".csv", "text/plain")),
        tags$p("or"),
        textAreaInput(ns("seq"), "Enter sequence(s), separated by commas or spaces",
                      placeholder = "e.g. ATCGATCGATCGATCGATCG, GTACGTACGTACGTACGTAC", rows = 4),
        textAreaInput(ns("orig"),
                      "Optional: existing annotation(s), one per sequence, separated by commas or spaces",
                      placeholder = "e.g. TP53, negative_control", rows = 4),
        textInput(ns("control_string"),
                  "Optional: existing annotation(s) that indicate controls",
                  placeholder = "e.g. negative_control"),
        tags$hr(),
        tags$h4(tags$strong("Reannotation policy")),
        textInput(ns("pam"), "PAM sequence to append", value = "NGG", placeholder = "e.g. NGG"),
        selectizeInput(ns("genome"), "Select reference genome", choices = names(EXORCISE_GENOMES)),
        selectizeInput(ns("mode"), "Select CRISPR chemistry", choices = EXORCISE_MODE_CHOICES),
        textInput(ns("mode_advanced"), "Advanced: custom chemistry mode string",
                  placeholder = "e.g. beCA0213"),
        actionButton(ns("example"), "Example", icon = icon("fire")),
        actionButton(ns("reset"), "Reset form", icon = icon("arrows-rotate"), class = "btn-danger"),
        tags$br(), tags$br(),
        downloadButton(ns("download_input"), "Download input", icon = icon("download")),
        tags$hr(),
        tags$p(tags$strong("Help")),
        tags$p("Choose", tags$strong("sequences"), "to Exorcise by uploading a CRISPick",
               "sgRNA designs file or by typing sequences into the box."),
        tags$p("Sequences must be in [ACGT] and must not contain long runs of the same",
               "nucleotide. The length of the sequence plus protospacer-adjacent motif",
               "(PAM) must be at least 21 excluding any Ns."),
        tags$p("Then choose the", tags$strong("reannotation policy"), "according to the",
               "genome to align to, the PAM to append, and the CRISPR chemistry."),
        tags$p("The PAM can be in [ACGTN] and will be appended to each sequence during",
               "sequence search."),
        tags$p("Knockout chemistry annotates guides with features occurring at the Cas9",
               "cut site. Interference and activation chemistry finds features up to 500",
               "nucleotides away in linear distance. Cytosine base editor chemistry edits",
               "C to G on both strands and assumes a base editing window of [2, 8].",
               "Adenine base editor chemistry edits A to T on both strands and assumes a",
               "base editing window of [4, 9]."),
        tags$p("Custom CRISPR chemistry can be a number indicating the number of",
               "nucleotides upstream and downstream of the alignment of the guide + PAM",
               "to look for features; or a string in the format `beNM0123`, where `N` and",
               "`M` indicate the base edited from and to, both in [ACTG], and `01` and",
               "`23` indicate the start and end of the base editing window, a one-based",
               "fully-closed range counting from the PAM-distal end."),
        tags$p("Examples: `1000`: CRISPRa/i affecting genes up to 1 kb upstream or",
               "downstream in linear distance. `beCA0110`: base editing genomic C to",
               "edited A within the base editing window [1, 10].")
      ),
      mainPanel(
        tags$h4("Exorcise"),
        tags$p("Reannotate CRISPR guides by aligning to other reference genomes using the",
               extLink(DOI_EXORCISE, "Exorcise"), "algorithm."),
        tags$p("This tool uses the implementation maintained at the",
               extLinkTight(URL_EXORCISE_REPO, "Exorcise GitHub repo"), "."),
        tags$br(),
        actionButton(ns("submit"), "Exorcise", icon = icon("ghost"), class = "btn-success"),
        tags$br(), tags$br(),
        textOutput(ns("console")),
        tags$br(), tags$br(),
        downloadButton(ns("download_result"), "Download result", icon = icon("download"))
      )
    )
  )
}

exorciseServer <- function(id, config) {
  moduleServer(id, function(input, output, session) {

    observeEvent(input$example, {
      updateTextAreaInput(session, "seq", value = paste(
        "AAGGAGCCAACATAACAGAT", "AGAAACCTACAACTCATGGA", "GGCTCAGGGTTACCGAAGAG",
        "TGGTACTGATTATGGCACTC", "CAACGTCGCGAACGTCGTAT", "ACCACCGTTCGTACCGGTCG",
        "AGACCGTCGTCGATCGATAC", "GTCGTACGGATTCGCGCGTA", sep = ", "))
      updateTextAreaInput(session, "orig", value = paste(
        "BRCA1", "BRCA1", "BRCA1", "BRCA1",
        "non-targeting", "non-targeting", "non-targeting", "non-targeting", sep = ", "))
      updateTextInput(session, "pam", value = "NGG")
      updateSelectizeInput(session, "genome", selected = names(EXORCISE_GENOMES)[1])
      updateSelectizeInput(session, "mode", selected = "Knockout")
      updateTextInput(session, "mode_advanced", value = "")
      updateTextInput(session, "control_string", value = "non-targeting")
    })

    observeEvent(input$reset, {
      updateTextAreaInput(session, "seq", value = "")
      updateTextAreaInput(session, "orig", value = "")
      updateTextInput(session, "pam", value = "NGG")
      updateTextInput(session, "mode_advanced", value = "")
      updateTextInput(session, "control_string", value = "")
    })

    result <- eventReactive(input$submit, {
      withProgress(message = "Running Exorcise, this can take a minute...", value = NULL, {
        runExorcise(
          params = list(
            upload         = input$upload,
            seq            = input$seq,
            orig           = input$orig,
            control_string = input$control_string,
            pam            = input$pam,
            genome         = input$genome,
            mode           = input$mode,
            mode_advanced  = input$mode_advanced
          ),
          config = config
        )
      })
    }, ignoreInit = TRUE)

    output$console <- renderText(result()$message)

    output$download_input <- downloadHandler(
      filename = function() paste0("exorcise-input_", result()$taskId %||% "input", ".tsv.gz"),
      content  = function(f) data.table::fwrite(result()$input %||% data.table(), f, sep = "\t")
    )

    output$download_result <- downloadHandler(
      filename = function() paste0("exorcise-output_", result()$taskId %||% "result", ".tsv.gz"),
      content  = function(f) data.table::fwrite(result()$result %||% data.table(), f, sep = "\t")
    )

    stateIds <- c("seq", "orig", "control_string", "pam", "genome", "mode", "mode_advanced")
    textStateIds <- c("seq", "orig", "control_string", "pam", "mode_advanced")

    list(
      state   = function() captureInputs(input, stateIds),
      restore = function(st) restoreInputs(session, st, stateIds, textStateIds)
    )
  })
}
