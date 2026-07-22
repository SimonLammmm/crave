run_exorcise <- function(input) {
  output <- NULL
  # File dictionaries
  exorcise_prefix <- paste0(exorcise_hostroot, "/")
  exorcise_temp <- paste0(exorcise_root, "/temp/")
  dir.create(exorcise_temp, showWarnings = F, recursive = T)
  
  exorcise_genomes <- c("GRCh38 Ensembl 111"    = "hg38.2020-09-22.2bit",
                        "GRCm39 Ensembl 111"    = "mm39.2020-07-30.2bit")
  exorcise_exomes <- c("GRCh38 Ensembl 111"     = "hsa.grch38.refseqall.tsv.gz",
                       "GRCm39 Ensembl 111"     = "mmu.grcm39.refseqall.tsv.gz")
  exorcise_priorities <- c("GRCh38 Ensembl 111" = "hsa.priorities.tsv.gz",
                           "GRCm39 Ensembl 111" = "mmu.priorities.tsv.gz")
  
  # Mode dictionary
  exorcise_modes <- c("Knockout" = "ko",
                      "Inhibition" = "i",
                      "Activation" = "a",
                      "Cytosine base editor" = "cbe",
                      "Adenine base editor" = "abe")
  
  # Create infile
  exorcise_taskid <- paste0(sample(c(LETTERS[1:6], 0:9), 8, F), collapse = "")
  exorcise_input <- list()
  opt <- list()
  if(is.null(input$exorcise_upload)) {
    # From form
    exorcise_input$seq <- unlist(strsplit(input$exorcise_seq, "[ \t\n\r,]"))
    exorcise_input$seq <- exorcise_input$seq[exorcise_input$seq != ""]
    if (length(input$exorcise_orig) == length(input$exorcise_seq)) {
      exorcise_input$orig <- unlist(strsplit(input$exorcise_orig, "[ \t\n\r,]"))
      exorcise_input$orig <- exorcise_input$orig[exorcise_input$orig != ""]
      exorcise_harm <- T
      opt$control <- input$exorcise_control_string
      opt$control_type <- NULL
    }
  } else {
    # From CRISPick
    crispick <- fread(input$exorcise_upload$datapath[1])
    exorcise_input$seq <- crispick$`sgRNA Sequence`
    exorcise_input$orig <- crispick$`Input`
    exorcise_harm <- T
    opt$control <- "\\(NEG_CONTROL\\)"
    opt$control_type <- NULL
  }
  # Fail if no sequences
  if(length(exorcise_input$seq) == 0) {
    output$exorcise_console <- renderText("No sequences submitted. Please revise your input.")
    # Fail if annotations submitted but not the same number as the number of sequences
  } else if(length(exorcise_input$orig) > 0 & length(exorcise_input$seq) != length(exorcise_input$orig)) {
    output$exorcise_console <- renderText("Number of submitted sequences not equal to number of submitted annotations. Please revise your input.")
    # Fail if any polyN >= 17
  } else if(any(grepl("A{17}|C{17}|G{17}|T{17}", exorcise_input$seq))) {
    output$exorcise_console <- renderText("Mononucleotide run detected, not allowed. Please revise your input.")
    # Fail if any N + PAM not N < 21
  } else if(any(nchar(exorcise_input$seq) + nchar(gsub("[^ACGT]", "", input$exorcise_pam)) < 21)) {
    output$exorcise_console <- renderText("Short sequence detected, not allowed. Please revise your input or PAM.")
    # Succeed
  } else {
    fwrite(exorcise_input, paste0(exorcise_temp, "exorcise-input_", exorcise_taskid, ".tsv"), sep = "\t")
    
    # Populate options
    opt$infile <-       paste0("exorcise-input_", exorcise_taskid, ".tsv")
    opt$outdir <-       paste0("exorcise-output_", exorcise_taskid)
    opt$seq <-          "1"
    if(exorcise_harm) { opt$harm <- "2" } else { opt$harm <- NULL}
    if(nchar(input$exorcise_pam) > 0) {
      opt$pam <-        input$exorcise_pam
    } else {
      opt$pam <-        "NGG"
    }
    opt$genome <-            exorcise_genomes[input$exorcise_genome]
    opt$exome <-             exorcise_exomes[input$exorcise_genome]
    opt$expression <-        NULL
    opt$exprcutoff <-        NULL
    if(nchar(input$exorcise_mode_advanced) > 0) {
      opt$mode <- input$exorcise_mode_advanced
    } else {
      opt$mode <- exorcise_modes[input$exorcise_mode]
    }
    # Run Exorcise Docker
    # Build command
    exorcise_image <- exorcise_docker
    command <- paste0("docker run --rm -v ", exorcise_prefix, ":/data/ -v ", exorcise_prefix, "/temp/:/temp/ ", exorcise_image, " exorcise",
                      " -i /temp/", opt$infile,
                      " -o /temp/", opt$outdir,
                      " -z ", opt$pam,
                      " -q ", opt$mode,
                      " -g ", opt$seq,
                      " -v /data/", opt$genome,
                      " -w /data/", opt$exome)
    if(!is.null(opt$expression) & !is.null(opt$exprcutoff)) {
      command <- paste0(command, " -x /data/", opt$expression, " -k ", opt$exprcutoff)
    }
    if(!is.null(opt$harm)) {
      command <- paste0(command, " -n ", opt$harm)
    }
    # Execute
    system(command)
    # Check the output to see if it was successful
    exorcise_outfile <- paste0(exorcise_temp, opt$outdir, "/exorcise.tsv")
    # Show error if fail
    if(!file.exists(exorcise_outfile)) {
      exorcise_log <- fread(Sys.glob(paste0(exorcise_temp, opt$outdir, "/logfile*.log")))
      output$exorcise_console <- renderText(exorcise_log[nrow(exorcise_log), 1])
    } else {
      # Offer output if success
      output$exorcise_console <- renderText("Exorcise done, results are ready for download.")
      exorcise_result <- fread(paste0(exorcise_temp, opt$outdir, "/exorcise.tsv"))
      output$exorcise_download_input <- downloadHandler(
        filename = paste0("exorcise-input_", exorcise_taskid, ".tsv.gz"),
        content = function(f) fwrite(exorcise_input, f, sep = "\t")
      )
      output$exorcise_download_result <- downloadHandler(
        filename = paste0("exorcise-output_", exorcise_taskid, ".tsv.gz"),
        content = function(f) fwrite(exorcise_result, f, sep = "\t")
      )
    }
    # Clean up
    file.remove(paste0(exorcise_temp, "exorcise-input_", exorcise_taskid, ".tsv"))
    unlink(paste0(exorcise_temp, opt$outdir), recursive = T)
  }
  return(output)
}