#### Shiny server ####
server <- function(input, output, session) {
  #### Privacy message ####
  showNotification(
    ui = privacy_toast,
    duration = 8, type = "message")
  
  #### Shiny server: init ####
  # Explore pickers
  initialiseExplorePickers <- function() {
    updateSelectizeInput(session, "comparisons_selected", server = TRUE, choices = comparisons$`FriendlyID`, options = list(maxOptions = 1000000))
    updateSelectizeInput(session, "genes_selected", server = TRUE, choices = genes, options = list(maxOptions = 1000))
  }
  # Correlate pickers
  initialiseCorrelatePickers <- function() {
    updateSelectizeInput(session, "genes_queried", server = TRUE, choices = genes, options = list(maxOptions = 1000))
    #updateCheckboxInput(session, 'gq_gene_filter_baseedit', value = F)
    #updateCheckboxInput(session, 'gq_gene_filter_array', value = F)
    #updateCheckboxInput(session, 'gq_gene_filter_proteincoding', value = T)
    #updateCheckboxInput(session, 'gq_gene_filter_pseudo', value = F)
    #updateCheckboxInput(session, 'gq_gene_filter_ncrna', value = F)
    #updateCheckboxInput(session, 'gq_gene_filter_snrna', value = F)
    #updateCheckboxInput(session, 'gq_gene_filter_snorna', value = F)
    #updateCheckboxInput(session, 'gq_gene_filter_scrna', value = F)
    #updateCheckboxInput(session, 'gq_gene_filter_rrna', value = F)
    #updateCheckboxInput(session, 'gq_gene_filter_trna', value = F)
    updateSelectizeInput(session, 'gq_filter_citation', server = T, choices = sort(gsub("<.+?>", "", levels(comparisons$`Citation`))))
    updateSelectizeInput(session, 'gq_filter_library', server = T, choices = sort(levels(comparisons$`Library`)))
    updateSelectizeInput(session, 'gq_filter_source', server = T, choices = sort(levels(comparisons$`Source`)), selected = sort(levels(comparisons$`Source`)))
    updateSelectizeInput(session, 'gq_filter_organism', server = T, choices = sort(levels(comparisons$`Organism`)))
    updateSelectizeInput(session, 'gq_filter_timepoint', server = T, choices = sort(levels(comparisons$`Timepoint`)))
    updateSelectizeInput(session, 'gq_filter_kind', server = T, choices = sort(levels(comparisons$`Kind`)))
    updateSelectizeInput(session, 'gq_filter_days_grown_diff', server = T, choices = sort(levels(comparisons$`Days grown (diff)`)))
    updateSelectizeInput(session, 'gq_filter_treatment_diff', server = T, choices = sort(levels(comparisons$`Treatment (diff)`)))
    updateSelectizeInput(session, 'gq_filter_dose_diff', server = T, choices = sort(levels(comparisons$`Dose (diff)`)))
    updateSelectizeInput(session, 'gq_filter_knockout_diff', server = T, choices = sort(levels(comparisons$`Knockout (diff)`)))
    updateSelectizeInput(session, 'gq_filter_cellline_diff', server = T, choices = sort(levels(comparisons$`Cell line (diff)`)))
    updateSelectizeInput(session, 'gq_filter_days_grown_ref', server = T, choices = sort(levels(comparisons$`Days grown (ref)`)))
    updateSelectizeInput(session, 'gq_filter_treatment_ref', server = T, choices = sort(levels(comparisons$`Treatment (ref)`)))
    updateSelectizeInput(session, 'gq_filter_dose_ref', server = T, choices = sort(levels(comparisons$`Dose (ref)`)))
    updateSelectizeInput(session, 'gq_filter_knockout_ref', server = T, choices = sort(levels(comparisons$`Knockout (ref)`)))
    updateSelectizeInput(session, 'gq_filter_cellline_ref', server = T, choices = sort(levels(comparisons$`Cell line (ref)`)))
    updateSelectizeInput(session, "gq_filter_contrast", server = TRUE, choices = comparisons$`FriendlyID`, options = list(maxOptions = 1000000))
    updateNumericInput(session, 'gq_cutoff', value = defaults_gq_cutoff_selected)
    updateSelectizeInput(session, 'gq_filter_timepoint', selected = c("Matched timepoints", "From experiment start"))
    updateCheckboxInput(session, 'gq_filter_exorcised', value = T)
    updateTextInput(session, 'gq_filter_custom', value = "")
    updateSelectizeInput(session, 'analysis_method_gq', selected = "DrugZ")
  }
  # Explore customisations
  initialiseExploreCustomise <- function() {
    updateNumericInput(session, 'explore_customise_height', value = defaults_explore_customise_height)
    updateNumericInput(session, 'explore_customise_width', value = defaults_explore_customise_width)
    updateCheckboxInput(session, 'explore_customise_width_auto', value = defaults_explore_customise_width_auto)
    updateColourInput(session, 'explore_customise_point_colour', value = defaults_explore_customise_point_colour)
    updateColourInput(session, 'explore_customise_goi_colour', value = defaults_explore_customise_goi_colour)
    updateSelectizeInput(session, 'explore_customise_x_axis_log', selected = defaults_explore_customise_x_axis_log)
    updateSelectizeInput(session, 'explore_customise_y_axis_log', selected = defaults_explore_customise_y_axis_log)
    updateCheckboxInput(session, 'explore_customise_y_equals_x', value = defaults_explore_customise_y_equals_x)
    updateCheckboxInput(session, 'explore_customise_density', value = defaults_explore_customise_density)
    updateCheckboxInput(session, 'explore_customise_rug_x', value = defaults_explore_customise_rug_x)
    updateCheckboxInput(session, 'explore_customise_rug_y', value = defaults_explore_customise_rug_y)
  }
  # Correlate customisations
  initialiseCorrelateCustomise <- function() {
    updateNumericInput(session, 'correlate_customise_height', value = defaults_correlate_customise_height)
    updateNumericInput(session, 'correlate_customise_width', value = defaults_correlate_customise_width)
    updateCheckboxInput(session, 'correlate_customise_width_auto', value = defaults_correlate_customise_width_auto)
    updateColourInput(session, 'correlate_customise_point_colour', value = defaults_correlate_customise_point_colour)
    updateSelectizeInput(session, 'correlate_customise_genequery_fillscheme', selected = defaults_correlate_customise_genequery_fillscheme)
    updateColourInput(session, 'correlate_customise_heatmap_high', value = defaults_correlate_customise_heatmap_high)
    updateColourInput(session, 'correlate_customise_heatmap_mid', value = defaults_correlate_customise_heatmap_mid)
    updateColourInput(session, 'correlate_customise_heatmap_low', value = defaults_correlate_customise_heatmap_low)
  }
  initialiseExplorePickers()
  initialiseCorrelatePickers()
  
  # Guide library pickers
  initialiseGuideLibraryPickers <- function() {
    updateSelectizeInput(session, 'guideLibrary_targets', selected = character(0), choices = guideLibraryTargetList, server = T, options = list(maxOptions = 1000))
    updateSelectizeInput(session, 'guideLibrary_assembly', selected = character(0), choices = guideLibraryAssemblyList, server = T, options = list(maxOptions = 1000))
    updateSelectizeInput(session, 'guideLibrary_pam', selected = character(0), choices = guideLibraryPAMList, server = T, options = list(maxOptions = 1000))
    updateSelectizeInput(session, 'guideLibrary_library', selected = character(0), choices = guideLibraryLibraryList, server = T, options = list(maxOptions = 1000))
    updateSelectizeInput(session, 'guideLibrary_organism', selected = character(0), choices = guideLibraryOrganismList, server = T, options = list(maxOptions = 1000))
    updateCheckboxGroupInput(session, 'guideLibrary_geneclasses', selected = character(0))
    updateCheckboxGroupInput(session, 'guideLibrary_chromosomes', selected = character(0))
    updateCheckboxGroupInput(session, 'guideLibrary_chemistry', selected = character(0))
  }
  
  # Guide library table
  initialiseGuideLibraryTable <- function() {
    if(!exists("file_guideLibrary")) return()   # no guides dataset loaded
    # SQL
    guideLibrary <<- dbConnect(SQLite(), file_guideLibrary)
    guideLibraryTargetList <<- sqldf(connection = guideLibrary, x = "SELECT DISTINCT Target FROM guideLibrary") %>% unlist() %>% setNames(NULL)
    guideLibraryAssemblyList <<- sqldf(connection = guideLibrary, x = "SELECT DISTINCT Assembly FROM guideLibrary") %>% unlist() %>% setNames(NULL)
    guideLibraryPAMList <<- sqldf(connection = guideLibrary, x = "SELECT DISTINCT PAM FROM guideLibrary") %>% unlist() %>% setNames(NULL)
    guideLibraryLibraryList <<- sqldf(connection = guideLibrary, x = "SELECT DISTINCT Library FROM guideLibrary") %>% unlist() %>% setNames(NULL)
    guideLibraryOrganismList <<- sqldf(connection = guideLibrary, x = "SELECT DISTINCT Organism FROM guideLibrary") %>% unlist() %>% setNames(NULL)
    
    # Draw table
    dummy <- tibble(V1 = "Fill the form in the sidebar and then click the Guides button.")
    output$guideLibrary_table <- renderDT(dummy, filter = "top", server = T,
                                          options = list(iDisplayLength = 100, scrollX = T), escape = F)
    output$guideLibrary_download <- downloadHandler(
      filename = "export.csv.gz",
      content = function(f) fwrite(guideLibrary, f, sep = ",")
    )
  }
  
  # Libraries table
  initialiseLibrariesTable <- function() {
    if(!exists("file_libraries")) return()   # no libraries file loaded
    libraries <- fread(file_libraries) %>%
      transmute(Library = case_when(grepl("^http", `Source link`) ~ paste0("<a href=", `Source link`, " target=\"_blank\">", `Library`, "</a>"),
                                    grepl("^http", `Reference`) ~ paste0("<a href=", Reference, " target=\"_blank\">", `Library`, "</a>"),
                                    T ~ Library),
                `Guide length`, `Number of guides`, `Number of targets`, Species, `CRISPR chemistry`, PAM,
                Citation = case_when(grepl("^http", Reference) ~ paste0("<a href=", Reference, " target=\"_blank\">", Citation, "</a>"),
                                     T ~ Citation),
                Lab)
    output$libraries_table <- renderDT(libraries, filter = "top", server = T,
                                       options = list(iDisplayLength = 100, scrollX = T), escape = F)
  }
  # Ontology table
  initialiseOntologyTable <- function() {
    if(is.null(ontology) || nrow(ontology) == 0) return()   # no ontology loaded
    output$ontology_table <- renderDT(ontology %>% mutate(class = factor(class)), filter = "top", server = T,
                                      options = list(iDisplayLength = 100, scrollX = T), escape = F)
  }
  
  # Experiments table
  initialiseExperimentsTables <- function() {
    # Experiments table
    output$experiments_table <- renderDT(experiments %>% select(-`Experiment ID`), filter = "top", server = TRUE,
                                         options = list(iDisplayLength = 100, scrollX = T), escape = F)
    experiments_table_proxy <<- dataTableProxy("experiments_table")
  }
  # Comparisons table
  initialiseComparisonsTables <- function() {
    # Comparisons table
    output$comparisons_table <- renderDT(comparisons %>% select(-`Comparison ID`, -`Experiment ID`, -`FriendlyID`), filter = "top", server = TRUE,
                                         options = list(iDisplayLength = 100, scrollX = T), escape = F)
    comparisons_table_proxy <<- dataTableProxy("comparisons_table")
  }
  initialiseComparisonSearch <- function() {
    # Comparisons search
    citations <- paste(gsub("<.+?>", "", experiments$`Citation`), collapse = "; ")
    citations <- case_when(nchar(citations) > 255 ~ paste0(sub("(^.{255}.+?);(.+)$", "\\1", citations),
                                                           "... and ",
                                                           length(gregexpr(";", sub("(^.{255}.+?);(.+)$", "\\2", citations))[[1]]),
                                                           " more."),
                           T ~ paste0(sub("^(.+);(.+?)$", "\\1", citations),
                                      "; and ",
                                      sub("^(.+);(.+?)$", "\\2", citations)))
    output$comparison_selector_experiments <- renderText(citations)
    updateSelectizeInput(session,     'comparison_selector_library',              selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Library`)))))
    updateSelectizeInput(session,     'comparison_selector_organism',             selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Organism`)))))
    updateSelectizeInput(session,     'comparison_selector_source',               selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Source`)))))
    updateSelectizeInput(session,     'comparison_selector_citation',             selected = "(any)", choices = c("(any)", sort(unique(as.character(gsub("<.+?>", "", comparisons$`Citation`))))))
    updateSelectizeInput(session,     'comparison_selector_kind',                 selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Kind`)))))
    updateSelectizeInput(session,     'comparison_selector_timepoint',            selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Timepoint`)))))
    updateSelectizeInput(session,     'comparison_selector_endpoint',             selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Endpoint`)))))
    updateSelectizeInput(session,     'comparison_selector_daysgrown_diff',       selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Days grown (diff)`)))))
    updateSelectizeInput(session,     'comparison_selector_daysgrown_ref',        selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Days grown (ref)`)))))
    updateSelectizeInput(session,     'comparison_selector_treatment_diff',       selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Treatment (diff)`)))))
    updateSelectizeInput(session,     'comparison_selector_treatment_ref',        selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Treatment (ref)`)))))
    updateSelectizeInput(session,     'comparison_selector_dose_diff',            selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Dose (diff)`)))))
    updateSelectizeInput(session,     'comparison_selector_dose_ref',             selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Dose (ref)`)))))
    updateSelectizeInput(session,     'comparison_selector_knockout_diff',        selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Knockout (diff)`)))))
    updateSelectizeInput(session,     'comparison_selector_knockout_ref',         selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Knockout (ref)`)))))
    updateSelectizeInput(session,     'comparison_selector_cell_line_diff',       selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Cell line (diff)`)))))
    updateSelectizeInput(session,     'comparison_selector_cell_line_ref',        selected = "(any)", choices = c("(any)", sort(unique(as.character(comparisons$`Cell line (ref)`)))))
  }
  
  # Hart 2017 Table S2 core essentials
  essentials <- 
    c("AARS", "ABCE1", "ABCF1", "ACTB", "ACTL6A", "ACTR10", "ACTR2", "ADSL", "ADSS",
      "AHCY", "ALG1", "ALG14", "ALG2", "ANAPC2", "ANAPC4", "ANAPC5", "AQR", "ARCN1",
      "ARIH1", "ARL2", "ATP2A2", "ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5J2-PTCD1", "ATP5L",
      "ATP5O", "ATP6V0B", "ATP6V0C", "ATP6V1A", "ATP6V1D", "ATP6V1E1", "ATR", "AURKB", "BANF1",
      "BIRC5", "BUB1B", "BUB3", "BUD31", "BYSL", "C10orf2", "C1orf109", "C21orf59", "C3orf17",
      "C9orf114", "CCDC84", "CCDC94", "CCNA2", "CCNH", "CCNK", "CCT2", "CCT3", "CCT4",
      "CCT5", "CCT6A", "CCT7", "CCT8", "CDC123", "CDC16", "CDC20", "CDC27", "CDC37",
      "CDC5L", "CDC73", "CDK1", "CDK7", "CDK9", "CDT1", "CEBPZ", "CENPA", "CENPC",
      "CFL1", "CHAF1A", "CHAF1B", "CHEK1", "CHERP", "CHMP2A", "CHMP6", "CIAO1", "CINP",
      "CIRH1A", "CKAP5", "CLNS1A", "CLP1", "CLTC", "CMPK1", "CMTR1", "CNOT3", "COA5",
      "COPA", "COPB1", "COPB2", "COPS3", "COPS6", "COPZ1", "COQ4", "COX10", "COX11",
      "COX15", "COX4I1", "COX5B", "COX6B1", "CPSF1", "CPSF2", "CPSF3", "CPSF4", "CRNKL1",
      "CSE1L", "CTDP1", "CTPS1", "CTR9", "CYCS", "DAD1", "DBR1", "DCTN5", "DDB1",
      "DDOST", "DDX10", "DDX18", "DDX20", "DDX21", "DDX27", "DDX41", "DDX47", "DDX49",
      "DDX55", "DDX56", "DGCR8", "DHODH", "DHPS", "DHX15", "DHX33", "DHX37", "DHX8",
      "DHX9", "DIEXF", "DIMT1", "DIS3", "DKC1", "DLST", "DMAP1", "DNAJA3", "DNAJC9",
      "DNM2", "DNMT1", "DOLK", "DONSON", "DPAGT1", "DTL", "DTYMK", "DYNC1I2", "ECD",
      "EEF2", "EFTUD2", "EIF2B1", "EIF2B3", "EIF2B5", "EIF2S1", "EIF2S2", "EIF2S3", "EIF3A",
      "EIF3B", "EIF3C", "EIF3D", "EIF3G", "EIF3I", "EIF4A3", "EIF5A", "EIF5B", "EIF6",
      "ELAC2", "ELL", "EPRS", "ERCC2", "ERCC3", "ERH", "EXOSC2", "EXOSC3", "EXOSC4",
      "EXOSC6", "EXOSC7", "EXOSC8", "FAM96B", "FARS2", "FARSA", "FARSB", "FAU", "FNTA",
      "FNTB", "FTSJ3", "GABPA", "GAPDH", "GART", "GEMIN5", "GEMIN8", "GFM1", "GGPS1",
      "GINS2", "GINS3", "GINS4", "GMPPB", "GMPS", "GNB2L1", "GNL3", "GPN3", "GPS1",
      "GRPEL1", "GRWD1", "GSPT1", "GTF2B", "GTF2H1", "GTF2H2C", "GTF2H4", "GTF3A", "GTF3C1",
      "GTF3C2", "GTF3C5", "GTPBP4", "GUK1", "HARS", "HAUS1", "HAUS5", "HCFC1", "HDAC3",
      "HEATR1", "HINFP", "HIST1H2AJ", "HIST2H2AA3", "HJURP", "HNRNPC", "HNRNPK", "HNRNPL", "HNRNPU",
      "HSD17B10", "HSPA9", "HSPD1", "HUWE1", "HYPK", "IARS", "IGBP1", "ILF3", "IMP3",
      "IMP4", "INTS1", "INTS3", "INTS8", "INTS9", "IPO13", "ISCU", "ISG20L2", "KANSL3",
      "KARS", "KAT8", "KIF11", "KIF23", "KPNB1", "KRI1", "KRR1", "LARS", "LAS1L",
      "LONP1", "LRR1", "LSG1", "LSM11", "LSM12", "LSM2", "LSM7", "LUC7L3", "MAD2L1",
      "MAGOH", "MAK16", "MARS", "MARS2", "MASTL", "MCM3", "MCM3AP", "MCM4", "MCM5",
      "MCM7", "MDN1", "MED11", "MED12", "MED18", "MED27", "MED30", "MEPCE", "METTL16",
      "MMS22L", "MPHOSPH10", "MRP63", "MRPL18", "MRPL28", "MRPL38", "MRPL4", "MRPL43", "MRPL45",
      "MRPL46", "MRPL53", "MRPS14", "MRPS24", "MRPS34", "MSTO1", "MTG2", "MVK", "MYBBP1A",
      "MYC", "NAA10", "NAA38", "NAA50", "NAMPT", "NAPA", "NARFL", "NARS", "NAT10",
      "NCBP1", "NCBP2", "NDC80", "NDUFA13", "NEDD8", "NELFB", "NHP2", "NHP2L1", "NIP7",
      "NKAP", "NLE1", "NMD3", "NMT1", "NOC4L", "NOL10", "NOL11", "NOL6", "NOL9",
      "NOP16", "NOP2", "NOP56", "NOP9", "NPLOC4", "NSA2", "NSF", "NUDC", "NUDCD3",
      "NUDT21", "NUDT4", "NUF2", "NUP133", "NUP155", "NUP160", "NUP214", "NUP85", "NUP88",
      "NUP93", "NUS1", "NUTF2", "NVL", "NXF1", "OGDH", "OGT", "ORAOV1", "ORC6",
      "OSGEP", "PABPC1", "PAFAH1B1", "PAICS", "PAK1IP1", "PCID2", "PCNA", "PFDN2", "PFN1",
      "PGAM1", "PGGT1B", "PGK1", "PHB", "PHB2", "PHF5A", "PKMYT1", "PLK1", "PLRG1",
      "PMPCA", "PMPCB", "PNKP", "POLA2", "POLR1A", "POLR1B", "POLR1C", "POLR2A", "POLR2B",
      "POLR2C", "POLR2D", "POLR2E", "POLR2G", "POLR2H", "POLR2I", "POLR2L", "POLR3A", "POLR3C",
      "POLR3H", "POLR3K", "POLRMT", "POP1", "POP5", "PPA1", "PPAN", "PPAT", "PPIL2",
      "PPP2CA", "PPP2R4", "PPP4C", "PPWD1", "PREB", "PRELID1", "PRIM1", "PRMT1", "PRMT5",
      "PRPF19", "PRPF31", "PRPF38A", "PRPF38B", "PRPF4", "PRPF8", "PSMA1", "PSMA2", "PSMA3",
      "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB7",
      "PSMC2", "PSMC3", "PSMC5", "PSMC6", "PSMD1", "PSMD11", "PSMD12", "PSMD13", "PSMD14",
      "PSMD3", "PSMD4", "PSMG3", "PTPN23", "PUF60", "PWP2", "QARS", "RABGGTB", "RACGAP1",
      "RAD21", "RAD51C", "RAD51D", "RAE1", "RAN", "RANGAP1", "RARS2", "RBBP6", "RBM14",
      "RBM17", "RBM8A", "RBMX", "RBX1", "RCC1", "RCL1", "RFC2", "RFC4", "RFC5",
      "RFK", "RHEB", "RIOK2", "RNF20", "RNGTT", "ROMO1", "RPA1", "RPA2", "RPF2",
      "RPL10A", "RPL11", "RPL12", "RPL13", "RPL14", "RPL18", "RPL18A", "RPL19", "RPL23",
      "RPL24", "RPL27", "RPL27A", "RPL3", "RPL30", "RPL35", "RPL35A", "RPL36", "RPL37A",
      "RPL4", "RPL6", "RPL8", "RPLP0", "RPLP1", "RPLP2", "RPP21", "RPP38", "RPS11",
      "RPS12", "RPS13", "RPS15A", "RPS16", "RPS18", "RPS19", "RPS2", "RPS20", "RPS21",
      "RPS23", "RPS3", "RPS4X", "RPS5", "RPS6", "RPS7", "RPS8", "RRM1", "RRP1",
      "RRP12", "RRS1", "RTCB", "RUVBL2", "SACM1L", "SAE1", "SAMM50", "SAP18", "SARS",
      "SARS2", "SART3", "SBNO1", "SDAD1", "SDHC", "SEC13", "SEH1L", "SF1", "SF3A2",
      "SF3A3", "SF3B1", "SF3B2", "SF3B3", "SF3B5", "SKP1", "SLC35B1", "SLMO2", "SLU7",
      "SMC1A", "SMC2", "SMC4", "SMU1", "SNAPC1", "SNAPC2", "SNAPC4", "SNRNP200", "SNRNP25",
      "SNRNP27", "SNRNP35", "SNRNP70", "SNRPA1", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPF", "SNW1",
      "SPATA5L1", "SPC24", "SPC25", "SRBD1", "SRP19", "SRRM1", "SRRT", "SRSF1", "SRSF2",
      "SRSF3", "SRSF7", "SS18L2", "SSU72", "SUPT5H", "SUPT6H", "SUPV3L1", "SYMPK", "SYS1",
      "TAF1B", "TAF6", "TANGO6", "TARS", "TBCD", "TBL3", "TCP1", "TELO2", "TFAM",
      "TFRC", "THOC2", "THOC3", "THOC5", "TICRR", "TIMM10", "TIMM13", "TIMM23", "TIMM44",
      "TMEM258", "TNPO3", "TOMM22", "TOMM40", "TONSL", "TOP1", "TOP2A", "TPT1", "TPX2",
      "TRAPPC1", "TRAPPC3", "TRIAP1", "TRMT112", "TRMT5", "TRNAU1AP", "TRRAP", "TSR1", "TTC1",
      "TTC27", "TTI1", "TTI2", "TUBB", "TUBG1", "TUBGCP2", "TUBGCP3", "TUBGCP6", "TUFM",
      "TUT1", "TXN", "TXNL4A", "U2AF1", "U2AF2", "UBA1", "UBA52", "UBE2L3", "UBE2M",
      "UBE2N", "UBL5", "UBTF", "UPF1", "UPF2", "UQCRC1", "UQCRFS1", "UROD", "USP39",
      "USP5", "USPL1", "UTP15", "UTP20", "UTP23", "UXT", "VARS", "VARS2", "VCP",
      "VPS25", "VPS28", "WARS", "WBSCR22", "WDR12", "WDR25", "WDR3", "WDR33", "WDR43",
      "WDR61", "WDR70", "WDR74", "WDR75", "WDR77", "WDR92", "WEE1", "XAB2", "XPO1",
      "XRCC6", "YARS", "YARS2", "YRDC", "ZBTB8OS", "ZMAT5", "ZNF131", "ZNF259", "ZNF574")
  
  #### Shiny server: react to input: home ####
  # When jump buttons are clicked on the home page
  observeEvent(input$jumpCorrelate, {
    updateTabsetPanel(session, inputId = "main", selected = "Correlate")
  } )
  observeEvent(input$jumpExplore, {
    updateTabsetPanel(session, inputId = "main", selected = "Explore")
  } )
  observeEvent(input$jumpGuideLibrary, {
    updateTabsetPanel(session, inputId = "main", selected = "Guides")
  } )
  observeEvent(input$jumpLibraries, {
    updateTabsetPanel(session, inputId = "main", selected = "Libraries")
  } )
  observeEvent(input$jumpOntology, {
    updateTabsetPanel(session, inputId = "main", selected = "Ontology")
  } )
  observeEvent(input$jumpExorcise, {
    updateTabsetPanel(session, inputId = "main", selected = "Exorcise")
  } )
  # When the data refresh button is clicked
  # Store the last tab viewed by the user, except for the Refresh data and Save/Load tabs. Defaults to "Home"
  lastTab <- reactiveVal("Home")
  observeEvent(input$main, {
    if(!(input$main %in% c("Refresh data", "Save/Load"))) {
      lastTab(input$main)
    }
  })

  observeEvent(input$main, {
    if (input$main == "Refresh data") {
      # Return the user to their original tab
      updateTabsetPanel(session, "main", selected = lastTab())
      # Refresh data
    load()
    initialiseExplorePickers()
    initialiseCorrelatePickers()
    initialiseExploreCustomise()
    initialiseCorrelateCustomise()
    initialiseLibrariesTable()
    initialiseOntologyTable()
    initialiseExperimentsTables()
    initialiseComparisonsTables()
    initialiseComparisonSearch()
  }})
  #### Shiny server: react to input: explore ####
  # When the validate genes button is pressed
  observeEvent(input$genes_selected_text_validate, {
    if(nchar(input$genes_selected_text) >= 3) {
      if(exists("genes_filtered")) {
        genes_check <- genes_filtered
      } else {
        genes_check <- genes
      }
      genes_typed <- unique(unlist(strsplit(input$genes_selected_text, split = "[,; ]")))
      genes_typed <- unique(genes_check[toupper(genes_check) %in% toupper(genes_typed)])
      if(length(genes_typed) <= 1000) {
        updateSelectizeInput(session, "genes_selected", selected = unique(c(input$genes_selected, genes_typed)), choices = genes_check, server = T)
      } else {
        updateSelectizeInput(session, "genes_selected", selected = unique(c(input$genes_selected, genes_typed)), choices = genes_check, server = F)
      }
    }
  })
  # When the copy genes from correlate button is pressed
  observeEvent(input$genes_copy_from_correlate, {
    if(exists("genes_filtered")) {
      genes_check <- genes_filtered
    } else {
      genes_check <- genes
    }
    updateSelectizeInput(session, "genes_selected", selected = input$genes_queried, choices = genes_check, server = T)
  })
  # When the core essentials link is pressed
  observeEvent(input$genes_select_essentialome, {
    if(exists("genes_filtered")) {
      genes_check <- genes_filtered
    } else {
      genes_check <- genes
    }
    updateTextAreaInput(session, "genes_selected_text", value = paste(essentials, collapse = ", "))
  })
  # When the clear genes button is pressed
  observeEvent(input$genes_clear, {
    if(exists("genes_filtered")) {
      genes_check <- genes_filtered
    } else {
      genes_check <- genes
    }
    updateSelectizeInput(session, "genes_selected", selected = NULL, choices = genes_check, server = T)
  })
  # When experiments are picked from the table
  observeEvent(input$experiments_table_rows_selected, ignoreNULL = F, {
    experiments_filtered <- experiments[input$experiments_table_rows_selected,]
    experiments_selected <- experiments_filtered$`Experiment ID`
    # Filter the comparisons on the picked experiments
    if(length(experiments_selected) > 0) {
      comparisons_filtered <- comparisons %>% filter(`Experiment ID` %in% experiments_selected)
      output$comparisons_table <- renderDT(comparisons_filtered %>% select(-`Comparison ID`, -`Experiment ID`, -`FriendlyID`), filter = "top", server = TRUE,
                                           options = list(iDisplayLength = 100, scrollX = T), escape = F)
      # Update the comparisons text picker choices
      updateSelectizeInput(session, "comparisons_selected", choices = comparisons_filtered$`FriendlyID`, server = TRUE)
      # Update the comparisons search choices
      if(length(experiments_filtered$`Citation`) > 1) {
        citations_validated <- paste(gsub("<.+?>", "", experiments_filtered$`Citation`), collapse = "; ")
        citations_validated <- case_when(nchar(citations_validated) > 255 ~ paste0(sub("(^.{255}.+?);(.+)$", "\\1", citations_validated),
                                                                                   "... and ",
                                                                                   length(gregexpr(";", sub("(^.{255}.+?);(.+)$", "\\2", citations_validated))[[1]]),
                                                                                   " more."),
                                         T ~ paste0(sub("^(.+);(.+?)$", "\\1", citations_validated),
                                                    "; and ",
                                                    sub("^(.+);(.+?)$", "\\2", citations_validated)))
      } else {
        citations_validated <- experiments_filtered$`Citation`
        citations_validated <- gsub("<.+?>", "", citations_validated) # remove html tags
      }
      output$comparison_selector_experiments <- renderText(citations_validated)
      updateSelectizeInput(session,     'comparison_selector_library'             , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Library`)))))
      updateSelectizeInput(session,     'comparison_selector_organism'            , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Organism`)))))
      updateSelectizeInput(session,     'comparison_selector_source'              , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Source`)))))
      updateSelectizeInput(session,     'comparison_selector_citation'            , selected = "(any)", choices = c("(any)", as.character(sort(unique(gsub("<.+?>", "", comparisons_filtered$`Citation`))))))
      updateSelectizeInput(session,     'comparison_selector_kind'                , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Kind`)))))
      updateSelectizeInput(session,     'comparison_selector_timepoint'           , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Timepoint`)))))
      updateSelectizeInput(session,     'comparison_selector_endpoint'            , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Endpoint`)))))
      updateSelectizeInput(session,     'comparison_selector_daysgrown_diff'      , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Days grown (diff)`)))))
      updateSelectizeInput(session,     'comparison_selector_daysgrown_ref'       , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Days grown (ref)`)))))
      updateSelectizeInput(session,     'comparison_selector_treatment_diff'      , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Treatment (diff)`)))))
      updateSelectizeInput(session,     'comparison_selector_treatment_ref'       , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Treatment (ref)`)))))
      updateSelectizeInput(session,     'comparison_selector_dose_diff'           , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Dose (diff)`)))))
      updateSelectizeInput(session,     'comparison_selector_dose_ref'            , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Dose (ref)`)))))
      updateSelectizeInput(session,     'comparison_selector_knockout_diff'       , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Knockout (diff)`)))))
      updateSelectizeInput(session,     'comparison_selector_knockout_ref'        , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Knockout (ref)`)))))
      updateSelectizeInput(session,     'comparison_selector_cell_line_diff'      , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Cell line (diff)`)))))
      updateSelectizeInput(session,     'comparison_selector_cell_line_ref'       , selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered$`Cell line (ref)`)))))
    } else {
      initialiseExplorePickers()
      initialiseComparisonsTables()
      initialiseComparisonSearch()
    }
  })
  
  initialiseExperimentsTables()
  initialiseComparisonsTables()
  initialiseComparisonSearch()
  # When the explore customise plot button is pressed
  observeEvent(input$explore_customise, {
    # When the explore customise plot button is pressed for the first time
    if(input$explore_customise == 1) {
      initialiseExploreCustomise()
    }
    showModal(modalDialog(
      title = "Customise plot",
      easyClose = T,
      footer = NULL,
      numericInput('explore_customise_height', "Plot height (pixels)", min = 0, max = Inf, value = input$explore_customise_height),
      numericInput('explore_customise_width', "Plot width (pixels)", min = 0, max = Inf, value = input$explore_customise_width),
      #checkboxInput('explore_customise_width_auto', "Fit plot to window width", value = input$explore_customise_width_auto),
      colourInput('explore_customise_point_colour', "Point colour", value = input$explore_customise_point_colour),
      colourInput('explore_customise_goi_colour', "Genes of interest colour", value = input$explore_customise_goi_colour),
      selectizeInput('explore_customise_x_axis_log', "x-axis scale", choices = c("Linear", "Logarithmic", "Automatic"), selected = input$explore_customise_x_axis_log),
      selectizeInput('explore_customise_y_axis_log', "y-axis scale", choices = c("Linear", "Logarithmic", "Automatic"), selected = input$explore_customise_y_axis_log),
      checkboxInput('explore_customise_y_equals_x', "Equal x and y axes (biplot only)", value = input$explore_customise_y_equals_x),
      #checkboxInput('explore_customise_rug_x', "Show x-axis rug plot (biplot only)", value = input$explore_customise_rug_x),
      #checkboxInput('explore_customise_rug_y', "Show y-axis rug plot (biplot only)", value = input$explore_customise_rug_y),
      tags$hr(),
      actionButton('explore_customise_reset', "Restore defaults", icon = icon("arrows-rotate")),
      modalButton("Save options", icon = icon("floppy-disk"))
    )
    )
  })
  # When the customise plot restore defaults button is pressed
  observeEvent(input$explore_customise_reset, {
    initialiseExploreCustomise()
  })
  # When the customise plot auto width option is selected
  observeEvent(input$explore_customise_width_auto, {
    if(input$explore_customise_width_auto) shinyjs::disable('explore_customise_width')
    else shinyjs::enable('explore_customise_width')
  })
  # When a filter is applied in Correlate
  correlate_filter_updated <- reactive({
    list(
      input$gq_filter_library,
      input$gq_filter_organism,
      input$gq_filter_source,
      input$gq_filter_citation,
      input$gq_filter_timepoint,
      input$gq_filter_kind,
      input$gq_filter_days_grown_diff,
      input$gq_filter_days_grown_ref,
      input$gq_filter_treatment_diff,
      input$gq_filter_treatment_ref,
      input$gq_filter_dose_diff,
      input$gq_filter_dose_ref,
      input$gq_filter_knockout_diff,
      input$gq_filter_knockout_ref,
      input$gq_filter_cellline_diff,
      input$gq_filter_cellline_ref,
      input$gq_filter_contrast,
      input$gq_filter_custom
    )
  })
  observeEvent(correlate_filter_updated(), ignoreNULL = F, {
    # Return the number of screens found
    comparisons_validated = gq_apply_filter(input)
    comparisons_validated = sort(comparisons$FriendlyID[comparisons$`Comparison ID` %in% comparisons_validated])
    if(length(comparisons_validated) > 1) {
      comparisons_validated <- paste(comparisons_validated, collapse = "; ")
      comparisons_validated <- case_when(nchar(comparisons_validated) > 255 ~ paste0(sub("(^.{255}.+?);(.+)$", "\\1", comparisons_validated),
                                                                                     "... and ",
                                                                                     length(gregexpr(";", sub("(^.{255}.+?);(.+)$", "\\2", comparisons_validated))[[1]]),
                                                                                     " more."),
                                         T ~ paste0(sub("^(.+);(.+?)$", "\\1", comparisons_validated),
                                                    "; and ",
                                                    sub("^(.+);(.+?)$", "\\2", comparisons_validated)))
      comparisons_validated <- paste0("Screens found: ", comparisons_validated)
    } else if (length(comparisons_validated) == 1) {
      comparisons_validated <- paste0("Screen found: ", comparisons_validated)
    } else {
      comparisons_validated <- "No screens found with the selected filters."
    }
    output$gq_filtered_screens <- renderText(comparisons_validated)
  })
  # When a filter is selected in the experiments searcher
  comparison_selector_updated <- reactive({
    list(
      input$comparison_selector_library,
      input$comparison_selector_organism,
      input$comparison_selector_source,
      input$comparison_selector_citation,
      input$comparison_selector_kind,
      input$comparison_selector_timepoint,
      input$comparison_selector_endpoint,
      input$comparison_selector_daysgrown_diff,
      input$comparison_selector_daysgrown_ref,
      input$comparison_selector_treatment_diff,
      input$comparison_selector_treatment_ref,
      input$comparison_selector_dose_diff,
      input$comparison_selector_dose_ref,
      input$comparison_selector_knockout_diff,
      input$comparison_selector_knockout_ref,
      input$comparison_selector_cell_line_diff,
      input$comparison_selector_cell_line_ref,
      input$comparison_selector_validate,
      input$experiments_table_rows_selected
    )
  })
  observeEvent(comparison_selector_updated(), ignoreNULL = F, {
    comparisons$CitationWithoutURL <- gsub("<.+?>", "", comparisons$`Citation`)
    if(length(input$experiments_table_rows_selected) > 0) experiments_filtered <- experiments[input$experiments_table_rows_selected,]
    else experiments_filtered <- experiments
    experiments_selected <- experiments_filtered$`Experiment ID`
    # Filter the comparisons on the picked experiments
    comparisons_filtered <- comparisons %>% filter(`Experiment ID` %in% experiments_selected)
    if (all(input$comparison_selector_library             != "(any)") & length(input$comparison_selector_library             ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Library`            %in% input$comparison_selector_library)
    if (all(input$comparison_selector_organism            != "(any)") & length(input$comparison_selector_organism            ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Organism`           %in% input$comparison_selector_organism)
    if (all(input$comparison_selector_source              != "(any)") & length(input$comparison_selector_source              ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Source`             %in% input$comparison_selector_source)
    if (all(input$comparison_selector_citation            != "(any)") & length(input$comparison_selector_citation            ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`CitationWithoutURL` %in% input$comparison_selector_citation)
    if (all(input$comparison_selector_kind                != "(any)") & length(input$comparison_selector_kind                ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Kind`               %in% input$comparison_selector_kind)
    if (all(input$comparison_selector_timepoint           != "(any)") & length(input$comparison_selector_timepoint           ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Timepoint`          %in% input$comparison_selector_timepoint)
    if (all(input$comparison_selector_endpoint            != "(any)") & length(input$comparison_selector_endpoint            ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Endpoint`           %in% input$comparison_selector_endpoint)
    if (all(input$comparison_selector_daysgrown_diff      != "(any)") & length(input$comparison_selector_daysgrown_diff      ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Days grown (diff)`  %in% input$comparison_selector_daysgrown_diff)
    if (all(input$comparison_selector_daysgrown_ref       != "(any)") & length(input$comparison_selector_daysgrown_ref       ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Days grown (ref)`   %in% input$comparison_selector_daysgrown_ref)
    if (all(input$comparison_selector_treatment_diff      != "(any)") & length(input$comparison_selector_treatment_diff      ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Treatment (diff)`   %in% input$comparison_selector_treatment_diff)
    if (all(input$comparison_selector_treatment_ref       != "(any)") & length(input$comparison_selector_treatment_ref       ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Treatment (ref)`    %in% input$comparison_selector_treatment_ref)
    if (all(input$comparison_selector_dose_diff           != "(any)") & length(input$comparison_selector_dose_diff           ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Dose (diff)`        %in% input$comparison_selector_dose_diff)
    if (all(input$comparison_selector_dose_ref            != "(any)") & length(input$comparison_selector_dose_ref            ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Dose (ref)`         %in% input$comparison_selector_dose_ref)
    if (all(input$comparison_selector_knockout_diff       != "(any)") & length(input$comparison_selector_knockout_diff       ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Knockout (diff)`    %in% input$comparison_selector_knockout_diff)
    if (all(input$comparison_selector_knockout_ref        != "(any)") & length(input$comparison_selector_knockout_ref        ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Knockout (ref)`     %in% input$comparison_selector_knockout_ref)
    if (all(input$comparison_selector_cell_line_diff      != "(any)") & length(input$comparison_selector_cell_line_diff      ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Cell line (diff)`   %in% input$comparison_selector_cell_line_diff)
    if (all(input$comparison_selector_cell_line_ref       != "(any)") & length(input$comparison_selector_cell_line_ref       ) >= 1) comparisons_filtered <- comparisons_filtered %>% filter(`Cell line (ref)`    %in% input$comparison_selector_cell_line_ref)
    
    fields <- tibble(
      input = c("comparison_selector_library",
                "comparison_selector_organism",
                "comparison_selector_source",
                "comparison_selector_citation",
                "comparison_selector_kind",
                "comparison_selector_timepoint",
                "comparison_selector_endpoint",
                "comparison_selector_daysgrown_diff",
                "comparison_selector_daysgrown_ref",
                "comparison_selector_treatment_diff",
                "comparison_selector_treatment_ref",
                "comparison_selector_dose_diff",
                "comparison_selector_dose_ref",
                "comparison_selector_knockout_diff",
                "comparison_selector_knockout_ref",
                "comparison_selector_cell_line_diff",
                "comparison_selector_cell_line_ref"
      ),
      column = c("Library",
                 "Organism",
                 "Source",
                 "CitationWithoutURL",
                 "Kind",
                 "Timepoint",
                 "Endpoint",
                 "Days grown (diff)",
                 "Days grown (ref)",
                 "Treatment (diff)",
                 "Treatment (ref)",
                 "Dose (diff)",
                 "Dose (ref)",
                 "Knockout (diff)",
                 "Knockout (ref)",
                 "Cell line (diff)",
                 "Cell line (ref)"))
    # Obtain indexes of each input field into the comparison table
    comparisons_table_colnames = names(comparisons %>% select(-`Comparison ID`, -`Experiment ID`, -`FriendlyID`))
    comparisons_table_columns = foreach(c = fields$column, .combine = "c") %do% which(sub("CitationWithoutURL", "Citation", c) == comparisons_table_colnames)+1
    # Initialise all comparison table filters to empty
    comparisons_table_filters = NULL
    for (i in 1:nrow(fields)) {
      # 1. Update the choices
      if (all(input[[fields$input[i]]] == "(any)")) {
        updateSelectizeInput(session,  fields$input[i] , server = F,
                             selected = input[[fields$input[i]]], choices = c("(any)", as.character(sort(unique(comparisons_filtered[[fields$column[i]]])))))
        # 2. Clear the "(any)" if something was picked
      } else if (input[[fields$input[i]]] [1] == "(any)" & length(input[[fields$input[i]]] ) > 1) {
        updateSelectizeInput(session,  fields$input[i] , server = F,
                             selected = input[[fields$input[i]]][2:length(input[[fields$input[i]]])], choices = c("(any)", as.character(sort(unique(comparisons_filtered[[fields$column[i]]])))))
        # 3. Clear choices if "(any)" was picked
      } else if (input[[fields$input[i]]][length(input[[fields$input[i]]] )] == "(any)") {
        updateSelectizeInput(session,  fields$input[i] , server = F, selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered[[fields$column[i]]]))))) }
      # 4. Restore the "(any)" if everything was cleared
      if (is.null(input[[fields$input[i]]])) {
        updateSelectizeInput(session,  fields$input[i] , server = F, selected = "(any)", choices = c("(any)", as.character(sort(unique(comparisons_filtered[[fields$column[i]]])))))
      }
      # 5. Update the comparisons table filter
      filterValue = input[[fields$input[i]]]
      filterValue = filterValue[filterValue != "(any)"]
      filterValue = paste0("[", paste0("\"", filterValue, "\"", collapse = ","), "]")
      comparisons_table_filters[comparisons_table_columns[i]] = filterValue
    }
    
    # Update the citations found
    if(length(experiments_filtered$`Citation`) > 1) {
      citations_validated <- paste(gsub("<.+?>", "", experiments_filtered$`Citation`), collapse = "; ")
      citations_validated <- case_when(nchar(citations_validated) > 255 ~ paste0(sub("(^.{255}.+?);(.+)$", "\\1", citations_validated),
                                                                                 "... and ",
                                                                                 length(gregexpr(";", sub("(^.{255}.+?);(.+)$", "\\2", citations_validated))[[1]]),
                                                                                 " more."),
                                       T ~ paste0(sub("^(.+);(.+?)$", "\\1", citations_validated),
                                                  "; and ",
                                                  sub("^(.+);(.+?)$", "\\2", citations_validated)))
    } else {
      citations_validated <- experiments_filtered$`Citation`
      citations_validated <- gsub("<.+?>", "", citations_validated) # remove html tags
    }
    output$comparison_selector_experiments <- renderText(citations_validated)
    # Update the list of screens found
    if(length(comparisons_filtered$FriendlyID) > 1) {
      comparisons_validated <- paste(sort(comparisons_filtered$FriendlyID), collapse = "; ")
      comparisons_validated <- case_when(nchar(comparisons_validated) > 255 ~ paste0(sub("(^.{255}.+?);(.+)$", "\\1", comparisons_validated),
                                                                                     "... and ",
                                                                                     length(gregexpr(";", sub("(^.{255}.+?);(.+)$", "\\2", comparisons_validated))[[1]]),
                                                                                     " more."),
                                         T ~ paste0(sub("^(.+);(.+?)$", "\\1", comparisons_validated),
                                                    "; and ",
                                                    sub("^(.+);(.+?)$", "\\2", comparisons_validated)))
    } else {
      comparisons_validated <- comparisons_filtered$`FriendlyID`
    }
    output$comparisons_selector_list <- renderText(comparisons_validated)
    # Export the screens found ready for accepting
    comparisons_exported <<- comparisons_filtered$`FriendlyID`
    # Allow users to reset the form
    comparison_selector_clearable(T)
    # Update filters in the table
    updateSearch(comparisons_table_proxy, keywords = list(global = NULL, columns = comparisons_table_filters))
  }
  
  
  )
  # When the Clear button is pressed in the comparisons searcher
  # Don't let users spam the button
  comparison_selector_clearable <- reactiveVal(T)
  observeEvent(input$comparisons_selector_clear, {
    if(comparison_selector_clearable()) {
      comparisons$CitationWithoutURL <- gsub("<.+?>", "", comparisons$`Citation`)
      if(length(input$experiments_table_rows_selected) > 0) experiments_filtered <- experiments[input$experiments_table_rows_selected,]
      else experiments_filtered <- experiments
      experiments_selected <- experiments_filtered$`Experiment ID`
      # Filter the comparisons on the picked experiments
      comparisons_filtered <- comparisons %>% filter(`Experiment ID` %in% experiments_selected) %>%
        mutate(`TreatmentDose (diff)` = paste(`Treatment (diff)`, `Dose (diff)`),
               `TreatmentDose (ref)` = paste(`Treatment (ref)`, `Dose (ref)`))
      if(length(experiments_filtered$`Citation`) > 1) {
        citations_validated <- paste(gsub("<.+?>", "", experiments_filtered$`Citation`), collapse = "; ")
        citations_validated <- case_when(nchar(citations_validated) > 255 ~ paste0(sub("(^.{255}.+?);(.+)$", "\\1", citations_validated),
                                                                                   "... and ",
                                                                                   length(gregexpr(";", sub("(^.{255}.+?);(.+)$", "\\2", citations_validated))[[1]]),
                                                                                   " more."),
                                         T ~ paste0(sub("^(.+);(.+?)$", "\\1", citations_validated),
                                                    "; and ",
                                                    sub("^(.+);(.+?)$", "\\2", citations_validated)))
      } else {
        citations_validated <- gsub("<.+?>", "", experiments_filtered$`Citation`)
      }
      output$comparison_selector_experiments <- renderText(citations_validated)
      
      updateSelectizeInput(session,     'comparison_selector_library'            , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Library`)))))
      updateSelectizeInput(session,     'comparison_selector_organism'           , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Organism`)))))
      updateSelectizeInput(session,     'comparison_selector_source'             , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Source`)))))
      updateSelectizeInput(session,     'comparison_selector_citation'           , selected = "(any)", choices = c("(any)",             as.character(sort(unique(gsub("<.+?>", "", comparisons_filtered$`Citation`))))))
      updateSelectizeInput(session,     'comparison_selector_kind'               , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Kind`)))))
      updateSelectizeInput(session,     'comparison_selector_timepoint'          , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Timepoint`)))))
      updateSelectizeInput(session,     'comparison_selector_endpoint'           , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Endpoint`)))))
      updateSelectizeInput(session,     'comparison_selector_daysgrown_diff'     , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Days grown (diff)`)))))
      updateSelectizeInput(session,     'comparison_selector_daysgrown_ref'      , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Days grown (ref)`)))))
      updateSelectizeInput(session,     'comparison_selector_treatment_diff'     , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Treatment (diff)`)))))
      updateSelectizeInput(session,     'comparison_selector_treatment_ref'      , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Treatment (ref)`)))))
      updateSelectizeInput(session,     'comparison_selector_dose_diff'          , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Dose (diff)`)))))
      updateSelectizeInput(session,     'comparison_selector_dose_ref'           , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Dose (ref)`)))))
      updateSelectizeInput(session,     'comparison_selector_knockout_diff'      , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Knockout (diff)`)))))
      updateSelectizeInput(session,     'comparison_selector_knockout_ref'       , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Knockout (ref)`)))))
      updateSelectizeInput(session,     'comparison_selector_cell_line_diff'     , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Cell line (diff)`)))))
      updateSelectizeInput(session,     'comparison_selector_cell_line_ref'      , selected = "(any)", choices = c("(any)",             as.character(sort(unique(comparisons_filtered$`Cell line (ref)`)))))
      comparison_selector_clearable(F)
    }
  })
  # When the Select screens button is pressed in the comparisons searcher
  observeEvent(input$comparisons_selector_select, {
    # If fewer than 100 screens are validated
    if(length(comparisons_exported) < 100) {
      # Update the comparisons table picker selection
      experiments_selected <- experiments$`Experiment ID`[input$experiments_table_rows_selected]
      if(length(experiments_selected) > 0) comparisons_filtered <- comparisons %>% filter(`Experiment ID` %in% experiments_selected)
      else comparisons_filtered <- comparisons
      comparisons_table_rows_selected <- which(comparisons_filtered$`FriendlyID` %in% comparisons_exported)
      selectRows(comparisons_table_proxy, c(comparisons_table_rows_selected, input$comparisons_table_rows_selected))
    } else {
      output$comparisons_selector_list <- renderText(paste0("Won't select ", length(comparisons_exported), " comparisons. Please refine your search."))
    }
  })
  # When the "select all experiments" button is pressed
  observeEvent(input$experiments_table_selectAll, {
    # Select all visible rows in the experiments table
    selectRows(experiments_table_proxy, input$experiments_table_rows_all)
  })
  # When the "select none experiments" button is pressed
  observeEvent(input$experiments_table_selectNone, {
    # Clear the selection in the experiments table
    selectRows(experiments_table_proxy, NULL)
    # Remove the filter on the comparisons table
    output$comparisons_table <- renderDT(comparisons %>% select(-`Comparison ID`, -`Experiment ID`, -`FriendlyID`), filter = "top", server = TRUE, options = list(scrollX = T), escape = F)
    # Clear the comparisons text picker
    updateSelectizeInput(session, "comparisons_selected", server = TRUE, selected = NA, choices = comparisons$`FriendlyID`)
  })
  # Picked comparisons race condition logic
  comparisons_selectable <- reactiveVal(FALSE)
  # When comparisons are picked from the table
  observeEvent(input$comparisons_table_rows_selected, ignoreNULL = F, {
    experiments_selected <- experiments$`Experiment ID`[input$experiments_table_rows_selected]
    genes_preselected <- input$genes_selected
    if(length(experiments_selected) > 0) comparisons_filtered <- comparisons %>% filter(`Experiment ID` %in% experiments_selected)
    # If no comparisons are picked
    else {
      comparisons_filtered <- comparisons
      # Reset the gene list
      updateSelectizeInput(session, "genes_selected", server = TRUE, choices = genes, selected = genes_preselected, options = list(maxOptions = 1000))
    }
    comparisons_selected <- comparisons_filtered$`Comparison ID`[input$comparisons_table_rows_selected]
    comparisons_displayed <- comparisons_filtered$`FriendlyID`[input$comparisons_table_rows_selected]
    # Update the comparisons text picker selection
    if(comparisons_selectable()) {
      comparisons_selectable(FALSE)
      updateSelectizeInput(session, "comparisons_selected", selected = comparisons_displayed)
    } else {
      comparisons_selectable(TRUE)
    }
    # Update the gene list for genes in the comparisons selected
    if(length(comparisons_selected) > 0 & length(comparisons_selected) < 120) {
      q <- paste0("SELECT [gene_id] from stat WHERE [comparison_id] = ", paste0("'", comparisons_selected, "'", collapse = " OR [comparison_id] = "))
      genes_filtered <- bind_rows(lapply(cons, function(con) sqldf(q, connection = con))) %>% unique()
      genes_filtered <- unique(genes[genes %in% genes_filtered$gene_id])
      genes_filtered <<- genes_filtered
      updateSelectizeInput(session, "genes_selected", server = TRUE, choices = genes_filtered, selected = genes_preselected[genes_preselected %in% genes_filtered], options = list(maxOptions = 1000))
      # Deal when too many comparisons are selected
    } else {
      rm(genes_filtered)
      updateSelectizeInput(session, "genes_selected", server = TRUE, choices = genes, selected = genes_preselected, options = list(maxOptions = 1000))
    }
  })
  # When comparisons are picked from the text picker
  observeEvent(input$comparisons_selected, ignoreNULL = F, {
    experiments_selected <- experiments$`Experiment ID`[input$experiments_table_rows_selected]
    if(length(experiments_selected) > 0) comparisons_filtered <- comparisons %>% filter(`Experiment ID` %in% experiments_selected)
    else comparisons_filtered <- comparisons
    comparisons_table_rows_selected <- which(comparisons_filtered$`FriendlyID` %in% input$comparisons_selected)
    # Select those rows in the comparisons table
    if(comparisons_selectable()) {
      comparisons_selectable(FALSE)
      selectRows(comparisons_table_proxy, comparisons_table_rows_selected)
    } else {
      comparisons_selectable(TRUE)
    }
    # Update the individual comparison pickers
    updateSelectizeInput(session, "comparison_volcano", server = TRUE, choices = input$comparisons_selected)
    updateSelectizeInput(session, "comparison_rank", server = TRUE, choices = input$comparisons_selected)
    updateSelectizeInput(session, "comparison_enrichment", server = TRUE, choices = input$comparisons_selected)
    updateSelectizeInput(session, "comparison_X", server = TRUE, choices = input$comparisons_selected)
    if(length(input$comparisons_selected) > 1) {
      updateSelectizeInput(session, "comparison_Y", server = TRUE, choices = input$comparisons_selected, selected = input$comparisons_selected[2])
      updateSelectizeInput(session, "rocauc_head", server = TRUE, choices = input$comparisons_selected, selected = input$comparisons_selected[2:min(length(input$comparisons_selected),6)])
    } else {
      updateSelectizeInput(session, "comparison_Y", server = TRUE, choices = input$comparisons_selected)
    }
    updateSelectizeInput(session, "comparison_overlap", server = TRUE, choices = input$comparisons_selected, selected = input$comparisons_selected[1:min(length(input$comparisons_selected),7)])
    updateSelectizeInput(session, "rocauc_base", server = TRUE, choices = input$comparisons_selected, selected = input$comparisons_selected[1])
  })
  # When the "select all comparisons" button is pressed
  observeEvent(input$comparisons_table_selectAll, {
    # Select all visible rows in the comparisons table
    selectRows(comparisons_table_proxy, input$comparisons_table_rows_all)
  })
  # When the "select none experiments" button is pressed
  observeEvent(input$comparisons_table_selectNone, {
    # Clear the selection in the experiments table
    selectRows(comparisons_table_proxy, NULL)
    # Clear the comparisons text picker
    updateSelectizeInput(session, "comparisons_selected", selected = NA)
  })
  # When the volcano plot button is pressed
  observeEvent(input$volcano_submit, {
    volcano <- plotVolcano(input)
    output$volcano_table <- renderDT(volcano$plotdata, server = TRUE, filter = "top", options = list(scrollX = T), escape = F)
    output$volcano_plot <- renderUI({
      p = renderPlotly({ volcano$p })
      attr(p, 'outputArgs') <- list(height = volcano$height)
      return(p)
    })
    output$volcano_download_plot <- downloadHandler(
      filename = "export.html",
      content = function(f) htmlwidgets::saveWidget(volcano$p, f)
    )
    output$volcano_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(volcano$plotdata, f, sep = ",")
    )
    brush_table <<- volcano$plotdata %>% transmute(x = Score, y = FDR, Gene)
  })
  # When the rank plot button is pressed
  observeEvent(input$rank_submit, {
    rank <- plotRank(input)
    output$rank_table <- renderDT(rank$plotdata, server = TRUE, filter = "top", options = list(scrollX = T), escape = F)
    output$rank_plot <- renderUI({
      p = renderPlotly({ rank$p })
      attr(p, 'outputArgs') <- list(height = rank$height)
      return(p)
    })
    output$rank_download_plot <- downloadHandler(
      filename = "export.html",
      content = function(f) htmlwidgets::saveWidget(rank$p, f)
    )
    output$rank_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(rank$plotdata, f, sep = ",")
    )
    brush_table <<- rank$plotdata %>% transmute(x = Rank, y = Score, Gene)
  })
  # When the biplot button is pressed
  observeEvent(input$bi_submit, {
    biplot <- plotBiplot(input)
    output$bi_table <- renderDT(biplot$plotdata, server = TRUE, filter = "top", options = list(scrollX = T), escape = F)
    output$bi_plot <- renderUI({
      p = renderPlotly({ biplot$p })
      attr(p, 'outputArgs') <- list(height = biplot$height)
      return(p)
    })
    output$bi_download_plot <- downloadHandler(
      filename = "export.html",
      content = function(f) htmlwidgets::saveWidget(biplot$p, f)
    )
    output$bi_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(biplot$plotdata, f, sep = ",")
    )
    brush_table <<- biplot$plotdata %>% transmute(x = biplot$plotdata[[2]], y = biplot$plotdata[[4]], Gene)
  })
  # When the overlap button is pressed
  observeEvent(input$overlap_submit, {
    if (input$overlap_style == "Venn diagram") {
      overlap <- plotOverlap(input)
    } else if (input$overlap_style == "Upset plot") {
      overlap <- plotUpset(input)
    }
    output$overlap_table <- renderDT(overlap$plotdata, server = TRUE, filter = "top", options = list(scrollX = T), escape = F)
    output$overlap_plot <- renderUI({
      if (input$overlap_style == "Venn diagram") {
        p = renderPlotly({ overlap$p })
      } else if (input$overlap_style == "Upset plot") {
        p = renderPlot({ overlap$p })
      }
      attr(p, 'outputArgs') <- list(height = overlap$height)
      return(p)
    })
    output$overlap_download_plot <- downloadHandler(
      filename = "export.html",
      content = function(f) htmlwidgets::saveWidget(overlap$p, f)
    )
    output$overlap_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(overlap$plotdata, f, sep = ",")
    )
  })
  # When the ROC-AUC button is pressed
  observeEvent(input$rocauc_submit, {
    rocauc <- plotRoc(input)
    output$rocauc_table <- renderDT(rocauc$plotdata, server = TRUE, filter = "top", options = list(scrollX = T), escape = F)
    output$rocauc_plot <- renderUI({
      p = renderPlotly({ rocauc$p })
      attr(p, 'outputArgs') <- list(height = rocauc$height)
      return(p)
    })
    output$rocauc_download_plot <- downloadHandler(
      filename = "export.html",
      content = function(f) htmlwidgets::saveWidget(rocauc$p, f)
    )
    output$rocauc_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(rocauc$plotdata, f, sep = ",")
    )
  })
  #### Shiny server: react to input: correlate ####
  # When the validate genes button is pressed
  observeEvent(input$genes_queried_text_validate, {
    if(nchar(input$genes_queried_text) >= 3) {
      genes_typed <- unique(unlist(strsplit(input$genes_queried_text, split = "[,; ]")))
      genes_typed <- unique(genes[toupper(genes) %in% toupper(genes_typed)])
      if(length(genes_typed) <= 1000) {
        updateSelectizeInput(session, "genes_queried", selected = unique(c(input$genes_queried, genes_typed)), choices = genes, server = T)
      } else {
        updateSelectizeInput(session, "genes_queried", selected = unique(c(input$genes_queried, genes_typed)), choices = genes, server = F)
      }
    }
  })
  # When the core essentials link is pressed
  observeEvent(input$genes_queried_essentials, {
    updateTextAreaInput(session, "genes_queried_text", value = paste(essentials, collapse = ", "))
  })
  # When the copy genes from explore button is pressed
  observeEvent(input$genes_copy_from_explore, {
    updateSelectizeInput(session, "genes_queried", selected = input$genes_selected, choices = genes, server = T)
  })
  # When the clear genes button is pressed
  observeEvent(input$gq_clear_genes, {
    updateSelectizeInput(session, "genes_queried", selected = NULL, choices = genes, server = T)
  })
  # When the copy contrasts from explore button is pressed
  observeEvent(input$gq_filter_contrast_copy_from_explore, {
    updateSelectizeInput(session, "gq_filter_contrast", selected = input$comparisons_selected, server = TRUE, choices = comparisons$`FriendlyID`, options = list(maxOptions = 1000000))
  })
  # When the reset form button is pressed
  observeEvent(input$gq_reset, {
    initialiseCorrelatePickers()
  })
  # When the example button is pressed
  observeEvent(input$gq_example, {
    updateTextAreaInput(session, "genes_queried_text", value = "ABCC1 AMBRA1 CIP2A CUL3 DCLRE1C ERCC6L2 H2AX LIG4 MCPH1 NBN NHEJ1 RAD54L2 TDP2 TIAL1 UBA3 UBE2M UBE2K ZNF451 TP53")
    updateSelectizeInput(session, "genes_queried", selected = c("ABCC1", "AMBRA1", "CIP2A", "CUL3", "DCLRE1C", "ERCC6L2", "H2AX", "LIG4", "MCPH1", "NBN", "NHEJ1", "RAD54L2", "TDP2", "TIAL1", "UBA3", "UBE2M", "UBE2K", "ZNF451", "TP53"), choices = genes, server = T)
    updateSelectizeInput(session, 'gq_filter_citation', selected = NULL)
    updateSelectizeInput(session, 'gq_filter_library', selected = NULL)
    updateSelectizeInput(session, 'gq_filter_source', selected = names(cons))
    updateSelectizeInput(session, 'gq_filter_timepoint', selected = NULL)
    updateSelectizeInput(session, 'gq_filter_organism', selected = "Human")
    updateSelectizeInput(session, 'gq_filter_kind', selected = "Treatment")
    updateSelectizeInput(session, 'gq_filter_days_grown_diff', selected = NULL)
    updateSelectizeInput(session, 'gq_filter_treatment_diff', selected = "Etoposide")
    updateSelectizeInput(session, 'gq_filter_dose_diff', selected = NULL)
    updateSelectizeInput(session, 'gq_filter_knockout_diff', selected = NULL)
    updateSelectizeInput(session, 'gq_filter_cellline_diff', selected = NULL)
    updateSelectizeInput(session, 'gq_filter_days_grown_ref', selected = NULL)
    updateSelectizeInput(session, 'gq_filter_treatment_ref', selected = NULL)
    updateSelectizeInput(session, 'gq_filter_dose_ref', selected = NULL)
    updateSelectizeInput(session, 'gq_filter_knockout_ref', selected = NULL)
    updateSelectizeInput(session, 'gq_filter_cellline_ref', selected = NULL)
    updateSelectizeInput(session, 'gq_stat', selected = "Score")
    updateNumericInput(session, 'gq_cutoff', value = 1e-15)
    updateSelectizeInput(session, 'gq_filter_contrast', selected = NULL)
    updateTextInput(session, 'gq_filter_custom', value = "")
    updateSelectizeInput(session, 'analysis_method_gq', selected = "DrugZ")
  })
  
  initialiseCorrelateCustomise()
  # When the correlate customise plot button is pressed
  observeEvent(input$correlate_customise, {
    showModal(modalDialog(
      title = "Customise plot",
      easyClose = T,
      footer = NULL,
      numericInput('correlate_customise_height', "Plot height", min = 0, max = Inf, value = input$correlate_customise_height),
      numericInput('correlate_customise_width', "Plot width", min = 0, max = Inf, value = input$correlate_customise_width),
      checkboxInput('correlate_customise_width_auto', "Fit plot to window width", value = input$correlate_customise_width_auto),
      colourInput('correlate_customise_point_colour', "Point colour", value = input$correlate_customise_point_colour),
      selectizeInput('correlate_customise_genequery_fillscheme', "Colour violins by (gene query only)", choices = c("Treatment", "Knockout", "Cell line", "Library"), selected = input$correlate_customise_genequery_fillscheme),
      colourInput('correlate_customise_heatmap_high', "Positive correlation colour (heatmap only)", value = input$correlate_customise_heatmap_high),
      colourInput('correlate_customise_heatmap_mid', "Zero correlation colour (heatmap only)", value = input$correlate_customise_heatmap_mid),
      colourInput('correlate_customise_heatmap_low', "Negative correlation colour (heatmap only)", value = input$correlate_customise_heatmap_low),
      tags$hr(),
      actionButton('correlate_customise_reset', "Restore defaults", icon = icon("arrows-rotate")),
      modalButton("Save options", icon = icon("floppy-disk"))
    )
    )
    if(is.null(input$correlate_customise_height)) initialiseCorrelateCustomise()
    else {
      if(input$correlate_customise_width_auto) shinyjs::disable('correlate_customise_width')
    }
  })
  # When the customise plot restore defaults button is pressed
  observeEvent(input$correlate_customise_reset, {
    initialiseCorrelateCustomise()
  })
  # When the customise plot auto width option is selected
  observeEvent(input$correlate_customise_width_auto, {
    if(input$correlate_customise_width_auto) shinyjs::disable('correlate_customise_width')
    else shinyjs::enable('correlate_customise_width')
  })
  # When the gene query button is pressed
  observeEvent(input$genequery_submit, {
    # Plot gene(s)' scores across screens
    genequery <- genequery(input)
    output$genequery_plot <- renderPlotly({ genequery$p })
    output$genequery_table <- renderDT(genequery$displaydata, server = TRUE, filter = "top", options = list(scrollX = T), escape = F)
    output$genequery_download_plot <- downloadHandler(
      filename = "export.html",
      content = function(f) htmlwidgets::saveWidget(genequery$p, f)
    )
    output$genequery_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(genequery$displaydata, f, sep = ",")
    )
    brush_table <<- tibble(x = as.numeric(factor(genequery$displaydata$Contrast)), y = genequery$displaydata$Score, Gene = genequery$displaydata$Gene)
  })
  # When the hitmap button is pressed
  observeEvent(input$hitmap_submit, {
    # Plot gene(s)' scores across screens
    hitmap <- hitmap(input)
    output$hitmap_plot <- renderPlotly({ hitmap$p })
    output$hitmap_table <- renderDT(hitmap$displaydata, server = TRUE, filter = "top", options = list(scrollX = T), escape = F)
    output$hitmap_download_plot <- downloadHandler(
      filename = "export.html",
      content = function(f) htmlwidgets::saveWidget(hitmap$p, f)
    )
    output$hitmap_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(hitmap$displaydata, f, sep = ",")
    )
  })
  # When the heatmap button is pressed
  observeEvent(input$heatmap_submit, {
    heatmap <- plot_heatmap(input)
    output$heatmap_plot <- renderPlotly({ heatmap$plot })
    output$heatmap_table <- renderDT(heatmap$table, server = TRUE, filter = "top", options = list(scrollX = T), escape = F)
    output$heatmap_download_plot <- downloadHandler(
      filename = "export.html",
      content = function(f) htmlwidgets::saveWidget(heatmap$plot, f)
    )
    output$heatmap_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(heatmap$table, f, sep = ",")
    )
  } 
  )
  # When the UMAP button is pressed
  observeEvent(input$umap_submit, {
    umap <- plot_umap(input)
    output$umap_plot <- renderPlotly({ umap$plot })
    output$umap_table <- renderDT(umap$table, server = TRUE, filter = "top", options = list(scrollX = T), escape = F)
    output$umap_download_plot <- downloadHandler(
      filename = "export.html",
      content = function(f) htmlwidgets::saveWidget(umap$plot, f)
    )
    output$umap_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(umap$table, f, sep = ",")
    )
    brush_table <<- tibble(x = umap$table$UMAP1, y = umap$table$UMAP2, Gene = umap$table$Gene)
  })
  # When the network button is pressed
  observeEvent(input$network_submit, {
    network <- plot_network(input)
    output$network_plot <- renderVisNetwork({ network$plot })
    output$network_table <- renderDT(network$table, server = TRUE, filter = "top", options = list(scrollX = T), escape = F)
    output$network_download_plot <- downloadHandler(
      filename = "export.html",
      content = function(f) htmlwidgets::saveWidget(network$plot, f)
    )
    output$network_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(network$table, f, sep = ",")
    )
  })
  # When the enrichment analysis button is pressed
  observeEvent(input$enrichment_submit, {
    enrichment <- plotEnrichment(input)
    output$enrichment_table <- renderDT(enrichment$plotdata, server = TRUE, filter = "top", options = list(scrollX = T), escape = F)
    output$enrichment_plot <- renderPlotly({ enrichment$p })
    output$enrichment_download_plot <- downloadHandler(
      filename = "export.html",
      content = function(f) htmlwidgets::saveWidget(enrichment$p, f)
    )
    output$enrichment_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(enrichment$plotdata, f, sep = ",")
    )
  })
  # When the Pendragonator button is pressed
  observeEvent(input$pendragonator_submit, {
    pendragonated <- pendragonate(input)
    output$pendragonator_table <- renderDT(pendragonated, server = T, filter = "top", options = list(scrollX = T), escape = F)
    output$pendragonator_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(pendragonated, f, sep = ",")
    )
  })
  # When the bulk download button is pressed
  observeEvent(input$download_submit, {
    download <- bulkDownload(input)
    output$download_console <- renderText(download$console)
    output$download_table <- renderDT(download$table %>% head(6), server = T, filter = "top", options = list(scrollX = T), escape = F)
    output$download_download <- downloadHandler(
      filename = download$filename,
      content = function(f) fwrite(download$table, f, sep = download$delim)
    )
  })
  #### Shiny server: react on manipulating the plot ####
  # When scatter plots are brushed
  observeEvent(input$`plotly_brushed-A`, {
    if(exists("brush_table")) {
      brushed_x <- as.numeric(unlist(strsplit(sub("^.+x\\\":\\[(.+?)\\],\\\"y.+$", "\\1", input$`plotly_brushed-A`), ",")))
      brushed_y <- as.numeric(unlist(strsplit(sub("^.+y\\\":\\[(.+?)\\].+$", "\\1", input$`plotly_brushed-A`), ",")))
      brushed <- brush_table %>% filter(x >= brushed_x[1] & x <= brushed_x[2] & y >= brushed_y[1] & y <= brushed_y[2])
      if(input$main == "Explore") {
        if(exists("genes_filtered")) {
          genes_check <- genes_filtered
        } else {
          genes_check <- genes
        }
        updateSelectizeInput(session, "genes_selected", selected = c(input$genes_selected, brushed$Gene), choices = genes_check, server = T)
      } else if(input$main == "Correlate") {
        updateSelectizeInput(session, "genes_queried", selected = c(input$genes_queried, brushed$Gene), choices = genes, server = T)
      }
    }
  })
  
  #### Shiny server: react to input: guide ####
  # When the guide tab is selected
  observeEvent(input$main, {
    if(input$main == "Guides") {
      if(!exists("guideLibrary")) {
        initialiseGuideLibraryTable()
        initialiseGuideLibraryPickers()
      }
    }
  })
  #### Shiny server: react to input: library ####
  # When the library tab is selected
  observeEvent(input$main, {
    if(input$main == "Libraries") {
      initialiseLibrariesTable()
    }
  })
  
  # When the ontology tab is selected
  observeEvent(input$main, {
    if(input$main == "Ontology") {
      initialiseOntologyTable()
    }
  })
  
  # When the guide library example button is clicked
  observeEvent(input$guideLibrary_example, {
    initialiseGuideLibraryPickers()
    updateCheckboxGroupInput(session, "guideLibrary_geneclasses", selected = c("Protein-coding"))
    updateCheckboxGroupInput(session, "guideLibrary_chromosomes", selected = c("10"))
    updateCheckboxGroupInput(session, "guideLibrary_chemistry", selected = c("Knockout"))
    updateCheckboxGroupInput(session, "guideLibrary_organism", selected = c("Human"))
  })
  # When the guide library reset button is clicked
  observeEvent(input$guideLibrary_reset, {
    initialiseGuideLibraryPickers()
  })
  # When the guide library button is clicked
  observeEvent(input$guideLibrary_submit, {
    if(!exists("guideLibrary")) return()   # no guides dataset loaded
    # Call main function
    guideLibraryFiltered <- guideLibraryinate(input)
    # Draw table
    output$guideLibrary_table <- renderDT(guideLibraryFiltered, filter = "top", server = T,
                                          options = list(iDisplayLength = 100, scrollX = T), escape = F)
    output$guideLibrary_download <- downloadHandler(
      filename = "export.csv",
      content = function(f) fwrite(guideLibraryFiltered, f, sep = ",")
    )
  })
  
  #### Shiny server: react to input: exorcise ####
  # When the exorcise example button is pressed
  observeEvent(input$exorcise_example, {
    updateTextAreaInput(session, 'exorcise_seq', value = "AAGGAGCCAACATAACAGAT, AGAAACCTACAACTCATGGA, GGCTCAGGGTTACCGAAGAG, TGGTACTGATTATGGCACTC, CAACGTCGCGAACGTCGTAT, ACCACCGTTCGTACCGGTCG, AGACCGTCGTCGATCGATAC, GTCGTACGGATTCGCGCGTA")
    updateTextInput(session, 'exorcise_pam', value = "NGG")
    updateSelectizeInput(session, 'exorcise_genome', selected = "GRCh38 Ensembl 111")
    updateSelectizeInput(session, 'exorcise_mode', selected = "Knockout")
    updateTextInput(session, 'exorcise_mode_advanced', value = NULL)
    updateTextAreaInput(session, 'exorcise_orig', value = "BRCA1, BRCA1, BRCA1, BRCA1, non-targeting, non-targeting, non-targeting, non-targeting")
    updateTextInput(session, 'exorcise_control_string', value = "non-targeting")
  })
  # When the exorcise reset button is pressed
  observeEvent(input$exorcise_reset, {
    updateTextAreaInput(session, 'exorcise_seq', value = "")
    updateTextInput(session, 'exorcise_pam', value = "NGG")
    updateTextInput(session, 'exorcise_mode_advanced', value = "")
    updateTextAreaInput(session, 'exorcise_orig', value = "")
    updateTextInput(session, 'exorcise_control_string', value = "")
  })
  # When the exorcise button is pressed
  observeEvent(input$exorcise_submit, {
    exorcised <- run_exorcise(input)
    if(length(exorcised$exorcise_console) > 0) output$exorcise_console <- exorcised$exorcise_console
    if(length(exorcised$exorcise_download_input) > 0) output$exorcise_download_input <- exorcised$exorcise_download_input
    if(length(exorcised$exorcise_download_result) > 0) output$exorcise_download_result <- exorcised$exorcise_download_result
  })
  
  #### Shiny server: react to input: save/load
  # When the save/load button is pressed
  saveload_reload <- reactiveVal(0)
  observeEvent(input$main, {
    if (input$main == "Save/Load") {
      # Return the user to their original tab
      updateTabsetPanel(session, "main", selected = lastTab())
      # Show the dialog box
      showModal(modalDialog(
        title = "Save/Load",
        easyClose = T,
        footer = NULL,
        tags$h4("Save"),
        tags$p("Download the state of your current CRAVE session as a key file."),
        downloadButton('saveload_key_download', label = "Download key", class = "btn-success"),
        tags$br(),
        tags$hr(),
        tags$h4("Load"),
        tags$p("Restore the state of a previous CRAVE session by uploading a key file."),
        fileInput('saveload_key_upload', label = "Upload key", accept = ".txt"),
        textOutput('saveload_console'),
        tags$br(),
        actionButton("saveload_load", label = "Load key", icon = icon("folder-open"), class = "btn-success"),
        tags$hr(),
        modalButton("Close", icon = icon("xmark"))
      ))
      if(!is.null(input$saveload_key_upload$datapath)) {
        consoleText <- paste0("Existing key already uploaded: ", input$saveload_key_upload$name, ".")
      } else {
        consoleText <- "No key uploaded."
      }
      output$saveload_console <- renderText( { consoleText } )
      saveload_reload(8)
      # Calculate key
      key <- list()
      # Explore
      key$genes_selected_text <- input$genes_selected_text
      key$genes_selected <- input$genes_selected
      key$comparisons_table_rows_selected <- input$comparisons_table_rows_selected
      key$experiments_table_rows_selected <- input$experiments_table_rows_selected
      key$comparisons_selected <- input$comparisons_selected
      key$comparison_selector_citation <- input$comparison_selector_citation
      key$comparison_selector_kind <- input$comparison_selector_kind
      key$comparison_selector_library <- input$comparison_selector_library
      key$comparison_selector_timepoint <- input$comparison_selector_timepoint
      key$comparison_selector_endpoint <- input$comparison_selector_endpoint
      key$comparison_selector_organism <- input$comparison_selector_organism
      key$comparison_selector_source <- input$comparison_selector_source
      key$comparison_selector_daysgrown_diff <- input$comparison_selector_daysgrown_diff
      key$comparison_selector_daysgrown_ref <- input$comparison_selector_daysgrown_ref
      key$comparison_selector_treatment_diff <- input$comparison_selector_treatment_diff
      key$comparison_selector_treatment_ref <- input$comparison_selector_treatment_ref
      key$comparison_selector_dose_diff <- input$comparison_selector_dose_diff
      key$comparison_selector_dose_ref <- input$comparison_selector_dose_ref
      key$comparison_selector_knockout_diff <- input$comparison_selector_knockout_diff
      key$comparison_selector_knockout_ref <- input$comparison_selector_knockout_ref
      key$comparison_selector_cell_line_diff <- input$comparison_selector_cell_line_diff
      key$comparison_selector_cell_line_ref <- input$comparison_selector_cell_line_ref
      key$comparison_volcano <- input$comparison_volcano
      key$analysis_method_se_v <- input$analysis_method_se_v
      key$comparison_rank <- input$comparison_rank
      key$analysis_method_se_r <- input$analysis_method_se_r
      key$comparison_X <- input$comparison_X
      key$analysis_method_se_b_x <- input$analysis_method_se_b_x
      key$comparison_Y <- input$comparison_Y
      key$analysis_method_se_b_y <- input$analysis_method_se_b_y
      key$comparison_overlap <- input$comparison_overlap
      key$overlap_cutoff_method <- input$overlap_cutoff_method
      key$overlap_cutoff_stat <- input$overlap_cutoff_stat
      key$overlap_cutoff_value <- input$overlap_cutoff_value
      key$overlap_cutoff_direction <- input$overlap_cutoff_direction
      key$rocauc_base <- input$rocauc_base
      key$rocauc_method <- input$rocauc_method
      key$rocauc_cutoff_stat <- input$rocauc_cutoff_stat
      key$rocauc_cutoff_value <- input$rocauc_cutoff_value
      key$rocauc_cutoff_direction <- input$rocauc_cutoff_direction
      key$rocauc_head <- input$rocauc_head
      # Correlate
      key$genes_queried_text <- input$genes_queried_text
      key$genes_queried <- input$genes_queried
      key$gq_stat <- input$gq_stat
      key$gq_cutoff <- input$gq_cutoff
      key$analysis_method_gq <- input$analysis_method_gq
      key$gq_filter_citation <- input$gq_filter_citation
      key$gq_filter_kind <- input$gq_filter_kind
      key$gq_filter_library <- input$gq_filter_library
      key$gq_filter_timepoint <- input$gq_filter_timepoint
      key$gq_filter_organism <- input$gq_filter_organism
      key$gq_filter_source <- input$gq_filter_source
      key$gq_filter_days_grown_diff <- input$gq_filter_days_grown_diff
      key$gq_filter_days_grown_ref <- input$gq_filter_days_grown_ref
      key$gq_filter_treatment_diff <- input$gq_filter_treatment_diff
      key$gq_filter_treatment_ref <- input$gq_filter_treatment_ref
      key$gq_filter_dose_diff <- input$gq_filter_dose_diff
      key$gq_filter_dose_ref <- input$gq_filter_dose_ref
      key$gq_filter_knockout_diff <- input$gq_filter_knockout_diff
      key$gq_filter_knockout_ref <- input$gq_filter_knockout_ref
      key$gq_filter_cellline_diff <- input$gq_filter_cellline_diff
      key$gq_filter_cellline_ref <- input$gq_filter_cellline_ref
      key$gq_filter_contrast <- input$gq_filter_contrast
      key$gq_filter_custom <- input$gq_filter_custom
      key$genequery_violin <- input$genequery_violin
      key$corr_cutoff <- input$corr_cutoff
      key$enrichment_cutoff_stat <- input$enrichment_cutoff_stat
      key$enrichment_cutoff_value <- input$enrichment_cutoff_value
      key$enrichment_cutoff_direction <- input$enrichment_cutoff_direction
      key$enrichment_cutoff_tail <- input$enrichment_cutoff_tail
      key$pendragonator_n <- input$pendragonator_n
      key$pendragonator_stat <- input$pendragonator_stat
      key$download_filename <- input$download_filename
      # Exorcise
      key$exorcise_seq <- input$exorcise_seq
      key$exorcise_orig <- input$exorcise_orig
      key$exorcise_control_string <- input$exorcise_control_string
      key$exorcise_pam <- input$exorcise_pam
      key$exorcise_genome <- input$exorcise_genome
      key$exorcise_mode <- input$exorcise_mode
      key$exorcise_mode_advanced <- input$exorcise_mode_advanced
      # Guides
      key$guideLibrary_targets <- input$guideLibrary_targets
      key$guideLibrary_geneclasses <- input$guideLibrary_geneclasses
      key$guideLibrary_chromosomes <- input$guideLibrary_chromosomes
      key$guideLibrary_chemistry <- input$guideLibrary_chemistry
      key$guideLibrary_assembly <- input$guideLibrary_assembly
      key$guideLibrary_pam <- input$guideLibrary_pam
      key$guideLibrary_library <- input$guideLibrary_library
      key$guideLibrary_organism <- input$guideLibrary_organism
      output$saveload_key_download <- downloadHandler(
        filename = paste0("crave-key.txt"),
        content = function(f) save(key, file = f)
      )
    }})
  # When a key is uploaded
  observeEvent(input$saveload_key_upload, ignoreNULL = F, {
    if(!is.null(input$saveload_key_upload$datapath)) {
      consoleText <- paste0("Key uploaded: ", input$saveload_key_upload$name, ".")
    } else {
      consoleText <- "No key uploaded."
    }
    output$saveload_console <- renderText( { consoleText } )
  })
  # When the load button is pressed
  observeEvent(input$saveload_load, {
    if(!is.null(input$saveload_key_upload)) {
      base::load(input$saveload_key_upload$datapath[1])
      # Clear values
      if (saveload_reload() == 8) {
        # Explore
        updateTextAreaInput(session, 'genes_selected_text', value = NULL)
        updateSelectizeInput(session, 'genes_selected', selected = NULL)
        selectRows(experiments_table_proxy, NULL)
        selectRows(comparisons_table_proxy, NULL)
        updateSelectizeInput(session, 'comparisons_selected', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_citation', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_kind', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_library', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_timepoint', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_endpoint', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_organism', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_source', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_daysgrown_diff', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_daysgrown_ref', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_treatment_diff', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_treatment_ref', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_dose_diff', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_dose_ref', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_knockout_diff', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_knockout_ref', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_cell_line_diff', selected = NULL)
        updateSelectizeInput(session, 'comparison_selector_cell_line_ref', selected = NULL)
        updateSelectizeInput(session, 'comparison_volcano', selected = NULL)
        updateSelectizeInput(session, 'analysis_method_se_v', selected = NULL)
        updateSelectizeInput(session, 'comparison_rank', selected = NULL)
        updateSelectizeInput(session, 'analysis_method_se_r', selected = NULL)
        updateSelectizeInput(session, 'comparison_X', selected = NULL)
        updateSelectizeInput(session, 'analysis_method_se_b_x', selected = NULL)
        updateSelectizeInput(session, 'comparison_Y', selected = NULL)
        updateSelectizeInput(session, 'analysis_method_se_b_y', selected = NULL)
        updateSelectizeInput(session, 'comparison_overlap', selected = NULL)
        updateSelectizeInput(session, 'overlap_cutoff_method', selected = NULL)
        updateSelectizeInput(session, 'overlap_cutoff_stat', selected = NULL)
        updateSelectizeInput(session, 'overlap_cutoff_value', selected = NULL)
        updateSelectizeInput(session, 'overlap_cutoff_direction', selected = NULL)
        updateSelectizeInput(session, 'rocauc_base', selected = NULL)
        updateSelectizeInput(session, 'rocauc_method', selected = NULL)
        updateSelectizeInput(session, 'rocauc_cutoff_stat', selected = NULL)
        updateSelectizeInput(session, 'rocauc_cutoff_value', selected = NULL)
        updateSelectizeInput(session, 'rocauc_cutoff_direction', selected = NULL)
        updateSelectizeInput(session, 'rocauc_head', selected = NULL)
        # Correlate
        updateSelectizeInput(session, 'genes_queried_text', selected = NULL)
        updateSelectizeInput(session, 'genes_queried', selected = NULL)
        updateSelectizeInput(session, 'gq_stat', selected = NULL)
        updateNumericInput(session, 'gq_cutoff', value = NULL)
        updateSelectizeInput(session, 'analysis_method_gq', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_citation', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_kind', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_library', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_timepoint', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_organism', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_source', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_days_grown_diff', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_days_grown_ref', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_treatment_diff', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_treatment_ref', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_dose_diff', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_dose_ref', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_knockout_diff', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_knockout_ref', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_cellline_diff', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_cellline_ref', selected = NULL)
        updateSelectizeInput(session, 'gq_filter_contrast', selected = NULL)
        updateTextInput(session, 'gq_filter_custom', value = NULL)
        updateSelectizeInput(session, 'genequery_violin', selected = NULL)
        updateSliderTextInput(session, 'corr_cutoff', selected = NULL)
        updateSelectizeInput(session, 'enrichment_cutoff_stat', selected = NULL)
        updateSelectizeInput(session, 'enrichment_cutoff_value', selected = NULL)
        updateSelectizeInput(session, 'enrichment_cutoff_direction', selected = NULL)
        updateSelectizeInput(session, 'enrichment_cutoff_tail', selected = NULL)
        updateSelectizeInput(session, 'pendragonator_n', selected = NULL)
        updateSelectizeInput(session, 'pendragonator_stat', selected = NULL)
        updateTextInput(session, 'download_filename', value = NULL)
        # Exorcise
        updateTextAreaInput(session, 'exorcise_seq', value = NULL)
        updateTextAreaInput(session, 'exorcise_orig', value = NULL)
        updateTextInput(session, 'exorcise_control_string', value = NULL)
        updateTextInput(session, 'exorcise_pam', value = NULL)
        updateSelectizeInput(session, 'exorcise_genome', selected = NULL)
        updateSelectizeInput(session, 'exorcise_mode', selected = NULL)
        updateTextInput(session, 'exorcise_mode_advanced', value = NULL)
        # Guides
        updateSelectizeInput(session, 'guideLibrary_targets', selected = NULL)
        updateCheckboxGroupInput(session, 'guideLibrary_geneclasses', selected = NULL)
        updateCheckboxGroupInput(session, 'guideLibrary_chromosomes', selected = NULL)
        updateCheckboxGroupInput(session, 'guideLibrary_chemistry', selected = NULL)
        updateSelectizeInput(session, 'guideLibrary_assembly', selected = NULL)
        updateSelectizeInput(session, 'guideLibrary_pam', selected = NULL)
        updateSelectizeInput(session, 'guideLibrary_library', selected = NULL)
        updateSelectizeInput(session, 'guideLibrary_organism', selected = NULL)
      }
      # Restore values
      comparisons_selectable(FALSE)
      # Explore
      updateTextAreaInput(session, 'genes_selected_text', value = key$genes_selected_text)
      if (saveload_reload() == 0) shinyjs::click('genes_selected_text_validate')
      updateSelectizeInput(session, 'genes_selected', selected = key$genes_selected)
      selectRows(experiments_table_proxy, key$experiments_table_rows_selected)
      selectRows(comparisons_table_proxy, key$comparisons_table_rows_selected)
      updateSelectizeInput(session, 'comparisons_selected', selected = key$comparisons_selected)
      updateSelectizeInput(session, 'comparison_selector_citation', selected = key$comparison_selector_citation)
      updateSelectizeInput(session, 'comparison_selector_kind', selected = key$comparison_selector_kind)
      updateSelectizeInput(session, 'comparison_selector_library', selected = key$comparison_selector_library)
      updateSelectizeInput(session, 'comparison_selector_timepoint', selected = key$comparison_selector_timepoint)
      updateSelectizeInput(session, 'comparison_selector_endpoint', selected = key$comparison_selector_end)
      updateSelectizeInput(session, 'comparison_selector_organism', selected = key$comparison_selector_organism)
      updateSelectizeInput(session, 'comparison_selector_source', selected = key$comparison_selector_source)
      updateSelectizeInput(session, 'comparison_selector_daysgrown_diff', selected = key$comparison_selector_daysgrown_diff)
      updateSelectizeInput(session, 'comparison_selector_daysgrown_ref', selected = key$comparison_selector_daysgrown_ref)
      updateSelectizeInput(session, 'comparison_selector_treatment_diff', selected = key$comparison_selector_treatment_diff)
      updateSelectizeInput(session, 'comparison_selector_treatment_ref', selected = key$comparison_selector_treatment_ref)
      updateSelectizeInput(session, 'comparison_selector_dose_diff', selected = key$comparison_selector_dose_diff)
      updateSelectizeInput(session, 'comparison_selector_dose_ref', selected = key$comparison_selector_dose_ref)
      updateSelectizeInput(session, 'comparison_selector_knockout_diff', selected = key$comparison_selector_knockout_diff)
      updateSelectizeInput(session, 'comparison_selector_knockout_ref', selected = key$comparison_selector_knockout_ref)
      updateSelectizeInput(session, 'comparison_selector_cell_line_diff', selected = key$comparison_selector_cell_line_diff)
      updateSelectizeInput(session, 'comparison_selector_cell_line_ref', selected = key$comparison_selector_cell_line_ref)
      updateSelectizeInput(session, 'comparison_volcano', selected = key$comparison_volcano)
      updateSelectizeInput(session, 'analysis_method_se_v', selected = key$analysis_method_se_v)
      updateSelectizeInput(session, 'comparison_rank', selected = key$comparison_rank)
      updateSelectizeInput(session, 'analysis_method_se_r', selected = key$analysis_method_se_r)
      updateSelectizeInput(session, 'comparison_X', selected = key$comparison_X)
      updateSelectizeInput(session, 'analysis_method_se_b_x', selected = key$analysis_method_se_b_x)
      updateSelectizeInput(session, 'comparison_Y', selected = key$comparison_Y)
      updateSelectizeInput(session, 'analysis_method_se_b_y', selected = key$analysis_method_se_b_y)
      updateSelectizeInput(session, 'comparison_overlap', selected = key$comparison_overlap)
      updateSelectizeInput(session, 'overlap_cutoff_method', selected = key$overlap_cutoff_method)
      updateSelectizeInput(session, 'overlap_cutoff_stat', selected = key$overlap_cutoff_stat)
      updateSelectizeInput(session, 'overlap_cutoff_value', selected = key$overlap_cutoff_value)
      updateSelectizeInput(session, 'overlap_cutoff_direction', selected = key$overlap_cutoff_direction)
      updateSelectizeInput(session, 'rocauc_base', selected = key$rocauc_base)
      updateSelectizeInput(session, 'rocauc_method', selected = key$rocauc_method)
      updateSelectizeInput(session, 'rocauc_cutoff_stat', selected = key$rocauc_cutoff_stat)
      updateSelectizeInput(session, 'rocauc_cutoff_value', selected = key$rocauc_cutoff_value)
      updateSelectizeInput(session, 'rocauc_cutoff_direction', selected = key$rocauc_cutoff_direction)
      updateSelectizeInput(session, 'rocauc_head', selected = key$rocauc_head)
      # Correlate
      updateSelectizeInput(session, 'genes_queried_text', selected = key$genes_queried_text)
      if (saveload_reload() == 0) shinyjs::click('genes_queried_text_validate')
      updateSelectizeInput(session, 'genes_queried', selected = key$genes_queried)
      updateSelectizeInput(session, 'gq_stat', selected = key$gq_stat)
      updateNumericInput(session, 'gq_cutoff', value = key$gq_cutoff)
      updateSelectizeInput(session, 'analysis_method_gq', selected = key$analysis_method_gq)
      updateSelectizeInput(session, 'gq_filter_citation', selected = key$gq_filter_citation)
      updateSelectizeInput(session, 'gq_filter_kind', selected = key$gq_filter_kind)
      updateSelectizeInput(session, 'gq_filter_library', selected = key$gq_filter_library)
      updateSelectizeInput(session, 'gq_filter_timepoint', selected = key$gq_filter_timepoint)
      updateSelectizeInput(session, 'gq_filter_organism', selected = key$gq_filter_organism)
      updateSelectizeInput(session, 'gq_filter_source', selected = key$gq_filter_source)
      updateSelectizeInput(session, 'gq_filter_days_grown_diff', selected = key$gq_filter_days_grown_diff)
      updateSelectizeInput(session, 'gq_filter_days_grown_ref', selected = key$gq_filter_days_grown_ref)
      updateSelectizeInput(session, 'gq_filter_treatment_diff', selected = key$gq_filter_treatment_diff)
      updateSelectizeInput(session, 'gq_filter_treatment_ref', selected = key$gq_filter_treatment_ref)
      updateSelectizeInput(session, 'gq_filter_dose_diff', selected = key$gq_filter_dose_diff)
      updateSelectizeInput(session, 'gq_filter_dose_ref', selected = key$gq_filter_dose_ref)
      updateSelectizeInput(session, 'gq_filter_knockout_diff', selected = key$gq_filter_knockout_diff)
      updateSelectizeInput(session, 'gq_filter_knockout_ref', selected = key$gq_filter_knockout_ref)
      updateSelectizeInput(session, 'gq_filter_cellline_diff', selected = key$gq_filter_cellline_diff)
      updateSelectizeInput(session, 'gq_filter_cellline_ref', selected = key$gq_filter_cellline_ref)
      updateSelectizeInput(session, 'gq_filter_contrast', selected = key$gq_filter_contrast)
      updateTextInput(session, 'gq_filter_custom', value = key$gq_filter_custom)
      updateSelectizeInput(session, 'genequery_violin', selected = key$genequery_violin)
      updateSliderTextInput(session, 'corr_cutoff', selected = key$corr_cutoff)
      updateSelectizeInput(session, 'enrichment_cutoff_stat', selected = key$enrichment_cutoff_stat)
      updateSelectizeInput(session, 'enrichment_cutoff_value', selected = key$enrichment_cutoff_value)
      updateSelectizeInput(session, 'enrichment_cutoff_direction', selected = key$enrichment_cutoff_direction)
      updateSelectizeInput(session, 'enrichment_cutoff_tail', selected = key$enrichment_cutoff_tail)
      updateSelectizeInput(session, 'pendragonator_n', selected = key$pendragonator_n)
      updateSelectizeInput(session, 'pendragonator_stat', selected = key$pendragonator_stat)
      updateTextInput(session, 'download_filename', value = key$download_filename)
      # Exorcise
      updateTextAreaInput(session, 'exorcise_seq', value = key$exorcise_seq)
      updateTextAreaInput(session, 'exorcise_orig', value = key$exorcise_orig)
      updateTextInput(session, 'exorcise_control_string', value = key$exorcise_control_string)
      updateTextInput(session, 'exorcise_pam', value = key$exorcise_pam)
      updateSelectizeInput(session, 'exorcise_genome', selected = key$exorcise_genome)
      updateSelectizeInput(session, 'exorcise_mode', selected = key$exorcise_mode)
      updateTextInput(session, 'exorcise_mode_advanced', value = key$exorcise_mode_advanced)
      # Guides
      updateSelectizeInput(session, 'guideLibrary_targets', selected = key$guideLibrary_targets)
      updateCheckboxGroupInput(session, 'guideLibrary_geneclasses', selected = key$guideLibrary_geneclasses)
      updateCheckboxGroupInput(session, 'guideLibrary_chromosomes', selected = key$guideLibrary_chromosomes)
      updateCheckboxGroupInput(session, 'guideLibrary_chemistry', selected = key$guideLibrary_chemistry)
      updateSelectizeInput(session, 'guideLibrary_assembly', selected = key$guideLibrary_assembly)
      updateSelectizeInput(session, 'guideLibrary_pam', selected = key$guideLibrary_pam)
      updateSelectizeInput(session, 'guideLibrary_library', selected = key$guideLibrary_library)
      updateSelectizeInput(session, 'guideLibrary_organism', selected = key$guideLibrary_organism)
      # Return
      # Reapply due to reactive values falling through
      if (saveload_reload() > 0) {
        shinyjs::click('saveload_load')
        saveload_reload(saveload_reload() - 1)
        output$saveload_console <- renderText( { paste0("Please wait... loading state from: ", input$saveload_key_upload$name, ".") } )
      } else {
        comparisons_selectable(TRUE)
        output$saveload_console <- renderText( { paste0("State successfully loaded from key: ", input$saveload_key_upload$name) } )
      }
    } else {
      output$saveload_console <- renderText( { "Please upload a key." } )
    }
  }) 
  
  #### Logging ####
  log_shiny_input_changes(input)
  
  #### MOTD ####
  # MOTDs
  motds <- c(
    "Vibe check passed.",
    "Just one more clonogenic assay.",
    "Don't be afraid to struggle with R alone.",
    "ERCCCCCCCCCCCC6L2.",
    "Pink is better than orange.",
    "Your MOTD here!",
    "*hides in the supplementaries*",
    "When the reviewers' report is longer than some papers.",
    "One of these days I'll write the methods section.",
    "I'm having a BLAST writing these.",
    "I hope UBE OK. Have the BEst day ever.",
    "Take time to unwind. Relieve topological stress.",
    "Bayesian is bae. Log2(fold change) my mind.",
    "Add some MAGeCK to your day today.",
    "Don't forget to enrich. Adopt a positive growth phenotype.",
    "IC50 reasons to drop out all hypersensitivities.",
    "You are statistically significant.",
    "Exorcise your off-targets away.",
    "That's hypergeometric!",
    "Science is like magic, but real.",
    "I didn't CRISPick these genes. I inherited them.",
    "All models are wrong, but some are useful.",
    "May these results validate your validations.",
    "My idea of an upper-tail t-test is the first sip of an Earl Grey.",
    "One way or ANOVA, I'll test the population standard deviation.",
    "Next slide, please.",
    "This better be worth the eight hours I spent in TC.",
    "Warning: this product might contain homoscedastic stochasticity.",
    "I map, UMAP, we all MAP kinase.",
    "Take time to relax. Take a well-deserved double-strand break.",
    "Unclench. Helicase your negative supercoils away.",
    "May the bonds we make be covalent.",
    "May the path be lit with GFP and luciferase.",
    "While you're SLFN11, I'm in POLE position TRAIPsing towards the sticky end.",
    "Stand at the precipice, this 5' overhang, and behold the genome-wide view around you.",
    "Everyone say, \"single-ended double-strand break\"!",
    "DDR but it stands for Dance Dance Revolution.",
    "Do androids DDREAMM of electric sheep?",
    "What the fucose kinase are these results.",
    "You deserve a break, but if there's not enough time, then take a WEE1.",
    "Integrate your stress response. Translate your energy.",
    "Do you take your coffee full fat or DCAF12L1?",
    "Gen AI is not a substitute for thinking clearly.",
    "Feel like you're on TOP1 of the world.",
    "Back in the carbonyl cyanide m-chlorophenyl hydrazone."
  )
  set.seed(round(as.numeric(Sys.time())/3600))
  motd <- sample(motds, 1)
  output$motd <- renderText({ motd })
}
output$motd <- renderText({ motd })

  manageServer(input, output, session)
}
