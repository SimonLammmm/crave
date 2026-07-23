#### I/O ####
loadOne <- function(ds) {
  # ds: one entry from `datasets` (config.R) —
  #     list(name, path, guides = FALSE, citation_from_id = FALSE)
  path_dataset <- ds$path
  
  # Optional guide library / libraries (only if the dataset ships them)
  if (isTRUE(ds$guides)) {
    f_guide <- file.path(path_dataset, "guideLibrary.db")
    f_libs  <- file.path(path_dataset, "libraries.txt.gz")
    if (file.exists(f_guide)) file_guideLibrary <<- f_guide
    if (file.exists(f_libs)) file_libraries <<- f_libs
  }
  
  # I/O
  file_experiments <- file.path(path_dataset, "experiments_metadata.csv.gz")
  file_comparisons <- file.path(path_dataset, "comparisons_metadata.csv.gz")
  file_data        <- file.path(path_dataset, "database.db")
  file_ontology    <- file.path(path_dataset, "ontology.tsv.gz")
  
  # Return without crashing if data do not exist
  if(!file.exists(file_data)) {
    warning("Dataset ", file_data, " doesn't exist.")
    return()
  } else {
    con <- dbConnect(SQLite(), file_data)
    experiments <- fread(file_experiments, colClasses = "character")
    comparisons <- fread(file_comparisons, colClasses = "character")
    # Optional ontology
    ontology <- if (file.exists(file_ontology))
      fread(file_ontology, colClasses = "character") else data.table()
    # Wrangle experiments metadata
    if (isTRUE(ds$citation_from_id)) {
      experiments <- experiments %>%
        transmute(`Experiment ID`, Organism, DOI,
                  Citation = paste0(sub("^(.+)_(NVS\\d+|HS\\d+|ST_J\\d+|FGC\\d+|SLX\\d+)(.*?)$", "\\1 (\\2\\3)", `Experiment ID`)))
    }
    experiments <- experiments %>%
      filter(Citation != "") %>%
      select(`Experiment ID`, Citation, Organism, DOI) %>%
      mutate(Citation = case_when(grepl("https://doi.org", DOI) ~ paste0("<a href=", sub("^.*?\\((.+)\\).*?$", "\\1", DOI), " target=\"_blank\">", Citation, "</a>"),
                                  T ~ Citation)) %>%
      left_join(comparisons, by = "Experiment ID") %>%
      summarise(Treatments = paste(sort(unique(`Treatment`)), collapse = ", "),
                Knockouts = paste(sort(unique(`KO`)), collapse = ", "),
                `Cell lines` = paste(sort(unique(`Cell line`)), collapse = ", "),
                Libraries = paste(sort(unique(Library)), collapse = ", "),
                Source = ds$name,
                .by = c(`Experiment ID`, Citation, Organism)) %>%
      ungroup()
    # Wrangle comparisons metadata
    comparisons <- comparisons %>%
      left_join(experiments %>% select(`Experiment ID`, Citation, Organism), by = "Experiment ID") %>%
      filter(Citation != "") %>%
      # Standardise unspecified levels
      mutate(`Days grown` = as.numeric(case_when(is.na(`Days grown`) | `Days grown` == "" ~ NA, T ~ `Days grown`)),
             `Dose` = case_when(is.na(`Dose`) | `Dose` == "" ~ "Unspecified", T ~ `Dose`),
             `Treatment` = case_when(is.na(`Treatment`) | `Treatment` == "" ~ "Unspecified", T ~ `Treatment`),
             `KO` = case_when(is.na(`KO`) | `KO` == "" ~ "Unspecified", T ~ `KO`),
             `Cell line` = case_when(is.na(`Cell line`) | `Cell line` == "" ~ "Unspecified", T ~ `Cell line`),
             `ControlDays grown` = as.numeric(case_when(is.na(`ControlDays grown`) | `ControlDays grown` == "" ~ NA, T ~ `ControlDays grown`)),
             `ControlDose` = case_when(is.na(`ControlDose`) | `ControlDose` == "" ~ "Unspecified", T ~ `ControlDose`),
             `ControlTreatment` = case_when(is.na(`ControlTreatment`) | `ControlTreatment` == "" ~ "Unspecified", T ~ `ControlTreatment`),
             `ControlKO` = case_when(is.na(`ControlKO`) | `ControlKO` == "" ~ "Unspecified", T ~ `ControlKO`),
             `ControlCell line` = case_when(is.na(`ControlCell line`) | `ControlCell line` == "" ~ "Unspecified", T ~ `ControlCell line`),
             `Library` = case_when(is.na(`Library`) | `Library` == "" ~ "Unspecified", T ~ `Library`)) %>%
      # Standardise inert treatments
      mutate(Treatment = case_when(Treatment %in% c("DMSO", "N/A") ~ "No treatment", T ~ Treatment),
             ControlTreatment = case_when(ControlTreatment %in% c("DMSO", "N/A") ~ "No treatment", T ~ ControlTreatment)) %>%
      # Standardise plasmid
      mutate(`ControlCell line` = case_when(ControlKO == "Plasmid" ~ "Plasmid", T ~ `ControlCell line`),
             ControlKO = case_when(`ControlCell line` == "Plasmid" ~ "Wildtype", T ~ ControlKO)) %>%
      # Calculate timepoint identity
      group_by(`Experiment ID`, KO, `Cell line`, ControlKO, `ControlCell line`) %>%
      mutate(Timepoint = case_when(`Days grown` == max(c(`Days grown`, `ControlDays grown`)) ~ "Final timepoint",
                                   `Days grown` == min(c(`Days grown`, `ControlDays grown`)) ~ "Initial timepoint",
                                   T ~ "Midpoint")) %>%
      ungroup() %>%
      # Factorise
      transmute(Citation = factor(Citation),
                `Experiment ID`,
                Contrast,
                Kind = factor(case_when(grepl("[Ee]ssentialome", Contrast) ~ "Essentialome", # Essentialome by annotated contrast
                                        `ControlCell line` == "Plasmid" & `Cell line` != "Plasmid" ~ "Essentialome", # Essentialome by comparison with plasmid
                                        Treatment != ControlTreatment ~ "Treatment",
                                        KO != ControlKO ~ "Knockout",
                                        `Cell line` != `ControlCell line` ~ "Cell line",
                                        Dose != ControlDose ~ "Dose",
                                        `Days grown` != `ControlDays grown` ~ "Essentialome", # Essentialome by unmatched timepoints
                                        grepl("[Ss]orted", Contrast) ~ "Cell sorting",
                                        grepl("^CRISPR", Contrast) ~ "CRISPR targeting")),
                Endpoint = factor(case_when(grepl(" on gH2A.X", Contrast) ~ sub("^.+ on (.+?)(, .+?| at .+?| in .+?| under .+?| from .+?| to .+?)*$", "\\1", Contrast),
                                            grepl(" on ", Contrast) ~ toupper_first_initial(sub("^.+ on (.+?)(, .+?| at .+?| in .+?| under .+?| from .+?| to .+?)*$", "\\1", Contrast)),
                                            T ~ "Default")),
                Timepoint = factor(Timepoint),
                `Days grown (diff)` = factor(`Days grown`),
                `Dose (diff)` = factor(Dose),
                `Treatment (diff)` = factor(Treatment),
                `Knockout (diff)` = factor(KO),
                `Cell line (diff)` = factor(`Cell line`),
                `Days grown (ref)` = factor(`ControlDays grown`),
                `Dose (ref)` = factor(ControlDose),
                `Treatment (ref)` = factor(ControlTreatment),
                `Knockout (ref)` = factor(ControlKO),
                `Cell line (ref)` = factor(`ControlCell line`),
                Library = factor(Library),
                Organism = factor(Organism),
                `Comparison ID`,
                FriendlyID = paste0(gsub("<.+?>", "", Citation), ": ", Contrast),
                Source = factor(ds$name))
    
    # Obtain genes list
    genes <- sqldf("SELECT [symbol] from gene", connection = con) %>% unlist() %>% setNames(NULL)
    genes <- genes[!grepl("Non-targeting", genes)] # exclude Non-targeting
    genes <- genes[!grepl("Cutting", genes)] # exclude cutting controls
    genes <- genes[!grepl("Olfactory", genes)] # exclude promiscuous olfactory controls
    genes <- genes[!grepl("[ATCG]{19}", genes)] # exclude sequences
    genes <- genes[!grepl("\\|", genes)] # exclude gene arrays
    genes <- genes[!grepl(";", genes)] # exclude gene arrays
    genes <- genes[!grepl(":", genes)] # exclude base edits
    genes <- genes[genes != "X"] # exclude "X"
    genes <- sort(genes)
    return(list(con = con, experiments = experiments, comparisons = comparisons, ontology = ontology, genes = genes))
  }
}

load <- function() {
  # Load every dataset declared in `datasets` (config.R)
  loaded <- lapply(datasets, loadOne)
  names(loaded) <- vapply(datasets, `[[`, character(1), "name")
  loaded <- loaded[!vapply(loaded, is.null, logical(1))]   # drop datasets whose database.db was missing
  
  # Globalise. `cons` is a named list keyed by dataset name
  cons        <<- lapply(loaded, `[[`, "con")
  experiments <<- bind_rows(lapply(loaded, `[[`, "experiments"))
  comparisons <<- bind_rows(lapply(loaded, `[[`, "comparisons"))
  ontology    <<- unique(bind_rows(lapply(loaded, `[[`, "ontology")))
  symbol      <<- unique(unlist(lapply(loaded, `[[`, "symbol")))
  genes       <<- sort(unique(unlist(lapply(loaded, `[[`, "genes"))))
  
  # Decide whether features should be enabled or not
  canGuideLibrary <<- ifelse(exists("file_guideLibrary"), T, F)
  canLibraries <<- ifelse(exists("file_libraries"), T, F)
  canOntology <<- ifelse(nrow(ontology) > 0, T, F)
  if(length(exorcise_root > 0) & length(exorcise_docker) > 0) {
    if(system(paste0("docker images | grep -E '", sub("^(.+):(.+)$", "\\1", exorcise_docker), "\\s*", sub("^(.+):(.+)$", "\\2", exorcise_docker), "'")) == 0) {
      canExorcise <<- T
    } else {
      canExorcise <<- F
    }
  } else {
    canExorcise <<- F
  }
  canExplore <<- ifelse(nrow(experiments) > 0, T, F)
}

Formaldehyde <- NULL # silly way to mask away the `Formaldehyde` built-in dataset so we can actually use "Formaldehyde" in queries