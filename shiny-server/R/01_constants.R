#### Version ####
# See CHANGELOG.md for the full history.
CRAVE_VERSION <- "5.1.2"
CRAVE_UPDATED <- "2026-08-18"

# Retained for backwards compatibility with any deployer customisations that
# referenced these names directly.
ver     <- CRAVE_VERSION
updated <- CRAVE_UPDATED

# Bump when the shape of the cached metadata changes, to invalidate stale caches.
CRAVE_CACHE_VERSION <- 3L

# Suppresses all whitespace around and inside a tag. Use it only where the very
# next thing is punctuation, e.g. tags$strong("Correlate", .noWS = noWS), ", where".
# For a link in the middle of a sentence use extLink(), which leaves the single
# space on each side that the surrounding prose needs.
noWS <- c("before", "after", "outside", "after-begin", "before-end")

#### External documentation links ####
DOI_DRUGZ    <- "https://doi.org/10.1186/s13073-019-0665-3"
DOI_MAGECK   <- "https://doi.org/10.1186/s13059-014-0554-4"
DOI_CHRONOS  <- "https://doi.org/10.1186/s13059-021-02540-7"
DOI_AGNES    <- "https://doi.org/10.1002/9780470316801"
DOI_UMAP     <- "https://doi.org/10.48550/arXiv.1802.03426"
DOI_TSNE     <- "https://www.jmlr.org/papers/v9/vandermaaten08a.html"
DOI_EXORCISE <- "https://doi.org/10.1186/s13073-024-01414-4"
URL_EXORCISE_REPO <- "https://github.com/SimonLammmm/exorcise"

#### Gene-gene similarity methods ####
# Spearman, Pearson and Kendall go straight to cor(). Cosine similarity is computed
# separately; it is a similarity in [-1, 1] like the others, so the divergent colour
# scale and the network's abs(x) >= cutoff filter apply unchanged.
CRAVE_SIMILARITY_METHODS <- c("Spearman", "Pearson", "Kendall", "Cosine similarity")
CRAVE_SIMILARITY_DEFAULT <- "Spearman"
# Kendall's tau is O(n^2) per pair; warn rather than refuse above this many genes.
CRAVE_KENDALL_WARN_GENES <- 300L

#### Clustering methods ####
# Agnes is hierarchical and yields a dendrogram directly. Clara partitions instead,
# so a dendrogram is built over its medoids and each cluster occupies a contiguous
# block of the axis.
CRAVE_DENDRO_METHODS <- c(
  "Agnes: Ward", "Agnes: average", "Agnes: single", "Agnes: complete",
  "Agnes: flexible", "Agnes: weighted", "Clara"
)
CRAVE_DENDRO_DEFAULT <- "Agnes: Ward"

# agnes() accepts euclidean and manhattan; clara() also accepts jaccard.
CRAVE_DISTANCE_AGNES <- c("euclidean", "manhattan")
CRAVE_DISTANCE_CLARA <- c("euclidean", "manhattan", "jaccard")
CRAVE_DISTANCE_DEFAULT <- "euclidean"

# Lance-Williams parameter for agnes(method = "flexible").
CRAVE_AGNES_FLEXIBLE_PAR <- 0.625

CRAVE_CLARA_K_DEFAULT <- 5L

#' Translate a dendrogram method label into the string agnes() expects.
#'
#' @return NULL for Clara, which is not an agnes method.
agnesMethod <- function(label) {
  switch(
    label %||% "",
    "Agnes: Ward"     = "ward",
    "Agnes: average"  = "average",
    "Agnes: single"   = "single",
    "Agnes: complete" = "complete",
    "Agnes: flexible" = "flexible",
    "Agnes: weighted" = "weighted",
    NULL
  )
}

#### Dimensionality reduction methods ####
CRAVE_REDUCE_METHODS <- c("UMAP", "PCA", "t-SNE")
CRAVE_REDUCE_DEFAULT <- "UMAP"

# The other Correlate analyses keep only genes that pass the sidebar cutoff, because
# they draw a gene per axis position and a matrix of hits is the point. Reduce is
# different: cluster structure in an embedding only shows up against the whole
# library, so by default it projects every gene measured in the selected screens and
# ignores the cutoff. Set TRUE to make Reduce respect the cutoff like the others.
CRAVE_REDUCE_USES_CUTOFF <- FALSE
# t-SNE requires 3 * perplexity < n - 1; this is the ceiling before that constraint.
CRAVE_TSNE_MAX_PERPLEXITY <- 30

#### Text filter modes ####
CRAVE_TEXT_FILTER_MODES <- c("Plain text", "Regular expression")

# Separator used to join screen names into a Venn compartment label. A control
# character, so it cannot collide with anything in a contrast label, and it is what
# ggupset::axis_combmatrix() splits on to draw the upset matrix.
CRAVE_COMPARTMENT_SEP <- "\u001f"

#### Analysis method codes ####
# analysis_type_id values used in the `stat` table of a CRAVE dataset. Defined in
# one place so the query builder and every results table decode them identically.
CRAVE_METHODS <- c(MAGeCK = 1L, DrugZ = 2L, Chronos = 3L, Manual = 4L)
CRAVE_METHOD_NAMES <- names(CRAVE_METHODS)

#' Decode analysis_type_id to a method name.
methodName <- function(id) {
  out <- names(CRAVE_METHODS)[match(id, unname(CRAVE_METHODS))]
  out
}

#### Limits ####
# Guard rails for the analyses, collected here rather than left inline.
CRAVE_MAX_SCREENS_PLOT     <- 120   # refuse to draw more screens than this
CRAVE_MAX_SCREENS_DOWNLOAD <- 1000  # refuse a bulk download bigger than this
CRAVE_MAX_AUTO_GOI         <- 250   # cap on auto-detected genes of interest
# Reduce plots one point per gene rather than a square matrix, so the whole library
# is legitimate here and is what reveals cluster structure. Inf means no cap; set a
# number if projection is too slow on your hardware. See also
# CRAVE_REDUCE_USES_CUTOFF below, which is what actually admits non-hit genes.
CRAVE_MAX_AUTO_GOI_REDUCE  <- Inf
CRAVE_MAX_VENN_SETS        <- 7     # more sets than this cannot be drawn as a Venn
CRAVE_MAX_SELECTIZE_SERVER <- 1000  # above this, push choices client-side
CRAVE_ROC_STEPS            <- 2000  # resolution of the ROC sweep
CRAVE_SUMMARY_CHARS        <- 255   # truncate "screens found" text beyond this

#### Correlate defaults ####
defaults_gq_cutoff_choices  <- c(1e-30, 1e-20, 1e-15, 1e-10, 1e-9, 1e-8, 1e-7,
                                 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.25)
defaults_gq_cutoff_selected <- 1e-4
defaults_corr_cutoff_min    <- 0.4
defaults_corr_cutoff_max    <- 1
defaults_corr_cutoff_step   <- 0.05
defaults_corr_cutoff_value  <- 0.7

#### Explore plot customisation defaults ####
EXPLORE_DEFAULTS <- list(
  height          = 900,
  width           = 0,
  width_auto      = TRUE,
  point_colour    = "#ff0087",
  goi_colour      = "#007cff",
  x_axis_log      = "Automatic",
  y_axis_log      = "Automatic",
  y_equals_x      = TRUE,
  density         = FALSE,
  rug_x           = FALSE,
  rug_y           = FALSE
)

#### Correlate plot customisation defaults ####
CORRELATE_DEFAULTS <- list(
  height              = 900,
  width               = 900,
  width_auto          = TRUE,
  point_colour        = "#ff0087",
  genequery_fillscheme = "Treatment",
  heatmap_high        = "#ff0087",
  heatmap_mid         = "#ffffff",
  heatmap_low         = "#00007e"
)

# Flat aliases, kept because config.R files in the wild may reference them.
defaults_explore_customise_height           <- EXPLORE_DEFAULTS$height
defaults_explore_customise_width            <- EXPLORE_DEFAULTS$width
defaults_explore_customise_width_auto       <- EXPLORE_DEFAULTS$width_auto
defaults_explore_customise_point_colour     <- EXPLORE_DEFAULTS$point_colour
defaults_explore_customise_goi_colour       <- EXPLORE_DEFAULTS$goi_colour
defaults_explore_customise_x_axis_log       <- EXPLORE_DEFAULTS$x_axis_log
defaults_explore_customise_y_axis_log       <- EXPLORE_DEFAULTS$y_axis_log
defaults_explore_customise_y_equals_x       <- EXPLORE_DEFAULTS$y_equals_x
defaults_explore_customise_density          <- EXPLORE_DEFAULTS$density
defaults_explore_customise_rug_x            <- EXPLORE_DEFAULTS$rug_x
defaults_explore_customise_rug_y            <- EXPLORE_DEFAULTS$rug_y

defaults_correlate_customise_height              <- CORRELATE_DEFAULTS$height
defaults_correlate_customise_width               <- CORRELATE_DEFAULTS$width
defaults_correlate_customise_width_auto          <- CORRELATE_DEFAULTS$width_auto
defaults_correlate_customise_point_colour        <- CORRELATE_DEFAULTS$point_colour
defaults_correlate_customise_genequery_fillscheme <- CORRELATE_DEFAULTS$genequery_fillscheme
defaults_correlate_customise_heatmap_high        <- CORRELATE_DEFAULTS$heatmap_high
defaults_correlate_customise_heatmap_mid         <- CORRELATE_DEFAULTS$heatmap_mid
defaults_correlate_customise_heatmap_low         <- CORRELATE_DEFAULTS$heatmap_low

#### Gene symbols excluded from the picker ####
# One alternation, matched in a single pass over the symbol table.
CRAVE_GENE_EXCLUDE <- paste(
  "Non-targeting",      # non-targeting controls
  "Cutting",            # cutting controls
  "Olfactory",          # promiscuous olfactory controls
  "[ACGT]{19}",         # bare guide sequences
  "\\|", ";", ":",      # gene arrays and base edits
  sep = "|"
)

#### Hart 2017 Table S2 core essential genes ####
CRAVE_ESSENTIALS <- c(
  "AARS", "ABCE1", "ABCF1", "ACTB", "ACTL6A", "ACTR10", "ACTR2", "ADSL", "ADSS",
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
  "XRCC6", "YARS", "YARS2", "YRDC", "ZBTB8OS", "ZMAT5", "ZNF131", "ZNF259", "ZNF574"
)

#### Example gene set for the Correlate "Example" button ####
CRAVE_EXAMPLE_GENES <- c(
  "ABCC1", "AMBRA1", "CIP2A", "CUL3", "DCLRE1C", "ERCC6L2", "H2AX", "LIG4",
  "MCPH1", "NBN", "NHEJ1", "RAD54L2", "TDP2", "TIAL1", "UBA3", "UBE2M",
  "UBE2K", "ZNF451", "TP53"
)

#### Messages of the day ####
CRAVE_MOTDS <- c(
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

#### Exorcise reference data ####
EXORCISE_GENOMES <- c(
  "GRCh38 Ensembl 111" = "hg38.2020-09-22.2bit",
  "GRCm39 Ensembl 111" = "mm39.2020-07-30.2bit"
)
EXORCISE_EXOMES <- c(
  "GRCh38 Ensembl 111" = "hsa.grch38.refseqall.tsv.gz",
  "GRCm39 Ensembl 111" = "mmu.grcm39.refseqall.tsv.gz"
)
EXORCISE_PRIORITIES <- c(
  "GRCh38 Ensembl 111" = "hsa.priorities.tsv.gz",
  "GRCm39 Ensembl 111" = "mmu.priorities.tsv.gz"
)
# The UI offers "Interference"; both that spelling and "Inhibition" map to "i".
EXORCISE_MODES <- c(
  "Knockout"             = "ko",
  "Interference"         = "i",
  "Inhibition"           = "i",
  "Activation"           = "a",
  "Cytosine base editor" = "cbe",
  "Adenine base editor"  = "abe"
)
EXORCISE_MODE_CHOICES <- c("Knockout", "Interference", "Activation",
                           "Cytosine base editor", "Adenine base editor")
