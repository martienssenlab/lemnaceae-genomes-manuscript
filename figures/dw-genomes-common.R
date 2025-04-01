# Package management ----
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pacman::p_load(devtools)

# Data support
pacman::p_load(
  treeio,             # phylogenetic trees
  rtracklayer,        # bigwigs
  vcfR                # VCFs
)

# Plotting/Graphical packages
pacman::p_load(
  ggtree,             # Treemap plots
  patchwork,          # Compositing for ggplots
  colorspace,         # Lighten and darken colors algorithmically
  viridis,            # Color schemes for heatmaps, etc.
  # "ggridges",       # Ridgeline plots
  # "RColorBrewer",   # Basic color palettes
  glue,               # Interpolated text formatting
  ggrepel,            # Labels for individual datapoints
  ggtext,             # For markdown formatting in ggplot2
  #ggforce,           # Plot zooming
  # "ggnewscale",     # Support multiple color scales in a single ggplot
  ggpubr,             # Publication-ready themes
  lemon,              # Cap axis lines
  ggh4x,              # Independent axes scales for faceting variables, truncated axis lines
  Rgraphviz,          # GraphViz interface (bioconductor)
  ComplexHeatmap      # Heatmaps of all sorts
)
pacman::p_load_gh("zeehio/facetscales@archived")   # Use different scales per facet
pacman::p_install_gh("krassowski/complex-upset")   # ComplexUpset Upset diagrams
pacman::p_load_gh("krassowski/ComplexUpset")

# Essential packages
pacman::p_load(
  optparse,           # Command line args processing
  gtools,             # Sundry useful programming functions
  openxlsx,           # Import/Export MS Excel files
  naturalsort,        # Sort factors naturally
  memoise,            # Cache expensive function calls
  furrr,              # Parallelized purrr map functions
  fuzzyjoin,          # Join tables on partial matches
  tidyverse           # Make life less painful
)


##
## Data loading and formatting ----
##

#' Read chromosome lengths from FASTA index files (.fai) for a list of accessions and compile them into a single data
#' frame.
#'
#' @param faiFiles a named list of fai file paths. Names indicate the accession.
#' @return a data frame with the columns "spec", "chrom", "chrom_length".
readChromLengthsFromFai <- function(faiFiles) {
  colNames <- c("chrom", "chrom_length", "offset", "linebases", "linewidth")
  colTypes <- "cn???"
  #chromLengthsDf <- read_tsv(file = file, col_names = colNames, col_types = colTypes)
  chromLengthsDf <- map_dfr(faiFiles, .id = "accession", read_tsv, col_names = colNames, col_types = colTypes) %>%
    dplyr::select(accession, chrom, chrom_length)
  return(chromLengthsDf)
}

addCumIndex <- function(.data, cols = c("start", "end")) {
  cumData <- .data %>%
    dplyr::left_join(chromNameToNumberAndSubgen(accDwChromLengths)) %>%
    dplyr::select(accession, subgen, chrom, chrom_number, dplyr::everything()) %>%
    dplyr::arrange(accession, subgen, chrom_number) %>%
    dplyr::mutate(dplyr::across(.cols = dplyr::all_of(cols), .names = "{.col}_cum", ~ .x + chrom_cum_index)) %>%
    dplyr::select(-chrom_length, -chrom_cum_index)
  return(cumData)
}

#' Normalize chromosome variables in a data frame by parsing out the chromosome number and subgenome letter if any.
#'
#' @param df a data frame or tibble with chromosome names stored in the variable "chrom"
#' @return the same data frame with the new columns "chrom_number" and "subgen" parsed from the chromosome name
chromNameToNumberAndSubgen <- function(df) {
  df <- df %>% dplyr::mutate(
    subgen = dplyr::case_when(
      stringr::str_detect(chrom, "^chr[0-9]+[a-zA-Z]$") ~
        stringr::str_replace(chrom, "^.*([a-zA-Z])$", "\\1"),
      stringr::str_detect(chrom, "^(HiC_scaffold|scaffold)") ~ "unplaced",
      TRUE ~ as.character(NA)
    ) %>% factor(),
    chrom_number = stringr::str_to_lower(chrom) %>%
      stringr::str_replace(
        pattern = "^[^0-9]+([0-9]+)([a-zA-Z])?$",
        replacement = "\\1"
      ) %>% as.integer()
  )
  return(df)
}

shortAccToLong <- c(
  "sp9509"      = "Spolyrhiza9509",
  "lg7742a"     = "Lgibba7742a",
  "lj7182"      = "Ljaponica7182",
  "lj7182_m"    = "Ljaponica7182_M",
  "lj7182_t"    = "Ljaponica7182_T",
  "lj8627"      = "Ljaponica8627",
  "lj8627_m"    = "Ljaponica8627_M",
  "lj8627_t"    = "Ljaponica8627_T",
  "lj9421"      = "Ljaponica9421",
  "lj9421_m"    = "Ljaponica9421_M",
  "lj9421_t"    = "Ljaponica9421_T",
  "lm7210"      = "Lminor7210",
  "lm9252"      = "Lminor9252",
  "lt9434"      = "Lturionifera9434",
  "wa8730"      = "Waustraliana8730",
  "sp7498"      = "Spolyrhiza",
  "si8410"      = "Sintermedia8410",
  "lu5633"      = "Lminuta5633",
  "lm5500"      = "Lminor5500",
  "Osativa"     = "Osativa",
  "Bdistachyon" = "Bdistachyon"
)

subgenToSpecies <- function(df) {
  df2 <- df %>% dplyr::mutate(
    accession = dplyr::if_else(
      is.na(subgen) | subgen == "unplaced",
      as.character(accession),
      paste0(accession, "_", subgen)
    ) %>% forcats::fct_drop()
  ) %>%
    dplyr::select(-subgen)
  return(df2)
}

#' Load human-readable functional annotations produced by AHRD.
#'
#' @param ahrdAnnots annotation file output of AHRD
#' @param tairTopHits top blast hits from the TAIR database
#' @param gene2Go
#' @return a list of character vectors, where each element is a gene ID mapped to GO terms
loadAnnots <- function(
  ahrdAnnots,
  tairTopHits,
  gene2Go
)
{
  df <- readr::read_tsv(
    ahrdAnnots,
    skip = 3,
    col_names = c("locus", "blast_hit_id", "ahrd_quality", "description", "interpro", "go_term"),
    col_types = "ccfccc",
    col_select = c("locus", "description")
  ) %>%
    # Only keep the primary transcript annotation
    dplyr::filter(stringr::str_detect(locus, "_T001")) %>%
    dplyr::mutate(locus = stringr::str_remove(locus, "_T[0-9]{3,3}$")) %>%
    dplyr::relocate(locus) %>%
    arrange(locus) %>%
    dplyr::left_join(
      by = "locus",
      readr::read_tsv(
        tairTopHits,
        col_names = c("locus", "at_bh_locus"),
        col_select = c("locus", "at_bh_locus")
      ) %>%
        dplyr::filter(stringr::str_detect(locus, "_T001")) %>%
        dplyr::mutate(locus = stringr::str_remove(locus, "_T[0-9]{3,3}$")) %>%
        dplyr::mutate(at_bh_locus = stringr::str_remove(at_bh_locus, "\\.\\d+$"))
    ) %>%
    dplyr::left_join(
      by = "at_bh_locus",
      readr::read_tsv(
        "annotations/tair10_gene_aliases_20210630.clean.txt",
        col_names = c("at_bh_locus", "at_bh_symbol"),
        col_select = c("at_bh_locus", "at_bh_symbol")
      ) %>% dplyr::group_by(at_bh_locus) %>%
        dplyr::summarize(
          at_bh_symbol = paste(at_bh_symbol, collapse="/")
        )
    ) %>%
    dplyr::left_join(tibble::tibble(locus = names(gene2Go), GOs = unname(gene2Go))) %>%
    dplyr::mutate(GOs = purrr::modify_if(GOs, is.null, as.character))

  return(df)
}

#' Load EggNog mapper GO-term associations
#'
#' @param eggnog_annot_file the *.annot file output from EggNog Mapper
#' @param g2g_file simplified locus -> GO ID association file, necessary for TopGO
#' @return a list of character vectors, where each element is a gene ID mapped to GO terms
loadEnGomapToTopGo <- function(
  eggnog_annot_file,
  g2g_file,
  force_overwrite = TRUE
) {
  if (force_overwrite || (! file.exists(g2g_file) || ! file.size(g2g_file) > 0)) {
    readr::read_tsv(
      eggnog_annot_file,
      comment = "##") %>%
      dplyr::rename("query_name" = `#query`) %>%
      dplyr::mutate(
        query_name = stringr::str_remove(query_name, "_T[0-9]{3,3}$"),
        GOs = GOs %>% stringr::str_split(stringr::fixed(",")) %>%
          # IIRC, filter through AnnotationDbi::Term to eliminated deprecated GO Terms
          purrr::map(.f = mAnnDbiTerm) %>%
          purrr:::map(na.omit)
      ) %>%
      dplyr::select(query_name, GOs) %>%
      dplyr::mutate(
        GOs = GOs %>% purrr::map(names) %>%
          purrr::map(.f = stringr::str_flatten, collapse = ",") %>%
          unlist()
      ) %>%
      readr::write_tsv(
        g2g_file,
        col_names = FALSE
      )
  }

  return(topGO::readMappings(g2g_file))
}


#' Load bigwigs for coverage across each chromosome
#'
#' @param bigWigFiles a list of bigWig files
#' @return a tibble with the bigWig coverage info and additional computed coverage averages and normalizations
bwCovFilesToTibble <- function(bigWigFiles) {
  names(bigWigFiles) <- basename(bigWigFiles) %>% stringr::str_replace("([^-_]+).*$", "\\1")
  map_dfr(bigWigFiles, ~ rtracklayer::import.bw(.x) %>% 
            dplyr::as_tibble(), .id = "accession") %>%
    tibble::as_tibble() %>%
    dplyr::rename_with(~ stringr::str_replace_all(.x, pattern = "^value.(.+)$", replacement = "\\1")) %>%
    dplyr::select(-strand) %>%
    dplyr::rename(
      chrom = seqnames,
      coverage = score
    ) %>%
    dplyr::filter(stringr::str_detect(chrom, "^chr[0-9]+")) %>%
    dplyr::mutate(accession = accession %>% shortAccToLong[.] %>% factor()) %>% 
    chromNameToNumberAndSubgen() %>% 
    dplyr::group_by(accession, subgen) %>% 
    dplyr::mutate(
      coverage_subgenome_med = median(coverage),
      coverage_subgenome_mode = statMode(coverage)
    ) %>% 
    dplyr::ungroup() %>%
    dplyr::group_by(accession) %>% 
    dplyr::mutate(coverage_norm = coverage/min(coverage_subgenome_mode)) %>% 
    dplyr::ungroup() %>% 
    addCumIndex()
}



##
## Utility functions ----
##
uniqueCaseInsensitive <- function(x) {
  x[!duplicated(tolower(x))]
}


##
## Math ----
##

isOutlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

# https://stackoverflow.com/a/8189441/329881
# Caution: only returns highest mode of multimodal data
statMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

##
## Pretty printing ----
##
dwHighlighter <- function(orig_labels, pal = pal4) {
  labels <- dplyr::left_join(tibble(accession = orig_labels), accTblFullExtra) %>%
    dplyr::mutate(sci_name_abbrev = dplyr::if_else(is.na(sci_name_abbrev), accession, sci_name_abbrev)) %>%
    dplyr::pull(sci_name_abbrev)
  labels <- stringr::str_replace(labels, "([^ ]+) ([^ ]+)(.*)", "*\\1 \\2*\\3")
  labels <- dplyr::if_else(orig_labels %in% names(pal),
                  glue::glue("<strong style='color:{pal[orig_labels]};'>{labels}</strong>"),
                  glue::glue("{labels}"))
  return(labels)
}

dwHighlighterMultiline <- function(orig_labels, pal = pal4) {
  labels <- dplyr::left_join(tibble(accession = orig_labels), accTblFullExtra) %>%
    dplyr::mutate(sci_name_abbrev = dplyr::if_else(is.na(sci_name_abbrev), accession, sci_name_abbrev)) %>%
    dplyr::pull(sci_name_abbrev)
  labels <- stringr::str_replace(labels, "([^ ]+) ([^ ]+)(.*)", "*\\1 \\2*<br>\\3")
  labels <- dplyr::if_else(orig_labels %in% names(pal),
                           glue::glue("<strong style='color:{pal[orig_labels]};'>{labels}</strong>"),
                           glue::glue("{labels}"))
  return(labels)
}

dwHighlighterMsa <- function(orig_labels, pal = pal4, additionalBolded = NULL) {
  # Set aside anything following a space to rejoin later
  labels <- tibble::tibble(accession = orig_labels) %>%
    tidyr::separate(col = accession, into = c("accession", "suffix"), sep = " ") %>%
    dplyr::left_join(accTblFullExtra) %>%
    dplyr::mutate(sci_name_abbrev = dplyr::if_else(is.na(sci_name_abbrev), accession, sci_name_abbrev)) %>%
    dplyr::mutate(sci_name_abbrev = sci_name_abbrev %>% stringr::str_replace("([^ ]+) ([^ ]+)(.*)", "*\\1 \\2*\\3")) %>%
    dplyr::mutate(
      sci_name_abbrev = dplyr::case_when(
        accession %in% names(pal) ~ glue::glue("<strong style='color:{pal[accession]};'>{sci_name_abbrev}</strong>"),
        accession %in% additionalBolded ~ glue::glue("<strong>{sci_name_abbrev}</strong>"),
        TRUE ~ glue::glue("{sci_name_abbrev}")
      )
    ) %>%
    dplyr::select(sci_name_abbrev, suffix) %>%
    tidyr::unite(col = "new_label", sep = " ", sci_name_abbrev, suffix, na.rm = TRUE) %>%
    dplyr::select(new_label)
  return(labels)
}

dwHighlighter2Subgen <- function(orig_labels, palette = palsDfV4) {
  labels <- left_join(tibble(accession = orig_labels), accTblSubgen)$sci_name_abbrev
  #labels <- str_replace(labels, "([^ ]+) ([^ ]+) ?([0-9a-zA-Z]*) ?([A-Z]*)", "~bolditalic('\\1 \\2')~bold('\\3 \\4')")
  labels <- ifelse(
    orig_labels %in% palette$acc,
    str_replace(labels, "([^ ]+) ([^ ]+) ?([0-9a-zA-Z]*) ?([A-Z]*)", "~bolditalic('\\1 \\2')~bold('\\3 \\4')"),
    str_replace(labels, "([^ ]+) ([^ ]+) ?([0-9a-zA-Z]*) ?([A-Z]*)", "~italic('\\1 \\2')~plain('\\3 \\4')")
  )
  return(labels)
}

sciNameLabeller <- function(labels, multi_line = TRUE) {
  labels <- dplyr::left_join(tibble(accession = labels), accTbl %>% dplyr::bind_rows(accTblFullExtra) %>% unique())$sci_name_abbrev
  labels <- stringr::str_replace(labels, "([^ ]+) ([^ ]+)(.*)", "*\\1 \\2*\\3")
  return(label_value(labels, multi_line))
}

sciNameLabeller2 <- function(labels, multi_line = TRUE) {
  labels <- dplyr::left_join(tibble(accession = labels), accTbl %>% dplyr::bind_rows(accTblFullExtra) %>% unique())$sci_name_abbrev
  labels <- stringr::str_replace(labels, "([^ ]+) ([^ ]+)( [^ ]+)?( [^ ]+)?", "*\\1 \\2*\\4")
  return(label_value(labels, multi_line))
}

sciNameLabeller3 <- function(labels) {
  labels <- dplyr::left_join(tibble(accession = labels), accTbl %>% dplyr::bind_rows(accTblFullExtra) %>% unique())$sci_name_abbrev
  labels <- stringr::str_replace(labels, "([^ ]+) ([^ ]+)( [^ ]+)?( [^ ]+)?", "bquote(italic(.(\\1 \\2)))\\4")
  return(labels)
}

sciNamePlotMathLabeller <- function(labels, multi_line = TRUE) {
  labels <- dplyr::left_join(tibble(accession = labels), accTbl %>% dplyr::bind_rows(accTblFullExtra) %>% unique())$sci_name_abbrev
  labels <- stringr::str_replace(labels, "([^ ]+) ([^ ]+)(.*)", "italic('\\1 \\2')\\3")
  labels <- stringr::str_replace_all(labels, " ", "~")
  return(label_parsed(as.expression(labels), multi_line))
}



##
## Species metadata ----
##
accDwAll <- c(
  "Spolyrhiza",
  "Spolyrhiza9509",
  "Sintermedia8410",
  "Waustraliana8730",
  "Lgibba7742a",
  "Lminor7210",
  "Lminor9252",
  "Ljaponica7182_M",
  "Ljaponica8627_M",
  "Ljaponica9421_M",
  "Ljaponica7182",
  "Ljaponica8627",
  "Ljaponica9421",
  "Ljaponica8627_T",
  "Ljaponica7182_T",
  "Ljaponica9421_T",
  "Lturionifera9434"
)
accDwSubgen <- accDwAll %>% purrr::discard(~ .x %in% c("Ljaponica7182", "Ljaponica8627", "Ljaponica9421"))
accDw <- accDwAll %>% purrr::discard(~ stringr::str_detect(.x, pattern = "_T$|_M$"))

accDwFais <- c(
  "../../assemblies/sp9509/asm3v1/sp9509.asm3v1.fasta.fai",
  "../../assemblies/lg7742a/asm202106v1/lg7742a.asm202106v1.fasta.fai",
  "../../assemblies/lj7182/asm202106v1/lj7182.asm202106v1.fasta.fai",
  "../../assemblies/lj8627/asm202106v1/lj8627.asm202106v1.fasta.fai",
  "../../assemblies/lj9421/asm202106v1/lj9421.asm202106v1.fasta.fai",
  "../../assemblies/lm7210/asm202106v1/lm7210.asm202106v1.fasta.fai",
  "../../assemblies/lm9252/asm202106v1/lm9252.asm202106v1.fasta.fai",
  "../../assemblies/lt9434/asm202106v1/lt9434.asm202106v1.fasta.fai",
  "../../assemblies/wa8730/asm202106v1/wa8730.asm202106v1.fasta.fai"
)
names(accDwFais) <- accDwFais %>% basename() %>% stringr::str_replace(pattern = "^([^.]+).*", replacement = "\\1") %>% shortAccToLong[.]
accDwChromLengths <- readChromLengthsFromFai(accDwFais) %>%
  chromNameToNumberAndSubgen() %>%
  dplyr::group_by(accession, subgen) %>%
  dplyr::mutate(chrom_cum_index = cumsum(chrom_length) - chrom_length) %>%
  dplyr::ungroup()
accDwSubgenLengths <- accDwChromLengths %>%
  dplyr::filter(chrom %>% stringr::str_detect(pattern = '^chr[0-9]+')) %>%
  dplyr::group_by(accession, subgen) %>%
  dplyr::summarize(length = sum(chrom_length)) %>%
  tidyr::unite(col = accSubgen, sep = "_", remove = FALSE, accession, subgen) %>%
  dplyr::mutate(accession = accSubgen %>% stringr::str_remove(pattern = "_NA") %>% factor()) %>%
  dplyr::select(accession, length)

accAll <- c(
  accDwAll,
  "Zmarina",
  "Sbicolor",
  "Zmays",
  "Sitalica",
  "Osativa",
  "Bdistachyon",
  "Acomosus",
  "Eguineensis",
  "Macuminata",
  # "Aofficinalis",
  "Athaliana",
  "Vvinifera",
  "Cdemersum",
  "Ncolorata",
  "Eferox",
  "Atrichopoda",
  "Gmontanum"
)
acc <- accAll %>% purrr::discard(~ stringr::str_detect(.x, pattern = "_T$|_M$"))
accSubgen <- accAll %>% purrr::discard(~ .x %in% c("Ljaponica7182", "Ljaponica8627", "Ljaponica9421"))
accAngio <- acc %>% purrr::discard(~ .x == "Gmontanum")
accAngioSubgen <- accSubgen %>% purrr::discard(~ .x == "Gmontanum")

### Phylogenetic groupings ----
###
### Accessions and accession intersections of interest to determine exclusively common and missing
### hierarchical orthogroups.
###
accSets <- list(
  "Spirodelas" = c("Sintermedia8410", "Spolyrhiza", "Spolyrhiza9509"),
  "Lemnas" = c("Lgibba7742a", "Ljaponica7182", "Ljaponica8627", "Ljaponica9421", "Lminor7210", "Lminor9252",
               "Lturionifera9434"),
  "Wolffias" = c("Waustraliana8730"),
  "Spolyrhizas" = c("Spolyrhiza", "Spolyrhiza9509"),
  "Sintermedias" = c("Sintermedia8410"),
  "Lgibbas" = c("Lgibba7742a"),
  "Ljaponicas" = c("Ljaponica7182", "Ljaponica8627", "Ljaponica9421"),
  "Lminors" = c("Lminor7210", "Lminor9252"),
  "Lturioniferas" = c("Lturionifera9434"),
  "Waustralianas" = c("Waustraliana8730"),
  "Lgibba_and_Waustraliana" = c("Lgibba7742a", "Waustraliana8730"),
  "Spoly_and_Lturionifera" = c("Spolyrhiza", "Spolyrhiza9509", "Lturionifera9434"),
  "sinking_dw" = c("Spolyrhiza", "Spolyrhiza9509", "Lturionifera9434", "Waustraliana8730"),
  "Lemnaceae" = accDw,
  "dicots" = c("Athaliana", "Vvinifera"),
  "monocots" = c(
    "Acomosus", "Bdistachyon", "Eguineensis", "Macuminata", "Osativa", "Sbicolor", "Sitalica", "Zmays", "Zmarina", accDw
  ),
  "Ceratophyllales" = c("Cdemersum"),
  "Nymphaeales" = c("Eferox", "Ncolorata"),
  "Amborellales" = c("Atrichopoda"),
  "gymnosperms" = c("Gmontanum"),
  "aquatic_emergent" = c("Osativa"),
  "aquatic_submerged" = c("Cdemersum", "Zmarina"),
  "aquatic_submerged_dw" = c("Cdemersum", "Zmarina", accDw),
  "cdemersum_dw" = c("Cdemersum", accDw),
  "zmarina_dw" = c("Zmarina", accDw),
  "aquatic_floating_dw" = c("Eferox", "Ncolorata", accDw),
  "aquatic_sinking" = c("Spolyrhiza", "Spolyrhiza9509", "Lturionifera9434", "Waustraliana8730", "Cdemersum", "Zmarina"),
  "freshwater_sinking" = c("Cdemersum", "Spolyrhiza", "Spolyrhiza9509", "Lturionifera9434", "Waustraliana8730"),
  "aquatic_rootless" = c("Cdemersum", "Waustraliana8730"),
  "aquatic_ESF_dw" = c("Osativa", "Cdemersum", "Zmarina", "Eferox", "Ncolorata", accDw),
  "aquatic_SF_dw" = c("Cdemersum", "Zmarina", "Eferox", "Ncolorata", accDw),
  "allAcc" = acc,
  "angiosperms" = accAngio,
  "basal_angiosperms" = c("Atrichopoda", "Eferox", "Ncolorata")
  # "terrestrial_ESF" = acc %>% discard(~ .x %in% c("Osativa", "Cdemersum", "Zmarina", "Eferox", "Ncolorata", accDw)),
  # "terrestrial_SF" = acc %>% discard(~ .x %in% c("Cdemersum", "Zmarina", "Eferox", "Ncolorata", accDw)),
  # "terrestrial_F" = acc %>% discard(~ .x %in% c("Eferox", "Ncolorata", accDw))
)
accSets <- c(accSets, list("core_angiosperms" = setdiff(accSets$angiosperms, accSets$basal_angiosperms)))
accSetsAngio <- accSets %>% purrr::map(~ .x %>% purrr::discard(~ .x %in% accSets$gymnosperms)) %>% purrr::compact()
accDwSubgenList <- as.list(accDwSubgen)
names(accDwSubgenList) <- accDwSubgen
accSetsDwSubgen <- c(
  accDwSubgenList,
  accSets[c("Spirodelas", "Spolyrhizas", "Lminors", "Lturioniferas", "Spoly_and_Lturionifera", "sinking_dw")],
  list(
    "Ljaponicas_subM" = c("Ljaponica7182_M", "Ljaponica8627_M", "Ljaponica9421_M"),
    "Ljaponicas_subT" = c("Ljaponica7182_T", "Ljaponica8627_T", "Ljaponica9421_T"),
    "Lminors_with_subs" = c("Lminor7210", "Lminor9252", "Ljaponica7182_M", "Ljaponica8627_M", "Ljaponica9421_M"),
    "Lturioniferas_with_subs" = c("Lturionifera9434", "Ljaponica7182_T", "Ljaponica8627_T", "Ljaponica9421_T"),
    "Lemnas" = c(accSets$Lgibbas, accSets$Lminors, "Ljaponica7182_M", "Ljaponica8627_M", "Ljaponica9421_M", 
                 "Ljaponica7182_T", "Ljaponica8627_T", "Ljaponica9421_T", accSets$Lturioniferas),
    "Lemnaceae" = accDwSubgen
  )
)

accTbl <- read_tsv(
  "../../accessions.txt",
  col_names = c("accession", "sci_name_abbrev", "sci_name_full", "common_name", "source"),
  col_types = "fcccc"
) %>%
  dplyr::select(-source) %>%
  dplyr::mutate(label = paste0("italic('", sci_name_abbrev, "')")) %>%
  dplyr::mutate(accession = forcats::fct_relevel(accession, acc)) %>%
  dplyr::mutate(genus = sci_name_full %>% stringr::str_extract(pattern = "^([A-Za-z]+)"))

accTblSubgen <- read_tsv(
  "../../accessions.subgen.txt",
  col_names = c("accession", "sci_name_abbrev", "sci_name_full", "common_name", "source"),
  col_types = "fcccc"
) %>%
  dplyr::select(-source) %>%
  dplyr::mutate(label = paste0("italic('", sci_name_abbrev, "')")) %>%
  dplyr::filter(accession %in% accSubgen) %>% 
  dplyr::mutate(accession = forcats::fct_relevel(accession, accSubgen) %>% droplevels()) %>%
  dplyr::mutate(genus = sci_name_full %>% stringr::str_extract(pattern = "^([A-Za-z]+)"))
accTblFull <- accTbl %>% dplyr::bind_rows(accTblSubgen) %>% unique()
accTblFullExtra <- accTblFull %>% dplyr::bind_rows(
  tibble::tribble(
    ~accession, ~sci_name_abbrev, ~sci_name_full, ~common_name, ~label, ~genus,
    "Lminor5500", "L. minor 5500", "Lemna minor 5500", "common duckweed", "italic('L. minor 5500')", "Lemna",
    "Lminuta5633", "L. minuta 5633", "Lemna minuta 5633", "least duckweed", "italic('L. minor 5633')", "Lemna"
  )
)

# Color Palettes ----

#pal <- c("#152451","#945784","#e87744","#5398be")
#names(pal) <- c("Spolyrhiza9509","Lgibba7742a","Ljaponica8627","Waustraliana8730")
# coolers.co palette "dw genomes 4 bright 30"
#palDwGen <- c("#b096bc", "#f2b868", "#dd977d", "#e14840", "#4d40a4")

# 4-color all-purpose palette
palAux <- c("#083d77","#f4d35e","#ee964b","#f95738")

palsDfV4 <- tibble::tribble(
  ## Yellows:
  # Mikado Yellow:     #FFC20A
  # Orange Yellow:     #f3b700
  ## Purples:
  # Fuschsia crayola:  #BD4CB2
  # Byzantine:         #B342A7
  # English lavendar:  #a6808c
  # Plum:              #95378C
  #4c5b5c, #2ab7ca, #dbad6a, #f3b700, #ff8811, #f63e02, #a6808c
  ~acc,                ~base_color,
  "all",               "#4c5b5c",
  "Bdistachyon",       "#4c5b5c",
  "Osativa",           "#4c5b5c",
  "Spolyrhiza9509",    "#2ab7ca",
  "Lgibba7742a",       "#dbad6a",
  "Lminor7210",        "#ffc71f",
  "Lminor9252",        "#ffc71f",
  "Ljaponica7182",     "#ff8811",
  "Ljaponica7182_M",   "#ff8811",
  "Ljaponica7182_T",   "#ff8811",
  "Ljaponica8627",     "#ff8811",
  "Ljaponica8627_M",   "#ff8811",
  "Ljaponica8627_T",   "#ff8811",
  "Ljaponica9421",     "#ff8811",
  "Ljaponica9421_M",   "#ff8811",
  "Ljaponica9421_T",   "#ff8811",
  "Lturionifera9434",  "#f63e02",
  "Waustraliana8730",  "#95378C"
) %>%
  dplyr::mutate(
    light30 = base_color %>% colorspace::lighten(0.30),
    dark30 = base_color %>% colorspace::darken(0.30))
pal4 <- palsDfV4 %>% dplyr::pull(base_color, name = acc)
pal4Light <- palsDfV4 %>% dplyr::pull(light30, name = acc)
pal4Dark <- palsDfV4 %>% dplyr::pull(dark30, name = acc)
pal4DarkerLabels <- palsDfV4 %>% 
  dplyr::mutate(
    base_color = dplyr::case_when(
      acc %>% stringr::str_detect(pattern = "Lminor") ~ "#f3b700",
      .default = base_color
    )
  ) %>% 
  dplyr::pull(base_color, name = acc)

palsDfV3 <- tibble::tribble(
  #4c5b5c, #3891a6, #db4c40, #f3b700, #ff8811, #f63e02, #79b473
  ~acc,                ~base_color,
  "all",               "#4c5b5c",
  "Spolyrhiza9509",    "#3891a6",
  "Lgibba7742a",       "#cc3f0c",
  "Lminor7210",        "#f3b700",
  "Lminor9252",        "#f3b700",
  "Ljaponica7182",     "#ff8811",
  "Ljaponica7182_M",   "#ff8811",
  "Ljaponica7182_T",   "#ff8811",
  "Ljaponica8627",     "#ff8811",
  "Ljaponica8627_M",   "#ff8811",
  "Ljaponica8627_T",   "#ff8811",
  "Ljaponica9421",     "#ff8811",
  "Ljaponica9421_M",   "#ff8811",
  "Ljaponica9421_T",   "#ff8811",
  "Lturionifera9434",  "#f63e02",
  "Waustraliana8730",  "#A63D40"
) %>%
  dplyr::mutate(
    light30 = base_color %>% colorspace::lighten(0.30),
    dark30 = base_color %>% colorspace::darken(0.30))
pal3 <- palsDfV3 %>% dplyr::pull(base_color, name = acc)
pal3Light <- palsDfV3 %>% dplyr::pull(light30, name = acc)
pal3Dark <- palsDfV3 %>% dplyr::pull(dark30, name = acc)

palsDfV2 <- tibble::tribble(
  ~acc,                ~base_color,
  "all",               "#4c5b5c",
  "Spolyrhiza9509",    "#3891a6",
  "Lgibba7742a",       "#e57c04",
  "Lminor7210",        "#f3b700",
  "Lminor9252",        "#f3b700",
  "Ljaponica7182",     "#ff6201",
  "Ljaponica7182_M",   "#ff6201",
  "Ljaponica7182_T",   "#ff6201",
  "Ljaponica8627",     "#ff6201",
  "Ljaponica8627_M",   "#ff6201",
  "Ljaponica8627_T",   "#ff6201",
  "Ljaponica9421",     "#ff6201",
  "Ljaponica9421_M",   "#ff6201",
  "Ljaponica9421_T",   "#ff6201",
  "Lturionifera9434",  "#f63e02",
  "Waustraliana8730",  "#6b2d5c"
) %>%
  dplyr::mutate(
    light30 = base_color %>% colorspace::lighten(0.30),
    dark30 = base_color %>% colorspace::darken(0.30))
pal <- palsDfV2 %>% dplyr::pull(base_color, name = acc)
palLight <- palsDfV2 %>% dplyr::pull(light30, name = acc)
palDark <- palsDfV2 %>% dplyr::pull(dark30, name = acc)

# pal <- c(
#   "Spolyrhiza9509"    = "#868686",
#   "Lgibba7742a"       = "#945784",
#   "Ljaponica7182"     = "#e87744",
#   "Ljaponica8627"     = "#e87744",
#   "Ljaponica9241"     = "#e87744",
#   "Lminor7210"        = "#f0ad4f",
#   "Lminor9252"        = "#f0ad4f",
#   "Lturionifera9434"  = "#e04138",
#   "Waustraliana8730"  = "#5398be"
# )
palSubgen <- c(
  pal[! names(pal) %in% c("lj7182", "Ljaponica8627", "lj9421")],
  "Ljaponica7182_M"  = "#e87744",
  "Ljaponica8627_M"  = "#e87744",
  "Ljaponica9421_M"  = "#e87744",
  "Ljaponica7182_T"  = "#e87744",
  "Ljaponica8627_T"  = "#e87744",
  "Ljaponica9421_T"  = "#e87744"
)
# palDark <- c(
#   "Spolyrhiza9509"  = "#444444",
#   "Lgibba7742a" = "#4a2b42",
#   "Lminor7210"  = "#945c0d",
#   "Lminor9252"  = "#945c0d",
#   "Ljaponica7182"  = "#873a10",
#   "Ljaponica8627"  = "#873a10",
#   "Ljaponica9241"  = "#873a10",
#   "Lturionifera9434"  = "#791813",
#   "Waustraliana8730"  = "#264e64"
# )
# palLight <- c(
#   "Spolyrhiza9509"  = "#939393",
#   "Lgibba7742a" = "#a36392",
#   "Lminor7210"  = "#f1b663",
#   "Lminor9252"  = "#f1b663",
#   "Ljaponica8627"  = "#ea8658",
#   "Ljaponica7182"  = "#ea8658",
#   "Ljaponica9241"  = "#ea8658",
#   "Lturionifera9434"  = "#e3544c",
#   "Waustraliana8730"  = "#66a4c5"
# )

project_theme <- theme_pubr() +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    text = element_text(
      family = "Helvetica",
      face = "plain",
      colour = "#2e2e2e",
      size = 6,
      lineheight = 1.25
    ),
    axis.line = element_line(linewidth = 0.25, color = "#2e2e2e"),
    axis.text.x = element_markdown(
      colour = "#2e2e2e",
      size = 6
    ),
    axis.text.y = element_markdown(
      colour = "#2e2e2e",
      size = 6
    ),
    axis.ticks = element_line(linewidth = 0.25, color = "#2e2e2e"),
    axis.title.x = element_textbox(
      color = "#2e2e2e",
      size = 6,
      padding = margin(4, 4, 4, 4),
      margin = margin(4, 0, 0, 0)
    ),
    axis.title.y = element_textbox(
      color = "#2e2e2e",
      size = 6,
      padding = margin(4, 4, 4, 4),
      margin = margin(0, 0, 4, 0),
      orientation = "left-rotated"
    ),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_markdown(
      colour = "#2e2e2e",
      size = 5,
      margin = margin(0, 4, 0, 1),
      padding = margin(0, 0, 0, 0)
    ),
    # This sometimes causes an error
    # legend.text = element_markdown(
    #   colour = "#2e2e2e",
    #   size = 6,
    #   margin = margin(0, 6, 0, 0),
    #   padding = margin(0, 0, 0, 0)
    # ),
    legend.key.size = unit(0.75, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.justification = "left",
    legend.box.margin = margin(-10, 0, 0, -10),
    legend.spacing.x = unit(0, "lines"),
    strip.text = element_markdown(
      colour = "#2e2e2e",
      size = 6,
      margin = margin(t = 4, r = 4, b = 4, l = 4)
    ),
    strip.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(
      color = "#2e2e2e",
      size = 6,
      face = "bold",
      margin = margin(10, 0, 10, 0)
    ),
    panel.border = element_blank()
  )
