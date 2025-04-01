# topGO and UpSetR have problematic conflicts (e.g. with tidyverse)
pacman::p_load(
  topGO,                 # Gene Ontology enrichment analysis
  org.At.tair.db,        # Needed for rrvgo
  rrvgo                  # reduce and vizualize GOs
)

# For importing RNA-seq expression estimates from salmon
pacman::p_load(tximport)

# Updated Gene Ontology sqlite DB. Also loads AnnotationDBI, which has a method called "select", 
# so always a good idea to prefix external function calls by package name.
pacman::p_load(GO.db)

# Common theme elements, globals and functions
source("../dw-genomes-common.R")
theme_set(project_theme)


# 
# Data loading ----
# 

ofResultsDir <- "of-gmont-outgroup-ljsubgens"
dir.create("plots", showWarnings = FALSE, recursive = TRUE)

## 
## Reference annotations ----
## 

# TAIR10 primary gene symbols
athalGeneTbl <- readr::read_tsv(
  lazy = TRUE,
  "reference/athaliana/gene_aliases_20210630.clean.txt",
  col_types = "ccc"
) %>% 
  dplyr::group_by(locus_name) %>% 
  dplyr::summarize(
    symbol = paste(symbol, collapse="/"), 
    full_name = full_name %>% na.omit() %>% unique() %>% paste(collapse=";")
  )

# RAP-DB (including oryzabase and CGSNL) gene symbols. Have to map across RAP-DB locus IDs to MSU 
# locus IDs from the phytozome v13 derived proteome. 
osatiGeneTbl <- readr::read_tsv(
  lazy = TRUE,
  "reference/oryza_sativa/rapdb/IRGSP-1.0_representative_annotation_2020-12-02.tsv",
  col_types = "cccccccccccccccci"
) %>% 
  dplyr::left_join(
    readr::read_tsv(
      lazy = TRUE,
      "reference/oryza_sativa/rapdb/RAP-MSU_2020-12-02.txt", 
      col_names = c("rap_id","msu_ids")
      ) %>% 
      dplyr::mutate(locus_name = msu_ids %>% stringr::str_extract("^LOC_[^.]*")) %>% 
      dplyr::select(rap_id, locus_name),
    by = c("Locus_ID" = "rap_id")
  ) %>% 
  dplyr::filter(! is.na(locus_name)) %>% 
  dplyr::mutate(
    symbol = dplyr::case_when(
        ! is.na(`CGSNL Gene Symbol`) ~ `CGSNL Gene Symbol`,
        ! is.na(`RAP-DB Gene Symbol Synonym(s)`) ~ `RAP-DB Gene Symbol Synonym(s)` %>% 
          stringr::str_extract("^[^,]*"),
        ! is.na(`Oryzabase Gene Symbol Synonym(s)`) ~ `Oryzabase Gene Symbol Synonym(s)` %>% 
          stringr::str_extract("^[^,]*")
      ),
    full_name = dplyr::case_when(
      ! is.na(`CGSNL Gene Name`) ~ `CGSNL Gene Name`,
      ! is.na(`RAP-DB Gene Name Synonym(s)`) ~ `RAP-DB Gene Name Synonym(s)` %>% 
        stringr::str_extract("^[^,]*"),
      ! is.na(`Oryzabase Gene Name Synonym(s)`) ~ `Oryzabase Gene Name Synonym(s)` %>% 
        stringr::str_extract("^[^,]*")
    )
  ) %>% 
  dplyr::rename(description = Description) %>% 
  dplyr::select(locus_name, symbol, full_name, description)

zmaysGeneTbl <- readr::read_tsv(
  lazy = TRUE,
  "reference/zea_mays/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.full_gene_data.txt.gz"
) %>% 
  dplyr::select(gene_model, locus_symbol, locus_name) %>% 
  dplyr::rename(full_name = locus_name) %>% 
  dplyr::rename(locus_name = gene_model) %>% 
  dplyr::group_by(locus_name) %>% 
  dplyr::summarize(
    symbol = paste(locus_symbol %>% na.omit() %>% stringi::stri_trans_toupper(), collapse="/"), 
    full_name = full_name %>% na.omit() %>% unique() %>% paste(collapse=";")
  ) %>% dplyr::filter(! symbol %>% stringi::stri_isempty())

## Functional annotations ----
## Annotations of the included proteomes by eggnog-mapper on the v5 database.

# Greatly speeds up GO term retrieval
mAnnDbiTerm <- memoise::memoise(AnnotationDbi::Term)

eggMembers <- readr::read_tsv(
  lazy = TRUE,
  "eggnog-mapper/33090_members.tsv.gz",
  col_names = c(
    "TaxonomicLevel",
    "GroupName",
    "ProteinCount",
    "SpeciesCount",
    "ProteinIDs",
    "COGFunctionalCategory"
  ),
  col_types = "ccddcc"
) %>% tidyr::unite("GroupName", c("GroupName", "TaxonomicLevel"), sep = "@", remove = TRUE)

#' Get all of the A. thaliana gene symbols for the given EggNog orthogroup ID
#'
#' @param ogIDs String containing one or more EggNog v5 orthogroup IDs including taxonomy ID, e.g. 
#' "37Q3J@33090,37P3K@33090"
#' @return String containing gene symbols, e.g. "IOS1,SIFT3,SIF4"
athSymbolsInEggNogOg <- function(ogIds) {
  ogIdsVec <- ogIds %>% stringr::str_split(",") %>% unlist()
  protIds <- eggMembers %>% 
    dplyr::filter(GroupName %in% ogIdsVec) %>% 
    dplyr::pull(ProteinIDs) %>% 
    stringr::str_split(",") %>% 
    unlist() %>% 
    stringr::str_subset(pattern = "^3702\\.") %>% 
    stringr::str_remove(pattern = "^3702\\.") %>% 
    stringr::str_remove(pattern = "\\.[0-9]+$")
  symbols <- athalGeneTbl %>% 
    dplyr::filter(locus_name %in% protIds) %>% 
    dplyr::pull(symbol) %>% 
    stringr::str_flatten(collapse = ",")
  return(symbols)
}


#' Load eggNOG mapper and AHRD annotations for each query protein in the OrthoFinder analysis.
#'
#' @return a dataframe of the annotations.
loadAnnots <- function(
)
{
  eggFiles <- paste0(ofResultsDir, "/eggnog/", accTblSubgen$accession, ".emapper.annotations")
  names(eggFiles) <- accTblSubgen$accession
  
  ahrdFiles <- paste0(ofResultsDir, "/ahrd/", accTblSubgen$accession, ".ahrd_output.tsv.gz")
  names(ahrdFiles) <- accTblSubgen$accession
  
  # Read in each annotation file and compile into a single dataframe
  # eggCols <- c(
  # "query", "seed_ortholog", "evalue", "score", "eN_OGs", "max_annot_lvl", "COG_category", "Description", 
  # "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", 
  # "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
  # )
  eggColsSelect <- c(
    "#query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "Description", "Preferred_name", "GOs", "PFAMs"
  )
  eggColsTypes <- "ccddccccc"
  eggs <- purrr::map_dfr(
    eggFiles, .id = "accession", 
    readr::read_tsv, 
    lazy = TRUE,
    comment = "##",
    col_select = dplyr::all_of(eggColsSelect),
    col_types = eggColsTypes
  ) %>%
    dplyr::mutate(
      accession = as.factor(accession),
      # EggNOG v5 includes obsolete GO terms in its annotations. These are purged by querying the 
      # annotated terms against the GO.db.
      GOs = GOs %>% stringr::str_split(fixed(",")) %>% purrr::map(mAnnDbiTerm) %>% 
        purrr:::map(na.omit)
    ) %>% 
    dplyr::rename(
      "query_name" = "#query",
      "eN_seed_olog_score" = "score",
      "eN_OGs" = "eggNOG_OGs",
      "eN_seed_olog_evalue" = "evalue",
      "eN_desc" = "Description",
      "eN_seed_olog" = "seed_ortholog"
    )
  print(paste0("Loaded ", nrow(eggs), " eggNOG annotations."))
  
  # Read in each annotation file and compile into a single dataframe
  ahrds <- purrr::map_dfr(
    ahrdFiles, .id = "accession", 
    readr::read_tsv, 
    lazy = TRUE,
    skip = 3,
    col_names = c("query_name", "blast_hit_id", "ahrd_qual", "ahrd_desc", "interpro", "go_term"),
    col_types = "ccfccc",
    col_select = c("query_name", "ahrd_qual", "ahrd_desc")
  ) %>%
    dplyr::mutate(ahrd_desc = ahrd_desc %>% stringr::str_replace("Unknown protein(,)?", as.character(NA))) %>% 
    dplyr::relocate(query_name) %>% 
    dplyr::arrange(query_name)
  
  print(paste0("Loaded ", nrow(ahrds), " AHRD annotations."))
  
  df <- ahrds %>% left_join(eggs, by = c("query_name", "accession"))
  
  return(df)
}

annots <- loadAnnots()

# topGO needs this simplified locus -> GO ID association file
annots %>%
  dplyr::filter(! is.na(eN_seed_olog)) %>% 
  dplyr::select(query_name, GOs) %>%
  dplyr::mutate(
    GOs = GOs %>% purrr::map(names) %>% 
      purrr::map(.f = stringr::str_flatten, collapse = ",") %>% 
      unlist()
  ) %>% 
  readr::write_tsv(paste0(ofResultsDir, "/eggnog/genes2go.reduced.txt"), col_names = FALSE)

# Load the list of putative duckweed transposon- and organelle-derived coding gene predictions, and remove them from the
# annotation list.
putativeDwTeAndOrganelleGenes <- readr::read_tsv(
  "annotations/putative-te-organellar-genes/all-putative-dw-te-organelle-genes.txt.gz", 
  col_select = 1
) %>% pull(1)
annotsFilt <- annots %>% 
  dplyr::filter(! query_name %>% stringr::str_detect(pattern = paste(putativeDwTeAndOrganelleGenes, collapse = "|")))


### 
### Proteome summaries ----
### 

accTblSubgenSummary2 <- annotsFilt %>% 
  dplyr::group_by(accession) %>%
  dplyr::add_tally(name = "num_seqs") %>% 
  dplyr::filter(!is.na(eN_desc)) %>%
  dplyr::group_by(accession, num_seqs) %>%
  dplyr::tally(name = "num_eggnog_annots", sort = TRUE) %>%
  dplyr::mutate(annots_per_seq = num_eggnog_annots / num_seqs)



cor(accTblSubgenSummary2$annots_per_seq, accTblSubgenSummary2$num_seqs)
linearMod <- lm(annots_per_seq ~ num_seqs, data = accTblSubgenSummary2)


accTblSubgenSummary2 <- accTblSubgenSummary2 %>%
  dplyr::left_join(
    readr::read_tsv(
      lazy = TRUE,
      paste0(ofResultsDir, "/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv")
    ) %>%
      dplyr::rename(observation = 1) %>%
      head(n = 10) %>%
      tidyr::pivot_longer(-observation, names_to = "accession") %>%
      dplyr::mutate(accession = factor(accession)) %>%
      tidyr::pivot_wider(names_from = observation, values_from = value) %>%
      dplyr::rename_with(make.names) %>%
      dplyr::rename_with(~ stringr::str_replace(.x, "Number.of.", "")) %>%
      dplyr::rename_with(~ stringr::str_replace(.x, "Percentage.of", "percent")) %>%
      dplyr::mutate(across(-accession, as.numeric)) %>%
      dplyr::mutate(
        genes.in.common.orthogroups = genes.in.orthogroups - genes.in.species.specific.orthogroups
      )
  )

### 
### OrthoFinder HOG table ----
###
### OrthoFinder was run with hydrid subgenomes separated, yet for many comparisons, we want to consider them as a whole
### constitutive gene set, so we merge those gene lists here.
### 

# Use these getters instead of accessing the HOG table directly
getOfHogs <- function(node = "N1", subgen = FALSE, accUniverse = accAll) {
  if(is.null(.ofHogsSubgenCache[[node]])) {
    .ofHogsSubgenCache[[node]] <<- readr::read_tsv(
      lazy = TRUE,
      paste0(ofResultsDir, "/Phylogenetic_Hierarchical_Orthogroups/", node, ".tsv")
    ) %>%
      dplyr::select(where(~ any(! is.na(.)))) %>% 
      dplyr::rename("parent_clade" = "Gene Tree Parent Clade") %>%
      dplyr::mutate(dplyr::across(! c("HOG", "OG", "parent_clade"), ~ stringr::str_split(.x, pattern = ", "))) %>%
      dplyr::mutate(Ljaponica7182 = mapply(c, Ljaponica7182_M, Ljaponica7182_T, SIMPLIFY = FALSE), .before = matches("Ljaponica7182_M")) %>%
      dplyr::mutate(Ljaponica8627 = mapply(c, Ljaponica8627_M, Ljaponica8627_T, SIMPLIFY = FALSE), .before = matches("Ljaponica8627_M")) %>%
      dplyr::mutate(Ljaponica9421 = mapply(c, Ljaponica9421_M, Ljaponica9421_T, SIMPLIFY = FALSE), .before = matches("Ljaponica9421_M")) %>%
      tidyr::pivot_longer(cols = ! c("HOG", "OG", "parent_clade"), names_to = "accession", values_to = "genes") %>%
      dplyr::mutate(across(HOG:accession, as.factor)) %>%
      dplyr::group_by(HOG, accession) %>%
      dplyr::mutate(gene_counts = sum(!is.na(unlist(genes)))) %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(accession %in% accUniverse) %>%
      dplyr::mutate(accession = accession %>% droplevels())
  }
  
  tbl <- NULL
  if(subgen) {
    tbl <- .ofHogsSubgenCache[[node]] %>% 
      dplyr::filter(
        accession %in% accUniverse,
        ! accession %in% c("Ljaponica7182", "Ljaponica8627", "Ljaponica9421")
      ) %>% 
      dplyr::mutate(accession = suppressWarnings(forcats::fct_relevel(accession, accAngioSubgen)))
  } else {
    tbl <- .ofHogsSubgenCache[[node]] %>% 
      dplyr::filter(
        accession %in% accUniverse,
        ! accession %in% c("Ljaponica7182_M", "Ljaponica7182_T", "Ljaponica8627_M", "Ljaponica8627_T", 
                           "Ljaponica9421_M", "Ljaponica9421_T")
      ) %>% 
      dplyr::mutate(accession = suppressWarnings(forcats::fct_relevel(accession, accAngio)))
  }
  return(tbl) # %>% dplyr::select(where(~ any(. != 0)))
}


# HOG x accession matrix of gene counts
getOfHogsMatrix <- function(node = "N1", subgen = FALSE) {
  tbl <- getOfHogs(node = node, subgen = subgen) %>%
    dplyr::select(HOG, accession, gene_counts) %>%
    tidyr::pivot_wider(names_from = accession, values_from = gene_counts) %>%
    dplyr::mutate(HOG = as.character(HOG)) %>%
    dplyr::as_tibble()
  return(tbl)
}


#' For the given HOG matrix (tibble or data frame), return the HOGs satisfying presence/absence criteria.
#'
#' @param requiredAcc all of these accessions must have at least one gene in each HOG.
#' @param prohibitedAcc none of these accessions may have any genes in each HOGs.
#' @param optionalAcc accessions for which minOptionalAcc and maxOptionalAcc criteria will apply.
#' @param minOptionalAcc limit returned HOGs to those with at least this many optionalAcc represented. Can be different
#' accessions in different HOGs.
#' @param maxOptionalAcc limit returned HOGs to those with at most this many optionalAcc represented. Can be different
#' accessions in different HOGs.
#' @param hogMatrix HOG x accessions matrix with accession values coercible to logical (e.g. integer gene
#' counts or boolean)
#' @return vector of HOGs satisfying these conditions.
filterHogs <- function(
    requiredAcc = NULL, 
    prohibitedAcc = NULL, 
    optionalAcc = NULL, 
    minOptionalAcc = 0,
    maxOptionalAcc = length(optionalAcc),
    hogMatrix) {
  
  filtHogs <- hogMatrix %>% 
    dplyr::mutate(across(
      c(all_of(requiredAcc), all_of(prohibitedAcc), all_of(optionalAcc)), 
      ~ 1L * (as.integer(.x) != 0))
    ) %>% 
    dplyr::mutate(
      sum_req = rowSums(across(all_of(requiredAcc))),
      sum_pro = rowSums(across(all_of(prohibitedAcc))),
      sum_opt = rowSums(across(all_of(optionalAcc)))
    ) %>% 
    dplyr::filter(
      sum_req == length(requiredAcc),
      sum_pro == 0,
      sum_opt >= minOptionalAcc & sum_opt <= maxOptionalAcc
    ) %>% 
    dplyr::select(-sum_req, -sum_pro, -sum_opt)
    
  return(filtHogs %>% dplyr::pull(HOG))
}


#' Return a vector of all gene IDs in a given HOG for the given accessions.
#'
#' @param selectedHogs vector of HOGs of interest.
#' @param accUniverse only genes belonging to the accessions in this vector will be returned.
#' @param hogTable a table such as that returned by @function{getOfHogs}.
#' @return vector of gene IDs.
genesFromHogs <- function(selectedHogs, accUniverse = accAngio, hogTable) {
  genes <- hogTable %>% 
    dplyr::filter(accession %in% accUniverse) %>%
    dplyr::filter(HOG %in% selectedHogs) %>% 
    dplyr::pull(genes) %>% 
    unlist()
  genes <- genes[!is.na(genes)]
  return(genes)
}


#' Create a table, one row per HOG, including geneIDs for each accession and eggNOG and AHRD annotations.
#'
#' @param accUniverse the subset of accessions to include in the table.
#' @param selectedHogs the subset of HOGs to output.
#' @param annTbl the annotation table imported from eggNOG mapper.
#' @return an annotated HOG table.
getAnnotatedHogs <- function(
    selectedHogs = Inf,
    node = "N1",
    subgen = FALSE,
    accUniverse = accAngioSubgen,
    annTbl = annots
) {
  if(is.null(selectedHogs)) {
    return("")
  }
  
  annotFile <- paste0(hogTableOutDir, "/hogs.", node, ".annotated.tsv.gz")
  
  if(!is.null(.ofHogsSubgenAnnotCache[[node]])) {
    # First check the table cache
    print(paste0("Loading annotations from the cache for node ", node))
    tbl <- .ofHogsSubgenAnnotCache[[node]]
  } else if(file.exists(annotFile)) {
    # No hit, so try to read from disk
    print(paste0("Loading annotations from disk for node ", node))
    .ofHogsSubgenAnnotCache[[node]] <<- readr::read_tsv(
      lazy = TRUE,
      annotFile
    )
    tbl <- .ofHogsSubgenAnnotCache[[node]]
  } else {
    # Can't find the full table for this node, so let's build it from scratch (will take many minutes!!!).
    print(paste0("Building annotations from scratch for node ", node))
    # Starting point: a table with one row per HOG per accession. Build this table with subgenomes separated as in the
    # original OrthoFinder run, and leave the responsibility for returning the correct representation to post-processing
    # code below.
    hogTbl <- getOfHogs(node = node, subgen = TRUE) %>% dplyr::select(HOG, OG, accession, genes)
    
    # Now lengthen the table to include one row per gene in each HOG and lookup the annotations for each gene ID in the 
    # annotation table. Create new columns to store the Arabidopsis, rice, and maize orthologues of each gene and 
    # retrieve their gene symbols from their respective annotations.
    longHogTbl <- hogTbl %>% 
      tidyr::unnest(genes) %>% 
      dplyr::left_join(
        annTbl %>% dplyr::select(
          query_name, 
          eN_seed_olog, 
          eN_seed_olog_evalue, 
          eN_OGs,
          Preferred_name, 
          GOs,
          eN_desc,
          ahrd_desc
        ), 
        by = c("genes" = "query_name")
      ) %>% 
      dplyr::mutate(
        athal_locus = dplyr::if_else(
          accession == "Athaliana", 
          stringr::str_extract(genes, "^AT[^.]*"), 
          as.character(NA)
        ),
        rice_locus = dplyr::if_else(
          accession == "Osativa", 
          stringr::str_extract(genes, "^LOC[^.]*"), 
          as.character(NA)
        ),
        maize_locus = dplyr::if_else(
          accession == "Zmays", 
          stringr::str_extract(genes, "^Zm00001d[^.]*") %>% stringr::str_remove("_P[0-9]{3,3}$"), 
          as.character(NA)
        )
      ) %>% 
      dplyr::left_join(
        athalGeneTbl %>% dplyr::select(-full_name) %>% dplyr::rename(athal_symbol = symbol),
        by = c("athal_locus" = "locus_name")
      ) %>% 
      dplyr::left_join(
        osatiGeneTbl %>% dplyr::select(-full_name) %>% dplyr::rename(osati_symbol = symbol),
        by = c("rice_locus" = "locus_name")
      ) %>% 
      dplyr::left_join(
        zmaysGeneTbl %>% dplyr::select(-full_name) %>% dplyr::rename(zmays_symbol = symbol),
        by = c("maize_locus" = "locus_name")
      ) %>% 
      dplyr::mutate(
        eN_seed_olog_evalue = eN_seed_olog_evalue %>% stringr::str_replace_na(replacement = "Inf") %>% as.numeric(),
        athal_symbol = if_else(is.na(athal_symbol), athal_locus, athal_symbol),
        osati_symbol = if_else(is.na(osati_symbol), rice_locus, osati_symbol),
        zmays_symbol = if_else(is.na(zmays_symbol), maize_locus, zmays_symbol)
      ) %>% 
      dplyr::select(! dplyr::all_of(c("athal_locus", "rice_locus", "maize_locus")))
    
    # Widen the data again, collapsing (by concatenation) the genes and their annotations.
    summaryTbl <- longHogTbl %>% 
      tidyr::unnest(GOs, keep_empty = TRUE) %>% 
      dplyr::group_by(HOG) %>% 
      dplyr::summarize(
        At_syms = athal_symbol %>% na.omit() %>% unique() %>% stringr::str_flatten(collapse = ","),
        Os_syms = osati_symbol %>% na.omit() %>% unique() %>% stringr::str_flatten(collapse = ","),
        Zm_syms = zmays_symbol %>% na.omit() %>% unique() %>% stringr::str_flatten(collapse = ","),
        EN_seeds = eN_seed_olog %>% na.omit() %>% unique() %>% stringr::str_flatten(collapse = ","),
        # Note: extracting only the plant-level eggNOG OGs (Tax ID 33090)
        EN_OGs = eN_OGs %>% stringr::str_split(pattern = ",") %>% unlist() %>% 
          stringr::str_extract(pattern = ".*@33090") %>% na.omit() %>% uniqueCaseInsensitive() %>% 
          stringr::str_flatten(collapse = ","),
        EN_OGs_At_syms = EN_OGs %>% athSymbolsInEggNogOg(),
        EN_descs = eN_desc %>% stringr::str_replace("^-$", NA_character_) %>% na.omit() %>% uniqueCaseInsensitive() %>% 
          stringr::str_flatten(collapse = ","),
        AHRD_descs = ahrd_desc %>% na.omit() %>% uniqueCaseInsensitive() %>% stringr::str_flatten(collapse = ","),
        EN_GOs = names(GOs) %>% unique() %>% stringi::stri_remove_empty_na() %>% stringr::str_flatten(collapse = ",")
      )
    
    # Where there are multiple annotated genes per HOG, pick the EggNog annotation with the best
    # e-value as the representative for that HOG.
    bestMatchTbl <- longHogTbl %>%
      dplyr::group_by(HOG) %>% 
      slice_min(n = 1, order_by = eN_seed_olog_evalue, with_ties = FALSE) %>% 
      droplevels() %>% 
      dplyr::select(HOG, OG, genes, eN_seed_olog_evalue, eN_seed_olog) %>% 
      dplyr::rename(
        EN_best_query = genes, 
        EN_best_seed = eN_seed_olog, 
        EN_best_seed_eval = eN_seed_olog_evalue
      ) %>% 
      dplyr::ungroup()
    
    # Finally, combine the best matches and the summaries for each HOG to form the final output.
    tbl <- hogTbl %>%
      dplyr::group_by(HOG, accession) %>% 
      dplyr::summarize(genes = genes %>% unlist() %>% stringr::str_flatten(collapse = ",")) %>% 
      tidyr::pivot_wider(names_from = accession, values_from = genes) %>% 
      dplyr::left_join(bestMatchTbl, by = c("HOG" = "HOG")) %>% 
      dplyr::left_join(summaryTbl, by = c("HOG" = "HOG")) %>%
      dplyr::mutate(
        dplyr::across(
          .cols = dplyr::any_of(accSubgen), 
          .names = "{.col}_n",
          ~ stringr::str_count(.x, "[^,]+(,)?") %>% na.replace(0)
        )
      ) %>% 
      dplyr::relocate(
        OG, 
        HOG, 
        dplyr::any_of(paste0(accAll, "_n")),
        At_syms, 
        Os_syms, 
        Zm_syms,
        EN_OGs_At_syms,
        EN_descs,
        AHRD_descs,
        EN_OGs,
        EN_best_query, 
        EN_best_seed, 
        EN_best_seed_eval, 
        EN_seeds,
        dplyr::any_of(accAll)
      ) %>% 
      dplyr::ungroup()
    
    # Save our work!!!
    readr::write_tsv(
      tbl,
      annotFile
    )
    
    .ofHogsSubgenAnnotCache[[node]] <- tbl
  }
  
  # Return only the selected HOG subset if specified.
  if(!identical(selectedHogs, Inf)) {
    tbl <- tbl %>% 
      dplyr::filter(HOG %in% selectedHogs)
  }
  
  if(subgen) {
    tbl <- tbl %>%
      dplyr::select(! matches(ignore.case = FALSE, c("^Ljaponica[0-9]{4}$", "^Ljaponica[0-9]{4}_n$")))
  } else {
    tbl <- tbl %>% 
      rowwise() %>% 
      dplyr::mutate(
        Ljaponica7182 = c(
          Ljaponica7182_M %>% stringr::str_split(pattern = ",") %>% unlist(), 
          Ljaponica7182_T %>% stringr::str_split(pattern = ",") %>% unlist()
          ) %>% naturalsort::naturalsort() %>% stringr::str_flatten(collapse = ",") %>% dplyr::na_if(""), 
        Ljaponica8627 = c(
          Ljaponica8627_M %>% stringr::str_split(pattern = ",") %>% unlist(), 
          Ljaponica8627_T %>% stringr::str_split(pattern = ",") %>% unlist()
        ) %>% naturalsort::naturalsort() %>% stringr::str_flatten(collapse = ",") %>% dplyr::na_if(""),
        Ljaponica9421 = c(
          Ljaponica9421_M %>% stringr::str_split(pattern = ",") %>% unlist(), 
          Ljaponica9421_T %>% stringr::str_split(pattern = ",") %>% unlist()
        ) %>% naturalsort::naturalsort() %>% stringr::str_flatten(collapse = ",") %>% dplyr::na_if(""),
        .before = matches("^Ljaponica7182_M$")
      ) %>%
      as_tibble() %>%
      dplyr::mutate(
        dplyr::across(
          .cols = dplyr::any_of(accSets$Ljaponicas), 
          .names = "{.col}_n",
          ~ stringr::str_count(.x, "[^,]+(,)?") %>% na.replace(0)
        ),
        .after = matches("^Lminor9252_n$")
      ) %>%
      dplyr::select(! matches(ignore.case = FALSE, c("^Ljaponica[0-9]{4}_[A-Z]$", "^Ljaponica[0-9]{4}_[A-Z]_n$")))
  }

  tbl <- tbl %>% 
    dplyr::arrange(
      At_syms == "", 
      Os_syms == "", 
      Zm_syms == "",
      EN_OGs_At_syms == "", 
      EN_descs == "", 
      AHRD_descs == "", 
      OG
    )
  
  return(tbl)
}

#' Add columns of GO IDs and GO Terms to an already annotated HOG table.
#'
#' @param annHogTbl the annotated HOG table to decorate.
#' @param goResultsBp a table of GO terms enriched in this subset of HOGs as returned by goTermEnrichment()
#' @param goResultsMf a table of GO terms enriched in this subset of HOGs as returned by goTermEnrichment()
#' @param goResultsCc a table of GO terms enriched in this subset of HOGs as returned by goTermEnrichment()
#' @return the input annotated HOG table decorated with the go results.
addGoResultsToHogTbl <- function(
    annHogTbl,
    goResultsBp = NA,
    goResultsMf = NA,
    goResultsCc = NA
) {
  # Get the char vectors from the tibbles, sort
  topGosBp <- if(identical(goResultsBp, NA)) NA else sort(goResultsBp %>% pull(GO.ID))
  topGosMf <- if(identical(goResultsMf, NA)) NA else sort(goResultsMf %>% pull(GO.ID))
  topGosCc <- if(identical(goResultsCc, NA)) NA else sort(goResultsCc %>% pull(GO.ID))
  
  tbl <- annHogTbl %>% left_join(
    by = "HOG", 
    annHogTbl %>% 
      dplyr::mutate(GOs_list = EN_GOs %>% stringr::str_split(pattern = ",")) %>% 
      tidyr::unnest(GOs_list, keep_empty = TRUE) %>%
      dplyr::group_by(HOG) %>% 
      dplyr::summarize(
        Top_GOs_BP = GOs_list %>% unlist() %>% unique() %>% naturalsort::naturalsort() %>% intersect(y = topGosBp) %>% 
          stringr::str_flatten(collapse = ","),
        Top_GOs_MF = GOs_list %>% unlist() %>% unique() %>% naturalsort::naturalsort() %>% intersect(y = topGosMf) %>% 
          stringr::str_flatten(collapse = ","),
        Top_GOs_CC = GOs_list %>% unlist() %>% unique() %>% naturalsort::naturalsort() %>% intersect(y = topGosCc) %>% 
          stringr::str_flatten(collapse = ","),
        Top_GO_terms_BP = GOs_list %>% unlist() %>% unique() %>% naturalsort::naturalsort() %>% intersect(y = topGosBp) %>% 
          na.replace("") %>% mAnnDbiTerm() %>% stringr::str_flatten(collapse = ","),
        Top_GO_terms_MF = GOs_list %>% unlist() %>% unique() %>% naturalsort::naturalsort() %>% intersect(y = topGosMf) %>% 
          na.replace("") %>% mAnnDbiTerm() %>% stringr::str_flatten(collapse = ","),
        Top_GO_terms_CC = GOs_list %>% unlist() %>% unique() %>% naturalsort::naturalsort() %>% intersect(y = topGosCc) %>% 
          na.replace("") %>% mAnnDbiTerm() %>% stringr::str_flatten(collapse = ",")
      )
  ) %>%       
    dplyr::relocate(
      .after = AHRD_descs,
      any_of("Top_GO_terms_BP"), 
      any_of("Top_GO_terms_MF"),
      any_of("Top_GO_terms_CC")
    )
  
  if("Top_GO_terms_BP" %in% colnames(tbl)) {
    tbl <- tbl %>% 
      dplyr::arrange(
        Top_GO_terms_BP == "", 
        At_syms == "", 
        Os_syms == "", 
        Zm_syms == "",
        EN_OGs_At_syms == "", 
        EN_descs == "", 
        AHRD_descs == "", 
        OG
      )
  } else {
    tbl <- tbl %>% 
      dplyr::arrange(
        At_syms == "", 
        Os_syms == "", 
        Zm_syms == "",
        EN_OGs_At_syms == "", 
        EN_descs == "", 
        AHRD_descs == "", 
        OG
      )
  }
  
  return(tbl)
}


#' Remove entries from an annotated HOG table.
#'
#' @param annHogTbl annotated HOG table.
#' @param geneIds character vector of gene IDs. If this is provided, HOGs containing any of these will be purged.
#' @param purgeAtOrganellar purge HOGs containing Athaliana plastid or mitochondrial gene IDs.
#' @return the input table without organellar or transposon HOGs.
purgeHogs <- function(annHogTbl, geneIds = NULL, purgeAtOrganellar = FALSE) {
  hogTblFilt <- annHogTbl
  geneIdsPtn <- paste(geneIds, collapse = "|")
  if (length(geneIds) > 0) {
    hogTblFilt <- hogTblFilt %>% 
      dplyr::filter(
        ! dplyr::if_any(
          dplyr::any_of(accAll), 
          ~ stringr::str_replace_na(.) %>% stringr::str_detect(pattern = geneIdsPtn)
        )
      )
  }
  if (purgeAtOrganellar) {
    hogTblFilt <- hogTblFilt %>% 
      dplyr::filter(
        ! EN_seeds %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "AT[MC]G[0-9]{5}\\.[0-9]"),
        ! Athaliana %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "AT[MC]G[0-9]{5}\\.[0-9]"),
      )
  }
  return(hogTblFilt)
}


getOgsFromHogs <- function(hogs, node, subgen = FALSE) {
  ogs <- getOfHogs(node = node, subgen = subgen) %>% dplyr::filter(HOG %in% hogs) %>% dplyr::pull(OG) %>% as.character() %>% unique()
  return(ogs)
}


getHogsInOgs <- function(ogs, node, subgen = FALSE) {
  hogs <- getOfHogs(node = node, subgen = subgen) %>% dplyr::filter(OG %in% ogs) %>% dplyr::pull(HOG) %>% as.character() %>% unique()
  return(hogs)
}


getHogsInOgsWithAccs <- function(ogs, queryAccs, node, subgen = FALSE) {
  hogs <- getOfHogs(node = node, subgen = subgen) %>% dplyr::filter(OG %in% ogs)
  hogsFilt <- hogs %>% 
    dplyr::group_by(HOG, accession) %>% 
    dplyr::summarize(.groups = "drop_last", gene_ct = length(genes %>% unlist())) %>% 
    dplyr::filter(accession %in% queryAccs & gene_ct > 0) %>% 
    dplyr::pull(HOG) %>% 
    as.character() %>% 
    unique()
  return(hogsFilt)
}


#' Output HOG tables in Excel and tsv formats. The Excel workbook will only contain those HOGs that have an EggNog 
#' match, an Arabidopsis thaliana symbol, or functional annotations from AHRD.
#' 
#' @param hogTables named list of HOG tables. Names will be used to separate and name Excel sheets, one per table.
#' @param tblName filename prefix for HOG tables.
#' @param outDir directory to write HOG tables to.
#' @param tblNameAppendDate if true append the date to the tblName.
writeHogTables <- function(hogTables, tblName = "hog-table", outDir = ".", tblNameAppendDate = FALSE) {
  if(tblNameAppendDate) {
    dateString <- format(Sys.time(), "%Y%m%d")
    tblName <- paste0(tblName, ".", dateString)
  }

  # Enforce some prettifying things
  hogTables <- purrr::map(
    hogTables,
    .f = ~ {
      .x %>% dplyr::rename_with(
        .cols = dplyr::matches("^Spolyrhiza([^0-9]+.*)?$"), 
        .fn = ~ .x %>% gsub(pattern = "Spolyrhiza", replacement = "Spolyrhiza7498")
      )
    }
  )
  
  tsvOutPath <- paste0(outDir, "/", tblName, ".txt.gz")
  xlsxOutPath <- paste0(outDir, "/", tblName, ".xlsx")
  
  # First we write the TSV
  readr::write_tsv(
    c("comparison", hogTables[[1]] %>% colnames()) %>% as_tibble_row(.name_repair = ~ paste0("c", seq(1, length(.)))),
    tsvOutPath,
    col_names = FALSE
  )
  purrr::imap(
    hogTables,
    .f = ~ {
      readr::write_tsv(
        append = TRUE,
        .x %>% mutate(comparison = .y) %>% dplyr::relocate(comparison),
        tsvOutPath
      )
    }
  )
  
  # And now the Excel workbook
  hogWb <- openxlsx::createWorkbook(creator = "Evan Ernst", title = tblName)
  
  mapply(
    function(hogTbl, sheetName, workBook) {
      # !!!
      # Remove rows from each HOG table that have none of functional annotations from AHRD, EggNOG orthogroup 
      # assignments, or symbols from A. thaliana, rice or maize
      # !!!
      hogTbl <- hogTbl %>% 
        dplyr::filter(! (is.na(At_syms) & is.na(Os_syms) & is.na(Zm_syms) & is.na(AHRD_descs) & is.na(EN_OGs)))
      
      workBook %>% openxlsx::addWorksheet(sheetName)
      workBook %>% openxlsx::writeData(sheet = sheetName, hogTbl)
      workBook %>% openxlsx::pageSetup(sheet = sheetName, paperSize = 1)
      
      hogTblNCol <- ncol(hogTbl)
      hogTblNRow <- nrow(hogTbl)
      
      headerStyle <- openxlsx::createStyle(
        halign = "CENTER", 
        textDecoration = "Bold",
        border = "Bottom"
      )
      
      verticalHeaderStyle <- openxlsx::createStyle(
        textRotation = -60
      )
      
      # Column widths
      openxlsx::setColWidths(
        workBook, 
        sheetName, 
        cols = 1:hogTblNCol,
        widths = 12
      )
      openxlsx::setColWidths(
        workBook, 
        sheetName, 
        cols = c(
          colnames(hogTbl) %>% grep(pattern = "^contrast$")
        ),
        widths = 16
      )
      openxlsx::setColWidths(
        workBook, 
        sheetName, 
        cols = c(
          colnames(hogTbl) %>% grep(pattern = "^OG$")
        ),
        widths = 10
      )
      openxlsx::setColWidths(
        workBook, 
        sheetName, 
        cols = c(
          colnames(hogTbl) %>% grep(pattern = "^HOG$")
        ),
        widths = 13.5
      )
      openxlsx::setColWidths(
        workBook, 
        sheetName, 
        cols = c(
          colnames(hogTbl) %>% grep(pattern = "^EN_OGs_At_syms$"),
          colnames(hogTbl) %>% grep(pattern = "^EN_descs$"),
          colnames(hogTbl) %>% grep(pattern = "^AHRD_descs$")
        ),
        widths = 16
      )
      openxlsx::setColWidths(
        workBook, 
        sheetName, 
        cols = c(
          colnames(hogTbl) %>% grep(pattern = "^Top_GO_terms_")
        ),
        widths = 15
      )
      
      # Headers
      openxlsx::addStyle(
        workBook,
        sheetName,
        style = headerStyle,
        rows = 1,
        cols = 1:hogTblNCol,
        stack = TRUE
      )
      openxlsx::addFilter(
        workBook, 
        sheetName, 
        rows = 1,
        cols = c(
          colnames(hogTbl) %>% grep(pattern = "^contrast$")
        )
      )
      openxlsx::freezePane(workBook, sheetName, firstRow = TRUE)
      
      # Special formatting for gene copy number columns
      openxlsx::setColWidths(
        workBook, 
        sheetName, 
        cols = c(
          colnames(hogTbl) %>% grep(pattern = "^[A-z][a-z0-9]+(_[A-Z])?_n$")
        ),
        widths = 4
      )
      openxlsx::addStyle(
        workBook,
        sheetName,
        style = verticalHeaderStyle,
        rows = 1,
        cols = c(colnames(hogTbl) %>% grep(pattern = "^[A-z][a-z0-9]+(_[A-Z])?_n$")),
        stack = TRUE
      )
      openxlsx::addStyle(
        workBook,
        sheetName,
        style = openxlsx::createStyle(fgFill = pal4["Spolyrhiza9509"]),
        rows = 1,
        cols = c(colnames(hogTbl) %>% grep(pattern = "^S(polyrhiza|intermedia)[0-9]*_n$")),
        stack = TRUE
      )
      openxlsx::addStyle(
        workBook,
        sheetName,
        style = openxlsx::createStyle(fgFill = pal4["Lgibba7742a"]),
        rows = 1,
        cols = c(colnames(hogTbl) %>% grep(pattern = "^Lgibba[0-9]{4}[a-z]?_n$")),
        stack = TRUE
      )
      openxlsx::addStyle(
        workBook,
        sheetName,
        style = openxlsx::createStyle(fgFill = pal4["Lminor7210"]),
        rows = 1,
        cols = c(colnames(hogTbl) %>% grep(pattern = "^Lminor[0-9]{4}_n$")),
        stack = TRUE
      )
      openxlsx::addStyle(
        workBook,
        sheetName,
        style = openxlsx::createStyle(fgFill = pal4["Lminor7210"]),
        rows = 1,
        cols = c(colnames(hogTbl) %>% grep(pattern = "^Ljaponica[0-9]{4}_M_n$")),
        stack = TRUE
      )
      openxlsx::addStyle(
        workBook,
        sheetName,
        style = openxlsx::createStyle(fgFill = pal4["Ljaponica7182"]),
        rows = 1,
        cols = c(colnames(hogTbl) %>% grep(pattern = "^Ljaponica[0-9]{4}_n$")),
        stack = TRUE
      )
      openxlsx::addStyle(
        workBook,
        sheetName,
        style = openxlsx::createStyle(fgFill = pal4["Lturionifera9434"]),
        rows = 1,
        cols = c(colnames(hogTbl) %>% grep(pattern = "^Lturionifera[0-9]{4}_n$")),
        stack = TRUE
      )
      openxlsx::addStyle(
        workBook,
        sheetName,
        style = openxlsx::createStyle(fgFill = pal4["Lturionifera9434"]),
        rows = 1,
        cols = c(colnames(hogTbl) %>% grep(pattern = "^Ljaponica[0-9]{4}_T_n$")),
        stack = TRUE
      )
      openxlsx::addStyle(
        workBook,
        sheetName,
        style = openxlsx::createStyle(fgFill = pal4["Waustraliana8730"]),
        rows = 1,
        cols = c(colnames(hogTbl) %>% grep(pattern = "^Waustraliana[0-9]{4}_n$")),
        stack = TRUE
      )
      openxlsx::addStyle(
        workBook,
        sheetName,
        style = openxlsx::createStyle(halign = "CENTER"),
        rows = 1:(hogTblNRow + 1),
        cols = c(colnames(hogTbl) %>% grep(pattern = "^[A-z][a-z0-9]+(_[A-Z])?_n$")),
        gridExpand = TRUE,
        stack = TRUE
      )
      openxlsx::conditionalFormatting(
        workBook,
        sheetName,
        cols = colnames(hogTbl) %>% grep(pattern = "^[A-z][a-z0-9]+_n$"),
        rows = 1 + (1:nrow(hogTbl)),
        style = c("#6699CC", "#FFFFFF", "#FF9999"),
        rule = c(0, 1, 3),
        type = "colorScale"
      )
      
      return(workBook)
    },
    hogTables,
    names(hogTables),
    MoreArgs = list(workBook = hogWb),
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE
  )

  openxlsx::saveWorkbook(hogWb, xlsxOutPath, overwrite = TRUE)
}



### 
### OrthoFinder Duplications table ----
###
.ofDupsCache <- NULL
getOfDups <- function(node = NULL, minSupport = 0) {
  if (is.null(.ofDupsCache)) {
    dupTbl <- readr::read_tsv(
      name_repair = "universal",
      lazy = TRUE,
      paste0(ofResultsDir, "/Gene_Duplication_Events/Duplications.tsv")
    )
    .ofDupsCache <<- dupTbl
  }
  if (is.null(node)) {
    filteredDups <- .ofDupsCache
  } else {
    filteredDups <- .ofDupsCache %>% dplyr::filter(Species.Tree.Node == node)
  }
  filteredDups <- filteredDups %>% dplyr::filter(Support >= minSupport)
  return(filteredDups)
}
getOfDupsLoci <- function(node = NULL, minSupport = 0) {
  dupTbl <- getOfDups(node, minSupport) %>% 
    .$Genes.2 %>% purrr::map(~ stringr::str_split(., pattern = ",") %>% unlist() %>% stringr::str_trim())
  return(dupTbl)
}



# 
# HOG intersections/subsets ----
#

# HOG table globals
hogTableOutDir <- "hog-tables"
.ofHogsSubgenCache <- list()
.ofHogsSubgenAnnotCache <- list()
.ofHogsToOg <- list()

# Make sure these have been generated and cached
getAnnotatedHogs(node = "N0")
getAnnotatedHogs(node = "N1")
getAnnotatedHogs(node = "N4")


## All Angiosperms ----
# Genes missing and specific across accession groups of Lemnaceae genera, species, and others, relative to the remaining
# angiosperms in the analysis. 
selectedSetsAngio <- accSetsAngio[c(
"Lemnaceae",
"Spirodelas",
"Lemnas",
"Wolffias",
"Spolyrhizas",
"Lgibbas",
"Lminors",
"Ljaponicas",
"Lturioniferas",
"Waustralianas",
"aquatic_floating_dw",
"aquatic_submerged_dw",
"aquatic_SF_dw",
"aquatic_rootless",
#"aquatic_sinking", # this resulted in 0 HOGs
#"freshwater_sinking", # this resulted in 0 HOGs
"cdemersum_dw",
"zmarina_dw",
"monocots",
"dicots"
)]
# Sets of accessions to compute specific HOGs over, yielding e.g. HOGs which are present in all 
# "spirodelas" but missing from all other angiosperm accessions.
specificSetsAngio <- selectedSetsAngio
names(specificSetsAngio) <- names(selectedSetsAngio) %>% paste0("_spec")
# The set complements of the specific sets, yielding e.g. HOGs which are absent in "Spirodelas" 
# but present in all other angiosperm accessions.
missingSetsAngio <- purrr::map(
  specificSetsAngio,
  ~ setdiff(accAngio, .x)
)
names(missingSetsAngio) <- names(selectedSetsAngio) %>% paste0("_miss")
# interleave the specific and missing sets
specificAndMissingSetsAngio <- (c(missingSetsAngio, specificSetsAngio))[order(
  c(seq_along(missingSetsAngio), seq_along(specificSetsAngio)))]

# Get the lists of HOGs exclusive to each predefined set
filteredHogSetsAngio <- purrr::map(
  specificAndMissingSetsAngio, 
  ~ filterHogs(
    requiredAcc = .x, 
    prohibitedAcc = setdiff(accAngio, .x),
    hogMatrix = getOfHogsMatrix()
  )
)
# Get the gene IDs for the exclusive HOGs
selectedHogsGeneIdsAngio <- lapply(
  filteredHogSetsAngio,
  genesFromHogs,
  accUniverse = accAngio,
  hogTable = getOfHogs()
)

## Ath1Mono Missing ----
# Genes missing and specific from lemnaceae but present in both athal and one other monocot, thus 
# a looser bound than requiring them to be present in all other accessions.
accUniverseAth1Mono <- c(accSets$monocots, "Athaliana")
prohibitedSetsAth1Mono <- accSetsAngio[c(
  "Lemnaceae",
  "Spirodelas",
  "Lemnas",
  "Wolffias",
  "Spolyrhizas",
  "Lgibbas",
  "Lminors",
  "Ljaponicas",
  "Lturioniferas",
  "Waustralianas",
  "aquatic_rootless",
  "aquatic_sinking",
  "freshwater_sinking",
  "cdemersum_dw",
  "zmarina_dw"
)]
names(prohibitedSetsAth1Mono) <- names(prohibitedSetsAth1Mono) %>% paste0("_miss_ath_1mono")
# Get the lists of HOGs exclusive to each predefined set
filteredHogSetsAth1Mono <- purrr::map(
  prohibitedSetsAth1Mono, 
  ~ filterHogs(
    requiredAcc = c("Athaliana", setdiff(accDw, .x)),
    prohibitedAcc = .x,
    optionalAcc = setdiff(accSets$monocots, c(accDw)),
    minOptionalAcc = 1,
    hogMatrix = getOfHogsMatrix(node = "N4")
  )
)
# Get the gene IDs for the exclusive HOGs
selectedHogsGeneIdsAth1Mono <- lapply(
  filteredHogSetsAth1Mono,
  genesFromHogs,
  accUniverse = accUniverseAth1Mono,
  hogTable = getOfHogs(node = "N4")
)

## Ath-Osa Missing ----
# Genes missing and specific from lemnaceae but present in both athal and osati, thus 
# a looser bound than requiring them to be present in all other accessions.
accUniverseAthOsa <- c(accDw, "Cdemersum", "Zmarina", "Athaliana", "Osativa")
selectedSetsAthOsa <- accSetsAngio[c(
  "Lemnaceae",
  "Spirodelas",
  "Lemnas",
  "Wolffias",
  "Spolyrhizas",
  "Lgibbas",
  "Lminors",
  "Ljaponicas",
  "Lturioniferas",
  "Waustralianas"
)]
specificSetsAthOsa <- selectedSetsAthOsa
missingSetsAthOsa <- purrr::map(
  specificSetsAthOsa,
  ~ setdiff(accUniverseAthOsa, .x)
)
names(missingSetsAthOsa) <- names(selectedSetsAthOsa) %>% paste0("_miss_ath_osa")
names(specificSetsAthOsa) <- names(selectedSetsAthOsa) %>% paste0("_spec_ath_osa")
specificAndMissingSetsAthOsa <- (c(missingSetsAthOsa, specificSetsAthOsa))[order(
  c(seq_along(missingSetsAthOsa), seq_along(specificSetsAthOsa)))]
# Get the lists of HOGs exclusive to each predefined set
filteredHogSetsAthOsa <- purrr::map(
  specificAndMissingSetsAthOsa,
  ~ filterHogs(
    requiredAcc = .x,
    prohibitedAcc = setdiff(accUniverseAthOsa, .x),
    hogMatrix = getOfHogsMatrix(node = "N4")
  )
)  %>% 
  append(list(
    "aquatic_rootless_miss_ath_osa" = filterHogs(
      requiredAcc = c("Athaliana", "Osativa"),
      prohibitedAcc = accSetsAngio$aquatic_rootless,
      hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
    "aquatic_rootless_spec_ath_osa" = filterHogs(
      requiredAcc = accSetsAngio$aquatic_rootless,
      prohibitedAcc = c("Athaliana", "Osativa"),
      hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
    "cdemersum_dw_miss_ath_osa" = filterHogs(
      requiredAcc = c("Athaliana", "Osativa"),
      prohibitedAcc = accSetsAngio$cdemersum_dw,
      hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
    "cdemersum_dw_spec_ath_osa" = filterHogs(
      requiredAcc = accSetsAngio$cdemersum_dw,
      prohibitedAcc = c("Athaliana", "Osativa"),
      hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
    "zmarina_dw_miss_ath_osa" = filterHogs(
      requiredAcc = c("Athaliana", "Osativa"),
      prohibitedAcc = accSetsAngio$zmarina_dw,
      hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
    "zmarina_dw_spec_ath_osa" = filterHogs(
      requiredAcc = accSetsAngio$zmarina_dw,
      prohibitedAcc = c("Athaliana", "Osativa"),
      hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4"))
  ))
# Get the gene IDs for the exclusive HOGs
selectedHogsGeneIdsAthOsa <- lapply(
  filteredHogSetsAthOsa,
  genesFromHogs,
  accUniverse = accUniverseAthOsa,
  hogTable = getOfHogs(node = "N4")
)

## Intra-Lemnaceae ----
# Intra-Lemnaceae comparison, hybrid subgenomes treated separately
accUniverseLemnaceae <- accDw
selectedSetsLemnaceae <- accSetsAngio[c(
  "Spirodelas",
  "Lemnas",
  "Wolffias",
  "Spolyrhizas",
  "Lgibbas",
  "Ljaponicas",
  "Lminors",
  "Lturioniferas",
  "Waustralianas",
  "Spoly_and_Lturionifera",
  "sinking_dw"
)]
specificSetsLemnaceae <- selectedSetsLemnaceae
missingSetsLemnaceae <- purrr::map(
  specificSetsLemnaceae,
  ~ setdiff(accUniverseLemnaceae, .x)
)
names(missingSetsLemnaceae) <- names(selectedSetsLemnaceae) %>% paste0("_miss")
names(specificSetsLemnaceae) <- names(selectedSetsLemnaceae) %>% paste0("_spec")
specificAndMissingSetsLemnaceae <- (c(missingSetsLemnaceae, specificSetsLemnaceae))[order(
  c(seq_along(missingSetsLemnaceae), seq_along(specificSetsLemnaceae)))]
# Get the lists of HOGs exclusive to each predefined set
filteredHogSetsLemnaceae <- c(
  list("Lemnaceae_variable" = filterHogs(
    optionalAcc = accDw,
    minOptionalAcc = 1,
    maxOptionalAcc = length(accDw) - 1,
    hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4"))),
  purrr::map(
    specificAndMissingSetsLemnaceae, 
    ~ filterHogs(
      requiredAcc = .x, 
      prohibitedAcc = setdiff(accUniverseLemnaceae, .x),
      hogMatrix = getOfHogsMatrix(node = "N4")
    )
  )
)
# Get the gene IDs for the exclusive HOGs
selectedHogsGeneIdsLemnaceae <- lapply(
  filteredHogSetsLemnaceae,
  genesFromHogs,
  accUniverse = accUniverseLemnaceae,
  hogTable = getOfHogs(subgen = TRUE, node = "N4")
)

## Hybrids ----
# Hybrid subgenomes independent comparisons vs. diploids
filteredHogSetsHybrids <-
    list(
      "Ljaponica_variable" = filterHogs(
        optionalAcc = c(accSets$Lminors, accSets$Lturioniferas, accSets$Ljaponicas),
        minOptionalAcc = 1,
        maxOptionalAcc = length(c(accSets$Lminors, accSets$Lturioniferas, accSets$Ljaponicas)) - 1,
        hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
      "Ljaponica7182_miss_lemnas" = filterHogs(
        requiredAcc = setdiff(accSets$Lemnas, "Ljaponica7182"),
        prohibitedAcc = "Ljaponica7182",
        hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
      "Ljaponica7182_spec_lemnas" = filterHogs(
        requiredAcc = "Ljaponica7182",
        prohibitedAcc = setdiff(accSets$Lemnas, "Ljaponica7182"), 
        hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
      "Ljaponica8627_miss_lemnas" = filterHogs(
        requiredAcc = setdiff(accSets$Lemnas, "Ljaponica8627"),
        prohibitedAcc = "Ljaponica8627",
        hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
      "Ljaponica8627_spec_lemnas" = filterHogs(
        requiredAcc = "Ljaponica8627",
        prohibitedAcc = setdiff(accSets$Lemnas, "Ljaponica8627"), 
        hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
      "Ljaponica9421_miss_lemnas" = filterHogs(
        requiredAcc = setdiff(accSets$Lemnas, "Ljaponica9421"),
        prohibitedAcc = "Ljaponica9421",
        hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
      "Ljaponica9421_spec_lemnas" = filterHogs(
        requiredAcc = "Ljaponica9421",
        prohibitedAcc = setdiff(accSets$Lemnas, "Ljaponica9421"), 
        hogMatrix = getOfHogsMatrix(subgen = FALSE, node = "N4")),
      "Ljaponica_subgen_variable" = filterHogs(
        optionalAcc = c(accSetsDwSubgen$Ljaponicas_subM, accSetsDwSubgen$Ljaponicas_subT),
        minOptionalAcc = 1,
        maxOptionalAcc = length(c(accSetsDwSubgen$Ljaponicas_subM, accSetsDwSubgen$Ljaponicas_subT)) - 1,
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica7182_subM_miss_dip" = filterHogs(
        requiredAcc = accSetsDwSubgen$Lminors,
        prohibitedAcc = "Ljaponica7182_M",
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica7182_subM_spec_dip" = filterHogs(
        requiredAcc = "Ljaponica7182_M",
        prohibitedAcc = accSetsDwSubgen$Lminors,
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica7182_subT_miss_dip" = filterHogs(
        requiredAcc = accSetsDwSubgen$Lturioniferas,
        prohibitedAcc = "Ljaponica7182_T",
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica7182_subT_spec_dip" = filterHogs(
        requiredAcc = "Ljaponica7182_T",
        prohibitedAcc = accSetsDwSubgen$Lturioniferas,
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica8627_subM_miss_dip" = filterHogs(
        requiredAcc = accSetsDwSubgen$Lminors,
        prohibitedAcc = "Ljaponica8627_M",
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica8627_subM_spec_dip" = filterHogs(
        requiredAcc = "Ljaponica8627_M",
        prohibitedAcc = accSetsDwSubgen$Lminors,
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica8627_subT_miss_dip" = filterHogs(
        requiredAcc = accSetsDwSubgen$Lturioniferas,
        prohibitedAcc = "Ljaponica8627_T",
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica8627_subT_spec_dip" = filterHogs(
        requiredAcc = "Ljaponica8627_T",
        prohibitedAcc = accSetsDwSubgen$Lturioniferas,
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica9421_subM_miss_dip" = filterHogs(
        requiredAcc = accSetsDwSubgen$Lminors,
        prohibitedAcc = "Ljaponica9421_M",
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica9421_subM_spec_dip" = filterHogs(
        requiredAcc = "Ljaponica9421_M",
        prohibitedAcc = accSetsDwSubgen$Lminors,
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica9421_subT_miss_dip" = filterHogs(
        requiredAcc = accSetsDwSubgen$Lturioniferas,
        prohibitedAcc = "Ljaponica9421_T",
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponica9421_subT_spec_dip" = filterHogs(
        requiredAcc = "Ljaponica9421_T",
        prohibitedAcc = accSetsDwSubgen$Lturioniferas,
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponicas_subM_miss_dip" = filterHogs(
        requiredAcc = accSetsDwSubgen$Lminors,
        prohibitedAcc = accSetsDwSubgen$Ljaponicas_subM, 
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponicas_subM_spec_dip" = filterHogs(
        requiredAcc = accSetsDwSubgen$Ljaponicas_subM,
        prohibitedAcc = accSetsDwSubgen$Lminors,
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponicas_subM_miss_Lm9252" = filterHogs(
        requiredAcc = accSetsDwSubgen$Lminors,
        prohibitedAcc = accSetsDwSubgen$Ljaponicas_subM, 
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponicas_subM_spec_Lm9252" = filterHogs(
        requiredAcc = accSetsDwSubgen$Ljaponicas_subM,
        prohibitedAcc = accSetsDwSubgen$Lminors,
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponicas_subT_miss_dip" = filterHogs(
        requiredAcc = "Lturionifera9434",
        prohibitedAcc = accSetsDwSubgen$Ljaponicas_subT, 
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponicas_subT_spec_dip" = filterHogs(
        requiredAcc = accSetsDwSubgen$Ljaponicas_subT,
        prohibitedAcc = "Lturionifera9434",
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponicas_subT_miss_subM" = filterHogs(
        requiredAcc = accSetsDwSubgen$Ljaponicas_subM,
        prohibitedAcc = accSetsDwSubgen$Ljaponicas_subT,
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4")),
      "Ljaponicas_subM_miss_subT" = filterHogs(
        requiredAcc = accSetsDwSubgen$Ljaponicas_subT,
        prohibitedAcc = accSetsDwSubgen$Ljaponicas_subM,
        hogMatrix = getOfHogsMatrix(subgen = TRUE, node = "N4"))
    )
# Get the gene IDs for the exclusive HOGs
selectedHogsGeneIdsHybrids <- lapply(
  filteredHogSetsHybrids[1:7],
  genesFromHogs,
  accUniverse = accSets$Ljaponicas,
  hogTable = getOfHogs(subgen = FALSE, node = "N4")
) %>% append(lapply(
  filteredHogSetsHybrids[8:length(filteredHogSetsHybrids)],
  genesFromHogs,
  accUniverse = c(accSets$Ljaponicas, accSetsDwSubgen$Ljaponicas_subM, accSetsDwSubgen$Ljaponicas_subT),
  hogTable = getOfHogs(subgen = TRUE, node = "N4")))


# 
# HOG tables ----
# 
dir.create(hogTableOutDir, showWarnings = FALSE)

## Annotate HOGs ----
hogTablesAngio <- purrr::map(
  filteredHogSetsAngio, ~ getAnnotatedHogs(., node = "N1") %>% purgeHogs(geneIds = putativeDwTeAndOrganelleGenes, TRUE)
)
hogTablesAth1Mono <- purrr::map(
  filteredHogSetsAth1Mono, ~ getAnnotatedHogs(., node = "N4") %>% purgeHogs(geneIds = putativeDwTeAndOrganelleGenes, TRUE)
)
hogTablesAthOsa <- purrr::map(
  filteredHogSetsAthOsa, ~ getAnnotatedHogs(., node = "N4") %>% purgeHogs(geneIds = putativeDwTeAndOrganelleGenes, TRUE)
)
hogTablesLemnaceae <- purrr::map(
  filteredHogSetsLemnaceae, ~ getAnnotatedHogs(., node = "N4") %>% purgeHogs(geneIds = putativeDwTeAndOrganelleGenes, TRUE)
)
hogTablesHybrids <- purrr::map(
  filteredHogSetsHybrids, ~ getAnnotatedHogs(., node = "N4", subgen = TRUE) %>% purgeHogs(geneIds = putativeDwTeAndOrganelleGenes, TRUE)
)


#
# GO term enrichment ----
# (HOG Presence/Absence)

#' Compute GO Term enrichment for selected genes in a gene universe, reduce with rrvgo, and output resulting tables and
#' plots.
#'
#' @param selectedGeneIds vector of gene IDs of interest.
#' @param allGeneIds vector of IDs representing the gene universe to compute enrichment over.
#' @param gene2Go mappings of gene IDs to go terms in the format expected by TopGO
#' @param outDir directory to write results to. A "plots" subdir will be created beneath.
#' @param resName output filename stem.
#' @param ontologies vector of names of ontologies over which enrichment will be computed.
#' @param pValueCutoff minimum pval resulting from topGO test to retain as a significantly enriched term.
#' @return list of the number of significant terms discovered, indexed by ontology name.
goTermEnrichment <- function(
  selectedGeneIds,
  allGeneIds,
  gene2Go,
  outDir = "results",
  resName = "topgo",
  ontologies = c("BP", "MF", "CC"),
  pValueCutoff = 0.01
) {
  outDir <- paste0(hogTableOutDir, "/topGO/", outDir)
  print(paste0("Running TopGO and rrvgo on ", resName))
  print(paste0("Writing output to ", outDir))
  dir.create(paste0(outDir, "/plots"), showWarnings = FALSE, recursive = TRUE)
  
  geneList <- factor(as.integer(allGeneIds %in% selectedGeneIds))
  names(geneList) <- allGeneIds
  
  # Write empty sig GO files to cover all of the bailout cases that follow.
  header <- tibble::tribble(
    ~"GO.ID", ~"Term", ~"Annotated", ~"Significant", ~"Expected", ~"weight01.fisher", ~"genes"
  )
  headerReduced <- tibble::tribble(
    ~"GO.ID", ~"Term", ~"Annotated", ~"Significant", ~"Expected", ~"weight01.fisher", ~"parent", ~"parentSimScore", 
    ~"score", ~"parentTerm", ~"genes"
  )
  for (ontology in ontologies) {
    outPrefix = paste0(outDir, "/", resName, ".", ontology)
    readr::write_tsv(header, paste0(outPrefix, ".sig-gos.txt"))
    readr::write_tsv(headerReduced, paste0(outPrefix, ".sig-gos.reduced.txt"))
  }
  
  # Handle case where no geneIds are selected
  if (length(levels(geneList)) < 2) {
    return(NULL)
  }
  
  results <- list()
  
  for (ontology in ontologies) {
    print(paste0("Ontology: ", ontology))
    goData <- new(
      "topGOdata",
      description = resName,
      ontology = ontology,
      allGenes = geneList,
      annot = topGO::annFUN.gene2GO,
      gene2GO = gene2Go
    )
    
    resultFisher <- topGO::runTest(goData, algorithm = "weight01", statistic = "fisher")
    resultFisherSummary <- summary(attributes(resultFisher)$score <= pValueCutoff)
    print(resultFisherSummary)
    nSigTerms <- 0
    if (length(resultFisherSummary) == 3) {
      nSigTerms <- as.integer(resultFisherSummary[[3]])
    }
    
    if (!purrr::is_null(goData) && !purrr::is_null(resultFisher) && nSigTerms > 0) {
      # problem: weight01.fisher is returned from GenTable as chr type, rather than numeric, where the most
      # significant values are expressed as "< 1e-30" or "<1e-30", which breaks numeric casting. It should not be
      # necessary to filter for sig terms, as we take topNodes = nSigTerms in GenTable().
      tgTable <- topGO::GenTable(
        goData,
        weight01.fisher = resultFisher,
        orderBy = "weight01.fisher",
        ranksOf = "weight01.fisher",
        topNodes = nSigTerms,
        numChar = 1000
      ) %>%
        # Sometimes this column is misnamed as the function call apply(...)
        dplyr::rename_with(
          .cols = tidyr::starts_with("apply"), 
          ~ if (length(.x) > 0) { "weight01.fisher" } else { as.character() } ) %>% 
        dplyr::mutate(
          weight01.fisher = weight01.fisher %>%
            stringr::str_replace(pattern = "< *1e-30", replacement = "1e-30") %>%
            as.numeric()
        )
      
      results[[ontology]] <- nSigTerms
      
      # Add the selected genes in each significant GO term to the table
      sigTerms <- tgTable$GO.ID
      genesInTerms <- topGO::genesInTerm(goData, sigTerms)
      genesInTerms2 <- purrr::map(
        genesInTerms,
        ~ intersect(.x, selectedGeneIds) %>% paste(collapse = ",")
      )
      tgTable2 <- tgTable %>% left_join(
        tibble::tibble(GO.ID = names(genesInTerms2), genes = genesInTerms2) %>% tidyr::unnest(genes),
        by = "GO.ID"
      )
      
      outPrefix = paste0(outDir, "/", resName, ".", ontology)
      readr::write_tsv(tgTable2, paste0(outPrefix, ".sig-gos.txt"))
      
      # rrvgo plots
      if (! is.null(tgTable) && nrow(tgTable) > 0) {
        
        # Need to transform the pValues into a named vector of scores with increasing significance for rrvgo
        scores <- setNames(-log10(tgTable$weight01.fisher), tgTable$GO.ID)
        
        # Default the reduced term table to the full sig term table. Overwrite below if reduction is successful.
        reducedTerms <- tgTable2 %>%
          dplyr::rename("go" = "GO.ID", "term" = "Term") %>%
          dplyr::mutate(parent = go, parentSimScore = 1.000, parentTerm = term) %>% 
          dplyr::left_join(scores %>% tibble::enframe(name = "go", value = "score"), by = "go")

        simMatrix <- matrix(0, 0, 0)
        if (nrow(tgTable) > 1) {
          simMatrix <- rrvgo::calculateSimMatrix(
            tgTable$GO.ID,
            orgdb = "org.At.tair.db",
            ont = ontology, method = "Rel"
          )
          if(is.matrix(simMatrix) && length(simMatrix > 0)) {
            reducedTerms <- rrvgo::reduceSimMatrix(
              simMatrix = simMatrix, scores = scores, threshold = 0.7, orgdb="org.At.tair.db"
            )
          }
        }
        
        pdf(paste0(outDir, "/plots/", paste(sep = ".", resName, ontology, "treemap.pdf")),
            width = 7.2, height = 7.2)
        rrvgo::treemapPlot(
          reducedTerms,
          size = "score",
          force.print.labels = TRUE,
        )
        dev.off()
      }

      readr::write_tsv(
        tgTable2 %>% 
          dplyr::filter(GO.ID %in% reducedTerms$go) %>% 
          dplyr::left_join(
            reducedTerms %>% select(go, parent, parentSimScore, score, parentTerm), 
            by = c("GO.ID" = "go")
          ) %>% 
          dplyr::relocate(.after = parentTerm, genes),
        paste0(outPrefix, ".sig-gos.reduced.txt"))
    }
  }
  
  return(results)
}

#' Read topGO results from a txt file. This avoids recomputing enrichment, important because @function{goTermEnrichment} 
#' is so costly.
#'
#' @param inDir directory to write results to. A "plots" subdir will be created beneath.
#' @param resName output filename stem.
#' @param ontology one of "BP", "MF", or "CC". Will be appended the resName for reading.
#' @param reduced if TRUE, load the reduced (e.g. rrvgo) results instead of full TopGO results.
#' @return vector of gene IDs.
loadTopGoResults <- function(inDir, resName, ontology, reduced = FALSE) {
  if (reduced) {
    sigGoFile <- paste0(hogTableOutDir, "/topGO/", inDir, "/", resName, ".", ontology, ".sig-gos.reduced.txt")
    df <- readr::read_tsv(
      sigGoFile, 
      col_type = "ccddddcddcc", 
      col_select = c("GO.ID", "Term", "weight01.fisher", "parent", "parentSimScore", "score", "parentTerm"))
  } else {
    sigGoFile <- paste0(hogTableOutDir, "/topGO/", inDir, "/", resName, ".", ontology, ".sig-gos.txt")
    df <- readr::read_tsv(sigGoFile, col_type = "ccddddc", col_select = c("GO.ID", "Term", "weight01.fisher"))
  }

  return(df)
}


# This is large in memory: ~1.1G
gene2Go <- topGO::readMappings(paste0(ofResultsDir, "/eggnog/genes2go.reduced.txt"))

### 
### Run TopGO and rrvgo ----
### 
allGeneIdsAngio <- as.character(
  annots %>% dplyr::filter(
    stringr::str_detect(accession, pattern = paste(accSets$angiosperms, collapse = "|"))
  ) %>% 
    dplyr::pull(query_name)
)
allGeneIdsAth1Mono <- as.character(
  annots %>% dplyr::filter(
    stringr::str_detect(accession, paste(accUniverseAth1Mono, collapse = "|"))
  ) %>% 
    dplyr::pull(query_name)
)
allGeneIdsAthOsa <- as.character(
  annots %>% dplyr::filter(
    stringr::str_detect(accession, paste(accUniverseAthOsa, collapse = "|"))
  ) %>%
    dplyr::pull(query_name)
)
allGeneIdsLemnaceae <- as.character(
  annots %>% dplyr::filter(
    stringr::str_detect(accession, paste(c(accSets$Lemnaceae), collapse = "|"))
  ) %>% 
    dplyr::pull(query_name)
)
allGeneIdsHybrids <- as.character(
  annots %>% dplyr::filter(
    stringr::str_detect(accession, paste(c(accSets$Ljaponicas), collapse = "|"))
  ) %>% 
    dplyr::pull(query_name)
)


### N.B. - these goTermEnrichment() calls take a loooong time, so if they've already been computed, skip this and
### load below!!!
topGoResultsAngioNSigTerms <- purrr::imap(
  selectedHogsGeneIdsAngio,
  ~ {
    goTermEnrichment(
      selectedGeneIds = .x,
      allGeneIds = allGeneIdsAngio,
      resName = .y,
      gene2Go = gene2Go,
      outDir = "all_angiosperms",
      ontologies = c("BP", "MF", "CC"),
      pValueCutoff = 0.01
    )
  }
)
topGoResultsAth1MonoNSigTerms <- purrr::imap(
  selectedHogsGeneIdsAth1Mono,
  ~ {
    goTermEnrichment(
      selectedGeneIds = .x,
      allGeneIds = allGeneIdsAth1Mono,
      resName = .y,
      gene2Go = gene2Go,
      outDir = "ath_1mono_miss",
      ontologies = c("BP", "MF", "CC"),
      pValueCutoff = 0.01
    )
  }
)
topGoResultsAthOsaNSigTerms <- purrr::imap(
  selectedHogsGeneIdsAthOsa,
  ~ {
    goTermEnrichment(
      selectedGeneIds = .x,
      allGeneIds = allGeneIdsAthOsa,
      resName = .y,
      gene2Go = gene2Go,
      outDir = "ath_osa",
      ontologies = c("BP", "MF", "CC"),
      pValueCutoff = 0.01
    )
  }
)
topGoResultsLemnaceaeNSigTerms <- purrr::imap(
  selectedHogsGeneIdsLemnaceae,
  ~ {
    goTermEnrichment(
      selectedGeneIds = .x,
      allGeneIds = allGeneIdsLemnaceae,
      resName = .y,
      gene2Go = gene2Go,
      outDir = "intra_lemnaceae",
      ontologies = c("BP", "MF", "CC"),
      pValueCutoff = 0.01
    )
  }
)
topGoResultsHybridsNSigTerms <- purrr::imap(
  selectedHogsGeneIdsHybrids,
  ~ {
    goTermEnrichment(
      selectedGeneIds = .x,
      allGeneIds = allGeneIdsHybrids,
      resName = .y,
      gene2Go = gene2Go,
      outDir = "hybrids",
      ontologies = c("BP", "MF", "CC"),
      pValueCutoff = 0.01
    )
  }
)


## 
## Read GO enrichment results ----
## 
topGoResultsAngio <- list()
for(ontology in c("BP", "MF", "CC")) {
  topGoResultsAngio[[ontology]] <- as.list(names(selectedHogsGeneIdsAngio))
  names(topGoResultsAngio[[ontology]]) <- names(selectedHogsGeneIdsAngio)
  topGoResultsAngio[[ontology]] <- purrr::map(
    topGoResultsAngio[[ontology]],
    loadTopGoResults,
    inDir = "all_angiosperms",
    ontology = ontology
  )
}
topGoResultsAth1Mono <- list()
for(ontology in c("BP", "MF", "CC")) {
  topGoResultsAth1Mono[[ontology]] <- as.list(names(selectedHogsGeneIdsAth1Mono))
  names(topGoResultsAth1Mono[[ontology]]) <- names(selectedHogsGeneIdsAth1Mono)
  topGoResultsAth1Mono[[ontology]] <- purrr::map(
    topGoResultsAth1Mono[[ontology]],
    loadTopGoResults,
    inDir = "ath_1mono_miss",
    ontology = ontology
  )
}
topGoResultsAthOsa <- list()
for(ontology in c("BP", "MF", "CC")) {
  topGoResultsAthOsa[[ontology]] <- as.list(names(selectedHogsGeneIdsAthOsa))
  names(topGoResultsAthOsa[[ontology]]) <- names(selectedHogsGeneIdsAthOsa)
  topGoResultsAthOsa[[ontology]] <- purrr::map(
    topGoResultsAthOsa[[ontology]],
    loadTopGoResults,
    inDir = "ath_osa",
    ontology = ontology
  )
}
topGoResultsLemnaceae <- list()
for(ontology in c("BP", "MF", "CC")) {
  topGoResultsLemnaceae[[ontology]] <- as.list(names(selectedHogsGeneIdsLemnaceae))
  names(topGoResultsLemnaceae[[ontology]]) <- names(selectedHogsGeneIdsLemnaceae)
  topGoResultsLemnaceae[[ontology]] <- purrr::map(
    topGoResultsLemnaceae[[ontology]],
    loadTopGoResults,
    inDir = "intra_lemnaceae",
    ontology = ontology
  )
}
topGoResultsHybrids <- list()
for(ontology in c("BP", "MF", "CC")) {
  topGoResultsHybrids[[ontology]] <- as.list(names(selectedHogsGeneIdsHybrids))
  names(topGoResultsHybrids[[ontology]]) <- names(selectedHogsGeneIdsHybrids)
  topGoResultsHybrids[[ontology]] <- purrr::map(
    topGoResultsHybrids[[ontology]],
    loadTopGoResults,
    inDir = "hybrids",
    ontology = ontology
  )
}

# Reduced GO terms
topGoResultsReducedAngio <- list()
for(ontology in c("BP", "MF", "CC")) {
  topGoResultsReducedAngio[[ontology]] <- as.list(names(selectedHogsGeneIdsAngio))
  names(topGoResultsReducedAngio[[ontology]]) <- names(selectedHogsGeneIdsAngio)
  topGoResultsReducedAngio[[ontology]] <- purrr::map(
    topGoResultsReducedAngio[[ontology]],
    loadTopGoResults,
    inDir = "all_angiosperms",
    ontology = ontology,
    reduced = TRUE
  )
}
topGoResultsReducedAth1Mono <- list()
for(ontology in c("BP", "MF", "CC")) {
  topGoResultsReducedAth1Mono[[ontology]] <- as.list(names(selectedHogsGeneIdsAth1Mono))
  names(topGoResultsReducedAth1Mono[[ontology]]) <- names(selectedHogsGeneIdsAth1Mono)
  topGoResultsReducedAth1Mono[[ontology]] <- purrr::map(
    topGoResultsReducedAth1Mono[[ontology]],
    loadTopGoResults,
    inDir = "ath_1mono_miss",
    ontology = ontology,
    reduced = TRUE
  )
}
topGoResultsReducedAthOsa <- list()
for(ontology in c("BP", "MF", "CC")) {
  topGoResultsReducedAthOsa[[ontology]] <- as.list(names(selectedHogsGeneIdsAthOsa))
  names(topGoResultsReducedAthOsa[[ontology]]) <- names(selectedHogsGeneIdsAthOsa)
  topGoResultsReducedAthOsa[[ontology]] <- purrr::map(
    topGoResultsReducedAthOsa[[ontology]],
    loadTopGoResults,
    inDir = "ath_osa",
    ontology = ontology,
    reduced = TRUE
  )
}
topGoResultsReducedLemnaceae <- list()
for(ontology in c("BP", "MF", "CC")) {
  topGoResultsReducedLemnaceae[[ontology]] <- as.list(names(selectedHogsGeneIdsLemnaceae))
  names(topGoResultsReducedLemnaceae[[ontology]]) <- names(selectedHogsGeneIdsLemnaceae)
  topGoResultsReducedLemnaceae[[ontology]] <- purrr::map(
    topGoResultsReducedLemnaceae[[ontology]],
    loadTopGoResults,
    inDir = "intra_lemnaceae",
    ontology = ontology,
    reduced = TRUE
  )
}
topGoResultsReducedHybrids <- list()
for(ontology in c("BP", "MF", "CC")) {
  topGoResultsReducedHybrids[[ontology]] <- as.list(names(selectedHogsGeneIdsHybrids))
  names(topGoResultsReducedHybrids[[ontology]]) <- names(selectedHogsGeneIdsHybrids)
  topGoResultsReducedHybrids[[ontology]] <- purrr::map(
    topGoResultsReducedHybrids[[ontology]],
    loadTopGoResults,
    inDir = "hybrids",
    ontology = ontology,
    reduced = TRUE
  )
}


hogTablesAngio <- purrr::pmap(
  list(
    annHogTbl = hogTablesAngio, 
    goResultsBp = topGoResultsAngio$BP,
    goResultsMf = topGoResultsAngio$MF,
    goResultsCc = topGoResultsAngio$CC
  ), 
  addGoResultsToHogTbl
)
hogTablesAth1Mono <- purrr::pmap(
  list(
    annHogTbl = hogTablesAth1Mono, 
    goResultsBp = topGoResultsAth1Mono$BP,
    goResultsMf = topGoResultsAth1Mono$MF,
    goResultsCc = topGoResultsAth1Mono$CC
  ), 
  addGoResultsToHogTbl
)
hogTablesAthOsa <- purrr::pmap(
  list(
    annHogTbl = hogTablesAthOsa,
    goResultsBp = topGoResultsAthOsa$BP,
    goResultsMf = topGoResultsAthOsa$MF,
    goResultsCc = topGoResultsAthOsa$CC
  ),
  addGoResultsToHogTbl
)
hogTablesLemnaceae <- purrr::pmap(
  list(
    annHogTbl = hogTablesLemnaceae, 
    goResultsBp = topGoResultsLemnaceae$BP,
    goResultsMf = topGoResultsLemnaceae$MF,
    goResultsCc = topGoResultsLemnaceae$CC
  ), 
  addGoResultsToHogTbl
)
hogTablesHybrids <- purrr::pmap(
  list(
    annHogTbl = hogTablesHybrids, 
    goResultsBp = topGoResultsHybrids$BP,
    goResultsMf = topGoResultsHybrids$MF,
    goResultsCc = topGoResultsHybrids$CC
  ), 
  addGoResultsToHogTbl
)



## 
## Write tables ----
## 
writeHogTables(
  hogTables = list(
    "all_angiosperms" = hogTablesAngio %>% dplyr::bind_rows(.id = "contrast"),
    "ath_1mono" = hogTablesAth1Mono %>% dplyr::bind_rows(.id = "contrast"),
    "ath_osa" = hogTablesAthOsa %>% dplyr::bind_rows(.id = "contrast"),
    "intra_lemnaceae" = hogTablesLemnaceae %>% dplyr::bind_rows(.id = "contrast"),
    "hybrids" = hogTablesHybrids %>% dplyr::bind_rows(.id = "contrast")
  ), 
  outDir = hogTableOutDir, 
  tblName = "hog-table", 
  tblNameAppendDate = TRUE
)


##
## Plot GO term enrichment summary ----
##
topTermsAthOsaBp <- topGoResultsReducedAthOsa$BP %>% 
  dplyr::bind_rows(.id = "contrast") %>% 
  dplyr::mutate(
    contrast_type = dplyr::case_when(
      contrast %>% stringr::str_detect(pattern = "_miss_") ~ "missing HOGs",
      contrast %>% stringr::str_detect(pattern = "_spec_") ~ "unique paralogs"
    ),
    contrast = factor(contrast),
    Term = factor(Term)
  ) %>% 
  dplyr::filter(
    score >= 3,
    contrast %in% c("Lemnaceae_miss_ath_osa", "Lemnaceae_spec_ath_osa", "cdemersum_dw_miss_ath_osa", "cdemersum_dw_spec_ath_osa")
  ) %>%
  droplevels() %>% 
  dplyr::group_by(contrast) %>% 
  dplyr::ungroup() %>%
  droplevels() %>% 
  dplyr::group_by(contrast, parent) %>% 
  dplyr::add_count(name = "group_size")


treemapColors <- c(
  "#41575d",
  "#0072a8",
  "#122129",
  "#367dab",
  "#0c202e",
  "#1e7fa5",
  "#032437",
  "#26809b",
  "#002d3f",
  "#5e7a82",
  "#004c6f",
  "#64717e",
  "#1c3d45",
  "#517d8f",
  "#184461",
  "#5a758e",
  "#004959",
  "#21638d",
  "#566670",
  "#01708e",
  "#3f5568",
  "#006b7c",
  "#1d5561",
  "#3f6b76",
  "#015b72",
  "#2e5e69"
)
treemapColorsLight <- treemapColors %>% colorspace::lighten(0.55)
treemapColorsDark <- treemapColors %>% colorspace::darken(0.40)

plotTopTermsAthOsaTreeMap <- topTermsAthOsaBp %>% 
  filter(contrast == "Lemnaceae_miss_ath_osa") %>% 
  dplyr::mutate(Term = dplyr::case_when(group_size == 1 ~ "", .default = Term)) %>% 
  ggplot(aes(area = score, label = Term, fill = parentTerm, subgroup = parentTerm)) + 
  geom_treemap(
    aes(color = parentTerm),
    layout = "squarified",
    start = "topleft",
  ) +
  geom_treemap_subgroup_border(
    layout = "squarified",
    start = "topleft",
    color = "#00384c",
    size = 2
  ) +
  geom_treemap_text(
    aes(color = parentTerm),
    alpha = 0.9,
    layout = "squarified",
    start = "topleft",
    family = "Helvetica",
    fontface = "plain", 
    place = "bottomleft", 
    grow = FALSE, 
    reflow = TRUE,
    min.size = 0,
    size = 5
  ) +
  geom_treemap_subgroup_text(
    alpha = 1, 
    layout = "squarified",
    start = "topleft",
    family = "Helvetica",
    fontface = "bold", 
    colour = "#dedede", 
    place = "topleft", 
    grow = FALSE, 
    reflow = TRUE,
    min.size = 0,
    size = 7
  ) +
  scale_fill_manual(values = treemapColors) +
  scale_color_manual(values = treemapColorsLight) +
  theme(
    legend.position = "none"
  )
pdf("plots/topgo-overview.treemap.pdf", width = 3.25, height = 3.25)
plotTopTermsAthOsaTreeMap
dev.off()



plotTopTermsAthOsaTiles <- topTermsAthOsaBp %>% 
  ggplot() +
  theme_bw(base_size = 5) +
  facet_grid(cols = vars(contrast_type), scales = "free") +
  geom_tile(aes(x = factor(contrast), y = forcats::fct_rev(Term), fill = negLog10Pval), width = 1, height = 1) +
  scale_fill_gradient(low = "white", high = "blue", limits = c(0, NA))
pdf("plots/topgo-overview.tiles.pdf", width = 3.5, height = 5)
plotTopTermsAthOsaTiles
dev.off()


#
# Genes/OG stats overview ----
#
geneSummaryCols <- c(
  "unassigned.genes", 
  "genes.in.species.specific.orthogroups", 
  "genes.in.common.orthogroups"
)
geneSummary <- accTblSubgenSummary2 %>%
  dplyr::select(accession, genes, all_of(geneSummaryCols)) %>%
  tidyr::pivot_longer(
    names_to = "og_assignment", 
    values_to = "gene_count", 
    cols = all_of(geneSummaryCols)
  ) %>%
  dplyr::mutate(og_assignment = forcats::fct_rev(og_assignment))

plotOgOverview <- geneSummary %>%
  ggplot(aes(
    y = fct_relevel(accession, accAll),
    alpha = og_assignment,
    fill = accession
  )) +
  geom_col(aes(x = gene_count), position = "fill", width = 0.85) +
  scale_fill_manual(
    values = pal4, 
    na.value = pal4["all"] %>% colorspace::lighten(0.4),
    guide = guide_none()
  ) +
  scale_alpha_manual(
    values = c(0.33, 0.66, 1),
    labels = c("unassigned", "accession-specific", "multi-accession"),
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, .05)), 
    labels = scales::label_percent(accuracy = 1), 
    breaks = seq(0, 1, .1)
  ) +
  scale_y_discrete(expand = expansion(mult = c(.04, .025)), labels = sciNameLabeller) +
  coord_capped_cart(bottom = "both", xlim = c(0, 1)) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  ) +
  labs(x = "Gene assignment to orthogroups", y = NULL)
pdf("plots/orthofinder_og_overview.pdf", width = 3.6, height = 5)
plotOgOverview
dev.off()



#
# EggNog annotation plots ----
#
plotPredsVsEggNogAnnsCol <- accTblSubgenSummary2 %>%
  ggplot(aes(
    y = fct_reorder(accession, annots_per_seq),
    fill = fct_reorder(accession, annots_per_seq)
  )) +
  geom_col(aes(x = num_seqs), alpha = 0.5) +
  geom_col(aes(x = num_eggnog_annots)) +
  geom_text(aes(x = num_seqs, label = scales::percent(annots_per_seq, accuracy = 1)),
            nudge_x = 2500, size = 2.25, family = "Helvetica"
  ) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = pal4, na.value = "lightgrey") +
  scale_x_continuous(
    expand = expansion(mult = c(0, .05)), 
    labels = scales::comma, 
    breaks = seq(0, 50000, 10000)
  ) +
  scale_y_discrete(expand = expansion(mult = c(.04, .025)), labels = sciNameLabeller) +
  coord_capped_cart(bottom = "both", xlim = c(0, 50000)) +
  labs(x = "EggNog annotated / total sequences", y = NULL)
pdf("plots/proteome_summary_preds_vs_eggnog_anns.column.pdf", width = 3.6, height = 4)
plotPredsVsEggNogAnnsCol
dev.off()


plotPredsVsEggNogAnns <- accTblSubgenSummary2 %>%
  ggplot(aes(
    y = annots_per_seq,
    x = num_seqs,
    color = fct_reorder(accession, annots_per_seq)
  )) +
  geom_smooth(
    formula = y ~ x, 
    method = "lm", 
    color = "#4d4d4d", 
    fill = "#e9e9e9", 
    size = 0.3, 
    linetype = 3
  ) +
  geom_point() +
  geom_text_repel(aes(label = accession), force = 5, size = 2.25, family = "Helvetica", parse = TRUE) +
  #scale_color_viridis(discrete = TRUE, option = "magma", end = 0.8) +
  scale_color_manual(values = pal, na.value = palLight["all"]) +
  scale_x_continuous(
    expand = expansion(mult = c(0.025, 0.05)),
    labels = scales::comma, 
    breaks = seq(15000, 45000, 5000)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.025, 0.05)), 
    labels = scales::label_percent(accuracy = 1), 
    breaks = seq(0.7, 1, 0.10)
  ) +
  coord_capped_cart(bottom = "both", left = "both", ylim = c(.7, 1), xlim = c(15000, 45000)) +
  labs(
    title = NULL,
    y = "Functionally annotated by EggNog",
    x = "Predicted proteins"
  ) +
  theme(
    legend.position = "none"
  )
pdf("plots/proteome_summary_preds_vs_eggnog_anns.pdf", width = 5, height = 3)
plotPredsVsEggNogAnns
dev.off()


annots %>%
  ggplot(
    aes(
      y = eN_seed_olog_score, 
      x = fct_reorder(accession, eN_seed_olog_score), 
      fill = fct_reorder(accession, eN_seed_olog_score)
    )
  ) +
  geom_violin() +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "none"
  ) +
  coord_capped_cart(bottom = "both", left = "both", ylim = c(0, 2500)) +
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  ggtitle("Seed ortholog scores by accession") +
  labs(x = NULL, y = "Seed ortholog score")


annots %>%
  dplyr::group_by(best_tax_level) %>%
  dplyr::tally(sort = TRUE) %>%
  head(n = 10) %>%
  ggplot(aes(y = n, x = fct_reorder(best_tax_level, n), fill = best_tax_level)) +
  geom_col() +
  theme_pubclean() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    legend.position = "none"
  ) +
  scale_fill_viridis(discrete = TRUE, option = "magma") +
  ggtitle("Annotation counts by best taxon level")


# Upset plots ----
pacman::p_load_current_gh("const-ae/ggupset@master")

hogSetsColors <- c(
  "Spolyrhizas_spec" = palsDfV4 %>% filter(acc == "Spolyrhiza9509") %>% pull(base_color),
  "Spolyrhizas_miss" = palsDfV4 %>% filter(acc == "Spolyrhiza9509") %>% pull(light30),
  "Wolffias_spec" = palsDfV4 %>% filter(acc == "Waustraliana8730") %>% pull(base_color),
  "Wolffias_miss" = palsDfV4 %>% filter(acc == "Waustraliana8730") %>% pull(light30),
  "Waustralianas_spec" = palsDfV4 %>% filter(acc == "Waustraliana8730") %>% pull(base_color),
  "Waustralianas_miss" = palsDfV4 %>% filter(acc == "Waustraliana8730") %>% pull(light30),
  "Lgibbas_spec" = palsDfV4 %>% filter(acc == "Lgibba7742a") %>% pull(base_color),
  "Lgibbas_miss" = palsDfV4 %>% filter(acc == "Lgibba7742a") %>% pull(light30),
  "Lminors_spec" = palsDfV4 %>% filter(acc == "Lminor7210") %>% pull(base_color),
  "Lminors_miss" = palsDfV4 %>% filter(acc == "Lminor7210") %>% pull(light30),
  "Ljaponicas_spec" = palsDfV4 %>% filter(acc == "Ljaponica7182") %>% pull(base_color),
  "Ljaponicas_miss" = palsDfV4 %>% filter(acc == "Ljaponica7182") %>% pull(light30),
  "Lturioniferas_spec" = palsDfV4 %>% filter(acc == "Lturionifera9434") %>% pull(base_color),
  "Lturioniferas_miss" = palsDfV4 %>% filter(acc == "Lturionifera9434") %>% pull(light30),
  "Lemnaceae_spec" = palsDfV4 %>% filter(acc == "all") %>% pull(base_color),
  "Lemnaceae_miss" = palsDfV4 %>% filter(acc == "all") %>% pull(light30),
  "Spirodelas_spec" = palsDfV4 %>% filter(acc == "Spolyrhiza9509") %>% pull(base_color),
  "Spirodelas_miss" = palsDfV4 %>% filter(acc == "Spolyrhiza9509") %>% pull(light30),
  "Lemnas_spec" = palsDfV4 %>% filter(acc == "all") %>% pull(base_color),
  "Lemnas_miss" = palsDfV4 %>% filter(acc == "all") %>% pull(light30),
  "aquatic_submerged_dw_spec" = palsDfV4 %>% filter(acc == "all") %>% pull(base_color),
  "aquatic_submerged_dw_miss" = palsDfV4 %>% filter(acc == "all") %>% pull(light30),
  "aquatic_rootless_spec" = palsDfV4 %>% filter(acc == "all") %>% pull(base_color),
  "aquatic_rootless_miss" = palsDfV4 %>% filter(acc == "all") %>% pull(light30),
  "cdemersum_dw_spec" = palsDfV4 %>% filter(acc == "all") %>% pull(base_color),
  "cdemersum_dw_miss" = palsDfV4 %>% filter(acc == "all") %>% pull(light30),
  "zmarina_dw_spec" = palsDfV4 %>% filter(acc == "all") %>% pull(base_color),
  "zmarina_dw_miss" = palsDfV4 %>% filter(acc == "all") %>% pull(light30),
  "monocots_spec" = palsDfV4 %>% filter(acc == "all") %>% pull(base_color),
  "monocots_miss" = palsDfV4 %>% filter(acc == "all") %>% pull(light30),
  "dicots_spec" = palsDfV4 %>% filter(acc == "all") %>% pull(base_color),
  "dicots_miss" = palsDfV4 %>% filter(acc == "all") %>% pull(light30),
  "turion_competent_spec" = palsDfV4 %>% filter(acc == "all") %>% pull(base_color),
  "turion_competent_miss" = palsDfV4 %>% filter(acc == "all") %>% pull(light30),
  "Spoly_and_Lturionifera_spec" = palsDfV4 %>% filter(acc == "all") %>% pull(base_color),
  "Spoly_and_Lturionifera_miss" = palsDfV4 %>% filter(acc == "all") %>% pull(light30),
  "sinking_dw_spec" = palsDfV4 %>% filter(acc == "all") %>% pull(base_color),
  "sinking_dw_miss" = palsDfV4 %>% filter(acc == "all") %>% pull(light30)
)


## Angiosperms upset ----

specificAndMissingSetsAngioMerged <- specificAndMissingSetsAngio %>% 
  purrr::map_chr(.f = ~ stringr::str_flatten(.x %>% naturalsort::naturalsort(), collapse = "-"))

specificAndMissingSetsAngioMergedDf <- specificAndMissingSetsAngioMerged %>% 
  dplyr::as_tibble(rownames = "intersection_name") %>% dplyr::rename("intersection" = "value") %>% 
  filter(! intersection_name %>% 
           grepl(pattern = "^(Spirodelas|Wolffias|aquatic_sinking|freshwater|aquatic_SF_|aquatic_floating)")
         ) %>% 
  droplevels()

hogSetsAngioColors <- hogSetsColors %>% as_tibble(rownames = "intersection_name") %>% 
  left_join(specificAndMissingSetsAngioMergedDf)

hogSetsAngioColorTbl <- hogSetsAngioColors$value %>% as.character()
names(hogSetsAngioColorTbl) <- hogSetsAngioColors$intersection

hogSetsAngioCounts <- hogTablesAngio %>% 
  purrr::map_dfr(.id = "intersection_name", .f = ~ .x %>% dplyr::summarize(nHOGs = n())) %>% 
  dplyr::left_join(specificAndMissingSetsAngioMergedDf) %>% 
  drop_na(intersection) %>% 
  mutate(intersection = factor(intersection, levels = specificAndMissingSetsAngioMergedDf$intersection %>% unique())) %>% 
  mutate(intersection_name = intersection_name %>% forcats::fct_relevel(hogSetsAngioColors$intersection_name)) %>%
  arrange(intersection_name) %>% 
  mutate(intersection = intersection %>% fct_inorder()) %>% 
  droplevels()

plotUpsetAngio <- hogSetsAngioCounts %>% 
  ggplot(aes(x = intersection, y = nHOGs)) +
  geom_col(aes(fill = intersection)) +
  ggupset::axis_combmatrix(
    sep = "-", 
    levels = intersect(accAngio, specificAndMissingSetsAngio %>% purrr::reduce(c) %>% unique()), 
    override_plotting_function = function(df) {
      df <- df %>%
        dplyr::mutate(labels_dots = dplyr::case_when(
          ! observed ~ "not observed",
          observed ~ as.character(labels)
        ))
      print("AFTER:")
      print(class(df))
      print(df, n = 100)
      df %>%
        ggplot(aes(x = at, y = single_label)) +
        geom_rect(aes(fill = index %% 2 == 0), ymin = df$index-0.5, ymax = df$index+0.5, xmin = 0, xmax = 1) +
        geom_point(aes(color = labels_dots), size = 1) +
        geom_line(data = function(dat) dat[dat$observed, , drop = FALSE], aes(group = labels, color = labels_dots), size = 1) +
        ylab("") + xlab("") +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        scale_y_discrete(labels = dwHighlighter) +
        scale_fill_manual(values = c(`TRUE` = "white", `FALSE` = "#F7F7F7")) +
        scale_color_manual(values = c("not observed" = alpha("white", 0), hogSetsAngioColorTbl)) +
        guides(fill = "none") +
        theme(
          legend.position = "none",
          panel.background = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_markdown(size = 5),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length = unit(0, "pt"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank()
        )
    }) +
  xlab("") +
  ylab("exclusive HOGs") +
  scale_fill_manual(values = hogSetsAngioColorTbl) +
  scale_color_manual(values = colorspace::darken(hogSetsAngioColorTbl, .30)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 215)) +
  guides(fill = "none") +
  theme(
    plot.margin = margin(2, 1, 0, 0),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_markdown(size = 5),
    axis.ticks.x = element_blank(),
    axis.line.x.bottom = element_blank(),
    panel.grid.major.y = element_line(size = 0.25),
    combmatrix.label.extra_spacing = 2
  )
pdf("plots/upset-angio.v5.pdf", width = 3.4, height = 4)
plotUpsetAngio
dev.off()



## Intra-Lemnaceae upset ----
specificAndMissingSetsLemnaceaeMerged <- specificAndMissingSetsLemnaceae %>% 
  purrr::map_chr(.f = ~ stringr::str_flatten(.x %>% naturalsort::naturalsort(), collapse = "-"))

specificAndMissingSetsLemnaceaeMergedDf <- specificAndMissingSetsLemnaceaeMerged %>% 
  dplyr::as_tibble(rownames = "intersection_name") %>% dplyr::rename("intersection" = "value") %>% 
  filter(! intersection_name %>% grepl(pattern = "^(Spirodelas|Wolffias)")) %>% droplevels()

hogSetsLemnaceaeColors <- hogSetsColors %>% as_tibble(rownames = "intersection_name") %>% 
  left_join(specificAndMissingSetsLemnaceaeMergedDf)

hogSetsLemnaceaeColorTbl <- hogSetsLemnaceaeColors$value %>% as.character()
names(hogSetsLemnaceaeColorTbl) <- hogSetsLemnaceaeColors$intersection

hogSetsLemnaceaeCounts <- hogTablesLemnaceae %>% 
  purrr::map_dfr(.id = "intersection_name", .f = ~ .x %>% dplyr::summarize(nHOGs = n())) %>% 
  dplyr::left_join(specificAndMissingSetsLemnaceaeMergedDf) %>% 
  drop_na(intersection) %>% 
  mutate(intersection = factor(intersection, levels = specificAndMissingSetsLemnaceaeMergedDf$intersection %>% unique())) %>% 
  mutate(intersection_name = intersection_name %>% forcats::fct_relevel(hogSetsAngioColors$intersection_name)) %>%
  arrange(intersection_name) %>% 
  mutate(intersection = intersection %>% fct_inorder()) %>% 
  droplevels()

plotUpsetLemnaceae <- hogSetsLemnaceaeCounts %>% 
  mutate(intersection_name = intersection_name %>% fct_relevel(hogSetsLemnaceaeColors$intersection_name)) %>%
  arrange(intersection_name) %>% 
  mutate(intersection = intersection %>% fct_inorder()) %>% 
  droplevels() %>% 
  ggplot(aes(x = intersection, y = nHOGs)) +
  geom_col(aes(fill = intersection)) +
  axis_combmatrix(
    sep = "-", 
    levels = intersect(accAngio, specificAndMissingSetsLemnaceae %>% reduce(c) %>% unique()), 
    override_plotting_function = function(df) {
      df <- df %>%
        droplevels() %>% 
        dplyr::mutate(labels_dots = dplyr::case_when(
          ! observed ~ "not observed",
          observed ~ as.character(labels)
        ))
      print(df$single_label %>% levels())
      print(df %>% filter(single_label == "Bdistachyon"))
      df %>%
        ggplot(aes(x = at, y = single_label)) +
        geom_rect(aes(fill = index %% 2 == 0), ymin=df$index-0.5, ymax=df$index+0.5, xmin=0, xmax=1) +
        geom_point(aes(color = labels_dots), size = 2) +
        geom_line(data = function(dat) dat[dat$observed, , drop = FALSE], aes(group = labels, color = labels_dots), size= 1.2) +
        ylab("") + xlab("") +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        scale_y_discrete(labels = dwHighlighter) +
        scale_fill_manual(values = c(`TRUE` = "white", `FALSE` = "#F7F7F7")) +
        scale_color_manual(values = c("not observed" = alpha("white", 0), hogSetsLemnaceaeColorTbl)) +
        guides(fill = "none") +
        theme(
          legend.position = "none",
          panel.background = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_markdown(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length = unit(0, "pt"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank()
        )
    }) +
  xlab("") +
  ylab("HOGs") +
  scale_fill_manual(values = hogSetsLemnaceaeColorTbl) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = "none") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x.bottom = element_blank(),
    panel.grid.major.y = element_line(size = 0.25)
  )
pdf("plots/upset-intra-lemnaceae.pdf", width = 3.75, height = 2.9)
plotUpsetLemnaceae
dev.off()


# 
# HOG copy number heatmaps ----
# 
addHogTablePlotLabels <- function(hogTable, boldStrings = "a^") {
  # Note: the default value for boldStrings is a regex that will never match.
  hogTable <- hogTable %>% 
    dplyr::mutate(
      plotLabel = dplyr::case_when(
        At_syms != ""        ~ At_syms %>%
          stringr::str_split(pattern = ",") %>% 
          purrr::map(
            .f = ~ .x %>%
              # Keep only the last synonymous symbol
              stringr::str_replace(pattern = "^([^/]+/)+([^/]+)$", replacement = "\\2") %>%
              # Remove any "AT" or "At" prefixes from gene symbols, but not ATG IDs
              stringr::str_replace(pattern = "^(At|AT|at|aT)([^,/]{3,6})$", replacement = "\\2") %>%
              # Remove any dashes between the symbol and a trailing number
              stringr::str_replace(pattern = "^([^,/]{3,6})(−|-)([0-9]+)$", replacement = "\\1\\3") %>% 
              stringr::str_trim(side = "both") %>% 
              # Replace specific Ath symbols
              stringr::str_replace(pattern = "^ICU9$", replacement = "AGO1") %>% 
              stringr::str_replace(pattern = "^ZIP$", replacement = "AGO7") %>% 
              stringr::str_replace(pattern = "^OCP11$", replacement = "AGO4") %>% 
              stringr::str_replace(pattern = "^ZLL$", replacement = "AGO10") %>% 
              # Hack to sort ATG IDs to the back
              stringr::str_replace(pattern = "^(AT.G[0-9]+)$", replacement = "ZZZ\\1") %>% 
              naturalsort::naturalsort() %>%
              stringr::str_replace(pattern = "^ZZZ(AT.G[0-9]+)$", replacement = "\\1") %>% 
              stringr::str_replace(
                pattern = paste0("(", paste0(boldStrings, collapse = "|"), ")"), 
                replacement = "<strong>\\1</strong>"
              ) %>% 
              stringr::str_flatten(collapse = ", ")
          ) %>% 
          unlist(),
        
        Os_syms != ""        ~ Os_syms %>% 
          stringr::str_split(pattern = ",") %>% 
          purrr::map(
            .f = ~ .x %>% 
              # Keep only the last synonymous symbol
              stringr::str_replace(pattern = "^([^/]+/)+([^/]+)$", replacement = "\\2") %>%
              # Remove any "Os" prefixes from gene symbols
              stringr::str_replace(pattern = "/?(Os|OS)[^/,]{3,}/?", replacement = "") %>% 
              # Add the "Os" prefix to all symbols to distinguish from Athaliana symbols
              stringr::str_replace(pattern = "^(.+)$", replacement = "Os\\1") %>%
              stringr::str_trim(side = "both") %>% 
              naturalsort::naturalsort() %>% 
              stringr::str_flatten(collapse = ", ")
          ) %>%
          unlist(),
        
        # If the HOG contains neither Athaliana nor Osativa genes, use the symbols in the mapped EggNog OG
        EN_OGs_At_syms != "" ~ EN_OGs_At_syms,
        TRUE                   ~ ""
      )
    )
  return(hogTable)
}


# ComplexHeatmap setup
set.seed(123)
heatmapColFun <- function(x) {
  pal <- scales::brewer_pal(palette = 7, type = "div", direction = -1)(4)
  colx <- if_else(x <= 3, pal[x + 1], pal[length(pal)])
  return(colx)
}
heatmapColFun2 <- function(x) {
  pal <- scales::brewer_pal(palette = 7, type = "div", direction = -1)(4)
  colx <- if_else(x <= 3, pal[x + 1], pal[length(pal)])
  return(colx %>% alpha(alpha = 0.7))
}
heatmapColFun3 <- function(x) {
  pal <- scales::brewer_pal(palette = "RdGy", type = "seq", direction = -1)(4)
  colx <- if_else(x <= 3, pal[x + 1], pal[length(pal)])
  return(colx)
}
at <- seq(0, 3, 1)
heatmapLegend <- function(colFun, title = NULL) {
  lgd <- ComplexHeatmap::Legend(
    at = rev(at),
    title = title, 
    title_gp = gpar(fontsize = 5),
    title_position = "topcenter",
    legend_gp = gpar(fill = rev(colFun(at))),
    labels = rev(c("0", "1", "2", "3+")),
    labels_gp = gpar(fontsize = 6),
    grid_width = unit(3, "mm"),
    grid_height = unit(3, "mm"),
    direction = "horizontal"
  )
  return(lgd)
}


## 
## Developmental genes ----
## 
devHogsIds <- tibble::tribble(
  ~HOG,            ~plotLabel,
  "N4.HOG0012277", "OsZFP",
  "N4.HOG0002613", "NAL2,3 (WOX3)",
  "N4.HOG0014003", "ORC3",
  "N4.HOG0001584", "RL9/SLL1",
  "N4.HOG0002331", "OsSNDP1,2,3",
  "N4.HOG0017274", "EXPA7,18",
  "N4.HOG0000031", "MYB53,92,<strong>93</strong>",
  "N4.HOG0002304", "CMI1, PBP1",
  "N4.HOG0006359", "SOC1, <strong>XAL2</strong>, FYF, FYL1,2, AGL19",
  "N4.HOG0016364", "WOX5",
  "N4.HOG0002248", "RHD6, RSL1",
  "N4.HOG0002069", "RGI1,2",
  "N4.HOG0001615", "DOT5/WIP6",
  "N4.HOG0000965", "URP7",
  
  "N4.HOG0017043", "ACL5",
  "N4.HOG0003546", "BUD2",
  "N4.HOG0003479", "ATHB-8, ATHB-15/CNA",
  "N4.HOG0003480", "REV",
  "N4.HOG0003481", "PHB, PHV",
  "N4.HOG0010006", "SAC51, SACL1,2",
  "N4.HOG0008645", "SACL3",
  "N4.HOG0009777", "PHIP1",
  "N4.HOG0002615", "WOX4",
  "N4.HOG0007628", "HUP17",
  
  "N4.HOG0006359", "<strong>SOC1</strong>, XAL2, FYF, AGL19,71,72",
  "N4.HOG0002683", "ALMT12, ... [N4.HOG0002683]",
  "N4.HOG0004134", "UGTs [N4.HOG0004134]",
  
  "N4.HOG0006882", "ESB1, AT1G07730",
  "N4.HOG0000006", "MYB10,58,63",
  "N4.HOG0001926", "CYP86A1",
  "N4.HOG0002231", "SKS5,6,7,8,9,10",
  "N4.HOG0000382", "germin-like [N4.HOG0000382]"
)

devHogs <- devHogsIds %>% dplyr::left_join(getAnnotatedHogs(node = "N4", selectedHogs = devHogsIds$HOG))

devHogsMat <- devHogs %>% 
  dplyr::select(paste0(
    c(
      accSets$Spolyrhizas, 
      accSets$Sintermedias, 
      accSets$Lgibbas,
      accSets$Lminors,
      accSets$Ljaponicas,
      accSets$Lturioniferas,
      accSets$Waustralianas, 
      "Zmarina", "Zmays", "Osativa", "Athaliana", "Cdemersum"
    ), "_n")) %>% as.matrix()
rownames(devHogsMat) <- devHogs$plotLabel
devHogsMat <- t(devHogsMat)

devHogsHeatmap <- ComplexHeatmap::Heatmap(
  devHogsMat,
  column_labels = ComplexHeatmap::gt_render(colnames(devHogsMat)),
  row_labels = ComplexHeatmap::gt_render(
    rownames(devHogsMat) %>% 
      stringr::str_remove(pattern = "_n$") %>% 
      dwHighlighter()
  ),
  col = heatmapColFun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = factor(
    c(rep("A_lemnaceae", 11), rep("B_monocots", 3), rep("C_At_and_Cd", 2))
  ),
  column_split = factor(
    c(rep("A_roots", 14), rep("B_gross_development", 10), rep("C_stomata", 3), rep("D_metabolism", 5))
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(devHogsMat[i, j] > 3)
      grid.text(sprintf("%d", devHogsMat[i, j]), x, y, gp = gpar(fontsize = 5, col = "white"))
  },
  show_heatmap_legend = FALSE,
  rect_gp = gpar(col = "white", lwd = 2),
  width = unit(4.5, "in"), 
  height = unit(2.250, "in"),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  row_dend_width = unit(0.75, "cm"),
  column_dend_height = unit(0.75, "cm"),
  column_gap = unit(0.25, "cm"),
  row_names_side = "left",
  column_names_side = "top",
  column_names_rot = 45,
  row_dend_gp = gpar(col = "#cccccc"),
  column_dend_gp = gpar(col = "#cccccc"),
  column_title = c("root development", "frond development", "stomatal\nclosure", "metabolism"),
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 5),
  row_title = NULL,
  border = FALSE
)
pdf("plots/og-matrix.development.pdf", width = 7, height = 3.6)
draw(devHogsHeatmap)
draw(heatmapLegend(colFun = heatmapColFun, title = "copy\nnumber"), x = unit(0.95, "npc"), y = unit(0.5, "npc"), 
     just = c("centre"))
dev.off()


## 
## RNAi/methylation genes ----
## 
agoHogsPtn <- "(AGO)"
silencingHogsAthNamesPtn2 <- "(ALP1|AGO|CLF|CLSY|CMT|DCL|DML|DTF|DRB|DRM|EMF2|FIS2|HTR|HIS|JMJ|KTF1|MAIL|MAIN|MBD|MED|MET|MOM|MORC|NRP[A-E]|RAN|RDR|ROS|RTL|SHH|SHL|SmD1|SWC|SWN|TEK|VRN2)[0-9]*[^A-Za-z]"
silencingHogsAthNamesPtn3 <- "(ALP1|AGO|CLF|CLSY|CMT|DCL|DML|DTF|DRB|DRD|DRM|EMF2|FIS2|HIS|JMJ|KTF1|MAIL|MAIN|MBD|MED|MET|MOM|MORC|NRP[A-E]|RAN|RDR|ROS|RTL|SHH|SHL|SmD1|SWC|SWN|TEK|VRN2)[0-9]*[^A-Za-z]"

silencingHogsAth1MonoFilt <- 
  getAnnotatedHogs(node = "N4", selectedHogs = hogTablesAth1Mono$Lemnaceae_miss_ath_1mono %>% 
                     dplyr::pull(HOG)) %>% 
  dplyr::filter(
    At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = silencingHogsAthNamesPtn3) |
      Os_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = silencingHogsAthNamesPtn3) |
      Zm_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = silencingHogsAthNamesPtn3),
    ! Athaliana %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "(ATCG|ATMG)"),
    ! At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "(^RPS2|ATPP2-A3)"),
    ! At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "^(AT[0-9]G[0-9]+,?)+$"),
    # We omit AGO10 for this figure for simplicity - the AGO10 OG is split into two HOGs, most monocots in one, and At
    # in the other, but all Lemnaceae have 1 copy like At.
    ! At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "AGO10")
  ) %>% 
  addHogTablePlotLabels() %>% 
  dplyr::mutate(
    plotLabel = plotLabel %>% stringr::str_replace(
      pattern = "(, AT.G[0-9]{5})+",
      replacement = paste0(", ... [", HOG, "]"))
  ) %>% 
  dplyr::mutate(Lemnaceae_n = rowSums(dplyr::select(., paste0(accDw, "_n"))))
  
silencingHogsAth1MonoFilt <- silencingHogsAth1MonoFilt[
  naturalsort::naturalorder(silencingHogsAth1MonoFilt$plotLabel), 
]

silencingHogsAth1MonoMat <- silencingHogsAth1MonoFilt %>% 
  dplyr::select(paste0(c("Lemnaceae", "Zmarina", "Cdemersum", "Zmays", "Osativa", "Athaliana"), "_n")) %>% as.matrix()
rownames(silencingHogsAth1MonoMat) <- silencingHogsAth1MonoFilt$plotLabel
silencingHogsAth1MonoMat <- t(silencingHogsAth1MonoMat)

silencingHogsHeatmap <- ComplexHeatmap::Heatmap(
  silencingHogsAth1MonoMat,
  column_labels = ComplexHeatmap::gt_render(colnames(silencingHogsAth1MonoMat)),
  row_labels = ComplexHeatmap::gt_render(
    rownames(silencingHogsAth1MonoMat) %>% 
      stringr::str_remove(pattern = "_n$") %>% 
      ifelse(. == "Lemnaceae", ., sciNameLabeller(.))
  ),
  col = heatmapColFun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(silencingHogsAth1MonoMat[i, j] > 3)
      grid.text(sprintf("%d", silencingHogsAth1MonoMat[i, j]), x, y, gp = gpar(fontsize = 5, col = "white"))
  },
  show_heatmap_legend = FALSE,
  rect_gp = gpar(col = "white", lwd = 2),
  width = unit(2.75, "in"), 
  height = unit(0.75, "in"),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  row_dend_width = unit(0.75, "cm"),
  column_dend_height = unit(0.75, "cm"),
  row_names_side = "left",
  column_names_side = "top",
  column_names_rot = 45,
  row_dend_gp = gpar(col = "#cccccc"),
  column_dend_gp = gpar(col = "#cccccc"),
  row_title = NULL,
  column_title = NULL,
  border = FALSE
)
pdf("plots/og-matrix.silencing-horiz5.pdf", width = 4.2, height = 1.5)
draw(silencingHogsHeatmap)
draw(heatmapLegend(colFun = heatmapColFun, title = "copy\nnumber"), x = unit(0.95, "npc"), y = unit(0.26, "npc"), 
     just = c("centre"))
dev.off()  


silencingHogsLemnaceaeVarFilt <- 
  bind_rows(.id = "contrast", hogTablesLemnaceae) %>% 
  dplyr::filter(contrast %in% c("Lemnaceae_variable")) %>% 
  addHogTablePlotLabels() %>% 
  dplyr::mutate(
    plotLabel = plotLabel %>% stringr::str_replace(
      pattern = "(, AT.G[0-9]{5})+",
      replacement = paste0(", ... ", HOG))
  ) %>% 
  dplyr::mutate(Lemnaceae_n = rowSums(dplyr::select(., paste0(accDw, "_n")))) %>% 
  dplyr::filter(
    At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = silencingHogsAthNamesPtn2),
    ! Athaliana %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "(ATCG|ATMG)"),
    ! At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "(^RPS2|ATPP2-A3)"),
    ! At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "^(AT[0-9]G[0-9]+,?)+$"),
  )
silencingHogsLemnaceaeVarFilt <- silencingHogsLemnaceaeVarFilt[
  naturalsort::naturalorder(silencingHogsLemnaceaeVarFilt$plotLabel), 
]


silencingHogsLemnaceaeFilt <- 
  bind_rows(.id = "contrast", hogTablesLemnaceae$Lemnaceae_variable) %>% 
  dplyr::filter(! contrast %>% stringr::str_detect(pattern = "^Waustralianas")) %>% 
  addHogTablePlotLabels() %>% 
  dplyr::mutate(
    plotLabel = plotLabel %>% stringr::str_replace(
      pattern = "(, AT.G[0-9]{5})+",
      replacement = paste0(", ... ", HOG))
  ) %>%
  dplyr::filter(
    At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = silencingHogsAthNamesPtn2),
    ! Athaliana %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "(ATCG|ATMG)"),
    ! At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "(^RPS2|ATPP2-A3)"),
    ! At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "^(AT[0-9]G[0-9]+,?)+$"),
  )
silencingHogsLemnaceaeFilt <- silencingHogsLemnaceaeFilt[
  naturalsort::naturalorder(silencingHogsLemnaceaeFilt$plotLabel), 
]

silencingHogsLemnaceaeMat <- silencingHogsLemnaceaeFilt %>% 
  dplyr::select(paste0(c(accDw, "Zmarina", "Cdemersum", "Zmays", "Osativa", "Athaliana"), "_n")) %>% as.matrix()
rownames(silencingHogsLemnaceaeMat) <- silencingHogsLemnaceaeFilt$plotLabel
silencingHogsLemnaceaeMat <- t(silencingHogsLemnaceaeMat)

silencingHogsLemnaceaeMat <- silencingHogsLemnaceaeFilt %>% 
  dplyr::select(paste0(accAngio[-c(24:26)], "_n")) %>% as.matrix()
rownames(silencingHogsLemnaceaeMat) <- silencingHogsLemnaceaeFilt$plotLabel
silencingHogsLemnaceaeMat <- t(silencingHogsLemnaceaeMat)

silencingHogsLemnaceaeHeatmap <- ComplexHeatmap::Heatmap(
  silencingHogsLemnaceaeMat,
  column_labels = ComplexHeatmap::gt_render(colnames(silencingHogsLemnaceaeMat)),
  row_labels = ComplexHeatmap::gt_render(
    rownames(silencingHogsLemnaceaeMat) %>% 
      stringr::str_remove(pattern = "_n$") %>% 
      ifelse(. == "Lemnaceae", ., sciNameLabeller(.))
  ),
  col = heatmapColFun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = FALSE,
  rect_gp = gpar(col = "white", lwd = 2),
  width = unit(1.4, "in"), 
  height = unit(1.8, "in"),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  row_dend_width = unit(0.75, "cm"),
  column_dend_height = unit(0.75, "cm"),
  row_names_side = "left",
  column_names_side = "top",
  column_names_rot = 45,
  row_dend_gp = gpar(col = "#cccccc"),
  column_dend_gp = gpar(col = "#cccccc"),
  row_title = NULL,
  column_title = NULL,
  border = FALSE
)
pdf("plots/og-matrix.silencing-lemnaceae-variable.pdf", width = 3, height = 3.5)
draw(silencingHogsLemnaceaeHeatmap)
draw(heatmapLegend(colFun = heatmapColFun, title = "copy\nnumber"), x = unit(0.95, "npc"), y = unit(0.26, "npc"), 
     just = c("centre"))
dev.off()  

## 
## Argonaute genes ----
## 
argonauteFilt <- 
  getAnnotatedHogs(node = "N1", subgen = TRUE) %>% 
  dplyr::filter(At_syms %>% stringr::str_detect(pattern = ".*AGO[0-9]+.*")) %>% 
  addHogTablePlotLabels() %>% 
  dplyr::mutate(
    plotLabel = plotLabel %>% stringr::str_replace(
      pattern = "(, AT.G[0-9]{5})+",
      replacement = paste0(", ... ", HOG))
  ) %>% 
  dplyr::filter(
    At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = agoHogsPtn),
    ! (Athaliana %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "(ATCG|ATMG)")),
    ! (At_syms %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "(^RPS2|ATPP2-A3)"))
  )
argonauteFilt <- argonauteFilt[
  naturalsort::naturalorder(argonauteFilt$plotLabel), 
]


argonauteMat <- argonauteFilt %>% 
  dplyr::select(paste0(accAngioSubgen, "_n")) %>% 
  as.matrix()
rownames(argonauteMat) <- argonauteFilt$plotLabel
argonauteMat <- t(argonauteMat)

argonauteHeatmap <- ComplexHeatmap::Heatmap(
  argonauteMat,
  column_labels = ComplexHeatmap::gt_render(colnames(argonauteMat)),
  row_labels = ComplexHeatmap::gt_render(
    rownames(argonauteMat) %>% stringr::str_remove(pattern = "_n$") %>% dwHighlighter()
  ),
  col = heatmapColFun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = FALSE,
  rect_gp = gpar(col = "white", lwd = 2),
  width = unit(1, "in"), 
  height = unit(3.75, "in"),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  row_dend_width = unit(0.75, "cm"),
  column_dend_height = unit(0.75, "cm"),
  row_names_side = "left",
  column_names_side = "top",
  column_names_rot = 45,
  row_dend_gp = gpar(col = "#cccccc"),
  column_dend_gp = gpar(col = "#cccccc"),
  row_title = NULL,
  column_title = NULL,
  border = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%d", argonauteMat[i, j]), x, y, gp = gpar(fontsize = 5, fontface = "bold", color = "white"))
  }
)
pdf("plots/og-matrix.silencing-ago.pdf", width = 2.5, height = 4.5)
draw(argonauteHeatmap)
draw(heatmapLegend(colFun = heatmapColFun, title = "copy\nnumber"), x = unit(0.92, "npc"), y = unit(0.26, "npc"), 
     just = c("center"))
dev.off()

##
## Hybrids ----
## 
subgenVarWithNames <- 
  dplyr::bind_rows(
    hogTablesHybrids$Ljaponica_variable, 
    getAnnotatedHogs(c("N4.HOG0007627", "N4.HOG0007629"), subgen = TRUE, node = "N4")
  ) %>% 
  addHogTablePlotLabels() %>% 
  dplyr::mutate(
    plotLabel = plotLabel %>% stringr::str_replace(
      pattern = "AGL29, AGL40, AGL57, AGL60, AGL62, AGL91, AGL100, DIA",
      replacement = "<strong>AGL62</strong>,29,40,57,60,91,100, DIA"),
    plotLabel = plotLabel %>% stringr::str_replace(
      pattern = "(, AT.G[0-9]{5})+",
      replacement = paste0(", ... ", HOG))
  )

# This filter ensures that only HOGs missing in one or more hybrids but present in one or more "parents" are shown.
subgenVarWithNamesFilt <- subgenVarWithNames %>% 
  dplyr::filter(
    ! At_syms %>% is.na(),
    ! Athaliana %>% stringr::str_replace_na() %>% stringr::str_detect(pattern = "(ATCG|ATMG)"),
    ! At_syms %>% stringr::str_detect(pattern = "(^RPS2|ATPP2-A3)"),
    ! At_syms %>% stringr::str_detect(pattern = "^(AT[0-9]G[0-9]+,?)+$"),
  ) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    subgen_m_zeros_n = sum(c_across(paste0(accSetsDwSubgen$Ljaponicas_subM, "_n")) == 0),
    subgen_t_zeros_n = sum(c_across(paste0(accSetsDwSubgen$Ljaponicas_subT, "_n")) == 0),
    parents_m_non_zeros_n = sum(c_across(paste0(accSetsDwSubgen$Lminors, "_n")) > 0),
    parents_t_non_zeros_n = sum(c_across(paste0(accSetsDwSubgen$Lturioniferas, "_n")) > 0)
  ) %>% 
  ungroup() %>% 
  dplyr::filter(
    (subgen_m_zeros_n >= 1 & parents_m_non_zeros_n >= 1) | (subgen_t_zeros_n >= 1 & parents_t_non_zeros_n >= 1)
  )


subgenVarWithNamesMat <- subgenVarWithNamesFilt %>% 
  dplyr::select(paste0(c(accDwSubgen, "Athaliana", "Zmarina", "Cdemersum"), "_n")) %>% 
  as.matrix()
rownames(subgenVarWithNamesMat) <- subgenVarWithNamesFilt$plotLabel

subgenVarWithNamesMatSubset1 <- subgenVarWithNamesMat[ , 6:14]
subgenVarWithNamesMatSubset2 <- subgenVarWithNamesMat[ , c(1:5, 15:17)]

hybridVarHeatmapSubset1 <- ComplexHeatmap::Heatmap(
  subgenVarWithNamesMatSubset1,
  column_labels = ComplexHeatmap::gt_render(
    colnames(subgenVarWithNamesMatSubset1) %>% stringr::str_remove(pattern = "_n$") %>% dwHighlighter()
  ),
  row_labels = ComplexHeatmap::gt_render(rownames(subgenVarWithNamesMatSubset1)),
  top_annotation = HeatmapAnnotation(
    species = anno_block(
      height = unit(0.1, "cm"),
      gp = gpar(
        fill = c(
          pal4["Lminor7210"],
          pal4["Lminor7210"],
          pal4["Lturionifera9434"],
          pal4["Lturionifera9434"]
        ), 
        col = "white"
      )
    )
  ),
  col = heatmapColFun,
  row_split = 2,
  cluster_columns = FALSE,
  column_split = factor(
    c("Lm", "Lm", "LjM", "LjM", "LjM","LjT", "LjT", "LjT", "Lt"),
    levels = c("Lm", "LjM", "LjT", "Lt")
  ),
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D2",
  show_heatmap_legend = FALSE,
  rect_gp = gpar(col = "white", lwd = 2),
  width = unit(1.1, "in"), 
  height = unit(1.8, "in"),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 6),
  row_dend_width = unit(0.5, "cm"),
  column_dend_height = unit(0.5, "cm"),
  row_dend_gp = gpar(col = "#cccccc"),
  column_dend_gp = gpar(col = "#cccccc"),
  row_title = NULL,
  column_title = NULL,
  border = FALSE
)
hybridVarHeatmapSubset1RowOrder <- draw(hybridVarHeatmapSubset1) %>% ComplexHeatmap::row_order()
hybridVarHeatmapSubset2 <- ComplexHeatmap::Heatmap(
  subgenVarWithNamesMatSubset2,
  column_labels = ComplexHeatmap::gt_render(
    colnames(subgenVarWithNamesMatSubset2) %>% stringr::str_remove(pattern = "_n$") %>% dwHighlighter()
  ),
  row_labels = ComplexHeatmap::gt_render(rownames(subgenVarWithNamesMatSubset2)),
  top_annotation = HeatmapAnnotation(
    species = anno_block(
      height = unit(0.1, "cm"),
      gp = gpar(
        fill = c(
          pal4["all"],
          pal4["all"], 
          pal4["all"],
          pal4["Spolyrhiza9509"],
          pal4["all"],
          pal4["Waustraliana8730"],
          pal4["Lgibba7742a"]
        ), 
        col = "white"
      )
    )
  ),
  col = heatmapColFun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = factor(
    c("Sp", "Sp", "Si", "Wa", "Lg", "At", "Zm", "Cd"),
    levels = rev(c("Lg", "Wa", "Si", "Sp", "Zm", "Cd", "At"))
  ),
  show_heatmap_legend = FALSE,
  rect_gp = gpar(col = "white", lwd = 2),
  width = unit(1.1, "in"), 
  height = unit(1.8, "in"),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 6),
  row_dend_width = unit(0.75, "cm"),
  column_dend_height = unit(0.75, "cm"),
  row_dend_gp = gpar(col = "#cccccc"),
  column_dend_gp = gpar(col = "#cccccc"),
  row_title = NULL,
  column_title = NULL,
  border = FALSE
)
pdf("plots/og-matrix.hybrid-variable-fractionated.pdf", width = 3.8, height = 2.8)
draw(
  hybridVarHeatmapSubset2 + hybridVarHeatmapSubset1, 
  row_order = hybridVarHeatmapSubset1RowOrder %>% unlist(),
  ht_gap = unit(0, "mm")
)
draw(heatmapLegend(colFun = heatmapColFun, title = "copy\nnumber"), x = unit(0.9, "npc"), y = unit(0.2, "npc"), 
     just = c("centre"))
dev.off()


##
## Carbonic anhydrases ----
## 
caHogsIds <- tibble::tribble(
  ~HOG,            ~plotLabel,
  # Alpha
  "N4.HOG0003568", "ACA1",
  "N4.HOG0005831", "ACA4,ACA6",
  "N4.HOG0005832", "ACA3",
  "N4.HOG0005834", "ACA2,ACA7,ACA5,ACA8",
  "N4.HOG0021280", "AT3G54680",
  "N4.HOG0005835", "OsALPHACA3-7",
  "N4.HOG0018456", "LOC_Os09g28130",
  "N4.HOG0019005", "OsALPHACA10",
  "N4.HOG0019945", "OsALPHACA9",
  "N4.HOG0005833", "Zm00001d005844",
  "N4.HOG0024265", "ACA?",
  "N4.HOG0026154", "ACA?",
  "N4.HOG0026631", "Zm00001d046784",
  "N4.HOG0032495", "ACA?",
  # Beta
  "N4.HOG0003542", "BCA6,BCA5",
  "N4.HOG0003543", "ATBCA1,BETA CA2,BCA3,BCA4",
  # Gamma
  "N4.HOG0005632", "GAMMA CA1,2",
  "N4.HOG0005631", "GAMMA CA3",
  "N4.HOG0012071", "GAMMA CAL1,2",
  "N4.HOG0005633", "GCA?"
)
caHogs <- caHogsIds %>% dplyr::left_join(getAnnotatedHogs(node = "N4", selectedHogs = caHogsIds$HOG))

caHogsMat <- caHogs %>% 
  dplyr::select(paste0(c(accDw, "Zmarina", "Zmays", "Osativa", "Athaliana", "Cdemersum"), "_n")) %>% as.matrix()
rownames(caHogsMat) <- caHogs$plotLabel
caHogsMat <- t(caHogsMat)

caHogsHeatmap <- ComplexHeatmap::Heatmap(
  caHogsMat,
  column_labels = ComplexHeatmap::gt_render(colnames(caHogsMat)),
  row_labels = ComplexHeatmap::gt_render(
    rownames(caHogsMat) %>% 
      stringr::str_remove(pattern = "_n$") %>% 
      dwHighlighter()
  ),
  col = heatmapColFun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = factor(
    c(rep("A_lemnaceae", 11), rep("B_monocots", 3), rep("C_At_and_Cd", 2))
  ),
  column_split = factor(
    c(rep("Alpha", 14), rep("Beta", 2), rep("Gamma", 4))
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(caHogsMat[i, j] > 2)
      grid.text(sprintf("%d", caHogsMat[i, j]), x, y, gp = gpar(fontsize = 5, col = "white"))
  },
  show_heatmap_legend = FALSE,
  rect_gp = gpar(col = "white", lwd = 2),
  width = unit(3.25, "in"), 
  height = unit(2.50, "in"),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  row_dend_width = unit(0.75, "cm"),
  column_dend_height = unit(0.75, "cm"),
  row_names_side = "left",
  column_names_side = "top",
  column_names_rot = 45,
  row_dend_gp = gpar(col = "#cccccc"),
  column_dend_gp = gpar(col = "#cccccc"),
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 5),
  row_title = NULL,
  border = FALSE
)
pdf("plots/og-matrix.carbonic-anhydrases.pdf", width = 7, height = 3.6)
draw(caHogsHeatmap)
draw(heatmapLegend(colFun = heatmapColFun, title = "copy\nnumber"), x = unit(0.95, "npc"), y = unit(0.5, "npc"), 
     just = c("centre"))
dev.off()  





# 
# Time tree ----
#
pacman::p_load(phytools)
pacman::p_load(paleotree)
pacman::p_load(deeptime)
pacman::p_load(phyloseq)

# WGD events from:
# Guo, X., Wang, F., Fang, D. et al. The genome of Acorus deciphers insights into early monocot evolution. 
# Nat Commun 14, 3662 (2023). https://doi.org/10.1038/s41467-023-38836-4
wgdEvents <- tibble::tribble(
  ~node, ~wgd_name, ~wgd_time,
  54, "τ", 120.32,
  35, "γ", 116.51,
  # "σ", 97.83,
  # "ρ", 87.95,
  # "ρ", 32.10,
  # "ρ-1", 35.84,
  # "ρ-2", 58.72
)

otuLemnaceae <- accSetsDwSubgen$Lemnaceae
otuZosteraceae <- "Zmarina"
otuOtherMonocots <- setdiff(accSets$monocots, accSets$zmarina_dw)
otuDicots <- accSets$dicots
otuCeratophyllales <- accSets$Ceratophyllales
otuWaterlilies <- accSets$Nymphaeales
otuAmborellales <- accSets$Amborellales
otuAna <- c(accSets$Nymphaeales, accSets$Amborellales)

# Gene duplications at each node (via OrthoFinder)
dupSummary <- getOfDups() %>% 
  filter(Support >= 0.5) %>% 
  dplyr::count(Species.Tree.Node) %>% 
  dplyr::rename(node = Species.Tree.Node, gene_duplications = n)

angioTreeNodeLabs <- treeio::read.newick(
  "of-gmont-outgroup-ljsubgens/Species_Tree/SpeciesTree_rooted_node_labels.txt"
) %>% phyloseq::prune_taxa(taxa = ! phyloseq::taxa_names(.) %in% c("Gmontanum", "Aofficinalis"), .)

angioTree <- treeio::read.newick(
  "of-gmont-outgroup-ljsubgens/Species_Tree/SpeciesTree_rooted_node_labels.txt.ultrametric.tre"
) %>% 
  tidytree::as.ultrametric() %>% 
  phyloseq::prune_taxa(taxa = ! phyloseq::taxa_names(.) %in% c("Gmontanum", "Aofficinalis"), .) %>% 
  # Relabel the time tree with the OrthoFinder internal node labels
  treeio::as_tibble() %>% 
  dplyr::select(-label) %>% 
  dplyr::left_join(wgdEvents) %>% 
  dplyr::left_join(
    treeio::as_tibble(angioTreeNodeLabs) %>% dplyr::select(-branch.length), 
    by = dplyr::join_by(parent, node)
  ) %>% 
  dplyr::left_join(dupSummary, by = dplyr::join_by(label == node)) %>% 
  dplyr::rename(dups = gene_duplications) %>% 
  tidytree::as.treedata() %>% 
  ggtree::groupOTU(
    list(
      Amborellales = otuAmborellales,
      Nymphaeales = otuWaterlilies,
      Ceratophyllales = otuCeratophyllales, 
      Eudicots = otuDicots, 
      commelinids = otuOtherMonocots,
      Zosteraceae = otuZosteraceae,
      Lemnaceae = otuLemnaceae
    )
  )

mrcaLminors <- ggtree::MRCA(angioTree, "Lminor7210", "Lminor9252")
mrcaLturioniferas <- ggtree::MRCA(angioTree, "Lturionifera9434", "Ljaponica9421_T")
mrcaDuckweeds <- ggtree::MRCA(angioTree, "Lminor7210", "Spolyrhiza9509")
mrcaMonocots <- ggtree::MRCA(angioTree, "Ljaponica9421_M", "Macuminata")
mrcaOtherMonocots <- ggtree::MRCA(angioTree, "Sbicolor", "Macuminata")
mrcaDicots <- ggtree::MRCA(angioTree, "Cdemersum", "Athaliana")
mrcaWaterlilies <- ggtree::MRCA(angioTree, "Eferox", "Ncolorata")



# Diagnostic plot to show inner node labels
ggtree(angioTree) + geom_text2(aes(subset = !isTip, label = node), hjust = -.3) + geom_tiplab()

periodsCust <- deeptime::periods
epocsCust <- deeptime::epochs
epocsCust <- epocsCust %>% 
  dplyr::mutate(box_fill = dplyr::case_when(
    name == "Eocene" ~ "#59A2D7" %>% alpha(0.1),
    .default = "#ffffff00"
  ))

timeTreePlot <- angioTree %>% 
  tidytree::as_tibble() %>% 
  dplyr::mutate(is_tip = treeio::isTip(angioTree, node)) %>% 
  dplyr::mutate(label = dplyr::case_when(
    is_tip ~ label,
    .default = as.character(NA)
  )) %>% 
  treeio::as.treedata() %>% 
  ggtree::ggtree(aes(x = edge.length), size = 0.4, color = pal4["all"], right = FALSE) +
  geom_rect(
    data = epocsCust,
    aes(xmin = -min_age, xmax = -max_age, ymin = -Inf, ymax = Inf),
    fill = epocsCust$box_fill,
    inherit.aes = FALSE
  ) +
  ggtree::geom_rootedge(20, color = pal4["all"], size = 0.4) +
  ggtree::geom_tiplab(
    aes(label = dwHighlighter2Subgen(label), color = label), 
    size = 5/ggplot2::.pt, 
    as_ylab = FALSE, 
    offset = 0.5,
    parse=T
  ) +
  ggtree::geom_nodelab(
    node = "internal", 
    size = 5/ggplot2::.pt, 
    nudge_x = -4,
    nudge_y = 0.5
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.005, -0.025))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.23)),
    breaks = seq(-200, 0, 20), labels = abs(seq(-200, 0, 20))
  ) +
  scale_color_manual(values = pal4DarkerLabels, na.value = "#2e2e2e") +
  ggtree::theme_tree2() +
  deeptime::coord_geo(
    dat = list("epochs", "periods"),
    pos = list("bottom", "bottom"),
    abbrv = list(TRUE, FALSE),
    bord = list(NULL, c("bottom")),
    fill = list(
      (RColorBrewer::brewer.pal("Blues", n = 9) %>% colorspace::lighten(0.15))[c(1,2,3,4,5,6,7,8,8,9,9,9)],
      c("#E1EEFA", "#A7D3EA", "#4585C9", "#3B69B1", "#374B80")
    ),
    lab_color = list(
      c(rep("#374B80", 6), rep("#ffffff", 6)),
      c(rep("#374B80", 2), rep("#ffffff", 3))
    ),
    height = list(
      unit(10.5, "pt"), 
      unit(8, "pt")
    ),
    rot = list(90, 0),
    size = list(5/ggplot2::.pt, 5/ggplot2::.pt),
    expand = TRUE,
    xlim = c(-200, 0),
    ylim = c(1, treeio::Ntip(angioTree)),
    neg = TRUE,
    color = alpha("white", 0),
    lwd = 0.5
  ) +
  theme(
    plot.margin = margin(0, 2, 2, 4),
    axis.line = element_line(size = 0.25, color = "#2e2e2e"),
    axis.line.x = element_blank(),
    axis.ticks = element_line(size = 0.25, color = "#2e2e2e"),
    axis.text.x = element_markdown(
      colour = "#2e2e2e",
      size = 6
    ),
    axis.title.x = element_textbox(
      color = "#2e2e2e",
      size = 6,
      padding = margin(0, 4, 0, 4),
      margin = margin(4, 0, 0, 0)
    )
  ) + 
  guides(color = "none") +
  labs(x = "Age (Ma)")

transformClades <- function(treePlot) {
  newPlot <- treePlot %>% 
    ggtree::scaleClade(node = 51, scale = 0.75) %>%
    ggtree::scaleClade(node = 35, scale = 0.75) %>%
    ggtree::scaleClade(node = 32, scale = 0.75) +
    ggtree::geom_cladelab(
      node = 38, 
      label = "Lemnaceae", 
      offset = 27, 
      offset.text = 1.5, 
      extend = 0.4, 
      barcolour = pal4Light["all"], 
      textcolour = pal4Light["all"], 
      barsize = 0.25,
      fontsize = 5/ggplot2::.pt, 
      family = "Helvetica",
      fontface = "bold"
    ) +
    ggtree::geom_cladelab(
      node = 7, 
      label = "Zosteraceae", 
      offset = 27, 
      offset.text = 1.5, 
      extend = 0.4, 
      barcolour = pal4Light["all"], 
      textcolour = pal4Light["all"], 
      barsize = 0.25,
      fontsize = 5/ggplot2::.pt, 
      family = "Helvetica"
    ) +
    ggtree::geom_cladelab(
      node = 51, 
      label = "commelinids", 
      offset = 27, 
      offset.text = 1.5, 
      extend = 0.4, 
      barcolour = pal4Light["all"], 
      textcolour = pal4Light["all"], 
      barsize = 0.25,
      fontsize = 5/ggplot2::.pt, 
      family = "Helvetica"
    ) +
    ggtree::geom_cladelab(
      node = 35, 
      label = "Eudicots", 
      offset = 27, 
      offset.text = 1.5, 
      extend = 0.4, 
      barcolour = pal4Light["all"], 
      textcolour = pal4Light["all"], 
      barsize = 0.25, 
      fontsize = 5/ggplot2::.pt, 
      family = "Helvetica"
    ) +
    ggtree::geom_cladelab(
      node = 4, 
      label = "Ceratophyllales", 
      offset = 27, 
      offset.text = 1.5, 
      extend = 0.4, 
      barcolour = pal4Light["all"], 
      textcolour = pal4Light["all"], 
      barsize = 0.25,
      fontsize = 5/ggplot2::.pt, 
      family = "Helvetica"
    ) +
    ggtree::geom_cladelab(
      node = 32, 
      label = "Nymphaeales", 
      offset = 27, 
      offset.text = 1.5, 
      extend = 0.4, 
      barcolour = pal4Light["all"], 
      textcolour = pal4Light["all"], 
      barsize = 0.25,
      fontsize = 5/ggplot2::.pt, 
      family = "Helvetica"
    ) +
    ggtree::geom_cladelab(
      node = 1, 
      label = "Amborellales", 
      offset = 27, 
      offset.text = 1.5, 
      extend = 0.4, 
      barcolour = pal4Light["all"], 
      textcolour = pal4Light["all"], 
      barsize = 0.25,
      fontsize = 5/ggplot2::.pt, 
      family = "Helvetica"
    )
  return(newPlot)
}


timeTreePlot <- revts(timeTreePlot)
timeTreePlot2 <- timeTreePlot %>% transformClades()
nLayers <- length(timeTreePlot2$layers)
timeTreePlot2$layers <- c(timeTreePlot2$layers[3], timeTreePlot2$layers[1:2], timeTreePlot2$layers[4:nLayers])
timeTreePlot3 <- timeTreePlot %>% ggtree::flip(37, 51) %>% transformClades()
timeTreePlot3$layers <- c(timeTreePlot3$layers[3], timeTreePlot3$layers[1:2], timeTreePlot3$layers[4:nLayers])

pdf("plots/timetree.wider.2.pdf", width = 7, height = 3.25)
timeTreePlot2
dev.off()
pdf("plots/timetree.wider.2.rotated.pdf", width = 7, height = 3.25)
timeTreePlot3
dev.off()

##
## Zoomed tree ----
##
angioTreeZoomLem <- phyloseq::prune_taxa(
  taxa = phyloseq::taxa_names(timeTree) %in% 
    accDwSubgen,
  timeTree
)

timeTreeZoomPlot <- angioTreeZoomLem %>% 
  ggtree::ggtree(aes(x = edge.length), size = 0.4, color = unname(pal[1]), right = FALSE) + 
  ggtree::geom_tiplab(
    aes(label = dwHighlighter2Subgen(label), color = label), 
    size = 5/ggplot2::.pt, 
    as_ylab = FALSE, 
    offset = 0.2,
    parse=T
  ) +
  ggtree::geom_rootedge(5, color = "#2e2e2e") +
  deeptime::coord_geo(
    dat = list("epochs"),
    abbrv = list(TRUE, TRUE),
    fill = list(
      (RColorBrewer::brewer.pal("Greys", n = 9) %>% colorspace::lighten(0.15))[2:3]
    ),
    height = unit(12, "pt"),
    expand = TRUE, 
    xlim = c(-60, 0),
    center_end_labels = TRUE,
    ylim = c(1, treeio::Ntip(angioTreeZoomLem)), 
    neg = TRUE,
    size = 5/ggplot2::.pt,
    color = alpha("white", 0)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.25)),
    breaks = seq(-60, 0, 5), labels = abs(seq(-60, 0, 5))
  ) +
  scale_color_manual(values = pal4, na.value = "#2e2e2e", labels = sciNameLabeller) +
  ggtree::theme_tree2() +
  theme(
    plot.margin = margin(0, 2, 2, 4),
    axis.line = element_line(size = 0.25, color = "#2e2e2e"),
    axis.ticks = element_line(size = 0.25, color = "#2e2e2e"),
    axis.text.x = element_markdown(
      colour = "#2e2e2e",
      size = 6
    ),
    axis.title.x = element_textbox(
      color = "#2e2e2e",
      size = 6,
      padding = margin(0, 4, 0, 4),
      margin = margin(4, 0, 0, 0)
    )
  ) + 
  guides(color = "none") +
  labs(x = "Age (Ma)")
pdf("plots/timetree.wider.pdf", width = 3.875, height = 4.2)
revts(timeTreeZoomPlot)
dev.off()




