# Common theme elements, globals and functions
source("../dw-genomes-common.R")
theme_set(project_theme)

dir.create("plots", showWarnings = FALSE, recursive = TRUE)

###
### Load mosDepth coverage summaries ----
### 
mdCovFiles <- list.files(path = "coverage-2022", pattern = ".mosdepth.summary.txt$", full.names = TRUE)
names(mdCovFiles) <- basename(mdCovFiles) %>% stringr::str_replace("-ont-gdna.202106.min1k.q30.mosdepth.summary.txt$", "")
mdCov <- purrr::map_dfr(
  mdCovFiles,
  .id = "acc",
  read_tsv
) %>% 
  dplyr::filter(chrom %>% stringr::str_detect("^chr[0-9]+[MT]?$")) %>% 
  dplyr::mutate(
    acc = shortAccToLong[acc],
    chrom_number = stringr::str_replace(chrom, pattern = "^chr(([0-9]+|[MC])).*", replacement = "\\1"),
    subgen = dplyr::if_else(
      stringr::str_detect(chrom, "chr[0-9]+[MT]"), 
      stringr::str_extract(chrom, "[MT]$"),
      ""
    )
  ) %>% 
  dplyr::group_by(acc) %>% 
  dplyr::mutate(
    tot_bases = sum(bases),
    tot_length = sum(length)
  ) %>% 
  dplyr::group_by(acc, subgen) %>% 
  dplyr::mutate(
    tot_length_subgen = sum(length),
    mean_subgen = sum(bases) / sum(length)
  ) %>% 
  dplyr::group_by(acc, chrom_number) %>% 
  dplyr::mutate(
    tot_bases_chr = sum(bases),
    tot_length_chr = sum(length),
    mean_chr = sum(bases / length),
    mean_norm = mean / min(mean_subgen)
  ) %>% 
  dplyr::ungroup()

covGdna1k <- bwCovFilesToTibble(list.files(path = 'coverage-2022', pattern = '*.1000bp.q30.bw$', full.names = TRUE))
covGdna10k <- bwCovFilesToTibble(list.files(path = 'coverage-2022', pattern = '.10000bp.q30.bw$', full.names = TRUE))
covGdna100k <- bwCovFilesToTibble(list.files(path = 'coverage-2022', pattern = '.100000bp.q30.bw$', full.names = TRUE))

###
### Load variant vcfs ----
###
vcfToDf <- function (vcfFile, formatFields = NULL) {
  vcfDf <- vcfR2tidy(
    read.vcfR(vcfFile, verbose = FALSE), single_frame = TRUE, info_types = TRUE,
    format_types = TRUE, format_fields = formatFields
  ) %>% 
    .$dat %>% 
    tibble::as_tibble() %>%
    rename(chrom = CHROM) %>%
    dplyr::filter(stringr::str_detect(chrom, "^chr[0-9]+[A-Z]?")) %>%
    chromNameToNumberAndSubgen() %>% 
    dplyr::relocate(subgen, .after = chrom) %>%
    addCumIndex(cols = "POS") %>%
    tidyr::unite(col = accSubgen, sep = "_", remove = FALSE, accession, subgen) %>% 
    dplyr::mutate(
      accession = accSubgen %>% stringr::str_remove(pattern = "_NA") %>% factor(),
      FILTER = factor(FILTER)
    ) %>% 
    dplyr::select(!c(accSubgen, subgen)) 
  return(vcfDf)
}

vcfToDfCounts <- function(acc, vcfFile, genesIntFile, intactTesIntFile, svCalls = FALSE) {
  vcfDf <- vcfR2tidy(
    read.vcfR(vcfFile, verbose = FALSE, cols = seq(1:9)), 
    single_frame = TRUE,
    format_fields = "GT",
    alleles = FALSE
  ) %>% 
    .$dat %>% 
    tibble::as_tibble() %>%
    dplyr::left_join(
      readr::read_tsv(
        genesIntFile,
        col_names = if (svCalls) "ID" else c("CHROM", "POS")
      ) %>% 
        dplyr::mutate(intersectsGenes = TRUE)
    ) %>% 
    dplyr::left_join(
      readr::read_tsv(
        intactTesIntFile,
        col_names = if (svCalls) "ID" else c("CHROM", "POS")
      ) %>% 
        dplyr::mutate(intersectsTes = TRUE)
    ) %>% 
    dplyr::rename(chrom = CHROM) %>%
    dplyr::filter(
      stringr::str_detect(chrom, "^chr[0-9]+[A-Z]?"),
      dplyr::if_all(.cols = dplyr::any_of("SVLEN"), .fns = ~ .x > 50)
    ) %>%
    dplyr::mutate(
      accession = acc
    ) %>% 
    chromNameToNumberAndSubgen() %>% 
    dplyr::relocate(subgen, .after = chrom) %>%
    tidyr::unite(col = accSubgen, sep = "_", remove = FALSE, accession, subgen) %>% 
    dplyr::mutate(accession = accSubgen %>% stringr::str_remove(pattern = "_NA") %>% factor()) %>% 
    dplyr::count(name = "n_het_vars", accession, intersectsGenes, intersectsTes) %>% 
    dplyr::mutate(
      intersection = dplyr::case_when(intersectsTes ~ "intact TEs", intersectsGenes ~ "genes", TRUE ~ "none") %>% 
        forcats::fct_relevel("none", "genes", "intact TEs")
    ) %>%
    dplyr::group_by(accession, intersection) %>% 
    dplyr::summarize(n_het_vars = sum(n_het_vars))
  return(vcfDf)
}

snpVcfFiles <- list.files(path = 'vcf', pattern = '.short-variants.filt.het.vcf.gz$', full.names = TRUE)
names(snpVcfFiles) <- shortAccToLong[basename(snpVcfFiles) %>% stringr::str_replace("([^-.]+).*$", "\\1")]
snpGenesIntFiles <- list.files(path = 'intersections', pattern = '.short-variants.filt.het.intersect.genes.txt.gz$', full.names = TRUE)
names(snpGenesIntFiles) <- shortAccToLong[basename(snpVcfFiles) %>% stringr::str_replace("([^-.]+).*$", "\\1")]
snpIntactTesIntFiles <- list.files(path = 'intersections', pattern = '.short-variants.filt.het.intersect.intact-tes.txt.gz$', full.names = TRUE)
names(snpIntactTesIntFiles) <- shortAccToLong[basename(snpVcfFiles) %>% stringr::str_replace("([^-.]+).*$", "\\1")]

snpsCounts <- purrr::imap_dfr(
  snpVcfFiles, 
  ~ vcfToDfCounts(
    acc = .y, 
    vcfFile = .x, 
    genesIntFile = snpGenesIntFiles[.y], 
    intactTesIntFile = snpIntactTesIntFiles[.y]
  )
) %>% 
  dplyr::mutate(varType = "SNVs")

svVcfFiles <- list.files(path = 'vcf', pattern = '.structural-variants.filt.het.vcf.gz$', full.names = TRUE)
names(svVcfFiles) <- shortAccToLong[basename(svVcfFiles) %>% stringr::str_replace("([^-.]+).*$", "\\1")]
svGenesIntFiles <- list.files(path = 'intersections', pattern = '.structural-variants.filt.het.intersect.genes.txt.gz$', full.names = TRUE)
names(svGenesIntFiles) <- shortAccToLong[basename(svVcfFiles) %>% stringr::str_replace("([^-.]+).*$", "\\1")]
svIntactTesIntFiles <- list.files(path = 'intersections', pattern = '.structural-variants.filt.het.intersect.intact-tes.txt.gz$', full.names = TRUE)
names(svIntactTesIntFiles) <- shortAccToLong[basename(svVcfFiles) %>% stringr::str_replace("([^-.]+).*$", "\\1")]

svsCounts <- purrr::imap_dfr(
  svVcfFiles, 
  ~ vcfToDfCounts(
    acc = .y, 
    vcfFile = .x, 
    genesIntFile = svGenesIntFiles[.y], 
    intactTesIntFile = svIntactTesIntFiles[.y],
    svCalls = TRUE
  )
) %>% 
  dplyr::mutate(varType = "SVs")

###
### Load het variant density bigwigs ----
###
bwVarDensityFilesToTibble <- function(bigWigFiles) {
  names(bigWigFiles) <- basename(bigWigFiles) %>% stringr::str_replace("([^-.]+).*$", "\\1")
  purrr::map_dfr(bigWigFiles, ~ import.bw(.x) %>% tibble::as_tibble(), .id = "accession") %>%
    tibble::as_tibble() %>%
    dplyr::rename_with(~ stringr::str_replace_all(.x, pattern = "^value.(.+)$", replacement = "\\1")) %>%
    dplyr::select(-strand) %>%
    dplyr::rename(
      chrom = seqnames,
      snv_ct = score
    ) %>%
    dplyr::filter(stringr::str_detect(chrom, "^chr[0-9]+")) %>%
    dplyr::mutate(accession = accession %>% shortAccToLong[.] %>% factor()) %>% 
    chromNameToNumberAndSubgen() %>%
    addCumIndex() %>% 
    subgenToSpecies() %>% 
    dplyr::mutate(accession = accession %>% forcats::fct_relevel(c("Lgibba7742a", "Lminor7210", "Lminor9252")))
}

hetVarDensity100k <- bwVarDensityFilesToTibble(
  list.files(path = 'bigwig/variants', pattern = '*short-variants.filt.het.density.100k.bw$', full.names = TRUE)
)


###
### Load Infernal/RFAM ncRNA annotations ----
### 
infernalFiles <- list.files(
  path = 'gff/infernal', pattern = 'infernal.no-olap.gff.gz$', full.names = TRUE)
names(infernalFiles) <- basename(infernalFiles) %>% stringr::str_replace("([^.]+).*$", "\\1")

infernalAnns <- purrr::map_dfr(infernalFiles, as.data.frame(rtracklayer::import.gff), .id = "accession") %>%
  tibble::as_tibble() %>%
  dplyr::rename_with(~ stringr::str_replace_all(.x, pattern = "^value.(.+)$", replacement = "\\1")) %>%
  dplyr::rename(chrom = seqnames) %>% 
  chromNameToNumberAndSubgen() %>% 
  dplyr::mutate(accession = shortAccToLong[accession]) %>%
  dplyr::select(accession, subgen, chrom, chrom_number, start, end, width, strand, type, score, evalue, desc) %>%
  dplyr::ungroup() %>% 
  addCumIndex() %>% 
  dplyr::left_join(
    accDwChromLengths %>% dplyr::select(accession, chrom, chrom_length)
  ) %>% 
  dplyr::mutate(
    accession = accession %>% factor(levels = accDw) %>% droplevels()
  )


###
### Load BED telomere annotations ----
### 

telo1xBedFiles <- list.files(
  path = 'bed/telomeres', pattern = '(asm202106v1|asm3v1).telomeres.1x.bed.gz$', full.names = TRUE)
names(telo1xBedFiles) <- basename(telo1xBedFiles) %>% stringr::str_replace("([^.]+).*$", "\\1")

telomeres1x <- purrr::map_dfr(telo1xBedFiles, as.data.frame(import.bed), .id = "accession") %>%
  dplyr::as_tibble() %>%
  dplyr::rename_with(~ stringr::str_replace_all(.x, pattern = "^value.(.+)$", replacement = "\\1")) %>%
  dplyr::rename(chrom = seqnames) %>% 
  chromNameToNumberAndSubgen() %>% 
  dplyr::mutate(accession = shortAccToLong[accession]) %>%
  dplyr::select(-width, -name) %>% 
  dplyr::select(accession, subgen, chrom, tidyselect::everything()) %>%
  dplyr::ungroup() %>% 
  addCumIndex() %>% 
  dplyr::left_join(
    accDwChromLengths %>% dplyr::select(accession, chrom, chrom_length)
  ) %>%
  subgenToSpecies() %>% 
  dplyr::left_join(
    accTblSubgen %>% dplyr::select(accession, genus)
  ) %>%
  dplyr::mutate(
    accession = accession %>% factor(levels = accDwSubgen) %>% droplevels(),
    score = 1
  )

# Some stats
telEndSummary <- telomeres1x %>%
  dplyr::filter(
    ! is.na(accession),
    chrom %>% stringr::str_detect("^chr[0-9]+")
  ) %>%
  dplyr::group_by(accession, chrom, chrom_number, genus) %>% 
  dplyr::summarize(
    tel_start = any(start <= 10000),
    tel_end = any(start >= chrom_length - 10000)
  )


### Pool Plot ----
pos = position_jitterdodge(seed = 3, jitter.width = 0.3, dodge.width = 0.7)
hybridAccToParents <- function(accSubgen) {
  parent <- dplyr::case_when(
    stringr::str_detect(accSubgen, "M$") ~ "Lminor9252",
    stringr::str_detect(accSubgen, "T$") ~ "Lturionifera9434",
    TRUE ~ accSubgen
  )
}
plotChrCovPool <- mdCov %>%
  dplyr::filter(acc %in% c("Lminor9252", "Lturionifera9434", accSets$Ljaponicas)) %>% 
  dplyr::mutate(
    accSubgen = dplyr::if_else(subgen == "", acc, paste(acc, subgen, sep = "_")),
    acc = acc %>% forcats::fct_relevel("Lminor9252", "Lturionifera9434")
  ) %>%
  ggplot(
    aes(
      y = mean_norm,
      x = acc, 
      color = pal4[hybridAccToParents(accSubgen)], 
      fill = pal4Light[hybridAccToParents(accSubgen)]
    )
  ) +
  geom_hline(yintercept = 1, size = 0.25, color = "lightgrey") +
  geom_hline(yintercept = 2, size = 0.25, color = "lightgrey") +
  geom_point(
    shape = 21,
    position = pos,
    size = 3.25,
    alpha = 0.3
  ) +
  geom_text(
    aes(label = chrom_number), 
    position = pos, 
    size = 5/ggplot2::.pt,
    alpha = 1,
    fontface = "bold",
    key_glyph = "blank"
  ) +
  stat_summary(
    fun = median, fun.min = median, fun.max = median,
    geom = "crossbar", width = 0.75, size = 0.25, key_glyph = "blank",
    position = position_dodge(width = 0.7)) +
  scale_y_continuous(
    limits = c(0.75, 2.5),
    expand = expansion(mult = c(0.025, 0.05)), 
    breaks = c(0.75, 1, 1.25, 1.5, 1.75, 2, 2.25),
    labels = c(NA, 1, NA, 1.5, NA, 2, NA)
  ) +
  scale_x_discrete(labels = ~ stringr::str_replace(.x, "_$", "") %>% dwHighlighter()) +
  coord_capped_cart(left = 'both', ylim = c(0.75, 2.25)) +
  scale_color_identity(
    guide = "none"
  ) +
  scale_fill_identity(guide = "none") +
  theme(
    plot.margin = unit(c(0.25, 0, 0.25, 0), "mm"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.title.y = element_textbox(halign = 0.5)
  ) +
  labs(title = NULL,
       x = "Pseudomolecules", 
       y = "normalized mean chromosome <br/>read depth")
pdf("plots/read-cov-pool-plot.v3.pdf", width = 2.5, height = 2)
plotChrCovPool
dev.off()



### Plot cov and density ----

xAxisLblDf = subgenToSpecies(covGdna10k) %>% 
  dplyr::group_by(accession, chrom_number) %>% 
  dplyr::summarize(
    left = min(start_cum),
    right = max(start_cum),
    width = max(start_cum) - min(start_cum),
    center = (max(start_cum) + min(start_cum) ) / 2
  )


plotReadCovBases4 <- 
  hetVarDensity100k %>% 
  ggplot() +
  facet_wrap(
    scales = "fixed",
    strip.position = "right",
    accession ~ .,
    ncol = 1,
    labeller = labeller(
      .multi_line = TRUE,
      accession = ~ dwHighlighterMultiline(., pal = pal4DarkerLabels)
    )
  ) +
  geom_hline(yintercept = c(0, 0.5, 1, 1.5, 2), size = 0.25, color = "#dedede") + 
  geom_segment(
    data = infernalAnns %>%
      dplyr::filter(
        chrom %>% stringr::str_detect("^chr[0-9]+"),
        type == "LSU_rRNA_eukarya",
        accession != "Spolyrhiza9509",
        score < 3000 & score >= 2000,
      ) %>%
      subgenToSpecies(),
    aes(
      x = start_cum + (width / 2),
      xend = start_cum + (width / 2),
      y = 0,
      yend = Inf
    ),
    color = pal4Light["all"],
    size = 0.5
  ) +
  geom_segment(
    data = infernalAnns %>%
      dplyr::filter(
        chrom %>% stringr::str_detect("^chr[0-9]+"),
        type == "LSU_rRNA_eukarya",
        accession != "Spolyrhiza9509",
        score >= 3000
      ) %>%
      subgenToSpecies(),
    aes(
      x = start_cum + (width / 2),
      xend = start_cum + (width / 2),
      y = 0,
      yend = Inf
    ),
    color = pal4["all"],
    size = 0.75
  ) +
  geom_point(
    data = hetVarDensity100k %>% dplyr::filter(snv_ct != 0),
    aes(
      x = start_cum,
      y = snv_ct / 1000,
      color = interaction(accession, chrom_number %% 2, lex.order = TRUE)
    ),
    size = 0.001,
    alpha = 1
  ) +
  geom_line(
    data = covGdna100k %>%
      subgenToSpecies() %>%
      dplyr::mutate(accession = accession %>% forcats::fct_relevel(c("Lgibba7742a", "Lminor7210", "Lminor9252"))),
    aes(
      x = start_cum,
      y = coverage / coverage_subgenome_mode,
    ),
    size = 0.1,
    #color = "#4e4e4e",
    color = "#2AB7CA",
    alpha = 1
  ) +
  geom_rect(
    data = xAxisLblDf,
    aes(
      xmin = left,
      xmax = right,
      color = interaction(accession, chrom_number %% 2, lex.order = TRUE)
    ),
    linewidth = 0.25,
    ymin = -0.4,
    ymax = -0.0,
    fill = "#FFFFFF"
  ) +
  geom_text(
    data = telomeres1x %>%
      dplyr::filter(
        chrom %>% stringr::str_detect("^chr[0-9]+"),
        start <= 10000 | start >= chrom_length - 10000,
        accession != "Spolyrhiza9509"
      ),
    aes(
      y = -0.18,
      x = dplyr::if_else(start <= 10000, start_cum + 2e6, end_cum - 2e6),
      color = interaction(accession, chrom_number %% 2, lex.order = TRUE),
      label = dplyr::if_else(start > 10000, "\u25b6", "\u25c0")
    ),
    family = "AppleGothic Regular",
    size = 1
  ) +
  geom_text(
    data = xAxisLblDf, 
    aes(
      x = center, 
      y = -0.85,
      #y = -1,
      label = chrom_number, 
      color = interaction(accession, chrom_number %% 2, lex.order = TRUE)
    ),
    size = 5/ggplot2::.pt,
    fontface = "bold"
  ) +
  labs(
    title = NULL,
    x = NULL, 
    y = "heterozygous short variants / 100 bp<br>read depth (100 Kbp bins) / read depth mode"
  ) +
  scale_color_manual(
    values = c(pal4 %>% setNames(paste0(names(pal4), ".0")), pal4Dark %>% setNames(paste0(names(pal4Dark), ".1"))), 
    guide = "none"
  ) +
  scale_fill_manual(
    values = c(pal4 %>% setNames(paste0(names(pal4), ".0")), pal4Dark %>% setNames(paste0(names(pal4Dark), ".1"))), 
    guide = "none"
  ) +
  scale_y_continuous(
    oob = ~ scales::oob_squish(.x, range = c(-1.25, 2)),
    expand = expansion(add = c(0, 0)), 
    position = "left",
    limits = c(-1.25, 2),
    breaks = c(0, 2)
  ) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 1, scale_cut = scales::cut_si("bp")), 
    position = "top",
    expand = expansion(mult = c(0, 0)),
    limits = c(0, 500e6),
    breaks = seq(0, 500e6, by = 100e6)
  ) +
  theme(
    panel.spacing.x = unit(0, "mm"),
    panel.spacing.y = unit(2, "mm"),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text.y = element_markdown(angle = 0, vjust = 1, lineheight = 4 / ggplot2::.pt),
    plot.margin = unit(c(0.25, 1, 2.25, 0), "mm"),
    axis.text.x = element_markdown(hjust = c(0.15, 0.5, 0.5, 0.5, 0.5, 0.9)),
    axis.line.y.left = element_blank(),
    axis.title.y.left = element_markdown(lineheight = 4 / ggplot2::.pt, vjust = 0.5, angle = 90)
  ) +
coord_cartesian(ylim = c(-1.25, 2), xlim = c(0, 500e6), clip = "off")
png("plots/read-cov2.png", width = 7, height = 5.5, units = "in", res = 600, type = "quartz")
plotReadCovBases4
dev.off()


#### coverage scatter ----
plotReadCovBases3 <- 
  covGdna10k %>% subgenToSpecies() %>% 
  dplyr::mutate(accession = accession %>% forcats::fct_relevel(c("Lgibba7742a", "Lminor7210", "Lminor9252"))) %>% 
  ggplot() +
  facet_wrap(
    scales = "fixed",
    strip.position = "bottom",
    accession ~ .,
    ncol = 1,
    labeller = labeller(
      .multi_line = FALSE,
      accession = sciNameLabeller
    )
  ) +
  geom_hline(yintercept = c(0, 0.5, 1, 1.5, 2), size = 0.25, color = "#dedede") + 
  geom_vline(
    data = infernalAnns %>%
      dplyr::filter(
        chrom %>% stringr::str_detect("^chr[0-9]+"),
        type == "LSU_rRNA_eukarya",
        accession != "Spolyrhiza9509",
        score < 3000 & score >= 2000,
      ) %>%
      subgenToSpecies(),
    aes(
      xintercept = start_cum + (width / 2)
    ),
    color = pal4Light["Spolyrhiza9509"],
    size = 0.5
  ) +
  geom_vline(
    data = infernalAnns %>%
      dplyr::filter(
        chrom %>% stringr::str_detect("^chr[0-9]+"),
        type == "LSU_rRNA_eukarya",
        accession != "Spolyrhiza9509",
        score >= 3000
      ) %>%
      subgenToSpecies(),
    aes(
      xintercept = start_cum + (width / 2)
    ),
    color = pal4Dark["Spolyrhiza9509"],
    size = 0.75
  ) +
  geom_point(
    aes(
      x = start_cum,
      y = coverage / coverage_subgenome_mode,
      color = interaction(accession, chrom_number %% 2, lex.order = TRUE),
    ),
    size = 0.001,
    alpha = 1
  ) +
  geom_text(
    data = telomeres1x %>%
      dplyr::filter(
        chrom %>% stringr::str_detect("^chr[0-9]+"),
        start <= 10000 | start >= chrom_length - 10000,
        accession != "Spolyrhiza9509"
      ),
    aes(
      y = -0.3,
      x = dplyr::if_else(start <= 10000, start_cum + 2e6, end_cum - 2e6),
      color = interaction(accession, chrom_number %% 2, lex.order = TRUE),
      label = dplyr::if_else(start > 10000, "\u25b6", "\u25c0")
    ),
    family = "AppleGothic Regular",
    size = 1
  ) +
  geom_text(
    data = xAxisLblDf, 
    aes(
      x = center, 
      y = 1, 
      label = chrom_number, 
      color = interaction(accession, (chrom_number + 1) %% 2, lex.order = TRUE)
    ),
    size = 5/ggplot2::.pt,
    fontface = "bold"
  ) +
  labs(
    title = NULL,
    x = NULL, 
    y = "read depth (10 Kbp bins) / read depth mode"
  ) +
  scale_color_manual(
    values = c(pal4 %>% setNames(paste0(names(pal4), ".0")), pal4Dark %>% setNames(paste0(names(pal4Dark), ".1"))), 
    guide = "none"
  ) +
  scale_fill_manual(
    values = c(pal4 %>% setNames(paste0(names(pal4), ".0")), pal4Dark %>% setNames(paste0(names(pal4Dark), ".1"))), 
    guide = "none"
  ) +
  scale_y_continuous(
    oob = scales::squish,
    expand = expansion(add = c(0.48, 0)), 
    position = "left",
    breaks = c(0, 2)
  ) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 1, scale_cut = scales::cut_si("bp")), 
    position = "top",
    expand = expansion(mult = c(0.015, 0)),
    limits = c(0, 500e6),
    breaks = seq(0, 500e6, by = 100e6)
  ) +
  theme(
    panel.spacing.x = unit(0, "mm"),
    panel.spacing.y = unit(0.25, "mm"),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(fill = "white"),
    plot.margin = unit(c(0.25, 1, 0.25, 0), "mm"),
    axis.text.x = element_markdown(hjust = c(0.15, 0.5, 0.5, 0.5, 0.5, 0.9)),
  ) +
  coord_cartesian(ylim = c(0, 2))

#### density ----
plotReadCovDensitySubgens <- 
  subgenToSpecies(covGdna1k) %>% 
  dplyr::mutate(accession = accession %>% forcats::fct_relevel(c("Lgibba7742a", "Lminor7210", "Lminor9252"))) %>% 
  ggplot() +
  geom_density(
    aes(
      x = coverage_norm,
      color = accession,
      fill = accession
    ),
    size = 0.3,
    adjust = 2
  ) +
  facet_wrap(
    strip.position = "bottom",
    accession ~ .,
    ncol = 1,
    labeller = labeller(
      .multi_line = FALSE,
      species = sciNameLabeller
    )
  ) +
  labs(title = NULL,
       x = NULL,
       y = "normalized read depth density (1 Kbp bins)") +
  scale_y_continuous(
    position = "right",
    expand = expansion(mult = c(0.15, 0.05)),
  ) +
  scale_x_continuous(
    position = "top",
    expand = expansion(mult = c(0.1, 0.1)),
    limits = c(0, 3),
    breaks = c(0, 1, 2, 3)
  ) +
  scale_color_manual(values = pal4Dark, guide = "none") +
  scale_fill_manual(
    values = pal4, 
    labels = sciNameLabeller,
    guide = "none"
  ) +
  theme(
    panel.spacing.x = unit(0, "mm"),
    panel.spacing.y = unit(0.25, "mm"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_markdown(hjust = c(0.5, 0.8)),
    strip.background = element_blank(),
    strip.text.x = element_markdown(colour = "white"),
    plot.margin = unit(c(0.25, 0, 0.25, 0), "cm")
  )


#### compose cov + density ----
png("plots/read-cov-with-density2.png", width = 7, height = 5.5, units = "in", res = 600, type = "quartz")
plotReadCovBases3 + plotReadCovDensitySubgens + plot_layout(widths = c(9,1))
dev.off()

### 
### Plot telomere summaries ----
### 
plotTelEndSummary <- function(g) {
  telEndSummary %>%
    dplyr::filter(genus == g) %>% 
  ggplot() +
  geom_tile(
    aes(
      x = chrom_number %>% as.character() %>% forcats::fct_inseq(), 
      y = forcats::fct_rev(accession), 
      fill = factor(as.numeric(tel_start) + as.numeric(tel_end), levels = c("one" = 1, "both" = 2))
    ),
    color = "#222222",
    size = 0.25
  ) +
  scale_y_discrete(labels = sciNameLabellerSubgen) +
  ggthemes::scale_fill_canva(palette = "Simple but bold", labels = c("one", "both")) +
  guides(fill = guide_legend(), color = "none") +
  labs(x = "chromosome") +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.box.margin = margin(4, 0, 1, 0),
    legend.box.spacing = unit(0, "lines"),
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_markdown(
      size = 6,
      margin = margin(0, 6, 0, 2),
    )
  ) +
  coord_equal()
}
pSpirodelaTelEndSummary <- plotTelEndSummary("Spirodela")
pLemnaTelEndSummary <- plotTelEndSummary("Lemna")
pWolffiaTelEndSummary <- plotTelEndSummary("Wolffia")
pdf("plots/telomere-ends-summary-matrix.pdf", width = 3.5, height = 3.5)
pSpirodelaTelEndSummary / pLemnaTelEndSummary / pWolffiaTelEndSummary + plot_layout(guides = 'collect') + theme(
  legend.position = "bottom",
  legend.box.margin = margin(4, 0, 1, 0),
  legend.box.spacing = unit(0, "lines"),
  legend.direction = "horizontal",
  legend.justification = "center",
  legend.text = element_markdown(
    size = 6,
    margin = margin(0, 6, 0, 2),
  )
)
dev.off()
