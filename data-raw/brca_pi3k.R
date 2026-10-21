# Builds data/brca_pi3k.rda: TCGA BRCA RNA-seq expression of the PI3K
# pathway genes, Normal vs Tumor.
#
# Adapted from dina_simulation/tcga_demo/pancandata.Rmd (M. Janiczek), which
# prepared the TCGA analysis of the DiNA paper. The steps and their order are
# kept exactly, so the result is identical to `pi3k_dat_mat_list` in that
# project's prepared_data_tcga_brca.Rdata.
#
# Steps:
#   1. curatedTCGAData 2.1.1, BRCA, assay RNASeqGene (BRCA_RNASeqGene-20160128,
#      raw counts), via ExperimentHub.
#   2. Tumor = TCGA sample type 01, Normal = 11; the 3 metastatic (06)
#      samples are dropped.
#   3. Genes with count >= 10 in fewer than 20% of samples are dropped.
#   4. DESeq2 variance-stabilising transformation, vst(blind = FALSE), with
#      sample type as the design.
#   5. Per-gene z-score across all Normal + Tumor samples (NOT within group).
#   6. Keep the Normal samples and the Tumor samples of patients without a
#      paired Normal sample (so that the groups are independent).
#   7. Keep the genes of the PI3K pathway in Sanchez-Vega et al. (2018)
#      (data-raw/canonical_genes.xlsx; 28 of its 29 genes are present after
#      step 3 and gene-name matching).
#   8. Keep the first 100 Tumor samples (in the order produced by the
#      pipeline) to balance the 100 Normal samples.
#
# Needs Bioconductor: curatedTCGAData, TCGAutils, MultiAssayExperiment,
# SummarizedExperiment, DESeq2, ExperimentHub; and CRAN: dplyr, tidyr, purrr,
# tibble, readxl.
#
# By default the script uses only the local ExperimentHub cache and fails if
# the BRCA data are not cached. Set BRCA_ALLOW_DOWNLOAD=true to allow a
# download from ExperimentHub instead.
#
# Run from the package root:  Rscript data-raw/brca_pi3k.R
# Optional check against the DiNA object: set DINA_PREPARED to the path of
# prepared_data_tcga_brca.Rdata.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(curatedTCGAData)
  library(TCGAutils)
  library(MultiAssayExperiment)
  library(DESeq2)
})

if (!identical(tolower(Sys.getenv("BRCA_ALLOW_DOWNLOAD")), "true")) {
  ExperimentHub::setExperimentHubOption("LOCAL", TRUE)
}

# 1-2. Data and sample types
brca <- curatedTCGAData(diseaseCode = "BRCA", version = "2.1.1",
                        assays = "RNASeqGene", dry.run = FALSE)
brca_mat <- brca@ExperimentList@listData$`BRCA_RNASeqGene-20160128`@assays@data@listData[[1]]
dat_split <- TCGAutils::TCGAsplitAssays(brca, c("01", "11"))
dat_tum   <- dat_split@ExperimentList@listData$`01_BRCA_RNASeqGene-20160128`@assays@data@listData[[1]]
dat_norm  <- dat_split@ExperimentList@listData$`11_BRCA_RNASeqGene-20160128`@assays@data@listData[[1]]
dat_norm_annot <- data.frame(id_short = colnames(dat_norm)) %>%
  tidyr::separate(id_short, into = c("ID", "Instance"), sep = 12)
sample_type <- ifelse(colnames(brca_mat) %in% colnames(dat_tum),  "Tumor",
               ifelse(colnames(brca_mat) %in% colnames(dat_norm), "Normal",
                                                                  "Other"))

# 3-4. Low-count filter and VST
counts <- round(brca_mat)
mode(counts) <- "integer"
keep <- rowSums(counts >= 10) >= ceiling(0.20 * ncol(counts))
counts <- counts[keep, ]
keep_samples <- sample_type %in% c("Tumor", "Normal")
counts       <- counts[, keep_samples]
sample_type2 <- factor(sample_type[keep_samples], levels = c("Normal", "Tumor"))
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = data.frame(sample_type = sample_type2, row.names = colnames(counts)),
  design    = ~ sample_type
)
vsd       <- vst(dds, blind = FALSE)
brca_norm <- SummarizedExperiment::assay(vsd)

# 5. Per-gene z-score across all samples
standardize <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
brca_std    <- t(apply(brca_norm, 1, standardize))

# 6. Long format with sample annotation
brca_dat <- tibble::as_tibble(t(brca_std), rownames = "ID")
dat_long <- brca_dat %>%
  pivot_longer(cols = -ID, names_to = "gene", values_to = "value",
               values_drop_na = TRUE) %>%
  mutate(sampletype = ifelse(ID %in% colnames(dat_tum),  "Tumor",
                      ifelse(ID %in% colnames(dat_norm), "Normal", "Metastatic"))) %>%
  separate(ID, into = c("ID", "Instance"), sep = 12) %>%
  mutate(hasboth = ifelse(ID %in% dat_norm_annot$ID, "Both", "Tumor only")) %>%
  unite(id, ID, Instance, remove = TRUE, sep = "")

# 7. Canonical pathway genes
canon <- readxl::read_excel("data-raw/canonical_genes.xlsx")
canon$og_tsg[is.na(canon$og_tsg)] <- "Other"
canon[canon == "Unknown"] <- "Other"
canon[canon == "NA"] <- "Other"
dat_long_sub <- dat_long %>%
  arrange(gene) %>%
  filter(gene %in% canon$Gene) %>%
  filter((hasboth == "Both" & sampletype == "Normal") |
         (hasboth == "Tumor only" & sampletype == "Tumor")) %>%
  select(-hasboth)
canon_sub <- canon %>%
  filter(Gene %in% unique(dat_long_sub$gene)) %>%
  group_by(Gene) %>%
  dplyr::slice(1) %>%
  ungroup()
dat_long_sub <- dat_long_sub %>%
  left_join(canon_sub %>% select(gene = Gene, og_tsg), by = "gene") %>%
  arrange(og_tsg)
sample_groups <- sort(unique(dat_long_sub$sampletype))

pi3k_genes <- canon_sub %>% filter(pathway == "PI3K")
pi3k_long <- dat_long_sub %>% filter(gene %in% pi3k_genes$Gene)
brca_pi3k <- purrr::map(sample_groups, ~ pi3k_long %>%
                          filter(sampletype == .x) %>%
                          select(-c(sampletype, og_tsg)) %>%
                          pivot_wider(names_from = gene, values_from = value,
                                      values_fill = 0) %>%
                          select(-id) %>%
                          as.matrix())
names(brca_pi3k) <- sample_groups

# 8. Balance: first 100 Tumor samples
brca_pi3k$Tumor <- brca_pi3k$Tumor[1:100, ]
stopifnot(identical(names(brca_pi3k), c("Normal", "Tumor")),
          all(vapply(brca_pi3k, nrow, integer(1)) == 100),
          identical(colnames(brca_pi3k$Normal), colnames(brca_pi3k$Tumor)),
          !anyNA(unlist(brca_pi3k)))

# Gene annotation, in column order
attr(brca_pi3k, "genes") <- as.data.frame(
  pi3k_genes[match(colnames(brca_pi3k$Normal), pi3k_genes$Gene),
             c("Gene", "og_tsg", "pathway")])

prep <- Sys.getenv("DINA_PREPARED")
if (nzchar(prep)) {
  e <- new.env(); load(prep, envir = e)
  stopifnot(isTRUE(all.equal(unclass(brca_pi3k)[c("Normal", "Tumor")],
                             e$pi3k_dat_mat_list[c("Normal", "Tumor")],
                             tolerance = 1e-12)))
  message("Identical to pi3k_dat_mat_list in ", prep)
}

save(brca_pi3k, file = "data/brca_pi3k.rda", compress = "xz", version = 2)
