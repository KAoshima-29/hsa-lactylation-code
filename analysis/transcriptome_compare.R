#!/usr/bin/env Rscript
# ============================================================
# Transcriptome comparison: HU-HSA cell lines vs public canine cells
# No network calls. Reads ONLY local files provided by the user.
# Data sources:
#   - RNA-seq FASTQs: GEO GSE304507 (download & process separately)
# ============================================================
# How to run (example):
# Rscript analysis/transcriptome_compare.R \
#   --counts_ours data/HU-HSA-cells.tsv \
#   --counts_public data/GeneCounts.tsv \
#   --metadata metadata/Completed_pheno_table.csv \
#   --outdir results
#
# Notes:
# - Column names in the count matrices must match the 'sample' column in metadata after
#   trimming suffix '.genes.results' (the script handles that trimming).
# - Metadata must include: sample, group, study, instrument, lib_prep.
# - Outputs: DE tables (HU_vs_Public, HU_vs_EC), PCA plot, and endothelial marker barplot.
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ---------- CLI ----------
opt <- OptionParser(option_list = list(
  make_option("--counts_ours",   type="character", help="Path to our HU-HSA counts (TSV)"),
  make_option("--counts_public", type="character", help="Path to public counts (TSV)"),
  make_option("--metadata",      type="character", help="Path to metadata CSV"),
  make_option("--outdir",        type="character", default="results", help="Output directory [default: %default]"),
  make_option("--min_cpm",       type="double", default=0.5, help="Minimum CPM threshold for filterByExpr [default: %default]")
))
args <- parse_args(opt)

required <- c("counts_ours","counts_public","metadata")
missing  <- required[!nzchar(trimws(unlist(args[required])))]
if (length(missing)) stop("Missing required option(s): ", paste(missing, collapse=", "))

if (!dir.exists(args$outdir)) dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
figdir <- file.path(args$outdir, "figures"); dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Load data ----------
message("Reading count matrices and metadata...")
ours   <- read.delim(args$counts_ours,   row.names = 1, check.names = FALSE)
public <- read.delim(args$counts_public, row.names = 1, check.names = FALSE)
meta.total <- read.csv(args$metadata, stringsAsFactors = FALSE)

# Intersect genes and combine
common <- intersect(rownames(ours), rownames(public))
if (length(common) < 1000)
  warning("Low gene intersection (", length(common), "). Check that gene IDs use the same reference build.")
counts.total <- cbind(ours[common, , drop=FALSE], public[common, , drop=FALSE])

# Harmonize sample names and align metadata
colnames(counts.total) <- sub("\\.genes\\.results$", "", colnames(counts.total))
meta.total <- meta.total[match(colnames(counts.total), meta.total$sample), ]
stopifnot(all(meta.total$sample == colnames(counts.total)))

# Coerce factors used downstream
meta.total$study      <- factor(meta.total$study)
meta.total$instrument <- factor(meta.total$instrument)
meta.total$lib_prep   <- factor(meta.total$lib_prep)
meta.total$group      <- factor(meta.total$group)
meta.total$study2     <- factor(ifelse(meta.total$study == "HU","HU","Public"),
                                levels = c("Public","HU"))

# ---------- PART A — HU vs Public (all samples) ----------
message("Differential expression: HU vs Public")
dge <- DGEList(counts.total)
keep <- filterByExpr(dge, group = meta.total$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~ study2 + instrument + lib_prep + group, data = meta.total)
colnames(design) <- make.names(colnames(design))
v <- voom(dge, design, plot = FALSE)

fit  <- lmFit(v, design)
contr <- makeContrasts(HU_vs_Public = study2HU, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))

de_pub <- topTable(fit2, coef = "HU_vs_Public", number = Inf)
write.csv(de_pub, file = file.path(args$outdir, "DE_HU_vs_Public.csv"), row.names = TRUE)

# ---------- PART B — HU vs EC (HSA2/3 vs EC_A/B) ----------
message("Differential expression: HU (HSA2/3) vs EC (EC_A/B)")
meta2 <- meta.total
meta2$group2 <- factor(ifelse(meta2$group %in% c("HSA2","HSA3"), "HU",
                              ifelse(meta2$group %in% c("EC_A","EC_B"),  "EC", NA)),
                       levels = c("EC","HU"))
keep.idx <- !is.na(meta2$group2)
meta2 <- droplevels(meta2[keep.idx, ])
counts2 <- counts.total[, keep.idx, drop=FALSE]

dge2 <- DGEList(counts2)
keep2 <- filterByExpr(dge2, group = meta2$group2)
dge2 <- dge2[keep2, , keep.lib.sizes = FALSE]
dge2 <- calcNormFactors(dge2)

design2 <- model.matrix(~ group2, data = meta2)
colnames(design2) <- make.names(colnames(design2))
v2 <- voom(dge2, design2, plot = FALSE)

fitB  <- lmFit(v2, design2)
contr2 <- makeContrasts(HU_vs_EC = group2HU, levels = design2)
fitB2 <- eBayes(contrasts.fit(fitB, contr2))

de_ec <- topTable(fitB2, coef = "HU_vs_EC", number = Inf)
write.csv(de_ec, file = file.path(args$outdir, "DE_HU_vs_EC.csv"), row.names = TRUE)

# ---------- PART C — PCA of all libraries (batch-corrected) ----------
message("PCA (batch-corrected)...")
dge.all <- DGEList(counts.total)
keep.all <- filterByExpr(dge.all, group = meta.total$group)
dge.all <- dge.all[keep.all, , keep.lib.sizes = FALSE]
dge.all <- calcNormFactors(dge.all)

# Technical-only design for voom; biological structure retained for plotting
design.pca <- model.matrix(~ study + instrument + lib_prep, data = meta.total)
v.all <- voom(dge.all, design.pca, plot = FALSE)

expr.bc <- removeBatchEffect(
  v.all$E,
  batch  = meta.total$study,
  batch2 = meta.total$instrument,
  design = model.matrix(~ group, data = meta.total)
)

pca <- prcomp(t(expr.bc), scale. = FALSE)
pca.df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                     group = meta.total$group,
                     sample = meta.total$sample)

# Palette/ordering (adjust if needed)
grp.order <- c("HSA2","HSA3","EC_A","EC_B",
               "Fibroblast_FBS","MDCK",
               "Osteosarcoma_A","Osteosarcoma_B",
               "Prostate_cancer_A","Prostate_cancer_B")
pca.df$group <- factor(pca.df$group, levels = grp.order)

p <- ggplot(pca.df, aes(PC1, PC2, colour = group, label = sample)) +
  geom_point(size = 3) +
  labs(title = "PCA (batch-corrected)", x = "PC1", y = "PC2") +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")
ggsave(file.path(args$outdir, "PCA_all_groups.png"), p, width = 8, height = 6, dpi = 300)

# ---------- PART D — Endothelial markers (mean ± SD logCPM) ----------
message("Summarising endothelial markers...")
markers <- c(
  KDR    = "ENSCAFG00845012083",
  PECAM1 = "ENSCAFG00845022888",
  ERG    = "ENSCAFG00845029042",
  CDH5   = "ENSCAFG00845012097",
  VWF    = "ENSCAFG00845029798",
  CLDN5  = "ENSCAFG00845029410",
  TEK    = "ENSCAFG00845018055",
  ENG    = "ENSCAFG00845005326",
  CD34   = "ENSCAFG00845002498"
)

expr.mat <- v.all$E
present <- markers %in% rownames(expr.mat)
if (!all(present)) warning("Missing markers in matrix: ",
                           paste(names(markers)[!present], collapse = ", "))
expr.sub <- expr.mat[markers[present], , drop = FALSE]

long <- as.data.frame(expr.sub) |>
  tibble::rownames_to_column("GeneID") |>
  tidyr::pivot_longer(-GeneID, names_to = "sample", values_to = "logCPM") |>
  dplyr::left_join(meta.total[, c("sample","group")], by = "sample") |>
  dplyr::mutate(Gene = names(markers)[match(GeneID, markers)])

summ <- long |>
  dplyr::group_by(Gene, group) |>
  dplyr::summarise(mean = mean(logCPM), sd = sd(logCPM), .groups = "drop") |>
  dplyr::mutate(group = factor(group, levels = grp.order))

bp <- ggplot(summ, aes(Gene, mean, fill = group)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.8), width = 0.2) +
  labs(title = "Endothelial markers (mean ± SD logCPM)", x = NULL, y = "logCPM") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(args$outdir, "Endothelial_markers_barplot.png"), bp, width = 9, height = 6, dpi = 300)

# ---------- Session info ----------
writeLines(capture.output(sessionInfo()), con = file.path(args$outdir, "sessionInfo.txt"))
message("Done. Outputs written to: ", normalizePath(args$outdir))
