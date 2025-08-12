#!/usr/bin/env Rscript
# ============================================================
# edgeR two-group differential expression (raw counts, not TPM)
# - Input counts: genes x samples TSV (first column = gene_id)
# - Metadata CSV: must contain at least: sample, group
# - Robust GLM (QL) workflow with filterByExpr + TMM
# - Outputs: DE table (case vs control), TMM-CPM matrix, sessionInfo
# Data source (RNA-seq): GEO GSE304507 (download & process separately)
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(edgeR)
})

# ---------- CLI ----------
opt <- OptionParser(option_list = list(
  make_option("--counts",        type="character", help="TSV counts matrix (genes x samples). First column must be gene_id."),
  make_option("--metadata",      type="character", help="CSV metadata with columns: sample, group (and optionally others)."),
  make_option("--group_col",     type="character", default="group", help="Metadata column name for the condition [default: %default]"),
  make_option("--control_level", type="character", help="Control level in group factor (e.g., '25mM')"),
  make_option("--case_level",    type="character", help="Case level in group factor (e.g., '0mM')"),
  make_option("--min_count",     type="integer",   default=10, help="Minimum total count for filterByExpr (soft threshold) [default: %default]"),
  make_option("--outdir",        type="character", default="results/edgeR", help="Output directory [default: %default]")
))
args <- parse_args(opt)
req <- c("counts","metadata","group_col","outdir")
if (any(is.na(args[req]))) stop("Missing required arguments. See --help.")
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# ---------- IO ----------
message("Reading counts and metadata...")
counts_in <- readr::read_tsv(args$counts, col_types = cols())
stopifnot(ncol(counts_in) >= 2)
if (names(counts_in)[1] %in% c("gene","gene_id","GeneID","Gene")) {
  names(counts_in)[1] <- "gene_id"
}
stopifnot(names(counts_in)[1] == "gene_id")
mat <- as.data.frame(counts_in)
rownames(mat) <- mat$gene_id; mat$gene_id <- NULL

meta <- read.csv(args$metadata, stringsAsFactors = FALSE)
stopifnot(all(c("sample", args$group_col) %in% names(meta)))

# keep only intersecting samples, and order to metadata
common <- intersect(colnames(mat), meta$sample)
if (length(common) < 2) stop("Fewer than 2 overlapping samples between counts and metadata.")
meta <- meta[match(common, meta$sample), , drop = FALSE]
mat  <- mat[, common, drop = FALSE]

# ------- basic checks -------
# ensure numeric
if (!all(sapply(mat, is.numeric))) stop("Counts matrix must be numeric (raw counts).")
# aggregate duplicated gene IDs if any
if (any(duplicated(rownames(mat)))) {
  warning("Duplicated gene IDs found. Aggregating by sum.")
  mat <- rowsum(mat, group = rownames(mat))
}

# ------- define groups (two levels only) -------
grp <- factor(meta[[args$group_col]])
# if case/control are provided, subset to those two; else take first two levels
if (!is.null(args$control_level) && !is.null(args$case_level)) {
  keep_idx <- grp %in% c(args$control_level, args$case_level)
  meta <- meta[keep_idx, , drop = FALSE]
  mat  <- mat[, keep_idx, drop = FALSE]
  grp  <- factor(meta[[args$group_col]], levels = c(args$control_level, args$case_level))
} else {
  levs <- levels(grp)
  if (length(levs) != 2) stop("Group has ", length(levs), " levels. Provide --control_level and --case_level to select two.")
  grp <- factor(grp, levels = levs) # control=levs[1], case=levs[2]
}

# edgeR pipeline
y <- DGEList(counts = mat, group = grp)

# filtering
keep <- filterByExpr(y, group = grp)
y <- y[keep, , keep.lib.sizes = FALSE]

# normalization
y <- calcNormFactors(y, method = "TMM")

# design: ~0+group for explicit contrast (case - control)
design <- model.matrix(~ 0 + grp)
colnames(design) <- sub("^grp", "", colnames(design))  # names become control/case
if (ncol(design) != 2) stop("Design must have exactly 2 columns for two-group contrast.")

# dispersion + QL fit
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

# contrast: case - control
cn <- colnames(design)
control <- cn[1]; case <- cn[2]
contrast_vec <- makeContrasts(contrasts = paste0(case, " - ", control), levels = design)

qlf <- glmQLFTest(fit, contrast = contrast_vec)

# outputs
tt <- topTags(qlf, n = Inf)$table
tt <- tibble::rownames_to_column(tt, var = "gene_id")

de_out <- file.path(args$outdir, paste0("edgeR_DE_", case, "_vs_", control, ".tsv"))
readr::write_tsv(tt, de_out)

cpm_mat <- cpm(y, normalized.lib.sizes = TRUE)
cpm_out <- file.path(args$outdir, "edgeR_TMM_CPM.tsv")
readr::write_tsv(tibble::as_tibble(cpm_mat, rownames = "gene_id"), cpm_out)

libs <- data.frame(
  sample = colnames(y),
  libsize_raw = y$samples$lib.size,
  norm_factor = y$samples$norm.factors,
  stringsAsFactors = FALSE
)
readr::write_tsv(libs, file.path(args$outdir, "edgeR_library_sizes.tsv"))

writeLines(capture.output(sessionInfo()), con = file.path(args$outdir, "sessionInfo.txt"))
message("Done. Wrote:")
message("  - ", de_out)
message("  - ", cpm_out)
message("  - ", file.path(args$outdir, "edgeR_library_sizes.tsv"))
message("  - ", file.path(args$outdir, "sessionInfo.txt"))
