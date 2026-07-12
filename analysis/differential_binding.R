#!/usr/bin/env Rscript
# Differential binding with DESeq2, normalised according to the peak mode.
#
# The normalisation choice is the substantive decision here.
#
# Narrow marks (TF, H3K27ac, H3K4me3): size factors from the peak counts by
# median-of-ratios. This assumes most peaks do not change between conditions,
# which for a punctate mark under a physiological perturbation is reasonable —
# a few thousand enhancers move, tens of thousands do not.
#
# Broad marks (H3K27me3, H3K9me3, H3K36me3): size factors from genome-wide
# background bins instead. The "most regions unchanged" assumption fails for
# repressive domains, which can shift globally; the textbook case is EZH2
# inhibition, where H3K27me3 drops almost everywhere. Median-of-ratios on peak
# counts would absorb that global loss into the size factors and report no
# change, which is precisely backwards. Estimating depth from background bins —
# regions that are mostly free of the mark — keeps the global component in the
# result where it belongs. csaw (Lun & Smyth, NAR 2016) makes the same argument.
#
# Where a genuinely global shift is expected, background normalisation is still a
# proxy; an exogenous spike-in is the only unbiased answer. Flagged in the summary.

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(yaml)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts", type = "character"),
  make_option("--background", type = "character"),
  make_option("--coldata", type = "character"),
  make_option("--regions", type = "character"),
  make_option("--target", type = "character"),
  make_option("--registry", type = "character"),
  make_option("--reference", type = "character"),
  make_option("--treatment", type = "character"),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--lfc", type = "double", default = 1.0),
  make_option("--outdir", type = "character")
)))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

registry <- yaml.load_file(opt$registry)
entry <- registry$marks[[opt$target]]
if (is.null(entry)) {
  message("[", opt$target, "] not in the mark registry; using defaults")
  entry <- list()
}
spec <- modifyList(registry$defaults, entry)

peak_mode <- spec$peak_mode
norm_method <- spec$diffbind$normalize

message(sprintf("[%s] peak_mode=%s  normalisation=%s", opt$target, peak_mode, norm_method))

# ------------------------------------------------------------------ load counts

counts <- read.delim(opt$counts, check.names = FALSE)

# A target whose ChIP did not enrich arrives here with no regions. There is
# nothing to test, but the run must not die: the other targets are unaffected and
# their results are still wanted. Write the full set of empty outputs so the DAG
# completes, and let the QC gate and the summary report why this target is empty.
if (nrow(counts) == 0) {
  message(sprintf(
    "[%s] no consensus peaks - the ChIP did not enrich. Writing empty results.",
    opt$target))
  empty <- data.frame(
    region = character(0), chrom = character(0), start = integer(0),
    end = integer(0), baseMean = numeric(0), log2FoldChange = numeric(0),
    lfcSE = numeric(0), pvalue = numeric(0), padj = numeric(0),
    direction = character(0))
  write.table(empty, file.path(opt$outdir, "results.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  for (d in c("up", "down")) {
    cat("", file = file.path(opt$outdir, paste0(d, ".bed")))
  }
  for (f in c("normalized_counts.tsv", "rlog.tsv", "size_factors.tsv")) {
    cat("", file = file.path(opt$outdir, f))
  }
  quit(save = "no", status = 0)
}

rownames(counts) <- counts$region
counts$region <- NULL

coldata <- read.csv(opt$coldata)
rownames(coldata) <- coldata$sample_id
coldata <- coldata[colnames(counts), , drop = FALSE]

if (any(is.na(coldata$sample_id))) {
  stop("Count matrix columns do not match the sample sheet")
}

coldata$condition <- factor(coldata$condition, levels = c(opt$reference, opt$treatment))
coldata$replicate <- factor(coldata$replicate)
coldata$layout <- factor(coldata$layout)

message(sprintf("  %d regions x %d samples", nrow(counts), ncol(counts)))
print(table(coldata$condition, coldata$layout))

# Drop regions with essentially no coverage anywhere; they only inflate the
# multiple-testing burden.
keep <- rowSums(counts) >= 10
message(sprintf("  %d regions retained (>=10 reads total)", sum(keep)))
counts <- counts[keep, , drop = FALSE]

# ------------------------------------------------------------------ the design
#
# Replicate 1 is single-end and replicate 2 is paired-end in this dataset, and
# replicate is perfectly confounded with layout. Both are present in both
# conditions, so ~ replicate + condition is estimable and removes the batch
# component (which here is library layout) from the condition effect. Modelling
# layout separately is impossible — it is the same variable.

design <- if (nlevels(coldata$replicate) > 1) ~ replicate + condition else ~ condition
message("  design: ", paste(deparse(design), collapse = ""))

dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), coldata, design)

# ----------------------------------------------------- mark-aware size factors

if (identical(norm_method, "lib")) {
  bg <- read.delim(opt$background, check.names = FALSE)
  rownames(bg) <- bg$region
  bg$region <- NULL
  bg <- bg[, colnames(counts), drop = FALSE]
  # Bins with no reads in some library carry no depth information.
  bg <- bg[rowSums(bg == 0) == 0, , drop = FALSE]
  message(sprintf("  size factors from %d background bins (broad mark)", nrow(bg)))
  sf <- estimateSizeFactorsForMatrix(as.matrix(bg))
  sizeFactors(dds) <- sf[colnames(dds)]
} else {
  message("  size factors from peak counts, median-of-ratios (narrow mark)")
  dds <- estimateSizeFactors(dds)
}

sf_out <- data.frame(sample = names(sizeFactors(dds)),
                     size_factor = as.numeric(sizeFactors(dds)),
                     method = norm_method)
write.table(sf_out, file.path(opt$outdir, "size_factors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(sf_out)

# ------------------------------------------------------------------- the test

dds <- DESeq(dds, quiet = TRUE)

res <- results(dds,
               contrast = c("condition", opt$treatment, opt$reference),
               alpha = opt$fdr)

# Shrink the log2 fold changes. Low-count peaks otherwise produce enormous,
# meaningless fold changes that dominate every volcano plot.
res <- tryCatch(
  lfcShrink(dds, coef = paste0("condition_", opt$treatment, "_vs_", opt$reference),
            res = res, type = "apeglm", quiet = TRUE),
  error = function(e) {
    message("  apeglm unavailable (", conditionMessage(e), "); falling back to ashr")
    tryCatch(
      lfcShrink(dds, contrast = c("condition", opt$treatment, opt$reference),
                res = res, type = "ashr", quiet = TRUE),
      error = function(e2) {
        message("  shrinkage unavailable; reporting unshrunk MLE fold changes")
        res
      })
  })

regions <- read.delim(opt$regions, header = FALSE,
                      col.names = c("chr", "start", "end", "region", "score", "strand"))
rownames(regions) <- regions$region

out <- as.data.frame(res)
out$region <- rownames(out)
out <- merge(regions[, c("region", "chr", "start", "end")], out, by = "region")

# The call needs an effect size as well as an FDR: with enough depth, ChIP-seq
# will return a significant q-value for a 1.1x change that means nothing.
out$significant <- !is.na(out$padj) & out$padj < opt$fdr & abs(out$log2FoldChange) >= opt$lfc
out$direction <- ifelse(!out$significant, "ns",
                 ifelse(out$log2FoldChange > 0, "up", "down"))

out <- out[order(out$padj, -abs(out$log2FoldChange)), ]
write.table(out, file.path(opt$outdir, "results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

bed_cols <- c("chr", "start", "end", "region", "log2FoldChange", "strand")
for (dir in c("up", "down")) {
  sel <- out[out$direction == dir, ]
  # rep(".", nrow(sel)), not ".": a scalar recycles fine against non-empty
  # columns but has length 1 against zero-length ones, and data.frame() then
  # refuses to build ("differing number of rows: 0, 1"). A direction with no
  # significant peaks is an ordinary result, not an error.
  bed <- data.frame(chr = sel$chr, start = sel$start, end = sel$end,
                    region = sel$region, score = round(sel$log2FoldChange, 3),
                    strand = rep(".", nrow(sel)))
  bed <- bed[order(bed$chr, bed$start), , drop = FALSE]
  write.table(bed, file.path(opt$outdir, paste0(dir, ".bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# ------------------------------------------------------- matrices for plotting

norm_counts <- counts(dds, normalized = TRUE)
write.table(data.frame(region = rownames(norm_counts), norm_counts, check.names = FALSE),
            file.path(opt$outdir, "normalized_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# rlog is better behaved than VST on the small n typical of ChIP-seq, but it is
# slow on large matrices; VST above 30k regions keeps this tractable.
trans <- if (nrow(dds) > 30000) vst(dds, blind = FALSE) else rlog(dds, blind = FALSE)
mat <- assay(trans)
write.table(data.frame(region = rownames(mat), mat, check.names = FALSE),
            file.path(opt$outdir, "rlog.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

n_up <- sum(out$direction == "up")
n_dn <- sum(out$direction == "down")
message(sprintf("[%s] %d up, %d down (FDR<%.2f, |log2FC|>=%.1f) of %d tested",
                opt$target, n_up, n_dn, opt$fdr, opt$lfc, nrow(out)))
