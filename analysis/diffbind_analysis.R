args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
    stop("Usage: Rscript diffbind_analysis.R sample_sheet.csv output_dir factor numerator reference narrow_summits")
}

sample_sheet <- args[1]
outdir <- args[2]
factor_name <- args[3]
numerator_condition <- args[4]
reference_condition <- args[5]
narrow_summits <- suppressWarnings(as.integer(args[6]))
if (is.na(narrow_summits) || narrow_summits < 1) stop("narrow_summits must be a positive integer")
if (!nzchar(numerator_condition) || !nzchar(reference_condition) || numerator_condition == reference_condition) {
    stop("Numerator and reference conditions must be distinct and non-empty")
}

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

library(DiffBind)

samples <- read.csv(sample_sheet, stringsAsFactors = FALSE)
required <- c("SampleID", "Factor", "Condition", "Replicate", "bamReads", "bamControl", "Peaks", "PeakCaller", "PeakMode")
missing <- setdiff(required, colnames(samples))
if (length(missing) > 0) {
    stop(paste("DiffBind sample sheet is missing columns:", paste(missing, collapse = ", ")))
}

samples <- samples[samples$Factor == factor_name, , drop = FALSE]
if (nrow(samples) == 0) {
    stop(paste("No samples found for factor", factor_name))
}

conditions <- unique(samples$Condition)
if (length(conditions) != 2) {
    stop("Differential binding requires exactly two conditions for one unambiguous contrast")
}
if (!setequal(conditions, c(reference_condition, numerator_condition))) {
    stop("Sample-sheet conditions do not match the configured numerator/reference contrast")
}

if (length(unique(samples$Factor)) != 1) {
    stop("Each DiffBind run must contain exactly one factor")
}

replicate_counts <- table(samples$Condition)
if (any(replicate_counts < 2)) {
    stop("Differential binding requires at least two ChIP-seq replicates per condition")
}

peak_modes <- unique(tolower(trimws(samples$PeakMode)))
if (length(peak_modes) != 1 || !peak_modes %in% c("narrow", "broad")) {
    stop("Each factor must contain one valid PeakMode: narrow or broad")
}

dba_samples <- samples[, setdiff(colnames(samples), "PeakMode"), drop = FALSE]
dba_obj <- dba(sampleSheet = dba_samples)
if (peak_modes == "broad") {
    dba_obj <- dba.count(dba_obj, minOverlap = 2, summits = FALSE)
} else {
    dba_obj <- dba.count(dba_obj, minOverlap = 2, summits = narrow_summits)
}
dba_obj <- dba.contrast(
    dba_obj,
    contrast = c("Condition", numerator_condition, reference_condition),
    minMembers = 2
)
dba_obj <- dba.analyze(dba_obj)
report <- dba.report(dba_obj, th = 0.05)

report_df <- as.data.frame(report)
report_df$ContrastNumerator <- rep(numerator_condition, nrow(report_df))
report_df$ContrastReference <- rep(reference_condition, nrow(report_df))
write.csv(report_df, file.path(outdir, "diffbind_summary.csv"), row.names = FALSE)

contrast_metadata <- data.frame(
    Factor = factor_name,
    Numerator = numerator_condition,
    Reference = reference_condition,
    PeakMode = peak_modes,
    NarrowSummits = if (peak_modes == "narrow") narrow_summits else NA_integer_
)
write.table(
    contrast_metadata,
    file.path(outdir, "contrast.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

pdf(file.path(outdir, "diffbind_plots.pdf"), width = 10, height = 8)
plot(dba_obj)
dba.plotPCA(dba_obj, attributes = DBA_CONDITION, label = DBA_ID)
dba.plotHeatmap(dba_obj, contrast = 1, correlations = FALSE, th = 1, maxSites = 1000)
dba.plotMA(dba_obj, contrast = 1)
dev.off()

saveRDS(dba_obj, file.path(outdir, "diffbind.rds"))
