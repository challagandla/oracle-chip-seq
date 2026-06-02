args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript diffbind_analysis.R sample_sheet.csv output_dir")

sample_sheet <- args[1]
outdir <- args[2]

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

library(DiffBind)

samples <- read.csv(sample_sheet, stringsAsFactors = FALSE)
dba_obj <- dba(sampleSheet = samples)
dba_obj <- dba.count(dba_obj, minOverlap = 2)
dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION)
dba_obj <- dba.analyze(dba_obj)
report <- dba.report(dba_obj, th = 0.05)

write.csv(as.data.frame(report), file.path(outdir, "diffbind_summary.csv"), row.names = FALSE)

pdf(file.path(outdir, "diffbind_plots.pdf"), width = 10, height = 8)
plot(dba_obj)
dba.plotPCA(dba_obj, DBA_TREATMENT, label = DBA_ID)
dba.plotHeatmap(dba_obj, contrast = 1, correlations = FALSE)
dba.plotMA(dba_obj, contrast = 1)
dev.off()

saveRDS(dba_obj, file.path(outdir, "diffbind.rds"))
