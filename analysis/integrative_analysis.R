args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript integrative_analysis.R diffbind.csv deseq2.tsv annotation.gtf outdir")

diffbind_file <- args[1]
deseq_file <- args[2]
annotation_gtf <- args[3]
outdir <- args[4]

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(ggplot2)
library(ChIPseeker)

peaks <- read.csv(diffbind_file, stringsAsFactors = FALSE)
if (!all(c("Chr", "Start", "End") %in% colnames(peaks))) {
    stop("DiffBind summary must contain Chr, Start, End columns")
}

peak_gr <- GRanges(seqnames = peaks$Chr, ranges = IRanges(peaks$Start, peaks$End))
txdb <- makeTxDbFromGFF(annotation_gtf, format = "gtf")
anno <- annotatePeak(peak_gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
anno_df <- as.data.frame(anno)

if (!"geneId" %in% colnames(anno_df)) {
    stop("Annotation did not produce geneId values")
}

deseq <- read.delim(deseq_file, stringsAsFactors = FALSE)

merged <- anno_df %>%
    select(seqnames, start, end, geneId, annotation, distanceToTSS) %>%
    rename(GeneID = geneId, PeakChr = seqnames, PeakStart = start, PeakEnd = end) %>%
    left_join(deseq %>% rename(GeneID = gene), by = "GeneID")

summary_file <- file.path(outdir, "integrative_summary.csv")
write.csv(merged, summary_file, row.names = FALSE)

pdf(file.path(outdir, "integrative_scatter.pdf"), width = 10, height = 8)
if (all(c("Fold", "log2FoldChange") %in% colnames(merged))) {
    ggplot(merged, aes(x = log2FoldChange, y = Fold)) +
        geom_point(alpha = 0.5) +
        theme_minimal() +
        xlab("RNA log2FC") +
        ylab("ChIP log2FC") +
        ggtitle("Integrated gene expression and differential binding")
} else {
    plot(1, 1, type = "n", xlab = "RNA log2FC", ylab = "ChIP log2FC", main = "Integration plot requires Fold and log2FoldChange columns")
}
dev.off()
