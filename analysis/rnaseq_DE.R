args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript rnaseq_DE.R sample_metadata.tsv output_dir result.tsv plot.pdf")

metadata_file <- args[1]
outdir <- args[2]
result_file <- args[3]
plot_file <- args[4]

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(pheatmap)

samples <- read_tsv(metadata_file, col_types = cols())
files <- file.path("results/rnaseq/salmon", samples$SampleID, "quant.sf")
names(files) <- samples$SampleID

txi <- tximport(files, type = "salmon", txOut = FALSE)
coldata <- samples %>% select(SampleID, Condition, Replicate) %>% mutate(Condition = factor(Condition))
rownames(coldata) <- coldata$SampleID

dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ Condition)
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
res_df <- as.data.frame(res) %>% rownames_to_column("gene") %>% arrange(padj)

write_tsv(res_df, result_file)

pdf(plot_file, width = 10, height = 8)
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "Condition")
plotMA(res, main = "RNA-seq DESeq2", ylim = c(-5, 5))
hist(res$pvalue, breaks = 50, main = "P-value distribution", xlab = "p-value")
dev.off()
