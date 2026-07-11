#!/usr/bin/env Rscript
# Peak annotation and GO enrichment.
#
# Gene assignment depends on the peak mode:
#
#   narrow  A peak is a discrete regulatory element, so the nearest TSS is a
#           defensible guess at its target gene. ChIPseeker's annotatePeak does this.
#
#   broad   A domain can span hundreds of kb and dozens of genes. "Nearest gene" for
#           a 400 kb H3K27me3 block is not wrong so much as meaningless. For broad
#           marks the genes reported are all genes the domain overlaps, which is the
#           question actually being asked of a repressive domain: what is inside it.

suppressPackageStartupMessages({
  library(optparse)
  library(ChIPseeker)
  library(GenomicRanges)
  library(clusterProfiler)
  library(ggplot2)
  library(yaml)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--consensus", type = "character"),
  make_option("--up", type = "character"),
  make_option("--down", type = "character"),
  make_option("--target", type = "character"),
  make_option("--registry", type = "character"),
  make_option("--txdb", type = "character"),
  make_option("--orgdb", type = "character"),
  make_option("--outdir", type = "character")
)))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(opt$txdb, character.only = TRUE)
  library(opt$orgdb, character.only = TRUE)
})
txdb <- get(opt$txdb)
orgdb <- opt$orgdb

registry <- yaml.load_file(opt$registry)
entry <- registry$marks[[opt$target]]
if (is.null(entry)) entry <- list()
spec <- modifyList(registry$defaults, entry)
peak_mode <- spec$peak_mode
is_broad <- identical(peak_mode, "broad")

message(sprintf("[%s] peak_mode=%s -> gene assignment by %s",
                opt$target, peak_mode,
                if (is_broad) "domain overlap" else "nearest TSS"))

read_bed <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(GRanges())
  df <- read.delim(path, header = FALSE)
  colnames(df)[1:3] <- c("chr", "start", "end")
  GRanges(df$chr, IRanges(df$start + 1, df$end),
          score = if (ncol(df) >= 5) df[[5]] else 0)
}

consensus <- read_bed(opt$consensus)
message(sprintf("  %d consensus peaks", length(consensus)))

# ------------------------------------------------------- feature distribution
#
# Where a mark sits relative to gene structure is the first sanity check of the
# whole experiment. H3K4me3 that is not overwhelmingly promoter-proximal means
# the antibody or the analysis is wrong, and no downstream result can be trusted.

anno <- annotatePeak(consensus, TxDb = txdb, annoDb = orgdb,
                     tssRegion = c(-3000, 3000), verbose = FALSE)
anno_df <- as.data.frame(anno)
write.table(anno_df, file.path(opt$outdir, "peak_annotation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(opt$outdir, "feature_distribution.pdf"), width = 9, height = 4)
print(plotAnnoBar(anno, title = sprintf("%s (%s peaks): genomic distribution",
                                        opt$target, peak_mode)))
print(plotDistToTSS(anno, title = sprintf("%s: distance to TSS", opt$target)))
dev.off()

# ------------------------------------------------------------ gene assignment

genes_for <- function(gr) {
  if (length(gr) == 0) return(character(0))
  if (is_broad) {
    # Every gene whose TSS falls inside the domain.
    tss <- GenomicFeatures::promoters(txdb, upstream = 0, downstream = 1)
    hits <- findOverlaps(tss, gr)
    ids <- unique(names(tss)[queryHits(hits)])
    if (is.null(ids) || length(ids) == 0) {
      ids <- unique(as.character(mcols(tss)$tx_name[queryHits(hits)]))
    }
    suppressMessages(
      sym <- AnnotationDbi::select(get(orgdb), keys = unique(ids),
                                   keytype = "ENTREZID", columns = "SYMBOL")$SYMBOL
    )
    return(unique(na.omit(sym)))
  }
  a <- as.data.frame(annotatePeak(gr, TxDb = txdb, annoDb = orgdb,
                                  tssRegion = c(-3000, 3000), verbose = FALSE))
  unique(na.omit(a$SYMBOL))
}

up <- read_bed(opt$up)
down <- read_bed(opt$down)
message(sprintf("  %d up, %d down differential peaks", length(up), length(down)))

gene_rows <- list()
for (dir in c("up", "down")) {
  gr <- if (dir == "up") up else down
  g <- genes_for(gr)
  if (length(g)) {
    gene_rows[[dir]] <- data.frame(direction = dir, gene = g)
  }
}
gene_tab <- if (length(gene_rows)) do.call(rbind, gene_rows) else
  data.frame(direction = character(0), gene = character(0))
write.table(gene_tab, file.path(opt$outdir, "differential_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------- GO enrichment
#
# The universe is the genes near the consensus peaks, not the whole genome.
# Testing differential H3K27ac genes against all ~20k genes would mostly
# rediscover that H3K27ac marks expressed genes, which we already knew.

universe <- unique(na.omit(anno_df$SYMBOL))
message(sprintf("  GO universe: %d genes near consensus peaks", length(universe)))

go_all <- list()
for (dir in unique(gene_tab$direction)) {
  g <- gene_tab$gene[gene_tab$direction == dir]
  if (length(g) < 10) {
    message(sprintf("  %s: only %d genes; skipping GO", dir, length(g)))
    next
  }
  ego <- tryCatch(
    enrichGO(gene = g, universe = universe, OrgDb = orgdb, keyType = "SYMBOL",
             ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.1,
             readable = FALSE),
    error = function(e) { message("  enrichGO failed: ", conditionMessage(e)); NULL }
  )
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    df <- as.data.frame(ego)
    df$direction <- dir
    go_all[[dir]] <- df
  }
}

go_tab <- if (length(go_all)) do.call(rbind, go_all) else
  data.frame(ID = character(0), Description = character(0), p.adjust = numeric(0),
             direction = character(0))
write.table(go_tab, file.path(opt$outdir, "go_enrichment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("[%s] %d enriched GO terms", opt$target, nrow(go_tab)))
