#!/usr/bin/env Rscript
# Motif enrichment in differential peaks, with monaLisa + JASPAR.
#
# monaLisa and JASPAR are permissively licensed and redistributable. HOMER is not
# (academic use only), which is why this pipeline does not use it — see
# THIRD_PARTY_LICENSES.md.
#
# What is compared matters as much as the tool. The test here is
#
#     peaks that CHANGED   vs   peaks of the same mark that did NOT change
#
# and not "peaks vs random genome". Against a random genomic background, any
# enhancer set is enriched for every enhancer-associated factor, and the answer
# carries no information. Against the mark's own unchanged peaks, an enrichment
# means the motif distinguishes the responsive elements from the static ones,
# which is the question actually being asked.
#
# Only invoked for targets whose registry entry enables motifs. Broad domains are
# excluded: scanning a 50 kb Polycomb block for 8-mers recovers its base
# composition, and the software will report that with confident p-values.

suppressPackageStartupMessages({
  library(optparse)
  library(monaLisa)
  library(JASPAR2020)
  library(TFBSTools)
  library(GenomicRanges)
  library(rtracklayer)
  library(Rsamtools)
  library(Biostrings)
  library(SummarizedExperiment)
  library(yaml)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--peaks", type = "character", help = "Differential peaks (up or down) BED"),
  make_option("--background", type = "character", help = "Consensus peaks BED"),
  make_option("--genome", type = "character"),
  make_option("--target", type = "character"),
  make_option("--direction", type = "character"),
  make_option("--registry", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--min_peaks", type = "integer", default = 50L)
)))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out_tsv <- file.path(opt$outdir, "motif_enrichment.tsv")

finish_empty <- function(reason) {
  message("[", opt$target, "/", opt$direction, "] ", reason)
  write.table(
    data.frame(motif = character(0), name = character(0), log2enr = numeric(0),
               negLog10P = numeric(0), negLog10Padj = numeric(0)),
    out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no", status = 0)
}

registry <- yaml.load_file(opt$registry)
entry <- registry$marks[[opt$target]]
if (is.null(entry)) entry <- list()
spec <- modifyList(registry$defaults, entry)

read_bed <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(GRanges())
  import(path, format = "BED")
}

fg <- read_bed(opt$peaks)
bg_all <- read_bed(opt$background)

if (length(fg) < opt$min_peaks) {
  finish_empty(sprintf("only %d %s peaks; too few for motif enrichment",
                       length(fg), opt$direction))
}

# Fixed-width windows on the peak centre. Motif enrichment is a comparison of
# sequence composition, so unequal foreground/background sequence lengths would be
# confounded with the very thing being measured.
w <- suppressWarnings(as.integer(spec$motifs$size))
if (is.na(w)) w <- 200L
fg <- resize(fg, width = w, fix = "center")
bg_all <- resize(bg_all, width = w, fix = "center")

# Background = consensus peaks that did NOT change. Leaving the changed peaks in
# the background would dilute precisely the signal being tested for.
bg <- bg_all[!overlapsAny(bg_all, fg)]
message(sprintf("[%s/%s] %d changed vs %d unchanged peaks, %d bp windows",
                opt$target, opt$direction, length(fg), length(bg), w))
if (length(bg) < opt$min_peaks) finish_empty("too few unchanged background peaks")

genome_fa <- opt$genome
if (!file.exists(paste0(genome_fa, ".fai"))) indexFa(genome_fa)
gen <- FaFile(genome_fa)
si <- seqinfo(gen)

# Windows can run off a contig end after resize(); getSeq would error.
constrain <- function(gr) {
  gr <- gr[as.character(seqnames(gr)) %in% seqnames(si)]
  seqlevels(gr) <- seqlevels(si)
  seqinfo(gr) <- si
  gr <- trim(gr)
  gr[width(gr) == w]
}
fg <- constrain(fg)
bg <- constrain(bg)
if (length(fg) < opt$min_peaks) finish_empty("too few foreground peaks after trimming")
if (length(bg) < opt$min_peaks) finish_empty("too few background peaks after trimming")

seqs <- c(getSeq(gen, fg), getSeq(gen, bg))
bins <- factor(c(rep("changed", length(fg)), rep("unchanged", length(bg))),
               levels = c("unchanged", "changed"))

pwms <- getMatrixSet(JASPAR2020, list(matrixtype = "PWM", tax_group = "vertebrates"))
message(sprintf("  scanning %d JASPAR vertebrate motifs", length(pwms)))

# G+C content differs systematically between regulatory element classes. monaLisa
# corrects for composition rather than reporting it as enrichment.
se <- calcBinnedMotifEnrR(seqs = seqs, bins = bins, pwmL = pwms,
                          background = "otherBins", verbose = FALSE)

res <- data.frame(
  motif        = rownames(se),
  name         = rowData(se)$motif.name,
  log2enr      = as.numeric(assay(se, "log2enr")[, "changed"]),
  negLog10P    = as.numeric(assay(se, "negLog10P")[, "changed"]),
  negLog10Padj = as.numeric(assay(se, "negLog10Padj")[, "changed"]),
  stringsAsFactors = FALSE
)
res <- res[order(-res$negLog10Padj, -abs(res$log2enr)), ]
write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

n_sig <- sum(res$negLog10Padj > -log10(0.05), na.rm = TRUE)
message(sprintf("[%s/%s] %d motifs at padj<0.05; top: %s",
                opt$target, opt$direction, n_sig,
                paste(utils::head(res$name, 5), collapse = ", ")))
