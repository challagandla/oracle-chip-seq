#!/usr/bin/env Rscript
# Open-source motif enrichment for narrow-factor ChIP-seq consensus peaks.
#
# The script tests JASPAR vertebrate CORE motifs in fixed-width foreground
# windows against a deterministic genomic background. Background windows come
# from the same chromosome distribution, do not overlap each other, and exclude
# blacklisted regions and the complete foreground consensus intervals.
#
# Usage: Rscript motif_enrichment.R peaks.bed genome.fa blacklist.bed outdir
#        window_bp background_multiplier seed [factor]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
    stop(
        "Usage: Rscript motif_enrichment.R peaks.bed genome.fa blacklist.bed ",
        "outdir window_bp background_multiplier seed [factor]"
    )
}

peaks_bed <- args[1]
genome_fa <- args[2]
blacklist_bed <- args[3]
outdir <- args[4]
window_bp <- suppressWarnings(as.integer(args[5]))
background_multiplier <- suppressWarnings(as.integer(args[6]))
seed <- suppressWarnings(as.integer(args[7]))
factor_label <- if (length(args) >= 8) args[8] else basename(outdir)

if (is.na(window_bp) || window_bp < 1L) stop("window_bp must be positive")
if (is.na(background_multiplier) || background_multiplier < 1L) {
    stop("background_multiplier must be positive")
}
if (is.na(seed) || seed < 0L) stop("seed must be non-negative")
for (path in c(peaks_bed, genome_fa, blacklist_bed)) {
    if (!file.exists(path)) stop("Required input does not exist: ", path)
}
if (!file.exists(paste0(genome_fa, ".fai"))) {
    stop("Genome FASTA index does not exist: ", genome_fa, ".fai")
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

suppressPackageStartupMessages({
    library(monaLisa)
    library(JASPAR2020)
    library(TFBSTools)
    library(GenomicRanges)
    library(rtracklayer)
    library(Rsamtools)
    library(Biostrings)
    library(SummarizedExperiment)
    library(BiocParallel)
    library(readr)
    library(dplyr)
    library(ggplot2)
})

max_fraction_n <- 0.1
gc_breaks <- c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8)
genome <- FaFile(genome_fa)
genome_info <- seqinfo(genome)

# Preserve the full consensus intervals for exclusion from the background, but
# use exact-width midpoint windows for motif testing. Consensus BEDs do not
# retain MACS summit offsets, so these midpoints must not be called summits.
consensus_regions <- import(peaks_bed)
if (length(consensus_regions) == 0L) stop("No peaks found in ", peaks_bed)
common_contigs <- intersect(seqlevels(consensus_regions), seqlevels(genome_info))
if (length(common_contigs) == 0L) {
    stop(
        "No overlapping sequence names between peaks and genome FASTA. ",
        "Check that both use the same assembly and chromosome naming."
    )
}

consensus_regions <- keepSeqlevels(
    consensus_regions, common_contigs, pruning.mode = "coarse"
)
seqlengths(consensus_regions) <- seqlengths(genome_info)[common_contigs]
strand(consensus_regions) <- "*"
consensus_regions <- trim(consensus_regions)
foreground_regions <- trim(
    resize(consensus_regions, width = window_bp, fix = "center")
)
foreground_regions <- unique(foreground_regions[width(foreground_regions) == window_bp])
if (length(foreground_regions) == 0L) {
    stop("No full-width foreground windows remain after trimming to genome bounds")
}

foreground_sequences <- scanFa(genome, param = foreground_regions)
foreground_n <- letterFrequency(
    foreground_sequences, letters = "N", as.prob = TRUE
)[, 1]
keep_foreground <- foreground_n <= max_fraction_n
foreground_regions <- foreground_regions[keep_foreground]
foreground_sequences <- foreground_sequences[keep_foreground]
if (length(foreground_sequences) == 0L) {
    stop("No foreground windows remain after filtering ambiguous sequence")
}

blacklist_regions <- GRanges()
blacklist_size <- file.info(blacklist_bed)$size
if (!is.na(blacklist_size) && blacklist_size > 0L) {
    blacklist_regions <- import(blacklist_bed)
    blacklist_contigs <- intersect(seqlevels(blacklist_regions), common_contigs)
    if (length(blacklist_contigs) > 0L) {
        blacklist_regions <- keepSeqlevels(
            blacklist_regions, blacklist_contigs, pruning.mode = "coarse"
        )
        seqlengths(blacklist_regions) <- seqlengths(genome_info)[blacklist_contigs]
        strand(blacklist_regions) <- "*"
        blacklist_regions <- trim(blacklist_regions)
    } else {
        blacklist_regions <- GRanges()
    }
}

# Draw non-overlapping windows from disjoint slots. Sampling is stratified by
# the observed foreground chromosome counts, preserving that distribution.
set.seed(seed)
foreground_chromosomes <- as.character(seqnames(foreground_regions))
chromosomes <- unique(foreground_chromosomes)
foreground_counts <- table(factor(foreground_chromosomes, levels = chromosomes))

sample_chromosome_background <- function(chromosome, target_count) {
    chromosome_length <- seqlengths(genome_info)[chromosome]
    excluded <- ranges(consensus_regions[seqnames(consensus_regions) == chromosome])
    if (length(blacklist_regions) > 0L) {
        excluded <- c(
            excluded,
            ranges(blacklist_regions[seqnames(blacklist_regions) == chromosome])
        )
    }
    excluded <- reduce(excluded)
    allowed <- IRanges::setdiff(IRanges::IRanges(1L, chromosome_length), excluded)
    allowed <- allowed[width(allowed) >= window_bp]
    if (length(allowed) == 0L) {
        stop("No background sequence space remains on ", chromosome)
    }

    slot_counts <- floor(width(allowed) / window_bp)
    total_slots <- sum(slot_counts)
    if (total_slots < target_count) {
        stop(
            "Not enough non-overlapping background windows on ", chromosome,
            ": need ", target_count, ", found ", total_slots
        )
    }

    pool_size <- min(total_slots, max(target_count * 10L, target_count + 1000L))
    slot_ids <- sample.int(total_slots, size = pool_size, replace = FALSE)
    cumulative_slots <- cumsum(slot_counts)
    interval_index <- findInterval(slot_ids - 1L, cumulative_slots) + 1L
    previous_slots <- c(0, cumulative_slots)[interval_index]
    within_interval <- slot_ids - previous_slots
    candidate_starts <- start(allowed)[interval_index] +
        (within_interval - 1L) * window_bp
    candidate_ranges <- GRanges(
        seqnames = chromosome,
        ranges = IRanges::IRanges(start = candidate_starts, width = window_bp)
    )
    candidate_sequences <- scanFa(genome, param = candidate_ranges)
    candidate_n <- letterFrequency(
        candidate_sequences, letters = "N", as.prob = TRUE
    )[, 1]
    keep_candidates <- which(candidate_n <= max_fraction_n)
    candidate_ranges <- candidate_ranges[keep_candidates]
    candidate_sequences <- candidate_sequences[keep_candidates]

    foreground_for_chromosome <- foreground_sequences[
        foreground_chromosomes == chromosome
    ]
    foreground_gc <- rowSums(
        letterFrequency(
            foreground_for_chromosome, letters = c("G", "C"), as.prob = TRUE
        )
    )
    candidate_gc <- rowSums(
        letterFrequency(
            candidate_sequences, letters = c("G", "C"), as.prob = TRUE
        )
    )
    foreground_gc_bins <- findInterval(
        foreground_gc, gc_breaks, all.inside = TRUE
    )
    candidate_gc_bins <- findInterval(candidate_gc, gc_breaks, all.inside = TRUE)
    selected <- integer()
    for (gc_bin in sort(unique(foreground_gc_bins))) {
        needed <- sum(foreground_gc_bins == gc_bin) * background_multiplier
        available <- which(candidate_gc_bins == gc_bin)
        if (length(available) < needed) {
            stop(
                "Too few GC-matched background windows on ", chromosome,
                " in GC bin ", gc_bin, ": need ", needed, ", found ",
                length(available), ". Review the assembly or lower ",
                "background_multiplier."
            )
        }
        selected <- c(selected, available[seq_len(needed)])
    }
    if (length(selected) != target_count) {
        stop(
            "Internal background-count mismatch on ", chromosome, ": expected ",
            target_count, ", selected ", length(selected)
        )
    }
    list(
        ranges = candidate_ranges[selected],
        sequences = candidate_sequences[selected]
    )
}

background_parts <- lapply(seq_along(chromosomes), function(index) {
    target <- as.integer(foreground_counts[index]) * background_multiplier
    sample_chromosome_background(chromosomes[index], target)
})
background_ranges <- do.call(c, lapply(background_parts, `[[`, "ranges"))
background_sequences <- do.call(c, lapply(background_parts, `[[`, "sequences"))
if (length(background_sequences) == 0L) stop("No background sequences were sampled")
if (anyDuplicated(background_ranges)) stop("Duplicate background windows were sampled")

names(foreground_sequences) <- paste0("foreground_", seq_along(foreground_sequences))
names(background_sequences) <- paste0("background_", seq_along(background_sequences))
all_sequences <- c(foreground_sequences, background_sequences)
bins <- factor(
    c(
        rep("foreground", length(foreground_sequences)),
        rep("background", length(background_sequences))
    ),
    levels = c("foreground", "background")
)

pwms <- getMatrixSet(
    JASPAR2020,
    list(
        collection = "CORE",
        matrixtype = "PWM",
        tax_group = "vertebrates",
        all_versions = FALSE
    )
)
if (length(pwms) == 0L) stop("JASPAR returned no vertebrate CORE motifs")

enrichment <- calcBinnedMotifEnrR(
    seqs = all_sequences,
    bins = bins,
    pwmL = pwms,
    background = "otherBins",
    maxFracN = max_fraction_n,
    GCbreaks = gc_breaks,
    BPPARAM = SerialParam(RNGseed = seed),
    verbose = TRUE
)
foreground_column <- match("foreground", colnames(enrichment))
if (is.na(foreground_column)) stop("monaLisa did not return a foreground column")

motif_ids <- rownames(enrichment)
motif_names <- as.character(rowData(enrichment)$motif.name)
missing_names <- is.na(motif_names) | motif_names == ""
motif_names[missing_names] <- motif_ids[missing_names]

neg_log10_p <- assay(enrichment, "negLog10P")[, foreground_column]
raw_p <- 10^(-neg_log10_p)
adjusted_p <- p.adjust(raw_p, method = "BH")
neg_log10_adjusted <- -log10(adjusted_p)
log2_enrichment <- assay(enrichment, "log2enr")[, foreground_column]
results <- tibble(
    Factor = factor_label,
    Motif.ID = motif_ids,
    Motif.Name = motif_names,
    log2enr = log2_enrichment,
    negLog10P = neg_log10_p,
    negLog10Padj = neg_log10_adjusted,
    P.value = raw_p,
    P.adjust = adjusted_p,
    significant = !is.na(adjusted_p) & adjusted_p <= 0.05 & log2_enrichment > 0
) %>%
    arrange(desc(negLog10Padj), desc(log2enr), desc(negLog10P))

message(
    "Motif background: ", length(foreground_sequences), " foreground and ",
    length(background_sequences), " background windows; seed=", seed
)
write_tsv(results, file.path(outdir, "motif_enrichment.tsv"))

# Keep infinite scores in the output; they represent P values below floating-
# point precision. Cap only the plotted value so the PDF remains readable.
top_motifs <- results %>%
    filter(!is.na(negLog10Padj), !is.na(Motif.Name)) %>%
    slice_head(n = 20)
write_csv(top_motifs, file.path(outdir, "motif_summary.csv"))

finite_scores <- top_motifs$negLog10Padj[is.finite(top_motifs$negLog10Padj)]
plot_cap <- if (length(finite_scores) > 0L) max(10, max(finite_scores) + 1) else 10
plot_data <- top_motifs %>%
    mutate(plot_score = if_else(is.infinite(negLog10Padj), plot_cap, negLog10Padj))

pdf(file.path(outdir, "motif_summary.pdf"), width = 10, height = 8)
print(
    ggplot(plot_data, aes(x = reorder(Motif.Name, plot_score), y = plot_score)) +
        geom_col(fill = "steelblue") +
        coord_flip() +
        theme_minimal() +
        labs(
            title = paste("Top enriched motifs:", factor_label),
            subtitle = "monaLisa with JASPAR vertebrate CORE motifs",
            x = "Motif",
            y = "-log10(adjusted P-value)"
        )
)
dev.off()
