#!/usr/bin/env Rscript
# Publication figure panels.
#
# Every panel that can be split by peak mode is, because the narrow/broad
# distinction is the claim this pipeline makes and a figure is the honest place to
# show whether it holds. Figure 2 in particular is the direct test: if the registry
# is right, narrow targets produce sub-kb peaks and broad targets produce
# multi-kb domains, and the width distributions must separate.

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(patchwork)
  library(ggrepel)
  library(yaml)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--results", type = "character", default = "results"),
  make_option("--registry", type = "character", default = "config/mark_registry.yaml"),
  make_option("--samples", type = "character", default = "samples.tsv"),
  make_option("--config", type = "character", default = "config.yaml"),
  make_option("--outdir", type = "character", default = "results/figures")
)))

source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "theme_pub.R"))

R <- opt$results
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

registry <- yaml.load_file(opt$registry)
cfg <- yaml.load_file(opt$config)
samples <- read_tsv(opt$samples, show_col_types = FALSE)

mode_of <- function(target) {
  e <- registry$marks[[target]]
  if (is.null(e$peak_mode)) registry$defaults$peak_mode else e$peak_mode
}
targets <- sort(unique(samples$target[samples$assay == "chip"]))
modes <- setNames(vapply(targets, mode_of, character(1)), targets)
# Narrow marks first so every legend and facet strip orders consistently.
target_levels <- names(sort(modes))

message("targets: ", paste(sprintf("%s (%s)", names(modes), modes), collapse = ", "))

theme_set(theme_pub())

# =====================================================================
# Figure 1 — QC. Enrichment judged against each mark's own threshold.
# =====================================================================

gate <- read_tsv(file.path(R, "qc", "qc_gate.tsv"), show_col_types = FALSE)

frip_thresh <- tibble(
  target = targets,
  min_frip = vapply(targets, function(t) {
    e <- registry$marks[[t]]
    v <- if (is.null(e$qc$min_frip)) registry$defaults$qc$min_frip else e$qc$min_frip
    as.numeric(v)
  }, numeric(1)),
  peak_mode = modes[targets]
)

chip <- gate %>%
  filter(assay == "chip") %>%
  mutate(target = factor(target, levels = target_levels))

p1a <- ggplot(chip, aes(target, FRiP, colour = peak_mode)) +
  geom_hline(data = frip_thresh %>% mutate(target = factor(target, levels = target_levels)),
             aes(yintercept = min_frip), linetype = "22", linewidth = 0.3, colour = "grey40") +
  geom_point(aes(shape = condition), size = 1.6, position = position_jitter(width = 0.12, height = 0)) +
  scale_colour_manual(values = MODE_COLS, name = "peak mode") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL, y = "FRiP",
       title = "a",
       subtitle = "Dashed line: threshold for that mark") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p1b <- ggplot(gate, aes(target, usable_reads / 1e6, fill = assay)) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.3, width = 0.6) +
  geom_hline(yintercept = cfg$qc$min_usable_reads_narrow / 1e6,
             linetype = "22", linewidth = 0.3, colour = "grey40") +
  scale_fill_manual(values = c(chip = OKABE_ITO[["skyblue"]], input = "grey75")) +
  labs(x = NULL, y = "usable reads (M)", title = "b",
       subtitle = "Dashed: ENCODE narrow-mark floor (20M)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p1c <- ggplot(gate, aes(NRF, PBC1, colour = status)) +
  geom_vline(xintercept = cfg$qc$min_nrf, linetype = "22", linewidth = 0.3, colour = "grey60") +
  geom_hline(yintercept = cfg$qc$min_pbc1, linetype = "22", linewidth = 0.3, colour = "grey60") +
  geom_point(size = 1.4) +
  scale_colour_manual(values = STATUS_COLS, name = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "NRF", y = "PBC1", title = "c", subtitle = "Library complexity")

p1d <- ggplot(gate, aes(fragment_length, reorder(sample, fragment_length))) +
  geom_col(aes(fill = assay), width = 0.7) +
  scale_fill_manual(values = c(chip = OKABE_ITO[["skyblue"]], input = "grey75"), guide = "none") +
  labs(x = "fragment length (bp)", y = NULL, title = "d") +
  theme(axis.text.y = element_text(size = 5))

fig1 <- (p1a | p1b) / (p1c | p1d)
save_pub(fig1, file.path(opt$outdir, "fig1_qc.pdf"), width_mm = 183, height_mm = 130)

# =====================================================================
# Figure 2 — the peak-mode rule, tested.
#
# If the registry classification is correct, peak widths must separate by mode.
# This is the panel that would expose a mis-classified mark: an H3K27ac called
# with --broad would show a multi-kb width distribution and land with H3K27me3.
# =====================================================================

read_peaks <- function(path, target) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  df <- read_tsv(path, col_names = FALSE, show_col_types = FALSE,
                 col_types = cols(.default = col_character()))
  tibble(target = target, width = as.numeric(df$X3) - as.numeric(df$X2))
}

peak_files <- samples %>%
  filter(assay == "chip") %>%
  mutate(ext = ifelse(modes[target] == "broad", "broadPeak", "narrowPeak"),
         path = file.path(R, "peaks", paste0(sample_id, "_peaks.", ext)))

widths <- bind_rows(lapply(seq_len(nrow(peak_files)), function(i) {
  read_peaks(peak_files$path[i], peak_files$target[i])
}))

widths <- widths %>%
  mutate(peak_mode = modes[target],
         target = factor(target, levels = target_levels))

p2a <- ggplot(widths, aes(width, target, fill = peak_mode)) +
  geom_violin(scale = "width", linewidth = 0.2, colour = "white", alpha = 0.9) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25, fill = "white") +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale()),
                breaks = c(100, 1e3, 1e4, 1e5, 1e6)) +
  scale_fill_manual(values = MODE_COLS, name = "peak mode") +
  annotation_logticks(sides = "b", linewidth = 0.2,
                      short = unit(0.5, "mm"), mid = unit(0.8, "mm"), long = unit(1.2, "mm")) +
  labs(x = "peak width (bp)", y = NULL, title = "a",
       subtitle = "Called peak widths separate by registry class")

counts <- widths %>% count(target, peak_mode)
p2b <- ggplot(counts, aes(n, target, fill = peak_mode)) +
  geom_col(width = 0.7, position = position_dodge2(preserve = "single")) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_fill_manual(values = MODE_COLS, guide = "none") +
  labs(x = "peaks per library", y = NULL, title = "b")

med <- widths %>%
  group_by(target, peak_mode) %>%
  summarise(median_width = median(width), .groups = "drop")
p2c <- ggplot(med, aes(median_width, reorder(target, median_width), fill = peak_mode)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = comma(round(median_width))), hjust = -0.15, size = 2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25)),
                     labels = label_number(scale_cut = cut_short_scale())) +
  scale_fill_manual(values = MODE_COLS, guide = "none") +
  labs(x = "median peak width (bp)", y = NULL, title = "c")

fig2 <- p2a / (p2b | p2c) + plot_layout(heights = c(1.4, 1))
save_pub(fig2, file.path(opt$outdir, "fig2_peak_landscape.pdf"), width_mm = 140, height_mm = 130)

# =====================================================================
# Figure 3 — sample relationships on log2(ChIP/Input) signal.
# =====================================================================

cor_path <- file.path(R, "qc", "correlation", "spearman_matrix.tsv")
if (file.exists(cor_path)) {
  cm <- read.delim(cor_path, row.names = 1, check.names = FALSE)
  colnames(cm) <- rownames(cm)
  long <- as.data.frame(as.table(as.matrix(cm)))
  names(long) <- c("a", "b", "rho")
  ord <- hclust(as.dist(1 - as.matrix(cm)))$order
  lv <- rownames(cm)[ord]
  long$a <- factor(long$a, levels = lv)
  long$b <- factor(long$b, levels = lv)

  p3a <- ggplot(long, aes(a, b, fill = rho)) +
    geom_tile() +
    scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2",
                         midpoint = 0, limits = c(-1, 1), name = expression(rho)) +
    coord_equal() +
    labs(x = NULL, y = NULL, title = "a",
         subtitle = "Spearman, log2(ChIP/Input), 10 kb bins") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
          axis.text.y = element_text(size = 5),
          panel.grid.major.y = element_blank())

  # PCA on the same binned signal matrix.
  sig <- read.delim(file.path(R, "qc", "correlation", "signal_matrix.tsv"),
                    check.names = FALSE, comment.char = "")
  mat <- as.matrix(sig[, -(1:3)])
  mat <- mat[complete.cases(mat), , drop = FALSE]
  mat <- mat[apply(mat, 1, var) > 0, , drop = FALSE]
  colnames(mat) <- gsub("^'|'$", "", colnames(mat))
  pca <- prcomp(t(mat), scale. = TRUE)
  ve <- round(100 * summary(pca)$importance[2, 1:2])
  pd <- as.data.frame(pca$x[, 1:2])
  pd$sample <- rownames(pd)
  pd <- pd %>%
    left_join(samples %>% select(sample_id, target, condition), by = c("sample" = "sample_id")) %>%
    mutate(peak_mode = modes[target])

  p3b <- ggplot(pd, aes(PC1, PC2, colour = target, shape = condition)) +
    geom_point(size = 2) +
    geom_text_repel(aes(label = target), size = 2, show.legend = FALSE,
                    max.overlaps = 20, segment.size = 0.2) +
    scale_colour_manual(values = unname(OKABE_ITO[c(6, 2, 4, 8, 3, 7)])[seq_along(targets)],
                        breaks = targets) +
    labs(x = sprintf("PC1 (%d%%)", ve[1]), y = sprintf("PC2 (%d%%)", ve[2]),
         title = "b", subtitle = "Samples separate by mark, then by condition")

  fig3 <- p3a | p3b
  save_pub(fig3, file.path(opt$outdir, "fig3_sample_relationships.pdf"),
           width_mm = 183, height_mm = 85)
}

# =====================================================================
# Figure 4 — differential binding: one volcano per target.
# =====================================================================

volc <- list()
stats <- list()
for (t in target_levels) {
  f <- file.path(R, "differential", t, "results.tsv")
  if (!file.exists(f)) next
  d <- read_tsv(f, show_col_types = FALSE)
  d$target <- t
  d$peak_mode <- modes[t]
  volc[[t]] <- d
  stats[[t]] <- tibble(
    target = t, peak_mode = modes[t],
    tested = nrow(d),
    up = sum(d$direction == "up", na.rm = TRUE),
    down = sum(d$direction == "down", na.rm = TRUE)
  )
}

if (length(volc)) {
  V <- bind_rows(volc) %>%
    filter(!is.na(padj)) %>%
    mutate(target = factor(target, levels = target_levels),
           negLog = -log10(pmax(padj, 1e-300)))

  p4 <- ggplot(V, aes(log2FoldChange, negLog, colour = direction)) +
    geom_point(size = 0.35, alpha = 0.6) +
    geom_vline(xintercept = c(-cfg$differential$min_lfc, cfg$differential$min_lfc),
               linetype = "22", linewidth = 0.25, colour = "grey50") +
    geom_hline(yintercept = -log10(cfg$differential$fdr),
               linetype = "22", linewidth = 0.25, colour = "grey50") +
    scale_colour_manual(values = DIR_COLS, name = NULL,
                        breaks = c("up", "down"),
                        labels = c(sprintf("up in %s", cfg$contrast$treatment),
                                   sprintf("up in %s", cfg$contrast$reference))) +
    facet_wrap(~ target, nrow = 1, scales = "free_y") +
    guides(colour = guide_legend(override.aes = list(size = 1.8))) +
    labs(x = sprintf("log2 fold change (%s / %s)",
                     cfg$contrast$treatment, cfg$contrast$reference),
         y = expression(-log[10]~italic(P)[adj]),
         title = "a", subtitle = "Differential binding per mark")

  S <- bind_rows(stats) %>%
    pivot_longer(c(up, down), names_to = "direction", values_to = "n") %>%
    mutate(target = factor(target, levels = target_levels),
           n_signed = ifelse(direction == "down", -n, n))

  p5 <- ggplot(S, aes(target, n_signed, fill = direction)) +
    geom_col(width = 0.65) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_text(aes(label = comma(n), vjust = ifelse(direction == "up", -0.4, 1.2)), size = 2) +
    scale_fill_manual(values = DIR_COLS, guide = "none") +
    scale_y_continuous(labels = function(x) comma(abs(x))) +
    labs(x = NULL, y = "differential peaks", title = "b") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  fig4 <- p4 / p5 + plot_layout(heights = c(1.5, 1))
  save_pub(fig4, file.path(opt$outdir, "fig4_differential.pdf"),
           width_mm = 183, height_mm = 120)

  write_tsv(bind_rows(stats), file.path(opt$outdir, "differential_summary.tsv"))
}

# =====================================================================
# Figure 5 — genomic distribution and GO enrichment.
# =====================================================================

anno <- list()
for (t in target_levels) {
  f <- file.path(R, "annotation", t, "peak_annotation.tsv")
  if (!file.exists(f)) next
  d <- read_tsv(f, show_col_types = FALSE)
  if (!"annotation" %in% names(d)) next
  anno[[t]] <- tibble(target = t, peak_mode = modes[t],
                      feature = sub(" \\(.*", "", d$annotation))
}

if (length(anno)) {
  A <- bind_rows(anno) %>%
    mutate(feature = recode(feature,
                            "5' UTR" = "5'UTR", "3' UTR" = "3'UTR",
                            "Distal Intergenic" = "Distal intergenic",
                            "Downstream" = "Downstream")) %>%
    count(target, peak_mode, feature) %>%
    group_by(target) %>% mutate(frac = n / sum(n)) %>% ungroup() %>%
    mutate(target = factor(target, levels = target_levels))

  p6 <- ggplot(A, aes(frac, target, fill = feature)) +
    geom_col(width = 0.7) +
    scale_x_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
    scale_fill_manual(values = unname(OKABE_ITO[c(6, 3, 2, 4, 8, 7, 5, 1)]), name = NULL) +
    labs(x = "peaks", y = NULL, title = "a",
         subtitle = "Genomic distribution of consensus peaks")

  go <- list()
  for (t in target_levels) {
    f <- file.path(R, "annotation", t, "go_enrichment.tsv")
    if (!file.exists(f) || file.info(f)$size < 10) next
    d <- try(read_tsv(f, show_col_types = FALSE), silent = TRUE)
    if (inherits(d, "try-error") || nrow(d) == 0 || !"Description" %in% names(d)) next
    d$target <- t
    go[[t]] <- d
  }

  if (length(go)) {
    G <- bind_rows(go) %>%
      group_by(target, direction) %>%
      slice_min(p.adjust, n = 6, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(label = paste0(Description),
             label = ifelse(nchar(label) > 45, paste0(substr(label, 1, 42), "..."), label))

    p7 <- ggplot(G, aes(-log10(p.adjust), reorder(label, -log10(p.adjust)),
                        colour = direction, size = Count)) +
      geom_point() +
      scale_colour_manual(values = DIR_COLS, name = NULL) +
      scale_size_continuous(range = c(1, 3.5), name = "genes") +
      facet_wrap(~ target, scales = "free_y", ncol = 2) +
      labs(x = expression(-log[10]~italic(q)), y = NULL, title = "b",
           subtitle = "GO biological process, genes near differential peaks") +
      theme(axis.text.y = element_text(size = 5.5))

    fig5 <- p6 / p7 + plot_layout(heights = c(1, 2.2))
    save_pub(fig5, file.path(opt$outdir, "fig5_annotation_go.pdf"),
             width_mm = 183, height_mm = 175)
  } else {
    save_pub(p6, file.path(opt$outdir, "fig5_annotation_go.pdf"),
             width_mm = 140, height_mm = 60)
  }
}

# =====================================================================
# Figure 6 — motif enrichment. Only exists for targets the registry allows.
# =====================================================================

mot <- list()
for (t in target_levels) {
  for (dir in c("up", "down")) {
    f <- file.path(R, "motifs", t, dir, "motif_enrichment.tsv")
    if (!file.exists(f)) next
    d <- try(read_tsv(f, show_col_types = FALSE), silent = TRUE)
    if (inherits(d, "try-error") || nrow(d) == 0) next
    d$target <- t
    d$direction <- dir
    # Rank by adjusted significance, but only among motifs actually enriched
    # (log2enr > 0). A depleted motif is a different statement and does not belong
    # on the same axis.
    mot[[paste(t, dir)]] <- d %>%
      filter(is.finite(negLog10Padj), log2enr > 0) %>%
      slice_max(negLog10Padj, n = 8, with_ties = FALSE)
  }
}

if (length(mot)) {
  M <- bind_rows(mot) %>%
    filter(negLog10Padj > 0) %>%
    mutate(facet = paste0(target, " — ", direction)) %>%
    arrange(facet, negLog10Padj)

  if (nrow(M)) {
    # facet_wrap(scales = "free") shares one y factor across panels, so a plain
    # reorder() would order every panel by the global ranking. Build a unique key
    # per row and relabel it with the motif name.
    M <- M %>% mutate(key = factor(seq_len(n()), levels = seq_len(n()), labels = name))

    p8 <- ggplot(M, aes(negLog10Padj, key, fill = direction)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = -log10(0.05), linetype = "22",
                 linewidth = 0.25, colour = "grey40") +
      scale_fill_manual(values = DIR_COLS, guide = "none") +
      facet_wrap(~ facet, scales = "free", ncol = 2) +
      labs(x = expression(-log[10]~italic(q)), y = NULL,
           title = "Motif enrichment in differential peaks (JASPAR / monaLisa)",
           subtitle = "Background: unchanged peaks of the same mark. Dashed: q = 0.05") +
      theme(axis.text.y = element_text(size = 5))
    save_pub(p8, file.path(opt$outdir, "fig6_motifs.pdf"),
             width_mm = 183, height_mm = 30 + 22 * ceiling(length(unique(M$facet)) / 2))
  }
}

message("figures written to ", opt$outdir)
