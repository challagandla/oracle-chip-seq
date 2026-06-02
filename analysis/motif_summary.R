args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript motif_summary.R homer_knownResults.txt outdir")

homer_file <- args[1]
outdir <- args[2]

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

library(readr)
library(dplyr)
library(ggplot2)

motifs <- read_tsv(homer_file, comment = "#", col_names = TRUE)
if (!"Motif.Name" %in% colnames(motifs)) {
    stop("Expected a Motif.Name column in HOMER results")
}

top_motifs <- motifs %>%
    mutate(logP = -log10(P.value)) %>%
    arrange(desc(logP)) %>%
    slice_head(n = 20)

write_csv(top_motifs, file.path(outdir, "motif_summary.csv"))

pdf(file.path(outdir, "motif_summary.pdf"), width = 10, height = 8)
ggplot(top_motifs, aes(x = reorder(Motif.Name, logP), y = logP)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top HOMER Motifs", x = "Motif", y = "-log10(P-value)")
dev.off()
