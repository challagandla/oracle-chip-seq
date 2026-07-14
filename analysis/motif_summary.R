args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript motif_summary.R homer_knownResults.txt outdir")

homer_file <- args[1]
outdir <- args[2]

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

library(readr)
library(dplyr)
library(ggplot2)

# HOMER column names themselves contain '#', so treating it as a comment marker
# truncates the header and silently corrupts the trailing columns.
motifs <- read_tsv(homer_file, col_names = TRUE, show_col_types = FALSE)
names(motifs) <- make.names(names(motifs), unique = TRUE)

if (!"Motif.Name" %in% colnames(motifs)) {
    stop("Expected a Motif Name column in HOMER results")
}
if (!"P.value" %in% colnames(motifs)) {
    stop("Expected a P-value column in HOMER results")
}

if ("Log.P.value" %in% colnames(motifs)) {
    motifs <- motifs %>%
        mutate(
            P.value = suppressWarnings(as.numeric(P.value)),
            Log.P.value = suppressWarnings(as.numeric(Log.P.value)),
            # HOMER reports the natural log; convert to the plotted -log10 scale.
            logP = -Log.P.value / log(10)
        )
} else {
    # Older HOMER tables may omit Log P-value. Keep underflowed values instead
    # of dropping the strongest motifs, while acknowledging double precision.
    motifs <- motifs %>%
        mutate(
            P.value = suppressWarnings(as.numeric(P.value)),
            logP = -log10(pmax(P.value, .Machine$double.xmin))
        )
}

top_motifs <- motifs %>%
    filter(!is.na(logP), is.finite(logP), logP >= 0) %>%
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
