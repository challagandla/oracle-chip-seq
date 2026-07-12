# Shared figure style.
#
# Journal figures are read at column width (~89 mm) in print, so everything here
# is sized for that: 7 pt base type, hairline rules, no decorative fill. Colour is
# Okabe–Ito, which stays distinguishable under all common colour-vision
# deficiencies — a red/green "up vs down" volcano is unreadable to ~8% of male
# readers, and journals increasingly reject it.
#
# Peak mode has one fixed colour throughout every figure: narrow is blue, broad is
# orange, everywhere, so the reader learns the encoding once.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

OKABE_ITO <- c(
  black    = "#000000",
  orange   = "#E69F00",
  skyblue  = "#56B4E9",
  green    = "#009E73",
  yellow   = "#F0E442",
  blue     = "#0072B2",
  vermilion= "#D55E00",
  purple   = "#CC79A7"
)

# One encoding for peak mode, used by every panel.
MODE_COLS <- c(narrow = "#0072B2", broad = "#E69F00")

# Direction of change: blue = gained, vermilion = lost, grey = not significant.
DIR_COLS <- c(up = "#0072B2", down = "#D55E00", ns = "grey80")

STATUS_COLS <- c(PASS = "#009E73", WARN = "#E69F00", FAIL = "#D55E00")

theme_pub <- function(base_size = 7, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.border      = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"),
      axis.line         = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
      axis.text         = element_text(colour = "black", size = base_size),
      axis.title        = element_text(colour = "black", size = base_size + 1),
      plot.title        = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle     = element_text(size = base_size, colour = "grey30"),
      plot.caption      = element_text(size = base_size - 1, colour = "grey40", hjust = 0),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size + 1, face = "bold"),
      legend.key        = element_blank(),
      legend.key.size   = unit(3, "mm"),
      legend.title      = element_text(size = base_size),
      legend.text       = element_text(size = base_size),
      legend.background = element_blank()
    )
}

# mm -> inches, since journals specify figure widths in mm.
mm <- function(x) x / 25.4

save_pub <- function(plot, path, width_mm = 89, height_mm = 70) {
  ggsave(path, plot, width = mm(width_mm), height = mm(height_mm),
         units = "in", device = cairo_pdf, bg = "white")
  # A raster copy at 300 dpi for slides and for reviewers who cannot open vectors.
  ggsave(sub("\\.pdf$", ".png", path), plot,
         width = mm(width_mm), height = mm(height_mm), units = "in", dpi = 300, bg = "white")
  message("  wrote ", path)
}
