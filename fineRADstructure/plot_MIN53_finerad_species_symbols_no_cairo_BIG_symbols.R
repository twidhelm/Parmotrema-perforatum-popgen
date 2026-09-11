#!/usr/bin/env Rscript

# MIN53-only fineRADstructure plot with PCA/DAPC-style species symbols
# -------------------------------------------------------------------
# This script is a stand-alone version of plot_finerad_Xtar_reduced.R.
# It plots only:
#   MIN53_Xtar_reduced_min53
#
# Instead of long sample names on the x/y axes, it replaces each label with
# a species symbol using the same species colors/shapes as the PCA/DAPC:
#   square = P. perforatum, P. subrigidum, P. preperforatum
#   circle = P. hypotropum, P. louisianae, P. hypoleucinum
#
# Input files expected in the working directory:
#   MIN53_Xtar_reduced_min53_chunks.out
#   MIN53_Xtar_reduced_min53_chunks.mcmcTree.xml
#   FinestructureLibrary.R

library(ape)
library(XML)

source("FinestructureLibrary.R")

prefix <- "MIN53_Xtar_reduced_min53"

# ------------------------------------------------------------
# Species parsing and PCA/DAPC-matched aesthetics
# ------------------------------------------------------------

get_species <- function(labels) {
  labs <- tolower(labels)

  species <- rep(NA_character_, length(labels))

  # Important: preperforatum must be matched before perforatum
  species[grepl("preperforatum", labs)] <- "P. preperforatum"
  species[grepl("hypoleucinum", labs)]  <- "P. hypoleucinum"
  species[grepl("hypotropum", labs)]    <- "P. hypotropum"
  species[grepl("subrigidum", labs)]    <- "P. subrigidum"
  species[grepl("louisianae", labs)]    <- "P. louisianae"

  # perforatum after preperforatum to avoid accidental matches
  species[grepl("perforatum", labs) & !grepl("preperforatum", labs)] <- "P. perforatum"

  species[is.na(species)] <- "Unknown"
  species
}

species_cols <- c(
  "P. perforatum"    = "#006400",  # dark green
  "P. hypotropum"    = "#90EE90",  # light green
  "P. subrigidum"    = "#8B0000",  # dark red
  "P. louisianae"    = "#F4A3A3",  # light red / salmon
  "P. preperforatum" = "#6A0DAD",  # dark purple
  "P. hypoleucinum"  = "#D8B7FF",  # light purple
  "Unknown"          = "black"
)

# Unicode symbols used as the axis labels.
# These are text labels, so they work with the existing plotFinestructure()
# label machinery without needing to rewrite FinestructureLibrary.R.
species_symbols <- c(
  "P. perforatum"    = "\u25A0", # square
  "P. hypotropum"    = "\u25CF", # circle
  "P. subrigidum"    = "\u25A0", # square
  "P. louisianae"    = "\u25CF", # circle
  "P. preperforatum" = "\u25A0", # square
  "P. hypoleucinum"  = "\u25CF", # circle
  "Unknown"          = "\u25CF"
)

species_pch <- c(
  "P. perforatum"    = 15,
  "P. hypotropum"    = 16,
  "P. subrigidum"    = 15,
  "P. louisianae"    = 16,
  "P. preperforatum" = 15,
  "P. hypoleucinum"  = 16,
  "Unknown"          = 16
)

species_order <- c(
  "P. perforatum",
  "P. hypotropum",
  "P. subrigidum",
  "P. louisianae",
  "P. preperforatum",
  "P. hypoleucinum"
)

get_symbol_labels <- function(labels) {
  sp <- get_species(labels)
  species_symbols[sp]
}

get_label_colors <- function(labels) {
  sp <- get_species(labels)
  species_cols[sp]
}

# ------------------------------------------------------------
# Outgroup detection
# ------------------------------------------------------------

is_outgroup <- function(labels) {
  labs <- tolower(labels)

  grepl("cristiferum", labs) |
    grepl("reticulatum", labs) |
    grepl("reticulated", labs) |
    grepl("outgroup", labs)
}

# ------------------------------------------------------------
# Helper for reading fineRADstructure chunks.out matrix
# ------------------------------------------------------------

is_numeric_vector <- function(x) {
  suppressWarnings(all(!is.na(as.numeric(x))))
}

read_chunks_matrix <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nchar(trimws(lines)) > 0]

  parts <- strsplit(trimws(lines), "\\s+")

  # Look for a square-ish matrix block:
  # header row with sample names, followed by n rows:
  # sample_name numeric numeric numeric ...
  for (i in seq_along(parts)) {
    header <- parts[[i]]

    # Try both possibilities:
    # 1. Header is exactly sample names
    # 2. Header has one leading non-sample field
    possible_headers <- list(header, header[-1])

    for (h in possible_headers) {
      n <- length(h)

      if (n < 5) next
      if ((i + n) > length(parts)) next

      candidate_rows <- parts[(i + 1):(i + n)]

      row_lengths_ok <- all(sapply(candidate_rows, length) >= n + 1)
      if (!row_lengths_ok) next

      numeric_ok <- all(sapply(candidate_rows, function(z) {
        vals <- z[2:(n + 1)]
        is_numeric_vector(vals)
      }))
      if (!numeric_ok) next

      row_names <- sapply(candidate_rows, function(z) z[1])

      mat <- do.call(rbind, lapply(candidate_rows, function(z) {
        as.numeric(z[2:(n + 1)])
      }))

      rownames(mat) <- row_names
      colnames(mat) <- h

      message("Read matrix from ", file, ": ", nrow(mat), " x ", ncol(mat))
      return(mat)
    }
  }

  stop("Could not identify a square coancestry matrix in: ", file)
}

# ------------------------------------------------------------
# Plot helper
# ------------------------------------------------------------

add_species_legend <- function() {
  present <- species_order[species_order %in% get_species(rownames(mat))]

  legend(
    "topleft",
    legend = present,
    col = species_cols[present],
    pch = species_pch[present],
    pt.cex = 2.8,
    cex = 1.5,
    bty = "n",
    title = "Species",
    xpd = NA
  )
}

plot_min53 <- function(pdf_file, png_file) {
  row_symbols <- get_symbol_labels(rownames(mat))
  col_symbols <- get_symbol_labels(colnames(mat))

  row_cols <- get_label_colors(rownames(mat))
  col_cols <- get_label_colors(colnames(mat))

  # In this matrix row and column order should match, so the colors should match.
  # plotFinestructure() accepts one text.col vector; use row colors.
  # If row/column order ever differs, reorder the matrix before plotting.
  if (!identical(rownames(mat), colnames(mat))) {
    warning("Row and column names are not identical. Reordering columns to match rows.")
    mat <<- mat[rownames(mat), rownames(mat), drop = FALSE]
    col_symbols <- get_symbol_labels(colnames(mat))
  }

  # PDF output. Use quartz on macOS when available so we do not depend on XQuartz/cairo.
  if (capabilities("aqua")) {
    grDevices::quartz(
      file = pdf_file,
      type = "pdf",
      width = 22,
      height = 22
    )
  } else {
    grDevices::pdf(
      file = pdf_file,
      width = 22,
      height = 22,
      onefile = FALSE
    )
  }

  plotFinestructure(
    tmpmat = mat,
    labelsx = col_symbols,
    labelsy = row_symbols,
    dend = dend,
    cex.axis = 2.2,
    xcrt = 90,
    ycrt = 0,
    labmargin = 12,
    labelsoff = c(1, 1),
    main = paste0(prefix, " - species symbols"),
    scalelabel = "Coancestry",
    text.col = row_cols
  )

  add_species_legend()
  dev.off()
  message("Wrote: ", pdf_file)

  # PNG output.
  # PNG output. Avoid type="cairo" because that requires XQuartz on some Macs.
  if (capabilities("aqua")) {
    png(
      filename = png_file,
      width = 6600,
      height = 6600,
      res = 300,
      type = "quartz"
    )
  } else {
    png(
      filename = png_file,
      width = 6600,
      height = 6600,
      res = 300
    )
  }

  plotFinestructure(
    tmpmat = mat,
    labelsx = col_symbols,
    labelsy = row_symbols,
    dend = dend,
    cex.axis = 2.2,
    xcrt = 90,
    ycrt = 0,
    labmargin = 12,
    labelsoff = c(1, 1),
    main = paste0(prefix, " - species symbols"),
    scalelabel = "Coancestry",
    text.col = row_cols
  )

  add_species_legend()
  dev.off()
  message("Wrote: ", png_file)
}

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

message("====================================")
message("Plotting ", prefix)
message("====================================")

chunks_file <- paste0(prefix, "_chunks.out")
tree_file   <- paste0(prefix, "_chunks.mcmcTree.xml")

if (!file.exists(chunks_file)) {
  stop("Missing chunks file: ", chunks_file)
}

if (!file.exists(tree_file)) {
  stop("Missing tree file: ", tree_file)
}

# Read coancestry matrix
mat <- read_chunks_matrix(chunks_file)

# Read fineSTRUCTURE tree
txml <- xmlTreeParse(tree_file)
tr <- extractTree(txml)

# Keep only samples shared between tree and matrix
shared <- intersect(tr$tip.label, rownames(mat))

if (length(shared) == 0) {
  stop("No shared sample labels between tree and matrix for ", prefix)
}

tr <- drop.tip(tr, setdiff(tr$tip.label, shared))
mat <- mat[tr$tip.label, tr$tip.label, drop = FALSE]

# Remove outgroups
outgroup_tips <- tr$tip.label[is_outgroup(tr$tip.label)]

if (length(outgroup_tips) > 0) {
  message("Dropping outgroups from ", prefix, ":")
  message(paste(outgroup_tips, collapse = ", "))

  tr <- drop.tip(tr, outgroup_tips)

  keep_samples <- tr$tip.label
  mat <- mat[keep_samples, keep_samples, drop = FALSE]
} else {
  message("No outgroups detected for ", prefix)
}

# Final check
if (any(is_outgroup(rownames(mat)))) {
  stop("Outgroups still present in matrix for ", prefix)
}

# Convert tree to dendrogram for plotting
dend <- myapetodend(tr)

# Output filenames
pdf_file <- paste0(prefix, "_fineRADstructure_heatmap_species_symbols_BIG.pdf")
png_file <- paste0(prefix, "_fineRADstructure_heatmap_species_symbols_BIG.png")

# RStudio sometimes leaves a tiny plotting device active; this prevents
# "figure margins too large" if a file device fails to open.
while (dev.cur() > 1) dev.off()

plot_min53(pdf_file, png_file)

message("Done.")
