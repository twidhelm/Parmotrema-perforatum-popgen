## PPG MIN85 site80 PCA: remove individuals with >=20% missing genotypes.
## Independent sensitivity check; preserves all previous analyses.
suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(ggplot2)
})

base_dir <- "/home/FM/twidhelm/Parmotrema_RADseq/PPG"
vcf_path <- file.path(base_dir, "VCF_Xtar_reduced_ingroup_MAF05",
                      "MIN85_Xtar_reduced_ingroup_maf05_site80.recode.vcf")
miss_path <- file.path(base_dir, "VCF_Xtar_reduced_ingroup_MAF05",
                       "MIN85_Xtar_reduced_ingroup_maf05_site80_final_indv.imiss")
meta_path <- file.path(base_dir, "PPG_metadata.csv")
out_dir <- file.path(base_dir, "PCA_MIN85_site80_lowmissing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c(vcf_path, miss_path, meta_path)) {
  if (!file.exists(f)) stop("Missing input file: ", f)
}

miss <- read.delim(miss_path, check.names = FALSE)
if (anyDuplicated(miss$INDV)) stop("Duplicate individuals in missingness file")
keep_ids <- miss$INDV[!is.na(miss$F_MISS) & miss$F_MISS < 0.20]

vcf <- read.vcfR(vcf_path, verbose = FALSE)
vcf_ids <- colnames(vcf@gt)[-1L]
if (!setequal(miss$INDV, vcf_ids)) {
  stop("Missingness file and VCF contain different sample names")
}
keep_ids <- vcf_ids[vcf_ids %in% keep_ids]
if (length(keep_ids) < 10L) stop("Too few retained individuals")

meta <- read.csv(meta_path, stringsAsFactors = FALSE)
meta$sample <- trimws(as.character(meta$sample))
if (anyDuplicated(meta$sample)) stop("Duplicate samples in metadata")
meta_sub <- meta[match(keep_ids, meta$sample), , drop = FALSE]
if (anyNA(meta_sub$sample)) stop("Missing metadata for retained samples")
meta_sub$species <- tolower(trimws(as.character(meta_sub$species)))
meta_sub$species <- sub("^p[._]", "", meta_sub$species)
meta_sub$group <- paste0("P.", meta_sub$species)

levels_sp <- c("P.perforatum", "P.hypotropum", "P.subrigidum",
               "P.louisianae", "P.preperforatum", "P.hypoleucinum")
if (!all(meta_sub$group %in% levels_sp)) stop("Unexpected species labels")
meta_sub$group <- factor(meta_sub$group, levels = levels_sp)

## Report missingness of actual VCF genotypes within the reduced sample set.
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
gt <- gt[, keep_ids, drop = FALSE]
missing_matrix <- is.na(gt) | grepl("\\.", gt)
site_miss <- rowMeans(missing_matrix)
cat("Individuals retained:", length(keep_ids), "of", length(vcf_ids), "\n")
cat("Individuals excluded:", length(vcf_ids) - length(keep_ids), "\n")
cat("Mean site missingness after removal:", round(mean(site_miss) * 100, 2), "%\n")
cat("Maximum site missingness after removal:", round(max(site_miss) * 100, 2), "%\n")
cat("Sites with >=20% missingness after removal:", sum(site_miss >= 0.20), "\n")

## Rebuild genotype matrix and refit PCA, rather than merely hiding points.
gl <- vcfR2genlight(vcf)
ploidy(gl) <- rep(1L, nInd(gl))
indNames(gl) <- trimws(indNames(gl))
gl <- gl[match(keep_ids, indNames(gl))]
if (!identical(indNames(gl), keep_ids)) stop("Genlight sample-order mismatch")
X <- tab(gl, NA.method = "mean", freq = FALSE)
valid <- apply(X, 2L, function(z) all(is.finite(z)) && stats::var(z) > 0)
X <- X[, valid, drop = FALSE]
if (ncol(X) < 4L) stop("Fewer than four informative SNP columns")
fit <- prcomp(X, center = TRUE, scale. = FALSE)
variance <- 100 * fit$sdev^2 / sum(fit$sdev^2)

scores <- as.data.frame(fit$x[, 1:4, drop = FALSE])
scores$sample <- keep_ids
scores$species <- meta_sub$species
scores$group <- meta_sub$group
scores$F_MISS <- miss$F_MISS[match(keep_ids, miss$INDV)]
write.table(scores, file.path(out_dir, "MIN85_lowmissing_PCA_scores.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(PC = seq_along(variance), percent_variance = variance),
            file.path(out_dir, "MIN85_lowmissing_PCA_variance.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

colors <- c("P.perforatum" = "#006400", "P.hypotropum" = "#90EE90",
            "P.subrigidum" = "#8B0000", "P.louisianae" = "#F08080",
            "P.preperforatum" = "#4B0082", "P.hypoleucinum" = "#D8BFD8")
shapes <- c("P.perforatum" = 15, "P.hypotropum" = 16,
            "P.subrigidum" = 15, "P.louisianae" = 16,
            "P.preperforatum" = 15, "P.hypoleucinum" = 16)
for (pair in list(c(1, 2), c(1, 3))) {
  x <- paste0("PC", pair[1]); y <- paste0("PC", pair[2])
  p <- ggplot(scores, aes(x = .data[[x]], y = .data[[y]],
                          color = group, shape = group)) +
    geom_point(size = 3) +
    scale_color_manual(values = colors, drop = FALSE) +
    scale_shape_manual(values = shapes, drop = FALSE) +
    theme_minimal() +
    labs(x = sprintf("%s (%.1f%%)", x, variance[pair[1]]),
         y = sprintf("%s (%.1f%%)", y, variance[pair[2]]),
         title = "MIN85 PCA: individuals with <20% missing genotypes",
         color = "Species", shape = "Species") +
    theme(plot.title = element_text(hjust = 0.5))
  tag <- paste0(x, "_", y)
  ggsave(file.path(out_dir, paste0("MIN85_lowmissing_", tag, ".pdf")),
         p, width = 7, height = 5)
  ggsave(file.path(out_dir, paste0("MIN85_lowmissing_", tag, ".png")),
         p, width = 7, height = 5, dpi = 300)
}
saveRDS(list(pca = fit, scores = scores, percent_var = variance,
             retained_samples = keep_ids, retained_SNP_columns = colnames(X)),
        file.path(out_dir, "MIN85_lowmissing_PCA_objects.rds"))
cat("Cleaned SNP columns:", ncol(X), "\n")
cat(sprintf("PC1 = %.2f%%, PC2 = %.2f%%, PC3 = %.2f%%\n",
            variance[1], variance[2], variance[3]))
cat("Species counts:\n")
print(table(scores$group))
cat("Results:", out_dir, "\n")
