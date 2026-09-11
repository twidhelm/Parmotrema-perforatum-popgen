## ================================================================
## Filter PPG Xtar VCF files for ingroup PCA and DAPC
##
## Workflow:
##   1. Begin with original ipyrad VCF files containing 106 samples.
##   2. Remove P. cristiferum and P. reticulatum.
##   3. Calculate MAF among the remaining 104 ingroup samples.
##   4. Retain variants with MAF >= 0.05.
##
## This script does not rebuild the ipyrad assemblies.
## ================================================================

suppressPackageStartupMessages({
  library(vcfR)
})

## ----------------------------------------------------------
## 1. Paths and settings
## ----------------------------------------------------------

base_dir <- Sys.getenv("PPG_DATA_DIR")

if (base_dir == "") {
  stop("Set PPG_DATA_DIR to the directory containing the PPG analysis files.")
}

out_dir <- file.path(
  base_dir,
  "VCF_Xtar_reduced_ingroup_MAF05"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

maf_filter <- 0.05

assemblies <- data.frame(
  asm = c(
    "MIN4_Xtar_reduced",
    "MIN27_Xtar_reduced",
    "MIN53_Xtar_reduced"
  ),
  vcf_path = c(
    file.path(
      base_dir,
      "MIN4_Xtar_reduced_outfiles",
      "MIN4_Xtar_reduced.vcf"
    ),
    file.path(
      base_dir,
      "MIN27_Xtar_reduced_outfiles",
      "MIN27_Xtar_reduced.vcf"
    ),
    file.path(
      base_dir,
      "MIN53_Xtar_reduced_outfiles",
      "MIN53_Xtar_reduced.vcf"
    )
  ),
  stringsAsFactors = FALSE
)

outgroup_samples <- c(
  "TW4688_P_cristiferum_Madagascar",
  "TW5198_P_reticulatum_USA_IL"
)

## ----------------------------------------------------------
## 2. Log file
## ----------------------------------------------------------

log_file <- file.path(
  out_dir,
  "filter_Xtar_ingroup_MAF05_log.txt"
)

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_msg <- function(...) {
  line <- paste(..., collapse = "")
  message(line)
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Starting ingroup MAF-filtering pipeline")
log_msg("MAF threshold: ", maf_filter)
log_msg(
  "Outgroups removed before MAF: ",
  paste(outgroup_samples, collapse = ", ")
)

## ----------------------------------------------------------
## 3. Filter each assembly
## ----------------------------------------------------------

for (i in seq_len(nrow(assemblies))) {

  asm <- assemblies$asm[i]
  vcf_path <- assemblies$vcf_path[i]

  log_msg("")
  log_msg("========================================")
  log_msg("Assembly: ", asm)
  log_msg("Input VCF: ", vcf_path)
  log_msg("========================================")

  if (!file.exists(vcf_path)) {
    stop("Input VCF not found: ", vcf_path)
  }

  output_prefix <- file.path(
    out_dir,
    paste0(asm, "_ingroup_maf05")
  )

  ## Construct one --remove-indv argument for each outgroup.
  remove_args <- paste(
    "--remove-indv",
    shQuote(outgroup_samples),
    collapse = " "
  )

  cmd <- paste(
    "vcftools",
    "--vcf", shQuote(vcf_path),
    remove_args,
    "--maf", maf_filter,
    "--recode",
    "--out", shQuote(output_prefix)
  )

  log_msg("Running:")
  log_msg(cmd)

  status <- system(cmd)

  if (status != 0) {
    stop(
      "VCFtools failed for ",
      asm,
      " with exit status ",
      status
    )
  }

  output_vcf <- paste0(output_prefix, ".recode.vcf")

  if (!file.exists(output_vcf)) {
    stop("Expected output VCF was not created: ", output_vcf)
  }

  ## --------------------------------------------------------
  ## Validate the new VCF
  ## --------------------------------------------------------

  log_msg("Validating filtered VCF...")

  vcf <- read.vcfR(output_vcf, verbose = FALSE)

  ## First gt column is FORMAT; remaining columns are samples.
  sample_names <- colnames(vcf@gt)[-1]

  n_samples <- length(sample_names)
  n_records <- nrow(vcf@fix)

  remaining_outgroups <- intersect(
    sample_names,
    outgroup_samples
  )

  if (length(remaining_outgroups) > 0) {
    stop(
      "Outgroups remain in filtered VCF: ",
      paste(remaining_outgroups, collapse = ", ")
    )
  }

  if (n_samples != 104) {
    stop(
      "Expected 104 ingroup samples, but found ",
      n_samples,
      " in ",
      asm
    )
  }

  log_msg("Output VCF: ", output_vcf)
  log_msg("Samples retained: ", n_samples)
  log_msg("Variant records retained: ", n_records)
  log_msg("Outgroups present: no")

  rm(vcf)
  gc()
}

log_msg("")
log_msg("Finished all three assemblies")