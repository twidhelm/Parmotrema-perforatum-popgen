# Original Phoebe scripts for MIN85 sensitivity analyses

This folder contains **exact UTF-8 source files** from the archive supplied
by the corresponding author on 2026-11-23, named
`MIN85_GitHub_sources.tar.gz`. They preserve original server paths,
comments, analysis settings and file names. They are retained unchanged
for provenance; scripts in `PCA/`, `DAPC/` and `filtering/` are portable
adaptations intended for reuse on other systems.

The blob SHAs below were validated byte-for-byte against Git blob SHA-1
hashes calculated from the supplied archive.

| File | Git blob SHA-1 |
| --- | --- |
| `filter_Xtar_ingroup_MAF05.R` | `2b05f92c4832f1a24207f70fb96752d0bb5575e6` |
| `PPG_PCA_MIN85_site80.R` | `d98c6d28fcfcb6c7e311566f433f549220e152c5` |
| `PPG_DAPC_MIN85_site80.R` | `35a1d2ed92d6a54c34ba11a469e37c8c86264264` |
| `PPG_PCA_MIN85_site80_lowmissing.R` | `f492ef4136bea0ff68e41913ea538da11912f60e` |
| `PPG_PCA_MIN85_site80_lowmissing_strict.R` | `d1faf5d9378129508ee8abc34a3e1851114d1763` |

The ipyrad parameter file is at
`../../ipyrad/params-MIN85_Xtar_reduced.txt`, exact Git blob
SHA-1 `014a299daa172fd54d499e887a85d03a4c626c51`.

**Important:** The original `filter_Xtar_ingroup_MAF05.R` in this
archive predates MIN85 and lists only MIN4/MIN27/MIN53. The MIN85 filtering
steps were executed with VCFtools and are documented in the
`filtering/filter_MIN85_ingroup_site80.sh` reproduction script.

The primary MIN85 PCA/DAPC used 104 ingroup individuals and 833
cleaned SNP columns. The final stricter robustness PCA retained 85
individuals and 827 cleaned SNP columns after reapplying SNP-site
missingness. Neither the original large input VCFs nor the RDS
analysis objects are included in this small source-code repository.
No new genomic analyses were executed when these source files were archived.
