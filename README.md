# Population genomics of the *Parmotrema perforatum* group

Analysis scripts, parameter files, metadata, and provenance records associated with the manuscript:

**Population genomics reveals contrasting genomic patterns among lichen species pairs in the *Parmotrema perforatum* group (Parmeliaceae, Ascomycota)**

This repository contains code and configuration files used for RADseq assembly, phylogenetic inference, principal component analysis (PCA), discriminant analysis of principal components (DAPC), and fineRADstructure analyses.

Raw sequencing reads and large phylogenetic alignment files are archived separately in public repositories and are not duplicated here.

## Repository structure

```text
.
├── ipyrad/
│   ├── params-MIN4_Xtar_reduced.txt
│   ├── params-MIN27_Xtar_reduced.txt
│   └── params-MIN53_Xtar_reduced.txt
├── filtering/
│   └── filter_Xtar_ingroup_MAF05.R
├── PCA/
│   └── PPG_PCA_Xtar_reduced_ingroup_MAF05_all_loci.R
├── DAPC/
│   └── PPG_DAPC_Xtar_reduced_ingroup_MAF05.R
├── fineRADstructure/
│   ├── run_finerad_all.sh
│   └── plot_MIN53_finerad_species_symbols_no_cairo_BIG_symbols.R
├── RAxML/
│   ├── run_raxml.sh
│   ├── RAxML_info.MIN4_Xtar_reduced
│   ├── RAxML_info.MIN27_Xtar_reduced
│   └── RAxML_info.MIN53_Xtar_reduced
└── metadata/
    └── PPG_metadata.csv
```

## RADseq assembly

Assemblies were generated with **ipyrad v0.9.102** using the *Xanthoparmelia taractica* reference genome (NCBI accession **GCA_033085395.1**).

Three final datasets were analyzed:

- **MIN4**: loci retained when present in at least 4 samples
- **MIN27**: loci retained when present in at least 27 samples
- **MIN53**: loci retained when present in at least 53 samples

The original ipyrad parameter files are provided in `ipyrad/`. These retain the original local paths used during analysis for provenance. Users reproducing the analyses will need to substitute paths appropriate to their own system.

## SNP filtering

For PCA and DAPC, the two outgroup samples

- `TW4688_P_cristiferum_Madagascar`
- `TW5198_P_reticulatum_USA_IL`

were removed before applying a minor allele frequency threshold of 0.05.

The final filtered VCFs contained 104 ingroup samples:

| Dataset | VCF variant records after MAF filtering |
|---|---:|
| MIN4 | 116,258 |
| MIN27 | 99,532 |
| MIN53 | 71,775 |

The filtering workflow is implemented in `filtering/filter_Xtar_ingroup_MAF05.R`.

The analysis scripts use the environment variable `PPG_DATA_DIR` to identify the directory containing the analysis inputs. For example:

```bash
export PPG_DATA_DIR=/path/to/PPG
Rscript filtering/filter_Xtar_ingroup_MAF05.R
```

## Principal component analysis

PCA was performed on the filtered ingroup datasets using all usable loci after genotype conversion and removal of non-finite or invariant loci.

| Dataset | Individuals | Cleaned loci | PC1 | PC2 |
|---|---:|---:|---:|---:|
| MIN4 | 104 | 113,349 | 21.74% | 12.19% |
| MIN27 | 104 | 96,990 | 22.75% | 12.63% |
| MIN53 | 104 | 69,954 | 23.41% | 12.84% |

The PCA workflow is implemented in `PCA/PPG_PCA_Xtar_reduced_ingroup_MAF05_all_loci.R`.

## DAPC

DAPC was performed using current species determinations as predefined groups.

For each dataset:

- up to 10,000 loci were retained using random seed `999`;
- PC retention was evaluated with `adegenet::xvalDapc()`;
- candidate PC values ranged from 5 to 50 in increments of 5;
- 100 cross-validation replicates were performed;
- five discriminant axes were retained.

| Dataset | Retained PCs |
|---|---:|
| MIN4 | 10 |
| MIN27 | 15 |
| MIN53 | 10 |

The workflow is implemented in `DAPC/PPG_DAPC_Xtar_reduced_ingroup_MAF05.R`.

## fineRADstructure

Coancestry analyses used RADpainter and fineRADstructure following Malinsky et al. (2018).

The installed analysis software reported:

- **RADpainter v0.3.3 r111**
- fineRADstructure source/build package version **0.3.1**

The final datasets used minimum-sample thresholds corresponding to their ipyrad assemblies:

- MIN4: `--minsample=4`
- MIN27: `--minsample=27`
- MIN53: `--minsample=53`

The analysis workflow used:

```text
RADpainter paint
finestructure -x 100000 -y 100000 -z 1000
finestructure -m T -x 10000
```

The analysis script is `fineRADstructure/run_finerad_all.sh`.

The manuscript fineRADstructure heatmap was based on the MIN53 analysis and plotted with `fineRADstructure/plot_MIN53_finerad_species_symbols_no_cairo_BIG_symbols.R`.

The shell script expects:

```bash
export PPG_DATA_DIR=/path/to/PPG
export FINERAD_TOOLS_DIR=/path/to/fineRADstructure-tools
export FINERAD_PROG_DIR=/path/to/fineRADstructure
```

## Maximum-likelihood phylogenetic analyses

Maximum-likelihood analyses were conducted with **RAxML v8.2.12** under the GTRGAMMA model.

Each analysis used rapid bootstrapping followed by a thorough maximum-likelihood search (`-f a`) with 100 bootstrap replicates.

Recorded settings were:

```text
Parsimony seed:       194955
Rapid-bootstrap seed: 12345
Bootstrap replicates: 100
Threads:              20
Rooting taxon:        TW4688_P_cristiferum_Madagascar
```

RAxML reported:

| Dataset | Distinct alignment patterns |
|---|---:|
| MIN4 | 633,013 |
| MIN27 | 558,601 |
| MIN53 | 414,393 |

The exact program-generated records from the completed analyses are retained as:

```text
RAxML/RAxML_info.MIN4_Xtar_reduced
RAxML/RAxML_info.MIN27_Xtar_reduced
RAxML/RAxML_info.MIN53_Xtar_reduced
```

These files are preserved unchanged and therefore contain the original server paths used during analysis.

`RAxML/run_raxml.sh` provides a portable reproduction of the commands recorded in those files.

## Metadata

`metadata/PPG_metadata.csv` contains the sample-to-species assignments and associated geographic, chemotype, reproductive-mode, species-pair, and outgroup metadata used by the PCA and DAPC scripts.

## Data availability

Raw Illumina reads are available in the NCBI Sequence Read Archive (SRA) under BioProject accession **PRJNA1527490** and SRA Study accession **SRP735742**. Individual BioSample and SRA Run accessions are provided in Supplementary Table S2.

Phylogenetic alignment and tree files are being deposited in TreeBASE.

Repository accession numbers and permanent links will be added after deposition is complete.

## Citation

Please cite the associated manuscript when using these data or scripts.

## Reference

Malinsky M, Trucchi E, Lawson DJ, Falush D. 2018. RADpainter and fineRADstructure: population inference from RADseq data. *Molecular Biology and Evolution* 35:1284–1290. https://doi.org/10.1093/molbev/msy023
