#!/usr/bin/env bash
set -euo pipefail

# Portable reproduction of the RAxML commands recorded in the
# RAxML_info files from the analyses reported in the manuscript.

RAXML_BIN="${RAXML_BIN:-raxmlHPC-PTHREADS-AVX}"
THREADS="${THREADS:-20}"

PARSIMONY_SEED=194955
BOOTSTRAP_SEED=12345
BOOTSTRAPS=100
OUTGROUP="TW4688_P_cristiferum_Madagascar"

DATASETS=(
  "MIN4_Xtar_reduced"
  "MIN27_Xtar_reduced"
  "MIN53_Xtar_reduced"
)

for ASM in "${DATASETS[@]}"; do
    ALIGNMENT="${ASM}.phy"

    if [[ ! -f "${ALIGNMENT}" ]]; then
        echo "ERROR: alignment not found: ${ALIGNMENT}" >&2
        exit 1
    fi

    echo "========================================"
    echo "Running RAxML: ${ASM}"
    echo "========================================"

    "${RAXML_BIN}" \
        -n "${ASM}" \
        -s "${ALIGNMENT}" \
        -m GTRGAMMA \
        -f a \
        -p "${PARSIMONY_SEED}" \
        -x "${BOOTSTRAP_SEED}" \
        -# "${BOOTSTRAPS}" \
        -T "${THREADS}" \
        -o "${OUTGROUP}"
done
