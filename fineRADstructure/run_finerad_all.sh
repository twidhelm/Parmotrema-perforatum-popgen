#!/usr/bin/env bash
set -euo pipefail

: "${PPG_DATA_DIR:?Set PPG_DATA_DIR to the PPG analysis directory}"
: "${FINERAD_TOOLS_DIR:?Set FINERAD_TOOLS_DIR to the directory containing finerad_input.py}"
: "${FINERAD_PROG_DIR:?Set FINERAD_PROG_DIR to the fineRADstructure installation directory}"

BASE_DIR="${PPG_DATA_DIR}"
FINE_DIR="${BASE_DIR}/FineRAD"
FINE_TOOLS="${FINERAD_TOOLS_DIR}"
FINE_PROG_DIR="${FINERAD_PROG_DIR}"

# Xtar-reduced datasets used for DAPC/PCA
DATASETS=("MIN4_Xtar_reduced" "MIN27_Xtar_reduced" "MIN53_Xtar_reduced")

# Minimum samples per locus for each dataset
declare -A MINSAMP
MINSAMP["MIN4_Xtar_reduced"]=4
MINSAMP["MIN27_Xtar_reduced"]=27
MINSAMP["MIN53_Xtar_reduced"]=53

cd "${FINE_DIR}"

for ASM in "${DATASETS[@]}"; do
    echo "==============================="
    echo "Processing dataset: ${ASM}"
    echo "==============================="

    OUTDIR="${BASE_DIR}/${ASM}_outfiles"
    ALLELES="${OUTDIR}/${ASM}.alleles"

    if [[ ! -f "${ALLELES}" ]]; then
        echo "ERROR: Alleles file not found: ${ALLELES}" >&2
        echo "Check the actual filename with:"
        echo "find ${BASE_DIR} -name '*${ASM}*.alleles'"
        continue
    fi

    MS=${MINSAMP[${ASM}]}
    PREFIX="${FINE_DIR}/${ASM}_min${MS}"

    echo "Using --minsample=${MS}"
    echo "Input alleles: ${ALLELES}"
    echo "Output prefix: ${PREFIX}"

    echo "[1/4] Running finerad_input.py ..."
    python "${FINE_TOOLS}/finerad_input.py" \
        --input "${ALLELES}" \
        --minsample "${MS}" \
        -o "${PREFIX}.finerad"

    echo "[2/4] Running RADpainter paint ..."
    "${FINE_PROG_DIR}/RADpainter" paint "${PREFIX}.finerad"

    CHUNKS_OUT="${PREFIX}_chunks.out"

    if [[ ! -f "${CHUNKS_OUT}" ]]; then
        echo "ERROR: Expected chunks file not found: ${CHUNKS_OUT}" >&2
        continue
    fi

    echo "[3/4] Running finestructure MCMC ..."
    "${FINE_PROG_DIR}/finestructure" \
        -x 100000 -y 100000 -z 1000 \
        "${CHUNKS_OUT}" \
        "${PREFIX}_chunks.mcmc.xml"

    echo "[4/4] Running finestructure tree building ..."
    "${FINE_PROG_DIR}/finestructure" \
        -m T -x 10000 \
        "${CHUNKS_OUT}" \
        "${PREFIX}_chunks.mcmc.xml" \
        "${PREFIX}_chunks.mcmcTree.xml"

    echo "Done with ${ASM}."
    echo
done

echo "All Xtar-reduced fineRADstructure datasets processed."