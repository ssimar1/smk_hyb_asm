#!/usr/bin/env bash
#
# build_rvdb.sh
#
# One-time setup: build the BLAST databases the tiered classification step
# needs. Run this ONCE per RVDB release, before running the pipeline. The
# output prefixes below must match the rvdb_nucl_db / rvdb_prot_db paths in
# config/config.yaml.
#
# Inputs you must have downloaded first:
#   - C-RVDBv32.0.fasta            (clustered nucleotide RVDB)     -> rvdb.dbi.udel.edu
#   - U-RVDBv32.0-prot_unique      (clustered protein RVDB-prot)   -> rvdb-prot.pasteur.fr
#   - RVDBv32_taxonomy.tab         (nucleotide taxonomy table)     -> rvdb.dbi.udel.edu
# The taxonomy table needs no build step; the pipeline reads it directly.
#
# NOTE: no -parse_seqids. RVDB headers contain multiple pipes, which breaks
# seqid parsing. The pipeline pulls accessions from the subject title (stitle)
# instead, so -parse_seqids is intentionally omitted.

set -euo pipefail

#conda init bash
#conda activate smk_hyb_asm

# --- Edit these to match your paths ---
RVDB_DIR="/data/projects/ssimar/db/RVDB"
NUCL_FASTA="${RVDB_DIR}/C-RVDBv32.0.fasta"
PROT_FASTA="${RVDB_DIR}/U-RVDBv32.0-prot_unique"
NUCL_OUT="${RVDB_DIR}/C-RVDBv32.0"      # must match config: rvdb_nucl_db
PROT_OUT="${RVDB_DIR}/RVDBv32-prot"     # must match config: rvdb_prot_db
# --------------------------------------

echo "[build_rvdb] Building nucleotide BLAST DB from ${NUCL_FASTA}"
makeblastdb \
    -in "${NUCL_FASTA}" \
    -dbtype nucl \
    -out "${NUCL_OUT}" \
    -title "C-RVDBv32.0"

echo "[build_rvdb] Building protein BLAST DB from ${PROT_FASTA}"
makeblastdb \
    -in "${PROT_FASTA}" \
    -dbtype prot \
    -out "${PROT_OUT}" \
    -title "RVDBv32-prot"

echo "[build_rvdb] Done. Index files written next to:"
echo "  ${NUCL_OUT}.*"
echo "  ${PROT_OUT}.*"
echo "[build_rvdb] Confirm these prefixes match rvdb_nucl_db / rvdb_prot_db in config/config.yaml"
