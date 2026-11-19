#!/usr/bin/env bash
# Full path: download a subset of raw runs for GSE127506 from SRA.
# Fill RUNS with the SRR accessions for the MOLM-13 (or RS4;11) VTP-50469 vs DMSO
# samples you want — read them off the GSE127506 "SRA Run Selector".
set -euo pipefail
mkdir -p data/fastq
cd data/fastq

# EDIT: replace with the SRR accessions from GSE127506 Run Selector
# (pick one cell line, both conditions, with replicates — e.g. 3x DMSO + 3x VTP-50469).
RUNS=(
  # SRRXXXXXXX   # MOLM13 DMSO rep1
  # SRRXXXXXXX   # MOLM13 DMSO rep2
  # SRRXXXXXXX   # MOLM13 VTP-50469 rep1
  # SRRXXXXXXX   # MOLM13 VTP-50469 rep2
)

if [ ${#RUNS[@]} -eq 0 ]; then
  echo "No SRR accessions set. Open the GSE127506 SRA Run Selector, pick your"
  echo "contrast samples, and paste their SRR IDs into the RUNS array above."
  exit 1
fi

for srr in "${RUNS[@]}"; do
  echo ">>> $srr"
  prefetch "$srr"
  fasterq-dump "$srr" --split-files --threads 4 -O .
  gzip -f "${srr}"*.fastq
done
echo "Done. FASTQs in data/fastq/ — now run: snakemake --cores 4"
