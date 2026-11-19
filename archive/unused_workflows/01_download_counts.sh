#!/usr/bin/env bash
# Fast path: pull the processed gene-count matrix for GSE127506 from GEO.
# GEO stores processed files under the series' supplementary archive.
set -euo pipefail
mkdir -p data
cd data

GSE=GSE127506
# GEO supplementary tarball for the series (contains processed count tables).
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/${GSE}/suppl/${GSE}_RAW.tar"

echo "Downloading ${GSE} supplementary archive..."
curl -fL -o "${GSE}_RAW.tar" "${URL}"
tar -xvf "${GSE}_RAW.tar"
gunzip -f ./*.gz 2>/dev/null || true

echo
echo "Done. Inspect the extracted per-sample count files in data/,"
echo "then build a genes x samples matrix (scripts/02 expects data/counts_matrix.tsv)."
echo "Also download the series matrix for the sample metadata:"
echo "  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/${GSE}/matrix/${GSE}_series_matrix.txt.gz"


