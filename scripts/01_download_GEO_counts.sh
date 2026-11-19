
#!/usr/bin/env bash
set -euo pipefail

GSE="GSE127506"

mkdir -p data/raw data/metadata

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/${GSE}"

echo "Downloading processed MOLM13/RS411 counts..."

curl -fL --retry 3 \
  -o "data/raw/GSE127506_MOLM13_RS411_EPZ5676_VTP50469_d2_d7_180927_htseq_counts_rmdup.txt.gz" \
  "${BASE}/suppl/GSE127506_MOLM13_RS411_EPZ5676_VTP50469_d2_d7_180927_htseq_counts_rmdup.txt.gz"

echo "Downloading GEO sample metadata..."

curl -fL --retry 3 \
  -o "data/metadata/GSE127506_series_matrix.txt.gz" \
  "${BASE}/matrix/GSE127506_series_matrix.txt.gz"

echo "Validating compressed files..."

gzip -t data/raw/*.gz
gzip -t data/metadata/*.gz

echo
echo "Download complete!"
echo "Processed counts: data/raw/"
echo "Sample metadata: data/metadata/"

