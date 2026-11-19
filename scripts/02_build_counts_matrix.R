# ============================================================
# 02_build_counts_matrix.R
# GSE127506: MOLM13, Day 2
# Comparison: VTP50469 vs DMSO
# ============================================================

# Locate project root
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) == 0) {
  stop("Please run this script using Rscript.")
}

script_path <- normalizePath(
  sub("^--file=", "", file_arg[1])
)

project_root <- dirname(dirname(script_path))
setwd(project_root)

# Input file
input_file <- file.path(
  "data", "raw",
  "GSE127506_MOLM13_RS411_EPZ5676_VTP50469_d2_d7_180927_htseq_counts_rmdup.txt.gz"
)

if (!file.exists(input_file)) {
  stop("Count file not found: ", input_file)
}

# Read original count matrix
message("Reading GEO count data...")

all_counts <- read.delim(
  gzfile(input_file),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Select six MOLM13 Day 2 samples
sample_names <- grep(
  "MOLM13_RNA_(DMSO|VTP50469).*_d2_",
  colnames(all_counts),
  value = TRUE
)

if (length(sample_names) != 6) {
  stop("Expected 6 samples, found ", length(sample_names))
}

# Extract gene symbols and counts
gene_ids <- all_counts[[1]]

count_matrix <- as.matrix(all_counts[, sample_names])
storage.mode(count_matrix) <- "integer"

rownames(count_matrix) <- make.unique(gene_ids)

# Create sample metadata
samples <- data.frame(
  treatment = ifelse(
    grepl("DMSO", sample_names),
    "DMSO",
    "VTP50469"
  ),
  row.names = sample_names
)

samples$treatment <- factor(
  samples$treatment,
  levels = c("DMSO", "VTP50469")
)

# Check sample counts
stopifnot(all(colnames(count_matrix) == rownames(samples)))
stopifnot(sum(samples$treatment == "DMSO") == 3)
stopifnot(sum(samples$treatment == "VTP50469") == 3)

# Create output directories
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)

# Save count matrix
write.table(
  count_matrix,
  "data/processed/GSE127506_MOLM13_Day2_counts.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Save metadata
write.table(
  samples,
  "data/metadata/GSE127506_MOLM13_Day2_samples.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

message("Count matrix dimensions: ", paste(dim(count_matrix), collapse = " x "))
print(table(samples$treatment))
message("Data preparation complete.")
