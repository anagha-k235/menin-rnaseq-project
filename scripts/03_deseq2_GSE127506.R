dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

packages1 <- c("ggplot2", "ggrepel")

sapply(packages1, requireNamespace, quietly = TRUE)

# Check required packages
packages <- c("DESeq2", "apeglm", "fgsea", "msigdbr")

missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  stop(
    "Please install these packages before running the script: ",
    paste(missing, collapse = ", ")
  )
}


# Load processed count matrix and sample metadata
count_matrix <- read.delim(
  "data/processed/GSE127506_MOLM13_Day2_counts.tsv",
  row.names = 1,
  check.names = FALSE
)

samples <- read.delim(
  "data/metadata/GSE127506_MOLM13_Day2_samples.tsv",
  row.names = 1,
  check.names = FALSE
)

# Ensure count data are integers
count_matrix <- as.matrix(count_matrix)
storage.mode(count_matrix) <- "integer"

# Ensure sample metadata treatment is a factor
samples$treatment <- factor(
  samples$treatment,
  levels = c("DMSO", "VTP50469")
)

# Verify that count columns match metadata rows
stopifnot(all(colnames(count_matrix) == rownames(samples)))

cat("Count matrix dimensions:", dim(count_matrix), "\n")
print(table(samples$treatment))


library(DESeq2)
library(msigdbr)
library(fgsea)
library(apeglm)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = samples,
  design = ~ treatment
)

# Remove genes with very low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run differential expression analysis
dds <- DESeq(dds)

# View results: VTP50469 compared with DMSO
results_table <- results(
  dds,
  contrast = c("treatment", "VTP50469", "DMSO")
)

# Sort by adjusted p-value
results_table <- results_table[
  order(results_table$padj, na.last = TRUE),
]

head(as.data.frame(results_table))

# Save the complete differential expression results
write.csv(
  as.data.frame(results_table),
  "results/GSE127506_MOLM13_Day2_DESeq2.csv",
  row.names = TRUE
)

summary(results_table)

sum(results_table$padj < 0.05, na.rm = TRUE)

# Significant genes at adjusted p-value < 0.05
sig_genes <- as.data.frame(results_table)
sig_genes <- sig_genes[
  !is.na(sig_genes$padj) & sig_genes$padj < 0.05,
]

cat("Significant upregulated:", 
    sum(sig_genes$log2FoldChange > 0), "\n")

cat("Significant downregulated:", 
    sum(sig_genes$log2FoldChange < 0), "\n")

# Show the 15 most significant genes
sig_genes <- sig_genes[order(sig_genes$padj), ]

head(
  sig_genes[, c("log2FoldChange", "pvalue", "padj")],
  15
)


# Variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plot
p <- plotPCA(
  vsd,
  intgroup = "treatment"
)

p <- p +
  ggplot2::ggtitle("MOLM13 Day 2: VTP50469 vs DMSO") +
  ggplot2::theme_classic()

print(p)

# Save PCA plot
ggplot2::ggsave(
  "results/GSE127506_MOLM13_Day2_PCA.png",
  plot = p,
  width = 7,
  height = 5,
  dpi = 300
)


## Further Investigation

library(ggrepel)

pca_data <- plotPCA(
  vsd,
  intgroup = "treatment",
  returnData = TRUE
)

pca_labeled <- ggplot(
  pca_data,
  aes(PC1, PC2, color = group, label = name)
) +
  geom_point(size = 4) +
  geom_text_repel(size = 3) +
  ggtitle("MOLM13 Day 2: VTP50469 vs DMSO") +
  xlab("PC1: 80% variance") +
  ylab("PC2: 15% variance") +
  theme_classic()

print(pca_labeled)

sample_cor <- cor(
  assay(vsd),
  method = "pearson"
)

round(sample_cor, 3)

# CHeck replicate

pca_data[, c("name", "group", "PC1", "PC2")]
sample_cor <- cor(
  assay(vsd),
  method = "pearson"
)

round(sample_cor, 3)

boxplot(
  log2(counts(dds, normalized = TRUE) + 1),
  las = 2,
  main = "Normalized counts by sample",
  ylab = "log2 normalized counts"
)

round(cor(assay(vsd), method = "pearson"), 3)

# Extract DESeq2 results
res_df <- as.data.frame(results_table)
res_df$gene <- rownames(res_df)

# Exclude genes without a test statistic
res_df <- res_df[!is.na(res_df$stat), ]

# Exclude duplicated original gene symbols
# (these were assigned unique row names during count preparation)
res_df <- res_df[!duplicated(res_df$gene), ]

# Rank genes by DESeq2 Wald statistic
res_df <- as.data.frame(results_table)
res_df$feature_id <- rownames(res_df)

# Keep genes with valid DESeq2 statistics
res_df <- res_df[!is.na(res_df$stat), ]

# Remove suffixes added by make.unique()
res_df$gene <- sub("\\.[0-9]+$", "", res_df$feature_id)

# Retain one row per gene symbol for GSEA
res_df <- res_df[!duplicated(res_df$gene), ]

# Create ranked gene list
ranked_genes <- res_df$stat
names(ranked_genes) <- res_df$gene

ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Verify
length(ranked_genes)
sum(duplicated(names(ranked_genes)))


# Retrieve human Hallmark gene sets
hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

# Convert to pathways format for fgsea
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

# Run GSEA
fgsea_results <- fgsea(
  pathways = pathways,
  stats = ranked_genes,
  minSize = 15,
  maxSize = 500,
  nperm = 10000
)

# Sort by adjusted p-value
fgsea_results <- fgsea_results[
  order(fgsea_results$padj),
]

# Display significant pathways
fgsea_significant <- fgsea_results[
  !is.na(fgsea_results$padj) &
    fgsea_results$padj < 0.05,
]

head(
  as.data.frame(fgsea_significant)[
    , c("pathway", "NES", "pval", "padj", "size")
  ],
  20
)



# Write

fgsea_export <- as.data.frame(fgsea_results)

# Convert the leadingEdge list into a semicolon-separated string
fgsea_export$leadingEdge <- vapply(
  fgsea_export$leadingEdge,
  function(x) paste(x, collapse = ";"),
  character(1)
)


write.csv(
  fgsea_export,
  "results/GSE127506_MOLM13_Day2_Hallmark_GSEA.csv",
  row.names = FALSE
)


# Inspect Leading Genes

key_pathways <- c(
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_HYPOXIA"
)

fgsea_results[
  pathway %in% key_pathways,
  .(pathway, NES, padj, leadingEdge)
]

leading_edge_list <- setNames(
  fgsea_results$leadingEdge[
    match(key_pathways, fgsea_results$pathway)
  ],
  key_pathways
)

# Genes shared by at least two pathways
all_leading_genes <- unlist(leading_edge_list)

gene_frequency <- sort(
  table(all_leading_genes),
  decreasing = TRUE
)

gene_frequency[gene_frequency >= 2]


shared_genes <- names(gene_frequency[gene_frequency >= 2])

res_df <- as.data.frame(results_table)
res_df$gene <- sub("\\.[0-9]+$", "", rownames(res_df))

res_df[
  res_df$gene %in% shared_genes,
  c("gene", "log2FoldChange", "pvalue", "padj")
]


# Volcano Plot

volcano_data <- as.data.frame(results_table)

volcano_data$gene <- rownames(volcano_data)

volcano_data$significance <- "Not significant"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange > 0
] <- "Upregulated"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange < 0
] <- "Downregulated"

ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significance
  )
) +
  geom_point(alpha = 0.7, size = 1.5, na.rm = TRUE) +
  scale_color_manual(
    values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not significant" = "grey"
    )
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential Gene Expression: VTP50469 vs DMSO",
    subtitle = "MOLM13 cells, Day 2",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Expression"
  ) +
  theme_classic()

ggsave(
  "results/GSE127506_MOLM13_Day2_volcano.png",
  width = 8,
  height = 6,
  dpi = 300
)

table(volcano_data$significance)

gene_summary <- as.data.frame(
  table(volcano_data$significance)
)

colnames(gene_summary) <- c("Category", "Number_of_genes")

write.csv(
  gene_summary,
  "results/GSE127506_volcano_gene_summary.csv",
  row.names = FALSE
)


# Heatmap of Top 20 DEG

library(pheatmap)

top_genes <- volcano_data[
  !is.na(volcano_data$padj),
]

top_genes <- top_genes[
  order(top_genes$padj),
]

top20 <- head(top_genes$gene, 20)

top20

normalized_counts <- counts(dds, normalized = TRUE)

heatmap_counts <- normalized_counts[
  top20,
]

heatmap_log <- log2(heatmap_counts + 1)

heatmap_scaled <- t(scale(t(heatmap_log)))

annotation_col <- data.frame(
  Treatment = colData(dds)$treatment
)

rownames(annotation_col) <- colnames(heatmap_scaled)

pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  main = "Top 20 Differentially Expressed Genes"
)

pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  main = "Top 20 Differentially Expressed Genes",
  filename = "results/GSE127506_MOLM13_Day2_heatmap.png",
  width = 9,
  height = 8
)


#GSEA Enrichment Plot

library(fgsea)

plotEnrichment(
  pathways[["HALLMARK_MYC_TARGETS_V2"]],
  ranked_genes
) +
  labs(title = "MYC Targets V2")

ggsave(
  "results/MYC_TARGETS_V2_GSEA.png",
  width = 7,
  height = 5,
  dpi = 300
)

plotEnrichment(
  pathways[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
  ranked_genes
) +
  labs(title = "TNFA Signaling via NFKB")


ggsave(
  "results/TNFA_NFKB_GSEA.png",
  width = 7,
  height = 5,
  dpi = 300
)

plotEnrichment(
  pathways[["HALLMARK_MYC_TARGETS_V2"]],
  ranked_genes
) +
  labs(title = "MYC Targets V2: NES = -2.23")
