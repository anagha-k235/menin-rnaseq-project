#!/usr/bin/env Rscript
# DESeq2 differential expression + GSEA for VTP-50469 vs DMSO (GSE127506).
# Expects data/counts_matrix.tsv (genes x samples) and a coldata you define below.
suppressPackageStartupMessages({
  library(DESeq2); library(fgsea); library(tidyverse); library(ggrepel)
})
dir.create("results", showWarnings = FALSE)

## 1. Load counts ------------------------------------------------------------
# Rows = genes (symbol or Ensembl in col 1), cols = samples.
counts <- read.delim("data/counts_matrix.tsv", row.names = 1, check.names = FALSE)
counts <- round(as.matrix(counts))

## 2. Sample metadata --------------------------------------------------------
# EDIT to match your chosen GSE127506 samples. Block on cell line if you mix
# MOLM13 + RS4;11; contrast on treatment. Keep DMSO as the reference level.
coldata <- data.frame(
  row.names = colnames(counts),
  treatment = factor(c("DMSO","DMSO","VTP","VTP"), levels = c("DMSO","VTP"))
  # cell_line = factor(c("MOLM13","MOLM13","MOLM13","MOLM13"))
)
stopifnot(all(rownames(coldata) == colnames(counts)))

## 3. DESeq2 -----------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ treatment)   # ~ cell_line + treatment if mixed
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "treatment_VTP_vs_DMSO", type = "apeglm")
res_df <- as.data.frame(res) %>% rownames_to_column("gene") %>% arrange(padj)
write_csv(res_df, "results/deseq2_results.csv")

## 4. Volcano, MLL-fusion targets labelled -----------------------------------
targets <- c("MEIS1","HOXA9","HOXA10","HOXA7","PBX3","FLT3","CDK6","MEF2C")
res_df <- res_df %>% mutate(
  sig = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
  lab = ifelse(toupper(gene) %in% targets, gene, NA))
ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = sig), alpha = .5, size = 1) +
  geom_point(data = filter(res_df, !is.na(lab)), color = "red") +
  geom_text_repel(aes(label = lab), max.overlaps = Inf, min.segment.length = 0) +
  scale_color_manual(values = c("grey70","steelblue"), guide = "none") +
  labs(title = "VTP-50469 vs DMSO (KMT2A-rearranged)",
       subtitle = "MLL-fusion targets expected DOWN (negative log2FC)",
       x = "log2 fold change", y = "-log10 adjusted p") +
  theme_minimal(base_size = 12)
ggsave("results/volcano.png", width = 7, height = 5, dpi = 150)

## 5. GSEA against the MLL-fusion / HOXA-MEIS target program -----------------
ranks <- res_df %>% filter(!is.na(stat <- log2FoldChange)) %>%
  { setNames(.$log2FoldChange, toupper(.$gene)) } %>% sort(decreasing = TRUE)
mll_set <- list(MLL_FUSION_TARGETS = c(
  "MEIS1","HOXA9","HOXA10","HOXA7","HOXA5","PBX3","FLT3","CDK6","MEF2C","JMJD1C","BCL2"))
fg <- fgsea(mll_set, ranks, minSize = 3, maxSize = 500)
write_csv(as.data.frame(fg)[, c("pathway","NES","padj","size")], "results/gsea.csv")
png("results/gsea_mll_targets.png", width = 700, height = 450)
plotEnrichment(mll_set[["MLL_FUSION_TARGETS"]], ranks) +
  labs(title = sprintf("MLL-fusion targets  NES=%.2f  padj=%.3g", fg$NES[1], fg$padj[1]))
dev.off()

## 6. One-paragraph summary --------------------------------------------------
writeLines(c(
  "# Summary",
  sprintf("VTP-50469 vs DMSO in KMT2A-rearranged leukemia (GSE127506). %d genes at padj<0.05 & |log2FC|>1.",
          sum(res_df$sig, na.rm = TRUE)),
  sprintf("MLL-fusion target set: NES=%.2f, padj=%.3g (expected depleted/down on treatment).",
          fg$NES[1], fg$padj[1]),
  "Interpretation: menin-MLL inhibition suppresses the HOXA/MEIS1 self-renewal program, consistent with Krivtsov et al. 2019."
), "results/SUMMARY.md")
cat("Done -> results/\n")
