# ============================================================
# 04_timecourse_and_signatures.R
# GSE127506 extension: MOLM13 VTP-50469 vs DMSO, Day 2 and Day 7
#   (1) Per-timepoint DESeq2 + a formal time x treatment interaction test
#   (2) Preranked GSEA (fgseaMultilevel) at each timepoint against:
#         - MSigDB Hallmark
#         - a curated menin-KMT2A / HOXA-MEIS target panel
#         - a curated myeloid-differentiation panel
# Shows the on-target program suppressed early AND the differentiation
# program emerging by Day 7. Uses the Day-7 columns already present in
# the raw matrix downloaded by scripts/01 (no extra download).
# Run from the project root, after scripts/01_download_GEO_counts.sh.
# ============================================================

# ---- Resolve project root (works under Rscript) ----
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 0) stop("Please run this script using Rscript.")
project_root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
setwd(project_root)

# ---- Package check ----
pkgs <- c("DESeq2", "fgsea", "msigdbr", "ggplot2")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install these packages first: ", paste(missing, collapse = ", "))
}
suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- Load the raw GEO count matrix ----
raw_file <- file.path(
  "data", "raw",
  "GSE127506_MOLM13_RS411_EPZ5676_VTP50469_d2_d7_180927_htseq_counts_rmdup.txt.gz"
)
if (!file.exists(raw_file)) {
  stop("Raw count file not found: ", raw_file,
       "\nRun scripts/01_download_GEO_counts.sh first.")
}

message("Reading raw GEO count matrix...")
all_counts <- read.delim(
  gzfile(raw_file), header = TRUE,
  check.names = FALSE, stringsAsFactors = FALSE
)
gene_ids <- all_counts[[1]]

# ---- Select MOLM13 DMSO / VTP-50469 at Day 2 and Day 7 ----
# Parse the design from the column names rather than hard-coding them.
sample_names <- grep(
  "MOLM13_RNA_(DMSO|VTP50469).*_d(2|7)_",
  colnames(all_counts), value = TRUE
)
if (length(sample_names) == 0) {
  stop("No MOLM13 Day2/Day7 DMSO/VTP50469 columns found; check the raw header.")
}

treatment <- ifelse(grepl("DMSO", sample_names), "DMSO", "VTP50469")
time      <- ifelse(grepl("_d2_", sample_names), "d2", "d7")

coldata <- data.frame(
  treatment = factor(treatment, levels = c("DMSO", "VTP50469")),
  time      = factor(time,      levels = c("d2", "d7")),
  row.names = sample_names
)

# Require both treatments present at both timepoints
tab <- table(coldata$time, coldata$treatment)
print(tab)
stopifnot(all(tab >= 2))

count_matrix <- as.matrix(all_counts[, sample_names, drop = FALSE])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- make.unique(gene_ids)

message("Selected ", length(sample_names), " samples across Day 2 / Day 7.")

# ---- Helper: strip make.unique() suffix for gene-level GSEA ----
clean_symbol <- function(x) sub("\\.[0-9]+$", "", x)

# ---- Helper: per-timepoint DESeq2 (VTP50469 vs DMSO) ----
run_timepoint <- function(tp) {
  keep_cols <- rownames(coldata)[coldata$time == tp]
  cm  <- count_matrix[, keep_cols, drop = FALSE]
  cd  <- droplevels(coldata[keep_cols, , drop = FALSE])

  dds <- DESeqDataSetFromMatrix(cm, cd, design = ~ treatment)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("treatment", "VTP50469", "DMSO"))
  res <- as.data.frame(res)
  res$feature_id <- rownames(res)
  res$gene <- clean_symbol(res$feature_id)
  res$timepoint <- tp
  res
}

res_d2 <- run_timepoint("d2")
res_d7 <- run_timepoint("d7")

n_sig <- function(r) sum(r$padj < 0.05, na.rm = TRUE)
message("Significant genes (padj<0.05) — Day 2: ", n_sig(res_d2),
        " | Day 7: ", n_sig(res_d7))

# ---- Formal time x treatment interaction test (LRT) ----
# Genes whose VTP-50469 response differs between Day 2 and Day 7.
dds_full <- DESeqDataSetFromMatrix(
  count_matrix, coldata, design = ~ time + treatment + time:treatment
)
dds_full <- dds_full[rowSums(counts(dds_full)) >= 10, ]
dds_full <- DESeq(dds_full, test = "LRT", reduced = ~ time + treatment)
res_int <- as.data.frame(results(dds_full))
res_int$gene <- clean_symbol(rownames(res_int))
res_int <- res_int[order(res_int$padj), ]
write.csv(res_int, "results/MOLM13_timecourse_interaction_LRT.csv", row.names = TRUE)
message("Interaction (time:treatment) genes at padj<0.05: ",
        sum(res_int$padj < 0.05, na.rm = TRUE))

# ---- Merged per-timepoint log2FC table ----
m <- merge(
  res_d2[, c("gene", "log2FoldChange", "padj")],
  res_d7[, c("gene", "log2FoldChange", "padj")],
  by = "gene", suffixes = c("_d2", "_d7")
)
write.csv(m, "results/MOLM13_timecourse_DESeq2_d2_vs_d7.csv", row.names = FALSE)

# ---- Curated signatures (canonical biology; hand-curated panels) ----
mll_targets <- c(
  "HOXA9", "HOXA10", "HOXA7", "HOXA5", "MEIS1", "PBX3",
  "MEF2C", "MEF2D", "FLT3", "CDK6", "JMJD1C", "MYB", "SIX1", "SIX4"
)
differentiation <- c(
  "ITGAM", "ITGAX", "CD14", "CSF1R", "CEBPA", "CEBPE", "CEBPD",
  "LYZ", "MPO", "ELANE", "FCGR1A", "CD86", "S100A8", "S100A9"
)
custom_sets <- list(
  MLL_KMT2A_TARGETS      = mll_targets,
  MYELOID_DIFFERENTIATION = differentiation
)

# ---- Hallmark sets ----
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(hallmark$gene_symbol, hallmark$gs_name)

# ---- Preranked GSEA (fgseaMultilevel) at a timepoint ----
ranked_stats <- function(res) {
  r <- res[!is.na(res$stat), ]
  r <- r[!duplicated(r$gene), ]
  v <- r$stat
  names(v) <- r$gene
  sort(v, decreasing = TRUE)
}

run_gsea <- function(res, pathways, min_size) {
  stats <- ranked_stats(res)
  set.seed(1)
  fg <- fgseaMultilevel(pathways = pathways, stats = stats,
                        minSize = min_size, maxSize = 500)
  fg$timepoint <- res$timepoint[1]
  fg[order(fg$padj), ]
}

# Hallmark at each timepoint (parallels script 03's Day-2 run)
hm_d2 <- run_gsea(res_d2, hallmark_sets, min_size = 15)
hm_d7 <- run_gsea(res_d7, hallmark_sets, min_size = 15)
flatten_le <- function(df) {
  df$leadingEdge <- vapply(df$leadingEdge, paste, character(1), collapse = ";")
  df
}
write.csv(flatten_le(as.data.frame(hm_d2)),
          "results/MOLM13_Day2_Hallmark_GSEA_multilevel.csv", row.names = FALSE)
write.csv(flatten_le(as.data.frame(hm_d7)),
          "results/MOLM13_Day7_Hallmark_GSEA_multilevel.csv", row.names = FALSE)

# Curated signatures at each timepoint (small sets -> lower minSize)
cs_d2 <- run_gsea(res_d2, custom_sets, min_size = 5)
cs_d7 <- run_gsea(res_d7, custom_sets, min_size = 5)
custom_gsea <- rbind(
  as.data.frame(cs_d2)[, c("pathway", "NES", "pval", "padj", "size", "timepoint")],
  as.data.frame(cs_d7)[, c("pathway", "NES", "pval", "padj", "size", "timepoint")]
)
write.csv(custom_gsea, "results/MOLM13_custom_signature_GSEA.csv", row.names = FALSE)

# ============================================================
# Figures
# ============================================================

# ---- (A) Target-panel log2FC, Day 2 vs Day 7 ----
panel_df <- m[m$gene %in% mll_targets, ]
panel_long <- rbind(
  data.frame(gene = panel_df$gene, log2FC = panel_df$log2FoldChange_d2, timepoint = "Day 2"),
  data.frame(gene = panel_df$gene, log2FC = panel_df$log2FoldChange_d7, timepoint = "Day 7")
)
panel_long$gene <- factor(panel_long$gene,
  levels = panel_df$gene[order(panel_df$log2FoldChange_d2)])

p_panel <- ggplot(panel_long, aes(x = gene, y = log2FC, fill = timepoint)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("Day 2" = "#4C78A8", "Day 7" = "#B03A2E")) +
  labs(title = "Menin-KMT2A target program: VTP-50469 vs DMSO",
       subtitle = "Day 2 vs Day 7, MOLM13",
       x = NULL, y = "log2 fold change", fill = NULL) +
  theme_classic()
ggsave("figures/MOLM13_target_panel_timecourse.png", p_panel,
       width = 7, height = 5, dpi = 300)

# ---- (B) Differentiation-panel log2FC, Day 2 vs Day 7 ----
diff_df <- m[m$gene %in% differentiation, ]
diff_long <- rbind(
  data.frame(gene = diff_df$gene, log2FC = diff_df$log2FoldChange_d2, timepoint = "Day 2"),
  data.frame(gene = diff_df$gene, log2FC = diff_df$log2FoldChange_d7, timepoint = "Day 7")
)
diff_long$gene <- factor(diff_long$gene,
  levels = diff_df$gene[order(diff_df$log2FoldChange_d7)])

p_diff <- ggplot(diff_long, aes(x = gene, y = log2FC, fill = timepoint)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("Day 2" = "#4C78A8", "Day 7" = "#B03A2E")) +
  labs(title = "Myeloid-differentiation panel: VTP-50469 vs DMSO",
       subtitle = "Day 2 vs Day 7, MOLM13",
       x = NULL, y = "log2 fold change", fill = NULL) +
  theme_classic()
ggsave("figures/MOLM13_differentiation_panel_timecourse.png", p_diff,
       width = 7, height = 5, dpi = 300)

# ---- (C) Custom-signature NES, Day 2 vs Day 7 ----
custom_gsea$tp_lab <- ifelse(custom_gsea$timepoint == "d2", "Day 2", "Day 7")
custom_gsea$sig <- ifelse(custom_gsea$padj < 0.05, "padj < 0.05", "n.s.")
p_nes <- ggplot(custom_gsea,
  aes(x = pathway, y = NES, fill = tp_lab, alpha = sig)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("Day 2" = "#4C78A8", "Day 7" = "#B03A2E")) +
  scale_alpha_manual(values = c("padj < 0.05" = 1, "n.s." = 0.4)) +
  labs(title = "Curated signature enrichment (NES)",
       subtitle = "VTP-50469 vs DMSO, MOLM13",
       x = NULL, y = "Normalized enrichment score", fill = NULL, alpha = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
ggsave("figures/MOLM13_custom_signature_NES.png", p_nes,
       width = 7, height = 5, dpi = 300)

# ---- (D) Genome-wide log2FC, Day 2 vs Day 7, targets highlighted ----
m$class <- ifelse(m$gene %in% mll_targets, "MLL/KMT2A target",
            ifelse(m$gene %in% differentiation, "Differentiation", "Other"))
m_ord <- m[order(m$class != "Other"), ]  # draw highlighted points on top
p_scatter <- ggplot(m_ord, aes(log2FoldChange_d2, log2FoldChange_d7)) +
  geom_point(data = subset(m_ord, class == "Other"),
             color = "grey80", size = 0.6, alpha = 0.6) +
  geom_point(data = subset(m_ord, class != "Other"),
             aes(color = class), size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  scale_color_manual(values = c("MLL/KMT2A target" = "#B03A2E",
                                "Differentiation" = "#2E7D32")) +
  labs(title = "VTP-50469 response: Day 2 vs Day 7",
       x = "log2FC (Day 2)", y = "log2FC (Day 7)", color = NULL) +
  theme_classic()
ggsave("figures/MOLM13_timecourse_scatter.png", p_scatter,
       width = 6.5, height = 5.5, dpi = 300)

message("Done. Tables -> results/, figures -> figures/.")
