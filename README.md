# Menin Inhibitor Response in KMT2A-Rearranged Leukemia: Bulk RNA-seq Reanalysis

## Project overview

This project reanalyzes publicly available bulk RNA-seq data to investigate transcriptional changes associated with pharmacological inhibition of the menin–KMT2A interaction in KMT2A-rearranged leukemia cells.

**Research question:** How does treatment with the menin inhibitor VTP-50469 alter gene expression and Hallmark biological pathways in MOLM13 leukemia cells?

The analysis demonstrates a reproducible bulk RNA-seq workflow using GEO data, R, DESeq2, and gene set enrichment analysis (GSEA).

## Dataset and experimental design

- **GEO accession:** [GSE127506](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127506)
- **Study:** Krivtsov et al., *Cancer Cell* (2019)
- **Cell line:** MOLM13 (KMT2A-rearranged leukemia)
- **Timepoint:** Day 2
- **Treatment:** VTP-50469, 330 nM
- **Control:** DMSO vehicle
- **Biological replicates:** 3 per condition

The analysis uses six samples:

| Condition | Number of samples |
|---|---:|
| DMSO | 3 |
| VTP-50469 | 3 |

The original GEO count matrix contains 23,710 gene features. The six selected samples were extracted to create the analysis count matrix.

## Analysis workflow

The project uses GEO's processed gene-level count matrix rather than starting from raw sequencing reads.

1. **Data acquisition:** Download the processed count matrix and GEO sample metadata.
2. **Data preparation:** Select the six MOLM13 Day 2 samples and generate a count matrix and sample metadata file.
3. **Differential expression:** Use DESeq2 to compare VTP-50469 against DMSO, filtering features with fewer than 10 total counts.
4. **Exploratory analysis:** Examine sample relationships using PCA and sample correlations.
5. **Visualization:** Generate a volcano plot and a heatmap of the top differentially expressed genes.
6. **Pathway analysis:** Perform preranked GSEA using MSigDB Hallmark gene sets and the DESeq2 test statistic.

## Key findings

### Differential expression

At an adjusted p-value threshold of 0.05:

- **68 differentially expressed genes**
- **29 genes with higher expression** in VTP-50469-treated cells
- **39 genes with lower expression** in VTP-50469-treated cells

Several genes associated with leukemia biology showed expression changes, including MEIS1 and PBX3.

Positive log2 fold changes indicate higher expression in VTP-50469 relative to DMSO; negative values indicate lower expression.

### Hallmark gene set enrichment

Preranked GSEA identified **14 significantly enriched Hallmark pathways** at an adjusted p-value below 0.05.

Notable results included:

- TNFA signaling via NF-κB: positive enrichment
- Hypoxia: positive enrichment
- Apoptosis: positive enrichment
- IL6/JAK/STAT3 signaling: positive enrichment
- MYC targets V2: negative enrichment

Positive normalized enrichment scores (NES) indicate enrichment toward genes with higher expression in VTP-50469-treated cells; negative scores indicate enrichment toward genes with lower expression.

These results describe transcriptional enrichment patterns and do not independently establish pathway activation, apoptosis, or a specific mechanism of drug response.

## Figures

The `figures/` directory contains:

| Figure | Description |
|---|---|
| `GSE127506_MOLM13_Day2_PCA.png` | PCA of the six samples |
| `GSE127506_MOLM13_Day2_volcano.png` | Differential expression volcano plot |
| `GSE127506_MOLM13_Day2_heatmap.png` | Heatmap of the top differentially expressed genes |
| `MYC_TARGETS_V2_GSEA.png` | GSEA enrichment plot for MYC Targets V2 |
| `TNFA_NFKB_GSEA.png` | GSEA enrichment plot for TNFA signaling via NF-κB |

## Repository structure

```text
menin_rnaseq_project/
├── README.md
├── environment.yml
├── .gitignore
├── archive/
│   └── unused_workflows/
├── data/
│   ├── raw/
│   ├── metadata/
│   └── processed/
├── figures/
├── results/
└── scripts/
    ├── 01_download_GEO_counts.sh
    ├── 02_build_counts_matrix.R
    └── 03_deseq2_GSE127506.R
```

The `archive/` directory contains earlier experimental or unfinished workflows that are not part of the completed analysis.

## Reproducing the analysis

Run commands from the project root directory.

### 1. Download the GEO data

```bash
bash scripts/01_download_GEO_counts.sh
```

### 2. Prepare the count matrix

```bash
Rscript scripts/02_build_counts_matrix.R
```

This creates:

- `data/processed/GSE127506_MOLM13_Day2_counts.tsv`
- `data/metadata/GSE127506_MOLM13_Day2_samples.tsv`

### 3. Run differential expression and GSEA

```bash
Rscript scripts/03_deseq2_GSE127506.R
```

The analysis requires R and the packages used in the script, including DESeq2, apeglm, fgsea, msigdbr, ggplot2, ggrepel, and pheatmap.

## Limitations and interpretation

- This is a reanalysis of publicly available data for bioinformatics training and portfolio demonstration.
- The analysis focuses on one cell line, one treatment dose, and one timepoint.
- The study uses 3′-tag gene-level RNA-seq data; transcript-level interpretations are outside the scope of this analysis.
- The analysis uses processed GEO counts and does not include raw FASTQ alignment or independent read-level quality control.
- Gene set enrichment indicates coordinated expression patterns, not direct evidence of pathway function.
- Findings should be interpreted in the context of the original study and the limited sample size.

## Acknowledgments

Data are from GEO accession GSE127506, associated with Krivtsov et al., *Cancer Cell* (2019).

Reanalysis by Anagha Kaduskar.
