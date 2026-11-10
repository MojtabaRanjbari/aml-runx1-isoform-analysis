# Transcript-Isoform and Alternative Promoter Analysis of RUNX1 in AML

This repository contains R workflows developed to evaluate relative exon-bin usage, promoter switching (Distal P1 vs. Proximal P2), and somatic mutation associations for *RUNX1* across Acute Myeloid Leukemia (AML) patient cohorts.

---

## 🔬 Project Background
- **Institution:** Department of Hematology, Erasmus MC, Rotterdam, Netherlands
- **Focus:** Dissecting the transcriptional architecture of *RUNX1* isoforms (RUNX1a, RUNX1b, and RUNX1c) and their clinical/mutational correlates in AML.
- **Reference Coordinates:** Human genome GRCh37/hg19.

---

## 📊 Key Methodological Features
- **Normalization & Quantification:** Coverage-adjusted exon counts transformed into Reads Per Kilobase (RPK) and relative usage ratios across 11 DEXSeq bins.
- **Correlation Profiling:** Pairwise Pearson and Spearman correlation matrices across exon bins to distinguish promoter-specific expression clusters.
- **Statistical Testing & Multiple Testing Correction:** Two-sided Fisher's exact tests coupled with Benjamini–Hochberg False Discovery Rate (FDR) adjustments across 45 mutational categories (including *FLT3-ITD*, *CEBPA_DM*, *DNMT3A*, and *GATA2*).
- **Genomic Architecture Visualization:** Custom `ggplot2` compressed intron models and multi-track `Gviz` genomic visualizations comparing patient subgroups.

---

## 🛠️ Stack & Dependencies

**Genomic & Plotting Libraries:** `Gviz`, `GenomicRanges`, `ggplot2`, `gplots`, `scales`, `reshape2`, `dplyr`

---

## 📂 Repository Organization
```text
├── scripts/
│   └── runx1_isoform_pipeline.R   # Core data normalization, testing, and plotting code
├── docs/                          # Gene model schematics and PDF outputs
└── README.md                      # Project summary and documentation
