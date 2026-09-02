# ==============================================================================
# Pipeline: RUNX1 Transcript-Isoform and Promoter Usage Analysis in AML
# Description: Normalization (RPK), correlation analysis, Fisher's exact tests 
#              with FDR correction, and genomic visualizations (Gviz & ggplot2).
# ==============================================================================

# ----- Load Required Packages
library(reshape2)
library(dplyr)
library(ggplot2)
library(gplots)
library(scales)
library(Gviz)
library(GenomicRanges)

# ----- Setup Project Directories (Relative paths for reproducibility)
data_dir    <- "data"
results_dir <- "results"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# ==============================================================================
# 1. Data Ingestion & Coverage Normalization
# ==============================================================================

# Define input subdirectories
exon_counts_dir <- file.path(data_dir, "RUNX1_exoncount")
coverage_dir    <- file.path(data_dir, "dex_counts")

# Process exon count tables across patient samples
files <- list.files(path = exon_counts_dir, pattern = "\\.txt$", full.names = FALSE)
table_list <- vector("list", length(files))

for (i in seq_along(files)) {
  r_data <- read.delim(file.path(exon_counts_dir, files[i]), header = FALSE)
  colnames(r_data) <- c("isoform", "count")
  r_data$sample <- gsub(".dex_counts_RUNX1.txt", "", files[i])
  table_list[[i]] <- r_data
}
table <- do.call(rbind, table_list)
table <- subset(table, isoform != "chr21" & isoform != "chr8")
table$count <- as.numeric(table$count)

# Coverage calculation per sample
files_coverage <- list.files(path = coverage_dir, pattern = "\\.dex_counts$", full.names = FALSE)
coverage_list <- vector("list", length(files_coverage))

for (i in seq_along(files_coverage)) {
  rc_data <- read.delim(file.path(coverage_dir, files_coverage[i]), header = FALSE)
  colnames(rc_data) <- c("isoform", "count")
  rc_data$sample <- gsub(".dex_counts", "", files_coverage[i])
  coverage_list[[i]] <- rc_data
}
table_coverage <- do.call(rbind, coverage_list)
table_coverage <- aggregate(count ~ sample, table_coverage, sum)
colnames(table_coverage)[2] <- "coverage"
table_coverage$coverage <- table_coverage$coverage / min(table_coverage$coverage)

# Normalization across sequencing depth
table <- merge(table, table_coverage, by = "sample")
table$norm <- table$count / table$coverage

table_RUNX1 <- reshape2::dcast(table, isoform ~ sample, value.var = "norm")
table_RUNX1 <- table_RUNX1[1:11, ]

# ==============================================================================
# 2. Exon Coordinate Mapping & RPK Calculation
# ==============================================================================

coord_file <- file.path(data_dir, "hg19_RUNX1_coordinates.txt")
runx1_exon <- read.delim(coord_file, header = FALSE)
runx1_exon <- runx1_exon[2:12, c(1, 4, 5)]
colnames(runx1_exon) <- c("chr", "start", "end")
runx1_exon$exon_number <- table_RUNX1$isoform
runx1_exon$length_bp <- runx1_exon$end - runx1_exon$start

# Reads Per Kilobase (RPK)
table_RUNX1_rpk <- table_RUNX1
table_RUNX1_rpk[, 2:ncol(table_RUNX1)] <- apply(
  table_RUNX1[, 2:ncol(table_RUNX1)], 2, 
  function(x) { (x / runx1_exon$length_bp) * 1000 }
)

# Relative Usage Fraction (Ratio Table)
table_RUNX1_ratio <- table_RUNX1_rpk
table_RUNX1_ratio[, 2:ncol(table_RUNX1_ratio)] <- apply(
  table_RUNX1_rpk[, 2:ncol(table_RUNX1_rpk)], 2, 
  function(x) { x / sum(x) }
)

RUNX1_ratio <- table_RUNX1_ratio
rownames(RUNX1_ratio) <- RUNX1_ratio$isoform
RUNX1_ratio$isoform <- NULL
RUNX1_ratio <- as.data.frame(t(RUNX1_ratio))
rownames(RUNX1_ratio) <- sub("\\_.*", "", rownames(RUNX1_ratio))

# ==============================================================================
# 3. Correlation Profiling Across Exon Bins
# ==============================================================================

RUNX1_pearson  <- cor(RUNX1_ratio, method = "pearson")
RUNX1_spearman <- cor(RUNX1_ratio, method = "spearman")

colheatmap <- colorRampPalette(c("blue", "gray95", "red"))
colors <- seq(-1, 1, length.out = 201)

pdf(file.path(results_dir, "RUNX1_exon_correlation_pearson.pdf"), width = 9, height = 9)
heatmap.2(as.matrix(RUNX1_pearson), trace = "none", margins = c(8, 8),
          col = colheatmap, breaks = colors,
          main = "RUNX1 Exon-Bin Pearson Correlation")
dev.off()

# ==============================================================================
# 4. Mutation Association & Fisher's Exact Tests with FDR
# ==============================================================================

mutation_file <- file.path(data_dir, "mut_data.txt")
mutation <- read.delim(mutation_file, header = TRUE)
merge_RM <- merge(RUNX1_ratio, mutation, by.x = "row.names", by.y = "patient")
colnames(merge_RM)[1] <- "patient"

# Fisher's test generator across median split
fisher_writer <- function(x) {
  num_muts <- 45
  fisher_table <- as.data.frame(matrix(nrow = 4 * num_muts, ncol = 5))
  colnames(fisher_table) <- c("expression", "status", "isoform", "mutation", "count")
  fisher_table$expression <- rep(rep(c("low", "high"), each = 2), num_muts)
  fisher_table$status     <- rep(rep(c("wt", "mut"), 2), num_muts)
  fisher_table$isoform    <- x
  
  median_val <- median(merge_RM[[x]], na.rm = TRUE)
  
  for (i in 1:num_muts) {
    mut_col <- colnames(merge_RM)[12 + i]
    idx <- (i - 1) * 4 + (1:4)
    fisher_table[idx, "mutation"] <- mut_col
    
    fisher_table[idx[1], "count"] <- sum(merge_RM[[x]] < median_val & merge_RM[[mut_col]] == 0, na.rm = TRUE)
    fisher_table[idx[2], "count"] <- sum(merge_RM[[x]] < median_val & merge_RM[[mut_col]] == 1, na.rm = TRUE)
    fisher_table[idx[3], "count"] <- sum(merge_RM[[x]] > median_val & merge_RM[[mut_col]] == 0, na.rm = TRUE)
    fisher_table[idx[4], "count"] <- sum(merge_RM[[x]] > median_val & merge_RM[[mut_col]] == 1, na.rm = TRUE)
  }
  fisher_table$count[is.na(fisher_table$count)] <- 0
  return(fisher_table)
}

exon_ids <- paste0("RUNX1:0", sprintf("%02d", 1:11))
names(exon_ids) <- paste0("RUNX1.", sprintf("%03d", 1:11))
RUNX1_all <- lapply(exon_ids, fisher_writer)

fisher_calculator <- function(df, gene) {
  df_sub <- subset(df, mutation == gene)
  mat <- matrix(df_sub$count, nrow = 2, byrow = FALSE)
  fisher.test(mat, alternative = "two.sided")$p.value
}

mutation_list <- unique(RUNX1_all[[1]]$mutation)
exon_mutation <- as.data.frame(matrix(nrow = 11, ncol = length(mutation_list)))
colnames(exon_mutation) <- mutation_list
rownames(exon_mutation) <- exon_ids

for (i in 1:11) {
  for (j in seq_along(mutation_list)) {
    exon_mutation[i, j] <- fisher_calculator(RUNX1_all[[i]], mutation_list[j])
  }
}

exon_mutation_long <- reshape2::melt(as.matrix(exon_mutation))
colnames(exon_mutation_long) <- c("isoform", "mutation", "p.value")
exon_mutation_long$fdr <- p.adjust(exon_mutation_long$p.value, method = "BH")

write.csv(exon_mutation_long, file.path(results_dir, "RUNX1_mutation_associations_fdr.csv"), row.names = FALSE)

# ==============================================================================
# 5. Promoter Balance Quantification (P1 vs P2)
# ==============================================================================

merge_RM$P1_score <- rowSums(merge_RM[, c("RUNX1:009", "RUNX1:010", "RUNX1:011")], na.rm = TRUE)
merge_RM$P2_score <- merge_RM[["RUNX1:008"]]
merge_RM$P1_P2_ratio <- log2((merge_RM$P1_score + 1e-6) / (merge_RM$P2_score + 1e-6))

# Boxplot function for specific mutations
plot_promoter_by_mut <- function(df, mut_col) {
  df_sub <- df[!is.na(df[[mut_col]]), ]
  df_sub$status <- factor(ifelse(df_sub[[mut_col]] == 1, "mut", "wt"), levels = c("wt", "mut"))
  
  ggplot(df_sub, aes(x = status, y = P1_P2_ratio, fill = status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
    scale_fill_manual(values = c("wt" = "cornflowerblue", "mut" = "firebrick3")) +
    theme_bw() +
    labs(title = paste("RUNX1 P1/P2 Balance by", mut_col),
         x = paste(mut_col, "Status"),
         y = "log2(P1 / P2 Ratio)")
}

p_cebpa <- plot_promoter_by_mut(merge_RM, "CEBPA_DM")
ggsave(file.path(results_dir, "RUNX1_P1_P2_CEBPA_DM.png"), p_cebpa, width = 6, height = 5, dpi = 300)

p_flt3 <- plot_promoter_by_mut(merge_RM, "FLT3_ITD")
ggsave(file.path(results_dir, "RUNX1_P1_P2_FLT3_ITD.png"), p_flt3, width = 6, height = 5, dpi = 300)
