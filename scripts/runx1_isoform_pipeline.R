# ==============================================================================
# Pipeline: RUNX1 Exon Usage, Mutation Association & Promoter Analysis in AML
# ==============================================================================

# ----- Load Required Packages
library(reshape2)
library(dplyr)
library(ggplot2)
library(gplots)

# ----- Define Relative Paths
data_dir    <- "data"
results_dir <- "results"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# ==============================================================================
# 1. Read Exon Counts and Coverage Normalization
# ==============================================================================

dir_counts <- file.path(data_dir, "RUNX1_exoncount")
file_counts <- list.files(pattern = "\\.txt$", path = dir_counts)
rna <- as.data.frame(matrix(nrow = 32 * length(file_counts), ncol = 3))
colnames(rna) <- c("isoform", "count", "sample")

for (i in seq_along(file_counts)) {
  data <- read.delim(file.path(dir_counts, file_counts[i]), header = FALSE)
  colnames(data) <- c("isoform", "count")
  data$sample <- file_counts[i]
  rna[(i * 32 - 31):(i * 32), ] <- as.data.frame(data)
}

rna$sample <- gsub(".dex_counts_RUNX1.txt", "", rna$sample)
rna <- subset(rna, isoform != "chr21" & isoform != "chr8")
rna$count <- as.numeric(rna$count)

# Coverage normalization
dir_coverage <- file.path(data_dir, "dex_counts")
file_coverage <- list.files(pattern = "\\.dex_counts$", path = dir_coverage)
rna_coverage <- as.data.frame(matrix(nrow = 226070 * length(file_coverage), ncol = 3))
colnames(rna_coverage) <- c("isoform", "coverage", "sample")

for (i in seq_along(file_coverage)) {
  data <- read.delim(file.path(dir_coverage, file_coverage[i]), header = FALSE)
  colnames(data) <- c("isoform", "count")
  data$sample <- file_coverage[i]
  rna_coverage[(i * 226070 - 226069):(i * 226070), ] <- as.data.frame(data)
}

rna_coverage$sample <- gsub(".dex_counts", "", rna_coverage$sample)
rna_coverage <- aggregate(coverage ~ sample, rna_coverage, sum)
rna_coverage$coverage <- rna_coverage$coverage / min(rna_coverage$coverage)

rna <- merge(rna, rna_coverage, by = "sample")
rna$norm <- rna$count / (rna$coverage)

rna_runx1 <- reshape2::dcast(rna, isoform ~ sample, value.var = "count")
rna_runx1 <- rna_runx1[1:11, ]

# ==============================================================================
# 2. Coordinate Annotation & RPK / Ratio Normalization
# ==============================================================================

coord_file <- file.path(data_dir, "hg19_RUNX1_coordinates_2023-05-22.txt")
runx1_exon <- read.delim(coord_file, header = FALSE)
runx1_exon <- runx1_exon[2:12, c(1, 4:5)]
runx1_exon$exon_number <- rna_runx1$isoform
colnames(runx1_exon)[1:3] <- c("chr", "start", "end")
runx1_exon$bp <- runx1_exon$end - runx1_exon$start

# RPK calculation
rna_runx1_rpk <- as.data.frame(matrix(nrow = nrow(rna_runx1), ncol = ncol(rna_runx1)))
colnames(rna_runx1_rpk) <- colnames(rna_runx1)
rna_runx1_rpk[, 1] <- rna_runx1[, 1]
rna_runx1_rpk[, 2:ncol(rna_runx1)] <- apply(rna_runx1[, 2:ncol(rna_runx1)], 2, function(x) { (x / runx1_exon$bp) * 1000 })

# Exon ratio calculation
rna_runx1_ratio <- as.data.frame(matrix(nrow = nrow(rna_runx1_rpk), ncol = ncol(rna_runx1_rpk)))
colnames(rna_runx1_ratio) <- colnames(rna_runx1_rpk)
rna_runx1_ratio[, 1] <- rna_runx1_rpk[, 1]
rna_runx1_ratio[, 2:ncol(rna_runx1_ratio)] <- apply(rna_runx1_rpk[, 2:ncol(rna_runx1_ratio)], 2, function(x) { x / sum(x) })

# Prepare table for downstream analysis
runx1_ratio <- rna_runx1_ratio
rownames(runx1_ratio) <- runx1_ratio$isoform
runx1_ratio$isoform <- NULL
runx1_ratio <- as.data.frame(t(runx1_ratio))
rownames(runx1_ratio) <- sub("\\_.*", "", rownames(runx1_ratio))

# ==============================================================================
# 3. Correlation Analysis (Pearson & Spearman Heatmaps)
# ==============================================================================

runx1_pearson <- cor(runx1_ratio, method = "pearson")
runx1_spearman <- cor(runx1_ratio, method = "spearman")

colheatmap <- colorRampPalette(c("blue", "gray95", "red"))
colors <- seq(-1, 1, length.out = 201)
side_cols <- c("black", "blue", "cornflowerblue", "cyan", "aquamarine3", "yellow", 
               "gold", "chocolate1", "darkorchid1", "lightpink1", "red")

pdf(file.path(results_dir, "RUNX1_11exon_correlation.pdf"), width = 9, height = 9)

# Pearson without dendrogram
heatmap.2(as.matrix(runx1_pearson), trace = "none", margins = c(8, 8), Rowv = FALSE, Colv = FALSE, 
          dendrogram = "none", col = colheatmap, breaks = colors, 
          ColSideColors = side_cols, RowSideColors = side_cols, 
          cellnote = round(runx1_pearson, digits = 2), notecol = "black",
          main = "Pearson (No Dendrogram)")

# Pearson with dendrogram
heatmap.2(as.matrix(runx1_pearson), trace = "none", margins = c(8, 8), Rowv = TRUE, Colv = TRUE, 
          dendrogram = "both", col = colheatmap, breaks = colors, 
          hclust = function(x) hclust(x, method = "average"), 
          distfun = function(x) as.dist((1 - cor(t(x)) / 2)), 
          ColSideColors = side_cols, RowSideColors = side_cols, 
          cellnote = round(runx1_pearson, digits = 2), notecol = "black",
          main = "Pearson (Clustered)")

# Spearman without dendrogram
heatmap.2(as.matrix(runx1_spearman), trace = "none", margins = c(8, 8), Rowv = FALSE, Colv = FALSE, 
          dendrogram = "none", col = colheatmap, breaks = colors, 
          ColSideColors = side_cols, RowSideColors = side_cols, 
          cellnote = round(runx1_spearman, digits = 2), notecol = "black",
          main = "Spearman (No Dendrogram)")

# Spearman with dendrogram
heatmap.2(as.matrix(runx1_spearman), trace = "none", margins = c(8, 8), Rowv = TRUE, Colv = TRUE, 
          dendrogram = "both", col = colheatmap, breaks = colors, 
          hclust = function(x) hclust(x, method = "average"), 
          distfun = function(x) as.dist((1 - cor(t(x)) / 2)), 
          ColSideColors = side_cols, RowSideColors = side_cols, 
          cellnote = round(runx1_spearman, digits = 2), notecol = "black",
          main = "Spearman (Clustered)")

dev.off()

# ==============================================================================
# 4. Exon Correlation Scatter Plot (Exon 8 vs Exon 9)
# ==============================================================================

# Load mutation data & merge
mutation_file <- file.path(data_dir, "mut_data_20230523.txt")
mutation <- read.delim(mutation_file, header = TRUE)
merge_RM <- merge(runx1_ratio, mutation, by.x = "row.names", by.y = "patient")
merge_RM <- merge_RM[order(merge_RM$`RUNX1:010`), ]
colnames(merge_RM)[1] <- "patient"

theme_dot <- theme(panel.border = element_blank(), axis.line = element_line(),
                   legend.title = element_blank(),
                   legend.text = element_text(colour = "black", size = 12),
                   legend.position = "none",
                   axis.text = element_text(colour = "black", face = "bold"),
                   axis.line.x = element_line(linewidth = 1),
                   axis.line.y = element_line(linewidth = 1),
                   axis.ticks = element_line(colour = "black"),
                   axis.ticks.x = element_line(linewidth = 1),
                   axis.ticks.y = element_line(linewidth = 1),
                   axis.ticks.length = unit(0.5, "cm"),
                   strip.background = element_blank(),
                   strip.text.x = element_text(colour = "black", size = 12, face = "bold"),
                   strip.text.y = element_text(colour = "black", size = 12, face = "bold"),
                   axis.title.x = element_text(colour = "black", size = 12, face = "bold"),
                   axis.title.y = element_text(colour = "black", size = 12, face = "bold"),
                   axis.text.x = element_text(size = 12),
                   axis.text.y = element_text(size = 12))

p_scatter <- ggplot(merge_RM, aes(`RUNX1:008`, `RUNX1:009`)) + 
  geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ x) + 
  xlab("Fraction of read (exon 4, RUNX1A/B)") + 
  ylab("Fraction of read (exon 3, RUNX1C)") + 
  theme_bw() + 
  theme_dot

ggsave(file.path(results_dir, "RUNX1_Exon8_vs_Exon9_Scatter.pdf"), p_scatter, width = 7, height = 6)

# ==============================================================================
# 5. Fisher's Exact Tests Across 45 Mutations and FDR Calculation
# ==============================================================================

fisher_writer <- function(runxcol) {
  fisher_table <- as.data.frame(matrix(nrow = 4 * 45, ncol = 5))
  colnames(fisher_table) <- c("expression", "status", "isoform", "mutation", "count")
  fisher_table$expression <- rep(rep(c("low", "high"), each = 2), 45)
  fisher_table$status <- rep(rep(c("wt", "mut"), 2), 45)
  fisher_table$isoform <- runxcol
  
  for (i in 1:45) {
    median_val <- median(merge_RM[[runxcol]], na.rm = TRUE)
    mut_col <- colnames(merge_RM)[12 + i]
    idx <- (i - 1) * 4 + (1:4)
    fisher_table[idx, 4] <- mut_col
    
    fisher_table[idx[1], 5] <- table(subset(merge_RM, merge_RM[[runxcol]] < median_val)[[mut_col]])[1]
    fisher_table[idx[2], 5] <- table(subset(merge_RM, merge_RM[[runxcol]] < median_val)[[mut_col]])[2]
    fisher_table[idx[3], 5] <- table(subset(merge_RM, merge_RM[[runxcol]] > median_val)[[mut_col]])[1]
    fisher_table[idx[4], 5] <- table(subset(merge_RM, merge_RM[[runxcol]] > median_val)[[mut_col]])[2]
  }
  fisher_table$count[is.na(fisher_table$count)] <- 0
  return(fisher_table)
}

exon_names <- paste0("RUNX1:0", sprintf("%02d", 1:11))
RUNX1_all <- lapply(exon_names, fisher_writer)
names(RUNX1_all) <- paste0("RUNX1.", sprintf("%03d", 1:11))

mutation_list <- RUNX1_all[[1]][!duplicated(RUNX1_all[[1]]$mutation), 4]

fisher_calculator <- function(df, gene) {
  df_sub <- subset(df, mutation == gene)
  mat <- matrix(df_sub$count, nrow = 2, byrow = FALSE)
  fisher.test(mat, alternative = "two.sided")$p.value
}

exon_mutation <- as.data.frame(matrix(nrow = 11, ncol = 45))
colnames(exon_mutation) <- mutation_list
rownames(exon_mutation) <- colnames(merge_RM)[2:12]

for (i in 1:11) {
  for (j in 1:45) {
    exon_mutation[i, j] <- fisher_calculator(RUNX1_all[[i]], mutation_list[j])
  }
}

exon_mutation_long <- exon_mutation
exon_mutation_long$isoform <- rownames(exon_mutation_long)
exon_mutation_long <- exon_mutation_long[, c(46, 1:45)]
exon_mutation_long <- reshape2::melt(exon_mutation_long)
colnames(exon_mutation_long) <- c("isoform", "mutation", "p.value")
exon_mutation_long$fdr <- p.adjust(exon_mutation_long$p.value, method = "BH", n = nrow(exon_mutation_long))

# ==============================================================================
# 6. Heatmaps for Significance (P-value & FDR with Stars)
# ==============================================================================

p_split <- function(df, n, newcol) {
  df0 <- subset(df, df[, n] >= 0.05)
  df1 <- subset(df, df[, n] >= 0.01 & df[, n] < 0.05)
  df2 <- subset(df, df[, n] >= 0.001 & df[, n] < 0.01)
  df3 <- subset(df, df[, n] < 0.001)
  
  df0[[newcol]] <- ""
  df1[[newcol]] <- "*"
  df2[[newcol]] <- "**"
  df3[[newcol]] <- "***"
  
  rbind(df0, df1, df2, df3)
}

exon_mutation_long <- p_split(exon_mutation_long, 3, "p_summary")
exon_mutation_long <- p_split(exon_mutation_long, 4, "fdr_summary")

exon_pv <- reshape2::dcast(exon_mutation_long, mutation ~ isoform, value.var = "p.value")
exon_ps <- reshape2::dcast(exon_mutation_long, mutation ~ isoform, value.var = "p_summary")

exon_fv <- reshape2::dcast(exon_mutation_long, mutation ~ isoform, value.var = "fdr")
exon_fs <- reshape2::dcast(exon_mutation_long, mutation ~ isoform, value.var = "fdr_summary")

colheatmap.p <- colorRampPalette(c("gray95", "black"))
colors.p <- seq(0, 15, length.out = 201)

colheatmap.f <- colorRampPalette(c("gray95", "black"))
colors.f <- seq(0, 12, length.out = 201)

pdf(file.path(results_dir, "RUNX1_Significance_Heatmaps.pdf"), width = 10, height = 10)

# P-value Heatmap
heatmap.2(as.matrix(-log10(exon_pv[, 2:12])), trace = "none", margins = c(8, 8), 
          Rowv = TRUE, Colv = FALSE, dendrogram = "row", breaks = colors.p, 
          col = colheatmap.p, ColSideColors = side_cols, 
          cellnote = exon_ps[, 2:12], notecol = "cornflowerblue", notecex = 1.5, 
          labRow = as.character(exon_pv$mutation), main = "-log10(p-value) by Mutation")

# FDR Heatmap
heatmap.2(as.matrix(-log10(exon_fv[, 2:12])), trace = "none", margins = c(8, 8), 
          Rowv = TRUE, Colv = FALSE, dendrogram = "row", breaks = colors.f, 
          col = colheatmap.f, ColSideColors = side_cols, 
          cellnote = exon_fs[, 2:12], notecol = "cornflowerblue", notecex = 1.5, 
          labRow = as.character(exon_fv$mutation), main = "-log10(FDR) by Mutation")

dev.off()

# ==============================================================================
# 7. Individual Mutation Bar Plots
# ==============================================================================

sel_muts <- c("RAS_all", "FLT3_all", "NPM1", "CEBPA_DM", "DNMT3A")

example_all <- rbind(
  cbind(subset(RUNX1_all[["RUNX1.001"]], mutation %in% sel_muts), exon = "exon.001"),
  cbind(subset(RUNX1_all[["RUNX1.004"]], mutation %in% sel_muts), exon = "exon.004"),
  cbind(subset(RUNX1_all[["RUNX1.005"]], mutation %in% sel_muts), exon = "exon.005"),
  cbind(subset(RUNX1_all[["RUNX1.008"]], mutation %in% sel_muts), exon = "exon.008"),
  cbind(subset(RUNX1_all[["RUNX1.009"]], mutation %in% sel_muts), exon = "exon.009"),
  cbind(subset(RUNX1_all[["RUNX1.010"]], mutation %in% sel_muts), exon = "exon.010")
)

example_all$mutation <- gsub("RAS_all", "1.RAS_all", example_all$mutation)
example_all$mutation <- gsub("CEBPA_DM", "2.CEBPA_DM", example_all$mutation)
example_all$mutation <- gsub("FLT3_all", "3.FLT3_all", example_all$mutation)
example_all$mutation <- gsub("NPM1", "4.NPM1", example_all$mutation)
example_all$mutation <- gsub("DNMT3A", "5.DNMT3A", example_all$mutation)

theme_ggplot <- theme(panel.border = element_blank(), axis.line = element_line(),
                      legend.title = element_blank(),
                      legend.text = element_text(colour = "black", size = 12),
                      legend.position = "none",
                      axis.text = element_text(colour = "black", face = "bold"),
                      axis.line.x = element_line(linewidth = 1),
                      axis.line.y = element_line(linewidth = 1),
                      axis.ticks = element_line(colour = "black"),
                      axis.ticks.x = element_line(linewidth = 1),
                      axis.ticks.y = element_line(linewidth = 1),
                      axis.ticks.length = unit(0.5, "cm"),
                      strip.background = element_blank(),
                      strip.text.x = element_text(colour = "black", size = 12, face = "bold"),
                      strip.text.y = element_text(colour = "black", size = 12, face = "bold"),
                      axis.title.x = element_blank(),
                      axis.title.y = element_text(colour = "black", size = 12, face = "bold"),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 12))

pdf(file.path(results_dir, "RUNX1_example_ExonMutation.pdf"), width = 9, height = 6)

# Exons 1, 4, 5
p_bar1 <- ggplot(subset(example_all, exon %in% c("exon.001", "exon.004", "exon.005")), 
                 aes(expression, count, fill = status)) + 
  facet_grid(rows = vars(exon), cols = vars(mutation)) + 
  geom_col() + 
  scale_fill_manual(values = c("wt" = "aquamarine3", "mut" = "firebrick3")) + 
  theme_bw() + 
  theme_ggplot
print(p_bar1)

# Exons 8, 9, 10
p_bar2 <- ggplot(subset(example_all, exon %in% c("exon.008", "exon.009", "exon.010")), 
                 aes(expression, count, fill = status)) + 
  facet_grid(rows = vars(exon), cols = vars(mutation)) + 
  geom_col() + 
  ylab("number of patients") + 
  scale_fill_manual(values = c("wt" = "aquamarine3", "mut" = "firebrick3")) + 
  theme_bw() + 
  theme_ggplot
print(p_bar2)

dev.off()

# ==============================================================================
# 8. Exon 8 & 9 vs Mutation Distribution Heatmaps (rex_89)
# ==============================================================================

rex_89 <- merge_RM[, c(16:17, 21, 24)]
rex_89 <- rex_89[order(rex_89$`RUNX1:008`), ]

colheatmap.exon <- colorRampPalette(c("blue", "gray95", "orange"))
colheatmap.mut  <- c("black", "red")
colors.exon     <- seq(-4, 4, length.out = 201)

pdf(file.path(results_dir, "RUNX1_Exon8_Exon9_Mutation_Heatmaps.pdf"), width = 8, height = 8)

# Exon fraction heatmap
heatmap.2(as.matrix(rex_89[, 1:2]), trace = "none", margins = c(8, 8), 
          Rowv = FALSE, Colv = FALSE, dendrogram = "none", 
          col = colheatmap.exon, breaks = colors.exon, scale = "column",
          main = "Exon 8 & 9 Fractions")

# Mutation status heatmap
heatmap.2(as.matrix(rex_89[, 3:4]), trace = "none", margins = c(8, 8), 
          Rowv = FALSE, Colv = FALSE, dendrogram = "none", 
          col = colheatmap.mut, scale = "column",
          main = "Mutation Status (rex_89)")

dev.off()
