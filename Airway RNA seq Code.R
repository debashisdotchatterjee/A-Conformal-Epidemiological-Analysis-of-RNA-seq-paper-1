# ==============================================================================
# Title: Dexamethasone Exposure in Airway Smooth Muscle Cells
# Methodology: Distribution-Free Framework (Conformal Inference Context)
# Dataset: Bioconductor 'airway'
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Setup and Library Installation
# ------------------------------------------------------------------------------
# Note: Ensure Bioconductor is installed:
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("airway", "DESeq2", "SummarizedExperiment"))
# install.packages(c("ggplot2", "dplyr", "tidyr", "pheatmap", "viridis", "RColorBrewer"))

library(airway)
library(DESeq2)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(gridExtra) # Added for grid arrangement if needed, though mostly using base/ggplot

# Set seed for reproducibility (crucial for Randomized Conformal P-values)
set.seed(42)

# ------------------------------------------------------------------------------
# 1. Dataset Loading and Application (Section 1.1)
# ------------------------------------------------------------------------------
data("airway")
se <- airway

# Inspect design: 8 samples, 4 cell lines, Treated (Dex) vs Untreated
print("Experimental Design (Table 1 Context):")
print(colData(se)[, c("cell", "dex")])

# ------------------------------------------------------------------------------
# 2. Data Preprocessing (Section 1.2)
# ------------------------------------------------------------------------------

# Create DESeqDataSet (The Working Model Framework)
# Design controls for 'cell' (baseline Z) and tests 'dex' (Exposure E)
dds <- DESeqDataSet(se, design = ~ cell + dex)

# A. Filtering
# "Filtered out genes with low overall expression (fewer than 10 total counts)"
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(paste0("\nRemaining Genes (G): ", nrow(dds), "\n"))

# B. Normalization
# "Library size normalization ... using the median-of-ratios method"
dds <- estimateSizeFactors(dds)
print("Size Factors (s_i):")
print(sizeFactors(dds))

# ------------------------------------------------------------------------------
# 3. Exploratory Visualization (The "Paired" Design)
# ------------------------------------------------------------------------------

# Variance Stabilizing Transformation for PCA
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot 1: PCA Plot showing Paired Design & Exposure Effect
p1 <- ggplot(pcaData, aes(PC1, PC2, color = cell, shape = dex)) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_viridis_d(option = "turbo", name = "Cell Line (Z)") +
  scale_shape_manual(values = c(16, 17), name = "Exposure (E)") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ggtitle("Principal Component Analysis: Airway Dataset",
          subtitle = "Demonstrating Paired Design (Cell Line) and Treatment Effect")

print(p1)

# ------------------------------------------------------------------------------
# 4. Working Model Specification & Non-Conformity Scoring (Section 2.2 & 2.3)
# ------------------------------------------------------------------------------

# Fit the GLM (Negative Binomial)
# This provides estimates for mu (expected count) and dispersion
dds <- DESeq(dds)

# Extract fitted values (mu_ig) and Dispersions
mu_hat <- assays(dds)[["mu"]]
phi_hat <- dispersions(dds) # Vector of dispersions per gene

# Replicate phi_hat for matrix math
phi_matrix <- matrix(rep(phi_hat, ncol(dds)), ncol = ncol(dds))

# Calculate Variance according to NB Model: Var = mu + phi * mu^2
sigma_sq_hat <- mu_hat + (phi_matrix * mu_hat^2)
sigma_hat <- sqrt(sigma_sq_hat)

# --- Calculating Non-Conformity Scores (Equation 3) ---
# V_ig = |Y_ig - mu_hat_ig| / sigma_hat_ig
Y_ig <- counts(dds, normalized = TRUE) # Using normalized counts for scale consistency
V_ig <- abs(Y_ig - mu_hat) / sigma_hat

# Reshape for plotting
v_df <- as.data.frame(V_ig)
v_df$Gene <- rownames(v_df)
v_long <- pivot_longer(v_df, cols = -Gene, names_to = "Sample", values_to = "Score")

# Plot 2: Distribution of Non-Conformity Scores
# This validates the "surprisingness" of data given the working model
p2 <- ggplot(v_long, aes(x = Score, fill = Sample)) +
  geom_density(alpha = 0.4) +
  scale_fill_viridis_d(option = "plasma") +
  xlim(0, 5) + # Focus on the bulk of the distribution
  theme_classic(base_size = 14) +
  labs(title = "Distribution of Non-Conformity Scores (V_ig)",
       subtitle = "Quantifying 'Surprise' relative to the Working Model (Eq. 3)",
       x = "Non-Conformity Score",
       y = "Density")

print(p2)

# ------------------------------------------------------------------------------
# 5. Inference and Results (Section 2.4)
# ------------------------------------------------------------------------------
# Extract results for Dexamethasone treatment
# We use standard DESeq2 p-values here as a proxy for the 'Working Model'
# in a real conformal loop, these would be the randomized p-values.
res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Convert to dataframe for plotting
res_df <- as.data.frame(res)
res_df <- na.omit(res_df) # Remove NA p-values

# Add Classification categories
res_df$diffexpressed <- "NO"
res_df$diffexpressed[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "UP"
res_df$diffexpressed[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "DOWN"

# Plot 3: Colorful Volcano Plot
p3 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "#440154FF", "NO" = "grey80", "UP" = "#FDE725FF"),
                     name = "Response") +
  labs(title = "Volcano Plot: Dexamethasone Response",
       subtitle = "Differential Expression (Working Model Results)",
       x = "Log2 Fold Change (Exposure Effect)",
       y = "-Log10 Adjusted P-value") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

print(p3)

# Plot 4: P-value Histogram (To check uniformity under Null)
p4 <- ggplot(res_df, aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "#3B528BFF", color = "white") +
  theme_minimal() +
  labs(title = "Histogram of P-values",
       subtitle = "Anticipating Uniformity under Null + Signal Spike near 0",
       x = "P-value", y = "Frequency")

print(p4)

# ------------------------------------------------------------------------------
# 6. Heatmap of Top Genes (Visualizing Y_ig)
# ------------------------------------------------------------------------------

# Select top 20 genes by adjusted p-value
top_genes <- head(order(res$padj), 20)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat) # Center the rows (Z-score like)

# Annotations for columns
df_col <- as.data.frame(colData(dds)[,c("cell", "dex")])

# Plot 5: Heatmap
# Using a diverging palette for expression differences
pheatmap(mat, 
         annotation_col = df_col,
         main = "Top 20 Differentially Expressed Genes",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize = 10)

# ------------------------------------------------------------------------------
# 7. Summary Tables
# ------------------------------------------------------------------------------

# Top 10 Genes Table
top_10_table <- res_df %>%
  arrange(padj) %>%
  head(10) %>%
  select(baseMean, log2FoldChange, pvalue, padj)

print("Top 10 Differentially Expressed Genes:")
print(top_10_table)

# Counts of Significant Genes
sig_summary <- table(res_df$diffexpressed)
print("Summary of Differential Expression:")
print(sig_summary)

# ------------------------------------------------------------------------------
# 8. Extended Visualizations & Diagnostics (New Additions)
# ------------------------------------------------------------------------------

# --- Plot 6: MA Plot (Log Fold Change vs Mean Expression) ---
# Highlights where the significant genes lie across the expression range
p6 <- ggplot(res_df, aes(x = log10(baseMean), y = log2FoldChange, color = diffexpressed)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_hline(yintercept = 0, color = "red", linetype = "solid", size=1) +
  scale_color_manual(values = c("DOWN" = "#21908CFF", "NO" = "grey90", "UP" = "#FDE725FF"), name="Sig") +
  theme_dark(base_size = 14) +
  labs(title = "MA Plot: Expression Level vs. Fold Change",
       subtitle = "Checking for bias in fold changes across expression levels",
       x = "Log10 Mean Expression",
       y = "Log2 Fold Change")

print(p6)

# --- Plot 7: Sample-to-Sample Distance Heatmap ---
# Confirms biological replicates cluster together (QC)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$cell, vsd$dex, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample-to-Sample Euclidean Distances")

# --- Plot 8: Paired Design Visualization (Top Gene) ---
# Directly showing the "Consistency" and "Paired" nature from Section 1.1
top_gene_id <- rownames(res_df)[which.min(res_df$padj)]
gene_counts <- plotCounts(dds, gene = top_gene_id, intgroup = c("dex", "cell"), returnData = TRUE)

p8 <- ggplot(gene_counts, aes(x = dex, y = count, color = cell, group = cell)) +
  geom_point(size = 4) +
  geom_line(size = 1, linetype = "dashed") + # Connects paired samples
  scale_y_log10() + 
  scale_color_viridis_d(option="turbo", name="Cell Line") +
  theme_bw(base_size=14) + 
  labs(title = paste0("Expression of Top Hit: ", top_gene_id),
       subtitle = "Dashed lines connect paired samples (Cell Line Control)",
       y = "Normalized Counts (Log Scale)", x = "Treatment")

print(p8)

# --- Plot 9: Dispersion Estimates (Working Model Parameters) ---
# Visualizing the phi_hat parameters from Section 2.2
# We extract raw data from DESeq object to plot manually for full control
disp_df <- data.frame(
  baseMean = mcols(dds)$baseMean,
  dispGeneEst = mcols(dds)$dispGeneEst,
  dispersion = dispersions(dds)
)
disp_df <- disp_df[disp_df$baseMean > 0 & !is.na(disp_df$dispGeneEst),]

p9 <- ggplot(disp_df, aes(x=baseMean, y=dispGeneEst)) +
  geom_point(alpha=0.1, color="black", size=0.5) +
  geom_point(aes(y=dispersion), color="red", size=0.5, alpha=0.2) + # The fitted trend
  scale_x_log10() + scale_y_log10() +
  theme_light() +
  labs(title = "Dispersion Estimates (Working Model)",
       subtitle = "Black: Gene-wise estimates | Red: Final fitted estimates",
       x = "Mean Expression", y = "Dispersion")

print(p9)

# --- Plot 10: Conformal Diagnostic (Score Homoscedasticity) ---
# Checking if Non-Conformity Scores (V_ig) depend on expression magnitude
# Calculate mean V_ig per gene
v_mean_per_gene <- rowMeans(V_ig)
diag_df <- data.frame(
  MeanExpr = mu_hat[,1], # Just taking first col as proxy for mean magnitude
  MeanScore = v_mean_per_gene
)

p10 <- ggplot(diag_df, aes(x=MeanExpr, y=MeanScore)) +
  geom_point(alpha=0.2, color="#440154FF") +
  geom_smooth(method="gam", color="#FDE725FF") + # Add trend line
  scale_x_log10() +
  ylim(0, 3) +
  theme_classic() +
  labs(title = "Conformal Diagnostic: Score vs Expression",
       subtitle = "Checking if V_ig is independent of signal intensity (Homoscedasticity)",
       x = "Predicted Mean Expression (Log Scale)",
       y = "Mean Non-Conformity Score")

print(p10)

# ------------------------------------------------------------------------------
# 9. Advanced Summary Tables
# ------------------------------------------------------------------------------

# Create a rich summary table of top significant genes with StdErr and Stat
rich_table <- res_df %>%
  arrange(padj) %>%
  head(15) %>%
  select(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  mutate(
    baseMean = round(baseMean, 1),
    log2FoldChange = round(log2FoldChange, 3),
    lfcSE = round(lfcSE, 3),
    stat = round(stat, 2),
    Significance = ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**", "*"))
  )

print("--- Advanced Top 15 Genes Table ---")
print(rich_table)

# ------------------------------------------------------------------------------
# 10. Even More Visualizations (Deep Dive)
# ------------------------------------------------------------------------------

# --- Plot 11: Violin Plot of Normalized Counts ---
# Checks if normalization harmonized the distributions across samples
norm_counts <- counts(dds, normalized = TRUE)
norm_long <- as.data.frame(norm_counts) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Count") %>%
  mutate(LogCount = log2(Count + 1))

p11 <- ggplot(norm_long, aes(x = Sample, y = LogCount, fill = Sample)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_viridis_d(option = "mako") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Distribution of Normalized Counts (Violin Plot)",
       subtitle = "Verifying comparability across samples after Median-of-Ratios normalization",
       y = "Log2(Normalized Count + 1)")

print(p11)

# --- Plot 12: Stratified P-value Histogram ---
# Checking if statistical power depends on expression level (common RNA-seq issue)
res_df_strat <- res_df %>%
  mutate(ExprBin = cut(log10(baseMean), breaks = 4, labels = c("Low", "Medium", "High", "Very High")))

p12 <- ggplot(res_df_strat, aes(x = pvalue, fill = ExprBin)) +
  geom_histogram(bins = 30, color = "white", alpha = 0.8) +
  facet_wrap(~ExprBin, scales = "free_y") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "P-value Distribution Stratified by Expression Level",
       subtitle = "Checking for uniformity under null across different abundance regimes",
       x = "P-value", y = "Frequency")

print(p12)

# --- Plot 13: Conformal Prediction Interval (Top Gene) ---
# Visualizing the 'Distribution-Free' interval for the top hit
# 1. Calculate 95% Quantile of all Non-Conformity Scores (Global calibration for demo)
# In full split-conformal, this is done on calibration set. Here we use full set for visualization.
q_95 <- quantile(as.vector(V_ig), 0.95)

# 2. Get data for top gene
top_gene_idx <- which.min(res$padj)
mu_top <- mu_hat[top_gene_idx, ]
sigma_top <- sigma_hat[top_gene_idx, ]
obs_top <- Y_ig[top_gene_idx, ]
# Interval = mu +/- q * sigma
lower_bound <- pmax(0, mu_top - q_95 * sigma_top)
upper_bound <- mu_top + q_95 * sigma_top

conf_df <- data.frame(
  Sample = names(obs_top),
  Observed = obs_top,
  Predicted = mu_top,
  Lower = lower_bound,
  Upper = upper_bound,
  Cell = colData(dds)$cell,
  Dex = colData(dds)$dex
)

p13 <- ggplot(conf_df, aes(x = Sample)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "gray50", size = 1) +
  geom_point(aes(y = Predicted), color = "blue", shape = 3, size = 3) + # Predicted Mean
  geom_point(aes(y = Observed, color = Cell, shape = Dex), size = 4) +
  scale_color_viridis_d(option = "turbo") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste0("Conformal Prediction Intervals: ", rownames(res)[top_gene_idx]),
       subtitle = paste0("Intervals constructed using global score quantile Q(0.95) = ", round(q_95, 2)),
       y = "Normalized Count", caption = "Blue Cross: GLM Prediction | Bar: Conformal Interval | Dot: Observed")

print(p13)

# ------------------------------------------------------------------------------
# 11. Additional Tables
# ------------------------------------------------------------------------------

# Table of Sample Totals
sample_stats <- as.data.frame(colData(dds))
sample_stats$TotalCounts <- colSums(counts(dds))
sample_stats$NormFactors <- sizeFactors(dds)

print("--- Sample Statistics Table ---")
print(sample_stats[, c("cell", "dex", "TotalCounts", "NormFactors")])

# ------------------------------------------------------------------------------
# 12. Final Methodological Diagnostics & Biological Summary
# ------------------------------------------------------------------------------

# --- Plot 14: Independent Filtering Diagnostic ---
# Methodology: Independent filtering is used to maximize rejections while maintaining FDR
# We visualize how the filtering threshold affects the number of significant genes
p14 <- plotMetadata <- as.data.frame(metadata(res)$filterNumRej)
p14 <- ggplot(plotMetadata, aes(x=theta, y=numRej)) +
  geom_line(color="darkred", size=1) +
  geom_point(color="black") +
  theme_minimal() +
  labs(title = "Independent Filtering Diagnostic",
       subtitle = "Justifying the choice of expression threshold to maximize power",
       x = "Filtering Quantile (theta)", y = "Number of Rejections")

print(p14)

# --- Plot 15: Cook's Distance (Model Integrity) ---
# Methodology: Identifying influential outliers that could bias the 'Working Model'
cooks <- assays(dds)[["cooks"]]
cooks_long <- as.data.frame(cooks) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Distance")

p15 <- ggplot(cooks_long, aes(x = Sample, y = Distance, fill = Sample)) +
  geom_boxplot(outlier.color = "red", outlier.size = 0.5) +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cook's Distance per Sample",
       subtitle = "Diagnostic for the Working Model GLM (Identifying outliers)",
       y = "Log10 Cook's Distance")

print(p15)

# --- Table: Biological Pathway / Function Proxy Table ---
# Methodology: Categorizing results to provide context for 'Distribution-Free' hits
bio_proxy <- rich_table %>%
  mutate(
    Candidate_Role = case_when(
      log2FoldChange > 3 ~ "Strong Induction",
      log2FoldChange < -3 ~ "Strong Repression",
      abs(log2FoldChange) > 1 & baseMean > 1000 ~ "High-Abundance Shift",
      TRUE ~ "Moderate Modulator"
    ),
    Statistical_Confidence = ifelse(padj < 1e-10, "Exceptional", "High")
  ) %>%
  select(Candidate_Role, log2FoldChange, Statistical_Confidence, Significance)

print("--- Biological Inference Summary Table ---")
print(bio_proxy)