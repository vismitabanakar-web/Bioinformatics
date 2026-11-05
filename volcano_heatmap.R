############################################################
# Project: FN1-mediated anoikis resistance in OSCC
# Dataset: GSE160395 (GEO)
# Analysis: Differential expression, Volcano Plot, Heatmap
############################################################

############################
# 1. Install and load packages
############################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c("ggplot2", "ggrepel", "pheatmap")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

bioc_packages <- c("GEOquery", "limma")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)

############################
# 2. Download GSE160395
############################

gse <- getGEO("GSE160395", GSEMatrix = TRUE)
gset <- gse[[1]]

expr <- exprs(gset)
pheno <- pData(gset)
fdata <- fData(gset)

# Create an output folder for figures and result tables
dir.create("results", showWarnings = FALSE)

############################
# 3. Inspect sample information and define groups
############################

cat("\nSample titles:\n")
print(pheno$title)

# This follows the grouping rule in the supplied analysis:
# titles containing "M3" are classified as Metastatic;
# all other samples are classified as Primary.
# IMPORTANT: Check the printed titles and confirm this rule is
# appropriate for every sample before interpreting the results.

group <- ifelse(grepl("M3", pheno$title, ignore.case = TRUE),
                "Metastatic", "Primary")
group <- factor(group, levels = c("Primary", "Metastatic"))

if (length(group) != ncol(expr)) {
  stop("The number of group labels does not match the expression samples.")
}

if (any(table(group) == 0)) {
  stop("One group has no samples. Check the sample titles and grouping rule.")
}

sample_info <- data.frame(
  Sample = colnames(expr),
  Title = pheno$title,
  Group = group
)
write.csv(sample_info, "results/sample_group_information.csv",
          row.names = FALSE)

cat("\nSamples in each group:\n")
print(table(group))

############################
# 4. Normalize expression data
############################

# Appropriate for the microarray expression matrix used in this workflow.
expr_norm <- normalizeBetweenArrays(expr)

############################
# 5. Differential expression analysis using limma
############################

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_norm, design)

contrast_matrix <- makeContrasts(
  Metastatic_vs_Primary = Metastatic - Primary,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

deg <- topTable(
  fit2,
  coef = "Metastatic_vs_Primary",
  adjust.method = "BH",
  number = Inf,
  sort.by = "P"
)

############################
# 6. Add gene symbols
############################

# Use the annotation columns available in the GEO platform data.
# The supplied analysis uses GENE_SYMBOL and ID.
if ("GENE_SYMBOL" %in% colnames(fdata)) {
  gene_symbols <- as.character(fdata$GENE_SYMBOL)
} else if ("Gene Symbol" %in% colnames(fdata)) {
  gene_symbols <- as.character(fdata$`Gene Symbol`)
} else {
  stop("Gene-symbol annotation column not found in fData(gset). ",
       "Inspect colnames(fData(gset)) and update this section.")
}

# Match annotation rows to expression rows.
if (!is.null(rownames(fdata)) &&
    all(rownames(expr_norm) %in% rownames(fdata))) {
  names(gene_symbols) <- rownames(fdata)
} else if ("ID" %in% colnames(fdata)) {
  names(gene_symbols) <- as.character(fdata$ID)
} else {
  stop("Could not match platform annotations to expression probes.")
}

deg$ProbeID <- rownames(deg)
deg$GeneSymbol <- gene_symbols[deg$ProbeID]

# Remove missing/blank symbols for gene-labelled output.
deg_annotated <- deg[
  !is.na(deg$GeneSymbol) & trimws(deg$GeneSymbol) != "",
]

write.csv(deg, "results/differential_expression_all_probes.csv",
          row.names = FALSE)
write.csv(deg_annotated, "results/differential_expression_annotated.csv",
          row.names = FALSE)

############################
# 7. Classify genes for the Volcano Plot
############################

logFC_cutoff <- 1
adjusted_p_cutoff <- 0.05

deg$status <- "Not significant"
deg$status[deg$logFC >= logFC_cutoff &
             deg$adj.P.Val < adjusted_p_cutoff] <- "Upregulated"
deg$status[deg$logFC <= -logFC_cutoff &
             deg$adj.P.Val < adjusted_p_cutoff] <- "Downregulated"

# Genes of interest from the supplied OSCC project analysis.
genes_interest <- c(
  "ITGB4", "ITGA6", "ITGA2", "ITGAV", "CD44",
  "PTK2", "ILK", "YAP1", "FN1", "SPP1",
  "TIMP1", "BCL2L1", "ATF4", "CALR", "CTNNB1",
  "PTGS2", "NQO1", "HMCN1", "CSPG4", "ANXA5"
)

highlight <- deg[
  !is.na(deg$GeneSymbol) & deg$GeneSymbol %in% genes_interest,
]

############################
# 8. Create and save Volcano Plot
############################

# Avoid infinite values if adjusted p-values are zero.
deg$plot_p <- pmax(deg$adj.P.Val, .Machine$double.xmin)
deg$minus_log10_adj_p <- -log10(deg$plot_p)

volcano_plot <- ggplot(
  deg,
  aes(x = logFC, y = minus_log10_adj_p, color = status)
) +
  geom_point(size = 1.8, alpha = 0.7) +
  scale_color_manual(values = c(
    "Downregulated" = "blue",
    "Not significant" = "grey70",
    "Upregulated" = "red"
  )) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(adjusted_p_cutoff),
             linetype = "dashed", color = "black") +
  geom_text_repel(
    data = highlight,
    aes(label = GeneSymbol),
    size = 3.5,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  labs(
    title = "Volcano Plot – GSE160395",
    subtitle = "Metastatic versus Primary samples",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Expression status"
  ) +
  theme_classic(base_size = 12)

print(volcano_plot)

ggsave(
  "results/volcano_plot_GSE160395.png",
  plot = volcano_plot,
  width = 9,
  height = 7,
  dpi = 300
)

############################
# 9. Prepare expression matrix for Heatmap
############################

# Select the project genes that are present in the platform annotation.
heatmap_idx <- which(
  !is.na(gene_symbols) & gene_symbols %in% genes_interest
)

if (length(heatmap_idx) == 0) {
  stop("None of the selected genes were found in the platform annotation.")
}

heat_expr <- expr_norm[heatmap_idx, , drop = FALSE]
heat_symbols <- gene_symbols[heatmap_idx]

# If multiple probes map to the same gene, average their expression.
heat_expr_df <- as.data.frame(heat_expr)
heat_expr_df$GeneSymbol <- heat_symbols

heat_expr_gene <- aggregate(
  . ~ GeneSymbol,
  data = heat_expr_df,
  FUN = mean
)

rownames(heat_expr_gene) <- heat_expr_gene$GeneSymbol
heat_expr_gene$GeneSymbol <- NULL
heat_expr_gene <- as.matrix(heat_expr_gene)

# Z-score each gene across samples for visualization.
heat_expr_scaled <- t(scale(t(heat_expr_gene)))

# Remove rows with non-finite values (e.g., genes with zero variance).
heat_expr_scaled <- heat_expr_scaled[
  apply(heat_expr_scaled, 1, function(x) all(is.finite(x))),
  ,
  drop = FALSE
]

if (nrow(heat_expr_scaled) == 0) {
  stop("No variable genes remain for the heatmap after scaling.")
}

############################
# 10. Add sample group annotation and create Heatmap
############################

annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(heat_expr_scaled)

heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

png(
  "results/heatmap_OSCC_genes_GSE160395.png",
  width = 1800,
  height = 1400,
  res = 200
)

pheatmap(
  heat_expr_scaled,
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 9,
  fontsize_col = 8,
  main = "Heatmap of Selected OSCC-Related Genes\nGSE160395"
)

dev.off()

############################
# 11. Save selected gene results
############################

selected_gene_results <- deg_annotated[
  deg_annotated$GeneSymbol %in% genes_interest,
]

write.csv(
  selected_gene_results,
  "results/selected_OSCC_gene_results.csv",
  row.names = FALSE
)

cat("\nAnalysis completed.\n")
cat("Check the 'results' folder for plots and CSV files.\n")
