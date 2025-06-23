
library(tidyverse)
library(Seurat)
library(dplyr)
library(patchwork)
library(org.Mm.eg.db)
library(clusterProfiler)
library(msigdbr)
library(dittoSeq)
library(EnhancedVolcano)
library(GOSemSim)
library(enrichplot)
library(ggupset)
library(DESeq2)
library(ggpubr)
library(ggrepel)
library(AnnotationDbi)
library(enrichR)
library(ggforce)
library(ggVennDiagram)
library(stringr)
library(glmGamPoi)

# Set seed for reproducibility
set.seed(1234)

cat("=== Starting Comprehensive MSN Analysis ===\n")

# Updated sample comparisons (swapped x to y and y to x)
sample_comparisons <- list(
  list(x = "WT_Cas9sgCTG",  y = "WT_Cas9",         title = "WT_Cas9sgCTG vs WT_Cas9"),
  list(x = "R61_Cas9sgCTG", y = "R61_Cas9",        title = "R61_Cas9sgCTG vs R61_Cas9"),
  list(x = "WT_Cas9sgCTG",  y = "R61_Cas9",        title = "WT_Cas9sgCTG vs R61_Cas9"),
  list(x = "R61_Cas9sgCTG", y = "WT_Cas9",         title = "R61_Cas9sgCTG vs WT_Cas9"),
  list(x = "R61_Cas9",      y = "WT_Cas9",         title = "R61_Cas9 vs WT_Cas9"),
  list(x = "R61_Cas9sgCTG", y = "WT_Cas9sgCTG",    title = "R61_Cas9sgCTG vs WT_Cas9sgCTG")
)

# get % mitochondria
merged_obj = PercentageFeatureSet(merged_obj, "^mt-", col.name = "percent_mito")
library(ggplot2)  # Ensure ggplot2 is loaded for VlnPlot() and NoLegend()
# Open the Cairo PostScript device with specified dimensions (in inches)
cairo_ps(file = "QC_feature_violin.eps", width = 7, height = 7, onefile = FALSE, bg = "transparent")

# Explicitly print the ggplot object to the EPS device
print(VlnPlot(merged_obj, features = c("nCount_RNA", "nFeature_RNA", "percent_mito"), pt.size = 0) + NoLegend())

# Close the device to save the EPS file
dev.off()

# Scatter QC plots.
# Open the EPS device using cairo_ps() for potentially improved rendering
cairo_ps(file = "QC_feature_scatter.eps", width = 5, height = 5,
         onefile = FALSE, bg = "transparent")

# Create your scatter plot
FeatureScatter(merged_obj, feature1 = "nCount_RNA", feature2 = "percent_mito")

# Close the EPS device to save the file
dev.off()

# filter RNA
merged_obj = subset(merged_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent_mito < 5)
merged_obj

sex_map <- c(
  "WT_Cas9_1" = "F", "WT_Cas9_2" = "F", "WT_Cas9_3" = "F",
  "WT_Cas9sgCTG_1" = "F", "WT_Cas9sgCTG_2" = "M", "WT_Cas9sgCTG_3" = "F",
  "R61_Cas9_1" = "M", "R61_Cas9_2" = "M", "R61_Cas9_3" = "F",
  "R61_Cas9sgCTG_1" = "F", "R61_Cas9sgCTG_2" = "M", "R61_Cas9sgCTG_3" = "F"
)
merged_obj@meta.data$sex <- sex_map[merged_obj@meta.data$orig.ident]
merged_obj$sex_numeric <- ifelse(merged_obj$sex == "M", 1, 0)

# SCT to remove sex based batch effect
merged_obj = SCTransform(merged_obj, 
                        vars.to.regress = c("percent_mito", "nCount_RNA", "sex_numeric"), 
                        verbose = TRUE)

# Run the PCA
merged_obj = RunPCA(object = merged_obj, assay ="SCT", reduction.name="pca_unintegrated")

# QC plots and parameter setting
dims_to_use = 20

## Reduce and cluster - using dims_to_use
# get the clusters - we do not optimisise, as we will integrate.
merged_obj = FindNeighbors(object = merged_obj, reduction = "pca_unintegrated", dims = 1:dims_to_use)
merged_obj = FindClusters(object = merged_obj, resolution = 0.5)
merged_obj = RunUMAP(merged_obj, reduction = "pca_unintegrated", reduction.name ="umap_unintegrated", dims = 1:dims_to_use, min.dist = 0.5, spread = 1)



# UMAP_unintegrated plot
postscript(file = "UMAP_unintegrated.eps", width = 10, height = 10,
           horizontal = FALSE, onefile = FALSE, bg = "transparent", paper = "special")
{
  p1 <- DimPlot(merged_obj, reduction = "umap_unintegrated", group.by = "orig.ident", label = TRUE, pt.size = 1)
  p2 <- DimPlot(merged_obj, reduction = "umap_unintegrated", group.by = "seurat_clusters", label = TRUE, label.size = 3)
  # For patchwork objects (p1 + p2), wrap in print()
  print(p1 + p2)
}
dev.off()

# QC_Heatmap plot
postscript(file = "QC_Heatmap.eps", width = 1500/300, height = 5000/300,
           horizontal = FALSE, onefile = FALSE, bg = "transparent", paper = "special")
print(DimHeatmap(merged_obj, reduction = "pca_unintegrated", dims = 1:50, cells = 500, balanced = TRUE))
dev.off()

# Harmony with orig.ident, sex, AND nCount_RNA
library(harmony)
merged_obj <- RunHarmony(
  merged_obj,
  assay.use = "SCT",
  group.by.vars = "orig.ident",
  theta = 1,
  reduction = "pca_unintegrated",
  reduction.save = "pca_integrated"
)


# re-cluster after integration
merged_obj = FindNeighbors(merged_obj, reduction = "pca_integrated", dims = 1:dims_to_use)
merged_obj = FindClusters(merged_obj, resolution = 0.5)
merged_obj = RunUMAP(merged_obj, dims = 1:dims_to_use, reduction = "pca_integrated", reduction.name ="umap_integrated")

# Create individual plots 
library(ggpubr)

p1 <- DimPlot(merged_obj, group.by = "orig.ident", reduction = "umap_integrated",
              label = TRUE, pt.size = 0.5) + NoLegend()
p2 <- DimPlot(merged_obj, group.by = "seurat_clusters", reduction = "umap_integrated",
              label = TRUE, label.size = 3) + NoLegend()
combined_plot <- ggarrange(p1, p2, ncol = 2, nrow = 1)

# Save integrated UMAP
postscript(file = "UMAP_integrated.eps", width = 10, height = 10,
           horizontal = FALSE, onefile = FALSE, bg = "transparent", paper = "special")
print(combined_plot)
dev.off()

## Task 4 - Identify cell types ##
library('presto')

# prep markers - needed to combine the 4 samples
merged_obj = PrepSCTFindMarkers(merged_obj)

# get the marker genes
merged_obj.markers = FindAllMarkers(merged_obj, only.pos = TRUE)

# get the top n
merged_obj.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 15) %>% ungroup() -> topn_markers

# heatmap
# Open the EPS device using postscript()
postscript(file = "Markers_heatmap.eps", width = 10, height = 10,
           horizontal = FALSE, onefile = FALSE, bg = "transparent", paper = "special")

# Explicitly print the DoHeatmap output to the EPS device
print(DoHeatmap(merged_obj, features = topn_markers$gene) + NoLegend())

# Close the device to save the file
dev.off()

celltype_markers <- list(
  MSN=c("Ppp1r1b","Bcl11b","Pde1b","Drd2","Penk","Grik3","Ttc12","Gpr6","Adora2a","Drd1","Tac1","Pdyn","Ebf1","Slc35d3","Foxp2", "Cacng5", "Olfm3","Dcx"),
  Astrocytes = c("Slc1a2", "Gfap", "Slc6a11", "Grm3"),
  Microglia = c("Csf1r", "Trem2", "Gpr34", "Gal3st4"),
  OPC = c("Cspg4", "Gpr17", "Neu4", "Zfp488", "Olig2", "A930009A15Rik", "Olig1", "Rlbp1"),
  Oligodendrocytes = c("Plp1", "Mobp", "Prr5l", "Cdh20", "Mag", "Mog"),
  "Cholinergic Interneuron" = "Chat",
  "PV/Th Interneuron" = c("Hs3st2", "Cntnap4"),
  "Sst/Npy Interneuron" = "Nos1",
  Endothelial = c("Atp13a5", "Gm5127", "Fbln5", "Zic3", "Ltbp4", "Slc2a1"),
  Mural = c("Car3", "Mylk", "Cald1", "Flna", "Lpl", "Pdgfrb"),
  "Ciliated Ependymal" = c("Hydin", "Scgb1a1", "Sec14l3", "Cdc20b", "Cyp2f2"),
  "Secretory Ependymal" = "Ttr",
  "Cycling cell" = "Top2a"
)
merged_object <- AddModuleScore(
  object = merged_obj,
  features = celltype_markers,
  name = "CellTypeScore"
)
library(dplyr)

cell_types <- names(celltype_markers)
score_cols <- paste0("CellTypeScore", seq_along(cell_types))

avg_scores <- merged_object@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean)) %>%
  rowwise() %>%
  mutate(AssignedCellType = cell_types[which.max(c_across(all_of(score_cols)))])

print(avg_scores)

merged_object$sample_group <- gsub("_\\d+$", "", merged_object$orig.ident)

write.csv(avg_scores, file = "average_module_scores_with_assigned_cell_types.csv", row.names = FALSE)
cluster2annotation <- setNames(avg_scores$AssignedCellType, avg_scores$seurat_clusters)

new_cell_type <- cluster2annotation[ as.character(merged_object$seurat_clusters) ]
names(new_cell_type) <- rownames(merged_object@meta.data)

merged_object <- AddMetaData(merged_object, metadata = new_cell_type, col.name = "cell_type")

# Create the UMAP plot 
umap_celltype <- DimPlot(
  merged_object,
  reduction = "umap_integrated",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("UMAP with Automatically Assigned Cell Type Annotations") +
  NoLegend()  # ggpubr function to remove legend

# Save the plot
ggsave(filename = "UMAP_cell_type_annotation_automatic_updated.eps",
       plot = umap_celltype,
       device = cairo_ps,  # Use cairo_ps as the EPS device
       width = 10, height = 10, units = "in", bg = "transparent")


# Print the plot to the default device
print(umap_celltype)

# Bar Plot of Frequencies using dittoBarPlot

library(dittoSeq)
library(ggpubr)

# Create the plot and store it in a variable
p <- dittoBarPlot(object = merged_object, var = "cell_type", group.by = "orig.ident") +
  theme_pubr() +
  theme(
    text = element_text(size = 10),                          # All text elements size 10
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # X-axis text size 10
    axis.text.y = element_text(size = 10),                   # Y-axis text size 10
    axis.title = element_text(size = 10),                    # Axis titles size 10
    legend.text = element_text(size = 10),                   # Legend text size 10
    legend.title = element_text(size = 10),                  # Legend title size 10
    plot.title = element_text(size = 10),                    # Plot title size 10
    legend.position = "right"
  )

# Save the plot as PDF
ggsave(filename = "Cells_Frequency.pdf",
       plot = p,
       device = "pdf",
       width = 10, height = 8, units = "in", bg = "white")


# Create a frequency table of cell types
celltype_freq <- table(merged_object@meta.data$cell_type)

# Convert to a data frame with proportions
celltype_df <- as.data.frame(celltype_freq)
names(celltype_df)[names(celltype_df) == "Var1"] <- "cell_type"
names(celltype_df)[names(celltype_df) == "Freq"] <- "count"

celltype_df <- celltype_df %>%
  mutate(proportion = count / sum(count) * 100)

# Inspect the table
print(celltype_df)

library(ggplot2)

donut <- ggplot(celltype_df, aes(x = 2, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  xlim(1, 2.5) +
  theme_pubr() +
  theme(
    legend.position = "right",
    text = element_text(size = 10),           
    plot.title = element_text(size = 10),     
    legend.text = element_text(size = 10),    
    legend.title = element_text(size = 10)    
  ) +
  geom_text(aes(label = paste0(round(proportion, 2), "%")),
            position = position_stack(vjust = 0.5), 
            size = 3.5) + 
  ggtitle("Cell-Type Composition (Donut Chart)")

ggsave(filename = "Donut_all_cells.pdf",
       plot = donut,
       device = "pdf",
       width = 10, height = 10, units = "in", bg = "white")


#Differential gene expression - MSNs

cell_types <- c("MSN")

all_results <- list()


for (comp in sample_comparisons) {
  for (ct in cell_types) {
    message("Processing: ", comp$title, " | Cell type: ", ct)
    
    # Subset data
    sub_obj <- subset(merged_object, subset = cell_type == ct & sample_group %in% c(comp$x, comp$y))
    
    # Transform if needed
    if (!"SCT" %in% names(sub_obj@assays)) {
      sub_obj <- SCTransform(sub_obj, vars.to.regress = c("percent_mito", "sex_numeric","nCount_RNA"), verbose = FALSE)
    }
    
    # Set identity and run DE
    sub_obj <- SetIdent(sub_obj, value = "sample_group")
    
    tryCatch({
      de_res <- FindMarkers(
        sub_obj,
        ident.1 = comp$x,
        ident.2 = comp$y,
        test.use = "MAST",
        verbose = FALSE, 
        min.pct = 0.1,
        logfc.threshold = 0
      )
      
      if (nrow(de_res) == 0) {
        message("  No DE genes found, skipping.")
        next
      }
      
      # Add metadata
      de_res$gene <- rownames(de_res)
      de_res$comparison <- comp$title
      de_res$cell_type <- ct
      
      # Store results
      result_name <- paste(comp$title, ct, sep = "_")
      all_results[[result_name]] <- de_res
      
    }, error = function(e) {
      message("  Error in DE analysis: ", e$message)
    })
  }
}

#Pseudobulk

# Set seed for reproducibility
set.seed(1234)
DefaultAssay(merged_object) <- "RNA"

# Output directory
output_dir <- "/scratch/c.mpmrrp/sc_analysis/Jan_2025_AMB_single_cell_nuclei/output/new_align_March/28_May_2025/DE_results_pseudobulk_final_no_covariates"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Comparisons list
sample_comparisons <- list(
  list(x = "WT_Cas9sgCTG",  y = "WT_Cas9",         title = "WT_Cas9sgCTG vs WT_Cas9"),
  list(x = "R61_Cas9sgCTG", y = "R61_Cas9",        title = "R61_Cas9sgCTG vs R61_Cas9"),
  list(x = "R61_Cas9",      y = "WT_Cas9",         title = "R61_Cas9 vs WT_Cas9"),
  list(x = "R61_Cas9sgCTG", y = "WT_Cas9sgCTG",    title = "R61_Cas9sgCTG vs WT_Cas9sgCTG")
)

# Create pseudobulk data
counts <- GetAssayData(merged_object, slot = "counts")
meta <- merged_object@meta.data

# Get all unique replicates
all_replicates <- unique(meta$orig.ident)

# Create pseudobulk count matrix
pseudobulk_counts <- matrix(0, nrow = nrow(counts), ncol = length(all_replicates))
rownames(pseudobulk_counts) <- rownames(counts)
colnames(pseudobulk_counts) <- all_replicates

for (replicate in all_replicates) {
  replicate_cells <- colnames(merged_object)[meta$orig.ident == replicate]
  if (length(replicate_cells) > 0) {
    pseudobulk_counts[, replicate] <- Matrix::rowSums(counts[, replicate_cells, drop = FALSE])
  }
}

# Create sample metadata
sample_meta <- meta %>%
  group_by(orig.ident) %>%
  summarise(
    sample_group = sample_group[1],
    n_cells = n(),
    .groups = 'drop'
  ) %>%
  as.data.frame()

rownames(sample_meta) <- sample_meta$orig.ident

# Filter genes with very low expression
min_samples <- max(3, ceiling(ncol(pseudobulk_counts) * 0.5))
keep_genes <- rowSums(pseudobulk_counts >= 10) >= min_samples
pseudobulk_counts_filtered <- pseudobulk_counts[keep_genes, ]

# Function to run DESeq2 analysis
run_deseq2_analysis <- function(comp) {
  # Filter to samples in this comparison
  comp_samples <- sample_meta$orig.ident[sample_meta$sample_group %in% c(comp$x, comp$y)]
  
  if (length(comp_samples) < 4) {
    return(NULL)
  }
  
  # Subset data
  comp_counts <- pseudobulk_counts_filtered[, comp_samples, drop = FALSE]
  comp_metadata <- sample_meta[comp_samples, ]
  
  # Check group sizes
  group_counts <- table(comp_metadata$sample_group)
  
  if (any(group_counts < 2)) {
    return(NULL)
  }
  
  tryCatch({
    comp_metadata$sample_group <- factor(comp_metadata$sample_group, levels = c(comp$y, comp$x))
    
    dds <- DESeqDataSetFromMatrix(
      countData = comp_counts,
      colData = comp_metadata,
      design = ~ sample_group
    )
    
    # Additional filtering within comparison
    keep <- rowSums(counts(dds) >= 10) >= max(2, min(group_counts))
    dds_filtered <- dds[keep, ]
    
    # Run DESeq2
    dds_filtered <- DESeq(dds_filtered, quiet = TRUE)
    # Contrast: comp$x (test) vs comp$y (reference)
    res <- results(dds_filtered, contrast = c("sample_group", comp$x, comp$y))
    
    # Convert to data frame and filter NA values
    de_res <- as.data.frame(res)
    de_res <- de_res[!is.na(de_res$pvalue), ]
    de_res$gene <- rownames(de_res)
    
    # Rename columns
    colnames(de_res)[colnames(de_res) == "log2FoldChange"] <- "avg_log2FC"
    colnames(de_res)[colnames(de_res) == "padj"] <- "p_val_adj"
    colnames(de_res)[colnames(de_res) == "pvalue"] <- "p_val"
    
    # Add normalized counts
    if (nrow(de_res) > 0) {
      normalized_counts <- counts(dds_filtered, normalized = TRUE)
      group1_samples <- comp_metadata$orig.ident[comp_metadata$sample_group == comp$x]
      group2_samples <- comp_metadata$orig.ident[comp_metadata$sample_group == comp$y]
      
      genes_present <- intersect(rownames(de_res), rownames(normalized_counts))
      
      de_res$mean_counts_group1 <- NA
      de_res$mean_counts_group2 <- NA
      
      if (length(genes_present) > 0) {
        de_res[genes_present, "mean_counts_group1"] <- rowMeans(normalized_counts[genes_present, group1_samples, drop = FALSE])
        de_res[genes_present, "mean_counts_group2"] <- rowMeans(normalized_counts[genes_present, group2_samples, drop = FALSE])
      }
    }
    
    return(list(
      results = de_res,
      dds = dds_filtered,
      comparison = comp
    ))
    
  }, error = function(e) {
    return(NULL)
  })
}

# Run analysis for each comparison
all_results <- list()

for (comp in sample_comparisons) {
  comparison_name <- gsub(" ", "_", comp$title)
  
  # File names
  all_genes_fname <- file.path(output_dir, paste0(comparison_name, "_pseudobulk_ALL_genes.csv"))
  sig_genes_fname <- file.path(output_dir, paste0(comparison_name, "_pseudobulk_SIGNIFICANT_genes.csv"))
  
  # Run analysis
  analysis_result <- run_deseq2_analysis(comp)
  if (is.null(analysis_result)) next
  
  de_res <- analysis_result$results
  all_results[[comparison_name]] <- analysis_result
  
  # Define significance
  de_res <- de_res %>%
    mutate(
      Significance = case_when(
        p_val_adj < 0.01 & avg_log2FC >= 0.585 ~ "Significantly Up",
        p_val_adj < 0.01 & avg_log2FC <= -0.585 ~ "Significantly Down",
        TRUE ~ "Not Significant"
      )
    )
  
  # Save results
  write.csv(de_res, all_genes_fname, row.names = FALSE)
  
  # Save significant genes
  sig_genes_df <- de_res %>%
    filter(Significance %in% c("Significantly Up", "Significantly Down"))
  
  if (nrow(sig_genes_df) > 0) {
    write.csv(sig_genes_df, sig_genes_fname, row.names = FALSE)
  }
}

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

#Volcano plot for both Pseudobulk & Findmarkers

# 1. DIRECTORY SETTINGS
INPUT_DIR <- "/scratch/c.mpmrrp/sc_analysis/Jan_2025_AMB_single_cell_nuclei/output/new_align_March/28_May_2025/DE_results/less_threshold/"  # Directory containing CSV files
OUTPUT_DIR <- file.path(INPUT_DIR, "volcano_plots")  
FILE_PATTERN <- "\\.csv$"                           

# 2. COLUMN NAME MAPPING (adjust to match your data)
GENE_COL <- "gene"           # Column with gene names
LFC_COL <- "avg_log2FC"      # Column with log2 fold change
PVAL_COL <- "p_val_adj"      # Column with adjusted p-values
SIGNIFICANCE_COL <- "Significance"  # Column with significance categories (optional)

# 3. SIGNIFICANCE THRESHOLDS
P_THRESHOLD <- 0.01          # P-value threshold for significance
LFC_THRESHOLD <- 0.585       # Log2FC threshold for significance

# 4. PLOT APPEARANCE
DOT_SIZE <- 2.4              # Size of points (increase for larger dots)
DOT_ALPHA <- 0.7             # Transparency of points (0-1)
LABEL_TOP_N <- 0            # Number of top genes to label

# 5. AXIS LIMITS AND SPACING
X_LIMITS <- c(-4, 4.2)         # X-axis limits [min, max]
Y_LIMIT_MAX <- 350         # Maximum Y-axis value
X_BREAKS_BY <- 1             # X-axis break intervals
Y_BREAKS_BY <- 50        # Y-axis break intervals (in 100s)

# 6. COLORS
COLORS <- list(
  up = "#D73027",            # Color for up-regulated genes
  down = "#1A9850",          # Color for down-regulated genes
  not_sig = "#CCCCCC"        # Color for non-significant genes
)

# 7. OUTPUT OPTIONS
CREATE_COMPACT_VERSION <- TRUE   # Create an additional compact version?
SAVE_PDF <- TRUE                # Save PDF version?
SAVE_PNG <- TRUE                # Save PNG version?
DPI <- 300                      # Resolution for PNG
PLOT_WIDTH <- 10                # Plot width in inches
PLOT_HEIGHT <- 7                # Plot height in inches

# 8. FILE FILTERING OPTIONS
EXCLUDE_PATTERNS <- c("correlation", "summary", "enrichment")  # Skip files containing these words
MIN_FILE_SIZE_KB <- 1           # Skip files smaller than this (in KB)

# Function to create volcano plot
create_volcano_plot <- function(input_file, output_dir, output_prefix,
                               gene_col = "gene", lfc_col = "avg_log2FC", 
                               pval_col = "p_val_adj", significance_col = NULL,
                               p_threshold = 0.01, lfc_threshold = 0.585,
                               x_limits = c(-4, 4), y_limit_max = 20,
                               x_breaks_by = 1, y_breaks_by = 100,
                               dot_size = 2.4, dot_alpha = 0.7, label_top_n = 15,
                               colors = list(up = "#E31A1C", down = "#1F78B4", not_sig = "#CCCCCC"),
                               create_compact = TRUE, save_pdf = TRUE, save_png = TRUE,
                               dpi = 300, plot_width = 10, plot_height = 7) {
  
  # Check if this is an ALL_genes file for special formatting
  is_all_genes <- grepl("ALL_genes", basename(input_file), ignore.case = TRUE)
  
  # Read data with error handling
  tryCatch({
    data <- read.csv(input_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    message("   ❌ Error reading file: ", e$message)
    return(NULL)
  })
  
  if (is.null(data) || nrow(data) == 0) {
    message("   ⚠️ Empty or invalid file, skipping")
    return(NULL)
  }
  
  # Validate required columns
  required_cols <- c(gene_col, lfc_col, pval_col)
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    message("   ⚠️ Missing columns: ", paste(missing_cols, collapse = ", "), " - skipping")
    return(NULL)
  }
  
  # Handle p_val_adj = 0 by setting to smallest representable number
  zero_genes_count <- sum(data[[pval_col]] == 0, na.rm = TRUE)
  data[[pval_col]][data[[pval_col]] == 0] <- .Machine$double.xmin
  
  message("   🔧 Imputed ", zero_genes_count, " genes with p_val_adj = 0")
  
  # Create significance categories if not provided
  if (is.null(significance_col) || !significance_col %in% colnames(data)) {
    data$Significance <- case_when(
      data[[pval_col]] < p_threshold & data[[lfc_col]] >= lfc_threshold ~ "Significantly Up",
      data[[pval_col]] < p_threshold & data[[lfc_col]] <= -lfc_threshold ~ "Significantly Down",
      TRUE ~ "Not Significant"
    )
    significance_col <- "Significance"
  }
  
  # Count significant genes
  sig_counts <- table(data[[significance_col]])
  sig_up <- ifelse("Significantly Up" %in% names(sig_counts), sig_counts["Significantly Up"], 0)
  sig_down <- ifelse("Significantly Down" %in% names(sig_counts), sig_counts["Significantly Down"], 0)
  total_sig <- sig_up + sig_down
  
  message("   📊 ", nrow(data), " genes | ", total_sig, " significant (", sig_up, " up, ", sig_down, " down)")
  
  if (total_sig == 0) {
    message("   ⚠️ No significant genes found, creating plot anyway")
  }
  
  # Prepare plot data
  plot_data <- data %>%
    filter(!is.na(.data[[lfc_col]]) & !is.na(.data[[pval_col]])) %>%
    mutate(
      neg_log10_p = ifelse(.data[[pval_col]] == 0, 
                          y_limit_max, 
                          pmin(-log10(.data[[pval_col]]), y_limit_max))
    )
  
  # Auto-detect title from filename
  plot_title <- gsub("_", " ", tools::file_path_sans_ext(basename(input_file)))
  plot_title <- gsub("ALL genes|all genes", "", plot_title)
  plot_title <- trimws(plot_title)
  
  # Set up colors
  color_mapping <- c(
    "Significantly Up" = colors$up,
    "Significantly Down" = colors$down, 
    "Not Significant" = colors$not_sig
  )
  
  # Create y-axis breaks - FIXED: NO EXTENSION
  y_breaks <- seq(0, y_limit_max, by = y_breaks_by)
  
  # Create base plot
  p_volcano <- ggplot(plot_data, aes(x = .data[[lfc_col]], y = neg_log10_p)) +
    geom_point(aes(color = .data[[significance_col]]), 
               size = dot_size, alpha = dot_alpha) +
    scale_color_manual(values = color_mapping, name = "Gene Regulation") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      plot.margin = margin(10, 10, 5, 5),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.title.y = element_text(margin = margin(r = 5)),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      panel.border = element_blank()
    ) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), 
               linetype = "dashed", color = "grey50", alpha = 0.7) +
    geom_hline(yintercept = -log10(p_threshold), 
               linetype = "dashed", color = "grey50", alpha = 0.7) +
    scale_x_continuous(
      limits = x_limits,
      breaks = seq(x_limits[1], x_limits[2], by = x_breaks_by),
      expand = expansion(mult = 0, add = 0)  # No expansion - axes meet exactly
    ) +
    scale_y_continuous(
      limits = c(0, y_limit_max),  # FIXED Y-LIMIT - no extension
      breaks = y_breaks,
      expand = expansion(mult = 0, add = 0)  # No expansion - axes meet exactly
    ) +
    labs(
      x = "Log2 Fold Change",
      y = "-log10 Adjusted P-value", 
      title = plot_title,
      subtitle = paste0("p_adj < ", p_threshold, ", |LFC| ≥ ", lfc_threshold)
    ) +
    # Add gene count annotations
    annotate("text", x = x_limits[1] + (x_limits[2] - x_limits[1]) * 0.1, 
             y = y_limit_max * 0.95,
             label = paste0("Down: ", sig_down), hjust = 0, 
             color = colors$down, size = 6, fontface = "bold") +
    annotate("text", x = x_limits[2] - (x_limits[2] - x_limits[1]) * 0.3, 
             y = y_limit_max * 0.95,
             label = paste0("Up: ", sig_up), hjust = 0,
             color = colors$up, size = 6, fontface = "bold")
  
  # Add top gene labels
  if (total_sig > 0 && label_top_n > 0) {
    top_genes <- plot_data %>%
      filter(.data[[significance_col]] %in% c("Significantly Up", "Significantly Down")) %>%
      arrange(.data[[pval_col]]) %>%
      head(label_top_n)
    
    if (nrow(top_genes) > 0) {
      p_volcano <- p_volcano +
        geom_text_repel(
          data = top_genes,
          aes(label = .data[[gene_col]]),
          size = 2.8,
          max.overlaps = label_top_n,
          force = 2,
          force_pull = 0.1,
          box.padding = 0.3,
          point.padding = 0.2,
          segment.color = "grey50",
          segment.size = 0.3,
          min.segment.length = 0.1
        )
    }
  }
  
  # Save plots
  files_created <- c()
  
  if (save_png) {
    png_file <- file.path(output_dir, paste0(output_prefix, "_volcano.png"))
    ggsave(png_file, plot = p_volcano, width = plot_width, height = plot_height, 
           dpi = dpi, bg = "white")
    files_created <- c(files_created, basename(png_file))
  }
  
  if (save_pdf) {
    pdf_file <- file.path(output_dir, paste0(output_prefix, "_volcano.pdf"))
    ggsave(pdf_file, plot = p_volcano, width = plot_width, height = plot_height, 
           bg = "white")
    files_created <- c(files_created, basename(pdf_file))
  }
  
  # Create compact version
  if (create_compact) {
    y_breaks_compact <- seq(0, y_limit_max, by = y_breaks_by * 2)
    
    p_compact <- p_volcano +
      theme(
        plot.margin = margin(5, 5, 2, 2),
        axis.title.x = element_text(margin = margin(t = 2)),
        axis.title.y = element_text(margin = margin(r = 2)),
        legend.margin = margin(0, 0, 0, 5)
      ) +
      scale_x_continuous(
        limits = x_limits,
        breaks = seq(x_limits[1], x_limits[2], by = x_breaks_by * 2),
        expand = expansion(mult = 0, add = 0)
      ) +
      scale_y_continuous(
        limits = c(0, y_limit_max),
        breaks = y_breaks_compact,
        expand = expansion(mult = 0, add = 0)
      )
    
    if (save_png) {
      compact_png <- file.path(output_dir, paste0(output_prefix, "_volcano_compact.png"))
      ggsave(compact_png, plot = p_compact, width = plot_width, height = plot_height, 
             dpi = dpi, bg = "white")
      files_created <- c(files_created, basename(compact_png))
    }
  }
  
  return(list(
    sig_up = sig_up,
    sig_down = sig_down,
    total_sig = total_sig,
    zero_genes_imputed = zero_genes_count,
    files_created = files_created
  ))
}

# Validate input directory
if (!dir.exists(INPUT_DIR)) {
  stop("❌ Input directory not found: ", INPUT_DIR)
}

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  message("📁 Created output directory: ", OUTPUT_DIR)
}

# Find all CSV files
message("🔍 Scanning for CSV files in: ", INPUT_DIR)
all_files <- list.files(INPUT_DIR, pattern = FILE_PATTERN, full.names = TRUE, recursive = FALSE)

# Filter files
valid_files <- c()
for (file in all_files) {
  # Check file size
  file_size_kb <- file.info(file)$size / 1024
  if (file_size_kb < MIN_FILE_SIZE_KB) {
    message("   ⏭️ Skipping small file: ", basename(file), " (", round(file_size_kb, 1), " KB)")
    next
  }
  
  # Check exclude patterns
  filename <- basename(file)
  skip_file <- FALSE
  for (pattern in EXCLUDE_PATTERNS) {
    if (grepl(pattern, filename, ignore.case = TRUE)) {
      message("   ⏭️ Skipping excluded file: ", filename, " (contains '", pattern, "')")
      skip_file <- TRUE
      break
    }
  }
  
  if (!skip_file) {
    valid_files <- c(valid_files, file)
  }
}

message("📋 Found ", length(valid_files), " files to process:")
for (file in valid_files) {
  message("   • ", basename(file))
}

if (length(valid_files) == 0) {
  stop("❌ No valid CSV files found to process")
}

# Process each file
message("\n🚀 Processing files...")
results_summary <- data.frame(
  File = character(),
  Genes = integer(),
  Significant = integer(),
  Up = integer(),
  Down = integer(),
  Zero_Genes_Imputed = integer(),
  Files_Created = character(),
  Status = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(valid_files)) {
  file <- valid_files[i]
  filename <- basename(file)
  output_prefix <- tools::file_path_sans_ext(filename)
  
  message("\n📄 [", i, "/", length(valid_files), "] Processing: ", filename)
  
  # Process file
  result <- tryCatch({
    create_volcano_plot(
      input_file = file,
      output_dir = OUTPUT_DIR,
      output_prefix = output_prefix,
      gene_col = GENE_COL,
      lfc_col = LFC_COL,
      pval_col = PVAL_COL,
      significance_col = SIGNIFICANCE_COL,
      p_threshold = P_THRESHOLD,
      lfc_threshold = LFC_THRESHOLD,
      x_limits = X_LIMITS,
      y_limit_max = Y_LIMIT_MAX,
      x_breaks_by = X_BREAKS_BY,
      y_breaks_by = Y_BREAKS_BY,
      dot_size = DOT_SIZE,
      dot_alpha = DOT_ALPHA,
      label_top_n = LABEL_TOP_N,
      colors = COLORS,
      create_compact = CREATE_COMPACT_VERSION,
      save_pdf = SAVE_PDF,
      save_png = SAVE_PNG,
      dpi = DPI,
      plot_width = PLOT_WIDTH,
      plot_height = PLOT_HEIGHT
    )
  }, error = function(e) {
    message("   ❌ Error processing file: ", e$message)
    return(NULL)
  })
  
  # Record results
  if (!is.null(result)) {
    data_size <- tryCatch({
      temp_data <- read.csv(file, stringsAsFactors = FALSE)
      nrow(temp_data)
    }, error = function(e) 0)
    
    results_summary <- rbind(results_summary, data.frame(
      File = filename,
      Genes = data_size,
      Significant = result$total_sig,
      Up = result$sig_up,
      Down = result$sig_down,
      Zero_Genes_Imputed = result$zero_genes_imputed,
      Files_Created = paste(result$files_created, collapse = ", "),
      Status = "✅ Success",
      stringsAsFactors = FALSE
    ))
    message("   ✅ Created: ", paste(result$files_created, collapse = ", "))
  } else {
    results_summary <- rbind(results_summary, data.frame(
      File = filename,
      Genes = 0,
      Significant = 0,
      Up = 0,
      Down = 0,
      Zero_Genes_Imputed = 0,
      Files_Created = "",
      Status = "❌ Failed",
      stringsAsFactors = FALSE
    ))
  }
}

# Save summary
summary_file <- file.path(OUTPUT_DIR, "volcano_plots_summary.csv")
write.csv(results_summary, summary_file, row.names = FALSE)

# Display summary table
print(results_summary)

#Enrichment plots for both Pseudobulk & Findmarkers

# === FOLDER PATHS - EDIT THESE ===
input_folder <- "/scratch/c.mpmrrp/sc_analysis/Jan_2025_AMB_single_cell_nuclei/output/new_align_March/28_May_2025/DE_results/less_threshold/"  # Change this to your CSV files folder
output_folder <- "/scratch/c.mpmrrp/sc_analysis/Jan_2025_AMB_single_cell_nuclei/output/new_align_March/28_May_2025/DE_results/less_threshold/enrichment_updated"  # Change this to your desired output folder

# === ANALYSIS PARAMETERS ===
padj_thresh <- 0.01           # Adjusted p-value threshold
lfc_thresh <- 0.585           # Log fold change threshold (0.585 = ~1.5 fold change)
gene_col <- "gene"            # Gene column name in your CSV
sig_col <- "Significance"     # Significance column name (if exists)
padj_col <- "p_val_adj"       # Adjusted p-value column name (updated for Seurat output)
lfc_col <- "avg_log2FC"       # Log fold change column name (updated for Seurat output)
top_terms <- 15              # Number of top terms to plot
create_plots <- TRUE         # Set to FALSE to skip plot generation
verbose <- TRUE              # Set to FALSE for less output

# === ENRICHMENT DATABASES ===
databases <- c(
  "GO_Biological_Process_2023",
  "GO_Cellular_Component_2023", 
  "GO_Molecular_Function_2023",
  "KEGG_2021_Human",
  "WikiPathway_2023_Human"
)


# Function to clean term names for plotting
clean_term <- function(term, max_length = 60) {
  # Remove database prefixes and clean up
  term <- gsub("^GO:\\d+~", "", term)
  term <- gsub("^KEGG_\\d+", "", term)
  term <- gsub("WP\\d+", "", term)
  term <- gsub("R-HSA-\\d+", "", term)
  
  # Truncate if too long
  if (nchar(term) > max_length) {
    term <- paste0(substr(term, 1, max_length-3), "...")
  }
  
  return(term)
}

# Function to detect and extract significant genes from CSV
extract_gene_lists <- function(df, gene_col, sig_col = NULL, padj_col = NULL, lfc_col = NULL, 
                              padj_thresh = 0.01, lfc_thresh = 0.585) {
  
  log_message(paste("Extracting gene lists from", nrow(df), "rows"))
  
  # Auto-detect gene column if not found
  if (!gene_col %in% colnames(df)) {
    log_message(paste("Gene column '", gene_col, "' not found. Available columns:", paste(colnames(df), collapse = ", ")), "WARNING")
    
    # Try common gene column names
    possible_gene_cols <- c("gene", "Gene", "GENE", "gene_name", "Gene_name", "symbol", "Symbol", "SYMBOL", "X")
    found_gene_col <- NULL
    
    for (possible_col in possible_gene_cols) {
      if (possible_col %in% colnames(df)) {
        found_gene_col <- possible_col
        break
      }
    }
    
    if (!is.null(found_gene_col)) {
      log_message(paste("Auto-detected gene column:", found_gene_col), "SUCCESS")
      gene_col <- found_gene_col
    } else {
      # Use first column as fallback
      gene_col <- colnames(df)[1]
      log_message(paste("Using first column as gene column:", gene_col), "WARNING")
    }
  }
  
  # Check if gene column exists now
  if (!gene_col %in% colnames(df)) {
    stop("No valid gene column found in data. Available columns: ", paste(colnames(df), collapse = ", "))
  }
  
  # Remove rows with missing gene names
  df <- df[!is.na(df[[gene_col]]) & df[[gene_col]] != "", ]
  log_message(paste("After removing missing genes:", nrow(df), "rows"))
  
  gene_lists <- list()
  
  # Method 1: Use existing significance column if available
  if (!is.null(sig_col) && sig_col %in% colnames(df)) {
    log_message(paste("Using existing significance column:", sig_col))
    
    # Get all significant genes
    sig_genes <- df[[gene_col]][!is.na(df[[sig_col]]) & 
                               df[[sig_col]] != "Not Significant" & 
                               df[[sig_col]] != ""]
    
    # Get up and down regulated genes
    up_genes <- df[[gene_col]][grepl("Up", df[[sig_col]], ignore.case = TRUE)]
    down_genes <- df[[gene_col]][grepl("Down", df[[sig_col]], ignore.case = TRUE)]
    
    if (length(sig_genes) > 0) gene_lists[["All_Significant"]] <- unique(sig_genes)
    if (length(up_genes) > 0) gene_lists[["Up_Regulated"]] <- unique(up_genes)
    if (length(down_genes) > 0) gene_lists[["Down_Regulated"]] <- unique(down_genes)
  }
  
  # Method 2: Use padj and log2FC thresholds if columns are available
  if (!is.null(padj_col) && padj_col %in% colnames(df) && 
      !is.null(lfc_col) && lfc_col %in% colnames(df)) {
    
    log_message(paste("Using thresholds: padj <", padj_thresh, ", |LFC| >=", lfc_thresh))
    
    # Create significance based on thresholds
    sig_idx <- !is.na(df[[padj_col]]) & !is.na(df[[lfc_col]]) & 
               df[[padj_col]] < padj_thresh & abs(df[[lfc_col]]) >= lfc_thresh
    
    up_idx <- sig_idx & df[[lfc_col]] > 0
    down_idx <- sig_idx & df[[lfc_col]] < 0
    
    sig_genes_thresh <- df[[gene_col]][sig_idx]
    up_genes_thresh <- df[[gene_col]][up_idx]
    down_genes_thresh <- df[[gene_col]][down_idx]
    
    # Add to gene lists (use threshold-based if no significance column found)
    if (length(gene_lists) == 0) {
      if (length(sig_genes_thresh) > 0) gene_lists[["All_Significant"]] <- unique(sig_genes_thresh)
      if (length(up_genes_thresh) > 0) gene_lists[["Up_Regulated"]] <- unique(up_genes_thresh)
      if (length(down_genes_thresh) > 0) gene_lists[["Down_Regulated"]] <- unique(down_genes_thresh)
    } else {
      # Add as additional lists for comparison
      if (length(sig_genes_thresh) > 0) gene_lists[["All_Significant_Threshold"]] <- unique(sig_genes_thresh)
      if (length(up_genes_thresh) > 0) gene_lists[["Up_Regulated_Threshold"]] <- unique(up_genes_thresh)
      if (length(down_genes_thresh) > 0) gene_lists[["Down_Regulated_Threshold"]] <- unique(down_genes_thresh)
    }
  }
  
  # Method 3: Fallback - use all genes if no significance info found
  if (length(gene_lists) == 0) {
    log_message("No significance information found, using all genes", "WARNING")
    all_genes <- unique(df[[gene_col]][!is.na(df[[gene_col]])])
    if (length(all_genes) > 0) {
      gene_lists[["All_Genes"]] <- all_genes
    }
  }
  
  # Log results
  for (list_name in names(gene_lists)) {
    log_message(paste("-", list_name, ":", length(gene_lists[[list_name]]), "genes"))
  }
  
  return(gene_lists)
}

# Function to run enrichment analysis
run_enrichment_analysis <- function(gene_lists, databases, output_dir, file_prefix, 
                                   top_terms = 15, create_plots = TRUE, verbose = FALSE) {
  
  log_message("Starting enrichment analysis", "PROCESS")
  
  enrichment_results <- list()
  all_results_combined <- data.frame()
  
  for (list_name in names(gene_lists)) {
    if (verbose) log_message(paste("Analyzing gene list:", list_name))
    
    enrichment_results[[list_name]] <- list()
    
    for (db_name in databases) {
      if (verbose) log_message(paste("  Database:", db_name))
      
      tryCatch({
        # Run enrichment
        enr_results <- enrichr(gene_lists[[list_name]], db_name)
        enr_df <- enr_results[[db_name]]
        
        if (!is.null(enr_df) && nrow(enr_df) > 0) {
          # Clean and process results
          enr_df <- enr_df[order(enr_df$P.value), ]
          enr_df$Term_Clean <- sapply(enr_df$Term, clean_term)
          enr_df$Gene_List <- list_name
          enr_df$Database <- db_name
          
          # Store results
          enrichment_results[[list_name]][[db_name]] <- enr_df
          
          # Save individual results
          enrich_fname <- file.path(output_dir, paste0(file_prefix, "_", list_name, "_", db_name, "_enrichment.csv"))
          write.csv(enr_df, enrich_fname, row.names = FALSE)
          
          # Add to combined results (top 10)
          top_results <- head(enr_df, 10)
          all_results_combined <- rbind(all_results_combined, top_results)
          
          # Create plots if requested
          if (create_plots && nrow(enr_df) > 0) {
            create_enrichment_plot(enr_df, list_name, db_name, output_dir, file_prefix, top_terms)
          }
          
          if (verbose) log_message(paste("    Found", nrow(enr_df), "terms"))
          
        } else {
          if (verbose) log_message("    No significant terms found")
        }
        
      }, error = function(e) {
        log_message(paste("Error in", db_name, ":", e$message), "ERROR")
      })
    }
  }
  
  # Save combined results
  if (nrow(all_results_combined) > 0) {
    combined_fname <- file.path(output_dir, paste0(file_prefix, "_ALL_enrichment_results.csv"))
    write.csv(all_results_combined, combined_fname, row.names = FALSE)
    
    # Create overview plot
    if (create_plots) {
      create_overview_plot(all_results_combined, output_dir, file_prefix)
    }
  }
  
  return(enrichment_results)
}

# Function to create enrichment plots (PDF OUTPUT)
create_enrichment_plot <- function(enr_df, list_name, db_name, output_dir, file_prefix, top_terms = 15) {
  
  top_n <- min(top_terms, nrow(enr_df))
  if (top_n == 0) return()
  
  # Determine colors based on gene list
  fill_color <- switch(list_name,
                      "All_Significant" = "#4292c6",
                      "Up_Regulated" = "#E31A1C", 
                      "Down_Regulated" = "#1F78B4",
                      "All_Genes" = "#666666",
                      "#4292c6")
  
  # Create bar plot
  p <- ggplot(head(enr_df, top_n),
              aes(x = reorder(Term_Clean, -log10(P.value)), y = -log10(P.value))) +
    geom_bar(stat = "identity", fill = fill_color, alpha = 0.8) +
    coord_flip() + 
    theme_pubr(base_size = 10) +
    labs(x = "Pathway", y = "-log10(P-value)", 
         title = paste0(db_name, " - ", gsub("_", " ", list_name)),
         subtitle = paste0("Analysis: ", gsub("_", " ", file_prefix))) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 8)
    )
  
  # Add significance line
  if (max(-log10(head(enr_df, top_n)$P.value)) > 1.3) {
    p <- p + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7)
  }
  
  # Save plot as PDF only
  plot_fname <- file.path(output_dir, paste0(file_prefix, "_", list_name, "_", db_name, "_barplot.pdf"))
  ggsave(plot_fname, p, width = 12, height = 8, device = "pdf", bg = "white")
}

# Function to create overview plot (PDF OUTPUT)
create_overview_plot <- function(all_results, output_dir, file_prefix, top_n = 20) {
  
  # Get top pathways (p < 0.01)
  top_pathways <- all_results %>%
    filter(P.value < 0.01) %>%
    arrange(P.value) %>%
    head(top_n)
  
  if (nrow(top_pathways) == 0) return()
  
  # Color by gene list
  colors_genelist <- c(
    "All_Significant" = "#4292c6", 
    "Up_Regulated" = "#E31A1C", 
    "Down_Regulated" = "#1F78B4",
    "All_Genes" = "#666666"
  )
  
  p <- ggplot(top_pathways, aes(x = reorder(Term_Clean, -log10(P.value)), 
                               y = -log10(P.value), 
                               fill = Gene_List)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = colors_genelist, name = "Gene List") +
    coord_flip() + 
    theme_pubr(base_size = 10) +
    labs(
      x = "Pathway", 
      y = "-log10(P-value)", 
      title = paste0("Top ", min(top_n, nrow(top_pathways)), " Enriched Pathways"),
      subtitle = paste0("Analysis: ", gsub("_", " ", file_prefix))
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom"
    ) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red", alpha = 0.7)
  
  # Save overview plot as PDF only
  overview_fname_pdf <- file.path(output_dir, paste0(file_prefix, "_overview.pdf"))
  ggsave(overview_fname_pdf, p, width = 14, height = 10, device = "pdf", bg = "white")
}

# Function to process a single CSV file
process_csv_file <- function(csv_file, output_dir) {
  
  log_message(paste("Processing file:", basename(csv_file)), "PROCESS")
  
  # Read CSV file
  tryCatch({
    df <- read.csv(csv_file, stringsAsFactors = FALSE)
    log_message(paste("Loaded", nrow(df), "rows,", ncol(df), "columns"))
  }, error = function(e) {
    log_message(paste("Failed to read", csv_file, ":", e$message), "ERROR")
    return(NULL)
  })
  
  # Extract gene lists
  gene_lists <- extract_gene_lists(
    df, 
    gene_col = gene_col,
    sig_col = if(sig_col %in% colnames(df)) sig_col else NULL,
    padj_col = if(padj_col %in% colnames(df)) padj_col else NULL,
    lfc_col = if(lfc_col %in% colnames(df)) lfc_col else NULL,
    padj_thresh = padj_thresh,
    lfc_thresh = lfc_thresh
  )
  
  if (length(gene_lists) == 0) {
    log_message("No gene lists extracted, skipping file", "WARNING")
    return(NULL)
  }
  
  # Create output subdirectory for this file
  file_prefix <- tools::file_path_sans_ext(basename(csv_file))
  file_output_dir <- file.path(output_dir, file_prefix)
  dir.create(file_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Run enrichment analysis
  results <- run_enrichment_analysis(
    gene_lists = gene_lists,
    databases = databases,
    output_dir = file_output_dir,
    file_prefix = file_prefix,
    top_terms = top_terms,
    create_plots = create_plots,
    verbose = verbose
  )
  
  # Generate summary report
  generate_report(file_output_dir, file_prefix, gene_lists)
  
  log_message(paste("Completed processing:", basename(csv_file)), "SUCCESS")
  return(results)
}

# Function to generate summary report
generate_report <- function(output_dir, file_prefix, gene_lists) {
  
  report_lines <- c(
    "PATHWAY ENRICHMENT ANALYSIS REPORT",
    "==================================",
    "",
    paste("File:", file_prefix),
    paste("Date:", Sys.Date()),
    paste("Time:", format(Sys.time(), "%H:%M:%S")),
    "",
    "PARAMETERS USED:",
    paste("- Adjusted p-value threshold:", padj_thresh),
    paste("- Log fold change threshold:", lfc_thresh),
    paste("- Gene column:", gene_col),
    paste("- Significance column:", sig_col),
    paste("- Top terms plotted:", top_terms),
    "",
    "GENE LISTS ANALYZED:",
    sapply(names(gene_lists), function(x) paste("-", x, ":", length(gene_lists[[x]]), "genes")),
    "",
    "DATABASES QUERIED:",
    paste("-", databases),
    "",
    "OUTPUT FILES:",
    "- Individual enrichment CSV files for each gene list × database",
    "- Bar plots for top enriched terms (PDF format)",
    "- Combined results file with top terms from all analyses", 
    "- Overview plot showing top pathways across databases (PDF format)",
    "- This summary report",
    "",
    paste("Total files created:", length(list.files(output_dir, recursive = TRUE)))
  )
  
  report_fname <- file.path(output_dir, paste0(file_prefix, "_enrichment_report.txt"))
  writeLines(report_lines, report_fname)
}

main <- function() {
  
  log_message("=== PATHWAY ENRICHMENT ANALYSIS (PDF OUTPUT) ===")
  log_message(paste("Version: RStudio 1.0 |", Sys.Date()))
  
  # Check if paths are set
  if (input_folder == "/scratch/c.mpmrrp/your_input_folder") {
    stop("❌ Please edit 'input_folder' path at the top of the script!")
  }
  
  if (output_folder == "/scratch/c.mpmrrp/your_output_folder") {
    stop("❌ Please edit 'output_folder' path at the top of the script!")
  }
  
  # Check if input directory exists
  if (!dir.exists(input_folder)) {
    stop("❌ Input directory not found: ", input_folder)
  }
  
  # Find CSV files with filtering - ONLY ALL_genes.csv files
  all_csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
  
  # Filter CSV files to only include those ending with "ALL_genes.csv"
  csv_files <- all_csv_files[grepl("ALL_genes\\.csv$", basename(all_csv_files), ignore.case = TRUE)]
  
  log_message(paste("Found", length(all_csv_files), "total CSV files in:", input_folder))
  log_message(paste("Filtered to", length(csv_files), "files ending with 'ALL_genes.csv'"))
  
  if (length(csv_files) < length(all_csv_files)) {
    skipped_files <- basename(all_csv_files[!all_csv_files %in% csv_files])
    log_message(paste("Skipped files (first 5):", paste(head(skipped_files, 5), collapse = ", ")))
    if (length(skipped_files) > 5) {
      log_message(paste("... and", length(skipped_files) - 5, "more files (focusing only on ALL_genes.csv files)"))
    }
  }
  
  if (length(csv_files) == 0) {
    stop("❌ No CSV files found in input directory")
  }
  
  # Create output directory
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  log_message(paste("Output directory:", output_folder))
  
  # Display configuration
  log_message("Configuration:")
  log_message(paste("- Input files:", length(csv_files)))
  log_message(paste("- Databases:", length(databases)))
  log_message(paste("- Significance thresholds: padj <", padj_thresh, ", |LFC| >=", lfc_thresh))
  log_message(paste("- Gene column:", gene_col))
  log_message(paste("- Create plots:", create_plots))
  log_message("- Output format: PDF only")
  
  # Process each file
  all_results <- list()
  success_count <- 0
  
  for (csv_file in csv_files) {
    result <- process_csv_file(csv_file, output_folder)
    if (!is.null(result)) {
      all_results[[basename(csv_file)]] <- result
      success_count <- success_count + 1
    }
  }
  
  # Final summary
  log_message("=== ANALYSIS COMPLETE ===", "SUCCESS")
  log_message(paste("Files processed successfully:", success_count, "/", length(csv_files)))
  log_message(paste("Results saved in:", output_folder))
  
  if (success_count > 0) {
    log_message("Key output files to check:")
    log_message("  • Individual enrichment CSV files")
    log_message("  • Combined results files (*_ALL_enrichment_results.csv)")
    log_message("  • Overview plots (*_overview.pdf)")
    log_message("  • Individual pathway plots (*_barplot.pdf)")
    log_message("  • Analysis reports (*_enrichment_report.txt)")
  }
  
  log_message("All plots saved as PDF files!", "SUCCESS")
}
# Execute the analysis
main()