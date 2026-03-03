# """
# Author: Julieth López, Ania Lorenc, adapted by Dámaris Esquén
# Date: 10/10/2025
# Project: JAGUAR
# Contact: julalopezcas@unal.edu.co, al16@sanger.ac.uk
# Description: This notebook receive as input a seurat object and some gene lists, then create module scores for each cell
# """
setwd("D:/HORAS NO LECTIVAS/PROYECTO JAGUAR/gene_module_score/")

# packages
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(readr)
  library(tidyr)
  library(stringr)
  library(dplyr) #database manipulation
  library(ggplot2) # for plotting
  library(patchwork)
  library(rlang) # to read character as symbol
  source("scripts/functions.R")
})


# Obtain input parameters
#args = commandArgs(trailingOnly=TRUE)
# Run in console
args <- c("params_jaguar.yaml", "D:/HORAS NO LECTIVAS/PROYECTO JAGUAR/gene_module_score/results/")

usage <- "Usage: It can be run by taking a parameters YAML file, or command line arguments as follows.\n
    Rscript scripts/27.gen_module_score.R <input_yaml> <output_directory>\n
    or\n
    Rscript scripts/27.gen_module_score.R --query_dir <query.rds> --save_toseurat <TRUE/FALSE> --gene_sets <gene_sests.tsv> --annotation_reference <Azimuth> --print_histograms <TRUE/FALSE> --print_statistics <TRUE/FALSE> --outdir <outdir>\n
    
    Input parameters:\n
    - query_dir: Query Seurat object (.rds/h5Seurat)\n
    - gene_sets: A list of gene sets to use for module scoring (default provided in script)\n
    - annotation_reference: Variable in metadata to use for subsetting (e.g., Azimuth)\n
    - outdir: Output directory\n

    Optional parameters:\n
    - save_toseurat: If set, save results back into query Seurat as *_mod.rds\n
    
    Outdir: Where results will be saved, e.g, results\n

    Expected outputs:\n
    - csv with modules results for each cell, contains the columns: barcodes, predicted labels, prediction scores, and mapping scores\n
    - if --save_toseurat, is retorned query_mod.rds with predictions added to metadata\n
    - if print_histograms==TRUE, histograms per module score per annotation_reference level\n
    - if print_statistics==TRUE, table with statistics per module score per annotation_reference level\n

    ** if run in sanger server use: module load HGI/softpack/users/hn4/seurat5/7.0"

# verification of arguments
if (length(args) < 2) {
    stop(usage)}

# Load parameters from command line or YAML file
paramsDict <- load_parameters(args)
outdir <- as.character(args[length(args)])

# Assign parsed arguments to variables
query_dir = as.character(paramsDict["query_dir"])
gene_sets_tsv = as.character(paramsDict["gene_sets_tsv"])
save_toseurat = as.character(paramsDict["save_toseurat"])
annotation_reference = as.character(paramsDict["annotation_reference"])

# Example for testing
query_dir <- "D:/HORAS NO LECTIVAS/PROYECTO JAGUAR/gene_module_score/jaguar_atlas_seurat.rds.rds"
gene_sets_tsv <- "D:/HORAS NO LECTIVAS/PROYECTO JAGUAR/gene_module_score/gene_set_list.tsv"
outdir <- "D:/HORAS NO LECTIVAS/PROYECTO JAGUAR/gene_module_score/results/"
annotation_reference <- "Azimuth.predicted.celltype.l1"
save_toseurat <- TRUE
print_histograms <- TRUE
print_statistics <- TRUE

# Check parameters
if (is.null(query_dir) ) stop("❌ You must provide  --query_dir \n\n", usage)

if (is.null(save_toseurat) || save_toseurat == "" || is.na(save_toseurat)) {
  save_toseurat <- FALSE
}

if (is.null(outdir)) {
  warning("⚠️ You did not provide --outdir. Using current directory instead.\n")
  outdir <- getwd()
}

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


# --- Load query Seurat object ---
# Function to read Seurat object based on file extension
read_seurat_file <-function(file) {
  if (grepl("\\.h5Seurat$", file, ignore.case = TRUE)) {
    LoadH5Seurat(file)
  } else if (grepl("\\.rds$", file, ignore.case = TRUE)) {
    readRDS(file)
  } else {
    stop("❌ Input must be a .h5Seurat or .rds file")
  }
}

# Read TSV (assuming first column = list name, second = gene, third = called gene sets)
read_gene_sets <-function(gene_sets_tsv, query) {
  gene_sets_df <- read.delim(gene_sets_tsv, sep = "\t")

  # organice gene sets from tsv in lists
  genes_list <- lapply(gene_sets_df[[2]], function(x) unlist(strsplit(as.character(x), ",")))
  names(genes_list) <- gene_sets_df[[1]]
  called_gene_sets <- lapply(gene_sets_df[[3]], function(x) {terms = unlist(strsplit(as.character(x), ","))
                                                  complete_genes_in_set <- c()
                                                  for (t in terms){
                                                    if (length(t)>0){
                                                      genes_in_set <- tryCatch({
                                                                        if (t %in% names(gene_sets)) gene_sets[[t]] else get(as.name(t))
                                                                      }, error = function(e) {NULL})}
                                                    else {
                                                      genes_in_set = NULL}
                                                    complete_genes_in_set <- c(complete_genes_in_set, genes_in_set)
                                                  }
                                                  return(complete_genes_in_set)
                                                  })

  # join both lists
  gene_sets <- Map(c,genes_list,called_gene_sets)

  # verify if genes are in the dataset
  for (name in names(gene_sets)) {
    if (any(!(gene_sets[[name]] %in% rownames(query)))){
      warning(paste("Some genes in gene set", name, "are not present in the dataset:", 
                    paste(gene_sets[[name]][!(gene_sets[[name]] %in% rownames(query))], collapse = ", ")))
      print(paste0("Absent genes were deleted from gene set ", name))
      gene_sets[[name]] <- gene_sets[[name]][gene_sets[[name]] %in% rownames(query)]
    }
  }

  # Check gene sets
  print("Gene sets loaded:")
  print(names(gene_sets))

  return(gene_sets)
}


# --- Calculate gene modules scores ---
# Try AddModuleScore with default settings, if it fails, decrease nbin by 1 until it works or nbin==1
calculate_module_score <- function(query, gene_sets, all_metadata = TRUE){

  query_mod <- NULL
  nbin_val <- min(24, ncol(query))
  success <- FALSE

    while (nbin_val >= 1 && !success) {
    try_result <- try({
      query_mod <- AddModuleScore(
        object   = query,
        features = gene_sets,
        name     = names(gene_sets),
        nbin     = nbin_val
      )
    }, silent = TRUE)
    if (!inherits(try_result, "try-error") && !is.null(query_mod)) {
      success <- TRUE
      message(paste("AddModuleScore succeeded with nbin =", nbin_val))
    } else {
      message(paste("AddModuleScore failed with nbin =", nbin_val, ". Trying nbin =", nbin_val - 1))
      nbin_val <- nbin_val - 1
    }
  }

  if (!success) stop("AddModuleScore failed for all attempted nbin values.")

  # selecting new columns added to meta.data
  module_cols <- setdiff(colnames(query_mod@meta.data), colnames(query@meta.data))
  
  # Remove trailing digits only from these columns
  added_columns <- sub("\\d+$", "", module_cols)
  
  # Rename in meta.data
  colnames(query_mod@meta.data)[match(module_cols, colnames(query_mod@meta.data))] <- added_columns

  # Return only new columns or all metadata
  if (all_metadata == FALSE)
    scores <- query_mod@meta.data[added_columns]
  else
    scores <- query_mod@meta.data

  return(scores)
}


# --- Graphics ---
# This could be done for clusters. And could be not necessary for azimuth levels
hist_module_score <- function(azimuth_module_scores, variable_group_by, added_columns, outdir){
  for (module in added_columns){
    p <- azimuth_module_scores %>%
        group_by(!!sym(variable_group_by)) %>%
        ggplot(aes(x = !!sym(module), fill = !!sym(variable_group_by))) +
          geom_histogram(bins = 30, color = "white", alpha = 0.7) +
          facet_wrap(vars(!!sym(variable_group_by)), scales = "free_y") +
          labs(title = paste0("Distribution module score:",module),
              x = "Module Score")

    # Save plot
    ggsave(paste0(outdir,"/hist_module_score_",module,".png"), plot = p, width = 10, height = 6, dpi = 300)
  }
  return("Plots saved")
}


# --- Table with statisics ---
# This could be done for clusters, not necessary if is done for azimuth levels
statistics_module_score <- function(query_mod, variable_group_by, added_columns){
  sum_table <- query_mod@meta.data %>%
    group_by(!!sym(variable_group_by)) %>%
    summarise(across(all_of(added_columns), 
              list(min = min, mean = mean, median = median, sd = sd, max = max), .names = "{col}_{fn}")) %>% 
    mutate_if(is.numeric, round, 3)

  return(sum_table)
}


# --- Apply functions ---
# Load RDS Seurat object
cat("Loading query:", query_dir, "\n")
query_all <- read_seurat_file(query_dir)
query <- subset(query_all, subset = orig.ident == "ATLAS09a" & anydoublet == "singlet")
print(table(query$orig.ident))

# Load gene sets
cat("Loading gene sets:", gene_sets_tsv, "\n")
gene_sets <- read_gene_sets(gene_sets_tsv, query)

# Only in this case, because we have few cells in some populations
table(query$Azimuth.predicted.celltype.l1)
query <- subset(query, subset = Azimuth.predicted.celltype.l1 %in% c("CD4 T", "CD8 T", "NK"))
azimuth_levels <- unique(query@meta.data[[annotation_reference]])

# Calculate module scores for groups of cells for specific variable like "Azimuth.predicted.celltype.l1"
cat("Calculating module score")
azimuth_levels = unique(query@meta.data[[annotation_reference]])
azimuth_module_scores = list()
for (level in azimuth_levels){
  print(level)
  subset_query <- subset(query, subset = `Azimuth.predicted.celltype.l1` == level)
  subset_module_scores <- calculate_module_score(subset_query, gene_sets[level], all_metadata = FALSE)
  azimuth_module_scores[[level]] <- subset_module_scores
}
names(azimuth_module_scores)

# Calculate module scores for CD4, CD8 and NK cells
filtered_gene_sets <- gene_sets
module_score_metadata <- calculate_module_score(query, filtered_gene_sets)
names(module_score_metadata)

# Scale data taking as reference subsetted cells (azimuth levels), for L1 level
# If not reference, module scores can be scaled per module and thresholds can be selected by module. However, modules will not be comparable between them
modules_to_run <- names(gene_sets)
for (module in modules_to_run){
  print(module)
  if (!is.null(azimuth_module_scores[[module]]) && is.numeric(azimuth_module_scores[[module]])) {
    max_value <- max(abs(azimuth_module_scores[[module]]))
    azimuth_module_scores[[module]] <- azimuth_module_scores[[module]]/max_value
    module_score_metadata[module] <- module_score_metadata[module]/max_value
  }
}

# Rules to select label
for (module in modules_to_run){
  print(module)
  if (!is.null(azimuth_module_scores[[module]])) {
    max_value <- max(abs(azimuth_module_scores[[module]]))
    azimuth_module_scores[[module]] <- azimuth_module_scores[[module]]/max_value
    module_score_metadata[module] <- module_score_metadata[module]/max_value
  }
}

# Add module scores to query metadata
added_columns <- setdiff(colnames(module_score_metadata), colnames(query@meta.data))
print(paste("Metadata now: ", colnames(query@meta.data)))
print(paste("Added columns: ", added_columns))
head(module_score_metadata[added_columns])

# 1. Identify what's new (You already did this)
added_columns <- setdiff(colnames(module_score_metadata), colnames(query@meta.data))
# 2. ACTUALLY ATTACH the data to the object
query@meta.data <- cbind(query@meta.data, module_score_metadata[added_columns])
# Convert list to table so the plotting function doesn't crash
if (is.list(azimuth_module_scores)) {
  azimuth_module_scores <- dplyr::bind_rows(azimuth_module_scores, .id = "Azimuth.predicted.celltype.l1")
}
# Identify only numeric columns for the functions
plot_columns <- added_columns[!grepl("_label$", added_columns)]

if (print_histograms == TRUE) {
  available_cols <- intersect(plot_columns, colnames(azimuth_module_scores))
  hist_module_score(azimuth_module_scores, "Azimuth.predicted.celltype.l1", available_cols, outdir)
  cat("✅ Saved histograms in:", outdir, "\n")
}

if (print_statistics==TRUE){
  stats_table <- statistics_module_score(query, "Azimuth.predicted.celltype.l1", plot_columns)
  write.csv(stats_table, file = paste0(outdir,"/module_score_statistics.csv"), row.names = FALSE)
  cat("Saved statistics table:", paste0(outdir,"/module_score_statistics.csv"), "\n")
}

# --- Save results ---
if (!is.null(save_toseurat) && save_toseurat) {
  out_file <- sub("\\.rds$", "_mod.rds", query_dir)
  saveRDS(query, out_file)
  cat("✅ Saved annotated Seurat object:", out_file, "\n")
} else {
  out_file <- sub("\\.rds$", "_gene_scores.csv.gz", query_dir)
  
  df <- cbind(Barcode = colnames(query), query@meta.data[, added_columns])
  write.csv(df, gzfile(out_file), row.names = TRUE)
  cat("✅ Saved predictions table:", out_file, "\n")
}

# Quick check of the new columns
head(query@meta.data[, added_columns])


# Visualizacion of gene set selected related activation
# UMAP plot
# 1. Config
base_modules <- c("resting_CD4_T", "early_activation_CD4_T", "activated_CD4_T")
cell_types   <- c("CD4 T", "CD8 T", "NK")
palette_3    <- c("#088082", "darkgray", "#ab158f")

# Identify actual metadata columns
actual_modules <- sapply(base_modules, function(x) grep(paste0("^", x), colnames(query@meta.data), value = TRUE)[1])

for (ct in cell_types) {
  # Subset by cell type
  sub_obj <- subset(query, subset = Azimuth.predicted.celltype.l1 == ct)
  
  # Calculate shared limits for consistent comparison
  limits <- range(sub_obj@meta.data[, actual_modules], na.rm = TRUE)
  
  plots <- lapply(1:3, function(i) {
    FeaturePlot(
      object = sub_obj, 
      features = actual_modules[i], 
      reduction = "umap", 
      raster = TRUE, 
      order = TRUE, 
      pt.size = 1.5
    ) + 
      theme_minimal() + 
      labs(title = base_modules[i]) +
      scale_color_gradientn(
        colors = palette_3, 
        limits = limits,
        name = "Module Score", # <--- Legend Title
        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
      ) + 
      theme(
        legend.position = "right", 
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.line = element_line(color = "gray30", linewidth = 0.5),
        axis.ticks = element_line(color = "gray30")
      )
  })
  
  # Layout and Export
  combined <- wrap_plots(plots, ncol = 3) + 
    plot_layout(guides = "collect") +
    plot_annotation(title = paste("Population:", ct), theme = theme(plot.title = element_text(size = 16, face = "bold")))
  
  ggsave(file.path(outdir, paste0("UMAP_", gsub(" ", "_", ct), "_CD4modules.png")), 
         plot = combined, width = 16, height = 6, dpi = 600)
  
  rm(sub_obj, combined); gc()
}


# Visualization of gene set selected related activation 
target_stages <- c("resting_CD4_T", "early_activation_CD4_T", "activated_CD4_T")
target_sites <- c("Argentina", "Chile", "Peru")

# Data preparation
query$overlap_count <- rowSums(query@meta.data[, target_stages] > 0)

df_overlap_viz <- query@meta.data %>%
  filter(overlap_count > 0) %>% 
  select(all_of(target_stages), overlap_count) %>%
  pivot_longer(cols = all_of(target_stages), 
               names_to = "Score_Type", 
               values_to = "Intensity")

# 2. Violin Plot of overlap level
plot_scores_overlap <- ggplot(df_overlap_viz, aes(x = Score_Type, y = Intensity, fill = Score_Type)) +
  geom_violin(alpha = 0.7, trim = FALSE, scale = "width") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  facet_wrap(~overlap_count, 
             labeller = as_labeller(c("1" = "1 State", "2" = "2 States", "3" = "3 States"))) +
  theme_bw() +
  scale_fill_manual(values = c("resting_CD4_T" = "#619CFF", 
                               "early_activation_CD4_T" = "#0fa348", 
                               "activated_CD4_T" = "#b82d14")) +
  labs(title = "Signature Intensity by Overlap Level",
       y = "Module Score", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none")

print(plot_scores_overlap)
ggsave(paste0(outdir, "/CD4_Overlap_Violin.png"), plot = plot_scores_overlap, width = 10, height = 7)

# cleaning by high status 
# clean regarding high score 
query$dominant_status <- target_stages[max.col(query@meta.data[, target_stages], ties.method = "first")]
# clean of cells with NA value
df_clean <- query@meta.data %>%
  filter(site %in% target_sites, !is.na(dominant_status))

# Boxplot of gene sets 
plot_data_intensity <- df_clean %>%
  select(Azimuth.predicted.celltype.l1, site, all_of(target_stages)) %>%
  pivot_longer(cols = all_of(target_stages),
               names_to = "Stage",
               values_to = "Score")

plot_data_intensity$Stage <- factor(plot_data_intensity$Stage, levels = target_stages)

final_boxplot <- ggplot(plot_data_intensity, aes(x = Stage, y = Score, fill = site)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, median.linewidth = 1) +
  facet_grid(Azimuth.predicted.celltype.l1 ~ site) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  scale_fill_manual(values = c("Argentina" = "#619CFF", "Chile" = "#0fa348", "Peru" = "#b82d14")) +
  labs(title = "CD4 T cell Activation (Intensity Score)", 
       y = "Module Score", x = "Activation Stage")
final_boxplot
ggsave(paste0(outdir, "/FINAL_intensity_boxplot.png"), plot = final_boxplot, width = 12, height = 8)

