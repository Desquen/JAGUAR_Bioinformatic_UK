## Exploring T and NK cells in Jaguar Data — `gene_module_scoring_scRNAseq.R`

In single-cell RNA-seq analysis, once cells have been sequenced, clustered, and annotated, a natural next question is: *how active is a particular biological program in each cell?* That's exactly what this script addresses.

After cell type annotation (performed upstream using Azimuth), this script calculates **gene module scores** for each cell in the dataset. A module score is essentially a summary statistic that reflects how strongly a cell expresses a predefined set of genes — in this case, gene sets related to T cell and NK cell activation states. It's particularly useful for capturing continuous biological signals that don't map cleanly onto discrete clusters.

### What it does

The script takes a Seurat object (`.rds` or `.h5Seurat`) and a TSV file listing the gene sets of interest, then runs the following steps:

- Loads and validates input parameters from a YAML file or command line arguments
- Subsets the data to the cell types of interest: CD4 T, CD8 T, and NK cells, as annotated by Azimuth (`Azimuth.predicted.celltype.l1`)
- Calculates module scores per cell using Seurat's `AddModuleScore`, with an automatic fallback that adjusts the `nbin` parameter if the function fails (useful for small datasets)
- Normalizes scores within each Azimuth cell type group, so scores are comparable across modules
- Generates histograms showing the distribution of each module score per cell type
- Custom UMAP plots showing exactly where activation signatures turn on for CD4, CD8, and NK cells using custom colors
- Produces a violin plot showing the overlap between activation states (resting, early activation, and activated CD cells)
- Produces a boxplot comparing activation intensity across collection sites: Argentina, Chile, and Peru
- Saves a summary statistics table (min, mean, median, SD, max) per module and cell type
- Optionally saves results back into the Seurat object or exports them as a compressed CSV

### Inputs

| Parameter | Description |
|---|---|
| `query_dir` | Path to the Seurat object (`.rds` or `.h5Seurat`) |
| `gene_sets_tsv` | TSV file with gene set names and gene lists |
| `annotation_reference` | Metadata column used for grouping (e.g., `Azimuth.predicted.celltype.l1`) |
| `outdir` | Output directory |
| `save_toseurat` | If `TRUE`, saves results back into the Seurat object |
| `print_histograms` | If `TRUE`, generates histogram plots per module |
| `print_statistics` | If `TRUE`, exports a statistics summary table |

### Outputs

- `module_score_statistics.csv` — summary statistics per module and cell type
- `hist_module_score_<module>.png` — score distribution histograms
- `UMAP_<cell_type>_CD4modules.png` — UMAP plots showing the activation level for each cell type
- `CD4_Overlap_Violin.png` — violin plot of activation state overlap
- `FINAL_intensity_boxplot.png` — boxplot of activation intensity by site and cell type
- `*_mod.rds` or `*_gene_scores.csv.gz` — annotated Seurat object or scores table, depending on `save_toseurat`

### Usage
The script was adapted for execution within RStudio in a Microsoft Windows environment. To this end, the required variables were specifically declared, allowing for an interactive and documented workflow.

```r
# Example for testing
query_dir            <- "D:/HORAS NO LECTIVAS/PROYECTO JAGUAR/gene_module_score/jaguar_atlas_seurat.rds"
gene_sets_tsv        <- "D:/HORAS NO LECTIVAS/PROYECTO JAGUAR/gene_module_score/gene_set_list.tsv"
outdir               <- "D:/HORAS NO LECTIVAS/PROYECTO JAGUAR/gene_module_score/results/"
annotation_reference <- "Azimuth.predicted.celltype.l1"
save_toseurat        <- TRUE
print_histograms     <- TRUE
print_statistics     <- TRUE
```

### Dependencies

`Seurat`, `SeuratDisk`, `readr`, `tidyr`, `stringr`, `dplyr`, `ggplot2`, `rlang`, `patchwork`

Also requires a local `scripts/functions.R` file containing helper functions used by this script (`load_parameters`, among others).

---
## My Contributions

* **Gene Set Selection:** Defined the specific gene sets used to evaluate immune activation within the Jaguar Atlas dataset.
* **Custom Visualizations:** Developed specific plots including site-specific boxplots (Argentina, Chile, Peru), activation state violin plots, and UMAPs showing exactly where signatures turn on for CD4, CD8, and NK cells.
* **Script Optimization:** Implemented automated cell-type subsetting and adaptive `nbin` parameter adjustment to ensure script stability across different dataset sizes.


## Authors

* **Original Authors:** Julieth López ([julalopezcas@unal.edu.co](mailto:julalopezcas@unal.edu.co)) and Ania Lorenc ([al16@sanger.ac.uk](mailto:al16@sanger.ac.uk)).
* **Adapted and Optimized by:** Dámaris Esquén ([esquendamaris@gmail.com](mailto:esquendamaris@gmail.com) / [db36@sanger.ac.uk](mailto:db36@sanger.ac.uk)).
* **Project:** JAGUAR (2025).

---
