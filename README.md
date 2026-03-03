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
- `CD4_Overlap_Violin.png` — violin plot of activation state overlap
- `FINAL_intensity_boxplot.png` — boxplot of activation intensity by site and cell type
- `*_mod.rds` or `*_gene_scores.csv.gz` — annotated Seurat object or scores table, depending on `save_toseurat`

### Usage
The script was adapted for execution within RStudio in a Microsoft Windows environment. To this end, the required variables were specifically declared, allowing for an interactive and documented workflow.

```r
# RStudio/Windows execution configuration
query_dir            <- "D:/HORAS NO LECTIVAS/PROYECTO JAGUAR/gene_module_score/jaguar_atlas_seurat.rds"
gene_sets_tsv        <- "D:/HORAS NO LECTIVAS/PROYECTO JAGUAR/gene_module_score/gene_set_list.tsv"
outdir               <- "D:/HORAS NO LECTIVAS/PROYECTO JAGUAR/gene_module_score/results/"
annotation_reference <- "Azimuth.predicted.celltype.l1"
save_toseurat        <- TRUE
print_histograms     <- TRUE
print_statistics     <- TRUE

```bash
Rscript scripts/gene_module_scoring_scRNAseq.R <input_yaml> <output_directory>
```

Or with explicit arguments:

```bash
Rscript scripts/gene_module_scoring_scRNAseq.R \
  --query_dir jaguar_atlas.rds \
  --gene_sets gene_set_list.tsv \
  --annotation_reference Azimuth.predicted.celltype.l1 \
  --save_toseurat TRUE \
  --print_histograms TRUE \
  --print_statistics TRUE \
  --outdir results/
```



> **Note:** If running on the Sanger Institute cluster, load the required environment first:
> `module load HGI/softpack/users/hn4/seurat5/7.0`

### Dependencies

`Seurat`, `SeuratDisk`, `readr`, `tidyr`, `stringr`, `dplyr`, `ggplot2`, `rlang`

Also requires a local `scripts/functions.R` file containing helper functions used by this script (`load_parameters`, among others).
