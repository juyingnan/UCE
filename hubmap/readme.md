# README.md for HubMAP Folder

## Overview

This repository contains a collection of Python scripts designed to process and analyze single-cell RNA sequencing (scRNA-seq) data using Universal Cell Embedding (UCE) and UMAP visualizations. This README provides an overview of each script and useful external resources related to the HuBMAP project.

## Scripts

### Data Preprocessing and Utilities
#### `data_fix.py`
This script processes `.h5ad` files and corresponding annotation CSV files within a specified directory. It ensures that NaN values in the 'hugo_symbol' column are replaced with a placeholder and that cell type labels are mapped correctly. This file works for Lung data only.

#### `label_finder.py`
This script scans a directory for annotation CSV files and extracts unique cell type labels (CL_Labels) for further analysis and categorization.

### UMAP Calculation and Visualization
#### `calculate_umap.py`
This script processes `.h5ad` files to calculate UMAP embeddings based on UCE embeddings if they do not already exist. The processed files are saved back to their original location to avoid duplicated calculation in the future to save time.

#### `preview_umap.py`
This script visualizes UMAP embeddings from a specified `.h5ad` file. It generates UMAP plots colored by predicted cell type labels and saves the visualization as a PNG file.

#### `preview_umap_cat.py`
This script generates UMAP visualizations categorized by cell type using different strategies: 'CL_Label', 'high_level', and 'detailed_level'. It also calculates the Calinski-Harabasz Index for each strategy and saves the visualizations.

#### `preview_umap_gtex.py`
This script is for visualizing GTEx data UMAP embeddings. It allows for the use of precomputed UMAP embeddings or recalculates them based on UCE embeddings.

#### `dendrogram.py`
This script computes hierarchical clustering on UMAP coordinates and generates dendrogram plots to visualize the clustering of UMAP embeddings based on cell type labels.

### H5AD File Preview
#### `preview_h5ad_dim.py`
This script provides a quick check of the dimensions of raw input data and UCE embeddings within a specified `.h5ad` file.

#### `preview_h5ad_format.py`
This script previews the structure and contents of `.h5ad` files in a specified directory, including data matrices and metadata.

### Command-Line Instructions
#### `commands.txt`
This file contains command-line instructions for running UCE evaluations and visualizations on various datasets, including lung and GTEx data.

## Useful Resources

- **UCE Paper**: [Universal Cell Embedding Paper](https://www.biorxiv.org/content/10.1101/2023.11.28.568918v1.full.pdf)
- **UCE Model Files**: [UCE Model Files on Figshare](https://figshare.com/articles/dataset/Universal_Cell_Embedding_Model_Files/24320806?file=43423236)
- **Human Brain Paper**: [Human Brain Single-cell Paper](https://www.biorxiv.org/content/10.1101/2022.10.12.511898v1.full.pdf)
- **Mouse Brain Paper**: [Mouse Brain Single-cell Paper](https://www.biorxiv.org/content/10.1101/2020.07.02.184051v1.full)
- **Lung Data Spreadsheet**: [Lung Data Google Spreadsheet](https://docs.google.com/spreadsheets/d/1TdCODl0UpsRZu-erH6t-noTyPeI-NKr7OmF3WTjNFeY/edit?gid=0#gid=0)
- **GTEx Data Portal**: [GTEx Data Portal](https://gtexportal.org/home/downloads/adult-gtex/single_cell)
- **GTEx Paper**: [GTEx Single-cell Paper](https://www.science.org/doi/10.1126/science.abl4290)
