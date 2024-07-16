import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

# Load the h5ad file
file_path = r'C:\Users\yiju\Desktop\chicken_heart.h5ad'
adata = sc.read_h5ad(file_path)

# Define the higher-level categories based on two strategies

# Strategy 1: Functional and Structural Roles
cell_type_to_functional = {
    'Cardiomyocytes-1': 'Cardiac Muscle Cells',
    'Cardiomyocytes-2': 'Cardiac Muscle Cells',
    'MT-enriched cardiomyocytes': 'Cardiac Muscle Cells',
    'Immature myocardial cells': 'Cardiac Muscle Cells',
    'Fibroblast cells': 'Structural Support Cells',
    'Mural cells': 'Structural Support Cells',
    'Vascular endothelial cells': 'Vascular Cells',
    'Endocardial cells': 'Vascular Cells',
    'Dendritic cells': 'Immune Cells',
    'Macrophages': 'Immune Cells',
    'Erythrocytes': 'Blood Cells',
    'Epi-epithelial cells': 'Epithelial Cells',
    'Epi-mesenchymal cells': 'Epithelial Cells',
    'Valve cells': 'Valve Cells',
    'TMSB4X high cells': 'Specialized or State-Specific Cells'
}

# Strategy 2: Cell Lineage and Type
cell_type_to_lineage = {
    'Cardiomyocytes-1': 'Cardiomyocytes (Heart Muscle Cells)',
    'Cardiomyocytes-2': 'Cardiomyocytes (Heart Muscle Cells)',
    'MT-enriched cardiomyocytes': 'Cardiomyocytes (Heart Muscle Cells)',
    'Immature myocardial cells': 'Cardiomyocytes (Heart Muscle Cells)',
    'Fibroblast cells': 'Supportive Cells (Stromal Cells)',
    'Mural cells': 'Supportive Cells (Stromal Cells)',
    'Vascular endothelial cells': 'Endothelial Cells',
    'Endocardial cells': 'Endothelial Cells',
    'Dendritic cells': 'Immune Cells',
    'Macrophages': 'Immune Cells',
    'Erythrocytes': 'Blood Cells',
    'Epi-epithelial cells': 'Epithelial Cells',
    'Epi-mesenchymal cells': 'Epithelial Cells',
    'Valve cells': 'Valve Cells',
    'TMSB4X high cells': 'Other Specific Cells'
}

# Map the original cell types to the new higher-level categories
adata.obs['functional_superclass'] = adata.obs['cell_type'].map(cell_type_to_functional)
adata.obs['lineage_superclass'] = adata.obs['cell_type'].map(cell_type_to_lineage)

# Check if 'X_umap' and 'X_uce' are in the obsm slot
embeddings = []
if 'X_umap' in adata.obsm:
    embeddings.append('X_umap')
if 'X_uce' in adata.obsm:
    embeddings.append('X_uce')

# Define the strategies
strategies = ['cell_type', 'functional_superclass', 'lineage_superclass']

for embedding in embeddings:
    print(f"Visualizing using embedding: {embedding}")
    for strategy in strategies:
        print(f"Visualizing based on {strategy}")

        # Set up the figure with a larger size
        fig, ax = plt.subplots(figsize=(12,8))

        if embedding == 'X_umap':
            # Plot the UMAP embeddings with custom figure and axis
            sc.pl.embedding(adata, basis=embedding, color=strategy, ax=ax, show=False)  # show=False to avoid immediate display
        else:
            # Run UMAP on the UCE embeddings
            sc.pp.neighbors(adata, use_rep='X_uce')  # Use UCE embeddings for neighborhood graph
            sc.tl.umap(adata)  # Run UMAP

            # Plot the UMAP
            sc.pl.umap(adata, color=strategy, ax=ax, show=False)  # Adjust 'color' based on your dataset

        # Adjust the layout to make space for the legend
        plt.tight_layout(rect=[0, 0, 0.8, 1])  # Adjust the right padding to make room for the legend

        # Save the plot
        plt.savefig(f'umap_{embedding}_{strategy}.png')

        # Show the plot
        # plt.show()
