import scanpy as sc
import matplotlib.pyplot as plt

# Load the primary AnnData file (with UCE embeddings and cell type labels)
file_path = r'C:\Users\yiju\Desktop\temp_data\HBM265.FQWZ.38410k_pbmcs_proc_fixed_uce_adata.h5ad'
adata = sc.read_h5ad(file_path)

# Check if 'X_uce' is in the obsm slot of the primary file
if 'X_uce' in adata.obsm:
    print("UCE embeddings found. Visualizing...")

    # Check if 'predicted_label' is in the obs slot
    if 'predicted_label' in adata.obs:
        print("Cell type labels found. Proceeding with UMAP visualization...")
        
        # Run UMAP on the UCE embeddings
        sc.pp.neighbors(adata, use_rep='X_uce')  # Use UCE embeddings for neighborhood graph
        sc.tl.umap(adata)  # Run UMAP

        # Set up the figure with a larger size
        fig, ax = plt.subplots(figsize=(12, 8))

        # Plot the UMAP with cell type information
        sc.pl.umap(adata, color='predicted_label', ax=ax, show=False)  # Adjust 'color' based on cell types

        # Adjust the layout to make space for the legend
        plt.tight_layout(rect=[0, 0, 0.8, 1])  # Adjust the right padding to make room for the legend

        # Show the plot
        plt.show()
    else:
        print("Cell type labels not found in the primary file.")
else:
    print("UCE embeddings not found in the primary file.")
