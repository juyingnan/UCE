import scanpy as sc
import matplotlib.pyplot as plt
import os
import argparse
from sklearn.metrics import calinski_harabasz_score

def visualize_umap(h5ad_file_path, original=False):
    # Load the AnnData file (with UCE embeddings and cell type labels)
    adata = sc.read_h5ad(h5ad_file_path)

    annotation_str = 'Broad cell type (numbers)'

    # Check if 'X_uce' is in the obsm slot of the file
    if 'X_uce' in adata.obsm:
        print("UCE embeddings found. Visualizing...")

        # Check if 'predicted_label' is in the obs slot
        if annotation_str in adata.obs:
            print("Cell type labels found. Proceeding with UMAP visualization...")

            if not original:
                # Run UMAP on the UCE embeddings
                print("Running UMAP calculation...")
                sc.pp.neighbors(adata, use_rep='X_uce')  # Use UCE embeddings for neighborhood graph
                sc.tl.umap(adata)  # Run UMAP
            else:
                print("Using original UMAP, skip calculating")

            # Set up the figure with a larger size
            fig, ax = plt.subplots(figsize=(12, 8))

            # Calculate Calinski-Harabasz Index
            calinski_harabasz_avg = calinski_harabasz_score(adata.obsm['X_umap'], adata.obs[annotation_str])
            print(f"Calinski-Harabasz Index for 'annotation': {calinski_harabasz_avg:.3f}")

            # Add Calinski-Harabasz Index to the plot
            plt.text(0.5, -0.1, f'Calinski-Harabasz Index: {calinski_harabasz_avg:.3f}', ha='center', va='center', transform=ax.transAxes, fontsize=12)

            # Plot the UMAP with cell type information
            sc.pl.umap(adata, color=annotation_str, ax=ax, show=False)  # Adjust 'color' based on cell types

            # Adjust the layout to make space for the legend
            # plt.tight_layout(rect=[0, 0, 0.25, 1])  # Adjust the right padding to make room for the legend
            ax.legend()

            # Save the plot to the same folder as the h5ad file
            output_file_path = os.path.join(os.path.dirname(h5ad_file_path), 'original_umap.png' if original else 'uce_umap.png')
            plt.savefig(output_file_path)

            # Show the plot
            plt.show()
        else:
            print("Cell type labels not found in the file.")
    else:
        print("UCE embeddings not found in the file.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize UMAP with UCE embeddings and cell type labels.')
    parser.add_argument('--h5ad_path', type=str, required=True, help='Full path to the h5ad file.')
    parser.add_argument('--original', type=bool, required=False, default=False, help='Full path to the h5ad file.')
    args = parser.parse_args()

    visualize_umap(args.h5ad_path, args.original)
