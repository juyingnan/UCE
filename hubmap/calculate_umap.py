import scanpy as sc
import os

def process_and_save_umap(h5ad_file_path):
    # Load the AnnData file (with UCE embeddings and cell type labels)
    adata = sc.read_h5ad(h5ad_file_path)

    # Check if 'X_uce' is in the obsm slot of the file
    if 'X_uce' in adata.obsm:
        if 'X_umap' in adata.obsm:
            print(f"UMAP embeddings already exist in {h5ad_file_path}. Skipping...")
            return  # Skip the UMAP computation
        else:
            print(f"UCE embeddings found in {h5ad_file_path}. Running UMAP...")

            # Run UMAP on the UCE embeddings
            sc.pp.neighbors(adata, use_rep='X_uce')  # Use UCE embeddings for neighborhood graph
            sc.tl.umap(adata)  # Run UMAP

            # Save the UMAP result back to the AnnData file
            adata.write(h5ad_file_path)
            print(f"UMAP result saved to {h5ad_file_path}.")
    else:
        print(f"UCE embeddings not found in {h5ad_file_path}.")

def process_all_subfolders(root_folder):
    for subdir, _, files in os.walk(root_folder):
        for file in files:
            if file.endswith('_uce_adata.h5ad'):
                h5ad_file_path = os.path.join(subdir, file)
                process_and_save_umap(h5ad_file_path)

if __name__ == "__main__":
    root_folder = r"C:/Users/yiju/Desktop/temp_data/"
    process_all_subfolders(root_folder)
