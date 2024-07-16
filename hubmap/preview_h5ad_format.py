import os
import scanpy as sc

# Define the root folder
root_folder = r'C:\Users\yiju\Desktop'

# List all files in the root folder
all_files = os.listdir(root_folder)

# Filter for h5ad files
h5ad_files = [f for f in all_files if f.endswith('.h5ad')]

# Iterate over each h5ad file and preview it
for file in h5ad_files:
    file_path = os.path.join(root_folder, file)
    print(f"Processing file: {file_path}")
    adata = sc.read_h5ad(file_path)
    
    # Display a summary of the AnnData object
    print(adata)
    
    # Preview the first few rows of the data matrix
    print("Data matrix preview:")
    print(adata.X[:5, :5])  # Adjust the slice as needed
    
    # Preview cell metadata (first few rows)
    print("Cell metadata preview:")
    print(adata.obs.head())
    
    # Preview gene metadata (first few rows)
    print("Gene metadata preview:")
    print(adata.var.head())
    print("\n" + "="*40 + "\n")
