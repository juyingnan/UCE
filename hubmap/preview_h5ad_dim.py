import scanpy as sc

# Load AnnData object
file_path = r'c:\users\yiju\desktop\10k_pbmcs_proc_uce_adata.h5ad'
adata = sc.read_h5ad(file_path)

# Check the dimensions of the raw input data
input_shape = adata.X.shape
print(f"Input data shape (cells x genes): {input_shape}")

# Check the dimensions of the UCE embeddings
uce_shape = adata.obsm["X_uce"].shape
print(f"UCE embeddings shape (cells x embedding dimensions): {uce_shape}")
