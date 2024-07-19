import scanpy as sc
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import pandas as pd
import preview_umap_cat

# Load the AnnData file (with UMAP and cluster information)
file_path = r'C:\Users\yiju\Desktop\temp_data\HBM762.RPDR.282\HBM762.RPDR.282_data_fixed_uce_adata.h5ad'
adata = sc.read_h5ad(file_path)

adata.obs['high_level'] = adata.obs['CL_Label'].map(preview_umap_cat.high_level)
adata.obs['detailed_level'] = adata.obs['CL_Label'].map(preview_umap_cat.detailed_level)

strategies = ['CL_Label', 'high_level', 'detailed_level']
selected = 0

# Ensure that UMAP and cluster information is present
if 'X_umap' not in adata.obsm:
    raise ValueError("UMAP embeddings not found in the AnnData object.")

# Compute the hierarchical clustering
# Using 'ward' linkage method on the UMAP coordinates
# single complete average weighted centroid median ward
Z = sch.linkage(adata.obsm['X_umap'], method='centroid')

# Ensure that CL_Label is present
if 'CL_Label' not in adata.obs:
    raise ValueError("CL_Label not found in the AnnData object.")

# Create a DataFrame for the clustering results
# No need to create new cluster labels, use existing CL_Label
labels = adata.obs[strategies[selected]].values

# Plot the dendrogram
fig, ax = plt.subplots(figsize=(20, 10))
dendrogram = sch.dendrogram(Z, labels=labels, ax=ax, truncate_mode="level", p=8)

plt.title('Dendrogram of UMAP Clusters')
plt.xlabel('Cell Index')
plt.ylabel('Distance')
plt.show()