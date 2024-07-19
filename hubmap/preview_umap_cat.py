import scanpy as sc
import matplotlib.pyplot as plt
import os
import argparse

high_level = {
    'club cell:nasal': 'Epithelial Cells',
    'multi-ciliated epithelial cell:nasal': 'Epithelial Cells',
    'nasal mucosa goblet cell': 'Epithelial Cells',
    'brush cell of trachebronchial tree': 'Epithelial Cells',
    'tracheobronchial goblet cell': 'Epithelial Cells',
    'multi-ciliated epithelial cell:non-nasal': 'Epithelial Cells',
    'type I pneumocyte': 'Epithelial Cells',
    'type II pneumocyte': 'Epithelial Cells',
    'type II pneumocyte:proliferating': 'Epithelial Cells',
    'alveolar type 1 fibroblast cell': 'Epithelial Cells',
    'alveolar type 2 fibroblast cell': 'Epithelial Cells',
    'serous secreting cell:activated': 'Epithelial Cells',
    'serous secreting cell of bronchus submucosal gland': 'Epithelial Cells',
    'serous secreting cell:nasal': 'Epithelial Cells',
    'airway submucosal gland collecting duct epithelial cell': 'Epithelial Cells',
    'respiratory basal cell': 'Epithelial Cells',
    'respiratory basal cell:resting': 'Epithelial Cells',
    'ionocyte': 'Epithelial Cells',
    'club cell:non-nasal': 'Epithelial Cells',
    'CD1c-positive myeloid dendritic cell': 'Immune Cells',
    'alveolar macrophage': 'Immune Cells',
    'Alveolar Mφ proliferating': 'Immune Cells',
    'monocyte': 'Immune Cells',
    'Non-classical monocytes': 'Immune Cells',
    'Monocyte-derived Mφ': 'Immune Cells',
    'Interstitial Mφ perivascular': 'Immune Cells',
    'natural killer cell': 'Immune Cells',
    'plasmacytoid dendritic cell, human': 'Immune Cells',
    'CD4-positive helper T cell': 'Immune Cells',
    'effector memory CD8-positive, alpha-beta T cell': 'Immune Cells',
    'T cell:proliferating': 'Immune Cells',
    'B cell': 'Immune Cells',
    'plasma cell': 'Immune Cells',
    'mast cell': 'Immune Cells',
    'endothelial cell of venule:pulmonary': 'Vascular Cells',
    'alveolar capillary type 2 endothelial cell': 'Vascular Cells',
    'endothelial cell of venule': 'Vascular Cells',
    'endothelial cell of artery': 'Vascular Cells',
    'capillary endothelial cell': 'Vascular Cells',
    'endothelial cell of lymphatic vessel:differentiating': 'Vascular Cells',
    'endothelial cell of lymphatic vessel:mature': 'Vascular Cells',
    'smooth muscle cell': 'Muscle Cells',
    'lung pericyte': 'Support Cells',
    'myofibroblast cell': 'Support Cells',
    'pulmonary interstitial fibroblast': 'Support Cells',
    'deuterosomal cell': 'Specialized Cells',
    'Transitional Club-AT2': 'Specialized Cells'
}


# Strategy 2: Cell Lineage and Type
detailed_level = {
    'club cell:nasal': 'Nasal Epithelium',
    'multi-ciliated epithelial cell:nasal': 'Nasal Epithelium',
    'nasal mucosa goblet cell': 'Nasal Epithelium',
    'brush cell of trachebronchial tree': 'Tracheobronchial Tree',
    'tracheobronchial goblet cell': 'Tracheobronchial Tree',
    'multi-ciliated epithelial cell:non-nasal': 'Tracheobronchial Tree',
    'type I pneumocyte': 'Alveolar Epithelium',
    'type II pneumocyte': 'Alveolar Epithelium',
    'type II pneumocyte:proliferating': 'Alveolar Epithelium',
    'alveolar type 1 fibroblast cell': 'Alveolar Epithelium',
    'alveolar type 2 fibroblast cell': 'Alveolar Epithelium',
    'serous secreting cell:activated': 'Secretory Cells',
    'serous secreting cell of bronchus submucosal gland': 'Secretory Cells',
    'serous secreting cell:nasal': 'Secretory Cells',
    'airway submucosal gland collecting duct epithelial cell': 'Secretory Cells',
    'respiratory basal cell': 'Basal Cells',
    'respiratory basal cell:resting': 'Basal Cells',
    'ionocyte': 'Specialized Epithelial Cells',
    'club cell:non-nasal': 'Specialized Epithelial Cells',
    'CD1c-positive myeloid dendritic cell': 'Myeloid Lineage',
    'alveolar macrophage': 'Myeloid Lineage',
    'Alveolar Mφ proliferating': 'Myeloid Lineage',
    'monocyte': 'Myeloid Lineage',
    'Non-classical monocytes': 'Myeloid Lineage',
    'Monocyte-derived Mφ': 'Myeloid Lineage',
    'Interstitial Mφ perivascular': 'Myeloid Lineage',
    'natural killer cell': 'Lymphoid Lineage',
    'plasmacytoid dendritic cell, human': 'Lymphoid Lineage',
    'CD4-positive helper T cell': 'Lymphoid Lineage',
    'effector memory CD8-positive, alpha-beta T cell': 'Lymphoid Lineage',
    'T cell:proliferating': 'Lymphoid Lineage',
    'B cell': 'Lymphoid Lineage',
    'plasma cell': 'Lymphoid Lineage',
    'mast cell': 'Mast Cells',
    'myofibroblast cell': 'Fibroblasts and Myofibroblasts',
    'pulmonary interstitial fibroblast': 'Fibroblasts and Myofibroblasts',
    'lung pericyte': 'Pericytes',
    'endothelial cell of venule:pulmonary': 'Vascular Endothelium',
    'endothelial cell of venule': 'Vascular Endothelium',
    'alveolar capillary type 2 endothelial cell': 'Vascular Endothelium',
    'endothelial cell of artery': 'Vascular Endothelium',
    'capillary endothelial cell': 'Vascular Endothelium',
    'endothelial cell of lymphatic vessel:differentiating': 'Lymphatic Endothelium',
    'endothelial cell of lymphatic vessel:mature': 'Lymphatic Endothelium',
    'Transitional Club-AT2': 'Transitional Cells',
    'alveolar type 1 fibroblast cell': 'Alveolar Cells',
    'alveolar type 2 fibroblast cell': 'Alveolar Cells',
    'endothelial cell of lymphatic vessel': 'Miscellaneous'
}


def visualize_umap(h5ad_file_path):
    # Load the AnnData file (with UCE embeddings and cell type labels)
    adata = sc.read_h5ad(h5ad_file_path)
    adata.obs['high_level'] = adata.obs['CL_Label'].map(high_level)
    adata.obs['detailed_level'] = adata.obs['CL_Label'].map(detailed_level)


    # Define the strategies
    strategies = ['CL_Label', 'high_level', 'detailed_level']

    # Check if 'X_umap' is in the obsm slot of the file
    if 'X_umap' in adata.obsm:
        print("UMAP embeddings found. Visualizing...")
        for strategy in strategies:
            print(f"Visualizing based on {strategy}")

            # Check if 'predicted_label' is in the obs slot
            if 'CL_Label' in adata.obs:

                # Set up the figure with a larger size
                fig, ax = plt.subplots(figsize=(20, 20))

                # Plot the UMAP with cell type information
                sc.pl.umap(adata, color=strategy, ax=ax, show=False)

                # Move the legend outside the plot
                #plt.tight_layout(rect=[0, 0, 0.6, 1])  # Adjust the right padding to make room for the legend
                ax.legend()

                # Save the plot to the same folder as the h5ad file
                output_file_path = os.path.join(os.path.dirname(h5ad_file_path), f'umap_{strategy}.png')
                plt.savefig(output_file_path, bbox_inches='tight')

                # Show the plot
                # plt.show()
            else:
                print("Cell type labels not found in the file.")
    else:
        print("UMAP embeddings not found in the file.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize UMAP with UCE embeddings and cell type labels.')
    parser.add_argument('--h5ad_path', type=str, required=True, help='Full path to the h5ad file.')
    args = parser.parse_args()

    visualize_umap(args.h5ad_path)