import os
import scanpy as sc
import pandas as pd

# Define the root folder
root_folder = r"C:/Users/yiju/Desktop/temp_data/"

# Loop over each subfolder in the root folder
for subfolder in os.listdir(root_folder):
    subfolder_path = os.path.join(root_folder, subfolder)
    
    if os.path.isdir(subfolder_path):
        # Define file paths
        h5ad_file = os.path.join(subfolder_path, f"{subfolder}_data.h5ad")
        csv_file = os.path.join(subfolder_path, f"{subfolder}_annotations.csv")
        
        # Load AnnData
        adata = sc.read(h5ad_file)
        
        # Replace NaN values in the 'hugo_symbol' column with a placeholder
        placeholder = "Unknown"
        
        # Ensure 'hugo_symbol' is treated as a categorical column
        adata.var['hugo_symbol'] = adata.var['hugo_symbol'].astype('category')
        
        # Add the placeholder as a new category
        adata.var['hugo_symbol'].cat.add_categories([placeholder], inplace=True)
        
        # Replace NaN values with the placeholder
        adata.var['hugo_symbol'].fillna(placeholder, inplace=True)
        
        # Load the CSV file
        annotations = pd.read_csv(csv_file)
        
        # Check if 'CL_Label' column exists in the CSV file
        if 'CL_Label' in annotations.columns:
            # Map the 'CL_Label' to the corresponding observations in AnnData
            adata.obs['CL_Label'] = annotations.set_index(adata.obs_names)['CL_Label']
        
        # Save the adjusted AnnData with _fixed added to the filename
        fixed_h5ad_file = os.path.join(subfolder_path, f"{subfolder}_data_fixed.h5ad")
        adata.write(fixed_h5ad_file)

print("Processing complete.")
