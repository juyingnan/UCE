import os
import pandas as pd

def find_cl_labels(root_folder):
    cl_labels = []

    # Walk through all subdirectories and files in the root folder
    for subdir, _, files in os.walk(root_folder):
        for file in files:
            if file.endswith('_annotations.csv'):
                csv_file_path = os.path.join(subdir, file)
                print(f"Processing file: {csv_file_path}")
                
                # Read the CSV file
                df = pd.read_csv(csv_file_path)
                
                # Check if 'CL_label' column exists and extract its values
                if 'CL_Label' in df.columns:
                    labels = df['CL_Label'].dropna().unique()
                    cl_labels.extend(labels)
                else:
                    print(f"'CL_Label' column not found in {csv_file_path}")

    # Remove duplicates by converting to a set and back to a list
    cl_labels = list(set(cl_labels))
    
    return cl_labels

if __name__ == "__main__":
    root_folder = r"C:/Users/yiju/Desktop/temp_data/"
    cl_labels = find_cl_labels(root_folder)
    
    # Print the found CL_labels
    print("Found CL_Labels:")
    for label in cl_labels:
        print(label)
