import scanpy as sc
import matplotlib.pyplot as plt

def run_clustering(input_file, output_file):
    # Read the data
    adata = sc.read_h5ad(input_file)

    # Perform UMAP
    sc.tl.umap(adata)

    # Cluster cells into subgroups using leiden algorithm
    sc.tl.leiden(adata, resolution=0.5) 

    # Plot the UMAP
    sc.pl.umap(adata, color=['leiden'], title="UMAP Leiden Clustering", save="_umap_plot.png")

    # Save the data
    adata.write(output_file)

    print(f"Clustering completed. Data saved to {output_file}.")

if __name__ == '__main__':
    run_clustering("reduced_data", "clustered_data")