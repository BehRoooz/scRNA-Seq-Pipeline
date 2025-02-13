import scanpy as sc
import matplotlib.pyplot as plt
import leidenalg

def run_dimensionality_reduction(input_file, output_file):
    # Read the data
    adata = sc.read_h5ad(input_file)

    #Scale each gene to the unit variance of that gene
    sc.pp.scale(adata, max_value=10)

    # Perform PCA
    sc.tl.pca(adata, svd_solver='arpack')

    # Plot the PCA
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save="_pca_variance_ratio.png")

    # Calculate Neighbors
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

    # Save the data
    adata.write(output_file)

    print(f"Dimensionality reduction completed. Data saved to {output_file}.")

if __name__ == '__main__':
    run_dimensionality_reduction("filtered_data", "reduced_data")