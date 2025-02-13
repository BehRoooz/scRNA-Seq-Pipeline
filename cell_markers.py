import scanpy as sc
import pandas as pd

def run_cell_markers(input_file, output_file):
    # Read the data
    adata = sc.read_h5ad(input_file)

    # Find marker genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

    # Get the marker genes in a dataframe
    markers = sc.get.rank_genes_groups_df(adata, None)

    # Filter the marker genes
    markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > 0.5)]
    
    # Save the marker genes
    markers.to_csv(output_file)

    print(f"Cell markers analysis completed. Marker genes saved to {output_file}.")

if __name__ == "__main__":
    run_cell_markers("clustered_data", "cell_markers.csv")