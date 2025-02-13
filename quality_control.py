import os
import scanpy as sc
import matplotlib.pyplot as plt

def run_qc(input_file, output_file):
    # Read the data
    adata = sc.read_10x_mtx(input_file, var_names='gene_symbols', cache=True)

    # Find mitochondrial genes
    adata.var["mito"] = adata.var_names.str.startswith("MT-")

    #  Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars="mito", inplace=True)

    # Create output directory (FIXED: Ensure it exists before saving)
 #   os.makedirs(output_dir, exist_ok=True)

    # Set figure directory
 #   sc.settings.figdir = f"{output_dir}/figures/"

    # Plot the QC metrics
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mito'],
                 jitter=0.4, multi_panel=True, save="_qc_violin_plot.png")
    
    # Filter cells based on QC metrics
    sc.pp.filter_cells(adata, min_genes=200, inplace=True)
    sc.pp.filter_genes(adata, min_cells=3, inplace=True)

    # Save the filtered data
 #   output_file = os.path.join(output_dir, "filtered_data")
    adata.write(output_file)

    print(f"Quality control completed. Filtered data saved to {output_file}.")

if __name__ == '__main__':
 #   input_dir = input("Enter the directory containing 10X data: ").strip()

    # Extract last folder name from input_dir and append "_results"
#    last_folder_name = os.path.basename(os.path.normpath(input_dir))
#    output_dir = f"{last_folder_name}_results"

    run_qc("raw_data", "filtered_data")