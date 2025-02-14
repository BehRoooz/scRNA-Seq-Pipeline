# scRNA-Seq-Pipeline

## 📌 Overview

This project provides a Nextflow-based pipeline for processing single-cell RNA sequencing (scRNA-Seq) data using Scanpy. It includes quality control, dimensionality reduction, clustering, and cell marker identification in a fully automated workflow.

The pipeline is containerized using Docker, ensuring reproducibility and ease of use.

## 🚀 Features

Automated scRNA-Seq Processing: From raw 10X Genomics data to processed results.

Modular Nextflow Workflow: Easily scalable and customizable.

Containerized Execution: Runs in a Docker container for reproducibility.

Customizable Parameters: Adjust settings as needed.

## 📂 Pipeline Structure

The Nextflow pipeline consists of the following processes:

Quality Control (QC) → Filters cells and genes based on QC metrics.

Dimensionality Reduction → Performs PCA and calculates neighbors.

Clustering → Identifies clusters in the data.

Cell Marker Identification → Determines marker genes for clusters.

## Prerequisites

- Docker (20.10.0 or later)
- Git
- At least 16GB RAM recommended
- Sufficient disk space for your dataset

## 🔧 Installation Steps

1. Clone the repository:
```bash
git clone https://github.com/BehRoooz/scRNA-Seq-Pipeline.git
cd scRNA-Seq-Pipeline
```

2. Build the Docker image:
```bash
docker build -t scrna-seq-pipeline .
```

## Running the Pipeline

1. Prepare your input data:
   - The pipeline expects 10x Genomics format data
   - Input files should include:
     - `matrix.mtx.gz`
     - `features.tsv.gz`
     - `barcodes.tsv.gz`

2. Run the pipeline:
```bash
docker run --rm \
  -v $(pwd):/app \
  -v /path/to/your/data:/app/raw_data \
  -v /path/to/output:/app/results \
  -w /app \
  scrna-seq-pipeline
```

Replace `/path/to/your/data` with the directory containing your input files and `/path/to/output` with your desired output directory.

## 📊 Output

Results will be saved in your specified output directory, including:
- Filtered data
- Quality control metrics
- Dimensionality reduction plots
- Clustering results
- Cell type marker analysis

## 🐛 Troubleshooting

Common issues:

1. Insufficient memory:
   - Increase Docker's memory allocation in Docker Desktop settings

2. File not found errors:
   - Verify that your input files are in the correct format
   - Check that volume mounting paths are correct

3. Permission issues:
   - Ensure write permissions for output directory

## 🤝 Contributing

Feel free to fork this repository, create a new branch, and submit a pull request with improvements!

## Support

For issues and feature requests, please create an issue in the GitHub repository.

📜 License

This project is licensed under the MIT License.
