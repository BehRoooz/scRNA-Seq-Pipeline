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

## 🔧 Installation

### 1️⃣ Install Docker

Ensure Docker is installed on your system:

sudo apt-get update && sudo apt-get install -y docker.io

### 2️⃣ Install Nextflow
Install Nextflow globally:

curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/

### 3️⃣ Clone the Repository

git clone https://github.com/BehRoooz/scRNA-Seq-Pipeline.git
cd scRNA-Seq-Pipeline

### 4️⃣ Build the Docker Image

docker build -t scrna-seq-pipeline .

## ▶️ Running the Pipeline

Run the pipeline using Docker:

docker run --rm -v $(pwd):/app -w /app scrna-seq-pipeline

Alternatively, if you have Nextflow installed:

nextflow run main.nf -with-docker scrna-seq-pipeline

## ⚙️ Customizing Parameters

Modify params in main.nf to adjust input/output directories:

params.input = "raw_data"
params.output = "results"

## 📊 Output

The results will be stored in the results/ directory:

filtered_data.h5ad → QC-filtered data

reduced_data.h5ad → PCA-reduced data

clustered_data.h5ad → Clustered data

cell_markers.csv → Identified cell markers

## 🐛 Troubleshooting

Check the Nextflow Logs

If a process fails, inspect Nextflow logs:

cat .nextflow.log

Debug with Interactive Mode

Enter the container:

docker run -it --rm -v $(pwd):/app -w /app scrna-seq-pipeline /bin/bash

Run Nextflow manually:

nextflow run main.nf

## 🤝 Contributing

Feel free to fork this repository, create a new branch, and submit a pull request with improvements!

📜 License

This project is licensed under the MIT License.