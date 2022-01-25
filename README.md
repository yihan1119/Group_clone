# Computational identification of clonal cells in single-cell CRISPR screens

## Overview
#### This repository described the scripts for identifying clones in single-cell CRISPR screens with multiplexed sgRNA barcodes. 

![Over-view](./MISC/overview.png "Overview")

## Requirement
* python (3.6 +)
* Numpy (1.16 +)
* Pandas (0.25 +)
* Scipy (1.1 +)
* Matplotlib (3.3 +)

## Original Fastq Files
GEO accession: [GSE185995](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185995)

## Pipelines 
### Step 1: Mapping and pre-processing 
* Map the transcriptome libraries.
* Map the sgRNA libraries. 
* Filter out the doublet cells.
* Generate a pandas dataframe of sgRNA library with singlets only (rows: sgRNA sequence / columns: cell IDs).
* Save as python pickle file.

### Step 2: Group clone analysis 
* Run the group clone script on shell.
* [Group clones](./Scripts/log.group_clones.sh "log.group_clones.sh")

### Step 3: Clonality visualization
* Read the clonal dictionary from jupyter lab and visualize clonality
* [Visualize clonality](./Notebooks/Visualize_clonality-Github.ipynb "Visualize_clonality")

## Contributors 
* First Author: Yihan Wang `Yihan.Wang@UTSouthwestern.edu`
* Corresponding Author: Gary Hon `Gary.Hon@UTSouthwestern.edu`
