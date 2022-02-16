# Computational identification of clonal cells in single-cell CRISPR screens

## Overview
#### This repository described the scripts for identifying clones in single-cell CRISPR screens with multiplexed sgRNA barcodes. 

![Over-view](./MISC/overview.png "Overview")

## Requirement
* Python (3.6 +)
* Numpy (1.16 +)
* Pandas (0.25 +)
* Scipy (1.1 +)
* Matplotlib (3.3 +)

## Original Fastq Files
GEO accession: [GSE185995](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185995)

## Pipelines 
### Step 1: Mapping and pre-processing 
* Map the [transcriptome](./Scripts/File_prep/log.10x_map.sh) libraries.
* Map the [sgRNA](./Scripts/File_prep/log_sgRNA_map.sh) libraries and [filter](./Scripts/File_prep/log.filter_umi.sh). 
* Filter out the [experimental doublet cells](./Scripts/File_prep/Filter_HTO_nova.ipynb).
* Generate a pandas dataframe of sgRNA library (rows: sgRNA sequence / columns: cell IDs) and save as [python pickle file](./Scripts/File_prep/generate_sgrna_df_regex.ipynb). 
* Example file can be found in the supplementary file session of [GSE185995](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185995): GSE185995_sgRNA_df_adj_regex.pkl.gz

### Step 2: Group clones analysis 
* Perform the [group clones analysis](./Scripts/Group_clonal_cells.py) on [shell](./Scripts/log.group_clones.sh).

### Step 3: Clonality visualization
* Read the clonal dictionary file with jupyter lab and [visualize clonality](./Notebooks/Visualize_clonality-Github.ipynb).

## How to cite
Wang, Y., Xie, S., Armendariz, D., Hon, GC. [Computational identification of clonal cells in single-cell CRISPR screens] BMC Genomics 23, 135 (2022)
(https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08359-1)

## Contributors 
* First Author: Yihan Wang `Yihan.Wang@UTSouthwestern.edu`
* Corresponding Author: Gary Hon `Gary.Hon@UTSouthwestern.edu`
