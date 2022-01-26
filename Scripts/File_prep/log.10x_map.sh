#!/bin/bash

#SBATCH --job-name=10x_HTO
#SBATCH --partition=256GB,256GBv1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=30-00:00:00
#SBATCH --output=job_log.%j.out
#SBATCH --error=job_log.%j.err
#SBATCH --mail-user=Yihan.Wang@utSouthwestern.edu
#SBATCH --mail-type=ALL

STAR=/project/GCRB/Hon_lab/s166631/00.bin/STAR-2.5.2a/bin/Linux_x86_64/STAR
GENOME_REFERENCE_DIR=/project/GCRB/Hon_lab/s160875/02.annotation/Reference/10X_Genomics/hg38_cellranger-3.0.1/GRCh38/

module load cellranger/3.1.0

echo $SLURM_CPUS_ON_NODE
export PATH=$PATH:$CELLRANGER

export libraries="TRANS_LIB"

RESULTS_DIR=$PWD
for library_name in $libraries; do
    echo $library_name

    LIBRARY_DIR=${library_name}_10x_HTO_nova

    cellranger count \
	       --id=$LIBRARY_DIR \
		--libraries=lib_info_$library_name.csv \
		--feature-ref=capture_feature.csv \
		--transcriptome=$GENOME_REFERENCE_DIR \
		--expect-cells=20000 \
		--localcores=$SLURM_CPUS_ON_NODE \
		--chemistry=SC3Pv3
done
