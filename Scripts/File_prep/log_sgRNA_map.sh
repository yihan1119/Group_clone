#!/bin/tcsh

#SBATCH --job-name=extract_barcode_ID                     # job name
#SBATCH --partition=256GB,256GBv1,384GB                   # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=60-00:00:00                                # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=User@email		                  # specify an email address
#SBATCH --mail-type=end                                   # send email when job status change (start, end, abortion and etc.)

setenv PATH ~/.conda/envs/py36/bin:$PATH

set sgRNA_ref = ./YWsg1.txt
set BARCODE_DIR=/project/GCRB/Hon_lab/s426305/Sequencing_data_analysis/10X/YW49-YW54_10x_mapping/nova_seq_data/10x_pipeline_transcriptome/10x_lib_10x_HTO_nova/outs/filtered_feature_bc_matrix
set FASTQ_DIR=/project/GCRB/Hon_lab/s426305/Sequencing_data_analysis/10X/YW49-YW54_10x_mapping/nova_seq_data/sgRNA_lib_regex/SGRNA.tmp_fastq_combined

zcat $BARCODE_DIR/barcodes.tsv.gz > 10x_lib.barcodes.tsv

./extract_sgRNA_from_reads.py\
    -i $FASTQ_DIR/SGRNA_R1.fastq.gz $FASTQ_DIR/SGRNA_R2.fastq.gz\
    -b 10x_lib.barcodes.tsv\
    -r $sgRNA_ref\
    -t $SLURM_CPUS_ON_NODE\
    -m 1\
    -o SGRNA.sgRNA_cell-bc.summary.txt

	
./summarize_sgRNA.ver2.py\
    -i SGRNA.sgRNA_cell-bc.summary.txt\
    -o SGRNA.sgRNA_count.txt


