#!/bin/tcsh

#SBATCH --job-name=filter_umi                     # job name
#SBATCH --partition=super                   # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=60-00:00:00                                # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=Yihan.Wang@utsouthwestern.edu           # specify an email address
#SBATCH --mail-type=end                                   # send email when job status change (start, end, abortion and etc.)

setenv PATH /home2/s426305/.conda/envs/py36/bin:$PATH
echo 'Program is running with the current python version:'
which python
python --version

set DATA_DIR=/project/GCRB/Hon_lab/s426305/Sequencing_data_analysis/10X/YW49-YW54_10x_mapping/nova_seq_data/sgRNA_lib_regex
set SCRIPT_DIR=/project/GCRB/Hon_lab/s426305/Script/sgrna_regex

foreach LIB(\
		    YW61\
		    YW62\
		    YW63\
		    YW64\
		    YW65\
		    YW66\
	    )

set NUM=`grep $LIB $DATA_DIR/lib_num.txt | cut -f2`
echo $LIB
echo $NUM

$SCRIPT_DIR/process_regex_sgrna_df.py\
    -i $DATA_DIR/$LIB.sgRNA_count.txt\
    -n $NUM\
    -o $DATA_DIR/$LIB.sgRNA_adj_df.pkl
end
