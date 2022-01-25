#!/bin/tcsh

#SBATCH --job-name=group_clones                           # job name
#SBATCH --partition=128GB,256GB,256GBv1                   # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=30-00:00:00                                # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=user@email                            # specify an email address
#SBATCH --mail-type=end                                   # send email when job status change (start, end, abortion and etc.)

echo 'Hello World'

#set python environment
setenv PATH ~/.conda/envs/py36/bin:$PATH
echo 'Program is running with the current python version:'
which python
python --version

set sgRNA_df=
set CUTOFF=0.05

./Group_clonal_cells.py\
    -s $sgRNA_df\
    -c $CUTOFF\
    -o ./Clone_tree.$CUTOFF.txt \
    -op ./Clone_pval.$CUTOFF.txt
