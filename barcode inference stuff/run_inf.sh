#!/bin/bash
# set working directory for this job
#SBATCH --chdir=/srv/gsfs0/projects/snyder/smartis
#
# set the maximum run time, in this example is set for 4 days which is the maximum primary time you can request, you will wait longer on the queue based on the time and the memeo
#SBATCH --time=00-24:00:00
#
# set number of cores
#SBATCH --ntasks-per-node=1
# 
#SBATCH --cpus-per-task=16
#
#This is optional, you can specify on which computer you want it run or you can let SLURM decide based on the availability. For example, here I mentioned the computer  specify to run on nih UV300 prcessor computer named nih_s10

# set the maximum memory usage
#SBATCH --mem=128G
#
# send mail when job ends or aborts
#SBATCH --mail-type=END 
#
# specify an email address
#SBATCH --mail-user smartis@berkeley.edu
#
# Account to charge
#SBATCH --account=mpsnyder

python calculate_total_model_snp_snp_associations.py $1 
