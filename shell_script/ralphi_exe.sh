#!/bin/bash
#PBS -N ralphi_chr
#PBS -l nodes=1:ppn=48
#PBS -q regular
#PBS -l walltime=24:00:00
#PBS -o ralphi_chr.out
#PBS -e ralphi_chr.err
 
# Directories to various paths
MYHOME=/home/anandita
PROJHOME=${MYHOME}/ralphi
WORKING=/home/anandita/ralphi
DATA=${WORKING}/s_data
 
# Software used in the job
PHASE=${PROJHOME}/engine/phase.py

# Input Files
YAML=${PROJHOME}/n_config/illumina.yaml  # OR ont.yaml for ONT reads

# Output Files

# Logs (including versions)
date

# Initialise environments
source $MYHOME/.bashrc
 
cd $WORKING
# Activate conda environment
conda activate ralphi_env
# Job execution
python ${PHASE} --config ${YAML}  # Edit the YAML file as per your needs

