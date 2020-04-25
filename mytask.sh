#!/bin/bash

#BSUB -J LUIS[1-2]
#BSUB -P -----
#BSUB -q normal
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -W 4:00
#BSUB -o ./logs/out.%J.%I
#BSUB -e ./logs/err.%J.%I

# =-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=
# Automatic LSF script generated via easypar
# 2020-04-24 17:10:26
# =-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=

# Required modules
module load R/3.6.0 
R_LIBS_USER="-----"

# Input file and R script
file_input=mytask_input_jobarray.csv
R_script=mytask_rscript.R
line=$LSB_JOBINDEX

# Data loading
pnad=$( awk -v line=$line 'BEGIN {FS="\t"}; FNR==line {print $1}' $file_input)
pnac=$( awk -v line=$line 'BEGIN {FS="\t"}; FNR==line {print $2}' $file_input)
pnae=$( awk -v line=$line 'BEGIN {FS="\t"}; FNR==line {print $3}' $file_input)

# Job run
Rscript $R_script $pnad $pnac $pnae
