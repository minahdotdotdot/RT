#!/bin/bash
#SBATCH --job-name dd2D8
#SBATCH --time 7-00:00:00
#SBATCH --ntasks 8
#SBATCH --qos blanca-igg
#SBATCH --partition blanca-igg
#SBATCH --account blanca-igg
#SBATCH --output out.txt
#SBATCH --nodes 1
#SBATCH --mem 24GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luya7574@colorado.edu
#echo "$(hostname --fqdn)"
module load matlab/R2019b
matlab -nodisplay < runsim.m > out.txt 2>&1