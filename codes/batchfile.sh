#!/bin/bash
#SBATCH --job-name 8_0.05
#SBATCH --qos=blanca-igg
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 24:00:00  
#SBATCH --mail-type=all
#SBATCH --mail-user=luya7574@colorado.edu

cd /projects/luya7574/RT/codes
#julia ./setupMMT.jl

#module purge
#module load matlab
#matlab -nodesktop -nodisplay -r "addpath(fullfile(cd,'chebfun'));genRationalApprox('IFRK3Rh=0.025d6');exit"

julia ./simplot.jl



