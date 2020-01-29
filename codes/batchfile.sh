#!/bin/bash
#SBATCH --job-name A:ARK3
#SBATCH --qos=blanca-igg
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 168:00:00  
#SBATCH --mail-type=all
#SBATCH --mail-user=luya7574@colorado.edu

cd /projects/luya7574/RT/codes
#julia ./setupMMT.jl

#module purge
#module load matlab
#matlab -nodesktop -nodisplay -r "addpath(fullfile(cd,'chebfun'));genRationalApprox('A');exit"

julia ./simplot.jl
