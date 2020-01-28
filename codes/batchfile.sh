#!/bin/bash
#SBATCH --job-name C:RK3
#SBATCH --qos=blanca-igg
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 72:00:00  
#SBATCH --mail-type=all
#SBATCH --mail-user=luya7574@colorado.edu

#cd /projects/luya7574/RT/codes

julia=/Users/luciaminahyang/julia/julia
matlab=/Applications/MATLAB_R2019a.app/bin/matlab
$julia ./setupMMT.jl
$matlab -nodesktop -nosplash -r "genRationalApprox('A');exit"
$julia ./simplot.jl
