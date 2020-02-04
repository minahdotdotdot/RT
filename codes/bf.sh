#!/bin/bash

for h in 0.06 0.05 0.025 0.01 0.005 0.0025 0.001 0.0005 0.00025 0.0001
do
	sed -i '3 a h='$h ./setupMMT.jl
	sed -i '5 d' ./setupMMT.jl
	sed -i '1 a #SBATCH --job-name I3'$h ./batchfile.sh
	sed -i '1 d' ./batchfile.sh
	sbatch batchfile.sh
done
