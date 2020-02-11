#!/bin/bash

#for h in 0.06 #0.05 0.025 0.01 0.005 0.0025 0.001 0.0005 0.00025 0.0001
for h in 0.05 #0.025 0.01 0.005 0.0025 #0.001 0.0005 0.00025 0.0001

do
	sed -i '3 a h='$h ./setupMMT.jl
	sed -i '5 d' ./setupMMT.jl
	sed -i '1 a #SBATCH --job-name I3R'$h ./batchfile.sh
	sed -i '3 d' ./batchfile.sh
	cm='"addpath(fullfile(cd,'\''chebfun'\''));genRationalApprox('\''IFRK3Rh='
        cm2='d6'\'')exit"'
	sed -i '14 a matlab -nodesktop -nodisplay -r '${cm}$h${cm2} ./batchfile.sh 
	#sed -i '14 a matlab -nodesktop -nodisplay -r "addpath(fullfile(cd,'chebfun'));genRationalApprox('IFRK3Rh='${h}'d6');exit"' ./batchfile.sh
	sed -i '16 d' ./batchfile.sh
	#sbatch batchfile.sh
	#sleep 45
done
