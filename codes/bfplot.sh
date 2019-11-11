#!/bin/bash
#SBATCH --job-name B
#SBATCH --qos=blanca-igg
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 72:00:00  
#SBATCH --mail-type=all
#SBATCH --mail-user=luya7574@colorado.edu
cd /projects/luya7574/RT/codes
julia ./plotimages.jl

ffmpeg="/projects/luya7574/ffmpeg/ffmpeg-3.3.4-64bit-static/ffmpeg"
name="B"
$ffmpeg -framerate 6 -i /projects/luya7574/RT/plots/${name}%03d.png -vb 20M ./${name}.avi

