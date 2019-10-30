#!/bin/bash
#SBATCH --job-name plot0125
#SBATCH --qos=blanca-igg
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 72:00:00  
#SBATCH --mail-type=all
#SBATCH --mail-user=luya7574@colorado.edu
cd /projects/luya7574/RT/codes
julia ./plotimages.jl

ffmpeg="/projects/luya7574/ffmpeg/ffmpeg-3.3.4-64bit-static/ffmpeg"
name="test"
$ffmpeg -framerate 6 -i /projects/luya7574/RT/plots/${name}%04d.png -vb 20M ./Dh0125.avi

