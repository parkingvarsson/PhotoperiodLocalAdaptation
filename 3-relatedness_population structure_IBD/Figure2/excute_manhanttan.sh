#! /bin/bash -l

#SBATCH -A b2011141
#SBATCH -p core
#SBATCH -o manhanttan.plot.out
#SBATCH -e manhanttan.plot.err
#SBATCH -J manhanttan.plot.job
#SBATCH -t 12:00:00

Rscript manhanttan.3methods.R


