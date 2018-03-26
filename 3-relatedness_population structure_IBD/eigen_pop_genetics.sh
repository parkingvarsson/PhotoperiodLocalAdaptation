#! /bin/bash -l

module load bioinfo-tools
module load eigensoft

smartpca -p parfile_smartpca_rm_outlier_ind > smartpca_rm_outlier_ind.log
