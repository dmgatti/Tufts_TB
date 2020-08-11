#!/bin/bash -l
#PBS -l nodes=1:ppn=1,mem=64GB,walltime=8:00:00
cd /hpcdata/dgatti/Tufts/scripts
module load R/3.3.1
R CMD BATCH --no-save tufts_qtl.R
