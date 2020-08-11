#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=4:00:00
module load R
cd /hpcdata/dgatti/PowerSim/
R CMD BATCH --no-save full.null.harvest.R

