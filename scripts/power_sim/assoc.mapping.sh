#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=72:00:00
module load R
cd ~/hpcdata/PowerSim
R CMD BATCH --no-save assoc.mapping.simulations.R

