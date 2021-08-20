#!/bin/bash
#PBS -N test
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -l select=1:ncpus=24:mem=1gb
#PBS -M of15641@bristol.ac.uk 

# Move to the directory the job was executed in
cd $PBS_O_WORKDIR
pwd
module load apps/gaussian/16

#Execute gaussian
g16 < test.com > test.log

