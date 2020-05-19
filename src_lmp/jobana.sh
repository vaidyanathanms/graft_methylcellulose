#!/bin/bash -l
#PBS -l walltime=2:00:00,nodes=1:ppn=1,pmem=2580mb
#PBS -m abe
#PBS -M vsethura@umn.edu
#PBS -q  mesabi
cd ${PBS_O_WORKDIR}
echo job_start
export OMP_NUM_THREADS=1
./ana.o anainp.txt
