#!/bin/bash -l
#PBS -l walltime=12:00:00,nodes=1:ppn=6,pmem=2580mb
#PBS -m abe
#PBS -N py_jobname
#PBS -M vsethura@umn.edu
#PBS -q  mesabi
cd ${PBS_O_WORKDIR}
echo job_start
mpirun -np 6 ./lmp_mesabi -in in.umbrella -e screen
