#!/bin/bash -l
#PBS -l walltime=02:00:00,nodes=1:ppn=12,pmem=2580mb
#PBS -m abe
#PBS -N py_jobname
#PBS -M vsethura@umn.edu
#PBS -q  mesabi
cd ${PBS_O_WORKDIR}
echo job_start
mpirun -np 12 ./lmp_mesabi -in in.makemovie -e screen
