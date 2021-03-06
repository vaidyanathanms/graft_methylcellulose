#!/bin/bash -l
#PBS -l walltime=01:00:00,nodes=1:ppn=1,pmem=2580mb
#PBS -m abe
#PBS -N py_jobname
#PBS -M vsethura@umn.edu
#PBS -q  mesabi
cd ${PBS_O_WORKDIR}
echo job_start
export OMP_NUM_THREADS=1
py_objfile py_anafile

wait
mkdir all_output_data
mv rgall_* all_output_data
mv rgavgall_* all_output_data
mv rgMConly_* all_output_data
mv rgsys_* all_output_data
mv mainpersistautocf_* all_output_data
mv indeig_* all_output_data
mv avgeigval_* all_output_data
mv eigMCavg_* all_output_data
mv persist2autocf_* all_output_data
mv compos_* all_output_data
mv distcom_* all_output_data
