#!/bin/bash
#PBS -o /data/docs/vikas/IT_QMC_H2SQ/run.out
#PBS -e /data/docs/vikas/IT_QMC_H2SQ/run.err
#PBS -l nodes=4:ppn=24


cd  $PBS_O_WORKDIR
cat $PBS_NODEFILE > pbs
#export I_MPI_FABRICS=tcp	

############################## MPI PART #####################################

mpicxx -O3 -fopenmp -ltbb  -O3 -std=c++11 -Wall  -pedantic-errors -g -msse2 -march=native -funroll-loops -ffast-math -o /data/docs/vikas/IT_QMC_H2SQ/main -ggdb /data/docs/vikas/IT_QMC_H2SQ/*.cpp   #$SRC/src/models/SUB_c.f90

time mpirun -np 96 \
/data/docs/vikas/IT_QMC_H2SQ/./main

rm /data/docs/vikas/IT_QMC_H2SQ/main 2>/dev/null
