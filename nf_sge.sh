#!/bin/bash
#$ -cwd
#$ -j y
#$ -V
#$ -o jianhua
#$ -l virtual_free=10G
#$ -q main.q@DNA
#$ -N impute-test
#$ -pe mpi 2
export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)
mpirun --pernode nextflow run impute_with_onephased_prephasing.nf -with-mpi -c impute_with_onephased_prephasing.config