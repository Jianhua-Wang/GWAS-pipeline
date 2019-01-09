#!/bin/bash
#$ -cwd
#$ -j y
#$ -V
#$ -o jianhua
#$ -q main.q
#$ -N impute-test
nextflow run impute_with_onephased_prephasing.nf -c impute_with_onephased_prephasing.config