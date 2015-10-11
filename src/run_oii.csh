#!/bin/tcsh
#$ -pe openmp 20
#$ -cwd
#$ -S /bin/tcsh

python get_oii.py
