#!/bin/tcsh
#$ -pe parallel 100
#$ -cwd
#$ -S /bin/tcsh
#$ -M aaorsi@cefca.es
#$ -m ea
set nproc = 100
set outfiles = '../out/lines_z'

./get_all_lines.py $nproc $outfiles

