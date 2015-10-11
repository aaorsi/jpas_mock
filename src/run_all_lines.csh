#!/bin/tcsh
#$ -pe parallel 20
#$ -cwd
#$ -S /bin/tcsh
#$ -M aaorsi@cefca.es
#$ -m ea
set nproc = 20
set outfiles = '../out/lines_z'

./get_all_lines.py $nproc $outfiles

