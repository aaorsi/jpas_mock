#!/bin/tcsh
#$ -pe parallel 20
#$ -cwd
#$ -S /bin/tcsh
#$ -M aaorsi@cefca.es
#$ -m ea
#$ -N lightcone_emlines
set nproc = 1
#set outfiles = '../out/lines_mags_sfr_mstellar_z'
set outfiles = '../out/test'

./get_all_lines.py $nproc $outfiles

