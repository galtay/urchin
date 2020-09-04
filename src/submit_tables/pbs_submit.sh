#!/bin/tcsh 
# 
#PBS -S /bin/tcsh
#PBS -N slab_tables
#PBS -l nodes=1:ppn=12 
#PBS -j oe 
#PBS -t 1-41


set EXEDIR = /gpfs/data/galtay/google-svn/sea-urchin/src 
set EXE = ${EXEDIR}/make_slab_lookup_table

setenv OMP_NUM_THREADS 24 

set REDSHIFT = (0.0, 0.25, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0, 2.25, 2.50, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.50, 6.75, 7.0, 7.25, 7.50, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.50, 9.75, 10.0)

cd ${EXEDIR} 
${EXE} 6.0 1000000 $REDSHIFT[${PBS_ARRAYID}] 4.0d0 128 64 32 8 >& ${EXEDIR}/submit_tables/screen.6.1000000-${PBS_ARRAYID}".out" 
