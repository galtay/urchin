#!/bin/tcsh
#
# Make a random sample of an Eagle snapshot
#
#BSUB -L /bin/tcsh
#BSUB -n 36
#BSUB -J sample_snapshot
#BSUB -oo sample.out
#BSUB -q cosma
#

module purge
module load python
module load platform_mpi/gnu_4.4.6/8.2.0

setenv INFILE  /gpfs/data/Eagle/mfTestRuns/EAGLE_RUNS/L0100N1440/ANARCHY_CSP_048/COOL_OWLS_LIM_1/AGN_HALO_11p60_SNe_T0_5p5_x2Fe_x2SNIA/data/snipshot_065_z009p266/snip_065_z009p266.11.hdf5
setenv OUTBASE /gpfs/data/Eagle/ttTestRuns/EAGLE_RUNS/L0100N1440/ANARCHY_CSP_048/COOL_OWLS_LIM_1/AGN_HALO_11p60_SNe_T0_5p5_x2Fe_x2SNIA/data/dsnip_002 
setenv sample 0.001

# Run the program
mpirun python ./sample_volume.py $INFILE $OUTBASE $sample 0.001

