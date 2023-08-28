#!/bin/csh -xf 
#SBATCH -A s0911
#PBS -l select=1:ncpus=2   
#PBS -l walltime=00:01:00 
#PBS -o SIMCONV.%j

# Script for creating CONV obs (excluding SATWINDS) .  needs 20m/synop time

# Specify modules and number of mpi processors
setenv NBRNODES 2  # must be on a single node
echo "NBRNODES = ${NBRNODES}"

  source ~/GOSS2_ENV
  source $GOSS_G5MODS
  setenv SIMHOME    $GOSS_PATH/Error/TEST  # dir for some files/executables
  setenv SIMWORK    $NOBACKUP/SimERRwork.$$
#

cp $SIMHOME/testMPI.x prog.x

    mpirun -np $NBRNODES ./prog.x 
# --------
  exit 0
# --------
