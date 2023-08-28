#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=SHMEMTEST
#SBATCH --time=00:05:00 
#SBATCH --nodes=1 --constraint=hasw 
#SBATCH --output=SHMEMTEST.%j 
#SBATCH --account=s0911 
# ------------------------------
setenv NBRPROCS 12  # must be on a single node
source ~/GOWASP_3/GOWASP_ENV
setenv SIMHOME $GOWASP_PATH/Shmem_test   # dir for some files/executables
setenv SIMWORK /discover/nobackup/rerrico/SHMEMTEST.$$
if (! -e $SIMWORK ) mkdir -p $SIMWORK
cd  $SIMWORK

cp $SIMHOME/shmem_reader.x prog.x
mpirun -np $NBRPROCS prog.x
date
exit
