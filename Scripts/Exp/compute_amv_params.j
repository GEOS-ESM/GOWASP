#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=SIMWINDP
#SBATCH --time=2:40:00 
#SBATCH --nodes=1 
#SBATCH --output=SIMWINDP.%j 
#SBATCH --partition=preops
#SBATCH --account=s0911
#SBATCH --qos=dastest
##SBATCH --qos=debug

# Script for computing a file of parameters to be used to define the 
# probability functions that will be used to determine if an observation 
# is present at each possible (thinned) viewing location. 

# Requires 30min CPU/ day for 0.5 hr NR files (but if temporal interpolation
# is required, it can double the time).  

# Input are 3 files: One is a version of the a kx-table file that 
# includes a valid table of obs counts distributed by lat, lev, sfctype, 
# kx and satellite. Another is the distribution of observation locations 
# that is a function of time, longitude, and N/S hemisphere.  These two 
# inputs are produced by count_amv.j. A third input is an rc file describing 
# the the NR fields and times to be used. 
#
# Output is a copied and updated version of the input kx-table file.           
# It now contains the table of parameters used to define the observation 
# probability function. Also the desired mean counts are replaced by the 
# modified values that have been inflated to account for a percentage of 
# later GSI QC rejections.
#
# Specify modules and number of mpi processors
setenv NBRPROCS 14  # must be on a single node
echo "NBRPROCS = ${NBRPROCS}"

  source ~/GOWASP3_ENV16
  setenv BUFRBIN    $GOWASP_BUFRBIN          # dir for executable "block"
  setenv SIMHOME    $GOWASP_PATH/Sim_satwind # dir for some files/executables
  setenv ADDTIME    $GOWASP_PATH/addtime.x   # addtime executable
  setenv RMSHMKEY   $GOWASP_RMSHMKEY         # dir for ESMF a shmem script
  setenv SIMWORK    $NOBACKUP/WORK/SimSWwork.$$
#
# RC_SATWIND_TABLES = path to satwind tables  
# RC_FILE_FIELDS = file describing NR fields to use 
# RC_FILE_KX_IN = file describing obs types and how they are to be used   
# RC_FILE_KX_OUT = copy RC_FILE_KX_IN and add table of prob params
# HISTOGRAM_FILE = Optional diagnostic OUTPUT file of field histograms (or set to "none")
  setenv RC_SATWIND_TABLES $GOWASP_PATH/Rcfiles/Satwind_tables 
  setenv RC_FILE_FIELDS    $GOWASP_PATH/Rcfiles/field_list_satwind.rc
  setenv RC_FILE_KX_IN     $RC_SATWIND_TABLES/kx_table_counts_10.txt
  setenv RC_FILE_KX_OUT    $RC_SATWIND_TABLES/kx_table_params_10.txt
  setenv HISTOGRAM_FILE    $RC_SATWIND_TABLES/f_hist_10.txt
#  
set test_print="F"  # T or F: print sample output for testing purposes          

# END OF USER SET VARIABLES
# ----------------------------------------------------------------------

# Make sure files to be written are acessible to others
umask 022

# Create working directory
if (! -e $SIMWORK            ) mkdir -p $SIMWORK

cd $SIMWORK
/bin/rm -f  *

cp $SIMHOME/compute_params.x prog.x
cp $ADDTIME addtime.x
#
# x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
  $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
  mpirun -np $NBRPROCS ./prog.x $RC_FILE_FIELDS $RC_FILE_KX_IN $RC_FILE_KX_OUT $HISTOGRAM_FILE  1111 $test_print 
    $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
# 
echo "end script"
# --------
  exit 0
# --------
