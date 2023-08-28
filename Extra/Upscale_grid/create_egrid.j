#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=EGRID
#SBATCH --time=00:05:00 
#SBATCH --nodes=1 --constraint=hasw 
#SBATCH --output=EGRID.%j 
#SBATCH --account=s0911 

# Script for putting NR data onto E grid using conservative mapping

# Specify modules and number of mpi processors
setenv NBRNODES 16  # must be on a single node
echo "NBRNODES = ${NBRNODES}"

  source ~/GOSS4_ENV
  setenv SIMHOME    /discover/home/rerrico/OSSE3/Egrid  # dir for some files/executables
  setenv ADDTIME    /discover/home/rerrico/OSSE3/addtime.x  # addtime executable
  setenv RMSHMKEY   /discover/home/rerrico/bin # directory for ESMF shmem cleaning scripts
  setenv SIMWORK    $NOBACKUP/WORK/EGRIDwork.$$
  setenv DIR_OUT   $NOBACKUP/EGRIDfiles    # dir for output
#
set datetimeN=2006072900  # starting time for nr field data
@ ntimes =  2           # number of times to process

set test_print="F"  # T or F: print sample output for testing purposes

# END OF USER SET VARIABLES
# ----------------------------------------------------------------------

date

# Make sure files to be written are acessible to others
umask 022

# Increment in hours between central synoptic input times
set add6h=+06

# Create working directory
if (! -e $SIMWORK            ) mkdir -p $SIMWORK

cd $SIMWORK
/bin/rm -f  *

cp $SIMHOME/NR2Egrid.x prog.x
cp $ADDTIME addtime.x

# x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
# loop over synoptic times

@ ierr_input  = 0  # initialize error counter
@ ierr_output = 0  # initialize error counter

@ i = 1
@ ilast =  $ntimes  
while ( $i <= $ntimes ) 
  
  if ($i == 1) then
    set addh=000
  else
    set addh=$add6h 
  endif
#

# set time for data to be created (must be same as NR time)
  set result=`./addtime.x $datetimeN $addh`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set datetimeN="${yyyy}${mm}${dd}${hh}"    # NR synoptic time
  set cdatetime=${datetimeN}0000
#
# set output data path to proper year and month
  set dir_output_path="$DIR_OUT/Y$yyyy/M$mm"
  if (! -e $dir_output_path ) mkdir -p $dir_output_path
#
  set file_out="$dir_output_path/egridNR.${datetimeN}.nc4"
#
    /bin/rm -f data_out 
    $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
    mpirun -np $NBRNODES ./prog.x $cdatetime data_out $test_print
    $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
    if (-e data_out ) then # output file created and can be saved
      /bin/cp data_out $file_out
    else
      echo 'NO DATA OUTPUT'
      @ ierr_output = $ierr_output + 1
    endif        
#  
  @ i++
end    # loop over times  

echo ' '
echo 'Number of output files not created =' $ierr_output 
#
echo "end script"
date

# --------
  exit 0
# --------
