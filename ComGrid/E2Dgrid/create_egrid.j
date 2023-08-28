#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=EGRID
#SBATCH --time=00:05:00 
#SBATCH --nodes=1 --constraint=hasw 
#SBATCH --output=EGRID.%j 
#SBATCH --account=s0911 
#SBATCH --qos=debug 

# Script for putting NR data onto E grid using conservative mapping

# Specify modules and number of mpi processors
setenv NBRNODES 16  # must be on a single node
echo "NBRNODES = ${NBRNODES}"

  source /home/rerrico/GOWASP_3/GOWASP_ENV
  setenv SIMHOME /discover/nobackup/rerrico/GRIDS/E2Dgrid  # dir for executable
  setenv ADDTIME $GOWASP_PATH/addtime.x # addtime executable
  setenv RMSHMKEY $GOWASP_RMSHMKEY # directory for ESMF shmem cleaning scripts
  setenv SIMWORK    $NOBACKUP/WORK/EGRIDwork.$$
  setenv DIR_OUT   $NOBACKUP/WORK/TestE2D
#
set datetimeN=2006072821  # starting time for nr field data
set prog_hours=27
@ ntimes =  2           # number of times to process

set test_print="T"  # T or F: print sample output for testing purposes

# END OF USER SET VARIABLES
# ----------------------------------------------------------------------

date

# Make sure files to be written are acessible to others
umask 022

# Increment in hours between central synoptic input times
set add6h=+24

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
  set cname1=${yyyy}${mm}${dd}_${hh}
  set dir_output_path="$DIR_OUT/Y$yyyy/M$mm"
#
# set prog time based on forecats length
  set result=`./addtime.x $datetimeN $prog_hours`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set datetimeP="${yyyy}${mm}${dd}${hh}"    # NR synoptic time
  set cdatetimeP=${datetimeP}0000
  set cname2=${yyyy}${mm}${dd}_${hh}
#
# set output data path to proper year and month
  if (! -e $dir_output_path ) mkdir -p $dir_output_path
#
  set file_out=$dir_output_path/E2D.prog.eta.${cname1}z+${cname2}z.nc4
#
    /bin/rm -f data_out 
    $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
    mpirun -np $NBRNODES ./prog.x $cdatetime $cdatetimeP data_out $test_print
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
