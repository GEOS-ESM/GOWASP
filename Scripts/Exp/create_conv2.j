#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=SIMCONV
#SBATCH --time=11:59:00 
#SBATCH --partition=preops
#SBATCH --qos=dastest
#SBATCH --nodes=1 --constraint=cas 
#SBATCH --output=SIMCONV.%j 
#SBATCH --account=s0911 
##SBATCH --qos=debug 

# Script for creating CONV obs (excluding SATWINDS) .  needs 20m/synop time

# Specify modules and number of mpi processors
setenv NBRPROCS 16  # must be on a single node
echo "NBRPROCS = ${NBRPROCS}"

  source ~/GOWASP3_ENV16
  setenv SIMHOME    $GOWASP_PATH/Sim_conv   # dir for some files/executables
  setenv ADDTIME    $GOWASP_PATH/addtime.x  # addtime executable
  setenv RMSHMKEY   $GOWASP_RMSHMKEY # dir for ESMF shmem cleaning scripts
  setenv SIMWORK    $NOBACKUP/WORK/SimCVwork.$$ #
  setenv RC_FILE    $GOWASP_PATH/Rcfiles/field_list_conv.rc # rc file
  setenv BUFR_IN    /discover/nobackup/projects/gmao/nwposse/develop/BUFRDATA  # dir for BUFR data to get bufr tab
  setenv BUFR_OUT   /discover/nobackup/projects/gmao/nwposse/develop/OSSEobs/N010   # dir for BUFR output             
#                                                                               
set datetimeN=2006101100  # starting time for nr field data  
set datetimeB=2020101100  # starting time for real obs data bufr files       
set ntimes =   48           # number of times to process                        

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

cp $SIMHOME/create_conv.x prog.x
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
  @ yy = $yyyy % 10 
  set yy=`printf "%2.2d\n" $yy`
  set datetimeNP="${yy}${mm}${dd}.t${hh}z"      # part of output BUFR file name
  set cdtime_new="${datetimeN}0000"             # 00 min and 00 sec indicated
#
# set time for bufr data file to be input
  set result=`./addtime.x $datetimeB $addh`
  set yyyyB=`echo $result | cut -d" " -f1`
  set yyB=`echo $yyyyB | cut -c 3-4`
  set mmB=`echo $result | cut -d" " -f2`
  set ddB=`echo $result | cut -d" " -f3`
  set hhB=`echo $result | cut -d" " -f4`
  set datetimeB="${yyyyB}${mmB}${ddB}${hhB}"    # time for BUFR input data
  set datetimeBP="${yyB}${mmB}${ddB}.t${hhB}z"  # part of input BUFR file name
  set cdtime_old="${datetimeB}0000"             # 00 min and 00 sec indicated
#
  echo 'Simulation i=' $i 'uses NR time' $datetimeN '& Obs time' $datetimeB
#
# set input and output data paths to proper year and month
  set BF_input_path="$BUFR_IN/Y$yyyyB/M$mmB"
  set BF_output_path="$BUFR_OUT/CONV/Y$yyyy/M$mm"
  if (! -e $BF_output_path ) mkdir -p $BF_output_path
#
  set BF_file_in="$BF_input_path/gdas1.${datetimeBP}.prepbufr"
  set BF_file_out="$BF_output_path/gdas1.${datetimeNP}.prepbufr"
#
  if (-e $BF_file_in ) then  # input data files exist and can be used   
    /bin/rm -f data_out 
    $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
    mpirun -np $NBRPROCS ./prog.x ALL $cdtime_old $cdtime_new $RC_FILE $BF_file_in data_out $test_print 
    $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
    if (-e data_out ) then # output file created and can be saved
      /bin/cp data_out $BF_file_out
    else
      echo 'NO DATA OUTPUT'
      @ ierr_output = $ierr_output + 1
    endif        
  else
    echo 'NONEXISTING INPUT FILE= ' $BF_file_in
    @ ierr_input = $ierr_input + 1
    @ ierr_output = $ierr_output + 1
  endif
#  
  @ i++
end    # loop over times  

echo ' '
echo 'Number of nonexisting input files=' $ierr_input 
echo 'Number of output files not created =' $ierr_output 
#
date
echo "end script"

# --------
  exit 0
# --------
