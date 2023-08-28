#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=SIMRADM
#SBATCH --time=02:59:00 
#SBATCH --nodes=1   --constraint=cas 
#SBATCH --output=SIMRADM.%j 
#SBATCH --account=s0911 
##SBATCH --qos=debug

# Script for creating binary files of possibly thinned header information 
# drawn from obs BUFR files for ingestion into create_rad_profs 
  source //home/nprive/GOWASP3_ENV16
  setenv SIMHOME    $GOWASP_PATH/Sim_rad    # dir for some files/executables
  setenv ADDTIME    $GOWASP_PATH/addtime.x  # addtime executable
  setenv SIMWORK    $NOBACKUP/WORK/SimRMwork.$$
#
  setenv PROF_IN1   /discover/nobackup/projects/gmao/nwposse/develop/OSSEobs/N010/ObsProf  # dir for profile input
  setenv PROF_IN2   /discover/nobackup/projects/gmao/nwposse/develop/OSSEobs/N010/ObsProf_QR  # dir for profile input
  setenv PROF_OUT   /discover/nobackup/projects/gmao/nwposse/develop/OSSEobs/N010/ObsProf  # dir for reordered output
#
set datetimeN=2006100100  # starting time to process
@ ntimes =   124            # number of times to process

# set which radiance data sets to create
set instr_set= (GMI AMSR2 MHS)

set test_print="F"  # T or F: print sample output for testing purposes

# END OF USER SET VARIABLES
# ----------------------------------------------------------------------
# Make sure files to be written are acessible to others
umask 022

# Increment in hours between central synoptic input times
set add6h=+06

# Create working directory
if (! -e $SIMWORK            ) mkdir -p $SIMWORK

cd $SIMWORK
/bin/rm -f  *

cp $SIMHOME/create_rad_merge.x prog.x
cp $ADDTIME addtime.x

# x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
# loop over synoptic times

@ i = 1
@ ilast =  $ntimes  
while ( $i <= $ntimes ) 
  
  if ($i == 1) then
    set addh=000
  else
    set addh=$add6h 
  endif
#
# set time for data to be created 
  set result=`./addtime.x $datetimeN $addh`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set datetimeN="${yyyy}${mm}${dd}${hh}"    # current synoptic time considered
  @ yy = $yyyy % 10 
  set yy=`printf "%2.2d\n" $yy`
  set datetimeNP="${yy}${mm}${dd}.t${hh}z"  # part of output BUFR file name
#
# Loop over instruments types to consider
  foreach instrument ($instr_set)
#
# Set intrument-dependent file name information
#
# Set name of input directory and file for this instrument
    set prof_dir_in1=$PROF_IN1/$instrument/Y$yyyy/M$mm
    set prof_dir_in2=$PROF_IN2/$instrument/Y$yyyy/M$mm
    set prof_dir_out=$PROF_OUT/$instrument/Y$yyyy/M$mm
    set prof_file_in1=${prof_dir_in1}/profile_${instrument}.${datetimeN}.bin 
    set prof_file_in2=${prof_dir_in2}/profile_${instrument}.${datetimeN}.bin 
    set prof_file_out=${prof_dir_out}/profile_${instrument}.${datetimeN}.bin 

    if (! -e $prof_dir_out ) mkdir -p $prof_dir_out    # create output DIR 

    if (-e $prof_file_in1 && -e $prof_file_in2 ) then  # input data files exist and can be used   
      /bin/rm data_out  
      ./prog.x $prof_file_in1 $prof_file_in2 data_out $test_print   
       if (-e data_out ) then # output file created and can be saved
         /bin/cp data_out $prof_file_out
       else
         echo 'NO DATA OUTPUT'
       endif        
     else
       echo 'NONEXISTING INPUT FILE= ' $prof_file_in1
       echo '         OR INPUT FILE= ' $prof_file_in2
       echo 'Skip processing this file'
     endif
#
  end  # loop over instrument
#
  @ i++
end    # loop over times  
# All done
echo "end script"
date

# --------
  exit 0
# --------
