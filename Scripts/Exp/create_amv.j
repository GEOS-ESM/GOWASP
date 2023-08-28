#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=SIMWINDC
#SBATCH --time=5:59:00 
#SBATCH --nodes=1  --constraint=cas 
#SBATCH --output=SIMWINDC.%j 
#SBATCH --account=s0911 
#SBATCH --partition=preops
#SBATCH --qos=dastest
##SBATCH --qos=debug

# Script for creating satwind data.  Requires xxx min per synoptic time 
# if using using half-hourly sub-windows in each 6-hour synoptic period

# Requires 30min CPU/ day for 0.5 hr NR files.   

# Input are 4 files: One is a version of the a kx-table file that includes 
# a valid table of inflated (to account for QC rejections) obs counts and 
# probability-function parameters distributed by lat, lev, sfctype, kx and 
# satellite (as produced by compute_amv_params.j). Another is the distribution 
# of observation locations that is a function of time, longitude, and N/S 
# hemisphere (as produced by count_amv.j). 
# A third input is an rc file describing the the NR fields and times to be used 
# A fourth contains the BUFR table to be copied.

# Output are 2-3 files: One is a kx-table copied from the input one but 
# with the table of counts replaced by the actual 6-hour mean counts produced 
# when applying theobserving probability function. Amother is a file of wind 
# obs in NCEP .prepbufr format. The third optional output file is a text 
# file of all the obs locations for the first 6-hour period that can be used 
# for testing purposes.

# Specify modules and number of mpi processors
setenv NBRPROCS 16  # must be on a single node
echo "NBRPROCS = ${NBRPROCS}"

  source /home/nprive/GOWASP3_ENV16
  setenv BUFRBIN    $GOWASP_BUFRBIN         # dir for executable "block"
  setenv SIMHOME    $GOWASP_PATH/Sim_satwind    # dir for some files/executables
  setenv ADDTIME    $GOWASP_PATH/addtime.x  # addtime executable
  setenv RMSHMKEY   $GOWASP_RMSHMKEY        # dir for ESMF shmem cleaning script
  setenv SIMWORK    $NOBACKUP/WORK/SimSWwork.$$
#
# RC_SATWIND_TABLES = path to satwind tables                                    
# RC_FILE_FIELDS = file describing NR fields to use                             
# RC_FILE_KX_IN = file describing obs types, desired counts, and prob params
# SATLOC_FILE = INPUT file of obs counts for lons/NS/times bins                 
# RC_FILE_KX_STATS = A file with this as the first part of the name will be 
#   created to hold the actual obs counts created in place of the target counts 
#   appearing in RC_FILE_KX_IN. The second part of the file name will be 
#   the date and hour of the last time processed in this run of the script. 
#
  setenv RC_SATWIND_TABLES $GOWASP_PATH/Rcfiles/Satwind_tables
  setenv RC_FILE_FIELDS    $GOWASP_PATH/Rcfiles/field_list_satwind.rc
  setenv RC_FILE_KX_IN     $RC_SATWIND_TABLES/kx_table_params_10.txt
  setenv RC_FILE_KX_STATS  $RC_SATWIND_TABLES/kx_table_stats_10
#
  setenv BUFR_FILE_IN  /discover/nobackup/rerrico/BUFRDATA/Y2020/M07/gdas1.200701.t00z.prepbufr  # only used to obtain BUFR table
  setenv BUFR_PATH_OUT /discover/nobackup/projects/gmao/nwposse/develop/OSSEobs/N010/SATWIND
  setenv TEXT_FILE_OUT none 
  setenv RAN_SEED 1111
#
#                                                                               
set datetimeN=2006103100  # starting time to process
set ntimes=4             # number of analysis periods to process  

set test_print="F"  # T or F: print sample output for testing purposes          

# END OF USER SET VARIABLES
# ----------------------------------------------------------------------
#
# Make sure files to be written are acessible to others
umask 022

# Create new directory for bufr output if required
# (Check all the times to be processed)
set datetimex=$datetimeN
set add6h="+00"
@ i = 1
while ( $i <= $ntimes ) 
#
  set result=`$ADDTIME $datetimex $add6h`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set datetimex="${yyyy}${mm}${dd}${hh}"
  set add6h="+06"
#
  set bufr_path=$BUFR_PATH_OUT/Y$yyyy/M$mm
  if (! -e $bufr_path) mkdir -p $bufr_path
  @ i++
end     # end loop over time
  setenv RC_FILE_KX_OUT   ${RC_FILE_KX_STATS}_${mm}${dd}${hh}.txt

# Create working directory
if (! -e $SIMWORK            ) mkdir -p $SIMWORK

cd $SIMWORK
/bin/rm -f  *

cp $SIMHOME/create_amv.x prog.x
#
# x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
  set cdtime="${datetimeN}0000"             # 00 min and 00 sec indicated
  $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
  mpirun -np $NBRPROCS ./prog.x $cdtime $ntimes $RC_FILE_FIELDS $RC_FILE_KX_IN $RC_FILE_KX_OUT $RAN_SEED $BUFR_FILE_IN $BUFR_PATH_OUT $TEXT_FILE_OUT $test_print 
    $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
# 
echo "end script"
# --------
  exit 0
# --------
