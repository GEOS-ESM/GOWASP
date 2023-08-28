#!/bin/csh -xf

#
# Compute new iterate of error parameters for conventional obs 
# by comparing statistics from real and OSSE exps.
# The only correlations considered here are horizontal.

source ~/GOWASP_3/GOWASP_ENV
setenv PROGHOME $GOWASP_PATH/Tuning/ErrProgs # home for executables

# "target" here means the (generally) real values the osse will try to match
set target_exp=517test
set target_dir=$target_exp
set osse_exp=517osse                
set osse_dir=$osse_exp
set STATS_PATH=/discover/nobackup/rerrico/ODSstats
set path_target_stats=${STATS_PATH}/${target_dir}/Counts_AMV_July
set path_osse_stats=${STATS_PATH}/${osse_dir}/Counts_AMV_July

# set filesin to "hzcor*.bin" if all data types within input directory 
# to be treated the same; otherwise set to, e.g., "hzcor*240*bin". 
# The " " is necessary in the set statement.  
set filesin="count_*.txt"  

# set corr_func_new to the desired shape of the correlation function
# for which primary set of correlation parameters will be determined
set corr_func_new=TOAR

# set location of files of error parameters used to produce the 
# previous OSSE whose statistics are being input here.
set RC_PATH=/discover/home/rerrico/GOWASP_3/Rcfiles
set ERROR_PATH=${RC_PATH}/ErrParams/V03.8
set version=03.8
set error_rc=${ERROR_PATH}/error.V${version}.rc

# set output path
set OUTPUT_PATH=${STATS_PATH}/${osse_dir}/Est_counts_July

#
# END OF USUAL USER-DEFINED VARIABLES
#
if (! -e $OUTPUT_PATH ) mkdir -p $OUTPUT_PATH
set prepbufr_table=prepobs_errtable
set hc_params=hc_params

set DIRWORK=${path_osse_stats}/FitWork
# Create working directory and move to it
if (! -e $DIRWORK ) mkdir -p $DIRWORK
cd $DIRWORK
/bin/rm -f  $DIRWORK/*

cp $PROGHOME/eparams_conv_hcorr.x prog.x
cp ${ERROR_PATH}/${prepbufr_table}*.txt .
cp ${ERROR_PATH}/${hc_params}*.bin .

foreach file1 (../$filesin)

  set cut_name_1 = `echo $file1 | cut -d/ -f 2- `      # remove "../" from begin
  set cut_name_2 = `echo $cut_name_1 | cut -d_ -f 4- ` # remove "hcor_expid_" 
  set instr = `echo $cut_name_2 | cut -d_ -f 1-1 `  
  set field = `echo $cut_name_2 | cut -d_ -f 2-2 `  
  set tail = `echo $cut_name_2 | cut -d_ -f 3- `   # includes region and ".txt"

set file_target_stats=${path_target_stats}/count_${target_exp}_conv_${instr}_${field}_${tail}

set file_osse_stats=${path_osse_stats}/count_${osse_exp}_conv_${instr}_${field}_${tail}

set file_output=${OUTPUT_PATH}/est_${osse_exp}_conv_${instr}_${field}.txt

if (-e $file_target_stats ) then  # then both stats files should exist
./prog.x $instr $field $error_rc $file_target_stats $file_osse_stats $file_output $corr_func_new >> logfile.txt
endif

end      # loop over input files 
#/bin/mv logfile.txt ${OUTPUT_PATH}/.
#/bin/rm -rf $DIRWORK/
exit 0
