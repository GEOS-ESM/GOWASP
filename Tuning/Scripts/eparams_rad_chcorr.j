#!/bin/csh -xf                                                                                

# Compute new iterate of error parameters for rad obs 
# by comparing statistics from real and OSSe exps.
# The only correlations considered here are inter-channel 

source /home/nprive/GOWASP3_ENV16
setenv PROGHOME $GOWASP_PATH/Tuning/ErrProgs # home for executables
set sat_info_file=$GOWASP_PATH/Rcfiles/sat_info.rc

# "target" here means the (generally) real values the osse will try to match
set target_exp=real529
set target_dir=$target_exp
set osse_exp=osse529          
set osse_dir=$osse_exp
set STATS_PATH=/discover/nobackup/projects/gmao/nwposse/develop/ODSstats
set path_target_stats=${STATS_PATH}/${target_dir}/Chcorr_June_5d
set path_osse_stats=${STATS_PATH}/${osse_dir}/Chcorr_June_5d

# set filesin to "chcor*.bin" if all data types within input directory 
# to be treated the same; otherwise set to, e.g., "chcor*WHEM.bin" for 
# only considering files for the region 'WHEM'.
# The " " is necessary in the set statement when * is in the generic name  
set filesin="chcor*.bin"  

# set location of files of error parameters used to produce the 
# previous OSSE whose statistics are being input here.
set RC_PATH=/discover/nobackup/nprive/GOWASP_3/Rcfiles
set ERROR_PATH=${RC_PATH}/ErrParams/V04.5
set version=04.5
set error_rc=${ERROR_PATH}/error.V${version}.rc

# set output path
set OUTPUT_PATH=${STATS_PATH}/${osse_dir}/Est_rad_June_5d

#
# END OF USUAL USER-DEFINED VARIABLES
#
if (! -e $OUTPUT_PATH ) mkdir -p $OUTPUT_PATH
set sat_table=sat_err_table
set hc_params=hc_params
set DIRWORK=${path_osse_stats}/FitWork

# Create working directory and move to it
if (! -e $DIRWORK ) mkdir -p $DIRWORK
cd $DIRWORK
/bin/rm -f  $DIRWORK/*

cp $PROGHOME/eparams_rad_chcorr.x prog.x
cp ${ERROR_PATH}/${sat_table}*.txt .
cp ${ERROR_PATH}/${hc_params}*.bin .

foreach file1 (../$filesin)

  set cut_name_1 = `echo $file1 | cut -d/ -f 2- `      # remove "../" from begin
  set cut_name_2 = `echo $cut_name_1 | cut -d_ -f 3- ` # remove "hcor_expid_" 
  set instr = `echo $cut_name_2 | cut -d_ -f 1-1 `  
  set sat = `echo $cut_name_2 | cut -d_ -f 2-2 `  
  set tail = `echo $cut_name_2 | cut -d_ -f 3- `   # includes region and ".txt"

set file_target_stats=${path_target_stats}/chcor_${target_exp}_${instr}_${sat}_${tail}

set file_osse_stats=${path_osse_stats}/chcor_${osse_exp}_${instr}_${sat}_${tail}

set file_output=${OUTPUT_PATH}/cov_${osse_exp}_${instr}_${sat}.txt
set file_cov_new=${OUTPUT_PATH}/cov_${osse_exp}_${instr}_${sat}.bin

echo $file_target_stats
if (-e $file_target_stats) then  # then both stats files should exist
  ./prog.x $instr $sat $error_rc $file_target_stats $file_osse_stats $sat_info_file $file_output $file_cov_new >> logfile

  if (-e corr_matrix_3) then  # this option for diagnostic purposes only
    /bin/mv corr_matrix_3 ${OUTPUT_PATH}/corr_matrix_3_${instr}_${sat}.bin
  endif
endif
end      # loop over input files 
#/bin/mv logfile ${OUTPUT_PATH}/logfile_chcorr.txt 
#/bin/rm -rf $DIRWORK/
exit 0
