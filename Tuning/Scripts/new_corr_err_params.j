#!/bin/tcsh -fx

# Create files of horiz correlation input parameters from param setup files
# No file is created for data type = MASS_

source /home/nprive/GOWASP3_ENV16
source $GOWASP_G5MODS
set CGHOME=$GOWASP_PATH/Tuning/ErrProgs
set VERSION_NUM=04.5
set INPUT_PATH=${GOWASP_PATH}/Rcfiles/ErrParams/V04.5
set OUTPUT_PATH=$INPUT_PATH
#set instr_set=(PREPBUFR AIRS IASI AMSUA ATMS SSMIS MHS WIND)
set instr_set=(PREPBUFR WIND )

set satinfo=$GOWASP_PATH/Rcfiles/sat_info.rc

# END OF USUAL USER-DEFINED VARIABLES

set CGWORK="/discover/nobackup/$user/WORK/CXwork.$$"
# Create working directory and move to it
if (! -e $CGWORK ) mkdir -p $CGWORK
cd $CGWORK
/bin/rm -f  $CGWORK/*

cp $CGHOME/new_corr_err_params.x corr_params.x

foreach dtype ($instr_set)
  set setup="${INPUT_PATH}/hc_setup_${dtype}_${VERSION_NUM}.txt"
  if (-e $setup) then
    set params="${OUTPUT_PATH}/hc_params_${dtype}_${VERSION_NUM}.bin"
    set printout="${OUTPUT_PATH}/hc_printout_${dtype}_${VERSION_NUM}.txt"
    rm -f $params
    ./corr_params.x $setup $params $satinfo > $printout
    if (-e corr_matrix_1) then  # this option for diagnostic purposes only
      /bin/mv corr_matrix_1 ${OUTPUT_PATH}/corr_matrix_1_${dtype}.bin
    endif  
  endif
end   # end loop over instrument types

rm -f corr_params.x

exit



