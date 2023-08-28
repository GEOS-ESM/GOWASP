#!/bin/tcsh -fx

# Create file of table of lat,lev gpsro error vstand deviations
#
# c_osse: first part of file name for osse stats file (upto the 'region' name)
# c_real: first part of file name for real stats file (upto the 'region' name)
# c_end:  last part of stats files names (after the 'region' name)
# c_names: list of 'region' names
# c_tab_old: name of previous file of error table to be updated (or 'none')
# c_tab_new: name of file of error table to be created

source ~/GOWASP3_ENV16
set PROG=$GOWASP_PATH/Tuning/ErrProgs/gps_vars.x
set c_osse=/discover/nobackup/projects/gmao/nwposse/develop/ODSstats/osse529/Vcorr_June_5d/vzcor_osse529_conv_
set c_real=/discover/nobackup/projects/gmao/nwposse/develop/ODSstats/real529/Vcorr_June_5d/vzcor_real529_conv_ 
set c_end=_omf_lats.txt
set c_names='LN75 LN50 LN30 LN10 LS10 LS30  LS50 LS75'
set c_tab_old=none
set c_tab_new=/discover/nobackup/nprive/GOWASP_3/Rcfiles/ErrParams/V04.4/gps_err_table_04.4.txt
$PROG $c_real $c_osse $c_end $c_tab_old $c_tab_new $c_names


