#!/bin/csh -fx

# Execute new_R_table_conv.f90 to take an old prepobs tables and results 
# est_* files to create a new table with up-dated R values.  Udates will be
# made for all kx and error-estimate files specified in a resource file
source ~/GOWASP3_ENV16
set exp_dir=ctl529
set path_new_R=/discover/nobackup/projects/gmao/nwposse/develop/ODSstats/${exp_dir}/Est_conv_June_5d
set file_old_R=$GOWASP_PATH/Rcfiles/ErrParams/V04.4/prepobs_errtable.V04.4.txt
set file_new_R=$GOWASP_PATH/Rcfiles/ErrParams/V04.6/prepobs_errtable.V04.6.txt

# END OF USUAL USER-DEFINED VARIABLES

set PROGHOME=$GOWASP_PATH/Tuning/ErrProgs
set DIRWORK=${path_new_R}/NewRWork
set rc_file=$GOWASP_PATH/Tuning/Scripts/new_R_table_conv.rc

# Create working directory and move to it
if (! -e $DIRWORK ) mkdir -p $DIRWORK
cd $DIRWORK
/bin/rm -f $DIRWORK/* 


/bin/cp $PROGHOME/new_R_table_conv.f90 prog.f90
/bin/cp $path_new_R/*.txt .
ifort -o prog.x prog.f90 

./prog.x $rc_file $file_old_R $file_new_R

/bin/rm -rf $DIRWORK

exit 0

