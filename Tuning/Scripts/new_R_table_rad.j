#!/bin/csh -fx

# Execute new_R_table.f90 to take an old sat_error_table and results from step 
# estimate_Rr.j to create a new table with up-dated R.  It will attempt to 
# update R for every file found in the directory "path_new_R"
#
# If there are multiple files in path_new_R fo each instrument (i.e., for 
# averages over different regions) then change "foreach file1 (../est*.txt)" 
# below to something that indicates which region, such as 
# "foreach file1 (../est*GLOB.txt)" 
#
# NOTE: The est* files are generally produced from software that computes 
# either channel or horizontal correlations, typically produced for 
# different observation types. If the R table is to be updated for all types, 
# then the user here must be sure either to run this code with est* files for 
# all types in the same directory or to run twice, on 2 different directories.
#
source ~/GOWASP3_ENV16
set exp_dir=osse529
set path_new_R=/discover/nobackup/projects/gmao/nwposse/develop/ODSstats/${exp_dir}/Est_rad_June_5d
set file_table_old=$GOWASP_PATH/Rcfiles/ErrParams/V04.4/sat_err_table_V04.4.txt
set file_table_new=$GOWASP_PATH/Rcfiles/ErrParams/V04.5/sat_err_table_V04.5.txt

# END OF USUAL USER-DEFINED VARIABLES

set PROGHOME=$GOWASP_PATH/Tuning/ErrProgs
set DIRWORK=${path_new_R}/NewRWork

# Create working directory and move to it
if (! -e $DIRWORK ) mkdir -p $DIRWORK
cd $DIRWORK
/bin/rm -f $DIRWORK/* 

/bin/cp ${PROGHOME}/new_R_table_rad.f90 prog.f90
/bin/cp $path_new_R/*.txt .
/bin/cp $file_table_old file_old

ifort -o prog.x prog.f90 

foreach file1 (est*.txt)

  ./prog.x $file1 file_old file_new 
  if (-e file_new) then
    /bin/mv file_new file_old # updates old file with successive
  endif

end      # loop over input files 
/bin/mv file_old $file_table_new
/bin/rm -rf $DIRWORK

exit 0

