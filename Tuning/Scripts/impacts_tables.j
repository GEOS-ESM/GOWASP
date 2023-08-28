#!/bin/csh -xf
# ------------------------------
# Execute countobs_sum.x and countobs_tables.x to create files of tables 
# of means, standard deviations, and counts of O-F for all radiance obs 
# for a selected region. Also, compute a distribution with latitude.
#
# If new instrument types are added, the number of channels must 
# be specified at the appropriate location in the script.
#

source ~/GOWASP_3/GOWASP_ENV

setenv PROGHOME $GOWASP_PATH/Tuning/StatsProgs  # home for executables

set rcfile=$GOWASP_PATH/Tuning/Scripts/impacts.rc
set outdir=$NOBACKUP/ODSstats/G514osse/Impact0
set outfile=impacts_table.txt

# END OF USUAL USER-DEFINED VARIABLES

setenv WORKDIR $NOBACKUP/WORK/IMPwork.$$ 
# Create and move to empty working directory
if (! -e $WORKDIR            ) mkdir -p $WORKDIR
cd $WORKDIR
/bin/rm -f  *

# Copy file to working directory
cp $PROGHOME/impacts_tables.x prog.x

./prog.x $rcfile > logfile
  /bin/mv logfile ${outdir}/impacts_table_printout
if (-e file_table_output) then 
  if (! -e $outdir ) mkdir -p $outdir
  cat file_table_output $rcfile > file_cat
  /bin/mv file_cat ${outdir}/${outfile}
endif

# --------
exit 0
# --------
