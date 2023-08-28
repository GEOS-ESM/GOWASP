#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=SIMWINDN
#SBATCH --time=05:59:00 
#SBATCH --nodes=1 
#SBATCH --output=SIMWINDN.%j 
#SBATCH --account=s0911 
##SBATCH --qos=debug

# Script for counting satwind data in ODS files and creating a table of 
# viewing locations.
#  
# Requires 20s CPU per day
#
# Input are 3 files: One is a kx_table but only the initial set of variables 
# and the kx_descriptor table are used. Another set of input files are the 
# ODS files that provide the data to be counted and located. A third input is an 
# rc file describing the the NR fields and times to be used (In particular, the 
# frquency of obs times to consider that is the same as the frquency to be used for
# the NR data in the other programs to be run; otherwise the only field used 
# here is the one that provides ocean fraction. 
# 
# Output are 2 files: One a copied and updated version of the input kx-table file. 
# It now contains a distribution of mean obs counts pulled from the ODS files. 
# The other is a file of obs location data that is a function of time, longitude, 
# and N/S hemisphere.  
#
  source /home/nprive/GOWASP3_ENV16
  setenv SIMHOME    $GOWASP_PATH/Sim_satwind  # dir for some files/executables
  setenv ADDTIME    $GOWASP_PATH/addtime.x  # addtime executable
  setenv RMSHMKEY   $GOWASP_RMSHMKEY        # dir for ESMF shmem cleaning script
  setenv SIMWORK    $NOBACKUP/WORK/SimSWwork.$$
#
# RC_SATWIND_TABLES = path to satwind tables 
# RC_FILE_FIELDS = file describing NR fields to use
# RC_FILE_KX_IN = file describing obs types and how they are to be used
# RC_FILE_KX_OUT = copy RC_FILE_KX_IN and add table of obs counts
  setenv RC_SATWIND_TABLES $GOWASP_PATH/Rcfiles/Satwind_tables 
  setenv RC_FILE_FIELDS    $GOWASP_PATH/Rcfiles/field_list_satwind.rc
  setenv RC_FILE_KX_IN     $RC_SATWIND_TABLES/kx_table_in_10.txt     
  setenv RC_FILE_KX_OUT    $RC_SATWIND_TABLES/kx_table_counts_10.txt 
  setenv DXMIN_FILE        none
  setenv ODS_TIME_OFFSET   0  # assumes times on file are relative to center time 

# Describe ODS files to input
  setenv ODS_PATH /discover/nobackup/projects/gmao/nwposse/develop/real529/obs/
  setenv ODS_TEMPLATE 'Y$yyyy#/M$mm#/D$dd#/H$hh#/real529.diag_conv.$yyyymmdd_hh#z.ods'

set datetimeNR=2006062200   # This can be any valid NR time since only used to get FROCEAN

# END OF USER SET VARIABLES
# ----------------------------------------------------------------------

# Make sure files to be written are acessible to others
umask 022

# Create working directory
if (! -e $SIMWORK            ) mkdir -p $SIMWORK

cd $SIMWORK
/bin/rm -f  *

cp $SIMHOME/count_amv.x prog.x
cp $ADDTIME addtime.x
#
# x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
  ./prog.x $datetimeNR  $RC_FILE_FIELDS $ODS_PATH $ODS_TEMPLATE $RC_FILE_KX_IN $RC_FILE_KX_OUT $ODS_TIME_OFFSET $DXMIN_FILE

# 
echo "end script"
# --------
  exit 0
# --------
