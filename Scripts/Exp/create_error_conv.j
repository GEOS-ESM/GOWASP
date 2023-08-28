#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=SIMERR2
#SBATCH --time=011:59:00 
#SBATCH --nodes=1
#SBATCH --partition=preops
#SBATCH --qos=dastest
#SBATCH --output=SIMERR.%j 
#SBATCH --account=s0911 
##SBATCH --qos=debug
#
# requires about 12 min/syn time for all types when some are hz or chan correl
# if extra printing for testing not set to true.
# ---------------------------------------------------------------------------
# Driver script for adding simulated observation errors to OSSE observations.
#
# Changes to make for different runs:
# 1. set SIMPARAMS to directory containing error parameters to use
# 2. set datetime1 to initial time
# 3. set ntimes to number of synoptic times to process
# 4. set proper names of sat_err_table, error.rc, and prepobs_errtable files
#    (These must match names in SIMPARAMS directory and in error.rc)
# 5. set BF_path_input to directory of input (no error) data files
# 6. set BF_path_output to directory of output (error added) data files
# 7. set which data types to process
#
#-------------------------
# Experiment environment

#-------------------------
  setenv NBRPROCS 16

  source /home/nprive/GOWASP3_ENV16
  setenv SIMPROGS   $GOWASP_PATH/Sim_error    # dir for executables 
  setenv RC_FILES   $GOWASP_PATH/Rcfiles      # dir for rc and descriptor files
  setenv SIMPARAMS  $RC_FILES/ErrParams/V04.4   # dir for error parameters
  setenv Err_file   error.V04.4.rc            # file in dir $SIMPARAMS
  setenv JCSDA_DATA $GOWASP_CRTM_COEF_DIR       # dir for CRTM files for IASI
  setenv AIRS_BUFR_TABLE $GOWASP_PATH/airs_bufr_table  # file for airs bufr tab
  setenv ADDTIME    $GOWASP_PATH/addtime.x   # file for addtime executable
  setenv RMSHMKEY   $GOWASP_RMSHMKEY         # dir for ESMF shmem clean scripts
  set sat_info_file=$RC_FILES/sat_info.rc    # file with satellite/instr info 
#
  set BF_path_input="/discover/nobackup/projects/gmao/nwposse/develop/OSSEobs/N010R"
#  set BF_path_input="/discover/nobackup/projects/gmao/nwposse/OSSEobs/TGPS3"
  set BF_path_output="/discover/nobackup/projects/gmao/nwposse/develop/OSSEobs/E013"
#  set BF_path_output="/discover/nobackup/projects/gmao/nwposse/OSSEobs/EGPS529"

  echo 'Error param file=' $SIMPARAMS $Err_file

  setenv SIMWORK $NOBACKUP/WORK/SimEwork.$$
  if (! -e $SIMWORK ) mkdir -p $SIMWORK

  umask 022
  cd  $SIMWORK

#
# USER SET VARIABLES
#---------------------------------------
# can only process within a single calendar month due to how the 
# input file paths are now set
#---------------------------------------
#
set datetime1="2006100100"   # set initial symoptic date time

@ ntimes = 124

# set which error data sets to create
#set instr_set=(PREPBUFR GPSRO AIRS AMSUAAQUA AMSUA MHS IASI ATMS GMI SSMIS CRISFSR AVCSAM AVCSPM AMSR2 HIRS4)
#set instr_set=(GPSRO)
set instr_set=(PREPBUFR)
#set instr_set=(MHS GMI AMSR2)

set test_print="F"  # T or F: print extra output for testing purposes

# ----------------------------------------------------------------------
# copy all required files to work directory
/bin/cp $ADDTIME addtime.x
/bin/cp $SIMPROGS/create_error.x prog.x
/bin/cp $SIMPARAMS/hc_params_*.bin .
/bin/cp $SIMPARAMS/sat_error_table_generic.txt .
/bin/cp $SIMPARAMS/sat_err_table_*.txt .
/bin/cp $SIMPARAMS/$Err_file error.rc    # always copy to error.rc
/bin/cp $SIMPARAMS/prepobs_errtable.*.txt .
/bin/cp $SIMPARAMS/gps_err_table_*.txt .

/bin/ls -l *

# write error.rc to log file 
echo 'cat error.rc'
cat error.rc

# set some constantsq
set sub6h=-06
set add6h=+06

# ----------------------------------------------------------------------
#                       Begin excutions of error simulator
# ----------------------------------------------------------------------

# initialize timing variables for loop

  set result=`./addtime.x $datetime1 $sub6h`
  echo 'processing data time1', $datetime1
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set datetime="${yyyy}${mm}${dd}${hh}"
  echo 'six hours before the processing data time', $datetime

# loop over synoptic times
  @ ierr_input  = 0  # initialize error counter
  @ ierr_size0  = 0  # initialize error counter
  @ ierr_output = 0  # initialize error counter
  @ ierr_other  = 0  # initialize error counter

  @ i = 1
  @ ilast =  $ntimes  
  while ( $i <= $ntimes ) 
  
# determine values for synoptic time to be processed

    set result=`./addtime.x $datetime $add6h`
    set yyyy=`echo $result | cut -d" " -f1`
    set mm=`echo $result | cut -d" " -f2`
    set dd=`echo $result | cut -d" " -f3`
    set hh=`echo $result | cut -d" " -f4`
    set datetime="${yyyy}${mm}${dd}${hh}"
    set datetime_syn="${datetime}"
    set yy=`echo $yyyy | cut -c 3-4`
#
    echo 'Add error for datetime =' $datetime
    set filetime="${yy}${mm}${dd}.t${hh}z"
#
# Loop over instruments types to consider
  foreach instrument ($instr_set)
#
# Set intrument-dependent file name information
    set d_type=$instrument
    set bufr_table_file="none "
    set file_path=$instrument/Y$yyyy/M$mm
    if ($instrument == "AIRS") then
      set bufr_table_file=$AIRS_BUFR_TABLE
      set file_path=AIRS/Y$yyyy/M$mm
      set file_tail=airsbufr_disc.$yyyy$mm$dd.t${hh}z.bufr
    else if ($instrument == "AMSUAAQUA") then
      set file_tail=amsua_disc.$yyyy$mm$dd.t${hh}z.bufr
    else if ($instrument == "IASI") then
      set file_tail=gdas1.${filetime}.mtiasi.tm00.bufr_d
    else if ($instrument == "AMSUA") then
      set file_tail=gdas1.${filetime}.1bamua.tm00.bufr_d
    else if ($instrument == "AMSUB") then
      set file_tail=gdas1.${filetime}.1bamub.tm00.bufr_d
    else if ($instrument == "HIRS2") then
      set file_tail=gdas1.${filetime}.1bhrs2.tm00.bufr_d
    else if ($instrument == "HIRS3") then
      set file_tail=gdas1.${filetime}.1bhrs3.tm00.bufr_d
    else if ($instrument == "HIRS4") then
      set file_tail=gdas1.${filetime}.1bhrs4.tm00.bufr_d
    else if ($instrument == "MSU") then
      set file_tail=gdas1.${filetime}.1bmsu.tm00.bufr_d
    else if ($instrument == "MHS") then
      set file_tail=gdas1.${filetime}.1bmhs.tm00.bufr_d
    else if ($instrument == "PREPBUFR") then
      set file_tail=gdas1.${filetime}.prepbufr
    else if ($instrument == "GPSRO") then
      set file_tail=gdas1.${filetime}.gpsro.tm00.bufr_d
    else if ($instrument == "ATMS") then
      set file_tail=gdas1.20${filetime}.atms.tm00.bufr_d
    else if ($instrument == "CRIS") then
      set file_tail=gdas1.20${filetime}.cris.tm00.bufr_d
    else if ($instrument == "CRISFSR") then
      set file_tail=gdas1.20${filetime}.crisf4.tm00.bufr_d
    else if ($instrument == "AVCSAM") then
      set file_tail=gdas1.${filetime}.avcsam.tm00.bufr_d
    else if ($instrument == "AVCSPM") then
      set file_tail=gdas1.${filetime}.avcspm.tm00.bufr_d
    else if ($instrument == "AMSR2") then
      set file_tail=gmao.amsr2_gw1_nrt.$yyyy$mm$dd.t${hh}z.bufr
    else if ($instrument == "GMI") then
      set file_tail=gmi_L1CR.$yyyy$mm$dd.t${hh}z.bufr
    else if ($instrument == "SSMIS") then
      set file_tail=gdas1.${filetime}.ssmisu.tm00.bufr_d
    else if ($instrument == "GENRADTXT") then
      set file_tail=gdas1.20${filetime}.genrad.txt
    else
      set d_type="NONE"
    endif 
#
# Create simulated obs for 1 data time and instrument
    if (d_type == "NONE") then
      echo $d_type 'MISSING FROM DATA TYPE SPECIFICATION OPTIONS'
      @ ierr_other = $ierr_other + 1
    else
      set filename="${BF_path_input}/${file_path}/${file_tail}"
      if (-e $filename ) then   # input exists and can therefore be processed
        if (! -z $filename ) then   # file has size nonzero
          rm -f file_in_bfr file_out.bfr
          ln -s $filename file_in_bfr 
          $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
          mpirun -np $NBRPROCS ./prog.x $d_type $datetime error.rc file_in_bfr file_out.bfr $bufr_table_file $JCSDA_DATA $sat_info_file $test_print
          $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
          if (-e file_out.bfr ) then # output file created, it therefore can be saved
            set path_output=${BF_path_output}/${file_path}
            if (! -e $path_output ) mkdir -p $path_output
            set filename="${path_output}/${file_tail}"
            if (-e $filename ) rm -f $filename # remove any previous file with same name
            /bin/cp file_out.bfr $filename
          else
            echo 'NO OUTPUT FILE CREATED'
            @ ierr_output = $ierr_output + 1
          endif   # test on existence of output file
        else
          echo 'REQUESTED INPUT FILE EXISTS BUT HAS SIZE 0:' $filename
          @ ierr_size0 = $ierr_size0 + 1
          @ ierr_output = $ierr_output + 1
        endif     # test on size of file
      else
        echo 'REQUESTED INPUT FILE DOES NOT EXIST:' $filename
        @ ierr_input = $ierr_input + 1
        @ ierr_output = $ierr_output + 1
      endif     # test on existence of input file
    endif       # test on d_type 
# 
  end   # end loop over instrument types

@ i++
end     # end loop over time
echo ' '
echo 'Number of nonexisting input files =' $ierr_input 
echo 'Number of existing input files with size 0 =' $ierr_size0 
echo 'Number of output files not created =' $ierr_output 
echo 'Number of other errors detected =' $ierr_other
# 
echo "end script"
exit

