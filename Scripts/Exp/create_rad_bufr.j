#!/bin/csh -xf 
#SBATCH --job-name=SIMRADB
#SBATCH --time=1:59:00 
##SBATCH --partition=preops
##SBATCH --qos=dastest
#SBATCH --nodes=1   --constraint=cas 
#SBATCH --output=SIMRADB.%j 
#SBATCH --account=s0818
#SBATCH --qos=obsdev

# 8 min/synoptic time for all 10 types
# Specify modules and number of mpi processors
setenv NBRPROCS 18    
echo "NBRPROCS = ${NBRPROCS}"

# Script for creating binary files of possibly thinned header information 
# drawn from obs BUFR files for ingestion into create_rad_profs 

  source ~/GOWASP3_ENV16
  setenv SIMHOME    $GOWASP_PATH/Sim_rad    # dir for some files/executables
  setenv ADDTIME    $GOWASP_PATH/addtime.x  # addtime executable
  setenv CRTM_COEF_DIR $GOWASP_CRTM_COEF_DIR
  setenv AIRS_BUFR_TABLE $GOWASP_PATH/airs_bufr_table  # file for air bufr table
  setenv SIMWORK    $NOBACKUP/WORK/SimRBwork.$$
  setenv RC_FILE    $GOWASP_PATH/Rcfiles/rad_prob_N010E.rc  # rc file 
  setenv SAT_PARAMS $GOWASP_PATH/Rcfiles/sat_info.rc # file sat/instr params
#  setenv BUFR_IN    $GOWASP_BUFR_DATA  # dir for real BUFR data to get buf tabl
  setenv BUFR_IN    /discover/nobackup/projects/gmao/nwposse/develop/BUFRDATA/   # dir for real BUFR data to get buf tabl
#  setenv BUFR_IN    /discover/nobackup/rerrico/BUFRDATA/   # dir for real BUFR data to get buf tabl
  setenv PROF_IN    /discover/nobackup/projects/gmao/obsdev/mkim1/OSSEobs/N010/ObsProf
  setenv BUFR_OUT   /discover/nobackup/projects/gmao/obsdev/mkim1/OSSEobs/N010   # dir for BUFR output
  setenv CLEAR_SKY  none # file name suffix if obs all-sky but clear-sky Tb witten .bin 
  setenv REALLOC    F  # T or F; T means use thinned and QCd locs in ods file
#
  if ($REALLOC == T ) then  
    setenv RC_FILE    $GOWASP_PATH/Rcfiles/rad_prob_none.rc  # rc file 
    setenv ODSEXP     517test                                                
    setenv ODSPATH    /discover/nobackup/rerrico/ODSfiles/$ODSEXP            
    setenv ODSTXT     $SIMHOME/Rad_ods
  endif
#
# For datetimeB, if REALLOC=F then this can be any time a BUFR file is available,
# since it only need the BUFR table; however, if REALLOC=T then the sequence of 
# times must match those used by create_rad_list.j since it will attempt to 
# match the BUFR and ods obs locations. 
#
set datetimeN=2006103100  # starting time to process
set datetimeB=2020103100  # set as explained above
@ ntimes =    4          # number of times to process

# set which radiance data sets to create
#set instr_set= (AIRS AMSUAAQUA AMSUA CRISFSR HIRS4 IASI ATMS SSMIS )
#set instr_set= (GMI AMSR2 MHS)
set instr_set= (AMSUA)
#set instr_set= (CRISFSR)
set test_print="F"  # T or F: print sample output for testing purposes

# END OF USER SET VARIABLES
# ----------------------------------------------------------------------

# Make sure files to be written are acessible to others
umask 022

# Increment in hours between central synoptic input times
set add6h=+06

# Create working directory
if (! -e $SIMWORK            ) mkdir -p $SIMWORK

cd $SIMWORK
/bin/rm -f  *

cp $SIMHOME/create_rad_bufr.x prog.x
cp $ADDTIME addtime.x
if ($REALLOC == T ) then  
  cp $ODSTXT/rad_ods2txt.x ods2txt.x 
endif

# x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
# loop over synoptic times
@ ierr_input  = 0  # initialize error counter
@ ierr_output = 0  # initialize error counter
@ ierr_other  = 0  # initialize error counter

@ i = 1
@ ilast =  $ntimes  
while ( $i <= $ntimes ) 
  
  if ($i == 1) then
    set addh=000
  else
    set addh=$add6h 
  endif
#
# set time for data to be created 
  set result=`./addtime.x $datetimeN $addh`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set datetimeN="${yyyy}${mm}${dd}${hh}"    # current synoptic time considered
  @ yy = $yyyy % 10 
  set yy=`printf "%2.2d\n" $yy`
  set datetimeNP="${yy}${mm}${dd}.t${hh}z"  # part of output BUFR file name
  set cdtime="${datetimeN}0000"             # 00 min and 00 sec indicated
#
# set time for bufr data file to be input
  set result=`./addtime.x $datetimeB $addh`
  set yyyyB=`echo $result | cut -d" " -f1`
  set yyB=`echo $yyyyB | cut -c 3-4`
  set mmB=`echo $result | cut -d" " -f2`
  set ddB=`echo $result | cut -d" " -f3`
  set hhB=`echo $result | cut -d" " -f4`
  set datetimeB="${yyyyB}${mmB}${ddB}${hhB}"    # time for BUFR input data
  set datetimeBP="${yyB}${mmB}${ddB}.t${hhB}z"  # part of input BUFR file name
#
#
# Loop over instruments types to consider
  foreach instrument ($instr_set)
#
# Set intrument-dependent file name information
#
# Set name of input directory and file for this instrument
    set prof_dir=$PROF_IN/$instrument/Y$yyyy/M$mm
    set prof_file=${prof_dir}/profile_${instrument}.${datetimeN}.bin 
#
# Set name of instrument in portion of the file name
    set bufr_table_file='none'
    if ($instrument == "AIRS") then
      set BF_file_instr=airsbufr_disc
      set bufr_table_file=$AIRS_BUFR_TABLE 
    else if ($instrument == "AMSUAAQUA") then
      set BF_file_instr=amsua_disc
    else if ($instrument == "IASI") then
      set BF_file_instr=mtiasi
    else if ($instrument == "AMSUA") then
      set BF_file_instr=1bamua
    else if ($instrument == "AMSUB") then
      set BF_file_instr=1bamub
    else if ($instrument == "HIRS2") then    
      set BF_file_instr=1bhrs2
    else if ($instrument == "HIRS3") then    
      set BF_file_instr=1bhrs3
    else if ($instrument == "HIRS4") then
      set BF_file_instr=1bhrs4
    else if ($instrument == "MSU") then
      set BF_file_instr=1bmsu
    else if ($instrument == "MHS") then
      set BF_file_instr=1bmhs
    else if ($instrument == "ATMS") then
      set BF_file_instr=atms
    else if ($instrument == "GMI") then
      set BF_file_instr=gmi_L1CR
    else if ($instrument == "CRIS") then
      set BF_file_instr=cris
    else if ($instrument == "CRISFSR") then
      set BF_file_instr=crisf4
    else if ($instrument == "AVCSAM") then
      set BF_file_instr=avcsam
    else if ($instrument == "AVCSPM") then
      set BF_file_instr=avcspm
    else if ($instrument == "AMSR2") then
      set BF_file_instr=gmao.amsr2_gw1_nrt
    else if ($instrument == "SSMIS") then
      set BF_file_instr=ssmisu
    else if ($instrument == "GENRADTXT") then
      set BF_file_instr=genrad
    else
      set BF_file_instr=NONE
    endif
#
# set file directory and name for input and output
    set BF_dir_in=${BUFR_IN}/Y${yyyyB}/M${mmB}
    set BF_dir_out=${BUFR_OUT}/${instrument}/Y${yyyy}/M${mm}
    if ($instrument == "AIRS" || $instrument == "AMSUAAQUA" || $instrument == "GMI" || $instrument == "AMSR2" ) then    
      set BF_file_in=${BF_dir_in}/${BF_file_instr}.${yyyyB}${mmB}${ddB}.t${hhB}z.bufr
      set BF_file_out=${BF_dir_out}/${BF_file_instr}.${yyyy}${mm}${dd}.t${hh}z.bufr
    else if ($instrument == "ATMS" || $instrument == "CRIS" || $instrument == "CRISFSR" ) then    
      set BF_file_in=${BF_dir_in}/gdas1.20${datetimeBP}.${BF_file_instr}.tm00.bufr_d
      set BF_file_out=${BF_dir_out}/gdas1.20${datetimeNP}.${BF_file_instr}.tm00.bufr_d
    else if ($instrument == "GENRADTXT") then
      set BF_file_in=${BF_dir_in}/gdas1.20${datetimeBP}.${BF_file_instr}.txt
      set BF_file_out=${BF_dir_out}/gdas1.20${datetimeNP}.${BF_file_instr}.txt
    else 
      set BF_file_in=${BF_dir_in}/gdas1.${datetimeBP}.${BF_file_instr}.tm00.bufr_d
      set BF_file_out=${BF_dir_out}/gdas1.${datetimeNP}.${BF_file_instr}.tm00.bufr_d     
    endif
#
# Create bufr file for simulated obs for 1 data time and instrument
# 
    if ($BF_file_instr == "NONE") then
      echo $instrument 'missing from file specifications'
      @ ierr_other = $ierr_other + 1
    else
      if (-e $prof_file && -e $BF_file_in ) then  # input data files exist and can be used   
        if ($REALLOC == T ) then  
          ./ods2txt.x ${datetimeB}0000 $SAT_PARAMS $ODSPATH $ODSEXP $instrument $test_print
        endif
        /bin/rm -f data_out 
# 
        mpirun -np $NBRPROCS ./prog.x $instrument $cdtime $RC_FILE $SAT_PARAMS $prof_file $BF_file_in $bufr_table_file data_out $CRTM_COEF_DIR $test_print $REALLOC $CLEAR_SKY
#
        if ($REALLOC == T ) then  
          set RADLOCS=radlocs_${instrument}.txt
          /bin/rm $RADLOCS
        endif
#
        if (-e data_out ) then # output file created and can be saved
          if (! -e $BF_dir_out ) mkdir -p $BF_dir_out  # create output DIR 
          /bin/cp data_out $BF_file_out
          if (-e data_out$CLEAR_SKY) then
            /bin/mv data_out${CLEAR_SKY} ${BF_file_out}${CLEAR_SKY}
          endif
        else
          echo 'NO DATA OUTPUT'
          @ ierr_output = $ierr_output + 1
        endif        
      else
        if (! -e $BF_file_in) then  
          echo 'NONEXISTING INPUT FILE= ' $BF_file_in
          echo 'Skip processing this file'
          @ ierr_input = $ierr_input + 1
        endif
        if (! -e $prof_file) then  
          echo 'NONEXISTING INPUT FILE= ' $prof_file
          echo 'Skip processing this file'
          @ ierr_input = $ierr_input + 1
        endif
        @ ierr_output = $ierr_output + 1
      endif
    endif
#  
  end  # loop over instrument
#
  @ i++
end    # loop over times  
# All done
echo "end script"
echo ' '
echo 'Number of nonexisting input files=' $ierr_input 
echo 'Number of output files not created =' $ierr_output 
echo 'Number of other errors detected=' $ierr_other
#

# --------
  exit 0
# --------
