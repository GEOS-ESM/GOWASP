#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=SIMRADL
#SBATCH --time=1:59:00 
#SBATCH --nodes=1  --constraint=cas  
#SBATCH --output=SIMRADL.%j 
#SBATCH --account=s0818
#SBATCH --qos=obsdev

# Script for creating binary files of possibly thinned header information 
# drawn from obs BUFR files for ingestion into create_rad_profs 
# Requires 12 min per synoptic time for 10 rad types

  source /home/mkim1/GOWASP3_ENV16
  setenv SIMHOME    $GOWASP_PATH/Sim_rad  # dir for some files/executables
  setenv ADDTIME    $GOWASP_PATH/addtime.x  # addtime executable
  setenv AIRS_BUFR_TABLE $GOWASP_PATH/airs_bufr_table  # file for air bufr table
  setenv SIMWORK    $NOBACKUP/WORK/SimRLwork.$$ 
  setenv RC_FILE    $GOWASP_PATH/Rcfiles/rad_thin_N010D.rc  # rc file
  setenv BUFR_IN    /discover/nobackup/projects/gmao/nwposse/develop/BUFRDATA/       # dir for real BUFR data 
#  setenv BUFR_IN    /discover/nobackup/rerrico/BUFRDATA/       # dir for real BUFR data 
#  setenv OUTPUT     /discover/nobackup/projects/gmao/nwposse/develop//OSSEobs/N010/ObsList # dir for output 
  setenv OUTPUT     /discover/nobackup/projects/gmao/obsdev/mkim1/OSSEobs/N010/ObsList # dir for output 
  setenv REALLOC    F # T or F; T if obs to use thinned and QCd real locations
#
  if ($REALLOC == T ) then  
    setenv RC_FILE    $GOWASP_PATH/Rcfiles/rad_nothin.rc  # rc file 
    setenv SAT_PARAMS $GOWASP_PATH/Rcfiles/sat_info.rc
    setenv ODSEXP     517test                                          
    setenv ODSPATH    /discover/nobackup/rerrico/ODSfiles/$ODSEXP     
    setenv ODSTXT     $SIMHOME/Rad_ods
  endif
#
  set HDR1="Exp"               # comments to add to file to identify run 
  set HDR2="rad_thin_N010.rc"     # comments to add to file to identify run 
#
# set NR and OBS (BUFR) starting date-time and number of times to process      
# Output simulated OBS files will have same time stamps as input OBS files     
# ADJUST datemine 1 FOR LEAP YEAR IN BUFR DATA IF NECESSARY!!                 
  
set datetimeN=2006103100  # starting time for NR           
set datetimeB=2020101100  # starting time for BUFR data    
@ ntimes = 4          # number of times to process                              
#   
# set which radiance data sets to create 
# GMI must be done separately as long as the dates differ (must use after 2014) 
#set instr_set= (AIRS AMSUAAQUA AMSUA CRISFSR HIRS4 MHS IASI ATMS GMI SSMIS  AMSR2)
set instr_set= ( AMSUA )
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

cp $SIMHOME/create_rad_obs_list.x prog.x
cp $ADDTIME addtime.x
if ($REALLOC == T ) then 
  cp $ODSTXT/rad_ods2txt.x ods2txt.x 
endif  
#
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
# set time for data to be created (must be same as NR time)
  set result=`./addtime.x $datetimeN $addh`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set datetimeN="${yyyy}${mm}${dd}${hh}"    # NR synoptic time
  @ yy = $yyyy % 10 
  set yy=`printf "%2.2d\n" $yy`
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
  echo 'Simulation i=' $i 'uses NR time' $datetimeN '& Obs time' $datetimeB
#
# Loop over instruments types to consider
  foreach instrument ($instr_set)
#
# Set intrument-dependent file name information
#
# Set name of output directory for this instrument
    set BF_dir_name=$instrument
    set obslist_dir=$OUTPUT/$BF_dir_name/Y$yyyy/M$mm
#
# Set name of instrument in portion of the file name
    set bufr_table_file='none'   # default
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
# set file names for input and output
    if ($instrument == "AIRS" || $instrument == "AMSUAAQUA" || $instrument == "GMI" || $instrument == "AMSR2" ) then    
      set BF_file_in=${BUFR_IN}/Y${yyyyB}/M${mmB}/${BF_file_instr}.${yyyyB}${mmB}${ddB}.t${hhB}z.bufr
    else if ($instrument == "ATMS" || $instrument == "CRIS" || $instrument == "CRISFSR"  ) then    
      set BF_file_in=${BUFR_IN}/Y${yyyyB}/M${mmB}/gdas1.20${datetimeBP}.${BF_file_instr}.tm00.bufr_d
    else if ($instrument == "GENRADTXT") then
      set BF_file_in=${BUFR_IN}/Y${yyyyB}/M${mmB}/gdas1.20${datetimeBP}.${BF_file_instr}.txt
    else
      set BF_file_in=${BUFR_IN}/Y${yyyyB}/M${mmB}/gdas1.${datetimeBP}.${BF_file_instr}.tm00.bufr_d
    endif
    set obslist_file_out=${obslist_dir}/obslist_${instrument}.${datetimeN}.bin 
#
# Create obs list of header and thiining info for 1 data time and instrument
# 
    if ($BF_file_instr == "NONE") then
      echo $instrument 'MISSING FROM FILE SPECIFICATION LIST'
      echo 'Skip processing this choice'
      @ ierr_other = $ierr_other + 1
    else
      if (-e $BF_file_in) then  # input data file exists and can be used   
#
        if ($REALLOC == T ) then  
          ./ods2txt.x ${datetimeB}0000 $SAT_PARAMS $ODSPATH $ODSEXP $instrument $test_print
        endif
#
        /bin/rm -f data_out 
#
        ./prog.x $instrument $datetimeB $cdtime $RC_FILE $BF_file_in data_out $bufr_table_file $HDR1 $HDR2 $test_print $REALLOC
#
        if ($REALLOC == T ) then  
          set RADLOCS=radlocs_${instrument}.txt
          /bin/rm $RADLOCS
        endif
#
        if (-e data_out ) then # output file created and can be saved
          if (! -e $obslist_dir ) mkdir -p ${obslist_dir}  # create output DIR 
          /bin/cp data_out $obslist_file_out
        else
          echo 'NO DATA OUTPUT'
          @ ierr_output = $ierr_output + 1
        endif        
      else
        echo 'NONEXISTING INPUT FILE= ' $BF_file_in
        @ ierr_input = $ierr_input + 1
        @ ierr_output = $ierr_output + 1
      endif
    endif
#  
  end  # loop over instrument
#
  @ i++
end    # loop over times  
# All done
echo ' '
echo 'Number of nonexisting input files=' $ierr_input 
echo 'Number of output files not created =' $ierr_output 
echo 'Number of other errors detected=' $ierr_other
# 
echo "end script"
date

# --------
  exit 0
# --------
