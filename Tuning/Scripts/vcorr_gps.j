#!/bin/csh -fx 
# ------------------------------
#SBATCH --job-name=VCORRG
#SBATCH --time=03:59:00 
#SBATCH --nodes=1 
#SBATCH --output=VCORRG.%j 
#SBATCH --account=s0911 
##SBATCH --qos=debug 
#

# Execute vcorr_sumods.x and vcorr_tables.x to create files of tables of means, 
# standard deviations, and vertical correlations of O-F 

# Set user variables. For every run, specify values indicated by "&&"
source /home/nprive/GOWASP3_ENV16
setenv PROGHOME $GOWASP_PATH/Tuning/StatsProgs # home for executables
setenv ADDTIME  $GOWASP_PATH  # home for addtime.x
#
#set exp="4dr_asky"
#set exp_dir=$exp                         
#set ods_path="/discover/nobackup/projects/gmao/nwposse/archive/${exp_dir}/obs" 
set exp="real529"
set exp_dir=$exp                         
set ods_path="/discover/nobackup/projects/gmao/nwposse/baselineruns/${exp_dir}/obs" 

set ods_name="${exp}.diag_conv."    # "${exp}.diag_conv."
set omf="omf"                       # omf, oma, amb, or obs 

setenv OUTPUT /discover/nobackup/projects/gmao/nwposse/develop/ODSstats/${exp_dir}/Vcorr_conv_5d # output dir 
setenv WORKDIR /discover/nobackup/$user/WORK/CORwork.$$ 

# Note that if addhours=+12 then this is appropriate for raobs also.
# In contrast, if addhours=+06, then the time averaging will include times 
# when almost no raobs are present, so that the time mean will be approx. 
# half what it is when addhours=+12. Using addhours=+12, however, excludes
# half of the available gpsro obs.

set region="GLOB"           # geo region to compute correls (see below)      &&
#set datetime1="2006062300"  # starting datetime yyyymmddhh                   &&
set datetime1="2020062300"  # starting datetime yyyymmddhh                   &&
set addhours="+12"     # hours between times to consider                    &&
@ ntimes =  10             # number of synoptic times to process          &

# If this is a restart intended to use files of previously calculated and saved sums
# to now add more times to the sum, several vlaues must be specified:
# RESTART_DIR is the direcotry containing the previously saved files of sums
# RESTART_I is the number of times previously considered 
# set ntimes above to the number of times considered in the combined datasets 
#     (i.e., for the previously plus now to be considered times)
# set datetime1 above to the beginning time for the new times to be considered
setenv RESTART_DIR /discover/nobackup/projects/gmao/nwposse/develop/ODSstats/${exp_dir}/Vcorr_conv_5d/Save_sums
setenv RESTART_I 0  # set to 0 unless this is a restart

#
# END OF USUAL USER-DEFINED VARIABLES
#
set ibins="90"         # number of bins (3 digits)

if ($region == 'GLOB' ) then      # Global
  set latN="  90.0"     # northern-most latitude in subdomain considered   
  set latS=" -90.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'GPS1' ) then 
  set latN="  40.0"     # northern-most latitude in subdomain considered   
  set latS="  39.0"     # southern-most latitude in subdomain considered   
  set lonW=" 116.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 117.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'PENN' ) then 
  set latN="  45.0"     # northern-most latitude in subdomain considered   
  set latS="  35.0"     # southern-most latitude in subdomain considered   
  set lonW=" 280.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 290.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'GPST' ) then 
  set latN=" -20.0"     # northern-most latitude in subdomain considered   
  set latS=" -40.0"     # southern-most latitude in subdomain considered   
  set lonW="  40.0"     # western-most longitude in subdomain considered 0-360
  set lonE="  90.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'TNPC' ) then # Tropical and North Pacific
  set latN="  55.0"     # northern-most latitude in subdomain considered   
  set latS=" -20.0"     # southern-most latitude in subdomain considered   
  set lonW=" 165.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 225.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'TPAC' ) then # Tropical Pacific
  set latN="  20.0"     # northern-most latitude in subdomain considered   
  set latS=" -20.0"     # southern-most latitude in subdomain considered   
  set lonW=" 165.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 225.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'NHET' ) then # Northern Hemisphere Extra-Tropics
  set latN="  90.0"     # northern-most latitude in subdomain considered   
  set latS="  20.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'N45A' ) then # Northern Hemisphere Extra-Tropics
  set latN="  90.0"     # northern-most latitude in subdomain considered   
  set latS="  45.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'NPAC' ) then # North Pacific
  set latN="  55.0"     # northern-most latitude in subdomain considered   
  set latS="  20.0"     # southern-most latitude in subdomain considered   
  set lonW=" 165.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 225.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'SPAC' ) then  # South Pacific
  set latN=" -20.0"     # northern-most latitude in subdomain considered   
  set latS=" -55.0"     # southern-most latitude in subdomain considered   
  set lonW=" 165.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 225.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'NAMR' ) then  # North America
  set latN="  60.0"     # northern-most latitude in subdomain considered   
  set latS="  30.0"     # southern-most latitude in subdomain considered   
  set lonW=" 240.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 290.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'TROP' ) then  # Tropics
  set latN="  20.0"     # northern-most latitude in subdomain considered   
  set latS=" -20.0"     # southern-most latitude in subdomain considered   
  set lonW="   0.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'SHET' ) then  # Southern Hemisphere Extra-Tropics
  set latN=" -20.0"     # northern-most latitude in subdomain considered   
  set latS=" -55.0"     # southern-most latitude in subdomain considered   
  set lonW="   0.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'RUSS' ) then  # Southern Hemisphere Extra-Tropics
  set latN="  70.0"     # northern-most latitude in subdomain considered   
  set latS="  50.0"     # southern-most latitude in subdomain considered   
  set lonW="  20.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 160.0"     # eastern-most longitude in subdomain considered 0-360
endif

# Create and move to empty working directory
if (! -e $WORKDIR            ) mkdir -p $WORKDIR
cd $WORKDIR
/bin/rm -f  *

# Copy file to working directory
cp $PROGHOME/vcorr_sumods.x vcorr_sumods.x  
cp $PROGHOME/vcorr_tables.x .
cp $PROGHOME/read_j_value.x .
cp $ADDTIME/addtime.x .

    set output_file="vzcor_${exp}_conv"
 # determine first set of date and time variables
    set result=`./addtime.x $datetime1 +00`
    set yyyy=`echo $result | cut -d" " -f1`
    set mm=`echo $result | cut -d" " -f2`
    set dd=`echo $result | cut -d" " -f3`
    set hh=`echo $result | cut -d" " -f4`
    set datetime="${yyyy}${mm}${dd}${hh}"

    rm -f input_stats
    rm -f correl_printout

# If not a restart just set $RESTART_J to 0, otherwise get previously saved 
# file of sums to provide the first input_stats file. Read the starting value 
# of j from it (This is the number of existing times for the particular data 
# type in the times previously considered).
    if ($RESTART_I == 0) then
      set RESTART_J=0
    else
      cp $RESTART_DIR/save_${output_file}_${region}_${omf}_${RESTART_I}.bin input_stats
      set RESTART_J=`./read_j_value.x input_stats`
    endif 

# initialize counters
    @ i = 1 + $RESTART_I
    @ j = $RESTART_J
    @ i_no = 0
    @ ilast =  $ntimes  

# loop over times 
    while ( $i <= $ntimes ) 

      echo "${datetime}"  

      set file_dir_in="${ods_path}/Y${yyyy}/M${mm}/D${dd}/H${hh}"
      set file_name="${file_dir_in}/${ods_name}${yyyy}${mm}${dd}_${hh}z.ods" 
      if (-e $file_name) then

        @ j = $j + 1
        ./vcorr_sumods.x ${j} ${ibins} ${yyyy}${mm}${dd} ${hh}0000 ${latN} ${latS} ${lonW} ${lonE} ${omf} ${file_name} >> correl_printout

        if (-e output_stats) then
          if (-e input_stats) then
            rm -f input_stats
          endif
          mv output_stats input_stats
        else
          @ j = $j - 1
        endif

      else
        @ i_no = $i_no + 1
        echo "FILE DOES NOT EXIST for datetime ${datetime}"
      endif

      set result=`./addtime.x $datetime $addhours`
      set yyyy=`echo $result | cut -d" " -f1`
      set mm=`echo $result | cut -d" " -f2`
      set dd=`echo $result | cut -d" " -f3`
      set hh=`echo $result | cut -d" " -f4`
      set datetime="${yyyy}${mm}${dd}${hh}"

      @ i++
    end  # loop over time

    @ i = $i - 1
    echo "${i} times considered"
    echo "${j} times processed"
    echo "${i_no} = number of files not existing"

    ./vcorr_tables.x input_stats $output_file >> correl_printout
    echo "Correlations program executed"

    if (! -e $OUTPUT ) mkdir -p $OUTPUT
    /bin/mv $output_file ${OUTPUT}/${output_file}_${region}_${omf}.txt
    /bin/mv correl_printout ${OUTPUT}/correl_printout.txt
    echo "Correlation file copied"
#
# Save restart files (files of sums over time)
    set savedir=$OUTPUT/Save_sums
    if (! -e $savedir ) mkdir -p $savedir
    set savefile=$savedir/save_${output_file}_${region}_${omf}_$i.bin
    /bin/mv input_stats $savefile

# --------
  echo "end script"
  exit 0
# --------
