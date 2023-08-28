#!/bin/csh -fx 
# ------------------------------
#SBATCH --job-name=VCORRG
#SBATCH --time=00:59:00 
#SBATCH --nodes=1   
#SBATCH --output=VCORRG.%j 
#SBATCH --account=s0911 
#SBATCH --qos=debug 
#

# Execute vcorr_sumods.x and vcorr_tables.x to create files of tables of means, 
# standard deviations, and vertical correlations of O-F 

# Set user variables. For every run, specify values indicated by "&&"
source /home/nprive//GOWASP3_ENV16
setenv PROGHOME $GOWASP_PATH/Tuning/StatsProgs # home for executables
setenv ADDTIME  $GOWASP_PATH  # home for addtime.x
#
#set exp="ctl529"
#set exp_dir=$exp
#set ods_path="/discover/nobackup/projects/gmao/nwposse/baselineruns/${exp_dir}/obs"
set exp="ctl529"
set exp_dir=$exp                         
set ods_path="/discover/nobackup/projects/gmao/nwposse/baselineruns/${exp_dir}/obs" 

set ods_name="${exp}.diag_conv."    # "${exp}.diag_conv."
set omf="obs"                       # omf, oma, amb, or obs 

setenv OUTPUT /discover/nobackup/projects/gmao/nwposse/develop/ODSstats/${exp_dir}/Vcorr_June_5d # output dir 
setenv WORKDIR /discover/nobackup/$user/WORK/CORwork.$$ 

# Note that if addhours=+12 then this is appropriate for raobs also.
# In contrast, if addhours=+06, then the time averaging will include times 
# when almost no raobs are present, so that the time mean will be approx. 
# half what it is when addhours=+12. Using addhours=+12, however, excludes
# half of the available gpsro obs.

set region="GLOB"           # geo region to compute correls (see below)      &&
#set datetime0="2006062300"  # starting datetime yyyymmddhh                   &&
set datetime0="2006062300"  # starting datetime yyyymmddhh                   &&
set addhours="+12"     # hours between times to consider                    &&
@ ntimes =  10             # number of synoptic times to process          &

set region_set= (LN75 LN50 LN30 LN10 LS75 LS50 LS30 LS10)

# If this is a restart intended to use files of previously calculated and saved sums
# to now add more times to the sum, several vlaues must be specified:
# RESTART_DIR is the direcotry containing the previously saved files of sums
# RESTART_I is the number of times previously considered 
# set ntimes above to the number of times considered in the combined datasets 
#     (i.e., for the previously plus now to be considered times)
# set datetime1 above to the beginning time for the new times to be considered
setenv RESTART_DIR /discover/nobackup/rerrico/ODSstats/${exp_dir}/Vcorr_June_5d/Save_sums
setenv RESTART_I 0  # set to 0 unless this is a restart

#
# END OF USUAL USER-DEFINED VARIABLES
#
set ibins="90"         # number of vertical bins (3 digits)

# Create and move to empty working directory
if (! -e $WORKDIR            ) mkdir -p $WORKDIR
cd $WORKDIR
/bin/rm -f  *

# Copy file to working directory
cp $PROGHOME/vcorr_sumods.x vcorr_sumods.x  
cp $PROGHOME/vcorr_tables.x .
cp $PROGHOME/read_j_value.x .
cp $ADDTIME/addtime.x .

foreach region ($region_set)
set datetime1=$datetime0

if ($region == 'GLOB' ) then      # Global
  set latN="  90.0"     # northern-most latitude in subdomain considered   
  set latS=" -90.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'LN75' ) then 
  set latN="  90.0"     # northern-most latitude in subdomain considered   
  set latS="  60.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'LN50' ) then 
  set latN="  60.0"     # northern-most latitude in subdomain considered   
  set latS="  40.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'LN30' ) then 
  set latN="  40.0"     # northern-most latitude in subdomain considered   
  set latS="  20.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'LN10' ) then 
  set latN="  20.0"     # northern-most latitude in subdomain considered   
  set latS="    .0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'LS75' ) then 
  set latN=" -60.0"     # northern-most latitude in subdomain considered   
  set latS=" -90.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'LS50' ) then 
  set latN=" -40.0"     # northern-most latitude in subdomain considered   
  set latS=" -60.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'LS30' ) then 
  set latN=" -20.0"     # northern-most latitude in subdomain considered   
  set latS=" -40.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'LS10' ) then 
  set latN="    .0"     # northern-most latitude in subdomain considered   
  set latS=" -20.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
endif

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
ls -l
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

ls -l

    if (! -e $OUTPUT ) mkdir -p $OUTPUT
    /bin/mv $output_file ${OUTPUT}/${output_file}_${region}_${omf}_lats.txt
    /bin/mv correl_printout ${OUTPUT}/correl_printout.txt
    echo "Correlation file copied"
#
# Save restart files (files of sums over time)
    set savedir=$OUTPUT/Save_sums
    if (! -e $savedir ) mkdir -p $savedir
    set savefile=$savedir/save_${output_file}_${region}_${omf}_$i.bin
    /bin/mv input_stats $savefile

end # loop over regions
# --------
  echo "end script"
  exit 0
# --------
