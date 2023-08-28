#!/bin/csh -xf
# ------------------------------
#SBATCH --job-name=COUNTR
#SBATCH --time=02:00:00 
#SBATCH --nodes=1 
#SBATCH --output=COUNTR.%j 
#SBATCH --account=s0911 
##SBATCH --qos=debug
# ------------------------------
# Execute countobs_sum.x and countobs_tables.x to create files of tables 
# of means, standard deviations, and counts of O-F for all radiance obs 
# for a selected region. 
#
source ~/GOWASP3_ENV16
setenv PROGHOME $GOWASP_PATH/Tuning/StatsProgs  # home for executables
setenv ADDTIME  $GOWASP_PATH                        # home for addtime.x
set sat_info_file=$GOWASP_PATH/Rcfiles/sat_info.rc
#
# Set the part of input path and .ods file name that is the experiment name
# $exp is the exp name that appears as part of the GSI ods file names
# $exp_dir is part of path for both input and output that varies with exp
# It is convienent to make $exp_dir=$exp, but they need not be, e.g., 
# as when multiple experiments use the same exp name. 
# It is assumed that the input file path will have sub-directories /Yyyyy/Mmm
# proceeding the ods_path (see file_dir_in below).
#
set exp="real529"
set exp_dir=$exp                         
set ods_path="/discover/nobackup/projects/gmao/nwposse/baselineruns/${exp_dir}/obs" 
#set exp="osse529"
#set exp_dir=$exp                         
#set ods_path="/discover/nobackup/projects/gmao/nwposse/develop/${exp_dir}/obs" 

# Set datetime (remember real and osse have different years) and region
set datetime1="2020062300"  # starting datetime yyyymmddhh                   
@ ntimes =  20          # number of synoptic times to process            

# Set the output directory
setenv OUTPUT /discover/nobackup/projects/gmao/nwposse/develop/ODSstats/${exp_dir}/Counts_rad_June_5d

# Set instr_platform for obs types to process
set kx_set= (amsua_n15 amsua_n18 cris-fsr_npp amsua_n19 amsua_metop-a amsua_metop-b hirs4_metop-a hirs4_n18 mhs_metop-a mhs_metop-b airs_aqua iasi_metop-a iasi_metop-b cris-fsr_n20 atms_npp atms_n20 ssmis_f17 amsr2_gcom-w1 avhrr_metop-a avhrr_n18 avhrr_n19 gmi_gpm)
#set kx_set= (airs_aqua)
#
# END OF USUAL USER-DEFINED VARIABLES
#
setenv WORKDIR /discover/nobackup/$user/WORK/CNTwork.$$ 
set addhours="+06"          # hours between times to consider (e.g., +06) 
set region="GLOB"           # geo region to compute correls (see below)    
set sfc_type="BOTH"         # SEA, LAND, or BOTH

if ($region == 'GLOB' ) then      # Global
  set latN="  90.0"     # northern-most latitude in subdomain considered   
  set latS=" -90.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'SSEA' ) then      # Southern ocean
  set latN=" -40.0"     # northern-most latitude in subdomain considered   
  set latS=" -60.0"     # southern-most latitude in subdomain considered   
  set lonW="   0.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'WHEM' ) then      # Western Hemisphere
  set latN="  90.0"     # northern-most latitude in subdomain considered   
  set latS=" -90.0"     # southern-most latitude in subdomain considered   
  set lonW=" 180.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
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
else if ($region == 'RUSS' ) then  # Russia
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
cp $PROGHOME/countobs_sum.x prog_sum.x
cp $PROGHOME/countobs_tables.x prog_tables.x
cp $ADDTIME/addtime.x .

foreach instrument ($kx_set)

  echo $instrument

  set sat_inst=`echo $instrument | cut -d_ -f 1-1`  # extract string before "_" denoting obs type
  set ods_name="${exp}.diag_${instrument}."
  set output_file="count_${exp}_${instrument}"  

# determine first set of date and time variables
  set result=`./addtime.x $datetime1 +00`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set datetime="${yyyy}${mm}${dd}${hh}"

  rm -f input_stats
  rm -f correl_printout

# initialize counters
  @ i = 1
  @ j = 0
  @ i_no = 0
  @ ilast =  $ntimes  

# loop over times 
  while ( $i <= $ntimes ) 

    echo "${datetime}"  

    set file_dir_in="${ods_path}/Y${yyyy}/M${mm}/D${dd}/H${hh}"
    set file_name="${file_dir_in}/${ods_name}${yyyy}${mm}${dd}_${hh}z.ods" 

    if (-e $file_name) then

      @ j = $j + 1
      ./prog_sum.x ${j} ${yyyy}${mm}${dd} ${hh}0000 ${latN} ${latS} ${lonW} ${lonE} 0 0 $sfc_type ${file_name} ${sat_inst} ${sat_info_file} >> count_printout

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

# end loop

    @ i++
  end

  @ i = $i - 1
  echo "${i} times considered"
  echo "${j} times processed"
  echo "${i_no} = number of files not existing"
  ./prog_tables.x input_stats $output_file  > count_printout
  echo "Correlations program executed"
  if (! -e $OUTPUT ) mkdir -p $OUTPUT
  /bin/mv ${output_file}.txt ${OUTPUT}/${output_file}_${region}_${sfc_type}.txt
  /bin/mv count_printout ${OUTPUT}/count_printout
  echo "Countobs file moved"
  echo "end consideration of instrument ${instrument}"
 
end  # loop over instruments

date
# --------
exit 0
# --------
