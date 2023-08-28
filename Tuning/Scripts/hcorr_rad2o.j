#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=HCRAD
#SBATCH --time=05:59:00 
#SBATCH --nodes=1 
#SBATCH --output=HCRAD.%j 
#SBATCH --account=s0911 
##SBATCH --qos=debug

# Execute hcorr_sumods.x and hcorr_tables.x to create files of tables 
# of means, standard deviations, and horizontal correlations of O-F 
# for all radiance obs for a selected region. 
#
# Requires about 40 min per day (assuming 4 times per day) if computed 
# for the entire GLOBE and without IASI, AIRS/ and CRIS
# For IASI, AIRS and CRIS computed for 1 hemisphere, 50min/day is required.
#
source ~/GOWASP3_ENV16
setenv PROGHOME $GOWASP_PATH/Tuning/StatsProgs # home for executables
setenv ADDTIME  $GOWASP_PATH  # home for addtime.x
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
#set exp="4dr_asky"
#set exp_dir=$exp
#set ods_path="/discover/nobackup/projects/gmao/nwposse/archive/${exp_dir}/obs"
set exp="osse529"
set exp_dir=$exp
set ods_path="/discover/nobackup/projects/gmao/nwposse/develop/${exp_dir}/obs"

# Set datetime (remember real and osse have different years) and region
set datetime1="2006062300"  # starting datetime yyyymmddhh                   
#set datetime1="2020062300"  # starting datetime yyyymmddhh                   
@ ntimes =20                # number of synoptic times to process            
#set region="GLOB"           # geo region to compute correls (see below)      
 set region="WHEM"     # for hyperspectral

# Set the output directory
setenv OUTPUT /discover/nobackup/projects/gmao/nwposse/develop/ODSstats/${exp_dir}/Corr_rad_June_5d

# If this is a restart intended to use files of previously calculated and saved sums
# to now add more times to the sum, several vlaues must be specified:
# RESTART_DIR is the direcotry containing the previously saved files of sums
# RESTART_I is the number of times previously considered 
# set ntimes above to the number of times considered in the combined datasets 
#     (i.e., for the previously plus now to be considered times)
# set datetime1 above to the beginning time for the new times to be considered
setenv RESTART_DIR /discover/nobackup/projects/gmao/nwposse/develop/ODSstats/${exp_dir}/Corr_rad_June_5d/Save_sums
setenv RESTART_I 0  # set to 0 unless this is a restart

# Set instr_platform for .ods files to process
#set kx_set= (amsua_n15 amsua_n18 amsua_n19 amsua_metop-a amsua_metop-b hirs4_metop-a hirs4_n18 mhs_metop-a mhs_metop-b airs_aqua iasi_metop-a iasi_metop-b cris_fsr_n20 atms_npp atms_n20 ssmis_f17 amsr2_gcom-w1 avhrr_metop-a avhrr_n18 avhrr_n19 gmi_gpm)
set kx_set= (airs_aqua iasi_metop-a iasi_metop-b)

#
# END OF USUAL USER-DEFINED VARIABLES
#
setenv WORKDIR /discover/nobackup/$user/WORK/CORwork.$$ 
set addhours="+06"          # hours between times to consider (e.g., +06) 
set comp_diag="omf"         # compute diagnostic: "omf", "oma", "obs", or "amb"
set ibins="025"         # number of separation distance bins (3 digits)
set x_max="960.0000"    # maximum distance to consider in km (8 digits)

if ($region == 'GLOB' ) then      # Global
  set latN="  90.0"     # northern-most latitude in subdomain considered   
  set latS=" -90.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
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
else if ($region == 'RUSS' ) then  # Southern Hemisphere Extra-Tropics
  set latN="  70.0"     # northern-most latitude in subdomain considered   
  set latS="  50.0"     # southern-most latitude in subdomain considered   
  set lonW="  20.0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 160.0"     # eastern-most longitude in subdomain considered 0-360
endif

set kx="  0"       # For sat obs, set kx=0 
set kt="  0"       # not used for sat obs
set delp="   5.0"  # not used for sat obs

# Create and move to empty working directory
if (! -e $WORKDIR            ) mkdir -p $WORKDIR
cd $WORKDIR
/bin/rm -f  *

# Copy file to working directory
cp $PROGHOME/hcorr_sumods.x .
cp $PROGHOME/hcorr_tables.x .
cp $PROGHOME/read_j_value.x .
cp $ADDTIME/addtime.x .

foreach instrument ($kx_set)

  echo $instrument

  set sat_inst=`echo $instrument | cut -d_ -f 1-1`  # extract string before "_" denoting obs type

  set ods_name="${exp}.diag_${instrument}."
  set output_file="hcor_${exp}_${instrument}"  

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
    cp $RESTART_DIR/save_${output_file}_${region}_$RESTART_I.bin input_stats
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

  if ($i == 2) set file_name=xxx  #QQQQQQQQ

    if (-e $file_name) then

      @ j = $j + 1
      ./hcorr_sumods.x ${j} ${ibins} ${x_max} ${yyyy}${mm}${dd} ${hh}0000 ${latN} ${latS} ${lonW} ${lonE} ${delp} ${kt} ${kx} ${comp_diag} ${file_name} ${sat_inst} ${sat_info_file} >> correl_printout

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
  ./hcorr_tables.x input_stats $output_file >> correl_printout
  echo "Correlations program executed"
  if (! -e $OUTPUT ) mkdir -p $OUTPUT
  /bin/mv $output_file ${OUTPUT}/${output_file}_${region}.txt
  /bin/mv correl_printout ${OUTPUT}/correl_printout
  /bin/mv corr.bin ${OUTPUT}/${output_file}_${region}.bin
  echo "Correlation file copied"
#
# Save restart files (files of sums over time)
  set savedir=$OUTPUT/Save_sums
  if (! -e $savedir ) mkdir -p $savedir
  set savefile=$savedir/save_${output_file}_${region}_$i.bin
  /bin/mv input_stats $savefile
#
  echo "end consideration of instrument ${instrument}"
 
end  # loop over instruments

# --------
exit 0
# --------
