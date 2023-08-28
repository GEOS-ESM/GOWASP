#!/bin/csh -xf
# ------------------------------
#SBATCH --job-name=COUNTC
#SBATCH --time=00:10:00 
#SBATCH --nodes=1 --constraint=hasw  
#SBATCH --output=COUNTC.%j 
#SBATCH --account=s0911 
#SBATCH --qos=debug
# ------------------------------
# Execute countobs_sum.x and countobs_tables.x to create files of tables 
# of means, standard deviations, and counts of O-F for conv obs for selected  
# kx values and region. 
#
# If new instrument types are added, the number of channels must 
# be specified at the appropriate location in the script.
#
source ~/GOWASP_3/GOWASP_ENV

setenv PROGHOME $GOWASP_PATH/Tuning/StatsProgs  # home for executables
setenv ADDTIME  $GOWASP_PATH                        # home for addtime.x
#
# Set the part of input path and .ods file name that is the experiment name
# $exp is the exp name that appears as part of the GSI ods file names
# $exp_dir is part of path for both input and output that varies with exp
# It is convienent to make $exp_dir=$exp, but they need not be, e.g., 
# as when multiple experiments use the same exp name. 
# It is assumed that the input file path will have sub-directories /Yyyyy/Mmm
# proceeding the ods_path (see file_dir_in below).
#
set exp=517test
set exp_dir=AIRT       
set ods_path="/discover/nobackup/rerrico/ODSfiles/Balance/${exp_dir}" 

# Set datetime (remember real and osse have different years) and region
set datetime1="2015070800"  # starting datetime yyyymmddhh                   
@ ntimes_other =   1      # number of synoptic times to process            
@ ntimes_raobs =   1     # number of synoptic times to process for raobs  

# Set the output directory
setenv OUTPUT /discover/nobackup/rerrico/ODSstats/Balance/${exp_dir}/Counts_conv_08

# Set kx values to process (=89 for GPSRO) 
set kx_set= (120 130 131 132 133 180 181 182 187 220 221 223 224 229 230 231 232 233 242 243 245 246 250 252 253 254 257 258 259 280 282 290 89)

#
# END OF USUAL USER-DEFINED VARIABLES
#
setenv WORKDIR /discover/nobackup/$user/WORK/CNTwork.$$ 
set addhours_other="+06"    # hours between times to consider (e.g., +06) 
set addhours_raobs="+12"    # hours between times to consider for raobs   
set region="NAMR"           # geo region to compute correls (see below)   
set sfc_type="BOTH"         # LAND, SEA, or BOTH

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
  set latS=" -90.0"     # southern-most latitude in subdomain considered   
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

foreach kx ($kx_set)

  if ($kx == "89") then 
    set field_set=(Z)   # for GPS
  else  
    set kx1=`echo $kx | cut -c 1-1`  # extract first character of kx
    if ($kx1 == "1") then
      if ($kx == "120" || $kx == "180" || $kx == "182"  ) then               
        set field_set=(T q p)      # temperature specific and relative humidity 
      else if ($kx == "181"  ) then
        set field_set=(p)
      else               
        set field_set=(T)
      endif
    else 
      set field_set=(u v)
    endif
  endif 

  if ($kx == "120" || $kx == "220") then 
    set addhours=$addhours_raobs
    set ntimes=$ntimes_raobs
  else
    set addhours="$addhours_other"
    set ntimes="$ntimes_other"
  endif 

  foreach field ($field_set)

    set ods_name="${exp}.diag_conv."
    set output_file="count_${exp}_conv_${kx}_${field}"  
    if ( $field == "T" ) then 
      set kt="44" 
    else if ( $field == "p" ) then 
      set kt="33"
    else if ( $field == "q" ) then 
      set kt="11"
    else if ( $field == "r" ) then 
      set kt="10"
    else if ( $field == "u" ) then 
      set kt=" 4"
    else if ( $field == "v" ) then 
      set kt=" 5"
    else if ( $field == "Z" ) then # GPS
      set kt="89"
    endif

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

      set file_dir_in="${ods_path}"
      set file_name="${file_dir_in}/${ods_name}${yyyy}${mm}${dd}_${hh}z.ods" 

      if (-e $file_name) then

        @ j = $j + 1
        ./prog_sum.x ${j} ${yyyy}${mm}${dd} ${hh}0000 ${latN} ${latS} ${lonW} ${lonE} ${kt} ${kx} ${sfc_type} ${file_name} conv none >> count_printout

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

# end time loop

      @ i++
    end

    @ i = $i - 1
    echo "${i} times considered"
    echo "${j} times processed"
    echo "${i_no} = number of files not existing"
    ./prog_tables.x input_stats $output_file  >> count_printout
    echo "Counting program executed"
    if (! -e $OUTPUT ) mkdir -p $OUTPUT
    /bin/mv ${output_file}.txt ${OUTPUT}/${output_file}_${region}.txt
    /bin/mv count_printout ${OUTPUT}/count_printout
    echo "Count table  file moved"
    echo "End consideration of kt= ${kt} kx=${kx}"
 
  end  # loop over field types
end  # loop over obs types

date
# --------
exit 0
# --------
