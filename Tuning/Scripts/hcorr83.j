#!/bin/csh -fx 
# ------------------------------
#SBATCH --job-name=HCCONV
#SBATCH --time=03:59:00 
#SBATCH --nodes=1 --constraint=hasw  
#SBATCH --output=HCCONV.%j 
#SBATCH --account=s0911 
##SBATCH --qos=debug
#
# Execute hcorr_sumods.x and hcorr_tables.x to create files of tables of means, 
# standard deviations, and horizontal correlations of O-F for conv obs for 
# a given region. 

# Set user variables. For every run, 
source ~/GOWASP_3/GOWASP_ENV
setenv PROGHOME $GOWASP_PATH/Tuning/StatsProgs # home for executables
setenv ADDTIME  $GOWASP_PATH                   # home for addtime.x
#
# Set the part of input path and .ods file name that is the experiment name
# $exp is the exp name that appears as part of the GSI ods file names
# $exp_dir is part of path for both input and output that varies with exp
# It is convienent to make $exp_dir=$exp, but they need not be, e.g., 
# as when multiple experiments use the same exp name. 
# It is assumed that the input file path will have sub-directories /Yyyyy/Mmm
# proceeding the ods_path (see file_dir_in below).
#
set exp="d5124_m2_jan79"
set exp=d517_osse_83
set exp_dir=${exp}
set ods_path="/discover/nobackup/rerrico/ODSfiles/${exp_dir}" 
#
# Set datetime (remember real and osse have different years) and region
set datetime1="2006062300"  # starting datetime yyyymmddhh                  
@ ntimes_other =   28       # number of synoptic times to process            
@ ntimes_raobs =   14       # number of synoptic times to process for raobs  
set region="GLOB"           # geo region to compute correls (see below)      

# Set the output directory
setenv OUTPUT /discover/nobackup/rerrico/ODSstats/${exp_dir}/Corr_conv_June
#
# Set kx of obs to process
set kx_set= (120 130 131 132 180 187 220 221 223 224 229 230 231 280)  


#
# END OF USUAL USER-DEFINED VARIABLES
#
setenv WORKDIR /discover/nobackup/$user/WORK/CORwork.$$ 
set ods_name="${exp}.diag_conv."     # This is general 
set comp_diag="omf"         # compute diagnostic: "omf", "oma", "obs", or "amb"
set addhours_other="+06"    # hours between times to consider (e.g., +06)    
set addhours_raobs="+12"    # hours between times to consider for raobs      

set ibins="025"         # number of bins (3 digits)
set x_max="960.0000"    # maximum distance to consider in km (8 digits)
set delp_r="2."      # consider pressures within +/- range of mb for raobs
set delp_o="20."     # consider pressures within +/- range of mb for other

if ($region == 'GLOB' ) then      # Global
  set latN="  90.0"     # northern-most latitude in subdomain considered   
  set latS=" -90.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
else if ($region == 'LB02' ) then      # Western Hemisphere
  set latN=" -45.0"     # northern-most latitude in subdomain considered   
  set latS=" -75.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
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
cp $PROGHOME/hcorr_sumods.x hcorr_sumods.x  
cp $PROGHOME/hcorr_tables.x .
cp $ADDTIME/addtime.x .

foreach kx ($kx_set)

  set kx1=`echo $kx | cut -c 1-1`  # extract first character of kx
  if ($kx1 == "1") then
    if ($kx == "120" || $kx == "180" ) then               
      set field_set=( T q r)      # temperature specific and relative humidity 
    else
      set field_set=( T)
    endif
  else 
    set field_set=( u v)
  endif

  if ($kx == "120" || $kx == "220") then 
    set addhours=$addhours_raobs
    set ntimes=$ntimes_raobs
    set delp=$delp_r
  else
    set addhours="$addhours_other"
    set ntimes="$ntimes_other"
    set delp=$delp_o
  endif 

  foreach field ($field_set)

    set output_file="hzcor_${exp}_conv_${kx}_${field}" 
    if ( $field == "T" ) then 
      set kt="44" 
    else if ( $field == "q" ) then 
      set kt="11"
    else if ( $field == "r" ) then 
      set kt="10"
    else if ( $field == "u" ) then 
      set kt=" 4"
    else if ( $field == "v" ) then 
      set kt=" 5"
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
      set file_dir_in="${ods_path}/Y${yyyy}/M${mm}"
      set file_name="${file_dir_in}/${ods_name}${yyyy}${mm}${dd}_${hh}z.ods" 

      if (-e $file_name) then

        @ j = $j + 1
        ./hcorr_sumods.x ${j} ${ibins} ${x_max} ${yyyy}${mm}${dd} ${hh}0000 ${latN} ${latS} ${lonW} ${lonE} ${delp} ${kt} ${kx} ${comp_diag} ${file_name} conv none >> correl_printout
    /bin/mv correl_printout ${OUTPUT}/correl_printout.xxx
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

    ./hcorr_tables.x input_stats $output_file >> correl_printout
    echo "Correlations program executed"

    if (! -e $OUTPUT ) mkdir -p $OUTPUT
    /bin/mv $output_file ${OUTPUT}/${output_file}_${region}.txt
    /bin/mv correl_printout ${OUTPUT}/correl_printout.txt
    /bin/mv corr.bin ${OUTPUT}/${output_file}_${region}.bin
    echo "Correlation file copied"
    echo "end script"
 
  end  # loop over kt
end  # loop over kx

# --------
  exit 0
# --------
