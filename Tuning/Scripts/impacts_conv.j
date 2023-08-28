#!/bin/csh -xf
# ------------------------------
# Execute countobs_sum.x and countobs_tables.x to create files of tables 
# of means, standard deviations, and counts of O-F for all radiance obs 
# for a selected region. Also, compute a distribution with latitude.
#
# If new instrument types are added, the number of channels must 
# be specified at the appropriate location in the script.
#

source ~/GOWASP_3/GOWASP_ENV

setenv PROGHOME $GOWASP_PATH/Tuning/StatsProgs  # home for executables
setenv ADDTIME  $GOWASP_PATH                        # home for addtime.x
set eof_file=$GOWASP_PATH/Tuning/Scripts/impacts_eof.txt

set exp="G514osse"                          # experiment name
set exp_dir=$exp                         # experiment subdir    
#set ods_path=/discover/nobackup/rerrico/FsensProg # path to ods files
set ods_path=/discover/nobackup/rerrico/ODSfiles/${exp_dir}/Impact0 # path to ods files

setenv OUTPUT /discover/nobackup/rerrico/ODSstats/${exp_dir}/Impact0 # output dir 
setenv WORKDIR /discover/nobackup/$user/WORK/IMPwork.$$ 
set archive="F"  # F if not archive directory, T otherwise

set datetime_bkg="2006070115"  # starting datetime for  bkg yyyymmddhh          
set datetime_fcs="2006070206"  # starting datetime for fcst yyyymmddhh 
set datetime_ana="2006070200"  # starting datetime for  ana yyyymmddhh
@ ntimes =   30              # number of synoptic times to process            
set addhours="+24"          # hours between times to consider  
set region="GLOB"           # geo region to compute correls (see below)    
set sfc_type="BOTH"         # SEA, LAND, or BOTH
set impact_metric="imp3_twe"

# END OF USUAL USER-DEFINED VARIABLES

if ($region == 'GLOB' ) then      # Global
  set latN="  90.0"     # northern-most latitude in subdomain considered   
  set latS=" -90.0"     # southern-most latitude in subdomain considered   
  set lonW="    .0"     # western-most longitude in subdomain considered 0-360
  set lonE=" 360.0"     # eastern-most longitude in subdomain considered 0-360
endif 

# Create and move to empty working directory
if (! -e $WORKDIR            ) mkdir -p $WORKDIR
cd $WORKDIR
/bin/rm -f  *

# Copy file to working directory
cp $PROGHOME/impacts_conv.x prog.x
cp $ADDTIME/addtime.x .

  set instrument=conv
  set ods_name="${exp}.${impact_metric}_${instrument}."
  set file_out="impacts_${exp}_${instrument}"  
  touch file_cat1

# determine first set of date and time variables
  set result=`./addtime.x $datetime_bkg +00`
  set yyyy1=`echo $result | cut -d" " -f1`
  set mm1=`echo $result | cut -d" " -f2`
  set dd1=`echo $result | cut -d" " -f3`
  set hh1=`echo $result | cut -d" " -f4`
  set datetime1="${yyyy1}${mm1}${dd1}${hh1}"

  set result=`./addtime.x $datetime_fcs +00`
  set yyyy2=`echo $result | cut -d" " -f1`
  set mm2=`echo $result | cut -d" " -f2`
  set dd2=`echo $result | cut -d" " -f3`
  set hh2=`echo $result | cut -d" " -f4`
  set datetime2="${yyyy2}${mm2}${dd2}${hh2}"

  set result=`./addtime.x $datetime_ana +00`
  set yyyy3=`echo $result | cut -d" " -f1`
  set mm3=`echo $result | cut -d" " -f2`
  set dd3=`echo $result | cut -d" " -f3`
  set hh3=`echo $result | cut -d" " -f4`
  set datetime3="${yyyy3}${mm3}${dd3}${hh3}"

  rm -f input_stats
  rm -f correl_printout

# initialize counters
  @ i = 1
  @ j = 0
  @ i_no = 0
  @ ilast =  $ntimes  

# loop over times 
  while ( $i <= $ntimes ) 

    echo "${datetime3}"  
    set file_dir_in=${ods_path}/Y${yyyy1}/M${mm1}
    set file_time1=${yyyy1}${mm1}${dd1}_${hh1}z
    set file_time2=${yyyy2}${mm2}${dd2}_${hh2}z
    set file_time3=${yyyy3}${mm3}${dd3}_${hh3}z
    set file_tail=obs.${file_time1}+${file_time2}-${file_time3}.ods
    set file_name=${file_dir_in}/${ods_name}${file_tail} 

    if (-s $file_name) then

      @ j = $j + 1
      ./prog.x ${yyyy3}${mm3}${dd3} ${hh3}0000 ${latN} ${latS} ${lonW} ${lonE} $sfc_type ${file_name} >> impact_printout

      if (-e output_stats) then
        cat file_cat1 output_stats > file_cat2
        /bin/rm/ file_cat1 output_stats 
        mv file_cat2 file_cat1 
      else 
        @ j = $j - 1
      endif

    else
      @ i_no = $i_no + 1
      echo "FILE DOES NOT EXIST for datetime ${datetime3}"
    endif

  set result=`./addtime.x $datetime1 $addhours`
  set yyyy1=`echo $result | cut -d" " -f1`
  set mm1=`echo $result | cut -d" " -f2`
  set dd1=`echo $result | cut -d" " -f3`
  set hh1=`echo $result | cut -d" " -f4`
  set datetime1="${yyyy1}${mm1}${dd1}${hh1}"

  set result=`./addtime.x $datetime2 $addhours`
  set yyyy2=`echo $result | cut -d" " -f1`
  set mm2=`echo $result | cut -d" " -f2`
  set dd2=`echo $result | cut -d" " -f3`
  set hh2=`echo $result | cut -d" " -f4`
  set datetime2="${yyyy2}${mm2}${dd2}${hh2}"

  set result=`./addtime.x $datetime3 $addhours`
  set yyyy3=`echo $result | cut -d" " -f1`
  set mm3=`echo $result | cut -d" " -f2`
  set dd3=`echo $result | cut -d" " -f3`
  set hh3=`echo $result | cut -d" " -f4`
  set datetime3="${yyyy3}${mm3}${dd3}${hh3}"

# end loop

    @ i++
  end

  @ i = $i - 1
  echo "${i} times considered"
  echo "${j} times processed"
  echo "${i_no} = number of files not existing"

  if (! -e $OUTPUT ) mkdir -p $OUTPUT
  cat file_cat1 $eof_file > file_cat2
  /bin/mv file_cat2 ${OUTPUT}/${file_out}.txt
  /bin/mv impact_printout ${OUTPUT}/impacts_conv_printout
  echo "Impact file moved"
  echo "end consideration of instrument ${instrument}"
 
date
# --------
exit 0
# --------
