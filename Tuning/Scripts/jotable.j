#!/bin/csh -xf

Script for computing Jo table for manuscript

set PROGHOME="${GOSS_PATH}/src/Aux_progs/ODS_Progs/Progs" # home for executables
set ADDTIME="${GOSS_PATH}/src/goss"                       # home for addtime.x
set RCFILE="${GOSS_PATH}/src/Aux_progs/ODS_Progs/Obs_Impacts/impacts_multi_e30.rc" # describes obs setss
set archive="F"  # F if not archive directory, T otherwise (use dmget 1st)
set exp="real591"                                   # experiment name      &&
#set exp="g571dctl"                                   # experiment name      &&
set ods_path="/discover/nobackup/rerrico/ODS2/${exp}" # path to ods files    &&
setenv OUTPUT /discover/nobackup/rerrico/OSSEpp/${exp}/Jo_july22to31 # output dir &&
setenv WORKDIR /discover/nobackup/$user/JOTwork.$$ 

set datetime1="2011072200"  # starting datetime yyyymmddhh                   &&
@ ntimes =     40            # number of synoptic times to process          &&
set addhours="+06"          # hours between times to consider (e.g., +06) 
set lat2 = +90
set lat1 = -90

if (! -e $OUTPUT ) mkdir -p $OUTPUT

# Create and move to empty working directory
if (! -e $WORKDIR ) mkdir -p $WORKDIR
cd $WORKDIR
/bin/rm -f  *

# Copy file to working directory
cp $PROGHOME/jo_tables.x .
cp $ADDTIME/addtime.x    .
cp $RCFILE file.rc

# determine first set of date and time variables
  set result=`./addtime.x $datetime1 +00`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set datetime="${yyyy}${mm}${dd}${hh}"

# loop over times 

  @ i = 1 
  while ( $i <= $ntimes ) 

    echo "${datetime}"  

    if ($archive == "T") then
      set file_dir_in="${ods_path}/Y${yyyy}/M${mm}/D${dd}/H${hh}"
    else
      set file_dir_in="${ods_path}/Y${yyyy}/M${mm}"
    endif
    set files_in="${file_dir_in}/*.${yyyy}${mm}${dd}_${hh}z.ods" 

    \rm -f odsout
    echo "Processing odsfiles in: " $file_dir_in
    jo_tables.x -rc file.rc -lat ${lat1}:${lat2} $files_in > odsout

#       Create stand-alone accumalted summaries for sum and numneg attributes
#       ---------------------------------------------------------------------
    if ( -e odsout ) then
      
      grep "${yyyy}${mm}${dd} ${hh}0000" odsout > accum_stats
      cat accum_stats | cut -c 1-60 >> $OUTPUT/accum.sum_odsstats.txt
      \rm odsout accum_stats 
    else
      echo "STOP... odsstats output file not found"
      exit
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
#
    exit(0)
