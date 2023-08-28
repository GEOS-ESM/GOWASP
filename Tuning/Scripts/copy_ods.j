#!/bin/csh -fx

# Copy ods files from several archive directories to a single nobackup 
# directory

source ~/GOWASP_3/GOWASP_ENV
setenv addtime  $GOWASP_PATH/addtime.x  # home for addtime.x

# the user should specify input and ouput paths, datatime1, and ntimes 

set input_path=/discover/nobackup/projects/gmao/nwposse/develop/allsky/obs
set output_path=/discover/nobackup/rerrico/ODSfiles/allsky


set datetime1="2006062200"    # starting datetime
@ ntimes = 46  # number of synoptic times to process
set addhours="+06" 

if (! -e $output_path ) mkdir -p $output_path

# determine first set of date and time variables
    set result=`$addtime $datetime1 +00`
    set yyyy=`echo $result | cut -d" " -f1`
    set mm=`echo $result | cut -d" " -f2`
    set dd=`echo $result | cut -d" " -f3`
    set hh=`echo $result | cut -d" " -f4`
    set datetime="${yyyy}${mm}${dd}${hh}"

# initialize counters
    @ i = 1

# loop over times 
    while ( $i <= $ntimes ) 
      set file_in="${input_path}/Y${yyyy}/M${mm}/D${dd}/H${hh}/*.ods"
      set output_dir="${output_path}/Y${yyyy}/M${mm}"

      dmget $file_in

      if (! -e $output_dir ) mkdir -p $output_dir
      cp ${file_in} ${output_dir}/.

      set result=`$addtime $datetime $addhours`
      set yyyy=`echo $result | cut -d" " -f1`
      set mm=`echo $result | cut -d" " -f2`
      set dd=`echo $result | cut -d" " -f3`
      set hh=`echo $result | cut -d" " -f4`
      set datetime="${yyyy}${mm}${dd}${hh}"

      @ i++
    end

exit 0
