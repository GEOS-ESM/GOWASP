#!/bin/csh -fx

# dmget ana and bkg files valid at same time for several synoptic times

cd /archive/u/rerrico/G514osse/ana

set addtime="/discover/home/rerrico/GOWASP_3/addtime.x"

# set expid name for GDAS prog files
set exp_id="G514osse"

# set initial symoptic date time
set datetime1="2006070100"   

# set number of times to process 
@ ntimes =   28

# set period (in hours) between times to process (for 6 hours, use +06 )
set period=+12


  set result=`$addtime $datetime1 +00`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  @ yy = $yyyy % 10 
  set yy=`printf "%2.2d\n" $yy`

set flist=""

# Loop over times

@ i = 1
@ ilast =  $ntimes  
while ( $i <= $ntimes ) 

# increment date and time
  if ($i == 1) then
    set add_period=0
    set datetime=$datetime1
  else
    set add_period=$period
  endif
 echo 'datetime=' ${datetime}
 echo 'add_period=' ${add_period}
  echo `addtime_2.x $datetime $add_period`
  set result=`$addtime $datetime $add_period`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set datetime="${yyyy}${mm}${dd}${hh}"
  @ yy = $yyyy % 10 
  set yy=`printf "%2.2d\n" $yy`
  set datetime2="${yyyy}${mm}${dd}"_"${hh}z"

# set input file names

 set flist="${flist} Y${yyyy}/M${mm}/${exp_id}.bkg.eta.${datetime2}.nc4 Y${yyyy}/M${mm}/${exp_id}.ana.eta.${datetime2}.nc4"  

echo $flist

# end loop over times
@ i++
end

dmget $flist


  echo " "
  echo "Script complete"
  date

exit 0
