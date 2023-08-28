#!/bin/csh -xf 

# cat obs files for CONV and SATWIND to create PREPBUFR obs

set yy=06     # last 2 digits of year
set mm=06    # 2-digit month

#set days_set = (01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31)   # 2 digit days

set days_set = (20 21 22 23 24 25 26 27 28 29 30)   

set hours_set = (00 06 12 18)     # 2 digit hours


set dir_sat=/discover/nobackup/projects/gmao/nwposse/develop/OSSEobs/N010
set dir_conv=/discover/nobackup/rerrico/OSSEobs/N010

set dir_out=$dir_sat/PREPBUFR/Y20$yy/M$mm
if (! -e $dir_out ) mkdir -p $dir_out

foreach dd ($days_set)
  foreach hh ($hours_set)

    set dtime=$yy$mm$dd.t${hh}z
    set satwind_file=$dir_sat/SATWIND/Y20$yy/M$mm/satwind.$dtime.prepbufr
    set windnosat_file=$dir_conv/CONV/Y20$yy/M$mm/gdas1.$dtime.prepbufr
    set windboth_file=$dir_out/gdas1.$dtime.prepbufr

    cat $windnosat_file $satwind_file > $windboth_file

  end
end
exit
