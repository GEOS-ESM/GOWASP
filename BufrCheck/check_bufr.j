#!/bin/csh -xf 

source ~/GOWASP3_ENV16
source $GOWASP_G5MODS
setenv SIMPROGS   $GOWASP_PATH/BufrCheck          # dir for executables 
setenv JCSDA_DATA $GOWASP_CRTM_COEF_DIR    # dir for CRTM files for IASI
setenv AIRS_BUFR_TABLE $GOWASP_PATH/airs_bufr_table  # file for airs bufr table 

set dtype='GPSRO'  

#set defaults 
set datetime='2006070100'
set ctest='F'
set bufr_table='none'
set scoef_file='none'

if ($dtype == 'AMSUA') then
set BF_file_1='/discover/nobackup/rerrico/OSSEobs/T015/AMSUA/Y2006/M06/gdas1.060610.t00z.1bamua.tm00.bufr_d'
set BF_file_2='/discover/nobackup/rerrico/OSSEobs/TE15/AMSUA/Y2006/M06/gdas1.060610.t00z.1bamua.tm00.bufr_d'
endif

if ($dtype == 'GENRADTXT') then
set BF_file_1='/discover/nobackup/rerrico/OSSEobs/T015/GENRADTXT/Y2006/M06/gdas1.20060610.t00z.genrad.txt'
set BF_file_2='/discover/nobackup/rerrico/OSSEobs/TE15/GENRADTXT/Y2006/M06/gdas1.20060610.t00z.genrad.txt'
endif


if ($dtype == 'PREPBUFR') then 
  set BF_file_1="/discover/nobackup/rerrico/OSSEobs/TEST1/PREPBUFR/Y2005/M07/gdas1.050701.t00z.prepbufr"
  set BF_file_2="/discover/nobackup/rerrico/OSSEobs/TESTERR2/PREPBUFR/Y2005/M07/gdas1.050701.t00z.prepbufr"
endif

if ($dtype == 'AIRS') then 
set BF_file_1='/discover/nobackup/rerrico/OSSEobs/T006A/AIRS/Y2006/M08/airsbufr_disc.20060801.t00z.bufr'
set BF_file_2='/discover/nobackup/rerrico/OSSEobs/T006B/AIRS/Y2006/M08/airsbufr_disc.20060801.t00z.bufr'
  set bufr_table=$AIRS_BUFR_TABLE
endif

if ($dtype == 'GPSRO') then 
  set BF_file_1='/discover/nobackup/projects/gmao/nwposse/OSSEobs/ROPP1D/GPSRO/Y2006/M07/gdas1.060701.t00z.gpsro.tm00.bufr_d'
  set BF_file_2='/discover/nobackup/projects/gmao/nwposse/OSSEobs/GSIMOD/GPSRO/Y2006/M07/gdas1.060701.t00z.gpsro.tm00.bufr_d'

endif


$SIMPROGS/check_bufr.x $dtype $datetime $BF_file_1 $BF_file_2 $bufr_table $scoef_file $ctest > CHECKBUF.out

exit
