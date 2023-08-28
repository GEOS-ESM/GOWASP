#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=ENORM
#SBATCH --time=00:59:00 
#SBATCH --nodes=1 
#SBATCH --output=ENORM.%j 
#SBATCH --account=s0911 
#SBATCH --qos=debug

# script for running diff_ana_phi

source ~/GOWASP_3/GOWASP_ENV
source $GOWASP_G5MODS 
setenv CGHOME $GOWASP_PATH/AnaAdj
setenv ADDTIME $GOWASP_PATH
setenv CGWORK /discover/nobackup/$user/WORK/CEwork.$$
setenv AKBK   $GOWASP_PATH/Rcfiles/eta_akbk_g5nr_1.txt 

set ab_path=/discover/nobackup/rerrico/ABfiles
set nr_path=/discover/nobackup/rerrico/ComGrid2/GMAONR_ana
set fsens_file_in=/discover/nobackup/rerrico/Fsens/fsens.tmpl.nc4 # file used as template
set output_path=/discover/nobackup/rerrico/$CGWORK
#set output_path=/discover/nobackup/rerrico/G514osse/asens
set expid=G514osse
# set initial symoptic date time
set datetime1="2006070100"

# set number of times to process 
@ ntimes =  31

# set period (in hours) between times to process (for 6 hours, use +06 )
set period=24

# set lat, lon p range
set latN=90.
set latS=-90.
set lonW=0.
set lonE=360.      # 
set pmax=110000.   # units Pa
set pmin=0.        # units Pa
#
# norm choices are:
# twe (total moist energy) includes T,q,u,v,ps 
# txe (total dry energy) includes T,u,v,ps 
# kxe (total kinetic energy) includes u,v 
# qxe (only moist 'energy') includes q
# ape (available potential energy) includes T,ps
# tps (ps part of ape only) includes ps
# tte (total thermal available energy) includes T 
set norm=twe

# set ctest=T if test mode; F otherwise
set ctest=F

# Create working directory and move to it
if (! -e $CGWORK            ) mkdir -p $CGWORK
cd $CGWORK
/bin/rm -f  *

cp ${CGHOME}/enorm.x      prog.x
cp ${ADDTIME}/addtime.x   addtime.x
cp ${AKBK}                akbk_file
cp ${fsens_file_in}       fsens_in
# separate datetime1 into components

  set result=`addtime.x $datetime1 00`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set yy=`echo $yyyy | cut -c 3-4`

# make output directory if it does not yet exist

if (! -e $output_path ) mkdir -p $output_path

# Loop over times

@ i = 1
@ ilast =  $ntimes  
while ( $i <= $ntimes ) 
  
# x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 

# increment date and time
  if ($i == 1) then
    set add_period=0
    set datetime=$datetime1
  else
    set add_period=$period
  endif
  set result=`addtime.x $datetime $add_period`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set yy=`echo $yyyy | cut -c 3-4`
  set datetime=${yyyy}${mm}${dd}${hh}
  set datetimeA=${yyyy}${mm}${dd}_${hh}z
#
  set result=`addtime.x $datetime -9`
  set yyyy2=`echo $result | cut -d" " -f1`
  set mm2=`echo $result | cut -d" " -f2`
  set dd2=`echo $result | cut -d" " -f3`
  set hh2=`echo $result | cut -d" " -f4`
  set yy2=`echo $yyyy2 | cut -c 3-4`
  set datetimeB=${yyyy2}${mm2}${dd2}_${hh2}z
#
  set result=`addtime.x $datetime -3`
  set yyyy3=`echo $result | cut -d" " -f1`
  set mm3=`echo $result | cut -d" " -f2`
  set dd3=`echo $result | cut -d" " -f3`
  set hh3=`echo $result | cut -d" " -f4`
  set yy3=`echo $yyyy3 | cut -c 3-4`
  set datetimeC=${yyyy3}${mm3}${dd3}_${hh3}z
#
  set result=`addtime.x $datetime  6`
  set yyyy4=`echo $result | cut -d" " -f1`
  set mm4=`echo $result | cut -d" " -f2`
  set dd4=`echo $result | cut -d" " -f3`
  set hh4=`echo $result | cut -d" " -f4`
  set yy4=`echo $yyyy3 | cut -c 3-4`
  set datetimeD=${yyyy4}${mm4}${dd4}_${hh4}z
#
  set time_ana=${datetimeC}+${datetimeD}-${datetimeA}
  set time_bkg=${datetimeB}+${datetimeD}-${datetimeA}

# set input file names

 set ana_file=${ab_path}/${expid}.ana.eta.${datetimeA}.nc4
 set bkg_file=${ab_path}/${expid}.bkg.eta.${datetimeA}.nc4
 set nr_file=${nr_path}/Y${yyyy}/M${mm}/comgrid.true.eta.${datetimeA}.nc4
 set fsens_out_ana=${output_path}/${expid}.fsens_${norm}.eta.${time_ana}.nc4
 set fsens_out_bkg=${output_path}/${expid}.fsens_${norm}.eta.${time_bkg}.nc4

 if ($ctest == "T") then 

 endif

 touch file_enorms1
 touch file_prnt1

# run executable
 ./prog.x $ana_file $nr_file fsens_in ${datetime}0000 $latN $latS $lonW $lonE $pmax $pmin $norm $ctest ana
 if ($ctest == "T") then 
   cat file_prnt1 enorm_prnt_test > file_prnt2
   /bin/mv file_prnt2 file_prnt1
 else
   cp fsens_in $fsens_out_ana
 endif  
 cat file_enorms1 file_norms.txt > file_enorms2
 /bin/mv file_enorms2 file_enorms1
#
 ./prog.x $bkg_file $nr_file fsens_in ${datetime}0000 $latN $latS $lonW $lonE $pmax $pmin $norm $ctest bkg
 if ($ctest == "T") then 
   cat file_prnt1 enorm_prnt_test > file_prnt2
   /bin/mv file_prnt2 file_prnt1
 else
   cp fsens_in $fsens_out_bkg
 endif
 cat file_enorms1 file_norms.txt > file_enorms2
 /bin/mv file_enorms2 file_enorms1
 
# end loop over times
@ i++
end

#
 if ($ctest == "T") then 
   /bin/mv file_prnt1 $CGHOME/file_prnt_test.txt     
 endif
 /bin/mv file_enorms1 $CGHOME/file_norms_$norm.txt     
#
  echo "Script complete"
exit 0

