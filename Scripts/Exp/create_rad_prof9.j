#!/bin/csh -xf 
# ------------------------------
#SBATCH --job-name=SIMRADP
#SBATCH --time=10:59:00 
#SBATCH --nodes=1   --constraint=cas 
#SBATCH --output=SIMRADP.%j 
#SBATCH --account=s0911 
##SBATCH --qos=debug

# Specify modules and number of mpi processors
# 25 min per synoptic time for all 10 rad types if 1/2 hourly NR used

setenv NBRPROCS 18  # must be on a single node
echo "NBRPROCS = ${NBRPROCS}"

# Experiment environment
# ----------------------
  source ~/GOWASP3_ENV16
  setenv SIMHOME    $GOWASP_PATH/Sim_rad    # dir for some files/executables
  setenv ADDTIME    $GOWASP_PATH/addtime.x  # addtime executable
  setenv RMSHMKEY   $GOWASP_RMSHMKEY        # dir for ESMF shmem cleaning script
  setenv SIMWORK    $NOBACKUP/WORK/SimRPwork.$$ 
# 
  setenv RC_FILE    $GOWASP_PATH/Rcfiles/field_list_rad.rc # rc file name
  setenv INPUT      /discover/nobackup/projects/gmao/nwposse/develop/OSSEobs/N010/ObsList # dir for output
  setenv OUTPUT     /discover/nobackup/projects/gmao/nwposse/develop/OSSEobs/N010/ObsProf # dir for output
#
set datetimeN=2006092712  # starting time to process
@ ntimes =  28           # number of times to process


# set which radiance data sets to create  
set instr_set= (AIRS AMSUAAQUA AMSUA HIRS4 IASI ATMS MHS GMI SSMIS CRISFSR AMSR2)
#set instr_set= (ATMS)
set instr_in1="T"   # T or F: compute all types in 1 execution per time 
set test_print="F"  # T or F: print sample output for testing purposes  
# END OF USER SET VARIABLES
# ----------------------------------------------------------------------

if ($instr_in1 == "T") then
  @ icount=0
  foreach instrument ($instr_set)
    @ icount = $icount + 1
    if ($icount == 1) then
      set instr_list=$instrument
    else 
      set instr_list="$instr_list/$instrument"
    endif
  end
else
  set instr_list= ($instr_set)
endif

# Make sure files to be written are acessible to others
umask 022

# Increment in hours between central synoptic input times
set add6h=+06

# Create working directory
if (! -e $SIMWORK            ) mkdir -p $SIMWORK

cd $SIMWORK
/bin/rm -f  *
mkdir Data_out

cp $SIMHOME/create_rad_profs.x prog.x
cp $SIMHOME/prof_io_tmpl.rc prof_io_file.tmpl
cp $SIMHOME/create_rad_reord.x reord.x
cp $ADDTIME addtime.x 

 # x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
# loop over synoptic times

@ i_files_made = 0     # counter for number of files produced
@ i_files_exist = 0    # counter for number of files produced
@ i_files_noexist = 0  # counter for number of files produced
@ i = 1                # counter for number of times processed
@ ilast =  $ntimes  
while ( $i <= $ntimes ) 
  
  if ($i == 1) then
    set addh=000
  else
    set addh=$add6h 
  endif
#
# set time for data to be created (must be same as NR time)
  set result=`./addtime.x $datetimeN $addh`
  set yyyy=`echo $result | cut -d" " -f1`
  set mm=`echo $result | cut -d" " -f2`
  set dd=`echo $result | cut -d" " -f3`
  set hh=`echo $result | cut -d" " -f4`
  set datetimeN="${yyyy}${mm}${dd}${hh}"    # NR synoptic time
  @ yy = $yyyy % 10 
  set yy=`printf "%2.2d\n" $yy`
  set cdtime="${datetimeN}0000"             # 00 min and 00 sec indicated
#
# Loop over instruments types to consider
  foreach instrument ($instr_list)
    echo 'Process instrument '  $instrument
#
# Check that all required input files exist
# The file name must be of the same template as in prof_io_file.tmpl
# 
    set some_files_exist="F"
    if ($instr_in1 == "T") then    
      @ icount=0
      foreach inst_1 ($instr_set)
        set input_dir=$INPUT/$inst_1/Y$yyyy/M$mm
        set input_file=$input_dir/obslist_$inst_1.$datetimeN.bin  
        if (-e $input_file) then
          @ i_files_exist = $i_files_exist + 1
          set some_files_exist="T"
          @ icount = $icount + 1
          if ($icount == 1) then
            set inst_2=$inst_1
          else 
            set inst_2="$inst_2/$inst_1"  # only include instruments whose file exist
          endif
        else 
          @ i_files_noexist = $i_files_noexist + 1
          echo "ERROR: nonexistent input file = ${input_file}"
        endif
      end
    else
      set input_dir=$INPUT/$instrument/Y$yyyy/M$mm
      set input_file=$input_dir/obslist_$instrument.$datetimeN.bin  
      if (-e $input_file) then
        @ i_files_exist = $i_files_exist + 1
        set some_files_exist="T"
        set inst_2=$instrument
      else 
        @ i_files_noexist = $i_files_noexist + 1
        echo "ERROR: nonexistent input file = ${input_file}"
      endif
    endif      

    if ($some_files_exist == "T") then  # some files exist
#
# Create profiles of obs for 1 data time and requested instrument(s)
# 
      $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
      mpirun -np $NBRPROCS ./prog.x  $inst_2 $cdtime $INPUT $RC_FILE Data_out $test_print
      $RMSHMKEY/RmShmKeys  # clean up any residual shared memory
#
# Reorder data in the output profile files
      foreach file1 (Data_out/*)     
        if (-e $file1) then # output file can be reordered
          set cut_file1=`echo $file1  | cut -c 10- ` 
          set cut_file2=`echo $cut_file1 | cut -d "." -f 1 `
          set output_dir="$OUTPUT/$cut_file2/Y$yyyy/M$mm"
          ./reord.x $file1 data_reordered    
          if (-e data_reordered) then # reordered file created and can be saved
            set output_file="$output_dir/profile_$cut_file1"  
            if (! -e $output_dir) mkdir -p $output_dir
            /bin/mv data_reordered $output_file
            if (-e $output_file) then 
              @ i_files_made = $i_files_made + 1
              echo 'data reordered and moved'
            endif 
          else
            echo 'NO DATA REORDERED'     
          endif
        else
          echo 'NO DATA OUTPUT FROM CREATE_RAD_PROF'     
        endif
      end  # loop over files in Data_out
#
    else
      echo 'NO INPUT FILES OF REQUESTED TYPE(s) EXIST FOR TIME = ' $cdtime
      echo 'Skip processing'
    endif
    /bin/rm -f Data_out/*
#  
  end  # loop over instrument
#
  @ i++
end    # loop over times  
# All done
echo "Number of reordered files made = $i_files_made"
echo "Number of input files existing = $i_files_exist"
echo "Number of input files not-existing = $i_files_noexist"
echo "end script"
date

# --------
  exit 0
# --------
