   program common_grid
!
! Take GEOS-5 NR low resolution 720x361 data sets and place them 
! on the common d-grid at 576x361 resolution. 
!
   use m_sft, only : sft_setup
   use m_sft, only : sft_apply
   use m_write_nc4, only : write_nc4_setup
   use m_write_nc4, only : write_nc4_2dfld
   use m_write_nc4, only : write_nc4_3dfld
   use m_write_nc4, only : write_nc4_close
!
   implicit none
!
   integer, parameter :: nlonNR=720
   integer, parameter :: nlonCG=576
   integer, parameter :: nlats=361
   integer, parameter :: nlevs=72
   integer, parameter :: nfields2d=2
   integer, parameter :: nfields3d=4
   integer, parameter :: nfzs=1
   integer, parameter :: nfps=2
   integer, parameter :: nfields=nfields2d+nfields3d
   integer :: ier
   integer :: nf,nx
   integer :: argc
   integer(4) :: iargc
!  
   real(4) :: fieldNR(nlonNR,nlats,nlevs)
   real(4) :: fields2d(nlonCG,nlats,nfields2d)
   real(4) :: fields3d(nlonCG,nlats,nlevs,nfields3d)
!
   character(len=240) :: common_path
   character(len=240) :: file_name_2d, file_name_3d, file_name_phis
   character(len=240) :: file_name
   character(len=240) :: file1
   character(len=240) :: file_name_out
   character(len=8)   :: field_name(nfields)
   character(len=14)  :: cdatetime    ! yyyymmddhhmmss
!
! Get program arguments
   argc=iargc()
   if (argc /= 2) then
     print *,' usage must be: prog.x cdtime fileout'
     stop
   endif
   call GetArg( 1_4, cdatetime)
   call GetArg( 2_4, file_name_out)
!    
! Get file names
  common_path='/discover/nobackup/projects/gmao/osse2/stage/c1440_NR_BETA9/DATA/0.5000_deg'
  file_name_2d='/inst/inst01hr_2d_met1_Cx/Y$yyyy#/M$mm#/D$dd#/c1440_NR.inst01hr_2d_met1_Cx.$yyyymmdd_hhmm#z.nc4'
  file_name_3d='/inst/inst01hr_3d_$ffff#_Cv/Y$yyyy#/M$mm#/D$dd#/c1440_NR.inst01hr_3d_$ffff#_Cv.$yyyymmdd_hhmm#z.nc4'    
  file_name_phis='/const/const_2d_asm_Cx/Y$yyyy#/M$mm#/D$dd#/c1440_NR.const_2d_asm_Cx.$yyyymmdd#.nc4'
  field_name(1)='PHIS'
  field_name(2)='PS'
  field_name(3)='T'
  field_name(4)='QV'
  field_name(5)='U'
  field_name(6)='V'
!
! Setup Fourier transform coefficients
  call sft_setup (nlonNR,nlonCG)
!
! Loop over input 2d fields
  do nf=1,nfields2d
    file1=file_name_2d
    if (trim(field_name(nf)) == 'PHIS') then
      file1=file_name_phis
    endif
    call set_field_file_name (field_name(nf),file1,common_path, &
                              cdatetime,file_name,ier)
    call read_nc4_2dfield (nlonNR,nlats,file_name,field_name(nf), &
                  fieldNR(:,:,1),.true.,.true.,ier)
    call sft_apply (nlats,1,fieldNR(:,:,1:1),fields2d(:,:,nf:nf)) 
  enddo
!
! Loop over input 3d fields 
  do nf=1,nfields3d
    nx=nf+nfields2d
    call set_field_file_name (field_name(nx),file_name_3d,common_path, &
                              cdatetime,file_name,ier)
    call read_nc4_3dfield (nlonNR,nlats,nlevs,file_name,field_name(nx), &
                  fieldNR,.true.,.true.,ier)  
    call sft_apply (nlats,nlevs,fieldNR,fields3d(:,:,:,nf))
  enddo
!
! Output fields
  call write_nc4_setup (nlonCG,nlats,nlevs,file_name_out)
  do nf=1,nfields2d
    call write_nc4_2dfld (nlonCG,nlats,fields2d(:,:,nf),field_name(nf))
  enddo
  do nf=1,nfields3d
    nx=nf+nfields2d
    call write_nc4_3dfld (nlonCG,nlats,nlevs,fields3d(:,:,:,nf),field_name(nx))
  enddo
  call write_nc4_close
!
  print *,'file_out written: ',file_name_out
  print *,'Program Done'
!
  end program common_grid
