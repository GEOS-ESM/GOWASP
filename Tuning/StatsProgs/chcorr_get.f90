      program chcorr_get
!
! Get subset of data for looking at stats of H(x_b)-H(x_t)
!
      use m_ods_RE
!
      implicit none
!
      logical :: lev_OK
!
      integer, parameter :: rkinds=8
      integer, parameter :: ikinds=8
      integer, parameter :: nlook=1000000
      integer(ikinds), allocatable :: iaccept(:,:)
      integer(ikinds), allocatable :: ibins(:,:,:)
      integer :: ksave(nlook)  ! array for saving sounding index
      integer :: ib
      integer :: iks,     iks2
      integer :: i
      integer :: i_list(15)
      integer :: i_list_last
      integer :: i_list_skip
      integer :: i_values(15)
      integer :: iqc
      integer :: ifield, itype
      integer :: ks1, ks2
      integer :: n, n1, n2
      integer :: nch1, nch2
      integer :: nsave
      integer :: n_bins,  n_bins2
      integer :: n_count, n_count2
      integer :: n_ksmax, n_ksmax2
      integer :: n_date,  n_date2
      integer :: n_time,  n_time2
      integer :: n_obs
      integer :: argc
      integer :: n_channel
!
      real(rkinds), parameter :: zero=0._rkinds
      real(rkinds) :: pifac
      real(rkinds) :: earthr2
      real(rkinds) :: d, d1, d2, d3
      real(rkinds) :: latN, latS, lonW, lonE
      real(rkinds) :: latN2, latS2, lonW2, lonE2
      real(rkinds) :: ods_lev
      real(rkinds) :: x
      real(rkinds) :: xlev, xlon
      real(rkinds) :: xm1, xm2, xs1, xs2, xa, xb, xc
      real(rkinds) :: x1,    x12
      real(rkinds) :: xbin,  xbin2
      real(rkinds) :: x_values(15)
      real(rkinds) :: xsave(nlook,5)
      real(rkinds), allocatable :: xsum(:,:)
      real(rkinds), allocatable :: xsumsq(:,:)
      real(rkinds), allocatable :: xmean(:)
      real(rkinds), allocatable :: xmsqr(:)
      real(rkinds), allocatable :: xstdv(:)
      real(rkinds), allocatable :: xbins(:,:,:,:)
!
      character(len=3) :: c_ksmax
      character(len=8) :: c_date
      character(len=6) :: c_time
      character(len=3) :: c_channel
      character(len=220) :: c_file_in, c_file_out
      character(len=11), parameter :: file_in='input_stats'
      character(len=12), parameter :: file_out='output_stats'
!
      type (ods_vect) :: ods
      character (len=30) :: type
      integer :: ierr
!       
! Read and check arguments
      argc = iargc()
      if (argc .ne. 5) then
        stop
      endif
      call GetArg( 1_4, c_channel)
      call GetArg( 2_4, c_date)
      call GetArg( 3_4, c_time)
      call GetArg( 4_4, c_file_in)
      call GetArg( 5_4, c_file_out)
      read (c_date, '(i8)')   n_date
      read (c_time, '(i6)')   n_time
      read (c_channel,'(i3)') n_channel
!
! Read ods file for 1 time
      call ODS_get_RE (.true.,trim(c_file_in),n_time,n_obs,ods,ierr)
      print *,'n_obs, ierr=',n_obs,ierr
      if (n_obs == 0) then
        print *,'Job stopping since no accepted obs found on file'
        stop
      endif 
!
      open (10,file=trim(c_file_out),form='formatted')
!
! save obs in regional subset
      nsave=0  ! counts # observations in desired region, channels, or levels
      do n=1,n_obs
!
! reset longitude so that it is in the range 0-360 
        xlon=mod(ods%data%lon(n),360.)
        if (xlon < zero)  xlon=xlon+360.
        if (xlon >= 360.) xlon=zero       ! accounts for round-off of xlon
!
! check for obs level set too large (nint function below otherwise fails)
        ods_lev=min(1.e7,ods%data%lev(n)) ! 1.e7 is just some large number 
!
        lev_OK=.false.
        if ( (nsave < nlook)                    .and. & 
             (ods%data%qcexcl(n) == 0)          .and. &
             (nint(ods_lev) == n_channel)       .and. &
             (abs(ods%data%omf(n)) < 1.e8) )    then
          lev_OK=.true.
        endif
!
        if (lev_OK) then
          nsave=nsave+1
          write (10,'(i6,2f8.2,f12.4,i9,i7,i5,2f8.3)') & 
            nsave,ods%data%lat(n),ods%data%lon(n),ods%data%time(n), &
            n_date,n_time,n_channel,ods%data%obs(n),ods%data%omf(n)
        endif
!
      enddo  ! loop over all obs
!
      close (10)
!
      call ods_clean_RE (ods)
      print *,'Program end for time=',n_count,n_date,n_time
!
      end program chcorr_get
