      program channel_sumods
!
! Augment pre-computed sums of ods data for previous times by adding 
! contributions for a new time for subsequent computation of satellite 
! channel correlations and covariances.
!
      use m_ods_RE
      use m_sat_info_table, only : sat_info_table_read
      use m_sat_info_table, only : sat_info_table_get_1i
!
      implicit none
!
      logical :: lev_OK
      logical :: lon_range1, lon_range2
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
      character(len=3) :: c_comp, c_comp2  ! "omf","oma","amb", or "obs"
      character(len=3) :: c_count
      character(len=3) :: c_nbins
      character(len=8) :: c_xmax
      character(len=6) :: c_latN, c_latS, c_lonW, c_lonE
      character(len=6) :: c_delp
      character(len=3) :: c_kt, c_kx
      character(len=8) :: c_date
      character(len=6) :: c_time
      character(len=20)  :: c_sat_inst
      character(len=220) :: c_sat_info_file
      character(len=220) :: c_file, c_file2
      character(len=11), parameter :: file_in='input_stats'
      character(len=12), parameter :: file_out='output_stats'
!
      type (ods_vect) :: ods
      character (len=30) :: type
      integer :: ierr
!       
! Read and check arguments
      argc = iargc()
      if (argc .ne. 11) then
        print *,' usage must be: channel_sumods.x ncount date time ', &
                ' latN latS lonW lonE c_comp filename sat_inst sat_info_file'
        stop
      endif
      call GetArg( 1_4, c_count)
      call GetArg( 2_4, c_date)
      call GetArg( 3_4, c_time)
      call GetArg( 4_4, c_latN)
      call GetArg( 5_4, c_latS)
      call GetArg( 6_4, c_lonW)
      call GetArg( 7_4, c_lonE)
      call GetArg( 8_4, c_comp)
      call GetArg( 9_4, c_file)
      call GetArg(10_4, c_sat_inst)
      call GetArg(11_4, c_sat_info_file)
!
      read (c_count,'(i3)')   n_count
      read (c_date, '(i8)')   n_date
      read (c_time, '(i6)')   n_time
      read (c_latN, '(f6.1)') latN  
      read (c_latS, '(f6.1)') latS  
      read (c_lonW, '(f6.1)') lonW
      read (c_lonE, '(f6.1)') lonE
!
      call sat_info_table_read (c_sat_info_file,.true.,ierr)
      if (ierr /= 0) stop
      call sat_info_table_get_1i ('instr','nchan',c_sat_inst,n_ksmax,ierr)
      if (ierr /= 0) stop
!
      if (lonW < 0.) then
        lonW=lonW+360.
      endif
      if (lonE < 0.) then
        lonE=lonE+360.
      endif
      if (lonW < lonE) then
        lon_range1=.true.
      else
        lon_range1=.false.
      endif
!
      print *,' '
      print ('(2i5,2i10)'),n_count,n_ksmax,n_date,n_time
      print ('(4f8.1,2x,a3)'),latN,latS,lonW,lonE,c_comp
      print ('(a)'),trim(c_file)
!
! Check that lats requested are OK
      if (latN <= latS) then
        print *,'ERROR: latN=',latN,' must be > latS=',latS
        stop
      endif
!
! Check that c_comp has acceptable value
      if ( (c_comp /= 'omf') .and. (c_comp /= 'oma') .and. &
           (c_comp /= 'amb') .and. (c_comp /= 'obs') ) then
        print *,'ERROR: c_comp is not one of omf, oma, amb, or obs'
        print *,'       c_comp=',c_comp
        stop
      endif
!
! Read ods file for 1 time
      call ods_get_RE (.true.,trim(c_file),n_time,n_obs,ods,ierr)
      if (n_obs == 0 .or. ierr /= 0) then
        print *,'Job stopping since either there are no obs on ods file ', &
                '(nobs=',n_obs,') or reading errors detected (ierr=',ierr,')'    
        stop
      else
        print *,'number of obs in file =',n_obs
      endif
!
! allocate all arrays
      allocate (xsum(n_ksmax,2))
      allocate (xsumsq(n_ksmax,2))
      allocate (xmean(n_ksmax))
      allocate (xmsqr(n_ksmax))
      allocate (xstdv(n_ksmax))
      allocate (xbins(n_ksmax,n_ksmax,5,2))
      allocate (iaccept(n_ksmax,2))
      allocate (ibins(n_ksmax,n_ksmax,2))
!
! initialize arrays
      ibins(:,:,1)=0
      iaccept(:,1)=0
      xsum(:,1)=zero
      xsumsq(:,1)=zero
      xmean(:)=zero
      xmsqr(:)=zero
      xstdv(:)=zero
      xbins(:,:,:,1)=zero
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
! Determine if obs longitude in desired range
        lon_range2=(lon_range1 .and. xlon >= lonW .and. xlon <= lonE) .or. & 
           ((.not. lon_range1) .and. (xlon >= lonW .or. xlon <= lonE)) 
!
        lev_OK=.false.
        if ( (nsave < nlook)                    .and. & 
             (ods%data%qcexcl(n) == 0)          .and. &
             (nint(ods_lev) <= n_ksmax)         .and. &
             (abs(ods%data%omf(n)) < 1.e8)      .and. &    
             (ods%data%lat(n) > latS)           .and. &
             (ods%data%lat(n) < latN)           .and. &
             lon_range2                        ) then
          xlev=ods%data%lev(n)  ! sat channel
          lev_OK=.true.
        endif
!
        if (lev_OK) then
          nsave=nsave+1
          ksave(nsave)=ods%data%ks(n)
          xsave(nsave,1)=ods%data%lat(n)*pifac
          xsave(nsave,2)=ods%data%lon(n)*pifac
          xsave(nsave,3)=xlev 
          xsave(nsave,4)=ods%data%time(n)
          if (c_comp == 'omf') then
            xsave(nsave,5)=ods%data%omf(n)
          elseif (c_comp == 'oma') then
            xsave(nsave,5)=ods%data%oma(n)
          elseif (c_comp == 'amb') then
            xsave(nsave,5)=ods%data%oma(n) - ods%data%omf(n)   ! a-b
          else 
            xsave(nsave,5)=ods%data%obs(n)
          endif 
        endif
!
      enddo  ! loop over all obs
!
      do n1=1,nsave                   ! loop over all saved obs
        ks1=ksave(n1)
        nch1=nint(xsave(n1,3))
        iaccept(nch1,1)=iaccept(nch1,1)+1
        xsum(nch1,1)=xsum(nch1,1)+xsave(n1,5)
        xsumsq(nch1,1)=xsumsq(nch1,1)+xsave(n1,5)*xsave(n1,5)
        do n2=1,nsave   
          ks2=ksave(n2)
          if (ks1 == ks2) then  ! obs from same sounding
            nch2=nint(xsave(n2,3))
            xbins(nch1,nch2,1,1)=xbins(nch1,nch2,1,1) + &
                                   xsave(n1,5)*xsave(n2,5)
            xbins(nch1,nch2,2,1)=xbins(nch1,nch2,2,1) + xsave(n1,5)
            xbins(nch1,nch2,3,1)=xbins(nch1,nch2,3,1) + &
                                   xsave(n1,5)*xsave(n1,5)
            xbins(nch1,nch2,4,1)=xbins(nch1,nch2,4,1) + xsave(n2,5)
            xbins(nch1,nch2,5,1)=xbins(nch1,nch2,5,1) + &
                                   xsave(n2,5)*xsave(n2,5)
            ibins(nch1,nch2,1)=ibins(nch1,nch2,1)+1
          endif
        enddo    
      enddo      
!
! compute means and standard deviations
      do iks=1,n_ksmax
        if (iaccept(iks,1) > 0) then
          xmean(iks)=xsum(iks,1)/iaccept(iks,1)
          xmsqr(iks)=xsumsq(iks,1)/iaccept(iks,1) - xmean(iks)*xmean(iks)
          if (xmsqr(iks) > zero) then
            xstdv(iks)=sqrt(xmsqr(iks))
          endif
        endif
      enddo 
!
! print table of means and standard deviations
      print *,' '
      print *,'Table 1: counts, means, stdvs'
      do iks=1,n_ksmax
        print ('(2i8,2f12.6)'),iks,iaccept(iks,1),xmean(iks),xstdv(iks)
      enddo
!
! Determine parameters for printing a subset of distributions
      i_list_skip=1+(n_ksmax-1)/15
      i=0
      do iks=1,n_ksmax,i_list_skip
        if (i < 15) then
          i=i+1
          i_list(i)=iks
        endif 
      enddo
      i_list_last=i
!
! Print a subset of distributions
      print *,' '
      print *,'Table 2: sample sizes (in hundreds)'
      print ('(6x,a,15i6)'),'channel=',i_list(1:i_list_last)
      do i=1,i_list_last
        i_values(i)=iaccept(i_list(i),1)
      enddo 
      print ('(6x,a,15i6)'),'iaccept=',i_values(1:i_list_last)
      print *,' '
      do n=1,n_ksmax
        do i=1,i_list_last
          i_values(i)=ibins(n,i_list(i),1)/100
        enddo 
        print ('(i4,10x,15i6)'),n,i_values(1:i_list_last)
      enddo
!
! Print a subset of distributions
      print *,' '
      print *,'Table 3: sample of correlations'
      print ('(6x,a,15i6)'),'channel=',i_list(1:i_list_last)
      do i=1,i_list_last
        i_values(i)=iaccept(i_list(i),1)
      enddo 
      print ('(6x,a,15i6)'),'iaccept=',i_values(1:i_list_last)
      print *,' '
      do n=1,n_ksmax
        do i=1,i_list_last
          if (ibins(n,i_list(i),1) > 0) then
            xm1=xbins(n,i_list(i),2,1)/ibins(n,i_list(i),1)
            xm2=xbins(n,i_list(i),4,1)/ibins(n,i_list(i),1)
            xs1=xbins(n,i_list(i),3,1)/ibins(n,i_list(i),1)
            xs2=xbins(n,i_list(i),5,1)/ibins(n,i_list(i),1)
            xa=(xs1-xm1*xm1)*(xs2-xm2*xm2)
            if (xa > zero) then   
              xb=xbins(n,i_list(i),1,1)/ibins(n,i_list(i),1)-xm1*xm2
              x_values(i)=xb/sqrt(xa)
            else
              x_values(i)=zero
            endif
          else
            x_values(i)=zero
          endif
        enddo 
        print ('(i4,10x,15f6.2)'),n,x_values(1:i_list_last)
      enddo
!
! read file pf previously accumulated values
      if (n_count > 1) then
        open (unit=10,file=file_in,form='unformatted')
        read (10) n_count2,n_ksmax2,n_date2,n_time2,       &
                  latN2,latS2,lonW2,lonE2,c_comp2,c_file2
        read (10) iaccept(:,2),xsum(:,2),xsumsq(:,2)
        read (10) ibins(:,:,2) 
        read (10) xbins(:,:,:,2) 
        close (10)     
        print *,'File read'      
!
! check that parameters on accumulation file to be updated agree with program
        if ( (n_ksmax2 /= n_ksmax)           .or. &
             (latN2 /= latN)                 .or. &
             (latS2 /= latS)                 .or. &
             (lonW2 /= lonW)                 .or. &
             (lonE2 /= lonE)                 .or. &
             (c_comp2 /= c_comp)           ) then
          print *,' '
          print *,'WRONG FILE ACCESSED X X X X X' 
          print *,'Parameters on file are: '
          print *,'n_ksmax,c_comp=',n_ksmax2,c_comp2
          print *,'latN,latS,lonW,lonE=',latN2,latS2,lonW2,lonE2
          print *,'Parameters specified in current program are: '
          print *,'n_ksmax,c_comp=',n_ksmax,c_comp
          print *,'latN,latS,lonW,lonE=',latN,latS,lonW,lonE
          stop
        endif
!
! update accumulations
        iaccept(:,1)=iaccept(:,1)+iaccept(:,2)
        xsum(:,1)=xsum(:,1)+xsum(:,2)
        xsumsq(:,1)=xsumsq(:,1)+xsumsq(:,2)
        ibins(:,:,1)=ibins(:,:,1)+ibins(:,:,2)
        xbins(:,:,:,1)=xbins(:,:,:,1)+xbins(:,:,:,2)
!
      endif
!
! write file of updated accumulated values
      open (unit=10,file=file_out,form='unformatted')
      write (10) n_count,n_ksmax,n_date,n_time, &
                 latN,latS,lonW,lonE,c_comp,c_file
      write (10) iaccept(:,1),xsum(:,1),xsumsq(:,1)
      write (10) ibins(:,:,1) 
      write (10) xbins(:,:,:,1) 
      close (10)     
!
      call ods_clean_RE (ods)
      print *,'Program end for time=',n_count,n_date,n_time
!
      end program channel_sumods
