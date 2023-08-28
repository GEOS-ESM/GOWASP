      program hcorr_sumods
!
! Augment pre-computed sums of ods data for previous times by adding 
! contributions for a new time for subsequent computation of horizontal 
! correlations and covariances.
!
      use m_ods_RE
      use m_sat_info_table, only : sat_info_table_read
      use m_sat_info_table, only : sat_info_table_get_1i
!
!
      implicit none
!
      logical :: lev_OK
      logical :: lookat_RH
      logical :: lon_range1, lon_range2
!
      integer, parameter :: rkinds=8
      integer, parameter :: ikinds=8
      integer, parameter :: nlook=1000000
      integer(ikinds), allocatable :: iaccept(:,:)
      integer(ikinds), allocatable :: ibins(:,:,:)
      integer :: ib
      integer :: iks,     iks2
      integer :: i
      integer :: i_list(15)
      integer :: i_list_last
      integer :: i_list_skip
      integer :: i_values(15)
      integer :: iqc
      integer :: ifield, itype
      integer :: kt_look, kt_look2, kt_ods
      integer :: k
      integer :: kx_look, kx_look2
      integer :: n, n1, n2
      integer :: nch1, nch2
      integer :: nsave
      integer :: n_bins,  n_bins2
      integer :: n_count, n_count2
      integer :: n_ksmax, n_ksmax2
      integer :: n_date,  n_date2
      integer :: n_time,  n_time2
      integer :: nlevs2
      integer :: n_obs
      integer :: argc
!
      real(rkinds), parameter :: zero=0._rkinds
      real(rkinds) :: pifac
      real(rkinds) :: earthr2
      real(rkinds) :: d, d1, d2, d3
      real(rkinds) :: latN, latS, lonW, lonE
      real(rkinds) :: latN2, latS2, lonW2, lonE2
      real(rkinds) :: delp, half_delp, delp2
      real(rkinds) :: ods_lev
      real(rkinds) :: x
      real(rkinds) :: xlev, xlon
      real(rkinds) :: xm1, xm2, xs1, xs2, xa, xb, xc
      real(rkinds) :: x_max, x_max2
      real(rkinds) :: x1,    x12
      real(rkinds) :: xbin,  xbin2
      real(rkinds) :: x_values(15)
      real(rkinds) :: xsave(nlook,5)
      real(rkinds), allocatable :: qsave(:)
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
      character(len=20) :: c_sat_inst
      character(len=220) :: c_sat_info_file
      character(len=220) :: c_file, c_file2
      character(len=11), parameter :: file_in='input_stats'
      character(len=12), parameter :: file_out='output_stats'
!
      type (ods_vect) :: ods
      character (len=30) :: type
      integer :: ierr
!
! The plevs are are conventional obs.  These should not be changed unless
! all the subsequent programs and applications are also altered.
      integer, parameter :: nlevs=10
      real(rkinds) :: plevs(nlevs)
      data plevs /1000., 925., 850., 700., 500., 400., 300., 200., 100., 50./
!       
! Read and check arguments
      argc = iargc()
      if (argc .ne. 16) then
        print *,' usage must be: hcorr_sumods.x ncount nbins xmax date',  &
                ' time latN latS lonW lonE c_delp c_kt c_kx c-comp', &
                ' filename sat_instr sat_info_file'
        stop
      endif
      call GetArg( 1_4, c_count)
      call GetArg( 2_4, c_nbins)
      call GetArg( 3_4, c_xmax)
      call GetArg( 4_4, c_date)
      call GetArg( 5_4, c_time)
      call GetArg( 6_4, c_latN)
      call GetArg( 7_4, c_latS)
      call GetArg( 8_4, c_lonW)
      call GetArg( 9_4, c_lonE)
      call GetArg(10_4, c_delp)
      call GetArg(11_4, c_kt)
      call GetArg(12_4, c_kx)
      call GetArg(13_4, c_comp)
      call GetArg(14_4, c_file)
      call GetArg(15_4, c_sat_inst)
      call GetArg(16_4, c_sat_info_file)
!
      read (c_count,'(i3)')   n_count
      read (c_nbins,'(i3)')   n_bins
      read (c_xmax, '(f8.4)') x_max
      read (c_date, '(i8)')   n_date
      read (c_time, '(i6)')   n_time
      read (c_latN, '(f6.1)') latN  
      read (c_latS, '(f6.1)') latS  
      read (c_lonW, '(f6.1)') lonW
      read (c_lonE, '(f6.1)') lonE
      read (c_delp, '(f6.1)') delp
      read (c_kt,   '(i3)')   kt_look
      read (c_kx,   '(i3)')   kx_look
!
      if ( (kx_look > 99) .and. (kx_look < 300 ) ) then ! conv obs
        half_delp=0.5*delp
        n_ksmax=nlevs               ! number of levels for conv obs 
        if (kt_look == 10) then     ! look at relative humidity
          lookat_RH=.true.
          allocate (qsave(nlook))
          kt_ods=11                 ! read q (specific humidity) to compute rh
        else
          lookat_RH=.false.
          kt_ods=kt_look
        endif
      else     ! sat obs
        call sat_info_table_read (c_sat_info_file,.true.,ierr)
        if (ierr /= 0) stop
        call sat_info_table_get_1i ('instr','nchan',c_sat_inst,n_ksmax,ierr)
        if (ierr /= 0) stop
        plevs(:)=0.
      endif
!
! Determine direction of longitude range
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
      print ('(3i5,f10.4,2i10)'),n_count,n_ksmax,n_bins,x_max,n_date,n_time
      print ('(4f8.1)'),latN,latS,lonW,lonE
      print ('(f8.1,2i5,2x,a3)'),delp,kt_look,kx_look,c_comp
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
      print *,'n_obs, ierr=',n_obs,ierr
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
      allocate (xbins(n_bins,5,n_ksmax,2))
      allocate (iaccept(n_ksmax,2))
      allocate (ibins(n_bins,n_ksmax,2))
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
      xbin=x_max/(n_bins-1)
      x1=xbin 
      pifac=atan(1._rkinds)/45._rkinds
      earthr2=2.*6371.
!
! if required, save temperatures so that RH can be computed.
      if (lookat_RH) then
        call q2rh (nlook,nlevs,n_obs,0,kx_look,   &
                   latN,latS,lonW,lonE,half_delp, &
                   plevs,xsave,qsave,ods,0)
      endif
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
        if ( (kx_look < 100) .or. (kx_look > 299 ) ) then    ! sat obs 
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
        else                                                 ! conv obs 
!       
          ifield=ods%data%kt(n)
          itype=ods%data%kx(n)
          if ( (nsave < nlook)               .and. & 
               (ods%data%qcexcl(n) == 0)     .and. &
               (ifield == kt_ods)            .and. &
               (itype == kx_look)            .and. &
               (abs(ods%data%omf(n)) < 1.e8) .and. &
               (ods%data%lat(n) > latS)      .and. &
               (ods%data%lat(n) < latN)      .and. &
               lon_range2                   ) then
!
! determine p-level of accepted obs
            do k=1,nlevs
              if  ( (ods_lev < plevs(k)+half_delp)  .and. &
                    (ods_lev > plevs(k)-half_delp) ) then
                lev_OK=.true.
                xlev=real(k)
                exit
              endif
            enddo
!
          endif  ! check on selected conv obs
        endif    ! check on whether sat or conv obs
!
        if (lev_OK) then
          nsave=nsave+1
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
          if (lookat_RH) then
            qsave(nsave)=ods%data%lev(n)
          endif
        endif
!
      enddo  ! loop over all obs
!
! if required, replace q by RH
      if (lookat_RH) then
        call q2rh (nlook,nlevs,n_obs,nsave,kx_look, &
                   latN,latS,lonW,lonE,half_delp,   &
                   plevs,xsave,qsave,ods,1)
        if (nsave == 0) then
          print *,'No T found for q obs to compute RH'
          stop
        endif
      endif
!
      do n1=1,nsave                   ! loop over all saved obs
        nch1=nint(xsave(n1,3))
        iaccept(nch1,1)=iaccept(nch1,1)+1
        xsum(nch1,1)=xsum(nch1,1)+xsave(n1,5)
        xsumsq(nch1,1)=xsumsq(nch1,1)+xsave(n1,5)*xsave(n1,5)
        do n2=n1,nsave ! include self product in first bin 
          nch2=nint(xsave(n2,3))
          if (nch1 == nch2) then
            if (n1 == n2) then
              ib=1
            else
              d1=sin(0.5*(xsave(n1,1)-xsave(n2,1)))
              d2=sin(0.5*(xsave(n1,2)-xsave(n2,2)))
              d3=cos(xsave(n1,1))*cos(xsave(n2,1))
              d=earthr2*asin(sqrt(d1*d1+d2*d2*d3))
              ib=2+d/xbin
            endif
            if (ib <= n_bins) then
              xbins(ib,1,nch1,1)=xbins(ib,1,nch1,1)+xsave(n1,5)*xsave(n2,5)
              xbins(ib,2,nch1,1)=xbins(ib,2,nch1,1)+xsave(n1,5)
              xbins(ib,3,nch1,1)=xbins(ib,3,nch1,1)+xsave(n1,5)*xsave(n1,5)
              xbins(ib,4,nch1,1)=xbins(ib,4,nch1,1)+xsave(n2,5)
              xbins(ib,5,nch1,1)=xbins(ib,5,nch1,1)+xsave(n2,5)*xsave(n2,5)
              ibins(ib,nch1,1)=ibins(ib,nch1,1)+1
            endif 
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
      print *,'Table 2: sample of separation distributions'
      print ('(6x,a,15i6)'),'channel=',i_list(1:i_list_last)
      do i=1,i_list_last
        i_values(i)=iaccept(i_list(i),1)
      enddo 
      print ('(6x,a,15i6)'),'iaccept=',i_values(1:i_list_last)
      print *,' '
      do n=1,n_bins
        do i=1,i_list_last
          i_values(i)=ibins(n,i_list(i),1)
        enddo 
        if (n == 1) then
          d=0.
        else  
          d=(n-1.5)*xbin
        endif
        print ('(i4,f10.0,15i6)'),n,d,i_values(1:i_list_last)
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
      do n=1,n_bins
        do i=1,i_list_last
          if (ibins(n,i_list(i),1) > 0) then
            xm1=xbins(n,2,i_list(i),1)/ibins(n,i_list(i),1)
            xm2=xbins(n,4,i_list(i),1)/ibins(n,i_list(i),1)
            xs1=xbins(n,3,i_list(i),1)/ibins(n,i_list(i),1)
            xs2=xbins(n,5,i_list(i),1)/ibins(n,i_list(i),1)
            xa=(xs1-xm1*xm1)*(xs2-xm2*xm2)
            if (xa > zero) then   
              xb=xbins(n,1,i_list(i),1)/ibins(n,i_list(i),1)-xm1*xm2
              x_values(i)=xb/sqrt(xa)
            else
              x_values(i)=zero
            endif
          endif
        enddo 
        if (n == 1) then
          d=0.
        else  
          d=(n-1.5)*xbin
        endif
        print ('(i4,f10.0,15f6.3)'),n,d,x_values(1:i_list_last)
      enddo
!
! read file pf previously accumulated values
      if (n_count > 1) then
        open (unit=10,file=file_in,form='unformatted')
        read (10) n_count2,n_ksmax2,n_date2,n_time2,n_bins2,x_max2, &
                  x12,xbin2,latN2,latS2,lonW2,lonE2,                &
                  delp2,kt_look2,kx_look2,nlevs2,c_comp2,c_file2
        read (10) plevs(:)
        read (10) iaccept(:,2),xsum(:,2),xsumsq(:,2)
        read (10) ibins(:,:,2) 
        read (10) xbins(:,:,:,2) 
        close (10)     
        print *,'File read'      
!
! check that parameters on accumulation file to be updated agree with program
        if ( (n_ksmax2 /= n_ksmax)           .or. &
             (n_bins2 /= n_bins)             .or. &
             (abs(x_max2-x_max)/x_max > .01) .or. & 
             (latN2 /= latN)                 .or. &
             (latS2 /= latS)                 .or. &
             (lonW2 /= lonW)                 .or. &
             (lonE2 /= lonE)                 .or. &
             (delp2 /= delp)                 .or. &
             (kt_look2 /= kt_look)           .or. &
             (kx_look2 /= kx_look)           .or. &
             (c_comp2 /= c_comp)           ) then
          print *,' '
          print *,'WRONG FILE ACCESSED X X X X X' 
          print *,'Parameters on file are: '
          print *,'n_ksmax,n_bins,x_max=',n_ksmax2,n_bins2,x_max2
          print *,'latN,latS,lonW,lonE=',latN2,latS2,lonW2,lonE2
          print *,'delp,kt_look,kx_look,c_comp=', &
                    delp2,kt_look2,kx_look2,c_comp2
          print *,'Parameters specified in current program are: '
          print *,'n_ksmax,n_bins,x_max=',n_ksmax,n_bins,x_max
          print *,'latN,latS,lonW,lonE=',latN,latS,lonW,lonE
          print *,'delp,kt_look,kx_look,c_comp=',delp,kt_look,kx_look,c_comp
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
      write (10) n_count,n_ksmax,n_date,n_time,n_bins,x_max, &
                 x1,xbin,latN,latS,lonW,lonE,                &
                 delp,kt_look,kx_look,nlevs,c_comp,c_file
      write (10) plevs(:)
      write (10) iaccept(:,1),xsum(:,1),xsumsq(:,1)
      write (10) ibins(:,:,1) 
      write (10) xbins(:,:,:,1) 
      close (10)     
!
      call ods_clean_RE (ods)
      print *,'Program end for time=',n_count,n_date,n_time
!
      end program hcorr_sumods
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
      subroutine q2rh   (nlook,nlevs,n_obs,nsave,kx_look, &
                         latN,latS,lonW,lonE,half_delp,   &
                         plevs,xsave,qsave,ods,ifunc)
!      
      use m_ods_RE
!
      implicit none
!
      integer, parameter :: rkinds=8
      integer :: nlook,nlevs,n_obs,ifunc
      integer :: nsave,kx_look
      real(rkinds) :: latN, latS, lonW, lonE, half_delp
      real(rkinds) :: plevs(nlevs)
      real(rkinds) :: xsave(nlook,5)
      real(rkinds) :: qsave(nlook)
      type (ods_vect) :: ods
!
      logical :: lev_OK
      integer, parameter :: nk_dim=100000
      integer  :: k, klev, n, n1, n2  
      integer  :: ifield, itype
      integer, allocatable, save :: nk_save(:)      
      real(rkinds), parameter :: zero=0._rkinds
      real(rkinds) :: pifac
      real(rkinds)  :: vs, qs
      real(rkinds)  :: xlon, ods_lev
      real(rkinds), allocatable, save :: t_save(:,:,:)
!
      if (ifunc==0) then   ! init
        allocate (t_save(nk_dim,nlevs,5))
        allocate (nk_save(nlevs))
        nk_save(:)=0       ! count of observations at each desired plev
        pifac=atan(1._rkinds)/45._rkinds
!
        do n=1,n_obs
!
! reset longitude so that it is in the range 0-360 
          xlon=mod(ods%data%lon(n),360.)
          if (xlon < zero)  xlon=xlon+360.
          if (xlon >= 360.) xlon=0.       ! accounts for round-off of xlon
!
! check for obs level set too large (nint function below otherwise fails)
          ods_lev=min(1.e7,ods%data%lev(n)) ! 1.e7 is just some large number 
!
          ifield=ods%data%kt(n)
          itype=ods%data%kx(n)
          if ( (itype        == kx_look)     .and. &
               (ifield       == 44     )     .and. & 
               (ods%data%qcexcl(n) == 0)     .and. &
               (abs(ods%data%obs(n)) < 400.) .and. &
                 (ods%data%lat(n) > latS)      .and. &
                 (ods%data%lat(n) < latN)      .and. &
                 (xlon > lonW)                 .and. &
                 (xlon < lonE)              )  then 
!
            lev_OK=.false.
            do k=1,nlevs
              if  ( (ods_lev < plevs(k)+half_delp) .and. &
                  (ods_lev > plevs(k)-half_delp) ) then
                lev_OK=.true.
                exit
              endif
            enddo 
!
            if (lev_OK .and. (nk_save(k) < nk_dim) ) then
              nk_save(k)=nk_save(k)+1
              t_save(nk_save(k),k,1)=ods%data%lat(n)*pifac
              t_save(nk_save(k),k,2)=ods%data%lon(n)*pifac
              t_save(nk_save(k),k,3)=ods%data%lev(n)
              t_save(nk_save(k),k,4)=ods%data%time(n)
              t_save(nk_save(k),k,5)=ods%data%obs(n)
            endif
!           
          endif  ! check on selected conv obs
        enddo    ! loop over obs
!
!
      else  ! find t obs corresponding to all q obs and replace q by rh         
!  
        n2=0
        do n1=1,nsave
          k=nint(xsave(n1,3))     
          do n=1,nk_save(k)
            if ( (xsave(n1,1) == t_save(n,k,1)) .and. &
                 (xsave(n1,2) == t_save(n,k,2)) .and. &
                 (qsave(n1)   == t_save(n,k,3)) .and. &
                 (xsave(n1,4) == t_save(n,k,4)) ) then
              n2=n2+1
              vs=6.112*exp(17.67*(t_save(n,k,5)-273.15)/ &
                                 (t_save(n,k,5)-29.65))  ! units mb            
              qs=622.*vs/(qsave(n1)-0.378*vs)            ! units g/kg 
              if (n2 < n1) then
                xsave(n2,1:4)=xsave(n1,1:4)
              endif 
              xsave(n2,5)=xsave(n1,5)/qs
              exit 
            endif
          enddo
        enddo
        nsave=n2
!
      endif
!              
      end subroutine q2rh 
