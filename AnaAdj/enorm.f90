   program enorm
!
!  program for computing enorm of ana or bkg error and file of fields of 
!  d(enorm)/d(field).
!  Assumes that all grids are ordered SP to NP, and have the same vertical 
!  eta coordinates.
!
   use m_prnt_test
   implicit none
!
   logical, parameter :: ltvq=.true.  ! adjust for use of Tv rather than T 
   integer, parameter :: ifield1=5    ! max # of fields defining norm
   logical :: ltest                   ! true if run to test (AD+BD)(AD-BD)
!
   integer, parameter :: iunit_norm=32
   integer :: imax,jmax,kmax
   integer :: kmaxp1
   integer :: ier
   integer :: nf, k
   integer :: iprnt
!
   real(4), parameter :: Tref=270.       ! referance temperature
   real(4), parameter :: pref=1.e5       ! referance surface pessure Pa
   real(4), parameter :: Rgas=287.04     ! gas constant for dry air 
   real(4), parameter :: Cpgas=1004.6    ! specific heat capacity for dry air
   real(4), parameter :: Lvapor=2.5104e6 ! latent heat of vaporization
   real(4), parameter :: tvfac=0.622     ! ratio of mole weight H20 : dry air
   real(4), parameter :: qfac=0.3        ! factor for reducing moist norm
   real(4) :: scaleF, scaleT, scaleQ, scaleW, scaleP
   real(4) :: enorm_sum
   real(4) :: lon_first
   real(4) :: latN, latS, lonW, lonE, pmax, pmin
   real(4), allocatable :: enorms(:)
   real(4), allocatable :: q(:,:,:,:)
   real(4), allocatable :: f2d(:,:,:)
   real(4), allocatable :: f1(:,:,:)
   real(4), allocatable :: f2(:,:,:)
   real(4), allocatable :: f3(:,:,:)
   real(4), allocatable :: akbk(:,:)
   real(4), allocatable :: vweights(:)
   real(4), allocatable :: hweights(:,:)
!
   character(len=1)   :: ctest
   character(len=3)   :: caorb
   character(len=4)   :: field_name(ifield1)
   character(len=14)  :: cdatetime
   character(len=220) :: cfile_ana
   character(len=220) :: cfile_nr
   character(len=220) :: cfile_adj
   character(len=220) :: file_names(3)
   character(len=8) :: clatN, clatS, clonW, clonE, cpmax, cpmin, cnorm
!
! read arguments
   call GetArg( 1_4, cfile_ana)
   call GetArg( 2_4, cfile_nr)
   call GetArg( 3_4, cfile_adj)
   call GetArg( 4_4, cdatetime)
   call GetArg( 5_4, clatN)
   call GetArg( 6_4, clatS)
   call GetArg( 7_4, clonW)
   call GetArg( 8_4, clonE)
   call GetArg( 9_4, cpmax)
   call GetArg(10_4, cpmin)
   call GetArg(11_4, cnorm)
   call GetArg(12_4, ctest)
   call GetArg(13_4, caorb)
!
   if (ctest == 'T') then
     ltest=.true.
   else
     ltest=.false.
   endif
!
   file_names(1)=cfile_nr
   file_names(2)=cfile_ana
   file_names(3)=cfile_adj
   read (clatN,('(f8.2)')) latN
   read (clatS,('(f8.2)')) latS
   read (clonW,('(f8.2)')) lonW
   read (clonE,('(f8.2)')) lonE
   read (cpmax,('(f8.1)')) pmax
   read (cpmin,('(f8.1)')) pmin
!
! Set field names; The order for 1-3 matters. 6-11 are 3d fields. Rest are 2d
! 
   field_name(1:ifield1)=(/'ps','sphu','tv','u','v'/)
! 
! Set scale factors depending on field and norm
   if (trim(cnorm) == 'txe' .or. trim(cnorm) == 'twe' .or. &
       trim(cnorm) == 'ape' .or. trim(cnorm) == 'tte') then
     scaleT=0.5*Cpgas/Tref
   else
     scaleT=0.
   endif
!
   if (trim(cnorm) == 'qxe' .or. trim(cnorm) == 'twe' ) then
     scaleQ=0.5*qfac*Lvapor**2/(Cpgas*Tref)
   else
     scaleQ=0.
   endif
!
   if (trim(cnorm) == 'txe' .or. trim(cnorm) == 'twe' .or. &
       trim(cnorm) == 'kxe' ) then
     scaleW=0.5
   else
     scaleW=0.
   endif
!
   if (trim(cnorm) == 'txe' .or. trim(cnorm) == 'twe' .or. &
       trim(cnorm) == 'tps' .or. trim(cnorm) == 'ape' ) then
     scaleP=0.5*Rgas*Tref/pref**2
   else
     scaleP=0.
   endif
!
   call get_grid_info (imax,jmax,kmax,3,lon_first,file_names,ier)
   if (ier > 0) then
     print *,'ERROR IN READING FILE INFO: ier=',ier
     print *,'Either file not read or grid parameter mismatch'
     stop
   endif
!
   print *,' '
   print ('(a,4f8.2)'),'latN,latS,lonW,lonE =',latN,latS,lonW,lonE
   print ('(a,2f10.2)'),'pmax,pmin =',pmax,pmin
   print ('(a,3i5,f12.4)'),'imax,jmax,kmax,lon_first =',imax,jmax,kmax,lon_first
   print ('(2a,1p4e14.4)'),'norm,scaleT,scaleQ,scaleW,scaleP = ', &
                           trim(cnorm),scaleT,scaleQ,scaleW,scaleP
!
! allocate all grid-dimension arrays
   kmaxp1=kmax+1
   allocate (q(imax,jmax,kmax,3))
   allocate (f2d(imax,jmax,3))
   allocate (f1(imax,jmax,kmax))
   allocate (f2(imax,jmax,kmax))
   allocate (f3(imax,jmax,kmax))
   allocate (akbk(kmaxp1,2))
   allocate (vweights(kmax))
   allocate (hweights(imax,jmax))
   allocate (enorms(ifield1))
!
! compute vertical and horizontal weights
   call read_akbk (kmaxp1,kmax,10,'akbk_file',.true.,akbk,ier)   
   call compute_vweights (kmaxp1,kmax,pmax,pmin,akbk,vweights)
   call compute_hweights (imax,jmax,latN,latS,lonW,lonE,lon_first,hweights)
!
! Open file if test output values printed to file
   if (ltest) then 
     call ptest_open (ltvq,caorb,scaleT,scaleQ,scaleW,scaleP)
     vweights(:)=1.
     hweights(:,:)=1.
   endif
! 
! Computations for ps field   
   nf=1  
   call read_nc4_2dfield (imax,jmax,cfile_ana,field_name(nf), &
                          f1(:,:,1),.true.,.true.,ier)
   call read_nc4_2dfield (imax,jmax,cfile_nr,field_name(nf), &
                          f2(:,:,1),.true.,.true.,ier)
   f2(:,:,1)=f2(:,:,1)*100.  ! change mb to Pa for NR data
   call compute_adj (imax,jmax,1,scaleP,hweights,vweights,enorms(nf), &
                     f1(:,:,1:1),f2(:,:,1:1),f3(:,:,1:1))
   call write_nc4_2dfld (imax,jmax,cdatetime,field_name(nf),cfile_adj, &
                         f3(:,:,1))
!
   if (ltest) then
     call ptest_fill (1,imax,jmax,1,'FIELD ',field_name(nf),f1(:,:,1:1))
     call ptest_fill (2,imax,jmax,1,'TRUTH ',field_name(nf),f2(:,:,1:1))
     call ptest_fill (3,imax,jmax,1,'DE/DF ',field_name(nf),f3(:,:,1:1))
   endif
!
! Set adjoint of delp to agree with adjoint ps
   do k=2,kmax
     f3(:,:,k)=f3(:,:,1)
   enddo
   call write_nc4_3dfld (imax,jmax,kmax,cdatetime,'delp',cfile_adj,f3)
! 
! Computations for 3-d fields   
   do nf=2,ifield1
     print *,'processing nf,field=',nf,field_name(nf)
     if (field_name(nf) == 'tv') then
       scaleF=scaleT
     elseif (field_name(nf) == 'sphu') then
       scaleF=scaleQ
     else
       scaleF=scaleW
     endif 
! 
     call read_nc4_3dfield (imax,jmax,kmax,cfile_ana,field_name(nf), &
                            f1,.true.,.true.,ier)
     call read_nc4_3dfield (imax,jmax,kmax,cfile_nr,field_name(nf), &
                            f2,.true.,.true.,ier)
!
     if (ltest) then
       iprnt=3*(nf-1)
       call ptest_fill (iprnt+1,imax,jmax,kmax,'FIELD ',field_name(nf),f1)
       call ptest_fill (iprnt+2,imax,jmax,kmax,'TRUTH ',field_name(nf),f2)
     endif
!
     if (ltvq .and. field_name(nf) == 'tv') then
       f1(:,:,:)=f1(:,:,:)/(1.+tvfac*q(:,:,:,1))  ! change Tv to T for field
       f2(:,:,:)=f2(:,:,:)/(1.+tvfac*q(:,:,:,2))  ! change Tv to T for truth
       if (ltest) then
         iprnt=3*ifield1
         call ptest_fill (iprnt+1,imax,jmax,kmax,'FIELD ','T',f1)
         call ptest_fill (iprnt+2,imax,jmax,kmax,'TRUTH ','T',f2)
       endif
     endif
!
     call compute_adj (imax,jmax,kmax,scaleF,hweights,vweights,enorms(nf), &
                       f1,f2,f3)
!
! Corrections for virtual temperature as control variable
! (dE/dTv)_q = (dE/dT)q / (dTv/dT)_q 
! (dE/dq)_Tv = (dE/dq)T - (dE/dTv)_q (dTv/dq)_T 
     if (ltvq .and. field_name(nf) == 'tv') then
       f3(:,:,:)=f3(:,:,:)/(1.+tvfac*q(:,:,:,1))
       q(:,:,:,3)=q(:,:,:,3)-f3(:,:,:)*tvfac*f1(:,:,:)
       if (ltest) then
         iprnt=3*(nf-1)
         call ptest_fill (6,imax,jmax,kmax,'DE/DF ',field_name(2),q(:,:,:,3))
       endif
       call write_nc4_3dfld (imax,jmax,kmax,cdatetime,'sphu', &
                             cfile_adj,q(:,:,:,3))
!
! if Tv control variable, save q to compute T and (dE/dF) at constant Tv, q
     else if (ltvq .and. field_name(nf) == 'sphu') then  
       q(:,:,:,1)=f1(:,:,:)
       q(:,:,:,2)=f2(:,:,:)
       q(:,:,:,3)=f3(:,:,:)
     endif
!
     if (.not. (ltvq .and. field_name(nf) == 'sphu')) then  
       if (ltest) then
         iprnt=3*nf
         call ptest_fill (iprnt,imax,jmax,kmax,'DE/DF ',field_name(nf),f3)
       endif
       call write_nc4_3dfld (imax,jmax,kmax,cdatetime,field_name(nf), &
                             cfile_adj,f3)
     endif
!  
   enddo
!
! Set other adjoint values to 0
   f3(:,:,:)=0.   
   call write_nc4_2dfld (imax,jmax,cdatetime,'ts',cfile_adj,f3(:,:,1))
   call write_nc4_3dfld (imax,jmax,kmax,cdatetime,'ozone',cfile_adj,f3)
   call write_nc4_3dfld (imax,jmax,kmax,cdatetime,'qitot',cfile_adj,f3)
   call write_nc4_3dfld (imax,jmax,kmax,cdatetime,'qltot',cfile_adj,f3)
   call write_nc4_3dfld (imax,jmax,kmax,cdatetime,'qrtot',cfile_adj,f3)
   call write_nc4_3dfld (imax,jmax,kmax,cdatetime,'qstot',cfile_adj,f3)
!   
   if (ltest) then
     call ptest_close
   endif
!
! print norms to log file
   print *,' '
   print *,'file_out written: ',trim(cfile_adj)
   enorm_sum=sum(enorms(:))
   print *,' '
   print ('(a,1p1e13.4)'),'enorm sum =',enorm_sum
   do nf=1,ifield1
     print ('(a,a4,1p1e13.4)'),'enorm ',field_name(nf),enorms(nf)
   enddo
!
! write norms to separate file
   open (iunit_norm,file='file_norms.txt')
   write (iunit_norm,'(5a,i2)') 'data_type= ',caorb,'  norm_type=',cnorm, &
                                '  nfields=',ifield1 
   write (iunit_norm,'(a,1p1e13.4)') 'enorm sum=',enorm_sum
   do nf=1,ifield1
     write (iunit_norm,'(a,a4,1p1e13.4)') 'enorm ',field_name(nf),enorms(nf)
   enddo
   close (iunit_norm)
!
   print *,' '
   print *,'Program Done'
!
   end program enorm
!
!
! x x x x x x x x x x x x x x  x x x x x x x x x x x x x x x x x x x x x  
!
   subroutine compute_vweights (kmaxp1,kmax,pmax,pmin,akbk,vweights)
   implicit none 
   integer :: kmaxp1,kmax
   real(4) :: pmax,pmin
   real(4) :: akbk(kmaxp1,2)
   real(4) :: vweights(kmaxp1)
   integer :: k
   real(4), parameter :: ps=1.e5
   real(4) :: pmid,delp
!
   delp=0.
   do k=1,kmax
     pmid=0.5*(akbk(k+1,1)+akbk(k,1)+ps*(akbk(k+1,2)+akbk(k,2)))
     if (pmid <= pmax .and. pmid >= pmin) then 
       vweights(k)=akbk(k+1,1)-akbk(k,1)+ps*(akbk(k+1,2)-akbk(k,2))
       delp=delp+vweights(k)
     else
       vweights(k)=0.
     endif
   enddo  
!
   if (delp > 0.) then
     vweights(:)=vweights(:)/delp
   endif
!
   end subroutine compute_vweights 
!
!
! x x x x x x x x x x x x x x  x x x x x x x x x x x x x x x x x x x x x  
!
   subroutine compute_hweights (imax,jmax,latN,latS,lonW,lonE,lon_first, &
                                hweights)
   implicit none 
   integer :: imax,jmax
   real(4) :: latN,latS,lonW,lonE
   real(4) :: lon_first
   real(4) :: hweights(imax,jmax)
!
   logical :: lon_case
   integer :: i,j
   integer :: ilon(imax)                ! =1 if pt included
   real(4), parameter :: epsilon=1.e-4  ! small value to account for round-off

   real(4) :: elatN,elatS,elonW,elonE
   real(4) :: pifac,dlat,dlon,latedge,longp,area
   real(4) :: latgp(jmax)
   real(4) :: slat(jmax+1)  ! sine of lats at edges
!
   elatN=latN+epsilon
   elatS=latS-epsilon
!
   if (lonW < 0.) then
     lonW=lonW+360.
   endif
   if (lonE < 0.) then
     lonE=lonE+360.
   endif
   elonW=lonW-epsilon
   elonE=lonE+epsilon
   if (lonW < lonE) then
     lon_case=.true.
   else
     lon_case=.false.
   endif
!
! determine lats of grid points and delta sine of lats at edges of grid boxes
   pifac=atan(1.)/45.
   dlat=180./(jmax-1)
   slat(1)=-1.
   latgp(1)=-90.
   do j=2,jmax   
     latedge=-90.+(j-1.5)*dlat
     latgp(j)=-90+(j-1)*dlat
     slat(j)=sin(latedge*pifac)
   enddo
   latgp(jmax)=90. 
   slat(jmax+1)=1.
!
! determine longitudes of grid points and if point falls within desired area
   dlon=360./real(imax)
   do i=1,imax
     longp=lon_first+(i-1)*dlon
     if (longp < 0.) then
       longp=360.+longp
     elseif (longp > 360.) then
       longp=longp-360.
     endif
!
     if ((lon_case .and. longp >= elonW .and. longp <= elonE) .or. & 
         ((.not. lon_case) .and. (longp >= elonW .or. longp <= elonE)) ) then
       ilon(i)=1
     else
       ilon(i)=0
     endif
   enddo
!
   area=0.           
   hweights(:,:)=0.  ! default value means do not include point
   do j=1,jmax
     if (latgp(j) >= elatS .and. latgp(j) <= elatN) then 
       do i=1,imax
         if (ilon(i) == 1) then 
           hweights(i,j)=0.5*(slat(j+1)-slat(j))
           area=area+hweights(i,j)
         endif
       enddo
     endif 
   enddo
!
   if (area > 0.) then 
     hweights(:,:)=hweights(:,:)/area
   endif
!
   end subroutine compute_hweights 
!
!
! x x x x x x x x x x x x x x  x x x x x x x x x x x x x x x x x x x x x  
!   
   subroutine compute_adj (imax,jmax,kmax,scaleF,hweights,vweights,enorm1, &
                          f1,f2,f3)
   implicit none 
   integer :: imax,jmax,kmax
   real(4) :: hweights(imax,jmax),vweights(kmax)
   real(4) :: scaleF, enorm1
   real(4) :: f1(imax,jmax,kmax),f2(imax,jmax,kmax),f3(imax,jmax,kmax)
   integer :: i,j,k
   real(4) :: w,vw,fd,w2
!
   enorm1=0.
   do k=1,kmax
     if (kmax == 1) then 
       vw=1.
     else
       vw=vweights(k)
     endif
     vw=vw*scaleF
     do j=1,jmax
       do i=1,imax
         w=vw*hweights(i,j)
         w2=w*2.
         fd=f1(i,j,k)-f2(i,j,k)
         f3(i,j,k)=fd*w2
         enorm1=enorm1+fd*fd*w
       enddo
     enddo
   enddo
!   
   end subroutine compute_adj

