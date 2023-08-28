   program diff_stats
   implicit none
   integer, parameter :: nbins=90
   integer, parameter :: nrmax=1000000
   integer :: nread,n,ib,nbad,ndatamax,ndataib,nc,nb,nf,ng
   integer :: idata(nrmax), idataBKG(nrmax), idataGSI(nrmax)
   real(4) :: data(7,nrmax), dataBKG(7,nrmax), dataGSI(6,nrmax), datadiff(7)
   real(4) :: omf(2,nbins)
   real(4) :: xdata1, latmin, latmax
   real(8) :: xbins(nbins,6), z, dz, xr, xx
   real(8) :: pbins(nbins,6), p, dplog, plog
   integer :: ibins(nbins), pibins(nbins)
   character(len=8) :: cnread, clatmin, clatmax
!
   xdata1=0.
   call GetArg( 1_4, clatmin)
   call GetArg( 2_4, clatmax)
   call GetArg( 3_4, cnread)
   read (cnread,'(i8)') nread
   read (clatmin,'(f8.2)') latmin
   read (clatmax,'(f8.2)') latmax
   print *,'nread=',nread,nrmax
   print *,'latmin,latmax=',latmin,latmax
!
   open (10,file='/discover/home/rerrico/GOWASP_3/GPSROTEST/omf.txt')
   do n=1,90 ! nbins
     read (10,*) ib,z,nc,omf(:,n)
   enddo
   close (10)
!
!   open (93,file='/discover/home/rerrico/GOWASP_3/GPSROTEST/xxgg45.txt',form='unformatted')
   open (93,file='/discover/nobackup/rerrico/TEST_GPS_517/allNR1p.txt',form='unformatted')
   open (94,file='/discover/nobackup/rerrico/TEST_GPS_517/allBKG0p.txt',form='unformatted')
   open (95,file='/discover/nobackup/rerrico/vcorr_OUT_all.txt')
   nbad=0
   do n=1,min(nread,nrmax)
     read (93) idata(n),data(:,n)
     read (94) idataBKG(n),dataBKG(:,n)
     read (95,'(i7,2f10.4,4e15.6 )'),idataGSI(n),dataGSI(:,n)
!ods%data%lat(n),ods%data%lon(n),ods%data%xm(n), &
!ods%data%obs(n),ods%data%omf(n),ods%data%obs(n)-ods%data%omf(n)



    if (abs(data(1,n)+50.2)<0.3 .and. abs(data(2,n)-159.4)<0.3) then
      print ('(a,i7,i8,1p7e13.4)'),'NNNR',n,idata(n),data(:,n) 
      print ('(a,i7,i8,1p7e13.4)'),'NNNB',n,idataBKG(n),dataBKG(:,n) 
    endif
!    if (n> 300    .and. n< 700   ) print ('(a,i7,i8,1p6e14.5)'),'NNNR',n,idata(n),data(:,n) 
!    if (n> 300    .and. n< 700   ) print ('(a,i7,i8,1p6e14.5)'),'NNNB',n,idataBKG(n),dataBKG(:,n) 
   
!if (idata(n) /= idataBKG(n) .or. data(1,n) /= dataBKG(1,n) .or. &
!    data(2,n) /= dataBKG(2,n) .or. data(3,n) /= dataBKG(3,n) ) then 
!  nbad=nbad+1
!   if (abs(idataBKG(n) - 4751) <5 .or. abs(idata(n) - 4751) <5) then
!  if (nbad < 100) then
!   print ('(a,i4,3i8,4f10.3,2f10.1)'),'HHH',nbad,n,idata(n),idataBKG(n),data(1,n),dataBKG(1,n), &
!         data(2,n),dataBKG(2,n),data(3,n),dataBKG(3,n)
!  endif

   enddo
!
   ibins(:)=0
   xbins(:,:)=0.
   pibins(:)=0
   pbins(:,:)=0.
   dz=6.0e4/real(nbins)
   dplog=log(1.e5)/real(nbins)
   ndataib=int(1+16.*1000./dz)
   nbad=0
   do n=1,min(nread,nrmax)
 
   if (data(1,n) >= latmin .and. data(1,n) <= latmax ) then

   nf=0
   ng=0
   do  nb=1,min(nread,nrmax)
     if (data(1,n) == dataBKG(1,nb) .and. data(2,n) == dataBKG(2,nb) .and. data(3,n) == dataBKG(3,nb)) then
       nf=nb
     endif
     if (abs(data(1,n)-dataGSI(1,nb))<0.1 .and. abs(data(2,n)-dataGSI(2,nb))<0.1 .and. &
         abs(data(3,n)-dataGSI(3,nb))<1.0) then
       ng=nb
     endif
   enddo

    if (nf /= 0 .and. ng /= 0) then 
      if (data(5,n) > .999 .or. data(5,n) < 0. .or. &
          data(4,n) > .999 .or. data(4,n) < 0. .or. data(3,n) < 0. .or. &
          dataBKG(5,nf) > .999 .or. dataBKG(5,nf) < 0. .or. &
          dataBKG(4,nf) > .999 .or. dataBKG(4,nf) < 0. .or. dataBKG(3,nf) < 0.) then 
        nf=0
      endif
    endif
    
    if (nf == 0 .or. ng==0) then    
       nbad=nbad+1
     else
!       datadiff(4)=data(4,n)-data(5,n)    ! ROPP_NR-GSI_NR
!       datadiff(5)=dataBKG(4,n)-data(5,n) ! ROPP_BKG-GSI_NR
!       datadiff(6)=data(4,n)-dataBKG(5,n) ! O-F  ROPP_NR-GSI_BKG
       datadiff(4)=data(4,n)-dataGSI(4,ng)    ! delta O = ROPP NR - ROPP ODS 
!       datadiff(5)=dataBKG(5,n)-data(5,n)    ! ROPP_NR-GSI_NR
       datadiff(5)=dataBKG(5,nf)-dataGSI(6,ng)    ! delta GSI = GSI_BKG - GSI_ODS
       datadiff(6)=dataGSI(5,ng) ! O-F  ROPP_ODAS-GSI_ODS
!
! sort into bins according to impact parameter - radius of curvature
       ib=1+int(data(3,n)/dz)
       ib=min(ib,nbins)    
       ibins(ib)=ibins(ib)+1
       xbins(ib,1)=xbins(ib,1)+datadiff(4)
       xbins(ib,2)=xbins(ib,2)+datadiff(5)
       xbins(ib,3)=xbins(ib,3)+datadiff(6)
       xbins(ib,4)=xbins(ib,4)+datadiff(4)*datadiff(4)
       xbins(ib,5)=xbins(ib,5)+datadiff(5)*datadiff(5)
       xbins(ib,6)=xbins(ib,6)+datadiff(6)*datadiff(6)

!QQQQ
       if (ib==6) print ('(a,1p8e14.5)'),'SSS',(data(nb,n),dataGSI(nb,ng),nb=1,4)
!
! sort into bins according to impact parameter - radius of curvature
       ib=1+int(log(data(7,n))/dplog)
       ib=min(ib,nbins)    
       pibins(ib)=pibins(ib)+1
       pbins(ib,1)=pbins(ib,1)+datadiff(4)
       pbins(ib,2)=pbins(ib,2)+datadiff(5)
       pbins(ib,3)=pbins(ib,3)+datadiff(6)
       pbins(ib,4)=pbins(ib,4)+datadiff(4)*datadiff(4)
       pbins(ib,5)=pbins(ib,5)+datadiff(5)*datadiff(5)
       pbins(ib,6)=pbins(ib,6)+datadiff(6)*datadiff(6)
!
     endif 
   endif  ! test on lat range
   enddo
!
   do ib=1,nbins
     if (ibins(ib) > 0) then
       xbins(ib,:)=xbins(ib,:)/ibins(ib)
     endif
     if (pibins(ib) > 0) then
       pbins(ib,:)=pbins(ib,:)/pibins(ib)
     endif
   enddo
!
   do ib=1,nbins
     if (ibins(ib) > 0) then
       do n=4,6
         xbins(ib,n)=xbins(ib,n)-xbins(ib,n-3)*xbins(ib,n-3)
         if (xbins(ib,n) > 0.) xbins(ib,n)=sqrt(xbins(ib,n))
       enddo
     endif
     if (pibins(ib) > 0) then
       do n=4,6
         pbins(ib,n)=pbins(ib,n)-pbins(ib,n-3)*pbins(ib,n-3)
         if (pbins(ib,n) > 0.) pbins(ib,n)=sqrt(pbins(ib,n))
       enddo
     endif
   enddo
!
   print *,' '
   print *,'STATS Z'
   do ib=nbins,1,-1
     z=(ib-1)*dz/1000.
     print ('(i3,f9.2,i6,1p6e12.3)'),ib,z,ibins(ib),xbins(ib,:)
   enddo
!
   print *,' '
   print *,'STATS P'
   do ib=1,nbins
     plog=(ib-1)*dplog
     p=exp(plog)
     print ('(i3,f9.2,i6,1p6e12.3)'),ib,p,pibins(ib),pbins(ib,:)
   enddo
!
   print *,' '
   print *,'nbad=',nbad
!
   print *,' '
   do ib=nbins,1,-1
     z=(ib-1)*dz/1000.
     if (omf(2,ib) > 0.) then
       xx=xbins(ib,6)/omf(2,ib)
     else
       xx=0.
     endif
     print ('(i3,f7.2,f8.4,1p2e12.3)'),ib,z,xx,omf(2,ib),xbins(ib,6)
   enddo


   end program diff_stats  


