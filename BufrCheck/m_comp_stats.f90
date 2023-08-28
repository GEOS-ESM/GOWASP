   module m_comp_stats
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   use m_kinds, only : rkind1, rkind2
   use m_saved_data, only : obs_ids,obs_sing,obs_multi,count_types
   use m_saved_data, only : obs_types_max, obs_multi_dim1, obs_types_max
   use m_saved_data, only : nobs_sing_2,nobs_multi_2,nobs_all_2
   use m_saved_data, only : count_types_num
!
   implicit none
!
   private
   public :: comp_stats_meanvar
   public :: comp_stats_chan_corr
   public :: comp_stats_hcorr_rad
   public :: comp_stats_hcorr_conv1
   public :: comp_stats_vcorr
!
   integer, parameter :: earthr=6.371e6
   integer, parameter :: obs_dim1=700
   integer, parameter :: obs_dim2=70
!
   integer :: ocount(obs_dim1,obs_dim2)
   real(rkind1) :: omean(obs_dim1,obs_dim2,2)
   real(rkind1) :: ovars(obs_dim1,obs_dim2,2)
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine comp_stats_meanvar (dtype)
!
   implicit none
!
   character(len=*), intent(in) :: dtype
!
   integer :: olevs
   integer :: id,nobs
   integer :: n,k,j,kid
   integer :: klevs
   real(rkind1) :: dlevs
!
   ocount(:,:)=0.
   omean(:,:,:)=0.
   ovars(:,:,:)=0.
!
! Compute sums of obs values
   if (trim(dtype) == 'PREPBUFR') then
!
     klevs=10
     do n=1,nobs_multi_2(1)+nobs_sing_2(1)
       olevs=int(obs_ids(5,n,1))
       id=int(obs_ids(7,n,1))
       nobs=int(obs_ids(6,n,1))
       if (olevs == 1) then
         if (obs_sing(2,nobs,1) < 1.e9 .and. &
             obs_sing(2,nobs,2) < 1.e9 .and. &
             obs_sing(1,nobs,1) < 1.e9) then
           k=max(klevs-int(obs_sing(1,nobs,1)/10000.),1)
           ocount(k,id)=ocount(k,id)+1
           do j=1,2
             omean(k,id,j)=omean(k,id,j)+obs_sing(2,nobs,j)
           enddo
         endif
       else  
         do kid=1,olevs
           if (obs_multi(kid,2,nobs,1) < 1.e9 .and. &
               obs_multi(kid,2,nobs,2) < 1.e9 .and. &
               obs_multi(kid,1,nobs,1) < 1.e9) then
             k=max(klevs-int(obs_multi(kid,1,nobs,1)/10000.),1)
             ocount(k,id)=ocount(k,id)+1
             do j=1,2
               omean(k,id,j)=omean(k,id,j)+obs_multi(kid,2,nobs,j)
             enddo
           endif
         enddo
       endif
     enddo
!
   elseif (trim(dtype) == 'GPSRO') then
!
     klevs=240
     dlevs=60000./(klevs-1)
     do n=1,nobs_multi_2(1) 
       olevs=int(obs_ids(5,n,1))
       id=int(obs_ids(7,n,1))
       nobs=int(obs_ids(6,n,1))
       do kid=1,olevs
         if (obs_multi(kid,2,nobs,1) < 1.e9 .and. &
             obs_multi(kid,2,nobs,2) < 1.e9 .and. &
             obs_multi(kid,1,nobs,1) < 1.e9) then
           k=min(1+int(obs_multi(kid,1,nobs,1)/dlevs),klevs) ! impact-loc_radius
           ocount(k,id)=ocount(k,id)+1
           do j=1,2
             omean(k,id,j)=omean(k,id,j)+obs_multi(kid,2,nobs,j)
           enddo
         endif
       enddo
     enddo    
!
   else   ! RAD TYPE
!
     olevs=int(obs_ids(5,1,1))
     klevs=olevs
     do n=1,nobs_multi_2(1) 
       id=int(obs_ids(7,n,1))
       nobs=int(obs_ids(6,n,1))
       do k=1,klevs 
         if (obs_multi(k,2,nobs,1) < 1.e9 .and. obs_multi(k,2,nobs,2) < 1.e9) then
           ocount(k,id)=ocount(k,id)+1
           do j=1,2
             omean(k,id,j)=omean(k,id,j)+obs_multi(k,2,nobs,j)
           enddo
         endif
       enddo
     enddo
!
   endif
!
! Replace sums of obs values by mean values
   do id=1,count_types_num
     do k=1,klevs
       if (ocount(k,id) > 0) then
         do j=1,2
           omean(k,id,j)=omean(k,id,j)/ocount(k,id)
         enddo
       endif
     enddo
   enddo
!
! Compute sums of squared deviations from mean 
   if (trim(dtype) == 'PREPBUFR') then
!
     do n=1,nobs_multi_2(1)+nobs_sing_2(1)
       olevs=int(obs_ids(5,n,1))
       id=int(obs_ids(7,n,1))
       nobs=int(obs_ids(6,n,1))
       if (olevs == 1) then
         if (obs_sing(2,nobs,1) < 1.e9 .and. &
             obs_sing(2,nobs,2) < 1.e9 .and. &
             obs_sing(1,nobs,1) < 1.e9) then
           k=max(10-int(obs_sing(1,nobs,1)/10000.),1)
           do j=1,2
             ovars(k,id,j)=ovars(k,id,j)+(obs_sing(2,nobs,j)-omean(k,id,j))**2
           enddo
         endif
       else  
         do kid=1,olevs
           if (obs_multi(kid,2,nobs,1) < 1.e9 .and. &
               obs_multi(kid,2,nobs,2) < 1.e9 .and. &
               obs_multi(kid,1,nobs,1) < 1.e9) then
             k=max(10-int(obs_multi(kid,1,nobs,1)/10000.),1)
             do j=1,2
               ovars(k,id,j)=ovars(k,id,j)+ &
                             (obs_multi(kid,2,nobs,j)-omean(k,id,j))**2
             enddo
           endif 
         enddo
       endif
     enddo
!
   elseif (trim(dtype) == 'GPSRO') then
!
     do n=1,nobs_multi_2(1) 
       olevs=int(obs_ids(5,n,1))
       id=int(obs_ids(7,n,1))  
       nobs=int(obs_ids(6,n,1))
       do kid=1,olevs
         if (obs_multi(kid,2,nobs,1) < 1.e9 .and. &
             obs_multi(kid,2,nobs,2) < 1.e9 .and. &
             obs_multi(kid,1,nobs,1) < 1.e9) then
           k=min(1+int(obs_multi(kid,1,nobs,1)/dlevs),klevs)
           do j=1,2
             ovars(k,id,j)=ovars(k,id,j)+ &
                           (obs_multi(kid,2,nobs,j)-omean(k,id,j))**2
           enddo
         endif 
       enddo
     enddo    
!
   else   ! RAD TYPE
!
     do n=1,nobs_multi_2(1) 
       id=int(obs_ids(7,n,1))
       nobs=int(obs_ids(6,n,1))
       do k=1,klevs
         if (obs_multi(k,2,nobs,1) < 1.e9 .and. &
             obs_multi(k,2,nobs,2) < 1.e9) then
           do j=1,2
             ovars(k,id,j)=ovars(k,id,j)+ &
                           (obs_multi(k,2,nobs,j)-omean(k,id,j))**2
           enddo
         endif
       enddo
     enddo
!
   endif
!
! Replace sums of squared deviations from mean with standard deviations
   do id=1,count_types_num
     do k=1,klevs
       if (ocount(k,id) > 0) then
         do j=1,2
           ovars(k,id,j)=ovars(k,id,j)/ocount(k,id)
           if (ovars(k,id,j) > 0.) then
             ovars(k,id,j)=sqrt(ovars(k,id,j))
           endif
         enddo
       endif
     enddo
   enddo
!
! Print results table
   do id=1,count_types_num
     print *,' '
     print ('(3a,i4,a,i5)'),'Means and stdv for type=',trim(dtype), &
             ' and subtype=',count_types(id,1),' with num=', &
             count_types(id,2)
     print ('(a,4(7x,a8))'),'lev  count','fld mean','fld stdv', &
                'dif mean','dif stdv' 
     do k=1,klevs
       print ('(i3,i7,1p4e15.3)'),k,ocount(k,id),omean(k,id,1), &
             ovars(k,id,1),omean(k,id,2),ovars(k,id,2)
     enddo
   enddo
!
   end subroutine comp_stats_meanvar
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine comp_stats_chan_corr (dtype,j)
!
   implicit none
!
   integer, intent(in) :: j
   character(len=*), intent(in) :: dtype
!
   integer :: olevs
   integer :: nprnt, kskip
   integer :: n,k,k1,k2,id,nobs
   integer :: kprnt(20)
   real(rkind1) :: ocorr(obs_multi_dim1,obs_multi_dim1,obs_types_max)
   real(rkind1) :: var(obs_multi_dim1)
   real(rkind1) :: var2
!
   ocorr(:,:,:)=0.
   olevs=int(obs_ids(5,1,1))
!
   do n=1,nobs_multi_2(1) 
     id=int(obs_ids(7,n,1))
     nobs=int(obs_ids(6,n,1))
     do k1=1,olevs 
       do k2=1,olevs 
        if (obs_multi(k1,2,nobs,j) < 1.e9 .and. obs_multi(k2,2,nobs,j) < 1.e9) then
          ocorr(k1,k2,id)=ocorr(k1,k2,id)+(obs_multi(k1,2,nobs,j)-omean(k1,id,j))* &
                                          (obs_multi(k2,2,nobs,j)-omean(k2,id,j))
        endif
       enddo
     enddo
   enddo
!
   do id=1,count_types_num
!
     do k=1,olevs
       var(k)=ocorr(k,k,id)  ! variance * number of accepted obs
     enddo
!
     do k1=1,olevs
       do k2=1,olevs 
         var2=var(k1)*var(k2) 
         if (var2 > 0.) then
           ocorr(k1,k2,id)=ocorr(k1,k2,id)/sqrt(var2)
         else
           ocorr(k1,k2,id)=0.
         endif
       enddo
     enddo
   enddo
!
   if (trim(dtype) == 'AIRS') then
     nprnt=20
     kprnt=(/10,21,31,47,67,76,85,87,90,95,100,105,110,115,121,166,172,176,193,202 /)
   elseif (trim(dtype) == 'IASI') then
     nprnt=20
     kprnt=(/11,32,57,68,76,88,95,104,125,129,131,145,151,155,161,166,180,186,200,211 /)
   else
     kskip=1+(olevs-1)/20
     nprnt=0
     do k=1,olevs,kskip
       nprnt=nprnt+1
       kprnt(nprnt)=k
     enddo
   endif
!       
   do id=1,count_types_num
     print *,' '
     print *,'Channel correlations for type=',trim(dtype),' and subtype=', &
              count_types(id,1),' with num=',count_types(id,2)
     print ('(3x,20i6)'),kprnt(1:nprnt)
     do k1=1,olevs
       print ('(i3,20f6.2)'),k1,(ocorr(k1,kprnt(k2),id),k2=1,nprnt)  
     enddo
   enddo
!
   end subroutine comp_stats_chan_corr
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine comp_stats_hcorr_rad (dtype,j)
!
   implicit none
!
   integer, intent(in) :: j
   character(len=*), intent(in) :: dtype
!
   integer :: olevs
   integer :: nprnt, kskip
   integer :: k,m,n1,n2,id,id1,id2,nobs1,nobs2,ib
   integer :: nbins
   integer, allocatable :: ibins(:,:,:)
!
   real(rkind1), parameter :: dbin=30.   
   real(rkind1), parameter :: dmax=901.   
   real(rkind1) :: earthr2
   real(rkind2) :: pifac
   real(rkind2) :: rlat1,rlat2,rlon1,rlon2
   real(rkind2) :: d1,d2,d3,d
   real(rkind2) :: obs1,obs2
   real(rkind2) :: xm(5), xa
   real(rkind2), allocatable :: xc(:,:,:)
   real(rkind2), allocatable :: xbins(:,:,:,:)
!
   pifac=atan(1.d0)/45.d0
   earthr2=2.*earthr/1000.   ! 2 x earth radius in km
!
   nbins=1+int(dmax/dbin)
   allocate (ibins(nbins,obs_dim1,obs_dim2))
   allocate (xbins(nbins,5,obs_dim1,obs_dim2))
   allocate (xc(nbins,obs_dim1,obs_dim2))
   ibins(:,:,:)=0
   xbins(:,:,:,:)=0.d0
   xc(:,:,:)=0.d0
!
   olevs=int(obs_ids(5,1,1))
!
   do n1=1,nobs_multi_2(1) 
     id1=int(obs_ids(7,n1,1))
     nobs1=int(obs_ids(6,n1,1))
     rlat1=obs_ids(2,n1,1)*pifac
     rlon1=obs_ids(3,n1,1)*pifac
     do n2=1,nobs_multi_2(1) 
       id2=int(obs_ids(7,n2,1))
       nobs2=int(obs_ids(6,n2,1))
       rlat2=obs_ids(2,n2,1)*pifac
       rlon2=obs_ids(3,n2,1)*pifac
       if (id1 == id2) then   ! same subtype 
!
         if (nobs1 == nobs2) then
           ib=1
         else
           d1=sin(0.5*(rlat1-rlat2))
           d2=sin(0.5*(rlon1-rlon2))
           d3=cos(rlat1)*cos(rlat2)
           d=earthr2*asin(sqrt(d1*d1+d2*d2*d3))
           ib=2+d/dbin
         endif
!        
         if (ib <= nbins) then
           do k=1,olevs
             obs1=obs_multi(k,2,nobs1,j)
             obs2=obs_multi(k,2,nobs2,j)
             if (obs1 < 1.e9 .and. obs2 < 1.e9) then
               xbins(ib,1,k,id1)=xbins(ib,1,k,id1)+obs1*obs2
               xbins(ib,2,k,id1)=xbins(ib,2,k,id1)+obs1
               xbins(ib,3,k,id1)=xbins(ib,3,k,id1)+obs1*obs1
               xbins(ib,4,k,id1)=xbins(ib,4,k,id1)+obs2
               xbins(ib,5,k,id1)=xbins(ib,5,k,id1)+obs2*obs2
               ibins(ib,k,id1)=ibins(ib,k,id1)+1               
             endif
           enddo
         endif
!
       endif  ! test on subtype same
     enddo    ! loop over n2
   enddo      ! loop over n1
!
   do id=1,count_types_num
     do ib=1,nbins
       do k=1,olevs
         if (ibins(ib,k,id) > 0) then 
           do m=1,5         
             xm(m)=xbins(ib,m,k,id)/real(ibins(ib,k,id),8)
           enddo
           xa=(xm(3)-xm(2)*xm(2))*(xm(5)-xm(4)*xm(4)) 
           if (xa > 0.d0) then
             xc(ib,k,id)=(xm(1)-xm(2)*xm(4))/sqrt(xa)
           else
             xc(ib,k,id)=0.d0  
           endif
         endif
       enddo
     enddo
   enddo
!
   do id=1,count_types_num
     print *,' '
     print ('(3a,i4,a,i6,a,f4.0)'),'Horizontal correlations for type=', &
            trim(dtype),' and subtype=',count_types(id,1), &
            ' with num=',count_types(id,2),' and dbin=',dbin
     do k=1,olevs
       print ('(i3,20f6.2)'),k,(xc(ib,k,id),ib=1,min(nbins,20))
       if (nbins > 20) then
         print ('(3x,20f6.2)'),k,(xc(ib,k,id),ib=21,nbins) 
       endif
     enddo
   enddo
!
   end subroutine comp_stats_hcorr_rad
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine comp_stats_hcorr_conv1 (dtype,j)
!
! compute horiz correl of single-level conv obs
!
   implicit none
!
   integer, intent(in) :: j
   character(len=*), intent(in) :: dtype
!
   integer, parameter :: klevs=21 
   integer :: olevs1,olevs2
   integer :: nprnt, kskip
   integer :: k,k1,k2,m,n1,n2,id,id1,id2,nobs1,nobs2,ib
   integer :: nbins
   integer, allocatable :: ibins(:,:,:)
!
   real(rkind1), parameter :: dbin=30.   
   real(rkind1), parameter :: dmax=901.
   real(rkind1), parameter :: dp_tol=500.  ! tolerance in Pa
   real(rkind1) :: delp,xp  
   real(rkind1) :: earthr2
   real(rkind2) :: pifac
   real(rkind2) :: rlat1,rlat2,rlon1,rlon2
   real(rkind2) :: d1,d2,d3,d
   real(rkind2) :: obs1,obs2
   real(rkind2) :: xm(5), xa
   real(rkind2), allocatable :: xc(:,:,:)
   real(rkind2), allocatable :: xbins(:,:,:,:)
!
   pifac=atan(1.d0)/45.d0
   earthr2=2.*earthr/1000.   ! 2 x earth radius in km
   delp=1.0e5/(klevs-1)
!
   nbins=1+int(dmax/dbin)
   allocate (ibins(nbins,klevs,obs_dim2))
   allocate (xbins(nbins,5,klevs,obs_dim2))
   allocate (xc(nbins,klevs,obs_dim2))
   ibins(:,:,:)=0
   xbins(:,:,:,:)=0.d0
   xc(:,:,:)=0.d0
!
   do n1=1,nobs_multi_2(1)+nobs_sing_2(1)
     id1=int(obs_ids(7,n1,1))
     nobs1=int(obs_ids(6,n1,1))
     rlat1=obs_ids(2,n1,1)*pifac
     rlon1=obs_ids(3,n1,1)*pifac
     olevs1=int(obs_ids(5,n1,1))
!
     if (olevs1 == 1) then 
       if (obs_sing(2,nobs1,2) < 1.e9 .and. obs_sing(1,nobs1,1) < 1.e9) then
         k1=1+int(obs_sing(1,nobs1,1)/delp)                
         xp=delp*(k1-1)
         if (abs(obs_sing(1,nobs1,1)-xp) <= dp_tol) then  
!
     do n2=1,nobs_multi_2(1)+nobs_sing_2(1) 
       id2=int(obs_ids(7,n2,1))
       nobs2=int(obs_ids(6,n2,1))
       rlat2=obs_ids(2,n2,1)*pifac
       rlon2=obs_ids(3,n2,1)*pifac
       olevs2=int(obs_ids(5,n2,1))
       if (id1 == id2) then     ! same subtype 
!
     if (olevs2 == 1) then 
       if (obs_sing(2,nobs1,2) < 1.e9 .and. obs_sing(1,nobs2,2) < 1.e9) then
         k2=1+int(obs_sing(1,nobs2,1)/delp)                
         xp=delp*(k2-1)
         if (abs(obs_sing(1,nobs2,1)-xp) <= dp_tol) then  
           if (k1 == k2) then   ! same level range
!
         if (nobs1 == nobs2) then
           ib=1
         else
           d1=sin(0.5*(rlat1-rlat2))
           d2=sin(0.5*(rlon1-rlon2))
           d3=cos(rlat1)*cos(rlat2)
           d=earthr2*asin(sqrt(d1*d1+d2*d2*d3))
           ib=2+d/dbin
         endif
!        
         if (ib <= nbins) then
           do k=k1,k1
             obs1=obs_sing(2,nobs1,j)
             obs2=obs_sing(2,nobs2,j)
             if (obs1 < 1.e9 .and. obs2 < 1.e9) then
               xbins(ib,1,k,id1)=xbins(ib,1,k,id1)+obs1*obs2
               xbins(ib,2,k,id1)=xbins(ib,2,k,id1)+obs1
               xbins(ib,3,k,id1)=xbins(ib,3,k,id1)+obs1*obs1
               xbins(ib,4,k,id1)=xbins(ib,4,k,id1)+obs2
               xbins(ib,5,k,id1)=xbins(ib,5,k,id1)+obs2*obs2
               ibins(ib,k,id1)=ibins(ib,k,id1)+1               
             endif
           enddo
         endif
!
                endif ! test on p values for both obs in same level range
              endif   ! test on p in some range of tolerance
            endif     ! test on obs < bmiss
          endif       ! test on olevs2
        endif         ! test on subtype same
      enddo           ! loop over n2
 
        endif   ! test on p in some range of tolerance
      endif     ! test on obs < bmiss
    endif       ! test on olevs1
  enddo         ! loop over n1
!
   do id=1,count_types_num
     do ib=1,nbins
       do k=1,klevs
         if (ibins(ib,k,id) > 0) then 
           do m=1,5         
             xm(m)=xbins(ib,m,k,id)/real(ibins(ib,k,id),8)
           enddo
           xa=(xm(3)-xm(2)*xm(2))*(xm(5)-xm(4)*xm(4)) 
           if (xa > 0.d0) then
             xc(ib,k,id)=(xm(1)-xm(2)*xm(4))/sqrt(xa)
           else
             xc(ib,k,id)=0.d0  
           endif
         endif
       enddo
     enddo
   enddo
!
   do id=1,count_types_num
     print *,' '
     print ('(3a,i4,a,i6,a,f4.0)'),'Horizontal correlations for type=', &
            trim(dtype),' and subtype=',count_types(id,1), &
            ' with num=',count_types(id,2),' and dbin=',dbin
     print *,'Note: The only non-zero values will be for single-level reports'
     do k=1,klevs
       xp=delp*(k-1)/100.  ! level p in mb
       print ('(i3,f6.0,i6,20f6.2)'),k,xp,ibins(1,k,id), &
                                  (xc(ib,k,id),ib=1,min(nbins,20))
       if (nbins > 20) then
         print ('(15x,20f6.2)'),(xc(ib,k,id),ib=21,nbins) 
       endif
     enddo
   enddo
!
   end subroutine comp_stats_hcorr_conv1
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine comp_stats_vcorr (dtype,j)
!
! compute vert correl of multi-level conv obs
!
   implicit none
!
   integer, intent(in) :: j
   character(len=*), intent(in) :: dtype
!
   integer, parameter :: klevs=51 
   integer :: olevs1,olevs2
   integer :: nprnt, kskip
   integer :: k,k1,k2,ka,kb,m,n1,id,id1,nobs1
   integer :: otypes_max,otypes_1,otypes_2
   integer, allocatable :: ibins(:,:,:)
!
   real(rkind1) :: dlev,xlev1,xlev2,slev,flev
   real(rkind2) :: obs1,obs2,obsx
   real(rkind2) :: xm(5), xa, xp
   real(rkind2), allocatable :: xc(:,:,:)
   real(rkind2), allocatable :: xbins(:,:,:,:)
!
   if (trim(dtype) == 'PREPBUFR') then
     dlev=1.0e5/(klevs-1)
     slev=0.
     flev=0.01  ! changes Pa to hPa
     otypes_max=count_types_num
   else
     dlev=6.0e4/(klevs-1)
     slev=0.
     flev=0.001 ! changes m to km
     otypes_max=1 
   endif
!
   allocate (ibins(klevs,klevs,otypes_max))
   allocate (xbins(klevs,klevs,5,otypes_max))
   allocate (xc(klevs,klevs,otypes_max))
   ibins(:,:,:)=0
   xbins(:,:,:,:)=0.d0
   xc(:,:,:)=0.d0
!
   do n1=1,nobs_all_2(1)
     if (trim(dtype) == 'PREPBUFR') then
       id1=int(obs_ids(7,n1,1))
     else
       id1=1   ! treat all subtypes together to increase counts
     endif
     nobs1=int(obs_ids(6,n1,1))
     olevs1=int(obs_ids(5,n1,1))
     if (olevs1 > 1) then

     do k1=1,olevs1
       obsx=obs_multi(k1,2,nobs1,1)
       obs1=obs_multi(k1,2,nobs1,2)
       xlev1=obs_multi(k1,1,nobs1,1)
       if ( obs1 < 1.e9 .and.  obsx < 1.e9 .and. xlev1 < 1.e9) then
         ka=1+int((xlev1-slev)/dlev)    
         ka=min(ka,klevs)            
         ka=max(ka,1)            
         do k2=1,olevs1
           obsx=obs_multi(k2,2,nobs1,1)
           obs2=obs_multi(k2,2,nobs1,2)
           xlev2=obs_multi(k2,1,nobs1,1)
           if ( obs2 < 1.e9 .and. obsx < 1.e9 .and. xlev2 < 1.e9) then
             kb=1+int((xlev2-slev)/dlev)                
             kb=min(kb,klevs)            
             kb=max(kb,1)            
!            
             xbins(kb,ka,1,id1)=xbins(kb,ka,1,id1)+obs1*obs2
             xbins(kb,ka,2,id1)=xbins(kb,ka,2,id1)+obs1
             xbins(kb,ka,3,id1)=xbins(kb,ka,3,id1)+obs1*obs1
             xbins(kb,ka,4,id1)=xbins(kb,ka,4,id1)+obs2
             xbins(kb,ka,5,id1)=xbins(kb,ka,5,id1)+obs2*obs2
             ibins(kb,ka,id1)=ibins(kb,ka,id1)+1               
           endif
         enddo   ! loop over k2
       endif    
     enddo       ! loop over k1
   endif         ! test on whther muti-level obs
   enddo         ! loop over n1
!
   do id=1,otypes_max
     do kb=1,klevs
       do ka=1,klevs
         if (ibins(kb,ka,id) > 0) then 
           do m=1,5         
             xm(m)=xbins(kb,ka,m,id)/real(ibins(kb,ka,id),8)
           enddo
           xa=(xm(3)-xm(2)*xm(2))*(xm(5)-xm(4)*xm(4)) 
           if (xa > 0.d0) then
             xc(kb,ka,id)=(xm(1)-xm(2)*xm(4))/sqrt(xa)
           else
             xc(kb,ka,id)=0.d0  
           endif
         endif
       enddo
     enddo
   enddo
!
   do id=1,otypes_max
     print *,' '
     if (trim(dtype) == 'PREPBUFR') then
       otypes_1=count_types(id,1) 
       otypes_2=count_types(id,2)
     else
       otypes_1=9999
       otypes_2=sum(count_types(:,2))
     endif
!
     print ('(3a,i4,a,i6,a,f6.0)'),'Vertical correlations for type=',   &
            trim(dtype),' and subtype=',otypes_1,' with num=',otypes_2, & 
            ' and dlev=',dlev
     print *,'Note: The only non-zero values will be for multi-level reports'
     do k=1,klevs
       xp=dlev*(k-1)*flev   ! level p in mb or h in km
       print ('(i3,f8.2,i7,20f6.2)'),k,xp,ibins(1,k,id), &
                                  (xc(kb,k,id),kb=1,min(klevs,20))
       if (klevs > 20) then
         print ('(18x,20f6.2)'),(xc(kb,k,id),kb=21,klevs) 
       endif
     enddo
   enddo
!
   end subroutine comp_stats_vcorr
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_comp_stats
