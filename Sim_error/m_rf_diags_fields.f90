   module m_rf_diags_fields
!
! Module used to compute diagnostics of correlated random fields, for 
! examination or code testing purposes
!
! Initial Code by Ronald Errico NASA/GMAO Sept. 2014
!
   use m_kinds, only : rkind1, rkind2
   use m_parameters, only : earthr
   use m_shtrans_rf, only : imax,jmax
   use m_random_fields, only : random_fields_get_1_value
   use m_random_fields, only : random_fields_get_1_lev
   use m_random_fields, only :ifields_corr, ilevs_corr 
   use m_obs_error_table, only : et_itypes_corr, et_hcorr_lengths
!
   implicit none
!
   private
   public :: rf_diags_fields_calc
   
!
   real(rkind1), parameter :: earthrkm=earthr*0.001   ! earth radius in km
   real(rkind1) :: hcorr_length1
!
   contains
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rf_diags_fields_calc (dtype,call_name)
!
! Call routine to calculate and print diagnostics of a sample of random fields
!
   implicit none
   integer :: id,kid,nf,n,k
   character(len=*), intent(in) :: dtype
   character(len=*), intent(in) :: call_name
   character(len=80) :: info
!
   print *,' '
   print *,'rf_diags_fields_calc: ',dtype,' ',call_name
   do id=1,et_itypes_corr
     do nf=1,ifields_corr
       kid=nf+(id-1)*ifields_corr
       do n=min(6,ilevs_corr),min(10,ilevs_corr),4
         k=n+(kid-1)*ilevs_corr
         hcorr_length1=et_hcorr_lengths(n,id)
         write(info,'(a,4i4,a,f7.2)'),                             &
                      'Fields_in_setup subtype,nfield,lev,index=', &
                      id,nf,n,k,'   hcorr_length=',hcorr_length1
         call rf_diags_fields_corrs (ilevs_corr,1,k,id,kid,'S',info)
       enddo
     enddo
   enddo 
!
   end subroutine rf_diags_fields_calc
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rf_diags_fields_corrs (nlevs,nfields,kpid,id,kid, &
                                     field_type,info)
!
! Print some sample of information for testing of random fields if desired
!
   implicit none
!
   integer, intent(in) :: nlevs
   integer, intent(in) :: nfields
   integer, intent(in) :: kpid, kid
   integer, intent(in) :: id
   character(len=*), intent(in) :: field_type
   character(len=*), intent(in) :: info
!
   logical :: hc_same
   integer :: nlons,nlats
   integer :: nbins
   integer :: nlevs10,k1,k2
   integer :: i_region(imax),j_region(2)
   integer, allocatable :: ibins(:,:)
!
   real(rkind1) :: max_bin=1500. ! max dx to consider for bins to determine corr
   real(rkind1) :: dx_bin
   real(rkind1) :: two_pi
   real(rkind1) :: xlon1, xlon2, xlat1, xlat2
   real(rkind1) :: dlon,dlat
   real(rkind1) :: fmean,fvar
   real(rkind1) :: fld1(9) 
   real(rkind1), allocatable :: fld2(:,:,:) 
   real(rkind2), allocatable :: xbins(:,:,:,:) 
   real(rkind1), allocatable :: hcorr(:,:,:)
!
   integer :: k 
   real(rkind1), allocatable :: hcorr_sum(:,:,:) 
!
   print *,' '
   print ('(a,i4)'),'Diagnostics for random horriz. correl. for kpid=', &
                     kpid
   print ('(a)'),'  ',info 
!
! Check pole values
   print *,' Check that fields at poles are specified properly: SP,NP,EQ'
   call rf_diags_fields_1_lat (9,kpid,-90.,0.,360.,fld1)
   print ('(9f12.6)'),fld1 
   call rf_diags_fields_1_lat (9,kpid, 90.,0.,360.,fld1)
   print ('(9f12.6)'),fld1
   call rf_diags_fields_1_lat (9,kpid, 0.,0.,360.,fld1)
   print ('(9f12.6)'),fld1
!  
! print portion of sample of fields
   nlons=24
   nlats=24
   dlon=360./imax
   dlat=180./(jmax-1)
   xlon1=dlon
   xlon2=dlon*nlons
   xlat1=0.
   xlat2=dlat*(nlats-1)
   allocate (fld2(nlons,nlats,1))
   call rf_diags_fields_get_part (nlons,nlats,kpid,xlon1,xlon2, &
                                  xlat1,xlat2,fld2(:,:,1))
   call rf_diags_fields_print_f (nlons,nlats,kpid,xlon1,xlon2,  &
                                  xlat1,xlat2,fld2(:,:,1))
   deallocate (fld2) 
!
! print sample of correlations for portion of field
   nlons=imax
   nlats=jmax
   dlon=360./imax
   dlat=180./(jmax-1)
   xlon1=0.
   xlon2=dlon*(nlons-1)
   xlat1=-90.
   xlat2=90.
   i_region(:)=1
   j_region(1)=jmax/3
   j_region(2)=2*jmax/3
   two_pi=8.*atan(1.) 
   dx_bin=dlat*earthrkm*two_pi/360.
   nbins=1+nint(max_bin/dx_bin)
   allocate (fld2(nlons,nlats,1))
   allocate (xbins(nbins,5,1,2))
   allocate (ibins(nbins,2))
   allocate (hcorr(nbins,1,2))
   allocate (hcorr_sum(nbins,1,2))
   xbins(:,:,:,:)=0.d0
   ibins(:,:)=0
   call random_fields_get_1_lev (kpid,fld2(:,:,1))
   call rf_diags_fields_int_var (nlons,nlats,xlon1,xlon2,xlat1,xlat2, &
                                 fld2(:,:,1),fmean,fvar) 
   print *,' '
   print ('(a)'),'Statistics for a field on a single level follow:'
   print *,' '
   print ('(a)'),'Mean and variance of field integrated over region:'
   print ('(a,4f10.4)'),'xlon1,xlon2,xlat1,xlat2=',xlon1,xlon2,xlat1,xlat2
   print ('(a,i3,a)'),'kpid=',kpid,' ',info
   print ('(2(a,f10.6))'),'mean=',fmean,'   var=',fvar
   call rf_diags_fields_sump (imax,jmax,1,nbins,i_region,j_region, &
                              ibins(:,1),ibins(:,2), &
                              dx_bin,fld2,xbins(:,:,:,1),xbins(:,:,:,2))
   call rf_diags_fields_compc (nbins,1,ibins(:,1),xbins(:,:,:,1),hcorr(:,:,1))
   call rf_diags_fields_compc (nbins,1,ibins(:,2),xbins(:,:,:,2),hcorr(:,:,2))
   call rf_diags_fields_prntc (nbins,1,kpid,dx_bin,ibins,hcorr,info)
!
!  Only compute the following average over the 1st 10 levels if all 
!  et_hcorr_lengths for these are the same
   nlevs10=min(nlevs,10)
   hc_same=.true.
   do k=2,nlevs10
     if (abs(et_hcorr_lengths(k,id)-et_hcorr_lengths(1,id)) > .01 ) then
       hc_same=.false.
     endif
   enddo 
!
   if (hc_same) then 
     hcorr_sum(:,:,:)=0.
     xbins(:,:,:,:)=0.d0
     ibins(:,:)=0
     k1=(kid-1)*nlevs+1 ! starts at first lev for this subtype group and field
     k2=k1+nlevs10-1
!
     do k=k1,k2
       call random_fields_get_1_lev (k,fld2(:,:,1)) 
       call rf_diags_fields_sump (imax,jmax,1,nbins,i_region,j_region, &
                                ibins(:,1),ibins(:,2), &
                                dx_bin,fld2,xbins(:,:,:,1),xbins(:,:,:,2))
       call rf_diags_fields_compc (nbins,1,ibins(:,1),xbins(:,:,:,1), &
                                   hcorr(:,:,1))
       call rf_diags_fields_compc (nbins,1,ibins(:,2),xbins(:,:,:,2), & 
                                   hcorr(:,:,2))
       hcorr_sum(:,:,:)=hcorr_sum(:,:,:)+hcorr(:,:,:)
     enddo
     hcorr_sum(:,:,:)=hcorr_sum(:,:,:)/nlevs10
     print *,' '
     print ('(a,i3,a,f7.2,a)'), &
             'Statistics averaged for a field over first',nlevs10, &
             ' levels present follow (hcorr_lengths=',             &
             et_hcorr_lengths(1,id),'):' 
     call rf_diags_fields_prntc (nbins,1,999,dx_bin,ibins,hcorr_sum,info)
   endif
!
   deallocate (fld2,xbins,ibins,hcorr,hcorr_sum)
!
! check interpolation at a sample of lats and lons
   nlons=5
   nlats=5
   allocate (fld2(nlons,nlats,1))
   dlon=360./imax
   dlat=180./(jmax-1)
!
   xlon1=0.
   xlon2=xlon1+(nlons-1)*dlon
   xlat1=-dlat
   xlat2=xlat1+(nlats-1)*dlat
   call rf_diags_fields_get_part (nlons,nlats,kpid,xlon1,xlon2, &
                                  xlat1,xlat2,fld2(:,:,1))
   call rf_diags_fields_print_f (nlons,nlats,kpid,xlon1,xlon2, &
                                 xlat1,xlat2,fld2(:,:,1))
   xlon1=dlon*0.2
   xlon2=xlon1+(nlons-1)*dlon
   xlat1=-dlat*0.8
   xlat2=xlat1+(nlats-1)*dlat
   call rf_diags_fields_get_part (nlons,nlats,kpid,xlon1,xlon2, &
                                  xlat1,xlat2,fld2(:,:,1))
   call rf_diags_fields_print_f (nlons,nlats,kpid,xlon1,xlon2,  &
                                 xlat1,xlat2,fld2(:,:,1))
!       
   deallocate (fld2)
!
   print ('(a)'),' Testing of horizontal correlations completed'
!
   end subroutine rf_diags_fields_corrs
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rf_diags_fields_int_var (nlons,nlats,xlon1,xlon2, &
                                       xlat1,xlat2,fld,fmean,fvar) 
!
! Integrate the square of a field over the area of the globe
! 
   integer, intent(in) :: nlons,nlats
   real(rkind1), intent(in) :: xlon1,xlon2,xlat1,xlat2
   real(rkind1), intent(in) :: fld(nlons,nlats)
   real(rkind1), intent(out) :: fmean,fvar
!
   integer :: i,j
   real (rkind1) :: dlat,dlon,area,fd
   real (rkind1) :: slat,slat1,slat2,pifac
   real (rkind1) :: gbaf(nlats), lats(nlats)
!
! Set grid lats
   dlat=(xlat2-xlat1)/(nlats-1)
   do j=1,nlats
     lats(j)=xlat1+(j-1)*dlat
   enddo
!
! Set grid box area factor (i.e., 'average' cosine(lat) factor)
   pifac=atan(1.)/45.
   do j=2,nlats-1
     slat1=sin(0.5*(lats(j)+lats(j-1))*pifac)
     slat2=sin(0.5*(lats(j)+lats(j+1))*pifac)
     gbaf(j)=(slat2-slat1)/dlat
     if (j==2) then
       slat=max(-90.,lats(1)-dlat*0.5)
       slat=sin(slat*pifac)
       gbaf(1)=(slat1-slat)/dlat
     endif
     if (j==nlats-1) then
       slat=min(90.,lats(nlats)+dlat*0.5)
       slat=sin(slat*pifac)
       gbaf(nlats)=(slat-slat2)/dlat
     endif
   enddo
   area=sum(gbaf(1:nlats))
   gbaf(:)=gbaf(:)/(area*nlons)
! 
   fmean=0.
   do j=1,nlats
     do i=1,nlons
       fmean=fmean+gbaf(j)*fld(i,j)
     enddo
   enddo
!
   fvar=0.
   do j=1,nlats
     do i=1,nlons
       fd=fld(i,j)-fmean
       fvar=fvar+gbaf(j)*fd*fd
     enddo
   enddo
!
   end subroutine rf_diags_fields_int_var        
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rf_diags_fields_get_part (nlons,nlats,kpid,xlon1,xlon2, &
                                        xlat1,xlat2,fld)
!
!  Get portion of field in desired range
!
   implicit none
   integer, intent(in) :: nlons,nlats,kpid
   real(rkind1), intent(in) :: xlon1,xlon2,xlat1,xlat2
   real(rkind1), intent(out) :: fld(nlons,nlats)
!
   integer :: nlat
   real(rkind1) :: dlat,xlat
!
   dlat=abs(xlat2-xlat1)/max(nlats-1,1)
   do nlat=1,nlats
     xlat=min(xlat2,xlat1)+(nlat-1)*dlat
     call rf_diags_fields_1_lat (nlons,kpid,xlat,xlon1,xlon2,fld(:,nlat)) 
   enddo  
!
   end subroutine rf_diags_fields_get_part 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rf_diags_fields_print_f (nlons,nlats,kpid,xlon1,xlon2, &
                                       xlat1,xlat2,fld)
!
! Print a portion of the random fields
!
   implicit none
!
   integer, intent(in) :: nlons,nlats,kpid
   real(rkind1), intent(in) :: xlon1,xlon2,xlat1,xlat2
   real(rkind1), intent(in) :: fld(nlons,nlats)
!
   integer :: i, j
   integer :: ifac
   real(rkind1) :: fmax,sfac
   integer :: ifield(nlons)
!
! determine scaling
   fmax=1.e-20
   do j=1,nlats
     do i=1,nlons
       fmax=max(fmax,abs(fld(i,j)))
     enddo 
   enddo  
   ifac=3-int(log10(fmax))
   sfac=10.**ifac
!
   print *,' '
   print *,'Sample field for nlons,nlats,kpid,xlon1,xlon2,xlat1,xlat2 ='
   print ('(3i5,4f8.2)'),nlons,nlats,kpid,xlon1,xlon2,xlat1,xlat2
!             
   do j=nlats,1,-1
     do i=1,nlons     
       ifield(i)=int(sfac*fld(i,j))
     enddo
     print ('(25i5)'),j,ifield
   enddo
!
   end subroutine rf_diags_fields_print_f
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rf_diags_fields_1_lat (nlons,kpid,xlat,xlon1,xlon2,f1lat)
! 
! Get values at a sample of points at a single latitude
! 
   integer, intent(in) :: nlons,kpid
   real(rkind1), intent(in) :: xlat,xlon1,xlon2
   real(rkind1), intent(out) :: f1lat(nlons)
!
   integer :: nlon
   real(rkind1) :: xlon,dlon   
!
   dlon=abs(xlon2-xlon1)/max(nlons-1,1)     
   do nlon=1,nlons
     xlon=min(xlon1,xlon2)+(nlon-1)*dlon
     call random_fields_get_1_value (xlon,xlat,kpid,f1lat(nlon))
   enddo
!
   end subroutine rf_diags_fields_1_lat 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rf_diags_fields_prntc (nbins,nlevs,kpid,dx_bin,ibins, &
                                     hcorr,info) 
!
!  Print binned NS or EW horizontal correlations as function of distance
!  Optionally, also print equivalent values for some standard functions
!  (Gaussian, exponential, third-order auto-regressive)
!
   implicit none
!
   integer, intent(in) :: nbins, nlevs
   integer, intent(in) :: kpid
   integer, intent(in) :: ibins(nbins,2)
   real(rkind1), intent(in) :: dx_bin
   real(rkind1), intent(in) :: hcorr(nbins,nlevs,2)
   character(len=*), intent(in) :: info
!
   logical, parameter :: lprint_shapes=.true.
   integer :: n
   integer :: ibins100(nbins)
   real(rkind1) :: x
   real(rkind1) :: ref_corr(nbins)
!
   print *,' ' 
   print ('(a)'),'Horizontal correlations as a function of distance'
   print ('(a,f6.0,a,i3,a)'),'dx_bin=',dx_bin,'  kpid=',kpid,' ',info
   print ('(a)'),'ibins/100 first, then NS corr,then EW' 
   ibins100(:)=ibins(:,1)/100
   print ('(20i6)'),ibins100(:)
   print ('(20f6.2)'),hcorr(:,1,1)
   ibins100(:)=ibins(:,2)/100
   print ('(20i6)'),ibins100(:)
   print ('(20f6.2)'),hcorr(:,1,2)
!
   if (lprint_shapes) then ! print sample of shapes for comparison 
     print ('(a,f7.2,a)'),'Gauss function with L=',hcorr_length1,'km' 
     ref_corr(1)=1.
     do n=2,nbins
       x=(n-1.5)*dx_bin/hcorr_length1
       ref_corr(n)=exp(-0.5*x*x)
     enddo
     print ('(20f6.2)'),ref_corr(:)
!
     print ('(a,f7.2,a)'),'Exp function with L=',hcorr_length1,'km' 
     ref_corr(1)=1.
     do n=2,nbins
       x=(n-1.5)*dx_bin/hcorr_length1
       ref_corr(n)=exp(-x)
     enddo
     print ('(20f6.2)'),ref_corr(:)
!
     print ('(a,f7.2,a)'),'TOAR function with L=',hcorr_length1,'km' 
     ref_corr(1)=1.
     do n=2,nbins
       x=(n-1.5)*dx_bin/hcorr_length1
       ref_corr(n)=(1.+x+x*x/3.)*exp(-x)
     enddo
     print ('(20f6.2)'),ref_corr(:)
!
   endif
!
   end subroutine rf_diags_fields_prntc
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine rf_diags_fields_sump (imax,jmax,nlevs,nbins, &
                               i_region,j_region,ns_count,ew_count, &
                               x_bin,fields,ns_corr,ew_corr)
!
!  Compute sums of cross products for points separated in N-S and, 
!  separately, E-W directions, for all pairs within a range of lats 
!  and lons.  The pairs are assigned to bins based on d=s/x_bin, where 
!  s is the separation and x_bin is the prescribed bin resolution:
!  bin 1 for d=0
!  bin 2 for 0<d<1.5
!  bin j+1 for j-0.5<d<j+0.5
!
      implicit none
!
      integer :: imax,jmax,nlevs,nbins
      integer :: i_region(imax),j_region(2)
      integer :: ns_count(nbins), ew_count(nbins)
      real(4) :: x_bin
      real(4) :: fields(imax,jmax,nlevs)
      real(8) :: ns_corr(nbins,5,nlevs), ew_corr(nbins,5,nlevs)
!
      integer :: i1, i2, j1, j2, k, ib, idiff
      real(4) :: alat, alon, dlat, dlon, xlat, xlon
      real(4) :: pifac, afac, d, alatxb, alonxb
!
      dlat=180./(jmax-1)
      dlon=360./imax
!
! sum NS cross products in region
      pifac=atan(1.)/45.
      alat=earthrkm*pifac*dlat
      alatxb=alat/x_bin
      do i1=1,imax
        if (i_region(i1) == 1) then ! consider this longitude
          do j1=j_region(1),j_region(2)
            do j2=j_region(1),j_region(2)
!
              if (j1==j2) then
                ib=1                 !
              else
                d=abs(j1-j2)*alatxb
                ib=max(int(d+1.5),2) 
              endif
!
              if (ib<=nbins) then              
                ns_count(ib)=ns_count(ib)+1
                do k=1,nlevs     
                  ns_corr(ib,1,k)=ns_corr(ib,1,k) + &
                                    fields(i1,j1,k)*fields(i1,j2,k)
                  ns_corr(ib,2,k)=ns_corr(ib,2,k) + fields(i1,j1,k)
                  ns_corr(ib,3,k)=ns_corr(ib,3,k) + &
                                    fields(i1,j1,k)*fields(i1,j1,k)
                  ns_corr(ib,4,k)=ns_corr(ib,4,k) + fields(i1,j2,k)
                  ns_corr(ib,5,k)=ns_corr(ib,5,k) + &
                                    fields(i1,j2,k)*fields(i1,j2,k)
                enddo
              endif
!
            enddo
          enddo
        endif
      enddo
! 
! sum EW cross products in region
      afac=earthrkm*pifac*dlon
      do j1=j_region(1),j_region(2)
        xlat=-90.+dlat*(j1-1)
        alon=afac*cos(xlat*pifac)
        alonxb=alon/x_bin
        do i1=1,imax
          if (i_region(i1) == 1) then ! consider this longitude
            do i2=1,imax
              if (i_region(i2) == 1) then ! consider this longitude
!
                if (i1==i2) then
                  ib=1
                else
                  if (i1 <= i2) then
                    idiff=min(i2-i1, i1+imax-i2)
                  else
                    idiff=min(i1-i2, i2+imax-i1)
                  endif
                  d=alonxb*idiff
                  ib=max(int(d+1.5),2) 
                endif
!
                if (ib<=nbins) then              
                  ew_count(ib)=ew_count(ib)+1
                  do k=1,nlevs          
                    ew_corr(ib,1,k)=ew_corr(ib,1,k) + &
                                    fields(i1,j1,k)*fields(i2,j1,k)
                    ew_corr(ib,2,k)=ew_corr(ib,2,k) + fields(i1,j1,k)
                    ew_corr(ib,3,k)=ew_corr(ib,3,k) + &
                                    fields(i1,j1,k)*fields(i1,j1,k)
                    ew_corr(ib,4,k)=ew_corr(ib,4,k) + fields(i2,j1,k)
                    ew_corr(ib,5,k)=ew_corr(ib,5,k) + &
                                    fields(i2,j1,k)*fields(i2,j1,k)
                  enddo
                endif
!   
              endif
            enddo
          endif
        enddo
      enddo
!
   end subroutine rf_diags_fields_sump
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine rf_diags_fields_compc (nbins,nlevs,ibins,xbins,hcorr)
!
! Compute correlations from sums
!
      implicit none
!
      integer    :: nbins,nlevs
      integer    :: ibins(nbins)
      real(8)    :: xbins(nbins,5,nlevs)
      real(4)    :: hcorr(nbins,nlevs)
!
      integer :: ib, k
      real(8) :: xm1, xm2, xs1, xs2, xa, xb, zero8
      real(4) :: zero       
!
      zero8=0.0_8
      zero=0.0_4
!
      do ib=1,nbins
        if (ibins(ib) > 0) then
          do k=1,nlevs
            xm1=xbins(ib,2,k)/ibins(ib)
            xm2=xbins(ib,4,k)/ibins(ib)
            xs1=xbins(ib,3,k)/ibins(ib)
            xs2=xbins(ib,5,k)/ibins(ib)
            xa=(xs1-xm1*xm1)*(xs2-xm2*xm2)
            if (xa > zero8) then   
              xb=(xbins(ib,1,k)/ibins(ib)-xm1*xm2)/sqrt(xa)
              hcorr(ib,k)=xb
            else
              hcorr(ib,k)=zero
            endif
          enddo
        endif
      enddo 
!
   end subroutine rf_diags_fields_compc
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   end module m_rf_diags_fields
