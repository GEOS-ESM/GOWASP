   module m_random_fields
!
! Module used to create spatially correlated random fields and to get 
! get values of those fields at requested locations.
!
! This is accomplished by:
! (0) reading desired parameters specifying the desired correlations from a
!   file, based on correlation shapes and lengths as functions of level
!   or channel and of sub-type  
! (1) projecting the desired correlation function shape (specified also by 
!   its length scale) onto Legendre polynomials P(m,n) for m=0.
!   (This is done in a separate module)
! (2) modifying the power spectrum from step 1 if random vector wind fields 
!   are desired
! (3) creating random coefficients for spherical harmonics so that the  
!   expected power spectrum is the one defined in step 2.
! (4) creating global random fields partitioned onto several processors by
!   their reconstruction from the spherical harmonic coefficients from step 3
! (5) pulling values for the random fields at desired locations using 
!   vertical and horizontal interpolation. 
! 
! Channel or vertical correlations are applied using EOFs as basis functions, 
! with the random fields definining random coefficients for each EOF.
!
! Initial Code by Ronald Errico NASA/GMAO Sept. 2015
!
   use m_kinds, only : rkind1, rkind2
   use m_shtrans_rf, only : sh_init_factors 
   use m_shtrans_rf, only : sh_init_lats 
   use m_shtrans_rf, only : sh_calc_power
   use m_shtrans_rf, only : sh_trans
   use m_shtrans_rf, only : sh_clean
   use m_shtrans_rf, only : rkindf, rkinds
   use m_shtrans_rf, only : imax,jmax,nspects
   use m_obs_error_table, only : et_itypes_corr, et_l_vc_corr, et_nlevels 
   use m_obs_error_table, only : et_rf_nmax, et_ran_fields_kinds
   use m_obs_error_table, only : et_ks_list, et_hcorr_ks
   use m_obs_error_table, only : et_hcorr_lengths, et_frac_corr, et_corr_shapes
   use m_obs_error_table, only : et_e_vects, et_e_v_sqrt 
   use m_obs_error_table, only : et_pmax, et_pmin
   use m_obs_error_table, only : 
   use m_random_gauss, only : random_gauss_r8 
   use MAPL_ShmemMod    ! The SHMEM infrastructure
!
   implicit none
!
   private
   public :: random_fields_clean 
   public :: random_fields_setup
   public :: random_fields_get_values 
   public :: random_fields_get_1_value 
   public :: random_fields_get_1_lev
!
   integer, public :: ifields_corr ! =et_ran_fields_kinds
   integer, public :: ilevs_corr   ! =et_nlevels (or number of channels)
   integer :: nlevs_types_kinds         ! number of levels*subtypes*datafields 
   integer :: np_used   ! number of procs containing portions of ran fields
   integer, allocatable :: jrange(:,:)  ! lat index range for each processor 
!
   real(rkind1), pointer :: ran_fields(:,:,:) => null()  ! SHMEM field
!
   contains
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine random_fields_clean (myid)
!
! Deallocate arrays allocated in random_fields_setup
!
   implicit none
!
   integer, intent(in) :: myid
   integer :: ierr
!
   if (allocated(jrange)) then
     deallocate (jrange)
     call MAPL_DeallocNodeArray (ran_fields,rc=ierr)
   endif
!
   end subroutine random_fields_clean 
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine random_fields_spectra (lprint,lev1,lev2,vars,corr_shape, &
                                     h_lengths,fnameSV,scoefs)
!
! Create random spectral coefficicients drawn from a distribution such that
! the expected power spectrum corresponds to that for random fields on the 
! sphere characterized by horizontal correlations functions described as 
! Gaussians in distance. 
!
   use m_parameters, only : earthr
   use m_shtrans_rf, only : mmax,nmax,kmax
   use m_shtrans_rf, only : jindex
   use m_random_power, only : random_power_compute
!
! Use module for printing info when aborting jobs   
   use m_die, only : mpi_die
!
   implicit none
!
   integer, parameter :: r4=4
   logical, intent(in) :: lprint
   integer, intent(in) :: lev1,lev2
   real(r4), intent(in) :: vars(et_nlevels)      ! expected variance of field
   real(r4), intent(in) :: h_lengths(et_nlevels) ! correlation length in km 
   complex(rkinds), intent(out) :: scoefs(nspects,et_nlevels) ! spectral coefs
   character(len=*), intent(in) :: corr_shape
   character(len=*), intent(in) :: fnameSV    ! field is 'S'=scalar, 'V'=vector
!
   integer, parameter :: rr8=8  ! precision required for ran function call
   integer :: n, k, m, j    
   integer :: ierr
   real(rr8), parameter :: rmax=5.d0     ! don't allow greater stdv
   real(rr8), parameter :: xmean=0.d0    ! use unbiased distributions
   real(rr8) :: gamma, cr, ci
   real(rkind1) :: x1
   real(rkind1) :: xc,xsum,earthr2
   real(rkind1) :: power(0:nmax,et_nlevels)
!
   character(len=*), parameter :: my_name='random_fields_spectra'
!
   earthr2=earthr**2
!
! Compute desired expectation power of random spectral coefs given corr. func. 
! Power here is normalized such that global mean variance is 1.
   call random_power_compute (nmax,et_nlevels,lev1,lev2,lprint, &
                              corr_shape,h_lengths,power,ierr)   
   if (ierr /= 0) then
     print *,'ERROR in call to random_power_compute'
     call mpi_die(my_name,44)
   endif
!
! Modify power if vector field since the random spectral coefficients will 
! be in terms of those for vorticity and divergence rather than u and v.
! First, change power in terms of u,v to in terms of vort, divg 
   if (fnameSV == 'V') then ! spectra for a vector field 
     power(0,:)=0.          ! reset global mean of vort and divg to 0  
     do n=1,nmax
       xc=real(n*n+n)/earthr2
       power(n,:)=power(n,:)*xc
     enddo 
!
! Re-normalize power now that power(0,:) has been reset to 0
     do k=lev1,lev2
       xsum=0.
       do n=1,nmax
         xc=real(2*n+1)*earthr2/real(n*n+n)  ! 2n+1 is number of m given n
         xsum=xsum+power(n,k)*xc
       enddo
       power(:,k)=power(:,k)/xsum
     enddo
   endif
!
! Rescale power by desired variance
   do k=lev1,lev2
     power(:,k)=power(:,k)*vars(k)  
   enddo
!
! In the following, the assumption of isotropy implies that the contribution 
! to the power for each wavenumber n is the same for all -n <= m <= n and 
! all possibly non-zero real and imaginary components. 
!
! First synchronize random numbers so that results are independent of number
!  of processors
   if (lev1 > 1) then
     do k=1,lev1-1 
       do n=0,nmax 
         do m=0,n
           call random_gauss_r8 (xmean,gamma,rmax,cr)
           if (.not.((m == 0) .or. (2*m==imax))) then
             call random_gauss_r8 (xmean,gamma,rmax,ci)
           endif  
         enddo  ! loop over m
       enddo    ! loop over n 
     enddo      ! loop over k
   endif        ! check on lev1>1
!
   do k=lev1,lev2
     do n=0,nmax 
       do m=0,n
         if (m == 0 .or. imax == 2*m) then
           x1=1.
         else
           x1=0.5
         endif
         gamma=sqrt(power(n,k)*x1) 
!
         call random_gauss_r8 (xmean,gamma,rmax,cr)
         if ((m == 0) .or. (2*m==imax)) then
           ci=0._rkind2          ! the imag. part for m=0 or m=imax/2 is zero
         else
           call random_gauss_r8 (xmean,gamma,rmax,ci)
         endif  
         j=jindex(m,n-m)
         scoefs(j,k)=cmplx(cr,ci,rkinds)
       enddo  ! loop over m
     enddo    ! loop over n 
   enddo      ! loop over k
!
! Complete Synchronization of random numbers
   if (lev2 < et_nlevels) then 
     do k=lev2+1,et_nlevels
       do n=0,nmax 
         do m=0,n
           call random_gauss_r8 (xmean,gamma,rmax,cr)
           if (.not.((m == 0) .or. (2*m==imax))) then
             call random_gauss_r8 (xmean,gamma,rmax,ci)
           endif  
         enddo  ! loop over m
       enddo    ! loop over n 
     enddo      ! loop over k
   endif        ! check on lev1>1
!
   end subroutine random_fields_spectra
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine random_fields_setup (lprint,dtype,myid,npet,ltest,iret)
!
! Setup and create the random, horizontally correlated fields for later
! extraction of spatially correlated perturbations.
!
   use m_rf_diags_power, only : rf_diags_power_calc
!
! Use module for printing info when aborting jobs
   use m_die, only : mpi_die
!
   implicit none
!
   include "mpif.h"
!
   logical, intent(in) :: ltest
   logical, intent(in) :: lprint
   integer, intent(in) :: myid  ! processor id number  
   integer, intent(in) :: npet  ! number of processors
   integer, intent(out) :: iret ! return code (0 if OK)
   character(len=*), intent(in) :: dtype     ! data type name
!
   integer :: dim3(3)
   integer :: jmd,jmdp,proc
   integer :: kmd
   integer :: id, kid, nf, n, n2, k, i, j, ns, k2
   integer :: ierr
   integer :: sdim3
   integer :: nscoefs
   integer, allocatable :: klev1(:), klev2(:)
   integer :: lev1, lev2
   real(rkind1), allocatable :: vars(:,:)
   real(rkindf), allocatable :: field1(:,:,:)
   real(rkind1), allocatable :: power(:,:,:)
   real(rkind1), allocatable :: evs(:,:)
   complex(rkinds), allocatable :: scoefs(:,:,:)
   complex(rkinds), allocatable :: scoefs1(:,:,:)
   character(len=1) :: fname,fnameSV
   character(len=21) :: info
!
   if (lprint) then 
     print *,' '
     print *,'Create spatially correlated random fields'
   endif
!
! Set variances of all random fields equal to 1. When referencing them later, 
! the random field values will be recaled by appropriate standard deviations.
   allocate (vars(et_nlevels,et_itypes_corr))
   vars(:,:)=1. 
!
! Set spectral truncation parameters and horizontal grid dimensions
   call sh_init_factors (et_rf_nmax)
!
! Allocate shared memory array to hold random fields 
! If PREPBUFR, then allow enough space for vector fields.
   if (trim(dtype) == 'PREPBUFR') then
     et_ran_fields_kinds=2
   else
     et_ran_fields_kinds=1
   endif
   ifields_corr=et_ran_fields_kinds
   ilevs_corr=et_nlevels
   nlevs_types_kinds=et_nlevels*et_itypes_corr*et_ran_fields_kinds
   dim3=(/imax,jmax,nlevs_types_kinds/)
   call MAPL_AllocNodeArray (ran_fields,dim3,rc=iret)
   if (iret /= 0) call mpi_die ('m_random_fields:AllocNodeArray',iret)
!
! Determine max number of lats in any single processor and the number 
! of processors actually used.
   jmd=1+(jmax-1)/npet
   if (jmd > 1) then
     np_used=npet
   else
     np_used=jmax
   endif
   allocate (jrange(2,0:np_used-1))
!
! Determine first and last latitudes assigned to each processor.
! Include overlap 
   jmdp=jmax-(jmd-1)*npet
   jrange(1,0)=1
   jrange(2,0)=jmd+1
   do proc=1,np_used-1    
     jrange(1,proc)=jrange(2,proc-1)-1
     if (proc < jmdp) then
       jrange(2,proc)=min(jrange(1,proc)+jmd+2,jmax)
     else
       jrange(2,proc)=min(jrange(1,proc)+jmd+1,jmax)
     endif
   enddo 
!
   do proc=0,np_used-1
     if (myid == proc) then
       allocate (field1(imax,jrange(1,proc):jrange(2,proc),et_nlevels))
     endif
   enddo
!
! Determine distribution of levels over processors for calc of spectral coefs
   allocate (klev1(0:npet-1))
   allocate (klev2(0:npet-1))
   kmd=(et_nlevels-1)/np_used
   klev1(0)=1
   klev2(0)=klev1(0)+kmd
   do proc=1,npet-1
     klev1(proc)=klev2(proc-1)+1
     klev2(proc)=min(klev1(proc)+kmd,et_nlevels)
   enddo
!
   if (et_l_vc_corr) then 
     allocate (evs(et_nlevels,et_nlevels))
   endif
!
! Loop over subsets of this data type. All subtypes in each subset will 
! be correlated similarly (i.e., using the same sets of correlation lengths)
   do id=1,et_itypes_corr           
!
! Determine if a random scalar or vector field is required for this subtype
     if (et_ran_fields_kinds == 2 .and. trim(dtype) == 'PREPBUFR' .and. &
         et_ks_list(1,id) > 199 .and. et_ks_list(1,id) < 300) then 
       fnameSV='V'  ! denotes vector fields
       sdim3=2
     else   
       fnameSV='S'  ! denotes scalar fields   
       sdim3=1  
     endif
!
! If required, compute product of eigenvectors of cov matrix with sqrt 
! of corresponding eigenvalues
     if (et_l_vc_corr) then 
       do n=1,et_nlevels
         do k=1,et_nlevels
           evs(n,k)=et_e_v_sqrt(k,id)*et_e_vects(n,k,id)
         enddo
       enddo
     endif 
!
     lev1=klev1(myid)
     lev2=klev2(myid)
     nscoefs=et_nlevels*sdim3
     allocate (scoefs(nspects,et_nlevels,sdim3))
     allocate (scoefs1(et_nlevels,sdim3,2))    
     do nf=1,et_ran_fields_kinds  ! loop over number of different obs data fields 
       kid=nf+(id-1)*et_ran_fields_kinds
!
! Make sure that random number seed the same for all processors.
       call random_fields_synch_ranf (myid,npet)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)      
!     
! Create spectral coefs of random fields, distributed over processors 
! according to the range of channels or levels lev1:lev2.
! For vector field, this is done only the first time since the same 2 sets of
! coefs (for vort and divg) are needed for both u and v fields.  
! Note that the same hcorr lengths are used for u and v 
       if (sdim3 == 2 .and. nf == 1 .and. fnameSV == 'V') then
         fname='U'       ! spectra for psi and chi both required
         call random_fields_spectra (lprint,lev1,lev2,vars(:,id),               &
                                     et_corr_shapes(id),et_hcorr_lengths(:,id), & 
                                     fnameSV,scoefs(:,:,1))
         call random_fields_spectra (lprint,lev1,lev2,vars(:,id),               &
                                     et_corr_shapes(id),et_hcorr_lengths(:,id), & 
                                     fnameSV,scoefs(:,:,2))
       elseif  (sdim3 == 2 .and. nf == 2 .and. fnameSV == 'V') then
         fname='V'      ! use same spectra as for U field
       else             ! default is for scalar fields
         fname='S'
         call random_fields_spectra (lprint,lev1,lev2,vars(:,id),               &
                                     et_corr_shapes(id),et_hcorr_lengths(:,id), & 
                                     fnameSV,scoefs(:,:,1))
       endif
!
! Test on whether random V field is to be produced. If so, spectral coefs 
! have already been copied when nf=1. If required, they also have already 
! been transformed vertically.
!
       if (trim(fname) /= 'V') then
!
         do n=1,nspects   
           scoefs1(:,:,1)=cmplx(0._rkinds,0._rkinds,rkinds)
           do proc=0,np_used-1
             if (myid == proc) then
!
! If required, apply partial transform from coefs of vertical EOFs to coefs for
! horizontal fields. Otherwise, simply copy coefs to arrray used to pass 
! coefs just produced to all other processors. The transform is only partial
! because it is applied only to spectral coefs for those EOFS considered on
! the current processor.  
               do i=1,sdim3
                 do k2=lev1,lev2
                   if (et_l_vc_corr) then 
                     do k=1,et_nlevels
                       scoefs1(k,i,1)=scoefs1(k,i,1)+evs(k,k2)*scoefs(n,k2,i)
                     enddo
                   else
                      scoefs1(k2,i,1)=scoefs(n,k2,i)
                   endif
! 
                 enddo  ! loop over levels available on processor
               enddo    ! loop over fields
             endif      ! test on proc id
           enddo        ! loop over processors used for this calculation
! 
! Replace scoefs for current value of n by sum of scoefs1 over all processors.
! This makes scoefs on all processors identical.
           call MPI_Barrier(MPI_COMM_WORLD,ierr)             
           call MPI_ALLreduce (scoefs1(:,:,1),scoefs1(:,:,2),nscoefs, &
                                   MPI_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
           call MPI_Barrier(MPI_COMM_WORLD,ierr)      
           scoefs(n,:,:)=scoefs1(:,:,2)
         enddo          ! loop over nspects
!
       endif  ! test on whether V field is being produced
!
! Print diagnsotics if requested for testing purposes
       if (lprint .and. ltest .and. &
           (sdim3 == 1 .or. (sdim3 == 2 .and. nf == 1))) then 
         write(info,'(a,i2,i2,2a)'),'Spects_in_setup',id,nf,' ',fnameSV
         call rf_diags_power_calc (et_nlevels,sdim3,kid,scoefs,fnameSV,info)
       endif
!
! Project spectral coefs on to lat range of fields
       do proc=0,np_used-1
         if (myid == proc) then
           call sh_trans (et_nlevels,jrange(1,proc),jrange(2,proc),scoefs, &
                          sdim3,field1(:,:,:),fname)
!
! Copy from portion of fields on specified processor to shared memory location
           do n=1,et_nlevels
             k=n+(kid-1)*et_nlevels 
             do j=jrange(1,proc),jrange(2,proc)
               ran_fields(:,j,k)=field1(:,j,n)
             enddo
           enddo
         endif 
       enddo                   
       call MPI_Barrier(MPI_COMM_WORLD,ierr)      
!
     enddo  ! loop over nf 
     deallocate (scoefs)
     deallocate (scoefs1)
   enddo    ! loop over id
!
! Make sure that random number seed the same for all processors.
   call random_fields_synch_ranf (myid,npet)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)      
!
   if (allocated(evs)) deallocate (evs)
   deallocate (klev1,klev2)
   deallocate (vars)
   deallocate (field1)
!
   call sh_clean 
!
   if (lprint) then
     print ('(i4,2a,i2,a)'),et_nlevels,' spatially correlated random ',&
                      'fields created for ',et_itypes_corr,' types'         
     print ('(a,2i5)'),' random field grid imax,jmax =',imax,jmax
   endif
!
  end subroutine random_fields_setup 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine random_fields_synch_ranf (myid,npet)
!
!  Synchronize random numbers so all processors create the same sequence
!  at this point
!
   use m_shtrans_rf, only : nmax
!
   implicit none
   include "mpif.h"
! 
   integer, intent(in) :: myid
   integer, intent(in) :: npet
!
   integer :: ierr
   integer :: nproc
   integer :: itag
   integer :: iseed
   integer :: k,n,m
   integer :: stat(MPI_STATUS_SIZE)
   integer :: i_random_seed(2)
   real(8) :: x
   real(8) :: dum
!
   if (myid == 0) then
     call random_number (x)
     iseed=int(1.e6*x)
     do nproc=1,npet-1
       itag=nproc
       call MPI_SEND (iseed,1,MPI_INTEGER,nproc,itag,MPI_COMM_WORLD,ierr)
     enddo
   endif
! 
   do nproc=1,npet-1
     if (myid == nproc) then
       itag=nproc
       call MPI_RECV (iseed,1,MPI_INTEGER,0,itag,MPI_COMM_WORLD,stat,ierr)
     endif
   enddo
!
   i_random_seed(1)=iseed
   i_random_seed(2)=1111    ! an arbitrary integer here
   call random_seed (put=i_random_seed(1:2))
!
   end subroutine random_fields_synch_ranf 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine random_fields_get_1_lev (kpid,f1lev)
!
! Get full horiz field for a single level
!                                    
   implicit none
!
   integer, intent(in) :: kpid
   real(rkind1), intent(out) :: f1lev(imax,jmax)
!
   f1lev(:,:)=ran_fields(:,:,kpid)
!
   end subroutine random_fields_get_1_lev 
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine random_fields_get_1_value (xlon,xlat,kpid,f1)
!
! Get single interpolated value at requested lat, lon, and third index
!                                    
   implicit none
!
   integer, intent(in) :: kpid  ! index for requeted level, type and field
   real(rkind1), intent(in) :: xlat, xlon   ! location 
   real(rkind1), intent(out) :: f1          ! point value
!
   integer :: lonE, lonW, latN, latS
   real(rkind1) :: rlon, dlon, rlat, dlat, lonx
   real(rkind1) :: weightE, weightW, weightN, weightS
   real(rkind1) :: wNE, wNW, wSE, wSW
!
! determine indexes and weights for E-W interpolation
   lonx=mod(xlon,360.)
   if (lonx < 0.) then     ! ensures values are in range (0.,360.)
     lonx=lonx+360.
   endif
   if (lonx > 360.) then   ! accounts for round off
     lonx=0.
   endif
   dlon=360./real(imax)
   rlon=lonx/dlon        ! longitude in units of dlon
   lonW=1+int(rlon)      ! grid lon to the west (1=first grid point; r=0)
   lonW=min(lonW,imax)
   if (lonW==imax) then   
     lonE=1
   else 
     lonE=lonW+1         ! grid lon to the east
   endif
   weightE=rlon-lonW+1.  
   weightW=1.-weightE
!
! determine indexes and weights for N-S interpolation
   dlat=180./(jmax-1)
   rlat=(xlat+90.)/dlat
   latS=1+int(rlat)
   latS=max(latS,1)
   latS=min(latS,jmax-1)
   latN=latS+1
   weightN=rlat-latS+1.
   weightS=1.-weightN     
!
! set weights for interpolations
   wNE=weightN*weightE
   wNW=weightN*weightW
   wSE=weightS*weightE
   wSW=weightS*weightW
!   
   f1=wSE*ran_fields(lonE,latS,kpid)+wSW*ran_fields(lonW,latS,kpid)+ &
      wNE*ran_fields(lonE,latN,kpid)+wNW*ran_fields(lonW,latN,kpid) 
!
   end subroutine random_fields_get_1_value
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine random_fields_get_values (ndim1,nobs_levs,ndum,nf,nobs,plevels, &
                                      xlat,xlon,stdv,fpoints,weights_cu,itype)
!
! Perform horizontal interpolation of random fields
! for single-level conventional data, also perform vertical interpolation
!                                    
   implicit none
!
   integer, intent(in) :: ndim1  ! dimension of plevels variable
   integer, intent(in) :: nobs_levs  ! number of obs levels or channels
   integer, intent(in) :: ndum   ! not used
   integer, intent(in) :: nf     ! index for data field (e.g., (1) u; (2) v)
   integer, intent(in) :: nobs   ! obs counter number   
   integer, intent(in) :: itype  ! index of group of data subtypes considered
   real(rkind2), intent(in) :: plevels(ndim1)          ! pressure levels of obs 
   real(rkind2), intent(in) :: xlat, xlon              ! location of obs
   real(rkind2), intent(in) :: stdv(nobs_levs)         ! standard deviations of errors 
   real(rkind2), intent(out) :: fpoints(nobs_levs)      ! point values of random corr fields
   real(rkind2), intent(out) :: weights_cu(nobs_levs,2) ! weights for correl and uncorr parts
!
   integer :: id3_offset   
   integer :: lonE, lonW, latN, latS
   integer :: k, k1, k2, kp, kpid, id
   integer :: iopt
   real(rkind1), parameter :: zero=0._rkind1
   real(rkind1), parameter :: one=1._rkind1
   real(rkind2) :: rlon, dlon, rlat, dlat, lonx, latx
   real(rkind1) :: weightE, weightW, weightN, weightS
   real(rkind1) :: wNE, wNW, wSE, wSW
   real(rkind1) :: wv, wv1, wv2
   real(rkind1) :: frac_corr_int
   real(rkind1) :: xtest1, xtest2
!
! Determine offset for 3rd dimension in ran_fields array 
   id3_offset=(nf+(itype-1)*et_ran_fields_kinds-1)*et_nlevels
   id=itype
!
! Augment lon and lat so that when more than one obs subtype uses the 
! same random field, each type will actually use a different region. 
! (This current implementation does not do so near the poles however.)
   if (1 == mod(itype,2)) then
     latx=xlat
   else
     latx=-xlat
   endif 
!
! Re-define lon so that it is in range 0.<= lon < 360.
   lonx=xlon+360.*(et_hcorr_ks(1)-1)/et_hcorr_ks(2)
   lonx=mod(lonx,360._rkind2)
   if (lonx < zero) lonx=lonx+360. 
   if (lonx >= 360. ) lonx=zero       ! accounts for round off
!
! Determine indexes and weights for E-W interpolation
   dlon=360./real(imax)
   rlon=lonx/dlon        ! longitude in units of dlon
   lonW=1+int(rlon)      ! grid lon to the west (1=first grid point; r=0)
   lonW=min(lonW,imax)
   if (lonW==imax) then   
     lonE=1
   else 
     lonE=lonW+1         ! grid lon to the east
   endif
   weightE=rlon-lonW+1.  
   weightW=1.-weightE
!
! Determine indexes and weights for N-S interpolation
   dlat=180./(jmax-1)
   rlat=(latx+90.)/dlat
   latS=1+int(rlat)
   latS=max(latS,1)
   latS=min(latS,jmax-1)
   latN=latS+1
   weightN=rlat-latS+1.
   weightS=1.-weightN     
!
! Set weights for interpolations
   wNE=weightN*weightE
   wNW=weightN*weightW
   wSE=weightS*weightE
   wSW=weightS*weightW
!
   if (.not. et_l_vc_corr) then  
     iopt=1
!
! Option 1:
! Perform only horizontal interpolation at each separate level or channel
! since there is no vertical or channel correlation. The assumption here is that
! the horizontally correlated random fields each have variance 1. 
! The standard deviations of the random field values output here are therefore 
! given by the values of stdv input to this routine.
!
       if (nobs_levs > et_nlevels) then
         print *,' '
         print *,'ERROR IN random_fields_get_values Option 1:'
         print *,' nobs_levs=',nobs_levs,' must not exceed et_nlevels=',et_nlevels
         return
       endif
       do kp=1,nobs_levs
         call random_fields_frac_weights (et_frac_corr(kp,id), &
                                       weights_cu(kp,1),weights_cu(kp,2))
         kpid=kp+id3_offset
         fpoints(kp)=stdv(kp) *                                         &
               ( wSE*ran_fields(lonE,latS,kpid)+ &
                 wSW*ran_fields(lonW,latS,kpid)+ &
                 wNE*ran_fields(lonE,latN,kpid)+ &
                 wNW*ran_fields(lonW,latN,kpid) ) 
       enddo 
!
     elseif (et_pmax <= et_pmin) then
       iopt=2
!
! Option 2:
! Determine values by summing contributions by all EOFs, but without 
! considering interpolation between levels or channels. Here, the assumption 
! is that the horizontally correlated random fields already include the 
! transformation from horizontal fields of EOF coefficients to horizontal
! fields of channel error. When there are many observation locations, 
! each having many levels, it is cheaper to apply this transformation once to
! all the random spectral coefficients than to each observation location.  
! The random fields here include the proper weighting functions of the 
! EOF error variances. 
!
       do kp=1,nobs_levs
         call random_fields_frac_weights (et_frac_corr(kp,id), &
                                       weights_cu(kp,1),weights_cu(kp,2))
         kpid=kp+id3_offset
         fpoints(kp)=wSE*ran_fields(lonE,latS,kpid)+ &
                     wSW*ran_fields(lonW,latS,kpid)+ &
                     wNE*ran_fields(lonE,latN,kpid)+ &
                     wNW*ran_fields(lonW,latN,kpid)  
       enddo
!
     else
       iopt=3
!
! Option 3:
! Determine values by interpolating between summed contributions by all 
! EOFs at levels above and below those sandwich the levels of the obs.
! The pressure levels at which the EOFs are defined aew assumed to be ordered
! from lowest to highest in elevation (i.e., with decreasing p).
! k1 is index for level above obs level (in elevation).    
! k2 is index for level above obs level (in elevation).    
! The p levels for the random fields are assumed here to be equally spaced 
! in distance between the p-levels et_pmaxX and et_pminX, assming an isothermal 
! column so that delta log p is proportional to distance. wv is then the 
! distance below the level et_pminX in units of the vertical spacing of the 
! fields.  wv1 and wv2 are the weights used for averaging between the 
! pressure levels using a linear function of distance. Here, the assumption 
! is that the horizontally correlated random fields already include the 
! transformation from horizontal fields of EOF coefficients to horizontal
! fields of channel error. When there are many observation locations, 
! each having many levels, it is cheaper to apply this transformation once to
! all the random spectral coefficients than to each observation location.  
! The random fields here include the proper weighting functions of the 
! EOF error variances. 
!
       do kp=1,nobs_levs
!
! First determine levels to interpolate vertically
         if (plevels(kp) <= et_pmin) then      ! obs p above p-range top
           k1=et_nlevels
           k2=k1
           wv=-1   ! used as flag here
           wv1=one
           wv2=zero
         elseif (plevels(kp) >= et_pmax) then  ! obs p below p-range bottom
           k1=1
           k2=k1
           wv=-2    ! used as flag here
           wv1=one
           wv2=zero
!
         else                               ! obs p within p-range
           wv=(et_nlevels-1)*log(et_pmax/plevels(kp))/log(et_pmax/et_pmin)
           k1=max(0,int(wv))
           k1=min(k1,et_nlevels-2)
           wv2=max(zero,wv-real(k1))
           wv2=min(wv2,one)
           wv1=max(zero,one-wv2)
           k1=k1+1  ! augment index since level indexes start at 1, not 0
           k2=k1+1  ! indicates 1st level above obs level (in terms of height)
         endif
!
         frac_corr_int=wv1*et_frac_corr(k1,id)+wv2*et_frac_corr(k2,id)
         call random_fields_frac_weights (frac_corr_int, &
                                       weights_cu(kp,1),weights_cu(kp,2))
         kpid=k1+id3_offset
         k=k2+id3_offset
         fpoints(kp)=wv1*                         &
               ( wSE*ran_fields(lonE,latS,kpid)+  &
                 wSW*ran_fields(lonW,latS,kpid)+  &
                 wNE*ran_fields(lonE,latN,kpid)+  &
                 wNW*ran_fields(lonW,latN,kpid) ) & 
                    +wv2*                      &
               ( wSE*ran_fields(lonE,latS,k)+  &
                 wSW*ran_fields(lonW,latS,k)+  &
                 wNE*ran_fields(lonE,latN,k)+  &
                 wNW*ran_fields(lonW,latN,k) )  
       enddo  ! loop over kp
     endif    ! test on which type of vertical or channel correlation
!        
   end subroutine random_fields_get_values
!
!
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   subroutine random_fields_frac_weights (frac,w_corr,w_uncorr)
!
! Determine weights for correlated and uncorrelated portions of 
! the perturbations based on the fraction of total variance that 
! is correlated
!
   real(rkind1), intent(in)  :: frac
   real(rkind2), intent(out) :: w_corr
   real(rkind2), intent(out) :: w_uncorr
!
   real(rkind1) :: var_uncorr   
!
   if (frac > 0.) then
     var_uncorr=1.-frac
     w_corr=sqrt(frac)
   else
     var_uncorr=1.
     w_corr=0.
   endif
!
   if (var_uncorr > 0.) then
     w_uncorr=sqrt(var_uncorr)
   else
     w_uncorr=0.
   endif
!
   end subroutine random_fields_frac_weights
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
   end module m_random_fields
