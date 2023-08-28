   module m_gpsro_names
!
!  Specify set of variable names and arrays used for GPSRO observational data
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   implicit none
!
   public
!
   integer, parameter :: obs_nhead=16        ! header size 
   integer, parameter :: obs_nvalues=5
   integer, parameter :: obs_max_levs=400
   integer, parameter :: obs_info_extra=6
   integer, parameter :: obs_info_num=obs_info_extra+obs_nhead
!
! declare indexes for various obs values read  
   integer :: bbsiid, bbbaz, bbyear, bblat, bblon, bbsaid
   integer :: bbnum, bbnlold, bbnlnew, bbnmsg, bbdhr, bbtslot
   integer :: bbgeodu, bbcurve 
   integer :: bbimpp, bbbang          ! impact parameter, bending angle
   integer :: bblatk, bblonk, bbbazk  ! lat, lon, bearing that vary with impact
   integer :: bbptid, bbpccf          ! used in the code merge_bufr
!
   character(len=*), parameter :: obs_hdstr1=        &
        'SIID CLATH CLONH YEAR MNTH DAYS HOUR MINU '
   character(len=*), parameter :: obs_hdstr2=        &
        'SECO SAID PTID ELRC BEARAZ GEODU PCCF QFRO '
   character(len=*), parameter :: info_extra_names=  &
        'obsnum nlevold nlevnew nmsg dhr ntslot'
   character(len=*), parameter :: obs_hdstr=obs_hdstr1//obs_hdstr2
   character(len=8) :: obs_value_names(obs_nvalues)
   character(len=8) :: obs_info_names(obs_info_num)
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine gpsro_names_setup (ier)
!
!  Specify indexes for GPSRO observational data array values corresponding 
!  to particular information.
!
   implicit none
!
   integer, intent(out) :: ier
!
   logical, parameter :: lstop=.false.
   integer :: nn
   integer :: bbb1, bbb2
   character(len=*), parameter :: cb='m_gpsro_names' 
   character(len=240) :: charstr
!
   ier=0
!
   obs_value_names(1)='imppara '
   bbimpp=1
   obs_value_names(2)='bendang '
   bbbang=2
   obs_value_names(3)='latimp  '
   bblatk=3
   obs_value_names(4)='lonimp  '
   bblonk=4
   obs_value_names(5)='bazimp  '
   bbbazk=5
!
   nn=obs_info_num
   charstr=obs_hdstr//info_extra_names
   read (charstr,*) obs_info_names(1:nn)
!
   call find_name (nn,obs_info_names,lstop,cb,'SIID',     bbsiid)
   call find_name (nn,obs_info_names,lstop,cb,'BEARAZ',   bbbaz)
   call find_name (nn,obs_info_names,lstop,cb,'YEAR',     bbyear)
   call find_name (nn,obs_info_names,lstop,cb,'CLATH',    bblat)
   call find_name (nn,obs_info_names,lstop,cb,'CLONH',    bblon)
   call find_name (nn,obs_info_names,lstop,cb,'SAID',     bbsaid)
   call find_name (nn,obs_info_names,lstop,cb,'obsnum',   bbnum)
   call find_name (nn,obs_info_names,lstop,cb,'nlevold',  bbnlold)
   call find_name (nn,obs_info_names,lstop,cb,'nlevnew',  bbnlnew)
   call find_name (nn,obs_info_names,lstop,cb,'nmsg',     bbnmsg)
   call find_name (nn,obs_info_names,lstop,cb,'dhr',      bbdhr)
   call find_name (nn,obs_info_names,lstop,cb,'ntslot',   bbtslot)
   call find_name (nn,obs_info_names,lstop,cb,'GEODU',    bbgeodu)
   call find_name (nn,obs_info_names,lstop,cb,'ELRC',     bbcurve)
   call find_name (nn,obs_info_names,lstop,cb,'PTID',     bbptid)
   call find_name (nn,obs_info_names,lstop,cb,'PCCF',     bbpccf)
!
! Simply check that no returned indexes are 0 (i.e., at least 1 not found)
   bbb1=bbsiid*bbbaz*bbyear*bblat*bblon*bbsaid*bbgeodu*bbpccf
   bbb2=bbnum*bbnlold*bbnlnew*bbnmsg*bbdhr*bbtslot*bbcurve*bbptid
  if (bbb1 == 0 .or. bbb2 == 0) ier=1
!
   end subroutine gpsro_names_setup 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_gpsro_names
