   module m_conv_names
!
!  Module containing variable names and indexes for obs information, including
!  bufr header variable names
!
!  Code History:
!  Ronald Errico      Oct 1 2016  Initial GOWASP-3 code
!
   implicit none
!
   public
!
   integer, parameter :: nheadc=7  ! header size common for wind and mass  
   integer, parameter :: nheadm=7  ! header size for mass
   integer, parameter :: nheadw=8  ! header size for wind
   integer, parameter :: conv_nhead=1+nheadw+nheadm-nheadc 
   integer, parameter :: conv_nfields=16
   integer, parameter :: conv_max_levs=255
   integer, parameter :: obs_info_extra=10
   integer, parameter :: obs_info_num=obs_info_extra+conv_nhead  
   integer :: bbx, bby, bbr, bbp, bbz, bbu, bbv, bbt, bbq
   integer :: bbc, bbf, bbpq, bbzq, bbwq, bbtq, bbqq
   integer :: bbnknd, bbnall, bbnlold, bbnlnew, bblints
   integer :: bbpmin, bbtslot, bbtmin, bbnmsg, bboslv
   integer :: ityp,isid,ielv,ilzf,idhr,ixob,iyob ! indexes in common part header
   integer :: iwsd ! index in wind header only  
!
! header string hdstrc is the string in common (must include SID, TYP, and ELV)
! header string hdstrm is the string for ps,t,q (MASS) reports 
! header string hdstrw is the string for u,v (WIND) reports 
   character(len=*), parameter :: hdstrc='SID XOB YOB DHR TYP ELV T29 '
   character(len=*), parameter :: hdstrm_extra=' '
   character(len=*), parameter :: hdstrw_extra='SAID'
   character(len=*), parameter :: hdstrm=hdstrc//hdstrm_extra
   character(len=*), parameter :: hdstrw=hdstrc//hdstrw_extra
   character(len=*), parameter :: hdstr_merge=hdstrc//hdstrw_extra//hdstrm_extra
   character(len=*), parameter :: fldstr1='lat lon rtime p z u v T q '
   character(len=*), parameter :: fldstr2='cat tvflag pqm zqm wqm tqm qqm' 
   character(len=*), parameter :: fldstr=fldstr1//fldstr2
   character(len=*), parameter ::info_extra= & 
     'nobkind noball nlevsold nlevsnew nlevsts pmin tmin tslot nmsg oldslev' 
   character(len=8) :: conv_info_names(conv_nhead)
   character(len=8) :: conv_value_names(conv_nfields)
   character(len=8) :: obs_info_names(obs_info_num)
!
   contains
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   subroutine conv_names_setup (ier)
!
!  Specify indexes for conventional observational data array values 
!  corresponding to particular information.
!
   implicit none
!
   integer, intent(out) :: ier
!
   logical, parameter :: lstop=.false.
   integer :: nn
   integer :: bbb1, bbb2
   character(len=*), parameter :: cb='m_conv_names' 
   character(len=240) :: charstr
!
   ier=0
!
   nn=conv_nfields
   charstr=fldstr
   read (charstr,*) conv_value_names(1:nn)
   call find_name (nn,conv_value_names,lstop,cb,'lon',bbx)
   call find_name (nn,conv_value_names,lstop,cb,'lat',bby)
   call find_name (nn,conv_value_names,lstop,cb,'rtime',bbr)
   call find_name (nn,conv_value_names,lstop,cb,'p',bbp)
   call find_name (nn,conv_value_names,lstop,cb,'z',bbz)
   call find_name (nn,conv_value_names,lstop,cb,'u',bbu)
   call find_name (nn,conv_value_names,lstop,cb,'v',bbv)
   call find_name (nn,conv_value_names,lstop,cb,'T',bbt)
   call find_name (nn,conv_value_names,lstop,cb,'q',bbq)
   call find_name (nn,conv_value_names,lstop,cb,'cat',bbc)
   call find_name (nn,conv_value_names,lstop,cb,'tvflag',bbf)
   call find_name (nn,conv_value_names,lstop,cb,'pqm',bbpq)
   call find_name (nn,conv_value_names,lstop,cb,'zqm',bbzq)
   call find_name (nn,conv_value_names,lstop,cb,'wqm',bbwq)
   call find_name (nn,conv_value_names,lstop,cb,'tqm',bbtq)
   call find_name (nn,conv_value_names,lstop,cb,'qqm',bbqq)
!
   bbb1=bbx*bby*bbr*bbp*bbz*bbu*bbv*bbt*bbq
   bbb2=bbc*bbf*bbpq*bbzq*bbwq*bbtq*bbqq
   if (bbb1 == 0 .or. bbb2 == 0) ier=1
!
   nn=conv_nhead
   charstr=hdstr_merge
   read (charstr,*) conv_info_names(1:nn-1)
   conv_info_names(nn)='lzflags '
   call find_name (nn,conv_info_names,lstop,cb,'TYP',ityp)
   call find_name (nn,conv_info_names,lstop,cb,'SID',isid)
   call find_name (nn,conv_info_names,lstop,cb,'ELV',ielv)
   call find_name (nn,conv_info_names,lstop,cb,'DHR',idhr)
   call find_name (nn,conv_info_names,lstop,cb,'XOB',ixob)
   call find_name (nn,conv_info_names,lstop,cb,'YOB',iyob)
   call find_name (nn,conv_info_names,lstop,cb,'lzflags',ilzf)
   call find_name (nn,conv_info_names,lstop,cb,'SAID',iwsd)
!
   bbb1=ityp*isid*ielv*ilzf*iwsd*idhr*ixob*iyob
   if (bbb1 == 0) ier=ier+10
!
   nn=conv_nhead
   obs_info_names(1:nn)=conv_info_names(1:nn)
   charstr=info_extra
   read (charstr,*) obs_info_names(nn+1:nn+obs_info_extra)
!
   nn=conv_nhead+obs_info_extra 
   call find_name (nn,obs_info_names,lstop,cb,'nobkind',bbnknd)
   call find_name (nn,obs_info_names,lstop,cb,'noball',bbnall)
   call find_name (nn,obs_info_names,lstop,cb,'nlevsold',bbnlold)
   call find_name (nn,obs_info_names,lstop,cb,'nlevsnew',bbnlnew)
   call find_name (nn,obs_info_names,lstop,cb,'nlevsts',bblints)
   call find_name (nn,obs_info_names,lstop,cb,'pmin',bbpmin)
   call find_name (nn,obs_info_names,lstop,cb,'tmin',bbtmin)
   call find_name (nn,obs_info_names,lstop,cb,'tslot',bbtslot)
   call find_name (nn,obs_info_names,lstop,cb,'nmsg',bbnmsg)
   call find_name (nn,obs_info_names,lstop,cb,'oldslev',bboslv)
!
   bbb1=bbnknd*bbnall*bbnlold*bbnlnew*bblints
   bbb2=bbpmin*bbtslot*bbtmin*bbnmsg*bboslv
   if (bbb1 == 0 .or. bbb2 == 0) ier=ier+100
!
   end subroutine conv_names_setup 
!
!
!  x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x 
!
   end module m_conv_names
