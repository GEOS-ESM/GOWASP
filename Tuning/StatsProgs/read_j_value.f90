  Program read_j_value
!
!  Read the counter from a binary accumulated sums file 
!
  integer, parameter :: iunit=10
  integer :: j
  character(len=240) :: filename
!
  call GetArg(1_4,filename)
  open (iunit,file=trim(filename),form='unformatted')
  read (iunit) j
  write (*,*) j
!
  end Program read_j_value
