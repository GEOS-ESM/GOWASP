  program xx
  integer :: n
  character(len=10) :: c
  n=99
  c='ABCD='
  open (10,file='testxxx')
  write (10,'(a,i4)') trim(c),n
  close (10)
  open (10,file='testxxx')
  read (10,*) c,n
  print *,c,n
  close (10)
  n=44
  open (10,file='testxxx2')
  write (10,'(a,i4)') trim(c),n
  close (10)

  end program xx
