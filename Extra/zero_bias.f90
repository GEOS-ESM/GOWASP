      Program zero_bias
!           
      integer, parameter :: nlines=2877
      integer, parameter :: nchars=116
      character(len=nchars) :: char_line 
! 
      integer, parameter :: ngroups=nlines
      integer, parameter :: nbang1=47
      integer, parameter :: nbang2=74
      character(len=nbang1) :: bang1 
      character(len=nbang2) :: bang2 
!
      open(unit=10,file='bias_in.txt',form='formatted')
      open(unit=11,file='bias_out.txt',form='formatted')
      do n=1,nlines
        read (10,'(a)') char_line
                          
        char_line(33:44)='    0.000000'
        char_line(45:56)='    0.000000'
        char_line(57:68)='    0.000000'
        char_line(69:80)='    0.000000'
        char_line(81:92)='    0.000000'
        char_line(93:104)='    0.000000'
        char_line(105:116)='    0.000000'
        write (11,'(a)') char_line  
      enddo
      close (10)
      close (11)
!
      open(unit=10,file='bang_in.txt',form='formatted')
      open(unit=11,file='bang_out.txt',form='formatted')
      do n=1,ngroups
        read (10,'(a)') bang1 
! keep this value as read        bang1(35:47)=' 0.000000E-00'
        write (11,'(a)') bang1 
        do m=1,10
          read (10,'(a1)') bang2(1:1)
        enddo
        bang2(1:74)= &
           '      0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000'
        do m=1,9 
          write (11,'(a)') bang2
        enddo 
        bang2(1:74)= &
     '                                                                          '
        write (11,'(a)') bang2                        
      enddo
      close (10)
      close (11)
!
      end Program zero_bias
