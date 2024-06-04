!-------------------------------------------------------------------------------
subroutine wait (ierr)

use iomodule,      only: stdin, stdout, ilu, lucmd
use errors,        only: ok, eof, stop

implicit none

!-pass
integer :: ierr
!-locl
character(len=1)c
!-init
ierr = OK

10    write(stdout,1,advance='no')

 1    format(1x,'press "return" to continue, "q" to quit > ')

read(stdin,'(A)',err=20,end=30) c

if (c.eq.'q') then
    ierr = STOP
    if (lucmd(ilu).ne.stdin) then
        ierr=EOF
    endif
endif

20    return
30    ierr = STOP
goto 20
end
