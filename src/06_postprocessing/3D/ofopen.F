subroutine OfOpen (fname,monitor,ierr)

use iomodule,        only: stderr, append, stdout, lumon
use errors,          only: ok, exec

!
! Opens monitor file if neccessary
implicit none

!-pass
 integer ierr
 logical monitor
 character(len=*) fname


!-locl
integer, external :: lnblk,ifopen, nextlu, lueof
integer :: l
 ierr = OK
!
 if (.not. monitor) return
!
 l = LEN(fname)
 if (l.le.0) return
 l = lnblk(fname)
 if(l.le.0) return
!
 call shvar(fname)
!    redirect monitor data to (ascii) file
 luMon = nextlu(0)
 ierr = ifopen(luMon,fname,'unknown')
 if (ierr.ne.0) goto 100
 if (APPEND) then
ierr = lueof(luMon)
if (ierr.ne.0) goto 100
 endif
 return
!
100     ierr = EXEC
 call closef(luMon)
 luMon = stdout
 return
!
 end
