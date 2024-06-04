 subroutine OfClose (file,monitor)

use iomodule, only: lumon, stdout
!
! Opens monitor file if neccessary

implicit none
!
!-pass
 logical monitor
 character(len=*) file

!-locl
 integer lnblk
!
 if (.not. monitor) return
 if (lnblk(file).le.0) return
!
 call closef(luMon)
 luMon = stdout
!
 return
 end
