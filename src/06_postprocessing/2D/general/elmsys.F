subroutine elmsys (ierr)

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: plotdatactx
use debugmodule,     only: iecho
use iomodule,        only: stderr
use modelTopology,   only: NEN,NSUR
use errors,          only: ok

implicit none
!-pass
integer ierr



!-locl
#include "lenchar.i"
logical makeij
character(len=LENCHAR)IJFILE
!-init
ierr = OK
IJFILE = 'IJOUT.DAT'

!    test modification times of FE-input file relative to IJFILE
if ( makeij(IJFILE) ) then
!   set up array of element neighbours
    call shared (meshdatactx%ien,plotdatactx%IJSUR,meshdatactx%neglobal,NEN,NSUR)
    call writij (plotdatactx%IJSUR,meshdatactx%neglobal,NSUR,IJFILE,ierr)
else
    if (iecho.eq.1) then
        write(stderr,1)
    endif
 1        format(1x,'reading "IJOUT.DAT"')
    call readij (plotdatactx%IJSUR,meshdatactx%neglobal,NSUR,IJFILE,ierr)
endif

end subroutine
