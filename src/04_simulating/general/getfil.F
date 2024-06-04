character(len=255) function getfil(string)

use debugmodule,     only: xit
use iomodule,        only: stderr

implicit none

!-pass
character(len=255) :: string

!-locl
character(len=255) :: HOME
integer            :: fnblnk
integer, external  :: lnblk

integer            :: nf, nl

call gethome(HOME)

nf = fnblnk(HOME)
nl = lnblk(HOME)

if (HOME(nl:nl).eq.'/') then
    nl=nl-1
endif

if (nf.gt.nl .or. nl.le.0) then
    write(stderr,1)
1   format(1x,'getfil: get_environment_variable(HOME) failed')
    call xit(1," ")
endif

if (string.eq.'temp0.F') then
    getfil = HOME(nf:nl)//'/src/TempSrc/temp0.F'
elseif (string.eq.'addan.F') then
    getfil = HOME(nf:nl)//'/src/TempSrc/addan.F'
else
    if (lnblk(string).gt.0) then
        write(stderr,2) string(1:lnblk(string))
 2            format(1x,'getfil: argument not recognized :',80a)
        call xit(1," ")
    else
        write(stderr,3)
 3            format(1x,'getfil: no argument')
        call xit(1," ")
    endif
endif

return
end
