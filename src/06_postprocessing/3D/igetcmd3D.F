integer function igetcmd (prgname,cmdline,tothetop)

use iomodule,    only: stdout, stdin, ilu, lucmd
use errors,      only: eof, ok, syntax

! gets a command line from unit lucmd
! INPUT returns  1: character(len=)string
!         2: real
!        -1: EOF
!        -2: input error

implicit none
!-pass
character(len=*) prgname,cmdline
logical tothetop
integer, external :: lnblk
!-local
integer :: itype
double precision :: rcmd
integer, external :: input

!-init
rcmd = 0d0

if (lucmd(ilu).eq.stdin) then
    write(stdout,1,advance='no') prgname(1:lnblk(prgname))
 1        format(1x,A)
    write(stdout,2,advance='no')
 2        format(1x,'> ')
endif

if (lucmd(ilu).ne.stdin .and. tothetop) then
    rewind(lucmd(ilu))
endif

itype = input (cmdline,rcmd,lucmd(ilu))

if (itype.eq.-2) then
    igetcmd = SYNTAX
elseif (itype.eq.-1) then
    igetcmd = EOF
else
    igetcmd = OK
endif

return
end
