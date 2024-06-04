subroutine bufopn ()

! opens command buffer file

use iomodule, only: lubuf

implicit none

!-locl
#include "lenchar.i"
 character(len=LENCHAR) tmpnm,fbuf
 logical opned

fbuf = tmpnm('Buf')
inquire (unit=lubuf,opened=opned)
if (opned) then
    close (lubuf,status='delete')
endif

call openf (lubuf,fbuf,'unknown')

return
end subroutine
!-------------------------------------------------------------------------------
subroutine bufcls ()

use iomodule,        only: stderr, stdin, stdout, lubuf
use errors,          only: ok, prgerr
! closes and removes command buffer file

implicit none

!-locl
#include "lenchar.i"
logical opned
integer :: ierr

ierr = OK
if (lubuf.eq.stdin.or.lubuf.eq.stdout.or.lubuf.eq.stderr) then
    write(stderr,1)
 1        format(1x,'bufcls: fatal logical unit number error')
    ierr = PRGERR
    goto 10
endif

inquire (unit=lubuf,opened=opned)

if (opned) then
    close (unit=lubuf,err=20,status='delete')
endif

10 if (ierr.eq.PRGERR) then
    call exitp(1)
endif

return

20    call prerr ('bufcls')
ierr = PRGERR
goto 10
end
!-------------------------------------------------------------------------------
subroutine bufwrt (lu,file,ierr)

! writes command buffer to a file

use iomodule,        only: stdin, stdout, stderr, lubuf
use errors,          only: ok, prgerr

implicit none

!-pass
integer :: lu, ierr
character file*(*)
!-locl
#include "lenchar.i"
character(len=LENCHAR) :: record
logical                :: opned
integer, external      :: lnblk
integer                :: n

!-init
ierr = OK

if (LEN(file).le.0 .or. lnblk(file).le.0) then
    return
endif

inquire (unit=lubuf,opened=opned)
if (.not.opned) then
    return
endif
if (lu.eq.stdout .or. &
    lu.eq.stdin  .or. &
    lu.eq.stderr) then
    return
endif

rewind(unit=lubuf,err=100)
inquire (unit=lu,opened=opned)

if (opned) then
    call closef(lu)
endif

call openf(lu,file,'unknown')
10    read(lubuf,'(a256)',err=100,end=20) record
n = lnblk(record)

if (n.gt.0) then
    write(lu,'(256a)',err=100) record(1:n)
else
    write(lu,'(a)',err=100)
endif

goto 10

20    call closef (lu)
return

100    call prerr('bufwrt')
ierr = PRGERR
goto 20
end subroutine
!-------------------------------------------------------------------------------
subroutine bufput (cmdline,ierr)

use iomodule,        only: stderr, lubuf
use errors,          only: prgerr, ok

! Puts a command line into buffer
implicit none

!-pass
integer           :: ierr
character(len=*)  :: cmdline
!-locl
logical           :: opned
integer, external :: lnblk
integer           :: n

!-init
ierr = ok

inquire(unit=lubuf,opened=opned)

if (.not.opned) then
    write(stderr,1)
 1        format(1x,'bufput: open error')
    ierr = PRGERR
else
    n = lnblk(cmdline)
    if (n.gt.0) then
        write(lubuf,'(256a)',err=20) cmdline(1:n)
    else    
        write(lubuf,'(a)',err=20)
    endif
endif

10    return
20    call prerr('bufput')
ierr = PRGERR
goto 10
end subroutine
!-------------------------------------------------------------------------------
subroutine cptobuf (lu,file,ierr)

use errors,     only: ok, exec
use iomodule,   only: lubuf

implicit none

!-pass
integer            :: lu, ierr
character(len=*)   :: file
!-locl
character(len=256) :: record
logical            :: opned
integer, external  :: lnblk
integer            :: n
!-init
ierr = OK

inquire(unit=lubuf,opened=opned)

if (opned) then
    rewind(unit=lubuf)
else
    call bufopn()
endif

inquire(unit=lu,opened=opned)

if (opned) then
    call closef(lu)
endif

call openf(lu,file,'old')

10    read(lu,'(a256)',err=100,end=20) record
n = lnblk(record)

if (n.gt.0) then
    write(lubuf,'(256a)',err=100) record(1:n)
else
    write(lubuf,'(a)',err=100)
endif

goto 10

20    return
100    call prerr('cptobuf')
ierr = EXEC
goto 20
end subroutine
!-------------------------------------------------------------------------------
subroutine bufcat (n,ierr)

use iomodule,        only: stderr, lubuf
use errors,          only: ok, exec
! catenates n last lines of buffer. If n < 0 do full buffer

implicit none

!-pass
integer n,ierr
!-locl
character(len=256) :: record
character(len=1)   :: c
logical            :: opned
integer, external  :: lnblk
integer            :: i, nhist, nlb, nr
!-init
ierr = OK

nhist = n

if (nhist.eq.0) then
    return
endif

inquire(unit=lubuf,opened=opned)

if (.not. opned) then
    return
endif

rewind(lubuf,err=50)

!    count number of records in buffer
nr = 0
10    read(lubuf,'(a1)',err=50,end=20) c
nr = nr + 1
goto 10

20  if (nhist.lt.0) then
    nhist=nr
endif

i = 0
rewind(lubuf,err=50)
30    read(lubuf,'(a256)',err=50,end=40) record
i = i + 1
nlb=lnblk(record)

if ((nr-i).lt.nhist) then
    if (nlb.gt.0) then
        write(stderr,'(i6,1x,256a)',err=50) &
               i,record(1:nlb)
    else
        write(stderr,'(i6,1x,a)',err=50) &
               i
    endif
endif
goto 30

40    return

50    call prerr('bufcat')
ierr = EXEC
goto 40
end subroutine
!-------------------------------------------------------------------------------
subroutine bufudo (n,ierr)

use errors,          only: ok, exec
use iomodule,        only: lubuf

! undo n last lines of buffer. If n < 0 undo full buffer

implicit none

!-pass
integer n,ierr



!-locl
#include "lenchar.i"
character(len=1)       :: record
character(len=LENCHAR) :: file
logical                :: opned
integer                :: i, nr, nundo, lueof
!-init
ierr = OK

nundo = n
if (nundo.eq.0) then
    return
endif

inquire(unit=lubuf,opened=opned)

if (.not. opned) then
    return
endif

rewind(lubuf,err=50)

!    count number of records in buffer
nr = 0
10    read(lubuf,'(a1)',err=50,end=20) record
nr = nr + 1
goto 10

20    if (nr.eq.0) then
    return
endif

inquire(unit=lubuf,name=file)

if (nundo.lt.0) then
    nundo=nr
endif

i = 0
rewind(lubuf,err=50)
30    read(lubuf,'(a256)',err=50,end=40) record
i = i + 1

if ((nr-i).gt.nundo) then
    goto 30
endif

40    endfile (unit=lubuf,err=50)
close(unit=lubuf,err=50)
open(unit=lubuf,file=file,err=50)
ierr = lueof(lubuf)

if (ierr.eq.0) then
    ierr = OK
else
    ierr = EXEC
endif

return

50    call prerr('bufudo')
ierr = EXEC
goto 40
end subroutine
