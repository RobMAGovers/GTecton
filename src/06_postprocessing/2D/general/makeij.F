logical function makeij (IJFILE)

use debugmodule,     only: iecho
use iomodule,        only: fein, stderr, stdout
use errors,          only: ok

implicit none
!-pass
character(len=*) IJFILE
!-locl
integer, parameter :: MAXFIN=20
character(len=255) FOUT(MAXFIN)
!double precision FOUT(MAXFIN)
logical there
integer statb(13),cstat,fnblnk
integer, external :: lnblk, isofil

integer :: i, j, k, l, m, n, ierr
integer :: nf
integer :: modOut, modIn

!-init

data FOUT(1),NF,statb(1)/' ',0,0/

makeij = .false.

n = fnblnk(IJFILE)
m = lnblk(IJFILE)

inquire (file=IJFILE,exist=there)
if (.not. there) then
    if (iecho.eq.1) then
        write(stdout,1) '"',IJFILE(n:m), &
       '" doesn''t exist'
 1        format(1x,100A)
    endif
    goto 300
endif

ierr = cstat(IJFILE(n:m),statb)

if (ierr.ne.0) then
    write(stderr,1) 'makeij: cstat error on "', &
      IJFILE(n:m),'"'
    goto 300
endif

MODOUT = statb(10)

ierr = iSoFil (FEIN,FOUT,NF,MAXFIN)

if (ierr.ne.OK) then
    write(stderr,1) 'makeij: iSoFil error'
    goto 300
endif

if (NF.le.0) then
    write(stderr,1) 'makeij: no files'
    goto 300
endif

do i=1,NF
    k = fnblnk(FOUT(i))
    l = lnblk(FOUT(i))
    inquire (file=FOUT(i),exist=there)
    if (.not.there) then
        write(stderr,1) 'makeij: file "', &
       FOUT(i)(k:l),'" doesn''t exist'
        goto 300
    endif

    ierr = cstat(FOUT(i)(k:l),statb)

    if (ierr.ne.0) then
        write(stderr,1) 'makeij: cstat error on "', &
       FOUT(i)(k:l),'"'
        goto 300
    endif

    MODIN = statb(10)

    if (MODOUT.lt.MODIN) then
        if (iecho.eq.1) then
              write(stderr,1) 'Creating "',IJFILE(n:m), &
           '" since its older than "',FOUT(i)(k:l), &
           '"'
        endif
        makeij = .true.
    endif
enddo

200    return
300    makeij = .true.
goto 200
end function
