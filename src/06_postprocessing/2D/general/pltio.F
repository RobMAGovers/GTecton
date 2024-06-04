subroutine OfClose (ierr,fname,monitor,lu,luMon)

use errors,          only: exec, ok

implicit none

!-pass
integer ierr,lu,luMon
logical monitor
character(len=*) fname


integer, external :: lnblk
!-local
integer :: luloc


luloc = lu
if (monitor) then
    luloc=luMon
endif

if (ierr.eq.OK) then
    if (monitor) then
        if (lnblk(fname).gt.0) then
            close(luMon,err=100)
            goto 200
100         call prerr('OfClose')
            ierr = EXEC
200         continue
        endif
    else
        close(luloc,err=300)
        goto 400
300     call prerr('OfClose')
        ierr = EXEC
400     continue
    endif
endif
return
end subroutine
!-------------------------------------------------------------------------------
subroutine OfOpen &
 (ierr,fname,ascii,monitor,append,triang,lu,luMon)

use iomodule,        only: stderr, stdin, stdout, luout
use errors,          only: exec, ok, prgerr, syntax

implicit none

!-pass
 integer ierr,lu,luMon
 logical ascii,monitor,append,triang
 character(len=*) fname


!-locl
 integer, external :: lnblk, lueof
integer :: l


 ierr = OK

 if (triang) then
goto 1000
 endif

 l = LEN(fname)
 if (l.gt.0) then
l = lnblk(fname)
 endif

 if (monitor) then
if (l.gt.0) then
    call shvar(fname)
!        redirect monitor data to (ascii) file
    luMon = luout(2)
    open(unit=luMon,file=fname,err=100,iostat=ierr)

    if (append) then
!                file pointer to end of file
        ierr = lueof(luMon)
    else
        rewind(luMon)
    endif

100            if (ierr.ne.0) then
        call prerr('OfOpen')
        ierr  = EXEC
        luMon = stdout
    endif

else
    luMon = stdout
endif

goto 1000
 endif

 if (lnblk(fname).le.0) then
write(stderr,1)
 1        format(1x,'file argument missing')
ierr = SYNTAX
goto 1000
 endif

 if (lu.eq.stdin.or.lu.eq.stdout.or.lu.eq.stderr) then
write(stderr,2) lu
 2        format(1x,'OfOpen: unit ',I1,' cannot be opened')
ierr = PRGERR
goto 1000
 endif

 call shvar(fname)

 if (ascii) then
open(unit=lu, &
     file=fname, &
     err=200, &
     status='unknown', &
     form='formatted', &
     iostat=ierr)
 else
open(unit=lu, &
     file=fname, &
     err=200, &
     status='unknown', &
     form='unformatted', &
     iostat=ierr)
 endif

 if (append) then
!        file pointer to end of file
ierr = lueof(lu)
 else
rewind(lu)
 endif

200    if (ierr.ne.0) then
call prerr('OfOpen')
ierr  = EXEC
luMon = stdout
 endif

1000    return
 end
!-------------------------------------------------------------------------------
integer function openform (lu,f,ascii)

use errors,          only: exec

implicit none

#include "lenchar.i"
!-pass
integer lu
character f*(*)
logical ascii
!-locl
logical AsciF

integer :: ierr
integer, external :: ifopen, ifbopen

ascii = AsciF(f)

!      write(*,*) "Input ", f
!      if (ascii) then
!          write(*,*) "is ASCII"
!      else
!          write(*,*) "not ASCII"
!      endif

if (ascii) then
    ierr = ifopen (lu,f,'old')
else
    ierr = ifbopen(lu,f,'old')
endif

if (ierr.ne.0) then
    ierr = EXEC
endif

50    openform = ierr
return
end function
!-------------------------------------------------------------------------------
 subroutine svdSolve (A,B,IDIAG,NEQ,FACT,BACK)

! dummy routine for loading "elastc.f"

implicit none
!-pass
integer :: NEQ
double precision :: A
dimension A(*)
double precision :: B(NEQ),IDIAG(NEQ)
logical FACT,BACK

return
end

!-------------------------------------------------------------------------------
integer function iSoFil (FIN,FOUT,NF,MAXF)

use iomodule,        only: stderr
use errors,          only: exec, ok, prgerr

! retrieves all files nested in FIN (by ".so")

implicit none
!-pass
integer MAXF,NF,ierr
character FIN*(*),FOUT*(*)
dimension FOUT(MAXF)
!-locl
integer, parameter :: maxnest=10
#include "lenchar.i"
logical there
integer :: fnblnk,lui(maxnest)
character(len=LENCHAR) tmp,file,dirname,home
character(len=255)record
integer, external :: lnblk,nextlu
integer :: m, n, l, ilevel

!-init
 lui(1) = nextlu(0)
 NF     = 0
 ierr   = OK
 ilevel = 1

 l = LEN(FIN)

 if (l.le.0) then
write(stderr,1)
 1        format(1x,'iSoFil: length of FIN error')
ierr = PRGERR
goto 1000
 endif

 m = LEN(FOUT(1))

 if (m.le.l) then
write(stderr,2)
 2        format(1x,'iSoFil: length of FOUT error')
ierr = PRGERR
goto 1000
 endif

 n = fnblnk(FIN)
 m = lnblk(FIN)

 if (n.ge.m) then
write(stderr,3)
 3        format(1x,'iSoFil: empty input file error')
ierr = EXEC
goto 1000
 endif

 file = FIN(n:m)

100    inquire (file=file,exist=there)

 if (.not.there) then
write(stderr,4) '"',file(n:m),'"'
 4        format(1x,'iSoFil: file does''nt exist: ',80A)
ierr = EXEC
goto 1000
 endif

 open(unit=lui(ilevel),status='old',file=file,err=2000)
 rewind(lui(ilevel),err=2000)
 NF = NF + 1

 if (NF.gt.MAXF) then
write(stderr,5) MAXF
 5       format(1x,'iSoFil: number of include files exceeds',1x, &
     'maximum = ',i3)
ierr = PRGERR
goto 1000
 endif

 FOUT(NF) = file(n:m)
200    if (ilevel.gt.0) then
300        read(lui(ilevel),6,err=2000,end=400) record
 6        format(A80)
if (record(1:4).ne.'.so ') goto 300
! ".so "
!        increment nesting level
ilevel = ilevel + 1

if (ilevel.gt.maxnest) then
    write(stderr,7) maxnest
 7        format(1x,'iSoFil: maximum ".so" nesting depth = ',i3)
    ierr = EXEC
    goto 1000
endif

!        read file name
l = lnblk(record)

if (l.le.4) then
    write(stderr,8) '"',file(n:m),'"'
 8        format(1x,'iSoFil: ".so "-error in ',80A)
    ierr = EXEC
    goto 1000
endif

lui(ilevel) = nextlu(0)
tmp = record(5:l)
n   = fnblnk(tmp)
m   = l - 4

if (tmp(n:n).eq.'/') then
    file = tmp(n:m)
else
    home = dirname(FIN)
    file = home(1:lnblk(home)) // '/' // tmp(n:m)
    n    = fnblnk(file)
    m    = lnblk(file)
endif

goto 100
! "EOF"
400        close(unit=lui(ilevel),err=2000)
ilevel = ilevel - 1
goto 200
 endif

1000    iSoFil = ierr
 return
2000    write(stderr,9) '"',file(n:m),'"'
 9    format(1x,'iSoFil: I/O error on file ',80a)
 call prerr('iSoFil')
 ierr = EXEC
 goto 1000
 end
