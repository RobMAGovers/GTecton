subroutine timinc (IREQ,increment,ierr)

! Increment or decrement NSTEP to next or previous time of output

USE MODELDEFINITION, only: MODE,ICVIS,ISTART,NELAST
use timestepmodule
use iomodule,        only: ifein, ifeout, stderr
use errors,          only: exec
use materialsmodule, only: lmat


implicit none

!-pass
integer IREQ,ierr
logical increment

! todo where does this stuff come from?
logical verbose,append,Mout,Tout

!-locl
logical Mech,Temp
integer :: i, iend, istp
!-init
IREQ  = 1

if (ifein.eq.0) then
    write(stderr,1)
 1        format(1x,'set itime: first read FE-input file')
    ierr = EXEC
    goto 300
endif

if (ifeout.eq.0) then
    write(stderr,2)
 2        format(1x,'set itime: first read FE-output file')
    ierr = EXEC
    goto 300
endif

if (NINTG.le.0) then
    if (increment) then
        write(stderr,3)
 3            format(1x,'set itime ++: cannot increment')
        ierr = EXEC
    else
        write(stderr,4)
 4            format(1x,'set itime --: cannot decrement')
        ierr = EXEC
    endif
    goto 300
endif

!    determine last integer output time IEND
call LSTSTP (NINTG,MODE,NMPRT,NTPRT,Iend)

if (increment) then
    I = NSTEP + 1
    Istp = 1
else
    I = NSTEP - 1
    Istp = -1
endif

100    if (I.ge.0 .and. I.le.Iend) then

!        time step I in output arrays ?
    call chktim(I,Mech,Temp)

    if (Mech .or. Temp) then
        goto 200
    endif

    I = I + Istp
    goto 100
endif

if (increment) then
    write(stderr,3)
    ierr = EXEC
else
    write(stderr,4)
    ierr = EXEC
endif

goto 300

200    Mout = Mech
Tout = Temp
IREQ = I

300    return
end subroutine
!-------------------------------------------------------------------------------
subroutine LSTSTP (NINTG,MODE,NMPRT, NTPRT,IEND)

use timestepmodule, only: outputcontroldatactx, &
                          timestepdatactx 

implicit none

! Determines last integer output time Iend

!-pass
integer NINTG,MODE,NMPRT,NTPRT,IEND
!-local
integer :: im, it, last, i
!-init
IEND = 0
IM = 0
IT = 0

if (NINTG.le.0) then
    return
endif

Last = 0
do i=1,NINTG
    Last = Last + timestepdatactx%MAXSTP(i)
enddo

if (Last.le.0) then
    return
endif

if ((MODE.ne.2 .and. MODE.ne.6).or.NMPRT.eq.0) then
    goto 300
endif

if (NMPRT.gt.0) then
!        loop IMPRINT until either LAST or IMPRINT(NMPRT) is reached
    I  = 0
100 if (I.eq.NMPRT) then
        goto 200
    endif
    I  = I + 1
    IM = outputcontroldatactx%IMPRINT(I)
    if (IM.eq.Last) then
        goto 200
    endif
    if (IM.gt.Last) then
        IM = 0
        I  = I - 1
        if (I.gt.0) then
            IM = outputcontroldatactx%IMPRINT(I)
        endif
        goto 200
    endif
    goto 100

200        continue

else
    IM = INT( REAL(Last) / REAL(NMPRT) ) * NMPRT
endif

300    if ((MODE.ne.5 .and. MODE.ne.6).or.NTPRT.eq.0) then
    goto 600
endif

if (NTPRT.gt.0) then
!        loop ITPRINT until either LAST or ITPRINT(NTPRT) is reached
    I  = 0
400 if (I.eq.NTPRT) then
        goto 500
    endif
    I  = I + 1
    IT = outputcontroldatactx%ITPRINT(I)
    if (IT.eq.Last) then
        goto 400
    endif
    if (IT.gt.Last) then
        IT = 0
        I  = I - 1
        if (I.gt.0) then
            IT = outputcontroldatactx%ITPRINT(I)
        endif
        goto 400
    endif
    goto 400

500        continue
else
    IT = INT( REAL(Last) / REAL(NTPRT) ) * NTPRT
endif

600    IEND = MAX0(IM,IT)
return
end subroutine
