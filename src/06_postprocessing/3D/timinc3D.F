subroutine timinc (IREQ,increment,ierr)

! Increment or decrement NSTEP to next or previous time of output

USE MODELDEFINITION, only: MODE,ICVIS,ISTART,NELAST
use timestepmodule
use iomodule,        only: ifein, ifeout, stderr
use errors,          only: exec

implicit none

!-pass
 integer IREQ,ierr
 logical increment

 logical verbose,append,Mout,Tout

!-locl
 logical Mech,Temp
integer :: i, istp, iend
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

if (Mech .or. Temp) goto 200
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
 end
!-------------------------------------------------------------------------------
 subroutine LSTSTP (NINTG,MODE,NMPRT, NTPRT,IEND)

! Determines last integer output time Iend

use timestepmodule, only: outputcontroldatactx,    &
                                   timestepdatactx

implicit none

!-pass
 integer NINTG,MODE,NMPRT,NTPRT,IEND

!-local
integer :: im, it, last, i


!-init
 IEND = 0
 IM = 0
 IT = 0

 if (NINTG.le.0) return

 Last = 0

do i=1,NINTG
    Last = Last + timestepdatactx%MAXSTP(i)
enddo

 if (Last.le.0) return

 if ((MODE.ne.2 .and. MODE.ne.6).or.NMPRT.eq.0) goto 300

 if (NMPRT.gt.0) then
!        loop IMPRINT until either LAST or IMPRINT(NMPRT) is reached
I  = 0
100        if (I.eq.NMPRT) goto 200
I  = I + 1
IM = outputcontroldatactx%IMPRINT(I)
if (IM.eq.Last) goto 200
if (IM.gt.Last) then
IM = 0
I  = I - 1
if (I.gt.0) IM = outputcontroldatactx%IMPRINT(I)
goto 200
endif
goto 100

200        continue
 else
IM = INT( REAL(Last) / REAL(NMPRT) ) * NMPRT
 endif

300    if ((MODE.ne.5 .and. MODE.ne.6).or.NTPRT.eq.0) goto 600

 if (NTPRT.gt.0) then
!        loop ITPRINT until either LAST or ITPRINT(NTPRT) is reached
I  = 0
400        if (I.eq.NTPRT) goto 500
I  = I + 1
IT = outputcontroldatactx%ITPRINT(I)
if (IT.eq.Last) goto 400
if (IT.gt.Last) then
IT = 0
I  = I - 1
if (I.gt.0) IT = outputcontroldatactx%ITPRINT(I)
goto 400
endif
goto 400

500        continue
 else
IT = INT( REAL(Last) / REAL(NTPRT) ) * NTPRT
 endif

600    IEND = MAX0(IM,IT)
 return
 end
