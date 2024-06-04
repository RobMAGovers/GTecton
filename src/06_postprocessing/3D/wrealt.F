!-------------------------------------------------------------------------------
subroutine wrealt (IMPRINT,ITPRINT,itime,lu)

USE CONSTANTS,      only: year
use timestepmodule, only: chktim, nintg, &
                          timestepdatactx
use iomodule,       only: stderr

implicit none
!-pass
integer          :: itime,lu
integer          :: IMPRINT(*),ITPRINT(*)
!-locl
integer          :: i, n, igrp
double precision :: x

logical          :: Mout,Tout
character(len=1) :: C

n    = 0
igrp = 1
i    = 1
x    = 0d0

if (NINTG.le.0) then
    if (itime.gt.0) then
        write(stderr,1)
    endif
 1        format(1x,'NINTG = 0: itime > 0 not available')
    return
endif

if (itime.le.0) then
    do igrp=1,NINTG
        do i=1,timestepdatactx%MAXSTP(igrp)
            n = n + 1
            x = x + timestepdatactx%DELT(igrp)
            if (n.eq.IABS(itime)) then
                C='<'
            else
                C=' '
            endif
            call chktim(n,Mout,Tout)
            if (Mout.and.Tout) then
                write(lu,2) n,x,x/YEAR,x/YEAR*1D-6,C
 2                   format(1x,'Mech+Temp',3x,i5,1x,1Pg12.3,1x,'sec', &
                 3x,1PG12.3,1x,'year',3x,1Pg15.3,1x,'My',1X,A1)
            else if (Mout) then
                write(lu,3) n,x,x/YEAR,x/YEAR*1D-6,C
 3                   format(1x,'Mech     ',3x,i5,1x,1Pg12.3,1x,'sec', &
                 3x,1PG12.3,1x,'year',3x,1Pg15.3,1x,'My',1X,A1)
            else if (Tout) then
                write(lu,4) n,x,x/YEAR,x/YEAR*1D-6,C
 4                   format(1x,'     Temp',3x,i5,1x,1Pg12.3,1x,'sec', &
                 3x,1PG12.3,1x,'year',3x,1Pg15.3,1x,'My',1X,A1)
            endif
        enddo
    enddo
else

300 if (igrp.gt.NINTG) then
        write(stderr,5) itime
 5        format(1x,'itime = ',I5,' not in data')
        goto 500
    endif

400 if (i.gt.timestepdatactx%MAXSTP(igrp)) then
        i    = 1
        igrp = igrp + 1
        goto 300
    endif

    n = n + 1
    x = x + timestepdatactx%DELT(igrp)
    if (n.eq.itime) then
        call chktim(n,Mout,Tout)
        if (Mout.and.Tout) then
            write(lu,6) n,x,x/YEAR,x/YEAR*1D-6
 6               format(1x,'M+T',3x,i5,1x,1Pg12.3,1x,'sec',3x, &
             1PG12.3,1x,'year',3x,1Pg12.3,1x,'My')
        else if (Mout) then
            write(lu,7) n,x,x/YEAR,x/YEAR*1D-6
 7               format(1x,'M  ',3x,i5,1x,1Pg12.3,1x,'sec',3x, &
             1PG12.3,1x,'year',3x,1Pg12.3,1x,'My')
        else if (Tout) then
            write(lu,8) n,x,x/YEAR,x/YEAR*1D-6
 8               format(1x,'  T',3x,i5,1x,1Pg12.3,1x,'sec',3x, &
             1PG12.3,1x,'year',3x,1Pg12.3,1x,'My')
        else
            write(stderr,5) itime
        endif
        goto 500
    endif
i = i + 1
goto 400
endif

500    call lstout(IMPRINT,ITPRINT,timestepdatactx%MAXstp)

return
end
