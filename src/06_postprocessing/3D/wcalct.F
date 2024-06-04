subroutine wcalct (ICALC,lu)

USE MODELDEFINITION, only: mode
USE CONSTANTS,       only: year
use timestepmodule,  only: timestepdatactx, &
                           nintg, &
                           ntcalc

implicit none
!-pass
integer lu
integer :: icalc
dimension ICALC(*)
!-local
integer :: m, n, i, igrp
double precision :: x

!-init
n    = 0
igrp = 1
i    = 1
x    = 0.d0
!
if (MODE.lt.3)   return
if (NTCALC.eq.0) return
if (NINTG.le.0)  return
!
m = 1
do igrp=1,NINTG
    do i=1,timestepdatactx%MAXSTP(igrp)
        n = n + 1
        x = x + timestepdatactx%DELT(igrp)
        if (NTCALC.gt.0 .and. ICALC(m).eq.n) then
            write(lu,1) n,x,x/YEAR,x/YEAR*1D-6
 1               format(1x,i5,1x,1Pg12.3,1x, &
             'sec',3x,1PG12.3,1x,'year',3x,1Pg15.3,1x,'My')
            m = m + 1
            if (m.gt.NTCALC) then
                return
            endif
        else if (NTCALC.lt.0 .and. MOD(n,-NTCALC).eq.0) then
            write(lu,1) n,x,x/YEAR,x/YEAR*1D-6
        endif
    enddo
enddo

return
end subroutine
