subroutine centrs (IEN,X,D,LMF,TFAULT,DXE, &
                   SKEW,NEN,NSD,NDOF,DEFORM,deffac,lu,ierr,ibase,istep, &
                   verbose,luMon)

USE MODELDEFINITION
use meshdatamodule,  only: meshdatactx
use iomodule,        only: ascout
use constants,       only: zero
use postProcessData, only: GMT,connout
use monitorVars,     only: elmon, monel, nnel
use errors,          only: ok
implicit none

!-pass
logical           :: DEFORM,verbose
integer           :: NEN,NSD,NDOF,lu,ierr,ibase,istep,luMon
double precision  :: deffac

integer           :: ien, lmf
double precision  :: x, d, tfault, dxe, skew
dimension         :: IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*), &
                     TFAULT(NDOF,*),DXE(NDOF,NEN,*),SKEW(2,*)
!-locl
double precision  :: def
double precision  :: CENTR(3)
real              :: D2REAL
integer           :: i, n, nn, npts
integer, external :: lnblk

!-init
ierr = OK
def  = ZERO

if (DEFORM) then
    def = deffac
endif

if (elmon) then
    do nn=1,nnel
        n = monel(nn)

        CENTR(1) = meshdatactx%X(1,meshdatactx%IEN(1,n)) + &
                   meshdatactx%X(1,meshdatactx%IEN(2,n)) + &
                   meshdatactx%X(1,meshdatactx%IEN(3,n))  
        CENTR(2) = meshdatactx%X(2,meshdatactx%IEN(1,n)) + &
                   meshdatactx%X(2,meshdatactx%IEN(2,n)) + &
                   meshdatactx%X(2,meshdatactx%IEN(3,n)) 

        if (NSD.eq.3) then
            CENTR(1) = CENTR(1) + meshdatactx%X(1,meshdatactx%IEN(4,n))
            CENTR(2) = CENTR(2) + meshdatactx%X(2,meshdatactx%IEN(4,n))

            CENTR(3) = meshdatactx%X(3,meshdatactx%IEN(1,n)) + &
                       meshdatactx%X(3,meshdatactx%IEN(2,n)) + &
                       meshdatactx%X(3,meshdatactx%IEN(3,n)) + &
                       meshdatactx%X(3,meshdatactx%IEN(4,n))

            CENTR(1) = 0.25 * CENTR(1)
            CENTR(2) = 0.25 * CENTR(2)
            CENTR(3) = 0.25 * CENTR(3)
        else
            CENTR(1) = 0.333333333333 * CENTR(1)
            CENTR(2) = 0.333333333333 * CENTR(2)
        endif 


        if (verbose) then
            write(luMon,2) n,(CENTR(i),i=1,NSD)
 2            format(1x,'center of element ',I5,': ',3(1PG12.4))
        else
            write(luMon,*) (CENTR(i),i=1,NSD)
        endif

    enddo
else

    npts = INT((NUMEL-ibase)/istep)+1

    if ((.not.GMT).and.ascout) then
        write(lu,*,err=100) 2,npts
    endif

    do nn=ibase,NUMEL,istep
        n = nn

! LW replaced by just center computation
        CENTR(1) = meshdatactx%X(1,meshdatactx%IEN(1,n)) + &
                   meshdatactx%X(1,meshdatactx%IEN(2,n)) + &
                   meshdatactx%X(1,meshdatactx%IEN(3,n))
        CENTR(2) = meshdatactx%X(2,meshdatactx%IEN(1,n)) + &
                   meshdatactx%X(2,meshdatactx%IEN(2,n)) + &
                   meshdatactx%X(2,meshdatactx%IEN(3,n)) 

        if (NSD.eq.3) then
            CENTR(1) = CENTR(1) + meshdatactx%X(1,meshdatactx%IEN(4,n))
            CENTR(2) = CENTR(2) + meshdatactx%X(2,meshdatactx%IEN(4,n))

            CENTR(3) = meshdatactx%X(3,meshdatactx%IEN(1,n)) + &
                       meshdatactx%X(3,meshdatactx%IEN(2,n)) + &
                       meshdatactx%X(3,meshdatactx%IEN(3,n)) + &
                       meshdatactx%X(3,meshdatactx%IEN(4,n))

            CENTR(1) = 0.25 * CENTR(1)
            CENTR(2) = 0.25 * CENTR(2)
            CENTR(3) = 0.25 * CENTR(3)
        else
            CENTR(1) = 0.333333333333 * CENTR(1)
            CENTR(2) = 0.333333333333 * CENTR(2)
        endif

        if (GMT) then
            if (ascout) then
                write(lu,*,err=100) (CENTR(i),i=1,NSD)
            else
                write(lu,err=100) (D2REAL(CENTR(i)),i=1,NSD)
            endif
        else
            stop "non-GMT plotting not implemented"
        endif
    enddo
endif
!
 return
100    call prerr('centrs')

call exitp(1)

end subroutine

!-------------------------------------------------------------------------------
subroutine echo (cmdln,cmd,ncmd)

use iomodule, only: stderr, stdout

implicit none

!-pass
character :: cmd*(*),cmdln*(*)
dimension :: cmd(*)
!-local
integer   :: ncmd, k, lnb, n, m

integer, external :: lnblk, indx

if (ncmd.gt.1) then
    k=2

    if (cmd(2)(1:3).eq.'-n ') then
        k=3
    endif

    if (ncmd.ge.k) then
        lnb = lnblk(cmd(k))
        if (lnb.le.0) then
            write(stderr,*) ' echo: parser error'
            call exitp(1)
        endif

        n = indx(cmdln,cmd(k)(1:lnb))

        if (n.gt.0) then
            m = lnblk(cmdln)

            if (cmdln(n:n).eq.'"' .or. cmdln(n:n).eq.'''') then
                n=n+1
            endif

            if (cmdln(m:m).eq.'"' .or. cmdln(m:m).eq.'''') then
                m=m-1
            endif

            if (n.le.m) then
                if (k.eq.2) then
                    write(stdout,1) cmdln(n:m)
 1                format(1x,256A)
                else
                    write(stdout,2,advance='no') cmdln(n:m)
 2                format(1x,A)
                endif
            else
                write(stdout,4)
            endif
        else
            write(stderr,3)
 3            format(1x,'echo: cmdln indx error')
        endif
    endif
else
    write(stdout,4)
 4        format(1X)
endif

return
end subroutine

!-------------------------------------------------------------------------------

integer function iprinfo(cmd,ncmd)

use documentation, only: plnplt_manual

implicit none

!     Print command info

!-pass
integer ncmd
character cmd*(*)
dimension cmd(ncmd)



if (ncmd.le.1) then
!       single "?": give list of commands
    call plnplt_manual('list')
else
!       information on a specific command requested
    call plnplt_manual(cmd(2))
endif

return
end function

!-------------------------------------------------------------------------------

subroutine nodals (IEN,X,D,NSLIP,DXE,NUMBRS,DEFORM,deffac,lufe, &
                   luMon,VERBOSE,TOTAL,TECTON)

!    USE MODELDEFINITION, only:  RADIUS,NUMSN
USE MODELDEFINITION
USE CONSTANTS,       only: deg2rad
USE MESHDATAMODULE
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: stderr
use postProcessData, only: GMT,connout
use monitorVars,     only: monnp, npmon, nnp
use plotControl,     only: jShell


!
implicit none

!-pass
integer          :: lufe,luMon
double precision :: deffac

integer          :: nslip, ien
double precision :: x, d, dxe
dimension           X(NSD,*),D(NDOF,*),NSLIP(5,*), &
                    DXE(NDOF,NEN,*), IEN(NEN,*)

logical          :: NUMBRS,DEFORM,VERBOSE,TOTAL,TECTON, deformed
!-locl
logical          :: slippery
double precision :: XL(24),DL(24),XX(3)
double precision :: def, ylat
integer          :: i, j, k, kk, ll, local, luo, n, ne, nout, npnum

!-init
XX = 0.
!
if (npmon) then
    npnum = nnp
    luo   = luMon
else
    npnum = meshdatactx%Nvglobal
    luo   = lufe
endif
!
if (DEFORM) then
    def = deffac
else
    def = 0d0
endif
!
if (NUMFN.gt.0 .and. DEFORM .and. TOTAL) then
    write(stderr,*) 'nodals: faulted nodes not implemented'
endif
!

if (.not.NUMBRS.and..not.VERBOSE.and..not.GMT) then
    if (NUMSN.gt.0 .and. DEFORM .and. TOTAL) then
!            count number of outputs
        nout = 0
        do i=1,npnum
            n = i

            if (npmon) then
                n = monnp(i)
            endif
 
            slippery = .false.

            do j=1,NUMSLP
                if (NSLIP(2,j).eq.n) then
                    slippery = .true.
                    nout = nout + 1
                endif
            enddo

            if (.not.slippery) then
                nout = nout + 1
            endif
        enddo

        write(luo,*) 2,nout
    else
        write(luo,*) 2,npnum
    endif
endif
!

if (NUMSN.gt.0 .and. DEFORM .and. TOTAL) then
    do i=1,npnum
        n = i

        if (npmon) then 
            n = monnp(i)
        endif

        slippery = .false.
!            is this node slippery?
        do j=1,NUMSLP
            if (NSLIP(2,j).eq.n) then
                slippery = .true.
                ne = NSLIP(1,j)

                call LCOORD (X,XL,IEN(1,ne))
                call LDISP (DL,D,IEN(1,ne),NDOF,NEN)
                call ADDSNE (DL,ne,NDOF,NEN,0)

                do kk=1,NEN
                    local = kk
                    if (IEN(local,ne).eq.n) then
                        goto 100
                    endif
                enddo

100             ll = (local-1)*NDOF

                if (JSHELL.eq.0) then
                    XX(1) = REAL(xl(ll+1) + def*dl(ll+1))
                    XX(2) = REAL(xl(ll+2) + def*dl(ll+2))
                else
                    ylat = xl(ll+2)*DEG2RAD
                    XX(1) = REAL(XL(ll+1) + def*dl(ll+1) / &
                            (DEG2RAD*RADIUS*COS(ylat)) )
                    XX(2) = REAL(XL(ll+2) + def*dl(ll+2) / &
                            (DEG2RAD*RADIUS) )
                endif

                if (NDOF.gt.2) then
                    XX(3) = REAL(def*dl(ll+3))
                endif

                if (NUMBRS) then
                    write(luo,1) XX(1),XX(2),n
 1                        format('pst ',2(1PG14.6),1x,I5)
                else
                    if (VERBOSE) then
                        write(luo,2) n,(XX(k),k=1,NDOF)

!                            ************ write actual file (np.dat) *******
 2                           format(' nodal point ',I5,'  ', &
                         3(1PG14.6))
                    else
                        write(luo,*) XX(1),XX(2)
                    endif
                endif
            endif
        enddo

        if (.not.slippery) then
            if (JSHELL.eq.0) then
                XX(1) = REAL( X(1,n)+def*D(1,n) )
                XX(2) = REAL( X(2,n)+def*D(2,n) )
            else    
!                    surface is rounded, compensate.
                ylat = X(2,n)*DEG2RAD
                XX(1) = REAL( X(1,n) + def*D(1,n) / &
                       (DEG2RAD*RADIUS*COS(ylat)) )
                XX(2) = REAL( X(2,n) + def*D(2,n) / &
                       (DEG2RAD*RADIUS) )
            endif

            if (NDOF.eq.3) then
                XX(3) = REAL(def*D(3,n))
            endif

            if (NUMBRS) then
                write(luo,1) XX(1),XX(2),n
            else
                if (VERBOSE) then
!           ************ write actual file (np.dat) *******                        
                    write(luo,2) n,(XX(k),k=1,NDOF)
                else
                    write(luo,*) XX(1),XX(2)
                endif
            endif
        endif
    enddo
else
    do i=1,npnum
        n = i
        if (npmon) then 
            n = monnp(i)
        endif

        if (JSHELL.eq.0) then
            XX(1) = REAL( X(1,n)+def*D(1,n) )
            XX(2) = REAL( X(2,n)+def*D(2,n) )
        else
            ylat = X(2,n)*DEG2RAD
            XX(1) = REAL( X(1,n) + def*D(1,n) / &
                      (DEG2RAD*RADIUS*COS(ylat)) )
            XX(2) = REAL( X(2,n) + def*D(2,n) / &
                      (DEG2RAD*RADIUS) )
        endif

        if (NDOF.eq.3) then 
            XX(3) = REAL(def*D(3,n))
        endif

        if (TECTON) then
            write(luo,3) n,0,(XX(k),k=1,NDOF)
 3                     format(2I5,1PG14.6,1PG14.6,1PG14.6)
        else
            if (NUMBRS) then
                write(luo,1) XX(1),XX(2),n
            else
                if (VERBOSE) then
!               ************ write actual file (np.dat) *******
                    write(luo,2) n,(XX(k),k=1,NDOF)
                else
                    write(luo,*) XX(1),XX(2)
                endif
            endif
        endif
    enddo

    if (TECTON) then
        write(luo,4)
 4        format('end deformed nodal coordinates')
    endif
endif

return
end subroutine

!-------------------------------------------------------------------------------

subroutine tstria (XPT,NPLT,NUMTR,f)

use meshdatamodule,  only: plotDOFperElem
use debugmodule,     only: iecho, &
                           parallelLog
use iomodule,        only: ascout, luout, stderr
use errors,          only: exec, ok
! writes triangle coordinates that are contained in XPT

implicit none

!-pass
integer           :: NPLT,NUMTR
double precision  :: XPT(plotDOFperElem,NPLT,NUMTR)
!character f*(*)
character(len=*)  :: f


!-locl
real              :: tst(4),D2REAL
integer           :: i, j, k
integer           :: ndata
integer           :: ierr, itype
integer, external :: ifopen, ifbopen, lnblk ! from tools

if (iecho.eq.8) then
    call ParallelLog("tstria", "Entered")
endif

!-init
tst(1) = 3.0
tst(2) = 2.0
tst(3) = 2.0
tst(4) = 2.0
ierr   = OK

if (ascout) then
    ierr = ifopen  (luout(1),f,"unknown")
else
    ierr = ifbopen (luout(1),f,"unknown")
endif

if (ierr.ne.0) then
    write(*,*) "rstria could not open output file; error: ", ierr
    stop "Exiting"
    ierr = EXEC
endif

if (NUMTR.gt.0 .and. ierr.eq.OK) then
    ndata = NUMTR*(NPLT+1)

    if (plotDOFperElem.eq.3) then
        itype = 4
    else if (plotDOFperElem.eq.4) then
        itype = 6
    else
        write(stderr,1)
 1        format(1x,'tstria: NPLDOF error')
        call exitp(1)
    endif

    if (ascout) then

        if (iecho.eq.8) then
            call ParallelLog("tstria", "writing ascii data")
        endif

        write(luout(1),*,err=200) REAL(itype),REAL(ndata)

          if (itype.eq.4) then
            write(luout(1),2,err=200) &
         (((D2REAL(XPT(i,j,k)),i=1,plotDOFperElem-1),tst(j), &
           j=1,NPLT), (D2REAL(XPT(i,1,k)),i=1,plotDOFperElem-1), &
           tst(4),k=1,NUMTR)
 2            format(1x,3(1PG15.6))
        else
            write(luout(1),3,err=200) &
         (((D2REAL(XPT(i,j,k)),i=1,plotDOFperElem-1),tst(j), &
           j=1,NPLT), (D2REAL(XPT(i,1,k)),i=1,plotDOFperElem-1), &
           tst(4),k=1,NUMTR)
 3            format(1x,4(1PG15.6))
        endif

        if (iecho.eq.8) then
            call ParallelLog("tstria", "Finished writing ascii data")
        endif


    else

        if (iecho.eq.8) then
            call ParallelLog("tstria", "Writing binary data")
        endif


        write(luout(1),err=200) REAL(itype),REAL(ndata)
        write(luout(1),err=200) &
     (((D2REAL(XPT(i,j,k)),i=1,plotDOFperElem-1),tst(j),j=1,NPLT), &
       (D2REAL(XPT(i,1,k)),i=1,plotDOFperElem-1),tst(4),k=1,NUMTR)

        if (iecho.eq.8) then
            call ParallelLog("tstria", "Finished writing binary data")
        endif
    endif
endif

call closef(luout(1))

write(stderr,4) '"',f(1:lnblk(f)),'"'
 4    format(1x,'EOF ',80a)

99 if (iecho.eq.8) then
    call ParallelLog("tstria", "Done")
endif

return

200    call prerr('tstria')

ierr = EXEC
goto 99

end subroutine

!-------------------------------------------------------------------------------

subroutine wait (ierr)

use iomodule,        only: stdin, stdout, ilu, lucmd
use errors,          only: exec, stop, eof, ok

implicit none

!-pass
integer          :: ierr



!-locl
character(len=1) :: c
!-init
ierr = OK

10    write(stdout,1,advance='no')
 1    format(1x,'press "return" to continue, "q" to quit > ')

read(stdin,'(A)',err=20,end=30) c

if (c.eq.'q') then
    ierr = STOP
    if (lucmd(ilu).ne.stdin) then
        ierr=EOF
    endif
endif

20    return
30    ierr = STOP

goto 20
end subroutine

!-------------------------------------------------------------------------------
subroutine wrealt (itime, lu)

USE TIMESTEPMODULE
USE CONSTANTS,     only: year
use debugmodule,   only: iecho
use iomodule,      only: stderr, stdout


implicit none
!-pass
integer itime,lu

!-locl
logical            :: Mout,Tout
character(len=256) :: ftmp
character(len=1)   :: C
character(len=350) :: cmd
integer            :: i, igrp, n, ierr, lutmp
double precision   :: x

integer, external  :: csystem, lnblk, nextlu
character(len=256), external :: tmpnm



if (NINTG.le.0) then
    if (itime.gt.0) then
        write(stderr,1)
 1        format(1x,'NINTG = 0: itime > 0 not available')
    endif
    return
endif

if (itime.le.0) then
    if (lu.eq.stdout) then
        lutmp = nextlu(10)
        ftmp = tmpnm('wrealt')
        if (iecho.gt.1) then
            write(stderr,99) 'call openf(', &
                              ftmp(1:lnblk(ftmp)),')'
        endif
        call openf(lutmp,ftmp,'unknown')
    else
        lutmp = lu
    endif

    n = 0
    x = 0d0

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
                write(lutmp,2) n,x,x/YEAR,x/YEAR*1D-6,C
 2                   format(1x,'Mech+Temp',3x,i5,1x,1Pg12.3,1x,'sec', &
                 3x,1PG12.3,1x,'year',3x,1Pg15.3,1x,'Ma',1X,A1)
            else if (Mout) then
                write(lutmp,3) n,x,x/YEAR,x/YEAR*1D-6,C
 3                   format(1x,'Mech     ',3x,i5,1x,1Pg12.3,1x,'sec', &
                 3x,1PG12.3,1x,'year',3x,1Pg15.3,1x,'Ma',1X,A1)
            else if (Tout) then
                write(lutmp,4) n,x,x/YEAR,x/YEAR*1D-6,C
 4                   format(1x,'     Temp',3x,i5,1x,1Pg12.3,1x,'sec', &
                 3x,1PG12.3,1x,'year',3x,1Pg15.3,1x,'Ma',1X,A1)
            endif

        enddo
    enddo

    if (lu.eq.stdout) then
        call closef(lutmp)
        if (n.gt.0) then
            cmd = 'more ' // ftmp(1:lnblk(ftmp))

            if (iecho.gt.1) then
                write(stderr,99) 'call system(', &
               cmd(1:lnblk(cmd)),')'
            endif

            ierr = csystem(cmd)

            if (ierr.ne.0) then
                write(stderr,5) ftmp(1:lnblk(ftmp)),'"'
 5                format(1x,'system call failed: "more ',80a)
            endif

            if (ierr.eq.0) then
                if (iecho.gt.1) then
                     write(stderr,99) 'call rm(', &
                               ftmp(1:lnblk(ftmp)),')'
                endif
                call rm(ftmp)
            endif
        endif
    endif
else
    n = 0
    x = 0d0
    igrp = 1
    do while (igrp.le.NINTG)
        i = 1
        do while (i.le.timestepdatactx%MAXSTP(igrp))
            n = n + 1
            x = x + timestepdatactx%DELT(igrp)
            if (n.eq.itime) then
               call chktim(n,Mout,Tout)

                if (Mout.and.Tout) then
                    write(lu,6) n,x,x/YEAR,x/YEAR*1D-6
 6                       format(1x,'M+T',3x,i5,1x,1Pg12.3,1x,'sec', &
                     3x,1PG12.3,1x,'year',3x,1Pg12.3,1x,'Ma')
                else if (Mout) then
                    write(lu,7) n,x,x/YEAR,x/YEAR*1D-6
 7                       format(1x,'M  ',3x,i5,1x,1Pg12.3,1x,'sec', &
                     3x,1PG12.3,1x,'year',3x,1Pg12.3,1x,'Ma')
                else if (Tout) then
                    write(lu,8) n,x,x/YEAR,x/YEAR*1D-6
 8                       format(1x,'  T',3x,i5,1x,1Pg12.3,1x,'sec', &
                     3x,1PG12.3,1x,'year',3x,1Pg12.3,1x,'Ma')
                else
                    write(stderr,5) itime
                endif

                goto 100
            endif
            i = i + 1
        enddo
        i    = 1
        igrp = igrp + 1
    enddo
    write(stderr,9) itime
 9        format(1x,'itime = ',I5,' not in data')
    goto 100
endif
!
100    call lstout(outputcontroldatactx%IMPRINT, &
                   outputcontroldatactx%ITPRINT, &
                   timestepdatactx%MAXSTP)
!
return
99    format(270a)
end subroutine
!-------------------------------------------------------------------------------

subroutine wcalct (ICALC,lu)

USE MODELDEFINITION, only: mode
USE CONSTANTS,       only: year
use timestepmodule,  only: timestepdatactx, &
                           nintg, ntcalc

implicit none
!-pass
integer          :: lu
integer          :: icalc
dimension        :: ICALC(*)
!-local
integer          :: m, n, i, igrp
double precision :: x

!-init
n    = 0
igrp = 1
i    = 1
x    = 0d0

if (  MODE.lt.3 .or. &
    NTCALC.eq.0 .or. &
     NINTG.le.0) then
    return
endif

m = 1

do  igrp=1,NINTG
    do i=1,timestepdatactx%MAXSTP(igrp)

        n = n + 1
        x = x + timestepdatactx%DELT(igrp)

        if (NTCALC.gt.0 .and. ICALC(m).eq.n) then

            write(lu,2) n,x,x/YEAR,x/YEAR*1D-6
 2               format(1x,i5,1x,1Pg12.3,1x, &
             'sec',3x,1PG12.3,1x,'year',3x,1Pg15.3,1x,'Ma')
            m = m + 1

            if (m.gt.NTCALC) then
                return
            endif

        else if (NTCALC.lt.0 .and. MOD(n,-NTCALC).eq.0) then
            write(lu,2) n,x,x/YEAR,x/YEAR*1D-6
        endif
    enddo
enddo

return
end subroutine

!-------------------------------------------------------------------------------

subroutine lstout(IMPRINT,ITPRINT,MAXSTP)

USE TIMESTEPMODULE,  only: NINTG, NMPRT,IMPRT,NTPRT,ITPRT
USE MODELCTX
use modeldefinition, only: MODE
use iomodule,        only: stderr
use materialsmodule, only: lmat

! Check if there are calculations after last result

implicit none

!-pass
integer   :: IMPRINT,ITPRINT,MAXSTP
dimension :: IMPRINT(*),ITPRINT(*),MAXSTP(*)
!-locl
integer   :: LstMch,LstTmp,LstStp,LstPrt
integer   :: i
!-init

LstMch = 0
LstTmp = 0
LstStp = 0

if (NINTG.eq.0) then
    return
endif

do i=1,NINTG
    LstStp = LstStp + MAXSTP(i)
enddo

if (LstStp.eq.0) then
    write(stderr,10) NINTG
10       format(1x,'Inconsistency in time step data: NINTG = ',I3,1x, &
     'and last step number = 0')
    return
endif

if (MODE.eq.2 .or. MODE.eq.6) then
    if (NMPRT.gt.0) then
        LstMch = IMPRINT(NMPRT)
    else if (NMPRT.lt.0) then
        LstMch = INT(LstStp/NMPRT) * NMPRT
    endif
endif

if (MODE.eq.5 .or. MODE.eq.6) then
    if (NTPRT.gt.0) then
        LstTmp = ITPRINT(NTPRT)
    else if (NTPRT.lt.0) then
        LstTmp = INT(LstStp/NTPRT) * NTPRT
    endif
endif

LstPrt = MAX(LstMch,LstTmp)

if (LstPrt.lt.LstStp) then
    if (MODE.eq.2) then
        write(stderr,20) LstStp,LstMch
20           format(1x,'Too much work: last calculation at time',1x, &
         'step ',I6/1x, &
         '               last mechanical output at     ',I6)
    else if (MODE.eq.5) then
        write(stderr,30) LstStp,LstTmp
30           format(1x,'Too much work: last calculation at time',1x, &
         'step ',I6/1x, &
         '               last thermal output at        ',I6)
    else
        write(stderr,40) LstStp,LstMch,LstTmp
40           format(1x,'Too much work: last calculation at time',1x, &
         'step ',I6/1x, &
         '               last mechanical output at     ', &
         I6/1x,'               last thermal output at        ', &
         I6)
    endif
endif

return
end

!-------------------------------------------------------------------------------
subroutine writbc (SYMBL,MAXSYM,lu)

use iomodule,           only: ascout

implicit none

!-pass
integer          :: lu, maxSym
double precision :: SYMBL
dimension        :: SYMBL(2,MAXSYM)
!-local
integer          :: l, n, i, j

l = 1
100    n = 0
if (l.lt.MAXSYM) then
    n = INT(SYMBL(2,l))
endif

if (n.le.0) then
    return
endif

if (ascout) then
    write(lu,1,err=200) (SYMBL(i,l),i=1,2)
    write(lu,1,err=200) ((SYMBL(i,l+j),i=1,2),j=1,n)
 1        format(1x,2(1PG15.6))
else
    write(lu,err=200) (REAL(SYMBL(i,l)),i=1,2)
    write(lu,err=200) ((REAL(SYMBL(i,l+j)),i=1,2),j=1,n)
endif

l = l + n + 1
goto 100

200    call prerr('writbc')
call exitp(1)

end subroutine

!-------------------------------------------------------------------------------

subroutine writij (IJSUR,NUMEL,NSUR,IJFILE,ierr)

use errors,          only: exec, ok

implicit none

!-pass
integer   :: NUMEL,NSUR

integer   :: ijsur
dimension :: IJSUR(NSUR,*)

character :: IJFILE*(*)



!-locl

integer   :: luij,i,j, ierr

integer, external :: nextlu
!-init
luij = nextlu(0)
ierr = OK

call fbopen (luij,IJFILE,'unknown')

write(luij,err=100) &
      ((IJSUR(i,j),i=1,NSUR),j=1,NUMEL)

call closef(luij)

50 return
100  call prerr('writij')
ierr = EXEC
goto 50
end subroutine

!-------------------------------------------------------------------------------

subroutine writria (XPT,NUMTR,lu,callerID,ierr)

use meshdatamodule,  only: plotDOFperElem, &
                           NPlotEntriesPerElem
use debugmodule,     only: iecho, &
                           parallellog
use iomodule,        only: ascout, stderr
use postProcessData, only: GMT,connout
use errors,          only: exec, ok

! writes triangle coordinates that are contained in XPT

implicit none

!-pass
integer          :: NPLT,NUMTR,lu, callerID, ierr

double precision :: XPT(plotDOFperElem,NPlotEntriesPerElem,NUMTR)


!-local
integer          :: itype
!-init

if (iecho.eq.8) then
    call ParallelLog("writria", "Entered")
endif

ierr = OK

if (NUMTR.le.0) then
    write(stderr,2)
 2        format(1x,'writria: no data')
    ierr = EXEC
    return
endif

!    NPLDOF is number of degrees of freedom for plotting a triangle
!    so for a value that is constant per triangle, it is 1
!    and for a value that is linear, it is 3.

if (plotDOFperElem.eq.3) then
    itype = 7
! itype determines on whether the file writing routines in tools
! (wGMTfile  and  wPfile) write a
! 'PSXY multisegment file' (itype == 4) which is ???
! or 
! 'contour file' (itype == 7)
else
    write(stderr,*) &
      'In pltinp; writria; plotDOFperElem should be 3, but is ', &
       plotDOFperElem, '. Exiting...'
    call exitp(1)
endif


call wGMTfile (XPT,NPlotEntriesPerElem*NUMTR, &
                   itype,lu,ascout,stderr,ierr)

if (ierr.ne.0) then
    ierr = EXEC
endif

if (iecho.eq.8) then
    call ParallelLog("writria", "Done")
endif


return
end
