subroutine comat (MAT,IJSUR,XPT,IEN,STN,PRPMAT,NSUR,NPLDOF,NPLT, &
                 NEN,NSTR,CST)

USE ALGEBRA
use materialsmodule, only: MaterialIDOfElement
use meshdatamodule,  only: meshdatactx
use constants,       only: half
use plotControl,     only: Vscale, viscli, vislog

! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element

implicit none

integer, parameter :: NTRIA=2
!-pass
logical            :: CST
integer            :: NSUR,NPLDOF,NPLT,NEN,NSTR

integer            :: mat, ijsur, itria, ien
double precision   :: xpt, stn, prpmat
dimension          :: MAT(*),IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      ITRIA(NTRIA,2),IEN(NEN,*),STN(NSTR,*),PRPMAT(6,*)
!-locl
logical            :: QUAD
integer            :: numtr
integer            :: n, nn, m, mm, ntr, i
double precision   :: emhu, anpwr, eff

double precision, external :: EFFSTN, etaLin

NUMTR = 0

do nn=1,meshdatactx%neglobal
    n = nn
    m = MaterialIDOfElement(n)

    QUAD=(IEN(3,n).ne.IEN(4,n))

    if (CST) then

        emhu  = PRPMAT(3,m)
        anpwr = PRPMAT(4,m)

        if (viscli) then
            eff = EFFSTN (STN(1,n),HALF,.false.)
            emhu = etaLin (eff,emhu,anpwr)
        endif

        if (Vscale.gt.0.) then
            emhu = Vscale*emhu
        endif

        if (vislog) then
            emhu = TENLOG(emhu)
        endif

        NUMTR = NUMTR + 1
        XPT(3,1,NUMTR) = emhu
        XPT(3,2,NUMTR) = emhu
        XPT(3,3,NUMTR) = emhu
        if (QUAD) then
            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = emhu
            XPT(3,2,NUMTR) = emhu
            XPT(3,3,NUMTR) = emhu
        endif

    else ! CST 
!       determine triangle element indices
        call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!       number of local triangles NTR
        if (ntr.gt.0) then

            emhu  = PRPMAT(3,m)
            anpwr = PRPMAT(4,m)

            if (viscli) then
                eff = EFFSTN (STN(1,n),HALF,.false.)
                emhu  = etaLin (eff,emhu,anpwr)
            endif

            if (Vscale.gt.0.) then
                emhu = Vscale*emhu
            endif

            if (vislog) then
                emhu = TENLOG(emhu)
            endif

            do mm=1,ntr
                XPT(3,1,NUMTR+mm) = emhu
            enddo

!           loop triangles
            do mm=1,ntr
!               loop surrounding elements
                do i=1,2
                    n  = ITRIA(i,mm)
                    m = MaterialIDOfElement(n)

                    emhu  = PRPMAT(3,m)
                    anpwr = PRPMAT(4,m)

                    if (viscli) then
                        eff  = EFFSTN (STN(1,n),HALF,.false.)
                        emhu = etaLin (eff,emhu,anpwr)
                    endif

                    if (Vscale.gt.0.) then
                        emhu = Vscale*emhu
                    endif

                    if (vislog) then
                        emhu = TENLOG(emhu)
                    endif

                    XPT(3,i+1,NUMTR+mm) = emhu

                enddo
            enddo
            NUMTR = NUMTR + ntr
        endif
    endif
enddo
!
return
end subroutine

!-------------------------------------------------------------------------------

subroutine cmat (IJSUR,XPT,PRP,LDP,IPRP,NSUR,NPLDOF, &
   NPLT,NEN,NSTR,CST)

use materialsmodule, only: matdatactx, &
                   MaterialIDOfElement
use meshdatamodule,  only: meshdatactx


!    USE MODELDEFINITION

! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element

! PRP is from the materialdatamodule and can be either 
! matdatactx%PRP[MAT,CRP,DIF,LAS,TEM]

! and the flag, scaled from 1 to max withing the block


! construction during subroutine xma
! 0 MAT           > PRP = matdatactx%PRPMAT   iprp 0
! 1-6 PRPMAT      > PRP = matdatactx%PRPMAT   iprp 1-6
! 7-12 PRPCRP     > PRP = matdatactx%PRPCRP   iprp 1-6
! 13-21 PRPLAS    > PRP = matdatactx%PRPLAS   iprp 1-9
! 22-27 PRPTEM    > PRP = matdatactx%PRPTEM   iprp 1-6


! IPRP is the material property requested in this specific group


implicit none

integer, parameter :: NTRIA=2
!-pass
logical            :: CST
integer            :: NSUR,NPLDOF,NPLT,NEN,NSTR,LDP,IPRP

integer            :: ijsur, itria
double precision   :: xpt, prp
dimension          :: IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      ITRIA(NTRIA,2),PRP(LDP,*)
!-local
logical            :: QUAD
integer            :: numtr
integer            :: n, nn, m, mm, ntr, i
double precision   :: emhu, anpwr, qty

numtr = 0

do nn=1,meshdatactx%neglobal
    n = nn
    ! determine material number of element 
    m = MaterialIDOfElement(n)

    if (IPRP.le.0) then
        ! give material number
        qty = FLOAT(m)
    else
        ! get a value from the properties array
        ! corresponding to the material number
        qty = PRP(IPRP,m)
    endif

    quad = (meshdatactx%IEN(3,n).ne.meshdatactx%IEN(4,n))

    if (CST) then
        numtr = numtr + 1
        XPT(3,1,numtr) = qty
        XPT(3,2,numtr) = qty
        XPT(3,3,numtr) = qty
        if (quad) then
            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = qty
            XPT(3,2,NUMTR) = qty
            XPT(3,3,NUMTR) = qty
        endif
    else
!       determine triangle element indices
        call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,quad)
!       number of local triangles NTR
        if (ntr.gt.0) then
            do mm=1,ntr
                XPT(3,1,NUMTR+mm) = qty
            enddo
!           loop triangles
            do mm=1,ntr
!               loop surrounding elements
                do i=1,2
                    n  = ITRIA(i,mm)
                    m = MaterialIDOfElement(n)

                    if (IPRP.le.0) then
                        qty = FLOAT(m)
                    else
                        qty = PRP(IPRP,m)
                    endif

                    XPT(3,i+1,NUMTR+mm) = qty
                enddo
            enddo
            NUMTR = NUMTR + ntr
        endif
    endif
enddo
!
return
end subroutine

!-------------------------------------------------------------------------------


subroutine triangles (IJSUR,NSUR,ITRIA,NTRIA,NEL,nt,QUAD)

use iomodule, only: stderr

! finds element indices of plot triangles for element NEL
! for rectangular elements triangles are created from two indices
! greater than NEL or two indices less than NEL. Triangular elements
! generate overlapping triangles with all bordering elements.

implicit none

!-pass
integer      :: NSUR, NTRIA, NEL, nt

integer      :: IJSUR,        ITRIA
dimension    :: IJSUR(NSUR,*),ITRIA(NTRIA,*)

integer      :: IPLUS(2),IMIN(2)
logical      :: QUAD

!-local
integer      :: i, nmin, nplus, ielm

i = 1

if (QUAD) then

    nplus = 0
    nmin  = 0

100 if (IJSUR(i,NEL).ne.0) then

        ielm = IJSUR(i,NEL)

        if (ielm.gt.NEL) then
            nplus = nplus + 1
           if (nplus.gt.2) then
                write(stderr,1)
 1                 format(1x,'triangles: nplus error')
                call exitp(1)
            endif
            IPLUS(nplus) = ielm
        else if (ielm.lt.NEL) then
            nmin = nmin + 1
            if (nmin.gt.2) then
                write(stderr,2)
 2                format(1x,'triangles: nmin error')
                call exitp(1)
            endif
            IMIN(nmin) = ielm
        endif

        i = i + 1

        if (i.le.NSUR) then
            goto 100
        endif

    endif

    nt = 0

    if (nplus.eq.2) then
        nt = nt + 1
        ITRIA(1,nt) = IPLUS(1)
        ITRIA(2,nt) = IPLUS(2)
    endif

    if (nmin.eq.2) then
        nt = nt + 1
        ITRIA(1,nt) = IMIN(1)
        ITRIA(2,nt) = IMIN(2)
    endif

else
!        count number of surrounding elms
200 if (IJSUR(i,NEL).ne.0) then
        i = i + 1
        if (i.le.NSUR) then
            goto 200
        endif
    endif

    i  = i - 1
    nt = 0

    if (i.ge.2) then
        nt = nt + 1
        ITRIA(1,nt) = IJSUR(1,NEL)
        ITRIA(2,nt) = IJSUR(2,NEL)
        if (i.ge.3) then
            nt = nt + 1
            ITRIA(1,nt) = IJSUR(1,NEL)
            ITRIA(2,nt) = IJSUR(3,NEL)
            nt = nt + 1
            ITRIA(1,nt) = IJSUR(2,NEL)
            ITRIA(2,nt) = IJSUR(3,NEL)
        endif
    endif
endif

return
end subroutine

!-------------------------------------------------------------------------------
subroutine tricrd (IEN,X,D,IJSUR,XPT,LMF,TFAULT,DXE, &
                   SKEW,DEFORM,DEFFAC,NUMTR,CST)

! subroutine sets element-based coordinates XPT and NUMTR

USE MODELDEFINITION, only:  RADIUS,NUMEL
USE CONSTANTS,       only: zero, deg2rad, useDXE
use debugmodule,     only: iecho
use modeltopology,   only: ndof, nen, nsd
use iomodule,        only: stderr
use plotControl,     only: jShell

!
! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element
!

implicit none

integer, parameter :: NPLT=3, NPLDOF=3, NSUR=4
!-pass
integer            :: ien, ijsur, lmf, numtr
double precision   :: x, d, xpt, skew, tfault, dxe
dimension          :: IEN(NEN,*),X(NSD,*),D(NDOF,*),IJSUR(NSUR,*), &
                      XPT(NPLDOF,NPLT,*),SKEW(2,*),LMF(NDOF,NEN,*), &
                      TFAULT(NDOF,*), DXE(NDOF,NEN,*)

logical            :: DEFORM,CST
double precision   :: deffac
!-locl
integer, parameter :: NTRIA=2
integer            :: ITRIA(NTRIA,2)
double precision   :: CENTR(NSD),XL(NSD,NEN),DL(NDOF,NEN)
logical            :: QUAD
double precision   :: def, ylat

integer            :: i, j, mm, n, nn, ntr, ierr

!-init

def  = ZERO

if (DEFORM) then
    def = deffac
endif

NUMTR = 0

do nn=1,NUMEL

    n = nn
    QUAD = (IEN(3,n).ne.IEN(4,n))

    if (.not.CST) then
!       determine triangle element indices

        call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)

        if (iecho.eq.5) then
            write(stderr,1) n,ntr
        endif

 1      format(1x,'element ',I12,': generates ',I1, &
         ' triangles')

!       number of local triangles NTR
        if (ntr.gt.0) then
            if (iecho.eq.5) then
                write(stderr,2) (j,n,(ITRIA(i,j),i=1,NTRIA), &
                 j=1,ntr)
 2              format(1x,'triangle ',I1,': elms ',I12,1X,I12,1X,I12)
            endif

!           calculate element center coordinates
            call ELCTR (CENTR,IEN(1,n),X,D,n,TFAULT, &
             SKEW, def,ierr)

            if (ierr.ne.0) then
                return
            endif

            do mm=1,ntr
               XPT(1,1,NUMTR+mm) = CENTR(1)
               XPT(2,1,NUMTR+mm) = CENTR(2)
            enddo

!           loop triangles
            do mm=1,ntr
!               loop surrounding elements
                do i=1,2
                    n = ITRIA(i,mm)
                    call ELCTR (CENTR,IEN(1,n),X,D,n, &
                    TFAULT, SKEW,def,ierr)
                    if (ierr.eq.0) then
                        XPT(1,i+1,NUMTR+mm) = CENTR(1)
                        XPT(2,i+1,NUMTR+mm) = CENTR(2)
                    endif
                enddo
            enddo

            NUMTR = NUMTR + ntr
        endif
    else

        call LCOORD (X,XL,IEN(1,n))

        if (DEFORM) then
            call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,165)
            call ADDSNE (DL,n,NDOF,NEN,useDXE)
            if (JSHELL.ne.0) then
                call SPHDISP (XL,DL)
            endif
            call REZONE(XL,DL,def)
        endif

        NUMTR = NUMTR + 1
        XPT(1,1,NUMTR) = XL(1,1)
        XPT(2,1,NUMTR) = XL(2,1)
        XPT(1,2,NUMTR) = XL(1,2)
        XPT(2,2,NUMTR) = XL(2,2)
        XPT(1,3,NUMTR) = XL(1,3)
        XPT(2,3,NUMTR) = XL(2,3)
        if (QUAD) then
            NUMTR = NUMTR + 1
            XPT(1,1,NUMTR) = XL(1,1)
            XPT(2,1,NUMTR) = XL(2,1)
            XPT(1,2,NUMTR) = XL(1,3)
            XPT(2,2,NUMTR) = XL(2,3)
            XPT(1,3,NUMTR) = XL(1,4)
            XPT(2,3,NUMTR) = XL(2,4)
        endif

    endif
enddo

return
end subroutine

!-------------------------------------------------------------------------------

subroutine SPHDISP (XL,DL)

! spherical shell routine: converts element-based displacements to degrees in
! easting and northing direction

USE MODELDEFINITION, only: RADIUS
USE CONSTANTS,       only: DEG2RAD
USE MODELTOPOLOGY,   only: NDOF, NEN, NSD

implicit none

integer          :: i
double precision :: XL(NSD,NEN),DL(NDOF,NEN),ylat

do i=1,NEN
    ylat = XL(2,i) * DEG2RAD
    DL(1,i) = DL(1,i) / (DEG2RAD*RADIUS*COS(ylat))
    DL(2,i) = DL(2,i) / (DEG2RAD*RADIUS)
enddo

return
end
!-------------------------------------------------------------------------------
subroutine INFROT (IJSUR,XPT,IEN,X,D,LMF,TFAULT, &
                   SKEW,verbose,luMon,CST)

USE MODELDEFINITION, only: lgdef
use meshdatamodule,  only: meshdatactx
use constants,       only: useDXE, half
use modeltopology,   only: ndof, nsd, nen
use monitorVars,     only: elmon, nnel, monel
use plotControl,     only: jShell

!
! displacement gradients or infinitesimal rotation

implicit none


integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3, NPLDOF=3, NSUR=4


!-pass
integer            :: luMon
logical            :: CST,verbose

integer            :: ijsur, ien, lmf
double precision   :: xpt, x, d, tfault, skew
dimension          :: IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      IEN(NEN,*),X(NSD,*),D(NDOF,*), &
                      LMF(NDOF,NEN,*), TFAULT(NDOF,*),SKEW(2,*)
!-locl
logical            :: QUAD,rotation
integer            :: ITRIA(NTRIA,2)
double precision   :: BP(NDOF,NSD),XL(NSD,NEN),DL(NDOF,NEN)

integer            :: n, nn, mm, numtr, i, j, ierr, ntr
double precision   :: rxy, val

!-init
rotation = .false.
NUMTR = 0
!
if (elmon) then
    do nn=1,nnel
        n = monel(nn)
        call LCOORD (X,XL,IEN(1,n))
        call LDISP (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,166)
        call ADDSNE (DL,n,NDOF,NEN,0)

        QUAD=(IEN(3,n).ne.IEN(4,n))

        if (JSHELL.ne.0) then
            call SHLCRD (XL,QUAD)
        endif

        call BPMATRIX (XL,DL,QUAD,BP,ierr)

        if (rotation) then
            rxy = HALF * (BP(2,1)-BP(1,2))

            if (verbose) then
                write(luMon,1) n,rxy
 1                   format(1x,'Wxy in element ',I12,' = ', &
                 1PG14.6)
            else
                write(luMon,*) rxy
            endif

        else

            if (verbose) then
                write(luMon,2) n,BP(1,1),BP(1,2),BP(2,1),BP(2,2)
 2                   format(1x,'element ',I12, &
                 ': dUx/dx, dUx/dy, dUy/dx, dUy,dy = ', &
                 4(1X,1PG16.6))
            else
                write(luMon,*) BP(1,1),BP(1,2),BP(2,1),BP(2,2)
            endif
        endif
    enddo
else
    do nn=1,meshdatactx%neglobal
        n = nn
        QUAD=(IEN(3,n).ne.IEN(4,n))
        if (CST) then

            call LCOORD (X,XL,IEN(1,n))
            call LDISP (DL,D,IEN(1,n),NDOF,NEN)
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,167)
            call ADDSNE (DL,n,NDOF,NEN,0)

            QUAD=(IEN(3,n).ne.IEN(4,n))

            if (JSHELL.ne.0) then
                call SHLCRD (XL,QUAD)
            endif

            call BPMATRIX (XL,DL,QUAD,BP,ierr)

            val = HALF * (BP(2,1)-BP(1,2))

            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = val
            XPT(3,2,NUMTR) = val
            XPT(3,3,NUMTR) = val
            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = val
                XPT(3,2,NUMTR) = val
                XPT(3,3,NUMTR) = val
            endif
        else
!                determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!                number of local triangles NTR
            if (ntr.gt.0) then

                call LCOORD (X,XL,IEN(1,n))
                call LDISP (DL,D,IEN(1,n),NDOF,NEN)
                call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,168)
                call ADDSNE (DL,n,NDOF,NEN,0)

                QUAD=(IEN(3,n).ne.IEN(4,n))

                if (JSHELL.ne.0) then
                    call SHLCRD (XL,QUAD)
                endif

                call BPMATRIX (XL,DL,QUAD,BP,ierr)

                val = HALF * (BP(2,1)-BP(1,2))

                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = val
                enddo
!                loop triangles
                do mm=1,ntr
!                    loop surrounding elements
                    do i=1,2

                        n  = ITRIA(i,mm)

                        call LCOORD (X,XL,IEN(1,n))
                        call LDISP (DL,D,IEN(1,n),NDOF,NEN)
                        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,169)
                        call ADDSNE (DL,n,NDOF,NEN,0)

                        QUAD=(IEN(3,n).ne.IEN(4,n))

                        if (JSHELL.ne.0) then
                            call SHLCRD (XL,QUAD)
                        endif

                        call BPMATRIX (XL,DL,QUAD,BP,ierr)

                        val = HALF * (BP(2,1)-BP(1,2))
                        XPT(3,i+1,NUMTR+mm) = val
                    enddo
                enddo
                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif

return
end subroutine

!-------------------------------------------------------------------------------

subroutine INFROTR (IJSUR,XPT,IEN,X,D,DELD,LMF,TFAULT, &
                    SKEW,verbose,luMon,CST)


USE MODELDEFINITION, only: lgdef
use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use constants,       only: useDXE, useDELDXE, half
use modeltopology,   only: ndof, nsd, nen
use monitorVars,     only: elmon, nnel, monel
use plotControl,     only: jShell

!
! velocity gradients or infinitesimal rotation rate
!
implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3, NPLDOF=3, NSUR=4

!-pass
integer          :: luMon
logical          :: CST,verbose

integer          :: ijsur, ien, lmf
double precision :: xpt, x, d, tfault, skew, deld
dimension        :: IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                    IEN(NEN,*),X(NSD,*),D(NDOF,meshdatactx%nvglobal),&
                    LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
                    SKEW(2,*),DELD(NDOF,meshdatactx%nvglobal)
!-locl
logical          :: QUAD
double precision :: BP(NDOF,NSD),XL(NSD,NEN),DL(NDOF,NEN)
integer          :: ITRIA(NTRIA,2)
integer          :: i, mm, n, nn, ntr, numtr
double precision :: rxy, val
integer          :: ierr

!-init
NUMTR = 0

if (elmon) then
    do nn=1,nnel
        n = monel(nn)

        if (n.le.0) then
            stop 'infrotr node < 0 error'
        endif

        call LCOORD (X,XL,IEN(1,n))
        QUAD=(IEN(3,n).ne.IEN(4,n))

        if (JSHELL.ne.0) then
            call SHLCRD (XL,QUAD)
        endif

        call LDISP (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,168)
        call ADDSNE (DL,n,NDOF,NEN,useDXE)
        if (lgdef.ne.0) call REZONE(XL,DL,1d0)

!! missing here are faulted node velocity contributions !!
        call LDISP (DL,DELD,IEN(1,n),NDOF,NEN)
        call ADDSNE (DL,n,NDOF,NEN,useDELDXE)
        call BPMATRIX (XL,DL,QUAD,BP,ierr)

        if (ierr.ne.0) then
            cycle
        endif

        rxy = HALF * (BP(2,1)-BP(1,2))

        if (verbose) then
            write(luMon,1) n,rxy
1           format(1x,'Wxy in element ',I12,' = ', 1PG14.6)
        else
            write(luMon,*) rxy

        endif
    enddo
else
    ! no monitoring, do all the elements
    do nn=1,meshdatactx%neglobal
        n = nn

        if (CST) then
            ! CST = constant stress mode, for interpolation.
            ! your guess is as good as mine...

            call LCOORD (X,XL,IEN(1,n))
            QUAD=(IEN(3,n).ne.IEN(4,n))

            if (JSHELL.ne.0) then
                call SHLCRD (XL,QUAD)
            endif

! first determine local coordinates: XL
            call LDISP (DL,D,IEN(1,n),NDOF,NEN)
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,168)
            call ADDSNE (DL,n,NDOF,NEN,useDXE)
            if (lgdef.ne.0) call REZONE(XL,DL,1d0)



!! missing here are faulted node velocity contributions !!

! after that, determine the velocity in that frame
! overwrite DL, that was the local dispalcement, but will now be the local velocity
            call LDISP (DL,DELD,IEN(1,n),NDOF,NEN)

! include differential velocities
            call ADDSNE (DL,n,NDOF,NEN,useDELDXE)

            call BPMATRIX (XL,DL,QUAD,BP,ierr)

            val = HALF * (BP(2,1)-BP(1,2))

            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = val
            XPT(3,2,NUMTR) = val
            XPT(3,3,NUMTR) = val

            if (QUAD) then
                NUMTR = NUMTR + 1
                XPT(3,1,NUMTR) = val
                XPT(3,2,NUMTR) = val
                XPT(3,3,NUMTR) = val
            endif

        else
!           determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!           number of local triangles NTR
            if (ntr.gt.0) then
                call LCOORD (X,XL,IEN(1,n))
                QUAD=(IEN(3,n).ne.IEN(4,n))

                if (JSHELL.ne.0) then
                    call SHLCRD (XL,QUAD)
                endif

                call LDISP (DL,D,IEN(1,n),NDOF,NEN)
                call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,168)
                call ADDSNE (DL,n,NDOF,NEN,useDXE)
                if (lgdef.ne.0) call REZONE(XL,DL,1d0)

!! todo missing here are faulted node velocity contributions !!
                call LDISP (DL,DELD,IEN(1,n),NDOF,NEN)
                call ADDSNE (DL,n,NDOF,NEN,useDELDXE)
                call BPMATRIX (XL,DL,QUAD,BP,ierr)

                val = HALF * (BP(2,1)-BP(1,2))
                do mm=1,ntr
                    XPT(3,1,NUMTR+mm) = val
                enddo

!               loop triangles
                do mm=1,ntr
!                   loop surrounding elements
                    do i=1,2
                        n  = ITRIA(i,mm)
                        call LCOORD (X,XL,IEN(1,n))
                        QUAD=(IEN(3,n).ne.IEN(4,n))

                        if (JSHELL.ne.0) then
                            call SHLCRD (XL,QUAD)
                        endif

                        call LDISP (DL,D,IEN(1,n),NDOF,NEN)
                        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,168)
                        call ADDSNE (DL,n,NDOF,NEN,useDXE)
                        if (lgdef.ne.0) call REZONE(XL,DL,1d0)

!! todo missing here are faulted node velocity contributions !!
                        call LDISP (DL,DELD,IEN(1,n),NDOF,NEN)
                        call ADDSNE (DL,n,NDOF,NEN,useDELDXE)
                        call BPMATRIX (XL,DL,QUAD,BP,ierr)

                        val = HALF * (BP(2,1)-BP(1,2))
                        XPT(3,i+1,NUMTR+mm) = val
                    enddo
                enddo
                NUMTR = NUMTR + ntr
            endif
        endif
    enddo
endif

return
end subroutine infrotr
