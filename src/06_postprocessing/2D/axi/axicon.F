subroutine comat (MAT,IJSUR,XPT,IEN,STN,PRPMAT,NSUR,NPLDOF,NPLT, &
                  NEN,NSTR,CST)

USE ALGEBRA
use materialsmodule, only: MaterialIDOfElement
USE MODELDEFINITION, only: numel
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
logical QUAD
integer            :: numtr, n, nn, m, mm, ntr, i
double precision   :: emhu, anpwr, eff

double precision, external :: EFFSTN
double precision, external :: etaLin

NUMTR = 0
do nn=1,NUMEL

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

    else
!            determine triangle element indices
        call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!            number of local triangles NTR
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
!            loop triangles
            do mm=1,ntr
!                loop surrounding elements
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

return
end subroutine
!-------------------------------------------------------------------------------
subroutine cmat (MAT,IJSUR,XPT,IEN,PRP,LDP,IPRP,NSUR,NPLDOF, &
                 NPLT,NEN,NSTR,CST)

use meshdatamodule,  only: meshdatactx
use materialsmodule, only: MaterialIDOfElement

! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element

implicit none

integer, parameter :: NTRIA=2

!-pass
logical            :: CST
integer            :: NSUR,NPLDOF,NPLT,NEN,NSTR,IPRP,LDP

integer            :: mat, ijsur, itria, ien
double precision   :: xpt, prp
dimension          :: MAT(*),IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      ITRIA(NTRIA,2),IEN(NEN,*),PRP(LDP,*)
!-locl
integer            :: numtr,nn,n,m, mm, i, ntr
logical            :: quad
double precision   :: qty
!
numtr = 0

do nn=1,meshdatactx%neglobal
 
    n = nn
    m = MaterialIDOfElement(n)

    if (IPRP.le.0) then
        qty = FLOAT(m)
    else
        qty = PRP(IPRP,m)
    endif

    quad = (IEN(3,n).ne.IEN(4,n))

    if (CST) then
        numtr = numtr + 1
        XPT(3,1,numtr) = qty
        XPT(3,2,numtr) = qty
        XPT(3,3,numtr) = qty
        if (IEN(3,n).ne.IEN(4,n)) then
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

! This sub computes NUMTR, which is the number of triangles
! that is plotted in a contour plot, which, i essence, comes down
! to num of mesh triangles + 2 * num of quads... probably


USE MODELDEFINITION, only: RADIUS,NUMEL
USE CONSTANTS,       only: zero, deg2rad
use debugmodule,     only: iecho
use iomodule,        only: stderr
use plotControl,     only: jShell

!
! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element
!

implicit none

integer NPLDOF,NPLT,NSUR,NEN,NSD,NDOF
parameter (NPLT=3, NPLDOF=3, NDOF=2, NSD=2, NEN=4, NSUR=4)
!-pass
integer          :: numtr
double precision :: deffac

integer          :: ien, ijsur, lmf
double precision :: x, d, xpt, skew, tfault, dxe

dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),IJSUR(NSUR,*), &
   XPT(NPLDOF,NPLT,*),SKEW(2,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   DXE(NDOF,NEN,*)

logical DEFORM,CST
!-locl
integer, parameter :: NTRIA=2
integer            :: ITRIA(NTRIA,2)
double precision   :: CENTR(NSD),XL(NSD,NEN),DL(NDOF,NEN)

logical            :: QUAD

integer            :: n, nn, m, mm, ntr, i, j, ierr
double precision   :: def, ylat

!-init
def  = ZERO
if (DEFORM) then
    def = deffac
endif

NUMTR = 0
!
!    call clear (XPT,NPLDOF*NPLT*NSUR*NUMEL)
!
do nn=1,NUMEL

    n = nn
    QUAD = (IEN(3,n).ne.IEN(4,n))

    if (.not.CST) then

!       determine triangle element indices
        call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)


        if (iecho.eq.5) then
            write(stderr,1) n,ntr
 1               format(1x,'element ',I12,': generates ',I1, &
             ' triangles')
        endif

!       number of local triangles NTR
        if (ntr.gt.0) then
            if (iecho.eq.5) then
                write(stderr,2) (j,n,(ITRIA(i,j),i=1,NTRIA), &
                   j=1,ntr)
 2                   format(1x,'triangle ',I1,': elms ',I12,1X,I12,1X,I12)
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
                    TFAULT, &
                    SKEW,def,ierr)
                    if (ierr.eq.0) then
                        XPT(1,i+1,NUMTR+mm) = CENTR(1)
                        XPT(2,i+1,NUMTR+mm) = CENTR(2)
                    endif
                enddo
            enddo

            NUMTR = NUMTR + ntr
        endif

    else  ! CST (????) is .true.

        call LCOORD (X,XL,IEN(1,n))

        if (DEFORM) then

            call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,136)
            call ADDSNE (DL,n,NDOF,NEN,0)
 
            do j=1,NEN
                if (JSHELL.eq.0) then
                    XL(1,j) = XL(1,j) + def*dl(1,j)
                    XL(2,j) = XL(2,j) + def*dl(2,j)
                else
                    ylat    = XL(2,j) * DEG2RAD
                    XL(1,j) = XL(1,j) + def*dl(1,j) / (DEG2RAD*RADIUS*COS(ylat))
                    XL(2,j) = XL(2,j) + def*dl(2,j) / (DEG2RAD*RADIUS)
                endif
            enddo

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
subroutine INFROT (IJSUR,XPT,IEN,X,D,LMF,TFAULT,DXE, &
                   SKEW,verbose,luMon,CST)

USE MODELDEFINITION, only: numel
use constants,       only: half, deg2rad
use monitorVars,     only: elmon, nnel, monel
use plotControl,     only: jShell

!
! displacement gradients or infinitesimal rotation
!
implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3, NPLDOF=3, NSUR=4
integer, parameter :: NDOF=2, NSD=2, NEN=4
!-pass
integer            :: luMon
logical            :: CST,verbose

integer            :: ijsur, ien, lmf
double precision   :: xpt, x, d, tfault, dxe, skew
dimension          :: IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*),&
                      TFAULT(NDOF,*),DXE(NDOF,NEN,*),SKEW(2,*)
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
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,137)
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
    do nn=1,NUMEL
        n = nn
        QUAD=(IEN(3,n).ne.IEN(4,n))
        if (CST) then

            call LCOORD (X,XL,IEN(1,n))
            call LDISP (DL,D,IEN(1,n),NDOF,NEN)
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,138)
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
                call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,139)
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
                        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,140)
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
subroutine INFROTR (IJSUR,XPT,IEN,X,D,DELD,LMF,TFAULT,DXE, &
                    SKEW,verbose,luMon,CST)

USE MODELDEFINITION, only: numel
use constants,       only: half
use monitorVars,     only: elmon, nnel, monel
use plotControl,     only: jShell
!
! velocty gradients or infinitesimal rotation rate
!
implicit none

integer, parameter :: NTRIA=2

integer, parameter :: NPLT=3, NPLDOF=3, NSUR=4
integer, parameter :: NDOF=2, NSD=2, NEN=4
!-pass
integer            :: luMon
logical            :: CST,verbose

integer            :: ijsur, ien, lmf
double precision   :: xpt, x, d, tfault, dxe, skew, deld
dimension          :: IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*),&
                      TFAULT(NDOF,*), DXE(NDOF,NEN,*),SKEW(2,*),DELD(NDOF,*)
!-locl
integer            :: n, nn, m, mm, ntr, i, j, ierr, numtr
double precision   :: rxy, val
logical            :: QUAD,rotation
integer            :: ITRIA(NTRIA,2)
double precision   :: BP(NDOF,NSD),XL(NSD,NEN),DL(NDOF,NEN)

!-init
rotation = .false.
NUMTR = 0
!
if (elmon) then
    do nn=1,nnel
        n = monel(nn)

        call LCOORD (X,XL,IEN(1,n))
        call LDISP (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,137)
        call ADDSNE (DL,n,NDOF,NEN,0)

        QUAD=(IEN(3,n).ne.IEN(4,n))

        if (JSHELL.ne.0) then
            call SHLCRD (XL,QUAD)
        endif

        call BPMATRIX (XL,DELD,QUAD,BP,ierr)

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
    do nn=1,NUMEL
        n = nn

        QUAD=(IEN(3,n).ne.IEN(4,n))

        if (CST) then

            call LCOORD (X,XL,IEN(1,n))
            call LDISP (DL,D,IEN(1,n),NDOF,NEN)
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,138)
            call ADDSNE (DL,n,NDOF,NEN,0)

            QUAD=(IEN(3,n).ne.IEN(4,n))

            if (JSHELL.ne.0) then
                call SHLCRD (XL,QUAD)
            endif

            call BPMATRIX (XL,DELD,QUAD,BP,ierr)

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
                call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,139)
                call ADDSNE (DL,n,NDOF,NEN,0)
                QUAD=(IEN(3,n).ne.IEN(4,n))
                if (JSHELL.ne.0) then
                    call SHLCRD (XL,QUAD)
                endif
                call BPMATRIX (XL,DELD,QUAD,BP,ierr)
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
                        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,140)
                        call ADDSNE (DL,n,NDOF,NEN,0)
                        QUAD=(IEN(3,n).ne.IEN(4,n))
                        if (JSHELL.ne.0) then
                            call SHLCRD (XL,QUAD)
                        endif
                        call BPMATRIX (XL,DELD,QUAD,BP,ierr)
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

!---------------------------------------------------------------------------

subroutine SHLCRD (XL,QUAD)
!
use modeldefinition, only: radius
use constants,       only: zero, deg2rad

implicit none

integer, parameter :: NEN=4,NSD=2
!-pass
double precision   :: XL
logical            :: QUAD
dimension          :: XL(NSD,*)
!-locl
integer            :: i
double precision   :: xlon,ylat,xc,yc,zc
double precision   :: v(3,NEN),north(3),east(3)
!
!       compute element center coordinate
if (QUAD) then
   xc = (XL(1,1)+XL(1,2)+XL(1,3)+XL(1,4))*2.5D-1
   yc = (XL(2,1)+XL(2,2)+XL(2,3)+XL(2,4))*2.5D-1
else
   xc = (XL(1,1)+XL(1,2)+2d0*XL(1,3))*2.5D-1
   yc = (XL(2,1)+XL(2,2)+2d0*XL(2,3))*2.5D-1
endif

xc = xc * DEG2RAD
yc = yc * DEG2RAD

east(1) = -DSIN(xc)
east(2) = DCOS(xc)
east(3) = ZERO

north(1) = -DSIN(yc)*DCOS(xc)
north(2) = -DSIN(yc)*DSIN(xc)
north(3) = DCOS(yc)

!       Derive local coordinate directions
xlon = XL(1,1)*DEG2RAD
ylat = XL(2,1)*DEG2RAD

v(1,1) = DCOS(ylat)
v(2,1) = v(1,1)*DSIN(xlon)
v(1,1) = v(1,1)*DCOS(xlon)
v(3,1) = DSIN(ylat)

xlon = XL(1,2)*DEG2RAD
ylat = XL(2,2)*DEG2RAD

v(1,2) = DCOS(ylat)
v(2,2) = v(1,2)*DSIN(xlon)
v(1,2) = v(1,2)*DCOS(xlon)
v(3,2) = DSIN(ylat)

xlon = XL(1,3)*DEG2RAD
ylat = XL(2,3)*DEG2RAD

v(1,3) = DCOS(ylat)
v(2,3) = v(1,3)*DSIN(xlon)
v(1,3) = v(1,3)*DCOS(xlon)
v(3,3) = DSIN(ylat)

if (QUAD) then
    xlon = XL(1,4)*DEG2RAD
    ylat = XL(2,4)*DEG2RAD
    v(1,4) = DCOS(ylat)
    v(2,4) = v(1,4)*DSIN(xlon)
    v(1,4) = v(1,4)*DCOS(xlon)
    v(3,4) = DSIN(ylat)
endif

!       Compute node coordinates relative to element center
if (QUAD) then
    xc = (v(1,1)+v(1,2)+v(1,3)+v(1,4))*2.5D-1
    yc = (v(2,1)+v(2,2)+v(2,3)+v(2,4))*2.5D-1
    zc = (v(3,1)+v(3,2)+v(3,3)+v(3,4))*2.5D-1
    do i=1,4
        v(1,i) = (v(1,i)-xc)*RADIUS
        v(2,i) = (v(2,i)-yc)*RADIUS
        v(3,i) = (v(3,i)-zc)*RADIUS
    enddo
else
    xc = (v(1,1)+v(1,2)+2d0*v(1,3))*2.5D-1
    yc = (v(2,1)+v(2,2)+2d0*v(2,3))*2.5D-1
    zc = (v(3,1)+v(3,2)+2d0*v(3,3))*2.5D-1
    do i=1,3
        v(1,i) = (v(1,i)-xc)*RADIUS
        v(2,i) = (v(2,i)-yc)*RADIUS
        v(3,i) = (v(3,i)-zc)*RADIUS
        write(*,'(1X,3F20.6)') v(1,i),v(2,i),v(3,i)
    enddo
endif
!
!       Project node coordinates on local east and north directions
do i=1,3
    XL(1,i) = v(1,i)*east(1)+v(2,i)*east(2)
    XL(2,i) = v(1,i)*north(1)+v(2,i)*north(2)+v(3,i)*north(3)
enddo

if (QUAD) then
    XL(1,4) = v(1,4)*east(1)+v(2,4)*east(2)
    XL(2,4) = v(1,4)*north(1)+v(2,4)*north(2)+v(3,4)*north(3)
else
    XL(1,4) = XL(1,3)
    XL(2,4) = XL(2,3)
endif

return
end
!-------------------------------------------------------------------------------
subroutine SPHDISP(XL,DL)
implicit none
double precision :: XL, DL
! potentially called by contur
return
end
