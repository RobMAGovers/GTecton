subroutine comat (MAT,IJSUR,XPT,IEN,STN,PRPMAT,NSUR,NPLDOF,NPLT, &
                  NEN,NSTR,CST)

USE ALGEBRA,         only: tenlog
use materialsmodule, only: MaterialIDOfElement
USE MODELDEFINITION, only: numel, numtr
use plotControl,     only: Vscale, viscli, vislog

! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element

implicit none

integer, parameter :: NTRIA=2

!-pass
logical            :: CST
integer            :: NSUR,NPLDOF,NPLT,NEN,NSTR

integer            :: mat, ijsur , itria, ien
double precision   :: xpt, stn, prpmat
dimension          :: MAT(*),IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      ITRIA(NTRIA,2),IEN(NEN,*),STN(NSTR,*),PRPMAT(6,*)
!-locl
logical            :: QUAD
integer            :: n, nn, m, mm, ntr, i
double precision   :: pois, eff, emhu, anpwr

double precision, external :: effstn, etaLin

NUMTR = 0

do nn=1,NUMEL

    n = nn
    m = MaterialIDOfElement(n)
    QUAD=(IEN(3,n).ne.IEN(4,n))

    if (CST) then

        EMHU  = PRPMAT(3,m)
        ANPWR = PRPMAT(4,m)

        if (viscli) then
            POIS = PRPMAT(2,m)
            eff = effstn (STN(1,n),POIS,.false.)
            EMHU  = etaLin (eff,EMHU,ANPWR)
        endif

        if (Vscale.gt.0.) then
            EMHU = Vscale*EMHU
        endif

        if (vislog) then
            EMHU = TENLOG(EMHU)
        endif

        NUMTR = NUMTR + 1
        XPT(3,1,NUMTR) = EMHU
        XPT(3,2,NUMTR) = EMHU
        XPT(3,3,NUMTR) = EMHU
        if (QUAD) then
            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = EMHU
            XPT(3,2,NUMTR) = EMHU
            XPT(3,3,NUMTR) = EMHU
        endif

    else
!            determine triangle element indices
        call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!            number of local triangles NTR
        if (ntr.gt.0) then
            EMHU  = PRPMAT(3,m)
            ANPWR = PRPMAT(4,m)
            if (viscli) then
                POIS = PRPMAT(2,m)
                eff = effstn (STN(1,n),POIS,.false.)
                EMHU  = etaLin (eff,EMHU,ANPWR)
            endif


            if (Vscale.gt.0.) then
                EMHU = Vscale*EMHU
            endif

            if (vislog) then
                EMHU = TENLOG(EMHU)
            endif

            do mm=1,ntr
                XPT(3,1,NUMTR+mm) = EMHU
            enddo
!            loop triangles
            do mm=1,ntr
!                loop surrounding elements
                do i=1,2
                    n  = ITRIA(i,mm)
                    m = MaterialIDOfElement(n)
                    EMHU  = PRPMAT(3,m)
                    ANPWR = PRPMAT(4,m)
                    if (viscli) then
                        POIS = PRPMAT(2,m)
                        eff  = effstn (STN(1,n),POIS,.false.)
                        EMHU = etaLin (eff,EMHU,ANPWR)
                    endif

                    if (Vscale.gt.0.) then
                        EMHU = Vscale*EMHU
                    endif

                    if (vislog) then
                        EMHU = TENLOG(EMHU)
                    endif

                    XPT(3,i+1,NUMTR+mm) = EMHU
                enddo
            enddo
            NUMTR = NUMTR + ntr
        endif
    endif
enddo

return
end subroutine

!-------------------------------------------------------------------------------
subroutine cmat (MAT,PRP,IJSUR,XPT,IEN,LDPRP,IPRP,NSUR,NPLDOF, &
                 NPLT,NEN,NSTR,CST)

use materialsmodule, only: MaterialIDOfElement
USE MODELDEFINITION, only: numel, numtr


! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element

implicit none

integer, parameter :: NTRIA=2

!-pass
logical            :: CST
integer            :: NSUR,NPLDOF,NPLT,NEN,NSTR,LDPRP,IPRP


integer            :: mat, ijsur, itria, ien
double precision   :: xpt, prp

dimension          :: MAT(*),IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      ITRIA(NTRIA,2),IEN(NEN,*),PRP(LDPRP,*)
!-locl
logical            :: quad
integer            :: n, nn, m, mm, ntr, i

NUMTR = 0

do nn=1,NUMEL
    n    = nn
    QUAD = (IEN(3,n).ne.IEN(4,n))
    m    = MaterialIDOfElement(n)
    if (CST) then

        NUMTR = NUMTR + 1
        XPT(3,1,NUMTR) = FLOAT(m)
        XPT(3,2,NUMTR) = FLOAT(m)
        XPT(3,3,NUMTR) = FLOAT(m)
        if (QUAD) then
            NUMTR = NUMTR + 1
            XPT(3,1,NUMTR) = FLOAT(m)
            XPT(3,2,NUMTR) = FLOAT(m)
            XPT(3,3,NUMTR) = FLOAT(m)
        endif

    else
!       determine triangle element indices
        call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!       number of local triangles NTR
        if (ntr.gt.0) then
            do mm=1,ntr
                XPT(3,1,NUMTR+mm) = FLOAT(m)
            enddo
!           loop triangles
            do mm=1,ntr
!               loop surrounding elements
                do i=1,2
                    n  = ITRIA(i,mm)
                    m = MaterialIDOfElement(n)
                    XPT(3,i+1,NUMTR+mm) = FLOAT(m)
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
integer            :: ijsur, itria, nsur, ntria
double precision   :: xpt, prp

dimension          :: IJSUR(NSUR,*), ITRIA(NTRIA,*)
integer            :: IPLUS(2), IMIN(2)
integer            :: nel, nt
logical            :: QUAD

!-locl
integer            :: n, nn, m, mm, ntr, i, ielm
integer            :: nmin, nplus

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
 1                format(1x,'triangles: nplus error')
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
!   count number of surrounding elms
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

!    USE MODELDEFINITION
use debugmodule,     only: iecho
USE MODELDEFINITION, only: numel
use iomodule,        only: stderr
use constants,       only: zero

! XPT (x-y-f, number of nodes on a triangle, number of triangles)
! NTRIA number of triangles that can be generated from one element
!
implicit none

integer, parameter :: NPLT=3, NPLDOF=3, NDOF=3, NSD=2, NEN=4, NSUR=4
!-pass
logical            :: DEFORM,CST
integer            :: numtr
double precision   :: deffac

integer            :: ien, ijsur, lmf
double precision   :: x, d, xpt, skew, tfault, dxe
dimension          :: IEN(NEN,*),X(NSD,*),D(NDOF,*),IJSUR(NSUR,*), &
                      XPT(NPLDOF,NPLT,*),SKEW(2,*),LMF(NDOF,NEN,*), &
                      TFAULT(NDOF,*), DXE(NDOF,NEN,*)

!-locl
integer, parameter :: NTRIA=2
double precision   :: CENTR(3),XL(2,4),DL(3,4)
integer            :: ITRIA(NTRIA,2)
logical            :: QUAD
double precision   :: def
integer            :: n, nn, mm, ntr, i, j, ierr

!-init
def  = ZERO
if (DEFORM) then
    def = DEFFAC
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
 1               format(1x,'element ',I5,': generates ',I1, &
             ' triangles')
        endif

!       number of local triangles NTR
        if (ntr.gt.0) then

            if (iecho.eq.5) then
                write(stderr,2) (j,n,(ITRIA(i,j),i=1,NTRIA), &
                               j=1,ntr)
 2                   format(1x,'triangle ',I1,': elms ',I5,1X,I5,1X, &
                 I5)
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
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,123)
            call ADDSNE (DL,n,NDOF,NEN,0)
            call REZONE (XL,DL,def)
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

subroutine INFROT (IJSUR,XPT,IEN,X,D,LMF,TFAULT, &
                   SKEW,verbose,luMon,CST)

USE MODELDEFINITION, only: numel, numtr
use constants,       only: half
use monitorVars,     only: elmon, nnel, monel

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

integer            :: ien, ijsur, lmf
double precision   :: x, d, xpt, skew, tfault, dxe
dimension          :: IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*), &
                      TFAULT(NDOF,*), SKEW(2,*)


!-locl
logical            :: QUAD,rotation
double precision   :: BP(NDOF,NSD),XL(NSD,NEN),DL(NDOF,NEN)
integer            :: ITRIA(NTRIA,2)
integer            :: n, nn, m, mm, ntr, i, ierr
double precision   :: val, rxy

!-init
rotation = .false.
NUMTR = 0
!
if (elmon) then
    do nn=1,nnel

        n = monel(nn)
 
        call LCOORD (X,XL,IEN(1,n))
        call LDISP (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,124)
        call ADDSNE (DL,n,NDOF,NEN,0)
        QUAD=(IEN(3,n).ne.IEN(4,n))
        call BPMATRIX (XL,DL,QUAD,BP,ierr)

        if (rotation) then
            rxy = HALF * (BP(2,1)-BP(1,2))
            if (verbose) then
                write(luMon,1) n,rxy
 1                   format(1x,'Wxy in element ',I5,' = ', &
                 1PG14.6)
            else
                write(luMon,*) rxy
            endif
        else
            if (verbose) then
                write(luMon,2) n,BP(1,1),BP(1,2),BP(2,1),BP(2,2)
 2                   format(1x,'element ',I5, &
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
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,125)
            call ADDSNE (DL,n,NDOF,NEN,0)
            QUAD=(IEN(3,n).ne.IEN(4,n))
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
                call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,126)
                call ADDSNE (DL,n,NDOF,NEN,0)
                QUAD=(IEN(3,n).ne.IEN(4,n))
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
                        call AddFaultDisplacement (DL,n,TFAULT,NDOF, &
                                                   NEN,127)
                        call ADDSNE (DL,n,NDOF,NEN,0)
                        QUAD=(IEN(3,n).ne.IEN(4,n))
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

USE MODELDEFINITION, only: numel, numtr
use constants,       only: half
use monitorVars,     only: elmon, nnel, monel

!
! velocity gradients or infinitesimal rotation rate
!
implicit none

integer, parameter :: NTRIA=2
integer, parameter :: NPLT=3, NPLDOF=3, NSUR=4
integer, parameter :: NDOF=2, NSD=2, NEN=4

!-pass
integer            :: luMon
logical            :: CST,verbose

integer            :: ien, ijsur, lmf
double precision   :: x, d, xpt, skew, tfault, deld
dimension          :: IJSUR(NSUR,*),XPT(NPLDOF,NPLT,*), &
                      IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*), &
                      TFAULT(NDOF,*), SKEW(2,*), DELD(NDOF,*)



!-locl
logical            :: QUAD,rotation
double precision   :: BP(NDOF,NSD),XL(NSD,NEN),DL(NDOF,NEN)
integer            :: ITRIA(NTRIA,2)
integer            :: n, nn, m, mm, ntr, i, ierr
double precision   :: val, rxy


!-init
rotation = .false.
NUMTR = 0

if (elmon) then
    do nn=1,nnel
        n = monel(nn)
        call LCOORD (X,XL,IEN(1,n))
        call LDISP (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,124)
        call ADDSNE (DL,n,NDOF,NEN,0)
        QUAD=(IEN(3,n).ne.IEN(4,n))
        call BPMATRIX (XL,DELD,QUAD,BP,ierr)

        if (rotation) then
            rxy = HALF * (BP(2,1)-BP(1,2))
            if (verbose) then
                write(luMon,1) n,rxy
 1              format(1x,'Wxy in element ',I5,' = ', &
                 1PG14.6)
            else
                write(luMon,*) rxy
            endif
        else
            if (verbose) then
                write(luMon,2) n,BP(1,1),BP(1,2),BP(2,1),BP(2,2)
 2              format(1x,'element ',I5, &
                 ': dVx/dx, dVx/dy, dVy/dx, dVy/dy = ', &
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
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,125)
            call ADDSNE (DL,n,NDOF,NEN,0)
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
!           determine triangle element indices
            call triangles (IJSUR,NSUR,ITRIA,NTRIA,n,ntr,QUAD)
!           number of local triangles NTR
            if (ntr.gt.0) then

                call LCOORD (X,XL,IEN(1,n))
                call LDISP (DL,D,IEN(1,n),NDOF,NEN)
                call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,126)
                call ADDSNE (DL,n,NDOF,NEN,0)
                QUAD=(IEN(3,n).ne.IEN(4,n))
                call BPMATRIX (XL,DELD,QUAD,BP,ierr)
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
                        call LDISP (DL,D,IEN(1,n),NDOF,NEN)
                        call AddFaultDisplacement (DL,n,TFAULT,NDOF, &
                                                   NEN,127)
                        call ADDSNE (DL,n,NDOF,NEN,0)
                        QUAD=(IEN(3,n).ne.IEN(4,n))
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

!-----------------------------------

subroutine SPHDISP(XL,DL)
! this sub is a stub.
implicit none
integer, parameter :: NDOF=2, NSD=2, NEN=4
double precision   :: XL(NSD,NEN),DL(NDOF,NEN)
return
end subroutine
