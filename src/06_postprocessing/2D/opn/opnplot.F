subroutine Eigen (ST,EVAL,EVEC)

USE ALGEBRA
use constants, only: dfmin

implicit none

integer, parameter :: NDIM=3, NSTR=5

double precision   :: ST(NSTR),EVAL(NDIM),EVEC(NDIM,NDIM),a(NDIM,NDIM)
logical            :: nonull

integer            :: i, N, nrot

call clear(a,NDIM*NDIM,"a")
nonull = .false.

do i=1,NSTR
    nonull = (nonull .or. ABS(ST(i)).gt.dfmin)
enddo

if (nonull) then
    a(1,1) = ST(1)
    a(2,2) = ST(2)
    a(1,2) = ST(3)
    a(2,1) = ST(3)

    N = NDIM
    call princip (a,N,NDIM,EVAL,EVEC,nrot)
else

    do i=1,NDIM
        EVEC(i,i) = 1d0
        EVAL(i)   = 0d0
    enddo

endif

return
end
!-------------------------------------------------------------------------------
subroutine EigVec (A,EVAL,EVEC,N)

USE ALGEBRA
use constants, only: dfmin

implicit none

!-pass
integer            :: N
double precision   :: A(N),EVAL(3),EVEC(3,3)
!-local
integer            :: i
double precision   :: sum


call clear (EVAL,3,"EVAL")
call clear (EVEC,9,"EVEC")

sum = 0d0

do i=1,N
    sum = sum + A(i)*A(i)
enddo

if (sum.lt.dfmin) then

    do i=1,3
        EVEC(i,i) = 1d0
    enddo

    return
endif

sum = SQRT(sum)

EVEC(1,1) = A(1)/sum
EVEC(2,1) = A(2)/sum
EVEC(1,2) = -EVEC(1,1)
EVEC(2,2) =  EVEC(2,1)

EVAL(1)   = sum

return
end
!-------------------------------------------------------------------------------
subroutine pltcon ()

use constants,       only: isOpnplt
use documentation,   only: plotVariety
USE MODELDEFINITION, only: NSUR, IWORD
use plotControl,     only: nplt, npldof

implicit none

NSUR   = 4
NPLT   = 3
NPLDOF = 3
plotVariety = isOpnplt

return
end
!-------------------------------------------------------------------------------
subroutine wmnmx (Amin,Amax,NSTR,modus)

use iomodule, only: stderr

implicit none

integer          :: NSTR,modus
double precision :: Amin(NSTR),Amax(NSTR)

if (modus.eq.0) then
    write(stderr,1) Amin(1),Amax(1),Amin(2),Amax(2), &
     Amin(3),Amax(3)
 1       format(1x,1PG15.6,' < SXX < ',1PG15.6/1x, &
     1PG15.6,' < SYY < ',1PG15.6/1x, &
     1PG15.6,' < SXY < ',1PG15.6)
elseif (modus.eq.1) then
    write(stderr,2) Amin(1),Amax(1),Amin(2),Amax(2), &
     Amin(3),Amax(3)
 2       format(1x,1PG15.6,' < EXX < ',1PG15.6/1x, &
     1PG15.6,' < EYY < ',1PG15.6/1x, &
     1PG15.6,' < EXY < ',1PG15.6)
elseif (modus.eq.2) then
    write(stderr,3) Amin(1),Amax(1),Amin(2),Amax(2), &
     Amin(3),Amax(3)
 3       format(1x,1PG15.6,' < RXX < ',1PG15.6/1x, &
     1PG15.6,' < RYY < ',1PG15.6/1x, &
     1PG15.6,' < RXY < ',1PG15.6)
else
    write(stderr,4)
 4        format(1x,'wmnmx: fatal mode error')
    call exitp(1)
endif

return
end
!-------------------------------------------------------------------------------
subroutine wrelms (XL,DL,XPT,QUAD,DEFORM,deffac,NXPT)

USE MODELDEFINITION, only: RADIUS
USE CONSTANTS,       only: deg2rad
use plotControl,     only: jShell, shrink, shrinkFactor

implicit none

integer,          parameter :: NSD=2,NDOF=3,NEN=4
double precision, parameter :: UP=3d0, DOWN=2d0
!-pass
integer                     :: NXPT
logical                     :: QUAD, DEFORM
double precision            :: deffac

double precision            :: xl, dl, xpt
dimension                   :: XL(NSD,*),DL(NDOF,*),XPT(3,*)
double precision            :: SH(3,NEN),XS(NSD,NSD),XDL(NSD,NEN)
integer                     :: ierr, m
double precision            :: detJ, f, def

if (DEFORM) then
    def = deffac
else
    def = 0d0
endif

NXPT = NXPT + 1
m = NXPT

if (JSHELL.eq.0) then
    XPT(1,NXPT) = XL(1,1)+def*DL(1,1)
    XPT(2,NXPT) = XL(2,1)+def*DL(2,1)
else
    XPT(1,NXPT) = XL(1,1)+def*DL(1,1)/(DEG2RAD*RADIUS* &
                  COS(DEG2RAD*XL(2,1)))
    XPT(2,NXPT) = XL(2,1)+def*DL(2,1)/(DEG2RAD*RADIUS)
endif

XDL(1,1) = XPT(1,NXPT)
XDL(2,1) = XPT(2,NXPT)

XPT(3,NXPT) = up
NXPT = NXPT + 1

if (JSHELL.eq.0) then
    XPT(1,NXPT) = XL(1,2)+def*DL(1,2)
    XPT(2,NXPT) = XL(2,2)+def*DL(2,2)
else
    XPT(1,NXPT) = XL(1,2)+def*DL(1,2)/(DEG2RAD*RADIUS* &
                  COS(DEG2RAD*XL(2,2)))
    XPT(2,NXPT) = XL(2,2)+def*DL(2,2)/(DEG2RAD*RADIUS)
endif

XDL(1,2) = XPT(1,NXPT)
XDL(2,2) = XPT(2,NXPT)

XPT(3,NXPT) = down
NXPT = NXPT + 1

if (JSHELL.eq.0) then
    XPT(1,NXPT) = XL(1,3)+def*DL(1,3)
    XPT(2,NXPT) = XL(2,3)+def*DL(2,3)
else
    XPT(1,NXPT) = XL(1,3)+def*DL(1,3)/(DEG2RAD*RADIUS* &
                  COS(DEG2RAD*XL(2,3)))
    XPT(2,NXPT) = XL(2,3)+def*DL(2,3)/(DEG2RAD*RADIUS)
endif

XDL(1,3) = XPT(1,NXPT)
XDL(2,3) = XPT(2,NXPT)

XPT(3,NXPT) = down

if (QUAD) then
    NXPT = NXPT + 1
    if (JSHELL.eq.0) then
        XPT(1,NXPT) = XL(1,4)+def*DL(1,4)
        XPT(2,NXPT) = XL(2,4)+def*DL(2,4)
    else
        XPT(1,NXPT) = XL(1,4)+def*DL(1,4)/(DEG2RAD*RADIUS* &
                      COS(DEG2RAD*XL(2,4)))
        XPT(2,NXPT) = XL(2,4)+def*DL(2,4)/(DEG2RAD*RADIUS)
    endif
    XDL(1,4) = XPT(1,NXPT)
    XDL(2,4) = XPT(2,NXPT)
    XPT(3,NXPT) = down
endif

NXPT = NXPT + 1

if (JSHELL.eq.0) then
    XPT(1,NXPT) = XL(1,1)+def*DL(1,1)
    XPT(2,NXPT) = XL(2,1)+def*DL(2,1)
else
    XPT(1,NXPT) = XL(1,1)+def*DL(1,1)/(DEG2RAD*RADIUS* &
                  COS(DEG2RAD*XL(2,1)))
    XPT(2,NXPT) = XL(2,1)+def*DL(2,1)/(DEG2RAD*RADIUS)
endif

XPT(3,NXPT) = down

if (shrink) then

    f = ABS(shrinkFactor)

    call SHAP20 (-f,-f,XDL,detJ,SH,XS,QUAD,ierr)

    XPT(1,m) = SH(3,1)*XDL(1,1) + &
               SH(3,2)*XDL(1,2) + &
               SH(3,3)*XDL(1,3) + &
               SH(3,4)*XDL(1,4)

    XPT(2,m) = SH(3,1)*XDL(2,1) + &
               SH(3,2)*XDL(2,2) + &
               SH(3,3)*XDL(2,3) + &
               SH(3,4)*XDL(2,4)

    XPT(1,NXPT)=XPT(1,m)
    XPT(2,NXPT)=XPT(2,m)

    m = m + 1

    call SHAP20 ( f,-f,XDL,detJ,SH,XS,QUAD,ierr)

    XPT(1,m) = SH(3,1)*XDL(1,1) + &
               SH(3,2)*XDL(1,2) + &
               SH(3,3)*XDL(1,3) + &
               SH(3,4)*XDL(1,4)

    XPT(2,m) = SH(3,1)*XDL(2,1) + &
               SH(3,2)*XDL(2,2) + &
               SH(3,3)*XDL(2,3) + &
               SH(3,4)*XDL(2,4)

    m = m + 1

    call SHAP20 ( f, f,XDL,detJ,SH,XS,QUAD,ierr)

    XPT(1,m) = SH(3,1)*XDL(1,1) + &
               SH(3,2)*XDL(1,2) + &
               SH(3,3)*XDL(1,3) + &
               SH(3,4)*XDL(1,4)

    XPT(2,m) = SH(3,1)*XDL(2,1) + &
               SH(3,2)*XDL(2,2) + &
               SH(3,3)*XDL(2,3) + &
               SH(3,4)*XDL(2,4)

    if (QUAD) then
        m = m + 1
        call SHAP20 (-f, f,XDL,detJ,SH,XS,QUAD,ierr)

        XPT(1,m) = SH(3,1)*XDL(1,1) + &
                   SH(3,2)*XDL(1,2) + &
                   SH(3,3)*XDL(1,3) + &
                   SH(3,4)*XDL(1,4)

        XPT(2,m) = SH(3,1)*XDL(2,1) + &
                   SH(3,2)*XDL(2,2) + &
                   SH(3,3)*XDL(2,3) + &
                   SH(3,4)*XDL(2,4)

 endif

endif

 return
 end
!-------------------------------------------------------------------------------
subroutine surflx (X,IEN,MAT,T,PRPTEM,D,LMF,TFAULT,DXE,&
                   TX,LMTX,SKEW,VERBOSE,luMon,iside,ierr)

USE MODELDEFINITION
use constants,       only: dfmin,dlogmin,dfmax,dlogmax, zero, one, useDXE
use materialsmodule, only: MaterialIDOfElement
use iomodule,        only: stderr
use monitorVars,     only: elmon, nnel, monel
use errors,          only: exec, ok

! Routine to calculate surface heat flux magnitudes

implicit none

integer, parameter :: NSD=2
integer, parameter :: NDOF=3
integer, parameter :: NEN=4

!-pass
integer            :: ierr,luMon,iside
logical            :: VERBOSE

integer            :: ien, mat
double precision   :: x, t, prptem, d, lmf, tfault, dxe, tx, lmtx, skew

dimension          :: X(NSD,*),IEN(NEN,*),MAT(*),T(*),PRPTEM(6,*), &
                      D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
                      DXE(NDOF,NEN,*), TX(*),LMTX(NEN,*),SKEW(2,*)


!-locl
logical            :: QUAD
double precision   :: GRAD(NSD),XL(NSD,NEN),DL(NDOF,NEN),TL(NEN),XELM(NSD)

integer            :: sn, tn, nn, i, m, n
double precision   :: flux

!-init
ierr  = OK

if (.not.elmon) then
    write(stderr,10)
10        format(1x,'surflx: no monitor elements set')
    ierr = EXEC
    return
endif

do nn=1,nnel
    n = monel(nn)
    m=MaterialIDOfElement(n)
    QUAD = (IEN(3,n).ne.IEN(4,n))

!   Set natural coordinates
    if (iside.le.1) then
        sn = zero
        tn = -one
    else if (iside.eq.2) then
        sn = one
        tn = zero
    else if (iside.eq.3) then
        sn = zero
        tn = one
    else if (iside.ge.4) then
        if (QUAD) then
            sn = -one
            tn = zero
        else
            sn = zero
            tn = one
        endif
    endif

    call LCOORD (X,XL,IEN(1,n))
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
    call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,159)
    call ADDSNE (DL,n,NDOF,NEN,useDXE)
    call REZONE (XL,DL,1d0)
    call LTEMP  (IEN(1,n),T,TL,NEN)

    if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
        call LTEMPX (LMTX(1,n),IEN(1,n),TX,TL,NEN)
    endif

    call thgrad (sn,tn,XL,TL,GRAD,QUAD,ierr)

    if (ierr.ne.0) then
        write(stderr,20) n
20           format(1x,'srflx: error in element ',I5)
        ierr = EXEC
        return
    endif

    FLUX = zero

    do i=1,NSD
        FLUX = FLUX + (PRPTEM(3+i,m)*GRAD(i))**2
    enddo

    if (FLUX.gt.dfmin) then
        FLUX=SQRT(FLUX)
    endif

    call ELCRD (sn,tn,XELM,IEN(1,n),X,D,n,TFAULT, &
                SKEW,one, ierr)

    if (VERBOSE) then
        write(luMon,30) iside,n,XELM(1),XELM(2),FLUX
30           format(1x,'Heat flux magnitude at side ',I1,' of',1x, &
         'element ',I5,' = ',3(1PG15.6,2X))
    else
        write(luMon,*) XELM(1),XELM(2),FLUX
    endif

enddo

return
end subroutine
!-----------------------------------------------------------------------
subroutine nsmpl (IEN,X,D,LMF,TFAULT, &
                 SKEW,xr,yr,QNTY,M,QS,deformed, &
                 lastSampledElement, samplePointFound)

USE MODELDEFINITION
use constants,       only: eps, useDXE
use modeltopology,   only: nen, nsd, ndof
use debugmodule,     only: iecho
use iomodule,        only: stderr
use errors,          only: ok
!
! Routine to sample nodal point quantities on a coordinate
!
implicit none

!-pass
integer          :: M
logical          :: deformed
double precision :: XR,YR,QS(*)

integer          :: ien, lmf
double precision :: x, qnty, d, tfault, skew
dimension        :: X(NSD,*),IEN(NEN,*),QNTY(M,*), &
                    D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
                    SKEW(2,*)
integer          :: lastSampledElement
logical          :: samplePointFound


!-locl
logical          :: QUAD
integer          :: i,inside, ierr, k, n
double precision :: XL(NSD,NEN),DL(NDOF,NEN),QL(NDOF,NEN)
double precision :: w1,w2,w3, area, dist

!-init
samplePointFound = .false.
n = 1

100    if (n.gt.NUMEL) then
    return
endif

QUAD = (IEN(3,n).ne.IEN(4,n))

if (QUAD) then
    write(stderr,1)
 1        format(1x,'sample: quad elements not implemented')
    return
endif

call LCOORD (X,XL,IEN(1,n))
if (deformed) then
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
    call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,160)
    call ADDSNE (DL,n,NDOF,NEN,useDXE)
    call REZONE (XL,DL,1d0)
endif


INSIDE = -1

w1 = 0d0
w2 = 0d0
w3 = 0d0

area = 1d0
dist = SQRT((XR-XL(1,1))**2+(YR-XL(2,1))**2)

if (dist.le.EPS) then
    INSIDE = 0
    w1 = 1d0
endif

if (INSIDE.eq.-1) then
    dist = SQRT((XR-XL(1,2))**2+(YR-XL(2,2))**2)
    if (dist.le.EPS) then
        INSIDE = 0
        w2 = 1d0
    endif
endif

if (INSIDE.eq.-1) then
    dist = SQRT((XR-XL(1,3))**2+(YR-XL(2,3))**2)
    if (dist.le.EPS) then
        INSIDE = 0
        w3 = 1d0
    endif
endif

if (INSIDE.eq.-1) then
    call PNPOLY (XL,XR,YR,area,w1,w2,w3,INSIDE,ierr,n)
    if (ierr.ne.OK) then
        write(stderr,10) n
10            format(1x,'sample: degenerate element ',I8)
        return
    endif
endif

if (INSIDE.eq.-1) then
    ! LWITODO: select an element in the direction of the sample point
    ! and not just increment by one.
    n = n + 1
    goto 100
endif

samplePointFound = .true.


w1 = w1 / area
w2 = w2 / area
w3 = w3 / area

if (iecho.gt.1) then
    if (INSIDE.eq.1) then
        write(stderr,20) n,area,w1,w2,w3
20   format(1x,'sample: coordinates in element ',I8/9x, &
 'element area = ',1PG16.6/9x,'weights are ',3(1PG16.6))
    elseif (INSIDE.eq.0) then
        write(stderr,30) n,area,w1,w2,w3
30   format(1x,'sample: coordinates along element ',I8/9x, &
 'element area = ',1PG16.6/9x,'weights are ',3(1PG16.6))
    endif
endif

do i=1,M
    k = IEN(1,n)
    QL(i,1) = QNTY(i,k)
    k = IEN(2,n)
    QL(i,2) = QNTY(i,k)
    k = IEN(3,n)
    QL(i,3) = QNTY(i,k)
enddo

do i=1,M
    QS(i) = w1*QL(i,1)+w2*QL(i,2)+w3*QL(i,3)
enddo

return
end subroutine
!-----------------------------------------------------------------------
subroutine esmpl (IEN,X,D,LMF,TFAULT, &
                  SKEW,XR,YR,QNTY,M,QS,deformed,proc, &
                  lastSampledElement, samplePointFound)

USE MODELDEFINITION
use constants,       only: eps, useDXE
use debugmodule,     only: iecho
use iomodule,        only: stderr
use errors,          only: ok

implicit none

integer, parameter :: NEN=4,NSD=2,NDOF=3
!-pass
integer            :: M
logical            :: deformed, proc
double precision   :: XR,YR,QS(M)

integer            :: ien, lmf
double precision   :: x, qnty, d, tfault, skew
dimension          :: X(NSD,*),IEN(NEN,*),QNTY(*), &
                      D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
                      SKEW(2,*)
integer            :: lastSampledElement
logical            :: samplePointFound

!-locl
logical            :: QUAD
double precision   :: XL(NSD,NEN),DL(NDOF,NEN),QL(NDOF,NEN),BP(NDOF,NSD)
double precision   :: w1, w2, w3, area, dist
integer            :: i, k, n, ierr
integer            :: INSIDE

!-init
samplePointFound = .false.
n = 1

100 if (n.gt.NUMEL) then
    return
endif

QUAD = (IEN(3,n).ne.IEN(4,n))

if (QUAD) then
    write(stderr,1)
 1        format(1x,'sample: quad elements not implemented')
    return
endif

call LCOORD (X,XL,IEN(1,n))

if (deformed.or.proc) then
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
    call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,161)
    call ADDSNE (DL,n,NDOF,NEN,useDXE)
    if (.not.proc) then
       call REZONE (XL,DL,1d0)
    endif
endif

INSIDE = -1
dist = SQRT((XR-XL(1,1))**2 + &
            (YR-XL(2,1))**2)

if (dist.le.EPS) then
    INSIDE = 0
endif

if (INSIDE.eq.-1) then
    dist = SQRT((XR-XL(1,2))**2 + &
                (YR-XL(2,2))**2)
    if (dist.le.EPS) then
        INSIDE = 0
    endif
endif

if (INSIDE.eq.-1) then
    dist = SQRT((XR-XL(1,3))**2 + &
                (YR-XL(2,3))**2)
    if (dist.le.EPS) then
        INSIDE = 0
    endif
endif

if (INSIDE.eq.-1) then
    call PNPOLY (XL,XR,YR,area,w1,w2,w3,INSIDE,ierr,n)
    if (ierr.ne.OK) then
        write(stderr,10) n
10            format(1x,'sample: degenerate element ',I8)
        return
    endif
endif

if (INSIDE.eq.-1) then
    ! LWITODO: select an element in the direction of the sample point
    ! and not just increment by one.
    n = n + 1
    goto 100
endif

samplePointFound = .true.


!if (iecho.gt.1) then
!    if (INSIDE.eq.1) then
!        write(stderr,20) n
!20           format(1x,'sample: coordinate in element ',I8)
!    else
!        write(stderr,30) n
!30           format(1x,'sample: coordinate along element ',I8)
!    endif
!endif

if (proc) then
!           post-process: determine displacement gradients (for rot)
    call BPMATRIX (XL,DL,QUAD,BP,ierr)

    if (iecho.gt.1) then
        write(stderr,40) BP(1,1),BP(1,2),BP(2,1),BP(2,2)
  40            format(1x,'dUx/dx, dUx/dy, dUy/dx, dUy/dy=', &
          4(1PG16.6))
    endif

    QS(1) = 5D-1 * (BP(2,1)-BP(1,2))
else
    do i=1,M
        k = i + (n-1)*M
        QS(i) = QNTY(k)
    enddo
endif
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine dsmpl (IEN,X,D,DELD,LMF,TFAULT,DFAULT, &
                  SKEW,XR,YR,QS,IVELO,deformed, &
                  lastSampledElement, samplePointFound)

USE MODELDEFINITION
use constants,       only: eps, useDXE, useDELDXE
use iomodule,        only: stderr
use errors,          only: ok
!
! Routine to sample displacements/velocities on coordinate XR, YR
! different from nsmpl because there can be multiple displacement/
! velocity components
!
! IVELO = 1 displ
!       = 2 velo
!
implicit none

integer, parameter :: NEN=4,NSD=2,NDOF=3
!-pass
logical            :: deformed
integer            :: IVELO
double precision   :: XR,YR

integer            :: ien, lmf
double precision   :: x, d, deld, tfault, dfault, skew, qs
dimension          :: X(NSD,*),IEN(NEN,*),D(NDOF,*),DELD(NDOF,*), &
                      LMF(NDOF,NEN,*),TFAULT(NDOF,*),DFAULT(NDOF,*), &
                      SKEW(2,*),QS(NDOF)
integer            :: lastSampledElement
logical            :: samplePointFound


!-locl
logical            :: QUAD
integer            :: INSIDE,ierr,n
double precision   :: dist,area,w1,w2,w3
double precision   :: XL(NSD,NEN),DL(NDOF,NEN),QL(NDOF,NEN)
integer            :: i

!-init
samplePointFound = .false.
n = 1

100     if (n.gt.NUMEL) then
    return
endif

QUAD = (IEN(3,n).ne.IEN(4,n))

if (QUAD) then
    write(stderr,1)
 1       format(1x,'sample: quad elements not implemented')
    return
endif

call LCOORD (X,XL,IEN(1,n))

if (deformed) then
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
    call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,162)
    call ADDSNE (DL,n,NDOF,NEN,useDXE)
    call REZONE (XL,DL,1d0)
endif

INSIDE = -1

! check whether sample point is close to vertex 1 of elem n
dist = SQRT((XR-XL(1,1))**2 + &
            (YR-XL(2,1))**2)
if (dist.le.EPS) then
    INSIDE = 0
endif

! check whether sample point is close to vertex 2 of elem n
if (INSIDE.eq.-1) then
    dist = SQRT((XR-XL(1,2))**2 + &
                (YR-XL(2,2))**2)
    if (dist.le.EPS) then
        INSIDE = 0
    endif
endif

! check whether sample point is close to vertex 3 of elem n
if (INSIDE.eq.-1) then
    dist = SQRT((XR-XL(1,3))**2 + &
                (YR-XL(2,3))**2)
    if (dist.le.EPS) then
        INSIDE = 0
    endif
endif

! check whether sample point is inside element n
if (INSIDE.eq.-1) then
    call PNPOLY (XL,XR,YR,area,w1,w2,w3,INSIDE,ierr,n)
    if (ierr.ne.OK) then
        write(stderr,10) n
10           format(1x,'sample: degenerate element ',I8)
        return
    endif
endif

! if not, may better luck with next elem.
if (INSIDE.eq.-1) then
    ! LWITODO: select an element in the direction of the sample point
    ! and not just increment by one.
    n = n + 1
    goto 100
endif

! element found. point in element n.

!if (iecho.gt.1) then
!    if (INSIDE.eq.1) then
!        write(stderr,20) n
!20           format(1x,'sample: coordinate in element ',I12)
!    else
!        write(stderr,30) n
!30           format(1x,'sample: coordinate along element ',I12)
!    endif
!endif

samplePointFound = .true.


w1 = w1 / area
w2 = w2 / area
w3 = w3 / area

if (IVELO.eq.1) then
    call LDISP (QL,D,IEN(1,n),NDOF,NEN)
    call ADDSNE (QL,n,NDOF,NEN,useDXE)
    call AddFaultDisplacement (QL,n,TFAULT,NDOF,NEN,163)
else
    call LDISP (QL,DELD,IEN(1,n),NDOF,NEN)
    call ADDSNE (QL,n,NDOF,NEN,useDELDXE)
    call AddFaultDisplacement (QL,n,DFAULT,NDOF,NEN,164)
endif

do i=1,NDOF
    QS(i) = w1*QL(i,1)+w2*QL(i,2)+w3*QL(i,3)
enddo
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine PNPOLY (X,XR,YR,D0,D1,D2,D3,INOUT,IERR, elemID)
!
use constants,       only: eps,zero,one
use errors,          only: exec, ok

implicit none

integer            :: elemID
integer, parameter :: NSD=2,NEN=3
!-pass
integer            :: INOUT,IERR
double precision   :: X(NSD,NEN),XR,YR,D0,D1,D2,D3


!-locl
integer            :: i0,i1,i2,i3
double precision   :: A(3,3)
integer, external  :: NINT

!-init
IERR = OK
INOUT = -1

A(1,1) = X(1,1)
A(1,2) = X(2,1)
A(1,3) = ONE
A(2,1) = X(1,2)
A(2,2) = X(2,2)
A(2,3) = ONE
A(3,1) = X(1,3)
A(3,2) = X(2,3)
A(3,3) = ONE

D0 = A(1,1) * (A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
   - A(2,1) * (A(1,2)*A(3,3)-A(1,3)*A(3,2)) &
   + A(3,1) * (A(1,2)*A(2,3)-A(1,3)*A(2,2))

if (ABS(D0).le.EPS) then
    ierr = EXEC
    goto 100
endif

A(1,1) = XR
A(1,2) = YR

D1 = A(1,1) * (A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
   - A(2,1) * (A(1,2)*A(3,3)-A(1,3)*A(3,2)) &
   + A(3,1) * (A(1,2)*A(2,3)-A(1,3)*A(2,2))

if (D1.eq.ZERO) then
    INOUT = 0
endif

A(1,1) = X(1,1)
A(1,2) = X(2,1)
A(2,1) = XR
A(2,2) = YR

D2 = A(1,1) * (A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
   - A(2,1) * (A(1,2)*A(3,3)-A(1,3)*A(3,2)) &
   + A(3,1) * (A(1,2)*A(2,3)-A(1,3)*A(2,2))

if (D2.eq.ZERO) then
    INOUT = 0
endif

A(2,1) = X(1,2)
A(2,2) = X(2,2)
A(3,1) = XR
A(3,2) = YR

D3 = A(1,1) * (A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
   - A(2,1) * (A(1,2)*A(3,3)-A(1,3)*A(3,2)) &
   + A(3,1) * (A(1,2)*A(2,3)-A(1,3)*A(2,2))

if (D3.eq.ZERO) then
    INOUT = 0
endif

i0 = INT(SIGN(1.1d0,D0))

if (D1.ne.ZERO) then
    i1 = INT(SIGN(1.1d0,D1))
else
    i1 = i0
endif

if (D2.ne.ZERO) then
    i2 = INT(SIGN(1.1d0,D2))
else
    i2 = i0
endif

if (D3.ne.ZERO) then
    i3 = INT(SIGN(1.1d0,D3))
else
    i3 = i0
endif

if (i0*i1.gt.0.and.i0*i2.gt.0.and.i0*i3.gt.0) then
    if (INOUT.eq.-1) then
        INOUT = 1
    endif
else
    INOUT = -1
endif

100    return
end subroutine
!-----------------------------------------------------------------------
logical function PosOrient (XL)

use constants, only: dfmin

! Determines if triangular element is positively oriented

implicit none

integer, parameter :: NSD=2
!-pass
double precision   :: XL
dimension          :: XL(NSD,*)
!-locl
double precision   :: V1x, V1y
double precision   :: V2x, V2y
double precision   :: outprd

!    Create vector from node 1 to 2
V1x = XL(1,2) - XL(1,1)
V1y = XL(2,2) - XL(2,1)

!    Create vector from node 2 to 3
V2x = XL(1,3) - XL(1,2)
V2y = XL(2,3) - XL(2,2)

outprd = V1x * V2y - V1y * V2x
PosOrient = (outprd.gt.-dfmin)

return
end function

!-----------------------------------------------------------------------

subroutine shslp (NSLIP,NDIM,IEN,X,NUMSLP,lumon,parallel)
!
implicit none

integer, parameter :: NEN=4,NSD=2
!-pass
logical            :: parallel
integer            :: NUMSLP, lumon, ndim

integer            :: nslip, ien
double precision   :: x
dimension          :: NSLIP(NDIM,*),IEN(NEN,*),X(NSD,*)
!-locl
double precision   :: XL(NSD,NEN)
logical            :: QUAD
integer            :: i, j, k
integer            :: ne, np
double precision   :: dx, dy


if (NUMSLP.le.0) then
    return
endif

do i=1,NUMSLP
    ne = NSLIP(1,i)
    np = NSLIP(2,i)
    k = 1
    !        It has been checked before that np occurs in IEN(.,ne)

! todo make a nice while loop.
100 if (IEN(k,ne).ne.np) then
        k = k + 1
        goto 100
    endif

    QUAD = (IEN(3,ne).ne.IEN(4,ne))
    call LCOORD (X,XL,IEN(1,ne))
    j = k - 1
    if (j.lt.1) then
        if (QUAD) then
            j = NEN
        else
            j = NEN-1
        endif
    endif

    dx = XL(1,k)-XL(1,j)
    dy = XL(2,k)-XL(2,j)

    if (parallel) then
        write(lumon,1) np,ne,j,dx,dy
 1           format(1x,'Node ',I6,': tangent to element ',I6, &
         ', side ',I1,' =',2(1x,1PG14.6))
    else
        write(lumon,2) np,ne,j,-dy,dx
 2           format(1x,'Node ',I6,': normal to element ',I6, &
         ', side ',I1,' =',2(1x,1PG14.6))
    endif

    j = k
    k = k + 1

    if (QUAD) then
        if (k.gt.NEN) then
            k = 1
        endif
    else
        if (k.gt.NEN-1) then
            k = 1
        endif
    endif

    dx = XL(1,k)-XL(1,j)
    dy = XL(2,k)-XL(2,j)

    if (parallel) then
        write(lumon,1) np,ne,j,dx,dy
    else
        write(lumon,2) np,ne,j,-dy,dx
    endif
enddo

return
end subroutine
!-------------------------------------------------------------------------------
