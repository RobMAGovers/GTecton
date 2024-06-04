!    These routines are the dimensionality-dependent portion of GTECTON
!    The routines here are for axisymmetric stress and strain
!-------------------------------------------------------------------------------
!
! DEFINITIONS
!
! Nodal points:                            Faces:
!            4             3                        3
!            o-------------o                 o-------------o
!            |             |                 |             |
!            |             |                 |             |
!            |             |                 |             |
!            |             |               4 |             | 2
!            |             |                 |             |
!            |             |                 |             |
!            |             |                 |             |
!            o-------------o                 o-------------o
!            1             2                        1
!
!
! Collapse of nodal points 4 onto 3 (i.e. assigning identical global node
! numbers to nodes 3 and 4) gives a triangle-element.
!-------------------------------------------------------------------------------
!    THE ROUTINES LISTED HERE ARE:
!
!    NAME                    FUNCTION
!
!
!    ADDFLX        NORMAL BOUNDARY FLUX CONTRIBUTIONS TO LOAD VECTOR
!    ADDPR        TRANSLATES APPLIED PRESSURES INTO LOAD VECTOR
!    ADDSTR        STRESS LOAD COMPUTATION
!    ADDWNK        ADDS WINKLER PRESSURES TO STIFFNESS MATRIX
!    AJ1        FIRST STRESS TENSOR INVARIANT
!    AJ2        SECOND STRESS TENSOR INVARIANT
!    AJ3        THIRD STRESS TENSOR INVARIANT
!    BDELD        STRAIN FROM DISPLACEMENTS
!    BMATRIX        STRAIN-DISPLACEMENT MATRIX
!    CPSTIFF        LOCAL HEAT CAPACITY STIFFNESS MATRIX
!    DEVSTN        STRESS DEVIATOR CALCULATION
!    EFORCE        EFFECTIVE LOAD VECTOR FOR VISCOUS FLOW
!    ELCRD        CALCULATES COORDINATES OF DEFORMED ELEMENT
!    ELCTR        CALCULATES CENTER COORDINATES OF DEFORMED ELEMENT
!    EIGMAX        CALCULATES MAXIMUM EIGENVALUE
!    FCGEOM        COMPUTES FACE NORMAL VECTOR AND FACE AREA
!    FLOWVP        VISCOPLASTIC STRAIN RATE
!    FORMBT        VISCOUS FLOW CORRECTION FOR STRAIN
!    FORMJB        JACOBEAN MATRIX FOR VISCOUS FLOW
!    GRAVLD        GRAVITATIONAL BODY FORCE COMPUTATION
!    INTCON        DIMENSIONAL AND TOPOLOGIC DATA
!    YIELDS        COMPUTES SCALAR STRESS FOR COMPARISON WITH YIELD STRESS
!    LCOORD        LOCALIZES COORDINATES
!    LFLUX        ELEMENT HEAT FLOW
!    LHEAT        HEAT PRODUCTION CONTRIBUTION TO LOAD VECTOR
!    FLUID        DETERMINES VISCOPLASTIC FLUIDITY
!    MATERL        MATERIAL MATRICES
!    REZONE        LOCAL REZONE FOR LARGE DEFORMATION
!    SDIMNT        SURFACE ELEVATION AND LOADS FROM SEDIMENTARY TRANSPORT
!    SRFHGT        DETERMINES SURFACE HIGHT AT A GIVEN LOCATION
!    SHAP20        SHAPE FUNCTION FOR 2-D ELEMENTS
!    SHAREL        IDENTIFIES NEIGHBOURING ELEMENTS
!    SIDENP        RETURNS LOCAL NODE NUMBERS ON A SPECIFIC SIDE
!    SHLCRD        COMPUTES COORDINATES IN LOCAL (SPHERICAL) FRAME
!    STIFF        STIFFNESS MATRIX
!    TANGVP        VISCOPLASTIC JACOBIAN MATRIX
!    TCENTR        COMPUTES ELEMENT CENTER TEMPERATURES
!    THGRAD        THERMAL GRADIENTS
!    TSTIFF        LOCAL CONDUCTION STIFFNESS MATRIX
!    WINKLP        WINKLER PRESSURES CONTRIBUTION TO GLOBAL FORCE VECTOR
!    WAPPLY        APPLIES RESIDUAL WINKLER FORCES
!    YIELDF        CONSTRUCTS PLASTIC FLOW VECTOR
!
!**********************************************************************

subroutine INTCON ()

USE LABELMODULE
USE MODELDEFINITION, only: NTYPE, NSUR, NEE, NEC, simVariety, &
                           allTemperaturesFixed
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use constants,       only: isAxi

! Set variables that depend on topology

 implicit none

!-pass

simVariety = isAxi

 NTYPE = 3
 NSD   = 2
 NDOF  = 2
!    write(*,*) 'INTCON set NDOF: ', NDOF
 NSTR  = 4
 NEN   = 4
 NSUR  = 4
 NEC   = NSD*NEN
 NEE   = NDOF*NEN
!      call labelsinit('pln')
   LABELC(1)='  R '
   LABELC(2)='  Z '
   LABELD(1)=LABELC(1)
   LABELD(2)=LABELC(2)
   LABELS(1)='SRR '
   LABELS(2)='SZZ '
   LABELS(3)='STT '
   LABELS(4)='SRZ '
   LABELS(5)='    '
   LABELS(6)='    '
   LABELR(1)='RRR '
   LABELR(2)='RZZ '
   LABELR(3)='RTT '
   LABELR(4)='RRZ '
   MATLBL(1)='YOUNGS MOD'
   MATLBL(2)='POIS RATIO'
   MATLBL(3)='VISCOSITY '
   MATLBL(4)='POWER     '
   MATLBL(5)='DENSITY   '
   MATLBL(6)='THICKNESS '
   CRPLBL(1)='CREEP ACT '
   CRPLBL(2)='CREEP PRE '
   CRPLBL(3)='CREEP PWR '
   CRPLBL(4)='BRITTLE MU'
   CRPLBL(5)='BRITTLE S0'
   TEMLBL(1)='HEAT PROD '
   TEMLBL(2)='SPEC HEAT '
   TEMLBL(3)='THERM EXP '
   TEMLBL(4)='R-CONDUCT'
   TEMLBL(5)='Z-CONDUCT'

allTemperaturesFixed = .false.

 return
 end
!-------------------------------------------------------------------------------
SUBROUTINE SHAP20 (S,T,X,DET,SH,XS,QUAD,ierr)

use constants, only: half, zero
use modeltopology,   only: nen, nsd
use iomodule,        only: stderr

! Program to compute shape functions for quadrilateral or triangle

!    S,T,                = natural coordinates
!    SH(1,NEN),SH(2,NEN) = r and z derivatives of shape functions
!    SH(3,NEN)           = shape functions
!    XS(NSD,NSD)         = inverse jacobian matrix
!    DET                 = Jacobian matrix determinant
!    X(2,NEN)            = global coordinates

implicit none

!-pass
double precision :: s,t
double precision :: SH(3,NEN),X(NSD,NEN),XS(NSD,NSD)
double precision :: DET
logical QUAD
integer ierr

!-locl
double precision :: SA(NEN),TA(NEN)
!save SA,TA
integer          :: i,j,k
double precision :: temp





!-init
!data SA/-0.5d0, 0.5d0,0.5d0,-0.5d0/
!data TA/-0.5d0,-0.5d0,0.5d0, 0.5d0/
SA=(/-0.5d0, 0.5d0,0.5d0,-0.5d0/)
TA=(/-0.5d0,-0.5d0,0.5d0, 0.5d0/)


ierr = 0

!    Calculate shape functions and s- and t-derivatives of shape functions
!    these derivatives are local.
do i=1,NEN
    SH(1,i) = SA(i)*(half+TA(i)*T)
    SH(2,i) = TA(i)*(half+SA(i)*S)
    SH(3,i) = (half+SA(i)*S)*(half+TA(i)*T)
enddo

if (.not.QUAD) then
!        adjust shape function of third node for triangular elements
    do i=1,3
        SH(i,3) = SH(i,3)+SH(i,4)
        SH(i,4) = zero
    enddo
endif

!    Jacobian matrix calculation
do i=1,2
    do j=1,2
        XS(i,j) = zero
        do  k=1,NEN
            XS(i,j) = XS(i,j) + SH(i,k)*X(j,k)
        enddo
    enddo
enddo

DET = XS(1,1)*XS(2,2)-XS(1,2)*XS(2,1)


if (DET.le.zero) then
    write(stderr,*) 'Shape function fails! Determinant is: ',DET
    write(stderr,*) 'due to matrix entries: '
    write(stderr,*) XS(1,1), XS(1,2)
    write(stderr,*) XS(2,1), XS(2,2)
    write(stderr,*) 'Cause: too much deformation. Solve by reducing force or timespan.'
    ierr = 1
    return
endif

!    Calculatate inverse Jacobian determinant (Cramers rule)
TEMP    =  XS(1,1)/DET
XS(1,1) =  XS(2,2)/DET
XS(2,2) =  TEMP
XS(1,2) = -XS(1,2)/DET
XS(2,1) = -XS(2,1)/DET

do i=1,NEN
    TEMP    = XS(1,1)*SH(1,i) + & 
              XS(1,2)*SH(2,i)
    SH(2,i) = XS(2,1)*SH(1,i) + & 
              XS(2,2)*SH(2,i)
    SH(1,i) = TEMP
enddo


return
end

!-------------------------------------------------------------------------------
SUBROUTINE ELMTEST (X,IEN)

USE MESHDATAMODULE
USE MODELDEFINITION
use constants, only: eps, zero, half
use debugmodule, only: iecho
use modeltopology,   only: nen, nsd
use iomodule

!
! tests element definitions
!
 implicit none

!-pass
 integer IEN
 double precision X
 dimension X(NSD,*),IEN(NEN,*)

!-locl
 logical quad
 integer i,n
 double precision xl,area1,area2,area,a,b
 dimension xl(NSD,NEN),a(NSD),b(NSD)
!
!    write (*,*) 'Entered ELMTEST'
 
if (meshdatactx%Neglobal.lt.1) then
!    if (NUMEL.lt.1) then
    write(stderr,10)
10       format(1x,'NUMEL < 1')
   call exitp(1)
endif
!
do n=1,meshdatactx%Nelocal

    do i=1,NEN
        if (IEN(i,n).lt.1 .or. &
            IEN(i,n).gt.meshdatactx%Nvglobal) then
            write(stderr,20) n,IEN(i,n)
20               format(1x,'element ',I8,': reference to node ', &
             'number ',I8)
            call exitp(1)
        endif
    enddo

    quad = (IEN(3,n).ne.IEN(4,n))

!       write (*,*) 'elmtest calling LCOORD'
    call LCOORD (X,xl,IEN(1,n),n)


    a(1) = xl(1,2) - xl(1,1)
    a(2) = xl(2,2) - xl(2,1)
    b(1) = xl(1,3) - xl(1,2)
    b(2) = xl(2,3) - xl(2,2)
    area1 = HALF*(a(1)*b(2)-a(2)*b(1))

    if (quad) then
        a(1) = xl(1,4) - xl(1,3)
        a(2) = xl(2,4) - xl(2,3)
        b(1) = xl(1,1) - xl(1,4)
        b(2) = xl(2,1) - xl(2,4)
        area2 = HALF*(a(1)*b(2)-a(2)*b(1))
        area = area1 + area2
    else
        area2 = 0
        area = area1
    endif

    if (ABS(area).le.EPS) then
        write(stderr,30) n
 !           write(*,*) a(1),a(2),b(1),b(2),ABS(area)
 !           write(*,*) 'xl(1,1) xl(2,1): ', xl(1,1), xl(2,1)
 !           write(*,*) 'xl(1,2) xl(2,2): ', xl(1,2), xl(2,2)
 !           write(*,*) 'xl(1,3) xl(2,3): ', xl(1,3), xl(2,3)
30       format(1x,'degenerate element ',I8)
        call exitp(1)
    endif

    if (area.lt.ZERO) then
        write(stderr,40) n
40           format(1x,'element ',I8,' is negatively oriented (node', &
         ' order is wrong)')
        call exitp(1)
    endif

    if (ABS(area1).le.EPS) then
        write(stderr,50) n
50           format(1x,'WARNING: degenerate triangle 123 in element', &
         1X,I8)
        call exitp(1)
    endif

    if (area1.lt.ZERO) then
        write(stderr,60) n
60           format(1x,'WARNING: negatively oriented triangle 123', &
         ' in element ',I8)
        call exitp(1)
    endif

    if (quad .and. ABS(area2).le.EPS) then
        write(stderr,70) n
70           format(1x,'WARNING: degenerate triangle 341 in element', &
         1X,I8)
        call exitp(1)
    endif

    if (quad .and. area2.lt.ZERO) then
        write(stderr,80) n
80           format(1x,'WARNING: negatively oriented triangle 341', &
         ' in element ',I8)
        call exitp(1)
    endif

enddo

if (iecho.eq.2) then
    write(stderr,90)
90   format(1x,'Elements tested ok')
endif

return
end

!-------------------------------------------------------------------------------
subroutine MATERL (DMAT,E,POIS)

USE MODELDEFINITION
USE constants, only: zero, one, two, three, half

implicit none

integer, parameter :: NSTR=4
double precision :: AM, AL
!
! Constructs the material matrix relating stress and strain.
! IOPT=0, plane strain
! IOPT=1, plane stress
!
!-pass
double precision :: DMAT(NSTR,NSTR)
double precision :: POIS, E
!
!            Compute Lame parameters AM(=MU) and AL(=LAMBDA)
AM=E/(ONE+POIS)
AL=AM*POIS/(ONE-TWO*POIS)
AM=HALF*AM

DMAT(1,1)=TWO*AM+AL
DMAT(1,2)=AL
DMAT(1,3)=AL

DMAT(1,4)=ZERO

DMAT(2,1)=AL
DMAT(2,2)=TWO*AM+AL
DMAT(2,3)=AL

DMAT(2,4)=ZERO

DMAT(3,1)=AL
DMAT(3,2)=AL 
DMAT(3,3)=TWO*AM+AL

DMAT(3,4)=ZERO

DMAT(4,1)=ZERO
DMAT(4,2)=ZERO
DMAT(4,3)=ZERO
DMAT(4,4)=AM

return
end
!-------------------------------------------------------------------------------
subroutine BDELD (XL,DL,EE,QUAD,ierr)

use constants, only: eps, zero
use modeltopology,   only: nen, ndof, nsd, nstr
use iomodule,        only: stderr

!
! Subroutine to compute strains from displacements in each element
!
implicit none

!-pass
integer ierr
double precision :: XL(NSD,NEN),DL(NDOF,NEN),EE(NSTR)
logical QUAD
!-locl
double precision :: R
double precision :: DET
double precision :: BP(NDOF,NSD)


!-locl
double precision :: SH(3,NEN),XS(NSD,NSD)
logical firstentry
save firstentry
data firstentry/.true./        
!
call SHAP20 (ZERO,ZERO,XL,DET,SH,XS,QUAD,ierr)
if (ierr.ne.0) then
    return
endif

R = SH(3,1)*XL(1,1) + &
    SH(3,2)*XL(1,2) + &
    SH(3,3)*XL(1,3) + &
    SH(3,4)*XL(1,4)

if (R.lt.EPS) then
    write(stderr,10) EPS
10  format(///1x,'BDELD: R < ',1PE20.6)
    ierr = 1
    return
endif

EE(1) = SH(1,1)*DL(1,1) + &
        SH(1,2)*DL(1,2) + &
        SH(1,3)*DL(1,3) + &
        SH(1,4)*DL(1,4)

EE(2) = SH(2,1)*DL(2,1) + &
        SH(2,2)*DL(2,2) + &
        SH(2,3)*DL(2,3) + &
        SH(2,4)*DL(2,4)

EE(3) =(SH(3,1)*DL(1,1) + & 
        SH(3,2)*DL(1,2) + &
        SH(3,3)*DL(1,3) + &
        SH(3,4)*DL(1,4) ) / R

EE(4) = SH(2,1)*DL(1,1) + &
        SH(2,2)*DL(1,2) + &
        SH(2,3)*DL(1,3) + &
        SH(2,4)*DL(1,4) + &
        SH(1,1)*DL(2,1) + &
        SH(1,2)*DL(2,2) + &
        SH(1,3)*DL(2,3) + &
        SH(1,4)*DL(2,4)

return
end

!-------------------------------------------------------------------------------
subroutine BMATRIX (B,SH,R)

use modeldefinition, only: nee
use modeltopology,   only: nen, ndof, nstr


!
! Computes the strain-displacement matrix B
!
implicit none

!-pass
double precision :: B(NSTR,NEE)
double precision :: SH(3,NEN)
double precision :: R

integer :: i,k

k = 1

do i=1,NEN

    B(1,k  ) = SH(1,i)
    B(1,k+1) = 0d0
    B(2,k  ) = 0d0
    B(2,k+1) = SH(2,i)
    B(3,k  ) = SH(3,i)/R
    B(3,k+1) = 0d0
    B(4,k  ) = SH(2,i)
    B(4,k+1) = SH(1,i)

    k = k + 2
enddo

return
end
!-------------------------------------------------------------------------------
subroutine GRAVLD (P,XL,GRAV,TH,WT,QUAD,ierr)

use modeldefinition, only: nee, igrav
USE CONSTANTS, only: G, twopi, zero, four
use modeltopology,   only: nen, ndof, nsd, nstr

!
! Computes the contribution to the load from gravity body forces
!
implicit none

integer, parameter :: NINT=4

!-pass
integer ierr
logical QUAD
double precision :: P(NEE),XL(NSD,NEN),GRAV(NDOF)
double precision :: TH, WT
!-locl
double precision :: SG(NEN),TG(NEN),XS(NSD,NSD),SH(NSTR,NEN)
!save TG,SG
double precision :: r
double precision :: GRAV1, GRAV2
double precision :: rc, zc
double precision :: DET
integer :: j
!-init

call SHAP20 (ZERO,ZERO,XL,DET,SH,XS,QUAD,ierr)

if (ierr.ne.0) then
    return
endif

! find center of element.
rc = ( XL(1,1) + XL(1,2) + XL(1,3) ) / 3.0
zc = ( XL(2,1) + XL(2,2) + XL(2,3) ) / 3.0


if (igrav.eq.1) then
    GRAV1 = GRAV(1) * 16.0d0 * DET * WT * TWOPI * rc
    GRAV2 = GRAV(2) * 16.0d0 * DET * WT * TWOPI * rc
else if (igrav.eq.2) then

    ! gravity of all elements points to the origin,
    ! with a magnitude of the first entry.
    ! assume constant gravity within an element. All nodes point in the same direction,
    ! while actually they should converge on the central point.

    ! both use grav(1), because there is just one magnitude
    GRAV1 = sqrt ( (GRAV(1)**2) / ( 1.0+(zc**2 / rc**2) )  )
    GRAV2 = sqrt ( (GRAV(1)**2) / ( 1.0+(rc**2 / zc**2) )  )

    ! add sign, based on whether points are left/right down/up from origin
    GRAV1 =    sign(GRAV1,-rc) * 16.0d0 * DET * WT * TWOPI * rc
    GRAV2 =    sign(GRAV2,-zc) * 16.0d0 * DET * WT * TWOPI * rc

else
    ! this is not supposed to happen
endif

do j=1,NEN
    P(2*j-1) = P(2*j-1) + GRAV1 * SH(3,j)
    P(2*j  ) = P(2*j  ) + GRAV2 * SH(3,j)
enddo

return
end
!-------------------------------------------------------------------------------
subroutine STIFF (DMAT,XL,S,TH,QUAD,ierr, nodeID)

!USE MODELDEFINITION
USE ALGEBRA,         only: IDFORM, CLEAR
USE CONSTANTS,       only: G, twopi, eps, zero
use modeltopology,   only: nen, nsd, ndof, nstr
use modeldefinition, only: nee
use iomodule,        only: stderr
!
! Computes the super-diagonal part of the local stiffness matrix.
!
implicit none
!    parameter (NEE=NEN*NDOF)
double precision,  parameter :: W=4d0


!-pass
integer          :: ierr
logical          :: QUAD
double precision :: XL(NSD,NEN),DMAT(NSTR,NSTR),S(NEE,NEE)
double precision :: th ! dummy
integer          :: nodeID
!-locl
double precision :: XS(NSD,NSD),SH(3,NEN),B(NSTR,NEE)

logical firstentry
save firstentry

double precision :: DB1, DB2, DB3, DB4
double precision :: r, c, det

integer          :: i, j

!-init
data firstentry /.true./
NEE=8
!
call SHAP20 (zero,zero,XL,DET,SH,XS,QUAD,ierr)
if (ierr.ne.0) then
    return
endif

 
R = SH(3,1)*XL(1,1) + & 
    SH(3,2)*XL(1,2) + &
    SH(3,3)*XL(1,3) + &
    SH(3,4)*XL(1,4)


if (R.lt.EPS) then
    write(stderr,10) EPS
10      format(///1x,'STIFF: R < ',1PE20.6)
    ierr = 1
    return
endif


call BMATRIX (B,SH,R)

C = DET * W * twopi * r


do j=1,NEE
!   Compute column J of DB-matrix
    DB1=DMAT(1,1)*B(1,j) + &
        DMAT(1,2)*B(2,j) + &
        DMAT(1,3)*B(3,j) + &
        DMAT(1,4)*B(4,j)

    DB2=DMAT(2,1)*B(1,j) + &
        DMAT(2,2)*B(2,j) + &
        DMAT(2,3)*B(3,j) + &
        DMAT(2,4)*B(4,j)

    DB3=DMAT(3,1)*B(1,j) + &
        DMAT(3,2)*B(2,j) + &
        DMAT(3,3)*B(3,j) + &
        DMAT(3,4)*B(4,j)

    DB4=DMAT(4,1)*B(1,j) + &
        DMAT(4,2)*B(2,j) + &
        DMAT(4,3)*B(3,j) + &
        DMAT(4,4)*B(4,j)


    do i=1,j
        S(i,j)=S(i,j)+C*(B(1,i)*DB1 + &
                         B(2,i)*DB2 + &
                         B(3,i)*DB3 + &
                         B(4,i)*DB4)
    enddo

enddo



return
end
!-------------------------------------------------------------------------------
subroutine ADDPR (B,X,D,IEN,LM,MAT,PRPMAT,IELNO,ISIDE,pressure, &
                  LMF,TFAULT,SKEW,DXE)

USE MODELDEFINITION, only: numpr, lgdef
USE ALGEBRA,         only: equate, vmprd, formrt
use constants,       only: half, one, useDXE
use constants,       only: eps, pi, twopi
use modeltopology,   only: nen, nsd, ndof
use debugmodule,     only: iecho
use iomodule,        only: stderr, stdout

! Subroutine to add pressures to load vector

 implicit none
!-pass
integer          :: ielno, mat, ien, iside, lm, lmf
double precision :: x, d, pressure, prpmat, b, tfault, skew, dxe
 dimension IELNO(*),MAT(*),IEN(NEN,*),X(NSD,*),D(NDOF,*), &
 ISIDE(*),pressure(*),PRPMAT(6,*),B(*),LM(NDOF,NEN,*), &
 LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
 DXE(NDOF,NEN,*)
!-locl
double precision :: XL(NSD,NEN),DL(NDOF,NEN),xn(NSD),tmp(NSD)
double precision :: rot(3,3)
integer          :: lien(2)
double precision :: forceVec(4)
double precision :: thisForce(2)

double precision :: w,r, angle, anglevec(2), area
double precision :: f(ndof)
double precision :: fl(ndof)

integer          :: offset
integer          :: globalnodenumber, globalelementnumber
integer          :: entrynumber

integer          :: ii, n, idof, nodecounter, numadd

double precision :: r1, r2, z1, z2

logical          :: quad

numadd = 0


do entrynumber = 1,IABS(NUMPR)

!	write(*,*) "------------------ Presssure entry ", entrynumber

    globalelementnumber = IELNO(entrynumber)

    if (globalelementnumber.le.0) then
        write(stderr,1) entrynumber,globalelementnumber
 1	format(1x,'ADDPR: fatal IELNO error. Entry ',I5,&
          ' has element number ',I8)
        call stoper()
    endif

    numadd = numadd + 1

    call LCOORD (X,XL,IEN(1,globalelementnumber),0)
    if (LGDEF.ne.0) then
       call LDISP  (DL,D,IEN(1,globalelementnumber),NDOF,NEN)
       call AddFaultDisplacement (DL,globalelementnumber,TFAULT,NDOF,NEN,42)
       call ADDSNE (DL,globalelementnumber,NDOF,NEN,useDXE)
       call REZONE(XL,DL,ONE)
    endif

!   fill lien with local side numbers for element side ISIDE(n)
    quad = IEN(3,globalelementnumber).ne.IEN(4,globalelementnumber)
    call SIDENP(ISIDE(entrynumber),quad,lien)

!   compute side length or surface area into variable AREA
    call FCGEOM (XL,lien,area,xn)

	! take average angle of the two nodes.
    angle =(-SKEW(1,IEN(lien(1),globalelementnumber))-&
             SKEW(1,IEN(lien(2),globalelementnumber)))*HALF

    anglevec(1) = angle
    anglevec(2) = 0d0

!    if (ABS(angle).gt.EPS) then
!        call FORMRT (anglevec,rot,NDOF)
!        call EQUATE (tmp,xn,NDOF)
!        call VMPRD (rot,tmp,xn,NDOF,3)
!    endif

    ! xn is the normal vector to the element side, pointing in the direction of the pressure

    ! the element side in the z,r system describes a (knotted) conical shell.
    ! the pressure works on the area of this shell.
    r1 = XL(1,lien(1))
    r2 = XL(1,lien(2))

    z1 = XL(2,lien(1))
    z2 = XL(2,lien(2))

!    write(*,*) "computing forces of pressure entry", entrynumber, &
!               "on element ", globalelementnumber, &
!               "having nodes: ",IEN(1:3,globalelementnumber), &
!               "with points ", lien(1), lien(2)
!    write(*,*) "and coords: ", 1, r1, z1
!    write(*,*) "and coords: ", 2, r2, z2

    ! compute forces corresponding to this pressure.
    ! takes into account that due to axisymmetry the edges are conical shells.
    ! Tested with a McTigue source (finite spherical pressure ball) and found OK.
    call forces(r1, z1, r2, z2, pressure(entrynumber), forceVec)

!	write(*,*) "gives unrotated force vec", 1, &
!                   IEN(lien(1),globalelementnumber), &
!                    forceVec(1:2)
!    write(*,*) "gives unrotated force vec", 2, &
!       	       	   IEN(lien(2),globalelementnumber), &
!       	       	    forceVec(3:4)



    do nodecounter = 1,2

!		write(*,*) 'processing node: ', nodecounter

        globalnodenumber = IEN(lien(nodecounter),globalelementnumber)
        anglevec(1) = -SKEW(1,globalnodenumber)
        anglevec(2) = 0d0

!		write(*,*) 'anglevec: ', anglevec

		! forceVec contains two vectors.
		! one for the first point on position 1 and 2
		! another for the seconde point on positions 3 and 4.
		! Put the proper one in vector thisForce:
		if (nodecounter.eq.1) then
			thisForce(1:2) = forceVec(1:2)
		else
            thisForce(1:2) = forceVec(3:4) 
		endif

!		write(*,*) 'unrotated thisForce: ', thisForce


        if (ABS(anglevec(1)).gt.EPS) then
            call FORMRT (anglevec,rot,NDOF)
            call EQUATE (tmp,thisForce,NDOF)
            call VMPRD (rot,tmp,fl,NDOF,3)
        else
            call EQUATE (fl,thisForce,NDOF)
        endif

!    	write(*,*) 'rotated thisForce: ', thisForce
!        write(*,*) 'rotated fg       : ', fl


#ifdef GIA
        do idof=2,2
#else
        do idof=1,NDOF
#endif
		    ii = LM(idof,lien(nodecounter),globalelementnumber)
		    if (ii.ne.0) then
		!        write(*,*) n,lien(1),"1 adding r pressure", forceVec(1)," to b entry", ii 

 !      	       	    write(*,*) "ADDPR adds force", &
  !                             globalelementnumber, &
!                               globalnodenumber, &
!                               nodecounter, &
!       	       	       	       ii, &
!       	       	       	       idof, fl(idof)



		        B(ii) = B(ii) + fl(idof)
		    endif
		enddo

	enddo

enddo

if (iecho.eq.1 .and. numadd.gt.0) then
    write(stdout,3) numadd
3   format(1x,I8,' pressure loads are being added')
endif


!
return
end

!-------------------------------------------------------------------------------
subroutine forces(r1, z1, r2, z2, pressure, fources)
! computes forces on two nodes, derived from pressure working on them.
! This takes into account that:
! 1: the exact area of the side forms a conical shell
! 2: the outer point gets more force, because it takes a bigger part of the surface

use constants, only: pi

implicit none

double precision            :: r1, z1, r2, z2, pressure
double precision            :: fources(4)  ! well, there are four of them, after all
logical                     :: echo

double precision            :: sideHeight, coneBaseRadius, coneTopRadius
double precision            :: surface, surfaceTop, surfaceFull
double precision            :: heightFull, heightTop
double precision            :: slantHeight, slantheightTop, slantHeightFull

double precision            :: normalLength, normalVecR, normalVecZ

double precision :: rs, rm, rb
double precision :: zs, zm, zb, dz
double precision :: ls, lm, lb
double precision :: as, am, ab
double precision :: surfaceSmall, surfaceLarge
double precision :: multiplier1, multiplier2

echo = .false.

if (echo) then
    write(*,*) "in forces with", r1, z1, r2, z2, pressure
endif

sideHeight = abs(z2 - z1)
coneBaseRadius = max(r1, r2)
coneTopRadius = min(r1, r2)


if (r1.eq.r2) then
    ! the surface is a not a conical shell but a cylindrical shell.
    ! Forces work only in the radial direction.
    multiplier1 = 0.5d0
    surface = 2d0 * pi * r1 * abs(z2 -z1)
else
    ! we have rs, rm, rb (radius small, medium and big)
    !         zs, zm, zb (height of the cones)
    !         ls, lmed, lb (lateral height of the cones)  (slant height, along the shell)
    !         as, am, ab (area of the conical shell)
    ! dz (delta z)
    rs = coneTopRadius
    rb = coneBaseRadius
    rm = 0.5d0 * (rs + rb)

    if (echo) then
        write(*,*) "rs, rm, rb", rs, rm, rb
    endif

    dz = sideHeight

    if (echo) then
        write(*,*) "dz: ", dz
    endif

    zs = rs * dz / (rb - rs)
!    zm = rs * dz * 0.5 / (rb - rm)
    zm = zs + 0.5 * dz
    zb = zs + dz

    if (echo) then
        write(*,*) "zs, zm, zb: ", zs, zm, zb
    endif

    ls = sqrt(zs**2 + rs**2)
    lm = sqrt(zm**2 + rm**2)
    lb = sqrt(zb**2 + rb**2)

if (echo) then
    write(*,*) "ls, lm, lb: ", ls, lm, lb
endif

    as = pi * ls * rs
    am = pi * lm * rm
    ab = pi * lb * rb

if (echo) then
    write(*,*) "as, am, ab: ", as, am, ab
endif

    surfaceSmall = am - as
    surfaceLarge = ab - am
    surface = ab - as

if (echo) then
    write(*,*) "surfaceSmall: ", surfaceSmall
    write(*,*) "surfaceLarge: ", surfaceLarge
endif

    if (r1.lt.r2) then
        ! r1 on the inside, so less force
        multiplier1 = surfaceSmall / (surfaceSmall + surfaceLarge)
    else
        ! r1 on the outside. so more force
        multiplier1 = surfaceLarge / (surfaceSmall + surfaceLarge)
    endif
endif

if (echo) then
    write(*,*) "multiplier 1 ", multiplier1
endif

multiplier2 = 1.0d0 - multiplier1

normalVecR = z2 - z1
normalVecZ = r2 - r1

normalLength = sqrt(normalVecR**2 + normalVecZ**2)

normalVecR = normalVecR / normalLength
normalVecZ = normalVecZ / normalLength

fources(1) =  normalVecR * pressure * surface * multiplier1
fources(2) = -normalVecZ * pressure * surface * multiplier1
fources(3) =  normalVecR * pressure * surface * multiplier2
fources(4) = -normalVecZ * pressure * surface * multiplier2


end subroutine



!-------------------------------------------------------------------------------
subroutine FCGEOM (XL,LIEN,AREA,XN)

use constants,     only: dfmin
use modeltopology, only: nen, nsd
use iomodule,      only: stderr
!
! Calculates area of element face and unit normal vector.
!
implicit none
!-pass
double precision :: XL(NSD,NEN),XN(NSD)
double precision :: area
integer          :: LIEN(2)

!-locl
XN(1) = XL(2,LIEN(2)) - XL(2,LIEN(1))
XN(2) = XL(1,LIEN(1)) - XL(1,LIEN(2))
AREA  = XN(1)*XN(1)+XN(2)*XN(2)

if (AREA.lt.dfmin) then
    write(stderr,1)
 1       format(1x,'FCGEOM: zero element face area')
    return
endif

AREA  = SQRT(AREA)
XN(1) = XN(1)/AREA
XN(2) = XN(2)/AREA

return
end
!-------------------------------------------------------------------------------
subroutine LFLUX (XL,TL,IEN,PRPTEM,FLUX,ierr)

use modeltopology, only: nsd, ndof, nen
!
! Routine to calculate element heat flux
!
 implicit none
!-pass
integer          :: IEN(NEN), ierr
double precision :: XL(NSD,NEN),TL(NEN),PRPTEM(6),FLUX(NSD)

!-locl

! todo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine EFORCE (EVP,P,XL,TH,QUAD,ierr)

USE MODELDEFINITION, only: IPOINT, nee
use constants,       only: twopi, eps
use modeltopology,   only: nsd, ndof, nen, nstr
use iomodule,        only: stderr

!
! Computes the effective forces at each node (modified 8/4/00)
!
implicit none

double precision, parameter :: GP=0.57735026918963d0

!-pass
integer ierr
logical QUAD
double precision :: XL(NSD,NEN),P(NEE),EVP(NSTR)
double precision :: th  ! dummy
!-locl
double precision :: SG(NEN),TG(NEN),XS(NSD,NSD),SH(3,NEN),B(NSTR,NEE)
! save SG,TG
integer          :: i
double precision :: c, det, w, r
!-init

double precision :: g
integer          :: nint, L

if (IPOINT.eq.4) then
    NINT = 4
    W    = 1D0
    G    = GP
else
    NINT = 1
    W    = 4D0
    G    = 0D0
endif

do L=1,NINT

    call SHAP20 (SG(L)*G,TG(L)*G,XL,DET,SH,XS,QUAD,ierr)

    if (ierr.ne.0) then
        return
    endif

    R = SH(3,1)*XL(1,1) + &
        SH(3,2)*XL(1,2) + &
        SH(3,3)*XL(1,3) + &
        SH(3,4)*XL(1,4)

    if (R.lt.EPS) then
        write(stderr,10) EPS
10      format(///1x,'EFORCE: R < ',1PE20.6)
        ierr = 1
        return
    endif

    call BMATRIX (B,SH,R)
    C = DET * W * twopi * r

    do i=1,8        
        P(i)=P(i) + C * (B(1,i)*EVP(1) + &
                         B(2,i)*EVP(2) + &
                         B(3,i)*EVP(3) + &
                         B(4,i)*EVP(4))
    enddo

enddo

return
end
!-------------------------------------------------------------------------------
subroutine DEVSTN (STN,SDEV)

USE MODELDEFINITION
use constants,       only: two, third
use modeltopology,   only: nstr

! Computes stress deviator SDEV

implicit none
! lwi be careful, newer version has NST3, which is incorrect. Should be 4

!-pass
double precision :: STN(NSTR),SDEV(NSTR)
!
SDEV(1)= (TWO*STN(1) - STN(2) - STN(3)) * THIRD
SDEV(2)= (TWO*STN(2) - STN(1) - STN(3)) * THIRD
SDEV(3)= (TWO*STN(3) - STN(1) - STN(2)) * THIRD
SDEV(4)= STN(4)

!
 return
 end
!-------------------------------------------------------------------------------
subroutine FORMBT (STN,BETA,BETB,PRPMAT,PRPLAS, elemID)

!USE MODELCTX
USE MODELDEFINITION,  only: icvis
use constants,        only: dfmin,dlogmin,dfmax,dlogmax,eps, third, one, two, six
use modeltopology,    only: nstr
use algebra,          only: idig

!
! Subroutine to compute strain rates from stress field. It defines BETA and
! evaluates the product BETA*STRESS. BETB is the viscoplastic strain rate.
!
implicit none
!-pass
double precision :: STN(NSTR),BETA(NSTR),BETB(NSTR),PRPMAT(6),PRPLAS(9)
integer          :: elemID
!-locl
double precision :: AVECT(NSTR)
double precision :: sigma, emhu, anpwr, streff, pstrs

!-init
sigma = 0d0
emhu  = PRPMAT(3)
anpwr = PRPMAT(4)

if (IDIG(ICVIS,4,2).ne.0) then
!   Calculate scalar stress PSTRS that will be compared to yield stress
    call YIELDS (STN,PRPLAS,PRPMAT,PSTRS, elemID)
!   Construct plastic flow vector AVECT
    call YIELDF (STN,PRPLAS,PRPMAT,AVECT, elemID)
!   Evaluate viscoplastic strain rate matrix BETB
    call FLOWVP (AVECT,STN,PRPLAS,PRPMAT,PSTRS,BETB, elemID)
endif


streff = ( STN(1)*STN(1) + &
           STN(2)*STN(2) + &
           STN(3)*STN(3) - &
           STN(1)*STN(2) - &
           STN(1)*STN(3) - &
           STN(2)*STN(3) )*THIRD + &
           STN(4)*STN(4)

if (streff.gt.dfmin) then
    streff = SQRT(streff)
    sigma = (streff/emhu)**(anpwr-ONE)/(SIX*emhu)
endif

! viscous stress-strain rate tensor
BETA(1) = sigma * (TWO*STN(1)-STN(2)-STN(3))
BETA(2) = sigma * (TWO*STN(2)-STN(1)-STN(3))
BETA(3) = sigma * (TWO*STN(3)-STN(1)-STN(2))
BETA(4) = sigma * STN(4) * SIX

if (IDIG(ICVIS,4,2).ne.0) then
    ! add visco-plastic contributions
    BETA(1) = BETA(1) + BETB(1)
    BETA(2) = BETA(2) + BETB(2)
    BETA(3) = BETA(3) + BETB(3)
    BETA(4) = BETA(4) + BETB(4)
endif

return
end subroutine
!-------------------------------------------------------------------------------
subroutine FORMJB (BETAJ,STN,PRPMAT,PRPLAS,BETB, elemID)

 USE MODELCTX
 USE MODELDEFINITION, only: icvis
 USE ALGEBRA,         only: CLEAR, idig
 use constants,       only: dfmin,dlogmin,dfmax,dlogmax,eps, one, two, six, third
use debugmodule,     only: delay
!
! Constructs the Jacobian of the stress-strain rate matrix BETA.
!
implicit none

!implicit none


integer, parameter :: NSTR=4
!-pass
double precision :: STN(NSTR),PRPLAS(9),PRPMAT(6),BETB(NSTR)
integer          :: elemID
!double precision :: BETAJ(6,6)
double precision :: BETAJ(NSTR,NSTR)
!-locl
double precision :: betajp(NSTR,NSTR)

double precision :: emhu, anpwr, avect
double precision :: a1, a2, a3
double precision :: strinv
double precision :: srr, szz, stt, sxy

integer          :: j
double precision :: sigma, streff, yield


!    compute effective stress
streff = ( STN(1)*STN(1) + &
           STN(2)*STN(2) + &
           STN(3)*STN(3) - &
           STN(1)*STN(2) - &
           STN(1)*STN(3) - &
           STN(2)*STN(3) )*THIRD + &
           STN(4)*STN(4) 

 if (streff.le.dfmin) then
BETAJ = 0d0
!call CLEAR (BETAJ,36,"BETAJ")
return
 endif



streff = SQRT(streff)


!   Plastic flow?
if (IDIG(ICVIS,4,2).ne.0) then
!        Calculate yield stress YIELD
    call YIELDS (STN,PRPLAS,prpmat,YIELD, elemID)
!        Construct Jacobian due to plastic yielding BETAJP
!!!!!
!!!!! do we not need a call to yieldf hgere?
!!!!!
!    call CLEAR (BETAJP,NSTR*NSTR)
    BETAJP = 0d0
! subroutine TANGVP (STN,PSTRS,PRPLAS,PRPMAT,BETB,AVECT,BETAJP)
    call TANGVP (STN,YIELD,PRPLAS,prpmat,betb,avect,BETAJP, elemID)

 endif
!

 emhu  = PRPMAT(3)
 anpwr = PRPMAT(4)
 a1 = anpwr-ONE
 a2 = a1/SIX
 a3 = a1*SIX
 strinv = ONE/streff
 srr = (TWO*STN(1)-STN(2)-STN(3))*strinv
 szz = (TWO*STN(2)-STN(1)-STN(3))*strinv
 stt = (TWO*STN(3)-STN(1)-STN(2))*strinv
 sxy = STN(4)*strinv
 BETAJ(1,1) =  TWO + a2*srr*srr
 BETAJ(1,2) = -ONE + a2*srr*szz
 BETAJ(1,3) = -ONE + a2*srr*stt
 BETAJ(1,4) =        a1*srr*sxy
 BETAJ(2,2) =  TWO + a2*szz*szz
 BETAJ(2,3) = -ONE + a2*szz*stt
 BETAJ(2,4) =        a1*szz*sxy
 BETAJ(3,3) =  TWO + a2*stt*stt
 BETAJ(3,4) =        a1*stt*sxy
 BETAJ(4,4) =  SIX + a3*sxy*sxy
 BETAJ(2,1) = BETAJ(1,2)
 BETAJ(3,1) = BETAJ(1,3)
 BETAJ(3,2) = BETAJ(2,3)
 BETAJ(4,1) = BETAJ(1,4)
 BETAJ(4,2) = BETAJ(2,4)
 BETAJ(4,3) = BETAJ(3,4)

 sigma=((streff/emhu)**a1)/(SIX*emhu)

 if (IDIG(ICVIS,4,2).ne.0) then
!        Add Jacobians of viscoplastic and viscous flow
! LWITODO? also apply loop to first index.
! is this correct?



do j=1,NSTR
    BETAJ(1,j) = sigma*BETAJ(1,j)+BETAJP(1,j)
    BETAJ(2,j) = sigma*BETAJ(2,j)+BETAJP(2,j)
    BETAJ(3,j) = sigma*BETAJ(3,j)+BETAJP(3,j)
    BETAJ(4,j) = sigma*BETAJ(4,j)+BETAJP(4,j)
enddo


else
! lwitodo same

do j=1,NSTR
    BETAJ(1,j) = sigma*BETAJ(1,j)
    BETAJ(2,j) = sigma*BETAJ(2,j)
    BETAJ(3,j) = sigma*BETAJ(3,j)
    BETAJ(4,j) = sigma*BETAJ(4,j)

enddo


 endif

 return
end subroutine
!-------------------------------------------------------------------------------
subroutine ADDSTR (B,X,D,IEN,LM,MAT,PRPMAT,IELSTR,ISSIDE, &
                   ISTR,LMF,TFAULT,SKEW,DXE)

USE MATERIALSMODULE,  only: LMAT
USE MODELdatamodule,  only: modeldatactx
USE MODELDEFINITION
USE ALGEBRA
use constants,        only: eps, twopi, half, useDXE
use debugmodule,      only: iecho
use modeltopology,    only: nsd, ndof, nstr, nen
use iomodule
use timestepmodule,   only: nstep

!
! Subroutine to add stresses to load vector. Routine variable GLOBAL controls
! whether stress boundary conditions are interpreted in a global or in a
! local coordinate frame which is defined by the normal to the plane.
!
 implicit none
!-pass
integer :: mat, ien, lm, ielstr, isside, istr, lmf
double precision :: x, d, prpmat, b, tfault, skew, dxe
 dimension MAT(*),IEN(NEN,*),X(NSD,*),D(NDOF,*),PRPMAT(6,*),B(*), &
   LM(NDOF,NEN,*),IELSTR(*),ISSIDE(*),ISTR(2,*), &
   LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN,*)
!-locl
double precision :: f(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),tmp(NSD)
double precision :: rot(3,3)
integer          :: lien(2)
 logical GLOBAL
 save GLOBAL
double precision, external :: FLOAT
double precision :: r,w, angle, anglevec(2)
double precision :: srr, szz, srz, stt, sn, st
integer          :: n, k, ii, numadd
integer          :: it0, it1
!-init
 data GLOBAL/.TRUE./
 numadd = 0
!

!write(*,*) "------------------------------------------------------------"
!write(*,*) "numstr", numstr, NUMSTRglobal, nstr
!write(*,*) "------------------------------------------------------------"


do n=1,NUMSTR

    k = IELSTR(n)
    if (k.le.0) then
        write(stderr,1) n
 1           format(1x,'ADDSTR: fatal IELSTR error. Index=',I5)
        call stoper()
    endif

    it0 = ISTR(1,n)
    it1 = ISTR(2,n)

!   add the load vector components when the increment time range has started...
    if (NSTEP.ge.it0 .and. NSTEP.le.it1) then

        numadd = numadd + 1
        call LCOORD (X,xl,IEN(1,k),0)

        if (LGDEF.ne.0) then
            call LDISP  (dl,D,IEN(1,k),NDOF,NEN)
            call AddFaultDisplacement (dl,k,TFAULT,NDOF,NEN,43)
            call ADDSNE (dl,k,NDOF,NEN,useDXE)
            call REZONE (xl,dl,1d0)
        endif

        call SIDENP(ISSIDE(n),IEN(3,k).ne.IEN(4,k),lien)


!        write(*,*) "read stress", modeldatactx%STRS(1:4,n)

        srr = modeldatactx%STRS(1,n) / FLOAT(it1-it0+1)
        szz = modeldatactx%STRS(2,n) / FLOAT(it1-it0+1)
!        srz = STRS(3,n) / FLOAT(it1-it0+1)
        stt = modeldatactx%STRS(3,n) / FLOAT(it1-it0+1)
        srz = modeldatactx%STRS(4,n) / FLOAT(it1-it0+1)

        xn(1) = XL(2,lien(2)) - XL(2,lien(1))
        xn(2) = XL(1,lien(1)) - XL(1,lien(2))

!            write(*,*) 'lien: ', lien(1), lien(2)
!            write(*,*) 'XL 1: ', XL(1,lien(1)), XL(1,lien(2))
!            write(*,*) 'XL 2: ', XL(2,lien(1)), XL(2,lien(2))

        angle = (-SKEW(1,IEN(lien(1),k))- &
                  SKEW(1,IEN(lien(2),k)))*HALF

        anglevec(1)    = angle
        anglevec(2)    = 0d0

        if (ABS(angle).gt.EPS) then
            call FORMRT (anglevec,rot,NDOF)
            call EQUATE (tmp,xn,NDOF)
            call VMPRD (rot,tmp,xn,NDOF,3)
        endif

        R = (XL(1,lien(1))+XL(1,lien(2)))*HALF
        w = TWOPI * R * HALF


!        write(*,*) "elem, side ID has normal", xn(1:2)
!        write(*,*) "and stresses rr, zz, tt, rz", srr, szz, stt, srz

        if (GLOBAL) then
            f(1) = -(xn(1)*srr + xn(2)*srz)*w   
            f(2) = -(xn(1)*srz + xn(2)*szz)*w
        else
            sn = srr
            st = srz
            f(1) = -(xn(1)*sn - xn(2)*st)*w
            f(2) = -(xn(1)*st + xn(2)*sn)*w
        endif

        ! first point
        ii = LM(1,lien(1),k)
        if (ii.gt.0) then
            B(ii) = B(ii) + f(1)
!            write(*,*) "added x", n, lien(1), f(1)
        endif

        ii = LM(2,lien(1),k)
        if (ii.gt.0) then
            B(ii) = B(ii) + f(2)
!            write(*,*) "added y", n, lien(1), f(2)
        endif

        ! second point
        ii = LM(1,lien(2),k)
        if (ii.gt.0) then
            B(ii) = B(ii) + f(1)
!            write(*,*) "added x", n, lien(2), f(1)
        endif

        ii = LM(2,lien(2),k)
        if (ii.gt.0) then
            B(ii) = B(ii) + f(2)
!            write(*,*) "added y", n, lien(2), f(2)
        endif

    endif
enddo

if (iecho.eq.1 .and. numadd.gt.0) then
    if (GLOBAL) then
        write(stdout,3)
 3           format(1x,'Global stress loads are being added')
    else
        write(stdout,4)
 4           format(1x,'Local stress loads are being added')
    endif
endif

return
end subroutine
!-------------------------------------------------------------------------------
subroutine ADDWNK (IWELM,IWSIDE,IWTIME,WPRES,IEN,MAT,X,D,PRPMAT, &
                   A,LM,IDIAG,LMF,TFAULT,SKEW,DXE,LMX,NSLIP)

USE MATERIALSMODULE,  only: LMAT
USE MODELDEFINITION,  only: numwnk, lgdef, numslp, neq
USE ALGEBRA
use constants,        only: eps, twopi, half, one, useDXE
use modeltopology,    only: nsd, ndof, nen
use petscksp,         only: add_values
use timestepmodule,   only: nstep
use modeldatamodule,  only: modeldatactx

#ifdef SPARSE
use spetscmodule,     only: m_a
#endif
use iomodule
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule,  only: elementSlipperyID
#endif
!
! Subroutine to add Winkler pressure contributions to the stiffness matrix
!
 implicit none
!-pass
integer          :: iwelm, iwside, iwtime, idiag(neq), lm, lmf, lmx, ien, nslip, mat
double precision :: wpres, x, prpmat, a, d, tfault, skew, dxe
 dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),A(*),LM(NDOF,NEN,*), &
   D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
   DXE(NDOF,NEN,*),LMX(NDOF,NEN,*),NSLIP(5,*)

#ifdef SPARSE
!#include "pde.i"
!#include "petsc/finclude/petscmatdef.h"
!#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petscvec.h"       
! for ADD_VALUES
!#include "petsc/finclude/petscsysdef.h"
!#include "petsc/finclude/petsckspdef.h"

   PetscMPIInt irank
#endif
!-locl
double precision :: w(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),rot(3,3),tmp(NSD)
integer          :: lien(2)
double precision :: angle, anglevec(2), area, r, wt
integer          :: n, LL, i, ii, j, ni, nj, nn, ierr
integer          :: idof, nodecounter
double precision :: r1, r2, z1, z2
double precision :: forceVec(4)
double precision :: thisForce(2)
double precision :: fl(ndof)
integer          :: globalnodenumber, globalelementnumber

#ifdef SPARSE
double precision :: v
#endif

if (NUMWNK.eq.0) then
    return
endif

do nn=1,IABS(NUMWNK)
    globalelementnumber = IABS(IWELM(nn))
	n = globalelementnumber
    LL= IWTIME(nn)


if (LL.lt.0) then
    if (LL.ne.-1.and.NSTEP.ge.-LL) cycle
elseif (LL.eq.0) then
    cycle
else
    if (NSTEP.lt.(LL-1)) cycle
endif

call LCOORD (X,XL,IEN(1,n),0)
 
if (LGDEF.ne.0) then
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
    call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,44)
    call ADDSNE (DL,n,NDOF,NEN,useDXE)
    call REZONE (XL,DL,ONE)
endif

call SIDENP(IABS(IWSIDE(nn)),IEN(3,n).ne.IEN(4,n),lien)
call FCGEOM (XL,lien,area,xn)

i = lien(1)
j = lien(2)

!write(*,*) "----------------- addwink wp to entry", nn, n, "nodes", i, j

angle = (-SKEW(1,IEN(i,n))-SKEW(1,IEN(j,n)))*HALF
anglevec(1)    = angle
anglevec(2)    = 0d0

!write(*,*) "angle: ", angle

!if (ABS(angle).gt.EPS) then
!    call FORMRT (anglevec,rot,NSD)
!    call EQUATE (tmp,xn,NSD)
!    call VMPRD (rot,tmp,xn,NSD,3)
!endif

r1 = XL(1,lien(1))
r2 = XL(1,lien(2))

z1 = XL(2,lien(1))
z2 = XL(2,lien(2))

! compute forces corresponding to this pressure.
! takes into account that due to axisymmetry the edges are conical shells.
! Particularly close to the symmetry axis this makes a noticable difference
call forces(r1, z1, r2, z2, WPRES(nn), forceVec)








!R = (XL(1,lien(1))+XL(1,lien(2)))*HALF
!wt = TWOPI * R * area * HALF
!w(1) = ABS(xn(1)) * WPRES(nn) * wt
!w(2) = ABS(xn(2)) * WPRES(nn) * wt

! nodes rotated, winkler pressure in y direction:
!w(1) = WPRES(nn) * wt
!w(2) = WPRES(nn) * wt

!write(*,*) "adding winklers: ", w

!-------------------------------------------------------------------
! explanation by Lukas
!
! For many models, the Winkler pressures are a mechanism to
! simulate the forces induced by contrasts in density.        
! On a round planet, these forces should only work in the radial direction.
! To accomplish this, nodes need to be rotated using an Euler angle,
! so that the y-axis is perpendicular to the surface.
! The magnitude provided in the input file is set to the 
! y-axis, and the x-axis will not provide Winkler contributions to the stiffness matrix.



do nodecounter=1,2


    globalnodenumber = IEN(lien(nodecounter),globalelementnumber)
    anglevec(1) = -SKEW(1,globalnodenumber)
    anglevec(2) = 0d0



!        write(*,*) 'processing node: ', nodecounter


    if (nodecounter.eq.1) then
        thisForce(1:2) = forceVec(1:2)
    else
        thisForce(1:2) = forceVec(3:4)
    endif

!	write(*,*) 'unrotated thisForce: ', thisForce


        if (ABS(anglevec(1)).gt.EPS) then
            call FORMRT (anglevec,rot,NDOF)
            call EQUATE (tmp,thisForce,NDOF)
            call VMPRD (rot,tmp,fl,NDOF,3)
        else
            call EQUATE (fl,thisForce,NDOF)
        endif

#ifdef GIA
            fl = abs(fl)
#endif


!	write(*,*) 'rotated thisForce: ', thisForce
!        write(*,*) 'rotated fg       : ', fl



#ifdef GIA
    do idof=2,2  ! assuming only the z direction has a non-zero Winkler contribution.
                 ! if the direction is different, Euler angles rotate the nodal point.
#else
    do idof=1,2
#endif

        ii = LM(idof,lien(nodecounter),n)

#ifdef SPARSE
        if (ii.gt.0) then

!		                    write(*,*) "ADDW adds force", &
!                               nn, idof, fl(idof)


            call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,fl(idof),ADD_VALUES, ierr)
            if (ierr.ne.0) then
                write(stderr,1) ierr
 1               format(1x,'ADDWNK: MatSetValues error ',I4)
                call stoper()
            endif
        endif
#else
        if (ii.gt.0) ii = IDIAG(ii)
        if (ii.gt.0) A(ii) = A(ii) + w(idof)
#endif

    enddo
enddo

!ii = LM(2,i,n)

!#ifdef SPARSE
!if (ii.gt.0) then
!call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,w(2),ADD_VALUES, ierr)
!        if (ierr.ne.0) then
!            write(stderr,1) ierr
!            call stoper()
!        endif
!endif
!#else
!if (ii.gt.0) ii = IDIAG(ii)
!if (ii.gt.0) A(ii) = A(ii) + w(2)
!#endif

!ii = LM(1,j,n)

!#ifdef SPARSE
!if (ii.gt.0) then
!call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,w(1),ADD_VALUES, ierr)
!        if (ierr.ne.0) then
!            write(stderr,1) ierr
!            call stoper()
!        endif
!endif
!#else
!if (ii.gt.0) ii = IDIAG(ii)
!if (ii.gt.0) A(ii) = A(ii) + w(1)
!#endif

!ii = LM(2,j,n)

!#ifdef SPARSE
!if (ii.gt.0) then
!call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,w(2),ADD_VALUES, ierr)
!        if (ierr.ne.0) then
!            write(stderr,1) ierr
!            call stoper()
!        endif
!endif
!#else
!if (ii.gt.0) ii = IDIAG(ii)
!if (ii.gt.0) A(ii) = A(ii) + w(2)
!#endif

if (NUMSLP.eq.0) cycle
if (IWELM(nn).lt.0) then
   ni = 1
   nj = 1
else
    ni = 0
    nj = 0
    do ii = 1,NUMSLP
    if (NSLIP(2,ii).eq.IEN(i,n)) ni = ni + 1
    if (NSLIP(2,ii).eq.IEN(j,n)) nj = nj + 1
    enddo
endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
if(elementSlipperyID(n,10).gt.0) then
    ii = IABS(LMX(1,i,elementSlipperyID(n,10)))
else
    ii = 0
endif
#else
ii = IABS(LMX(1,i,n))
#endif


if (ii.gt.0) then

#ifdef SPARSE
v = w(1)/DBLE(ni)
call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,v,ADD_VALUES, ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
            call stoper()
        endif
#else
ii = IDIAG(ii)
A(ii) = A(ii) + w(1)/DBLE(ni)
#endif

endif


#ifdef EXPERIMENTAL_ReducedSlipAlloc
if(elementSlipperyID(n,11).gt.0) then
    ii = IABS(LMX(1,j,elementSlipperyID(n,11)))
else
    ii = 0
endif
#else
ii = IABS(LMX(1,j,n))
#endif


if (ii.gt.0) then

#ifdef SPARSE
v = w(1)/DBLE(nj)
call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,v,ADD_VALUES, ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
            call stoper()
        endif
#else
ii = IDIAG(ii)
    A(ii) = A(ii) + w(1)/DBLE(nj)
#endif

endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
if(elementSlipperyID(n,12).gt.0) then
    ii = IABS(LMX(2,i,elementSlipperyID(n,12)))
else
    ii = 0
endif
#else
ii = IABS(LMX(2,i,n))
#endif


if (ii.gt.0) then

#ifdef SPARSE
v = w(2)/DBLE(ni)
call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,v,ADD_VALUES, ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
            call stoper()
        endif
#else
ii = IDIAG(ii)
A(ii) = A(ii) + w(2)/DBLE(ni)
#endif

endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
if(elementSlipperyID(n,13).gt.0) then
    ii = IABS(LMX(2,j,elementSlipperyID(n,13)))
else
    ii = 0
endif
#else
ii = IABS(LMX(2,j,n))
#endif


if (ii.gt.0) then

#ifdef SPARSE
    v = w(2)/DBLE(nj)
    call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,v,ADD_VALUES, ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
            call stoper()
        endif
#else
ii = IDIAG(ii)
A(ii) = A(ii) + w(2)/DBLE(nj)
#endif

endif
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine ELCRD (S,T,XELM,IEN,X,D,elemID,TFAULT,SKEW, &
                  FACTOR,ierr)

USE MODELDEFINITION
USE ALGEBRA,         only: CLEAR
use constants,       only: useDXE
use modeltopology,   only: nsd, ndof, nen
!
! Calculates (deformed) coordinates XELM at natural coordinates (R,S,T)
!
implicit none

!-pass
integer ierr
integer :: elemID
double precision :: factor

integer          :: ien
double precision :: s, t, xelm, x, d, tfault, skew

dimension        :: XELM(NSD),IEN(NEN),X(NSD,*),D(NDOF,*), &
                    TFAULT(NDOF,*), SKEW(2,*)
!-locl
double precision :: XL(NSD,NEN),DL(NDOF,NEN),SH(3,NEN),XS(NSD,NSD)
double precision :: det
logical          :: QUAD
integer          :: i,j
!-init
call CLEAR (XELM,NSD,"XELM")
!
QUAD = (IEN(3).ne.IEN(4))

call LCOORD (X,XL,IEN,0)

if (FACTOR.gt.0.) then
    call LDISP  (DL,D,IEN,NDOF,NEN)
    call AddFaultDisplacement (DL,elemID,TFAULT,NDOF,NEN,45)
    call ADDSNE (DL,elemID,NDOF,NEN,useDXE)
    call REZONE (XL,DL,FACTOR)
endif

call SHAP20 (S,T,XL,DET,SH,XS,QUAD,ierr)

if (ierr.ne.0) then
    return
endif

do j=1,NEN
    do i=1,NSD
        XELM(i) = XELM(i) + SH(3,j)*XL(i,j)
    enddo
enddo

return
end
!-------------------------------------------------------------------------------
subroutine ELCTR (XELM,IEN,X,D,elemID,TFAULT, &
 SKEW,FACTOR,ierr)

!      USE MODELDEFINITION
USE ALGEBRA,       only: CLEAR
USE CONSTANTS,     only: zero, deg2rad, third, useDXE
use modeltopology, only: nsd, ndof, nen
!
! Calculates element center coordinates XELM
!
implicit none
double precision, parameter :: THIRDM=-third
!-pass
integer :: ierr
integer :: elemID

integer          :: ien
double precision :: xelm, x, d, tfault, skew, factor
dimension XELM(NSD),IEN(NEN),X(NSD,*),D(NDOF,*), &
 TFAULT(NDOF,*), SKEW(2,*)
!-locl
double precision :: xl(NSD,NEN),dl(NDOF,NEN),sh(3,NEN),xs(NSD,NSD)
logical quad
double precision :: t, det
integer          :: j
!-init
call CLEAR (XELM,NSD,"XELM")
!
quad = (IEN(3).ne.IEN(4))
t = ZERO
if (.not.quad) then
    t = THIRDM
endif
!
call LCOORD (X,xl,IEN,0)

if (FACTOR.gt.0.) then
    call LDISP  (dl,D,IEN,NDOF,NEN)
    call AddFaultDisplacement (dl,elemID,TFAULT,NDOF,NEN,46)
    call ADDSNE (dl,elemID,NDOF,NEN,useDXE)

    do j=1,NEN
        XL(1,j) = XL(1,j) + FACTOR*dl(1,j)
        XL(2,j) = XL(2,j) + FACTOR*dl(2,j)
    enddo

endif
!
call SHAP20 (ZERO,t,XL,DET,SH,XS,quad,ierr)

if (ierr.ne.0) return
!
do j=1,NEN
    XELM(1) = XELM(1) + SH(3,j)*xl(1,j)
    XELM(2) = XELM(2) + SH(3,j)*xl(2,j)
enddo
!
return
end
!-------------------------------------------------------------------------------

subroutine TSTIFF (XL,QUAD,S,Conduction,TH,ierr)

! todo this function was a stub
! first attempt by direct copy from plnlib


USE CONSTANTS,      only: G, twopi
use modeltopology,  only: nsd
!
! computes the local conduction stiffness matrix from integration at Barlow
! points.
!
implicit none

integer, parameter :: NINT=4, NEN=4
double precision,  parameter :: W=4d0

!-pass
integer          :: ierr
logical          :: QUAD
double precision :: XL(2,NEN)
double precision :: S(NEN,NEN)
double precision :: Conduction(NSD)
double precision :: TH
!-locl
double precision :: sh(3,NEN)
double precision, save :: sg(NEN)
double precision, save :: tg(NEN)
double precision :: db(NSD,NEN)
double precision :: xs(NSD,NSD)
data tg/-1d0,-1d0, 1d0, 1d0/
data sg/-1d0, 1d0, 1d0,-1d0/

double precision :: detJ, dV, r, c
integer          :: i,j,k,l


do l=1,NINT
    call shap20 (sg(l)*G,tg(l)*G,XL,detJ,sh,xs,QUAD,ierr)

    r = (XL(1,1) + XL(1,2) + XL(1,3)) / 3d0
    C = detJ * W * twopi * r


    if (ierr.ne.0) then
        return
    endif

    dV = detJ * TH * C

    do j=1,NEN
        do i=1,NSD
            db(i,j) = Conduction(i) * sh(i,j) * dV
        enddo
    enddo

    do j=1,NEN
        do i=1,j
            do k=1,NSD
                S(i,j) = S(i,j) + sh(k,i) * db(k,j)
            enddo
        enddo
    enddo

enddo

return
end

!-------------------------------------------------------------------------------


subroutine ADVSTF (S,XL,VL,PL,COND,RHOCP,TH,HEAT,QUAD,ierr)
!
! convection (advection) update of local conduction stiffness matrix for
! steady convection-conduction solution
!

! todo: this sub was a stub.
! First try by Lukas, direct copy from pln

USE ALGEBRA
use constants, only: dfmin,dlogmin,dfmax,dlogmax,eps, &
                     zero, third, one, two, three
use debugmodule, only: iecho
use iomodule,        only: stderr
!
! NUMOPT = 0: no upwinding
! NUMOPT = 1: classical upwind scheme
! NUMOPT = 2: optimal upwind scheme
! NUMOPT = 3: double asymptotic approximation scheme
! NUMOPT = 4: critical approximation scheme
! For triangle elements, the Mizukami upwind scheme is used
!
implicit none
integer, parameter :: NUMOPT=2
integer,  parameter :: NEN=4, NSD=2, NINT=1
double precision, parameter :: G=ZERO, WGHT=4d0
!-pass
logical QUAD
integer ierr
double precision RHOCP,TH, heat
double precision :: XL(2,NEN),VL(2,NEN),S(NEN,NEN),COND(NSD),PL(NEN)
!-locl
double precision :: sh(3,NEN),sg(NEN),tg(NEN),xs(NSD,NSD),size(NSD)
double precision :: SM(NEN,NEN),SA(NEN,NEN),SD(NEN,NEN)
logical          :: firstentry
double precision :: vx, vy, vmag, xix, xiy, tau
integer          :: i, j, l
double precision :: alphax, alphay, cplump, detj, qi, samax, w
double precision, external :: apprx
save             :: tg,sg,firstentry
!-init
data tg/-1d0,-1d0, 1d0, 1d0/
data sg/-1d0, 1d0, 1d0,-1d0/
data firstentry /.true./

if (firstentry) then
    if (iecho.ne.0) then
        if (NUMOPT.eq.0) then
            write(stderr,5)
 5               format(1x,'No upwinding')
        else
            if (QUAD) then
                if (NUMOPT.eq.1) then
                    write(stderr,10)
10                       format(1x,'Upwind: classical upwind scheme')
                else if (NUMOPT.eq.2) then
                    write(stderr,20)
20                       format(1x,'Upwind: optimal upwind scheme')
                else if (NUMOPT.eq.3) then
                    write(stderr,30)
30                       format(1x,'Upwind: double assymptotic',1x, &
                     'approximation')
                else if (NUMOPT.eq.4) then
                    write(stderr,40)
40                       format(1x,'Upwind: critical approximation', &
                     ' scheme')
                else
                    write(stderr,50) NUMOPT
50                       format(1x,'ADVSTF: illegal NUMOPT value ', &
                     I2)
                    call stoper()
                endif
            else
                write(stderr,60)
60                   format(1x,'Upwind: Mizukami scheme')
            endif
        endif
    endif
    firstentry = .false.
endif

call CLEAR(SA,NEN*NEN,"sa")
size(1) = ZERO
size(2) = ZERO

if (NUMOPT.ne.0) then
    if (QUAD) then
        call ELMSIZE (XL,QUAD,size)
    else
        call CLEAR(SM,NEN*NEN,"sm")
        call CLEAR(SD,NEN*NEN,"sd")
!            consistent mass matrix SM
        call CPSTIFF (XL,QUAD,SM,RHOCP,TH,ierr)
!            conduction stiffness matrix SD
        call TSTIFF (XL,QUAD,SD,COND,TH,ierr)
    endif
endif
!
do l=1,NINT
    call SHAP20 (sg(l)*G,tg(l)*G,XL,detJ,sh,xs,QUAD,ierr)
    if (ierr.ne.0) then
        return
    endif

    Vx = VL(1,1)*sh(3,1) + &
         VL(1,2)*sh(3,2) + &
         VL(1,3)*sh(3,3) + &
         VL(1,4)*SH(3,4)

    Vy = VL(2,1)*sh(3,1) + &
         VL(2,2)*sh(3,2) + &
         VL(2,3)*sh(3,3) + &
         VL(2,4)*SH(3,4)

    Vmag = SQRT(Vx*Vx+Vy*Vy)

    if (Vmag.gt.1d-35) then
        if (QUAD .and. NUMOPT.ne.0) then
!            Directional grid Peclet numbers
            alphax = RHOCP*Vx*size(1)/(TWO*COND(1))
            alphay = RHOCP*Vy*size(2)/(TWO*COND(2))
            if (NUMOPT.eq.1) then
!                classical upwind
                xix = SIGN(ONE,alphax)
                xiy = SIGN(ONE,alphay)
            else if (NUMOPT.eq.2) then
!                optimal upwind
                xix = APPRX(alphax)
                xiy = APPRX(alphay)
            else if (NUMOPT.eq.3) then
!                double assymptotic approximation
                if (ABS(alphax).gt.THREE) then
                    xix = SIGN(ONE,alphax)
                else
                    xix = alphax*THIRD
                endif

                if (ABS(alphay).gt.THREE) then
                    xiy = SIGN(ONE,alphay)
                else
                    xiy = alphay*THIRD
                endif
            else if (NUMOPT.eq.4) then
!                critical approximation
                if (alphax.lt.-ONE) then
                    xix = -ONE-ONE/alphax
                elseif (alphax.gt.ONE) then
                    xix =  ONE-ONE/alphax
                else
                    xix = ZERO
                endif

                if (alphay.lt.-ONE) then
                    xiy = -ONE-ONE/alphay
                elseif (alphay.gt.ONE) then
                    xiy =  ONE-ONE/alphay
                else
                    xiy = ZERO
                endif
            endif
            tau=(size(1)*xix*Vx+size(2)*xiy*Vy)/(TWO*Vmag*Vmag)
        endif

        ! QUAD and NUMOPT>0
        w  = detJ * TH * WGHT
        do j=1,NEN
            do i=1,NEN
                SA(i,j) = &
                   RHOCP*w*sh(3,i)*(Vx*sh(1,j)+Vy*sh(2,j))
            enddo
        enddo

        do i=1,NEN
            do j=1,NEN
                S(i,j)  = S(i,j) + SA(i,j)
            enddo
        enddo

        if (NUMOPT.ne.0) then
            if (.not.QUAD) then
                call SYMMET(SM,NEN)
            endif
            do i=1,NEN
                if (.not.QUAD) then
                    if (i.lt.NEN) then
                        samax = ABS(SA(i,1))
                        do j=2,3
                            samax = MAX(ABS(SA(i,j)),samax)
                        enddo
!                    element ii of lumped mass matrix
                        cplump = SM(i,1)+SM(i,2)+SM(i,3)
                        tau = 5D-1*(cplump/(samax+SD(i,i)))
                    else
                    tau = ZERO
                    endif
                endif
                Qi = tau*(Vx*sh(1,i)+Vy*sh(2,i))
                PL(i) = PL(i) + HEAT*w*Qi
                do j=1,NEN
                    S(i,j) = S(i,j) + &
                       RHOCP*w*Qi*(Vx*sh(1,j)+Vy*sh(2,j))
                enddo
            enddo
        endif
! NUMOPT > 0
    endif
! Vmag > 0
enddo
!
return
end
!-------------------------------------------------------------------------------
double precision function APPRX(X)

use constants, only: half, one

implicit none
!-pass
double precision X
double precision a,b,y

!double precision, intrinsic :: exp

y = DABS(X)

if (y.le.1D-4) then
    APPRX = X*HALF
else if (y.ge.10.0) then
    APPRX = SIGN(ONE,X) - ONE/X
else
    a = EXP(X)
    b = ONE/a
    APPRX = (a+b)/(a-b) - ONE/X
endif

return
end

!-------------------------------------------------------------------------------
subroutine ELMSIZE (XL,QUAD,SIZEofElm)
!
! computes characteristic sizes of element

use modeltopology, only: nen, nsd
use constants,     only: zero, third, half, one, two
!
implicit none
!-pass
double precision :: XL(NSD,NEN),SIZEofElm(NSD)
logical QUAD
!-locl
double precision :: sh(NSD,NEN),xs(NSD,NSD),sa(NEN),ta(NEN)
!save sa,ta
integer :: i
!-init
!data sa/-0.5d0, 0.5d0,0.5d0,-0.5d0/
!data ta/-0.5d0,-0.5d0,0.5d0, 0.5d0/

sa= (/-0.5d0, 0.5d0,0.5d0,-0.5d0/) 
ta= (/-0.5d0,-0.5d0,0.5d0, 0.5d0/)

!
if (QUAD) then
    do i=1,NEN
        sh(1,i) = sa(i)*HALF
        sh(2,i) = ta(i)*HALF
    enddo
else
    do i=1,NEN
        sh(1,i) = sa(i)*(HALF-ta(i)*THIRD)
        sh(2,i) = ta(i)*HALF
    enddo
    sh(1,3) = sh(1,3)+sh(1,4)
    sh(2,3) = sh(2,3)+sh(2,4)
    sh(1,4) = ZERO
    sh(2,4) = ZERO
endif

xs(1,1) = sh(1,1)*XL(1,1) + &
          sh(1,2)*XL(1,2) + &
          sh(1,3)*XL(1,3) + &
          sh(1,4)*XL(1,4)

xs(2,1) = sh(2,1)*XL(1,1) + &
          sh(2,2)*XL(1,2) + &
          sh(2,3)*XL(1,3) + &
          sh(2,4)*XL(1,4)

xs(1,2) = sh(1,1)*XL(2,1) + &
          sh(1,2)*XL(2,2) + &
          sh(1,3)*XL(2,3) + &
          sh(1,4)*XL(2,4)

xs(2,2) = sh(2,1)*XL(2,1) + &
          sh(2,2)*XL(2,2) + &
          sh(2,3)*XL(2,3) + &
          sh(2,4)*XL(2,4)



SIZEofElm(1) = (ABS(xs(1,1))+ABS(xs(2,1)))*TWO
SIZEofElm(2) = (ABS(xs(1,2))+ABS(xs(2,2)))*TWO
!
return
end
!-------------------------------------------------------------------------------
subroutine LHEAT (XL,QUAD,PL,HEAT,TH,ierr)

use modeltopology, only: nen, nsd

!
! computes the local load vector contribution of heat production.
!
implicit none

integer,  parameter :: NINT=4
double precision,  parameter :: G=0.57735026918963d0
!-pass
double precision :: heat, th
integer ierr
logical QUAD
double precision :: XL(2,NEN),PL(NEN)

 return
 end subroutine
!-------------------------------------------------------------------------------
double precision function EigMax (STN,E,v,emhu,Npl,alpha)

! USE MODELDEFINITION
 use constants, only: dfmin,dlogmin,dfmax,dlogmax,eps, one, two, three, four, eight
!
! Gives an upper limit estimate of the maximum local eigenvalue
!
 implicit none
!-pass
! SNGL    real E,v,emhu,Npl,alpha
double precision :: E,v,emhu,Npl,alpha
double precision :: STN(4)


 EigMax = 0d0
 return
 end function
!-------------------------------------------------------------------------------
#ifdef EXPERIMENTALTRANSIENTDISCREEP
subroutine EffectiveDislCreepRate(InvEffViscDisl,STN,EffectiveRate)
! computes beta vector (viscous strain rate) for dislocation creep only
! return effective beta, i.e. square root of second invariant
! WARNING: the current equations apply for 3D, and should be rewritten for axisymmetric
! currently an error is returned

use debugmodule, only: xit
use iomodule, only: stderr
use modeltopology, only: nstr
use constants, only: six, third

implicit none

!-pass 
double precision :: STN(nstr) ! stress vector
double precision :: InvEffViscDisl ! inverse viscosity for dislocation creep
!-local
double precision :: Beta(nstr) ! viscous strain rate vector
double precision :: sigma,EffectiveRate

write(stderr,1) 
 1   format(1x,'WARNING: effective dislocation creep rate not implemented for axisymmetric')
call xit(1," ")

! inverse viscosity divided by 6
!sigma = InvEffViscDisl/(SIX)

!    viscous stress-strain rate vector
!BETA(1) = sigma * (TWO*STN(1)-STN(2)-STN(3))
!BETA(2) = sigma * (TWO*STN(2)-STN(1)-STN(3))
!BETA(3) = sigma * (TWO*STN(3)-STN(1)-STN(2))
!BETA(4) = sigma * STN(4) * SIX
!BETA(5) = sigma * STN(5) * SIX
!BETA(6) = sigma * STN(6) * SIX

! convert into effective strain rate
!EffectiveRate = sqrt( Beta(1)*Beta(1) + Beta(2)*Beta(2) + Beta(3)*Beta(3) &
!         - Beta(1)*Beta(2) - Beta(1)*Beta(3) - Beta(2)*Beta(3) )*THIRD &
!         + Beta(4)*Beta(4) + Beta(5)*Beta(5) + Beta(6)*Beta(6)

EffectiveRate=0
return
end subroutine
#endif
!-----------------------------------------------------------------------
 subroutine ADDFLX (BT,IEN,X,D,MAT,PRPMAT,IFLS,IFLX,BFLX,LMT, &
 LMTX,LMF,TFAULT,SKEW,DXE)!,NFLX)

USE MATERIALSMODULE, only : LMAT
USE MODELDEFINITION
use modeltopology,   only: nen, nsd, ndof
!
! Routine to add normal boundary fluxes to load vector
!
implicit none

integer, parameter :: NET=3
!-pass
!    integer NFLX

integer          :: ien, mat, ifls, iflx, lmt, lmf, lmtx
double precision :: bt, x, bflx, d, tfault, skew, dxe, prpmat
 dimension BT(*),IEN(NEN,*),X(NSD,*),MAT(*),PRPMAT(6,*),IFLS(*), &
 IFLX(*),BFLX(*),LMT(NEN,*),D(NDOF,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),LMTX(NEN,*),SKEW(2,*),DXE(NDOF,NEN,*)


 return
 end
!-------------------------------------------------------------------------------
subroutine TEMWNK (AA,LMT,ITDIAG,ITWINK,TWINK,IEN,MAT,PRPMAT, &
                   X,D,LMF,TFAULT,SKEW,DXE)

USE MATERIALSMODULE, only: LMAT
USE MODELDEFINITION
use constants,       only: half, one
use modeltopology,   only: nsd, ndof, nen
use timestepmodule,  only: nstep

!
! Subroutine to add Winkler fluxes to diagonal elements of thermal stiffness
! matrix
!
 implicit none
!-pass
integer          :: itwink, ien, mat, itdiag, lmt, lmf
double precision :: twink, x, prpmat, aa, d, tfault, skew, dxe
 dimension ITWINK(3,*),TWINK(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),AA(*),ITDIAG(*),LMT(NEN,*), &
   D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*), &
   DXE(NDOF,NEN,*)

if (NTWINK.eq.0) then
    return
endif

! todo fill empty sub
stop "TEMWNK not yet implemented in axi. Please contact model support. Leaving GTECTON..." 

return
end subroutine
!-------------------------------------------------------------------------------
subroutine THGRAD (S,T,XL,TL,GRAD,QUAD,ierr)

use constants,     only: zero, quart, one
use modeltopology, only: nsd, nen
!
! Routine for calculating thermal gradients at natural coordinates s,t
!
implicit none
!-pass
integer NPE
integer ierr
logical QUAD
double precision :: s, t
double precision :: XL(NSD,NEN),TL(NEN),GRAD(NSD)

NPE = NSD+1

stop "THGRAD not yet implemented in axi. Please contact model support. Leaving GTECTON..."

return
end
!-------------------------------------------------------------------------------
double precision function Tcentr (T,TX,XL,IEN,LMTX,ierr)



USE MODELDEFINITION
use constants,      only: zero, third
use iomodule,       only: stderr

implicit none
integer, parameter :: NSD=2, NDOF=2, NSTR=3, NEN=4
double precision, parameter :: sn=ZERO
!-pass
integer :: ierr
integer          :: IEN(NEN),LMTX(NEN)
double precision :: XL(NSD,NEN), T, TX
dimension :: T(*),TX(*)
!-locl
integer :: i
logical :: QUAD
double precision :: TL(NEN),SH(3,NEN),XS(NSD,NSD)
double precision :: detJ, Tc, tn


QUAD = (IEN(3).ne.IEN(4))
tn = ZERO

if (.not.QUAD) then
    tn = -third
endif

call LTEMP (IEN,T,TL,NEN)

if (NUMSLP+NUMFN.gt.0.and.IDIFT.eq.1) then
    call LTEMPX (LMTX,IEN,TX,TL,NEN)
endif

call SHAP20 (sn,tn,XL,detJ,SH,XS,QUAD,ierr)

if (ierr.ne.0) then
    write(stderr,1)
 1  format(1x,'Tcentr: error in calculation of shape functions')
    Tcentr = ZERO
    return
endif

Tc = ZERO

do i=1,NEN
    Tc = Tc + SH(3,i)*TL(i)
enddo

Tcentr = Tc

return
end
!-------------------------------------------------------------------------------
subroutine CPSTIFF (XL,QUAD,S,RhoCp,TH,ierr)
!
! computes the local heat capacity stiffness matrix, or consistent mass
! matrix, from integration at Barlow points.
!
use modeltopology, only: nen, nsd


implicit none

integer,  parameter :: NINT=4
double precision,  parameter :: G=0.57735026918963d0
!-pass
integer ierr
logical QUAD
double precision :: RhoCP, Th
double precision :: XL(2,NEN),S(NEN,NEN)



!-locl
double precision :: sh(3,NEN),xs(NSD,NSD)
double precision :: sg(4),tg(4)
save tg,sg
data tg/-1d0,-1d0, 1d0, 1d0/
data sg/-1d0, 1d0, 1d0,-1d0/
integer :: i,j,l
double precision :: dV, detJ

! todo copied from pln, may need axi modifications

do l=1,NINT
    ! shap20 should take the axisymmetric one...
    call shap20 (sg(l)*G,tg(l)*G,XL,detJ,sh,xs,QUAD,ierr)
    if (ierr.ne.0) then
        return
    endif
    dV = detJ * TH
    do j=1,NEN
        do i=1,j
            S(i,j) = S(i,j) + RhoCp * sh(3,i) * sh(3,j) * dV
        enddo
    enddo
enddo

return

end
!-------------------------------------------------------------------------------
double precision function AJ1 (ST,POIS)

!
! Calculates first invariant of stress tensor
!
implicit none
!-pass
double precision :: ST(3)
double precision :: pois

! todo sure this is correct?
! aj1 in plnlib takes the poisson ration into account...

! adjustements for nocompr?

AJ1 = (ST(1)+ST(2)+ST(3))/3d0

return
end function
!-------------------------------------------------------------------------------
double precision function AJ2 (ST,STLength,POIS)

USE MODELDEFINITION
use constants, only: third
!
! Calculates second invariant of deviatoric stress tensor
!
implicit none

integer, parameter :: NSTR=4
!-pass
integer          :: STLength
double precision :: ST(STLength)
double precision :: pois

! todo sure this is correct?
! aj1 in plnlib takes the poisson ration into account...


!write(*,*) "doing aj2"

AJ2 = (ST(1)*ST(1) + &
       ST(2)*ST(2) + &
       ST(3)*ST(3) - &
       ST(1)*ST(2) - &
       ST(1)*ST(3) - &
       ST(2)*ST(3))*THIRD + & ! - or + ?
       ST(4)*ST(4)

return
end
!-------------------------------------------------------------------------------
double precision function AJ3 (ST,POIS)

USE MODELDEFINITION
USE CONSTANTS, only: one, two, third
!
! Calculates third invariant of deviatoric stress tensor
!
implicit none

!-pass
double precision :: ST(3)
double precision :: pois
!-init
!

stop "AJ3 not yet implemented in axi. Please contact model support. Leaving GTECTON..."

! todo... ahum...

AJ3 = 0d0
!
return
end
!-------------------------------------------------------------------------------
 subroutine REZONE (XL,DL,FACTOR)

use modeltopology,   only: nen, ndof, nsd

!
! Routine to update the local coordinate XL by displacements DL for
! large deformation update
!
 implicit none

!-pass
 double precision ::  XL(NSD,NEN),DL(NDOF,NEN)
 double precision ::  FACTOR
!-locl
 integer :: i
!
 do i=1,NEN
XL(1,i) = XL(1,i) + DL(1,i)*FACTOR
XL(2,i) = XL(2,i) + DL(2,i)*FACTOR
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine TSTRN (EE,DT,VEXP)

! calculates strain due to thermal expansion

 implicit none

 integer, parameter :: NSTR=4
!-pass
 double precision :: EE(NSTR)
 double precision :: DT, VEXP

 EE(1) = DT*VEXP/3d0
 EE(2) = EE(1)
 EE(3) = EE(1)
 EE(4) = 0d0

 return
 end
!-------------------------------------------------------------------------------
subroutine YIELDS (STN,PRPLAS,PRPMAT,PSTRS, elemID)

USE CONSTANTS, only: root3, deg2rad, zero, half, one, two,  &
    three, dfmin,dlogmin,dfmax,dlogmax,eps

use modeltopology, only: nstr

! Routine to calculate the scalar stress PSTRS that will be compared to
! the (equivalent) uniaxial yield stress. See GTECTON notes, section
! on yield criteria.

implicit none
!-pass
double precision :: STN(NSTR),PRPLAS(9),PRPMAT(6)
double precision :: pstrs
integer          :: elemID
!-------------------------------------------------
! sub was an empty stub.
! first attempt by Lukas, direct copy from plnlib

!-locl
double precision :: i1
double precision, external :: EFFSTN,AJ1,AJ3

double precision :: stneff, sint3
double precision :: theta, sinphi, alpha, frangl
integer          :: nptype
!
nptype = INT(PRPLAS(1))

if (nptype.eq.0) then
    return
endif
!
!       compute (incompressible) effective stress and pressure
stneff = EFFSTN(STN,HALF,.false.)
i1 = AJ1(STN,HALF)

if (stneff.le.dfmin) then
   sint3 = ZERO
else
   sint3 = -THREE*root3*AJ3(STN,HALF)/(TWO*stneff**3)
   sint3 = MIN(ONE,sint3)
   sint3 = MAX(-ONE,sint3)
endif

theta = ASIN(sint3)/THREE
!
if (nptype.eq.1) then
!           Tresca
    PSTRS = TWO * stneff * cos(theta)
else if (nptype.eq.2) then
!           Von Mises
    PSTRS  = root3 * stneff
else if (nptype.eq.3) then
!           Mohr-Coulomb
    frangl = PRPLAS(2)
    sinphi = SIN(frangl*deg2rad)
    PSTRS  = ABS(i1) * sinphi / THREE &
           + stneff * (COS(theta) - ((SIN(theta)*sinphi)/root3))
else if (nptype.eq.4) then
!           Drucker-Prager
    frangl = PRPLAS(2)
    sinphi = SIN(frangl*deg2rad)
!           Outer apices:
    alpha = (TWO * sinphi)/(root3 * (THREE - sinphi))
!           Inner apices:
!           alpha = (TWO * sinphi)/(root3 * (THREE + sinphi))
    PSTRS = alpha * i1 + stneff
endif


return
end
!-------------------------------------------------------------------------------
subroutine YIELDF (STN,PRPLAS,PRPMAT,AVECT, elemID)
USE MODELDEFINITION
USE CONSTANTS, only: ZERO,HALF,ONE,TWO,THREE, root3, deg2rad, dfmin,dlogmin,dfmax,dlogmax,eps
!
!       Compute visco-plastic flow vector AVECT (see GTECTON notes).
!
implicit none
!
integer NSTR
parameter (NSTR=4)
!-pass
double precision STN,PRPLAS,AVECT,PRPMAT
dimension STN(NSTR),PRPLAS(9),AVECT(NSTR),PRPMAT(6)
integer          :: elemID

!-------------------------------------------------
! sub was an empty stub.
! first attempt by Lukas, direct copy from plnlib
! changed
! sxx -> srr
! syy -> szz
! szz -> stt
! stt = -sp  ->   stt = STN(3) - sp



!-locl
integer i,nptype
double precision stneff,EFFSTN,i2,sp,srr,szz,stt,veca1(NSTR), &
 veca2(NSTR),veca3(NSTR),const1,const2,const3,sint3,theta,cost, &
 cos3t,sint,tant,tan3t,frangl,sinp,AJ1,AJ3,alpha
external EFFSTN,AJ1,AJ3
!


! from subroutine FORMJB.
! Useful here?

! emhu  = PRPMAT(3)
! anpwr = PRPMAT(4)
! a1 = anpwr-ONE
! a2 = a1/SIX
! a3 = a1*SIX
! strinv = ONE/streff
! srr = (TWO*STN(1)-STN(2)-STN(3))*strinv
! szz = (TWO*STN(2)-STN(1)-STN(3))*strinv
! stt = (TWO*STN(3)-STN(1)-STN(2))*strinv
! sxy = STN(4)*strinv


sint = 0d0

do i=1,NSTR
    AVECT(i) = ZERO
enddo
!
nptype = INT(PRPLAS(1))
if (nptype.eq.0) then
    return
endif
!
stneff = EFFSTN(STN,HALF,.false.)
if (stneff.le.dfmin) then
    return
endif
i2 = stneff*stneff
!
!       calculate incompressible stress deviator
sp = AJ1(STN,HALF)/THREE
srr = STN(1) - sp
szz = STN(2) - sp
!
!       set up vectors A1, A2 and A3
if (IOPT.eq.0) then
!           plane strain
    veca1(1) = 1.5D0
    veca1(2) = 1.5D0
    veca1(3) = ZERO
    veca3(1) = ZERO
    veca3(2) = ZERO
    veca3(3) = ZERO
else
!           plane stress
    veca1(1) = ONE
    veca1(2) = ONE
    veca1(3) = ZERO
    stt = STN(3) - sp
    veca3(1) = szz*stt + i2/THREE
    veca3(2) = srr*stt + i2/THREE
!    veca3(3) = -TWO * stt * STN(3)
    veca3(3) = -TWO * stt * STN(4)

endif

veca2(1) = srr/(TWO*stneff)
veca2(2) = szz/(TWO*stneff)
! in pln this would be shear stress...
!veca2(3) = STN(3)/stneff
veca2(3) = STN(4)/stneff



!
if (nptype.eq.1) then
!           Tresca
    const1 = ZERO
    sint3 = -THREE*root3*AJ3(STN,HALF)/(TWO*(stneff**THREE))
    sint3 = MIN(ONE,sint3)
    sint3 = MAX(-ONE,sint3)
    theta = ASIN(sint3)/THREE
    if (ABS(theta/deg2rad-30.).le.1.) then
        const2 = root3
        const3 = ZERO
    else
        cost  = COS(theta)
        cos3t = COS(THREE*theta)
        sint  = SIN(theta)
        tant  = TAN(theta)
        tan3t = TAN(THREE*theta)
        const2= TWO*cost*(ONE+(tant*tan3t))
        const3= (root3*sint)/(i2*cos3t)
    endif
else if (nptype.eq.2) then
!           Von Mises
    const1 = ZERO
    const2 = root3
    const3 = ZERO
else if (nptype.eq.3) then
!           Mohr-Coulomb
    frangl = PRPLAS(2)
    sinp   = SIN(frangl*deg2rad)
    const1 = sinp/THREE
    sint3 = -THREE*root3*AJ3(STN,HALF)/(TWO*(stneff**THREE))
    sint3 = MIN(ONE,sint3)
    sint3 = MAX(-ONE,sint3)
    theta = ASIN(sint3)/THREE
    if (ABS(theta/deg2rad-30.).le.1.) then
        const2 = (root3+(SIGN(ONE,-theta)*(sint/ROOT3)))/TWO
        const3 = ZERO
    else
        cost  = COS(theta)
        cos3t = COS(THREE*theta)
        sint  = SIN(theta)
        tant  = TAN(theta)
        tan3t = TAN(THREE*theta)
        const2 = cost*(ONE+tant*tan3t+sint*(tan3t-tant)/root3)
        if (ABS(i2*cos3t).gt.dfmin) then
            const3 = (root3*sint + cost*sinp) /(TWO*i2*cos3t)
        else
            const3 = ZERO
        endif
    endif
else if (nptype.eq.4) then
!           Drucker-Prager
    frangl = PRPLAS(2)
    sinp   = SIN(frangl*deg2rad)
!           Outer apices:
    alpha = (TWO * sinp)/(root3 * (THREE - sinp))
!           Inner apices:
!           alpha = (TWO * sinp)/(root3 * (THREE + sinp))
    const1 = alpha
    const2 = ONE
    const3 = ZERO
else
    write(*,*) "Axilib subroutine YIELDF encountered and error"
    write(*,*) "Unknown nptype: ", nptype, " but has to be 2, 3 or 4"
    stop "Please contact model support"
endif
!
!       calculate A-vector
do i=1,NSTR
    AVECT(i) = const1*veca1(i)+const2*veca2(i)+const3*veca3(i)
enddo
!

return
end
!-------------------------------------------------------------------------------
!call TANGVP (STN,YIELD,PRPLAS,BETAJP)
subroutine TANGVP (STN,PSTRS,PRPLAS,PRPMAT,BETB,AVECT,BETAJP, elemID)
   USE MODELDEFINITION
   USE ALGEBRA, only: CLEAR
   USE CONSTANTS, only: ZERO,HALF,ONE,TWO,THREE,FOUR,root3, &
    THIRD,SIX,deg2rad,dfmin,dlogmin,dfmax,dlogmax,eps

!
! Compute BETAJP, which is the Jacobian matrix for viscoplastic flow
! See GTECTON notes for details.
!
implicit none

integer NSTR
parameter (NSTR=3)

!-pass
double precision :: STN(NSTR),BETAJP(NSTR,NSTR),PRPLAS(9)
double precision :: PRPMAT(6),BETB(NSTR),AVECT,PSTRS
integer          :: elemID

return
end
!-------------------------------------------------------------------------------
subroutine INYLD (STN,PRPLAS,PRPMAT,PSTRS,BETB,FNORM,UNIYLD)
USE CONSTANTS, only: ONE,TWO,THREE,SIX,root3,deg2rad,HALF,dfmin,dlogmin,dfmax,dlogmax,eps
!
! Computes FNORM = (F - F_0)/F_0 and uniaxial yield stress UNIYLD
! (see GTECTON notes). PSTRS is the scalar stress for specific yield
! criteria that is compared to UNIYLD
!

!-------------------------------------------------
! sub was an empty stub.
! first attempt by Lukas, direct copy from plnlib
! Is there nothing for ntype 2, von Mises?

implicit none

integer, parameter :: NSTR=4
!-pass
double precision STN(NSTR),PRPLAS(9),BETB(NSTR),PRPMAT(6), &
                 PSTRS,FNORM,UNIYLD
!-locl
integer nptype
double precision allow,frangl,stnhrd,debar,stneff,EFFSTN
external EFFSTN

!-init
FNORM  = -1d0
nptype = INT(PRPLAS(1))

if (nptype.eq.0) return

stneff = EFFSTN(STN,HALF,.false.)

if (stneff.le.dfmin) return

allow  = 1d-2
frangl = PRPLAS(2) * deg2rad
UNIYLD = PRPLAS(3)
stnhrd = PRPLAS(7)

if (nptype.eq.3) then
!           Mohr-Coulomb
    UNIYLD = UNIYLD * COS(frangl)
else if (nptype.eq.4) then
!           Drucker-Prager outer apices
    UNIYLD = SIX*UNIYLD*COS(frangl)/(root3*(THREE-SIN(frangl)))
!           Drucker-Prager inner apices
!           UNIYLD = SIX*UNIYLD*COS(frangl)/(root3*(THREE+SIN(frangl)))
endif

!       effective viscoplastic strain rate for hardening
if (stnhrd.gt.dfmin) then
    debar = BETB(1)**TWO + BETB(2)**TWO + HALF*BETB(3)**TWO
    debar = TWO * debar / THREE
    if (debar.gt.dfmin) then
        debar = SQRT(debar)
    endif
    UNIYLD = UNIYLD + stnhrd*debar
endif

UNIYLD = MAX(UNIYLD,1d-3)
FNORM = PSTRS/UNIYLD - ONE

if (FNORM.le.0d0) then
    FNORM = -1d0
endif

return
end
!-------------------------------------------------------------------------------
subroutine FLOWVP (AVECT,STN,PRPLAS,PRPMAT,PSTRS,BETB, elemID)

USE ALGEBRA,       only: CLEAR
use modeltopology, only: nstr
use constants,     only: one

!
! Computes BETB by multiplying AVECT with fluidity and PHI
!
implicit none
!-pass
double precision :: pstrs
double precision :: AVECT(NSTR),STN(NSTR),PRPLAS(9),BETB(NSTR),PRPMAT(6)
integer          :: elemID

!-------------------------------------------------
! sub was an empty stub.
! first attempt by Lukas, direct copy from plnlib
! which obviously failed...
!-------------------------------------------------

!local
double precision :: fluidy, yldpwr, cmult, fnorm, uniyld
integer          :: nflow, i
!-init
fluidy = PRPLAS(4)
yldpwr = PRPLAS(5)
nflow  = INT(PRPLAS(6))
!
!       Compute FNORM
call INYLD (STN,PRPLAS,PRPMAT,PSTRS,BETB,FNORM,uniyld)
call CLEAR (BETB,NSTR,"BETB")
!
if (INT(PRPLAS(1)).eq.0) then
    return
endif

if (FNORM.le.0d0) then
    return
endif
!
if (NFLOW.eq.1) then
    cmult = fluidy * FNORM ** yldpwr
else
    cmult = fluidy * ( EXP(yldpwr*FNORM) - ONE )
endif
!
do i=1,NSTR
    BETB(i) = cmult * AVECT(i)
enddo



!
return
end
!-------------------------------------------------------------------------------
subroutine WINKLP (BTOT,DELD,IWELM,IWSIDE,IWTIME,WPRES,IEN, &
                   MAT,X,D,PRPMAT,LM,LMF,TFAULT,SKEW, &
                   DXE,LMX,NSLIP,DX)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
USE MODELDEFINITION, only: numwnk, lgdef, numslp, &
                           nElemsWithSlipperyNodes
use modeldatamodule,   only: elementSlipperyID
#else
USE MODELDEFINITION, only: numwnk, lgdef, numslp
#endif

USE ALGEBRA,         only: equate, vmprd, formrt
use constants,       only: eps, twopi, half, one, useDXE
use meshdatamodule,  only: meshdatactx
use modeltopology,   only: nsd, ndof, nen
use timestepmodule,  only: nstep


!! Subroutine to compute Winkler restoring pressures from displacements
! and add them to the global force vector BTOT
!
 implicit none
!-pass

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer          :: LMX(NDOF,NEN,nElemsWithSlipperyNodes)
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
integer          :: LMX(NDOF,NEN,meshdatactx%neglobal)
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif


integer          :: iwelm, iwside, iwtime, ien, mat, lm, lmf, nslip
double precision :: wpres, x, prpmat, btot, deld, d, tfault, skew, dx
dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),LM(NDOF,NEN,*),BTOT(*), &
   DELD(NDOF,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),NSLIP(5,*), DX(NDOF,*)

!-locl
double precision :: w(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),rot(3,3)
double precision :: tmp(NSD)
integer :: lien(2)
double precision :: r, wt, angle, anglevec(2), area
integer :: slipElemID
integer :: n, LL, i, j, ni, nj, nn, ii


!
if (NUMWNK.eq.0) then
    return
endif
!
do nn=1,IABS(NUMWNK)
    n = IABS(IWELM(nn))
    LL= IWTIME(nn)

    if (LL.lt.0) then
        if (LL.ne.-1.and.NSTEP.ge.-LL) cycle
    elseif (LL.eq.0) then
        cycle
    else
        if (NSTEP.lt.(LL-1)) cycle
    endif
call LCOORD (X,XL,IEN(1,n),0)
if (LGDEF.ne.0) then
call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,47)
    call ADDSNE (DL,n,NDOF,NEN,useDXE)
call REZONE (XL,DL,ONE)
endif
call SIDENP(IABS(IWSIDE(nn)),IEN(3,n).ne.IEN(4,n),lien)
    call FCGEOM (XL,lien,area,xn)
i = lien(1)
j = lien(2)
    angle = (-SKEW(1,IEN(i,n))-SKEW(1,IEN(j,n)))*HALF
    anglevec(1)    = angle
    anglevec(2)    = 0d0


    if (ABS(angle).gt.EPS) then
        call FORMRT (anglevec,rot,NSD)
        call EQUATE (tmp,xn,NSD)
        call VMPRD (rot,tmp,xn,NSD,3)
    endif
    R = (XL(1,i)+XL(1,j))*HALF
    wt = TWOPI * R * area * HALF
    w(1) = ABS(xn(1)) * WPRES(nn) * wt
    w(2) = ABS(xn(2)) * WPRES(nn) * wt
    ii = LM(1,i,n)
    if (ii.gt.0) BTOT(ii) = BTOT(ii) - w(1)*DELD(1,IEN(i,n))
    ii = LM(2,i,n)
    if (ii.gt.0) BTOT(ii) = BTOT(ii) - w(2)*DELD(2,IEN(i,n))
    ii = LM(1,j,n)
    if (ii.gt.0) BTOT(ii) = BTOT(ii) - w(1)*DELD(1,IEN(j,n))
    ii = LM(2,j,n)
    if (ii.gt.0) BTOT(ii) = BTOT(ii) - w(2)*DELD(2,IEN(j,n))
  300       if (NUMSLP.eq.0) cycle
if (IWELM(nn).lt.0) then
    ni = 1
    nj = 1
else
        ni = 0
        nj = 0
        do ii = 1,NUMSLP
            if (NSLIP(2,ii).eq.IEN(i,n)) ni = ni + 1
            if (NSLIP(2,ii).eq.IEN(j,n)) nj = nj + 1
        enddo
endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc

        slipElemID = elementSlipperyID(n,14)

        if(slipElemID.gt.0) then
            ii = IABS(LMX(1,i,slipElemID))

            if (ii.gt.0) then
                BTOT(ii) = BTOT(ii) - w(1)*DX(1,IEN(i,slipElemID))/DBLE(ni)
            endif

            ii = IABS(LMX(2,i,slipElemID))

            if (ii.gt.0) then
                BTOT(ii) = BTOT(ii) - w(2)*DX(2,IEN(i,slipElemID))/DBLE(ni)
            endif

            ii = IABS(LMX(1,j,slipElemID))

            if (ii.gt.0) then
                BTOT(ii) = BTOT(ii) - w(1)*DX(1,IEN(j,slipElemID))/DBLE(nj)
            endif

            ii = IABS(LMX(2,j,slipElemID))

            if (ii.gt.0) then
                BTOT(ii) = BTOT(ii) - w(2)*DX(2,IEN(j,slipElemID))/DBLE(nj)
            endif

        endif
#else

    ii = IABS(LMX(1,i,n))
    if (ii.gt.0) BTOT(ii) = BTOT(ii) - &
       w(1)*DX(1,IEN(i,n))/DBLE(ni)
    ii = IABS(LMX(2,i,n))
    if (ii.gt.0) BTOT(ii) = BTOT(ii) - &
       w(2)*DX(2,IEN(i,n))/DBLE(ni)
    ii = IABS(LMX(1,j,n))
    if (ii.gt.0) BTOT(ii) = BTOT(ii) - &
       w(1)*DX(1,IEN(j,n))/DBLE(nj)
    ii = IABS(LMX(2,j,n))
    if (ii.gt.0) BTOT(ii) = BTOT(ii) - &
       w(2)*DX(2,IEN(j,n))/DBLE(nj)


#endif

enddo

!
return
end subroutine
!-------------------------------------------------------------------------------
 subroutine SHAREL (IEN,IJSUR,NUMEL,NEL)

!     Finds elements neighbouring to element NEL.
!    Neighbouring elements are identified by two shared nodes.

!use modeldefinition, only: nsur
use modeltopology,   only: nen
use debugmodule,     only: xit
use iomodule,        only: stderr

implicit none

integer, parameter :: NSUR=4
!-pass
integer NEL,NUMEL
integer :: ien, ijsur
dimension IEN(NEN,*),IJSUR(NSUR,*)

!-locl
integer, external ::  NINDX
integer :: n,m,nn,mm,nshare,isur, ielm, ii


nn   = NEN
!    set to last IJSUR(isur,NEL) <> 0 (last index that was written)
isur = NINDX (IJSUR,NEL)

if (IEN(3,NEL).eq.IEN(4,NEL)) then
    nn = 3
endif
!    loop other elements
do ielm=NEL+1,NUMEL
!        compare nodal points of element NEL with those of IELM > NEL
    mm = NEN
    nshare = 0
    if (IEN(3,ielm).eq.IEN(4,ielm)) then
        mm = 3
    endif

    do n=1,nn
        do m=1,mm
            if (IEN(n,NEL).eq.IEN(m,ielm)) then
                nshare=nshare+1
            endif
        enddo
    enddo

    if (nshare.gt.2) then
        write(stderr,1) NEL,ielm,nshare
 1           format(1x,'sharel: elements ',I5,' and ',I5,' share ', &
         I1,' nodal points')
        call xit(1," ")
    else if (nshare.eq.2) then
        isur = isur + 1
        if (isur.le.NSUR) then
            IJSUR(isur,NEL) = ielm
!            copy info NEL into IJSUR(*,ielm)
!            find free index
            ii = NINDX (IJSUR,ielm) + 1
            IJSUR(ii,ielm) = NEL
        else
            write(stderr,2) NEL,NSUR
 2               format(1x,'sharel: element ',I5, &
             ' surrounded by more than ',I1,' elements')
            call xit(1," ")
        endif

        if (isur.eq.NSUR) then
            goto 500
        endif
    endif
enddo

500    return
end subroutine
!-------------------------------------------------------------------------------
integer function NINDX (IJSUR,NEL)

use debugmodule,     only: xit
use iomodule,        only: stderr

implicit none

integer, parameter :: NSUR=4
!-pass
integer NEL
integer :: ijsur
dimension IJSUR(NSUR,*)
!-local
integer :: i

!    find free index
 i = 1
100    if (i.le.NSUR) then
if (IJSUR(i,NEL).ne.0) then
i = i + 1
goto 100
else
i = i - 1
endif
 else
write(stderr,1) NEL,NSUR
 1       format(1x,'NINDX: element ',I5,' surrounded by more than ', &
     I5,' elements')
write(stderr,2) (IJSUR(i,NEL),i=1,NSUR)
 2       format(1x,8(I5,1x))
call xit(1," ")
 endif
 NINDX = i
 return
 end
!-------------------------------------------------------------------------------
 subroutine FLUID (BETB,PRPMAT,PRPLAS)

!USE TIMESTEPMODULE
USE CONSTANTS, only: one, two, three, four, six, &
    root3, deg2rad, half,dfmin,dlogmin,dfmax,dlogmax,eps
use modeltopology, only: nstr

!    Routine to determine visco-plastic fluidity from
!    explicit time stepping limit (Cormeau) (LINEAR only)
!     (Ref.: Owen & Hinton, sec. 8.9)
!    ONLY WORKS WITH NEWTONIAN RELAXATION!

implicit none
!-pass
double precision :: PRPMAT(6),PRPLAS(9),BETB(NSTR)

 return
 end
!-------------------------------------------------------------------------------
subroutine BPMATRIX (XL,DL,QUAD,BP,ierr)

use algebra,   only: dot
use constants, only: zero, third
!
! Sets up derivatives matrix BP(i,j) = dUi/dXj
!
implicit none

integer, parameter :: NDOF=2, NSD=2, NEN=4
!-pass
logical QUAD
double precision :: XL(NSD,NEN),DL(NDOF,NEN),BP(NDOF,NSD)
integer :: ierr
!-locl
double precision :: SH(3,NEN),XS(NSD,NSD)
double precision :: det

integer          :: i, j

call SHAP20 (ZERO,ZERO,XL,DET,SH,XS,QUAD,ierr)

if (IERR.ne.0) then
    write(0,*) "SHAP20 returned with error value:", ierr
    stop "Leaving GTECTON..."
endif

do i=1,NDOF
    do j=1,NSD
        BP(i,j) = dot(SH(j,:), DL(i,:),nen)
    enddo
enddo

return
end
!-------------------------------------------------------------------------------
subroutine sdimnt (ISURF,SURF,B,ID,GRAV,PRPSED,DELTP)

USE MODELDEFINITION, only: nsurf, nsed
use modeltopology,   only: nsd, ndof

! Routine to calculate load increment resulting from sedimentary
! transport in following time step. The surface uplift SURF is
! updated. It is assumed that the first derivative of the surface
! elevation at the first and last points is zero. Other assumption
! is that the NDOF spatial component is vertical.

implicit none
!-pass
double precision :: surf, b, grav, prpsed
integer :: isurf, id
dimension :: ISURF(*),SURF(NSD,*),B(*),ID(NDOF,*),GRAV(*), PRPSED(2)

double precision :: deltp

if (NSURF.le.0 .or. NSED.eq.0) then
    return
endif


stop "SDIMNT not yet implemented in axi. Please contact model support. Leaving GTECTON..."

return
end
!-----------------------------------------------------------------------
subroutine SIDENP (ISS,QUAD,LIEN)

use debugmodule,     only: xit
use iomodule,        only: stderr
!
! returns local node point numbers in LIEN, belonging to side IS
! nodal points in LIEN are positively oriented.
!
! LWI: replaced IS by ISS, because IS has become a reserved word by PETSc (Index Set).
!-pass
 integer ISS,LIEN(2)
 logical QUAD

!
if (ISS.eq.1 .or. ISS.eq.2) then
    LIEN(1) = ISS
    LIEN(2) = ISS + 1
else if (ISS.eq.3 .or. ISS.eq.4) then
    if (QUAD) then
        LIEN(1) = ISS
        LIEN(2) = ISS + 1
        if (LIEN(2).gt.4) then
            LIEN(2) = 1
        endif
    else
        LIEN(1) = 3
        LIEN(2) = 1
    endif
else
    write(stderr,1) ISS
 1       format(1x,'SIDENP: quadrilaterals have 4 sides, not ',I2)
    call xit(1," ")
endif
!
return
end
!-----------------------------------------------------------------------
subroutine NPSIDE (NP,IEN,ISINT)

! returns side number 

integer, parameter:: NEN=4
integer :: NP(3),IEN(NEN),ISINT
if (NP(1).eq.IEN(1)) then
    if (NP(2).eq.IEN(2)) then
        ISINT = 1
    elseif (NP(2).eq.IEN(3)) then
        ISINT = 3
    endif
elseif (NP(1).eq.IEN(2)) then
    if (NP(2).eq.IEN(1)) then
        ISINT = 1
    elseif (NP(2).eq.IEN(3)) then
        ISINT = 2
    endif
elseif (NP(1).eq.IEN(3)) then
    if (NP(2).eq.IEN(1)) then
        ISINT = 3
    elseif (NP(2).eq.IEN(2)) then
        ISINT = 2
    endif
endif

end subroutine
!-----------------------------------------------------------------------
double precision function SRFHGT (SURF,NSURF,X)
!
! Determines surface height at X from surface nodal point data by search and
! interpolate.
!
USE ALGEBRA
use modeltopology, only: nsd
use iomodule

implicit none

!-pass
integer :: nsurf
double precision :: surf, x
dimension SURF(NSD,*)

#include "coord.i"
!-locl
double precision :: dx1, dx2, a, b
integer          :: i, k, m

 logical FIRST
 save FIRST
 data FIRST /.true./
srfhgt = 0d0
!
if (NSURF.le.0) return
!
 if (FIRST) then
!        Check whether coordinates are uniformly in- or decreasing
if (NSURF.gt.2) then
    dx2 = SURF(IX,2)-SURF(IX,1)
    do i=2,NSURF-1
    dx1 = dx2
    dx2 = SURF(IX,i+1)-SURF(IX,i)
    if (dx1*dx2.le.0.) then
        write(stderr,1)
 1                   format(1x,'srfhgt: non-uniform input coords')
    do k=1,NSURF
            write(stderr,2) k,SURF(1,k),SURF(2,k)
 2                       format(1x,'SURF(',I3,')= ',2(1PG14.6))
    enddo
        call stoper()
    endif
    enddo
endif
FIRST = .false.
 endif
!
i = 1
100     if (i.ge.NSURF) then
    write(stderr,3) X
 3       format(1x,'x = ',1PG14.5,' does not occur in surface data')
    return
endif
if ((X-SURF(IX,i))*(X-SURF(IX,i+1)).le.0.) goto 200
i = i + 1
goto 100

200     call linlin(SURF(IX,i),SURF(IY,i),SURF(IX,i+1),SURF(IY,i+1),a,b,m)

if (m.eq.0) then
    srfhgt = a*X+b
else
!           X = B
    write(stderr,4) X
 4       format(1x,'surface height indeterminate at x = ',1PG14.5)
    return
endif
!
return
end function

!-----------------------------------------------------------------------

subroutine SLipperySKEW (IDSLP,SKEW,NSELSD,IEN,X,D,TFAULT,DXE,LMF)

USE MODELDEFINITION, only: lgdef, nslskew, numrot, numsn
USE CONSTANTS,       only: zero, pi, rad2deg
use debugmodule,     only: debug, iecho
use modeltopology,   only: nen, ndof, nsd
use iomodule
use timestepmodule,  only: nstep

!
! Updates skew angles for slippery nodes
!
 implicit none
!-pass
integer :: idslp, nselsd, ien, lmf
double precision :: skew, x, d, tfault, dxe
 dimension IDSLP(*),SKEW(2,*),NSELSD(2,*),IEN(NEN,*),X(NSD,*), &
   D(NDOF,*),TFAULT(NDOF,*),DXE(NDOF,NEN,*),LMF(NDOF,NEN,*)

!-locl
integer :: nwait, lutmp
logical firstentry
save firstentry,lutmp,nwait

integer :: lien(2)
integer :: k, kk, n, nn
integer :: nels
double precision :: a, angle, oldskew
integer, external :: nextlu

!-init
data firstentry,nwait /.true.,50/
!
if (LGDEF.eq.0 .or. NSLSKEW.eq.0 .or. NUMROT.eq.0) then
    return
endif
!
if (firstentry) then
    if (iecho.ne.0) then
        write(stdout,*) 'Updating Euler angle'
    endif
    lutmp = nextlu(10)
    call openf(lutmp,'slskew.tmp','unknown')
    firstentry = .false.
endif
!
do kk=1,NUMSN
    k = IDSLP(kk)
    if (SKEW(1,k).eq.ZERO) then
        cycle
    endif

    angle = ZERO
    nels  = 0

    do nn=1,NSLSKEW
        n = NSELSD(1,nn)
        call SIDENP (NSELSD(2,nn),IEN(3,n).ne.IEN(4,n),lien)
        if (IEN(lien(1),n).eq.k .or. IEN(lien(2),n).eq.k) then
            call TANGENT (X,D,TFAULT,DXE(1,1,n), &
                           IEN(1,n),n,SKEW,lien,a)
            nels = nels + 1
            if (nels.gt.1 .and. ABS(angle-a).gt.0.9*PI) then
                if (angle.gt.a) then
                    a = a + PI
                else
                    a = a - PI
                endif
            endif
            angle = (angle*DBLE(nels-1)+a)/DBLE(nels)
        endif
    enddo

    if (nels.eq.0) then
        cycle
    endif

    if (DEBUG) then
        oldskew = SKEW(1,k)
    endif

    SKEW(1,k) = angle
!
    if (DEBUG) then
        if (NSTEP.eq.0) then
            write(stderr,1) k,SKEW(1,k)*rad2DEG
 1               format(1x,'SLSKEW: SKEW(',I5,') initialized at ', &
             F16.12,' deg')
        else
            write(stderr,2) k,oldskew*rad2DEG,SKEW(1,k)*rad2DEG
 2               format(1x,'SLSKEW: SKEW(',I5,') updated from ', &
             F16.12,' to ',F16.12,' deg')
        endif
    endif

    if (MOD(NSTEP,nwait).eq.0) then
        write(lutmp,3) NSTEP,k, SKEW(1,k)*rad2DEG
 3       format(1x,'NSTEP = ',I6,'   SKEW(1,',I6,') = ',F16.12, &
     ' deg')
    endif
enddo
!
if (MOD(NSTEP,nwait).eq.0) then
    call flushf(lutmp)
endif
!
return
end
!-----------------------------------------------------------------------
   subroutine TANGENT (X,D,TFAULT,DXE,IEN,elemID,SKEW,LIEN,ANGLE)

   USE MODELDEFINITION
   USE CONSTANTS, only: halfpi, one,eps, useDXE
use modeltopology,   only: nen, ndof, nsd
use iomodule, only: stderr

!
! returns angle of tangent of deformed element side
!
 implicit none
!-pass
 integer :: elemID

integer :: ien
double precision :: x, d, skew, tfault, dxe
 dimension IEN(*),X(NSD,*),D(NDOF,*),SKEW(2,*), &
  TFAULT(NDOF,*), DXE(NDOF,*)

integer :: LIEN(2)
double precision :: angle, area


!-locl
double precision :: xl(NSD,NEN),dl(NDOF,NEN),xn(NSD)
!
!       localize coordinate data
call LCOORD (X,xl,IEN,0)
!
!       deformed coordinates XL
call LDISP  (dl,D,IEN,NDOF,NEN)
call AddFaultDisplacement (dl,elemID,TFAULT,NDOF,NEN,48)
 call ADDSNE (dl,elemID,NDOF,NEN,useDXE)
call REZONE (xl,dl,ONE)
 call FCGEOM (xl,LIEN,area,xn)
!
 if (ABS(xn(2)).gt.EPS) then
!        returns angle in range -pi/2 to pi/2
    ANGLE = ATAN(xn(1)/xn(2))
 else
if (xn(1).gt.EPS) then
ANGLE = HALFPI
else if (xn(1).lt.-EPS) then
ANGLE = -HALFPI
else
write(stderr,1)
 1           format(1x,'TANGENT: collapsing nodes error')
call stoper()
endif
endif
!
return
end
!-----------------------------------------------------------------------
subroutine LCOORD (X,XL,IEN,elemID)

use modeltopology,   only: nen, nsd

! Program to localize coordinate data
!
implicit none
!-pass
integer :: ien
double precision :: x, xl
dimension IEN(NEN),X(NSD,*),XL(NSD,NEN)

integer :: elemID

XL(1,1) = X(1,IEN(1))
XL(2,1) = X(2,IEN(1))
XL(1,2) = X(1,IEN(2))
XL(2,2) = X(2,IEN(2))
XL(1,3) = X(1,IEN(3))
XL(2,3) = X(2,IEN(3))
XL(1,4) = X(1,IEN(4))
XL(2,4) = X(2,IEN(4))


return
end subroutine
!-----------------------------------------------------------------------

double precision function dilate (strain,POIS)

USE MODELDEFINITION, only: iopt, NOCOMPR
use constants,       only: third
!
! computes volume change or volume change rate, depending on what S is.
!
implicit none
!-pass
double precision :: strain(3)
double precision :: pois

if (NOCOMPR.ne.0) then
   dilate = 0d0
else
   dilate = (strain(1)+strain(2)+strain(3))*THIRD
endif

return
end function
!-------------------------------------------------------------------------------
double precision function EFFEE (S,POIS,SIGNINFO)

USE MODELDEFINITION
use constants, only: dfmin,eps,zero, half
!
! Compute effective deviatoric strain
!
implicit none
!-pass
logical          :: SIGNINFO
double precision :: S(3)
double precision :: pois
!-locl
double precision :: J2
double precision, external   :: AJ2 ! external function in this file
!-init
effee = ZERO

J2 = AJ2(S,3,POIS)

if (J2.gt.dfmin) then
    effee = SQRT(J2)
    if (SIGNINFO .and. (S(1)+S(2)+S(3)).lt.-EPS) then
        effee = -effee
    endif
endif
!
return
end function
!-------------------------------------------------------------------------------
subroutine STRDEV (S,POIS)

use constants, only: one
!
! computes strain (rate) deviator
!
implicit none

!-pass
double precision :: S(3)
!-locl
double precision :: dil
double precision, external :: dilate
double precision :: pois

dil = dilate(S,ONE)

S(1) = S(1)-dil
S(2) = S(2)-dil
S(3) = S(3)-dil  ! perhaps S(3) = S(2) - dil

!
return
end subroutine
!-------------------------------------------------------------------------------
double precision function EFFSTN (S,POIS,SIGNINFO)

use constants, only: dfmin,eps, zero

!
! computes effective deviatoric stress. S is the full stress tensor.
!
implicit none

integer, parameter :: NSTR=4
!-pass
logical            :: SIGNINFO
double precision   :: S(NSTR)
!local
double precision, external   :: AJ2
double precision   :: J2
double precision   :: pois

J2 = AJ2(S,4,POIS)
if (J2.gt.dfmin) then
    EFFSTN = SQRT(J2)
    if (SIGNINFO .and. (S(1)+S(2)+S(3)).lt.-EPS) then
        effstn = -effstn
    endif
else
    EFFSTN = ZERO
endif
!
return
end function
!-------------------------------------------------------------------------------
subroutine STNTEN (STN,POIS,SL)

!USE MODELDEFINITION
USE ALGEBRA, only: CLEAR
!
implicit none

integer, parameter :: NSTR=4
!-pass
double precision :: STN(NSTR)
double precision :: SL(3,3)
double precision :: pois
!
call CLEAR(SL,9,"SL")

SL(1,1) = STN(1)
SL(2,2) = STN(2)
SL(3,3) = STN(3)
SL(1,2) = STN(4)
SL(2,1) = STN(4)

!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine STNDEV (S,POIS)

USE MODELDEFINITION
use constants, only: third

!
implicit none
!-pass
double precision :: S(3)
double precision :: pois
!-locl
double precision :: p


p = (S(1)+S(2)+S(3))*THIRD

S(1) = S(1) - p
S(2) = S(2) - p
S(3) = S(3) - p


return
end subroutine
!-------------------------------------------------------------------------------
double precision function PRESSURE (S,POIS)

USE MODELDEFINITION
use constants, only: third
!
implicit none
!-pass
double precision :: S(3)
double precision :: pois
!-locl

pressure = (S(1)+S(2)+S(3))*THIRD

return
end function
!-------------------------------------------------------------------------------

#ifdef SHEARZONE_TRACTIONS
subroutine SHEARTRAX (IFLAG,NSTEP)
implicit none
integer :: iflag, nstep
end subroutine
#endif

!-------------------------------------------------------------------------------
                                                                                          
