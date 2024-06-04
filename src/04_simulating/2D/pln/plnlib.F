!    These routines are the dimensionality-dependent portion of GTECTON
!    The routines here are for plane strain and plane stress displacements
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
!            |             |               4 |     5/6     | 2
!            |             |                 |             |
!            |             |                 |             |
!            |             |                 |             |
!            o-------------o                 o-------------o
!            1             2                        1
!
!
! Collapse of nodal points 4 onto 3 (i.e. assigning identical global node
! numbers to nodes 3 and 4) gives a triangle-element. Side 5 refers to the top
! face and side 6 to the bottom face of the element.
!-------------------------------------------------------------------------------
!    THE ROUTINES LISTED HERE ARE:
!
!    NAME                    FUNCTION
!
!
!    ADDFLX       NORMAL BOUNDARY FLUX CONTRIBUTIONS TO LOAD VECTOR
!    ADDPR        TRANSLATES APPLIED PRESSURES INTO LOAD VECTOR
!    ADDSTR       STRESS LOAD COMPUTATION
!    ADDWNK       ADDS WINKLER PRESSURES TO STIFFNESS MATRIX
!    AJ1          FIRST STRESS TENSOR INVARIANT
!    AJ2          SECOND STRESS TENSOR INVARIANT
!    AJ3          THIRD STRESS TENSOR INVARIANT
!    BDELD        STRAIN FROM DISPLACEMENTS
!    BMATRIX      STRAIN-DISPLACEMENT MATRIX
!    CPSTIFF      LOCAL HEAT CAPACITY STIFFNESS MATRIX
!    DEVSTN       STRESS DEVIATOR CALCULATION
!    EFORCE       EFFECTIVE LOAD VECTOR FOR VISCOUS FLOW
!    ELCRD        CALCULATES COORDINATES OF DEFORMED ELEMENT
!    ELCTR        CALCULATES CENTER COORDINATES OF DEFORMED ELEMENT
!    EIGMAX       CALCULATES MAXIMUM EIGENVALUE
!    FCGEOM       COMPUTES FACE NORMAL VECTOR AND FACE AREA
!    FLOWVP       VISCOPLASTIC STRAIN RATE
!    FORMBT       VISCOUS FLOW CORRECTION FOR STRAIN
!    FORMJB       JACOBEAN MATRIX FOR VISCOUS FLOW
!    GRAVLD       GRAVITATIONAL BODY FORCE COMPUTATION
!    INTCON       DIMENSIONAL AND TOPOLOGIC DATA
!    YIELDS       COMPUTES SCALAR STRESS FOR COMPARISON WITH YIELD STRESS
!    LCOORD       LOCALIZES COORDINATES
!    LFLUX        ELEMENT HEAT FLOW
!    LHEAT        HEAT PRODUCTION CONTRIBUTION TO LOAD VECTOR
!    FLUID        DETERMINES VISCOPLASTIC FLUIDITY
!    MATERL       MATERIAL MATRICES
!    REZONE       LOCAL REZONE FOR LARGE DEFORMATION
!    SDIMNT       SURFACE ELEVATION AND LOADS FROM SEDIMENTARY TRANSPORT
!    SRFHGT       DETERMINES SURFACE HIGHT AT A GIVEN LOCATION
!    SHAP20       SHAPE FUNCTION FOR 2-D ELEMENTS
!    SHAREL       IDENTIFIES NEIGHBOURING ELEMENTS
!    SHEARTRAX    APPLIES VELOCITY-DEPENDENT SHEAR-ZONE TRACTION (NEEDS SHEARZONE_TRACTIONS FLAG)
!    SIDENP       RETURNS LOCAL NODE NUMBERS ON A SPECIFIC SIDE
!    SHLCRD       COMPUTES COORDINATES IN LOCAL (SPHERICAL) FRAME
!    STIFF        STIFFNESS MATRIX
!    TANGVP       VISCOPLASTIC JACOBIAN MATRIX
!    TCENTR       COMPUTES ELEMENT CENTER TEMPERATURES
!    THGRAD       THERMAL GRADIENTS
!    TSTIFF       LOCAL CONDUCTION STIFFNESS MATRIX
!    WINKLP       WINKLER PRESSURES CONTRIBUTION TO GLOBAL FORCE VECTOR
!    WAPPLY       APPLIES RESIDUAL WINKLER FORCES
!    YIELDF       CONSTRUCTS PLASTIC FLOW VECTOR
!
!**********************************************************************

subroutine INTCON ()

USE LABELMODULE,     only : LABELC, &
                            LABELD, &
                            LABELS, &
                            LABELR, &
                            MATLBL, &
                            CRPLBL, &
                            TEMLBL
USE MODELDEFINITION, only : NTYPE, NSUR, NEE, NEC, simVariety, &
                            allTemperaturesFixed
use modeltopology,   only : nSidesPerElem, NSD,NSTR,NEN, NDOF
use constants,       only : isPln

! Set variables that depend on topology

implicit none

simVariety = isPln

NTYPE = 1
NSD   = 2
NDOF  = 2
nSidesPerElem = 3
NSTR  = 3
NEN   = 4
NSUR  = 4
NEC   = NSD*NEN
NEE   = NDOF*NEN
!      call labelsinit('pln')
LABELC(1)='  X '
LABELC(2)='  Y '
LABELD(1)=LABELC(1)
LABELD(2)=LABELC(2)
LABELS(1)='SXX '
LABELS(2)='SYY '
LABELS(3)='SXY '
! LABELS(4)='SXZ ' ERROR INTRODUCED BY?
! LABELS(5)='SYZ ' ERROR INTRODUCED BY?
LABELS(6)='    '
LABELR(1)='RXX '
LABELR(2)='RYY '
LABELR(3)='RXY '
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
TEMLBL(4)='X-CONDUCT'
TEMLBL(5)='Y-CONDUCT'

allTemperaturesFixed = .false.

return
end
!-------------------------------------------------------------------------------
SUBROUTINE SHAP20 (S,T,X,DET,SH,XS,QUAD,ierr)

use constants,     only: zero, half
use iomodule,        only: stderr

implicit none
! Program to compute shape functions for quadrilateral or triangle

!    S,T,                = natural coordinates
!    SH(1,NEN),SH(2,NEN) = x and y derivatives of shape functions
!    SH(3,NEN)           = shape functions
!    XS(NSD,NSD)         = inverse jacobian matrix
!    DET                 = Jacobian matrix determinant
!    X(2,NEN)            = global coordinates

integer, parameter :: NEN=4, NSD=2

!-pass
logical QUAD
integer ierr
double precision :: SH(3,NEN),X(NSD,NEN),XS(NSD,NSD)
double precision :: S, T
double precision :: DET
!-locl
double precision :: TEMP
integer          :: i,j,k
double precision :: SA(NEN),TA(NEN)
save SA,TA
!-init

! s contains x coordinate of natural points, 
! t contains y coordinates.
! the points descrive the square with corners:
! (-1/2, -1/2)  (1/2, -1/2)   (1/2, 1/2)  and (-1/2, 1/2)
! for triangles, the last two points merge to
! (0, 1/2)
!
! see equation manual for more details

data SA/-0.5d0, 0.5d0,0.5d0,-0.5d0/
data TA/-0.5d0,-0.5d0,0.5d0, 0.5d0/

ierr = 0

!    Calculate shape functions and s- and t-derivatives of shape functions
!    these derivatives are local.
do i=1,NEN ! loop over the nodes of the element
    SH(1,i) = SA(i)*(half+TA(i)*T)
    SH(2,i) = TA(i)*(half+SA(i)*S)
    SH(3,i) = (half+SA(i)*S)*(half+TA(i)*T)
enddo

if (.not.QUAD) then
!   adjust shape function of third node for triangular elements
    do i=1,3
        SH(i,3) = SH(i,3)+SH(i,4)
        SH(i,4) = zero
    enddo
endif

! Jacobian matrix calculation
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
    write(stderr,*) "Element has inverted itself!"
    write(stderr,*) "Shape function fails! Determinant is: ",DET
    write(stderr,*) "due to matrix entries: "
    write(stderr,*) XS(1,1), XS(1,2)
    write(stderr,*) XS(2,1), XS(2,2)
    write(stderr,*) "Likely cause: too much deformation. "
    ierr = 1
    return
endif

!    Calculatate inverse Jacobian matrix (Cramers rule)
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
use constants,       only: eps, zero, half
use debugmodule,     only: iecho
use modelctx,        only: getrank
use iomodule,        only: stderr
!
! tests element definitions
!
implicit none
integer, parameter :: NEN=4, NSD=2

!-pass
integer IEN
double precision X
dimension X(NSD,*),IEN(NEN,*)
!-locl
logical quad
integer i,j,n
double precision xl,area1,area2,area,a,b
dimension xl(NSD,NEN),a(NSD),b(NSD)
!
area1 = 0d0
area2 = 0d0

if (meshdatactx%Neglobal.lt.1) then
    write(stderr,10)
10  format(1x,'NUMEL < 1')
    call exitp(1)
endif
!
do n=1,meshdatactx%Nelocal

    do i=1,4
        if (IEN(i,n).lt.1 .or. &
            IEN(i,n).gt.meshdatactx%Nvglobal) then

            write(stderr,*) 'plnlib has a problem with element ', n, &
                            'which has vertex ', i, ':', IEN(i,n), &
                            'which is outside range [1,',meshdatactx%Nvglobal,']'
            call exitp(1)
        endif
    enddo

    quad = (IEN(3,n).ne.IEN(4,n))

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
        area = area1
    endif

    if (ABS(area).le.EPS) then
        write(stderr,30) n
30      format(1x,'degenerate element ',I8)
        call exitp(1)
    endif

    if (area.lt.ZERO) then
        write(stderr,40) i,meshdatactx%Nelocal
40      format(1x,'elmtest: local element ',I6,' is negatively ',&
         'oriented (node order is wrong), ',I12,' elements total')
        if (quad) then
            write(stderr,42) (IEN(j,n),(xl(i,j),i=1,NSD),j=1,4)
42          format('X(',I12,')=(',1PG15.6,',',1PG15.6,')'/)
        else
            write(stderr,42) (IEN(j,n),(xl(i,j),i=1,NSD),j=1,3)
        endif
        call exitp(1)
    endif

    if (ABS(area1).le.EPS) then
        write(stderr,50) n
50      format(1x,'WARNING: degenerate triangle 123 in element', &
          1X,I12)
        call exitp(1)
    endif

    if (area1.lt.ZERO) then
        write(stderr,60) n
60      format(1x,'WARNING: negatively oriented triangle 123', &
          ' in element ',I12)
        call exitp(1)
    endif

    if (quad .and. ABS(area2).le.EPS) then
        write(stderr,70) n
70      format(1x,'WARNING: degenerate triangle 341 in element', &
          1X,I12)
        call exitp(1)
    endif

    if (quad .and. area2.lt.ZERO) then
        write(stderr,80) n
80      format(1x,'WARNING: negatively oriented triangle 341', &
          ' in element ',I12)
        call exitp(1)
    endif

enddo
!
if (iecho.eq.2) then
     write(*,*) 'Rank ', getrank(), ' tested elements OK'
endif

return
end
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
 1   format(1x,'WARNING: effective dislocation creep rate not implemented for planar 2D')
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

!-------------------------------------------------------------------------------


subroutine MATERL (DMAT,E,POIS)

USE MODELDEFINITION
USE constants, only: zero, one, two, three, half

implicit none

integer, parameter :: NSTR=3
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
if (NOCOMPR.ne.0) then
    DMAT(1,1)=TWO*E/THREE
    DMAT(1,2)=ZERO
    DMAT(1,3)=ZERO
    DMAT(2,1)=ZERO
    DMAT(2,2)=DMAT(1,1)
    DMAT(2,3)=ZERO
    DMAT(3,1)=ZERO
    DMAT(3,2)=ZERO
    DMAT(3,3)=DMAT(1,1)
else
    if (IOPT.eq.0) then
!            Compute Lame parameters AM(=MU) and AL(=LAMBDA)
        AM=E/(ONE+POIS)
        AL=AM*POIS/(ONE-TWO*POIS)
        AM=HALF*AM
        DMAT(1,1)=TWO*AM+AL
        DMAT(1,2)=AL
        DMAT(1,3)=ZERO
        DMAT(2,1)=AL
        DMAT(2,2)=TWO*AM+AL
        DMAT(2,3)=ZERO
        DMAT(3,1)=ZERO
        DMAT(3,2)=ZERO
        DMAT(3,3)=AM
    else
        DMAT(1,1)=E/(ONE-POIS*POIS)
        DMAT(1,2)=POIS*DMAT(1,1)
        DMAT(1,3)=ZERO
        DMAT(2,1)=DMAT(1,2)
        DMAT(2,2)=DMAT(1,1)
        DMAT(2,3)=ZERO
        DMAT(3,1)=ZERO
        DMAT(3,2)=ZERO
        DMAT(3,3)=(ONE-POIS)*DMAT(1,1)*HALF
    endif
endif

return
end
!-------------------------------------------------------------------------------
subroutine BDELD (XL,DL,EE,QUAD,ierr)

use plotControl,     only: jShell

!
! Subroutine to compute strains from displacements in each element
!
implicit none
integer, parameter :: NDOF=2, NSD=2, NSTR=3, NEN=4
!-pass
integer          :: ierr
double precision :: XL(NSD,NEN),DL(NDOF,NEN),EE(NSTR)
logical          :: QUAD
!-locl
double precision :: BP(NDOF,NSD)
!
!    Calculate matrix of partial displacement derivatives
if (JSHELL.ne.0) then
    call SHLCRD (XL,QUAD)
endif

call BPMATRIX (XL,DL,QUAD,BP,ierr)

if (ierr.ne.0) then
    return
endif

EE(1)=BP(1,1)
EE(2)=BP(2,2)
EE(3)=BP(1,2)+BP(2,1)    ! fluid dynamical definition of shear strain

return
end subroutine
!-------------------------------------------------------------------------------
subroutine BMATRIX (B,SH)
!
! Computes the strain-displacement matrix B
!
implicit none

integer, parameter :: NEN=4, NSTR=3, NDOF=2
integer, parameter :: NEE=NDOF*NEN

!-pass
double precision :: B(NSTR,NEE)
double precision :: SH(3,NEN)

integer :: i,k

k = 1

do i=1,NEN
    B(1,k  ) = SH(1,i)
    B(1,k+1) = 0d0
    B(2,k  ) = 0d0
    B(2,k+1) = SH(2,i)
    B(3,k  ) = SH(2,i)
    B(3,k+1) = SH(1,i)
    k = k + 2
enddo

return
end
!-------------------------------------------------------------------------------
subroutine GRAVLD (P,XL,GRAV,TH,WT,QUAD,ierr)
!
! Computes the contribution to the load from gravity body forces of a single element.
! igrav 1 = in one direction
! igrav 2 = aimed at the domain origina
!
USE CONSTANTS,       only: G
use modeldefinition, only: igrav

implicit none
integer, parameter :: NSD=2, NDOF=2, NEN=4, NEE=NEN*NDOF
integer, parameter :: NINT=4
!-pass
integer ierr
logical QUAD
double precision :: P(NEE),XL(NSD,NEN),GRAV(NDOF)
double precision :: TH  ! thickness of the element 
double precision :: WT  ! density
!-locl
double precision :: XS(NSD,NSD),SH(3,NEN)

double precision :: sg, tg
dimension        :: SG(NEN),TG(NEN) 
save    tg, sg

integer :: i, j
double precision :: DET, grav1, grav2
double precision :: xc, yc ! center of element.

!-init
data TG/-1d0,-1d0,1d0,1d0/
data SG/-1d0,1d0,1d0,-1d0/

P = 0d0

do i=1,NINT

    call SHAP20 (SG(i)*G,TG(i)*G,XL,DET,SH,XS,QUAD,ierr)
    if (ierr.ne.0) then
        return
    endif

    if (igrav.eq.1) then
        ! gravity points in one direction for every element
        GRAV1 = GRAV(1) * DET * WT * TH
        GRAV2 = GRAV(2) * DET * WT * TH

    else if (igrav.eq.2) then
        ! gravity of all elements points to the origin,
        ! with a magnitude of the first entry.
        ! assume constant gravity. Alle nodes point in the same direction,
        ! while actually they should converge on the central point.
        xc = ( XL(1,1) + XL(1,2) + XL(1,3) ) / 3.0
        yc = ( XL(2,1) + XL(2,2) + XL(2,3) ) / 3.0
        
        ! both grav(1), because there is just one magnitude
        GRAV1 = sqrt ( (GRAV(1)**2) / ( 1.0+(yc**2 / xc**2) )  )
        GRAV2 = sqrt ( (GRAV(1)**2) / ( 1.0+(xc**2 / yc**2) )  )

        ! add sign, based on whether points are left/right down/up from origin
        GRAV1 =    sign(GRAV1,-xc) * DET * WT * TH 
        GRAV2 =    sign(GRAV2,-yc) * DET * WT * TH 

    else
        ! this is not supposed to happen
    endif

    do j=1,NEN
        P(2*j-1) = P(2*j-1) + GRAV1 * SH(3,j)
        P(2*j  ) = P(2*j  ) + GRAV2 * SH(3,j)
    enddo

enddo

return
end
!-------------------------------------------------------------------------------
subroutine STIFF (DMAT,XL,S,TH,QUAD,ierr, nodeID)
!
! Computes the super-diagonal part of the local stiffness matrix.
!
USE MODELDEFINITION, only: INCOMP
USE ALGEBRA,         only: IDFORM, CLEAR
USE CONSTANTS,       only: G

implicit none
integer, parameter :: NDOF=2, NSD=2, NSTR=3, NEN=4
integer, parameter :: NINT=4, NMODES=2,NEE=NDOF*NEN
integer, parameter :: NEX=NMODES*NDOF
!-pass
integer          :: ierr
logical          :: QUAD
double precision :: XL(NSD,NEN)
double precision :: DMAT(NSTR,NSTR)
double precision :: S(NEE,NEE)
double precision :: TH
integer          :: nodeID
!-locl
double precision :: SG(NEN)
double precision :: TG(NEN)
double precision :: XS(NSD,NSD)
double precision :: SH(3,NEN)
double precision :: B(NSTR,NEE)
double precision :: SDAL(NEE,NEX)
double precision :: SAAL(NEX,NEX)
double precision :: SDA(NEE,NEX)
double precision :: SAA(NEX,NEX)
integer          :: IDSAA(NEX)
double precision :: SA(NEE,NEE)
logical :: firstentry

double precision :: db1, db2, db3
double precision :: sa1, sa2, sa3, sa4
integer          :: i,j,k,l
double precision :: C, DET

save             :: SG,TG,firstentry,IDSAA
!-init
data firstentry /.true./
data TG/-1d0,-1d0,1d0,1d0/
data SG/-1d0,1d0,1d0,-1d0/

if (INCOMP.ne.0 .and. QUAD) then
    if (firstentry) then
!       set up list of diagonal addresses of square array
        call IDFORM (IDSAA,NEX)
        firstentry = .false.
    endif
    call CLEAR (SDA,32,"SDA")
    call CLEAR (SAA,16,"SAA")
endif

loop_element_nodes: do l=1,NINT
    call SHAP20 (SG(l)*G,TG(l)*G,XL,DET,SH,XS,QUAD,ierr)
    if (ierr.ne.0) then
        return
    endif

    call BMATRIX (B,SH)

    C = DET*TH
    do j=1,NEE
!       Compute column J of DB-matrix
        DB1=DMAT(1,1)*B(1,j) + &
            DMAT(1,2)*B(2,j) + &
            DMAT(1,3)*B(3,j)
        DB2=DMAT(2,1)*B(1,j) + &
            DMAT(2,2)*B(2,j) + &
            DMAT(2,3)*B(3,j)
        DB3=DMAT(3,1)*B(1,j) + &
            DMAT(3,2)*B(2,j) + &
            DMAT(3,3)*B(3,j)
        do i=1,j
            S(i,j)=S(i,j)+C*(B(1,i)*DB1 + B(2,i)*DB2 + B(3,i)*DB3)
        enddo
    enddo

!   Incompatible modes
    if (INCOMP.ne.0.and.QUAD) then
        call INCSTFF(SG(l)*G,TG(l)*G,DMAT,XL,B,SAAL,SDAL,QUAD,ierr)
        if (ierr.ne.0) then
            return
        endif

        do j=1,NEX
            do k=1,8
                SDA(k,j)=SDA(k,j)+C*SDAL(k,j)
            enddo
            do k=1,4
                SAA(k,j)=SAA(k,j)+C*SAAL(k,j)
            enddo
        enddo

    endif
enddo loop_element_nodes

if (INCOMP.ne.0.and.QUAD) then
!        Replace SAA by its inverse
    call INVMAT (SAA,IDSAA,NEX)
    do j=1,NEE
!        column j of inv(SAA)*SAD
        SA1 = SAA(1,1)*SDA(j,1) + &
              SAA(1,2)*SDA(j,2) + &
              SAA(1,3)*SDA(j,3) + &
              SAA(1,4)*SDA(j,4)

        SA2 = SAA(2,1)*SDA(j,1) + &
              SAA(2,2)*SDA(j,2) + &
              SAA(2,3)*SDA(j,3) + &
              SAA(2,4)*SDA(j,4)

        SA3 = SAA(3,1)*SDA(j,1) + &
              SAA(3,2)*SDA(j,2) + &
              SAA(3,3)*SDA(j,3) + &
              SAA(3,4)*SDA(j,4)

        SA4 = SAA(4,1)*SDA(j,1) + &
              SAA(4,2)*SDA(j,2) + &
              SAA(4,3)*SDA(j,3) + &
              SAA(4,4)*SDA(j,4)

!           SDA*inv(SAA)*SAD: matrix to subtract from S
        do i=1,j
            SA(i,j) = SDA(i,1)*SA1 + &
                      SDA(i,2)*SA2 + &
                      SDA(i,3)*SA3 + &
                      SDA(i,4)*SA4
        enddo
!       adjust S
        do i=1,j
            S(i,j) = S(i,j) - SA(i,j)
        enddo
    enddo
endif

return
end
!-------------------------------------------------------------------------------
subroutine INCSTFF (SI,TI,DMAT,XL,B,SAAL,SDAL,QUAD,ierr)
!
!       Computes those parts of the local stiffness matrix which
!       are due to incompatible shape functions. Defines SDAL and SAAL
!
use constants, only: zero, one, two

implicit none
integer, parameter :: NDOF=2, NEN=4, NSD=2, NSTR=3
integer, parameter :: NMODES=2
integer, parameter :: NEX=NMODES*NDOF, NEE=NEN*NDOF
!-pass
integer ierr
logical QUAD
double precision :: XL(NSD,NEN)
double precision :: DMAT(NSTR,NSTR)
double precision :: B(NSTR,NEE)
double precision :: SAAL(NEX,NEX)
double precision :: SDAL(NEE,NEX)
double precision :: SI, TI
!-locl
double precision :: XS(NSD,NSD)
double precision :: SH(3,NEN)
double precision :: SHE(3,2)
double precision :: BE(NSTR,NEX)

double precision :: DBE1, DBE2, DBE3
double precision :: DET

double precision :: dSHdS, dSHdT

integer          :: j
!
!       two extra shape functions and s- and t-derivatives
SHE(3,1) = ONE - SI*SI
SHE(1,1) = - TWO*SI
SHE(2,1) = ZERO
SHE(3,2) = ONE - TI*TI
SHE(1,2) = ZERO
SHE(2,2) = - TWO*TI
!
!       x- and y-derivatives
!       inverse jacobian evaluated in centre element
call SHAP20 (ZERO,ZERO,XL,DET,SH,XS,QUAD,ierr)

if (ierr.ne.0) then
    return
endif

dSHdS = SHE(1,1)
dSHdT = SHE(2,1)

SHE(1,1) = XS(1,1)*dSHdS + XS(1,2)*dSHdT
SHE(2,1) = XS(2,1)*dSHdS + XS(2,2)*dSHdT

dSHdS = SHE(1,2)
dSHdT = SHE(2,2)

SHE(1,2) = XS(1,1)*dSHdS + XS(1,2)*dSHdT
SHE(2,2) = XS(2,1)*dSHdS + XS(2,2)*dSHdT
!
!       extra part strain-displacement matrix
BE(1,1) = SHE(1,1)
BE(1,2) = ZERO
BE(2,1) = ZERO
BE(2,2) = SHE(2,1)
BE(3,1) = SHE(2,1)
BE(3,2) = SHE(1,1)
BE(1,3) = SHE(1,2)
BE(1,4) = ZERO
BE(2,3) = ZERO
BE(2,4) = SHE(2,2)
BE(3,3) = SHE(2,2)
BE(3,4) = SHE(1,2)
!
do j=1,NEX

    DBE1=DMAT(1,1)*BE(1,j)+DMAT(1,2)*BE(2,j)+DMAT(1,3)*BE(3,j)
    DBE2=DMAT(2,1)*BE(1,j)+DMAT(2,2)*BE(2,j)+DMAT(2,3)*BE(3,j)
    DBE3=DMAT(3,1)*BE(1,j)+DMAT(3,2)*BE(2,j)+DMAT(3,3)*BE(3,j)

    SDAL(1,j) = B(1,1)*DBE1+B(2,1)*DBE2+B(3,1)*DBE3
    SDAL(2,j) = B(1,2)*DBE1+B(2,2)*DBE2+B(3,2)*DBE3
    SDAL(3,j) = B(1,3)*DBE1+B(2,3)*DBE2+B(3,3)*DBE3
    SDAL(4,j) = B(1,4)*DBE1+B(2,4)*DBE2+B(3,4)*DBE3
    SDAL(5,j) = B(1,5)*DBE1+B(2,5)*DBE2+B(3,5)*DBE3
    SDAL(6,j) = B(1,6)*DBE1+B(2,6)*DBE2+B(3,6)*DBE3
    SDAL(7,j) = B(1,7)*DBE1+B(2,7)*DBE2+B(3,7)*DBE3
    SDAL(8,j) = B(1,8)*DBE1+B(2,8)*DBE2+B(3,8)*DBE3

    SAAL(1,j) = BE(1,1)*DBE1+BE(2,1)*DBE2+BE(3,1)*DBE3
    SAAL(2,j) = BE(1,2)*DBE1+BE(2,2)*DBE2+BE(3,2)*DBE3
    SAAL(3,j) = BE(1,3)*DBE1+BE(2,3)*DBE2+BE(3,3)*DBE3
    SAAL(4,j) = BE(1,4)*DBE1+BE(2,4)*DBE2+BE(3,4)*DBE3
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
subroutine ADDPR (B,X,D,IEN,LM,MAT,PRPMAT,IELNO,ISIDE,PRES, &
   LMF,TFAULT,SKEW,DXE)

USE MATERIALSMODULE, only: LMAT,MaterialIDOfElement
USE MODELDEFINITION
use modeldatamodule, only: modeldatactx
USE ALGEBRA,         only: formrt, equate, vmprd
use constants,       only: eps, zero, half, useDXE
use modelctx,        only: getrank ! only for debug purposes. Can be removed
use debugmodule,     only: iecho
use iomodule,        only: stderr, stdout
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact

!
! Subroutine to add (incremental) tractions derived from stresses to load vector.
! That, and pressures
!
implicit none
integer, parameter :: NEN=4, NSTR=3, NDOF=2, NSD=2

!-pass
integer :: IEN, LM, IELNO, ISIDE, LMF, MAT
double precision :: X, D, PRPMAT, B, PRES(NUMPR), TFAULT, SKEW, DXE
dimension MAT(*),IEN(NEN,*),X(NSD,*),D(NDOF,*),PRPMAT(6,*),B(*), &
   LM(NDOF,NEN,*),IELNO(*),ISIDE(*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN,*)

!-locl
double precision :: f, fl, xl, dl, xn, tmp, rot
integer          :: lien(4)
dimension f(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),tmp(NDOF), &
   rot(3,3),fl(NDOF)
double precision,  external :: FLOAT
integer :: globalelementnumber, entrynumber, materialnumber, ieqn, nodecounter, idof, globalnodenumber, numadd
double precision :: sxx,syy,anglevec(2),area,w
logical :: quad

numadd = 0

do entrynumber = 1,NUMPR

	write(*,*) "-------------- adding pressure", entrynumber, "----------"


    globalelementnumber = IELNO(entrynumber)

    if (globalelementnumber.le.0) then
        write(stderr,1) entrynumber,globalelementnumber
 1      format(1x,'ADDPR: fatal IELNO error. Entry ',I5,&
          ' has element number ',I8)
        call stoper()
    endif

    numadd = numadd + 1

    call LCOORD (X,xl,IEN(1,globalelementnumber),0)
    if (LGDEF.ne.0) then
        call LDISP  (dl,D,IEN(1,globalelementnumber),NDOF,NEN)
        call AddFaultDisplacement (dl,globalelementnumber,TFAULT,NDOF,NEN,25)
        call ADDSNE (dl,globalelementnumber,NDOF,NEN,useDXE)
        call REZONE (xl,dl,1d0)
    endif

!   fill lien with local side numbers for element side ISIDE(n)
    quad = IEN(3,globalelementnumber).ne.IEN(4,globalelementnumber)
    call SIDENP(ISIDE(entrynumber),quad,lien)

!   compute side length or surface area into variable AREA
    call FCGEOM (xl,lien,area,xn,"addpr")

    sxx = PRES(entrynumber)
    syy = PRES(entrynumber)

    if (ISIDE(entrynumber).ge.1 .and. ISIDE(entrynumber).le.4) then
        materialnumber = MaterialIDOfElement(globalelementnumber)
        w = PRPMAT(6,materialnumber)*area/2D0    ! element thickness * side length / 2
        f(1) = xn(1)*sxx*w
        f(2) = xn(2)*syy*w

		write(*,*) "f direct", f(1:2)

!       compute local nodal forces FL from global nodal forces F for first two nodes
        do nodecounter = 1, 2

            globalnodenumber = IEN(lien(nodecounter),globalelementnumber)

    	    write(*,*) "applyling to node ", nodecounter, globalnodenumber


            anglevec(1) = -SKEW(1,globalnodenumber)
            anglevec(2) = 0d0

            write(*,*) "rotating with vec: ", anglevec(1:2)

            if (ABS(anglevec(1)).gt.EPS) then
                call FORMRT (anglevec,rot,NDOF)
                call EQUATE (tmp,f,NDOF)
                call VMPRD (rot,tmp,fl,NDOF,3)
            else
                call EQUATE (fl,f,NDOF)
            endif

			write(*,*) "rotated to: ", fl(1:2)

!           add local forces to load vector 


#ifdef GIA
            do idof=2,2
#else
            do idof=1,NDOF
#endif
                ieqn = LM(idof,lien(nodecounter),globalelementnumber)
                if (ieqn.gt.0) then
                    write(*,*) "adding force", &
                               globalelementnumber, &
       	       	       	       globalnodenumber, &
                               nodecounter, &
                               ieqn, &
                               idof, fl(idof)
                    B(ieqn) = B(ieqn) + fl(idof)
                endif
            enddo
        enddo
    endif
enddo

if (iecho.eq.1 .and. numadd.gt.0) then
    write(stdout,3) numadd
3   format(1x,I8,' pressure loads are being added')
endif

return
end
!-------------------------------------------------------------------------------
subroutine FCGEOM (XL,LIEN,AREA,XN,callerName)
!
! Calculates area of element face and unit normal vector.
!
use constants, only: dfmin, half
use iomodule,        only: stderr

implicit none
integer, parameter :: NSD=2, NEN=4
!-pass
double precision :: XL(NSD,NEN), XN(NSD)
integer          :: LIEN(4)
double precision :: AREA, a(2), b(2)
character(len=*) :: callerName


if (LIEN(3).eq.0) then
!   Two local nodes only so the area is distance between these nodes (later
!   multiplied with element thickness)
    XN(1) = XL(2,LIEN(2)) - XL(2,LIEN(1))
    XN(2) = XL(1,LIEN(1)) - XL(1,LIEN(2))

    AREA  = XN(1)*XN(1)+XN(2)*XN(2)

    if (AREA.lt.dfmin) then
        write(stderr,1)
 1      format(1x,"FCGEOM: zero element face area")
        write(stderr,*) "Element spanned by points:"
        write(stderr,*) LIEN(1), LIEN(2)
        write(stderr,*) "Element with spanning vectors:"
        write(stderr,*) XL(2,LIEN(2)), XL(2,LIEN(1))
        write(stderr,*) XL(1,LIEN(1)), XL(1,LIEN(2))
        return
    endif

    AREA  = SQRT(AREA)
    XN(1) = XN(1)/AREA
    XN(2) = XN(2)/AREA
else
!   three/four nodes, so here we're dealing with element sides 5 or 6
    a(1) = xl(1,2) - xl(1,1)
    a(2) = xl(2,2) - xl(2,1)
    b(1) = xl(1,3) - xl(1,2)
    b(2) = xl(2,3) - xl(2,2)
    AREA = HALF*(a(1)*b(2)-a(2)*b(1))

    if (LIEN(4).ne.0) then
!       quad element
        a(1) = xl(1,4) - xl(1,3)
        a(2) = xl(2,4) - xl(2,3)
        b(1) = xl(1,1) - xl(1,4)
        b(2) = xl(2,1) - xl(2,4)
        AREA = AREA + HALF*(a(1)*b(2)-a(2)*b(1))
    endif
endif

return
end
!-------------------------------------------------------------------------------
subroutine LFLUX (XL,TL,IEN,PRPTEM,FLUX,ierr)
!
! Routine to calculate element heat flux
!
implicit none
integer, parameter :: NEN=4, NSD=2, NDOF=2
!-pass
integer          :: IEN(NEN)
double precision :: XL(NSD,NEN),TL(NEN),PRPTEM(6),FLUX(NSD)
!-locl
double precision :: GRAD(NDOF)
logical          :: QUAD
double precision :: sn, tn
integer ierr

sn = 0d0
tn = 0d0
!
QUAD = (IEN(3).ne.IEN(4))

if (.not.QUAD) then
    tn = -1d0/3d0
endif

!    calculate thermal gradients in the center of the element
call THGRAD (sn,tn,XL,TL,GRAD,QUAD,ierr)

if (ierr.ne.0) then
    return
endif

FLUX(1) = -PRPTEM(4)*GRAD(1)
FLUX(2) = -PRPTEM(5)*GRAD(2)
!
return
end
!-------------------------------------------------------------------------------
subroutine EFORCE (EVP,P,XL,TH,QUAD,ierr)
!
! Computes the effective forces at each node 
!
USE MODELDEFINITION, only : IPOINT
use algebra, only: idig

implicit none
integer, parameter :: NSD=2, NDOF=2, NSTR=3, NEN=4
integer, parameter :: NEE = NEN*NDOF

double precision, parameter :: GP=0.57735026918963d0
!-pass
integer          :: ierr
logical          :: QUAD
double precision :: XL(NSD,NEN),P(NEE),EVP(NSTR)
double precision :: TH
!-locl
double precision :: SG(NEN),TG(NEN),XS(NSD,NSD),SH(3,NEN),B(NSTR,NEE)
save SG,TG
integer          :: nint, L
double precision :: C, DET, W, G
!-init
 data SG/-1d0, 1d0, 1d0,-1d0/
 data TG/-1d0,-1d0, 1d0, 1d0/
!
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
        write(*,*) 'EFORCE: error in shape function'
        return
    endif

    call BMATRIX (B,SH)
    C = DET * W * TH

    ! EVP = local element stress change matrix
    ! P(:) = P(:) + C * B(:,:) * EVP(:)

    P(1)=P(1)+C*(B(1,1)*EVP(1)+B(2,1)*EVP(2)+B(3,1)*EVP(3))
    P(2)=P(2)+C*(B(1,2)*EVP(1)+B(2,2)*EVP(2)+B(3,2)*EVP(3))
    P(3)=P(3)+C*(B(1,3)*EVP(1)+B(2,3)*EVP(2)+B(3,3)*EVP(3))
    P(4)=P(4)+C*(B(1,4)*EVP(1)+B(2,4)*EVP(2)+B(3,4)*EVP(3))
    P(5)=P(5)+C*(B(1,5)*EVP(1)+B(2,5)*EVP(2)+B(3,5)*EVP(3))
    P(6)=P(6)+C*(B(1,6)*EVP(1)+B(2,6)*EVP(2)+B(3,6)*EVP(3))
    P(7)=P(7)+C*(B(1,7)*EVP(1)+B(2,7)*EVP(2)+B(3,7)*EVP(3))
    P(8)=P(8)+C*(B(1,8)*EVP(1)+B(2,8)*EVP(2)+B(3,8)*EVP(3))

enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine DEVSTN (STN,SDEV)

USE MODELDEFINITION, only: nocompr

! Computes stress deviator SDEV

implicit none

integer, parameter :: NSTR=3
!-pass
double precision   :: STN(NSTR),SDEV(NSTR)
!-local
!
if (NOCOMPR.eq.0) then
    SDEV(1) = 0.5d0 * (STN(1) - STN(2))
    SDEV(2) = -SDEV(1)
    SDEV(3) = STN(3)
else
    SDEV(1) = STN(1)
    SDEV(2) = STN(2)
    SDEV(3) = STN(3)
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine FORMBT (STN,BETA,BETB,PRPMAT,PRPLAS, elemID)

USE MODELCTX
USE MODELDEFINITION
use constants, only: dfmin,dlogmin,dfmax,dlogmax, &
                     eps, one, two, four, six, third, half
use algebra,   only: idig
!
! Subroutine to compute strain rates from stress field. It defines BETA and
! evaluates the product BETA*STRESS. BETB is the viscoplastic strain rate.
!
implicit none

integer, parameter :: NSTR=3
!-pass
double precision :: STN(NSTR),BETA(NSTR),BETB(NSTR),PRPMAT(6),PRPLAS(9)
integer          :: elemID
!-locl
double precision :: AVECT(NSTR)
double precision :: sigma, emhu, anpwr, streff
double precision :: sxx, syy, sxy

double precision, external :: effstn
double precision :: PSTRS

!-init
sigma = 0d0
emhu  = PRPMAT(3)  ! effective viscosity \mhu_{eff}
anpwr = PRPMAT(4)  ! 'n' in gtecton equations manual
!
if (IDIG(ICVIS,4,2).ne.0) then
    ! Calculate scalar stress PSTRS that will be compared to yield stress
    call YIELDS (STN,PRPLAS,PRPMAT,PSTRS, elemID)
    ! Construct plastic flow vector AVECT
    call YIELDF (STN,PRPLAS,PRPMAT,AVECT, elemID)
    ! Evaluate viscoplastic strain rate matrix BETB
    call FLOWVP (AVECT,STN,PRPLAS,PRPMAT,PSTRS,BETB, elemID)
endif

streff = EFFSTN(STN,HALF,.false.)

if (IOPT.eq.0) then
    ! plane strain
    if (streff.gt.dfmin) then
        sigma=((streff/emhu)**(anpwr-ONE))/(FOUR*emhu)
    endif
    SXX = sigma * (STN(1)-STN(2))
    SYY = sigma * (STN(2)-STN(1))
    SXY = sigma * FOUR * STN(3)
else
    ! plane stress
    if (streff.gt.dfmin) then
        sigma=((streff/emhu)**(anpwr-ONE))/(SIX*emhu)
    endif
    SXX = sigma * (TWO*STN(1)-STN(2))
    SYY = sigma * (TWO*STN(2)-STN(1))
    SXY = sigma * SIX * STN(3)
endif

if (IDIG(ICVIS,4,2).ne.0) then
    BETA(1) = SXX + BETB(1)
    BETA(2) = SYY + BETB(2)
    BETA(3) = SXY + BETB(3)
else
    BETA(1) = SXX
    BETA(2) = SYY
    BETA(3) = SXY
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine FORMJB (BETAJ,STN,PRPMAT,PRPLAS,BETB, elemID)

USE MODELCTX
USE MODELDEFINITION
USE ALGEBRA, only: CLEAR, idig
use constants, only: dfmin,dlogmin,dfmax,dlogmax,eps,third,half,one,two,four,six
!
! Constructs the Jacobian of the stress-strain rate matrix BETAJ.
!
implicit none

integer, parameter :: NSTR=3
!-pass
double precision :: BETAJ(6,6),STN(NSTR),PRPLAS(9),PRPMAT(6),BETB(NSTR)
integer          :: elemID
!-locl
double precision :: betajp(NSTR,NSTR),avect(NSTR)
double precision :: emhu, anpwr
double precision :: streff
double precision :: sxx,syy,sxy
double precision :: a1,b1,d1,a,b,d
double precision :: sigma, dnum, pstrs
integer          :: j
!
emhu  = PRPMAT(3)
anpwr = PRPMAT(4)

if (IOPT.eq.0) then
    ! strain
!   Melosh-Raefsky scaling (incompressible stress)
    streff = HALF*(STN(1)-STN(2))
    streff = SQRT(streff*streff+STN(3)*STN(3))
else
    ! stress
!   Scaling by incompressible stress
    streff = THIRD*(STN(1)*STN(1)-STN(1)*STN(2)+ &
     STN(2)*STN(2)) + STN(3)*STN(3)
    streff = SQRT(streff)
endif

if (streff.lt.dfmin) then
    call CLEAR (BETAJ,36,"BETAJ")
    return
endif
!
!       Plastic flow?
if (IDIG(ICVIS,4,2).ne.0) then
    call CLEAR (betajp,9,"betajp")
!           calculate scalar plastic stress PSTRS
    call YIELDS (STN,PRPLAS,PRPMAT,PSTRS, elemID)
!           Construct plastic flow vector AVECT
    call YIELDF (STN,PRPLAS,PRPMAT,avect, elemID)
!           construct Jacobian due to plastic yielding BETAJP
    call TANGVP (STN,PSTRS,PRPLAS,PRPMAT,BETB,avect,betajp, elemID)
endif
!
if (IOPT.eq.0) then
    ! plane strain
    sxx = STN(1)
    syy = STN(2)
    sxy = STN(3)

    a1 = (sxx-syy)/(TWO*streff)
    a1 = a1*a1

    b1 = sxy/streff
    b1 = b1*b1

    d1 = (sxx-syy)/(streff*streff)

    dnum = anpwr-ONE

    a = ONE+dnum*a1
    b = ONE+dnum*b1
    d = dnum*d1*sxy

    BETAJ(1,1) =  a
    BETAJ(1,2) = -a
    BETAJ(1,3) =  d
    BETAJ(2,1) = -a
    BETAJ(2,2) =  a
    BETAJ(2,3) = -d
    BETAJ(3,1) =  d
    BETAJ(3,2) = -d
    BETAJ(3,3) =  FOUR*b

    sigma = ((streff/emhu)**dnum)/(FOUR*emhu)

else
    ! IOPT = 1; plane stress
    sxx = STN(1)
    syy = STN(2)
    sxy = STN(3)

    a = (TWO*sxx-syy)/streff
    b = (TWO*syy-sxx)/streff
    d = sxy/streff

    dnum = anpwr-ONE

    BETAJ(1,1) =  TWO + dnum*a*a/SIX
    BETAJ(1,2) = -ONE + dnum*a*b/SIX
    BETAJ(1,3) =  dnum*a*d
    BETAJ(2,1) =  BETAJ(1,2)
    BETAJ(2,2) =  TWO + dnum*b*b/SIX
    BETAJ(2,3) =  dnum*b*d
    BETAJ(3,1) =  dnum*a*d
    BETAJ(3,2) =  dnum*b*d
    BETAJ(3,3) =  SIX+SIX*dnum*d*d
    sigma = ((streff/emhu)**dnum)/(SIX*emhu)
endif

if (IDIG(ICVIS,4,2).ne.0) then
    ! Add Jacobians of viscoplastic yield and viscous flow
    do j=1,3
        BETAJ(1,j) = sigma*BETAJ(1,j) + betajp(1,j)
        BETAJ(2,j) = sigma*BETAJ(2,j) + betajp(2,j)
        BETAJ(3,j) = sigma*BETAJ(3,j) + betajp(3,j)
    enddo
else
    ! Jacobian of viscous flow
    do j=1,3
        BETAJ(1,j) = sigma*BETAJ(1,j)
        BETAJ(2,j) = sigma*BETAJ(2,j)
        BETAJ(3,j) = sigma*BETAJ(3,j)
    enddo
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine ADDSTR (B,X,D,IEN,LM,MAT,PRPMAT,IELSTR,ISSIDE, &
                   ISTR,LMF,TFAULT,SKEW,DXE)

USE MATERIALSMODULE, only: LMAT, &
                           MaterialIDOfElement
USE MODELDEFINITION
use modeldatamodule, only: modeldatactx
USE ALGEBRA,         only: formrt, equate, vmprd
use constants,       only: eps, half, useDXE
use modelctx,        only: getrank ! only for debug purposes. Can be removed
use debugmodule,     only: iecho
use iomodule,        only: stderr, stdout
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact

!
! Subroutine to add (incremental) tractions derived from stresses to load vector.
!
implicit none

!-pass
integer, parameter :: NEN=4, NSTR=3, NDOF=2, NSD=2
integer :: ien, lm, ielstr, isside, istr, lmf, mat
double precision :: x, d, prpmat, b, tfault, skew, dxe
dimension MAT(*),IEN(NEN,*),X(NSD,*),D(NDOF,*),PRPMAT(6,*),B(*), &
   LM(NDOF,NEN,*),IELSTR(*),ISSIDE(*),ISTR(2,*), &
   LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN,*)

!-locl
integer          :: lien(4)
double precision :: f(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),tmp(NDOF), &
                    rot(3,3)

double precision,  external :: FLOAT

integer          :: it0, it1, globalelementnumber
integer          :: n, m, ii, nodecounter, idof, globalnodenumber
double precision :: sxx, sxy, syy, sxz, syz
double precision :: anglevec(2), area, w
integer          :: numadd
double precision :: fl(ndof)

logical :: quad

numadd = 0

do n=1,NUMSTR

    globalelementnumber = IELSTR(n)
    if (globalelementnumber.le.0) then
        write(stderr,1) n,globalelementnumber
 1      format(1x,'ADDSTR: fatal IELSTR error. Entry ',I5,&
          ' has element number ',I8)
        call stoper()
    endif

    it0 = modeldatactx%ISTR(1,n)
    it1 = modeldatactx%ISTR(2,n)

!   add the load vector components when the increment time range has started...
    if (NSTEP.ge.it0 .and. NSTEP.le.it1) then

        numadd = numadd + 1

        call LCOORD (X,xl,IEN(1,globalelementnumber),0)
        if (LGDEF.ne.0) then
            call LDISP  (dl,D,IEN(1,globalelementnumber),NDOF,NEN)
            call AddFaultDisplacement (dl,globalelementnumber,TFAULT,NDOF,NEN,25)
            call ADDSNE (dl,globalelementnumber,NDOF,NEN,useDXE)
            call REZONE (xl,dl,1d0)
        endif

!       fill lien with local side numbers for element side ISSIDE(n)
        quad = IEN(3,globalelementnumber).ne.IEN(4,globalelementnumber)
        call SIDENP(ISSIDE(n),quad,lien)

!       compute side length or surface area into variable AREA
        call FCGEOM (xl,lien,area,xn,"addstr")

        sxx = modeldatactx%STRS(1,n) / FLOAT(it1-it0+1)
        syy = modeldatactx%STRS(2,n) / FLOAT(it1-it0+1)
        sxy = modeldatactx%STRS(4,n) / FLOAT(it1-it0+1)
        sxz = modeldatactx%STRS(5,n) / FLOAT(it1-it0+1)
        syz = modeldatactx%STRS(6,n) / FLOAT(it1-it0+1)

        if (ISSIDE(n).eq.5 .or. ISSIDE(n).eq.6) then
!           traction on top (5) or underside (6) of elements
            if (quad) then
                w = area/4d0
            else
                w = area/3d0
            endif
            if (ISSIDE(n).eq.6) then
                w = -w
            endif
            f(1) = w*sxz
            f(2) = w*syz
        else
!           traction along one of element sides 1-4
            m = MaterialIDOfElement(globalelementnumber)
            w = PRPMAT(6,m)*area/2D0    ! element thickness * side length / 2
            f(1) = -(xn(1)*sxx + xn(2)*sxy)*w
            f(2) = -(xn(1)*sxy + xn(2)*syy)*w

        endif

! compute local nodal forces FL from global nodal forces F
! for first two nodes
        do nodecounter=1,2
! todo, why are 3 and 4 different than 1 and 2? me blind, or just stoopid?
            globalnodenumber = IEN(lien(nodecounter),globalelementnumber)
            anglevec(1) = -SKEW(1,globalnodenumber)
            anglevec(2) = 0d0

            if (ABS(anglevec(1)).gt.EPS) then
!               rotate node force info local coordinates (because of -skew(1))
                call FORMRT (anglevec,rot,NDOF)
                call EQUATE (tmp,f,NDOF)
                call VMPRD (rot,tmp,fl,NDOF,3)
            else
                call EQUATE (fl,f,NDOF)
            endif

!           add local forces to load vector 
            do idof=1,NDOF
                ii = LM(idof,lien(nodecounter),globalelementnumber)
                if (ii.gt.0) then
                    B(ii) = B(ii) + fl(idof)
                endif
            enddo
        enddo

        nodecounter = 3
        if (lien(nodecounter).ne.0) then
! compute local nodal forces fl from global nodal forces f
! for third node
            globalnodenumber = IEN(lien(nodecounter),globalelementnumber)
            anglevec(1) = -SKEW(1,globalnodenumber)
            anglevec(2) = 0d0

            if (ABS(anglevec(1)).gt.EPS) then
                call FORMRT (anglevec,rot,NDOF)
                call EQUATE (tmp,f,NDOF)
                call VMPRD (rot,tmp,fl,NDOF,3)
            else
                call EQUATE (fl,f,NDOF)
            endif

!           add local forces to load vector 
            do idof=1,NDOF
                ii = LM(idof,lien(nodecounter),globalelementnumber)
                if (ii.gt.0) then
                    B(ii) = B(ii) + fl(idof)
                endif
            enddo
        endif

        nodecounter = 4
        if (lien(nodecounter).ne.0) then
! compute local nodal forces fl from global nodal forces f
! for fourth node
            globalnodenumber = IEN(lien(nodecounter),globalelementnumber)
            anglevec(1) = -SKEW(1,globalnodenumber)
            anglevec(2) = 0d0

            if (ABS(anglevec(1)).gt.EPS) then
                call FORMRT (anglevec,rot,NDOF)
                call EQUATE (tmp,f,NDOF)
                call VMPRD (rot,tmp,fl,NDOF,3)
            else
                call EQUATE (fl,f,NDOF)
            endif

!           add local forces to load vector 
            do idof=1,NDOF
                ii = LM(idof,lien(nodecounter),globalelementnumber)
                if (ii.gt.0) then
                    B(ii) = B(ii) + fl(idof)
                endif
            enddo
        endif

     endif
enddo

if (iecho.eq.1 .and. numadd.gt.0) then
    write(stdout,3) numadd
 3  format(1x,I5,' stress loads are being added')
endif

return
end subroutine
!-------------------------------------------------------------------------------

#ifdef SPARSE
subroutine ADDWNK (IWELM,IWSIDE,IWTIME,WPRES,IEN,MAT,X,D,PRPMAT, &
                   A,LM,IDIAG,LMF,TFAULT,SKEW,DXE,LMX,NSLIP)

use spetscmodule,    only: PETSc_MatSetValuesLocal, M_A

USE MATERIALSMODULE, only: LMAT, &
                           MaterialIDOfElement
USE MODELDEFINITION
USE ALGEBRA
use constants,       only: eps, zero, half, one, useDXE
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact

#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: modeldatactx, &
                           elementSlipperyID
#else
use modeldatamodule, only: modeldatactx
#endif
use iomodule,        only: stderr

!
! Subroutine to add Winkler pressure contributions to the stiffness matrix
!
implicit none

integer, parameter :: NEN=4, NSTR=3, NDOF=2, NSD=2

!-pass
integer :: iwelm, iwside, iwtime, ien, idiag, lm, lmf, lmx, nslip, mat
double precision :: wpres, x, prpmat, a, d, tfault, skew, dxe
dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),A(*),IDIAG(NEQ),LM(NDOF,NEN,*), &
   D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
   DXE(NDOF,NEN,*),LMX(NDOF,NEN,*),NSLIP(5,*)



#include "petsc/finclude/petscvec.h"       
   PetscMPIInt irank

!-locl
integer          :: lien(4)
double precision :: xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),rot(3,3), &
   tmp(NDOF), f(NDOF), fl(NDOF)

logical quad
integer :: WinklerEntry, globalelementnumber, materialnumber, timecode, eqno, &
   globalnodenumber, idof, nodecounter
double precision :: angle, anglevec(2), wt, area
integer :: ierr



if (NUMWNK.eq.0) then
    return
endif

!write(*,*) "#############################################################"
!write(*,*) "#####   adding winkler pressures ############################"
!write(*,*) "#############################################################"

do WinklerEntry = 1, IABS(NUMWNK)

    globalelementnumber = IABS(IWELM(WinklerEntry))

!    write(*,*) "Winkler entry ", WinklerEntry, globalelementnumber


    timecode = IWTIME(WinklerEntry)
    if (timecode.lt.0) then
        if (timecode.ne.-1.and.NSTEP.ge.-timecode) then
            cycle
        endif
    elseif (timecode.eq.0) then
        cycle
    else
        if (NSTEP.lt.(timecode-1)) then
            cycle
        endif
    endif

    call LCOORD (X,XL,IEN(1,globalelementnumber),0)
 
    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,globalelementnumber),NDOF,NEN)
        call AddFaultDisplacement (DL,globalelementnumber,TFAULT,NDOF,NEN,15)
        call ADDSNE (DL,globalelementnumber,NDOF,NEN,useDXE)
        call REZONE (XL,DL,ONE)
    endif

!   fill LIEN with local side numbers of element side IWSIDE(WinklerEntry)
    quad = IEN(3,globalelementnumber).ne.IEN(4,globalelementnumber)
    call SIDENP(IABS(IWSIDE(WinklerEntry)),quad,lien)
    
!   compute side length or surface area into variable AREA
    call FCGEOM (XL,lien,area,xn,"pln, addwnk")

!   Winkler pressure on top (5) or underside (6) of elements has no effect in pln
    if (IWSIDE(WinklerEntry).ge.5 .or. IWSIDE(WinklerEntry).eq.6) then
        cycle
    endif

!   Winkler pressure perpendicular to one of element sides 1-4
    materialnumber = MaterialIDOfElement(globalelementnumber)
    wt = PRPMAT(6,materialnumber)*area/2D0    ! element thickness * side length / 2

!    f(1) = ABS(xn(1)) * WPRES(WinklerEntry) * wt
!    f(2) = ABS(xn(2)) * WPRES(WinklerEntry) * wt

!	write(*,*) "wink",WinklerEntry,"xn", xn, "stiffness: ", f

!    fl(1) = WPRES(WinklerEntry) * wt
    fl(2) = WPRES(WinklerEntry) * wt



!   write(*,*) "wink",WinklerEntry,"xn", xn, "stiffness: ", fl

!-------------------------------------------------------------------
! explanation by Lukas
!
! For many models, the Winkler pressures are a mechanism to
! simulate the forces induced by contrasts in density.
! On a round planet, these forces should only work in the radial direction.
! To accomplish this, nodes need to be rotated using an Euler angle,
! so that the y-axis is perpendicular to the surface.
! The magnitude provided in the input file is set to the 
! y-axis, and the x-axis will not provide Winkler contributions to the stifness matrix.


! compute local nodal forces FL from global nodal forces F for first two nodes
    do nodecounter = 1, 2
        globalnodenumber = IEN(lien(nodecounter),globalelementnumber)

!		write(*,*) "applyling to node ", nodecounter, globalnodenumber


!        anglevec(1) = -SKEW(1,globalnodenumber)
!        anglevec(2) = 0d0

!        if (ABS(anglevec(1)).gt.EPS) then
!            call FORMRT (anglevec,rot,NDOF)
!            call EQUATE (tmp,f,NDOF)
!            call VMPRD (rot,tmp,fl,NDOF,3)
!        else
!            call EQUATE (fl,f,NDOF)
!        endif


!       add local forces to load vector 



!        do idof=1,NDOF
        do idof=2,2  ! there is only a Winkler component in the Y-direction
            eqno = LM(idof,lien(nodecounter),globalelementnumber)

			write(*,*) "adding stiffness", idof, fl(idof)


            if (eqno.gt.0) then
                call PETSc_MatSetValuesLocal(M_A,eqno-1,eqno-1,fl(idof),20,globalelementnumber,0,0, ierr)
                if (ierr.ne.0) then
                    write(stderr,10) nodecounter,ierr
10                  format('Error: subroutine ADDWNK(,',I3,',) received and error fMatSetValues'/ &
                     'Error code', i5/'This is most likely a problem in GTecton.'/&
                     'Contact model support'/'Leaving GTecton...')
                    call exitp(1)
                endif
            endif
        enddo
    enddo
enddo
!
return

end
#endif
!-------------------------------------------------------------------------------
subroutine ELCRD (S,T,XELM,IEN,X,D,elemID,TFAULT,SKEW,FACTOR,ierr)
!
! Calculates (deformed) coordinates XELM at natural coordinates (R,S,T)
!
use MODELDEFINITION
use ALGEBRA, only: CLEAR
use constants, only: useDXE

implicit none
integer, parameter :: NSD=2, NDOF=2, NSTR=3, NEN=4
!-pass
integer ierr
integer elemID
integer :: ien
double precision :: factor
double precision :: s, t
double precision :: xelm, x, d, tfault, skew
dimension XELM(NSD),IEN(NEN),X(NSD,*),D(NDOF,*),TFAULT(NDOF,*), SKEW(2,*)
!-locl
double precision :: XL(NSD,NEN),DL(NDOF,NEN),SH(3,NEN),XS(NSD,NSD)
logical QUAD
double precision :: DET
integer :: i, j
!-init
call CLEAR (XELM,NSD,"XELM")

call LCOORD (X,XL,IEN,0)
if (FACTOR.gt.0.) then
    call LDISP  (DL,D,IEN,NDOF,NEN)
    call AddFaultDisplacement (DL,elemID,TFAULT,NDOF,NEN,16)
    call ADDSNE (DL,elemID,NDOF,NEN,useDXE)
    call REZONE (XL,DL,FACTOR)
endif

QUAD = (IEN(3).ne.IEN(4))
call SHAP20 (S,T,XL,DET,SH,XS,QUAD,ierr)

if (ierr.ne.0) then
    return
endif

do  j=1,NEN
    do  i=1,NSD
        XELM(i) = XELM(i) + SH(3,j)*XL(i,j)
    enddo
enddo

return
end
!-------------------------------------------------------------------------------
subroutine ELCTR (XELM,IEN,X,D,elemID,TFAULT, SKEW,FACTOR,ierr)

! Calculates element center coordinates XELM
!
USE MODELDEFINITION, only: ISHELL, radius
USE ALGEBRA,         only: CLEAR
USE CONSTANTS,       only: zero, deg2rad, third, useDXE
use plotControl,     only: jShell

implicit none

integer, parameter :: NSD=2, NDOF=2, NSTR=3, NEN=4
!-pass
integer          :: ierr
integer          :: elemID
integer          :: ien
double precision :: xelm, x, d, tfault, skew
dimension        :: XELM(NSD),IEN(NEN),X(NSD,*),D(NDOF,*), &
                    TFAULT(NDOF,*), SKEW(2,*)
double precision :: factor
!-locl
double precision :: xl(NSD,NEN),dl(NDOF,NEN),sh(3,NEN),xs(NSD,NSD)
logical          :: quad
double precision :: t, det
integer          :: i, j
double precision :: ylat

!-init
call CLEAR (XELM,NSD,"XELM")
!
quad = (IEN(3).ne.IEN(4))
t = ZERO
if (.not.quad) then
    t = -THIRD
endif
!
call LCOORD (X,xl,IEN,0)


if (FACTOR.gt.0.) then
    call LDISP  (dl,D,IEN,NDOF,NEN)
    call AddFaultDisplacement (dl,elemID,TFAULT,NDOF,NEN,17)
    call ADDSNE (dl,elemID,NDOF,NEN,useDXE)

    do j=1,NEN
        if (JSHELL.eq.0) then
            XL(1,j) = XL(1,j) + FACTOR*dl(1,j)
            XL(2,j) = XL(2,j) + FACTOR*dl(2,j)
        else
            ylat = XL(2,j)*DEG2RAD
            XL(1,j) = XL(1,j) + FACTOR*dl(1,j) / &
             (DEG2RAD*RADIUS*COS(ylat))
            XL(2,j) = XL(2,j) + FACTOR*dl(2,j) / &
             (DEG2RAD*RADIUS)
        endif
    enddo
endif
!
call SHAP20 (ZERO,t,XL,DET,SH,XS,quad,ierr)

if (ierr.ne.0) then
    return
endif
!
do j=1,NEN
    XELM(1) = XELM(1) + SH(3,j)*xl(1,j)
    XELM(2) = XELM(2) + SH(3,j)*xl(2,j)
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine TSTIFF (XL,QUAD,S,C,TH,ierr)

USE CONSTANTS, only: G
!
! computes the local conduction stiffness matrix from integration at Barlow
! points.
!
implicit none

integer, parameter :: NEN=4, NINT=4, NSD=2

!-pass
integer          :: ierr
logical          :: QUAD
double precision :: XL(2,NEN)
double precision :: S(NEN,NEN)
double precision :: C(NSD)
double precision :: TH
!-locl
double precision :: sh(3,NEN)
double precision, save :: sg(NEN)
double precision, save :: tg(NEN)
double precision :: db(NSD,NEN)
double precision :: xs(NSD,NSD)
data tg/-1d0,-1d0, 1d0, 1d0/
data sg/-1d0, 1d0, 1d0,-1d0/

double precision :: detJ, dV
integer          :: i,j,k,l

do l=1,NINT

    call shap20 (sg(l)*G,tg(l)*G,XL,detJ,sh,xs,QUAD,ierr)

    if (ierr.ne.0) then
        return
    endif

    dV = detJ * TH

    do j=1,NEN
        do i=1,NSD
            db(i,j) = C(i) * sh(i,j) * dV
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
subroutine ELMSIZE (XL,QUAD,SIZE)
!
! computes characteristic sizes of element
!
use constants, only: zero, third, half, two

implicit none
integer, parameter :: NSD=2, NDOF=2, NSTR=3, NEN=4
!-pass
double precision :: XL(NSD,NEN),SIZE(NSD)
logical          :: QUAD
!-locl
double precision :: sh(NSD,NEN)
double precision :: xs(NSD,NSD)
double precision :: sa(NEN),ta(NEN)
save sa,ta
integer          :: i
!-init
data sa/-0.5d0, 0.5d0,0.5d0,-0.5d0/
data ta/-0.5d0,-0.5d0,0.5d0, 0.5d0/
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

SIZE(1) = (ABS(xs(1,1))+ABS(xs(2,1)))*TWO
SIZE(2) = (ABS(xs(1,2))+ABS(xs(2,2)))*TWO
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
 subroutine LHEAT (XL,QUAD,PL,HEAT,TH,ierr)
!
! computes the local load vector contribution of heat production.
!
implicit none
integer, parameter :: NEN=4, NINT=4, NSD=2
double precision,  parameter :: G=0.57735026918963d0
!-pass
 integer ierr
 logical QUAD
double precision :: XL(2,NEN),PL(NEN)
double precision :: heat
!-locl
double precision :: sh(3,NEN),sg(NEN),tg(NEN),xs(NSD,NSD)
save tg,sg
 data tg/-1d0,-1d0, 1d0, 1d0/
 data sg/-1d0, 1d0, 1d0,-1d0/
integer :: l, j
double precision :: dV, detJ, TH

do l=1,NINT
    call shap20 (sg(l)*G,tg(l)*G,XL,detJ,sh,xs,QUAD,ierr)
    if (ierr.ne.0) then
        return
    endif
    dV = detJ * TH
    do j=1,NEN
        PL(j) = PL(j) + HEAT*dV*sh(3,j)
    enddo
enddo

return
end
!-------------------------------------------------------------------------------
double precision function EigMax (STN,E,v,emhu,Npl,alpha)

USE MODELDEFINITION
use constants, only: dfmin,dlogmin,dfmax,dlogmax,eps, &
                     one, two, three, four, eight
use iomodule,        only: stderr
!
! Gives an upper limit estimate of the maximum local eigenvalue (of what?!)
!
 implicit none
!-pass
! SNGL    real E,v,emhu,Npl,alpha
double precision :: E,v,emhu,Npl,alpha
double precision ::  STN(3)


!-locl
double precision :: a,b,c,d, se, nplmin1, tmp
!
NplMin1 = Npl - one
SE      = (STN(1)-STN(2))*(STN(1)-STN(2))*2.5D-1+STN(3)*STN(3)

if (SE.gt.dfmin) then
    SE = SQRT(SE)
endif

tmp     = eight * (one+v)*(one-two*v) * emhu**Npl
tmp     = E * SE**NplMin1 / tmp

if (alpha.le.dfmin) then
    !        explicit stepping
    if (SE.gt.dfmin) then
        A = one + NplMin1 * ((STN(1)-STN(2))/(two * SE) ) ** two
        B = NplMin1 * (STN(1) - STN(2) ) / SE**2 * STN(3)
        C = four * ( one + NplMin1 * (STN(3)/SE)**two )
        D = (two*A - C)**two + eight*B**two
        if (D.ge.0.) then
            D = two*A + C + SQRT(D)
        else
            write(stderr,1)
 1               format(1x,'EigMax: determinant less than zero')
            D = one
        endif
    else
        D = one
    endif
else
!        implicit differencing
    D = two * (three*Npl + SQRT(Npl**two + two*NplMin1**two))
endif

EigMax = tmp * D
return
end
!-----------------------------------------------------------------------
subroutine ADDFLX (BT,IEN,X,D,MAT,PRPMAT,IFLS,IFLX,BFLX,LMT, &
       LMTX,LMF,TFAULT,SKEW,DXE)!,NFLX)
!
! Routine to add normal boundary fluxes to load vector
!
use MATERIALSMODULE, only : LMAT
use MODELDEFINITION
use constants, only: useDXE

implicit none
integer, parameter :: NEN=4, NSD=2, NDOF=2
!-pass
integer          :: IEN, IFLS, IFLX, LMT, LMF, LMTX, MAT
double precision :: BT, X, PRPMAT, BFLX, D, TFAULT, SKEW, DXE
dimension BT(*),IEN(NEN,*),X(NSD,*),MAT(*),PRPMAT(6,*),IFLS(*), &
   IFLX(*),BFLX(*),LMT(NEN,*),D(NDOF,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),LMTX(NEN,*),SKEW(2,*),DXE(NDOF,NEN,*)
!-locl
double precision :: XL(NSD,NEN),DL(NDOF,NEN)
integer sidenumber,eqno,lien(4),indexvalue,materialnumber,&
   globalelementnumber,nodeindex
double precision :: f,thick,area,xn(2)
logical quad

if (NFLX.le.0) then
    return
endif

do indexvalue=1,NFLX

    globalelementnumber = IFLX(indexvalue)
    if (IFLS(indexvalue).ge.1 .and. IFLS(indexvalue).le.4) then

!       define element coordinates
        call LCOORD (X,XL,IEN(1,globalelementnumber),0)
        if (LGDEF.ne.0) then
            call LDISP  (DL,D,IEN(1,globalelementnumber),NDOF,NEN)
            call AddFaultDisplacement (DL,globalelementnumber,TFAULT,NDOF,NEN,18)
            call ADDSNE (DL,globalelementnumber,NDOF,NEN,useDXE)
            call REZONE (XL,DL,1d0)
        endif

!       fill lien with local side numbers for element side ISSIDE(n)
        quad = IEN(3,globalelementnumber).ne.IEN(4,globalelementnumber)
        call SIDENP(IFLS(indexvalue),quad,lien)

!       compute side length or surface area into variable AREA
        call FCGEOM (xl,lien,area,xn,"addflx")

        if (LMAT) then
            materialnumber = globalelementnumber
        else
            materialnumber = MAT(globalelementnumber)
        endif
        thick = PRPMAT(6,materialnumber)

        f = -area * thick * BFLX(indexvalue)/2d0

        do nodeindex=1,2
            eqno = LMT(lien(nodeindex),globalelementnumber)
            if (eqno.gt.0) then
                BT(eqno) = BT(eqno) + f
            endif
        enddo
    endif
enddo

return
end
!-------------------------------------------------------------------------------
subroutine TEMWNK (AA,LMT,ITDIAG,ITWINK,TWINK,IEN,MAT,PRPMAT, &
                   X,D,LMF,TFAULT,SKEW,DXE)
!
! Subroutine to add Winkler fluxes to diagonal elements of thermal stiffness
! matrix
!
USE MATERIALSMODULE, only: LMAT
USE MODELDEFINITION
use debugmodule,     only: iecho
use constants,       only: useDXE, half, one
#ifdef SPARSE
use spetscmodule,     only: T_A
#endif
use petscksp,        only: add_values
use iomodule,        only: stderr, stdout
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact

implicit none
integer, parameter :: NSD=2, NDOF=2, NSTR=3, NEN=4
!-pass
integer :: itwink, ien, itdiag, lmt, lmf, mat
double precision :: twink, x, prpmat, aa, d, tfault, skew, dxe
 dimension ITWINK(3,*),TWINK(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),AA(*),ITDIAG(*),LMT(NEN,*), &
   D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), SKEW(2,*), &
   DXE(NDOF,NEN,*)

#ifdef SPARSE

#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscvec.h" 
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscksp.h"
   PetscMPIInt irank
#endif
!-locl
double precision :: xl(NSD,NEN),dl(NDOF,NEN),xn(NSD)
integer          :: lien(4)
integer          :: n, nn, ii, ll, ierr
double precision :: area, wt
!
if (NTWINK.eq.0) then
    return
endif
!
if (iecho.ne.0) then
    write(stdout,*) 'Adding Winkler fluxes'
endif
!
do nn=1,NTWINK
    n = ITWINK(1,nn)
    LL= ITWINK(3,nn)

    if (LL.lt.0) then
        if (LL.ne.-1.and.NSTEP.ge.-LL) then
            cycle
        endif
    elseif (LL.eq.0) then
        cycle
    else
        if (NSTEP.lt.(LL-1)) then
            cycle
        endif
    endif

    call LCOORD (X,XL,IEN(1,n),0)

    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,19)
        call ADDSNE (DL,n,NDOF,NEN,useDXE)
        call REZONE (XL,DL,ONE)
    endif

    call SIDENP (ITWINK(2,nn),IEN(3,n).ne.IEN(4,n),lien)

    call FCGEOM (XL,lien,area,xn,"pln, temwnk")    ! compute face area

    if (LMAT) then
        wt = PRPMAT(6,n) * area * TWINK(nn) * HALF
    else
        wt = PRPMAT(6,MAT(n)) * area * TWINK(nn) * HALF
    endif

    ii = LMT(lien(1),n)            ! equation number

    if (ii.gt.0) then
#ifdef SPARSE
        call MatSetValuesLocal(T_A,1,ii-1,1,ii-1,wt,ADD_VALUES,ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
 1          format(1x,'TEMWNK: MatSetValuesLocal error ',I4)
            call stoper()
        endif
#else
        ii = ITDIAG(ii)    ! matrix address of eqno
        AA(ii) = AA(ii) + wt
#endif
! SPARSE
    endif

    ii = LMT(lien(2),n)

    if (ii.gt.0) then
#ifdef SPARSE
        call MatSetValuesLocal(T_A,1,ii-1,1,ii-1,wt,ADD_VALUES,ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
            call stoper()
        endif
#else
        ii = ITDIAG(ii)
        AA(ii) = AA(ii) + wt
#endif
    endif
enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
 subroutine THGRAD (S,T,XL,TL,GRAD,QUAD,ierr)
!
! Routine for calculating thermal gradients at natural coordinates s,t
!
use constants, only: zero, one, quart

implicit none

integer, parameter :: NSD=2, NEN=4 

integer, parameter :: NPE=NSD+1

!-pass
integer ierr
logical QUAD
double precision :: s, t
double precision :: XL(NSD,NEN),TL(NEN),GRAD(NSD)
!-locl
double precision :: SA(NEN),TA(NEN),SH(NPE,NEN),XS(NSD,NSD)
double precision :: det
double precision :: dTds, dTdt
 save SA,TA
integer :: i
!-init
 data SA/-1d0, 1d0, 1d0, -1d0/
 data TA/-1d0,-1d0, 1d0,  1d0/
!
!    calculate inverse Jacobian matrix XS
call SHAP20 (s,t,XL,DET,SH,XS,QUAD,ierr)
if (ierr.ne.0) then
    return
endif
!    Calculate s- and t-derivatives of shape functions
do i=1,NEN
    SH(1,i) = QUART * SA(i) * (ONE + TA(i)*t)
    SH(2,i) = QUART * TA(i) * (ONE + SA(i)*s)
enddo

if (.not.QUAD) then
!        adjust shape function of third node for triangular elements
    do i=1,2
        SH(i,3) = SH(i,3)+SH(i,4)
        SH(i,4) = ZERO
    enddo
endif
!    calculate s- and t- thermal gradients
dTds = ZERO
dTdt = ZERO

do i=1,NEN
    dTds = dTds + SH(1,i)*TL(i)
    dTdt = dTdt + SH(2,i)*TL(i)
enddo

GRAD(1) = dTds*XS(1,1) + dTdt*XS(1,2)
GRAD(2) = dTds*XS(2,1) + dTdt*XS(2,2)
!
return
end
!-------------------------------------------------------------------------------
double precision function Tcentr (T,TX,XL,IEN,LMTX,ierr)
!
! Calculates the temperature in the center of the element
!
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
 1       format(1x,'Tcentr: error in calculation of shape functions')
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
implicit none

integer, parameter :: NEN=4, NINT=4, NSD=2
double precision, parameter :: G=0.57735026918963d0
!-pass
integer          :: ierr
logical          :: QUAD
double precision :: XL(2,NEN),S(NEN,NEN)
double precision :: RhoCp, Th
!-locl
double precision :: sh(3,NEN),sg(NEN),tg(NEN),xs(NSD,NSD)
save tg,sg
data tg/-1d0,-1d0, 1d0, 1d0/
data sg/-1d0, 1d0, 1d0,-1d0/
integer :: i,j,l
double precision :: dV, detJ

do l=1,NINT
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
function AJ1 (ST,POIS)

USE MODELDEFINITION, only: nocompr, &
                           iopt
use constants,       only: one
!
! Calculates first invariant of stress tensor
!
implicit none
!-pass
double precision :: AJ1
double precision :: ST(3)
double precision :: pois
!
if (NOCOMPR.eq.0) then
    if (IOPT.eq.0) then
        AJ1 = (ONE+POIS)*(ST(1)+ST(2))
    else
        AJ1 = ST(1)+ST(2)
    endif
else
    AJ1 = 0d0
endif
!
return
end
!-------------------------------------------------------------------------------
function AJ2 (ST,POIS)

USE MODELDEFINITION, only: iopt, &
                           nocompr
use constants,       only: half, one, third
!
! Calculates second invariant of deviatoric stress tensor
!
implicit none
!-pass
double precision :: AJ2
double precision :: ST(3)
double precision :: pois
!-local
double precision :: sp
double precision :: sxx, syy, szz

!write(*,*) "aj2 says: ", iopt, nocompr

if (IOPT.eq.0) then
    if (NOCOMPR.eq.0) then
        sp = (ONE+POIS)*(ST(1)+ST(2))*THIRD
        szz = POIS*(ST(1)+ST(2)) - sp
        sxx = ST(1)-sp
        syy = ST(2)-sp
    else
        sxx = ST(1)
        syy = ST(2)
        szz = 0d0
    endif
else
    if (NOCOMPR.eq.0) THEN
        sp = (ST(1)+ST(2))*THIRD
        sxx = ST(1)-sp
        syy = ST(2)-sp
        szz = -sp
    else
        sxx = ST(1)
        syy = ST(2)
        szz = -ST(1)-ST(2)
    endif
endif
!
AJ2 = (sxx*sxx+syy*syy+szz*szz)*HALF+ST(3)*ST(3)
!AJ2 = (sxx*sxx+syy*syy+szz*szz)-ST(3)*ST(3)

!
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
double precision :: ST(3), POIS
!-local
double precision :: sp, sxx, syy, szz
!-init
!
if (IOPT.eq.0) then
!           plane strain
    sp = (ONE+POIS)*(ST(1)+ST(2))*THIRD
    szz = THIRD*(TWO*POIS - ONE)*(ST(1)+ST(2))
    sxx = ST(1)-sp
    syy = ST(2)-sp
else
!           plane stress
    sp = (ST(1)+ST(2))*THIRD
    sxx = ST(1)-sp
    syy = ST(2)-sp
    szz = -sp
endif
!
AJ3 = szz * ( (sxx*syy) - (ST(3)*ST(3)) )
!
return
end
!-------------------------------------------------------------------------------
subroutine REZONE (XL,DL,FACTOR)
!
! Routine to update the local coordinate XL by displacements DL for
! large deformation update
!
implicit none

integer, parameter :: NSD=2, NDOF=2, NEN=4

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

 integer, parameter :: NSTR=3
!-pass
 double precision :: EE(NSTR)
 double precision :: DT, VEXP

 EE(1) = DT*VEXP/3d0
 EE(2) = EE(1)
 EE(3) = 0d0

 return
 end
!-------------------------------------------------------------------------------
subroutine YIELDS (STN,PRPLAS,PRPMAT,PSTRS, elemID)

USE CONSTANTS, only: root3, deg2rad, zero, half, one, two,  &
                     three, dfmin

! Routine to calculate the scalar stress PSTRS that will be compared to
! the (equivalent) uniaxial yield stress. See GTECTON notes, section
! on yield criteria.

implicit none

integer, parameter :: NSTR=3
!-pass
double precision :: STN(NSTR)
double precision :: PRPLAS(9)
double precision :: PRPMAT(6)
double precision :: PSTRS
integer          :: elemID

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
!
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
parameter (NSTR=3)
!-pass
integer          :: elemID
double precision STN,PRPLAS,AVECT,PRPMAT
dimension STN(NSTR),PRPLAS(9),AVECT(NSTR),PRPMAT(6)

!-locl
integer i,nptype
double precision stneff,EFFSTN,i2,sp,sxx,syy,szz,veca1(NSTR), &
 veca2(NSTR),veca3(NSTR),const1,const2,const3,sint3,theta,cost, &
 cos3t,sint,tant,tan3t,frangl,sinp,AJ1,AJ3,alpha
external EFFSTN,AJ1,AJ3
!

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
sxx = STN(1) - sp
syy = STN(2) - sp
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
    szz = -sp
    veca3(1) = syy*szz + i2/THREE
    veca3(2) = sxx*szz + i2/THREE
    veca3(3) = -TWO * szz * STN(3)
endif

veca2(1) = sxx/(TWO*stneff)
veca2(2) = syy/(TWO*stneff)
veca2(3) = STN(3)/stneff
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
    write(*,*) "Plnlib subroutine YIELDF encountered and error"
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
subroutine TANGVP (STN,PSTRS,PRPLAS,PRPMAT,BETB,AVECT,BETAJP, elemID)
!
! Compute BETAJP, which is the Jacobian matrix for viscoplastic flow
! See GTECTON notes for details.
!

USE MODELDEFINITION
USE ALGEBRA, only: CLEAR
USE CONSTANTS, only: ZERO,HALF,ONE,TWO,THREE,FOUR,root3, &
    THIRD,SIX,deg2rad,dfmin,dlogmin,dfmax,dlogmax,eps

implicit none

integer NSTR
parameter (NSTR=3)
!-pass
double precision STN(NSTR),BETAJP(NSTR,NSTR),PRPLAS(9), &
 PRPMAT(6),BETB(NSTR),AVECT(NSTR),PSTRS
integer          :: elemID
!-locl
integer nptype,nflow,j
double precision stneff,AJ1,AJ2,AJ3,fnorm,uniyld,i1,sxx,syy,sxy, &
 yldpwr,cmult,gradp,sint3,theta,p1,p2,p3,frangl,fluidy,sinp
double precision m1(NSTR,NSTR),m2(NSTR,NSTR),m3(NSTR,NSTR)
external AJ1,AJ2,AJ3
!-init
call CLEAR(BETAJP,NSTR*NSTR,"betajp")
nptype = INT(PRPLAS(1))
if (nptype.eq.0) return
stneff = AJ2(STN,HALF)
if (stneff.le.dfmin) return
stneff = SQRT(stneff)
call INYLD (STN,PRPLAS,PRPMAT,PSTRS,BETB,fnorm,uniyld)
if (fnorm.le.0.) return
i1 = AJ1(STN,HALF)
!       deviatoric stress tensor components
sxx = STN(1)-i1*THIRD
syy = STN(2)-i1*THIRD
sxy = STN(3)
!
fluidy = PRPLAS(4)
yldpwr = PRPLAS(5)
nflow  = INT(PRPLAS(6))
!
!       Determine PHI
if (nflow.eq.1) then
!           powerlaw
    cmult = fnorm**yldpwr
    gradp = yldpwr * (FNORM**(yldpwr-ONE)) / uniyld
else
!           exponential
    cmult = EXP(fnorm*yldpwr) - ONE
    gradp = yldpwr * EXP(fnorm*yldpwr) / uniyld
endif
!
if (nptype.eq.1) then
!           Tresca
    sint3 = -THREE*root3*AJ3(STN,HALF)/(TWO*stneff**THREE)
    sint3 = MIN(sint3,ONE)
    sint3 = MAX(sint3,-ONE)
    theta = asin(sint3)*THIRD
elseif (nptype.eq.2) then
!           von Mises
    p1 =  CMULT * root3*fluidy/(TWO*stneff)
    p2 = -CMULT * root3*fluidy/(FOUR*stneff**3)
    p3 = fluidy * gradp
    if (IOPT.eq.0) then
!               plane strain
        m1(1,1) = HALF
        m1(1,2) = -HALF
        m1(1,3) = ZERO
        m1(2,2) = HALF
        m1(2,3) = ZERO
        m1(3,3) = TWO
        m1(2,1) = m1(1,2)
        m1(3,1) = m1(1,3)
        m1(3,2) = m1(2,3)
    else
!               plane stress
        m1(1,1) = TWO*THIRD
        m1(1,2) = -THIRD
        m1(1,3) = ZERO
        m1(2,2) = m1(1,1)
        m1(2,3) = ZERO
        m1(3,3) = TWO
        m1(2,1) = m1(1,2)
        m1(3,1) = m1(1,3)
        m1(3,2) = m1(2,3)
    endif
    m2(1,1) = sxx*sxx
    m2(1,2) = sxx*syy
    m2(1,3) = TWO*sxx*sxy
    m2(2,2) = syy*syy
    m2(2,3) = TWO*syy*sxy
    m2(3,3) = FOUR*sxy*sxy
    m2(2,1) = m2(1,2)
    m2(3,1) = m2(1,3)
    m2(3,2) = m2(2,3)
    m3(1,1) = AVECT(1)*AVECT(1)
    m3(1,2) = AVECT(1)*AVECT(2)
    m3(1,3) = AVECT(1)*AVECT(3)
    m3(2,2) = AVECT(2)*AVECT(2)
    m3(2,3) = AVECT(2)*AVECT(3)
    m3(3,3) = AVECT(3)*AVECT(3)
    m3(2,1) = m3(1,2)
    m3(3,1) = m3(1,3)
    m3(3,2) = m3(2,3)
    do j=1,NSTR
        BETAJP(1,j) = P1*m1(1,j) + P2*m2(1,j) + P3*m3(1,j)
        BETAJP(2,j) = P1*m1(2,j) + P2*m2(2,j) + P3*m3(2,j)
        BETAJP(3,j) = P1*m1(3,j) + P2*m2(3,j) + P3*m3(3,j)
    enddo
elseif (nptype.eq.3) then
!           Mohr-Coulomb: incomplete
    frangl = PRPLAS(2) * deg2rad
    sinp  = SIN(frangl)
    sint3 = -THREE*root3*AJ3(STN,HALF)/(TWO*stneff**THREE)
    sint3 = MIN(sint3,ONE)
    sint3 = MAX(sint3,-ONE)
    theta = asin(sint3)*THIRD
elseif (nptype.eq.4) then
!           Drucker-Prager
    p1 =  CMULT * fluidy/(TWO*stneff)
    p2 = -CMULT * fluidy/(FOUR*stneff**3)
    p3 = fluidy * gradp
    if (IOPT.eq.0) then
!               plane strain
        m1(1,1) = HALF
        m1(1,2) = -HALF
        m1(1,3) = ZERO
        m1(2,2) = HALF
        m1(2,3) = ZERO
        m1(3,3) = TWO
        m1(2,1) = m1(1,2)
        m1(3,1) = m1(1,3)
        m1(3,2) = m1(2,3)
    else
!               plane stress
        m1(1,1) = TWO*THIRD
        m1(1,2) = -THIRD
        m1(1,3) = ZERO
        m1(2,2) = m1(1,1)
        m1(2,3) = ZERO
        m1(3,3) = TWO
        m1(2,1) = m1(1,2)
        m1(3,1) = m1(1,3)
        m1(3,2) = m1(2,3)
    endif
    m2(1,1) = sxx*sxx
    m2(1,2) = sxx*syy
    m2(1,3) = TWO*sxx*sxy
    m2(2,2) = syy*syy
    m2(2,3) = TWO*syy*sxy
    m2(3,3) = FOUR*sxy*sxy
    m2(2,1) = m2(1,2)
    m2(3,1) = m2(1,3)
    m2(3,2) = m2(2,3)
    m3(1,1) = AVECT(1)*AVECT(1)
    m3(1,2) = AVECT(1)*AVECT(2)
    m3(1,3) = AVECT(1)*AVECT(3)
    m3(2,2) = AVECT(2)*AVECT(2)
    m3(2,3) = AVECT(2)*AVECT(3)
    m3(3,3) = AVECT(3)*AVECT(3)
    m3(2,1) = m3(1,2)
    m3(3,1) = m3(1,3)
    m3(3,2) = m3(2,3)
    do j=1,NSTR
        BETAJP(1,j) = p1*m1(1,j) + p2*m2(1,j) + p3*m3(1,j)
        BETAJP(2,j) = p1*m1(2,j) + p2*m2(2,j) + p3*m3(2,j)
        BETAJP(3,j) = p1*m1(3,j) + p2*m2(3,j) + p3*m3(3,j)
    enddo
endif
!
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
implicit none
integer NSTR
parameter (NSTR=3)
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
!
allow  = 1d-2
frangl = PRPLAS(2) * deg2rad
UNIYLD = PRPLAS(3)
stnhrd = PRPLAS(7)
!
if (nptype.eq.3) then
!           Mohr-Coulomb
    UNIYLD = UNIYLD * COS(frangl)
else if (nptype.eq.4) then
!           Drucker-Prager outer apices
    UNIYLD = SIX*UNIYLD*COS(frangl)/(root3*(THREE-SIN(frangl)))
!           Drucker-Prager inner apices
!           UNIYLD = SIX*UNIYLD*COS(frangl)/(root3*(THREE+SIN(frangl)))
endif
!
!       effective viscoplastic strain rate for hardening
if (stnhrd.gt.dfmin) then
    debar = BETB(1)**TWO + BETB(2)**TWO + HALF*BETB(3)**TWO
    debar = TWO * debar / THREE
    if (debar.gt.dfmin) debar = SQRT(debar)
    UNIYLD = UNIYLD + stnhrd*debar
endif
!
UNIYLD = MAX(UNIYLD,1d-3)
FNORM = PSTRS/UNIYLD - ONE
if (FNORM.le.0d0) FNORM = -1d0
!
return
end
!-------------------------------------------------------------------------------
subroutine FLOWVP (AVECT,STN,PRPLAS,PRPMAT,PSTRS,BETB, elemID)

USE ALGEBRA, only: CLEAR
use constants, only: one
!
! Computes BETB by multiplying AVECT with fluidity and PHI
!
implicit none

integer, parameter :: NSTR=3
!-pass
double precision :: pstrs
double precision :: AVECT(NSTR),STN(NSTR),PRPLAS(9),BETB(NSTR),PRPMAT(6)
integer          :: elemID
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
subroutine WINKLP (BTOT,DELD,IWELM,IWSIDE,IWTIME,WPRES,IEN,MAT,X,D,PRPMAT,LM, &
  LMF,TFAULT,SKEW,DXE,LMX,NSLIP,DX)
!
! Subroutine to compute Winkler restoring pressures from displacements (in local
! coordinate frame), and add them to the global force vector BTOT
!
USE MATERIALSMODULE, only : LMAT, MaterialIDOfElement
USE MODELDEFINITION
USE ALGEBRA
use constants,       only : zero, eps, half, one, useDXE
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact

implicit none
integer, parameter :: NEN=4, NSTR=3, NDOF=2, NSD=2

!-pass
integer          :: iwelm, iwside, iwtime, ien, mat, lm, lmf, lmx, nslip
double precision :: wpres, x, prpmat, btot, deld, d, tfault, skew, dxe, dx
dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*),MAT(*),X(NSD,*),PRPMAT(6,*), &
   LM(NDOF,NEN,*),BTOT(*),DELD(NDOF,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN,*),LMX(NDOF,NEN,*),NSLIP(5,*),DX(NDOF,*)

!-locl
double precision :: f(NDOF),fl(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),rot(3,3),tmp(NDOF), &
  anglevec(2),wt,area
integer :: lien(4),WinklerEntry,globalelementnumber,timecode,materialnumber,nodecounter, &
  globalnodenumber,idof,eqno
logical :: quad
!
if (NUMWNK.eq.0) then
    return
endif

do WinklerEntry=1,IABS(NUMWNK)

    globalelementnumber = IABS(IWELM(WinklerEntry))
    timecode = IWTIME(WinklerEntry)
    if (timecode.lt.0) then
        if (timecode.ne.-1.and.NSTEP.ge.-timecode) then
            cycle
        endif
    elseif (timecode.eq.0) then
        cycle
    else
        if (NSTEP.lt.(timecode-1)) then
            cycle
        endif
    endif

    call LCOORD (X,XL,IEN(1,globalelementnumber),0)
 
    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,globalelementnumber),NDOF,NEN)
        call AddFaultDisplacement (DL,globalelementnumber,TFAULT,NDOF,NEN,15)
        call ADDSNE (DL,globalelementnumber,NDOF,NEN,useDXE)
        call REZONE (XL,DL,ONE)
    endif

!   fill LIEN with local side numbers of element side IWSIDE(WinklerEntry)
    quad = IEN(3,globalelementnumber).ne.IEN(4,globalelementnumber)
    call SIDENP(IABS(IWSIDE(WinklerEntry)),quad,lien)
    
!   compute side length or surface area into variable AREA
    call FCGEOM (XL,lien,area,xn,"pln, addwnk")

    if (IWSIDE(WinklerEntry).ge.1 .and. IWSIDE(WinklerEntry).le.4) then
!       Winkler pressure perpendicular to one of element sides 1-4
        materialnumber = MaterialIDOfElement(globalelementnumber)
        wt = PRPMAT(6,materialnumber)*area/2D0    ! element thickness * side length / 2
        f(1) = ABS(xn(1)) * WPRES(WinklerEntry) * wt
        f(2) = ABS(xn(2)) * WPRES(WinklerEntry) * wt

!       compute local nodal forces FL from global nodal forces F for first two nodes
        do nodecounter=1,2
            globalnodenumber = IEN(lien(nodecounter),globalelementnumber)
            anglevec(1) = -SKEW(1,globalnodenumber)
            anglevec(2) = 0d0
    
            if (ABS(anglevec(1)).gt.EPS) then
                call FORMRT (anglevec,rot,NDOF)
                call EQUATE (tmp,f,NDOF)
                call VMPRD (rot,tmp,fl,NDOF,3)
            else
                call EQUATE (fl,f,NDOF)
            endif

!           add local forces to load vector 
            do idof=1,NDOF
                eqno = LM(idof,lien(nodecounter),globalelementnumber)
                if (eqno.gt.0) then
                    BTOT(eqno) = BTOT(eqno) - fl(idof)*DELD(idof,globalnodenumber)
                endif
            enddo
        enddo
    endif
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine SHAREL (IEN,IJSUR,NUMEL,NEL)

!     Finds elements neighbouring to element NEL.
!    Neighbouring elements are identified by two shared nodes.

use debugmodule, only: xit
use iomodule,        only: stderr

implicit none

integer, parameter :: NSUR=4,NEN=4
!-pass
integer NEL,NUMEL, ijsur, IEN
dimension IEN(NEN,*),IJSUR(NSUR,*)
!-locl
integer, external :: NINDX
integer :: n, nn, ielm, m, ii, mm, isur
integer :: nshare

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

use debugmodule, only: xit
use iomodule,        only: stderr

implicit none

integer, parameter :: NSUR=4
!-pass
integer            :: NEL, ijsur
dimension IJSUR(NSUR,*)
!-loc
integer            :: i

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
! this is subroutine is never called.

USE TIMESTEPMODULE

USE CONSTANTS, only: one, two, three, four, six, &
    root3, deg2rad, half,dfmin,dlogmin,dfmax,dlogmax,eps

!    Routine to determine visco-plastic fluidity from
!    explicit time stepping limit (Cormeau) (LINEAR only)
!     (Ref.: Owen & Hinton, sec. 8.9)
!    ONLY WORKS WITH NEWTONIAN RELAXATION!

 implicit none
integer, parameter :: NSTR=3
!-pass
double precision :: PRPMAT(6),PRPLAS(9),BETB(NSTR)
!-local
integer :: nptype, nflow
double precision :: yldpwr, frangl, uniyld, stnhrd, fluidy, e, debar
double precision :: sinphi

!-init
NPTYPE = INT(PRPLAS(1))
if (NPTYPE.eq.0) return
 YLDPWR = PRPLAS(5)
 NFLOW  = INT(PRPLAS(6))
 if (NFLOW.ne.1 .or. ABS(YLDPWR-1d0).gt.dfmin) return
 frangl = PRPLAS(2) * deg2rad
uniyld = PRPLAS(3)
 STNHRD = PRPLAS(7)
 FLUIDY = PRPLAS(9)
 E      = PRPMAT(1)

 DEBAR = BETB(1)**TWO + BETB(2)**TWO + HALF*BETB(3)**TWO
 DEBAR = TWO * DEBAR / THREE
 if (DEBAR.gt.dfmin) DEBAR = SQRT(DEBAR)

 if (NPTYPE.eq.1) then
!        Tresca
continue
 else if (NPTYPE.eq.2) then
!        Von Mises
continue
 else if (NPTYPE.eq.3) then
!        Mohr-Coulomb
uniyld = uniyld * COS(frangl)
 else if (NPTYPE.eq.4) then
!        Drucker-Prager
uniyld = SIX * uniyld * COS(frangl) / &
       ( root3 * (THREE - SIN(frangl)) )
 endif
 if (STNHRD.gt.dfmin) uniyld = uniyld + STNHRD*DEBAR

 if (NPTYPE.eq.1) then
!        Tresca
FLUIDY = 1.5D0 * uniyld / (E * DELTP)
 else if (NPTYPE.eq.2) then
!        Von Mises
FLUIDY = FOUR * 1.5D0 * uniyld / (THREE * E * DELTP)
 else if (NPTYPE.eq.3) then
!        Mohr-Coulomb
frangl = PRPLAS(2)
SINPHI = SIN(frangl*deg2rad)
SINPHI = SINPHI*SINPHI
FLUIDY = FOUR*uniyld /(E*DELTP)
 else if (NPTYPE.eq.4) then
!        Drucker-Prager
continue
 endif
 if (FLUIDY.gt.dfmin) PRPLAS(4) = FLUIDY

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


! in case of tiangles, the shape function has 
! 0-entries in SH(:,4), so they can freely be
! added to the BP matrix without effectively
! counting the last node twice.

do i=1,NDOF
    do j=1,NSD
        BP(i,j) = dot(SH(j,:), DL(i,:),nen)
    enddo
enddo

 return
 end
!-------------------------------------------------------------------------------
subroutine sdimnt (ISURF,SURF,B,ID,GRAV,PRPSED,DELTP)

USE MODELDEFINITION
use constants, only: two
use debugmodule, only: xit
use iomodule,        only: stderr

! Routine to calculate load increment resulting from sedimentary
! transport in following time step. The surface uplift SURF is
! updated. It is assumed that the first derivative of the surface
! elevation at the first and last points is zero. Other assumption
! is that the NDOF spatial component is vertical.

implicit none

integer, parameter :: NSD=2,NDOF=2
!-pass
double precision :: deltp
integer :: isurf, id
double precision :: surf, b, grav, prpsed
 dimension ISURF(*),SURF(NSD,*),B(*),ID(NDOF,*),GRAV(*), &
   PRPSED(2)
!-locl
integer, parameter :: IX=1, IY=2
integer i,node,k
double precision RhoSed,Diffus,Him1,Dxim1,Hi,Hip1,Dxi, &
   Qim1,Qi,Qip1,Dhi,Area,Dfi

#ifdef DEBUG
write(stderr,10) NSURF
10   format(1x,'sdimnt: NSURF = ',I3)
#endif

if (NSURF.le.0 .or. NSED.eq.0) then
    return
endif

!    If NSURF <= 2 first derivative assumptions are violated
if (NSURF.le.2) then
    write(stderr,20)
20       format(1x,'sdimnt: NSURF needs to be greater than 2')
    call xit(1," ")
endif

RhoSed = PRPSED(1)
Diffus = PRPSED(2)

#ifdef DEBUG 
    write(stderr,30) RhoSed,Diffus
30   format(1x,'sdimnt: RhoSed = ',F6.1,'  Diffus = ',1PE12.2)
#endif

Him1  = two*SURF(IY,1) - SURF(IY,2)
Dxim1 = SURF(IX,2) - SURF(IX,1)

do i=1,NSURF

    Hi = SURF(IY,i)

    if (i.eq.NSURF) then
        Hip1 = two * SURF(IY,i) - Him1
        Dxi  = Dxim1
    else
        Hip1 = SURF(IY,i+1)
        Dxi  = SURF(IX,i+1) - SURF(IX,i)
    endif

    Qim1 = two/(Dxim1*(Dxi+Dxim1))
    Qi   = -two/(Dxi*Dxim1)
    Qip1 = two/(Dxi*(Dxi+Dxim1))
    Dhi = DELTP*Diffus*(Qim1*Him1+Qi*Hi+Qip1*Hip1)
    !        Update surface
    SURF(NSD,i) = SURF(NSD,i) + Dhi
    !        Add loads
    Area = ABS(Dxim1+Dxi)/two
    Dfi = Area * GRAV(IY) * RhoSed * Dhi
    node = ISURF(i)
    k = ID(IY,node)
    if (k.gt.0) then
        B(k) = B(k) + Dfi
    endif

#ifdef DEBUG
    write(stderr,40) i,Dxim1*1d-3,Dxi*1d-3
40           format(/1x,'sdimnt: i = ',I3/1x,'sdimnt: Dx(i-1)= ', &
         F6.2,' km',4x,'Dx(i)= ',F6.2,' km')
    write(stderr,50) Him1,Hi,Hip1
50           format(1x,'sdimnt: H(i-1)= ',F6.1,' m', &
         4x,'H(i)= ',F6.1,' m',5x,'H(i+1)= ',F6.1,' m')
    write(stderr,60) Qim1,Qi,Qip1
60           format(1x,'sdimnt: Q(i-1)= ',1PG12.4, &
         4x,'Q(i)= ',1PG12.4,5x,'Q(i+1)= ',1PG12.4)
    write(stderr,70) Dhi,Dfi
70           format(1x,'sdimnt: dHi = ',F7.2,' m',5x,'dFi = ', &
         1PG12.4,' N')
#endif

!        save for next loop
    Him1 = Hi
    Dxim1 = Dxi
enddo

return
end subroutine
!-----------------------------------------------------------------------
subroutine SIDENP (ISint,QUAD,LIEN)

use debugmodule, only: xit
use iomodule,        only: stderr

!
! returns local node point numbers in LIEN, belonging to side IS
! nodal points in LIEN are positively oriented.
!
implicit none
!-pass
integer ISint,LIEN(4)
logical QUAD

LIEN(3) = 0
LIEN(4) = 0

if (ISint.eq.1 .or. ISint.eq.2) then
    LIEN(1) = ISint
    LIEN(2) = ISint + 1
else if (ISint.eq.3 .or. ISint.eq.4) then
    if (QUAD) then
        LIEN(1) = ISint
        LIEN(2) = ISint + 1
        if (LIEN(2).gt.4) then
            LIEN(2) = 1
        endif
    else
        LIEN(1) = 3
        LIEN(2) = 1
    endif
elseif (ISint.eq.5 .or. ISint.eq.6) then
    LIEN(1) = 1
    LIEN(2) = 2
    LIEN(3) = 3
    if (QUAD) then
        LIEN(4) = 4
    endif
else
    write(stderr,1) ISint
    1       format(1x,'SIDENP: quadrilaterals have 6 sides, not ',I2)
    call xit(1," ")
endif

return
end
!-----------------------------------------------------------------------
subroutine NPSIDE (NP,IEN,ISINT)

! returns side number 

integer, parameter:: NEN=4
integer :: NP(3),IEN(NEN),ISINT

ISINT = 0
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

return

end subroutine
!-----------------------------------------------------------------------
function SRFHGT (SURF,NSURF,X)
!
! Determines surface height at X from surface nodal point data by search and
! interpolate.
!
USE ALGEBRA
use iomodule,        only: stderr

implicit none
integer, parameter :: NSD=2
!-pass
double precision :: x
integer          :: nsurf
double precision :: surf
dimension SURF(NSD,*)

#include "coord.i"
!-locl
 logical FIRST
 save FIRST
 data FIRST /.true./
double precision :: srfhgt
double precision :: dx1, dx2
integer          :: i,k,m
double precision :: a,b
!-init
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
end
!-----------------------------------------------------------------------
! routine, and routine that calls this one (elwinx) is never used
subroutine towinkx (NSLIP,IEN,WINX,IWINX,X,nel,iside,mode,f, &
                    NWINKX)

USE MODELDEFINITION, only: NUMEL, NUMSLP
use constants,       only: eps
use debugmodule,     only: xit
use iomodule,        only: stderr

!
!    Geometry-dependent part of elwinx, where stresses are converted
!    to forces. Updates IWINX, WINX and NWINKX.
!
 implicit none
integer, parameter :: NEN=4,NDOF=2,NSD=2
!-pass
integer :: iside, nwinkx
integer :: nslip, ien, iwinx, nel, mode
double precision :: winx, x, f
 dimension NSLIP(5,*),IEN(NEN,*),WINX(NDOF,*),IWINX(NDOF,*), &
   X(NSD,*),f(*)
!-locl
double precision :: XL(NSD,NEN),DW(2)
integer          :: IIEN(2), islp(2), np, i, k
double precision :: dx, dy
!
!    Some preliminary checks
if (iside.lt.1 .or. iside.gt.4) then
    write(stderr,1) nel,iside
 1       format(1x,'elwinx: illegal side number for element ',I12, &
     ': iside=',I5)
    call xit(1," ")
endif

if (nel.lt.1 .or. nel.gt.NUMEL) then
    write(stderr,2) NUMEL,nel
 2       format(1x,'elwinx: element number outside range [1, ',I6, &
     ']: ',I6)
    call xit(1," ")
endif
!
!    Check that both face nodes of side iside occur in NSLIP
IIEN(1) = iside
IIEN(2) = iside + 1

if (IEN(3,nel).eq.IEN(4,nel)) then
    if (IIEN(2).gt.3) then
        IIEN(2)=1
    endif
else
    if (IIEN(2).gt.NEN) then
        IIEN(2)=1
    endif
endif

do i=1,2
    np = IEN(IIEN(i),nel)
    islp(i) = 1

100 if (NSLIP(1,islp(i)).eq.nel.and.NSLIP(2,islp(i)).eq.np) then
        cycle
    endif

    islp(i) = islp(i) + 1

    if (islp(i).gt.NUMSLP) then
        write(stderr,3) nel,np,iside
 3           format(1x,'elwinx: element-node ',I5,'-',I5,1x, &
         'not in NSLIP array (iside=',i1,')')
        call xit(1," ")
    endif

    goto 100
enddo

call LCOORD (X,XL,IEN(1,nel),0)

dx = (XL(1,IIEN(1))-XL(1,IIEN(2)))*0.5d0
dy = (XL(2,IIEN(1))-XL(2,IIEN(2)))*0.5d0

DW(1) = ABS(dy*f(1)) + ABS(dx*f(2))
DW(2) = ABS(dx*f(1)) + ABS(dy*f(2))

do i=1,2
    np = IEN(IIEN(i),nel)
    do k=1,NSD

        if (DW(k).le.EPS) then
            cycle
        endif

        if (NSLIP(2+k,islp(i)).eq.0) then
            cycle
        endif

        WINX(k,np) = WINX(k,np) + DW(k)

        if (IWINX(k,np).eq.0) then
            NWINKX = NWINKX + 1
            IWINX(k,np) = mode
            cycle
        else if (IWINX(k,np).eq.mode) then
            cycle
        endif

        write(stderr,4) nel,np,iside,k
 4           format(1x,'elwinx: conflicting MODE for element ', &
         I5,', node ',I5,' (iside=',I1,', ',I1, &
         '-component)')
        call xit(1," ")
    enddo
enddo
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine SlipperySkew (IDSLP,SKEW,NSELSD,IEN,X,D,TFAULT,DXE,LMF)

! handles fault parallel elements

USE MODELDEFINITION  ! only: NUMSN and probably more
use meshdatamodule,  only: meshdatactx
USE CONSTANTS,       only: zero, &
                           pi, &
                           rad2deg
use debugmodule,     only: lustat, &
                           itest, &
                           iecho, &
                           parallellog, &
                           delay, &
                           debug, &
                           debugFileID
use modelctx,        only: getrank
use iomodule,        only: stderr
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact

#ifdef SPARSE
use AOmodule,        only: vertices, &
                           nvertices
#endif

!
! Updates skew angles for slippery nodes
!
implicit none
integer, parameter :: NSD=2, NDOF=2, NSTR=3, NEN=4

!-pass

double precision :: SKEW, D, tfault, dxe
integer :: IDSLP, NSELSD, IEN, lmf

dimension :: IDSLP(*),SKEW(2,meshdatactx%Nvglobal),NSELSD(2,NSLSKEW)
dimension :: IEN(NEN,meshdatactx%nelocal)
dimension :: D(NDOF,*),TFAULT(NDOF,*),DXE(NDOF,NEN,*),LMF(NDOF,NEN,*)

#ifdef SPARSE
    ! nvertices only available in SPARSE mode.
    double precision :: X(NSD,nvertices)
#else
    ! no SPARSE = single processor run; all vertices in the X array
    double precision :: X(NSD,meshdatactx%Nvglobal)
#endif

!-locl
logical :: firstentry
save    :: firstentry,lutmp,nwait
integer :: lien(4)
integer :: nels
double precision :: angle
integer :: nwait
integer :: k, kk, n, nn, lutmp
integer, external :: nextlu
double precision :: oldskew
double precision :: a
logical          :: quad
logical          :: ElementInThisPartition
integer          :: checkIndex
!-init
data firstentry,nwait /.true.,50/
a = 0d0
!

if (LGDEF.eq.0 .or. NSLSKEW.eq.0 .or. NUMROT.eq.0) then
    return
endif

if (iecho.eq.8) then
    call parallelLog("SlipperySkew","entered")
endif

if (firstentry) then
    lutmp = nextlu(10)
    call openf(lutmp,'slskew.tmp','unknown')
    firstentry = .false.
endif

! loop over the global number of slippery nodes
! to make sure to also get the ghost points.

! This is *not* the number of slippery node entries
! in the input file, as every slippery node will have 
! multiple entries, one for each of the elements
! bordering on it.

do kk=1,NUMSN

#ifdef DEBUG

!    call delay()

    write(*,*) "**********************************************************"
    write(*,*) "## rank",getrank(),"doing slippery node",kk,"of total", NUMSN, "########"
    write(*,*) "**********************************************************"


    write(*,*) "*** rank ", getrank(), "get kk entry", kk, "of IDSLP ", IDSLP(1:NUMSN)
#endif

    ! k is the local number of the slippery node
    k = IDSLP(kk)



    if (k.eq.0) then
        ! k 0, this node is not relevant for this partition
#ifdef DEBUG
        write(*,*) "k = 0, returning"
#endif
        cycle
#ifdef DEBUG
    else
        write(*,*) "k = ", k, "continuing"
#endif
    endif


    if (SKEW(1,k).eq.ZERO) then
#ifdef DEBUG
        write(*,*) 'Skew 0, returning'
#endif
        cycle
#ifdef DEBUG
    else
        ! the node has an Euler angle.
        ! continue and possibly adapt change it, depending on below.
        write(*,*) "rank", getrank(),"skew: (1,",k,") is ", skew(1,k), "continuing"
#endif
    endif


    angle = ZERO
    nels  = 0

    ! loop over the fault parallel entries in TECIN.DAT
    ! for every element bordering on this node, adapt the angle


    do nn=1,NSLSKEW ! global number of fault parallel entries, so that foreign fault parallel elements
                    ! are also included in the contribution to the slippery node

#ifdef DEBUG
        write(*,*) "-----------------------------------------------------------------------------"
        write(*,*) "## rank",getrank(),"doing parallel surface entry",nn,"of total", NSLSKEW, "########"
        ! n = element number of this entry
        ! NSELSD contains fault par elems in local element numbering 
        ! (to match numbering of IEN which is local)
#endif
        n = NSELSD(1,nn)

#ifdef DEBUG
        ! we need to find the points of this element.
        ! There are several options:

        ! 1) This element is in this partition, and all its points are 
        ! in his partition also.

        ! 2) This element is in this partition, but it has points outside
        ! this partition

        ! 3) This element is outside this partition, but it has a point,
        ! (At least the slippery point, but potentially others, that are in this partition)

        ! First we need to determine whether the element is in this partition
        ! The index must appear in locelt, then
    
        ! TODO, this can be sped up by a lookup table, at the cost of memory. 
        ! Maybe not necessary    
!        ElementInThisPartition = .false.
!        checkIndex = 0
!        do while (.not. ElementInThisPartition .and. checkIndex .lt. meshdatactx%nelocal) 
!            checkIndex = checkIndex + 1
!            if (meshdatactx%gloElt(checkIndex) .eq. n) then
!                ElementInThisPartition = .true.
!            endif
!        enddo

        ! possiblity 1 and 2 are essentially the same, because 
        ! the partitions can access the coordinates of their ghost points.
        ! which are included in the X array.
        ! However, we must find out IEN entry that belongs to this element.
        ! The elementID is local index, which runs from 1 to the global
        ! number of elements. However, the IEN arrays only store vertices
        ! for elements of this particular partition.

!        if (ElementInThisPartition) then
!            n = checkIndex ! our previous search stopped there.
!        else
!            ! we are in for a world of pain
!        endif
        

!        if (n.gt.0) then

!        write(*,*) "checks entry",nn, "of NSELSD", NSELSD, "which is n", n

        

        !n = n - meshdatactx%gloElt(1) + 1

!        write(*,*) "for this partitions count, that becomes ", n
#endif
        ! sidenp finds the local indices of the vertices that
        ! are the endpoint of this slippery interface, and puts them in the LIEN array
        ! They will be in [1,2,3], according to PLN convention.

        quad = IEN(3,n).ne.IEN(4,n)
        quad = .false.  ! because WAAAAAA
        call SIDENP (NSELSD(2,nn),quad,lien)



#ifdef DEBUG
!        write(*,*) "side ID: ", NSELSD(2,nn)

!        write(*,*) "sideNP says element ", n, &
!                   "lien", lien(1), lien(2), &
!                   "and IEN", meshdatactx%IEN(lien(1),n), meshdatactx%IEN(lien(2),n), &
!                   "must match", k

!        write(*,*) "1 try entry ", n, " out of ", size(meshdatactx%IEN,2)

        ! compare this to local number:

!        write(*,*) "2 try entry ", n - meshdatactx%lowestLocalVertexID, &
!                    " out of ", size(meshdatactx%IEN,2)


        ! if one of those end points is the slippery node k, of the outer loop.
        ! 2D, so we only have to validate two points to constrain the slippery line.

!!!
!        if (debug) then
!        write(*,*) "rank",getrank()," checks points ", lien(1), "and" ,lien(2), "of elem", k

!        write(*,*) "rank",getrank(), &
!                   "slippery node ", k, &
!                   "could be pos", n, &
!                   "which is", IEN(lien(1),n), &
!                   "or", IEN(lien(2),n)
!        endif
#endif

        ! IEN contains the numbering according to the partition itself,
        ! indices running from 1 to the number of relevant nodes for this partition.

        ! n is the global number as given in the TECIN.DAT
        ! we can translate

        if (IEN(lien(1),n).eq.k .or. &
            IEN(lien(2),n).eq.k) then


!            if (LIEN(1).eq.LIEN(2)) then
!                write(stderr,*) "Vulnerability detected; SlipperySkew calling tangent called"
!               write(stderr,*) "with two identical vertices", LIEN(1)
!                write(stderr,*) "Contact model support"
!                stop "Halting executing"
!            else if (X(1,LIEN(1)).eq.X(1,LIEN(2))  .and.  X(2,LIEN(1)).eq.X(2,LIEN(2))) then
!                write(stderr,*) "Vulnerability detected; SlipperySkew calling tangent with "
!                write(stderr,*) "two points with the same coordinates"
!                write(stderr,*) "Point",LIEN(1), "has coordinates", X(1,LIEN(1)), X(2,LIEN(1))
!                write(stderr,*) "Point",LIEN(2), "has coordinates", X(1,LIEN(2)), X(2,LIEN(2))
!                write(stderr,*) "Contact model support"
!                stop "Halting executing"
!            endif

!            write(*,*) "calling tangent", lien(1), lien(2)
!            write(*,*) "calling tangent", n, "of", meshdatactx%nvlocal


            call TANGENT (X,D,TFAULT,DXE(1,1,n), &
                          IEN(1,n),n,SKEW,lien,a)
#ifdef DEBUG

            write(*,*) "rank",getrank(),"angle after tangent", a
#endif
            nels = nels + 1

            if (nels.gt.1 .and. ABS(angle-a).gt.0.9*PI) then
                if (angle.gt.a) then
                    a = a + PI
                else
                    a = a - PI
                endif
            endif

#ifdef DEBUG
            write(*,*) "rank",getrank(),"angle after pi correction", a
#endif

            ! what is this? averaging angle over the elements involved?
            angle = (angle*DBLE(nels-1) + a) / DBLE(nels)

#ifdef DEBUG
            write(*,*) "rank",getrank(),"angle after averaging", a
#endif

!!!
        endif ! if vertix is points on fault parallel interface 
!!!

    enddo  ! loop over NSLSKEW

    if (nels.eq.0) then
        cycle
    endif

    if (DEBUG) then
        oldskew = SKEW(1,k)
    endif

#ifdef DEBUG
   write(*,*) "rank", getrank(), "updates angle ", k, "from ", SKEW(1,k), "to", angle
#endif

    SKEW(1,k) = angle
!

     if (DEBUG) then
         if (NSTEP.eq.0) then
             write(stderr,1) k,SKEW(1,k)*rad2DEG
 1              format(1x,'SLSKEW: SKEW(',I5,') initialized at ', &
            F16.12,' deg')
        else
               write(stderr,2) k,oldskew*rad2DEG,SKEW(1,k)*rad2DEG
 2               format(1x,'SLSKEW: SKEW(',I5,') updated from ', &
             F16.12,' to ',F16.12,' deg')
        endif
    endif


    if (MOD(NSTEP,nwait).eq.0) then
        write(lutmp,3) NSTEP,k, SKEW(1,k)*rad2DEG
 3    format(1x,'NSTEP = ',I6,'   SKEW(1,',I6,') = ',F16.12, ' deg')
    endif

enddo

!
if (MOD(NSTEP,nwait).eq.0) then 
    call flushf(lutmp)
endif


if (iecho.eq.8) then
    write(*,*) "rank", getrank(), "leaving SLSKEW"
endif


return
end
!-----------------------------------------------------------------------
subroutine TANGENT (X,D,TFAULT,DXE,IEN,elemID,SKEW,LIEN,ANGLE)
!
! returns angle of tangent of deformed element side
!
USE MODELDEFINITION
USE CONSTANTS,       only: halfpi, one,eps, useDXE
use modelctx,        only: getrank
use iomodule,        only: stderr

implicit none
integer, parameter :: NSD=2, NDOF=2, NSTR=3, NEN=4
!-pass
double precision :: X, D, SKEW, TFAULT, DXE
dimension        :: X(NSD,*),D(NDOF,*),SKEW(2,*), &
                    TFAULT(NDOF,*), DXE(NDOF,*)
double precision :: ANGLE
integer          :: elemID
integer          :: IEN(NEN), LIEN(4)
!-locl
double precision :: xl(NSD,NEN),dl(NDOF,NEN),xn(NSD)
double precision :: area


! localize coordinate data
call LCOORD (X,xl,IEN,0)

! deformed coordinates XL
call LDISP  (dl,D,IEN,NDOF,NEN)
call AddFaultDisplacement (dl,elemID,TFAULT,NDOF,NEN,21)
call ADDSNE (dl,elemID,NDOF,NEN,useDXE)
call REZONE (xl,dl,ONE)
call FCGEOM (xl,LIEN,area,xn,"pln, tangent")

if (ABS(xn(2)).gt.EPS) then
!   returns angle in range -pi/2 to pi/2
    ANGLE = ATAN(xn(1)/xn(2))
else
    ! wrap angle
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
!
! Program to localize coordinate data
!
USE MODELDEFINITION, only : ISHELL

implicit none

integer, parameter :: NSD=2, NDOF=2, NSTR=3, NEN=4
!-pass
double precision   :: X, XL
dimension          :: X(NSD,*),XL(NSD,NEN)
integer            :: IEN(NEN)
integer            :: elemID
logical            :: quad

XL(1,1) = X(1,IEN(1))
XL(2,1) = X(2,IEN(1))
XL(1,2) = X(1,IEN(2))
XL(2,2) = X(2,IEN(2))
XL(1,3) = X(1,IEN(3))
XL(2,3) = X(2,IEN(3))
XL(1,4) = X(1,IEN(4))
XL(2,4) = X(2,IEN(4))

if (ISHELL.ne.0) then
    quad = IEN(3).ne.IEN(4)
    call SHLCRD (XL, quad)
endif

return
end subroutine
!-----------------------------------------------------------------------
subroutine SHLCRD (XL,QUAD)

USE MODELDEFINITION, only : ISHELL, RADIUS
USE CONSTANTS, only: deg2rad, zero

implicit none

integer, parameter :: NEN=4,NSD=2
!-pass
double precision   :: XL
dimension          :: XL(NSD,*)

logical            :: QUAD

!-locl
integer i
double precision   :: xlon,ylat,xc,yc,zc
double precision   :: v(3,NEN),north(3),east(3)

! compute easting and northing vectors in element center (in global cartesian
! frame)
if (QUAD) then
    xlon = (XL(1,1)+XL(1,2)+XL(1,3)+XL(1,4))*2.5D-1
    ylat = (XL(2,1)+XL(2,2)+XL(2,3)+XL(2,4))*2.5D-1
else
    xlon = (XL(1,1)+XL(1,2)+XL(1,3)) * 0.333333333
    ylat = (XL(2,1)+XL(2,2)+XL(2,3)) * 0.333333333
endif

xlon = xlon * DEG2RAD
ylat = ylat * DEG2RAD

east(1) = -DSIN(xlon)
east(2) = DCOS(xlon)
east(3) = ZERO

north(1) = -DSIN(ylat)*DCOS(xlon)
north(2) = -DSIN(ylat)*DSIN(xlon)
north(3) = DCOS(ylat)

! compute vectors to nodes in global cartesian frame
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

! Compute node coordinates relative to element center
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
    xc = (v(1,1)+v(1,2)+v(1,3))*0.33333333
    yc = (v(2,1)+v(2,2)+v(2,3))*0.33333333
    zc = (v(3,1)+v(3,2)+v(3,3))*0.33333333
    do i=1,3
        v(1,i) = (v(1,i)-xc)*RADIUS
        v(2,i) = (v(2,i)-yc)*RADIUS
        v(3,i) = (v(3,i)-zc)*RADIUS
    enddo
endif

! Project node coordinates on easting and northing directions
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
end subroutine

!---------------------------------------------------------------

double precision function dilate (strain,POIS)

USE MODELDEFINITION, only: iopt, NOCOMPR
use constants,       only: one, two, third
!
! computes volume change or volume change rate, depending on what S is.
!
implicit none
integer, parameter :: NSTR=3
!-pass
double precision :: strain(NSTR)
double precision :: pois

if (NOCOMPR.ne.0) then
   dilate = 0d0
else
    if (IOPT.eq.0) then
!       plane strain
        dilate = (strain(1)+strain(2))*THIRD
    else
!       plane stress
        dilate = (ONE-TWO*POIS)/(ONE-POIS)*(strain(1)+strain(2))*THIRD
    endif
endif

return
end
!-------------------------------------------------------------------------------
double precision function EFFEE (S,POIS,SIGNINFO)

USE MODELDEFINITION, only: iopt, nocompr
use constants,       only: dfmin,zero,third,one,two, eps, half
!
! Compute effective deviatoric strain
!
implicit none
integer, parameter :: NSTR=3

!-pass
logical SIGNINFO
double precision :: S(NSTR)
double precision :: pois
!-locl
double precision :: J2
double precision :: dil
double precision :: exx, eyy, ezz
!-init
 effee = ZERO
!
if (IOPT.eq.0) then
!        plane strain
   if (NOCOMPR.eq.0) then
       dil = (S(1)+S(2))*THIRD
   else
       dil = 0d0
   endif

   exx = S(1)-dil
   eyy = S(2)-dil
   ezz =     -dil

else
!        plane stress
   if (NOCOMPR.eq.0) then
        dil = (ONE+TWO*POIS)/(ONE-POIS)*(S(1)+S(2))*THIRD
        exx = S(1)-dil
        eyy = S(2)-dil
        ezz = -THIRD*(ONE+POIS)/(ONE-POIS)*(S(1)+S(2))
   else
        exx = S(1)
        eyy = S(2)
        ezz = 0d0
   endif
endif

J2 = HALF*(exx*exx+eyy*eyy+ezz*ezz)+S(3)*S(3)

if (J2.gt.dfmin) then
    effee = SQRT(J2)
    if (SIGNINFO .and. (S(1)+S(2)).lt.-EPS) then
        effee = -effee
    endif
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine STRDEV (S,POIS)

USE MODELDEFINITION, only: iopt, nocompr
use constants,       only: third, one, two
!
! computes strain (rate) deviator
!
implicit none
integer, parameter :: NSTR=3
!-pass
double precision :: S(NSTR)
double precision :: pois
!-locl
double precision :: dil
!
if (NOCOMPR.ne.0) then
    return
endif
!
if (IOPT.eq.0) then
    !        plane strain
    dil = (S(1)+S(2))*THIRD
    S(1) = S(1)-dil
    S(2) = S(2)-dil
else
    !        plane stress
    dil = (ONE+TWO*POIS)/(ONE-POIS)*(S(1)+S(2))*THIRD
    S(1) = S(1)-dil
    S(2) = S(2)-dil
endif
!
return
end
!-------------------------------------------------------------------------------
double precision function EFFSTN (S,POIS,SIGNINFO)

use constants, only: zero, dfmin, eps

!
! computes effective deviatoric stress. S is the full stress tensor.
!
implicit none
integer, parameter :: NSTR=3
!-pass
logical          :: SIGNINFO
double precision :: S(NSTR)
double precision :: pois
!-locl
double precision :: J2
double precision, external :: AJ2
!
J2 = AJ2(S,POIS)

if (J2.gt.dfmin) then
    EFFSTN = SQRT(J2)
    if (SIGNINFO .and. (S(1)+S(2)).lt.-EPS) then
        effstn = -effstn
    endif
else
    EFFSTN = ZERO
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine STNTEN (STN,POIS,SL)

USE MODELDEFINITION, only: IOPT, nocompr
USE ALGEBRA,         only: CLEAR

implicit none
integer, parameter :: NSTR=3
!-pass
double precision :: STN(NSTR)
double precision :: SL(3,3)
double precision :: pois

call CLEAR(SL,9,"SL")

if (IOPT.eq.0) then
    !        plane strain
    SL(1,1) = STN(1)
    SL(2,2) = STN(2)
    if (NOCOMPR.eq.0) then
        SL(3,3) = POIS*(STN(1)+STN(2))
    endif
    SL(1,2) = STN(3)
    SL(2,1) = STN(3)
else
    !        plane stress
    SL(1,1) = STN(1)
    SL(2,2) = STN(2)
    SL(1,2) = STN(3)
    SL(2,1) = STN(3)
    if (NOCOMPR.ne.0) then
        SL(3,3) = -STN(1)-STN(2)
    endif
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine STNDEV (S,POIS)

USE MODELDEFINITION, only: iopt, nocompr
use constants,       only: third, two, one

!
implicit none
integer, parameter :: NSTR=3
!-pas2s
double precision :: S(NSTR)
double precision :: pois
!-locl
double precision :: sxx,syy,p

if (NOCOMPR.ne.0) then
    return
endif

!
if (IOPT.eq.0) then
   !        plane strain
   p = (ONE+POIS)*(S(1)+S(2))*THIRD
   S(1) = S(1) - p
   S(2) = S(2) - p
else
   !        plane stress
   sxx = (TWO*S(1)-S(2))*THIRD
   syy = (TWO*S(2)-S(1))*THIRD
   S(1) = sxx
   S(2) = syy
endif
!
return
end
!-------------------------------------------------------------------------------
double precision function PRESSURE (S,POIS)

USE MODELDEFINITION
use constants, only: third, one
!
implicit none
integer, parameter :: NSTR=3
!-pass
double precision :: S(NSTR)
double precision :: pois
!-locl

if (NOCOMPR.ne.0) then
    pressure = 0d0
else
    if (IOPT.eq.0) then
!            plane strain
        pressure = (ONE+POIS)*(S(1)+S(2))*THIRD
    else
!            plane stress
        pressure = (S(1)+S(2))*THIRD
    endif
endif
!
return
end
!-------------------------------------------------------------------------------
#ifdef SHEARZONE_TRACTIONS_MARIO_VERSION

subroutine SHEARTRAX (B,BTOT,X,D,IDX,IEN,TFAULT,NSLIP,DX,DT,shearETA,shearTHICK,ELNOSHEARZONE,PRPMAT,shearnpl,shearQpl,shearTEMP)

! Updates the load vector B for shear tractions on a shear zone with viscosity shearETA and thickness shearTHICK.
! The shear zone is parameterized by a slippery node interface. 
! *** ASSUMPTION: the interface is parallel to the local x-direction (see Euler angles) ***
! The current implementation is based on a linear (Newtonian) dependence of shear traction on differential velocities DX.
! Citation: D'Acquisto et al. (in prep)

USE MODELDEFINITION
use constants,      only: useDXE, ZERO, ONE, TWO, THREE, FOUR
USE MATERIALSMODULE,only: MaterialIDOfElement
use modelctx,       only: getrank
use iomodule,       only: stderr
use debugmodule

implicit none
    
integer, parameter          ::  NEN=4, NDOF=2, NSD=2
                                                                    ! Below, partition-local includes ghost numbers:
! pass
integer                     ::  IDX(NDOF,*),&                   ! partition-local equation numbers for each partition-local node number
                                                                ! with length NUMNP=partition-local number of nodes. However,
                                                                ! modeldatamodule.F allocates length numnpGlobal!
                                IEN(NEN,*),&                    ! partition-local node numbers for each partition-local element number.
                                NSLIP(5,*),&                    ! partition-local element and node numbers and associated slippery dofs.
                                                                ! Length NUMSLP is the partition-local number of slippery node entries, i.e.,
                                                                ! elements are local but nodes may not be. 
                                ELNOSHEARZONE(NUMELSHEARZONE)   ! list of partition-local element numbers along the shear zone(s)
double precision            ::  B(*),&                          ! incremental load for each partition-local equation number
                                BTOT(*),&                       ! cumulative load for each partition-local equation number
                                DX(NDOF,*),&                    ! differential velocity for each partition-local slippery node number
                                DT,&                            ! Timestep size
                                shearETA(NUMELSHEARZONE),&      ! shear zone viscosity for given local element numbers
                                shearTHICK(NUMELSHEARZONE),&    ! shear zone thickness for given local element numbers
                                X(NSD,*),&                      ! coordinates of partition-local node numbers
                                D(NDOF,*),&                     ! displacements of partition-local node numbers
                                TFAULT(NDOF,*),&                ! faulted node displacements of partition-local node numbers. Length NUMFN
                                                                ! is partition local NUMFN
                                PRPMAT(6,*)                     ! material properties of partition-local element numbers
! local
integer                     ::  ielm,ElmPlocal, islpn, counter, nodesnp(2), lien(4), &
                                NPloc1, NPloc2, LocalEq, inodeno, numnodes, nodenoarr(NUMELSHEARZONE*2), NParr(NUMELSHEARZONE)
!integer                    ::  ielm2, nodesnp_all_asc(2,NUMELSHEARZONE)
logical                     ::  newnodeno
double precision            ::  xl(NSD,NEN), dl(NDOF,NEN), xn(NSD), area, dF
integer                     ::  shearnpl                        ! power-law coefficient
double precision            ::  shearQpl(NUMELSHEARZONE)        ! power-law activation energy (J/mol)
double precision            ::  shearTEMP(NUMELSHEARZONE)       ! Shear-zone temperature for each element
double precision, parameter ::  Rpl=8.31446261815324            ! Universal gas constant (for power law)

if (NUMELSHEARZONE.le.0) then
    return
endif

if (debug) then
    if (getrank().eq.0) write(*,*) getrank(),shearnpl
5   format("Rank",i3," SHEARTRAX: shearnpl=",I1)
end if

! Create a list of shear zone nodes. Loop over shear zone elements
numnodes = 0 ! Reset counter of number of unique node numbers

loop_shear_zone_elements: do ielm = 1,NUMELSHEARZONE

    ElmPlocal = ELNOSHEARZONE(ielm)
    islpn = 1
    counter = 1
    nodesnp = 0
    do while (counter.le.2)                   ! Set interface node numbers in the current element
        if (NSLIP(1,islpn).eq.ElmPlocal) then
            nodesnp(counter) = NSLIP(2,islpn)   ! set partition-local node number of the overall mechanical problem
            counter = counter+1
        endif
        islpn = islpn+1
    enddo

    ! Get element-local node numbers defining the interface
    lien = 0
    do counter = 1,2
        do islpn = 1,NEN
            if (nodesnp(counter).eq.IEN(islpn,ElmPLocal)) then
                lien(counter) = islpn
                exit
            endif
        enddo
    enddo
                    
    call LCOORD (X,xl,IEN(1,ElmPLocal),ElmPlocal) ! Copy coordinates of the current element into xl
    if (LGDEF.ne.0) then                          ! Lagrangian update of node coordinates
        call LDISP  (dl,D,IEN(1,ElmPLocal),NDOF,NEN)
        call AddFaultDisplacement (dl,ElmPlocal,TFAULT,NDOF,NEN,1005)
        call ADDSNE (dl,ElmPlocal,NDOF,NEN,useDXE)
        call REZONE(xl,dl,ONE)
    endif
        
    ! Calculate surface area and unit normal vector of the interface
    call FCGEOM (xl,lien,area,xn,"SHEARTRAX")
    area = area * PRPMAT(6,MaterialIDOfElement(ElmPlocal))

    NPloc1 = nodesnp(1)
    NPloc2 = nodesnp(2)
    if (DT.gt.0d0) then
        if (debug) then
            write(*,30) getrank(),DT,NPloc1,DX(1,NPloc1)/DT,NPloc2,DX(1,NPloc2)/DT
30          format("Rank",I3," SHEARTRAX: DT=",1PG12.2,"   DX(1,node=",I12,")/DT=",1PG12.4, &
             " and DX(1,node=",i12,")=",1PG12.4)
        end if
        if (shearnpl.gt.0) then
            dF = -(TWO*shearETA(ielm) / shearTHICK(ielm))**(ONE/shearnpl) * (THREE*SIGN(ABS(DX(1,NPloc1)/DT)**&
                (ONE/shearnpl),DX(1,NPloc1)/DT)+SIGN(ABS(DX(1,NPloc2)/DT)**(ONE/shearnpl),DX(1,NPloc2)/DT)) / FOUR * &
                exp(shearQpl(ielm)/Rpl/shearTEMP(ielm)/shearnpl) * area / TWO
        else
            dF = -(THREE*DX(1,NPloc1)/DT+DX(1,NPloc2)/DT) * area * shearETA(ielm) / (FOUR*shearTHICK(ielm))
        end if
    else
        dF = ZERO
    end if

    if (debug) then
        write(*,40) getrank(),ElmPlocal,NPloc1,dF
40      format("Rank",i3," SHEARTRAX: shear force on element=",i12," node=",i12," DF=",1PG12.4)
    endif

! notyet
!   dF = SIGN(1d0,dF) * MIN(dFmax,dF)

    if (debug) then
        write(*,42) getrank(),dF
42      format("Rank",i3," SHEARTRAX: maxed DF=",1PG12.4)
    endif

    LocalEq = IDX(1,NPloc1)   ! equation number that is local to the parallel partition
    if (LocalEq.gt.0) then
        B(LocalEq) =  B(LocalEq) + dF
        if (debug) then
            write(*,50) getrank(),NPloc1,LocalEq,LocalEq,B(LocalEq)
50          format("Rank",i3," SHEARTRAX: IDX(1,node=",i12,")=",i12,". Updated B(",i12,")=",1PG12.4)
        endif
        ! Save list of (unique) node numbers used for forces
        newnodeno = .true.
        do inodeno = 1,numnodes
            if (nodenoarr(inodeno).eq.NPloc1) then
                newnodeno = .false.
            end if
        end do
        if (newnodeno) then
            numnodes = numnodes+1
            nodenoarr(numnodes) = NPloc1
        end if
    else
        write(stderr,*) "SHEARTRAX SURPRISE: no equation corresponding with local slippery node number ",NPloc2
    end if
    if (DT.gt.0) then
        if (shearnpl.gt.0) then
            dF = -(TWO*shearETA(ielm) / shearTHICK(ielm))**(ONE/shearnpl) * (SIGN(ABS(DX(1,NPloc1)/DT)**&
                (ONE/shearnpl),DX(1,NPloc1)/DT)+THREE*SIGN(ABS(DX(1,NPloc2)/DT)**(ONE/shearnpl),DX(1,NPloc2)/DT)) / &
                FOUR * exp(shearQpl(ielm)/Rpl/shearTEMP(ielm)/shearnpl) * area/TWO
        else
            dF = -(DX(1,NPloc1)/DT+THREE*DX(1,NPloc2)/DT) * area * shearETA(ielm) / (FOUR*shearTHICK(ielm))
        end if
    else
        dF = ZERO
    end if
    if (debug) then
        write(*,40) getrank(),ElmPlocal,NPloc2,dF
    end if
    LocalEq = IDX(1,NPloc2)   ! equation number that is local to the parallel partition
    if (LocalEq.gt.0) then
        B(LocalEq) =  B(LocalEq) + dF
        if (debug) then
            write(*,50) getrank(),NPloc2,LocalEq,LocalEq,B(LocalEq)
        endif
        ! Save list of (unique) node numbers used for forces
        newnodeno = .true.
        do inodeno = 1,numnodes
            if (nodenoarr(inodeno).eq.NPloc2) then
                newnodeno = .false.
            end if
        end do
        if (newnodeno) then
            numnodes = numnodes+1
            nodenoarr(numnodes) = NPloc2
        end if
    else
        write(stderr,*) "SHEARTRAX SURPRISE: no equation corresponding with local slippery node number ",NPloc2
    end if
enddo loop_shear_zone_elements

! Subtract BTOT entries from shear forces B entries

loop_shear_zone_nodes: do inodeno = 1,numnodes

        LocalEq = IDX(1,nodenoarr(inodeno))
        B(LocalEq) = B(LocalEq) - BTOT(LocalEq)
        if (debug) then
            write(*,60) getrank(),nodenoarr(inodeno),LocalEq,BTOT(LocalEq),localEq,B(LocalEq)
60          format("Rank",i3," SHEARTRAX: For dof=1 of node ",I12,", subtracting BTOT(",i12,")=", &
                1PG12.4," yields B(",i12,")=",1PG12.4) 
        endif

enddo loop_shear_zone_nodes

if (debug) then
    write(*,*) "Rank", getrank()," End of the SHEARTRAX subroutine"
end if

end subroutine

#endif
!--------------------------------------------------------------------------------
#ifdef SHEARZONE_TRACTIONS

subroutine SHEARTRAX (IFLAG,NSTEP)

! IFLAG=0: set differential Winkler forces to facilitate the calculation of the maximum shear forces
! IFLAG=1: compute maximum shear forces from the differential Winkler forces (from IFLAG=0)
! IFLAG=2: set nodal forces equivalent to shear traction but make them lower than the maximum shear forces

use iomodule          

! Routine for checking if implicit code needs to be invoked because
! Differential Winkler forces are being changed.

use CONSTANTS,      only    :   useDXE,ZERO,ONE,TWO,THREE,FOUR
use IOMODULE,       only    :   luin, stderr
use MODELDATAMODULE,only    :   modeldatactx
use MESHDATAMODULE, only    :   meshdatactx
use MATERIALSMODULE,only    :   matdatactx,MaterialIDOfElement
use DEBUGMODULE,    only    :   iecho, debug
use MODELCTX,       only    :   getrank,getsize
use MODELDEFINITION,only    :   NUMELSHEARZONE, NUMNPSHEARZONE,NTYPE,LGDEF
use ALGEBRA,        only    :   arraycontains
use TIMESTEPMODULE, only    :   DELTP

implicit none
integer,parameter :: NEN=4,NSD=2,NDOF=2
!pass
integer           :: IFLAG,NSTEP
!locl
logical           :: quad,applied,firstentry,Node1InPart,Node2InPart
integer           :: n,elno,m,side,lien(4),i,node,idof,ientry
integer           :: mysize,myrank,maxdof,lOCALeQ
double precision  :: E,pois,ShearModulus,xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),area, &
                     HalfWinklerConstant,shearmag,Bmag,velo(NDOF,2),dF(NDOF,2)
save              :: firstentry
!init
data firstentry /.true./

if (NUMELSHEARZONE.eq.0) return 

#ifdef ZISFAULTNORMAL
    maxdof=NDOF-1
#else
    maxdof=NDOF
#endif

mysize = 1
myrank = 1
#ifdef SPARSE
    mysize = getsize()
    myrank = getrank()
#endif
    
applied = .false.       ! just for echoing purposes

#ifdef SPARSE
if (firstentry) then
    ! convert global node numbers in the shear zone list to partition-local numbers
    do n=1,NUMNPSHEARZONE
        if (mysize.gt.1) then
            node = modeldatactx%NPNOSHEARZONE(n)        ! global node number
            if (meshdatactx%v2p(node).ne.myrank) then
                node = -1
            else
                node = meshdatactx%Glo2LocVertex(node)  ! partition-local node number
            endif
            modeldatactx%NPNOSHEARZONE(n) = node
        else
            exit
        endif
    enddo
endif
#endif
firstentry = .false.

if (IFLAG.eq.0) then
    ! set differential Winkler forces to facilitate the calculation of the maximum shear forces

    loop_shear_zone_elements: do n=1,NUMELSHEARZONE
        elno = modeldatactx%ELNOSHEARZONE(n)            ! global element number
#ifdef SPARSE
        if (mysize.gt.1) then
            elno = meshdatactx%Glo2LocElement(elno)     ! set ELNO to partition-local element number
            if (elno.eq.-1) cycle loop_shear_zone_elements
        endif
#endif
        m = MaterialIDOfElement(elno)
    
        call LCOORD (meshdatactx%X,xl,meshdatactx%IEN(1,elno),0)
        if (LGDEF.ne.0) then
            call LDISP  (dl,modeldatactx%D,meshdatactx%IEN(1,elno),NDOF,NEN)
            call AddFaultDisplacement (dl,elno,modeldatactx%TFAULT,NDOF,NEN,25)
            call ADDSNE (dl,elno,NDOF,NEN,useDXE)
            call REZONE (xl,dl,ONE)
        endif

        quad = meshdatactx%IEN(3,elno).ne.meshdatactx%IEN(4,elno)
        side = modeldatactx%SIDESHEARZONE(n)                    
        call SIDENP (side,quad,lien)                            ! get nodal point numbers of both nodes along the interface
        call FCGEOM (xl,lien,area,xn,'MaxShearTractionForces')  ! compute interface area and normal vector
        E = matdatactx%PRPMAT(1,m)                              ! Young's modulus
        pois = matdatactx%PRPMAT(2,m)                           ! Poisson's ratio
        ShearModulus = E/(TWO+TWO*pois)                         ! Use the elastic shear modulus as the differential Winkler constant
        HalfWinklerConstant = matdatactx%PRPMAT(6,m)*area*ShearModulus/TWO
HalfWinklerConstant=HalfWinklerConstant * 1d-5

        do i=1,2
            node = meshdatactx%IEN(lien(i),elno) ! note that this node can be outside the partition,
                                                 ! but that will be handled thru the equation numbers
            ! set timing and magnitude of the differential Winkler force
            ! on all dofs, will be handled by equation no's later 
            do idof=1,maxdof
#ifdef ZISFAULTNORMAL
                modeldatactx%IWINX(NDOF,node) = 0       ! phase shift
#endif
                modeldatactx%IWINX(idof,node) = NSTEP
                modeldatactx%WINX(idof,node) = modeldatactx%WINX(idof,node) + HalfWinklerConstant
                applied = .true.
            enddo
            if (debug) then
                write(*,10) myrank,IFLAG,elno,side,node,modeldatactx%IWINX(1,node), &
                    modeldatactx%WINX(1,node)
10              format('rank ',I3,' Sheartrax(IFLAG=',I1,'): elem=',I12,' side=',I1, &
                    ' node=',I12,' IWINX=',I5,' updated WINX=',1PG12.4)
            endif

        enddo
    enddo loop_shear_zone_elements

else if (IFLAG.eq.1) then

    ! (re)compute maximum equivalent nodal forces from differential Winkler ! forces
    modeldatactx%DFMAXSHEARZONE = ZERO

    loop_shear_zone_nodes: do n=1,NUMNPSHEARZONE
        node = modeldatactx%NPNOSHEARZONE(n)
        if (node.eq.-1) then
            cycle loop_shear_zone_nodes         ! not in this partition
        endif
        shearmag = ZERO
        do idof=1,maxdof
            if (modeldatactx%IWINX(idof,node).ne.0) then
                shearmag = shearmag + (modeldatactx%DX(idof,node)*modeldatactx%WINX(idof,node))**2
                applied = .true.
            endif
            modeldatactx%IWINX(idof,node) = 0
            modeldatactx%WINX(idof,node) = ZERO
        enddo
        modeldatactx%DFMAXSHEARZONE(n) = modeldatactx%DFMAXSHEARZONE(n) + SQRT(shearmag)
    enddo loop_shear_zone_nodes

    if (debug) then
        do n=1,NUMNPSHEARZONE
            node = modeldatactx%NPNOSHEARZONE(n)
            if (node.eq.-1) then
                cycle
            endif
            write(*,20) myrank,IFLAG,node,modeldatactx%DFMAXSHEARZONE(n)
20          format('rank ',I3,' Sheartrax(IFLAG=',I1,'): node=',I12, &
             ' DFMAX=',1PG12.4)
        enddo
    endif

else if (IFLAG.eq.2 .and. DELTP.gt.ZERO) then

    ! impose nodal forces equivalent to shear traction but make them lower than the maximum shear forces

    loop_shear_elements: do n=1,NUMELSHEARZONE
        elno = modeldatactx%ELNOSHEARZONE(n)
#ifdef SPARSE
        if (mysize.gt.1) then
            elno = meshdatactx%Glo2LocElement(elno)
            if (elno.eq.-1) cycle loop_shear_elements
        endif
#endif
        m = MaterialIDOfElement(elno)
        quad = meshdatactx%IEN(3,elno).ne.meshdatactx%IEN(4,elno)
        side = modeldatactx%SIDESHEARZONE(n)
        call SIDENP (side,quad,lien)
    
        call LCOORD (meshdatactx%X,xl,meshdatactx%IEN(1,elno),0)
        if (LGDEF.ne.0) then
            call LDISP  (dl,modeldatactx%D,meshdatactx%IEN(1,elno),NDOF,NEN)
            call AddFaultDisplacement (dl,elno,modeldatactx%TFAULT,NDOF,NEN,25)
            call ADDSNE (dl,elno,NDOF,NEN,useDXE)
            call REZONE (xl,dl,ONE)
        endif

        call FCGEOM (xl,lien,area,xn,'MaxShearTractionForces')
        area = area*matdatactx%PRPMAT(6,m)

        ! compute velocity magnitude of first node (if in this partition)
        node = meshdatactx%IEN(lien(1),elno)
        call arraycontains(modeldatactx%NPNOSHEARZONE,NUMNPSHEARZONE,node,i)
        Node1InPart = (i.le.NUMNPSHEARZONE) 
        if (Node1InPart) then
            do idof=1,maxdof
                velo(idof,1) = modeldatactx%DX(idof,node)/DELTP
                if (debug) write(*,43) myrank,IFLAG,idof,node,velo(idof,1)
43              format('Rank',i3,' SHEARTRAX(IFLAG=',i1,'): velo(idof=',i1,',node=',i12,')=',1PG12.4)
            enddo
        endif

        ! compute velocity magnitude of second node (if in this partition)
        node = meshdatactx%IEN(lien(2),elno)
        call arraycontains(modeldatactx%NPNOSHEARZONE,NUMNPSHEARZONE,node,i)
        Node2InPart = (i.le.NUMNPSHEARZONE) 
        if (Node2InPart) then
            do idof=1,maxdof
                velo(idof,2) = modeldatactx%DX(idof,node)/DELTP
                if (debug) write(*,43) myrank,IFLAG,idof,node,velo(idof,2)
            enddo
        endif

        ! equivalent shear force on the nodes that constitute the interface
        do idof=1,maxdof
            if (Node1InPart .and. Node2InPart) then
                dF(idof,1) = -(velo(idof,1)*THREE + velo(idof,2))/FOUR * &
                               area * modeldatactx%shearETA(n) / modeldatactx%shearTHICK(n)
                dF(idof,2) = -(velo(idof,2)*THREE + velo(idof,1))/FOUR * &
                               area * modeldatactx%shearETA(n) / modeldatactx%shearTHICK(n)
            else if ((.not.Node1InPart) .and. Node2InPart) then
                dF(idof,2) = -velo(idof,2) * area * modeldatactx%shearETA(n) / modeldatactx%shearTHICK(n)
                dF(idof,1) = dF(idof,2) ! won't be applied because not in partition
            else if (Node1InPart .and. (.not.Node2InPart)) then
                dF(idof,1) = -velo(idof,1) * area * modeldatactx%shearETA(n) / modeldatactx%shearTHICK(n)
                dF(idof,2) = dF(idof,1) ! won't be applied because not in partition
            else
                dF(idof,1) = ZERO
                dF(idof,2) = ZERO
            endif
            if (debug) then
                write(*,45) myrank,idof,dF(idof,1),idof,dF(idof,2)
45              format('Rank',i3,' SHEARTRAX: dF(idof=',I1,',node1)=',1PG12.4/ &
                                   4x,3x,12x,'dF(idof=',I1,',node2)=',1PG12.4)
            endif
        enddo

        ! add equivalent forces to load vector
        node = meshdatactx%IEN(lien(1),elno)
        do idof=1,maxdof
            LocalEq = modeldatactx%IDX(idof,node)
            if (LocalEq.gt.0) then
                modeldatactx%B(LocalEq) =  modeldatactx%B(LocalEq) + dF(idof,1)
                if (debug) then
                    write(*,50) myrank,idof,node,LocalEq,LocalEq,modeldatactx%B(LocalEq)
50                  format('Rank',i3,' SHEARTRAX: IDX(idof=',I1,',node=',i12,')=',i12, &
                     '. Updated B(',i12,')=',1PG12.4)
                endif
            endif
        enddo
        node = meshdatactx%IEN(lien(2),elno)
        do idof=1,maxdof
            LocalEq = modeldatactx%IDX(idof,node)
            if (LocalEq.gt.0) then
                modeldatactx%B(LocalEq) =  modeldatactx%B(LocalEq) + dF(idof,2)
                if (debug) then
                    write(*,50) myrank,idof,node,LocalEq,LocalEq,modeldatactx%B(LocalEq)
                endif
            endif
        enddo

    enddo loop_shear_elements

    ! compute force increment
    compute_force_increment: do n=1,NUMNPSHEARZONE
        node = modeldatactx%NPNOSHEARZONE(n)
        do idof=1,maxdof
            LocalEq = modeldatactx%IDX(idof,node)
            if (LocalEq.gt.0) then
                modeldatactx%B(LocalEq) = modeldatactx%B(LocalEq) - modeldatactx%BTOT(LocalEq)
                if (debug) then
                    write(*,60) myrank,node,LocalEq,modeldatactx%BTOT(LocalEq),localEq, &
                        modeldatactx%B(LocalEq)
60                  format('Rank',i3,' SHEARTRAX: of node ',I12,', subtracting BTOT(', &
                        i12,')=',1PG12.4,' yields B(',i12,')=',1PG12.4) 
                endif
            endif
        enddo
    enddo compute_force_increment

    ! impose upper limit on force increment
    curb_force_update: do n=1,NUMNPSHEARZONE
        node = modeldatactx%NPNOSHEARZONE(n)
        Bmag = ZERO
        do idof=1,maxdof
            LocalEq = modeldatactx%IDX(idof,node)
            if (LocalEq.gt.0) then
                Bmag = Bmag + modeldatactx%B(LocalEq)**2
            endif
        enddo
        Bmag = SQRT(Bmag)
        do idof=1,maxdof
            LocalEq = modeldatactx%IDX(idof,node)
            if (LocalEq.gt.0 .and. Bmag.gt.ZERO) then
                if (Bmag.gt.modeldatactx%DFMAXSHEARZONE(n)) then
                    modeldatactx%B(LocalEq) = modeldatactx%B(LocalEq) * &
                        modeldatactx%DFMAXSHEARZONE(n) / Bmag
                    if (debug) then
                        write(*,70) myrank,node,LocalEq,modeldatactx%B(LocalEq)
70                      format('Rank',i3,' SHEARTRAX: of node ',I12,', capped B(', &
                        i12,')=',1PG12.4) 
                    endif
                endif
            endif
        enddo
    enddo curb_force_update
    applied = .true.
endif

if (iecho.eq.1 .and. applied) then
    if (IFLAG.eq.0) then
        write(*,100) getrank(),NSTEP
100     format('Rank ',I3,': Computing maximum tractions on shear zone in time step ',I6)
    else if (IFLAG.eq.1) then
        write(*,110) getrank(),NSTEP
110     format('Rank ',I3,': Maximum shear tractions determined in time step ',i6)
    else if (IFLAG.eq.2) then
        write(*,120) getrank(),NSTEP
120     format('Rank ',I3,': Shear tractions applied for time step ',i6)
    endif
endif

return
end subroutine

#endif /* SHEARZONE_TRACTIONS */
