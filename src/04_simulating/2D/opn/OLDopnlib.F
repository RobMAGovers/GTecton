!       Version 2009.2
!
!       These routines are the dimensionality-dependent portion of GTECTON
!       The routines here are for out-of-plane displacements
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
! Collapse of nodal points 4 onto 3 gives a triangle-element.
!-------------------------------------------------------------------------------
!       THE ROUTINES LISTED HERE ARE:
!
!       NAME                    FUNCTION
!
!
!       ADDFLX          NORMAL BOUNDARY FLUX CONTRIBUTIONS TO LOAD VECTOR
!       ADDPR           TRANSLATES APPLIED PRESSURES INTO LOAD VECTOR
!       ADDSTR          STRESS LOAD COMPUTATION
!       ADDWNK          ADDS WINKLER PRESSURES TO STIFFNESS MATRIX
!       AJ1             FIRST STRESS TENSOR INVARIANT
!       AJ2             SECOND STRESS TENSOR INVARIANT
!       AJ3             THIRD STRESS TENSOR INVARIANT
!       BDELD           STRAIN FROM DISPLACEMENTS
!       BMATRIX         STRAIN-DISPLACEMENT MATRIX
!       CPSTIFF         LOCAL HEAT CAPACITY STIFFNESS MATRIX
!       DEVSTN          STRESS DEVIATOR CALCULATION
!       EFORCE          EFFECTIVE LOAD VECTOR FOR VISCOUS FLOW
!       ELCRD           CALCULATES COORDINATES OF DEFORMED ELEMENT
!       ELCTR           CALCULATES CENTER COORDINATES OF DEFORMED ELEMENT
!       EIGMAX          CALCULATES MAXIMUM EIGENVALUE
!       FLOWVP          VISCOPLASTIC STRAIN RATE
!       FORMBT          VISCOUS FLOW CORRECTION FOR STRAIN
!       FORMJB          JACOBEAN MATRIX FOR VISCOUS FLOW
!       GRAVLD          GRAVITATIONAL BODY FORCE COMPUTATION
!       INTCON          DIMENSIONAL AND TOPOLOGIC DATA
!       YIELDS          COMPUTES SCALAR STRESS FOR COMPARISON WITH YIELD STRESS
!       LCOORD          LOCALIZES COORDINATES
!       LFLUX           ELEMENT HEAT FLOW
!       LHEAT           HEAT PRODUCTION CONTRIBUTION TO LOAD VECTOR
!       FLUID           DETERMINES VISCOPLASTIC FLUIDITY
!       MATERL          MATERIAL MATRICES
!       REZONE          LOCAL REZONE FOR LARGE DEFORMATION
!       SDIMNT          SURFACE ELEVATION AND LOADS FROM SEDIMENTARY TRANSPORT
!       SHAP20          SHAPE FUNCTION FOR 2-D ELEMENTS
!       SHAREL          IDENTIFIES NEIGHBOURING ELEMENTS
!       SIDENP          RETURNS LOCAL NODE NUMBERS ON A SPECIFIC SIDE
!       STIFF           STIFFNESS MATRIX
!       TANGVP          VISCOPLASTIC JACOBIAN MATRIX
!       TCENTR          COMPUTES ELEMENT CENTER TEMPERATURES
!       THGRAD          THERMAL GRADIENTS
!       TSTIFF          LOCAL CONDUCTION STIFFNESS MATRIX
!       WINKLP          WINKLER PRESSURES CONTRIBUTION TO GLOBAL FORCE VECTOR
!       WAPPLY          APPLIES RESIDUAL WINKLER FORCES
!       YIELDF          CONSTRUCTS PLASTIC FLOW VECTOR
!-------------------------------------------------------------------------------
subroutine INTCON ()

USE LABELMODULE
USE MODELDEFINITION, only: NTYPE, NSUR, NEE, NEC, simVariety
use modeltopology,   only: nSidesPerElem, NSD,NDOF,NSTR,NEN
use constants,       only: isOpn

!
! Set variables that depend on topology
!
! SNGL  implicit real (a-h,o-z)
implicit none

simVariety = isOpn

nSidesPerElem = 3
NTYPE = 2
NSD   = 2
NDOF  = 3
NSTR  = 5
NEN   = 4
NEC   = NSD * NEN
NEE   = NDOF * NEN
LABELC(1)='  X '
LABELC(2)='  Y '
LABELD(1)=LABELC(1)
LABELD(2)=LABELC(2)
LABELD(3)='  Z '
LABELS(1)='SXX '
LABELS(2)='SYY '
LABELS(3)='SXY '
LABELS(4)='SXZ '
LABELS(5)='SYZ '
LABELS(6)='    '
LABELR(1)='RXX '
LABELR(2)='RYY '
LABELR(3)='RXY '
LABELR(4)='RXZ '
LABELR(5)='RYZ '
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
!


 return
 end
!------------------------------------------------------------------------
 subroutine SHAP20 (S,T,X,DET,SH,XS,QUAD,IERR,IFLAG)

use constants, only: zero, half
!
! Program to compute shape functions for quadrilateral or triangle
!
! IFLAG                = 1, calculate s-, t- derivatives of shape
!                                 functions and shape functions
!                               = 2, calculate s-, t- derivatives,
!                                 shape functions and Jacobian
!                               = 3, calculate x-, y- derivatives of shape
!                                 functions, shape functions and Jacobian &
! S,T                 = natural coordinates
! X(NSD,NEN)                    = global coordinates
! QUAD                          = TRUE if quadrilateral element
!                               = FALSE if triangle element
! SH(1,NEN)                     = x- or s- derivatives of shape functions
! SH(2,NEN)                     = y- or t- derivatives of shape functions
! SH(3,NEN)                     = shape functions
! XS(NSD,NSD)                    = Jacobian matrix
! DET                 = Jacobian matrix determinant
!
!
implicit none
integer, parameter :: NEN=4, NSD=2
!-pass
logical QUAD
integer ierr, iflag
double precision :: SH(3,NEN),X(NSD,NEN),XS(NSD,NSD)
double precision :: S,T,DET
!-locl
double precision :: temp
double precision :: SA(NEN),TA(NEN)
save SA,TA
integer :: i

!-init
data SA/-0.5d0, 0.5d0,0.5d0,-0.5d0/
data TA/-0.5d0,-0.5d0,0.5d0, 0.5d0/
ierr = 0
!
!       Calculate shape functions and s- and t-derivatives of shape functions
do i=1,NEN
    SH(3,i) = (HALF+SA(i)*S)*(HALF+TA(i)*T)
    SH(1,i) = SA(i)*(HALF+TA(i)*T)
    SH(2,i) = TA(i)*(HALF+SA(i)*S)
enddo

if (.not.QUAD) then
!           adjust shape function of third node for triangular elements
    do i=1,3
        SH(i,3) = SH(i,3)+SH(i,4)
        SH(i,4) = ZERO
    enddo
endif
!
!       Jacobian matrix calculation
XS(1,1) = SH(1,1)*X(1,1) + &
          SH(1,2)*X(1,2) + &
          SH(1,3)*X(1,3) + &
          SH(1,4)*X(1,4)

XS(1,2) = SH(1,1)*X(2,1) + &
          SH(1,2)*X(2,2) + &
          SH(1,3)*X(2,3) + &
          SH(1,4)*X(2,4)

XS(2,1) = SH(2,1)*X(1,1) + &
          SH(2,2)*X(1,2) + &
          SH(2,3)*X(1,3) + &
          SH(2,4)*X(1,4)

XS(2,2) = SH(2,1)*X(2,1) + &
          SH(2,2)*X(2,2) + &
          SH(2,3)*X(2,3) + &
          SH(2,4)*X(2,4)

DET = XS(1,1)*XS(2,2)-XS(1,2)*XS(2,1)

if (DET.le.zero) then
    go to 1000
endif

!
if (IFLAG.eq.1) then
    return
endif
!
!       Calculatate inverse Jacobian determinant (Cramers rule)
TEMP    =  XS(1,1)/DET
XS(1,1) =  XS(2,2)/DET
XS(2,2) =  TEMP
XS(1,2) = -XS(1,2)/DET
XS(2,1) = -XS(2,1)/DET
!
if (IFLAG.eq.2) then
    return
endif
!
!    compute x- and y-derivatives of shape functions
do i=1,NEN
    TEMP    = XS(1,1)*SH(1,i) + XS(1,2)*SH(2,i)
    SH(2,i) = XS(2,1)*SH(1,i) + XS(2,2)*SH(2,i)
    SH(1,i) = TEMP
enddo
!
return
!
1000   write(stderr,1001) DET

1001   format(///1x,'Shape function fails! Determinant is ',1PE20.4)

ierr = 1

return
end subroutine
!-------------------------------------------------------------------------------
SUBROUTINE ELMTEST (X,IEN)

USE MESHDATAMODULE
USE MODELDEFINITION
use constants, only: eps,zero, half
use debugmodule, only: iecho

!
! tests element definitions
!
implicit none
integer NEN,NSD
parameter (NEN=4,NSD=2)
!-pass
integer IEN
double precision X
dimension X(NSD,*),IEN(NEN,*)

!-locl
logical quad
integer i,k
double precision xl,area1,area2,area,a,b
dimension xl(NSD,NEN),a(NSD),b(NSD)
!
if (NUMEL.lt.1) then
    write(stderr,10)
10        format(1x,'NUMEL < 1')
    call exitp(1)
endif
!
do i=1,NUMEL
    k = IEN(1,i)

    if (k.lt.1 .or. k.gt.NUMNP) then
        write(stderr,20) i,k
20           format(1x,'element ',I6,': reference to node number ', &
         I6)
        call exitp(1)
    endif

    XL(1,1) = X(1,k)
    XL(2,1) = X(2,k)
    k = IEN(2,i)

    if (k.lt.1 .or. k.gt.NUMNP) then
        write(stderr,20) i,k
        call exitp(1)
    endif

    XL(1,2) = X(1,k)
    XL(2,2) = X(2,k)
    k = IEN(3,i)

    if (k.lt.1 .or. k.gt.NUMNP) then
        write(stderr,20) i,k
        call exitp(1)
    endif

    XL(1,3) = X(1,k)
    XL(2,3) = X(2,k)
    k = IEN(4,i)

    if (k.lt.1 .or. k.gt.NUMNP) then
        write(stderr,20) i,k
        call exitp(1)
    endif

    XL(1,4) = X(1,k)
    XL(2,4) = X(2,k)
    quad = (IEN(3,i).ne.IEN(4,i))
!
    a(1) = XL(1,2) - XL(1,1)
    a(2) = XL(2,2) - XL(2,1)
    b(1) = XL(1,3) - XL(1,2)
    b(2) = XL(2,3) - XL(2,2)
    area1 = HALF*(a(1)*b(2)-a(2)*b(1))

    if (quad) then
        a(1) = XL(1,4) - XL(1,3)
        a(2) = XL(2,4) - XL(2,3)
        b(1) = XL(1,1) - XL(1,4)
        b(2) = XL(2,1) - XL(2,4)
        area2 = HALF*(a(1)*b(2)-a(2)*b(1))
        area = area1 + area2
    else
        area = area1
    endif

    if (ABS(area).le.EPS) then
        write(stderr,30) i
30        format(1x,'degenerate element ',I6)
        call exitp(1)
    endif

    if (area.lt.ZERO) then
        write(stderr,40) i
40           format(1x,'element ',I6,' is negatively oriented (node', &
         ' order is wrong)')
        call exitp(1)
    endif

    if (ABS(area1).le.EPS) then
        write(stderr,50) i
50           format(1x,'WARNING: degenerate triangle 123 in element', &
         1X,I6)
        call exitp(1)
    endif

    if (area1.lt.ZERO) then
        write(stderr,60) i
60           format(1x,'WARNING: negatively oriented triangle 123', &
         ' in element ',I6)
        call exitp(1)
    endif

    if (quad .and. ABS(area2).le.EPS) then
        write(stderr,70) i
70           format(1x,'WARNING: degenerate triangle 341 in element', &
         1X,I6)
        call exitp(1)
    endif

    if (quad .and. area2.lt.ZERO) then
        write(stderr,80) i
80           format(1x,'WARNING: negatively oriented triangle 341', &
         ' in element ',I6)
        call exitp(1)
    endif
enddo
!
if (iecho.eq.2) then
    write(stderr,90)
90    format(1x,'Elements tested ok')
endif

return
end
!-----------------------------------------------------------------------
subroutine MATERL (DMAT,E,POIS)

USE ALGEBRA,    only: clear
use constants,  only: half, one, two
!
! Constructs the material matrix relating stress and strain.
!
implicit none

integer, parameter :: NSTR=5
!-pass
double precision :: DMAT(NSTR,NSTR)
double precision :: E, pois
!-local
double precision :: am, al

call CLEAR(DMAT,25,"DMAT")
!
!       Compute Lame parameters AM(=MU) and AL(=LAMBDA)
AM=E/(ONE+POIS)
AL=AM*POIS/(ONE-TWO*POIS)
AM=HALF*AM
DMAT(1,1)=TWO*AM+AL
DMAT(1,2)=AL
DMAT(2,1)=AL
DMAT(2,2)=DMAT(1,1)
DMAT(3,3)=AM
DMAT(4,4)=AM
DMAT(5,5)=AM
!
 return
 end
!-------------------------------------------------------------------------------
subroutine BDELD (XL,DL,EE,QUAD,ierr)

use constants, only: half

!
! Subroutine to compute strains from displacements in each element
!
implicit none
integer, parameter :: NDOF=3, NSD=2, NSTR=5, NEN=4
!-pass
integer ierr
double precision :: XL(NSD,NEN),DL(NDOF,NEN),EE(NSTR)
logical QUAD
!-locl
double precision :: BP(NDOF,NSD)
!
!       Calculate matrix of partial displacement derivatives
call BPMATRIX (XL,DL,QUAD,BP,ierr)

if (ierr.ne.0) then
    return
endif

EE(1)=BP(1,1)
EE(2)=BP(2,2)
EE(3)=BP(1,2)+BP(2,1)
EE(4)=BP(3,1)
EE(5)=BP(3,2)
!
return
end
!-------------------------------------------------------------------------------
subroutine BMATRIX (B,SH)

use constants, only: zero
!
! Computes the strain-displacement matrix B
!
implicit none
integer, parameter :: NEN=4, NSTR=5, NDOF=3
integer, parameter :: NEE=NDOF*NEN
!-pass
double precision :: B(NSTR,NEE),SH(3,NEN)
!-local
integer :: i,k

k = 0
do i=1,NEN
  k = k + 1
    B(1,k) = SH(1,i)
    B(2,k) = ZERO
    B(3,k) = SH(2,i)
    B(4,k) = ZERO
    B(5,k) = ZERO
  k = k + 1
    B(1,k) = ZERO
    B(2,k) = SH(2,i)
    B(3,k) = SH(1,i)
    B(4,k) = ZERO
    B(5,k) = ZERO
  k = k + 1
    B(1,k) = ZERO
    B(2,k) = ZERO
    B(3,k) = ZERO
    B(4,k) = SH(1,i)
    B(5,k) = SH(2,i)
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine GRAVLD (P,XL,GRAV,TH,WT,QUAD,ierr)

USE CONSTANTS, only: G

!
! Computes the contribution to the load from gravity body forces
!
implicit none

integer, parameter :: NINT=4, NEN=4, NSD=2, NDOF=3
integer, parameter :: NEE=NEN*NDOF
!-pass
integer ierr
logical QUAD
double precision :: P(NEE),XL(NSD,NEN),GRAV(NDOF)
double precision :: TH, WT
!-locl
double precision :: GRAV1, GRAV2, GRAV3, DET
double precision :: SG(NEN),TG(NEN),XS(NSD,NSD),SH(3,NEN)
save TG,SG
integer :: i, j

!-init
data TG/-1d0,-1d0,1d0,1d0/
data SG/-1d0,1d0,1d0,-1d0/
!
do i=1,NINT
    call SHAP20 (SG(i)*G,TG(i)*G,XL,DET,SH,XS,QUAD,ierr,1)
    if (ierr.ne.0) then
        return
    endif

    GRAV1 = GRAV(1) * DET * WT * TH
    GRAV2 = GRAV(2) * DET * WT * TH
    GRAV3 = GRAV(3) * DET * WT * TH

    do j=1,NEN
        P(3*j-2) = P(3*j-2) + GRAV1 * SH(3,j)
        P(3*j-1) = P(3*j-1) + GRAV2 * SH(3,j)
        P(3*j  ) = P(3*j  ) + GRAV3 * SH(3,j)
    enddo
enddo
!
return
end
!-----------------------------------------------------------------------
subroutine STIFF (DMAT,XL,S,TH,QUAD,ierr)
USE MODELDEFINITION
USE ALGEBRA
USE CONSTANTS, only: G
!
! Computes the super-diagonal part of the local stiffness matrix.
!
implicit none

integer, parameter :: NINT=4, NEN=4, NSD=2, NDOF=3, NSTR=5, NMODES=2
integer, parameter :: NEX=NMODES*NDOF
!-pass
integer ierr
logical QUAD
double precision ::  XL(NSD,NEN),DMAT(NSTR,NSTR),S(NEE,NEE), th
!-local
integer          :: IDSAA(NEX)
double precision :: SG(NEN),TG(NEN),XS(NSD,NSD),SH(3,NEN),B(NSTR,NEE), &
 SDAL(NEE,NEX),SAAL(NEX,NEX),SDA(NEE,NEX),SAA(NEX,NEX), &
 SA(NEE,NEE)
integer :: i,j,m

double precision :: db1, db2, db3, db4, db5
double precision :: sa1, sa2, sa3, sa4, sa5, sa6
double precision :: c, det

logical firstentry
save SG,TG,firstentry,IDSAA
!-init
data firstentry /.true./
data TG/-1d0,-1d0,1d0,1d0/
data SG/-1d0,1d0,1d0,-1d0/
!
if (INCOMP.ne.0 .and. QUAD) then
    if (firstentry) then
!               set up list of diagonal addresses of square array
        call IDFORM (IDSAA,NEX)
        firstentry = .false.
    endif
    call CLEAR (SDA,72,"SDA")
    call CLEAR (SAA,36,"SAA")
endif
!
do m=1,NINT
    call SHAP20 (SG(m)*G,TG(m)*G,XL,DET,SH,XS,QUAD,ierr,3)
    if (ierr.ne.0) then
        return
    endif
    call BMATRIX (B,SH)
    C = DET*TH
    do j=1,NEE
!               Compute column J of DB-matrix
        DB1 = DMAT(1,1)*B(1,j) + &
              DMAT(1,2)*B(2,j) + &
              DMAT(1,3)*B(3,j) + &
              DMAT(1,4)*B(4,j) + &
              DMAT(1,5)*B(5,j)

        DB2 = DMAT(2,1)*B(1,j) + &
              DMAT(2,2)*B(2,j) + &
              DMAT(2,3)*B(3,j) + &
              DMAT(2,4)*B(4,j) + &
              DMAT(2,5)*B(5,j)

        DB3 = DMAT(3,1)*B(1,j) + &
              DMAT(3,2)*B(2,j) + &
              DMAT(3,3)*B(3,j) + &
              DMAT(3,4)*B(4,j) + &
              DMAT(3,5)*B(5,j)

        DB4 = DMAT(4,1)*B(1,j) + &
              DMAT(4,2)*B(2,j) + &
              DMAT(4,3)*B(3,j) + &
              DMAT(4,4)*B(4,j) + &
              DMAT(4,5)*B(5,j)

        DB5 = DMAT(5,1)*B(1,j) + &
              DMAT(5,2)*B(2,j) + &
              DMAT(5,3)*B(3,j) + &
              DMAT(5,4)*B(4,j) + &
              DMAT(5,5)*B(5,j)

        do i=1,j
            S(i,j) = S(i,j)+ DET*(B(1,i)*DB1 + &
                                  B(2,i)*DB2 + &
                                  B(3,i)*DB3 + &
                                  B(4,i)*DB4 + &
                                  B(5,i)*DB5)
        enddo
    enddo
!           Incompatible modes
    if (INCOMP.ne.0.and.QUAD) then
        call INCSTFF(SG(m)*G,TG(m)*G,DMAT,XL,B,SAAL,SDAL,QUAD,ierr)
        if (ierr.ne.0) then
            return
        endif
        do j=1,NEX
            SDA(1,j) =SDA(1,j) +C*SDAL(1,j)
            SDA(2,j) =SDA(2,j) +C*SDAL(2,j)
            SDA(3,j) =SDA(3,j) +C*SDAL(3,j)
            SDA(4,j) =SDA(4,j) +C*SDAL(4,j)
            SDA(5,j) =SDA(5,j) +C*SDAL(5,j)
            SDA(6,j) =SDA(6,j) +C*SDAL(6,j)
            SDA(7,j) =SDA(7,j) +C*SDAL(7,j)
            SDA(8,j) =SDA(8,j) +C*SDAL(8,j)
            SDA(9,j) =SDA(9,j) +C*SDAL(9,j)
            SDA(10,j)=SDA(10,j)+C*SDAL(10,j)
            SDA(11,j)=SDA(11,j)+C*SDAL(11,j)
            SDA(12,j)=SDA(12,j)+C*SDAL(12,j)
            SAA(1,j) =SAA(1,j) +C*SAAL(1,j)
            SAA(2,j) =SAA(2,j) +C*SAAL(2,j)
            SAA(3,j) =SAA(3,j) +C*SAAL(3,j)
            SAA(5,j) =SAA(5,j) +C*SAAL(5,j)
            SAA(6,j) =SAA(6,j) +C*SAAL(6,j)
        enddo
    endif
enddo
!
if (INCOMP.ne.0.and.QUAD) then
!           Replace SAA by its inverse
    call INVMAT (SAA,IDSAA,NEX)
    do j=1,NEE
!               column j of inv(SAA)*SAD
        SA1 = &
          SAA(1,1)*SDA(j,1)+SAA(1,2)*SDA(j,2)+SAA(1,3)*SDA(j,3) &
         +SAA(1,4)*SDA(j,4)+SAA(1,5)*SDA(j,5)+SAA(1,6)*SDA(j,6)
        SA2 = &
          SAA(2,1)*SDA(j,1)+SAA(2,2)*SDA(j,2)+SAA(2,3)*SDA(j,3) &
         +SAA(2,4)*SDA(j,4)+SAA(2,5)*SDA(j,5)+SAA(2,6)*SDA(j,6)
        SA3 = &
          SAA(3,1)*SDA(j,1)+SAA(3,2)*SDA(j,2)+SAA(3,3)*SDA(j,3) &
         +SAA(3,4)*SDA(j,4)+SAA(3,5)*SDA(j,5)+SAA(3,6)*SDA(j,6)
        SA4 = &
          SAA(4,1)*SDA(j,1)+SAA(4,2)*SDA(j,2)+SAA(4,3)*SDA(j,3) &
         +SAA(4,4)*SDA(j,4)+SAA(4,5)*SDA(j,5)+SAA(4,6)*SDA(j,6)
        SA5 = &
          SAA(5,1)*SDA(j,1)+SAA(5,2)*SDA(j,2)+SAA(5,3)*SDA(j,3) &
         +SAA(5,4)*SDA(j,4)+SAA(5,5)*SDA(j,5)+SAA(5,6)*SDA(j,6)
        SA6 = &
          SAA(6,1)*SDA(j,1)+SAA(6,2)*SDA(j,2)+SAA(6,3)*SDA(j,3) &
         +SAA(6,4)*SDA(j,4)+SAA(6,5)*SDA(j,5)+SAA(6,6)*SDA(j,6)
!               SDA*inv(SAA)*SAD: matrix to subtract from S
        do i=1,j
            SA(i,j) = &
              SDA(i,1)*SA1+SDA(i,2)*SA2+SDA(i,3)*SA3 &
             +SDA(i,4)*SA4+SDA(i,5)*SA5+SDA(i,6)*SA6
        enddo
!               adjust S
        do i=1,j
            S(i,j) = S(i,j) - SA(i,j)
        enddo
    enddo
endif
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine INCSTFF (SI,TI,DMAT,XL,B,SAAL,SDAL,QUAD,ierr)

use constants, only: zero, one, two
!
!       Computes those parts of the local stiffness matrix which
!       are due to incompatible shape functions. Defines SDAL and SAAL
!
implicit none
integer, parameter :: NEN=4, NSD=2, NDOF=3, NSTR=5, NMODES=2
integer, parameter :: NEE=NEN*NDOF
integer, parameter :: NEX=NMODES*NDOF

!-pass
integer ierr
logical QUAD
double precision :: si, ti
double precision :: XL(NSD,NEN),DMAT(NSTR,NSTR),B(NSTR,NEE),SAAL(NEX,NEX)
double precision :: SDAL(NEE,NEX)

!-locl
integer          :: i,j,k
double precision :: XS(NSD,NSD),SH(3,NEN),SHE(3,2),BE(NSTR,NEX)
double precision :: dSHdS, dSHdT
double precision :: dbe1, dbe2, dbe3, dbe4, dbe5
double precision :: det
!
!       two extra shape functions and s- and t-derivatives
SHE(3,1) = ONE - SI*SI
SHE(1,1) = - TWO*SI
SHE(2,1) = ZERO
SHE(3,2) = ONE - TI*TI
SHE(1,2) = ZERO
SHE(2,2) = - TWO*TI
!
!    Compute Jacobian
call SHAP20 (ZERO,ZERO,XL,DET,SH,XS,QUAD,ierr,2)

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
k = 0
do i=1,NMODES
    k = k + 1
    BE(1,k) = SHE(1,i)
    BE(2,k) = ZERO
    BE(3,k) = SHE(2,i)
    BE(4,k) = ZERO
    BE(5,k) = ZERO
    k = k + 1
    BE(1,k) = ZERO
    BE(2,k) = SHE(2,i)
    BE(3,k) = SHE(1,i)
    BE(4,k) = ZERO
    BE(5,k) = ZERO
    k = k + 1
    BE(1,k) = ZERO
    BE(2,k) = ZERO
    BE(3,k) = ZERO
    BE(4,k) = SHE(1,i)
    BE(5,k) = SHE(2,i)
enddo
!
do j=1,NEX
   DBE1 = DMAT(1,1)*BE(1,j) + &
          DMAT(1,2)*BE(2,j) + &
          DMAT(1,3)*BE(3,j) + &
          DMAT(1,4)*BE(4,j) + &
          DMAT(1,5)*BE(5,j)

   DBE2 = DMAT(2,1)*BE(1,j) + &
          DMAT(2,2)*BE(2,j) + &
          DMAT(2,3)*BE(3,j) + &
          DMAT(2,4)*BE(4,j) + &
          DMAT(2,5)*BE(5,j)

   DBE3 = DMAT(3,1)*BE(1,j) + &
          DMAT(3,2)*BE(2,j) + &
          DMAT(3,3)*BE(3,j) + &
          DMAT(3,4)*BE(4,j) + &
          DMAT(3,5)*BE(5,j)

   DBE4 = DMAT(4,1)*BE(1,j) + &
          DMAT(4,2)*BE(2,j) + &
          DMAT(4,3)*BE(3,j) + &
          DMAT(4,4)*BE(4,j) + &
          DMAT(4,5)*BE(5,j)

   DBE5 = DMAT(5,1)*BE(1,j) + &
          DMAT(5,2)*BE(2,j) + &
          DMAT(5,3)*BE(3,j) + &
          DMAT(5,4)*BE(4,j) + &
          DMAT(5,5)*BE(5,j)

   do i=1,NEE
      SDAL(i,j) = B(1,i)*DBE1 + &
                  B(2,i)*DBE2 + &
                  B(3,i)*DBE3 + &
                  B(4,i)*DBE4 + &
                  B(5,i)*DBE5
   enddo

   do i=1,NEX
      SAAL(i,j) = BE(1,i)*DBE1 + &
                  BE(2,i)*DBE2 + &
                  BE(3,i)*DBE3 + &
                  BE(4,i)*DBE4 + &
                  BE(5,i)*DBE5
   enddo
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine ADDPR (B,X,D,IEN,LM,MAT,PRPMAT,IELNO,ISIDE,PRES, &
   LMF,TFAULT,SKEW,DXE)

USE MATERIALSMODULE, only : LMAT, &
                          MaterialIDOfElement
USE MODELDEFINITION
USE ALGEBRA,         only: formrt, equate, vmprd
use constants,       only: eps, half, one, useDXE

! Subroutine to add pressures to load vector

 implicit none

integer, parameter :: NEN=4, NSD=2, NDOF=3

!-pass
integer          :: ielno, mat, ien, iside, lm, lmf
double precision :: x, d, press, prpmat, b, tfault, skew, dxe
dimension IELNO(*),MAT(*),IEN(NEN,*),X(NSD,*),D(NDOF,*), &
 ISIDE(*),PRES(*),PRPMAT(6,*),B(*),LM(NDOF,NEN,*), &
 LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
 DXE(NDOF,NEN,*)

!-locl
double precision :: anglevec(2)
double precision :: F(NDOF),XL(NSD,NEN),DL(NDOF,NEN),xn(NSD),tmp(NSD), &
   rot(3,3)
integer          :: n,lien(2), m,ii,k
double precision :: w, angle, area, pres

do k=1,IABS(NUMPR)
    n = IELNO(k)
    call LCOORD (X,XL,IEN(1,n),0)

    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,24)
        call ADDSNE (DL,n,NDOF,NEN, useDXE)
        call REZONE(XL,DL,ONE)
    endif

    call SIDENP(ISIDE(k),IEN(3,n).ne.IEN(4,n),lien)
    call FCGEOM (XL,lien,area,xn)
    angle =(-SKEW(1,IEN(lien(1),n))-SKEW(1,IEN(lien(2),n)))*HALF

    anglevec(1) = angle
    anglevec(2) = 0d0

     if (ABS(angle).gt.EPS) then
        call FORMRT (anglevec,rot,NDOF)
           call EQUATE (tmp,xn,NDOF)
        call VMPRD (rot,tmp,xn,NDOF,3)
    endif

    m = MaterialIDOfElement(n)
    w = PRPMAT(6,m) * area * HALF

    F(1) = xn(1) * PRES(k) * w
    F(2) = xn(2) * PRES(k) * w

    ii = LM(1,lien(1),n)

    if (ii.ne.0) then
        B(ii) = B(ii) + F(1)
    endif

    ii = LM(2,lien(1),n)

    if (ii.ne.0) then
        B(ii) = B(ii) + F(2)
    endif

    ii = LM(1,lien(2),n)

    if (ii.ne.0) then
        B(ii) = B(ii) + F(1)
    endif

    ii = LM(2,lien(2),n)

    if (ii.ne.0) then
        B(ii) = B(ii) + F(2)
    endif

enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine FCGEOM (XL,LIEN,AREA,XN)

use constants, only: dfmin,dlogmin,dfmax,dlogmax, two

!
! Calculates area of element face and unit normal vector.
!
implicit none
integer, parameter :: NSD=2,NEN=4
!-pass
integer          :: LIEN(2)
double precision :: XL(NSD,NEN),XN(NSD)
!-locl
double precision :: area

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
!
return
end
!-------------------------------------------------------------------------------
subroutine LFLUX (XL,TL,IEN,PRPTEM,FLUX,ierr)
!
! Routine to calculate element heat flux
!
implicit none

integer, parameter :: NSD=2,NEN=4
!-pass
integer          :: IEN(NEN), ierr
double precision :: XL(NSD,NEN),TL(NEN),PRPTEM(6),FLUX(NSD)

!-locl
double precision :: sn, tn
double precision :: GRAD(NSD)
logical QUAD

sn = 0d0
tn = 0d0
!
QUAD = (IEN(3).ne.IEN(4))

if (.not.QUAD) then
    tn = -1d0/3d0
endif

! calculate thermal gradients in the center of the element
call thgrad (sn,tn,XL,TL,GRAD,QUAD,ierr)

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
USE MODELDEFINITION, only : IPOINT
!
! Computes the effective forces at each node
!
implicit none

integer, parameter :: NEN=4, NSD=2, NDOF=3, NSTR=5
integer, parameter :: NEE=NEN*NDOF
double precision, parameter :: GP=0.57735026918963d0
!-pass
integer ierr
logical QUAD
double precision :: th
double precision :: XL(NSD,NEN),P(NEE),EVP(NSTR)
!-locl
integer   :: NINT
double precision :: SG(NEN),TG(NEN),XS(NSD,NSD),SH(3,NEN),B(NSTR,NEE)
save SG,TG
integer :: i, m
double precision :: W, G,DET
!-init
data SG/-1d0, 1d0, 1d0,-1d0/
data TG/-1d0,-1d0, 1d0, 1d0/
!
if (IPOINT.eq.4) then
    NINT = 4
    W    = 1d0
    G    = GP
else
    NINT = 1
    W    = 4d0
    G    = 0d0
endif
!
do m=1,NINT
    call SHAP20 (SG(m)*G,TG(m)*G,XL,DET,SH,XS,QUAD,ierr,3)
    if (ierr.ne.0) then
        return
    endif
    call BMATRIX (B,SH)
    do i=1,NEE
        P(i) = P(i)+(B(1,i)*EVP(1)+B(2,i)*EVP(2)+B(3,i)*EVP(3)+ &
        B(4,i)*EVP(4)+B(5,i)*EVP(5))*DET*W*TH
    enddo
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine DEVSTN (STN,SDEV)
!
! Computes deviatior SDEV
!
! SNGL  implicit real (a-h,o-z)
implicit none

integer, parameter :: NSTR=5
!-pass
double precision :: STN(NSTR),SDEV(NSTR)
!
SDEV(1) = 0.5d0 * (STN(1) - STN(2))
SDEV(2) = -SDEV(1)
SDEV(3) = STN(3)
SDEV(4) = STN(4)
SDEV(5) = STN(5)
!
return
end
!-------------------------------------------------------------------------------
subroutine FORMBT (STN,BETA,BETB,PRPMAT,PRPLAS)

USE MODELCTX
use constants,       only: dfmin,half,one,four
use modeldefinition, only: icvis
use algebra,         only: idig

!
! Subroutine to compute strain rates from stress field. It defines BETA and
! evaluates the product BETA*STRESS. BETB is the viscoplastic strain rate.
!
! SNGL  implicit real (a-h,o-z)
implicit none
integer, parameter :: NSTR=5
!-pass
double precision :: STN(NSTR),BETA(NSTR),BETB(NSTR),PRPMAT(6),PRPLAS(9)
!-locl
double precision :: AVECT(NSTR)
double precision :: sigma, emhu, anpwr, streff, PSTRS
!-init
sigma = 0d0
emhu  = PRPMAT(3)
anpwr = PRPMAT(4)
!
if (IDIG(ICVIS,4,2).ne.0) then
!           Calculate scalar stress PSTRS that will be compared to yield stress
    call YIELDS (STN,PRPLAS,PRPMAT,PSTRS)
!           Construct plastic flow vector AVECT
    call YIELDF (STN,PRPLAS,PRPMAT,AVECT)
!           Evaluate viscoplastic strain rate matrix BETB
    call FLOWVP (AVECT,STN,PRPLAS,PRPMAT,PSTRS,BETB)
endif
!
streff = HALF*(STN(1)-STN(2))
streff = streff*streff+STN(3)*STN(3)+STN(4)*STN(4)+STN(5)*STN(5)

if (streff.gt.dfmin) then
    streff = SQRT(streff)
    sigma  = (streff/emhu)**(anpwr-ONE)/(FOUR*emhu)
endif

BETA(1) = sigma*(STN(1)-STN(2))
BETA(2) = -BETA(1)
BETA(3) = sigma*FOUR*STN(3)
BETA(4) = sigma*FOUR*STN(4)
BETA(5) = sigma*FOUR*STN(5)

if (IDIG(ICVIS,4,2).ne.0) then
    BETA(1) = BETA(1) + BETB(1)
    BETA(2) = BETA(2) + BETB(2)
    BETA(3) = BETA(3) + BETB(3)
    BETA(4) = BETA(4) + BETB(4)
    BETA(5) = BETA(5) + BETB(5)
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine FORMJB (BETAJ,STN,PRPMAT,PRPLAS,BETB)

USE MODELCTX
USE ALGEBRA
use constants,       only: dfmin, half, one, four
use modeldefinition, only: icvis

!
! Constructs the Jacobean of the stress-strain rate matrix BETAJ.
!
implicit none
integer, parameter :: NSTR=5
!-pass
double precision   :: BETAJ(6,6),STN(NSTR),PRPLAS(9),PRPMAT(6),BETB(NSTR)
!-locl
integer            :: j
double precision   :: betajp(NSTR,NSTR),avect(NSTR)
double precision   :: streff
double precision   :: emhu, anpwr, strinv, PSTRS, sigma
double precision   :: sxx, syy, sxy, sxz, syz
double precision   :: a1,a2,a3
!
streff = HALF*(STN(1)-STN(2))
streff = streff*streff+STN(3)*STN(3)+STN(4)*STN(4)+STN(5)*STN(5)

if (streff.lt.dfmin) then
    call CLEAR (BETAJ,36,"BETAJ")
    return
endif

streff = SQRT(streff)
!
!       Plastic flow?
if (IDIG(ICVIS,4,2).ne.0) then
    call CLEAR (betajp,9,"betajp")
!           calculate scalar plastic stress PSTRS
    call YIELDS (STN,PRPLAS,PRPMAT,PSTRS)
!           Construct plastic flow vector AVECT
    call YIELDF (STN,PRPLAS,PRPMAT,avect)
!           construct Jacobian due to plastic yielding BETAJP
    call TANGVP (STN,PSTRS,PRPLAS,PRPMAT,BETB,avect,betajp)
endif
!
emhu  = PRPMAT(3)
anpwr = PRPMAT(4)
strinv = ONE/streff
sxx = (STN(1)-STN(2))*strinv
syy = -sxx
sxy = STN(3)*strinv
sxz = STN(4)*strinv
syz = STN(5)*strinv
 a1 = anpwr-ONE
 a2 = a1/FOUR
 a3 = a1*FOUR

BETAJ(1,1) =  ONE + a2*sxx*sxx
BETAJ(1,2) = -ONE + a2*sxx*syy
BETAJ(1,3) =        a1*sxx*sxy
BETAJ(1,4) =        a1*sxx*sxz
BETAJ(1,5) =        a1*sxx*syz
BETAJ(2,2) =  ONE + a2*syy*syy
BETAJ(2,3) =        a2*syy*sxy
BETAJ(2,4) =        a2*syy*sxz
BETAJ(2,5) =        a2*syy*syz
BETAJ(3,3) = FOUR + a3*sxy*sxy
BETAJ(3,4) =        a3*sxy*sxz
BETAJ(3,5) =        a3*sxy*syz
BETAJ(4,4) = FOUR + a3*sxz*sxz
BETAJ(4,5) =        a3*sxz*syz
BETAJ(5,5) = FOUR + a3*syz*syz

BETAJ(2,1) = BETAJ(1,2)
BETAJ(3,1) = BETAJ(1,3)
BETAJ(3,2) = BETAJ(2,3)
BETAJ(4,1) = BETAJ(1,4)
BETAJ(4,2) = BETAJ(2,4)
BETAJ(4,3) = BETAJ(3,4)
BETAJ(5,1) = BETAJ(1,5)
BETAJ(5,2) = BETAJ(2,5)
BETAJ(5,3) = BETAJ(3,5)
BETAJ(5,4) = BETAJ(4,5)

sigma = ((streff/emhu)**a1)/(FOUR*emhu)
!
if (IDIG(ICVIS,4,2).ne.0) then
!           Add Jacobians of viscoplastic yield and viscous flow
    do j=1,NSTR
        BETAJ(1,j) = sigma*BETAJ(1,j) + BETAJP(1,j)
        BETAJ(2,j) = sigma*BETAJ(2,j) + BETAJP(2,j)
        BETAJ(3,j) = sigma*BETAJ(3,j) + BETAJP(3,j)
        BETAJ(4,j) = sigma*BETAJ(4,j) + BETAJP(3,j)
        BETAJ(5,j) = sigma*BETAJ(5,j) + BETAJP(3,j)
    enddo
else
!           Jacobian of viscous flow
    do j=1,NSTR
        BETAJ(1,j) = sigma*BETAJ(1,j)
        BETAJ(2,j) = sigma*BETAJ(2,j)
        BETAJ(3,j) = sigma*BETAJ(3,j)
        BETAJ(4,j) = sigma*BETAJ(4,j)
        BETAJ(5,j) = sigma*BETAJ(5,j)
    enddo
endif
!
return
end
!-----------------------------------------------------------------------
subroutine ADDSTR (B,X,D,IEN,LM,MAT,PRPMAT,IELSTR,ISSIDE,STRS, &
 ISTR,LMF,TFAULT,SKEW,DXE)

USE MATERIALSMODULE, only : LMAT, &
                          MaterialIDOfElement
USE MODELDEFINITION
USE ALGEBRA
use constants, only: eps, useDXE, half, two
use debugmodule,     only: iecho
use timestepmodule,  only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT
! Subroutine to add stresses to load vector

implicit none
integer, parameter :: NSD=2,NDOF=3,NSTR=5,NEN=4
!-pass
real denominator

integer          :: mat, ien, lm, ielstr, isside, istr, lmf
double precision :: x, d, prpmat, b, strs, tfault, skew, dxe
dimension MAT(*),IEN(NEN,*),X(NSD,*),D(NDOF,*),PRPMAT(6,*),B(*), &
 LM(NDOF,NEN,*),STRS(NSTR,*),IELSTR(*),ISSIDE(*),ISTR(2,*), &
 LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
 SKEW(2,*),DXE(NDOF,NEN,*)
!-locl
double precision :: anglevec(2)
double precision :: f(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),tmp(NSD)
double precision :: rot(3,3)
!double precision, external :: FLOAT
integer          :: i,ii,j,jj,k,m,n, it0, it1, numadd
double precision :: sxx, syy, sxy, sxz, syz, thick, angle
!-init
numadd = 0

do 100 n=1,NUMSTR
    k = IELSTR(n)

    if (k.le.0) then
        write(stderr,3) n
 3           format(1x,'ADDSTR: fatal IELSTR error. Index=',I5)
        call stoper()
    endif

    it0 = ISTR(1,n)
    it1 = ISTR(2,n)

    if (DBLE(NSTEP-it0)*DBLE(NSTEP-it1).gt.0d0) then
        goto 100
    endif

    numadd = numadd + 1
    call LCOORD (X,XL,IEN(1,k),0)

    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,k),NDOF,NEN)
        call AddFaultDisplacement (DL,k,TFAULT,NDOF,NEN,25)
        call ADDSNE (DL,k,NDOF,NEN,useDXE)
        call REZONE (XL,DL,1d0)
    endif

    i = ISSIDE(n)

    if (i.le.0 .or. i.gt.NEN) then
        write(stderr,4) k
 4           format(1x,'ADDSTR: fatal ISSIDE error on element ',I5)
        call stoper()
    endif

    j = i+1
    if (IEN(3,k).ne.IEN(4,k)) then
        if (j.gt.NEN) j=1
    else
        if (j.ge.NEN) j=1
    endif

    denominator = 1.0D0 / FLOAT(it1-it0+1)
!            sxx = STRS(1,n) / FLOAT(it1-it0+1)
!            syy = STRS(2,n) / FLOAT(it1-it0+1)
!            sxy = STRS(3,n) / FLOAT(it1-it0+1)
!            sxz = STRS(4,n) / FLOAT(it1-it0+1)
!            syz = STRS(5,n) / FLOAT(it1-it0+1)
    sxx = STRS(1,n) * denominator
    syy = STRS(2,n) * denominator
    sxy = STRS(3,n) * denominator
    sxz = STRS(4,n) * denominator
    syz = STRS(5,n) * denominator

    m = MaterialIDOfElement(k)

    thick = PRPMAT(6,m)
    xn(1) = xl(2,j) - xl(2,i)
    xn(2) = xl(1,i) - xl(1,j)
    angle = (-SKEW(1,IEN(i,k))-SKEW(1,IEN(j,k)))*HALF
    anglevec(1) = angle
    anglevec(2) = 0d0

    if (ABS(angle).gt.EPS) then
        call FORMRT (anglevec,rot,2)
        call EQUATE (tmp,xn,2)
        call VMPRD (rot,tmp,xn,2,3)
    endif

    f(1) = (xn(1)*SXX + xn(2)*SXY)*THICK
    f(2) = (xn(1)*SXY + xn(2)*SYY)*THICK
    f(3) = (xn(1)*SXZ + xn(2)*SYZ)*THICK

    do m=1,3
        ii=LM(m,i,k)
        jj=LM(m,j,k)
        if(ii.ne.0) B(ii)=B(ii)+F(m)/TWO
        if(jj.ne.0) B(jj)=B(jj)+F(m)/TWO
    enddo
100    continue
!
if (iecho.eq.1 .and. numadd.gt.0) then
    write(stdout,5)
 5   format(1x,'Stress loads are being added')
endif
!
return
end
!-------------------------------------------------------------------------------
 subroutine ADDWNK (IWELM,IWSIDE,IWTIME,WPRES,IEN,MAT,X,D,PRPMAT, &
   A,LM,IDIAG,LMF,TFAULT,SKEW,DXE,LMX,NSLIP)

USE MATERIALSMODULE, only : LMAT, &
                            MaterialIDOfElement
USE MODELDEFINITION
USE ALGEBRA
use constants,       only: eps, useDXE, half, one
use modeltopology,   only: nsd, ndof, nen
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: elementSlipperyID
#endif
use meshdatamodule,  only: meshdatactx
use petscmodule,     only: M_A
use petscksp,        only: add_values
use timestepmodule,  only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT
!
! Subroutine to add Winkler pressures to diagonal elements of stiffness
! matrix
!
! SNGL    implicit real (a-h,o-z)
implicit none

!-pass
integer          :: elemID

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer          :: LMX(NDOF,NEN,nElemsWithSlipperyNodes)
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
integer          :: LMX(NDOF,NEN,meshdatactx%neglobal)
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif

double precision :: anglevec(2)
integer          :: iwelm, iwside, iwtime, ien, mat, idiag, lm, lmf, nslip
double precision :: wpres, x, prpmat, a, d, tfault, skew

 dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),A(*),IDIAG(*),LM(NDOF,NEN,*), &
   D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), NSLIP(5,*)

#ifdef SPARSE
! for petscMPIInt
#include "petsc/finclude/petscsys.h"
   PetscMPIInt irank
#endif
!-locl
double precision :: w(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),rot(3,3)
double precision :: tmp(NSD)
integer          :: lien(2)

integer          :: i, ierr, ii, j, LL, m, n, ni, nj, nn, v, wt
double precision :: angle, area

if (NUMWNK.eq.0) then
    return
endif

 do 400 nn=1,IABS(NUMWNK)
n = IABS(IWELM(nn))
LL= IWTIME(nn)
    if (LL.lt.0) then
        if (LL.ne.-1 .and. NSTEP.ge.-LL) goto 400
    else if (LL.eq.0) then
        goto 400
    else
    if (NSTEP.lt.(LL-1)) goto 400
    endif
call LCOORD (X,XL,IEN(1,n),0)
if (LGDEF.ne.0) then
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,26)
    call ADDSNE (DL,n,NDOF,NEN,useDXE)
call REZONE (XL,DL,ONE)
endif
call SIDENP(IABS(IWSIDE(nn)),IEN(3,n).ne.IEN(4,n),lien)
    call FCGEOM (XL,lien,area,xn)
i = lien(1)
j = lien(2)
    angle = (-SKEW(1,IEN(i,n))-SKEW(1,IEN(j,n)))*HALF
anglevec(1)=angle
anglevec(2)=0d0


    if (ABS(angle).gt.EPS) then
        call FORMRT (anglevec,rot,NSD)
        call EQUATE (tmp,xn,NSD)
        call VMPRD (rot,tmp,xn,NSD,3)
    endif
m = MaterialIDOfElement(n)
wt = PRPMAT(6,m) * area * HALF
    w(1) = ABS(xn(1)) * WPRES(nn) * wt
    w(2) = ABS(xn(2)) * WPRES(nn) * wt
ii = LM(1,i,n)
#ifdef SPARSE
if (ii.gt.0) then
call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,w(1),ADD_VALUES, &
       ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
 1               format(1x,'ADDWNK: MatSetValuesLocal error ',I4)
            call stoper()
        endif
endif
#else
if (ii.gt.0) ii = IDIAG(ii)
if (ii.gt.0) A(ii) = A(ii) + w(1)
#endif
ii = LM(1,j,n)
#ifdef SPARSE
if (ii.gt.0) then
call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,w(1),ADD_VALUES, &
       ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
            call stoper()
        endif
endif
#else
if (ii.gt.0) ii = IDIAG(ii)
if (ii.gt.0) A(ii) = A(ii) + w(1)
#endif
ii = LM(2,i,n)
#ifdef SPARSE
if (ii.gt.0) then
call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,w(2),ADD_VALUES, &
       ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
            call stoper()
        endif
endif
#else
if (ii.gt.0) ii = IDIAG(ii)
if (ii.gt.0) A(ii) = A(ii) + w(2)
#endif
ii = LM(2,j,n)
#ifdef SPARSE
if (ii.gt.0) then
call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,w(2),ADD_VALUES, &
       ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
            call stoper()
        endif
endif
#else
if (ii.gt.0) ii = IDIAG(ii)
if (ii.gt.0) A(ii) = A(ii) + w(2)
#endif
   if (NUMSLP.eq.0) goto 400
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
if(elementSlipperyID(n,61).gt.0) then
    ii = IABS(LMX(1,i,elementSlipperyID(n,61)))
else
    ii = 0
endif
#else
ii = IABS(LMX(1,i,n))
#endif



if (ii.gt.0) then
#ifdef SPARSE
v = w(1)/DBLE(ni)
call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,v,ADD_VALUES, &
       ierr)
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
if(elementSlipperyID(n,62).gt.0) then
    ii = IABS(LMX(1,j,elementSlipperyID(n,62)))
else
    ii = 0
endif
#else
ii = IABS(LMX(1,j,n))
#endif


if (ii.gt.0) then
#ifdef SPARSE
   if (irank.eq.0) then
v = w(1)/DBLE(nj)
call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,v,ADD_VALUES, &
       ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
            call stoper()
        endif
   endif
#else
ii = IDIAG(ii)
    A(ii) = A(ii) + w(1)/DBLE(nj)
#endif
endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
if(elementSlipperyID(n,63).gt.0) then
    ii = IABS(LMX(2,i,elementSlipperyID(n,63)))
else
    ii = 0
endif
#else
ii = IABS(LMX(2,i,n))
#endif


if (ii.gt.0) then
#ifdef SPARSE
v = w(2)/DBLE(ni)
call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,v,ADD_VALUES, &
       ierr)
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
if(elementSlipperyID(n,64).gt.0) then
    ii = IABS(LMX(2,j,elementSlipperyID(n,64)))
else
    ii = 0
endif
#else
ii = IABS(LMX(2,j,n))
#endif



if (ii.gt.0) then
#ifdef SPARSE
v = w(2)/DBLE(nj)
call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,v,ADD_VALUES, &
       ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
            call stoper()
        endif
#else
ii = IDIAG(ii)
A(ii) = A(ii) + w(2)/DBLE(nj)
#endif
endif


400    continue
!
 return
 end
!-------------------------------------------------------------------------------
subroutine ELCRD (S,T,XELM,IEN,X,D,elemID,TFAULT, &
 SKEW,FACTOR,ierr)

USE MODELDEFINITION
USE ALGEBRA
use constants, only: useDXE
!
! Calculates (deformed) coordinates XELM at natural coordinates (S,T)
!
! SNGL  implicit real (a-h,o-z)
implicit none
integer, parameter :: NSD=2, NDOF=3, NEN=4
!-pass
integer ierr
integer :: elemID
integer :: ien
double precision :: xelm, x, d, tfault, skew,s,t,factor
dimension XELM(NSD),IEN(NEN),X(NSD,*),D(NDOF,*), &
   TFAULT(NDOF,*),SKEW(2,*)
!-locl
double precision :: XL(NSD,NEN),DL(NDOF,NEN),SH(3,NEN),XS(NSD,NSD)
double precision :: det
logical QUAD
integer :: j
!-init
call CLEAR (XELM,NSD,"XELM")
!
QUAD = (IEN(3).ne.IEN(4))
call LCOORD (X,XL,IEN,0)

if (FACTOR.gt.0.) then
    call LDISP  (DL,D,IEN,NDOF,NEN)
    call AddFaultDisplacement (DL,elemID,TFAULT,NDOF,NEN,27)
    call ADDSNE (DL,elemID,NDOF,NEN,useDXE)
    call REZONE (XL,DL,FACTOR)
endif
!
call SHAP20 (S,T,XL,DET,SH,XS,QUAD,ierr,1)

if (ierr.ne.0) then
    return
endif
!
do j=1,NEN
    XELM(1) = XELM(1) + SH(3,j)*XL(1,j)
    XELM(2) = XELM(2) + SH(3,j)*XL(2,j)
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine ELCTR (XELM,IEN,X,D,elemID,TFAULT, &
 SKEW,FACTOR,ierr)

use constants, only: third

!
! Calculates (deformed) element center coordinates XELM
!
implicit none
integer, parameter :: NSD=2, NDOF=3, NEN=4
!-pass
integer :: ierr
integer :: elemID

integer          :: ien
double precision :: xelm, x, d, tfault, skew, factor
dimension XELM(NSD),IEN(NEN),X(NSD,*),D(NDOF,*), &
 TFAULT(NDOF,*), SKEW(2,*)

double precision :: s, t
!-locl
save s
!-init
data s/0d0/
t = 0d0
!
if (IEN(3).eq.IEN(4)) then
    t = -THIRD
endif

call ELCRD (s,t,XELM,IEN,X,D,elemID,TFAULT, SKEW,FACTOR,ierr)
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
integer ierr
logical QUAD
double precision :: XL(2,NEN),S(NEN,NEN),C(NSD)
!-locl
integer          :: i, j, m
double precision :: SH(3,NEN),SG(NEN),TG(NEN),XS(NSD,NSD)
double precision :: volume, detJ,TH
save TG,SG
data TG/-1d0,-1d0, 1d0, 1d0/
data SG/-1d0, 1d0, 1d0,-1d0/
!
do m=1,NINT
!        compute spatial derivatives of shape functions
    call SHAP20 (SG(m)*G,TG(m)*G,XL,detJ,SH,XS,QUAD,ierr,3)

    if (ierr.ne.0) then
        return
    endif

    volume = detJ * TH

    do j=1,NEN
        do i=1,j
            S(i,j)=S(i,j)+volume*(C(1)*SH(1,i)*SH(1,j) &
             +C(2)*SH(2,i)*SH(2,j))
        enddo
    enddo
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine ADVSTF (XL,VL,QUAD,S,RhoCp,TH,ierr)

USE CONSTANTS, only: G,dfmin,dlogmin,dfmax,dlogmax

!
! convection (advection) update of local conduction stiffness matrix for
! steady convection-conduction solution
!
! SNGL  implicit real (a-h,o-z)
implicit none
integer, parameter :: NEN=4, NINT=4, NSD=2
!-pass
logical QUAD
integer ierr
double precision RhoCp,TH
double precision :: XL(2,NEN),VL(2,NEN),S(NEN,NEN)
!-locl
integer          :: i,j,l
double precision :: sh(3,NEN),sg(NEN),tg(NEN),xs(NSD,NSD),q(NEN)
double precision :: Vx, Vy,v,w,detJ
save tg,sg
data tg/-1d0,-1d0, 1d0, 1d0/
data sg/-1d0, 1d0, 1d0,-1d0/
!
do l=1,NINT
    call SHAP20 (sg(l)*G,tg(l)*G,XL,detJ,sh,xs,QUAD,ierr,3)

    if (ierr.ne.0) then
        return
    endif

    Vx = sh(3,1)*VL(1,1) + &
         sh(3,2)*VL(1,2) + &
         sh(3,3)*VL(1,3) + &
         sh(3,4)*VL(1,4)

    Vy = sh(3,1)*VL(2,1) + &
         sh(3,2)*VL(2,2) + &
         sh(3,3)*VL(2,3) + &
         sh(3,4)*VL(2,4)

    do i=1,NEN
        v = SQRT(VL(1,i)*VL(1,i) + &
                 VL(2,i)*VL(2,i))
        if (v.gt.dfmin) then
            q(i) = 50*(VL(1,i)*sh(1,i) + &
                       VL(2,i)*sh(2,i))/(2d0*v)
        else
            q(i) = 0d0
        endif
    enddo

    w = detJ * TH * RhoCp

    do i=1,NEN
        do j=1,i
            S(i,j) = S(i,j) &
             + w*(sh(3,i)+q(i))*(Vx*sh(1,j)+Vy*sh(2,j))
        enddo
    enddo

enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine LHEAT (XL,QUAD,PL,HEAT,TH,ierr)

USE CONSTANTS, only: G
!
! computes the local load vector contribution of heat production.
!
implicit none

integer, parameter :: NEN=4, NINT=4, NSD=2
!-pass
integer ierr
logical QUAD
double precision :: XL(NSD,NEN),PL(NEN),TH,HEAT
!-locl
double precision :: SH(3,NEN),SG(NEN),TG(NEN),XS(NSD,NSD)
integer          :: m
double precision :: volume, detJ
save SG,TG
data SG/-1.0d0, 1.0d0,1.0d0,-1.0d0/
data TG/-1.0d0,-1.0d0,1.0d0, 1.0d0/
!
do m=1,NINT
!           Compute shape functions
    call SHAP20 (SG(m)*G,TG(m)*G,XL,detJ,SH,XS,QUAD, ierr,1)

    if (ierr.ne.0) then
        return
    endif

    volume = detJ * TH

    PL(1) = PL(1) + HEAT*volume*SH(3,1)
    PL(2) = PL(2) + HEAT*volume*SH(3,2)
    PL(3) = PL(3) + HEAT*volume*SH(3,3)
    PL(4) = PL(4) + HEAT*volume*SH(3,4)
enddo
!
return
end
!-------------------------------------------------------------------------------
double precision function EigMax (STN,E,v,emhu,Npl,alpha)

use constants, only: dfmin, one, two, three, four, eight

!
! Gives an upper limit estimate of the maximum local eigenvalue
!
! SNGL  implicit real (a-h,o-z)
implicit none
!-pass
! SNGL  real E,v,emhu,Npl,alpha
double precision E,v,emhu,Npl,alpha
double precision :: STN(5)

!-locl
double precision, external :: AJ2
double precision           :: NplMin1
double precision           :: A,B,C,D,tmp,SE
!
NplMin1 = Npl - one
SE      = AJ2 (STN,5d-1)

if (SE.gt.dfmin) then
    SE = SQRT(SE)
endif

tmp     = eight * (one+v)*(one-two*v) * emhu**Npl
tmp     = E * SE**NplMin1 / tmp

if (alpha.le.dfmin) then
    !  explicit stepping
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
    ! implicit differencing
    D = two * (three*Npl + SQRT(Npl**two + two*NplMin1**two))
endif
EigMax = tmp * D
!
return
end
!-----------------------------------------------------------------------
subroutine ADDFLX (BT,IEN,X,D,MAT,PRPMAT,IFLS,IFLX,BFLX,LMT, &
 LMTX,LMF,TFAULT,SKEW,DXE)

   USE MATERIALSMODULE, only : LMAT
   USE MODELDEFINITION
use constants, only: useDXE
!
! Routine to add normal boundary fluxes to load vector
!
implicit none
integer, parameter :: NEN=4, NET=3, NSD=2, NDOF=3
!-pass
!        integer NFLX
integer            :: ien, mat, ifls, iflx, lmt, lmf, lmtx
double precision   :: bt, x, prpmat, bflx, d, tfault, skew, dxe

dimension BT(*),IEN(NEN,*),X(NSD,*),MAT(*),PRPMAT(6,*),IFLS(*), &
 IFLX(*),BFLX(*),LMT(NEN,*),D(NDOF,*),LMF(NDOF,NEN,*), &
 TFAULT(NDOF,*),LMTX(NEN,*),SKEW(2,*),DXE(NDOF,NEN,*)
!-locl
double precision   :: XL(NSD,NEN),DL(NDOF,NEN)
double precision   :: f, deltx, delty, thick, size
integer            :: i,j,m,ii,jj,n
!
if (NFLX.le.0) then
    return
endif
!
do n=1,NFLX

    i = IFLS(n)
    j = i + 1

    if (IEN(3,n).ne.IEN(4,n)) then
        if (j.gt.NEN) then
            j = 1
        endif
    else
        if (j.gt.NET) then
            j = 1
        endif
    endif

    m = IFLX(n)
    call LCOORD (X,XL,IEN(1,m),0)

    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,m),NDOF,NEN)
        call AddFaultDisplacement (DL,m,TFAULT,NDOF,NEN,28)
        call ADDSNE (DL,m,NDOF,NEN,useDXE)
        call REZONE (XL,DL,1d0)
    endif

    deltx = XL(1,i) - XL(1,j)
    delty = XL(2,i) - XL(2,j)
    size = SQRT(deltx*deltx + delty*delty)

    if (LMAT) then
        THICK = PRPMAT(6,m)
    else
        THICK = PRPMAT(6,MAT(m))
    endif

    f = -size*THICK*BFLX(n)/2d0
    ii = LMT(i,m)

    if (ii.eq.0.and.NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
       ii = LMTX(i,m)
    endif

    jj = LMT(j,m)

    if (jj.eq.0.and.NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
       jj = LMTX(j,m)
    endif

    if (ii.gt.0) then
        BT(ii) = BT(ii) + f
    endif

    if (jj.gt.0) then
        BT(jj) = BT(jj) + f
    endif
enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine TEMWNK (AA,LMT,ITDIAG,ITWINK,TWINK,IEN,MAT,PRPMAT, &
   X,D,LMF,TFAULT,SKEW,DXE)

USE MATERIALSMODULE, only: LMAT
USE MODELDEFINITION
use constants,       only: useDXE, half, one
use debugmodule,     only: iecho
use timestepmodule,  only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

!
! Subroutine to add Winkler fluxes to diagonal elements of thermal stiffness
! matrix
!
! SNGL    implicit real (a-h,o-z)
implicit none

integer, parameter :: NSD=2,NDOF=2,NEN=4
!-pass
integer            :: itwink, ien, mat, lmt, lmf
double precision   :: twink, x, prpmat, aa, itdiag, d, tfault, skew, dxe
 dimension ITWINK(3,*),TWINK(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),AA(*),ITDIAG(*),LMT(NEN,*), &
   D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN,*)

integer :: n,nn,LL
!-locl
double precision   :: xl(NSD,NEN),dl(NDOF,NEN),xn(NSD)
double precision   :: wt, area
integer            :: lien(2),ii
!
if (NTWINK.eq.0) then
    return
endif
!
if (iecho.ne.0) then
    write(stdout,*) 'Adding Winkler fluxes'
endif
!
do 400 nn=1,NTWINK
    n = ITWINK(1,nn)
    LL= ITWINK(3,nn)

    if (LL.lt.0) then
        if (LL.ne.-1 .and. NSTEP.ge.-LL) then
            goto 400
        endif
    else if (LL.eq.0) then
        goto 400
    else
        if (NSTEP.lt.(LL-1)) then
            goto 400
        endif
    endif

    call LCOORD (X,XL,IEN(1,n),0)

    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,29)
        call ADDSNE (DL,n,NDOF,NEN,useDXE)

!rezone
!double precision :: XL(NSD,NEN),DL(NDOF,NEN)
!double precision :: factor

        call REZONE (XL,DL,ONE)
    endif

    call SIDENP (ITWINK(2,nn),IEN(3,n).ne.IEN(4,n),lien)
    call FCGEOM (XL,lien,area,xn)    ! compute face area

    if (LMAT) then
        wt = PRPMAT(6,n) * area * TWINK(nn) * HALF
    else
        wt = PRPMAT(6,MAT(n)) * area * TWINK(nn) * HALF
    endif

    ii = LMT(lien(1),n)            ! equation number
    if (ii.gt.0) ii = ITDIAG(ii)    ! matrix address of eqno
    if (ii.gt.0) AA(ii) = AA(ii) + wt
    ii = LMT(lien(2),n)
    if (ii.gt.0) ii = ITDIAG(ii)
    if (ii.gt.0) AA(ii) = AA(ii) + wt
400    continue
!
 return
 end
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
double precision :: XL(NSD,NEN),TL(NEN),GRAD(NSD)
double precision :: s,t
!-locl
double precision :: SA(NEN),TA(NEN),SH(NPE,NEN),XS(NSD,NSD)
double precision :: dTds, dTdt, det
integer          :: i

save SA,TA
!-init
data SA/-1d0, 1d0, 1d0, -1d0/
data TA/-1d0,-1d0, 1d0,  1d0/

!       calculate inverse Jacobian matrix XS
call SHAP20 (s,t,XL,DET,SH,XS,QUAD,ierr,2)

if (ierr.ne.0) then
    return
endif
!       Calculate s- and t-derivatives of shape functions
do i=1,NEN
    SH(1,i) = quart * SA(i) * (one + TA(i)*t)
    SH(2,i) = quart * TA(i) * (one + SA(i)*s)
enddo

if (.not.QUAD) then
!           adjust shape function of third node for triangular elements
    do i=1,2
        SH(i,3) = SH(i,3)+SH(i,4)
        SH(i,4) = zero
    enddo
endif
!       calculate s- and t- thermal gradients
dTds = zero
dTdt = zero

do i=1,NEN
    dTds = dTds + SH(1,i)*TL(i)
    dTdt = dTdt + SH(2,i)*TL(i)
enddo

GRAD(1) = dTds*XS(1,1) + dTdt*XS(1,2)
GRAD(2) = dTds*XS(2,1) + dTdt*XS(2,2)
!
return
end
!-----------------------------------------------------------------------
double precision function Tcentr (T,TX,XL,IEN,LMTX,ierr)

USE MODELDEFINITION, only: NUMFN, NUMSLP, IDIFT
use constants,       only: zero, one, third
!
! Calculates the temperature in the center of the element
!
implicit none

integer, parameter :: NEN=4, NSD=2
double precision, parameter :: THIRDM=-third
!-pass
integer          :: IEN, LMTX, ierr
double precision :: T, TX, XL
dimension T(*),TX(*),IEN(NEN),LMTX(NEN),XL(NSD,NEN)
!-locl
logical QUAD
double precision :: TL(NEN),SH(3,NEN),XS(NSD,NSD)
double precision :: tn, det
!
QUAD=(IEN(3).ne.IEN(4))
!
!       Localize temperatures
call LTEMP (IEN,T,TL,NEN)
!
!       Modify for discontinuities
if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
   call LTEMPX (LMTX,IEN,TX,TL,NEN)
endif
!
!       Compute shape functions
tn = ZERO

if (.not.QUAD) then
    tn=THIRDM
endif

call SHAP20 (ZERO,tn,XL,DET,SH,XS,QUAD,ierr,1)
!
!       Calculate center temperature
Tcentr = SH(3,1)*TL(1)+SH(3,2)*TL(2)+SH(3,3)*TL(3)+SH(3,4)*TL(4)
!
return
end
!-------------------------------------------------------------------------------
subroutine CPSTIFF (XL,QUAD,S,RhoCp,TH,ierr)

USE CONSTANTS, only: G
!
! computes the local heat capacity stiffness matrix from integration at Barlow
! points.
!
implicit none

integer, parameter :: NEN=4, NINT=4, NSD=2
!-pass
integer ierr
logical QUAD
double precision :: XL(2,NEN),S(NEN,NEN),TH,RhoCp
!-locl
double precision :: sh(3,NEN),sg(NEN),tg(NEN),xs(NSD,NSD)
double precision :: detJ, volume
integer :: i, j, m

save tg,sg
data tg/-1d0,-1d0, 1d0, 1d0/
data sg/-1d0, 1d0, 1d0,-1d0/
!
do m=1,NINT
    call shap20 (sg(m)*G,tg(m)*G,XL,detJ,sh,xs,QUAD,ierr,1)
    if (ierr.ne.0) then
        return
    endif
    volume = detJ * TH
    do j=1,NEN
        do i=1,j
            S(i,j) = S(i,j) + RhoCp * sh(3,i) * sh(3,j) * volume
        enddo
    enddo
enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
double precision function AJ1 (ST,POIS)
!
! Calculates first invariant of stress tensor
!
! SNGL  implicit real (a-h,o-z)
implicit none
!-pass
double precision :: ST(5), pois
!-init
AJ1 = (ST(1)+ST(2))*(1d0+POIS)
!
return
end function
!-------------------------------------------------------------------------------
double precision function AJ2 (ST,POIS)

! Calculates second invariant of deviatoric stress tensor

implicit none
!-pass
double precision :: ST(5), pois
!-local
double precision :: sp, sum
double precision, external :: AJ1

!-init
sp = AJ1(ST,POIS)/3d0
sum = (ST(1)-sp)*(ST(1)-sp)+(ST(2)-sp)*(ST(2)-sp)
AJ2 = sum*0.5d0 + ST(3)*ST(3) + ST(4)*ST(4) + ST(5)*ST(5)
!
return
end function
!-------------------------------------------------------------------------------
double precision function AJ3 (ST,POIS)
!
! Calculates third invariant of deviatoric stress tensor
!
! SNGL  implicit real (a-h,o-z)
implicit none

!-pass
double precision :: ST(5), pois

!-local
double precision :: p, sxx,syy,szz,sxy,sxz,syz

double precision, external :: pressure

!-init
p = pressure(ST,POIS)
sxx = ST(1)-p
syy = ST(2)-p
szz = POIS*(ST(1)+ST(2)) - p
sxy = ST(3)
sxz = ST(4)
syz = ST(5)
AJ3 = sxx*(syy*szz-syz*syz) - &
      sxy*(sxy*szz-syz*syz) + &
      syz*(sxy*syz-sxz*syy)
return
end function
!-------------------------------------------------------------------------------
subroutine REZONE (XL,DL,FACTOR)
!
! Routine to update local coordinate XL by displacements DL
!
implicit none
integer, parameter :: NSD=2, NDOF=3, NEN=4
!-pass
double precision :: XL(NSD,NEN),DL(NDOF,NEN)
double precision :: factor
!
XL(1,1) = XL(1,1) + DL(1,1)*FACTOR
XL(2,1) = XL(2,1) + DL(2,1)*FACTOR
XL(1,2) = XL(1,2) + DL(1,2)*FACTOR
XL(2,2) = XL(2,2) + DL(2,2)*FACTOR
XL(1,3) = XL(1,3) + DL(1,3)*FACTOR
XL(2,3) = XL(2,3) + DL(2,3)*FACTOR
XL(1,4) = XL(1,4) + DL(1,4)*FACTOR
XL(2,4) = XL(2,4) + DL(2,4)*FACTOR
!
return
end
!-------------------------------------------------------------------------------
subroutine TSTRN (EE,DT,VEXP)

use constants, only: zero, one, third
!
! calculates strain due to thermal expansion
!
implicit none
integer, parameter :: NSTR=5
!-pass
double precision :: EE(NSTR)
double precision :: DT, VEXP
!
EE(1) = DT*VEXP*THIRD
EE(2) = EE(1)
EE(3) = ZERO
EE(4) = ZERO
EE(5) = ZERO
!
return
end
!-----------------------------------------------------------------------
subroutine YIELDS (STN,PRPLAS,PRPMAT,PSTRS)

! Routine to calculate the scalar stress PSTRS that will be compared to
! the (equivalent) uniaxial yield stress. See GTECTON notes, section
! on yield criteria.

USE CONSTANTS, only: dfmin, zero, half, one, two, three, root3, deg2rad

implicit none
integer, parameter :: NSTR=5
!-pass
double precision :: PSTRS
double precision :: STN(NSTR),PRPLAS(9),PRPMAT(6)
!-locl
double precision :: i1
double precision, external :: EFFSTN,AJ1,AJ3

double precision :: stneff, sint3, frangl, sinphi, alpha, theta

integer          :: nptype

!
nptype = INT(PRPLAS(1))

if (nptype.eq.0) then
    return
endif
!
!       compute effective stress and pressure
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
    sinphi = SIN(frangl*DEG2RAD)
    PSTRS  = ABS(i1) * sinphi / THREE &
     + stneff * (COS(theta) - ((SIN(theta)*sinphi)/root3))
else if (nptype.eq.4) then
!           Drucker-Prager
    frangl = PRPLAS(2)
    sinphi = SIN(frangl*DEG2RAD)
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
subroutine YIELDF (STN,PRPLAS,PRPMAT,AVECT)
USE MODELDEFINITION
USE CONSTANTS, only: deg2rad, root3, zero, half, one, two, three, dfmin,dlogmin,dfmax,dlogmax

!
!       Compute visco-plastic flow vector AVECT (see GTECTON notes).
!
implicit none
!
integer, parameter :: NSTR=5
!-pass
double precision STN,PRPLAS,AVECT,PRPMAT
dimension STN(NSTR),PRPLAS(9),AVECT(NSTR),PRPMAT(6)
!-locl
integer i,nptype
double precision stneff,i2,sp,sxx,syy,szz,veca1(NSTR), &
 veca2(NSTR),veca3(NSTR),const1,const2,const3,sint3,theta,cost, &
 cos3t,sint,tant,tan3t,frangl,sinp,alpha

double precision, external :: EFFSTN,AJ1,AJ3

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
veca1(1) = 1.5D0
veca1(2) = 1.5D0
veca1(3) = ZERO
veca1(4) = ZERO
veca1(5) = ZERO
veca2(1) = sxx/(TWO*stneff)
veca2(2) = syy/(TWO*stneff)
veca2(3) = STN(3)/stneff
veca2(4) = STN(4)/stneff
veca2(5) = STN(5)/stneff
veca3(1) = HALF*(STN(4)*STN(4)-STN(5)*STN(5))
veca3(2) = HALF*(STN(5)*STN(5)-STN(4)*STN(4))
veca3(3) = ZERO
veca3(4) = STN(4)*(STN(1)-STN(2))
veca3(5) = STN(5)*(STN(1)-STN(2))
!
if (nptype.eq.1) then
!           Tresca
    const1 = ZERO
    sint3 = -THREE*root3*AJ3(STN,HALF)/(TWO*(stneff**THREE))
    sint3 = MIN(ONE,sint3)
    sint3 = MAX(-ONE,sint3)
    theta = ASIN(sint3)/THREE
    if (ABS(theta/DEG2RAD-30.).le.1.) then
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
    sinp   = SIN(frangl*DEG2RAD)
    const1 = sinp/THREE
    sint3 = -THREE*root3*AJ3(STN,HALF)/(TWO*(stneff**THREE))
    sint3 = MIN(ONE,sint3)
    sint3 = MAX(-ONE,sint3)
    theta = ASIN(sint3)/THREE
    if (ABS(theta/DEG2RAD-30.).le.1.) then
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
    sinp   = SIN(frangl*DEG2RAD)
!           Outer apices:
    alpha = (TWO * sinp)/(root3 * (THREE - sinp))
!           Inner apices:
!           alpha = (TWO * sinp)/(root3 * (THREE + sinp))
    const1 = alpha
    const2 = ONE
    const3 = ZERO
else
    write(*,*) "Opnlib subroutine YIELDF encountered and error"
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
end subroutine
!-------------------------------------------------------------------------------
subroutine TANGVP (STN,PSTRS,PRPLAS,PRPMAT,BETB,AVECT,BETAJP)
USE MODELDEFINITION
USE ALGEBRA
USE CONSTANTS, only: ZERO,HALF,ONE,TWO,THREE,FOUR,root3, &
                     THIRD,SIX,deg2rad, dfmin,dlogmin,dfmax,dlogmax


!
! Compute BETAJP, which is the Jacobian matrix for viscoplastic flow
! See GTECTON notes for details.
!
implicit none
integer, parameter :: NSTR=5
!-pass
double precision STN(NSTR),BETAJP(NSTR,NSTR),PRPLAS(9), &
 PRPMAT(6),BETB(NSTR),AVECT(NSTR),PSTRS
!-locl
integer nptype,nflow,j
double precision stneff,AJ1,AJ2,AJ3,fnorm,uniyld,i1,sxx,syy,sxy, &
 yldpwr,cmult,gradp,sint3,theta,p1,p2,p3,frangl,fluidy,sinp,sxz, &
 syz
double precision m1(NSTR,NSTR),m2(NSTR,NSTR),m3(NSTR,NSTR)
external AJ1,AJ2,AJ3
!-init

call CLEAR(BETAJP,NSTR*NSTR)
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
sxz = STN(4)
syz = STN(5)
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
    m1(1,1) = HALF
    m1(1,2) = -HALF
    m1(1,3) = ZERO
    m1(1,4) = ZERO
    m1(1,5) = ZERO
    m1(2,2) = HALF
    m1(2,3) = ZERO
    m1(2,4) = ZERO
    m1(2,5) = ZERO
    m1(3,3) = TWO
    m1(3,4) = ZERO
    m1(3,5) = ZERO
    m1(4,4) = TWO
    m1(4,5) = ZERO
    m1(5,5) = TWO
    m1(2,1) = m1(1,2)
    m1(3,1) = m1(1,3)
    m1(4,1) = m1(1,4)
    m1(5,1) = m1(1,5)
    m1(3,2) = m1(2,3)
    m1(4,2) = m1(2,4)
    m1(5,2) = m1(2,5)
    m1(4,3) = m1(3,4)
    m1(5,3) = m1(3,5)
    m1(5,4) = m1(4,5)
    m2(1,1) = sxx*sxx
    m2(1,2) = sxx*syy
    m2(1,3) = TWO*sxx*sxy
    m2(1,4) = TWO*sxx*sxz
    m2(1,5) = TWO*sxx*syz
    m2(2,2) = syy*syy
    m2(2,3) = TWO*syy*sxy
    m2(2,4) = TWO*syy*sxz
    m2(2,5) = TWO*syy*syz
    m2(3,3) = FOUR*sxy*sxy
    m2(3,4) = FOUR*sxy*sxz
    m2(3,5) = FOUR*sxy*syz
    m2(4,4) = FOUR*sxz*sxz
    m2(4,5) = FOUR*sxz*syz
    m2(5,5) = FOUR*syz*syz
    m2(2,1) = m2(1,2)
    m2(3,1) = m2(1,3)
    m2(4,1) = m2(1,4)
    m2(5,1) = m2(1,5)
    m2(3,2) = m2(2,3)
    m2(4,2) = m2(2,4)
    m2(5,2) = m2(2,5)
    m2(4,3) = m2(3,4)
    m2(5,3) = m2(3,5)
    m2(5,4) = m2(4,5)
    m3(1,1) = AVECT(1)*AVECT(1)
    m3(1,2) = AVECT(1)*AVECT(2)
    m3(1,3) = AVECT(1)*AVECT(3)
    m3(1,4) = AVECT(1)*AVECT(4)
    m3(1,5) = AVECT(1)*AVECT(5)
    m3(2,2) = AVECT(2)*AVECT(2)
    m3(2,3) = AVECT(2)*AVECT(3)
    m3(2,4) = AVECT(2)*AVECT(4)
    m3(2,5) = AVECT(2)*AVECT(5)
    m3(3,3) = AVECT(3)*AVECT(3)
    m3(3,4) = AVECT(3)*AVECT(4)
    m3(3,5) = AVECT(3)*AVECT(5)
    m3(4,4) = AVECT(4)*AVECT(4)
    m3(4,5) = AVECT(4)*AVECT(5)
    m3(5,5) = AVECT(5)*AVECT(5)
    m3(2,1) = m3(1,2)
    m3(3,1) = m3(1,3)
    m3(4,1) = m3(1,4)
    m3(5,1) = m3(1,5)
    m3(3,2) = m3(2,3)
    m3(4,2) = m3(2,4)
    m3(5,2) = m3(2,5)
    m3(4,3) = m3(3,4)
    m3(5,3) = m3(3,5)
    m3(5,4) = m3(4,5)
    do j=1,NSTR
        BETAJP(1,j) = P1*m1(1,j) + P2*m2(1,j) + P3*m3(1,j)
        BETAJP(2,j) = P1*m1(2,j) + P2*m2(2,j) + P3*m3(2,j)
        BETAJP(3,j) = P1*m1(3,j) + P2*m2(3,j) + P3*m3(3,j)
    enddo
elseif (nptype.eq.3) then
!           Mohr-Coulomb
    frangl = PRPLAS(2) * DEG2RAD
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
    m1(1,1) = HALF
    m1(1,2) = -HALF
    m1(1,3) = ZERO
    m1(1,4) = ZERO
    m1(1,5) = ZERO
    m1(2,2) = HALF
    m1(2,3) = ZERO
    m1(2,4) = ZERO
    m1(2,5) = ZERO
    m1(3,3) = TWO
    m1(3,4) = ZERO
    m1(3,5) = ZERO
    m1(4,4) = TWO
    m1(4,5) = ZERO
    m1(5,5) = TWO
    m1(2,1) = m1(1,2)
    m1(3,1) = m1(1,3)
    m1(4,1) = m1(1,4)
    m1(5,1) = m1(1,5)
    m1(3,2) = m1(2,3)
    m1(4,2) = m1(2,4)
    m1(5,2) = m1(2,5)
    m1(4,3) = m1(3,4)
    m1(5,3) = m1(3,5)
    m1(5,4) = m1(4,5)
    m2(1,1) = sxx*sxx
    m2(1,2) = sxx*syy
    m2(1,3) = TWO*sxx*sxy
    m2(1,4) = TWO*sxx*sxz
    m2(1,5) = TWO*sxx*syz
    m2(2,2) = syy*syy
    m2(2,3) = TWO*syy*sxy
    m2(2,4) = TWO*syy*sxz
    m2(2,5) = TWO*syy*syz
    m2(3,3) = FOUR*sxy*sxy
    m2(3,4) = FOUR*sxy*sxz
    m2(3,5) = FOUR*sxy*syz
    m2(4,4) = FOUR*sxz*sxz
    m2(4,5) = FOUR*sxz*syz
    m2(5,5) = FOUR*syz*syz
    m2(2,1) = m2(1,2)
    m2(3,1) = m2(1,3)
    m2(4,1) = m2(1,4)
    m2(5,1) = m2(1,5)
    m2(3,2) = m2(2,3)
    m2(4,2) = m2(2,4)
    m2(5,2) = m2(2,5)
    m2(4,3) = m2(3,4)
    m2(5,3) = m2(3,5)
    m2(5,4) = m2(4,5)
    m3(1,1) = AVECT(1)*AVECT(1)
    m3(1,2) = AVECT(1)*AVECT(2)
    m3(1,3) = AVECT(1)*AVECT(3)
    m3(1,4) = AVECT(1)*AVECT(4)
    m3(1,5) = AVECT(1)*AVECT(5)
    m3(2,2) = AVECT(2)*AVECT(2)
    m3(2,3) = AVECT(2)*AVECT(3)
    m3(2,4) = AVECT(2)*AVECT(4)
    m3(2,5) = AVECT(2)*AVECT(5)
    m3(3,3) = AVECT(3)*AVECT(3)
    m3(3,4) = AVECT(3)*AVECT(4)
    m3(3,5) = AVECT(3)*AVECT(5)
    m3(4,4) = AVECT(4)*AVECT(4)
    m3(4,5) = AVECT(4)*AVECT(5)
    m3(5,5) = AVECT(5)*AVECT(5)
    m3(2,1) = m3(1,2)
    m3(3,1) = m3(1,3)
    m3(4,1) = m3(1,4)
    m3(5,1) = m3(1,5)
    m3(3,2) = m3(2,3)
    m3(4,2) = m3(2,4)
    m3(5,2) = m3(2,5)
    m3(4,3) = m3(3,4)
    m3(5,3) = m3(3,5)
    m3(5,4) = m3(4,5)
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
!
! Computes FNORM = (F - F_0)/F_0 and uniaxial yield stress UNIYLD
! (see GTECTON notes). PSTRS is the scalar stress for specific yield
! criteria that is compared to UNIYLD
!
USE CONSTANTS, only: ONE,TWO,THREE,SIX,root3,HALF, deg2rad, dfmin,dlogmin,dfmax,dlogmax


implicit none
integer NSTR
parameter (NSTR=5)
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
frangl = PRPLAS(2) * DEG2RAD
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
    debar = BETB(1)**TWO + BETB(2)**TWO + HALF*BETB(3)**TWO + &
     HALF*BETB(4)**TWO + HALF*BETB(5)**TWO
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
subroutine FLOWVP (AVECT,STN,PRPLAS,PRPMAT,PSTRS,BETB)

USE ALGEBRA,   only: clear
use constants, only: one
!
! Computes BETB by multiplying AVECT with fluidity and PHI
!
implicit none
integer, parameter :: NSTR=5
!-pass
double precision :: AVECT(NSTR),STN(NSTR),PRPLAS(9),BETB(NSTR),PRPMAT(6)
!-local
double precision :: cmult, fluidy, yldpwr, nflow, FNORM, uniyld, PSTRS
integer :: i
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

USE MATERIALSMODULE, only: LMAT, &
                           MaterialIDOfElement
USE MODELDEFINITION, only: numwnk, lgdef, numslp
USE ALGEBRA
use constants,       only: eps, useDXE, half
use timestepmodule,  only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: elementSlipperyID
#endif

!
! Subroutine to compute Winkler restoring pressures from displacements
! and add them to the global force vector BTOT
!
implicit none
integer, parameter :: NSD=2,NDOF=3,NEN=4
!-pass
integer          :: iwelm, iwside, iwtime, ien, mat, lm, lmf, lmx
integer          :: nslip
double precision :: wpres, x, prpmat, btot, deld,d
double precision :: tfault, skew, dxe, dx

dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
 MAT(*),X(NSD,*),PRPMAT(6,*),LM(NDOF,NEN,*),BTOT(*), &
 DELD(NDOF,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN,*),LMX(NDOF,NEN,*),NSLIP(5,*), &
 DX(NDOF,*)
!-locl
double precision :: angle, thick
double precision :: anglevec(2)
integer :: elemID
double precision :: w(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),rot(3,3)
double precision :: tmp(NSD)

integer i,ii,j,n,m,nn,ni,nj,ll

logical DIFWNK
save DIFWNK
data DIFWNK/.true./
!
if (NUMWNK.eq.0) then
    return
endif
!
do 400 nn=1,IABS(NUMWNK)

    n = IABS(IWELM(nn))
    m = MaterialIDOfElement(n)

    LL= IWTIME(nn)
    if (LL.lt.0) then
        if (LL.ne.-1 .and. NSTEP.ge.-LL) then
            goto 400
        endif
    else if (LL.eq.0) then
        goto 400
    else
        if (NSTEP.lt.(LL-1)) then
            goto 400
        endif
    endif

    call LCOORD (X,XL,IEN(1,n),0)

    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,30)
        call ADDSNE (DL,n,NDOF,NEN,useDXE)
        call REZONE (XL,DL,1d0)
    endif

    i = IWSIDE(nn)
    j = i+1

    if (IEN(3,n).ne.IEN(4,n)) then
        if (j.gt.NEN) then
            j=1
        endif
    else
        if (j.ge.NEN) then
            j=1
        endif
    endif

    xn(1) = XL(2,j) - XL(2,i)
    xn(2) = XL(1,i) - XL(1,j)

    angle = (-SKEW(1,IEN(i,n))-SKEW(1,IEN(j,n)))*HALF

    anglevec(1)=angle
    anglevec(2)=0d0

    if (ABS(angle).gt.EPS) then
        call FORMRT (anglevec,rot,NSD)
        call EQUATE (tmp,xn,NSD)
        call VMPRD (rot,tmp,xn,NSD,3)
    endif

    thick = PRPMAT(6,m)

    w(1) = ABS(xn(1)) * WPRES(nn) * thick * HALF
    w(2) = ABS(xn(2)) * WPRES(nn) * thick * HALF

    ii = LM(1,i,n)

    if (ii.gt.0) then
        BTOT(ii) = BTOT(ii) - w(1)*DELD(1,IEN(i,n))
    endif

    ii = LM(1,j,n)

    if (ii.gt.0) then
        BTOT(ii) = BTOT(ii) - w(1)*DELD(1,IEN(j,n))
    endif

    ii = LM(2,i,n)

    if (ii.gt.0) then
        BTOT(ii) = BTOT(ii) - w(2)*DELD(2,IEN(i,n))
    endif

    ii = LM(2,j,n)

    if (ii.gt.0) then
        BTOT(ii) = BTOT(ii) - w(2)*DELD(2,IEN(j,n))
    endif

    if (NUMSLP.eq.0 .or. (.not.DIFWNK)) then
        goto 400
    endif

    ni = 0
    nj = 0

    do ii = 1,NUMSLP
        if (NSLIP(2,ii).eq.IEN(i,n)) then
            ni = ni + 1
        endif
        if (NSLIP(2,ii).eq.IEN(j,n)) then
            nj = nj + 1
        endif
    enddo


#ifdef EXPERIMENTAL_ReducedSlipAlloc
    elemID = elementSlipperyID(n,65)
    if (elemID.gt.0) then

        ii = IABS(LMX(1,i,elemID))

        if (ii.gt.0) then
            BTOT(ii) = BTOT(ii) - w(1)*DX(1,IEN(i,elemID))/DBLE(ni)
        endif

        ii = IABS(LMX(1,j,elemID))

        if (ii.gt.0) then
            BTOT(ii) = BTOT(ii) - w(1)*DX(1,IEN(j,elemID))/DBLE(nj)
        endif

        ii = IABS(LMX(2,i,elemID))

        if (ii.gt.0) then
            BTOT(ii) = BTOT(ii) - w(2)*DX(2,IEN(i,elemID))/DBLE(ni)
        endif

        ii = IABS(LMX(2,j,elemID))

        if (ii.gt.0) then
            BTOT(ii) = BTOT(ii) - w(2)*DX(2,IEN(j,elemID))/DBLE(nj)
        endif

    endif

#else

    ii = IABS(LMX(1,i,n))

    if (ii.gt.0) then
        BTOT(ii) = BTOT(ii) - w(1)*DX(1,IEN(i,n))/DBLE(ni)
    endif

    ii = IABS(LMX(1,j,n))

    if (ii.gt.0) then
        BTOT(ii) = BTOT(ii) - w(1)*DX(1,IEN(j,n))/DBLE(nj)
    endif

    ii = IABS(LMX(2,i,n))

    if (ii.gt.0) then
        BTOT(ii) = BTOT(ii) - w(2)*DX(2,IEN(i,n))/DBLE(ni)
    endif

    ii = IABS(LMX(2,j,n))

    if (ii.gt.0) then
        BTOT(ii) = BTOT(ii) - w(2)*DX(2,IEN(j,n))/DBLE(nj)
    endif
#endif


! original:

!    ii = IABS(LMX(1,i,n))
!    if (ii.gt.0) BTOT(ii)=BTOT(ii)-w(1)*DX(1,IEN(i,n))/ &
!     DBLE(ni)
!    ii = IABS(LMX(1,j,n))
!    if (ii.gt.0) BTOT(ii)=BTOT(ii)-w(1)*DX(1,IEN(j,n))/ &
!     DBLE(nj)
!    ii = IABS(LMX(2,i,n))
!    if (ii.gt.0) BTOT(ii)=BTOT(ii)-w(2)*DX(2,IEN(i,n))/ &
!     DBLE(ni)
!    ii = IABS(LMX(2,j,n))
!    if (ii.gt.0) BTOT(ii)=BTOT(ii)-w(2)*DX(2,IEN(j,n))/ &
!     DBLE(nj)

400     continue
!
return
end
!-------------------------------------------------------------------------------
subroutine SHAREL (IEN,IJSUR,NUMEL,NEL)
use debugmodule,     only: xit

!
!       Finds elements neighbouring to element NEL.
!       Neighbouring elements are identified by two shared nodes.
!
implicit none

integer, parameter :: NSUR=4,NEN=4
!-pass
integer NEL,NUMEL
integer :: ijsur, ien
dimension IEN(NEN,*),IJSUR(NSUR,*)

!-locl
integer,external :: NINDX
integer          :: nn, isur, mm, n, m, ielm, ii
integer          :: nshare

!
nn   = NEN
!       set to last IJSUR(isur,NEL) <> 0 (last index that was written)
isur = NINDX (IJSUR,NEL)
if (IEN(3,NEL).eq.IEN(4,NEL)) then
    nn = 3
endif
!       loop other elements
do ielm=NEL+1,NUMEL
!           compare nodal points of element NEL with those of IELM > NEL
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
!                   copy info NEL into IJSUR(*,ielm)
!                   find free index
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

500     return
end
!-------------------------------------------------------------------------------
integer function NINDX (IJSUR,NEL)

use debugmodule,     only: xit

implicit none

integer, parameter :: NSUR=4
!-pass
integer NEL
integer :: ijsur, i
dimension IJSUR(NSUR,*)

!
!       find free index
i = 1
100     if (i.le.NSUR) then
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
!
return
end
!-------------------------------------------------------------------------------
subroutine FLUID (BETB,PRPMAT,PRPLAS)

USE TIMESTEPMODULE
USE CONSTANTS, only: one, two, three, four, six,  &
    root3, half, deg2rad, dfmin,dlogmin,dfmax,dlogmax


!
!       Routine to determine visco-plastic fluidity from
!       explicit time stepping limit (Cormeau) (LINEAR only)
!       (Ref.: Owen & Hinton, sec. 8.9)
!       ONLY WORKS WITH NEWTONIAN RELAXATION!
!
implicit none
integer, parameter :: NSTR=5
!-pass
double precision :: PRPMAT(6),PRPLAS(9),BETB(NSTR)
double precision :: e, pois, frangl, cohesn, stnhrd, fluidy
double precision :: debar
double precision :: YLDPWR
integer          :: nflow, nptype, sinphi

!-init
NPTYPE = INT(PRPLAS(1))

if (NPTYPE.eq.0) then
    return
endif

YLDPWR = PRPLAS(5)
NFLOW  = INT(PRPLAS(6))

if (NFLOW.ne.1 .or. ABS(YLDPWR-1d0).gt.dfmin) then
    return
endif

FRANGL = PRPLAS(2) * deg2rad
COHESN = PRPLAS(3)
STNHRD = PRPLAS(7)
FLUIDY = PRPLAS(9)
E      = PRPMAT(1)
POIS   = PRPMAT(2)
!
DEBAR = BETB(1)*BETB(1) + BETB(2)*BETB(2) + HALF*BETB(3)*BETB(3) &
   + HALF*BETB(4)*BETB(4) + HALF*BETB(5)*BETB(5)
DEBAR = TWO * DEBAR / THREE

if (DEBAR.gt.dfmin) then
    DEBAR = SQRT(DEBAR)
endif
!
if (NPTYPE.eq.1) then
!           Tresca
    continue
else if (NPTYPE.eq.2) then
!           Von Mises
    continue
else if (NPTYPE.eq.3) then
!           Mohr-Coulomb
    COHESN = COHESN * COS(FRANGL)
else if (NPTYPE.eq.4) then
!           Drucker-Prager
    COHESN = SIX * COHESN * COS(FRANGL) / &
     ( root3 * (THREE - SIN(FRANGL)) )
endif

if (STNHRD.gt.dfmin) then
    COHESN = COHESN + STNHRD*DEBAR
endif
!
if (NPTYPE.eq.1) then
!           Tresca
    FLUIDY = (ONE+POIS) * COHESN / (E * DELTP)
else if (NPTYPE.eq.2) then
!           Von Mises
    FLUIDY = FOUR * (ONE+POIS) * COHESN / (THREE * E * DELTP)
   FLUIDY = FLUIDY*1D-2
else if (NPTYPE.eq.3) then
!           Mohr-Coulomb
    FRANGL = PRPLAS(2)
    SINPHI = SIN(FRANGL*deg2rad)
    SINPHI = SINPHI*SINPHI
    FLUIDY = FOUR*(ONE+POIS)*(ONE-TWO*POIS)*COHESN / &
     ((ONE-TWO*POIS+SINPHI)*E*DELTP)
else if (NPTYPE.eq.4) then
!           Drucker-Prager
    continue
endif

if (FLUIDY.gt.dfmin) PRPLAS(4) = FLUIDY
!
return
end
!-------------------------------------------------------------------------------
subroutine BPMATRIX (XL,DL,QUAD,BP,ierr)

use constants, only: zero
!
! Sets up derivatives matrix BP(i,j) = dUi/dXj
!
implicit none
integer, parameter :: NDOF=3, NSD=2, NEN=4
!-pass
integer          :: ierr
logical          :: QUAD
double precision :: XL(NSD,NEN),DL(NDOF,NEN),BP(NDOF,NSD)
!-locl
double precision :: SH(3,NEN),XS(NSD,NSD),det
!
call SHAP20 (ZERO,ZERO,XL,DET,SH,XS,QUAD,ierr,3)
if (ierr.ne.0) then
    return
endif
!
BP(1,1) = SH(1,1)*DL(1,1)+SH(1,2)*DL(1,2) + &
          SH(1,3)*DL(1,3)+SH(1,4)*DL(1,4)
BP(1,2) = SH(2,1)*DL(1,1)+SH(2,2)*DL(1,2) + &
          SH(2,3)*DL(1,3)+SH(2,4)*DL(1,4)
BP(2,1) = SH(1,1)*DL(2,1)+SH(1,2)*DL(2,2) + &
          SH(1,3)*DL(2,3)+SH(1,4)*DL(2,4)
BP(2,2) = SH(2,1)*DL(2,1)+SH(2,2)*DL(2,2) + &
          SH(2,3)*DL(2,3)+SH(2,4)*DL(2,4)
BP(3,1) = SH(1,1)*DL(3,1)+SH(1,2)*DL(3,2) + &
          SH(1,3)*DL(3,3)+SH(1,4)*DL(3,4)
BP(3,2) = SH(2,1)*DL(3,1)+SH(2,2)*DL(3,2) + &
          SH(2,3)*DL(3,3)+SH(2,4)*DL(3,4)
!
return
end
!-------------------------------------------------------------------------------
subroutine sdimnt (ISURF,SURF,B,ID,GRAV,PRPSED,DELTP)
USE MODELDEFINITION
!
implicit none

integer, parameter :: NSD=2, NDOF=3
!-pass
integer          :: isurf, id
double precision :: surf, b, grav, prpsed, deltp

dimension ISURF(*),SURF(NSD,*),B(*),ID(NDOF,*),GRAV(*), &
 PRPSED(2)


if (NSURF.le.0 .or. NSED.eq.0) then
    return
endif
!
write(stderr,10)
10   format(1x,'"SDIMNT" NOT IMPLEMENTED IN OPN VERSION')
!
return
end
!-----------------------------------------------------------------------
subroutine SIDENP (ISS,QUAD,LIEN)
use debugmodule,     only: xit

!
! returns local node point numbers in LIEN, belonging to side IS
! nodal points in LIEN are positively oriented.
!
! LWI replaced IS by ISS, because IS has become reserved

implicit none
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
        if (LIEN(2).gt.4) LIEN(2) = 1
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
function SRFHGT (SURF,NSURF,X)

USE ALGEBRA, only: linlin
!
! Routine to interpolate surface height at (X,Y) from surface nodal point data.
!
implicit none
integer, parameter :: NSD=2
!-pass
integer          :: nsurf
double precision :: surf, x
dimension SURF(NSD,*)

#include "coord.i"
!
!-locl
integer          :: i,k,m
double precision :: dx1, dx2,srfhgt, a,b
logical FIRST
save FIRST
data FIRST /.true./
srfhgt = 0d0
!
if (NSURF.le.0) return
!
if (FIRST) then
!           Check whether coordinates are uniformly in- or decreasing
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

200 call linlin(SURF(IX,i),SURF(IY,i),SURF(IX,i+1),SURF(IY,i+1),a,b,m)
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

 USE MODELDEFINITION, only: NUMEL, numslp
 use constants, only: eps
use debugmodule,     only: xit

!
!       Geometry-dependent part of elwinx, where stresses are converted
!       to forces. Updates IWINX, WINX and NWINKX.
!
implicit none
integer, parameter :: NEN=4,NDOF=2,NSD=2
!-pass
integer          :: nslip, ien, iwinx, mode, nwinkx
double precision :: winx, x, f
dimension NSLIP(5,*),IEN(NEN,*),WINX(NDOF,*),IWINX(NDOF,*), &
 X(NSD,*),f(*)

!-locl
integer, external :: nextlu
integer :: iside, nel, np, i,k
double precision :: dx, dy
character(len=1)axes
integer          :: islp(2), IIEN(2)
double precision :: XL(NSD,NEN),DW(2)
dimension :: axes(2)
save axes
data axes/'x','y'/
!
!       Some preliminary checks
if (iside.lt.1 .or. iside.gt.4) then
    write(stderr,1) nel,iside
 1       format(1x,'elwinx: illegal side number for element ',I5, &
     ': iside=',I5)
    call xit(1," ")
endif

if (nel.lt.1 .or. nel.gt.NUMEL) then
    write(stderr,1) NUMEL,nel
 2       format(1x,'elwinx: element number outside range [1, ',I6, &
     ']: ',I6)
    call xit(1," ")
endif
!
!       Check that both face nodes of side iside occur in NSLIP
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

do 200 i=1,2
    np = IEN(IIEN(i),nel)
    islp(i) = 1
100         if (NSLIP(1,islp(i)).eq.nel.and.NSLIP(2,islp(i)).eq.np) &
     goto 200
    islp(i) = islp(i) + 1
    if (islp(i).gt.NUMSLP) then
        write(stderr,3) nel,np,iside
 3           format(1x,'elwinx: element-node ',I5,'-',I5,1x, &
         'not in NSLIP array (iside=',i1,')')
        call xit(1," ")
    endif
    goto 100
200     continue

call LCOORD (X,XL,IEN(1,nel),0)

dx = (XL(1,IIEN(1))-XL(1,IIEN(2)))*0.5d0
dy = (XL(2,IIEN(1))-XL(2,IIEN(2)))*0.5d0

DW(1) = ABS(dy*f(1)) + ABS(dx*f(2))
DW(2) = ABS(dx*f(1)) + ABS(dy*f(2))

do i=1,2
    np = IEN(IIEN(i),nel)
    do 300 k=1,NSD
        if (DW(k).le.EPS) then
            goto 300
        endif
        if (NSLIP(2+k,islp(i)).eq.0) then
            goto 300
        endif

        WINX(k,np) = WINX(k,np) + DW(k)

        if (IWINX(k,np).eq.0) then
            NWINKX = NWINKX + 1
            IWINX(k,np) = mode
            goto 300
        else if (IWINX(k,np).eq.mode) then
            goto 300
        endif
        write(stderr,4) nel,np,iside,axes(k)
 4           format(1x,'elwinx: conflicting MODE for element ', &
         I5,', node ',I5,' (iside=',I1,', ',A1, &
         '-component)')
        call xit(1," ")
300         continue
enddo
!
return
end
!-----------------------------------------------------------------------
subroutine SLipperySKEW (IDSLP,SKEW,NSELSD,IEN,X,D,TFAULT,DXE, &
    LMF)

USE MODELDEFINITION
USE CONSTANTS,       only: zero, pi, rad2deg
use debugmodule,     only: lustat, itest, idebug, iecho
use timestepmodule,  only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

!
! Updates skew angles for slippery nodes
!
implicit none
integer, parameter :: NSD=2, NEN=4, NDOF=3
!-pass
integer          :: idslp, nselsd, ien, lmf
double precision :: skew, x, d, tfault, dxe

dimension IDSLP(*),SKEW(2,*),NSELSD(2,*),IEN(NEN,*),X(NSD,*), &
   D(NDOF,*),TFAULT(NDOF,*),DXE(NDOF,NEN,*), &
   LMF(NDOF,NEN,*)
!-locl
integer, external :: nextlu
integer :: k, kk, n, nn, nels
double precision :: angle, a, aveangle, oldskew

logical firstentry,warning
integer nwait, lutmp
save firstentry,lutmp,nwait
integer :: lien(2)
!-init
data firstentry, nwait /.true.,50/
!
 if (LGDEF.eq.0 .or. NSLSKEW.eq.0 .or. NUMROT.eq.0) return
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
do 100 kk=1,NUMSN
    warning = .false.
    k = IDSLP(kk)

    if (SKEW(1,k).eq.ZERO) then
        goto 100
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
            if (nels.gt.0) then
                if (ABS(aveangle-a).ge.PI/2d0) then
                    if (aveangle.lt.ZERO .and. a.gt.ZERO) then
                        a = a-PI
                    elseif (aveangle.gt.ZERO .and. a.lt.ZERO) then
                        a = a+PI
                    else
                        write(stderr,10) NSELSD(2,nn),n,a*rad2DEG, &
                       aveangle*rad2DEG
10                           format(1x,'SLSKEW warning: side ',I1, &
                         ' of element ',I5,' has angle of ', &
                         F7.2,' deg. Average thusfar is ', &
                         F7.2,' deg.')
                        warning = .true.
                   endif
                endif
            endif

            angle = angle + a
            aveangle = angle / DBLE(nels)
        endif
    enddo

    if (nels.eq.0) then
        goto 100
    endif

    if (IDEBUG.ne.0) then
        oldskew = SKEW(1,k)
    endif

    SKEW(1,k) = angle / DBLE(nels)
!
    if (IDEBUG.ne.0) then
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

    if (MOD(NSTEP,nwait).eq.0 .or. warning) then
       write(lutmp,3) NSTEP,k,SKEW(1,k)*rad2DEG
 3       format(1x,'NSTEP = ',I6,'   SKEW(1,',I6,') = ',F16.12, &
     ' deg')
    endif

100    continue
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
USE CONSTANTS, only: pi, halfpi, one,eps, useDXE
!
! returns angle of tangent of deformed element side
!
implicit none
integer, parameter :: NSD=2, NEN=4, NDOF=3
!-pass
 integer :: elemID

integer          :: ien, lien
double precision :: x, d, skew, tfault, dxe
 dimension IEN(*),X(NSD,*),D(NDOF,*),SKEW(2,*), &
  TFAULT(NDOF,*),LIEN(2),DXE(NDOF,*)

!-locl
double precision :: angle, area
double precision :: xl(NSD,NEN),dl(NDOF,NEN),xn(NSD)
!
!       localize coordinate data
call LCOORD (X,xl,IEN,0)
!
!       deformed coordinates XL
call LDISP  (dl,D,IEN,NDOF,NEN)
call AddFaultDisplacement (dl,elemID,TFAULT,NDOF,NEN,31)
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
 1        format(1x,'TANGENT: collapsing nodes error')
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
implicit none
integer, parameter :: NEN=4,NSD=2
!-pass
integer :: ien
double precision :: x, xl
dimension IEN(NEN),X(NSD,*),XL(NSD,NEN)
integer :: elemID, j, k
!
do j=1,NEN
    k = IEN(j)
    XL(1,j) = X(1,k)
    XL(2,j) = X(2,k)
enddo
!
return
end
!-----------------------------------------------------------------------


double precision function dilate (S,POIS)

use constants, only: third
!
! computes volume change or volume change rate, depending on what S is.
!
! SNGL  implicit real (a-h,o-z)
implicit none
!-pass
double precision :: S, pois
dimension S(*)
!-locl
!
dilate = (S(1)+S(2))*THIRD
!
return
end
!-------------------------------------------------------------------------------
function effee (S,POIS,SIGNINFO)

use constants, only: dfmin, zero, third, half, eps

implicit none
!-pass
logical SIGNINFO
double precision :: s
dimension S(*)
double precision :: pois
!-locl
double precision :: dil, exx, eyy, ezz,effee, j2
!-init
effee = ZERO
!
dil = (S(1)+S(2))*THIRD
exx = S(1)-dil
eyy = S(2)-dil
ezz =     -dil

J2 = HALF*(exx*exx+eyy*eyy+ezz*ezz) + &
   S(3)*S(3)+S(4)*S(4)+S(5)*S(5)

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
subroutine strdev (S,POIS)

USE constants, only: third

!
! computes strain (rate) deviator
!
! SNGL  implicit real (a-h,o-z)
implicit none
!-pass
double precision :: S, pois
dimension S(*)
!-locl
double precision :: dil
!
dil = (S(1)+S(2))*THIRD
S(1) = S(1)-dil
S(2) = S(2)-dil
!
return
end
!-------------------------------------------------------------------------------
double precision function effstn (S,POIS,SIGNINFO)

 use constants, only: dfmin,zero,eps

!
implicit none
integer, parameter :: NSTR=5
!-pass
logical SIGNINFO
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
subroutine stndev (S,POIS)

USE MODELDEFINITION

implicit none
integer, parameter :: NSTR=5
!-pass
double precision :: S(NSTR)
double precision :: pois
!-locl
double precision, external :: pressure
double precision :: p

p = pressure(S,POIS)
S(1) = S(1) - p
S(2) = S(2) - p
!
return
end
!-------------------------------------------------------------------------------
double precision function pressure (S,POIS)
!
implicit none
integer, parameter :: NSTR=5
!-pass
double precision :: pois
double precision :: S(NSTR)
!
pressure = (S(1)+S(2))*(1D0+POIS)/3D0
!
return
end
