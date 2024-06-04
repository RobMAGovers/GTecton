!       Version 2009.2
!
!       These routines are the rheology-dependent portion of GTECTON
!    The routines here are for 3-dimensional displacements
!-------------------------------------------------------------------------------
! This library is for isoparametric linear tetrahedra (simplex) elements
!
!  Y
!  |
!  |
!  |
!  3
!  |\
!  | \            face 1: nodes 1 2 4
!  |__\2_____X         2:       1 3 2
!  1\ /                3:       1 4 3
!    \4                4:       2 3 4
!     \
!      Z
!
!
!-------------------------------------------------------------------------------
!    THE ROUTINES LISTED HERE ARE:
!
!    NAME                    FUNCTION
!
!
!    ADDFLX        NORMAL BOUNDARY FLUX CONTRIBUTIONS TO LOAD VECTOR
!    ADDPR        TRANSLATES APPLIED PRESSURES INTO LOAD VECTOR
!    ADDSTR        ADDS STRESSES TO LOAD VECTOR
!    ADDWNK        ADDS WINKLER PRESSURES TO STIFFNESS MATRIX
!    ADVSTF        CONVECTION UPDATE OF LOCAL CONDUCTION MATRIX
!    AJ1        FIRST INVARIANT OF STRESS OR STRAIN MATRIX
!    AJ2        SECOND INVARIANT OF DEVIATORIC MATRIX
!    AJ3        THIRD INVARIANT OF DEVIATORIC MATRIX
!    BDELD        STRAIN FROM DISPLACEMENTS
!    BMATRIX        STRAIN-DISPLACEMENT MATRIX
!    BPMATRIX    PARTIAL DISPLACEMENT DERIVATIVES MATRIX
!    CPSTIFF        LOCAL HEAT CAPACITY STIFFNESS MATRIX
!    DEVSTN        STRESS DEVIATOR
!    DILATE        VOLUME CHANGE
!    EFFEE        EFFECTIVE DEVIATORIC STRAIN
!    EFFSTN        EFFECTIVE DEVIATORIC STRESS
!    EFORCE        STRESS EQUIVALENT NODAL FORCES
!    EIGMAX        UPPER LIMIT ESTIMATE OF EIGENVALUE
!    ELCTR        ELEMENT CENTER COORDINATES
!    ELMTEST        TEST ELEMENT ORIENTATION
!    FCGEOM        GEOMETRICAL PROPERTIES OF AN ELEMENT FACE
!    FLOWVP        VISCOPLASTIC STRAIN RATE
!    FLUID        VISCOPLASTIC FLUIDITY
!    FORMBT        VISCOUS FLOW CORRECTION FOR STRAIN
!    FORMJB        JACOBEAN MATRIX FOR VISCOUS FLOW
!    GRAVLD        GRAVITATIONAL BODY FORCE
!    INTCON        DIMENSIONAL AND TOPOLOGIC DATA
!    INYLD        YIELD CRITERION
!    LCOORD        LOCALIZE COORDINATE DATA
!    LFLUX        ELEMENT HEAT FLUX
!    LHEAT        HEAT PRODUCTION CONTRIBUTION TO LOAD VECTOR
!    MATERL        MATERIAL MATRICES
!    NINDX        ELEMENT SHARING INDEX
!    PRESSURE    FROM FIRST INVARIANT OF STRESS MATRIX
!    REZONE        LOCAL REZONE FOR LARGE DEFORMATION
!    SDIMNT        SURFACE ELEVATION AND LOADS FROM SEDIMENTARY TRANSPORT
!    SHAREL        IDENTIFIES NEIGHBOURING ELEMENTS
!    SHPTET        SHAPE FUNCTION FOR 3-D SIMPLEX ELEMENTS
!    SIDENP        LOCAL NODE NUMBERS ON A SPECIFIC SIDE
!    SLSKEW        UPDATE SKEW ANGLES FOR SLIPPERY NODES
!    SRFHGT        HEIGHT TO THE FREE SURFACE
!    STIFF        MECHANICAL STIFFNESS MATRIX
!    STNDEV        DEVIATORIC STRESS
!    STNTEN        FULL STRESS TENSOR
!    STRDEV        DEVIATORIC STRAIN
!    TANGVP        VISCOPLASTIC JACOBIAN MATRIX
!    TCENTR        ELEMENT CENTER TEMPERATURES
!    TEMWNK        ADDS WINKLER FLUXES TO THERMAL STIFFNESS MATRIX
!    TSTIFF        LOCAL CONDUCTION STIFFNESS MATRIX
!    TSTRN        STRAIN DUE TO THERMAL EXPANSION
!    WINKLP        LOADS RESULTING FROM WINKLER PRESSURES
!    YIELDF        PLASTIC FLOW VECTOR
!    YIELDS        PLASTIC YIELD STRESS
!
!**********************************************************************
subroutine INTCON (NSD,NDOF,NSTR,NEN)

USE LABELMODULE
USE MODELDEFINITION, only : NTYPE, NSUR, NEE, NEC
!
! Set variables that depend on topology
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
!-pass
 integer NSD,NDOF,NSTR,NEN

 NTYPE = 6
 NSD   = 3
 NDOF  = 3
 NSTR  = 6
 NEN   = 4
 NSUR  = 4
 NEC   = NSD*NEN
 NEE   = NDOF*NEN
 LABELC(1)='  X '
 LABELC(2)='  Y '
 LABELC(3)='  Z '
 LABELD(1)=LABELC(1)
 LABELD(2)=LABELC(2)
 LABELD(3)=LABELC(3)
 LABELS(1)='SXX '
 LABELS(2)='SYY '
 LABELS(3)='SZZ '
 LABELS(4)='SXY '
 LABELS(5)='SXZ '
 LABELS(6)='SYZ '
 LABELR(1)='RXX '
 LABELR(2)='RYY '
 LABELR(3)='RZZ '
 LABELR(4)='RXY '
 LABELR(5)='RXZ '
 LABELR(6)='RYZ '
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
 TEMLBL(6)='Z-CONDUCT'
!
return
end subroutine
!-------------------------------------------------------------------------------
! tested 25/6/2003
 subroutine SHPTET (X,Y,Z,XL,DET,SH,IERR,IFLAG)
 USE CONSTANTS, only: sixth, one
!
! Computes shape functions and their derivatives for linear tetrahedra (simplex
! element).
!
 implicit none
 integer NEN,NSD
 parameter (NEN=4,NSD=3)
!-pass
 integer IERR,IFLAG
 double precision X,Y,Z,XL,SH,DET
 dimension XL(NSD,NEN),SH(4,NEN)

!-locl
 double precision A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3, &
   D4,RDET
!-init
 IERR = 0
!
 A1 = XL(1,2)*(XL(2,3)*XL(3,4)-XL(3,3)*XL(2,4))- &
       XL(1,3)*(XL(2,2)*XL(3,4)-XL(3,2)*XL(2,4))+ &
       XL(1,4)*(XL(2,2)*XL(3,3)-XL(3,2)*XL(2,3))
 A2 = XL(1,1)*(XL(3,3)*XL(2,4)-XL(2,3)*XL(3,4))- &
       XL(1,3)*(XL(3,1)*XL(2,4)-XL(2,1)*XL(3,4))+ &
       XL(1,4)*(XL(3,1)*XL(2,3)-XL(2,1)*XL(3,3))
 A3 = XL(1,1)*(XL(2,2)*XL(3,4)-XL(3,2)*XL(2,4))- &
       XL(1,2)*(XL(2,1)*XL(3,4)-XL(3,1)*XL(2,4))+ &
       XL(1,4)*(XL(2,1)*XL(3,2)-XL(3,1)*XL(2,2))
 A4 = XL(1,1)*(XL(3,2)*XL(2,3)-XL(2,2)*XL(3,3))- &
       XL(1,2)*(XL(3,1)*XL(2,3)-XL(2,1)*XL(3,3))+ &
       XL(1,3)*(XL(3,1)*XL(2,2)-XL(2,1)*XL(3,2))
 B1 = XL(3,3)*XL(2,4)-XL(2,3)*XL(3,4) - &
       XL(3,2)*XL(2,4)+XL(2,2)*XL(3,4) + &
       XL(3,2)*XL(2,3)-XL(2,2)*XL(3,3)
 B2 = XL(2,3)*XL(3,4)-XL(3,3)*XL(2,4) - &
       XL(2,1)*XL(3,4)+XL(3,1)*XL(2,4) + &
       XL(2,1)*XL(3,3)-XL(3,1)*XL(2,3)
 B3 = XL(3,2)*XL(2,4)-XL(2,2)*XL(3,4) - &
       XL(3,1)*XL(2,4)+XL(2,1)*XL(3,4) + &
       XL(3,1)*XL(2,2)-XL(2,1)*XL(3,2)
 B4 = XL(2,2)*XL(3,3)-XL(3,2)*XL(2,3) - &
       XL(2,1)*XL(3,3)+XL(3,1)*XL(2,3) + &
       XL(2,1)*XL(3,2)-XL(3,1)*XL(2,2)
 C1 = XL(1,3)*XL(3,4)-XL(3,3)*XL(1,4) - &
       XL(1,2)*XL(3,4)+XL(3,2)*XL(1,4) + &
       XL(1,2)*XL(3,3)-XL(3,2)*XL(1,3)
 C2 = XL(3,3)*XL(1,4)-XL(1,3)*XL(3,4) - &
       XL(3,1)*XL(1,4)+XL(1,1)*XL(3,4) + &
       XL(3,1)*XL(1,3)-XL(1,1)*XL(3,3)
 C3 = XL(1,2)*XL(3,4)-XL(3,2)*XL(1,4) - &
       XL(1,1)*XL(3,4)+XL(3,1)*XL(1,4) + &
       XL(1,1)*XL(3,2)-XL(3,1)*XL(1,2)
 C4 = XL(3,2)*XL(1,3)-XL(1,2)*XL(3,3) - &
       XL(3,1)*XL(1,3)+XL(1,1)*XL(3,3) + &
       XL(3,1)*XL(1,2)-XL(1,1)*XL(3,2)
 D1 = XL(2,3)*XL(1,4)-XL(1,3)*XL(2,4) - &
       XL(2,2)*XL(1,4)+XL(1,2)*XL(2,4) + &
       XL(2,2)*XL(1,3)-XL(1,2)*XL(2,3)
 D2 = XL(1,3)*XL(2,4)-XL(2,3)*XL(1,4) - &
       XL(1,1)*XL(2,4)+XL(2,1)*XL(1,4) + &
       XL(1,1)*XL(2,3)-XL(2,1)*XL(1,3)
 D3 = XL(2,2)*XL(1,4)-XL(1,2)*XL(2,4) - &
       XL(2,1)*XL(1,4)+XL(1,1)*XL(2,4) + &
       XL(2,1)*XL(1,2)-XL(1,1)*XL(2,2)
 D4 = XL(1,1)*XL(2,2)-XL(2,1)*XL(1,2) - &
       XL(1,1)*XL(2,3)+XL(2,1)*XL(1,3) + &
       XL(1,2)*XL(2,3)-XL(2,2)*XL(1,3)
!
 DET = A1+A2+A3+A4
 if (DET.le.0d0) goto 1000
 RDET = ONE/DET
!
 if (IFLAG.eq.1) then
!        Compute x-, y-, and z-derivatives of shape functions
SH(1,1) = B1*RDET
SH(1,2) = B2*RDET
SH(1,3) = B3*RDET
SH(1,4) = B4*RDET
SH(2,1) = C1*RDET
SH(2,2) = C2*RDET
SH(2,3) = C3*RDET
SH(2,4) = C4*RDET
SH(3,1) = D1*RDET
SH(3,2) = D2*RDET
SH(3,3) = D3*RDET
SH(3,4) = D4*RDET
 else
SH(4,1) = (A1+B1*X+C1*Y+D1*Z)*RDET
SH(4,2) = (A2+B2*X+C2*Y+D2*Z)*RDET
SH(4,3) = (A3+B3*X+C3*Y+D3*Z)*RDET
SH(4,4) = (A4+B4*X+C4*Y+D4*Z)*RDET
 endif
!
!    compute volume
 DET = DET*SIXTH
!
 return
!
1000    write(stderr,1) DET
 1    format(///1x,'Shape function fails! Determinant is ',1PE20.4)
 IERR = 1
 return
 end
!-------------------------------------------------------------------------------
   subroutine elmtest (X,IEN)
   USE MESHDATAMODULE
   USE MODELDEFINITION
   use constants, only: eps

!
! tests element definitions
!
 implicit none
 integer NEN,NSD
 parameter (NEN=4,NSD=3)
!-pass
 integer IEN
 double precision X
 dimension X(NSD,*),IEN(NEN,*)

!-locl
 integer i,k
 double precision xl,a1,a2,a3,a4,det
 dimension xl(NSD,NEN)
!
 if (NUMEL.lt.1) then
write(stderr,10)
10        format(1x,'NUMEL < 1')
   call exitp(1)
 endif
!

 write(*,*) 'tet; rank', getrank(), 'says: ', NUMEL, NUMNP, meshdatactx%nvlocal, meshdatactx%nelocal



 do i=1,NUMEL
k = IEN(1,i)
if (k.lt.1 .or. k.gt.NUMNP) then
write(stderr,20) i,k
20        format(1x,'element ',I8,': reference to node number ',I8)
call exitp(1)
endif
XL(1,1) = X(1,k)
XL(2,1) = X(2,k)
XL(3,1) = X(3,k)
k = IEN(2,i)
if (k.lt.1 .or. k.gt.NUMNP) then
write(stderr,20) i,k
call exitp(1)
endif
XL(1,2) = X(1,k)
XL(2,2) = X(2,k)
XL(3,2) = X(3,k)
k = IEN(3,i)
if (k.lt.1 .or. k.gt.NUMNP) then
write(stderr,20) i,k
call exitp(1)
endif
XL(1,3) = X(1,k)
XL(2,3) = X(2,k)
XL(3,3) = X(3,k)
k = IEN(4,i)
if (k.lt.1 .or. k.gt.NUMNP) then
write(stderr,20) i,k
call exitp(1)
endif
XL(1,4) = X(1,k)
XL(2,4) = X(2,k)
XL(3,4) = X(3,k)
a1 = xl(1,2)*(xl(2,3)*xl(3,4)-xl(3,3)*xl(2,4))- &
       xl(1,3)*(xl(2,2)*xl(3,4)-xl(3,2)*xl(2,4))+ &
       xl(1,4)*(xl(2,2)*xl(3,3)-xl(3,2)*xl(2,3))
a2 = xl(1,1)*(xl(3,3)*xl(2,4)-xl(2,3)*xl(3,4))- &
       xl(1,3)*(xl(3,1)*xl(2,4)-xl(2,1)*xl(3,4))+ &
       xl(1,4)*(xl(3,1)*xl(2,3)-xl(2,1)*xl(3,3))
a3 = xl(1,1)*(xl(2,2)*xl(3,4)-xl(3,2)*xl(2,4))- &
       xl(1,2)*(xl(2,1)*xl(3,4)-xl(3,1)*xl(2,4))+ &
       xl(1,4)*(xl(2,1)*xl(3,2)-xl(3,1)*xl(2,2))
a4 = xl(1,1)*(xl(3,2)*xl(2,3)-xl(2,2)*xl(3,3))- &
       xl(1,2)*(xl(3,1)*xl(2,3)-xl(2,1)*xl(3,3))+ &
       xl(1,3)*(xl(3,1)*xl(2,2)-xl(2,1)*xl(3,2))
     det = (a1+a2+a3+a4)/6d0
if (ABS(det).le.EPS) then
write(stderr,30) i
30        format(1x,'degenerate element ',I8)
call exitp(1)
endif
if (det.lt.0d0) then
write(stderr,40) i
40           format(1x,'element ',I8,' is negatively oriented (node', &
         ' order is wrong)')
call exitp(1)
endif
!
 enddo
!
 if (iecho.eq.2) write(stderr,50)
50    format(1x,'Elements tested ok')
 return
 end
!-------------------------------------------------------------------------------
! tested 26/6/2003
 subroutine MATERL (DMAT,E,POIS)
!
! Constructs the material matrix relating stress and strain.
!
 implicit none
 integer NSTR
 parameter (NSTR=6)
 double precision ONE,TWO,HALF
 parameter (ONE=1d0,TWO=2d0,HALF=.5d0)
!-pass
 double precision DMAT,E,POIS
 dimension DMAT(NSTR,NSTR)
!-locl
 double precision AM,AL
!
 call CLEAR(DMAT,36,"DMAT")
!    Compute Lame parameters AM(=MU) and AL(=LAMBDA)
 AM=E/(ONE+POIS)
 AL=AM*POIS/(ONE-TWO*POIS)
 AM=HALF*AM
 DMAT(1,1)=TWO*AM+AL
 DMAT(1,2)=AL
 DMAT(1,3)=AL
 DMAT(2,1)=AL
 DMAT(2,2)=DMAT(1,1)
 DMAT(2,3)=AL
 DMAT(3,1)=AL
 DMAT(3,2)=AL
 DMAT(3,3)=DMAT(2,2)
 DMAT(4,4)=AM
 DMAT(5,5)=AM
 DMAT(6,6)=AM
!
 return
 end
!-------------------------------------------------------------------------------
! no need to test
 subroutine BDELD (XL,DL,EE,QUAD,ierr)
!
! Subroutine to compute strains from displacements in each element
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NDOF=3, NSD=3, NSTR=6, NEN=4, NEE=12)
!-pass
 integer ierr
 dimension XL(NSD,NEN),DL(NEE),EE(NSTR)
 logical QUAD
!-locl
 dimension BP(NDOF,NSD)
!
!       Calculate partial displacement derivatives
call BPMATRIX (XL,DL,QUAD,BP,ierr)
if (ierr.ne.0) return
EE(1) = BP(1,1)
EE(2) = BP(2,2)
EE(3) = BP(3,3)
EE(4) = BP(1,2)+BP(2,1) ! fluid dynamical definition of shear strain
EE(5) = BP(1,3)+BP(3,1) ! fluid dynamical definition of shear strain
EE(6) = BP(3,2)+BP(2,3) ! fluid dynamical definition of shear strain
!
 return
 end
!-------------------------------------------------------------------------------
! tested 26/6/2003
 subroutine BMATRIX (B,SH)
!
! Defines the strain-displacement matrix B
!
 implicit none
 integer NEN,NSTR,NDOF,NEE
 parameter (NEN=4, NSTR=6, NDOF=3, NEE=NEN*NDOF)
 double precision ZERO
 parameter (ZERO=0d0)
!-pass
 double precision B,SH
 dimension B(NSTR,NEE),SH(4,NEN)
!-locl
 integer k,i
!
!    Column number k
 k = 0
 do i=1,NEN
k = k + 1
B(1,k) = SH(1,i)
B(2,k) = ZERO
B(3,k) = ZERO
B(4,k) = SH(2,i)
B(5,k) = SH(3,i)
B(6,k) = ZERO
k = k + 1
B(1,k) = ZERO
B(2,k) = SH(2,i)
B(3,k) = ZERO
B(4,k) = SH(1,i)
B(5,k) = ZERO
B(6,k) = SH(3,i)
k = k + 1
B(1,k) = ZERO
B(2,k) = ZERO
B(3,k) = SH(3,i)
B(4,k) = ZERO
B(5,k) = SH(1,i)
B(6,k) = SH(2,i)
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
! tested 25/6/2003
 subroutine GRAVLD (PL,XL,GRAV,DUM1,DENS,QUAD,IERR)
!
! Computes the contribution to the load from gravity body forces
!
 implicit none
 integer NEN,NSD,NDOF
 parameter (NEN=4, NSD=3, NDOF=3)
 double precision ZERO,QUART
 parameter (ZERO=0D0,QUART=0.25D0)
!-pass
 logical QUAD
 integer IERR
 double precision PL,XL,GRAV,DUM1,DENS
 dimension PL(NDOF,NEN),XL(NSD,NEN),GRAV(NDOF)
!-locl
 double precision DET,GRAV1,GRAV2,GRAV3,SH
 dimension SH(4,NEN)
!
!    Compute Jacobian determinant
 call SHPTET (ZERO,ZERO,ZERO,XL,DET,SH,IERR,1)
 if (IERR.ne.0) return
 GRAV1 = GRAV(1)*DENS*DET*QUART
 GRAV2 = GRAV(2)*DENS*DET*QUART
 GRAV3 = GRAV(3)*DENS*DET*QUART
 PL(1,1) = PL(1,1) + GRAV1
 PL(1,2) = PL(1,2) + GRAV1
 PL(1,3) = PL(1,3) + GRAV1
 PL(1,4) = PL(1,4) + GRAV1
 PL(2,1) = PL(2,1) + GRAV2
 PL(2,2) = PL(2,2) + GRAV2
 PL(2,3) = PL(2,3) + GRAV2
 PL(2,4) = PL(2,4) + GRAV2
 PL(3,1) = PL(3,1) + GRAV3
 PL(3,2) = PL(3,2) + GRAV3
 PL(3,3) = PL(3,3) + GRAV3
 PL(3,4) = PL(3,4) + GRAV3
!
 return
 end
!-------------------------------------------------------------------------------
! tested 26/6/2003
 subroutine STIFF (DMAT,XL,SL,DUM1,QUAD,IERR)
!
! Computes the super-diagonal part of the local stiffness matrix.
!
 implicit none
 integer NEN,NSD,NDOF,NSTR,NEE
 parameter (NEN=4, NSD=3, NDOF=3, NSTR=6, NEE=NEN*NDOF)
!-pass
 integer IERR
 logical QUAD
 double precision DMAT,XL,SL,DUM1
 dimension XL(NSD,NEN),DMAT(NSTR,NSTR),SL(NEE,NEE)
!-locl
 double precision ZERO
 parameter (ZERO=0d0)
 integer j,i
 double precision SH,B,DET,DB1,DB2,DB3,DB4,DB5,DB6
 dimension SH(4,NEN),B(NSTR,NEE)
!
!    Compute spatial derivatives of shape functions and determinant
 call SHPTET (ZERO,ZERO,ZERO,XL,DET,SH,IERR,1)
 if (IERR.ne.0) return
!    Assemble B-matrix
 call BMATRIX (B,SH)
 do j=1,NEE
!        Calculate column J of DB-matrix
DB1=DMAT(1,1)*B(1,j)+DMAT(1,2)*B(2,j)+DMAT(1,3)*B(3,j) &
       +DMAT(1,4)*B(4,j)+DMAT(1,5)*B(5,j)+DMAT(1,6)*B(6,j)
DB2=DMAT(2,1)*B(1,j)+DMAT(2,2)*B(2,j)+DMAT(2,3)*B(3,j) &
       +DMAT(2,4)*B(4,j)+DMAT(2,5)*B(5,j)+DMAT(2,6)*B(6,j)
DB3=DMAT(3,1)*B(1,j)+DMAT(3,2)*B(2,j)+DMAT(3,3)*B(3,j) &
       +DMAT(3,4)*B(4,j)+DMAT(3,5)*B(5,j)+DMAT(3,6)*B(6,j)
DB4=DMAT(4,1)*B(1,j)+DMAT(4,2)*B(2,j)+DMAT(4,3)*B(3,j) &
       +DMAT(4,4)*B(4,j)+DMAT(4,5)*B(5,j)+DMAT(4,6)*B(6,j)
DB5=DMAT(5,1)*B(1,j)+DMAT(5,2)*B(2,j)+DMAT(5,3)*B(3,j) &
       +DMAT(5,4)*B(4,j)+DMAT(5,5)*B(5,j)+DMAT(5,6)*B(6,j)
DB6=DMAT(6,1)*B(1,j)+DMAT(6,2)*B(2,j)+DMAT(6,3)*B(3,j) &
       +DMAT(6,4)*B(4,j)+DMAT(6,5)*B(5,j)+DMAT(6,6)*B(6,j)
do i=1,j
SL(i,j)=SL(i,j)+DET*(B(1,i)*DB1+B(2,i)*DB2+B(3,i)*DB3 &
       +B(4,i)*DB4+B(5,i)*DB5+B(6,i)*DB6)
     enddo
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
! not tested yet
 subroutine ADDPR (B,X,D,IEN,LM,DUM1,DUM2,IELNO,ISIDE,PRES, &
   LMF,TFAULT,SKEW,DXE)
   USE MATERIALSMODULE, only : LMAT
   USE MODELDEFINITION
   use constants, only: eps

!
! Subroutine to add pressures to load vector
!
 implicit none
 integer NEN,NSD,NDOF
 parameter (NEN=4, NSD=3, NDOF=3)
 double precision ONE,THIRD,HALF
 parameter (THIRD=1d0/3d0, HALF=.5d0, ONE=1d0)
!-pass
 integer IEN,LM,IELNO,ISIDE,LMF
 double precision B,X,D,DUM1,DUM2,PRES,TFAULT,SKEW,DXE
 dimension IELNO(*),DUM1(*),IEN(NEN,*),X(NSD,*),D(NDOF,*), &
   ISIDE(*),PRES(*),DUM2(*),B(*),LM(NDOF,NEN,*), &
 LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
 DXE(NDOF,NEN,*)

!-locl
 integer k,n,lien,l,ii
 double precision xl,dl,area,xn,f,angle,tmp,rot
 dimension f(NDOF),xl(NSD,NEN),dl(NDOF,NEN),lien(3),xn(NSD), &
   angle(2),tmp(NSD),rot(3,3)
!
 do k=1,IABS(NUMPR)
!        Set element number N
n = IELNO(k)
!        Localize coordinates
call LCOORD (X,xl,IEN(1,n))
    if (LGDEF.ne.0) then
       call LDISP  (dl,D,IEN(1,n),NDOF,NEN)
       call AddFaultDisplacement (dl,LMF(1,1,n),TFAULT,NDOF,NEN)
       call ADDSNE (dl,DXE(1,1,n),NDOF,NEN)
       call REZONE(xl,dl,ONE)
    endif
!        Get local nodal point numbers of the side on which pressure works
call SIDENP(ISIDE(k),.true.,lien)
!        Calculate surface area and unit normal vector
call FCGEOM (xl,lien,area,xn)
!        Define surface mean of SKEW angles
angle(1) = ( -SKEW(1,IEN(lien(1),n)) - &
       SKEW(1,IEN(lien(2),n))-SKEW(1,IEN(lien(3),n)) ) * HALF
angle(2) = ( -SKEW(2,IEN(lien(1),n)) - &
       SKEW(2,IEN(lien(2),n))-SKEW(2,IEN(lien(3),n)) ) * HALF
if (ABS(angle(1)).gt.EPS .or. ABS(angle(2)).gt.EPS) then
       call FORMRT (angle,rot,NDOF)
       call EQUATE (tmp,xn,NDOF)
       call VMPRD (rot,tmp,xn,NDOF,3)
    endif
f(1) = xn(1) * PRES(k) * area * THIRD
    f(2) = xn(2) * PRES(k) * area * THIRD
    f(3) = xn(3) * PRES(k) * area * THIRD
do l=1,NDOF
ii = LM(l,lien(1),n)
if (ii.ne.0) B(ii) = B(ii) + f(l)
ii = LM(l,lien(2),n)
if (ii.ne.0) B(ii) = B(ii) + f(l)
ii = LM(l,lien(3),n)
if (ii.ne.0) B(ii) = B(ii) + f(l)
enddo
 enddo
!
if (iecho.eq.1) write(stdout,1)
 1   format(1x,'Pressure loads are being added')
!
 return
 end
!-----------------------------------------------------------------------
 subroutine SIDENP (is,QUAD,lien)
!
! routine sets (local) nodal point numbers belonging to side IS.
! NPs in lien are positively oriented.
!
 implicit none
!-pass
 integer is,lien
 logical QUAD
 dimension lien(3)
!
 if (is.eq.1) then
lien(1) = 1
lien(2) = 2
lien(3) = 4
 elseif (is.eq.2) then
lien(1) = 1
lien(2) = 3
lien(3) = 2
 elseif (is.eq.3) then
lien(1) = 1
lien(2) = 4
lien(3) = 3
 elseif (is.eq.4) then
lien(1) = 2
lien(2) = 3
lien(3) = 4
 else
write(stderr,1) is
 1        format(1x,'SIDENP: a tetrahedron has 4 sides, not ',I2)
call xit(1," ")
 endif
!
 return
 end
!-----------------------------------------------------------------------
! tested 26/5/2003
 subroutine FCGEOM (XL,LIEN,AREA,XN)
 use constant, only: dfmin,dlogmin,dfmax,dlogmax

!
! Calculates area of an element face, perpendicular to face normal,
! and the unit normal vector.
!
 implicit none
 integer NSD,NEN
 parameter (NSD=3,NEN=4)
 double precision TWO
 parameter (TWO=2d0)
!-pass
 integer LIEN
 double precision XL,AREA,XN
 dimension XL(NSD,NEN),LIEN(3),XN(NSD)

!-locl
 integer i
 double precision XLL
 dimension XLL(NSD,4)
!
!    Localize face coordinates
 do i=1,3
XLL(1,i) = XL(1,LIEN(i))
XLL(2,i) = XL(2,LIEN(i))
XLL(3,i) = XL(3,LIEN(i))
 enddo
!
 XN(1)=(XLL(2,2)-XLL(2,1))*(XLL(3,3)-XLL(3,1)) &
    -(XLL(3,2)-XLL(3,1))*(XLL(2,3)-XLL(2,1))
 XN(2)=(XLL(3,2)-XLL(3,1))*(XLL(1,3)-XLL(1,1)) &
    -(XLL(1,2)-XLL(1,1))*(XLL(3,3)-XLL(3,1))
 XN(3)=(XLL(1,2)-XLL(1,1))*(XLL(2,3)-XLL(2,1)) &
    -(XLL(2,2)-XLL(2,1))*(XLL(1,3)-XLL(1,1))
 AREA=SQRT(XN(1)*XN(1)+XN(2)*XN(2)+XN(3)*XN(3))
 if (ABS(AREA).lt.dfmin) then
write(stderr,1)
 1        format(1x,'FCGEOM: zero element face area')
  return
 endif
 XN(1)=XN(1)/AREA
 XN(2)=XN(2)/AREA
 XN(3)=XN(3)/AREA
 AREA=AREA/TWO
!
 return
 end
!-------------------------------------------------------------------------------
! not tested
 subroutine LFLUX (XL,TL,IEN,PRPTEM,FLUX,IERR)
!
! Routine to calculate element heat flux
!
 implicit none
 integer NSD,NEN,NDOF
 parameter (NSD=3,NEN=4,NDOF=3)
!-pass
 integer IEN,IERR
 double precision XL,TL,PRPTEM,FLUX
 dimension XL(NSD,NEN),TL(NEN),PRPTEM(6),FLUX(NSD),IEN(NEN)

!-locl
 double precision DET,SH,GRAD
 dimension GRAD(NDOF),SH(4,NEN)
!
!    Compute temperature gradient
 call SHPTET (0d0,0d0,0d0,XL,DET,SH,IERR,1)
 if (IERR.ne.0) return
!
 GRAD(1)=TL(1)*SH(1,1)+TL(2)*SH(1,2)+TL(3)*SH(1,3)+TL(4)*SH(1,4)
 GRAD(2)=TL(1)*SH(2,1)+TL(2)*SH(2,2)+TL(3)*SH(2,3)+TL(4)*SH(2,4)
 GRAD(3)=TL(1)*SH(3,1)+TL(2)*SH(3,2)+TL(3)*SH(3,3)+TL(4)*SH(3,4)
 FLUX(1) = -PRPTEM(4)*GRAD(1)
 FLUX(2) = -PRPTEM(5)*GRAD(2)
 FLUX(3) = -PRPTEM(6)*GRAD(3)
!
 return
 end
!-------------------------------------------------------------------------------
! tested 26/5/2003
 subroutine EFORCE (EVP,PL,XL,DUM1,QUAD,IERR)
!
! Computes stress-equivalent nodal forces PL
!
 implicit none
 integer NEN,NSD,NDOF,NSTR,NEE
 parameter (NEN=4,NSD=3,NDOF=3,NSTR=6,NEE=NEN*NDOF)
!-pass
 integer IERR
 logical QUAD
 double precision EVP,PL,XL,DUM1
 dimension XL(NSD,NEN),PL(NEE),EVP(NSTR)

!-locl
 integer i
 double precision SH,B,ZERO,DET
 parameter (ZERO=0d0)
 dimension SH(4,NEN),B(NSTR,NEE)
!
!    Compute Jacobian determinant and shape function derivatives
 call SHPTET (ZERO,ZERO,ZERO,XL,DET,SH,IERR,1)
 if (IERR.ne.0) return
 call BMATRIX (B,SH)
 do i=1,NEE
PL(i) = PL(i) + DET*(B(1,i)*EVP(1)+B(2,i)*EVP(2)+ &
       B(3,i)*EVP(3)+B(4,i)*EVP(4)+B(5,i)*EVP(5)+B(6,i)*EVP(6))
 enddo
!
 return
 end

!-------------------------------------------------------------------------------
! no need to test
 subroutine DEVSTN (STN,SDEV)
!
! Computes stress deviatior SDEV
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=6)
 parameter (TWO=2d0, THIRD=1d0/3d0)
!-pass
 dimension STN(NSTR),SDEV(NSTR)
!
 SDEV(1)=(TWO*STN(1)-STN(2)-STN(3))*THIRD
 SDEV(2)=(TWO*STN(2)-STN(1)-STN(3))*THIRD
 SDEV(3)=(TWO*STN(3)-STN(1)-STN(2))*THIRD
 SDEV(4)=STN(4)
 SDEV(5)=STN(5)
 SDEV(6)=STN(6)
!
 return
 end
!-------------------------------------------------------------------------------
! no need to test
 subroutine FORMBT (STN,BETA,BETB,PRPMAT,PRPLAS)
   USE MODELCTX
 use constant, only: dfmin,dlogmin,dfmax,dlogmax

!
! Subroutine to compute strain rates from stress field. It defines BETA and
! evaluates the product BETA*STRESS. BETB is the viscoplastic strain rate.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NSTR=6)
 parameter (ZERO=0d0,THIRD=1D0/3D0,ONE=1d0,TWO=2d0,SIX=6D0)
!-pass
 dimension STN(NSTR),BETA(NSTR),BETB(NSTR),PRPMAT(6),PRPLAS(9)
!-locl
 dimension AVECT(NSTR)
!-init
 sigma = ZERO
 emhu  = PRPMAT(3)
 anpwr = PRPMAT(4)
!
 if (IDIG(ICVIS,4,2).ne.0) then
!        Calculate yield stress YIELD
call YIELDS (STN,PRPLAS,YIELD)
!        Construct plastic flow vector AVECT
call YIELDF (STN,PRPLAS,AVECT)
!        Evaluate viscoplastic strain rate matrix BETB
call FLOWVP (AVECT,STN,PRPLAS,YIELD,BETB)
 endif
!
 streff = ( STN(1)*STN(1) + STN(2)*STN(2) + STN(3)*STN(3) &
   - STN(1)*STN(2) - STN(1)*STN(3) - STN(2)*STN(3) )*THIRD &
   + STN(4)*STN(4) + STN(5)*STN(5) + STN(6)*STN(6)
 if (streff.gt.dfmin) then
streff = SQRT(streff)
sigma = (streff/emhu)**(anpwr-ONE)/(SIX*emhu)
 endif
!    viscous stress-strain rate tensor
 BETA(1) = sigma * (TWO*STN(1)-STN(2)-STN(3))
 BETA(2) = sigma * (TWO*STN(2)-STN(1)-STN(3))
 BETA(3) = sigma * (TWO*STN(3)-STN(1)-STN(2))
 BETA(4) = sigma * STN(4) * SIX
 BETA(5) = sigma * STN(5) * SIX
 BETA(6) = sigma * STN(6) * SIX
!
 if (IDIG(ICVIS,4,2).ne.0) then
!        add visco-plastic contributions
BETA(1) = BETA(1) + BETB(1)
BETA(2) = BETA(2) + BETB(2)
BETA(3) = BETA(3) + BETB(3)
BETA(4) = BETA(4) + BETB(4)
BETA(5) = BETA(5) + BETB(5)
BETA(6) = BETA(6) + BETB(6)
 endif
!
 return
 end
!-------------------------------------------------------------------------------
! no need to test
 subroutine FORMJB (BETAJ,STN,PRPMAT,PRPLAS)
   USE MODELCTX
 use constant, only: dfmin,dlogmin,dfmax,dlogmax

!
! Constructs the Jacobean of the stress-strain rate matrix BETA.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NSTR=6)
 parameter(THIRD=1d0/3d0,ONE=1d0,TWO=2d0,SIX=6d0)
!-pass
 dimension BETAJ(6,6),STN(NSTR),PRPLAS(9),PRPMAT(6)
!-locl
 dimension BETAJP(NSTR,NSTR)
!
!    compute effective stress
 streff = ( STN(1)*STN(1) + STN(2)*STN(2) + STN(3)*STN(3) &
 - STN(1)*STN(2) - STN(1)*STN(3) - STN(2)*STN(3) )*THIRD &
 + STN(4)*STN(4) + STN(5)*STN(5) + STN(6)*STN(6)
 if (streff.le.dfmin) then
call CLEAR (BETAJ,36,"BETAJ")
return
 endif
 streff = SQRT(streff)
!
!       Plastic flow?
if (IDIG(ICVIS,4,2).ne.0) then
!        Calculate yield stress YIELD
call YIELDS (STN,PRPLAS,YIELD)
!        Construct Jacobian due to plastic yielding BETAJP
call CLEAR (BETAJP,36,"BETAJP")
call TANGVP (STN,YIELD,PRPLAS,BETAJP)
 endif
!
 emhu  = PRPMAT(3)
 anpwr = PRPMAT(4)
 a1 = anpwr-ONE
 a2 = a1/SIX
 a3 = a1*SIX
 strinv = ONE/streff
 sxx = (TWO*STN(1)-STN(2)-STN(3))*strinv
 syy = (TWO*STN(2)-STN(1)-STN(3))*strinv
 szz = (TWO*STN(3)-STN(1)-STN(2))*strinv
 sxy = STN(4)*strinv
 sxz = STN(5)*strinv
 syz = STN(6)*strinv
 BETAJ(1,1) =  TWO + a2*sxx*sxx
 BETAJ(1,2) = -ONE + a2*sxx*syy
 BETAJ(1,3) = -ONE + a2*sxx*szz
 BETAJ(1,4) =        a1*sxx*sxy
 BETAJ(1,5) =        a1*sxx*sxz
 BETAJ(1,6) =        a1*sxx*syz
 BETAJ(2,2) =  TWO + a2*syy*syy
 BETAJ(2,3) = -ONE + a2*syy*szz
 BETAJ(2,4) =        a1*syy*sxy
 BETAJ(2,5) =        a1*syy*sxz
 BETAJ(2,6) =        a1*syy*syz
 BETAJ(3,3) =  TWO + a2*szz*szz
 BETAJ(3,4) =        a1*szz*sxy
 BETAJ(3,5) =        a1*szz*sxz
 BETAJ(3,6) =        a1*szz*syz
 BETAJ(4,4) =  SIX + a3*sxy*sxy
 BETAJ(4,5) =        a3*sxy*sxz
 BETAJ(4,6) =        a3*sxy*syz
 BETAJ(5,5) =  SIX + a3*sxz*sxz
 BETAJ(5,6) =        a3*sxz*syz
 BETAJ(6,6) =  SIX + a3*syz*syz
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
 BETAJ(6,1) = BETAJ(1,6)
 BETAJ(6,2) = BETAJ(2,6)
 BETAJ(6,3) = BETAJ(3,6)
 BETAJ(6,4) = BETAJ(4,6)
 BETAJ(6,5) = BETAJ(5,6)
 sigma=((streff/emhu)**a1)/(SIX*emhu)
!
 if (IDIG(ICVIS,4,2).ne.0) then
!        Add Jacobians of viscoplastic and viscous flow
do j=1,NSTR
    BETAJ(1,j) = sigma*BETAJ(1,j)+BETAJP(1,j)
    BETAJ(2,j) = sigma*BETAJ(2,j)+BETAJP(2,j)
    BETAJ(3,j) = sigma*BETAJ(3,j)+BETAJP(3,j)
    BETAJ(4,j) = sigma*BETAJ(4,j)+BETAJP(4,j)
    BETAJ(5,j) = sigma*BETAJ(5,j)+BETAJP(5,j)
    BETAJ(6,j) = sigma*BETAJ(6,j)+BETAJP(6,j)
enddo
 else
do j=1,NSTR
    BETAJ(1,j) = sigma*BETAJ(1,j)
    BETAJ(2,j) = sigma*BETAJ(2,j)
    BETAJ(3,j) = sigma*BETAJ(3,j)
    BETAJ(4,j) = sigma*BETAJ(4,j)
    BETAJ(5,j) = sigma*BETAJ(5,j)
    BETAJ(6,j) = sigma*BETAJ(6,j)
enddo
 endif
!
 return
 end
!-------------------------------------------------------------------------------
! not checked yet
 subroutine ADDSTR (B,X,D,IEN,LM,DUM1,DUM2,IELSTR,ISSIDE, &
   ISTR,LMF,TFAULT,SKEW,DXE)


USE MATERIALSMODULE, only: LMAT
use modeldatamodule, only: modeldatactx
USE MODELDEFINITION
use constants,       only: eps
use timestepmodule,  only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

!
! Subroutine to add stresses to load vector
!
 implicit none
 integer NSD,NDOF,NSTR,NEN
 parameter (NSD=3,NDOF=3,NSTR=6,NEN=4)
 double precision THIRD
 parameter (THIRD=0.33333333333333333d0)
!-pass
 integer IEN,LM,IELSTR,ISSIDE,ISTR,LMF
 double precision B,X,D,DUM1,DUM2,TFAULT,SKEW,DXE
 dimension DUM1(*),IEN(NEN,*),X(NSD,*),D(NDOF,*),DUM2(*),B(*), &
   LM(NDOF,NEN,*),IELSTR(*),ISSIDE(*),ISTR(2,*), &
   LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*),DXE(NDOF,NEN,*)
!-locl
 integer numadd,n,k,it0,it1,lien,i,np,l,ii
 double precision xl,dl,area,xn,angle,rot,tmp,sxx,syy,szz,sxy, &
   sxz,syz,f
 dimension xl(NSD,NEN),dl(NDOF,NEN),lien(3),xn(NSD),angle(2), &
   rot(3,3),tmp(NSD),f(NDOF)
!-init
 numadd = 0
!
 do 20 n=1,NUMSTR
k = IELSTR(n)
if (k.le.0) then
write(stderr,1) n
 1        format(1x,'ADDSTR: fatal IELSTR error. Index=',I5)
call stoper()
endif
it0 = ISTR(1,n)
it1 = ISTR(2,n)
if (DBLE(NSTEP-it0)*DBLE(NSTEP-it1).gt.0d0) goto 20
numadd = numadd + 1
call LCOORD (X,xl,IEN(1,k))
if (LGDEF.ne.0) then
    call LDISP  (dl,D,IEN(1,k),NDOF,NEN)
call AddFaultDisplacement (dl,LMF(1,1,k),TFAULT,NDOF,NEN)
    call ADDSNE (dl,DXE(1,1,k),NDOF,NEN)
call REZONE (xl,dl,1d0)
endif
!        Get nodal point numbers of the side on which pressure works
call SIDENP(ISSIDE(n),.true.,lien)
!        Calculate surface area and normal vector
call FCGEOM (xl,lien,area,xn)
!        Rotate face into local direction
angle(1) = 0d0
angle(2) = 0d0
do i=1,3
np = IEN(lien(i),k)
angle(1) = angle(1) - SKEW(1,np)
angle(2) = angle(2) - SKEW(2,np)
enddo
angle(1) = angle(1)*THIRD
angle(2) = angle(2)*THIRD
if (ABS(angle(1)).gt.EPS.or.ABS(angle(2)).gt.EPS) then
call FORMRT(angle,rot,NDOF)
call EQUATE (tmp,xn,NDOF)
call VMPRD (rot,tmp,xn,NDOF,3)
endif
 sxx = modeldatactx%STRS(1,n)
 syy = modeldatactx%STRS(2,n)
 szz = modeldatactx%STRS(3,n)
 sxy = modeldatactx%STRS(4,n)
 sxz = modeldatactx%STRS(5,n)
 syz = modeldatactx%STRS(6,n)
    f(1)=(sxx*xn(1)+sxy*xn(2)+sxz*xn(3))*area*THIRD
    f(2)=(sxy*xn(1)+syy*xn(2)+syz*xn(3))*area*THIRD
    f(3)=(sxz*xn(1)+syz*xn(2)+szz*xn(3))*area*THIRD
do l=1,NSD
ii=LM(l,lien(1),k)
      if(ii.ne.0) B(ii)=B(ii)+f(l)
      ii=LM(l,lien(2),k)
      if(ii.ne.0) B(ii)=B(ii)+f(l)
      ii=LM(l,lien(3),k)
      if(ii.ne.0) B(ii)=B(ii)+f(l)
enddo
20    continue
!
 if (iecho.eq.1 .and. numadd.gt.0) write(stdout,2)
 2    format(1x,'Global stress loads are being added')
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine ADDWNK (IWELM,IWSIDE,IWTIME,WPRES,IEN,MAT,X,D,PRPMAT, &
   A,LM,IDIAG,LMF,TFAULT,SKEW,DXE,LMX,NSLIP)
   USE MATERIALSMODULE, only : LMAT
   USE MODELDEFINITION
   use constants, only: eps
use timestepmodule, only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

!
! Subroutine to add Winkler pressures to diagonal elements of stiffness
! matrix
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NSD=3,NDOF=3,NEN=4)
 parameter (QUART=.25d0,THIRD=1d0/3d0,HALF=5D-1,ONE=1D0)
!-pass
 dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),A(*),IDIAG(NEQ),LM(NDOF,NEN,*), &
   D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
   DXE(NDOF,NEN,*),LMX(NDOF,NEN,*),NSLIP(5,*)
#ifdef SPARSE


   PetscMPIInt irank
#endif
!-locl
 dimension w(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),rot(3,3), &
   tmp(NSD),lien(3),angle(2)

 if (NUMWNK.eq.0) return

 do 400 nn=1,IABS(NUMWNK)
n = IABS(IWELM(nn))
LL= IWTIME(nn)
    if (LL.lt.0) then
    if (LL.ne.-1 .and. NSTEP.ge.-LL) goto 400
    elseif (LL.eq.0) then
        goto 400
    else
    if (NSTEP.lt.(LL-1)) goto 400
    endif
call LCOORD (X,XL,IEN(1,n))
if (LGDEF.ne.0) then
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
call AddFaultDisplacement (DL,LMF(1,1,n),TFAULT,NDOF,NEN)
    call ADDSNE (DL,DXE(1,1,n),NDOF,NEN)
call REZONE (XL,DL,ONE)
endif
!           Get nodal point numbers of the side on which pressure works
    call SIDENP (IABS(IWSIDE(nn)),.true.,lien)
!           Calculate surface area and unit normal vector
    call FCGEOM (XL,lien,area,xn)
if (NUMROT.ne.0) then
!               Rotate face normal into local direction
    np = IEN(lien(1),n)
    angle(1) = -SKEW(1,np)
    angle(2) = -SKEW(2,np)
    np = IEN(lien(2),n)
    angle(1) = angle(1)-SKEW(1,np)
    angle(2) = angle(2)-SKEW(2,np)
    np = IEN(lien(3),n)
    angle(1) = angle(1)-SKEW(1,np)
    angle(2) = angle(2)-SKEW(2,np)
    angle(1) = angle(1)*THIRD
    angle(2) = angle(2)*THIRD
    if (ABS(angle(1)).gt.EPS.or.ABS(angle(2)).gt.EPS) then
            call FORMRT (angle,rot,NSD)
            call EQUATE (tmp,xn,NSD)
            call VMPRD (rot,tmp,xn,NSD,3)
        endif
endif
    w(1) = ABS( xn(1) * area ) * WPRES(nn) * THIRD
    w(2) = ABS( xn(2) * area ) * WPRES(nn) * THIRD
    w(3) = ABS( xn(3) * area ) * WPRES(nn) * THIRD
    if (NUMSLP.gt.0) then
if (IWELM(nn).gt.0) then
    ni = 0
    nj = 0
    nk = 0
    nl = 0
            do ii = 1,NUMSLP
                if (NSLIP(2,ii).eq.IEN(lien(1),n)) ni = ni + 1
                if (NSLIP(2,ii).eq.IEN(lien(2),n)) nj = nj + 1
                if (NSLIP(2,ii).eq.IEN(lien(3),n)) nk = nk + 1
            enddo
else
    ni = 1
    nj = 1
    nk = 1
endif
    endif
!        add contributions to global stiffness matrix
do 300 idof=1,NDOF
if (ABS(w(idof)).lt.EPS) goto 300
do li=1,3
    ii = LM(idof,lien(li),n)
    if (ii.gt.0) then
#ifdef SPARSE
        call MatSetValuesLocal (M_A,1,ii-1,1,ii-1,w(idof), &
               ADD_VALUES,ierr)
        if (ierr.ne.0) then
        write(stderr,1) ierr
 1                format(1x,'ADDWNK: MatSetValuesLocal error ',I4)
        call stoper()
        endif
#else
        ii = IDIAG(ii)
        A(ii) = A(ii) + w(idof)
#endif
    endif
enddo
!
if (NUMSLP.eq.0) goto 300
ii = IABS(LMX(idof,lien(1),n))
if (ii.gt.0) then
    v = w(idof)/DBLE(ni)
#ifdef SPARSE
    call MatSetValuesLocal (M_A,1,ii-1,1,ii-1,v,ADD_VALUES, &
           ierr)
    if (ierr.ne.0) then
    write(stderr,1) ierr
    call stoper()
    endif
#else
    ii = IDIAG(ii)
    A(ii) = A(ii) + v
#endif
endif
ii = IABS(LMX(idof,lien(2),n))
if (ii.gt.0) then
    v = w(idof)/DBLE(nj)
#ifdef SPARSE
    call MatSetValuesLocal (M_A,1,ii-1,1,ii-1,v,ADD_VALUES, &
           ierr)
    if (ierr.ne.0) then
    write(stderr,1) ierr
    call stoper()
    endif
#else
    ii = IDIAG(ii)
    A(ii) = A(ii) + v
#endif
endif
ii = IABS(LMX(idof,lien(3),n))
if (ii.gt.0) then
    v = w(idof)/DBLE(nk)
#ifdef SPARSE
    call MatSetValuesLocal (M_A,1,ii-1,1,ii-1,v,ADD_VALUES, &
           ierr)
    if (ierr.ne.0) then
    write(stderr,1) ierr
    call stoper()
    endif
#else
    ii = IDIAG(ii)
    A(ii) = A(ii) + v
#endif
endif
300        continue
400    continue
!
 return
 end
!-------------------------------------------------------------------------------
! not tested yet
 subroutine ELCTR (XELM,IEN,X,D,LMF,TFAULT, &
 SKEW,DXE,FACTOR,IERR)
   USE MODELDEFINITION
   USE MATERIALSMODULE, only : RADIUS
!
! Calculates (deformed) element center coordinates XELM
!
 implicit none
 integer NSD,NDOF,NEN
 parameter (NSD=3, NDOF=3, NEN=4)
 double precision QUART
 parameter (QUART=0.25d0)
!-pass
 integer IEN,LMF,IERR
 double precision XELM,X,D,TFAULT,SKEW,DXE,FACTOR
 dimension XELM(NSD),IEN(NEN),X(NSD,*),D(NDOF,*), &
   LMF(NDOF,NEN),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN)
!-locl
 double precision xl,dl
 dimension xl(NSD,NEN),dl(NDOF,NEN)
!
!    localize coordinates
 call LCOORD (X,xl,IEN)
 if (FACTOR.gt.0.) then
call LDISP  (dl,D,IEN,NDOF,NEN)
call AddFaultDisplacement (dl,LMF,TFAULT,NDOF,NEN)
call ADDSNE (dl,DXE,NDOF,NEN)
call REZONE (xl,dl,FACTOR)
 endif
 XELM(1) = (xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4))*QUART
 XELM(2) = (xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4))*QUART
 XELM(3) = (xl(3,1)+xl(3,2)+xl(3,3)+xl(3,4))*QUART
!
 return
 end
!-------------------------------------------------------------------------------
! not tested yet
 subroutine TSTIFF (XL,QUAD,SL,C,DUM1,IERR)
!
! computes the local conduction stiffness matrix from integration at Barlow
! points. Solely calculates upper-diagonal matrix.
!
 implicit none
 integer NEN,NSD
 parameter (NEN=4, NSD=3)
!-pass
 integer IERR
 logical QUAD
 double precision C,XL,SL,DUM1
 dimension XL(NSD,NEN),C(NSD),SL(NEN,NEN)
!-locl
 double precision ZERO
 parameter (ZERO=0d0)
 integer j,i
 double precision SH,DET
 dimension SH(4,NEN)
!
!    Compute spatial derivatives of shape functions and determinant
 call SHPTET (ZERO,ZERO,ZERO,XL,DET,SH,IERR,1)
 if (IERR.ne.0) return
 do j=1,NEN
do i=1,j
SL(i,j) = SL(i,j) + DET*(C(1)*SH(1,i)*SH(1,j) &
       +C(2)*SH(2,i)*SH(2,j)+C(3)*SH(3,i)*SH(3,j))
     enddo
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
! not tested yet
 subroutine ADVSTF (S,XL,VL,PL,COND,RHOCP,TH,HEAT,QUAD,IERR)
!
! convection (advection) update of local conduction stiffness matrix for
! steady convection-conduction solution
!
 implicit none
 integer NUMOPT
 parameter (NUMOPT=1)
!
! NUMOPT = 0: no upwinding
! NUMOPT = 1: Mizukami upwind scheme
!
 integer NEN,NDOF,NSD
  parameter (NEN=4, NDOF=3, NSD=3)
 double precision ZERO
 parameter (ZERO=0d0)
!-pass
 logical QUAD
 integer IERR
 double precision S,XL,VL,PL,COND,RHOCP,TH,HEAT
 dimension XL(NSD,NEN),VL(NDOF,NEN),S(NEN,NEN),COND(NSD),PL(NEN)

!-locl
 integer j,i
 double precision sh,sm,sa,sd,det,w,vx,vy,vz,samax,cplump, &
   tau,qi
 dimension sh(4,NEN),sm(NEN,NEN),sa(NEN,NEN),sd(NEN,NEN)
 logical firstentry
save firstentry
!-init
 data firstentry /.true./
!
 if (firstentry) then
if (iecho.ne.0) then
    if (NUMOPT.eq.0) then
    write(stderr,5)
 5            format(1x,'No upwinding')
else
    write(stderr,10)
10            format(1x,'Upwind: Mizukami scheme')
endif
endif
firstentry = .false.
 endif
!
 call CLEAR(sa,16,"sa")
 if (NUMOPT.ne.0) then
!        consistent mass matrix SM
call CLEAR(sm,16,"sm")
call CPSTIFF (XL,QUAD,sm,RHOCP,TH,IERR)
if (IERR.ne.0) return
!        conduction stiffness matrix SD
call CLEAR(sd,16,"sd")
call TSTIFF (XL,QUAD,sd,COND,TH,IERR)
 endif
!
 call SHPTET (ZERO,ZERO,ZERO,XL,det,sh,IERR,1)
 if (IERR.ne.0) return
 w = det * 5.0D-2 * RHOCP
 vx = VL(1,1)+VL(1,2)+VL(1,3)+VL(1,4)
 vy = VL(2,1)+VL(2,2)+VL(2,3)+VL(2,4)
 vz = VL(3,1)+VL(3,2)+VL(3,3)+VL(3,4)
 do j=1,NEN
do i=1,NEN
sa(i,j) = w*( (vx+VL(1,i))*sh(1,j) + &
       (vy+VL(2,i))*sh(2,j) + (vz+VL(3,i))*sh(3,j) )
enddo
 enddo
 do i=1,NEN
do j=1,NEN
S(i,j)  = S(i,j) + sa(i,j)
sa(i,j) = ABS(sa(i,j))
enddo
 enddo
 if (NUMOPT.ne.0) then
call SYMMET(sm,NEN)
do i=1,NEN
samax = MAX(sa(i,1),sa(i,2),sa(i,3),sa(i,4))
!        element ii of lumped mass matrix
cplump = sm(i,1)+sm(i,2)+sm(i,3)+sm(i,4)
tau = 5.0D-1*(cplump/(samax+sd(i,i)))
qi = tau*(vx*sh(1,i)+vy*sh(2,i)+vz*sh(3,i))
PL(i) = PL(i) + HEAT*det*qi
w = det * 0.25D0 * RHOCP *qi
    do j=1,NEN
    S(i,j)=S(i,j)+w*(vx*sh(1,j)+vy*sh(2,j)+vz*sh(3,j))
enddo
enddo
 endif
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine LHEAT (XL,QUAD,PL,HEAT,DUM1,IERR)
!
! computes the local load vector contribution of heat production.
!
 implicit none
 integer NEN,NSD
  parameter (NEN=4,NSD=3)
 double precision ZERO
 parameter (ZERO=0d0)
!-pass
 integer IERR
 logical QUAD
 double precision XL,PL,HEAT,DUM1
 dimension XL(NSD,NEN),PL(NEN)
!-locl
 double precision det,sh,w
 dimension sh(4,NEN)
!
 call SHPTET (ZERO,ZERO,ZERO,XL,det,sh,IERR,1)
 if (IERR.ne.0) return
 w = det * HEAT * 0.25D0
 PL(1) = PL(1) + w
 PL(2) = PL(2) + w
 PL(3) = PL(3) + w
 PL(4) = PL(4) + w
!
 return
 end
!-------------------------------------------------------------------------------
 function EigMax (STN,E,V,EMHU,NPL,ALPHA)

 ! LWItodo: what does this do? eigmax is not argument. 

!
! Gives an upper limit estimate of the maximum local eigenvalue
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=6)
!-pass
 dimension STN(NSTR)

!
 EigMax = 1.d0
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine ADDFLX (BT,IEN,X,D,DUM1,DUM2,IFLS,IFLX,BFLX,LMT, &
 LMTX,LMF,TFAULT,SKEW,DXE)!,NFLX)
   USE MATERIALSMODULE, only : LMAT
   USE MODELDEFINITION
!
! Routine to add normal boundary fluxes to load vector
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NEN=4, NSD=3, NDOF=3)
 parameter (THIRD=1d0/3d0)
!-pass
!    integer NFLX
 dimension BT(*),IEN(NEN,*),X(NSD,*),DUM1(*),DUM2(*),IFLS(*), &
 IFLX(*),BFLX(*),LMT(NEN,*),D(NDOF,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),LMTX(NEN,*),SKEW(2,*),DXE(NDOF,NEN,*)
!-locl
 dimension XL(NSD,NEN),DL(NDOF,NEN),lien(3),XN(NSD)

!
 if (NFLX.le.0) return
!
 do n=1,NFLX
m = IFLX(n)
iside = IFLS(n)
call LCOORD (X,XL,IEN(1,m))
if (LGDEF.ne.0) then
    call LDISP  (DL,D,IEN(1,m),NDOF,NEN)
    call AddFaultDisplacement (DL,LMF(1,1,m),TFAULT,NDOF,NEN)
    call ADDSNE (DL,DXE(1,1,m),NDOF,NEN)
    call REZONE (XL,DL,1d0)
endif
    call SIDENP (iside,.true.,lien)
    call FCGEOM (XL,lien,AREA,XN)
    F=-AREA*BFLX(n)*THIRD
ii = LMT(lien(1),m)
if (ii.eq.0.and.NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) &
       ii = LMTX(lien(1),m)
if (ii.gt.0) BT(ii) = BT(ii) + F
ii = LMT(lien(2),m)
if (ii.eq.0.and.NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) &
       ii = LMTX(lien(2),m)
if (ii.gt.0) BT(ii) = BT(ii) + F
ii = LMT(lien(3),m)
if (ii.eq.0.and.NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) &
       ii = LMTX(lien(3),m)
if (ii.gt.0) BT(ii) = BT(ii) + F
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine TEMWNK (AA,LMT,ITDIAG,ITWINK,TWINK,IEN,MAT,PRPMAT, &
   X,D,LMF,TFAULT,SKEW,DXE)
   USE MATERIALSMODULE, only : LMAT
   USE MODELDEFINITION
use timestepmodule, only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

!
! Subroutine to add Winkler fluxes to diagonal elements of thermal stiffness
! matrix
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NSD=3,NDOF=3,NEN=4,ONE=1D0)
!-pass
 dimension ITWINK(3,*),TWINK(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),AA(*),ITDIAG(*),LMT(NEN,*), &
   D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
   DXE(NDOF,NEN,*)

#ifdef SPARSE
   PetscMPIInt irank
#endif
!-locl
 dimension xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),lien(3)
!
 if (NTWINK.eq.0) return
!
if (iecho.ne.0) write(stdout,*) 'Adding Winkler fluxes'
!
 do 400 nn=1,NTWINK
n = ITWINK(1,nn)
LL= ITWINK(3,nn)
    if (LL.lt.0) then
    if (LL.ne.-1 .and. NSTEP.ge.-LL) goto 400
    elseif (LL.eq.0) then
        goto 400
    else
    if (NSTEP.lt.(LL-1)) goto 400
    endif
call LCOORD (X,XL,IEN(1,n))
if (LGDEF.ne.0) then
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
call AddFaultDisplacement (DL,LMF(1,1,n),TFAULT,NDOF,NEN)
    call ADDSNE (DL,DXE(1,1,n),NDOF,NEN)
call REZONE (XL,DL,ONE)
endif
call SIDENP (ITWINK(2,nn),.true.,lien)
    call FCGEOM (XL,lien,area,xn)
wt = area * TWINK(nn)/3d0
ii = LMT(lien(1),n)
if (ii.gt.0) then
#ifdef SPARSE
call MatSetValuesLocal(T_A,1,ii-1,1,ii-1,wt,ADD_VALUES,ierr)
if (ierr.ne.0) then
    write(stderr,1) ierr
 1            format(1x,'TEMWNK: MatSetValuesLocal error ',I4)
    call stoper()
endif
#else
    ii = ITDIAG(ii)    ! matrix address of eqno
AA(ii) = AA(ii) + wt
#endif
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
ii = LMT(lien(3),n)
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
400    continue
!
 return
 end
!-----------------------------------------------------------------------
! not tested
 function Tcentr (T,TX,XL,IEN,LMTX,IDUM1)
   USE MODELDEFINITION
!
! Calculates the temperature in the center of the element
!
 implicit none
 integer NEN,NSD
 parameter (NEN=4, NSD=3)
 double precision QUART
 parameter (QUART=0.25d0)
!-pass
 integer IEN,LMTX,IDUM1
 double precision T,TX,XL,Tcentr
 dimension T(*),TX(*),IEN(NEN),LMTX(NEN),XL(NSD,NEN)
!-locl
 double precision TL
 dimension TL(NEN)
!
!    Localize temperatures
 call LTEMP (IEN,T,TL,NEN)
!
!    Modify for discontinuities
 if (NUMSLP+NUMFN.gt.0.and.IDIFT.eq.1) &
   call LTEMPX (LMTX,IEN,TX,TL,NEN)
!
!    Calculate center temperature
 Tcentr = (TL(1)+TL(2)+TL(3)+TL(4))*QUART
!
 return
 end
!-----------------------------------------------------------------------
! not tested
 subroutine CPSTIFF (XL,QUAD,S,RHOCP,DUM1,IERR)
!
! computes the upper local heat capacity stiffness matrix
!
 implicit none
 integer NEN,NSD
  parameter (NEN=4,NSD=3)
 double precision ZERO,TWO
 parameter (ZERO=0d0,TWO=2D0)
!-pass
 integer IERR
 logical QUAD
 double precision XL,S,RHOCP,DUM1
 dimension XL(NSD,NEN),S(NEN,NEN)
!-locl
 double precision det,sh,w
 dimension sh(4,NEN)
!
 call SHPTET (ZERO,ZERO,ZERO,XL,det,sh,IERR,1)
 if (IERR.ne.0) return
 w = det * 0.05D0 * RHOCP
 S(1,1) = S(1,1) + w*TWO
 S(1,2) = S(1,2) + w
 S(1,3) = S(1,3) + w
 S(1,4) = S(1,4) + w
 S(2,2) = S(2,2) + w*TWO
 S(2,3) = S(2,3) + w
 S(2,4) = S(2,4) + w
 S(3,3) = S(3,3) + w*TWO
 S(3,4) = S(3,4) + w
 S(4,4) = S(4,4) + w*TWO
!
 return
 end
!-------------------------------------------------------------------------------
 function AJ1 (ST)
!
! Calculates first invariant of stress or strain matrix
!
! SNGL  implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
!-pass
 dimension ST(6)
 AJ1 = ST(1)+ST(2)+ST(3)
!
 return
 end
!-------------------------------------------------------------------------------
 function AJ2 (ST)
!
! Calculates second deviatoric invariant of matrix ST
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (SIX=6.d0)
!-pass
 dimension ST(6)
!
 SXY = ST(1)-ST(2)
 SXZ = ST(1)-ST(3)
 SYZ = ST(2)-ST(3)
 AJ2 = (SXY*SXY+SXZ*SXZ+SYZ*SYZ)/SIX &
   +ST(4)*ST(4)+ST(5)*ST(5)+ST(6)*ST(6)
!
 return
 end
!-------------------------------------------------------------------------------
 function AJ3 (ST)
!
! Calculates third deviatoric invariant of matrix ST
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (TWO=2.d0,THIRD=1d0/3d0)
!-pass
 dimension ST(6)
!
 SDEV1 = (TWO*ST(1)-ST(2)-ST(3))*THIRD
 SDEV2 = (TWO*ST(2)-ST(1)-ST(3))*THIRD
 SDEV3 = (TWO*ST(3)-ST(1)-ST(2))*THIRD
 AJ3 = SDEV1*SDEV2*SDEV3-SDEV1*ST(6)*ST(6)-SDEV2*ST(5)*ST(5) &
   -SDEV3*ST(4)*ST(4)+TWO*ST(4)*ST(5)*ST(6)
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine REZONE (XL,DL,FACTOR)
!
! Routine to update local coordinate XL by displacements DL
!
 parameter (NSD=3, NDOF=3, NEN=4)
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
!-pass
 dimension XL(NSD,NEN),DL(NDOF,NEN)
!
XL(1,1) = XL(1,1) + DL(1,1)*FACTOR
XL(2,1) = XL(2,1) + DL(2,1)*FACTOR
XL(3,1) = XL(3,1) + DL(3,1)*FACTOR
XL(1,2) = XL(1,2) + DL(1,2)*FACTOR
XL(2,2) = XL(2,2) + DL(2,2)*FACTOR
XL(3,2) = XL(3,2) + DL(3,2)*FACTOR
XL(1,3) = XL(1,3) + DL(1,3)*FACTOR
XL(2,3) = XL(2,3) + DL(2,3)*FACTOR
XL(3,3) = XL(3,3) + DL(3,3)*FACTOR
XL(1,4) = XL(1,4) + DL(1,4)*FACTOR
XL(2,4) = XL(2,4) + DL(2,4)*FACTOR
XL(3,4) = XL(3,4) + DL(3,4)*FACTOR
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine TSTRN (EE,DT,VEXP)
!
! calculates strain due to thermal expansion
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=6)
 parameter (ZERO=0.d0, THIRD=1d0/3d0)
!-pass
 dimension EE(NSTR)
!
 EE(1) = DT*VEXP*THIRD
 EE(2) = EE(1)
 EE(3) = EE(1)
 EE(4) = ZERO
 EE(5) = ZERO
 EE(6) = ZERO
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine YIELDS (STN,PRPLAS,YIELD)

 USE CONSTANTS, only third, two, root3, deg2rad, dfmin,dlogmin,dfmax,dlogmax

!
! Routine to calculate the yield stress for various types of plasticity.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=6)
!-pass
 dimension STN(NSTR),PRPLAS(9)

!-locl
 external AJ1,AJ2
!
 STNEFF = AJ2(STN)
 if (STNEFF.gt.dfmin) STNEFF = SQRT(STNEFF)
 STMEAN = AJ1(STN) * THIRD
!
 NPTYPE = INT(PRPLAS(1))
 if (NPTYPE.eq.1) then
!        Tresca
YIELD  = TWO * STNEFF
 else if (NPTYPE.eq.2) then
!        Von Mises
YIELD  = root3 * STNEFF
 else if (NPTYPE.eq.3) then
!        Mohr-Coulomb
FRANGL = PRPLAS(2)
SINPHI = SIN(FRANGL*deg2rad)
YIELD  = STMEAN * SINPHI + STNEFF
 else if (NPTYPE.eq.4) then
!        Drucker-Prager
FRANGL = PRPLAS(2)
SINPHI = SIN(FRANGL*deg2rad)
YIELD  = 6.0d0 * STMEAN * SINPHI / (root3*(3.d0-SINPHI)) + &
       STNEFF
 else
write(stderr,1) NPTYPE
 1        format(1x,'YIELDS: unknown plasticity type ',I5)
call stoper()
 endif
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine YIELDF (STN,PRPLAS,AVECT)
!
 USE CONSTANTS, only: zero, one, two, three, third, eight,  &
    root3, deg2rad, THF, dfmin,dlogmin,dfmax,dlogmax

! Evaluate the plastic flow vector AVECT.(Owen & Hinton, sec. 7.4)
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=6)
!-pass
 dimension STN(NSTR),PRPLAS(9),AVECT(NSTR)
!-locl
 dimension SDEV(NSTR),A1(NSTR),A2(NSTR),A3(NSTR)
 external AJ2,AJ3
 save A1
!-init
 data A1/1.,1.,1.,0.,0.,0./
!
!    Compute second deviatoric invariant of stress matrix
 DJ2 = AJ2(STN)
 if (DJ2.gt.dfmin) goto 100
!    Zero stress: set flow vector to zero and return
 call CLEAR(AVECT,NSTR,"AVECT")
 return
!
100    TWOJ2I = ONE/(TWO*SQRT(DJ2))
!
!    Calculate stress deviator
 call DEVSTN (STN,SDEV)
!
!    Setup flow vectors A2 and A3
 A2(1) = SDEV(1)*TWOJ2I
 A2(2) = SDEV(2)*TWOJ2I
 A2(3) = SDEV(3)*TWOJ2I
 A2(4) = TWO*SDEV(4)*TWOJ2I
 A2(5) = TWO*SDEV(5)*TWOJ2I
 A2(6) = TWO*SDEV(6)*TWOJ2I
 A3(1) = SDEV(2)*SDEV(3)-SDEV(6)*SDEV(6)+DJ2*THIRD
 A3(2) = SDEV(1)*SDEV(3)-SDEV(5)*SDEV(5)+DJ2*THIRD
 A3(3) = SDEV(1)*SDEV(2)-SDEV(4)*SDEV(4)+DJ2*THIRD
 A3(4) = TWO*(SDEV(6)*SDEV(5)-SDEV(3)*SDEV(4))
 A3(5) = TWO*(SDEV(4)*SDEV(6)-SDEV(2)*SDEV(5))
 A3(6) = TWO*(SDEV(5)*SDEV(4)-SDEV(1)*SDEV(6))
!
!    Compute multipliers
 NPTYPE = INT(PRPLAS(1))
 if (NPTYPE.eq.1) then
!        Tresca
SIN3TH = THF*AJ3(STN)*EIGHT*TWOJ2I*TWOJ2I*TWOJ2I
if (SIN3TH.lt.-ONE) SIN3TH=-ONE
if (SIN3TH.gt. ONE) SIN3TH= ONE
THETA = THIRD*ASIN(SIN3TH)
C1 = ZERO
C2 = TWO*COS(THETA)*(ONE+TAN(THETA)*TAN(THREE*THETA))
C3 = root3/DJ2*SIN(THETA)/COS(THREE*THETA)
 else if (NPTYPE.eq.2) then
!        Von Mises
C1 = ZERO
C2 = root3
C3 = ZERO
 else if (NPTYPE.eq.3) then
!        Mohr-Coulomb
SIN3TH = THF*AJ3(STN)*EIGHT*TWOJ2I*TWOJ2I*TWOJ2I
if (SIN3TH.lt.-ONE) SIN3TH=-ONE
if (SIN3TH.gt. ONE) SIN3TH= ONE
THETA  = THIRD*ASIN(SIN3TH)
SINTH  = SIN(THETA)
COSTH  = COS(THETA)
TANTH  = TAN(THETA)
COS3TH = COS(THREE*THETA)
TAN3TH = TAN(THREE*THETA)
FRANGL = PRPLAS(2)
SINPHI = SIN(FRANGL*deg2rad)
C1 = SINPHI*THIRD
C2 = COSTH*(ONE+TANTH*TAN3TH+SINPHI/root3*(TAN3TH-TANTH))
C3 = (root3*SINTH+COSTH*SINPHI)/(TWO*DJ2*COS3TH)
 else if (NPTYPE.eq.4) then
!        Drucker-Prager
FRANGL = PRPLAS(2)
SINPHI = SIN(FRANGL*deg2rad)
C1 = TWO*SINPHI/(root3*(THREE-SINPHI))
C2 = ONE
C3 = ZERO
 endif
!
!    calculate A-vector
 AVECT(1)=C1*A1(1)+C2*A2(1)+C3*A3(1)
 AVECT(2)=C1*A1(2)+C2*A2(2)+C3*A3(2)
 AVECT(3)=C1*A1(3)+C2*A2(3)+C3*A3(3)
 AVECT(4)=C1*A1(4)+C2*A2(4)+C3*A3(4)
 AVECT(5)=C1*A1(5)+C2*A2(5)+C3*A3(5)
 AVECT(6)=C1*A1(6)+C2*A2(6)+C3*A3(6)
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine TANGVP (STN,YIELD,PRPLAS,BETAJP)

 USE CONSTANTS, only: one, two, four, threeq, root3, dfmin,dlogmin,dfmax,dlogmax,eps
!
! Compute the Jacobian matrix for viscoplastic flow BETAJP
! (Ref.: Owen & Hinton, sec. 8.7.2)
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=6)
!-pass
 dimension STN(NSTR),BETAJP(NSTR,NSTR),PRPLAS(9)

!-locl
 dimension SDEV(NSTR),AM1(NSTR,NSTR),AM2(NSTR,NSTR)
 save AM1
 external AJ2
!-init
 data AM1/0.666666666666666667d0,-0.333333333333333333d0, &
   -0.333333333333333333d0,0.d0,0.d0,0.d0,-0.333333333333333333d0, &
   0.666666666666666667d0,-0.333333333333333333d0,0.d0,0.d0,0.d0, &
   -0.333333333333333333d0,-0.333333333333333333d0, &
   0.666666666666666667d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,2.d0,0.d0, &
 0.d0,0.d0,0.d0,0.d0,0.d0,2.d0, &
   0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,2.d0/
 ALLOW  = 0.01d0
 NPTYPE = INT(PRPLAS(1))
 COHESN = PRPLAS(3)
 FLUIDY = PRPLAS(4)
 YLDPWR = PRPLAS(5)
 NFLOW  = INT(PRPLAS(6))
!
 if (ABS(COHESN).lt.EPS) then
write(stderr,1)
 1        format(1x,'TANGVP: zero plastic cohesion error')
call stoper()
 endif
 if (NPTYPE.eq.1) then
write(stderr,2)
 2        format(1x,'TANGVP: Tresca''s criterion not implemented')
call stoper()
 else if (NPTYPE.eq.3) then
write(stderr,3)
 3        format(1x,'TANGVP: Mohr-Coulomb criterion not implemented')
call stoper()
 else if (NPTYPE.eq.4) then
write(stderr,4)
 4       format(1x,'TANGVP: Drucker-Prager criterion not',1x, &
     'implemented')
call stoper()
 endif
!
 FNORM = (YIELD-COHESN)/COHESN
 if (FNORM.lt.ALLOW) return
!
 xJ2 = AJ2(STN)
 if (xJ2.le.dfmin) return
!
!    Compute PHI and dPHI/dF (GRADPHI)
 if (NFLOW.eq.1) then
if (ABS(YLDPWR-ONE).lt.EPS) then
PHI=FNORM
GRADPHI=FNORM/COHESN
else
PHI=FNORM**YLDPWR
GRADPHI=YLDPWR/COHESN*FNORM**(YLDPWR-ONE)
endif
 else
PHI = EXP(FNORM*YLDPWR)-ONE
GRADPHI = YLDPWR/COHESN*EXP(FNORM*YLDPWR)
 endif
!
!    Calculate effective stress
 STNEFF = SQRT(xJ2)
!
!    Matrix multipliers
 FACT1 = FLUIDY*root3*PHI/(TWO*STNEFF)
 FACT2 = THREEQ*GRADPHI/xJ2 - root3*PHI/(FOUR*STNEFF**3)
 FACT2 = FACT2*FLUIDY
!
!    Compute M2-matrix
 call DEVSTN (STN,SDEV)
 AM2(1,1) = SDEV(1)*SDEV(1)
 AM2(1,2) = SDEV(1)*SDEV(2)
 AM2(1,3) = SDEV(1)*SDEV(3)
 AM2(1,4) = TWO*SDEV(1)*SDEV(4)
 AM2(1,5) = TWO*SDEV(1)*SDEV(5)
 AM2(1,6) = TWO*SDEV(1)*SDEV(6)
 AM2(2,2) = SDEV(2)*SDEV(2)
 AM2(2,3) = SDEV(2)*SDEV(3)
 AM2(2,4) = TWO*SDEV(2)*SDEV(4)
 AM2(2,5) = TWO*SDEV(2)*SDEV(5)
 AM2(2,6) = TWO*SDEV(2)*SDEV(6)
 AM2(3,3) = SDEV(3)*SDEV(3)
 AM2(3,4) = TWO*SDEV(3)*SDEV(4)
 AM2(3,5) = TWO*SDEV(3)*SDEV(5)
 AM2(3,6) = TWO*SDEV(3)*SDEV(6)
 AM2(4,4) = FOUR*SDEV(4)*SDEV(4)
 AM2(4,5) = FOUR*SDEV(4)*SDEV(5)
 AM2(4,6) = FOUR*SDEV(4)*SDEV(6)
 AM2(5,5) = FOUR*SDEV(5)*SDEV(5)
 AM2(5,6) = FOUR*SDEV(5)*SDEV(6)
 AM2(6,6) = FOUR*SDEV(6)*SDEV(6)
 call SYMMET (AM2,NSTR)
!
 do 100 j=1,NSTR
BETAJP(1,j) = FACT1*AM1(1,j)+FACT2*AM2(1,j)
BETAJP(2,j) = FACT1*AM1(2,j)+FACT2*AM2(2,j)
BETAJP(3,j) = FACT1*AM1(3,j)+FACT2*AM2(3,j)
BETAJP(4,j) = FACT1*AM1(4,j)+FACT2*AM2(4,j)
BETAJP(5,j) = FACT1*AM1(5,j)+FACT2*AM2(5,j)
BETAJP(6,j) = FACT1*AM1(6,j)+FACT2*AM2(6,j)
100    continue
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine INYLD (STN,PRPLAS,YIELD,BETB,FNORM)
 USE CONSTANTS, only: half, one, hree, sixe, root3, deg2rad, dfmin,dlogmin,dfmax,dlogmax

!
! determines FNORM
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=6)
!-pass
 dimension STN(NSTR),PRPLAS(9),BETB(NSTR)

!-locl
 external AJ2
!-init
 ALLOW  = 0.01d0
 NPTYPE = INT(PRPLAS(1))
 FRANGL = PRPLAS(2)*deg2rad
 COHESN = PRPLAS(3)
 STNHRD = PRPLAS(7)
 FNORM  = -ONE
!
 DEBAR = half*AJ2(BETB)
 if (DEBAR.gt.dfmin) DEBAR = SQRT(DEBAR)
!
 STNEFF = AJ2(STN)
 if (STNEFF.lt.dfmin) return
 STNEFF = SQRT(STNEFF)
!
 if (NPTYPE.eq.1) then
!        Tresca
continue
 else if (NPTYPE.eq.2) then
!        Von Mises
continue
 else if (NPTYPE.eq.3) then
!        Mohr-Coulomb
COHESN = COHESN*COS(FRANGL)
 else if (NPTYPE.eq.4) then
!        Drucker-Prager
COHESN = SIX*COHESN*COS(FRANGL)/(root3*(THREE-SIN(FRANGL)))
 endif
 if (STNHRD.gt.dfmin) COHESN = COHESN+STNHRD*DEBAR
 if (COHESN.lt.0.001d0) COHESN = ONE
!
 FNORM = (YIELD-COHESN)/COHESN
 if (FNORM.lt.ALLOW) FNORM = -1.0d0
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine FLOWVP (AVECT,STN,PRPLAS,YIELD,BETB)

   use constants, only: eps


!
! Compute the viscoplastic strain rate matrix BETB
! (Ref.: Owen & Hinton, sec. 8.9)
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=6)
 parameter (ONE=1.d0)
!-pass
 dimension AVECT(NSTR),STN(NSTR),PRPLAS(9),BETB(NSTR)
!-init
 FLUIDY = PRPLAS(4)
 YLDPWR = PRPLAS(5)
 NFLOW  = INT(PRPLAS(6))
!
!    calculate FNORM, the normalized stress overshoot (uses BETB)
 call INYLD (STN,PRPLAS,YIELD,BETB,FNORM)
 call CLEAR (BETB,NSTR,"BETB")
!
 if (FNORM.le.0d0) return
!
 if (NFLOW.eq.1) then
if (ABS(YLDPWR-ONE).le.EPS) then
PHI = FLUIDY*FNORM
else
PHI = FLUIDY*FNORM**YLDPWR
endif
 else
PHI = FLUIDY*(EXP(YLDPWR*FNORM)-ONE)
 endif
!
 BETB(1) = PHI*AVECT(1)
 BETB(2) = PHI*AVECT(2)
 BETB(3) = PHI*AVECT(3)
 BETB(4) = PHI*AVECT(4)
 BETB(5) = PHI*AVECT(5)
 BETB(6) = PHI*AVECT(6)
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine WINKLP (BTOT,DELD,IWELM,IWSIDE,IWTIME,WPRES,IEN, &
   DUM1,X,D,DUM2,LM,LMF,TFAULT,SKEW, &
   DXE,LMX,NSLIP,DX)
   USE MATERIALSMODULE, only : LMAT
   USE MODELDEFINITION
   use constants, only: eps
use timestepmodule, only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

!
! Subroutine to compute Winkler restoring pressures from displacements
! and add them to the global force vector BTOT

! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NSD=3,NDOF=3,NEN=4)
 parameter (THIRD=1d0/3d0)
!-pass
 dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
   DUM1(*),X(NSD,*),DUM2(6,*),LM(NDOF,NEN,*),BTOT(*), &
   DELD(NDOF,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN,*),LMX(NDOF,NEN,*),NSLIP(5,*), &
 DX(NDOF,*)
!-locl
 dimension w(NDOF),xl(NSD,NEN),dl(NDOF,NEN),lien(3),xn(NSD), &
   tmp(NSD),rot(3,3),angle(2)
!
 if (NUMWNK.eq.0) return
!
 do 400 nn=1,IABS(NUMWNK)
n = IABS(IWELM(nn))
LL= IWTIME(nn)
    if (LL.lt.0) then
    if (LL.ne.-1 .and. NSTEP.ge.-LL) goto 400
    elseif (LL.eq.0) then
        goto 400
    else
    if (NSTEP.lt.(LL-1)) goto 400
    endif
call LCOORD (X,XL,IEN(1,n))
if (LGDEF.ne.0) then
call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
call AddFaultDisplacement (DL,LMF(1,1,n),TFAULT,NDOF,NEN)
    call ADDSNE (DL,DXE(1,1,n),NDOF,NEN)
call REZONE (XL,DL,1d0)
endif
    call SIDENP (IABS(IWSIDE(nn)),.true.,lien)
i = lien(1)
j = lien(2)
k = lien(3)
!           Calculate surface area
    call FCGEOM (XL,lien,area,xn)
    if (NUMROT.ne.0) then
!               Rotate face normal into local direction
        np = IEN(i,n)
        angle(1) = -SKEW(1,np)
        angle(2) = -SKEW(2,np)
        np = IEN(j,n)
        angle(1) = angle(1)-SKEW(1,np)
        angle(2) = angle(2)-SKEW(2,np)
        np = IEN(k,n)
        angle(1) = angle(1)-SKEW(1,np)
        angle(2) = angle(2)-SKEW(2,np)
angle(1) = angle(1)*THIRD
angle(2) = angle(2)*THIRD
        if (ABS(angle(1)).gt.EPS.or.ABS(angle(2)).gt.EPS) then
            call FORMRT (angle,rot,NSD)
            call EQUATE (tmp,xn,NSD)
            call VMPRD (rot,tmp,xn,NSD,3)
        endif
    endif
w(1) = WPRES(nn) * ABS( area * xn(1) ) * THIRD
w(2) = WPRES(nn) * ABS( area * xn(2) ) * THIRD
w(3) = WPRES(nn) * ABS( area * xn(3) ) * THIRD
    if (NUMSLP.gt.0) then
        if (IWELM(nn).gt.0) then
            ni = 0
            nj = 0
            nk = 0
            do ii = 1,NUMSLP
                if (NSLIP(2,ii).eq.IEN(i,n)) ni = ni + 1
                if (NSLIP(2,ii).eq.IEN(j,n)) nj = nj + 1
                if (NSLIP(2,ii).eq.IEN(k,n)) nk = nk + 1
            enddo
        else
            ni = 1
            nj = 1
            nk = 1
        endif
    endif
do 300 idof=1,NDOF
if (ABS(w(idof)).lt.EPS) goto 300
ii = LM(idof,i,n)
      if (ii.gt.0) BTOT(ii)=BTOT(ii)- &
       w(idof)*DELD(idof,IEN(i,n))
ii = LM(idof,j,n)
      if (ii.gt.0) BTOT(ii)=BTOT(ii)- &
       w(idof)*DELD(idof,IEN(j,n))
ii = LM(idof,k,n)
      if (ii.gt.0) BTOT(ii)=BTOT(ii)- &
       w(idof)*DELD(idof,IEN(k,n))
!
if (NUMSLP.eq.0) goto 300
ii = IABS(LMX(idof,i,n))
      if (ii.gt.0) BTOT(ii) = BTOT(ii) - &
       w(idof)*DX(idof,IEN(i,n))/DBLE(ni)
ii = IABS(LMX(idof,j,n))
      if (ii.gt.0) BTOT(ii) = BTOT(ii) - &
       w(idof)*DX(idof,IEN(j,n))/DBLE(nj)
ii = IABS(LMX(idof,k,n))
      if (ii.gt.0) BTOT(ii) = BTOT(ii) - &
       w(idof)*DX(idof,IEN(k,n))/DBLE(nk)
300        continue
400    continue
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine SHAREL (IEN,IJSUR,NUMEL,NEL)
!
!     Finds elements neighbouring to element NEL.
!    Neighbouring elements are identified by three nodes.
!
 parameter (NSUR=4,NEN=4)
!-pass
 integer NEL,NUMEL
 dimension IEN(NEN,*),IJSUR(NSUR,*)

!-locl
 integer NINDX
 external NINDX
!
!    set to last IJSUR(isur,NEL) <> 0 (last index that was written)
 isur = NINDX (IJSUR,NEL)
!    loop other elements
 do 200 iel=NEL+1,NUMEL
!        compare nodal points of element NEL with those of IEL > NEL
nshare = 0
do 100 n=1,NEN
if (IEN(n,NEL).eq.IEN(1,iel)) nshare = nshare+1
if (IEN(n,NEL).eq.IEN(2,iel)) nshare = nshare+1
if (IEN(n,NEL).eq.IEN(3,iel)) nshare = nshare+1
if (IEN(n,NEL).eq.IEN(4,iel)) nshare = nshare+1
100        continue
if (nshare.gt.3) then
write(stderr,1) NEL,iel,nshare
 1           format(1x,'sharel: elements ',I5,' and ',I5,' share ', &
         I1,' nodal points')
call xit(1," ")
else if (nshare.eq.3) then
isur = isur + 1
if (isur.le.NSUR) then
    IJSUR(isur,NEL) = iel
!            copy info NEL into IJSUR(*,iel)
!            find free index
    ii = NINDX (IJSUR,iel) + 1
    IJSUR(ii,iel) = NEL
else
    write(stderr,2) NEL,NSUR
 2               format(1x,'sharel: element ',I5, &
             ' surrounded by more than ',I1,' elements')
    call xit(1," ")
endif
if (isur.eq.NSUR) goto 300
endif
200    continue
!
300    return
 end
!-------------------------------------------------------------------------------
 function NINDX (IJSUR,NEL)
!
 parameter (NSUR=4)
!-pass
 integer NEL
 dimension IJSUR(NSUR,*)

!
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
 2        format(1x,8(I5,1x))
call xit(1," ")
 endif
 NINDX = i
 return
 end
!-------------------------------------------------------------------------------
! RG not adapted
!
 subroutine FLUID (BETB,PRPMAT,PRPLAS)
   USE TIMESTEPMODULE
 USE CONSTANTS, only: half, one, two, three, four, six, root3, deg2rad, dfmin,dlogmin,dfmax,dlogmax
!
!    Routine to determine visco-plastic fluidity from
!    explicit time stepping limit (Cormeau) (LINEAR only)
!     (Ref.: Owen & Hinton, sec. 8.9)
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=6)
!-pass
 dimension PRPMAT(6),PRPLAS(9),BETB(NSTR)

!-init
 NPTYPE = INT(PRPLAS(1))
 FRANGL = PRPLAS(2) * deg2rad
 COHESN = PRPLAS(3)
 STNHRD = PRPLAS(7)
 FLUIDY = PRPLAS(9)
 E      = PRPMAT(1)
 POIS   = PRPMAT(2)

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
COHESN = COHESN * COS(FRANGL)
 else if (NPTYPE.eq.4) then
!        Drucker-Prager
COHESN = SIX * COHESN * COS(FRANGL) / &
       ( root3 * (THREE - SIN(FRANGL)) )
 else
write(stderr,1) NPTYPE
 1        format(1x,'FLUID: unknown plasticity type ',I5)
call stoper()
 endif
 if (STNHRD.gt.dfmin) COHESN = COHESN + STNHRD*DEBAR

 if (NPTYPE.eq.1) then
!        Tresca
FLUIDY = (ONE+POIS) * COHESN / (E * DELTP)
 else if (NPTYPE.eq.2) then
!        Von Mises
FLUIDY = FOUR*(ONE+POIS)*COHESN/(THREE*E*DELTP)
 else if (NPTYPE.eq.3) then
!        Mohr-Coulomb
FRANGL = PRPLAS(2)
SINPHI = SIN(FRANGL*deg2rad)
SINPHI = SINPHI*SINPHI
FLUIDY = FOUR*(ONE+POIS)*(ONE-TWO*POIS)*COHESN / &
     ((ONE-TWO*POIS+SINPHI)*E*DELTP)
 else if (NPTYPE.eq.4) then
!        Drucker-Prager
continue
 else
write(stderr,1) NPTYPE
call stoper()
 endif
 if (FLUIDY.gt.dfmin) PRPLAS(4) = FLUIDY

 return
 end
!-------------------------------------------------------------------------------
 subroutine BPMATRIX (XL,DL,QUAD,BP,IERR)
!
! Sets up derivatives matrix BP
!
 implicit none
 integer NDOF,NSD,NEN
 parameter (NDOF=3, NSD=3, NEN=4)
 double precision ZERO
 parameter (ZERO=0d0)
!-pass
 logical QUAD
 integer IERR
 double precision XL,DL,BP
 dimension XL(NSD,NEN),DL(NDOF,NEN),BP(NDOF,NSD)
!-locl
 double precision det,sh
 dimension sh(4,NEN)
!
!    Compute spatial derivatives of shape functions
 call SHPTET (ZERO,ZERO,ZERO,XL,det,sh,IERR,1)
 if (IERR.ne.0) return
!
 BP(1,1) = SH(1,1)*DL(1,1)+SH(1,2)*DL(1,2)+SH(1,3)*DL(1,3)+ &
   SH(1,4)*DL(1,4)
 BP(2,1) = SH(1,1)*DL(2,1)+SH(1,2)*DL(2,2)+SH(1,3)*DL(2,3)+ &
   SH(1,4)*DL(2,4)
 BP(3,1) = SH(1,1)*DL(3,1)+SH(1,2)*DL(3,2)+SH(1,3)*DL(3,3)+ &
   SH(1,4)*DL(3,4)
 BP(1,2) = SH(2,1)*DL(1,1)+SH(2,2)*DL(1,2)+SH(2,3)*DL(1,3)+ &
   SH(2,4)*DL(1,4)
 BP(2,2) = SH(2,1)*DL(2,1)+SH(2,2)*DL(2,2)+SH(2,3)*DL(2,3)+ &
   SH(2,4)*DL(2,4)
 BP(3,2) = SH(2,1)*DL(3,1)+SH(2,2)*DL(3,2)+SH(2,3)*DL(3,3)+ &
   SH(2,4)*DL(3,4)
 BP(1,3) = SH(3,1)*DL(1,1)+SH(3,2)*DL(1,2)+SH(3,3)*DL(1,3)+ &
   SH(3,4)*DL(1,4)
 BP(2,3) = SH(3,1)*DL(2,1)+SH(3,2)*DL(2,2)+SH(3,3)*DL(2,3)+ &
   SH(3,4)*DL(2,4)
 BP(3,3) = SH(3,1)*DL(3,1)+SH(3,2)*DL(3,2)+SH(3,3)*DL(3,3)+ &
   SH(3,4)*DL(3,4)
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine sdimnt (ISURF,SURF,B,ID,GRAV,PRPSED,DELTP)
   USE MODELDEFINITION
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSD=3,NDOF=3)
!-pass
 dimension ISURF(*),SURF(NSD,*),B(*),ID(NDOF,*),GRAV(*), &
   PRPSED(2)

!-locl
 if (NSURF.le.0 .or. NSED.eq.0) return
!
 write(stderr,10)
10    format(1x,'"SDIMNT" NOT IMPLEMENTED IN 3-D VERSION')
!
 return
 end
!-----------------------------------------------------------------------
 function SRFHGT ()
 implicit double precision (a-z)
 SRFHGT = 0d0
 return
 end
!-----------------------------------------------------------------------
subroutine SLSKEW (IDSLP,SKEW,NSELSD,IEN,X,D,TFAULT,DXE,LMF)
USE CONSTANTS, only: rad2deg, zero, two
!
! Updates skew angles for slippery nodes
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (A-H,O-Z)
parameter (NSD=3, NEN=4, NDOF=3)
!-pass
dimension IDSLP(*),SKEW(2,*),NSELSD(2,*),IEN(NEN,*),X(NSD,*), &
 D(NDOF,*),TFAULT(NDOF,*),DXE(NDOF,NEN,*),LMF(NDOF,NEN,*)
!
 return
 end
!-----------------------------------------------------------------------
subroutine LCOORD (X,XL,IEN)
USE MODELDEFINITION
!
! Program to localize coordinate data
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (A-H,O-Z)
parameter (NSD=3)
!-pass
dimension IEN(*),X(NSD,*),XL(NSD,*)
!
k = IEN(1)
XL(1,1) = X(1,k)
XL(2,1) = X(2,k)
XL(3,1) = X(3,k)
k = IEN(2)
XL(1,2) = X(1,k)
XL(2,2) = X(2,k)
XL(3,2) = X(3,k)
k = IEN(3)
XL(1,3) = X(1,k)
XL(2,3) = X(2,k)
XL(3,3) = X(3,k)
k = IEN(4)
XL(1,4) = X(1,k)
XL(2,4) = X(2,k)
XL(3,4) = X(3,k)
!
return
end
!-------------------------------------------------------------------------------
function dilate (S)
!
! computes volume change or volume change rate, depending on what S is.
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (a-h,o-z)
!-pass
dimension S(*)
!-locl
 parameter (THIRD=1D0/3D0)
!-init
!
 dilate = AJ1(S)*THIRD
!
return
end
!-------------------------------------------------------------------------------
function EFFEE (S)
use constants, only: dfmin

!
! Compute effective deviatoric strain
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (a-h,o-z)
!-pass
dimension S(*)
!-locl
 double precision J2
!-init
 effee = ZERO
!
 J2 = AJ2(S)
 if (J2.gt.dfmin) effee = SQRT(J2)
!
return
end
!-------------------------------------------------------------------------------
subroutine STRDEV (S)

use constants, only: third
!
! computes strain (rate) deviator
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (a-h,o-z)
!-pass
dimension S(*)
!
 dil = AJ1(S)*THIRD
 S(1) = S(1)-dil
 S(2) = S(2)-dil
 S(3) = S(3)-dil
!
return
end
!-------------------------------------------------------------------------------
function EFFSTN (S)

use constants, only: dfmin
!
! computes effective deviatoric stress. S is the full stress tensor.
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (a-h,o-z)
!-pass
dimension S(*)
!-locl
 double precision J2
!
 J2 = AJ2(S)
 if (J2.gt.dfmin) EFFSTN = SQRT(J2)
!
return
end
!-------------------------------------------------------------------------------
subroutine STNTEN (STN,SL)
USE MODELDEFINITION
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (a-h,o-z)
 parameter (NSTR=6)
!-pass
dimension STN(NSTR),SL(3,3)

 SL(1,1) = STN(1)
 SL(2,2) = STN(2)
 SL(3,3) = STN(3)
 SL(1,2) = STN(4)
 SL(2,1) = STN(4)
 SL(1,3) = STN(5)
 SL(3,1) = STN(5)
 SL(2,3) = STN(6)
 SL(3,2) = STN(6)
!
return
end
!-------------------------------------------------------------------------------
subroutine STNDEV (S)
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (a-h,o-z)
!-pass
dimension S(*)
!-locl
 parameter (THIRD=1D0/3D0)
!
 p = AJ1(S)*THIRD
 S(1) = S(1) - p
 S(2) = S(2) - p
 S(3) = S(3) - p
!
return
end
!-------------------------------------------------------------------------------
function PRESSURE (S,POIS)
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (a-h,o-z)
!-pass
dimension S(*)
!-locl
 parameter (THIRD=1D0/3D0)
!
 pressure = AJ1(S)*THIRD
!
return
end
