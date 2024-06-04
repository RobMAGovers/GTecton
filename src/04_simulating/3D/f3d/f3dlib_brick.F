!       Version 2005.1
!
!       These routines are the rheology-dependent portion of GTECTON
!    The routines here are for 3-dimensional displacements
!-------------------------------------------------------------------------------
!
! DEFINITIONS
!
! Nodal points:    8             7         Faces:
!                  o-------------o                 o-------------o
!                 /|            /|                /|            /|
!                / |           / |               / |           / |
!               /  |          /  |              /  |   6      /  |
!              /   |         /   |             /   |      3  /   |
!            5/    |       6/    |            /    |        /    |
!            o-------------o     |           o-------------o     |
!            |     |       |     |           |  4  |       |  2  |
!            |     o-------|-----o           |     o-------|-----o
!            |    /4       |    /3           |    /        |    /
!            |   /         |   /             |   /  1      |   /
!            |  /          |  /              |  /      5   |  /
!            | /           | /               | /           | /
!            |/            |/                |/            |/
!            o-------------o                 o-------------o
!            1             2
!
!
! Collapse of nodal points 4 onto 3 AND 8 onto 7 is allowed; this will give
! a prism-element.
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
!    AJ1        FIRST TENSOR INVARIANT
!    AJ2        SECOND TENSOR INVARIANT
!    BDELD        STRAIN FROM DISPLACEMENTS
!    BMATRIX        STRAIN-DISPLACEMENT MATRIX
!    CPSTIFF        LOCAL HEAT CAPACITY STIFFNESS MATRIX
!    DEVSTN        STRESS DEVIATOR CALCULATION
!    EFORCE        EFFECTIVE LOAD VECTOR FOR VISCOUS FLOW
!    ELCRD        CALCULATES COORDINATES OF DEFORMED ELEMENT
!    ELCTR        CALCULATES CENTER COORDINATES OF DEFORMED ELEMENT
!    EIGMAX        CALCULATES MAXIMUM EIGENVALUE
!    FCGEOM        PERFORMS CALCULATIONS ON FACES OF ELEMENT
!    FLOWVP        VISCOPLASTIC STRAIN RATE
!    FORMBT        VISCOUS FLOW CORRECTION FOR STRAIN
!    FORMJB        JACOBEAN MATRIX FOR VISCOUS FLOW
!    GRAVLD        GRAVITATIONAL BODY FORCE COMPUTATION
!    INTCON        DIMENSIONAL AND TOPOLOGIC DATA
!    YIELDS        COMPUTES PLASTIC YIELD STRESS
!    LHEAT        HEAT PRODUCTION CONTRIBUTION TO LOAD VECTOR
!    LFLUX        CALCULATES ELEMENT HEAT FLUX
!    FLUID        DETERMINES VISCOPLASTIC FLUIDITY
!    MATERL        MATERIAL MATRICES
!    REZONE        LOCAL REZONE FOR LARGE DEFORMATION
!    SDIMNT        SURFACE ELEVATION AND LOADS FROM SEDIMENTARY TRANSPORT
!    SHAP30        SHAPE FUNCTION FOR 3-D ELEMENTS
!    SHAREL        IDENTIFIES NEIGHBOURING ELEMENTS
!       SIDENP          RETURNS LOCAL NODE NUMBERS ON A SPECIFIC SIDE
!    STIFF        STIFFNESS MATRIX
!    TANGVP        VISCOPLASTIC JACOBIAN MATRIX
!    TCENTR        COMPUTES ELEMENT CENTER TEMPERATURES
!    THGRAD        THERMAL GRADIENTS
!    TSTIFF        LOCAL CONDUCTION STIFFNESS MATRIX
!    WINKLP        WINKLER PRESSURES CONTRIBUTION TO GLOBAL FORCE VECTOR
!    YIELDF        CONSTRUCTS PLASTIC FLOW VECTOR
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
 integer NTYPE,NSD,NDOF,NSTR,NEN


 NTYPE = 5
 NSD   = 3
 NDOF  = 3
 NSTR  = 6
 NEN   = 8
 NSUR  = 6
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
 end
!-------------------------------------------------------------------------------
 SUBROUTINE SHAP30 (R,S,T,X,DET,SH,XS,QUAD,IERR,IFLAG)
!
! Program to compute shape functions for quadrilateral or triangle
!
!                = 1, calculate r-, s-, t- derivatives
!                  and shape functions only
!                = 2, calculate r-, s-, t- derivatives
!                  and shape functions and Jacobian
! IFLAG                = 3, full calculation &
! R,S T,            = natural coordinates
! X(NSD,NEN)            = global coordinates
! QUAD                = TRUE if brick element
!                = FALSE if prism element
! SH(1,NEN)            = x- or r- derivatives of shape functions
! SH(2,NEN)            = y- or s- derivatives of shape functions
! SH(3,NEN)            = z- or t- derivatives of shape functions
! SH(4,NEN)            = shape functions
! XS(NSD,NSD)            = inverse jacobian matrix
! DET                = Jacobian matrix determinant
! IERR                = 0, succesful calculation
!                = 1, Jacobian determinant error
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NEN=8, NSD=3)
 parameter (ZERO=0d0, ONE=1d0, F=.125d0)
!-pass
 logical QUAD
 integer ierr
 dimension SH(4,NEN),X(NSD,NEN),XS(NSD,NSD)
!-locl
 dimension RA(NEN),SA(NEN),TA(NEN)
 save RA,SA,TA

!-init
 data RA/-1., 1., 1.,-1.,-1., 1., 1.,-1./
 data SA/-1.,-1., 1., 1.,-1.,-1., 1., 1./
 data TA/-1.,-1.,-1.,-1., 1., 1., 1., 1./
 ierr = 0
!
!    Calculate shape functions and r-, s- and t-derivatives of shape functions
 do 10 i=1,NEN
RX = (ONE+RA(i)*R)
SX = (ONE+SA(i)*S)
TX = (ONE+TA(i)*T)
SH(4,i) = F*RX*SX*TX
SH(1,i) = F*RA(i)*SX*TX
SH(2,i) = F*SA(i)*RX*TX
SH(3,i) = F*TA(i)*RX*SX
10    continue
 if (.not.QUAD) then
!        Collapse of brick element on prism.
do 20 i=1,4
SH(i,3) = SH(i,3)+SH(i,4)
SH(i,7) = SH(i,7)+SH(i,8)
SH(i,4) = ZERO
SH(i,8) = ZERO
20        continue
 endif
!    Jacobian matrix calculation
 do 40 i=1,NSD
do 30 j=1,NSD
XS(i,j) = SH(i,1)*X(j,1)+SH(i,2)*X(j,2)+SH(i,3)*X(j,3) &
       +SH(i,4)*X(j,4)+SH(i,5)*X(j,5)+SH(i,6)*X(j,6) &
       +SH(i,7)*X(j,7)+SH(i,8)*X(j,8)
30        continue
40    continue
 DET = XS(1,1)*XS(2,2)*XS(3,3)-XS(1,1)*XS(2,3)*XS(3,2) &
       -XS(2,1)*XS(1,2)*XS(3,3)+XS(2,1)*XS(1,3)*XS(3,2) &
       +XS(3,1)*XS(1,2)*XS(2,3)-XS(3,1)*XS(1,3)*XS(2,2)
 if (DET.le.zero) go to 1000
!
 if (IFLAG.eq.1) return
!
!    Calculatate inverse Jacobian determinant and store it back in XS
 RDET=ONE/DET
 Xi11=(XS(2,2)*XS(3,3)-XS(3,2)*XS(2,3))*RDET
 Xi21=(XS(2,3)*XS(3,1)-XS(2,1)*XS(3,3))*RDET
 Xi31=(XS(2,1)*XS(3,2)-XS(2,2)*XS(3,1))*RDET
 Xi12=(XS(1,3)*XS(3,2)-XS(1,2)*XS(3,3))*RDET
 Xi22=(XS(1,1)*XS(3,3)-XS(1,3)*XS(3,1))*RDET
 Xi32=(XS(1,2)*XS(3,1)-XS(1,1)*XS(3,2))*RDET
 Xi13=(XS(1,2)*XS(2,3)-XS(1,3)*XS(2,2))*RDET
 Xi23=(XS(1,3)*XS(2,1)-XS(1,1)*XS(2,3))*RDET
 Xi33=(XS(1,1)*XS(2,2)-XS(1,2)*XS(2,1))*RDET
 XS(1,1)=Xi11
 XS(2,1)=Xi21
 XS(3,1)=Xi31
 XS(1,2)=Xi12
 XS(2,2)=Xi22
 XS(3,2)=Xi32
 XS(1,3)=Xi13
 XS(2,3)=Xi23
 XS(3,3)=Xi33
!
 if (IFLAG.eq.2) return
!
!    Calculate x-, y- and z-derivatives of shape functions
 do 50 i=1,NEN
dSHdR = SH(1,i)
dSHdS = SH(2,i)
dSHdT = SH(3,i)
SH(1,i) = XS(1,1)*dSHdR+XS(1,2)*dSHdS+XS(1,3)*dSHdT
SH(2,i) = XS(2,1)*dSHdR+XS(2,2)*dSHdS+XS(2,3)*dSHdT
SH(3,i) = XS(3,1)*dSHdR+XS(3,2)*dSHdS+XS(3,3)*dSHdT
50    continue
!
 return
!
 1000    write(stderr,2000) DET
 2000    format(///1x,'Shape function fails! Determinant is ',1PE20.4)
 ierr = 1
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
 parameter (NEN=8,NSD=3)
!-pass
 integer IEN
 double precision X
 dimension X(NSD,*),IEN(NEN,*)

!-locl
 logical quad
 integer i,j,n
 double precision ra,sa,ta,xl,sh,xs,det
 dimension xl(NSD,NEN),ra(NEN),sa(NEN),ta(NEN),sh(3,NEN), &
   xs(NSD,NSD)
 save ra,sa,ta
!-init
 data ra/-0.125, 0.125, 0.125,-0.125,-0.125, 0.125, 0.125,-0.125/
 data sa/-0.125,-0.125, 0.125, 0.125,-0.125,-0.125, 0.125, 0.125/
 data ta/-0.125,-0.125,-0.125,-0.125, 0.125, 0.125, 0.125, 0.125/
!
 if (NUMEL.lt.1) then
write(stderr,10)
10        format(1x,'NUMEL < 1')
   call exitp(1)
 endif
!
 do n=1,NUMEL
do j=1,NEN
    i = IEN(j,n)
    if (i.lt.1 .or. i.gt.NUMNP) then
    write(stderr,20) n,i
20               format(1x,'element ',I6,': reference to node ', &
             'number ',I6)
    call exitp(1)
endif
xl(1,j) = X(1,i)
xl(2,j) = X(2,i)
xl(3,j) = X(3,i)
enddo
do j=1,NEN
sh(1,j) = ra(j)
sh(2,j) = sa(j)
sh(3,j) = ta(j)
enddo
if (IEN(3,n).eq.IEN(4,n)) then
    if (IEN(7,n).eq.IEN(8,n)) then
    quad = .false.
else
    write(stderr,22) n,3,4,7,8
22               format(1x,'element ',I6,': nodes ',I1,' and ',I1, &
             ' are identical, nodes ',I1,' and ',I1,' are not')
    call exitp(1)
endif
else
    if (IEN(7,n).eq.IEN(8,n)) then
    write(stderr,22) n,7,8,3,4
    call exitp(1)
else
    quad = .true.
endif
endif
if (.not.quad) then
    do j=1,3
    sh(j,3) = sh(j,3)+sh(j,4)
    sh(j,7) = sh(j,7)+sh(j,8)
    sh(j,4) = 0d0
    sh(j,8) = 0d0
enddo
endif
do i=1,NSD
do j=1,NSD
    xs(i,j) = sh(i,1)*xl(j,1)+sh(i,2)*xl(j,2)+ &
           sh(i,3)*xl(j,3)+sh(i,4)*xl(j,4)+sh(i,5)*xl(j,5)+ &
           sh(i,6)*xl(j,6)+sh(i,7)*xl(j,7)+sh(i,8)*xl(j,8)
enddo
enddo
det = xs(1,1)*xs(2,2)*xs(3,3)-xs(1,1)*xs(2,3)*xs(3,2) &
       -xs(2,1)*xs(1,2)*xs(3,3)+xs(2,1)*xs(1,3)*xs(3,2) &
       +xs(3,1)*xs(1,2)*xs(2,3)-xs(3,1)*xs(1,3)*xs(2,2)
if (ABS(det).le.EPS) then
write(stderr,30) i
30        format(1x,'degenerate element ',I6)
call exitp(1)
endif
if (det.lt.0d0) then
write(stderr,40) i
40           format(1x,'element ',I6,' is negatively oriented (node', &
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
 subroutine MATERL (DMAT,E,POIS)
 USE CONSTANTS, only: one, two, half
!
! Constructs the material matrix relating stress and strain.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NSTR=6)
!-pass
 dimension DMAT(NSTR,NSTR)
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
 subroutine BDELD (XL,DL,EE,QUAD,ierr)
 USE CONSTANTS, only: half
!
! Subroutine to compute strains from displacements in each element
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NDOF=3, NSD=3, NSTR=6, NEN=8,NEE=24)
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
 subroutine BMATRIX (B,SH)
!
! Defines the strain-displacement matrix B
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NEN=8, NSTR=6, NDOF=3)
 parameter (NEE=NDOF*NEN)
 parameter (ZERO=0d0)
!-pass
 dimension B(NSTR,NEE),SH(4,NEN)
!
 k = 0
 do 100 i=1,NEN
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
100    continue
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine GRAVLD (P,XL,GRAV,DUM1,WT,QUAD,ierr)
 USE CONSTANTS, only: G
!
! Computes the contribution to the load from gravity body forces
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NINT=8, NEN=8, NSD=3, NDOF=3)
 parameter (NEE=NEN*NDOF)
!-pass
 integer ierr
 logical QUAD
 dimension P(NEE),XL(NSD,NEN),GRAV(NDOF)
!-locl
 dimension RG(NEN),SG(NEN),TG(NEN),XS(NSD,NSD),SH(4,NEN)
save RG,TG,SG
!-init
 data RG/-1., 1., 1.,-1.,-1., 1., 1.,-1./
 data SG/-1.,-1., 1., 1.,-1.,-1., 1., 1./
 data TG/-1.,-1.,-1.,-1., 1., 1., 1., 1./
!
 do i=1,NINT
!        Compute shape functions
call SHAP30(RG(i)*G,SG(i)*G,TG(i)*G,XL,DET,SH,XS,QUAD, &
       ierr,1)
if (ierr.ne.0) return
GRAV1=GRAV(1)*DET*WT
GRAV2=GRAV(2)*DET*WT
GRAV3=GRAV(3)*DET*WT
P( 1)=P( 1)+GRAV1*SH(4,1)
P( 2)=P( 2)+GRAV2*SH(4,1)
P( 3)=P( 3)+GRAV3*SH(4,1)
P( 4)=P( 4)+GRAV1*SH(4,2)
P( 5)=P( 5)+GRAV2*SH(4,2)
P( 6)=P( 6)+GRAV3*SH(4,2)
P( 7)=P( 7)+GRAV1*SH(4,3)
P( 8)=P( 8)+GRAV2*SH(4,3)
P( 9)=P( 9)+GRAV3*SH(4,3)
P(10)=P(10)+GRAV1*SH(4,4)
P(11)=P(11)+GRAV2*SH(4,4)
P(12)=P(12)+GRAV3*SH(4,4)
P(13)=P(13)+GRAV1*SH(4,5)
P(14)=P(14)+GRAV2*SH(4,5)
P(15)=P(15)+GRAV3*SH(4,5)
P(16)=P(16)+GRAV1*SH(4,6)
P(17)=P(17)+GRAV2*SH(4,6)
P(18)=P(18)+GRAV3*SH(4,6)
P(19)=P(19)+GRAV1*SH(4,7)
P(20)=P(20)+GRAV2*SH(4,7)
P(21)=P(21)+GRAV3*SH(4,7)
P(22)=P(22)+GRAV1*SH(4,8)
P(23)=P(23)+GRAV2*SH(4,8)
P(24)=P(24)+GRAV3*SH(4,8)
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine STIFF (DMAT,XL,S,DUM1,QUAD,ierr)
 USE MODELDEFINITION
 USE ALGEBRA, only: IDFORM
 USE CONSTANTS, only: G
!
! Computes the super-diagonal part of the local stiffness matrix.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NINT=8, NEN=8, NSD=3, NDOF=3, NSTR=6, NMODES=3)
 parameter (NEE=NEN*NDOF)
 parameter (NEX=NMODES*NDOF)
!-pass
 integer ierr
 logical QUAD
 dimension XL(NSD,NEN),DMAT(NSTR,NSTR),S(NEE,NEE)
!-locl
 dimension RG(NEN),SG(NEN),TG(NEN),XS(NSD,NSD),SH(4,NEN), &
   B(NSTR,NEE),SDAL(NEE,NEX),SAAL(NEX,NEX),SDA(NEE,NEX), &
   SAA(NEX,NEX),IDSAA(NEX),SA(NEE,NEE)
 logical firstentry
save RG,SG,TG,firstentry,IDSAA
!-init
 data firstentry /.true./
 data RG/-1., 1., 1.,-1.,-1., 1., 1.,-1./
data SG/-1.,-1., 1., 1.,-1.,-1., 1., 1./
data TG/-1.,-1.,-1.,-1., 1., 1., 1., 1./
!
if (INCOMP.ne.0 .and. QUAD) then
    if (firstentry) then
!               set up list of diagonal addresses of square array
        call IDFORM (IDSAA,NEX)
        firstentry = .false.
    endif
    call CLEAR (SDA,216,"SDA")
    call CLEAR (SAA,81,"SAA")
endif
!
 do l=1,NINT
!        calculate spatial derivatives of shape functions
call SHAP30 (RG(l)*G,SG(l)*G,TG(l)*G,XL,DET,SH,XS,QUAD, &
       ierr,3)
if (ierr.ne.0) return
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
    S(i,j)=S(i,j)+DET*(B(1,i)*DB1+B(2,i)*DB2+B(3,i)*DB3 &
           +B(4,i)*DB4+B(5,i)*DB5+B(6,i)*DB6)
enddo
enddo
!           Incompatible modes
if (INCOMP.ne.0.and.QUAD) then
call INCSTFF(RG(l)*G,SG(l)*G,TG(l)*G,DMAT,XL,B,SAAL, &
       SDAL,QUAD,ierr)
if (ierr.ne.0) return
do j=1,NEX
    do i=1,NEE
    SDA(i,j)=SDA(i,j)+DET*SDAL(i,j)
    enddo
    do i=1,NEX
    SAA(i,j)=SAA(i,j)+DET*SAAL(i,j)
    enddo
enddo
endif
 enddo
!
 if (INCOMP.ne.0.and.QUAD) then
!           Replace SAA by its inverse
    call INVMAT (SAA,IDSAA,NEX)
do j=1,NEE
!        column j of inv(SAA)*SAD
        SA1 = &
        SAA(1,1)*SDA(j,1)+SAA(1,2)*SDA(j,2)+SAA(1,3)*SDA(j,3) &
         +SAA(1,4)*SDA(j,4)+SAA(1,5)*SDA(j,5)+SAA(1,6)*SDA(j,6) &
         +SAA(1,7)*SDA(j,7)+SAA(1,8)*SDA(j,8)+SAA(1,9)*SDA(j,9)
        SA2 = &
        SAA(2,1)*SDA(j,1)+SAA(2,2)*SDA(j,2)+SAA(2,3)*SDA(j,3) &
         +SAA(2,4)*SDA(j,4)+SAA(2,5)*SDA(j,5)+SAA(2,6)*SDA(j,6) &
         +SAA(2,7)*SDA(j,7)+SAA(2,8)*SDA(j,8)+SAA(2,9)*SDA(j,9)
        SA3 = &
        SAA(3,1)*SDA(j,1)+SAA(3,2)*SDA(j,2)+SAA(3,3)*SDA(j,3) &
         +SAA(3,4)*SDA(j,4)+SAA(3,5)*SDA(j,5)+SAA(3,6)*SDA(j,6) &
         +SAA(3,7)*SDA(j,7)+SAA(3,8)*SDA(j,8)+SAA(3,9)*SDA(j,9)
        SA4 = &
        SAA(4,1)*SDA(j,1)+SAA(4,2)*SDA(j,2)+SAA(4,3)*SDA(j,3) &
         +SAA(4,4)*SDA(j,4)+SAA(4,5)*SDA(j,5)+SAA(4,6)*SDA(j,6) &
         +SAA(4,7)*SDA(j,7)+SAA(4,8)*SDA(j,8)+SAA(4,9)*SDA(j,9)
        SA5 = &
        SAA(5,1)*SDA(j,1)+SAA(5,2)*SDA(j,2)+SAA(5,3)*SDA(j,3) &
         +SAA(5,4)*SDA(j,4)+SAA(5,5)*SDA(j,5)+SAA(5,6)*SDA(j,6) &
         +SAA(5,7)*SDA(j,7)+SAA(5,8)*SDA(j,8)+SAA(5,9)*SDA(j,9)
        SA6 = &
        SAA(6,1)*SDA(j,1)+SAA(6,2)*SDA(j,2)+SAA(6,3)*SDA(j,3) &
         +SAA(6,4)*SDA(j,4)+SAA(6,5)*SDA(j,5)+SAA(6,6)*SDA(j,6) &
         +SAA(6,7)*SDA(j,7)+SAA(6,8)*SDA(j,8)+SAA(6,9)*SDA(j,9)
        SA7 = &
        SAA(7,1)*SDA(j,1)+SAA(7,2)*SDA(j,2)+SAA(7,3)*SDA(j,3) &
         +SAA(7,4)*SDA(j,4)+SAA(7,5)*SDA(j,5)+SAA(7,6)*SDA(j,6) &
         +SAA(7,7)*SDA(j,7)+SAA(7,8)*SDA(j,8)+SAA(7,9)*SDA(j,9)
        SA8 = &
        SAA(8,1)*SDA(j,1)+SAA(8,2)*SDA(j,2)+SAA(8,3)*SDA(j,3) &
         +SAA(8,4)*SDA(j,4)+SAA(8,5)*SDA(j,5)+SAA(8,6)*SDA(j,6) &
         +SAA(8,7)*SDA(j,7)+SAA(8,8)*SDA(j,8)+SAA(8,9)*SDA(j,9)
        SA9 = &
        SAA(9,1)*SDA(j,1)+SAA(9,2)*SDA(j,2)+SAA(9,3)*SDA(j,3) &
         +SAA(9,4)*SDA(j,4)+SAA(9,5)*SDA(j,5)+SAA(9,6)*SDA(j,6) &
         +SAA(9,7)*SDA(j,7)+SAA(9,8)*SDA(j,8)+SAA(9,9)*SDA(j,9)
!        SDA*inv(SAA)*SAD: matrix to subtract from S
do i=1,j
    SA(i,j) = &
            SDA(i,1)*SA1+SDA(i,2)*SA2+SDA(i,3)*SA3 &
           +SDA(i,4)*SA4+SDA(i,5)*SA5+SDA(i,6)*SA6 &
           +SDA(i,7)*SA7+SDA(i,8)*SA8+SDA(i,9)*SA9
enddo
!        adjust S
do i=1,j
    S(i,j) = S(i,j) - SA(i,j)
enddo
enddo
 endif
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine INCSTFF (RI,SI,TI,DMAT,XL,B,SAAL,SDAL,QUAD,ierr)
 USE CONSTANTS, only: zero, one, two
!
!       Computes those parts of the local stiffness matrix which
!       are due to incompatible shape functions.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NEN=8, NSD=3, NDOF=3, NSTR=6, NMODES=3)
 parameter (NEE=NEN*NDOF)
 parameter (NEX=NMODES*NDOF)
!-pass
 integer ierr
 logical QUAD
 dimension XL(NSD,NEN),DMAT(NSTR,NSTR)
!-locl
 dimension XS(NSD,NSD),SH(4,NEN),B(NSTR,NEE),SHE(4,3), &
   BE(NSTR,NEX),SDAL(NEE,NEX),SAAL(NEX,NEX)
!
!       three extra shape functions and r,- s- and t-derivatives
SHE(4,1) = ONE - RI*RI
    SHE(4,2) = ONE - SI*SI
    SHE(4,3) = ONE - TI*TI
 SHE(1,1) = - TWO*RI
 SHE(1,2) = ZERO
 SHE(1,3) = ZERO
 SHE(2,1) = ZERO
 SHE(2,2) = - TWO*SI
 SHE(2,3) = ZERO
 SHE(3,1) = ZERO
 SHE(3,2) = ZERO
 SHE(3,3) = - TWO*TI
!
!       x- and y-derivatives
!       inverse jacobian evaluated in centre element
 call SHAP30 (ZERO,ZERO,ZERO,XL,DET,SH,XS,QUAD,ierr,2)
 if (ierr.ne.0) return
do i=1,NMODES
    dSHdR = SHE(1,i)
    dSHdS = SHE(2,i)
    dSHdT = SHE(3,i)
    SHE(1,i) = XS(1,1)*dSHdR+XS(1,2)*dSHdS+XS(1,3)*dSHdT
    SHE(2,i) = XS(2,1)*dSHdR+XS(2,2)*dSHdS+XS(2,3)*dSHdT
    SHE(3,i) = XS(3,1)*dSHdR+XS(3,2)*dSHdS+XS(3,3)*dSHdT
 enddo
!
!       extra part strain-displacement matrix
k = 0
do i=1,NMODES
k = k + 1
BE(1,k) = SHE(1,i)
BE(2,k) = zero
BE(3,k) = zero
BE(4,k) = SHE(2,i)
BE(5,k) = SHE(3,i)
BE(6,k) = zero
k = k + 1
BE(1,k) = zero
BE(2,k) = SHE(2,i)
BE(3,k) = zero
BE(4,k) = SHE(1,i)
BE(5,k) = zero
BE(6,k) = SHE(3,i)
k = k + 1
BE(1,k) = zero
BE(2,k) = zero
BE(3,k) = SHE(3,i)
BE(4,k) = zero
BE(5,k) = SHE(1,i)
BE(6,k) = SHE(2,i)
 enddo
!
do j=1,NEX
!          column j of DBE
    DBE1=DMAT(1,1)*BE(1,j)+DMAT(1,2)*BE(2,j)+DMAT(1,3)*BE(3,j) &
    +DMAT(1,4)*BE(4,j)+DMAT(1,5)*BE(5,j)+DMAT(1,6)*BE(6,j)
    DBE2=DMAT(2,1)*BE(1,j)+DMAT(2,2)*BE(2,j)+DMAT(2,3)*BE(3,j) &
    +DMAT(2,4)*BE(4,j)+DMAT(2,5)*BE(5,j)+DMAT(2,6)*BE(6,j)
    DBE3=DMAT(3,1)*BE(1,j)+DMAT(3,2)*BE(2,j)+DMAT(3,3)*BE(3,j) &
    +DMAT(3,4)*BE(4,j)+DMAT(3,5)*BE(5,j)+DMAT(3,6)*BE(6,j)
    DBE4=DMAT(4,1)*BE(1,j)+DMAT(4,2)*BE(2,j)+DMAT(4,3)*BE(3,j) &
    +DMAT(4,4)*BE(4,j)+DMAT(4,5)*BE(5,j)+DMAT(4,6)*BE(6,j)
    DBE5=DMAT(5,1)*BE(1,j)+DMAT(5,2)*BE(2,j)+DMAT(5,3)*BE(3,j) &
    +DMAT(5,4)*BE(4,j)+DMAT(5,5)*BE(5,j)+DMAT(5,6)*BE(6,j)
    DBE6=DMAT(6,1)*BE(1,j)+DMAT(6,2)*BE(2,j)+DMAT(6,3)*BE(3,j) &
    +DMAT(6,4)*BE(4,j)+DMAT(6,5)*BE(5,j)+DMAT(6,6)*BE(6,j)
    do i=1,NEE
  SDAL(i,j) = B(1,i)*DBE1+B(2,i)*DBE2+B(3,i)*DBE3 &
       +B(4,i)*DBE4+B(5,i)*DBE5+B(6,i)*DBE6
    enddo
    do i=1,NEX
  SAAL(i,j) = BE(1,i)*DBE1+BE(2,i)*DBE2+BE(3,i)*DBE3 &
       +BE(4,i)*DBE4+BE(5,i)*DBE5+BE(6,i)*DBE6
    enddo
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine ADDPR (B,X,D,IEN,LM,DUM1,DUM2,IELNO,ISIDE,PRES, &
   LMF,TFAULT,SKEW,DXE)
 USE MODELDEFINITION
 USE CONSTANTS, only: quart, third, half, one,eps
!
! Subroutine to add pressures to load vector
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NEN=8, NSD=3, NDOF=3)
!-pass
 dimension IELNO(*),DUM1(*),IEN(NEN,*),X(NSD,*),D(NDOF,*), &
   ISIDE(*),PRES(*),DUM2(*),B(*),LM(NDOF,NEN,*), &
 LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
 DXE(NDOF,NEN,*)

!-locl
 logical quad
 dimension F(NDOF),XL(NSD,NEN),DL(NDOF,NEN),lien(4),xn(NSD), &
   angle(2),tmp(NSD),rot(3,3)
!-init
 numadd = 0
!
 do k=1,IABS(NUMPR)
!        Set element number N
n = IELNO(k)
quad=(.not.(IEN(3,n).eq.IEN(4,n).and.IEN(7,n).eq.IEN(8,n)))
!        Get local nodal point numbers of the side on which pressure works
call SIDENP(ISIDE(k),quad,lien)
!        Localize coordinates
call LCOORD (X,XL,IEN(1,n))
    if (LGDEF.ne.0) then
       call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
       call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN)
       call ADDSNE (DL,DXE(1,1,n),NDOF,NEN)
       call REZONE(XL,DL,ONE)
    endif
!        Calculate surface area and unit normal vector
call FCGEOM (XL,lien,AREA,XN)
quad=(lien(3).ne.lien(4))
!        Define surface mean of SKEW angles
if (quad) then
angle(1) = ( -SKEW(1,IEN(lien(1),n)) - &
       SKEW(1,IEN(lien(2),n)) - SKEW(1,IEN(lien(3),n)) - &
           SKEW(1,IEN(lien(4),n)) ) * HALF
angle(2) = ( -SKEW(2,IEN(lien(1),n)) - &
       SKEW(2,IEN(lien(2),n)) - SKEW(2,IEN(lien(3),n)) - &
           SKEW(2,IEN(lien(4),n)) ) * HALF
else
angle(1) = ( -SKEW(1,IEN(lien(1),n)) - &
       SKEW(1,IEN(lien(2),n))-SKEW(1,IEN(lien(3),n)) ) * HALF
angle(2) = ( -SKEW(2,IEN(lien(1),n)) - &
       SKEW(2,IEN(lien(2),n))-SKEW(2,IEN(lien(3),n)) ) * HALF
endif
    if (ABS(angle(1)).gt.EPS .or. ABS(angle(2)).gt.EPS) then
       call FORMRT (angle,rot,NDOF)
       call EQUATE (tmp,xn,NDOF)
       call VMPRD (rot,tmp,xn,NDOF,3)
    endif
    if (quad) then
    F(1) = XN(1) * PRES(k) * AREA * QUART
    F(2) = XN(2) * PRES(k) * AREA * QUART
    F(3) = XN(3) * PRES(k) * AREA * QUART
    else
    F(1) = XN(1) * PRES(k) * AREA * THIRD
    F(2) = XN(2) * PRES(k) * AREA * THIRD
    F(3) = XN(3) * PRES(k) * AREA * THIRD
    endif
do l=1,NDOF
ii = LM(l,lien(1),n)
if (ii.ne.0) B(ii) = B(ii) + F(l)
ii = LM(l,lien(2),n)
if (ii.ne.0) B(ii) = B(ii) + F(l)
ii = LM(l,lien(3),n)
if (ii.ne.0) B(ii) = B(ii) + F(l)
ii = LM(l,lien(4),n)
    if (ii.ne.0.and.quad) B(ii) = B(ii) + F(l)
enddo
 enddo
!
if (iecho.eq.1 .and. numadd.gt.0) write(stdout,1)
 1   format(1x,'Pressure loads are being added')
!
 return
 end
!-----------------------------------------------------------------------
 subroutine FCGEOM (XL,lien,AREA,XN)
 USE CONSTANTS, only: two,dfmin
!
! Calculates area of an element face, perpendicular to face normal,
! and the unit normal vector.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NSD=3,NEN=8)
!-pass
 dimension XL(NSD,NEN),lien(4),XN(NSD)

!-locl
 logical quad
 dimension XLL(NSD,4)
!
!    Localize face coordinates
 do 100 i=1,4
XLL(1,i) = XL(1,lien(i))
XLL(2,i) = XL(2,lien(i))
XLL(3,i) = XL(3,lien(i))
100    continue
 quad = (lien(3).ne.lien(4))
!
 if (quad) then
!        Determine (average) normal vector Y
XN(1)=(XLL(2,3)-XLL(2,1))*(XLL(3,4)-XLL(3,2)) &
     +(XLL(2,4)-XLL(2,2))*(XLL(3,1)-XLL(3,3))
XN(2)=(XLL(3,3)-XLL(3,1))*(XLL(1,4)-XLL(1,2)) &
     +(XLL(3,4)-XLL(3,2))*(XLL(1,1)-XLL(1,3))
XN(3)=(XLL(1,3)-XLL(1,1))*(XLL(2,4)-XLL(2,2)) &
     +(XLL(1,4)-XLL(1,2))*(XLL(2,1)-XLL(2,3))
 else
XN(1)=(XLL(2,2)-XLL(2,1))*(XLL(3,3)-XLL(3,1)) &
     -(XLL(3,2)-XLL(3,1))*(XLL(2,3)-XLL(2,1))
XN(2)=(XLL(3,2)-XLL(3,1))*(XLL(1,3)-XLL(1,1)) &
     -(XLL(1,2)-XLL(1,1))*(XLL(3,3)-XLL(3,1))
XN(3)=(XLL(1,2)-XLL(1,1))*(XLL(2,3)-XLL(2,1)) &
     -(XLL(2,2)-XLL(2,1))*(XLL(1,3)-XLL(1,1))
 endif
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
 subroutine LFLUX (XL,TL,IEN,PRPTEM,FLUX,ierr)
!
! Routine to calculate element heat flux
!
 implicit double precision (a-h,o-z)
! SNGL    implicit real (a-h,o-z)
 parameter (NSD=3,NEN=8,NDOF=3,THIRD=-1d0/3d0)
!-pass
 dimension XL(NSD,NEN),TL(NEN),PRPTEM(6),FLUX(NSD),IEN(NEN)

!-locl
 dimension GRAD(NDOF)
 logical quad
 save rn,tn
 data rn,tn/0.,0./
 sn = 0.d0
!
quad=(.not.(IEN(3).eq.IEN(4).and.IEN(7).eq.IEN(8)))
 if (.not.quad) sn = THIRD
!    calculate thermal gradients in the center of the element
 call thgrad (rn,sn,tn,XL,TL,GRAD,quad,ierr)
 if (ierr.ne.0) return
 FLUX(1) = -PRPTEM(4)*GRAD(1)
 FLUX(2) = -PRPTEM(5)*GRAD(2)
 FLUX(3) = -PRPTEM(6)*GRAD(3)
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine EFORCE (EVP,P,XL,DUM1,QUAD,ierr)
   USE MODELDEFINITION
 USE CONSTANTS, only: G
!
! Computes the effective forces P at each node
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NEN=8, NSD=3, NDOF=3, NSTR=6)
 parameter (NEE=NEN*NDOF)
!-pass
 integer ierr
 logical QUAD
 dimension XL(NSD,NEN),P(NEE),EVP(NSTR)

!-locl
 dimension RG(NEN),SG(NEN),TG(NEN),XS(NSD,NSD),SH(4,NEN), &
   B(NSTR,NEE)
 save RG,SG,TG
!-init
 data RG/-1., 1., 1.,-1.,-1., 1., 1.,-1./
data SG/-1.,-1., 1., 1.,-1.,-1., 1., 1./
data TG/-1.,-1.,-1.,-1., 1., 1., 1., 1./
!
 if (IPOINT.eq.1) then
NINT = 1
W    = 8d0
G_local    = 0d0
 else
NINT = 8
W    = 1d0
!        G_local    = GP
 endif
!
 do l=1,NINT
!        Compute spatial derivatives of shape functions
call SHAP30 (RG(l)*G_local, &
             SG(l)*G_local, &
             TG(l)*G_local, &
             XL,DET,SH,XS,QUAD, &
       ierr,3)
if (ierr.ne.0) return
call BMATRIX (B,SH)
do i=1,NEE
EFFE=B(1,i)*EVP(1)+B(2,i)*EVP(2)+B(3,i)*EVP(3) &
       +B(4,i)*EVP(4)+B(5,i)*EVP(5)+B(6,i)*EVP(6)
P(i)=P(i)+W*DET*EFFE
enddo
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine DEVSTN (STN,SDEV)
 USE CONSTANTS, only: two, third
!
! Computes stress deviatior SDEV
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=6)
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
 subroutine FORMBT (STN,BETA,BETB,PRPMAT,PRPLAS)
   USE MODELCTX
 USE CONSTANTS, only: zero, third, one, two, six,dfmin,dlogmin,dfmax,dlogmax
!
! Subroutine to compute strain rates from stress field. It defines BETA and
! evaluates the product BETA*STRESS. BETB is the viscoplastic strain rate.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NSTR=6)
!-pass
 dimension STN(NSTR),BETA(NSTR),BETB(NSTR),PRPMAT(6),PRPLAS(9)
!-locl
 dimension AVECT(NSTR)
!-init
 sigma = ZERO
 emhu  = PRPMAT(3)
 anpwr = PRPMAT(4)
!
 if (ICVIS.ge.2) then
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
 if (ICVIS.ge.2) then
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
 subroutine FORMJB (BETAJ,STN,PRPMAT,PRPLAS)
   USE MODELCTX
 USE CONSTANTS, only: third, one, two, six,dfmin,dlogmin,dfmax,dlogmax
!
! Constructs the Jacobean of the stress-strain rate matrix BETA.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NSTR=6)
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
if (ICVIS.ge.2) then
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
 if (ICVIS.ge.2) then
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
 subroutine ADDSTR (B,X,D,IEN,LM,DUM1,DUM2,IELSTR,ISSIDE, &
   ISTR,LMF,TFAULT,SKEW,DXE)

USE MATERIALSMODULE, only : LMAT
USE MODELDEFINITION
use modeldatamodule, only: modeldatactx
 USE CONSTANTS, only: quart, third,eps
use timestepmodule, only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

!
! Subroutine to add stresses to load vector
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NSD=3,NDOF=3,NSTR=6,NEN=8)
!-pass
 dimension DUM1(*),IEN(NEN,*),X(NSD,*),D(NDOF,*),DUM2(*),B(*), &
   LM(NDOF,NEN,*),IELSTR(*),ISSIDE(*),ISTR(2,*), &
   LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN,*)
!-locl
 logical quad
 dimension f(NDOF),xl(NSD,NEN),dl(NDOF,NEN),lien(4),xn(NSD), &
   tmp(NSD),rot(3,3),angle(2)
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
call AddFaultDisplacement (dl,k,TFAULT,NDOF,NEN)
    call ADDSNE (dl,DXE(1,1,k),NDOF,NEN)
call REZONE (xl,dl,1d0)
endif
!        Get nodal point numbers of the side on which pressure works
quad=(.not.(IEN(3,n).eq.IEN(4,n).and.IEN(7,n).eq.IEN(8,n)))
call SIDENP(ISSIDE(n),quad,lien)
!        Calculate surface area and normal vector
call FCGEOM (xl,lien,area,xn)
quad=(lien(3).ne.lien(4))
!        Rotate face into local direction
angle(1) = 0d0
angle(2) = 0d0
if (quad) then
do i=1,4
    np = IEN(lien(i),k)
    angle(1) = angle(1) - SKEW(1,np)
    angle(2) = angle(2) - SKEW(2,np)
enddo
angle(1) = angle(1)*QUART
angle(2) = angle(2)*QUART
else
do i=1,3
    np = IEN(lien(i),k)
    angle(1) = angle(1) - SKEW(1,np)
    angle(2) = angle(2) - SKEW(2,np)
enddo
angle(1) = angle(1)*THIRD
angle(2) = angle(2)*THIRD
endif
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
    if (quad) then
        f(1)=(sxx*xn(1)+sxy*xn(2)+sxz*xn(3))*area*QUART
        f(2)=(sxy*xn(1)+syy*xn(2)+syz*xn(3))*area*QUART
        f(3)=(sxz*xn(1)+syz*xn(2)+szz*xn(3))*area*QUART
    else
        f(1)=(sxx*xn(1)+sxy*xn(2)+sxz*xn(3))*area*THIRD
        f(2)=(sxy*xn(1)+syy*xn(2)+syz*xn(3))*area*THIRD
        f(3)=(sxz*xn(1)+syz*xn(2)+szz*xn(3))*area*THIRD
    endif
do l=1,NSD
ii=LM(l,lien(1),k)
      if(ii.ne.0) B(ii)=B(ii)+f(l)
      ii=LM(l,lien(2),k)
      if(ii.ne.0) B(ii)=B(ii)+f(l)
      ii=LM(l,lien(3),k)
      if(ii.ne.0) B(ii)=B(ii)+f(l)
      ii=LM(l,lien(4),k)
      if(ii.ne.0.and.quad) B(ii)=B(ii)+f(l)
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
 USE CONSTANTS, only: quart, third, half, one,eps
use timestepmodule, only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

!
! Subroutine to add Winkler pressures to diagonal elements of stiffness
! matrix
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NSD=3,NDOF=3,NEN=8)
!-pass
 dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),A(*),IDIAG(NEQ),LM(NDOF,NEN,*), &
   D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
   DXE(NDOF,NEN,*),LMX(NDOF,NEN,*),NSLIP(5,*)
#ifdef SPARSE


   PetscMPIInt irank
#endif
!-locl
 logical quad
 dimension w(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),rot(3,3), &
   tmp(NSD),lien(4),angle(2)

 if (NUMWNK.eq.0) return

 do 400 nn=1,IABS(NUMWNK)
n = IABS(IWELM(nn))
LL= IWTIME(nn)
if (LL.eq.-1) goto 250
if (LL) 100,400,200
100        if (NSTEP.ge.-LL) goto 400
200        if (NSTEP.lt.(LL-1)) goto 400
250        call LCOORD (X,XL,IEN(1,n))
if (LGDEF.ne.0) then
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN)
    call ADDSNE (DL,DXE(1,1,n),NDOF,NEN)
call REZONE (XL,DL,ONE)
endif
    quad=(.not.(IEN(3,n).eq.IEN(4,n).and.IEN(7,n).eq.IEN(8,n)))
!           Get nodal point numbers of the side on which pressure works
    call SIDENP (IABS(IWSIDE(nn)),quad,lien)
!           Calculate surface area and unit normal vector
    call FCGEOM (XL,lien,area,xn)
    quad=(lien(3).ne.lien(4))
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
    if (quad) then
        np = IEN(lien(4),n)
        angle(1) = (angle(1)-SKEW(1,np))*QUART
        angle(2) = (angle(2)-SKEW(2,np))*QUART
    else
        angle(1) = angle(1)*THIRD
        angle(2) = angle(2)*THIRD
    endif
        if (ABS(angle(1)).gt.EPS.or.ABS(angle(2)).gt.EPS) then
            call FORMRT (angle,rot,NSD)
            call EQUATE (tmp,xn,NSD)
            call VMPRD (rot,tmp,xn,NSD,3)
        endif
endif
    if (quad) then
        w(1) = ABS( xn(1) * area ) * WPRES(nn) * QUART
        w(2) = ABS( xn(2) * area ) * WPRES(nn) * QUART
        w(3) = ABS( xn(3) * area ) * WPRES(nn) * QUART
    else
        w(1) = ABS( xn(1) * area ) * WPRES(nn) * THIRD
        w(2) = ABS( xn(2) * area ) * WPRES(nn) * THIRD
        w(3) = ABS( xn(3) * area ) * WPRES(nn) * THIRD
    endif
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
                if (NSLIP(2,ii).eq.IEN(lien(4),n)) nl = nl + 1
            enddo
else
    ni = 1
    nj = 1
    nk = 1
    nl = 1
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
if (quad) then
    ii = LM(idof,lien(4),n)
        if (ii.gt.0) then
#ifdef SPARSE
        call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,w(idof), &
               ADD_VALUES,ierr)
                if (ierr.ne.0) then
                    write(stderr,1) ierr
                    call stoper()
                endif
#else
        ii = IDIAG(ii)
        A(ii) = A(ii) + w(idof)
#endif
        endif
endif
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
    call MatSetValuesLocal(M_A,1,ii-1,1,ii-1,v,ADD_VALUES, &
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
ii = IABS(LMX(idof,lien(4),n))
if (ii.gt.0 .and. quad) then
    v = w(idof)/DBLE(nl)
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
 subroutine ELCRD (R,S,T,XELM,IEN,X,D,elemID,TFAULT,SKEW,FACTOR,ierr)
 USE CONSTANTS, only: zero, useDXE
 f3d/f3dlib_brick.F
!
! Calculates (deformed) coordinates XELM at natural coordinates (R,S,T)
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSD=3, NDOF=3, NEN=8)
!-pass
 integer ierr
 dimension XELM(NSD),IEN(NEN),X(NSD,*),D(NDOF,*), &
   LMF(NDOF,NEN),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN)
!-locl
 logical quad
 dimension XL(NSD,NEN),DL(NDOF,NEN),SH(4,NEN),XS(NSD,NSD)

!-init
 quad = (.not.(IEN(3).eq.IEN(4).and.IEN(7).eq.IEN(8)))
!
!    Set local coordinates
 call LCOORD (X,XL,IEN)
 if (FACTOR.gt.0.) then
call LDISP  (DL,D,IEN,NDOF,NEN)
call AddFaultDisplacement (DL,elemID,TFAULT,NDOF,NEN)
call ADDSNE (DL,elemID,NDOF,NEN,useDXE)
call REZONE (XL,DL,FACTOR)
 endif
!
!    Compute shape functions
 call SHAP30 (R,S,T,XL,DET,SH,XS,quad,ierr,1)
!
 XELM(1)=ZERO
 XELM(2)=ZERO
 XELM(3)=ZERO
 do 20 j=1,NEN
XELM(1) = XELM(1)+SH(4,j)*XL(1,j)
XELM(2) = XELM(2)+SH(4,j)*XL(2,j)
XELM(3) = XELM(3)+SH(4,j)*XL(3,j)
20    continue
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine ELCTR (XELM,IEN,X,D,LMF,TFAULT, &
 SKEW,DXE,FACTOR,ierr)
   USE MODELDEFINITION
   USE MATERIALSMODULE, only : RADIUS
 USE CONSTANTS, only: third
!
! Calculates (deformed) element center coordinates XELM
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSD=3, NDOF=3, NEN=8)
!    parameter (THIRD=-1d0/3d0) really? minus sign?
!-pass
 integer ierr
 dimension XELM(NSD),IEN(NEN),X(NSD,*),D(NDOF,*), &
   LMF(NDOF,NEN),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN)
!-locl
 save rn,tn
!-init
 data rn,tn/0d0,0d0/
 sn = 0d0
!
!    if (IEN(3).eq.IEN(4).and.IEN(7).eq.IEN(8)) sn=THIRD
 if (IEN(3).eq.IEN(4).and.IEN(7).eq.IEN(8)) sn=-third
 call ELCRD (rn,sn,tn,XELM,IEN,X,D,LMF,TFAULT, &
   SKEW,DXE,FACTOR,ierr)
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine TSTIFF (XL,QUAD,S,C,DUM1,ierr)
 USE CONSTANTS, only: G
!
! computes the local conduction stiffness matrix from integration at Barlow
! points. Solely calculates upper-diagonal matrix.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
  parameter (NEN=8, NINT=8, NSD=3)
!-pass
 integer ierr
 logical QUAD
 dimension XL(NSD,NEN),S(NEN,NEN),C(NSD)
!-locl
 dimension SH(4,NEN),RG(NEN),SG(NEN),TG(NEN),XS(NSD,NSD)
save RG,SG,TG
 data RG/-1., 1., 1.,-1.,-1., 1., 1.,-1./
data SG/-1.,-1., 1., 1.,-1.,-1., 1., 1./
data TG/-1.,-1.,-1.,-1., 1., 1., 1., 1./
!
 do L=1,NINT
!        Calculate spatial derivatives of shape functions
call SHAP30 (RG(L)*G,SG(L)*G,TG(L)*G,XL,DET,SH,XS,QUAD, &
       ierr,3)
if (ierr.ne.0) return
do j=1,NEN
do i=1,j
    S(i,j)=S(i,j)+DET*(C(1)*SH(1,i)*SH(1,j) &
           +C(2)*SH(2,i)*SH(2,j)+C(3)*SH(3,i)*SH(3,j))
enddo
enddo
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine ADVSTF (S,XL,VL,PL,COND,RHOCP,TH,HEAT,QUAD,ierr)

 use constants, only: dfmin,dlogmin,dfmax,dlogmax,zero,third,one,two,three
!
! convection (advection) update of local conduction stiffness matrix for
! steady convection-conduction solution
!
 parameter (NUMOPT=2)
!
! NUMOPT = 0: no upwinding
! NUMOPT = 1: classical upwind scheme
! NUMOPT = 2: optimal upwind scheme
! NUMOPT = 3: double asymptotic approximation scheme
! NUMOPT = 4: critical approximation scheme
! For triangle elements, the Mizukami upwind scheme is used
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
  parameter (NEN=8, NDOF=3, NSD=3)
 parameter (NINT=1,G=ZERO, WGHT=8d0)
!-pass
 logical QUAD
 integer ierr
 double precision RHOCP,TH
 dimension XL(NSD,NEN),VL(NDOF,NEN),S(NEN,NEN),COND(NSD),PL(NEN)

!-locl
 dimension sh(4,NEN),rg(NEN),sg(NEN),tg(NEN),xs(NSD,NSD), &
   size(NSD),SM(NEN,NEN),SA(NEN,NEN),SD(NEN,NEN)
 logical firstentry
save rg,tg,sg,firstentry
!-init
 data rg/-1., 1., 1.,-1.,-1., 1., 1.,-1./
data sg/-1.,-1., 1., 1.,-1.,-1., 1., 1./
data tg/-1.,-1.,-1.,-1., 1., 1., 1., 1./
 data firstentry /.true./
!
 if (firstentry) then
if (iecho.ne.0) then
    if (NUMOPT.eq.0) then
    write(stderr,5)
 5            format(1x,'No upwinding')
else
    if (QUAD) then
        if (NUMOPT.eq.1) then
            write(stderr,10)
10                    format(1x,'Upwind: classical upwind scheme')
        else if (NUMOPT.eq.2) then
            write(stderr,20)
20                    format(1x,'Upwind: optimal upwind scheme')
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
60            format(1x,'Upwind: Mizukami scheme')
    endif
endif
endif
firstentry = .false.
 endif
!
 call CLEAR(SA,NEN*NEN)
 size(1) = ZERO
 size(2) = ZERO
 size(3) = ZERO
 if (NUMOPT.ne.0) then
if (QUAD) then
    call ELMSIZE (XL,QUAD,size)
else
    call CLEAR(SM,NEN*NEN)
    call CLEAR(SD,NEN*NEN)
!            consistent mass matrix SM
    call CPSTIFF (XL,QUAD,SM,RHOCP,TH,ierr)
!            conduction stiffness matrix SD
    call TSTIFF (XL,QUAD,SD,COND,TH,ierr)
endif
 endif
!
do l=1,NINT
call SHAP30 (rg(L)*G,sg(L)*G,tg(L)*G,XL,detJ,sh,xs,QUAD, &
       ierr,3)
if (ierr.ne.0) return
Vx = ZERO
Vy = ZERO
Vz = ZERO
do j=1,NEN
    Vx = Vx + VL(1,j)*sh(4,j)
    Vy = Vy + VL(2,j)*sh(4,j)
    Vz = Vz + VL(3,j)*sh(4,j)
enddo
    Vmag = SQRT(Vx*Vx+Vy*Vy+Vz*Vz)
    if (Vmag.gt.1d-35) then
if (QUAD .and. NUMOPT.ne.0) then
!            Directional grid Peclet numbers
            alphax = RHOCP*Vx*size(1)/(TWO*COND(1))
            alphay = RHOCP*Vy*size(2)/(TWO*COND(2))
            alphaz = RHOCP*Vz*size(3)/(TWO*COND(3))
    if (NUMOPT.eq.1) then
!                classical upwind
        xix = SIGN(ONE,alphax)
        xiy = SIGN(ONE,alphay)
    xiz = SIGN(ONE,alphaz)
    else if (NUMOPT.eq.2) then
!                optimal upwind
        xix = APPRX(alphax)
        xiy = APPRX(alphay)
        xiz = APPRX(alphaz)
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
        if (ABS(alphaz).gt.THREE) then
        xiz = SIGN(ONE,alphaz)
        else
        xiz = alphaz*THIRD
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
        if (alphaz.lt.-ONE) then
        xiz = -ONE-ONE/alphaz
        elseif (alphaz.gt.ONE) then
            xiz =  ONE-ONE/alphaz
        else
        xiz = ZERO
        endif
    endif
    tau=(size(1)*xix*Vx+size(2)*xiy*Vy+size(3)*xiz*Vz)/ &
           (TWO*Vmag*Vmag)
endif ! QUAD and NUMOPT>0
w  = detJ * TH * WGHT
do j=1,NEN
            do i=1,NEN
    SA(i,j) = &
           RHOCP*w*sh(4,i)*(Vx*sh(1,j)+Vy*sh(2,j)+ &
           Vz*sh(3,j))
       enddo
        enddo
do i=1,NEN
    do j=1,NEN
    S(i,j)  = S(i,j) + SA(i,j)
    enddo
enddo
if (NUMOPT.ne.0) then
    if (.not.QUAD) call SYMMET(SM,NEN)
    do i=1,NEN
    if (.not.QUAD) then
        if (i.lt.NEN) then
                samax = ABS(SA(i,1))
            do j=2,NEN
                samax = MAX(ABS(SA(i,j)),samax)
            enddo
!                    element ii of lumped mass matrix
            cplump = SM(i,1)+SM(i,2)+SM(i,3)+ &
               SM(i,5)+SM(i,6)+SM(i,7)
            tau = 5D-1*(cplump/(samax+SD(i,i)))
        else
        tau = ZERO
        endif
    endif
    Qi = tau*(Vx*sh(1,i)+Vy*sh(2,i)+Vz*sh(3,i))
    PL(i) = PL(i) + HEAT*w*Qi
        do j=1,NEN
        S(i,j) = S(i,j) + &
               RHOCP*w*Qi*(Vx*sh(1,j)+Vy*sh(2,j)+ &
               Vz*sh(3,j))
    enddo
    enddo
endif ! NUMOPT > 0
    endif ! Vmag > 0
enddo
!
 return
 end
!-------------------------------------------------------------------------------
subroutine ELMSIZE (XL,QUAD,SIZE)
!
! computes characteristic sizes of element
!
implicit double precision (A-H,O-Z)
parameter (NEN=8, NSD=3)
parameter (ZERO=0d0, THREE=3D0)
!-pass
dimension XL(NSD,NEN),SIZE(NSD)
logical QUAD
!-locl
dimension sh(NSD,NEN),xs(NSD,NSD)
!
 if (QUAD) then
sn = -1d0/THREE
 else
sn = ZERO
 endif
 call SHAP30 (ZERO,sn,ZERO,XL,det,sh,xs,QUAD,ierr,1)
 if (ierr.ne.0) return
!
SIZE(1) = (ABS(xs(1,1))+ABS(xs(2,1))+ABS(xs(3,1)))*THREE ! not sure
SIZE(2) = (ABS(xs(1,2))+ABS(xs(2,2))+ABS(xs(3,2)))*THREE
SIZE(3) = (ABS(xs(1,3))+ABS(xs(2,3))+ABS(xs(3,3)))*THREE
!
return
end
!-------------------------------------------------------------------------------
function APPRX(X)
!
 implicit double precision (a-h,o-z)
 parameter (HALF=5D-1,ONE=1D0)
!-pass
double precision APPRX,X
double precision a,b,y
!
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
!
return
end
!-------------------------------------------------------------------------------
 subroutine LHEAT (XL,QUAD,PL,HEAT,DUM1,ierr)
 USE CONSTANTS, only: G
!
! computes the local load vector contribution of heat production.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
  parameter (NEN=8, NINT=8, NSD=3)
!-pass
 integer ierr
 logical QUAD
 dimension XL(NSD,NEN),PL(NEN)
!-locl
 dimension SH(4,NEN),RG(NEN),SG(NEN),TG(NEN),XS(NSD,NSD)
save RG,SG,TG
 data RG/-1., 1., 1.,-1.,-1., 1., 1.,-1./
data SG/-1.,-1., 1., 1.,-1.,-1., 1., 1./
data TG/-1.,-1.,-1.,-1., 1., 1., 1., 1./
!
 do L=1,NINT
!        Compute shape functions
call SHAP30 (RG(L)*G,SG(L)*G,TG(L)*G,XL,DET,SH,XS,QUAD, &
       ierr,1)
if (ierr.ne.0) return
PL(1) = PL(1) + HEAT*DET*SH(4,1)
PL(2) = PL(2) + HEAT*DET*SH(4,2)
PL(3) = PL(3) + HEAT*DET*SH(4,3)
PL(4) = PL(4) + HEAT*DET*SH(4,4)
PL(5) = PL(5) + HEAT*DET*SH(4,5)
PL(6) = PL(6) + HEAT*DET*SH(4,6)
PL(7) = PL(7) + HEAT*DET*SH(4,7)
PL(8) = PL(8) + HEAT*DET*SH(4,8)
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 function EigMax (STN,E,V,EMHU,NPL,ALPHA)
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
 parameter (NEN=8, NSD=3, NDOF=3)
 parameter (QUART=.25d0, THIRD=1d0/3d0)
!-pass
!    integer NFLX
 dimension BT(*),IEN(NEN,*),X(NSD,*),DUM1(*),DUM2(*),IFLS(*), &
 IFLX(*),BFLX(*),LMT(NEN,*),D(NDOF,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),LMTX(NEN,*),SKEW(2,*),DXE(NDOF,NEN,*)
!-locl
 logical quad
 dimension XL(NSD,NEN),DL(NDOF,NEN),lien(4),XN(NSD)

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
!           Get nodal point numbers of the side
    quad=(.not.(IEN(3,m).eq.IEN(4,m).and.IEN(7,m).eq.IEN(8,m)))
    call SIDENP (iside,quad,lien)
!           Calculate surface area and normal vector
    call FCGEOM (XL,lien,AREA,XN)
quad=(lien(3).ne.lien(4))
    if (quad) then
        F=-AREA*BFLX(n)*QUART
    else
        F=-AREA*BFLX(n)*THIRD
    endif
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
ii = LMT(lien(4),m)
if (ii.eq.0.and.NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) &
       ii = LMTX(lien(4),m)
if (ii.gt.0 .and. quad) BT(ii) = BT(ii) + F
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
 parameter (NSD=3,NDOF=3,NEN=8,ONE=1D0)
!-pass
 dimension ITWINK(3,*),TWINK(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),AA(*),ITDIAG(*),LMT(NEN,*), &
   D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
   DXE(NDOF,NEN,*)
#ifdef SPARSE


   PetscMPIInt irank
#endif
!-locl
 dimension xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),lien(4)
!
 if (NTWINK.eq.0) return
!
if (iecho.ne.0) write(stdout,*) 'Adding Winkler fluxes'
!
 do 400 nn=1,NTWINK
n = ITWINK(1,nn)
LL= ITWINK(3,nn)
if (LL.eq.-1) goto 300
if (LL) 100,400,200
100        if (NSTEP.ge.-LL) goto 400
200        if (NSTEP.lt.(LL-1)) goto 400
300        call LCOORD (X,XL,IEN(1,n))
if (LGDEF.ne.0) then
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
call AddFaultDisplacement (DL,LMF(1,1,n),TFAULT,NDOF,NEN)
    call ADDSNE (DL,DXE(1,1,n),NDOF,NEN)
call REZONE (XL,DL,ONE)
endif
call SIDENP (ITWINK(2,nn),IEN(3,n).ne.IEN(4,n),lien)
    call FCGEOM (XL,lien,area,xn)    ! compute face area
if (lien(3).eq.lien(4)) then
    wt = area * TWINK(nn)/3d0
else
    wt = area * TWINK(nn)/4d0
endif
ii = LMT(lien(1),n)            ! equation number
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
if (lien(3).eq.lien(4)) goto 400
ii = LMT(lien(4),n)
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
!-------------------------------------------------------------------------------
 subroutine THGRAD (R,S,T,XL,TL,GRAD,QUAD,ierr)
!
! Routine for calculating thermal gradients at natural coordinates r,s,t
!
 implicit double precision (a-h,o-z)
! SNGL    implicit real (a-h,o-z)
 parameter (NSD=3, NEN=8)
!-pass
 integer ierr
 logical QUAD
 dimension XL(NSD,NEN),TL(NEN),GRAD(NSD)
!-locl
 dimension SH(4,NEN),XS(NSD,NSD)
!
!    Compute r-, s- and t-derivatives of shape functions and
!    inverse Jacobian matrix XS
 call SHAP30 (R,S,T,XL,DET,SH,XS,QUAD,ierr,2)
 if (ierr.ne.0) return
!
!    calculate r-, s- and t- thermal gradients
 dTdr = SH(1,1)*TL(1)+SH(1,2)*TL(2)+SH(1,3)*TL(3)+SH(1,4)*TL(4) &
   +SH(1,5)*TL(5)+SH(1,6)*TL(6)+SH(1,7)*TL(7)+SH(1,8)*TL(8)
 dTds = SH(2,1)*TL(1)+SH(2,2)*TL(2)+SH(2,3)*TL(3)+SH(2,4)*TL(4) &
   +SH(2,5)*TL(5)+SH(2,6)*TL(6)+SH(2,7)*TL(7)+SH(2,8)*TL(8)
 dTdt = SH(3,1)*TL(1)+SH(3,2)*TL(2)+SH(3,3)*TL(3)+SH(3,4)*TL(4) &
   +SH(3,5)*TL(5)+SH(3,6)*TL(6)+SH(3,7)*TL(7)+SH(3,8)*TL(8)
!
!    Compute thermal gradient vector
 GRAD(1) = dTdr*XS(1,1)+dTds*XS(1,2)+dTdt*XS(1,3)
 GRAD(2) = dTdr*XS(2,1)+dTds*XS(2,2)+dTdt*XS(2,3)
 GRAD(3) = dTdr*XS(3,1)+dTds*XS(3,2)+dTdt*XS(3,3)
!
 return
 end
!-----------------------------------------------------------------------
 function Tcentr (T,TX,XL,IEN,LMTX,IDUM1)
   USE MODELDEFINITION
!
! Calculates the temperature in the center of the element
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NEN=8, NSD=3)
 parameter (ZERO=0.d0, THIRDM=-1d0/3d0)
!-pass
 dimension T(*),TX(*),IEN(NEN),LMTX(NEN),XL(NSD,NEN)

!-locl
 logical quad
 dimension TL(NEN),SH(4,NEN),XS(NSD,NSD)
!-init
!
 quad=(.not.(IEN(3).eq.IEN(4).and.IEN(7).eq.IEN(8)))
!
!    Localize temperatures
 call LTEMP (IEN,T,TL,NEN)
!
!    Modify for discontinuities
 if (NUMSLP+NUMFN.gt.0.and.IDIFT.eq.1) &
   call LTEMPX (LMTX,IEN,TX,TL,NEN)
!
!    Compute shape functions
 sn = ZERO
 if (.not.quad) sn=THIRDM
 call SHAP30 (ZERO,sn,ZERO,XL,DET,SH,XS,quad,ierr,1)
!
!    Calculate center temperature
 Tcentr = SH(4,1)*TL(1)+SH(4,2)*TL(2)+SH(4,3)*TL(3)+SH(4,4)*TL(4) &
   +SH(4,5)*TL(5)+SH(4,6)*TL(6)+SH(4,7)*TL(7)+SH(4,8)*TL(8)
!
 return
 end
!-----------------------------------------------------------------------
 subroutine CPSTIFF (XL,QUAD,S,RhoCp,DUM1,ierr)
 USE CONSTANTS, only: G
!
! computes the upper local heat capacity stiffness matrix from
! integration at Barlow points.
!
 implicit double precision (a-h,o-z)
! SNGL    implicit real (a-h,o-z)
  parameter (NEN=8, NINT=8, NSD=3)
!-pass
 integer ierr
 logical QUAD
 dimension XL(NSD,NEN),S(NEN,NEN)
!-locl
 dimension SH(4,NEN),RG(NEN),SG(NEN),TG(NEN),XS(NSD,NSD)
save RG,SG,TG
!-init
 data RG/-1., 1., 1.,-1.,-1., 1., 1.,-1./
data SG/-1.,-1., 1., 1.,-1.,-1., 1., 1./
data TG/-1.,-1.,-1.,-1., 1., 1., 1., 1./
!
 do L=1,NINT
!        Compute shape functions
call SHAP30(RG(L)*G,SG(L)*G,TG(L)*G,XL,DET,SH,XS,QUAD, &
       ierr,1)
if (ierr.ne.0) return
do j=1,NEN
do i=1,j
    S(i,j)=S(i,j)+RhoCp*SH(4,i)*SH(4,j)*DET
enddo
enddo
 enddo
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
 parameter (NSD=3, NDOF=3, NEN=8)
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
XL(1,5) = XL(1,5) + DL(1,5)*FACTOR
XL(2,5) = XL(2,5) + DL(2,5)*FACTOR
XL(3,5) = XL(3,5) + DL(3,5)*FACTOR
XL(1,6) = XL(1,6) + DL(1,6)*FACTOR
XL(2,6) = XL(2,6) + DL(2,6)*FACTOR
XL(3,6) = XL(3,6) + DL(3,6)*FACTOR
XL(1,7) = XL(1,7) + DL(1,7)*FACTOR
XL(2,7) = XL(2,7) + DL(2,7)*FACTOR
XL(3,7) = XL(3,7) + DL(3,7)*FACTOR
XL(1,8) = XL(1,8) + DL(1,8)*FACTOR
XL(2,8) = XL(2,8) + DL(2,8)*FACTOR
XL(3,8) = XL(3,8) + DL(3,8)*FACTOR
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine TSTRN (EE,DT,VEXP)
!
 USE CONSTANTS, only: zero, third

! calculates strain due to thermal expansion
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
!-pass
 dimension EE(6)
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
!
 USE CONSTANTS, six_int, two, third, deg2rad, root3,dfmin,dlogmin,dfmax,dlogmax

! Routine to calculate the yield stress for various types of plasticity.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
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
 USE CONSTANTS, only: zero, one, two, three, eight, third, root3,  &
    deg2rad, six_int, THF,dfmin,dlogmin,dfmax,dlogmax

! Evaluate the plastic flow vector AVECT.(Owen & Hinton, sec. 7.4)
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=six_int)
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
 USE CONSTANTS, only: one, two, three, four, threeq, root3,dfmin,dlogmin,dfmax,dlogmax,eps
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

 USE CONSTANTS, only: deg2rad, root3, one three, six, half,dfmin,dlogmin,dfmax,dlogmax

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

   use constants, only: eps,one

!
! Compute the viscoplastic strain rate matrix BETB
! (Ref.: Owen & Hinton, sec. 8.9)
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (a-h,o-z)
 parameter (NSTR=6)
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
 parameter (NSD=3,NDOF=3,NEN=8)
 parameter (THIRD=1d0/3d0, QUART=.25d0)
!-pass
 dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
   DUM1(*),X(NSD,*),DUM2(6,*),LM(NDOF,NEN,*),BTOT(*), &
   DELD(NDOF,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN,*),LMX(NDOF,NEN,*),NSLIP(5,*), &
 DX(NDOF,*)
!-locl
 logical quad
 dimension w(NDOF),xl(NSD,NEN),dl(NDOF,NEN),lien(4),xn(NSD), &
   tmp(NSD),rot(3,3),angle(2)
!
 if (NUMWNK.eq.0) return
!
 do 400 nn=1,IABS(NUMWNK)
n = IABS(IWELM(nn))
LL= IWTIME(nn)
if (LL.eq.-1) goto 250
if (LL) 100,400,200
100        if (NSTEP.ge.-LL) goto 400
200        if (NSTEP.lt.(LL-1)) goto 400
250        call LCOORD (X,XL,IEN(1,n))
if (LGDEF.ne.0) then
call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
call AddFaultDisplacement (DL,LMF(1,1,n),TFAULT,NDOF,NEN)
    call ADDSNE (DL,DXE(1,1,n),NDOF,NEN)
call REZONE (XL,DL,1d0)
endif
quad=(.not.(IEN(3,n).eq.IEN(4,n).and.IEN(7,n).eq.IEN(8,n)))
!           Get nodal point numbers of the side on which pressure works
    call SIDENP (IABS(IWSIDE(nn)),quad,lien)
i = lien(1)
j = lien(2)
k = lien(3)
l = lien(4)
!           Calculate surface area
    call FCGEOM (XL,lien,area,xn)
    quad=(lien(3).ne.lien(4))
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
        if (quad) then
            np = IEN(l,n)
            angle(1) = (angle(1)-SKEW(1,np))*QUART
            angle(2) = (angle(2)-SKEW(2,np))*QUART
        else
            angle(1) = angle(1)*THIRD
            angle(2) = angle(2)*THIRD
        endif
        if (ABS(angle(1)).gt.EPS.or.ABS(angle(2)).gt.EPS) then
            call FORMRT (angle,rot,NSD)
            call EQUATE (tmp,xn,NSD)
            call VMPRD (rot,tmp,xn,NSD,3)
        endif
    endif
    if (quad) then
        w(1) = WPRES(nn) * ABS( area * xn(1) ) * QUART
        w(2) = WPRES(nn) * ABS( area * xn(2) ) * QUART
        w(3) = WPRES(nn) * ABS( area * xn(3) ) * QUART
    else
        w(1) = WPRES(nn) * ABS( area * xn(1) ) * THIRD
        w(2) = WPRES(nn) * ABS( area * xn(2) ) * THIRD
        w(3) = WPRES(nn) * ABS( area * xn(3) ) * THIRD
    endif
    if (NUMSLP.gt.0) then
        if (IWELM(nn).gt.0) then
            ni = 0
            nj = 0
            nk = 0
            nl = 0
            do ii = 1,NUMSLP
                if (NSLIP(2,ii).eq.IEN(i,n)) ni = ni + 1
                if (NSLIP(2,ii).eq.IEN(j,n)) nj = nj + 1
                if (NSLIP(2,ii).eq.IEN(k,n)) nk = nk + 1
                if (NSLIP(2,ii).eq.IEN(l,n)) nl = nl + 1
            enddo
        else
            ni = 1
            nj = 1
            nk = 1
            nl = 1
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
if (quad) then
     ii = LM(idof,l,n)
     if (ii.gt.0) BTOT(ii)=BTOT(ii)- &
            w(idof)*DELD(idof,IEN(l,n))
endif
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
ii = IABS(LMX(idof,l,n))
      if (ii.gt.0.and.quad) BTOT(ii) = BTOT(ii) - &
       w(idof)*DX(idof,IEN(l,n))/DBLE(nl)
300        continue
400    continue
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine SHAREL (IEN,IJSUR,NUMEL,NEL)
!
!     Finds elements neighbouring to element NEL.
!    Neighbouring elements are identified by three or more shared nodes.
!
 parameter (NSUR=6,NEN=8)
!-pass
 integer NEL,NUMEL
 dimension IEN(NEN,*),IJSUR(NSUR,*)

!-locl
 logical PRISM,quad
 integer NINDX
 external NINDX
!
!    set to last IJSUR(isur,NEL) <> 0 (last index that was written)
 PRISM=(IEN(3,NEL).eq.IEN(4,NEL).and.IEN(7,NEL).eq.IEN(8,NEL))
 isur = NINDX (IJSUR,NEL)
!    loop other elements
 do 200 iel=NEL+1,NUMEL
!        compare nodal points of element NEL with those of IEL > NEL
quad=(.not.(IEN(3,iel).eq.IEN(4,iel).and. &
       IEN(7,iel).eq.IEN(8,iel)))
nshare = 0
do 100 n=1,NEN
if (PRISM.and.(n.eq.4 .or. n.eq.8)) goto 100
if (IEN(n,NEL).eq.IEN(1,iel)) nshare = nshare+1
if (IEN(n,NEL).eq.IEN(2,iel)) nshare = nshare+1
if (IEN(n,NEL).eq.IEN(3,iel)) nshare = nshare+1
if (quad.and.IEN(n,NEL).eq.IEN(4,iel)) nshare = nshare+1
if (IEN(n,NEL).eq.IEN(5,iel)) nshare = nshare+1
if (IEN(n,NEL).eq.IEN(6,iel)) nshare = nshare+1
if (IEN(n,NEL).eq.IEN(7,iel)) nshare = nshare+1
if (quad.and.IEN(n,NEL).eq.IEN(8,iel)) nshare = nshare+1
100        continue
if (nshare.gt.4) then
write(stderr,1) NEL,iel,nshare
 1           format(1x,'sharel: elements ',I5,' and ',I5,' share ', &
         I1,' nodal points')
call xit(1," ")
else if (nshare.eq.3 .or. nshare.eq.4) then
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
integer, parameter (NSUR=6)
!-pass
integer :: NEL
 dimension IJSUR(NSUR,*)

!local
integer :: i
!
!    find free index
i = 1
100 if (i.le.NSUR) then
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
USE CONSTANTS, only: half, one, two, three, four, six, root3, deg2rad,dfmin,dlogmin,dfmax,dlogmax

 
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

if (DEBAR.gt.dfmin) then
    DEBAR = SQRT(DEBAR)
endif

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

if (STNHRD.gt.dfmin) then
    COHESN = COHESN + STNHRD*DEBAR
endif

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

if (FLUIDY.gt.dfmin) then
    PRPLAS(4) = FLUIDY
endif

return
end
!-------------------------------------------------------------------------------
subroutine BPMATRIX (XL,DL,QUAD,BP,ierr)

use constants, only: zero
!
! Sets up derivatives matrix BP
!
! SNGL    implicit real (a-h,o-z)
implicit none

integer, parameter :: NDOF=3
integer, parameter :: NSD=3
integer, parameter :: NEN=8

!-pass
logical :: QUAD
double precision :: XL(NSD,NEN)
double precision :: DL(NDOF,NEN)
double precision :: BP(NDOF,NSD)
!-locl
double precision :: SH(4,NEN)
double precision :: XS(NSD,NSD)
integer          :: ierr
!
!    Compute spatial derivatives of shape functions
call SHAP30(ZERO,ZERO,ZERO,XL,DET,SH,XS,QUAD,ierr,3)

if (ierr.ne.0) then
    return
endif
!
!    BP(1,1) = SH(1,1)*DL(1,1)+SH(1,2)*DL(1,2)+SH(1,3)*DL(1,3)+
!     >     SH(1,4)*DL(1,4)+SH(1,5)*DL(1,5)+SH(1,6)*DL(1,6)+
!     >     SH(1,7)*DL(1,7)+SH(1,8)*DL(1,8)
!    BP(2,1) = SH(1,1)*DL(2,1)+SH(1,2)*DL(2,2)+SH(1,3)*DL(2,3)+
!     >     SH(1,4)*DL(2,4)+SH(1,5)*DL(2,5)+SH(1,6)*DL(2,6)+
!     >     SH(1,7)*DL(2,7)+SH(1,8)*DL(2,8)
!    BP(3,1) = SH(1,1)*DL(3,1)+SH(1,2)*DL(3,2)+SH(1,3)*DL(3,3)+
!     >     SH(1,4)*DL(3,4)+SH(1,5)*DL(3,5)+SH(1,6)*DL(3,6)+
!     >     SH(1,7)*DL(3,7)+SH(1,8)*DL(3,8)
!    BP(1,2) = SH(2,1)*DL(1,1)+SH(2,2)*DL(1,2)+SH(2,3)*DL(1,3)+
!     >     SH(2,4)*DL(1,4)+SH(2,5)*DL(1,5)+SH(2,6)*DL(1,6)+
!     >     SH(2,7)*DL(1,7)+SH(2,8)*DL(1,8)
!    BP(2,2) = SH(2,1)*DL(2,1)+SH(2,2)*DL(2,2)+SH(2,3)*DL(2,3)+
!     >     SH(2,4)*DL(2,4)+SH(2,5)*DL(2,5)+SH(2,6)*DL(2,6)+
!     >     SH(2,7)*DL(2,7)+SH(2,8)*DL(2,8)
!    BP(3,2) = SH(2,1)*DL(3,1)+SH(2,2)*DL(3,2)+SH(2,3)*DL(3,3)+
!     >     SH(2,4)*DL(3,4)+SH(2,5)*DL(3,5)+SH(2,6)*DL(3,6)+
!     >     SH(2,7)*DL(3,7)+SH(2,8)*DL(3,8)
!    BP(1,3) = SH(3,1)*DL(1,1)+SH(3,2)*DL(1,2)+SH(3,3)*DL(1,3)+
!     >     SH(3,4)*DL(1,4)+SH(3,5)*DL(1,5)+SH(3,6)*DL(1,6)+
!     >     SH(3,7)*DL(1,7)+SH(3,8)*DL(1,8)
!    BP(2,3) = SH(3,1)*DL(2,1)+SH(3,2)*DL(2,2)+SH(3,3)*DL(2,3)+
!     >     SH(3,4)*DL(2,4)+SH(3,5)*DL(2,5)+SH(3,6)*DL(2,6)+
!     >     SH(3,7)*DL(2,7)+SH(3,8)*DL(2,8)
!    BP(3,3) = SH(3,1)*DL(3,1)+SH(3,2)*DL(3,2)+SH(3,3)*DL(3,3)+
!     >     SH(3,4)*DL(3,4)+SH(3,5)*DL(3,5)+SH(3,6)*DL(3,6)+
!     >     SH(3,7)*DL(3,7)+SH(3,8)*DL(3,8)
!
call CLEAR(BP,NDOF*NSD)

do k=1,NEN
    do j=1,NSD
        BP(1,j) = BP(1,j)+SH(j,k)*DL(1,k)
        BP(2,j) = BP(2,j)+SH(j,k)*DL(2,k)
        BP(3,j) = BP(3,j)+SH(j,k)*DL(3,k)
    enddo
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine sdimnt (ISURF,SURF,B,ID,GRAV,PRPSED, DELTP)

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
 subroutine SIDENP (is,QUAD,lien)
!
! routine sets (local) nodal point numbers belonging to side IS.
! NPs in lien are positively oriented.
!
!-pass
 logical QUAD
 dimension lien(4)
!
 if (is.eq.1) then
!        frontal side
lien(1)=1
lien(2)=2
lien(3)=6
lien(4)=5
 elseif (is.eq.2) then
!        right side
lien(1)=2
lien(2)=3
lien(3)=7
lien(4)=6
 elseif (is.eq.3) then
if (QUAD) then
!        back side
lien(1)=3
lien(2)=4
lien(3)=8
lien(4)=7
else
lien(1)=3
lien(2)=1
lien(3)=5
lien(4)=7
endif
 else if (is.eq.4) then
if (QUAD) then
!        left side
lien(1)=4
lien(2)=1
lien(3)=5
lien(4)=8
else
lien(1)=3
lien(2)=1
lien(3)=5
lien(4)=7
endif
 else if (is.eq.5) then
if (QUAD) then
!        bottom side
lien(1)=4
lien(2)=3
lien(3)=2
lien(4)=1
else
lien(1)=3
lien(2)=2
lien(3)=1
lien(4)=1
endif
 else if (is.eq.6) then
if (QUAD) then
!        top side
lien(1)=5
lien(2)=6
lien(3)=7
lien(4)=8
else
lien(1)=5
lien(2)=6
lien(3)=7
lien(4)=7
endif
 else
write(stderr,1) is
 1        format(1x,'SIDENP: a brick only has 6 sides, not ',I2)
call xit(1," ")
 endif
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
! routine, and routine that calls this one (elwinx) is never used
 subroutine towinkx (NSLIP,IEN,WINX,IWINX,X,nel,iside,mode,f, &
   NWINKX)
 USE MODELDEFINITION
 use constants, only: eps

!
!    Geometry-dependent part of elwinx, where stresses are converted
!    to forces. Updates IWINX, WINX and NWINKX.
!
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
 parameter (NEN=8,NDOF=3,NSD=3)
!-pass
 dimension NSLIP(5,*),IEN(NEN,*),WINX(NDOF,*),IWINX(NDOF,*), &
   X(NSD,*),f(*)

!-locl
 parameter (QUART=0.25d0, THIRD=0.333333333333d0)
 logical quad
 character(len=1) axes
 dimension XL(NSD,NEN),lien(4),W(3),XN(NSD),XS(NSD),XT(NSD), &
   DW(NSD),axes(NSD),islp(4)
 save axes
 data axes/'x','y','z'/
!
!    Some preliminary checks
 if (iside.lt.1 .or. iside.gt.6) then
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
!       Get nodal point numbers of the side
quad=(.not.(IEN(3,nel).eq.IEN(4,nel).and. &
   IEN(7,nel).eq.IEN(8,nel)))
call SIDENP (iside,quad,lien)
quad=(lien(3).ne.lien(4))
!    Check that each of the side face nodes occur in NSLIP
 if (quad) then
m = 4
 else
m = 3
 endif
!    Set up entry numbers in NSLIP array into ISLP
 do 200 i=1,m
np = IEN(lien(i),nel)
islp(i) = 1
100        if (NSLIP(1,islp(i)).eq.nel.and.NSLIP(2,islp(i)).eq.np) &
       goto 200
islp(i) = islp(i) + 1
if (islp(i).gt.NUMSLP) then
write(stderr,3) nel,np,iside
 3           format(1x,'elwinx: element-node ',I5,'-',I5,1x, &
         'not in NSLIP array (iside=',i1,')')
call xit(1," ")
endif
goto 100
200    continue
 call LCOORD (X,XL,IEN(1,nel))
!       Calculate surface area
call FCGEOM (XL,lien,AREA,XN)
if (quad) then
    W(1) = AREA * QUART
    W(2) = AREA * QUART
    W(3) = AREA * QUART
else
    W(1) = AREA * THIRD
    W(2) = AREA * THIRD
    W(3) = AREA * THIRD
endif
 if (ABS(ABS(XN(3))-1d0).le.1d-3) then
!        Near-vertical plane normal, f(2) is tangential in x-dir, f(3) in y
DW(1) = ABS(XN(1)*W(1)*f(1)) + ABS(W(1)*f(2))
DW(2) = ABS(XN(2)*W(2)*f(1)) + ABS(W(2)*f(3))
DW(3) = ABS(XN(3)*W(3)*f(1))
do i=1,m
np = IEN(lien(i),nel)
do 300 k=1,NSD
    if (DW(k).le.EPS) goto 300
    if (NSLIP(2+k,islp(i)).eq.0) goto 300
    WINX(k,np) = WINX(k,np) + DW(k)
    if (IWINX(k,np).eq.0) then
    NWINKX = NWINKX + 1
    IWINX(k,np) = mode
    goto 300
    else if (IWINX(k,np).eq.mode) then
    goto 300
    endif
    write(stderr,4) nel,np,iside,axes(k)
 4               format(1x,'elwinx: conflicting MODE for element ', &
             I5,', node ',I5,' (iside=',I1,', ',A1, &
             '-component)')
    call xit(1," ")
300        continue
enddo
 else
!        XS results from outer product between plane normal and z
s = SQRT(XN(1)*XN(1)+XN(2)*XN(2))
XS(1) =  XN(2)/s
XS(2) = -XN(1)/s
!        XT results from outer product between plane normal and XS
XT(1) = XN(2)*XS(3)-XN(3)*XS(2)
XT(2) = XN(3)*XS(1)-XN(1)*XS(3)
XT(3) = XN(1)*XS(2)-XN(2)*XS(1)
!        f(2) is horizontal tangential, f(3) is perpendicular to this
do k=1,NSD
DW(k) = ABS(XN(k)*W(k)*f(1)) + ABS(XS(k)*W(k)*f(2)) &
       + ABS(XT(k)*W(k)*f(3))
enddo
do i=1,m
np = IEN(lien(i),nel)
do 400 k=1,NSD
    if (DW(k).le.EPS) goto 400
    if (NSLIP(2+k,islp(i)).eq.0) goto 400
    WINX(k,np) = WINX(k,np) + DW(k)
    if (IWINX(k,np).eq.0) then
    NWINKX = NWINKX + 1
    IWINX(k,np) = mode
    goto 400
    else if (IWINX(k,np).eq.mode) then
    goto 400
    endif
    write(stderr,4) nel,np,iside,axes(k)
    call xit(1," ")
400        continue
enddo
 endif
!
 return
 end
!-----------------------------------------------------------------------
   subroutine TANGENT (X,D,TFAULT,DXE,IEN,LMF,SKEW,lien,ANGLE)
   USE MODELDEFINITION
   USE CONSTANTS, only: pi, halfpi, one, eps
!
! Returns the angle of the tangent of a deformed element side.
! Assumption here is that only the first Euler angle changes.
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (A-H,O-Z)
parameter (NSD=3, NEN=8, NDOF=3)
!-pass
dimension IEN(*),X(NSD,*),D(NDOF,*),SKEW(2,*),ANGLE(2), &
  LMF(NDOF,*),TFAULT(NDOF,*),LIEN(4),DXE(NDOF,*)

#include "coord.i"
!-locl
dimension XL(NSD,NEN),DL(NDOF,NEN),xn(NSD)
!
!       localize coordinate data
call LCOORD (X,XL,IEN)
!
!       deformed coordinates XL
call LDISP  (DL,D,IEN,NDOF,NEN)
call AddFaultDisplacement (DL,LMF,TFAULT,NDOF,NEN)
call ADDSNE (DL,DXE,NDOF,NEN)
call REZONE (XL,DL,ONE)
 call FCGEOM (XL,lien,area,xn)
!
 if (ABS(xn(IY)).gt.EPS) then
ANGLE(1) = ATAN(xn(IX)/xn(IY))
 else
    if (xn(IX).gt.EPS) then
        ANGLE(1) = HALFPI
    else if (xn(IX).lt.-EPS) then
        ANGLE(1) = -HALFPI
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
subroutine LCOORD (X,XL,IEN)
USE MODELDEFINITION, only : ISHELL
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
k = IEN(5)
XL(1,5) = X(1,k)
XL(2,5) = X(2,k)
XL(3,5) = X(3,k)
k = IEN(6)
XL(1,6) = X(1,k)
XL(2,6) = X(2,k)
XL(3,6) = X(3,k)
k = IEN(7)
XL(1,7) = X(1,k)
XL(2,7) = X(2,k)
XL(3,7) = X(3,k)
k = IEN(8)
XL(1,8) = X(1,k)
XL(2,8) = X(2,k)
XL(3,8) = X(3,k)
!
 if (ISHELL.ne.0) call SHLCRD (XL,(IEN(3).ne.IEN(4)))
!
return
end
!-----------------------------------------------------------------------
   subroutine SHLCRD (XL,QUAD)
   USE MODELDEFINITION, only : ISHELL
   USE MATERIALSMODULE, only : RADIUS
   USE ALGEBRA
   USE CONSTANTS, only: deg2rad, zero, half, one,eps
!
! SNGL  implicit real (a-h,o-z)
implicit double precision (A-H,O-Z)
 parameter (NEN=8,NSD=3)
!-pass
dimension XL(NSD,*)
 logical QUAD

#include "coord.i"
!-locl
 double precision north
 dimension v(3,NEN),a(3),b(3),xn(3),north(3),east(3)
!
!    Derive local coordinate directions
 xlon = XL(IX,1)*deg2rad
 ylat = XL(IY,1)*deg2rad
 v(1,1) = RADIUS*COS(ylat)
 v(2,1) = v(1,1)*SIN(xlon)
 v(1,1) = v(1,1)*COS(xlon)
 v(3,1) = RADIUS*SIN(ylat)
 xlon = XL(IX,2)*deg2rad
 ylat = XL(IY,2)*deg2rad
 v(1,2) = RADIUS*COS(ylat)
 v(2,2) = v(1,2)*SIN(xlon)
 v(1,2) = v(1,2)*COS(xlon)
 v(3,2) = RADIUS*SIN(ylat)
 xlon = XL(IX,3)*deg2rad
 ylat = XL(IY,3)*deg2rad
 v(1,3) = RADIUS*COS(ylat)
 v(2,3) = v(1,3)*SIN(xlon)
 v(1,3) = v(1,3)*COS(xlon)
 v(3,3) = RADIUS*SIN(ylat)
 if (QUAD) then
xlon = XL(IX,4)*deg2rad
ylat = XL(IY,4)*deg2rad
v(1,4) = RADIUS*COS(ylat)
v(2,4) = v(1,4)*SIN(xlon)
v(1,4) = v(1,4)*COS(xlon)
v(3,4) = RADIUS*SIN(ylat)
 endif
 a(1) = v(1,2) - v(1,1)
 a(2) = v(2,2) - v(2,1)
 a(3) = v(3,2) - v(3,1)
 b(1) = v(1,3) - v(1,1)
 b(2) = v(2,3) - v(2,1)
 b(3) = v(3,3) - v(3,1)
 call CROSSP(a,b,xn)
 if (QUAD) then
a(1) = v(1,4) - v(1,1)
a(2) = v(2,4) - v(2,1)
a(3) = v(3,4) - v(3,1)
call CROSSP(b,a,a)
xn(1) = (xn(1) + a(1))*HALF
xn(2) = (xn(2) + a(2))*HALF
xn(3) = (xn(3) + a(3))*HALF
 endif
 s = SQRT(DOT(xn,xn,3))
 if (s.lt.EPS) then
write(stderr,*) 'SHCRD: zero normal'
call stoper()
 endif
 call BMULT (xn,3,ONE/s)
 east(1) = -xn(2)
 east(2) =  xn(1)
 east(3) = ZERO
 s = SQRT(DOT(east,east,2))
 if (s.le.EPS) then
write(stderr,1)
 1        format(1X,'SHCRD: east direction cannot be derived')
call stoper()
 endif
 call BMULT (east,2,ONE/s)
 call CROSSP(xn,east,north)
 s = SQRT(DOT(north,north,3))
 if (s.le.EPS) then
write(stderr,2)
 2        format(1x,'SHCRD: north direction cannot be derived')
call stoper()
 endif
 call BMULT (north,3,ONE/s)
!
!    Project node coordinates on local east and north directions
 XL(IX,1) = DOT(v(1,1),east,2)
 XL(IY,1) = DOT(v(1,1),north,3)
 XL(IX,2) = DOT(v(1,2),east,2)
 XL(IY,2) = DOT(v(1,2),north,3)
 XL(IX,3) = DOT(v(1,3),east,2)
 XL(IY,3) = DOT(v(1,3),north,3)
 if (QUAD) then
XL(IX,4) = DOT(v(1,4),east,2)
XL(IY,4) = DOT(v(1,4),north,3)
 else
XL(IX,4) = XL(IX,3)
XL(IY,4) = XL(IY,3)
 endif
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
 implicit none

 parameter (NSTR=6)
!-pass
 double precision :: STN(NSTR)
 double precision :: SL(3,3)

!
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
 USE CONSTANTS, only: third

 implicit none

 double precision, dimension (:) ::  S
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
