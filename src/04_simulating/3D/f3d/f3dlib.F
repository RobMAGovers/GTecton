!    These routines are the rheology-dependent portion of GTECTON
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
!    EffectiveDislCreepRate   Effective strain rate for dislocation creep only
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
subroutine INTCON ()

USE LABELMODULE,     only: LABELC, &
                           LABELD, &
                           LABELS, &
                           LABELR, &
                           MATLBL, &
                           CRPLBL, &
                           TEMLBL

USE MODELDEFINITION, only: nType, NSUR, NEE, NEC, simVariety, &
                           allTemperaturesFixed
use modeltopology,   only: NSD,NSTR,NEN, NDOF
use constants,       only: isF3d
!
! Set variables that depend on topology
!
implicit none
!-pass
! integer NSD,NSTR,NEN, NDOF

simVariety = isF3d

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

allTemperaturesFixed = .false.

!
 return
 end
!-------------------------------------------------------------------------------
! tested 25/6/2003   
subroutine SHPTET (X,Y,Z,XL,DET,SH,IERR,IFLAG)

USE CONSTANTS, only: sixth, one
use iomodule, only: stderr
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
double precision :: A1,A2,A3,A4
double precision :: B1,B2,B3,B4
double precision :: C1,C2,C3,C4
double precision :: D1,D2,D3,D4
double precision :: RDET
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

if (DET.le.0d0) then
    goto 1000
endif

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
write(stderr,*) 'Cause: too much deformation. Solve by reducing force or timespan.'

IERR = 1
return
end
!-------------------------------------------------------------------------------

subroutine elmtest (X,IEN)

USE MESHDATAMODULE
USE MODELDEFINITION
use modelctx,         only: getrank
use constants,        only: eps
use debugmodule,      only: iecho
use meshdatamodule,   only: meshdatactx
use algebra,          only: volumeOfTetrahedron
use iomodule
!
! tests element definitions
! whether the sequence of the vertices is counter clockwise -> positive determinant
!
! Note that is usually redundant.
! Gmsh2tecton had already executed an elmtest!
! Only when a mesh has been made without using gmsh,
! this should be called.
! However, because it does not take that much time, 
! it is left in for robustness' sake.

implicit none

integer, parameter :: NEN=4,NSD=3
!-pass
integer            :: IEN
double precision   :: X
dimension          :: X(NSD,*),IEN(NEN,*)

!-locl
integer            :: i
integer            :: elemID, tempnode
integer            :: p(4)  ! indices of 4 points 
double precision   :: volume
double precision   :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4


if (NUMEL.lt.1) then
    write(stderr,10)
10        format(1x,'NUMEL < 1')
    call exitp(1)
endif
!


do elemID=1,NUMEL

     ! determine nodal points of this particular element
     do i=1,4
        p(i) = IEN(i,elemID)
     enddo

     ! makes sure that p1...4 are sensible numbers

     do i=1,4
        if (p(i).lt.1 ) then
            write(stderr,*) 'Error: node number ', p(i), 'less than 1'
            call exitp(1)
            return
        else if (p(i).gt.meshdatactx%neglobal) then
            write(stderr,*) 'Error: node number ', p(i), &
                    'greater than max of ', meshdatactx%neglobal
            call exitp(1)
            return
        endif
    enddo

    x1 = X(1,p(1))
    y1 = X(2,p(1))
    z1 = X(3,p(1))

    x2 = X(1,p(2))
    y2 = X(2,p(2))
    z2 = X(3,p(2))

    x3 = X(1,p(3))
    y3 = X(2,p(3))
    z3 = X(3,p(3))

    x4 = X(1,p(4))
    y4 = X(2,p(4))
    z4 = X(3,p(4))

    
    volume = VolumeOfTetrahedron(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)


    if (ABS(volume).le.eps) then
        write(*,*) 'Error: element ', elemID, ' is flat. Zero volume'
        write(*,*) 'All the nodes are in a single plane.'
        write(*,*) 'the .msh file contains an error'
        write(*,*) 'coordinates are: '
        write(*,*) "point 1: ", p(1), "has coords:", x1, y1, z1
        write(*,*) "point 2: ", p(2), "has coords:", x2, y2, z2
        write(*,*) "point 3: ", p(3), "has coords:", x3, y3, z3
        write(*,*) "point 4: ", p(4), "has coords:", x4, y4, z4
        write(*,*) "has volume: ", ABS(volume)
        call exitp(1)
        return
    endif

     if (volume.lt.0.0) then
        ! element has the wrong order of nodes, so that elements becomes negative.
        ! exchange two nodes of the element, and all should be well.
        write(*,*) 'Warning: element ', elemID, ' is inverted.'
        write(*,*) 'Switching vertex 1 and 2 to fix it'

        tempnode = IEN(1,elemID)
        IEN(1,elemID) = IEN(2,elemID)
        IEN(2,elemID) = tempnode

        return
    endif

enddo
!
if (iecho.eq.2) then
    write(stderr,50)
50    format(1x,'Elements tested ok')
endif

return

end subroutine
!-------------------------------------------------------------------------------
! tested 26/6/2003
 subroutine MATERL (DMAT,E,POIS)

 USE ALGEBRA
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
implicit none

integer, parameter :: NDOF=3, NSD=3, NSTR=6, NEN=4, NEE=12
!-pass
integer            :: ierr
double precision   :: XL(NSD,NEN),DL(NEE),EE(NSTR)
logical            :: QUAD
!-locl
double precision   ::  BP(NDOF,NSD)

!  Calculate partial displacement derivatives
call BPMATRIX (XL,DL,QUAD,BP,ierr)

if (ierr.ne.0) then
    return
endif

EE(1) = BP(1,1)
EE(2) = BP(2,2)
EE(3) = BP(3,3)

EE(4) = BP(1,2)+BP(2,1) ! fluid dynamical definition of shear strain
EE(5) = BP(1,3)+BP(3,1) ! fluid dynamical definition of shear strain
EE(6) = BP(3,2)+BP(2,3) ! fluid dynamical definition of shear strain

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
use modeldefinition, only: igrav

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
double precision :: xc, yc, zc
!
!    Compute Jacobian determinant
call SHPTET (ZERO,ZERO,ZERO,XL,DET,SH,IERR,1)



if (IERR.ne.0) return

    if (igrav.eq.1) then
        GRAV1 = GRAV(1)*DENS*DET*QUART
        GRAV2 = GRAV(2)*DENS*DET*QUART
        GRAV3 = GRAV(3)*DENS*DET*QUART
    else if (igrav.eq.2) then
        ! get center of element
        xc = ( XL(1,1) + XL(1,2) + XL(1,3) + XL(1,4)) * 0.25
        yc = ( XL(2,1) + XL(2,2) + XL(2,3) + XL(2,4)) * 0.25
        zc = ( XL(3,1) + XL(3,2) + XL(3,3) + XL(3,4)) * 0.25

        ! (total grav)^2 = Gx^2 + Gy^2 + Gz^2
        ! xc/yc = Gx/Gy
        ! xc/yz = Gx/Gz
        ! from these equations we compute Gx, Gy and Gz.
        
        GRAV1 = sqrt ( (GRAV(1)**2) / ( 1.0+(yc**2 / xc**2)+(zc**2 / xc**2) )  )
        GRAV2 = sqrt ( (GRAV(1)**2) / ( 1.0+(xc**2 / yc**2)+(zc**2 / yc**2) )  )
        GRAV3 = sqrt ( (GRAV(1)**2) / ( 1.0+(xc**2 / zc**2)+(yc**2 / zc**2) )  )

        ! add sign, based on whether points are left/right down/up from origin
        GRAV1 = sign(GRAV1,-xc)*DENS*DET * 0.25
        GRAV2 = sign(GRAV2,-yc)*DENS*DET * 0.25
        GRAV3 = sign(GRAV3,-zc)*DENS*DET * 0.25

    else
        ! this is not supposed to happen
    endif

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
    
return
end


!-------------------------------------------------------------------------------
! tested 26/6/2003
subroutine STIFF (DMAT,XL,SL,DUM1,QUAD,IERR, nodeID)
!
! Computes the super-diagonal part of the local stiffness matrix.
!
implicit none

integer NEN,NSD,NDOF,NSTR,NEE
parameter (NEN=4, NSD=3, NDOF=3, NSTR=6, NEE=NEN*NDOF)
!-pass
integer          :: IERR
logical          :: QUAD
double precision :: DMAT,XL,SL,DUM1
dimension        :: XL(NSD,NEN),DMAT(NSTR,NSTR),SL(NEE,NEE)
integer                :: nodeID

!-locl
double precision, parameter :: ZERO=0d0
integer          :: i, j, k
double precision :: DET,DB1,DB2,DB3,DB4,DB5,DB6
double precision :: SH(4,NEN),B(NSTR,NEE)
!
!    Compute spatial derivatives of shape functions and determinant
call SHPTET (ZERO,ZERO,ZERO,XL,DET,SH,IERR,1)


if (IERR.ne.0) then
    return
endif

! Assemble B-matrix
call BMATRIX (B,SH)

do j=1,NEE
!        Calculate column J of DB-matrix


DB1=0
DB2=0 
DB3=0 
DB4=0 
DB5=0 
DB6=0 

do k=1,6
    DB1 = DB1 + DMAT(1,k)*B(k,j)    
    DB2 = DB2 + DMAT(2,k)*B(k,j)
    DB3 = DB3 + DMAT(3,k)*B(k,j)
    DB4 = DB4 + DMAT(4,k)*B(k,j)
    DB5 = DB5 + DMAT(5,k)*B(k,j)
    DB6 = DB6 + DMAT(6,k)*B(k,j)
enddo

!        DB1=DMAT(1,1)*B(1,j)+
!     >      DMAT(1,2)*B(2,j)+
!     >      DMAT(1,3)*B(3,j)+
!     >      DMAT(1,4)*B(4,j)+
!     >      DMAT(1,5)*B(5,j)+
!     >      DMAT(1,6)*B(6,j)

!        DB2=DMAT(2,1)*B(1,j)+DMAT(2,2)*B(2,j)+DMAT(2,3)*B(3,j)
!     >         +DMAT(2,4)*B(4,j)+DMAT(2,5)*B(5,j)+DMAT(2,6)*B(6,j)

!        DB3=DMAT(3,1)*B(1,j)+DMAT(3,2)*B(2,j)+DMAT(3,3)*B(3,j)
!     >         +DMAT(3,4)*B(4,j)+DMAT(3,5)*B(5,j)+DMAT(3,6)*B(6,j)

!        DB4=DMAT(4,1)*B(1,j)+DMAT(4,2)*B(2,j)+DMAT(4,3)*B(3,j)
!     >         +DMAT(4,4)*B(4,j)+DMAT(4,5)*B(5,j)+DMAT(4,6)*B(6,j)

!        DB5=DMAT(5,1)*B(1,j)+DMAT(5,2)*B(2,j)+DMAT(5,3)*B(3,j)
!     >         +DMAT(5,4)*B(4,j)+DMAT(5,5)*B(5,j)+DMAT(5,6)*B(6,j)

!        DB6=DMAT(6,1)*B(1,j)+DMAT(6,2)*B(2,j)+DMAT(6,3)*B(3,j)
!     >         +DMAT(6,4)*B(4,j)+DMAT(6,5)*B(5,j)+DMAT(6,6)*B(6,j)


    do i=1,j
        SL(i,j)=SL(i,j) +DET*( B(1,i)*DB1 + &
                               B(2,i)*DB2 + &
                               B(3,i)*DB3 + &
                               B(4,i)*DB4 + &
                               B(5,i)*DB5 + &
                               B(6,i)*DB6)
    enddo
enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
! not tested yet
subroutine ADDPR (B,X,D,IEN,LM,DUM1,DUM2,IELNO,ISIDE,PRES, &
                  LMF,TFAULT,SKEW,DXE)

USE MATERIALSMODULE, only: LMAT
USE MODELDEFINITION
USE ALGEBRA
!use modeltopology,   only: ndof, nen, nsd
use constants,       only: eps, useDXE, one, third, half
use debugmodule,     only: iecho
use iomodule,        only: stdout
!
! Subroutine to add pressures to load vector
!
implicit none

integer, parameter :: NEN=4, NSD=3, NDOF=3

!-pass
integer            :: IEN,LM,IELNO,ISIDE,LMF
double precision   :: B,X,D,DUM1,DUM2,PRES,TFAULT,SKEW,DXE
dimension          :: IELNO(*),DUM1(*),IEN(NEN,*),X(NSD,*),D(NDOF,*), &
                      ISIDE(*),PRES(*),DUM2(*),B(*),LM(NDOF,NEN,*), &
                      LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
                      DXE(NDOF,NEN,*)
!-locl
integer            :: k,n,ii, idof
integer            :: lien(3)
double precision   :: area
double precision   :: f(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD), &
                      angle(2),tmp(NSD),rot(3,3)


do k=1,IABS(NUMPR)
!   Set element number N
    n = IELNO(k)

!   Localize coordinates
    call LCOORD (X,xl,IEN(1,n),n)

    if (LGDEF.ne.0) then
       call LDISP  (dl,D,IEN(1,n),NDOF,NEN)
       call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,35)
       call ADDSNE (dl,n,NDOF,NEN,useDXE)
       call REZONE(xl,dl,ONE)
    endif

!   Get local nodal point numbers of the side on which pressure works
    call SIDENP(ISIDE(k),.true.,lien)

!   Calculate surface area and unit normal vector
    call FCGEOM (xl,lien,area,xn)

!   Define surface mean of SKEW angles
    angle(1) = ( -SKEW(1,IEN(lien(1),n)) - &
                  SKEW(1,IEN(lien(2),n)) - &
                  SKEW(1,IEN(lien(3),n)) ) * HALF

    angle(2) = ( -SKEW(2,IEN(lien(1),n)) - &
                  SKEW(2,IEN(lien(2),n)) - &
                  SKEW(2,IEN(lien(3),n)) ) * HALF

    if (ABS(angle(1)).gt.EPS .or. &
        ABS(angle(2)).gt.EPS) then
        call FORMRT (angle,rot,NDOF)
        call EQUATE (tmp,xn,NDOF)
        call VMPRD (rot,tmp,xn,NDOF,3)
    endif

    f(1) = xn(1) * PRES(k) * area * THIRD
    f(2) = xn(2) * PRES(k) * area * THIRD
    f(3) = xn(3) * PRES(k) * area * THIRD

#ifdef GIA
	do idof = 3,3
#else
  	do idof=1,NDOF
#endif
        ii = LM(idof,lien(1),n)
        if (ii.ne.0) then
            B(ii) = B(ii) + f(idof)
        endif
        ii = LM(idof,lien(2),n)
        if (ii.ne.0) then
            B(ii) = B(ii) + f(idof)
        endif
        ii = LM(idof,lien(3),n)
        if (ii.ne.0) then
            B(ii) = B(ii) + f(idof)
        endif
    enddo
enddo
!
if (iecho.eq.1) then
    write(stdout,1)
 1   format(1x,'Pressure loads are being added')
endif

return
end
!-----------------------------------------------------------------------
subroutine SIDENP (iss,QUAD,lien)

use debugmodule,     only: xit
use iomodule,        only: stderr
!
! routine sets (local) nodal point numbers belonging to side IS.
! NPs in lien are positively oriented.
!
implicit none
!-pass
integer iss,lien
logical QUAD
dimension lien(3)
!
if (iss.eq.1) then
    lien(1) = 1
    lien(2) = 2
    lien(3) = 4
elseif (iss.eq.2) then
    lien(1) = 1
    lien(2) = 3
    lien(3) = 2
elseif (iss.eq.3) then
    lien(1) = 1
    lien(2) = 4
    lien(3) = 3
elseif (iss.eq.4) then
    lien(1) = 2
    lien(2) = 3
    lien(3) = 4
else
    write(stderr,1) iss
 1        format(1x,'SIDENP: a tetrahedron has 4 sides, not ',I2)
    call xit(1," ")
endif
!
return
end
!-----------------------------------------------------------------------
subroutine NPSIDE (NP,IEN,ISINT)

! returns side number 
! face 1: nodes 1 2 4
!      2:       1 3 2
!      3:       1 4 3
!      4:       2 3 4

integer, parameter:: NEN=4
integer :: NP(3),IEN(NEN),ISINT
! locl
logical :: N1_in_IEN=.false., N2_in_IEN=.false., N3_in_IEN=.false.
integer :: i

i=1
find_np1_in_IEN: do while (i.le.NEN)
    N1_in_IEN = NP(1).eq.IEN(i)
    if (N1_in_IEN) exit
    i = i + 1
enddo find_np1_in_IEN
if (.not.N1_in_IEN) then
    ISINT = 4
    return
endif

i=1
find_np2_in_IEN: do while (i.le.NEN)
    N2_in_IEN = NP(2).eq.IEN(i)
    if (N2_in_IEN) exit
    i = i + 1
enddo find_np2_in_IEN
if (.not.N2_in_IEN) then
    ISINT = 3
    return
endif

i=1
find_np3_in_IEN: do while (i.le.NEN)
    N3_in_IEN = NP(3).eq.IEN(i)
    if (N3_in_IEN) exit
    i = i + 1
enddo find_np3_in_IEN

if (.not.N3_in_IEN) then
    ISINT = 1
else
    ISINT = 2
endif

return

end subroutine
!-----------------------------------------------------------------------
! tested 26/5/2003
 subroutine FCGEOM (XL,LIEN,AREA,XN)

 use constants, only: dfmin
use iomodule, only: stderr
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
 GRAD(1)=TL(1)*SH(1,1) + &
    TL(2)*SH(1,2) + &
    TL(3)*SH(1,3) + &
    TL(4)*SH(1,4)

 GRAD(2)=TL(1)*SH(2,1) + &
    TL(2)*SH(2,2) + &
    TL(3)*SH(2,3) + &
    TL(4)*SH(2,4)

 GRAD(3)=TL(1)*SH(3,1) + &
    TL(2)*SH(3,2) + &
    TL(3)*SH(3,3) + &
    TL(4)*SH(3,4)

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
PL(i) = PL(i) + DET*(B(1,i)*EVP(1)+ &
                     B(2,i)*EVP(2)+ &
                     B(3,i)*EVP(3)+ &
                     B(4,i)*EVP(4)+ &
                     B(5,i)*EVP(5)+ &
                     B(6,i)*EVP(6))
 enddo
!
 return
 end

!-------------------------------------------------------------------------------
! no need to test
subroutine DEVSTN (STN,SDEV)

use constants, only: two, third
!
! Computes stress deviatior SDEV
!
implicit none
integer, parameter :: NSTR=6
!-pass
double precision   :: STN(NSTR),SDEV(NSTR)
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
subroutine FORMBT (STN,BETA,BETB,PRPMAT,PRPLAS, elemID)

USE MODELCTX
use constants, only: dfmin,dlogmin,dfmax,dlogmax, zero, third, one, two, six
use algebra,   only: idig
use modeldefinition, only: icvis
!
! Subroutine to compute strain rates from stress field. It defines BETA and
! evaluates the product BETA*STRESS. BETB is the viscoplastic strain rate.
!
implicit none

integer, parameter :: NSTR=6
!-pass
double precision   :: STN(NSTR),BETA(NSTR),BETB(NSTR),PRPMAT(6),PRPLAS(9)
integer            :: elemID
!-locl
double precision   ::  AVECT(NSTR)
double precision   :: sigma, emhu, anpwr
double precision   :: streff
double precision   :: yield

!-init
sigma = ZERO
emhu  = PRPMAT(3)
anpwr = PRPMAT(4)
!
if (IDIG(ICVIS,4,2).ne.0) then
!   Calculate yield stress YIELD
    call YIELDS (STN,PRPLAS,YIELD, elemID)
!   Construct plastic flow vector AVECT
    call YIELDF (STN,PRPLAS,AVECT, elemID)
!   Evaluate viscoplastic strain rate matrix BETB
    call FLOWVP (AVECT,STN,PRPLAS,YIELD,BETB, elemID)
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
!   add visco-plastic contributions
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
! no need to test  ... WHAT? Why TF not?
subroutine FORMJB (BETAJ,STN,PRPMAT,PRPLAS, elemID)

USE MODELCTX
USE ALGEBRA
use constants, only: dfmin,dlogmin,dfmax,dlogmax, one, two, third, six
use MODELDEFINITION, only: icvis
!
! Constructs the Jacobean of the stress-strain rate matrix BETA.
!
implicit none

integer, parameter :: NSTR=6
!-pass
double precision :: BETAJ(6,6),STN(NSTR),PRPLAS(9),PRPMAT(6)
integer          :: elemID
!-locl
double precision :: BETAJP(NSTR,NSTR)
integer          :: ICVISentry
double precision :: yield
double precision :: emhu, anpwr
double precision :: a1,a2,a3
double precision :: sxx,syy,sxy,szz,sxz,syz
double precision :: strinv
double precision :: streff
double precision :: sigma
integer          :: j
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
ICVISentry = IDIG(ICVIS,4,2)
!       Plastic flow?
if (ICVISentry.ne.0) then
    !        Calculate yield stress YIELD
    call YIELDS (STN,PRPLAS,YIELD, elemID)
!!!!!!!
!!! do we not need a call to yieldf here?
!!!!!!!!!
    !        Construct Jacobian due to plastic yielding BETAJP
    call CLEAR (BETAJP,36,"BETAJP")
    call TANGVP (STN,YIELD,PRPLAS,BETAJP, elemID)
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
if (ICVISentry.ne.0) then
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
USE ALGEBRA
use constants,       only: eps, useDXE
use debugmodule,     only: iecho
use iomodule
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact

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
do n=1,NUMSTR
    k = IELSTR(n)

    if (k.le.0) then
        write(stderr,1) n
 1        format(1x,'ADDSTR: fatal IELSTR error. Index=',I5)
        call stoper()
    endif

    it0 = ISTR(1,n)
    it1 = ISTR(2,n)

    if (DBLE(NSTEP-it0)*DBLE(NSTEP-it1).gt.0d0) then
        cycle
    endif

    numadd = numadd + 1
    call LCOORD (X,xl,IEN(1,k),0)

    if (LGDEF.ne.0) then
        call LDISP  (dl,D,IEN(1,k),NDOF,NEN)
        call AddFaultDisplacement (dl,k,TFAULT,NDOF,NEN,36)
        call ADDSNE (dl,k,NDOF,NEN,useDXE)
        call REZONE (xl,dl,1d0)
    endif

!       Get nodal point numbers of the side on which pressure works
    call SIDENP(ISSIDE(n),.true.,lien)

!       Calculate surface area and normal vector
    call FCGEOM (xl,lien,area,xn)

!       Rotate face into local direction
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

enddo
!
 if (iecho.eq.1 .and. numadd.gt.0) write(stdout,2)
 2    format(1x,'Global stress loads are being added')
!
 return
 end subroutine
!-------------------------------------------------------------------------------
 subroutine ADDWNK (IWELM,IWSIDE,IWTIME,WPRES,IEN,MAT,X,D,PRPMAT, &
   A,LM,IDIAG,LMF,TFAULT,SKEW,DXE,LMX,NSLIP)

USE MATERIALSMODULE, only : LMAT
USE ALGEBRA
use constants,       only : eps, quart, third, half, one, useDXE
!use meshdatamodule,  only : meshdatactx
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only : elementSlipperyID
USE MODELDEFINITION, only : numwnk, numslp, numrot, lgdef, neq, nElemsWithSlipperyNodes
#else
use meshdatamodule,  only : meshdatactx
USE MODELDEFINITION, only : numwnk, numslp, numrot, lgdef, neq
#endif
use modeltopology,   only : nsd, ndof, nen
use petscksp

#ifdef SPARSE
use spetscmodule,      only: m_a
#endif
use iomodule
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact

!
! Subroutine to add Winkler pressures to diagonal elements of stiffness
! matrix
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


integer          :: IWELM, IWSIDE, IWTIME, ien, idiag(neq), lm, lmf, nslip, mat
double precision :: wpres, x, prpmat, a, d, tfault, skew

dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
          MAT(*),X(NSD,*),PRPMAT(6,*),A(*),LM(NDOF,NEN,*), &
          D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
          NSLIP(5,*)

#ifdef SPARSE
#include "petsc/finclude/petscsys.h"
#endif
!-locl
double precision :: w(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD),rot(3,3), &
                    tmp(NSD),angle(2)
integer          :: lien(3)
double precision :: area
integer          :: ni, nj, nk, nl, ii, idof, li, ll, nn, np
integer          :: globalelementnumber, globalnodenumber, nodecounter
double precision :: v
integer          :: ierr

!init
ni = 0
nj = 0
nk = 0

if (NUMWNK.eq.0) then
    return
endif

do nn=1,IABS(NUMWNK)
    globalelementnumber = IABS(IWELM(nn))
    LL= IWTIME(nn)
    if (LL.lt.0) then
        if (LL.ne.-1 .and. NSTEP.ge.-LL) then
            cycle
        endif
    elseif (LL.eq.0) then
        cycle
    else
        if (NSTEP.lt.(LL-1)) then
            cycle
        endif
    endif

    call LCOORD (X,XL,IEN(1,globalelementnumber),0)
    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,globalelementnumber),NDOF,NEN)
        call AddFaultDisplacement (DL,globalelementnumber,TFAULT,NDOF,NEN,37)
        call ADDSNE (DL,globalelementnumber,NDOF,NEN,useDXE)
        call REZONE (XL,DL,ONE)
    endif
!           Get nodal point numbers of the side on which pressure works
    call SIDENP (IABS(IWSIDE(nn)),.true.,lien)
!           Calculate surface area and unit normal vector
    call FCGEOM (XL,lien,area,xn)

    if (NUMROT.ne.0) then
!               Rotate face normal into local direction
        np = IEN(lien(1),globalelementnumber)
        angle(1) = -SKEW(1,np)
        angle(2) = -SKEW(2,np)
        np = IEN(lien(2),globalelementnumber)
        angle(1) = angle(1)-SKEW(1,np)
        angle(2) = angle(2)-SKEW(2,np)
        np = IEN(lien(3),globalelementnumber)
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

!    w(1) = ABS( xn(1) * area ) * WPRES(nn) * THIRD
!    w(2) = ABS( xn(2) * area ) * WPRES(nn) * THIRD
!    w(3) = ABS( xn(3) * area ) * WPRES(nn) * THIRD

!    w(1) = area * WPRES(nn) * THIRD
!    w(2) = area * WPRES(nn) * THIRD
    w(3) = area * WPRES(nn) * THIRD



    if (NUMSLP.gt.0) then
        if (IWELM(nn).gt.0) then
            ni = 0
            nj = 0
            nk = 0
            nl = 0
            do ii = 1,NUMSLP
                if (NSLIP(2,ii).eq.IEN(lien(1),globalelementnumber)) then
                    ni = ni + 1
                endif
                if (NSLIP(2,ii).eq.IEN(lien(2),globalelementnumber)) then
                    nj = nj + 1
                endif
                if (NSLIP(2,ii).eq.IEN(lien(3),globalelementnumber)) then
                    nk = nk + 1
                endif
            enddo
        else
            ni = 1
            nj = 1
            nk = 1
        endif
    endif
!        add contributions to global stiffness matrix
#ifdef GIA
    do idof=3,3  ! assuming only the z direction has a non-zero Winkler contribution.
                 ! if the direction is different, Euler angles rotate the nodal point.
#else
    do idof=1,NDOF
#endif

        if (ABS(w(idof)).lt.EPS) then
            cycle
        endif
        do nodecounter = 1, 3
            globalnodenumber = IEN(lien(nodecounter),globalelementnumber)

!           write(*,*) "applyling to node ", nodecounter, globalnodenumber


            ii = LM(idof,lien(nodecounter),globalelementnumber)
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
        if (NUMSLP.eq.0) then
            cycle
        endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
if (elementSlipperyID(n,40).gt.0) then
    ii = IABS(LMX(idof,lien(1),elementSlipperyID(globalelementnumber,40)))
else
    ii = 0
endif
#else 
ii = IABS(LMX(idof,lien(1),globalelementnumber))
#endif


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

#ifdef EXPERIMENTAL_ReducedSlipAlloc
if (elementSlipperyID(globalelementnumber,41).gt.0) then
    ii = IABS(LMX(idof,lien(2),elementSlipperyID(globalelementnumber,41)))
else
    ii = 0
endif
#else 
ii = IABS(LMX(idof,lien(2),globalelementnumber))
#endif


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


#ifdef EXPERIMENTAL_ReducedSlipAlloc
if (elementSlipperyID(globalelementnumber,42).gt.0) then
    ii = IABS(LMX(idof,lien(3),elementSlipperyID(globalelementnumber,42)))
else
    ii = 0
endif
#else
ii = IABS(LMX(idof,lien(3),globalelementnumber))
#endif

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
    enddo
enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
! not tested yet
subroutine ELCTR (XELM,IEN,X,D,elemID,TFAULT, &
                  SKEW,FACTOR,IERR)

USE MODELDEFINITION
use modeltopology,   only: nen, ndof, nsd
use modeldatamodule, only: modeldatactx
use constants,       only: useDXE, quart
!
! Calculates (deformed) element center coordinates XELM
!
implicit none

!-pass
 integer IEN,LMF,IERR
 integer :: elemID
 double precision XELM,X,D,TFAULT,SKEW,FACTOR
 dimension XELM(NSD),IEN(NEN),X(NSD,*),D(NDOF,*), &
   TFAULT(NDOF,*), SKEW(2,*)
!-locl
 double precision xl,dl
 dimension xl(NSD,NEN),dl(NDOF,NEN)

ierr=0
!
!    localize coordinates
call LCOORD (X,xl,IEN,elemID)

if (FACTOR.gt.0.) then
    call LDISP  (dl,D,IEN,NDOF,NEN)
    call AddFaultDisplacement (dl,elemID,TFAULT,NDOF,NEN,38)
    call ADDSNE (DL,elemID,NDOF,NEN,useDXE)
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
    SL(i,j) = SL(i,j) + DET*(C(1)*SH(1,i)*SH(1,j) + &
                             C(2)*SH(2,i)*SH(2,j) + &
                             C(3)*SH(3,i)*SH(3,j))
enddo
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
! not tested yet
 subroutine ADVSTF (S,XL,VL,PL,COND,RHOCP,TH,HEAT,QUAD,IERR)
 USE ALGEBRA
 use debugmodule, only: iecho
use iomodule
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
                 (vy+VL(2,i))*sh(2,j) + &
                 (vz+VL(3,i))*sh(3,j) )
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
        S(i,j)=S(i,j)+w*(vx*sh(1,j)+ &
                         vy*sh(2,j)+ &
                         vz*sh(3,j))
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
double precision function EigMax (STN,E,V,EMHU,NPL,ALPHA)
!
! Gives an upper limit estimate of the maximum local eigenvalue
!
implicit none

integer, parameter :: NSTR=6
!-pass
double precision :: E, V, EMHU, NPL, ALPHA
double precision :: STN(NSTR)

!
EigMax = 1.d0
!
return
end
!-------------------------------------------------------------------------------
#ifdef EXPERIMENTALTRANSIENTDISCREEP
subroutine EffectiveDislCreepRate(InvEffViscDisl,STN,EffectiveRate)
! computes beta vector (viscous strain rate) for dislocation creep only
! return effective beta, i.e. square root of second invariant

use modeltopology, only: nstr
use constants, only: two, six, third, half

implicit none

!-pass 
double precision :: STN(nstr) ! stress vector
double precision :: InvEffViscDisl ! inverse viscosity for dislocation creep
!-local
double precision :: StrainTensor(nstr) ! viscous strain rate vector
double precision :: sigma,EffectiveRate

! inverse viscosity divided by 6
sigma = InvEffViscDisl/(SIX)

!    viscous stress-strain rate vector,
!  where StrainTensor(4:6) are tensor strain rate and not engineering strain rates
StrainTensor(1) = sigma * (TWO*STN(1)-STN(2)-STN(3))
StrainTensor(2) = sigma * (TWO*STN(2)-STN(1)-STN(3))
StrainTensor(3) = sigma * (TWO*STN(3)-STN(1)-STN(2))
StrainTensor(4) = sigma * STN(4) * SIX * half ! not engineering shear strain rates
StrainTensor(5) = sigma * STN(5) * SIX * half
StrainTensor(6) = sigma * STN(6) * SIX * half


! equation from Ranalli - Rheology of the Earth, eq 4.20
EffectiveRate = sqrt(half*(StrainTensor(1)*StrainTensor(1) &
                          + StrainTensor(2)*StrainTensor(2) &
                          + StrainTensor(3)*StrainTensor(3) ) &
                          + StrainTensor(4)*StrainTensor(4) &
                          + StrainTensor(5)*StrainTensor(5) &
                          + StrainTensor(6)*StrainTensor(6) )


return
end
#endif
!-------------------------------------------------------------------------------
subroutine ADDFLX (BT,IEN,X,D,DUM1,DUM2,IFLS,IFLX,BFLX,LMT, &
                   LMTX,LMF,TFAULT,SKEW,DXE)!,NFLX)

USE MATERIALSMODULE, only: LMAT
USE MODELDEFINITION, only: lgdef, nflx, numfn, numslp, idift
use constants,       only: third, useDXE
!
! Routine to add normal boundary fluxes to load vector
!
implicit none
integer, parameter :: NEN=4, NSD=3, NDOF=3
!-pass
!    integer NFLX

integer          :: ien, ifls, iflx, lmt, lmf, lmtx
double precision :: BT, x, dum1, dum2, bflx, d, tfault, skew, dxe

 dimension BT(*),IEN(NEN,*),X(NSD,*),DUM1(*),DUM2(*),IFLS(*), &
 IFLX(*),BFLX(*),LMT(NEN,*),D(NDOF,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),LMTX(NEN,*),SKEW(2,*),DXE(NDOF,NEN,*)
!-locl
double precision   :: XL(NSD,NEN),DL(NDOF,NEN),XN(NSD)
integer            :: lien(3)
integer            :: ii, n, m
double precision   :: F, area
integer            :: iside
!
if (NFLX.le.0) then
    return
endif

do n=1,NFLX
    m = IFLX(n)
    iside = IFLS(n)

    call LCOORD (X,XL,IEN(1,m),0)

    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,m),NDOF,NEN)
        call AddFaultDisplacement (DL,m,TFAULT,NDOF,NEN,39)
        call ADDSNE (DL,m,NDOF,NEN, useDXE)
        call REZONE (XL,DL,1d0)
    endif

    call SIDENP (iside,.true.,lien)
    call FCGEOM (XL,lien,AREA,XN)

    F=-AREA*BFLX(n)*THIRD
    ii = LMT(lien(1),m)

    if (ii.eq.0.and.NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
           ii = LMTX(lien(1),m)
    endif

    if (ii.gt.0) then
        BT(ii) = BT(ii) + F
    endif

    ii = LMT(lien(2),m)

    if (ii.eq.0.and.NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
       ii = LMTX(lien(2),m)
    endif

    if (ii.gt.0) then
        BT(ii) = BT(ii) + F
    endif

    ii = LMT(lien(3),m)

    if (ii.eq.0.and.NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
       ii = LMTX(lien(3),m)
    endif

    if (ii.gt.0) then
        BT(ii) = BT(ii) + F
    endif

enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine TEMWNK (AA,LMT,ITDIAG,ITWINK,TWINK,IEN,MAT,PRPMAT, &
                   X,D,LMF,TFAULT,SKEW,DXE)

USE MATERIALSMODULE, only: LMAT
USE MODELDEFINITION, only: ntwink, lgdef
use constants,       only: one, useDXE
use debugmodule,     only: debug, iecho
#ifdef SPARSE
use spetscmodule,      only: t_a
#endif
use petscksp,        only: add_values
use iomodule
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact

!
! Subroutine to add Winkler fluxes to diagonal elements of thermal stiffness
! matrix
!
implicit none
integer, parameter :: NSD=3,NDOF=3,NEN=4
!-pass

integer :: itwink, ien, itdiag, lmt, lmf
double precision :: twink, mat, x, prpmat, aa, d, tfault, skew, dxe

 dimension ITWINK(3,*),TWINK(*),IEN(NEN,*), &
   MAT(*),X(NSD,*),PRPMAT(6,*),AA(*),ITDIAG(*),LMT(NEN,*), &
   D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
   DXE(NDOF,NEN,*)
!-locl
double precision :: xl(NSD,NEN),dl(NDOF,NEN),xn(NSD)
integer          :: lien(3)
integer          :: n, nn, LL, ii
double precision :: wt, area
integer          :: ierr
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
        if (LL.ne.-1 .and. NSTEP.ge.-LL) then
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
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,40)
        call ADDSNE (DL,n,NDOF,NEN, useDXE)
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
    ! SPARSE
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
! SPARSE
    endif
enddo
!
 return
 end
!-----------------------------------------------------------------------
! not tested
function Tcentr (T,TX,XL,IEN,LMTX,IDUM1)

USE MODELDEFINITION, only: numslp, numfn, idift
use constants, only: quart
!
! Calculates the temperature in the center of the element
!
implicit none

integer, parameter :: NEN=4, NSD=3
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
if (NUMSLP+NUMFN.gt.0.and.IDIFT.eq.1) then
    call LTEMPX (LMTX,IEN,TX,TL,NEN)
endif
!
!    Calculate center temperature
Tcentr = (TL(1)+TL(2)+TL(3)+TL(4))*QUART
!
 return
 end
!-----------------------------------------------------------------------
! not tested
subroutine CPSTIFF (XL,QUAD,S,RHOCP,DUM1,IERR)

use constants, only: zero, two
!
! computes the upper local heat capacity stiffness matrix
!
implicit none

integer, parameter :: NEN=4,NSD=3
!-pass
integer            :: IERR
logical            :: QUAD
double precision   :: RHOCP,DUM1
double precision   :: XL(NSD,NEN),S(NEN,NEN)
!-locl
double precision det,sh,w
dimension sh(4,NEN)
!
call SHPTET (ZERO,ZERO,ZERO,XL,det,sh,IERR,1)

if (IERR.ne.0) then
    return
endif

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

return
end subroutine
!-------------------------------------------------------------------------------
double precision function AJ1 (ST)
!
! Calculates first invariant of stress or strain matrix
!
implicit none

!-pass
double precision :: ST(3)

AJ1 = ST(1)+ST(2)+ST(3)

return
end function
!-------------------------------------------------------------------------------
double precision function AJ2 (ST)
!
! Calculates second deviatoric invariant of matrix ST
!
use constants, only: six

implicit none
!-pass
double precision :: ST(6)
double precision :: SXY, SXZ, SYZ
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
double precision function AJ3 (ST)
!
! Calculates third deviatoric invariant of matrix ST
!
use constants, only: third, two

implicit none
!-pass
double precision :: ST(6)
!-local
double precision :: SDEV1, SDEV2, SDEV3


SDEV1 = (TWO*ST(1)-ST(2)-ST(3))*THIRD
SDEV2 = (TWO*ST(2)-ST(1)-ST(3))*THIRD
SDEV3 = (TWO*ST(3)-ST(1)-ST(2))*THIRD

AJ3 = SDEV1*SDEV2*SDEV3 - SDEV1*ST(6)*ST(6) - &
                          SDEV2*ST(5)*ST(5) - &
                          SDEV3*ST(4)*ST(4) + &
                            TWO*ST(4)*ST(5)*ST(6)

return
end function
!-------------------------------------------------------------------------------
subroutine REZONE (XL,DL,FACTOR)
!
! Routine to update local coordinate XL by displacements DL
!
implicit none

integer, parameter :: NSD=3
integer, parameter :: NDOF=3
integer, parameter :: NEN=4

!-pass
double precision :: XL(NSD,NEN),DL(NDOF,NEN)
double precision :: factor
!
integer :: i, j

do i=1,3
    do j=1,4
        XL(i,j) = XL(i,j) + DL(i,j) * factor
    enddo
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine TSTRN (EE,DT,VEXP)
!
! calculates strain due to thermal expansion
!
use constants, only: zero, third

implicit none

integer, parameter :: NSTR=6
!-pass
double precision   :: EE(NSTR)
double precision   :: DT, VEXP
!
EE(1) = DT*VEXP*THIRD
EE(2) = EE(1)
EE(3) = EE(1)
EE(4) = ZERO
EE(5) = ZERO
EE(6) = ZERO
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine YIELDS (STN,PRPLAS,YIELD, elemID)

USE CONSTANTS, only: third, two, root3, deg2rad, dfmin,dlogmin,dfmax,dlogmax
use iomodule,  only: stderr
!
! Routine to calculate the yield stress for various types of plasticity.
!
implicit none

integer, parameter :: NSTR=6
!-pass
double precision :: STN(NSTR),PRPLAS(9)
double precision :: yield
integer          :: elemID
!-locl
double precision, external :: AJ1,AJ2
double precision :: frangl, sinphi, stmean
integer          :: nptype
double precision :: stnEff
!
STNEFF = AJ2(STN)

if (STNEFF.gt.dfmin) then
    STNEFF = SQRT(STNEFF)
endif

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
end subroutine
!-------------------------------------------------------------------------------
subroutine YIELDF (STN,PRPLAS,AVECT, elemID)

USE ALGEBRA, only: clear
USE CONSTANTS, only: zero, one, two, three, third, &
    eight, root3, deg2rad, THF, dfmin,dlogmin,dfmax,dlogmax

!
! Evaluate the plastic flow vector AVECT.(Owen & Hinton, sec. 7.4)
!
implicit none

integer, parameter :: NSTR=6
!-pass
double precision :: STN(NSTR),PRPLAS(9),AVECT(NSTR)
integer          :: elemID
!-locl
double precision :: SDEV(NSTR),A1(NSTR),A2(NSTR),A3(NSTR)
double precision, external :: AJ2,AJ3
double precision :: DJ2
double precision :: sin3th, theta, sinth, costh,tanth,cos3th,tan3th,frangl,sinphi
double precision :: c1,c2,c3
double precision :: twoj2i

integer :: nptype


save A1
!-init
data A1/1.,1.,1.,0.,0.,0./
!
!    Compute second deviatoric invariant of stress matrix
DJ2 = AJ2(STN)

if (DJ2.gt.dfmin) then
    goto 100
endif
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
    if (SIN3TH.lt.-ONE) then
        SIN3TH=-ONE
    endif

    if (SIN3TH.gt. ONE) then
        SIN3TH= ONE
    endif

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
    if (SIN3TH.lt.-ONE) then
        SIN3TH=-ONE
    endif
    if (SIN3TH.gt. ONE) then
        SIN3TH= ONE
    endif

    THETA  = THIRD*ASIN(SIN3TH)  ! <- huh?
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
else
    write(*,*) "f3dlib; subroutine YIELDF says:"
    write(*,*) "Unknown NPTYPE", NPTYPE, "Must be 1,2,3 or 4"
    write(*,*) "Could not determine c[123]."
    stop "please contact model support"
endif
!
if (NPTYPE.eq.2 .or. NPTYPE.eq.3 .or. NPTYPE.eq.4) then
!    calculate A-vector
    AVECT(1)=C1*A1(1)+C2*A2(1)+C3*A3(1)
    AVECT(2)=C1*A1(2)+C2*A2(2)+C3*A3(2)
    AVECT(3)=C1*A1(3)+C2*A2(3)+C3*A3(3)
    AVECT(4)=C1*A1(4)+C2*A2(4)+C3*A3(4)
    AVECT(5)=C1*A1(5)+C2*A2(5)+C3*A3(5)
    AVECT(6)=C1*A1(6)+C2*A2(6)+C3*A3(6)
else
    write(*,*) "f3dlib; subroutine YIELDF says:"
    write(*,*) "Unknown NPTYPE", NPTYPE, "Must be 1,2,3 or 4"
    write(*,*) "Could not set AVECT"
    stop "please contact model support"
endif
!
 return
 end
!-------------------------------------------------------------------------------
subroutine TANGVP (STN,YIELD,PRPLAS,BETAJP, elemID)

USE ALGEBRA
USE CONSTANTS, only: one, two, four, threeq, root3, dfmin,dlogmin,dfmax,dlogmax,eps
use iomodule, only: stderr
!
! Compute the Jacobian matrix for viscoplastic flow BETAJP
! (Ref.: Owen & Hinton, sec. 8.7.2)
!
implicit none

integer, parameter :: NSTR=6
!-pass
double precision :: yield
double precision :: STN(NSTR),BETAJP(NSTR,NSTR),PRPLAS(9)
integer          :: elemID
!-locl
double precision :: SDEV(NSTR),AM1(NSTR,NSTR),AM2(NSTR,NSTR)

double precision :: allow
integer          :: nptype
double precision :: COHESN, FLUIDY, YLDPWR
integer          :: NFLOW
double precision :: fnorm
double precision :: xj2
double precision :: phi, gradphi, stneff
double precision :: fact1, fact2
integer          :: j

save AM1
double precision, external :: AJ2
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

if (FNORM.lt.ALLOW) then
    return
endif
!
xJ2 = AJ2(STN)

if (xJ2.le.dfmin) then
    return
endif
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
do j=1,NSTR
    BETAJP(1,j) = FACT1*AM1(1,j)+FACT2*AM2(1,j)
    BETAJP(2,j) = FACT1*AM1(2,j)+FACT2*AM2(2,j)
    BETAJP(3,j) = FACT1*AM1(3,j)+FACT2*AM2(3,j)
    BETAJP(4,j) = FACT1*AM1(4,j)+FACT2*AM2(4,j)
    BETAJP(5,j) = FACT1*AM1(5,j)+FACT2*AM2(5,j)
    BETAJP(6,j) = FACT1*AM1(6,j)+FACT2*AM2(6,j)
enddo
!
 return
 end
!-------------------------------------------------------------------------------
subroutine INYLD (STN,PRPLAS,YIELD,BETB,FNORM)

USE CONSTANTS, only: half, one, three, six, root3, deg2rad, dfmin,dlogmin,dfmax,dlogmax

!
! determines FNORM
!
implicit none

integer, parameter :: NSTR=6
!-pass
double precision :: STN(NSTR),PRPLAS(9),BETB(NSTR)
double precision :: yield
!-locl
double precision, external :: AJ2
double precision :: allow
integer          :: nptype
double precision :: frangl, cohesn, stnhrd, fnorm, debar, stneff
!-init
ALLOW  = 0.01d0
NPTYPE = INT(PRPLAS(1))
FRANGL = PRPLAS(2)*deg2rad
COHESN = PRPLAS(3)
STNHRD = PRPLAS(7)
FNORM  = -ONE
!
DEBAR = half*AJ2(BETB)

if (DEBAR.gt.dfmin) then
    DEBAR = SQRT(DEBAR)
endif
!
STNEFF = AJ2(STN)

if (STNEFF.lt.dfmin) then
    return
endif

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

if (STNHRD.gt.dfmin) then
    COHESN = COHESN+STNHRD*DEBAR
endif

if (COHESN.lt.0.001d0) then
    COHESN = ONE
endif
!
FNORM = (YIELD-COHESN)/COHESN

if (FNORM.lt.ALLOW) then
    FNORM = -1.0d0
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine FLOWVP (AVECT,STN,PRPLAS,YIELD,BETB, elemID)

USE ALGEBRA, only: clear
use constants, only: eps, one
!
! Compute the viscoplastic strain rate matrix BETB
! (Ref.: Owen & Hinton, sec. 8.9)
!
implicit none

integer, parameter :: NSTR=6
!-pass
double precision :: AVECT(NSTR),STN(NSTR),PRPLAS(9),BETB(NSTR)
double precision :: yield
integer          :: elemID
!-init
double precision :: fluidy, yldpwr
integer          :: nflow
double precision :: fnorm, phi

FLUIDY = PRPLAS(4)
YLDPWR = PRPLAS(5)
NFLOW  = INT(PRPLAS(6))
!
!    calculate FNORM, the normalized stress overshoot (uses BETB)
call INYLD (STN,PRPLAS,YIELD,BETB,FNORM)
call CLEAR (BETB,NSTR,"BETB")
!
if (FNORM.le.0d0) then
    return
endif
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
end subroutine
!-------------------------------------------------------------------------------
subroutine WINKLP (BTOT,DELD,IWELM,IWSIDE,IWTIME,WPRES,IEN, &
                   DUM1,X,D,DUM2,LM,LMF,TFAULT,SKEW, &
                   DXE,LMX,NSLIP,DX)

USE MATERIALSMODULE, only: LMAT
USE MODELDEFINITION, only: numwnk, lgdef, numrot, numslp
USE ALGEBRA
use constants,       only: eps, third, useDXE
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact

#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: elementSlipperyID
#endif

!
! Subroutine to compute Winkler restoring pressures from displacements
! and add them to the global force vector BTOT

implicit none

integer, parameter :: NSD=3,NDOF=3,NEN=4
!-pass
integer :: iwelm, iwside, iwtime, ien, dum1, lm, lmf, lmx, nslip
double precision :: wpres, x, dum2, btot, deld, d, tfault, skew, dxe, dx


 dimension IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*),IEN(NEN,*), &
   DUM1(*),X(NSD,*),DUM2(6,*),LM(NDOF,NEN,*),BTOT(*), &
   DELD(NDOF,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*),DXE(NDOF,NEN,*),LMX(NDOF,NEN,*),NSLIP(5,*), &
 DX(NDOF,*)

!-locl
double precision :: w(NDOF),xl(NSD,NEN),dl(NDOF,NEN),xn(NSD), tmp(NSD),rot(3,3),angle(2)
integer :: lien(3)
double precision :: area
integer :: i,j,k, ni,nj,nk, np, idof, ii, ll, n, nn
integer :: slipElemID

ni = 0
nj = 0
nk = 0
!
if (NUMWNK.eq.0) then
    return
endif
!
do nn=1,IABS(NUMWNK)
    n = IABS(IWELM(nn))
    LL= IWTIME(nn)
    if (LL.lt.0) then
        if (LL.ne.-1 .and. NSTEP.ge.-LL) then
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
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,41)
        call ADDSNE (DL,n,NDOF,NEN,useDXE)
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

    do idof=1,NDOF

        if (ABS(w(idof)).lt.EPS) then
            cycle
        endif

        ii = LM(idof,i,n)

        if (ii.gt.0) then
            BTOT(ii)=BTOT(ii)- w(idof)*DELD(idof,IEN(i,n))
        endif

        ii = LM(idof,j,n)

        if (ii.gt.0) then
            BTOT(ii)=BTOT(ii)- w(idof)*DELD(idof,IEN(j,n))
        endif

        ii = LM(idof,k,n)

           if (ii.gt.0) then
            BTOT(ii)=BTOT(ii)- w(idof)*DELD(idof,IEN(k,n))
         endif

        if (NUMSLP.eq.0) then
            cycle
        endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc

        slipElemID = elementSlipperyID(n,43)

        if(slipElemID.gt.0) then
            ii = IABS(LMX(idof,i,slipElemID))

            if (ii.gt.0) then
                BTOT(ii) = BTOT(ii) - w(idof)*DX(idof,IEN(i,slipElemID))/DBLE(ni)
            endif
    
            ii = IABS(LMX(idof,j,slipElemID))

            if (ii.gt.0) then
                BTOT(ii) = BTOT(ii) - w(idof)*DX(idof,IEN(j,slipElemID))/DBLE(nj)
            endif

            ii = IABS(LMX(idof,k,slipElemID))

            if (ii.gt.0) then
                BTOT(ii) = BTOT(ii) - w(idof)*DX(idof,IEN(k,slipElemID))/DBLE(nk)
            endif
        else
            ! this element contains no slippery nodes. skip it.
        endif
#else

        ii = IABS(LMX(idof,i,n))

          if (ii.gt.0) then
            BTOT(ii) = BTOT(ii) - w(idof)*DX(idof,IEN(i,n))/DBLE(ni)
        endif

        ii = IABS(LMX(idof,j,n))

          if (ii.gt.0) then
            BTOT(ii) = BTOT(ii) - w(idof)*DX(idof,IEN(j,n))/DBLE(nj)
        endif

        ii = IABS(LMX(idof,k,n))

          if (ii.gt.0) then
            BTOT(ii) = BTOT(ii) - w(idof)*DX(idof,IEN(k,n))/DBLE(nk)
        endif

#endif
    enddo
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine SHAREL (IEN,IJSUR,NUMEL,NEL)

use debugmodule,     only: xit
use iomodule,        only: stderr
!
!     Finds elements neighbouring to element NEL.
!    Neighbouring elements are identified by three nodes.
!
implicit none

integer, parameter :: NSUR=4,NEN=4
!-pass
integer            :: NEL,NUMEL
integer            :: IEN(NEN,*),IJSUR(NSUR,*)

!-locl
integer            :: isur, iel, ii, n, nshare
integer, external  :: NINDX
!
!    set to last IJSUR(isur,NEL) <> 0 (last index that was written)
isur = NINDX (IJSUR,NEL)
!    loop other elements
do iel=NEL+1,NUMEL
!        compare nodal points of element NEL with those of IEL > NEL
    nshare = 0

    do n=1,NEN
        if (IEN(n,NEL).eq.IEN(1,iel)) nshare = nshare+1
        if (IEN(n,NEL).eq.IEN(2,iel)) nshare = nshare+1
        if (IEN(n,NEL).eq.IEN(3,iel)) nshare = nshare+1
        if (IEN(n,NEL).eq.IEN(4,iel)) nshare = nshare+1
    enddo

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

        if (isur.eq.NSUR) then
            goto 300
        endif

    endif
enddo
!
300    return
end subroutine

!-------------------------------------------------------------------------------

integer function NINDX (IJSUR,NEL)

use debugmodule,     only: xit
use iomodule,        only: stderr
!
implicit none

integer, parameter :: NSUR=4
!-pass
integer :: NEL
integer :: IJSUR
dimension IJSUR(NSUR,*)
!local
integer :: i

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
! RG not adapted   RG?
!
 subroutine FLUID (BETB,PRPMAT,PRPLAS)

USE TIMESTEPMODULE
USE CONSTANTS, only: half, one, two, three, four, six, &
    root3, deg2rad, dfmin,dlogmin,dfmax,dlogmax
use iomodule, only: stderr
!
!    Routine to determine visco-plastic fluidity from
!    explicit time stepping limit (Cormeau) (LINEAR only)
!     (Ref.: Owen & Hinton, sec. 8.9)
!
implicit none

integer, parameter :: NSTR=6
!-pass
double precision   :: PRPMAT(6),PRPLAS(9),BETB(NSTR)

!-local
integer          :: nptype
double precision :: frangl, cohesn, stnhrd, fluidy, e, pois, debar
double precision :: sinphi
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

use modeltopology, only: ndof, nsd, nen
use constants,     only: zero
use algebra,       only: dot
!
! Sets up derivatives matrix BP
!
implicit none

!-pass
logical          :: QUAD
integer          :: IERR
double precision :: XL(NSD,NEN),DL(NDOF,NEN),BP(NDOF,NSD)
!-locl
double precision :: det, sh(4,NEN)
integer          :: i,j

!    Compute spatial derivatives of shape functions
call SHPTET (ZERO,ZERO,ZERO,XL,det,sh,IERR,1)

if (IERR.ne.0) then
    write(0,*) "SHPTET returned with error value:", ierr
    stop "Leaving GTECTON..."
endif

do i=1,3
    do j=1,3
        BP(i,j) = dot(SH(j,:), DL(i,:),nen)
    enddo
enddo

return
end
!-------------------------------------------------------------------------------
subroutine sdimnt (ISURF,SURF,B,ID,GRAV,PRPSED,DELTP)

USE MODELDEFINITION, only: nsurf, nsed
use iomodule,        only: stderr
!
implicit none

integer, parameter :: NSD=3,NDOF=3
!-pass
integer          :: isurf, id
double precision :: surf, b, grav, prpsed
double precision :: deltp

dimension ISURF(*),SURF(NSD,*),B(*),ID(NDOF,*),GRAV(*), &
          PRPSED(2)

!-locl
if (NSURF.le.0 .or. NSED.eq.0) then
    return
endif
!
 write(stderr,10)
10    format(1x,'"SDIMNT" NOT IMPLEMENTED IN 3-D VERSION')
!
 return
 end
!-----------------------------------------------------------------------

double precision function SRFHGT ()
 ! your guess is as good as mine... 'surface height'?
implicit none

SRFHGT = 0d0
return
end
!-----------------------------------------------------------------------

! TODO, does this thing do anything at all?
subroutine SlipperySkew (IDSLP,SKEW,NSELSD,IEN,X,D,TFAULT,DXE,LMF)

!
! Updates skew angles for slippery nodes
!
implicit none

integer, parameter          :: NSD=3, NEN=4, NDOF=3
double precision, parameter :: ZERO=0d0, TWO=2d0
double precision, parameter :: DEG=1.8d2/3.141592653589793d0
!-pass
integer          :: idslp, nselsd, ien, lmf
double precision :: skew, x, d, tfault, dxe
dimension IDSLP(*),SKEW(2,*),NSELSD(2,*),IEN(NEN,*),X(NSD,*), &
 D(NDOF,*),TFAULT(NDOF,*),DXE(NDOF,NEN,*),LMF(NDOF,NEN,*)
!
 return
 end
!-----------------------------------------------------------------------
subroutine LCOORD (X,XL,IEN,elemID)

use modeltopology,  only: nsd, nen
!
! localize coordinate data
!
implicit none
!-pass
integer          :: IEN(nen)
double precision :: X(NSD,*),XL(NSD,nen)
integer          :: elemID
!local
integer          :: k, iPoint

do iPoint = 1, 4
    k = IEN(iPoint)
    XL(:,iPoint) = X(:,k)
enddo

return
end
!-------------------------------------------------------------------------------
double precision function dilate (strain)

use constants, only: third
!
! computes volume change or volume change rate, depending on what S is.
!
implicit none
!-pass
double precision :: strain
dimension strain(*)
double precision, external :: AJ1
!-init
!
dilate = AJ1(strain)*THIRD
!
return
end
!-------------------------------------------------------------------------------
double precision function EFFEE (S)

use constants, only: dfmin, zero

!
! Compute effective deviatoric strain
!
implicit none
!-pass
double precision :: S
dimension S(*)
!-locl
double precision, external :: AJ2
double precision J2

!-init
effee = ZERO
!
J2 = AJ2(S)

if (J2.gt.dfmin) then
    effee = SQRT(J2)
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine STRDEV (S)

 use constants, only: third
!
! computes strain (rate) deviator
!
implicit none
!-pass
double precision :: S
dimension S(*)
!-local
double precision :: dil
double precision, external :: AJ1
!
dil = AJ1(S)*THIRD
S(1) = S(1)-dil
S(2) = S(2)-dil
S(3) = S(3)-dil
!
return
end
!-------------------------------------------------------------------------------
double precision function EFFSTN (S)

use constants, only: dfmin
!
! computes effective deviatoric stress. S is the full stress tensor.
!
implicit none
!-pass
double precision :: S(6)
!-locl
double precision J2
!-external
double precision, external :: AJ2

!
J2 = AJ2(S)
if (J2.gt.dfmin) then
    EFFSTN = SQRT(J2)
else
    write(*,*) "J2 is too small. Do not know what to do"
    stop "please contact model support"
endif
!
return
end
!-------------------------------------------------------------------------------
 subroutine STNTEN (STN,SL)

 USE MODELDEFINITION

implicit none

 integer, parameter :: NSTR=6
 double precision   :: STN(NSTR)
 double precision   :: SL(3,3)

 SL(1,1) = STN(1)
 SL(2,2) = STN(2)
 SL(3,3) = STN(3)

 SL(1,2) = STN(4)
 SL(2,1) = STN(4)

 SL(1,3) = STN(5)
 SL(3,1) = STN(5)

 SL(2,3) = STN(6)
 SL(3,2) = STN(6)

 return
 end subroutine
!-------------------------------------------------------------------------------
subroutine STNDEV (S)

use constants, only: third
!
implicit none
!-pass
double precision :: S(3)
!-locl
double precision :: p
!-external
double precision, external :: AJ1
!
 p = AJ1(S)*THIRD
 S(1) = S(1) - p
 S(2) = S(2) - p
 S(3) = S(3) - p
!
return
end
!-------------------------------------------------------------------------------
double precision function PRESSURE (S,POIS)
!
use constants, only: third

implicit none
!-pass
double precision S(3), POIS
!dimension S(*)
!-locl
!-external
double precision, external :: AJ1

!
pressure = AJ1(S)*THIRD
!
return
end

!-------------------------------------------------------------------------------


integer function sideIDofTriangle(iTriangle, iElem)

use algebra, only: arraycontains
use meshdatamodule, only: meshdatactx

implicit none

integer :: iTriangle, iElem
integer :: pos1, pos2, pos3, pos4

            ! Compute on which site of the element a triangle is.
            ! To refresh our memory:

!  | \            face 1: nodes 1 2 4
!  |__\2_____X         2:       1 3 2
!  1\ /                3:       1 4 3
!    \4                4:       2 3 4


            ! check the positions of the the IEN values in the triangle.
            ! The only one that is not in the triangle determines the value.
            call arraycontains(meshdatactx%triangleConnectivity(:,iTriangle), &
                               3, &
                               meshdatactx%ien(1,iElem), &
                               pos1)
            call arraycontains(meshdatactx%triangleConnectivity(:,iTriangle), &
                               3, &
                               meshdatactx%ien(2,iElem), &
                               pos2)
            call arraycontains(meshdatactx%triangleConnectivity(:,iTriangle), &
                               3, &
                               meshdatactx%ien(3,iElem), &
                               pos3)
            call arraycontains(meshdatactx%triangleConnectivity(:,iTriangle), &
                               3, &
                               meshdatactx%ien(4,iElem), &
                               pos4)

            ! array length is 3, so if a point is not in it, it will will return a 4

            if (pos1.eq.4) then
                sideIDofTriangle = 4
            else if (pos2.eq.4) then
                sideIDofTriangle = 3
            else if (pos3.eq.4) then
                sideIDofTriangle = 1
            else if (pos4.eq.4) then
                sideIDofTriangle = 2
            else
                sideIDofTriangle = 0
            endif

end function

!-------------------------------------------------------------------------------
#ifdef SHEARZONE_TRACTIONS_MARIO_VERSION
subroutine SHEARTRAX (B,BTOT,X,D,IDX,IEN,TFAULT,NSLIP,DX,DT,shearETA,shearTHICK,ELNOSHEARZONE,shearnpl,shearQpl,shearT)
    ! Updates the load vector B for shear tractions on a shear zone bounded by elements with numbers ELNOSHEARZONE with viscosity
    ! shearETA and thickness shearTHICK, and optionally with power-law exponent shearnpl and activation energy shearQpl and 
    ! temperature shearT.
    ! Note: if the power law variables are not provided, the formulaiton for Newtonian viscosity is used.
    !       If they are provided, the viscosity is interpreted as the reference viscosity eta_0, i.e. half of the inverse of
    !       the pre-exponential coefficient A_D.
    ! The shear zone is parameterized by a slippery node interface. 
    ! *** ASSUMPTION: the interface is perpendicular to the local z-direction (see Euler angles) ***
    ! The current implementation is based on a linear (Newtonian) dependence of shear traction on differential velocities DX and DY.
    ! Note: NUMELSHEARZONE, the number of shear traction entries, is read from MODELDEFINITION to avoid conflicts between 
    ! variables with the same name or differently-named variables defined in different places but describing the same thing.
    ! Citation: D'Acquisto et al. (in prep)
    
    USE MODELDEFINITION
    use constants,      only: useDXE, ONE, TWO, THREE, SIX
    use modelctx,       only: getrank
    use iomodule,       only: stderr
    use debugmodule,    only: iecho

    implicit none
    
    integer, parameter          ::  NEN=4, NSD=3, NDOF=3
    integer                     ::  IDX(NDOF,*), &                  ! partition-local equation numbers for each partition-local node number
                                    IEN(NEN,*), &                   ! partition-local node numbers for each partition-local element number.
                                    NSLIP(5,*), &                   ! partition-local element and node numbers and associated slippery dofs.
                                                                    ! Length NUMSLP is the partition-local number of slippery node entries, i.e.,
                                    ELNOSHEARZONE(NUMELSHEARZONE)   ! list of partition-local element numbers along the shear zone(s)
    double precision            ::  B(*), &                         ! incremental load for each partition-local equation 
                                    BTOT(*),&                       ! cumulative load for each partition-local equation number
                                    DX(NDOF,*), &                   ! differential velocity for each partition-local slippery node number
                                    DT, &                           ! Timestep size
                                    shearETA(NUMELSHEARZONE), &     ! shear zone viscosity for given local element numbers
                                    shearTHICK(NUMELSHEARZONE), &   ! shear zone thickness for given local element numbers
                                    X(NSD,*), &                     ! coordinates of partition-local node numbers
                                    D(NDOF,*), &                    ! displacements of partition-local node numbers
                                    TFAULT(NDOF,*)                  ! faulted node displacements of partition-local node numbers. Length NUMFN
    INTEGER                     ::  ielm, ielm2, islpn, counter, counter2, nodesnp(3), nodesnp_all_asc(3,NUMELSHEARZONE), lien(3), &
                                    idim, NPloc, ElmPlocal, LocalEq, ieqno, numeq, eqnoarr(NUMELSHEARZONE*3*2)
    logical                     ::  neweqno
    double precision            ::  xl(NSD,NEN), dl(NDOF,NEN), xn(NSD), area, dF
    INTEGER                     ::  shearnpl                        ! power-law coefficient
    double precision            ::  shearQpl(NUMELSHEARZONE)        ! power-law activation energy (J/mol)
    double precision            ::  shearT(NUMELSHEARZONE)          ! Shear-zone temperature for each element
    double precision, parameter ::  Rpl=8.31446261815324            ! Universal gas constant (for power law)

    if (NUMELSHEARZONE.le.0) then
        return
    endif
    if (NUMSLP.le.0) then
        write(stderr,10) getrank(),NUMSLP,NUMELSHEARZONE
10      format(1x,"FATAL error in sheartrax: rank ",I5,": NUMSLP = ",I8," while NUMELSHEARZONE = ",I8)
#ifdef UnitTestTractions
        call exit(1)
#else
        call exitp(1)
#endif
    endif

    numeq = 0 ! Reset counter of number of unique equation numbers

    ! Loop over shear zone elements
    do ielm = 1,NUMELSHEARZONE
        ElmPlocal = ELNOSHEARZONE(ielm)
        islpn = 1
        counter = 1
        nodesnp = 0
        do while (counter.le.3)                       ! Set interface node numbers in the current element
            if (NSLIP(1,islpn).eq.ElmPlocal) then
                nodesnp(counter) = NSLIP(2,islpn)       ! set partition-local node number of the overall mechanical problem
                counter = counter+1
            endif
            islpn = islpn+1
        enddo

        ! Check that the same group of nodes is not considered twice (by considering the two faces of the same surface)
        nodesnp_all_asc(1,ielm) = nodesnp(1)
        do counter = 2,3
            do counter2 = counter-1,1,-1
                if (nodesnp(counter).lt.nodesnp_all_asc(counter2,ielm)) then
                    nodesnp_all_asc(counter2+1,ielm) = nodesnp_all_asc(counter2,ielm)
                    nodesnp_all_asc(counter2,ielm) = nodesnp(counter)
                else
                    nodesnp_all_asc(counter2+1,ielm) = nodesnp(counter)
                    exit
                end if
            end do
        end do
        do ielm2 = 1,ielm-1
            if (nodesnp_all_asc(1,ielm2).eq.nodesnp_all_asc(1,ielm) .and. nodesnp_all_asc(2,ielm2).eq.nodesnp_all_asc(2,ielm) .and.&
                nodesnp_all_asc(3,ielm2).eq.nodesnp_all_asc(3,ielm)) then
                write(stderr,*) "FATAL error in sheartrax: two elements with the same slippery nodes were provided"
#ifdef UnitTestTractions
                call exit(1)
#else
                call exitp(1)
#endif
            end if
        end do

        ! Get element-local node numbers defining the interface
        lien = 0
        do counter = 1,3
            do islpn = 1,NEN
                if (nodesnp(counter).eq.IEN(islpn,ElmPlocal)) then
                    lien(counter) = islpn
                    exit
                end if
            end do
        enddo

        call LCOORD (X,xl,IEN(1,ElmPlocal),ElmPlocal)   ! Copy coordinates of the current element into xl
        if (LGDEF.ne.0) then
            call LDISP  (dl,D,IEN(1,ElmPlocal),NDOF,NEN)
            call AddFaultDisplacement (dl,ElmPlocal,TFAULT,NDOF,NEN,1001)
            call ADDSNE (dl,ElmPlocal,NDOF,NEN,useDXE)
            call REZONE(xl,dl,ONE)
        endif
                
        ! Calculate surface area (needed) and unit normal vector (not needed for now) of the element side
        call FCGEOM (xl,lien,area,xn)
        if (debug) then
            write(*,*) "Rank", getrank(), ": the area of the side of element", ElmPlocal, "delimited by slippery nodes is", area
        end if

        do counter = 1,3                        ! Loop over slippery nodes in the current element
            NPloc = nodesnp(counter)
            if (debug) then
                write(*,*) "Rank", getrank(), ": the x- and y-velocities of node", NPloc, "are", &
                           DX(1,NPloc)/DT, "and", DX(2,NPloc)/DT, ", respectively"
            end if
            do idim = 1,2                       ! Loop over local shear zone dimensions (x and y)
                if (IDX(idim,NPloc).gt.0) then
                    if (shearnpl.gt.0) then
                        if (debug) then
                            write(*,*) "Rank", getrank(), ": using power-law creep"
                        end if
                        dF = - (THREE*(DX(1,NPloc)/DT**TWO+DX(2,NPloc)/DT**TWO)**((ONE-shearnpl)/TWO/shearnpl)*DX(idim,NPloc)/DT &
                             + (DX(1,nodesnp(1))/DT**TWO+DX(2,nodesnp(1))/DT**TWO)**((ONE-shearnpl)/TWO/shearnpl) &
                             *DX(idim,nodesnp(1))/DT + (DX(1,nodesnp(2))/DT**TWO+DX(2,nodesnp(2))/DT**TWO)**((ONE-shearnpl) &
                             /TWO/shearnpl)*DX(idim,nodesnp(2))/DT + (DX(1,nodesnp(3))/DT**TWO+DX(2,nodesnp(3))/DT**TWO) &
                             **((ONE-shearnpl)/TWO/shearnpl)*DX(idim,nodesnp(3))/DT)/SIX*(TWO*shearETA(ielm)/shearTHICK(ielm)) &
                             **(ONE/shearnpl)*exp(shearQpl(ielm)/Rpl/shearT(ielm) / shearnpl) * area / THREE
                    else
                        dF = -(DX(idim,nodesnp(1))/DT+DX(idim,nodesnp(2))/DT+DX(idim,nodesnp(3))/DT+THREE*DX(idim,NPloc))/DT &
                             / THREE * shearETA(ielm) / shearTHICK(ielm) * area / THREE
                    end if
                    if (debug) then
                        write(*,*) "Rank", getrank(), ": the velocity-dependent shear force contribution of element", ElmPlocal, &
                                   "on node", NPloc, "is", dF
                    end if
                    LocalEq = IDX(idim,NPloc)
                    if (LocalEq.gt.0) then
                        B(LocalEq) = B(LocalEq) + dF
                        ! Save list of (unique) equation numbers used for forces
                        neweqno = .true.
                        do ieqno = 1,numeq
                            if (eqnoarr(ieqno).eq.LocalEq) then
                                neweqno = .false.
                                exit
                            end if
                        end do
                        if (neweqno) then
                            numeq = numeq+1
                            eqnoarr(numeq) = LocalEq
                        end if
                    end if
                endif
            enddo
        enddo
    enddo

    do ieqno = 1,numeq ! Subtract BTOT entries from shear forces B entries
        LocalEq = eqnoarr(ieqno)
#ifdef UnitTestTractions
        write(*,*) "B(", LocalEq, ") =", B(LocalEq), ", BTOT(", LocalEq, ") =", BTOT(LocalEq), ", new B(", LocalEq, ") =", &
                   B(LocalEq)-BTOT(LocalEq)
#endif
        B(LocalEq) = B(LocalEq)-BTOT(LocalEq)
    end do
end subroutine
#endif
!-------------------------------------------------------------------------------
#ifdef SHEARZONE_TRACTIONS
subroutine SHEARTRAX (IFLAG,NSTEP)
implicit none
integer :: IFLAG,NSTEP
! todo why you no do anything?
end subroutine
#endif
