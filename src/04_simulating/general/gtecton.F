program main
!
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! *                                                                           *
! *                             G - T E C T O N                               *
! *                                                                           *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! This is a finite element code, originally written by Jay Melosh and
! Arthur Raefsky, especially to solve tectonic problems. It solves the
! mechanical equilibrium and energy equations independently.
! G-TECTON computes solutions to static elastic problems as well as
! time-dependent stresses and displacements due to faulting in a
! -elastic, non-newtonian earth. The thermal finite element code
! was written by Rob Govers.
!
! The code can run in parallel on different cores ("CPU's"). The input data
! needs to be distributed over these cores by partitioning. Each compute thread
! is called a partition.
!
! Please type:
!       > man pln
! or
!       > man 7 pln
! for further information about the input data and the variables below:

!-------------------------------------------------------------------------------
! NAMES OF VARIABLES AND ARRAYS THAT ARE USED THROUGHOUT THE CODE:
!-------------------------------------------------------------------------------
!
! I. GLOBAL DATA

!    NUMAT                      Number of material types
!    NUMEL=nelocal              Number of elements in the partition, in meshdatactx
!          neglobal             Total number of elements, in meshdatactx
!    NUMFN                      Number of split node entries
!    NUMNP=nvlocal              Number of nodes ("vertices") in the partition, in meshdatactx
!          nvglobal             Total number of nodes, in meshdatactx
!    nvertices                  node counter used to create the local numbering
!    NEQ=NEQlocal               Number of equations in this partition for the mechanical problem
!        NEQglobal              Global number of equations for the mechanical problem
!    nequations                 Global number of equations for the mechanical problem
!    nequationsT                Global number of equations for thermal problem
 

!    IDIAG(NEQ)                 Addresses of As diagonal elements (classical)
!    ITDIAG(NTEQ)               Addresses of AAs diagonal elements (classical)
!
!    BTOT(nequations)           Accumulated mechanical force vector.
!    BRES(nequations)           Residual mechanical force vector.
!    BT(nequationsT)            Global thermal force vector.

!    ID(NDOF,nvlocal)           Global equation numbers for normal degrees of freedom ("dofs") in the mechanical problem for a partition
!    IDglobal(NDOF,nvglobal)    Global equation numbers in the mechanical problem.
!    IDX(NDOF,nvglobal)         Global equation numbers for slippery node degrees of freedom. Only nodes relevant for own partition
!    IDXglobal(NDOF,nvglobal)   Global equation numbers for slippery node degrees of freedom. Identical for every partition
!    IDT(NUMNP)                 Global equation # of each degree of freedom in the thermal problem.

!    X(NSD,NUMNP)               Nodal coordinates.
!    D(NDOF,NUMNP)              Nodal displacements.
!    DELD(NDOF,NUMNP)           Nodal displacement increments 
!    T(NUMNP)                   Nodal temperatures.
!    TLST(NUMNP)                Nodal temperatures in last time step.

!    GRAV(NDOF)                 Gravity vector (uniform over grid).
!    SKEW(2,NUMNP)              Nodal rotation ("Euler") angles for skew boundary conditions.
!    OSKEW(2,NUMNP)             Nodal Euler angles of previous time step

!    ISURF(NSURF)               Surface nodal point numbers.
!    SURF(NSD,NSURF)            Surface nodal coordinates.
!    SURF0(NSD,NSURF)           Initial surface nodal coordinates.
!
!
! II.    ELEMENT DATA
!
!    A.    ELEMENT SYSTEM
!
!    IEN(NEN,NUMEL)             Node numbers in each element.
!    LM(NDOF,NEN,NUMEL)         Mechanical equation number of element nodes
!    LMT(NEN,NUMEL)             Thermal equation numbers element nodes
!    MAT(NUMAT)                 Element material number 
!    IJSUR(NSUR,NUMEL)          Element neighbours array
!
! B.    ELEMENT PROPERTIES
!
! (i) GENERAL MATERIAL PROPERTIES
!
!    PRPMAT(1,NUMAT)        Youngs modulus.
!    PRPMAT(2,NUMAT)        Poissons ratio.
!    PRPMAT(3,NUMAT)        Viscosity.
!    PRPMAT(4,NUMAT)        Power.
!    PRPMAT(5,NUMAT)        Density.
!    PRPMAT(6,NUMAT)        Thickness.
!
! (ii) POWERLAW CREEP MATERIAL PROPERTIES
!
!    PRPCRP(1,NUMAT)        Power law activation energy.
!    PRPCRP(2,NUMAT)        Power law scaling constant.
!    PRPCRP(3,NUMAT)        Power law power to stress.
!    PRPCRP(4,NUMAT)        Viscosity minimum
!    PRPCRP(5,NUMAT)        Viscosity maximum
!    PRPCRP(6,NUMAT)        Power law activation volume.
!
! (iii) STRAIN (RATE) WEAKENING PROPERTIES
!
!    PRPDIF(1,NUMAT)         Viscosity weakening factor
!    PRPDIF(2,NUMAT)         Viscosity weakening interval
!    PRPDIF(3,NUMAT)         Yield stress weakening factor
!    PRPDIF(4,NUMAT)         Yield stress weakening interval
!    PRPDIF(5,NUMAT)         not used
!    PRPDIF(6,NUMAT)         not used
!    PRPDIF(7,NUMAT)         not used
!    PRPDIF(8,NUMAT)         backup Yield stress
!    PRPDIF(9,NUMAT)         backup viscosity
!
!#ifdef EXPERIMENTALTRANSIENTDISCREEP
! (iv) PLASTICITY MATERIAL PROPERTIES
!
!    PRPLAS(1,NUMAT)        Plasticity type.
!    PRPLAS(2,NUMAT)        Friction angle in degrees.
!    PRPLAS(3,NUMAT)        Uniaxial yield stress.
!    PRPLAS(4,NUMAT)        Fluidity.
!    PRPLAS(5,NUMAT)        Yield power.
!    PRPLAS(6,NUMAT)        Yielding law.
!    PRPLAS(7,NUMAT)        Strain hardening.
!    PRPLAS(8,NUMEL)        Uniaxial yield stress backup
!    PRPLAS(9,NUMEL)        Fluidity backup
!
! (v) TRANSIENT DISLOCATION CREEP MATERIAL PROPERTIES
!
!    PRPDISLC(1,NUMAT)        Activation energy.
!    PRPDISLC(2,NUMAT)        Pre-exponent constant.
!    PRPDISLC(3,NUMAT)        Peierls stress.
!    PRPDISLC(4,NUMAT)        Burgers length.
!    PRPDISLC(5,NUMAT)        Beta Threshold Stress Coefficient.
!    PRPDISLC(6,NUMAT)        M coefficient Taylor stress.
!    PRPDISLC(7,NUMAT)        Rate coefficient dislocation change due to grain boundary diffusion.
!    PRPDISLC(8,NUMAT)        Rate coefficient dislocation change due to pipe diffusion
!    PRPDISLC(9,NUMAT)        Activation volume
!    PRPDISLC(10,NUMAT)       Initial ratio of Taylor stress to elastic stress
!
! (vi) DIFFUSION CREEP MATERIAL PROPERTIES
!
!    PRPDIFFC(1,NUMAT)        Diffusion creep activation energy.
!    PRPDIFFC(2,NUMAT)        Diffusion creep scaling constant.
!    PRPDIFFC(3,NUMAT)        Grain size exponent.
!    PRPDIFFC(4,NUMAT)        Diffusion creep activation volume.
!    PRPDIFFC(5,NUMAT)        Grain size.
!endif
!
!
! (vii) THERMAL MATERIAL PROPERTIES
!
!    PRPTEM(1,NUMAT)        Heat production per unit volume.
!    PRPTEM(2,NUMAT)        Heat capacity at constant pressure.
!    PRPTEM(3,NUMAT)        Thermal expansivity.
!    PRPTEM(4,NUMAT)        Thermal X-conductivity.
!    PRPTEM(5,NUMAT)        Thermal Y-conductivity.
!    PRPTEM(6,NUMAT)        Thermal Z-conductivity.
!
! C.    POST-ELEMENT PROCESSING DATA
!
!    BETA(NSTR,NUMEL)       Viscous strain rate matrix storage.
!    BETB(NSTR,NUMEL)       Plastic strain rate matrix storage.
!    VPSTR(NSTR,NUMEL)      Plastic strain.
!    DMAT(NSTR,NSTR,NUMEL)  Material matrix for each element.
!    STN(NSTR,NUMEL)        Stresses
!    STR(NSTR,NUMEL)        Strains
!    RATE(NSTR,NUMEL)       Total strain rate saved from last step.
!    FLUX(NSD,NUMEL)        Element heat flux.
!
!#ifdef EXPERIMENTALTRANSIENTDISCREEP
!    STNTAY(NSTR,NUMEL)     Taylor stress
!    STNTAYRATE(NSTR,NUMEL) Taylor stress time derivative
!    STNTHRES(NSTR,NUMEL)   Threshold stress
!    RATEDISLEFF(NUMEL)     Effective scalar dislocation creep strain rate
!    INVEFFVISCDISL(NUMEL)  Effective inverse viscosity dislocation creep
!    INVEFFVISCDIFF(NUMEL)  Effective inverse viscosity diffusion creep
!#endif
!
! D.    SPLIT NODE DATA
!
!    NFAULT(3,NUMFN)        Element # and node # of split node and
!                           mode of application of split node
!    FAULT(NDOF,NUMFN)      Split vector for 1/2 node.
!    DFAULT(NDOF,NUMFN)     Fault displacement change
!    TFAULT(NDOF,NUMFN)     Total fault displacement at any given time
!    FAULTL(NDOF,0)         Total fault displ. at last refactor and reform
!    LMF(NDOF,NEN,NUMEL)    LM array for faulted degrees of freedom
!
! E.    SLIPPERY NODE FREE SLIP INTERFACE DATA
!
!    NSLIP(5,NUMSLP)        Element #, node #, and weight (+/-) on dof,
!                           depending on side of fault.
!    IDSLP(NUMSN)           node # of each slippery node.
!    IDSLE(NUMSE)           elements with slippery nodes.
!    DIFORC(NDOF,NUMNP)     Differential force applied across interface.
!    DX(NDOF,NUMNP)         Differential displacements across interfaces.
!    DXE(NDOF,NEN,NUMEL)    Element based total differential displacements
!                           in a global coordinate frame
!    DELX(NDOF,NEN,NUMEL)   Element based incremental differential
!                           displacements in a global coordinate frame
!    WINX(NDOF,NUMNP)       Winkler restoring forces for slippery DOF
!    IWINX(NDOF,NUMNP)      Application mode for Winkler forces
!    LMX(NDOF,NEN,NUMEL)    LM array for slippery degrees of freedom.
!    IDX(NDOF,NUMNP)        ID array for extra degrees of freedom in the
!                           mechanical problem.
!
!  III. NODAL BOUNDARY CONDITIONS
!
!    IBOND(NDOF,NUMNP)      Boundary condition code of mechanical problem.
!    BOND(NDOF,NUMNP)       Mechanical boundary conditions.
!    LINK(3,NLINK)          IDOF, NODE1 links to NODE2, i.e. displacement
!                           of NODE1 is equal to that of NODE2.
!    IWINK(NDOF,NUMNP)      Winkler application mode
!    WINK(NDOF,NUMNP)       Winkler restoring force coefficients
!    FMAXW(NDOF,NUMNP)      Maxwell bc magnitude
!    FTOT(NDOF,NUMNP)       Maxwell force vector
!    ITBC(NUMNP)            Boundary condition code of thermal problem.
!    TBC(NUMNP)             Thermal boundary conditions.
!    ITANOM(3,NTANOM)       Thermal anomaly node, mode and time AFTER which
!                           it is applied.
!    TANOM(NTANOM)          Nodal thermal anomaly
!
!
!   IV.    FLUX-TYPE BOUNDARY CONDITIONS
!
!    A.    PRESSURE BOUNDARY CONDITIONS
!
!    IELNO(NUMPR)           Element # where pressure is applied.
!    ISIDE(NUMPR)           Side of element where pressure is applied.
!    PRES(NUMPR)            Applied pressure.
!
!    B.    STRESS BOUNDARY CONDITIONS
!
!    IELSTR(NUMSTR)         Element # where stress is applied.
!    ISSIDE(NUMSTR)         Side of element where stress is applied.
!    ISTR(2,NUMSTR)         Application times (from/to)
!    STRS(NUMSTR,NSTR)      Applied stress.
!
!    C.    HEAT FLUX BOUNDARY CONDITIONS
!
!    IFLX(NFLX)             Element # where heat flux is applied.
!    IFLS(NFLX)             Side of element where heat flux is applied.
!    BFLX(NFLX)             Normal boundary heat flux.
!
!    D. WINKLER RESTORING PRESSURE BOUNDARY CONDITIONS
!
!    IWELM(NUMWNK)          Element # where Winkler pressure is applied
!    IWSIDE(NUMWNK)         Side of element where pressure is applied.
!    IWTIME(NUMWNK)         Application mode for Winkler pressures
!    WPRES(NUMWNK)          Winkler pressure spring constants.
!
! 
!   V  BODY-FORCE BOUNDARY CONDITIONS
!
!        (PRE-) STRESS BOUNDARY CONDITIONS
!
!    ISELM(NPRE)            Element # where stress is applied
!    ISTIME(NPRE)           Application mode for stresses
!    STN0(NSTR,NPRE)        Stresses.
!
!
!   VI.    OUTPUT STORAGE AND CONTROL
!
!    IMPRINT(NMPRT)         Time step number where full outputs of the
!                           mechanical solution are produced.
!    ITPRINT(NTPRT)         Time step number where full outputs of the
!                           thermal solution are produced.
!    IMATPR(NPRMAT)         Time step where matrix diagonal outputs
!                           are produced
!
!
!  VII.    TIME STEP DATA
!
!    MAXSTP(NINTG)          Number of steps in each step size group
!    DELT(NINTG)            Time step size in each group.
!    ALFA(NINTG)            Implicitness parameter for mechanical problem.
!    TALF(NINTG)            Implicitness parameter for thermal problem.
!    ICALC(NTCALC)          Time step number when thermal calculation is
!                           required.
!
!-------------------------------------------------------------------------------
!    LOCAL ARRAYS. PRESET DIMENSION IS GIVEN IN SECOND COLUMN.
!
!    S(NEE,NEE)        = S(576)      Local mechanical stiffness matrix
!    ST(NEN,NEN)       = ST(64)      Local thermal stiffness matrix
!    STEMP(NEE,NEE)    = STEMP(576)  Scratch storage for stiffness
!    P(NEE)            = P(24)       Local mechanical load vector.
!    PT(NEN)           = PT(24)      Local thermal load vector.
!    IDMAT(NSTR)       = IDMAT(6)    Addresses of square matrix diagonals
!    DL(NDOF,NEN)      = DL(24)      Local element displacements
!    TL(NEN)           = TL(8)       Local element temperatures
!    XL(NSD,NEN)       = XL(24)      Local element nodal coordinateS
!    ST(NSTR)          = ST(6)       Local element stress
!    DB(NSTR)          = DB(6)       Updated local element stress
!    EVP(NSTR)         = EVP(6)      Local element stress changes
!    EE(NSTR)          = EE(6)       Local element strains
!    DELSTR(NSTR)      = DELSTR(6)   Local element strain changes
!    DMATL(NSTR,NSTR)  = DMATL(36)   Local material matrix
!    AI(NSTR)          = AI(6)       Scratch vector for matrix inverse
!    ALIN(NALIN)       = ALIN(21)    Scratch array for matrix inverse
!    BETAJ(NSTR,NSTR)  = BETAJ(6,6)  Local jacobean matrix
!    ROT(NDOF,NDOF)    = ROT(3,3)    Rotation matrix for skew boundaries
!    DTEMP(NDOF)       = DTEMP(3)    Scratch storage for displacements
!
!
!-------------------------------------------------------------------------------
!
!        SCALERS USED IN FINITE ELEMENT CODE G-TECTON
!
!-------------------------------------------------------------------------------
!
!    ALFAP              = alpha parameter (implicitness of time stepping) mechanical problem
!    ASCIN              = logical saying if input is ascii or not
!    ASCOUT             = logical saying if output is ascii or not
!    UNIYLD             = Uniaxial yield stress for viscoplasticity
!    DELTP              = Scalar time step size mechanical problem
!    FEDSK              = Compact FE output file
!    FEIN               = Input file containing geometry and mechanical input data
!    FEOUT              = Verbose output file
!    FLUIDY             = Plastic fluidity
!    FRANGL             = Plastic friction angle
!    IADVEC             = advecting grid parameter for thermal calculation
!    ICVIS              = material properties switch
!    IMPRT              = Number of last full mechanical output
!    IMWORK             = (0,1) do not or do make a mechanical implicit calculation.
!    IMTWRK             = (0,1) do not or do calculate thermo-mechanical coupling loads.
!    IDIFT              = (0,1) do not or do compute differential temperatures
!    IRESDU             = Flag of addition of residual stress or force
!    ISHELL             = (0,1) flat (0) or spherical (1) shell
!    ITECT              = (0,1) do or do not add tectonic contribution from thermal anomalies.
!    ITMODE             = 0, calculate initial steady-state temperatures
!    ITPRES             = (0,1) do or do not add tectonic contributions resulting from temperature changes.
!    ITPRT              = Number of last full thermal output
!    ITWORK             = (0,1) do not or do make a thermal LU decomp.
!    IVLIM              = (0-2) do not or do put a lower limit on viscosities, determined by time step size
!    IWORD              = (1,2) for single or double precision
!    NPRE               = Number of pre-stress entries
!    LGDEF              = switch for use of large deformation stiffness
!    MAATOT             = Total length of AA-array 4-byte words
!    MAXIT              = Max number of iterations between mech. rff
!    MAXTIT             = Max number of iterations between therm. rff
!    MODE               = Solution mode
!    MTOT               = Total length of A-array 4-byte words
!    NAA                = Number of non-zero entries in the stiffness matrix of the thermal problem.
!    NALIN              = Number of entries in a triangular matrix of the mechanical problem
!    NDOF               = Number of degrees of freedom at each node
!    NEC = NSD*NEN      = Number of element coordinates
!    NEE = NDOF*NEN     = Number of mechanical element equations
!    NEN                = Number of nodes in an element
!    NEQ                = Number of equations in the mechanical stiffness matrix
!    NFLOW              = Plastic flow law: 0. Exponential, 1. powerlaw.
!    NINTG              = Number of time step groups mechanical problem
!    NMAXW              = Number of nodes with Maxwell bc
!    NMPRT              = Number of full stress, displacement outputs
!    NLINK              = Number of linked nodes
!    NPRMAT             = Number of matrix diagonal outputs
!    NPTYPE             = Plasticity type
!    NSD                = Number of space dimensions
!    NSED               = (0,1) switch for including sedimentary transport loads or not
!    NSLSKEW            = Number of elements with skew angles parallel to slippery fault
!    NSURF              = Number of surface nodal points
!    NSIZEA             = Number of non-zero entries in the stiffness matrix of the mechanical problem.
!    NSTEP              = Time step number mechanical problem
!    NSTR               = Number of stress tensor components
!    NTANOM             = Number of thermal anomaly inputs
!    NTCALC             = Number of thermal solutions.
!    NTEQ               = Number of equations in the thermal stiffness matrix
!    NTPRT              = Number of thermal outputs
!    NUMPR              = Number of pressure load entries
!    NUMROT             = Number of nodes with skew boundaries
!    NUMSLP             = Number of slippery node entries
!    NUMSE              = Number of elements with slippery nodes
!    NUMSN              = Number of slippery nodes
!    NWINK              = Number of Winkler forces
!    NWINKX             = Number of Winkler forces on slippery nodes
!    OUTFIL(2)          = Logicals for generating output files or not
!    RADIUS             = Spherical shell radius
!    REFACT             = Logical; Refactor mechanical stiffness
!    STNHRD             = Plastic strain hardening
!    TALFP              = Scalar alpha parameter thermal problem
!    TEMPIN             = Thermal input file
!    YLDPWR             = Plastic yield power
!    VERSION            = FE-program version
!
!-------------------------------------------------------------------------------
!   
!       FILE POINTERS USED THROUGHOUT GTECTON
!   
!-------------------------------------------------------------------------------
!
!       NAME            FILE DESCRIPTION
!
!       luin            can either point to:
!                         modeldata.dat  (in AOmodule.F)
!                         partitioninfo file (in mergers.F)
!       luout(1)        feout file
!       luout(2)        fedsk file
!                       These FEDSK files need a more thorough explanation
!                       When written in ascii mode ('as' as command line argument)
!                       the data is clearly labeled, indexed and ordered, this is
!                       not obvious in the binary files. The data is entered by 
!                       'record'. It is also to be read by record, otherwise the
!                       read command will fail. It is not possible to read it as 
!                       a byte stream.
!  The FEDSK file contains the following records:
!                        
!                        
!
!-------------------------------------------------------------------------------
!
!          SUBROUTINES THAT DEPEND ON TOPOLOGY ARE LISTED HERE
!            THEY ARE CONTAINED IN THE APPROPRIATE LIBRARIES
!
!
!       NAME            FUNCTION
!
!       ADDFLX          NORMAL BOUNDARY FLUX CONTRIBUTIONS TO LOAD VECTOR
!       ADDPR           TRANSLATES APPLIED PRESSURES INTO LOAD VECTOR
!       ADDSTR          STRESS LOAD COMPUTATION
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
!       FLOWVP          VISCOPLASTIC STRAIN RATE
!       FORMBT          VISCOUS FLOW CORRECTION FOR STRAIN
!       FORMJB          JACOBEAN MATRIX FOR VISCOUS FLOW
!       GRAVLD          GRAVITATIONAL BODY FORCE COMPUTATION
!       INTCON          DIMENSIONAL AND TOPOLOGIC DATA
!       INVAR           COMPUTES STRESS INVARIANTS
!       LHEAT           HEAT PRODUCTION CONTRIBUTION TO LOAD VECTOR
!       MATERL          MATERIAL MATRICES
!       PRESUR          ADD HYDROSTATIC PRESSURE RHO*G*H TO STRESS
!       REZONE          LOCAL REZONE FOR LARGE DEFORMATION
!       SHAP20          SHAPE FUNCTION FOR 2-D ELEMENTS
!       STIFF           STIFFNESS MATRIX
!       TANGVP          VISCOPLASTIC JACOBIAN MATRIX
!       TCENTR          COMPUTES ELEMENT CENTER TEMPERATURES
!       THGRAD          THERMAL GRADIENTS
!       TSTIFF          LOCAL CONDUCTION STIFFNESS MATRIX
!       YIELDF          CONSTRUCTS PLASTIC FLOW VECTOR
!
!-------------------------------------------------------------------------------
!
!     CURRENTLY AVAILABLE AND PLANNED LIBRARIES ARE:
!
!    NTYPE #     FILENAME             DESCRIPTION             STATUS
!
!       1       PLNLIB .......PLANE STRESS-STRAIN             A
!       2       OPNLIB .......OUT-OF-PLANE DISPLACEMENTS      A
!       3       AXILIB .......AXISYMMETRIC                   N/A
!       4       SPHLIB .......SPHERICAL SYMMETRY             N/A
!       5       F3DLIB .......3D BRICKS AND WEDGES            A
!       6       F3DLIB .......3D TETRAHEDRA                   A
!
!*************************************************************
 
 
USE MODELDEFINITION
USE MODELDATAMODULE, only: modeldatactx, &
                           idiag, itdiag, &
                           showfbc
USE MODELTOPOLOGY
USE MESHDATAMODULE,  only: meshdatactx, &
                           elementNeighborTableComplete
USE MATERIALSMODULE, only: LMAT
USE TIMESTEPMODULE
USE FILEMODULE
use debugmodule
use formatsmodule
use constants,       only: readRestart
#ifdef SPARSE
use spetscmodule,    only: PETSc_openlogfile, & 
                           PETSc_closelogfile
USE AOMODULE
#endif
USE MODELCTX
use iomodule,        only: ascout, luin, optijd
use plotControl,     only: jShell

implicit none

!-local
integer :: ierr
integer :: irank
integer :: n

#ifdef SPARSE
!#include "petsc/finclude/petscsysdef.h"
!#include "petsc/finclude/petscmatdef.h"
!#include "petsc/finclude/petscvecdef.h"
!#include "petsc/finclude/petsckspdef.h"



#ifdef TIMING
   PetscMPIInt irank
   PetscLogDouble time_0,time_start, time_finish, time_diff
   PetscErrorCode pierr
#endif 
#endif
!SPARSE

integer :: datetime(8)

!-init
feplt  = .false.

! this setting is only read when gtecton is run, not 
! that plnplt needs to have this variable set separately.
! It is set in 
#ifdef SPARSE
IWORD  = 2 ! cannot be anything else
#else
! SNGL    IWORD  = 1
IWORD  = 2
#endif

ierr   = 0
NSTEP  = 0
TIME   = 0d0
ISTART = 0
MAXITIME = 0
JSHELL = 0
SHOWFBC  = .false.
! SHOWFBC  = .true.
LMAT = .false.



   optijd = .false.


elementNeighborTableComplete = .false.



!    initialize mpi context
!write(*,*) "Initializing MPI"
call mpi_init(ierr)
!write(*,*) "Initialized MPI"

#ifdef SPARSE


if (ierr.ne.0) then
    write(*,*) 'failed MPI initialize, with error ', ierr
    call xit(1," ")
endif
#endif


#ifdef TEST
! this is a direct link to the debugmodule subroutine testIt
! This is not a full automatic tesing suite. That is left to
! the external automatic benchmarking. It is a way to test 
! subroutines and the like, without needing to write external
! test program. Moreover, all the functionality in gtecton is 
! present here, and we have had an MPI init, so we can do more
! or less anything we want to do.

call testIt

#ifdef SPARSE
call MPI_Finalize(ierr)
#endif

stop "Finished testing"
#endif


!    establish I/O unit numbers and open in- and output files
!    and read command line arguments
  call setio()

#ifdef SPARSE
 call PETSc_openlogfile()

!write(*,*) "Initializing petsc"
!  first, establish PETSc - MPI communicators
   call PetscInit()
!write(*,*) "Initialized petsc"

#ifdef TIMING
!      call PetscGetTime(time_0, pierr)
!      call PetscGetTime(time_start, pierr)
#endif
#endif


!    establish I/O unit numbers and open in- and output files
!    call setio()

call initialiseDebug()

!    establish the basic geometric dimensions and element type
call INTCON ()

if (getrank().eq.0) then
!           write(*,*) 'finished intcon'
endif


if (iecho.eq.6) then
    call startClock()
endif

!    read in system data, establish the grid and define mechanical problem.
call READFE (NSD,NDOF,NSTR,NEN,ierr)


call MPI_barrier(0,ierr)


if (iecho.eq.8) then
    call ParallelLog("main", "Finished readfe")
endif


if (ierr.ne.0) then
    call xit(1,"Reading TECIN failed")
endif

if (iecho.eq.6) then
    call LogTimeStamp("reading TECIN")
endif


!    read data relating to the thermal problem
!      call RDTHRM (NSD,NDOF,NEN,.true.,ierr)

call RDTHRM (.true.,ierr)


if (iecho.eq.8) then
    call ParallelLog("main", "Finished readthrm")
endif

if (getrank().eq.0) then
    call openf (lustat,statusFileName,'unknown')
    write(lustat,*) 'Initialising'
    close(lustat)
endif

if (ierr.ne.0) then
   call xit(1," ")
endif

#ifdef SPARSE
!#ifdef TIMING
!       call PetscGetTime(time_finish, pierr)
!      write(0,*) 'setrank'

    call MPI_Comm_rank (PETSC_COMM_WORLD,irank,pierr)


    if (irank.eq.0) then
!           time_diff = time_finish - time_start
!           write(*,*) "Time spent in very init phase:  ", time_diff
    endif

!       call PetscGetTime(time_start, pierr)
!#endif
#endif

    NUMNP = meshdatactx%Nvlocal
    NUMEL = meshdatactx%Nelocal
!    write(0,*) 'set NUMNP and NUMEL'
!   after reading all input continue with setting up parallel data
!   structures, petsc and per processor numbering, etc
#ifdef SPARSE


if (iecho.eq.8) then
     call ParallelLog("main", "writing modeldata")
endif


call writemodeldata(meshdatactx, &
                    modeldatactx, &
                    outputcontroldatactx, &
                    timestepdatactx, &
                    luin)

if (iecho.eq.6) then
    call startClock()
endif

if (iecho.eq.8) then
    call ParallelLog("main", "creating equation numbering")
endif

call createEquationNumbering(meshdatactx,modeldatactx,nen,ndof,nsd)

if (iecho.eq.8) then
    write(*,*) "rank", getrank(), "scatters data"
    call ParallelLog("main", "scattering data")
endif

call scatterdata(meshdatactx,modeldatactx,getrank(),FILE_outputf)

if (iecho.eq.8) then
    write(*,*) "rank", getrank(), "scattered data"
    call ParallelLog("main", "scattered data")
endif


if (iecho.eq.6) then
    call LogTimeStamp("Creating parallel numbering")
endif


! we have the read the mesh, and set everything in local numbering
! From here on, we can let each partition search for its ghost elements.
!#ifdef EXPERIMENTAL_GHOSTELEMS 
!call getGhostElements()
!#endif



if (ITMODE.ne.0) then
    do n=1,meshdatactx%Nvlocal!NUMNP
        if (modeldatactx%ITBC(1,n).eq.1) then
            modeldatactx%T(n) = modeldatactx%TBC(n)
        endif
    enddo
endif

!      ndofVecSet = .false.
!      call scatterNSDdata(meshdatactx,modeldatactx,getrank(),
!     .  meshdatactx%X,FILE_outputf)
!      call scatterNDOFdata(meshdatactx,modeldatactx,getrank(),
!     .  modeldatactx%BOND,FILE_outputf)
#endif

if (ITEST.eq.1) then
    if (iecho.eq.8) then
        call ParallelLog("main", "testing elements")
    endif
    call ELMTEST (meshdatactx%X, meshdatactx%IEN)
    if (iecho.eq.8) then
        call ParallelLog("main", "elements tested")
    endif

endif


!    mechanical data check only?
if (MODE.eq.0) then
    goto 200
endif

!    establish the topological information of A
if (MODE.le.2.or.MODE.eq.6) then

    if (iecho.eq.8) then
        call ParallelLog("main", "mode 2 or 6; calling sydata")
    endif

    call SYDATA ()

    if (iecho.eq.8) then
        call ParallelLog("main", "mode 2 or 6; finished sydata")
    endif

endif


#ifdef SPARSE
#ifdef TIMING
!       call PetscGetTime(time_finish, pierr)

    if (iecho.eq.8) then
        call ParallelLog("main", "determining Petsc rank")
    endif

    call MPI_Comm_rank (PETSC_COMM_WORLD,irank,pierr)

    if (iecho.eq.8) then
        call ParallelLog("main", "determined Petsc rank")
    endif


    if (irank.eq.0) then
        time_diff = time_finish - time_start
        write(0,*) "Time spent in sydata phase:     ", time_diff
    endif
!       call PetscGetTime(time_start, pierr)
#endif
#endif


if (iecho.eq.8) then
     call ParallelLog("main", "calling tempdat")
endif

!    establish the diagonal addresses of AA
!    call TEMPDAT (NEN)
call TEMPDAT ()

if (iecho.eq.8) then
     call ParallelLog("main", "finished tempdat")
endif


#ifdef SPARSE
#ifdef TIMING
!       call PetscGetTime(time_finish, pierr)
    call MPI_Comm_rank (PETSC_COMM_WORLD,irank,pierr)
    if (irank.eq.0) then
    time_diff = time_finish - time_start
   write(0,*) "Time spent in tempdt phase:     ", time_diff
    endif
!       call PetscGetTime(time_start, pierr)
#endif
#endif

!    open compacted output file and write file header

if (iecho.eq.8) then
    call ParallelLog("main", "OPNDSK; temporarily removed")
endif


 call OPNDSK(.true.,ierr)


   if (iecho.eq.8) then
        call ParallelLog("main", "OPNDSK done")
   endif


!    thermal data check only ?
 if (MODE.eq.3) go to 200

!    transfer material properties
call ELMAT()


! since elmat has the proper references
!
if (ISTART.lt.0) then
!            optijd = .true.
   if (iecho.eq.8) then
        call ParallelLog("main", "calling RESTIO")
   endif

    call RESTIO (modeldatactx%D, &
     modeldatactx%DXE,modeldatactx%stn,modeldatactx%str,modeldatactx%BTOT, &
       modeldatactx%T,modeldatactx%FTOT,modeldatactx%VPSTR,NDOF,NEN,NSTR,readRestart)

   if (iecho.eq.8) then
        call ParallelLog("main", "completed RESTIO")
   endif


endif

!call CheckArrayForNans(modeldatactx%BTOT,"btot in main")



!
#ifdef SPARSE
#ifdef TIMING
!       call PetscGetTime(time_finish, pierr)
    call MPI_Comm_rank (PETSC_COMM_WORLD,irank,pierr)
    if (irank.eq.0) then
    time_diff = time_finish - time_start
   write(*,*) "Time spent in setup phase:      ", time_diff
    endif
!       call PetscGetTime(time_start, pierr)
#endif
#endif


    if (irank.eq.0) then
!           write(*,*) 'starting to compute steady state solution '
    endif

   if (iecho.eq.6) then
  call startClock()
   endif


    if (iecho.eq.8) then
        call ParallelLog("gtecton", "Computing steady state solution")
    endif



! before we start to solve, make room by clearing



!    construct the initial solution to the thermal problem.
 call STEADY ()


    if (iecho.eq.8) then
        call ParallelLog("gtecton", "Finished Computing steady state solution")
    endif


   if (iecho.eq.6) then
  call LogTimeStamp("Computing steady state thermal solution")
   endif



!    thermal rank check only?
 if (MODE.eq.4) then
    call PRINTP (modeldatactx%AA,itdiag, modeldatactx%IDT,modeldatactx%IDTX,1,1)
    goto 200
 endif


if (iecho.eq.5 .or. iecho.eq.55) then
    write(debugFileID,*) "locind  : ", meshdatactx%locind
    write(debugFileID,*) "gloind  : ", meshdatactx%gloind
    write(debugFileID,*) "locelt  : ", meshdatactx%locelt
    write(debugFileID,*) "gloelt  : ", meshdatactx%gloelt
    write(debugFileID,*) "ID      : ", modeldatactx%ID
    write(debugFileID,*) "IDX     : ", modeldatactx%IDX
#ifdef SPARSE
    write(debugFileID,*) "vertices: ", vertices
    write(debugFileID,*) "vertmask: ", verticesmask    
#endif
endif


if (getrank().eq.0) then
    call openf (lustat,statusFileName,'unknown')
    write(lustat,*) 'Busy: 0'
    close(lustat)
endif

! construct the elastic solution to the mechanical fe problem.

if (iecho.eq.8) then
    call ParallelLog("main", "computing elastic solution")
endif

call INITAL()

!    mechanical rank check only?
if (MODE.eq.1) then
    call PRINTP (modeldatactx%A,IDIAG, modeldatactx%ID,modeldatactx%IDX,NDOF,0)
    goto 200
endif

if (iecho.eq.8) then
    call ParallelLog("gtecton", "Calling timsol.")
endif

!    perform time dependent integration.
if (NINTG.gt.0) then
    call TIMSOL ()
endif


if (getrank().eq.0) then
    call openf (lustat,statusFileName,'unknown')
    write(lustat,*) 'Done'
    close(lustat)
endif

if (iecho.eq.8) then
    call ParallelLog("gtecton", "Finished timsol.")
endif


#ifdef SPARSE
#ifdef TIMING
!       call PetscGetTime(time_finish, pierr)
call MPI_Comm_rank (PETSC_COMM_WORLD,irank,pierr)
if (irank.eq.0) then
    time_diff = time_finish - time_start
    write(*,*) "Time spent in solution fase 2: ", time_diff
    write(*,*) "Time spent in itsolve routine: ", time_soltot
    time_diff = time_finish - time_0
    write(*,*) "Time spent in gtecton (total): ", time_diff
endif
#endif
#endif

if (iecho.eq.8) then
    call ParallelLog("gtecton", "Finished sim timing.")
endif


if (debug) then
    close(FILE_outputf)
endif

#ifdef SPARSE
   call PETSc_closelogfile()
#endif

#ifdef SPARSE
if (getrank().eq.0) then
!    ! we do not want all partitions to write to the same file
    call writeGTectonrc()
endif
#endif



200 call date_and_time(values=datetime)

if (irank.eq.0) then
    if (datetime(2).eq.12) then
        if (datetime(3).gt.20 .and. datetime(3).le.28) then
            call MerryChristmas()
        endif
    else if (datetime(2).eq.1) then
        if (datetime(3).eq.1) then
            call HappyNewYear()
        endif
    endif
endif


if (iecho.eq.8) then
    call ParallelLog("gtecton", "Waves goodbye.")
endif

call xit(0,"Simulation completed succesfully.")
end program
