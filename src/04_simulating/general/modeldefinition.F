MODULE MODELDEFINITION
implicit none
!     this module contains variables related to the definition and
!     nature of the model. it does not contain data related to nodes,
!     boundary conditions, etc, but contains related to the proces such
!     as the variables MODE, ITMODE

integer :: NSUR

integer :: MODE,ICVIS,ISTART,NELAST

integer :: IWORD   ! 1 for single prec; 2 for double prec.
integer :: MTOT    ! total memory allowed; set in size.i
integer :: NUMTR   ! Number of triangles
integer :: NUMNP   ! HORRIBLE, should be replaced by meshdatactx%nv[local/global]
integer :: NEQ     ! HORRIBLE, should be replaced by meshdatactx%NEQ[local/global]
                 ! see below 'some auxiliary variables'
integer :: NSIZEA,NEE,NEC,NUMROT, NUMSN,NUMSE,NWINK,NWINKX,NLINK,NMAXW

integer :: NTYPE
integer :: NUMEL  ! HORRIBLE, should be replaced by meshdatactx%ne[local/global]
integer :: NUMAT  ! HORRIBLE, should be replaced by a variable in materials module
integer :: IOPT,IPOINT,NUMPR,NUMSLP, &
           NUMFN,NPRE,LGDEF,IRESDU,NUMSTR,IGRAV,IVLIM, &
           NUMWNK,NSURF,NSED,INCOMP,NSLSKEW
#ifdef SHEARZONE_TRACTIONS
integer :: NUMELSHEARZONE
integer :: NUMNPSHEARZONE
#endif

logical :: feplt

integer :: LinPL

integer, parameter :: NOCOMPR = 0

integer :: ISHELL
double precision :: RADIUS

integer :: ITMODE, IADVEC, IDIFT

integer :: MAATOT,NTEQ,NAA,MBTOT

! NTANOM = Number of thermal anomalies
integer :: NFLX,NTANOM,ITPRES,ITECT,NTWINK

! some auxiliary variables
integer :: NEQlocal, NEQglobal
integer :: NTEQlocal, NTEQglobal
integer :: NEQext
integer :: NTEQext

! this are the global numbers of pressure bcs etc as read from the
! input files. the variables above, like NUMPR, denote the local
! number of pressure bcs
integer :: NUMPRglobal    ! pressure?
integer :: NUMSLPglobal   ! slippery nodes?
integer :: NUMFNglobal    ! faulted nodes
integer :: NPREglobal     ! prestress
integer :: NUMSTRglobal   ! ?
integer :: NUMWNKglobal   ! Winkler... pressure of force?
integer :: NSURFglobal    ! surface nodes
integer :: NSLSKEWglobal  ! fault parallel surface

#ifdef EXPERIMENTAL_FRICTION
! integer :: nFrictionEntriesGlobal  ! to be deletec when obsolete
integer :: nFaultSegmentEntriesGlobal
#endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer :: nElemsWithSlipperyNodes
integer :: nLocalElemsWithSLipperyNodes
integer, allocatable :: nLocalElemsWithSLipperyNodesOfEachPartition(:)
#endif


! number of slippery nodes/elements
integer :: NUMSNglobal ! nodes
integer :: NUMSEglobal ! elements

integer :: NFLXglobal, NTANOMglobal, NTWINKglobal    

! and their local equivalents, containing the amount for onl this partition.
! For runs on a single partition, these will have the same value
integer :: NUMPRlocal   
integer :: NUMSTRlocal  
integer :: NUMWNKlocal 
integer :: NUMSLPlocal 
integer :: NSLSKEWlocal
integer :: NUMFNlocal  
integer :: NPRElocal   
integer :: NSURFlocal

integer :: simVariety ! is set to isPln, or isOpn, or so, from constants module


logical :: allTemperaturesFixed


! set to 1 when pln uses cyclic unlocking of differential Winkler forces
! and set to 0 when not.
#ifdef EXPERIMENTAL_CYCLIC
logical :: UseCyclicSlip
#endif

!********* TECIN TYPE, which is used to store the data of the TECIN.xml *******
!**** not yet used; for future development ***************

type type_ICVIS_flags
    logical :: powerlaw
    logical :: viscoplastic
    logical :: viscous_weakening
    logical :: visc__true_strain__false_strain_rate
    logical :: cohesion_weakening
    logical :: coh__true_strain__false_strain_rate
end type

type type_tecin_prelims
    character(len=255)      :: model_name
    integer                 :: n_points
    integer                 :: n_elements
    integer                 :: mode
    integer                 :: n_intg
    integer                 :: nmprt
    integer                 :: max_iter
    integer                 :: n_link
    integer                 :: i_start
    integer                 :: ne_last
    integer                 :: npr_mat
    integer                 :: FaultVelo
    type (type_ICVIS_flags) :: ICVIS_flags
    character(len=255)      :: output_times
    character(len=255)      :: mat_diag_output_time
end type

type type_nodal_bcs
    character(len=255) :: filename_nodal_coords
    character(len=255) :: filename_nodal_boundary_ID
    character(len=255) :: filename_nodal_boundary_magn
    character(len=255) :: filename_linked_coordinates
    character(len=255) :: filename_nodal_winkler_forces_ID
    character(len=255) :: filename_nodal_winkler_forces_magn
    character(len=255) :: filename_nodal_maxwell_forces
    character(len=255) :: filename_local_DOF_rotations
end type

type type_time_steps
    integer       :: n_time_steps
    integer       :: time_step_size
    character(len=255) :: time_step_unit
    integer       :: mech_alpha
end type

type type_tecin_switches  
    integer       :: iopt
    integer       :: ipoint
    integer       :: lgdef
    integer       :: iresdu
    integer       :: igrav
    integer       :: ivlim
    integer       :: incomp
    integer       :: nocompr
    integer       :: nsed
    integer       :: ishell
end type

type type_element_property
    integer       :: n_lines
    character(len=255) :: filename
end type

type type_element_props
    type (type_element_property) :: elements
    type (type_element_property) :: surfaces
    type (type_element_property) :: pressure_bcs
    type (type_element_property) :: stress_bcs
    type (type_element_property) :: winkler_restoring
    type (type_element_property) :: slippery
    type (type_element_property) :: fault_parallel_elems
    type (type_element_property) :: faulted_nodes
    type (type_element_property) :: initial_stresses
    type (type_element_property) :: triple_junctions
end type

   type type_materials
  integer          :: index
  integer          :: nr_increment
  double precision :: young
  double precision :: poisson
  double precision :: viscosity
  double precision :: power
  double precision :: density     
  double precision :: thickness
  character(len=255)    :: description
   end type

   type type_powerlaw
  integer          :: index 
  integer          :: nr_increment
  double precision :: qpl
  double precision :: apl
  double precision :: npl
  double precision :: vmin
  double precision :: vmax
  double precision :: vpl
  character(len=255)    :: description
   end type

   type type_visco_elasticity
  integer          :: index 
  integer          :: nr_increment
  double precision :: np_type
  double precision :: n_flow
  double precision :: frangl
  double precision :: unit_id
  double precision :: fluidy
  double precision :: yieldpower
  double precision :: stnhrd
  character(len=255)    :: description
   end type

   type type_weakening
  integer          :: index
  integer          :: nr_increment
  double precision :: v_weak
  double precision :: v_sigma
  double precision :: c_weak
  double precision :: c_sigma
  character(len=255)    :: description
   end type

   type type_sediment   
  double precision :: rho_sediment
  double precision :: k_erosion
  double precision :: radius
  double precision :: grav_x, grav_y, grav_z
   end type

   type type_tecin
  type (type_tecin_prelims)      :: tecin_prelims
  type (type_nodal_bcs)          :: nodal_bcs
  type (type_time_steps)         :: time_steps
  type (type_tecin_switches)     :: tecin_switches
  type (type_element_props)      :: element_props
  type (type_materials)          :: materials
  type (type_powerlaw)           :: powerlaw
  type (type_visco_elasticity)   :: visco_elasticity
  type (type_weakening)          :: weakening
  type (type_sediment)           :: sediment
   end type

   END MODULE MODELDEFINITION


