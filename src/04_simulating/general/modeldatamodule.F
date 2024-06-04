MODULE MODELDATAMODULE

USE MODELDEFINITION, only: NEQglobal

implicit none
#define PETSC_AVOID_DECLARATIONS
#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscao.h"
!#include "petsc/finclude/petscis.h"
#undef PETSC_AVOID_DECLARATIONS

integer, allocatable :: IDIAG(:)
integer, allocatable :: ITDIAG(:)


   ! modeldata type contains arrays associated with the model
TYPE modeldata
!#define PETSC_AVOID_DECLARATIONS
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscvec.h"
!#undef PETSC_AVOID_DECLARATIONS
!          PetscReal, ALLOCATABLE :: disp(:)
!         nodal displacements
  PetscReal, ALLOCATABLE :: D(:,:)
#ifdef EXPERIMENTAL_FRICTION
double precision, allocatable :: cumulativeDisp(:,:)
double precision, allocatable :: cumulativeDX(:,:)
#endif

!         differential displacement (DELta D)  (for velocity)
  PetscReal, ALLOCATABLE :: DELD(:,:)
!         split nodes displacement... or slippery nodes?
  PetscReal, ALLOCATABLE :: DX(:,:) ! differential displacement
  PetscReal, ALLOCATABLE :: DXE(:,:,:) ! differential displacement per element
  PetscReal, ALLOCATABLE :: DELX(:,:,:) ! differential velocity
  PetscReal, ALLOCATABLE :: DELDXE(:,:,:) ! differential velocity per element

#ifdef EXPERIMENTAL_ReducedSlipAlloc
  ! look at position local index, 
  ! get position in the DXE,DELX,DELDXE,LMX array
  PetscInt,  ALLOCATABLE :: SlipElemLookup(:)
#endif

!         nodal temperatures
  PetscReal, ALLOCATABLE :: AA(:)
  PetscReal, ALLOCATABLE :: T(:)
  PetscReal, ALLOCATABLE :: T0(:)
  PetscReal, ALLOCATABLE :: TLST(:)
  PetscReal, ALLOCATABLE :: TX(:)
  PetscReal, ALLOCATABLE :: TLSTX(:)
!         fluxes
  PetscReal, ALLOCATABLE :: FLUX(:,:)
  PetscReal, ALLOCATABLE :: ThermalFlux(:,:)

  PetscReal, ALLOCATABLE :: DMAT(:,:,:)

!         A-vector used for non-sparse solving.
  PetscReal, ALLOCATABLE :: A(:)

!         load vector
  PetscReal, ALLOCATABLE :: B(:)
!         residual after the solve
  PetscReal, ALLOCATABLE :: BRES(:)
!         total displacement
  PetscReal, ALLOCATABLE :: BTOT(:)
!         ?
  PetscReal, ALLOCATABLE :: BT(:)


! s does something with boundary conditions.
  PetscReal, ALLOCATABLE :: S(:)


!         stress
  PetscReal, ALLOCATABLE :: STN(:,:)
!         strain
  PetscReal, ALLOCATABLE :: STR(:,:)   ! xx, yy, xy  in 2D
  PetscReal, ALLOCATABLE :: VPSTR(:,:)
  double precision, ALLOCATABLE :: BETA(:,:)
  PetscReal, ALLOCATABLE :: BETB(:,:)


#ifdef EXPERIMENTALTRANSIENTDISCREEP
! taylor stress
  PetscReal, ALLOCATABLE :: STNTAY(:,:)
! taylor stress time derivative
  PetscReal, ALLOCATABLE :: STNTAYRATE(:,:)
! threshold stress, scalar
  PetscReal, ALLOCATABLE :: STNTHRES(:)
! dislocation creep effective strain rate (scalar)
  PetscReal, ALLOCATABLE :: RATEDISLEFF(:) 
! dislocation creep effective inverse viscosity (scalar)
  PetscReal, ALLOCATABLE :: INVEFFVISCDISL(:) 
! diffusion creep effective inverse viscosity (scalar)
  PetscReal, ALLOCATABLE :: INVEFFVISCDIFF(:) 

#endif  


#ifndef EXPERIMENTAL_NORATESTORAGE
!         strain rate
  PetscReal, ALLOCATABLE :: RATE(:,:)
#endif

!         global equation ID of each degree of freedom; mechnical
!         ID(i,j) is the equation number of the i-th d.o.f. of vertex j
  PetscInt, ALLOCATABLE :: ID(:,:)

!         extra equations numbers for slippery nodes, for each degree of freedom
  PetscInt, ALLOCATABLE :: IDX(:,:)

!         The equation numbers of the neighbouring vertices
!         AjdID(i,j,k) has equation number belonging to 
!         D.O.F.  i, of neighbor j, of vertex k   
!         (max neighbours is 50, by allocation. In 3D, sometimes it reaches slighty over 30, with unlucky meshes)
!  PetscInt, pointer :: AdjID(:,:,:)
!  PetscInt, pointer :: AdjIDX(:,:,:)
!  PetscInt, pointer :: IDglobal(:,:)
!  PetscInt, pointer :: IDXglobal(:,:)

  PetscInt, ALLOCATABLE :: AdjID(:,:,:)
  PetscInt, allocatable :: AdjIDX(:,:,:)
  PetscInt, allocatable :: IDglobal(:,:)
  PetscInt, allocatable :: IDXglobal(:,:)



!          PetscInt NEQlocal
!          PetscInt NEQglobal
!         equation numbers per dof per node per element
!         LM(i,j,k) is the eq nr of d.o.f. i of vertex j of element k
!         actual allocation takes place in the AOmodule.F
  PetscInt, ALLOCATABLE :: LM(:,:,:)

!         equation numbers per dof per faulted node per element
  PetscInt, ALLOCATABLE :: LMF(:,:,:)

!         equation numbers per dof per slippery node per element
  PetscInt, ALLOCATABLE :: LMX(:,:,:)

!         equation numbers per node per element thermal problem
  PetscInt, ALLOCATABLE :: LMT(:,:)
  PetscInt, ALLOCATABLE :: LMTX(:,:)

!         global equation of each degree of freedom thermal
  PetscInt, ALLOCATABLE :: IDT(:)
  PetscInt, ALLOCATABLE :: IDTX(:)
  PetscInt, pointer :: AdjIDT(:,:)
  PetscInt, pointer :: AdjIDTX(:,:)
  PetscInt, ALLOCATABLE :: &
                          IDTglobal(:)
  PetscInt, ALLOCATABLE :: &
                          IDTXglobal(:)
!          PetscInt :: NTEQlocal
!          PetscInt :: NTEQglobal
!         boundary condition identifiers and values mech
  PetscInt, ALLOCATABLE :: &
                      IBONDGLOBAL(:,:)
  PetscInt, ALLOCATABLE :: IBOND(:,:)
  PetscReal, ALLOCATABLE :: BOND(:,:)
!         skew array for rotated nodes
  PetscReal, ALLOCATABLE :: SKEW(:,:)
  PetscReal, ALLOCATABLE :: OSKEW(:,:)
!          PetscInt :: NUMROT
!         boundary condition identifiers and values thermal
  PetscInt, ALLOCATABLE :: &
                          ITBCGLOBAL(:)
  ! ITBC is the thermal bounday cond per point,
  ! and should only one-dimensional.
  ! However, because all the integer boundary conditions
  ! are read by the same routine (IGEN; in tecin.f)
  ! a two-dimensional array is expected.
  ! First dimension is allocated to 1.
  PetscInt, ALLOCATABLE :: ITBC(:,:)
  PetscReal, ALLOCATABLE :: TBC(:)
!         linked nodes
  PetscInt, ALLOCATABLE :: LINK(:,:)
  PetscInt, ALLOCATABLE :: LINKID(:,:)
!          PetscInt :: NLINK
!         faulted (split) nodes
  PetscInt, ALLOCATABLE :: &
         NFAULTglobal(:,:)
  PetscInt, ALLOCATABLE :: NFAULT(:,:)
  PetscReal, ALLOCATABLE :: FAULT(:,:)
  PetscReal, ALLOCATABLE :: DFAULT(:,:) ! delta fault; incremental
  PetscReal, ALLOCATABLE :: TFAULT(:,:) ! total fault
  PetscReal, ALLOCATABLE :: FAULTL(:,:)

!          PetscInt :: NUMFN
!         slippery nodes
  ! temporary global array needed to set equations numbers
  ! correctly after input has been read
  PetscInt,  ALLOCATABLE :: NSLIPglobal(:,:)
  PetscInt,  ALLOCATABLE :: NSLIP(:,:)
  PetscReal, ALLOCATABLE :: DIFORC(:,:)
  PetscInt,  ALLOCATABLE :: IDSLP(:)
  PetscInt,  ALLOCATABLE :: IDSLE(:)
#ifdef EXPERIMENTAL_FRICTION
integer, allocatable :: maxDiffForceMagnitude(:)
#endif
!          PetscInt :: NUMSLP
!          PetscInt :: NUMSLPglobal
!          PetscInt :: NUMSN, NUMSE
!          PetscInt :: NUMSNglobal, NUMSEglobal
  ! fault parallel slippery element (experimental)
  PetscInt, ALLOCATABLE :: NSELSD(:,:)
!          PetscInt :: NSLSKEW
  ! differential winkler forces,
  PetscInt, ALLOCATABLE :: IWINX(:,:)
  PetscReal, ALLOCATABLE :: WINX(:,:)
!          PetscInt :: NWINKX
!         winklers forces on nodes
  PetscInt, ALLOCATABLE :: IWINK(:,:)
  PetscReal, ALLOCATABLE :: WINK(:,:)
!          PetscInt :: NWINK
!         maxwells
  PetscReal, ALLOCATABLE :: FMAXW(:,:)
  PetscReal, ALLOCATABLE :: FTOT(:,:)
!          PetscInt :: NMAXW
!         pre-stresses
  ! contains the stress tensor entries
  PetscReal, ALLOCATABLE :: STN0(:,:)
  ! element number for which this is the stress
  PetscInt, ALLOCATABLE :: ISELM(:)
  PetscInt, ALLOCATABLE :: ISTIME(:)
!          PetscInt :: NPRE
!         element pressure
  PetscReal, ALLOCATABLE :: PRES(:)
  PetscInt, ALLOCATABLE :: IELNO(:)
  PetscInt, ALLOCATABLE :: ISIDE(:)
!          PetscInt :: NUMPR
!         stress bc
  PetscReal, ALLOCATABLE :: STRS(:,:)
  PetscInt, ALLOCATABLE :: IELSTR(:)
  PetscInt, ALLOCATABLE :: ISSIDE(:)
  PetscInt, ALLOCATABLE :: ISTR(:,:)
!          PetscInt :: NUMSTR
!         winkler restoring pressures
  PetscReal, ALLOCATABLE :: WPRES(:)
  PetscInt, ALLOCATABLE :: IWELM(:)
  PetscInt, ALLOCATABLE :: IWSIDE(:)
  PetscInt, ALLOCATABLE :: IWTIME(:)
!         shear zone tractions
  PetscInt, ALLOCATABLE :: ELNOSHEARZONE(:)
  PetscInt, ALLOCATABLE :: NPNOSHEARZONE(:)
  PetscInt, ALLOCATABLE :: SIDESHEARZONE(:)
  PetscReal, ALLOCATABLE :: DFMAXSHEARZONE(:)
  PetscReal, ALLOCATABLE :: shearETA(:)
  PetscReal, ALLOCATABLE :: shearTHICK(:)
  PetscReal, ALLOCATABLE :: shearTEMP(:)
  PetscInt shearnpl
  PetscReal, ALLOCATABLE :: shearQpl(:)
  PetscReal, ALLOCATABLE :: shearF(:) ! To replace B (containing only shear forces)
!          PetscInt :: NUMWNK
!         surface nodes
  PetscInt, ALLOCATABLE :: &
                     ISURFORDER(:)
  PetscInt, ALLOCATABLE :: ISURF(:)
  PetscReal, ALLOCATABLE :: SURF(:,:)
  PetscReal, ALLOCATABLE :: SURF0(:,:)
!          PetscInt :: NSURF
!         heat flux boundary conditions
  PetscInt, ALLOCATABLE :: IFLX(:)
  PetscInt, ALLOCATABLE :: IFLS(:)
  PetscReal, ALLOCATABLE :: BFLX(:)
!          PetscInt :: NFLX
!         thermal anomaly
  PetscInt, ALLOCATABLE :: ITANOM(:,:)
  PetscReal, ALLOCATABLE :: TANOM(:)
!          PetscInt :: NTANOM
!         thermal winkler flux boundary conditions
  PetscInt, ALLOCATABLE :: ITWINK(:,:)
  PetscReal, ALLOCATABLE :: TWINK(:)
!          PetscInt :: NTWINK

  PetscReal, allocatable :: &
                    displacement(:)

    integer, allocatable :: checkSignDuality(:)

! two mystery arrays for showfbc
    integer, allocatable :: forShowfbcIAeqnr(:,:)
    integer, allocatable :: forShowfbcIAXeqnr(:,:)


! for plt 3d, to check the sum [BPmatrix].n
! to check the error caused mby the mesh size.
  double precision, allocatable :: elemError(:,:)


#ifdef EXPERIMENTAL_FRICTION
! the shear strength on a fault, global array
    double precision, allocatable :: DIFFSTRENGTH(:)
#endif

END TYPE modeldata

! used to be the showfbc common block.
logical :: showfbc
integer :: MEQ



type (modeldata) :: modeldatactx

PetscInt, ALLOCATABLE :: equationreference(:)

type plotdata

   double precision, allocatable ::   XPT(:,:)
   integer,          allocatable ::   IPT(:,:)
   double precision, allocatable ::  SRF1(:,:)
   double precision, allocatable ::  XINP(:,:)
   integer,          allocatable :: IJSUR(:,:)

end type plotdata

type(plotdata) :: plotdatactx


contains

!-----mechmodeldata------------------------------------------
!--------------------------------------------------------------------
subroutine allocatePlotData(numnp, numel, nsurface, nsd)
!--------------------------------------------------------------------
use DebugModule,   only: AllocateError
use modeltopology, only: nen
!--------------------------------------------------------------------
implicit none

integer :: numnp, numel, nsurface, nsd
integer, parameter :: nplt=3

integer :: ERROR


allocate(plotdatactx%XPT(3         , nplt*nen*numel )  , STAT=ERROR)
allocate(plotdatactx%IPT(nen       , numel)  , STAT=ERROR)
allocate(plotdatactx%SRF1(nsurface , nsd  )  , STAT=ERROR)
allocate(plotdatactx%XINP(numnp    , nsd  )  , STAT=ERROR)
allocate(plotdatactx%IJSUR(nsurface, numel)  , STAT=ERROR)

plotdatactx%XPT = 0d0
plotdatactx%IPT = 0
plotdatactx%SRF1 = 0d0
plotdatactx%XINP = 0d0
plotdatactx%IJSUR = 0

if (ERROR.ne.0) then
    call AllocateError("modeldatamodule -> plotdata", ERROR)
endif

end subroutine allocatePlotData
!--------------------------------------------------------------------

!-----mechmodeldata------------------------------------------
!--------------------------------------------------------------------
subroutine allocatemechmodeldata()

use debugmodule,    only: xit
use meshdatamodule, only: meshdatactx
use modeltopology,  only: ndof
use iomodule
use debugmodule
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
PetscInt :: ERROR
!--------------------------------------------------------------------
if (allocated(modeldatactx%IDglobal)) then
    deallocate(modeldatactx%IDglobal)
endif
ALLOCATE(modeldatactx%IDglobal(ndof, meshdatactx%nvglobal),  STAT=ERROR)
if (ERROR.ne.0) then
   write(stderr,1) ERROR,'IDglobal'
1  format(' **FATAL**: allocatemechmodeldata: error ',i8, &
          'while allocating memory for ',a50)
   call xit(1," ")
endif
!if (debug) write(*,*) 'Allocated IDglobal(nvglobal=',meshdatactx%nvglobal,')'

if (allocated(modeldatactx%ID)) then
    deallocate(modeldatactx%ID)
endif
! should actually be local. TODO,, fix this
ALLOCATE(modeldatactx%ID(ndof, meshdatactx%nvglobal),  STAT=ERROR)
!   ALLOCATE(modeldat%ID(ndof, meshdatactx%nvlocal),  STAT=ERROR)
if (ERROR.ne.0) then
   write(stderr,1) ERROR,'ID'
   call xit(1," ")
endif
if (debug) write(*,*) '**WASTEFUL?** Allocated ID(ndof=',ndof,',nvglobal=',meshdatactx%nvglobal,')'

if (allocated(modeldatactx%D)) then
    deallocate(modeldatactx%D)
endif
ALLOCATE(modeldatactx%D(ndof, meshdatactx%nvglobal),         STAT=ERROR)
if (ERROR.ne.0) then
   write(stderr,1) ERROR,'D'
   call xit(1," ")
endif
if (debug) write(*,*) '**WASTEFUL?** Allocated D(ndof=',ndof,',nvglobal=',meshdatactx%nvglobal,')'


#ifdef EXPERIMENTAL_FRICTION
if (allocated(modeldatactx%cumulativeDisp)) then
    deallocate(modeldatactx%cumulativeDisp)
endif
ALLOCATE(modeldatactx%cumulativeDisp(ndof, meshdatactx%nvglobal), &
         STAT=ERROR)
if (ERROR.ne.0) then
   write(stderr,1) ERROR,'cumulativeDisp'
   call xit(1," ")
endif
modeldatactx%cumulativeDisp = 0d0
if (debug) write(*,*) '**WASTEFUL?** Allocated cumulativeDisp(ndof=',ndof,',nvglobal=',meshdatactx%nvglobal,')'

if (allocated(modeldatactx%cumulativeDX)) then
    deallocate(modeldatactx%cumulativeDX)
endif
ALLOCATE(modeldatactx%cumulativeDX(ndof,meshdatactx%nvglobal), &
         STAT=ERROR)
if (ERROR.ne.0) then
   write(stderr,1) ERROR,'cumulativeDX'
   call xit(1," ")
endif
modeldatactx%cumulativeDX = 0d0
if (debug) write(*,*) '**WASTEFUL?** Allocated cumulativeDX(ndof=',ndof,',nvglobal=',meshdatactx%nvglobal,')'

#endif

if (allocated(modeldatactx%DELD)) then
    deallocate(modeldatactx%DELD)
endif
ALLOCATE(modeldatactx%DELD(ndof, meshdatactx%nvglobal),      STAT=ERROR)
if (ERROR.ne.0) then
   write(stderr,1) ERROR,'DELD'
   call xit(1," ")
endif
if (debug) write(*,*) '**WASTEFUL?** Allocated DELD(ndof=',ndof,',nvglobal=',meshdatactx%nvglobal,')'

!       if (allocated(modeldatactx%displacement)) then
!           deallocate(modeldatactx%displacement)
!       endif
!       ALLOCATE(modeldatactx%displacement(NEQglobal),STAT=ERROR)
!       if (ERROR.ne.0) then
!          write(stderr,1) ERROR,'displacement'
!          call xit(1," ")
!       endif
!       if (debug) write(*,*) '**WASTEFUL?** Allocated displacement(NEQglobal=',NEQglobal,')'
!       modeldatactx%displacement = 0d0

modeldatactx%IDglobal     = 0
modeldatactx%ID           = 0
modeldatactx%D            = 0d0
modeldatactx%DELD         = 0d0

end subroutine allocatemechmodeldata
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine deallocatemechmodeldata(modeldat)

use debugmodule, only: xit
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
PetscInt :: ERROR
!--------------------------------------------------------------------

error = 0

if (allocated(modeldat%ID)) then
    deallocate(modeldat%ID, stat=error)
endif

if (allocated(modeldat%D)) then
    deallocate(modeldat%D, stat=error)
endif

if (error.ne.0) then
    write(*,*) 'Could not allocate mech model data. Error ', ERROR
    call xit(1," ")
endif

end subroutine deallocatemechmodeldata
!--------------------------------------------------------------------



!-----initial stresses bcs data--------------------------------------------
!--------------------------------------------------------------------
subroutine allocateinitstress(modeldat,nstr,npre)
use debugmodule, only: xit
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: nstr,npre,error
!--------------------------------------------------------------------
error = 0

if (allocated(modeldat%STN0))then
    deallocate(modeldat%STN0)
endif
ALLOCATE(modeldat%STN0(NSTR,NPRE), STAT=ERROR)

if (allocated(modeldat%ISELM))then
    deallocate(modeldat%ISELM)
endif
ALLOCATE(modeldat%ISELM(NPRE), STAT=ERROR)

if (allocated(modeldat%ISTIME))then
    deallocate(modeldat%ISTIME)
endif
ALLOCATE(modeldat%ISTIME(NPRE), STAT=ERROR)

if (ERROR.ne.0) then
    write(*,*) 'Could not allocate ini stress. Error ', ERROR
    call xit(1," ")
endif

modeldat%STN0   = 0d0
modeldat%ISELM  = 0
modeldat%ISTIME = 0

end subroutine allocateinitstress
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine deallocateinitstress(modeldat)
use debugmodule, only: xit
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ERROR
!--------------------------------------------------------------------
error = 0

if (allocated(modeldat%STN0))then
   DEALLOCATE(modeldat%STN0, STAT=ERROR)
endif

if (allocated(modeldat%ISELM))then
   DEALLOCATE(modeldat%ISELM, STAT=ERROR)
endif

if (allocated(modeldat%ISTIME))then
   DEALLOCATE(modeldat%ISTIME, STAT=ERROR)
endif

if (ERROR.ne.0) then
    write(*,*) 'Could not deallocate ini stress. Error ', ERROR
    call xit(1," ")
endif

!--------------------------------------------------------------------
end subroutine deallocateinitstress
!--------------------------------------------------------------------

!-----pressure boundary conditions-----------------------------------------
!--------------------------------------------------------------------
subroutine allocatepressurebcs(modeldat,numpr)

use debugmodule, only: xit
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: numpr,ERROR
!--------------------------------------------------------------------

if (allocated(modeldat%IELNO))then
    deallocate(modeldat%IELNO)
endif
ALLOCATE(modeldat%IELNO(NUMPR), STAT=ERROR)

if (allocated(modeldat%ISIDE))then
    deallocate(modeldat%ISIDE)
endif
ALLOCATE(modeldat%ISIDE(NUMPR), STAT=ERROR)

if (allocated(modeldat%PRES))then
    deallocate(modeldat%PRES)
endif
ALLOCATE(modeldat%PRES(NUMPR), STAT=ERROR)

if (ERROR.ne.0) then
    write(*,*) 'Could not allocate pressure bc data. Error ', ERROR
    call xit(1," ")
endif

modeldat%IELNO = 0
modeldat%ISIDE = 0
modeldat%PRES  = 0d0

!--------------------------------------------------------------------
end subroutine allocatepressurebcs
!--------------------------------------------------------------------


!--------------------------------------------------------------------
subroutine deallocatepressurebcs(modeldat)
use debugmodule, only: xit
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ERROR
!--------------------------------------------------------------------

error = 0

if (allocated(modeldat%IELNO))then
   DEALLOCATE(modeldat%IELNO, STAT=ERROR)
endif

if (allocated(modeldat%ISIDE))then
   DEALLOCATE(modeldat%ISIDE, STAT=ERROR)
endif

if (allocated(modeldat%PRES))then
   DEALLOCATE(modeldat%PRES, STAT=ERROR)
endif

if (ERROR.ne.0) then
    write(*,*) 'Could not deallocate pressure bc data. Error ', ERROR
    call xit(1," ")
endif


!--------------------------------------------------------------------
end subroutine deallocatepressurebcs
!--------------------------------------------------------------------


!-----stress boundary conditions-------------------------------------------
!--------------------------------------------------------------------
subroutine allocatestressbcs(modeldat,nstr,numstr)
use debugmodule, only: xit

!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer nstr,numstr,ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%STRS))then
    deallocate(modeldat%STRS)
endif
! todo perhaps NSTR incorrect here?
!   ALLOCATE(modeldat%STRS(NSTR, NUMSTR), STAT=ERROR)
ALLOCATE(modeldat%STRS(6, NUMSTR), STAT=ERROR) !overwrite nstr, to deal with thickness


if (allocated(modeldat%IELSTR))then
    deallocate(modeldat%IELSTR)
endif
ALLOCATE(modeldat%IELSTR(NUMSTR), STAT=ERROR)

if (allocated(modeldat%ISSIDE))then
    deallocate(modeldat%ISSIDE)
endif
ALLOCATE(modeldat%ISSIDE(NUMSTR), STAT=ERROR)

if (allocated(modeldat%ISTR))then
    deallocate(modeldat%ISTR)
endif
ALLOCATE(modeldat%ISTR(2,NUMSTR), STAT=ERROR)

if (ERROR.ne.0) then
    write(*,*) 'Could not allocate stress bcs. Error ', ERROR
    call xit(1," ")
endif

modeldat%STRS = 0d0
modeldat%IELSTR = 0
modeldat%ISSIDE = 0
modeldat%ISTR = 0

!--------------------------------------------------------------------
end subroutine allocatestressbcs
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine deallocatestressbcs(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer ERROR
!--------------------------------------------------------------------

if (allocated(modeldat%STRS))then
   DEALLOCATE(modeldat%STRS, STAT=ERROR)
endif

if (allocated(modeldat%IELSTR))then
   DEALLOCATE(modeldat%IELSTR, STAT=ERROR)
endif

if (allocated(modeldat%ISSIDE))then
   DEALLOCATE(modeldat%ISSIDE, STAT=ERROR)
endif

if (allocated(modeldat%ISTR))then
   DEALLOCATE(modeldat%ISTR, STAT=ERROR)
endif

!--------------------------------------------------------------------
end subroutine deallocatestressbcs
!--------------------------------------------------------------------

!-----winkler restoring pressures-------------------------------------------
!--------------------------------------------------------------------
subroutine allocatewinklerrestoring(modeldat,numwnk)

use debugmodule, only: xit

!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer numwnk, ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%IWELM))then
    deallocate(modeldat%IWELM)
endif
ALLOCATE(modeldat%IWELM(NUMWNK), STAT=ERROR)

if (allocated(modeldat%IWSIDE))then
    deallocate(modeldat%IWSIDE)
endif
ALLOCATE(modeldat%IWSIDE(NUMWNK), STAT=ERROR)

if (allocated(modeldat%IWTIME))then
    deallocate(modeldat%IWTIME)
endif
ALLOCATE(modeldat%IWTIME(NUMWNK), STAT=ERROR)

if (allocated(modeldat%WPRES))then
    deallocate(modeldat%WPRES)
endif
ALLOCATE(modeldat%WPRES(NUMWNK), STAT=ERROR)

if (ERROR.ne.0) then
    write(*,*) 'Could not allocate wink press. Error ', ERROR
    call xit(1," ")
endif

modeldat%IWELM  = 0
modeldat%IWSIDE = 0
modeldat%IWTIME = 0
modeldat%WPRES  = 0d0

!--------------------------------------------------------------------
end subroutine allocatewinklerrestoring
!********************************************************************


!********************************************************************
subroutine deallocatewinklerrestoring(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%IWELM))then
   DEALLOCATE(modeldat%IWELM,STAT=ERROR)
endif

if (allocated(modeldat%IWSIDE))then
   DEALLOCATE(modeldat%IWSIDE,STAT=ERROR)
endif

if (allocated(modeldat%IWTIME))then
   DEALLOCATE(modeldat%IWTIME,STAT=ERROR)
endif

if (allocated(modeldat%WPRES))then
   DEALLOCATE(modeldat%WPRES,STAT=ERROR)
endif

!--------------------------------------------------------------------
end subroutine deallocatewinklerrestoring
!********************************************************************

!-----faulted nodes---------------------------------------------------------
!********************************************************************
subroutine allocatefaultednodes(modeldat,ndof,nen,numfn,numel)
!--------------------------------------------------------------------
use debugmodule, only: xit
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ndof,nen,numfn,numel, ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%NFAULTglobal))then
    deallocate(modeldat%NFAULTglobal)
endif
ALLOCATE(modeldat%NFAULTglobal(3,NUMFN), STAT=ERROR)
modeldat%NFAULTglobal = 0

if (allocated(modeldat%NFAULT))then
    deallocate(modeldat%NFAULT)
endif
ALLOCATE(modeldat%NFAULT(3,NUMFN), STAT=ERROR)
modeldat%NFAULT = 0

if (allocated(modeldat%FAULT))then
    deallocate(modeldat%FAULT)
endif
ALLOCATE(modeldat%FAULT(NDOF,NUMFN), STAT=ERROR)
modeldat%FAULT = 0d0

if (allocated(modeldat%DFAULT))then
    deallocate(modeldat%DFAULT)
endif
ALLOCATE(modeldat%DFAULT(NDOF,NUMFN), STAT=ERROR)
modeldat%DFAULT = 0d0

if (allocated(modeldat%TFAULT))then
    deallocate(modeldat%TFAULT)
endif
ALLOCATE(modeldat%TFAULT(NDOF,NUMFN), STAT=ERROR)
modeldat%TFAULT = 0d0

if (allocated(modeldat%FAULTL))then
    deallocate(modeldat%FAULTL)
endif
ALLOCATE(modeldat%FAULTL(NDOF,0), STAT=ERROR)
modeldat%FAULTL = 0d0

if (allocated(modeldat%LMF))then
    deallocate(modeldat%LMF)
endif
ALLOCATE(modeldat%LMF(NDOF,NEN, numel), STAT = ERROR)  ! this is BAAAAAAAD, too large
modeldat%LMF = 0

if (ERROR.ne.0) then
    write(*,*) 'Could not allocate faulted nodes. Error ', ERROR
    call xit(1," ")
endif

end subroutine allocatefaultednodes
!********************************************************************

!********************************************************************
subroutine deallocatefaultednodes(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%NFAULTglobal))then
   DEALLOCATE(modeldat%NFAULTglobal, STAT=ERROR)
endif

if (allocated(modeldat%NFAULT))then
   DEALLOCATE(modeldat%NFAULT, STAT=ERROR)
endif

if (allocated(modeldat%FAULT))then
   DEALLOCATE(modeldat%FAULT, STAT=ERROR)
endif

if (allocated(modeldat%DFAULT))then
   DEALLOCATE(modeldat%DFAULT, STAT=ERROR)
endif

if (allocated(modeldat%TFAULT))then
   DEALLOCATE(modeldat%TFAULT, STAT=ERROR)
endif

if (allocated(modeldat%FAULTL))then
   DEALLOCATE(modeldat%FAULTL, STAT=ERROR)
endif

if (allocated(modeldat%LMF))then
   DEALLOCATE(modeldat%LMF, STAT = ERROR)
endif
end subroutine deallocatefaultednodes
!********************************************************************




!-----slippery nodes--------------------------------------------------------
!********************************************************************
subroutine allocateslipperynodes(modeldat,ndof,nen, &
        numslp,numnpLocal,numnpglobal,numel)

use meshdatamodule,  only: meshdatactx
use modelctx,        only: getrank
use debugmodule,     only: allocateError
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ndof,nen,numslp,numnpLocal,numnpglobal,numel, ERROR
!--------------------------------------------------------------------
! numel is the local nunber of equations.
! numslp from tecin.dat

if (allocated(modeldat%NSLIPglobal))then
    deallocate(modeldat%NSLIPglobal)
endif
ALLOCATE(modeldat%NSLIPglobal(5,NUMSLP), STAT=ERROR)
modeldat%NSLIPglobal = 0

if (allocated(modeldat%NSLIP))then
    deallocate(modeldat%NSLIP)
endif
ALLOCATE(modeldat%NSLIP(5,NUMSLP), STAT=ERROR)
modeldat%NSLIP = 0

if (allocated(modeldat%IDSLP))then
    deallocate(modeldat%IDSLP)
endif
ALLOCATE(modeldat%IDSLP(NUMSLP),STAT=ERROR)
modeldat%IDSLP = 0

if (allocated(modeldat%IDSLE))then
    deallocate(modeldat%IDSLE)
endif
ALLOCATE(modeldat%IDSLE(NUMSLP),STAT=ERROR)
modeldat%IDSLE = 0

!ALLOCATE(modeldat%DIFORC(NDOF,numnp), STAT=ERROR)
if (allocated(modeldat%DIFORC))then
    deallocate(modeldat%DIFORC)
endif
ALLOCATE(modeldat%DIFORC(NDOF,meshdatactx%nvglobal), STAT=ERROR)
! write(*,*) "allocating diforc", NDOF,meshdatactx%nvlocal
!ALLOCATE(modeldat%DIFORC(NDOF,meshdatactx%nvlocal), STAT=ERROR)

modeldat%DIFORC = 0d0

#ifdef EXPERIMENTAL_FRICTION
if (allocated(modeldat%DIFFSTRENGTH))then
    deallocate(modeldat%DIFFSTRENGTH)
end if    
ALLOCATE(modeldat%DIFFSTRENGTH(meshdatactx%nvglobal))
modeldat%DIFFSTRENGTH=0
#endif


if (allocated(modeldat%checkSignDuality))then
    deallocate(modeldat%checkSignDuality)
end if
allocate(modeldat%checkSignDuality(meshdatactx%nvglobal))
modeldat%checkSignDuality = 0


if (allocated(modeldat%IWINX))then
    deallocate(modeldat%IWINX)
endif
ALLOCATE(modeldat%IWINX(NDOF,meshdatactx%nvglobal),STAT=ERROR)

if (allocated(modeldat%WINX))then
    deallocate(modeldat%WINX)
endif
ALLOCATE(modeldat%WINX(NDOF,meshdatactx%nvglobal), STAT=ERROR)

if (allocated(modeldat%DX))then
    deallocate(modeldat%DX)
endif
ALLOCATE(modeldat%DX(NDOF,numnpLocal),STAT=ERROR)
!ALLOCATE(modeldat%DX(NDOF,numnpGlobal),STAT=ERROR)


if (allocated(modeldat%IDX))then
    deallocate(modeldat%IDX)
endif
!ALLOCATE(modeldat%IDX(NDOF,numnpLocal), STAT=ERROR)
ALLOCATE(modeldat%IDX(NDOF,numnpGlobal),STAT=ERROR)


if (allocated(modeldat%IDXglobal))then
    deallocate(modeldat%IDXglobal)
endif
ALLOCATE(modeldat%IDXglobal(NDOF,numnpglobal),STAT=ERROR)

! todo: allocation hell, this. Should be allocated only for the number of slippery nodes and elements.
! FOR SLIPPERY NODES, THIS IS CURRENTLY THE LIMITING FACTOR FOR RESOLUTION on models with slippery nodes.


modeldat%IWINX = 0
modeldat%WINX = 0d0
modeldat%DX = 0d0
modeldat%IDX = 0
modeldat%IDXglobal = 0


call allocateError("Slippery node data", ERROR)

end subroutine allocateslipperynodes
!********************************************************************
#ifdef SHEARZONE_TRACTIONS

subroutine allocateShearTraxArrays()

use debugmodule,     only:  AllocateError
use modeldefinition, only:  NUMELSHEARZONE

implicit none

integer :: ERROR

if (allocated(modeldatactx%ELNOSHEARZONE)) then
    deallocate(modeldatactx%ELNOSHEARZONE)
endif
allocate(modeldatactx%ELNOSHEARZONE(NUMELSHEARZONE), STAT=ERROR)
call AllocateError("ELNOSHEARZONE", ERROR)

if (allocated(modeldatactx%SIDESHEARZONE)) then
    deallocate(modeldatactx%SIDESHEARZONE)
endif
allocate(modeldatactx%SIDESHEARZONE(NUMELSHEARZONE), STAT=ERROR)
call AllocateError("SIDESHEARZONE", ERROR)

if (allocated(modeldatactx%shearETA)) then
    deallocate(modeldatactx%shearETA)
endif
allocate(modeldatactx%shearETA(NUMELSHEARZONE), STAT=ERROR)
call AllocateError("shearETA", ERROR)

if (allocated(modeldatactx%shearTHICK)) then
    deallocate(modeldatactx%shearTHICK)
endif
allocate(modeldatactx%shearTHICK(NUMELSHEARZONE), STAT=ERROR)
call AllocateError("shearTHICK", ERROR)

if (allocated(modeldatactx%shearTEMP)) then
    deallocate(modeldatactx%shearTEMP)
endif
allocate(modeldatactx%shearTEMP(NUMELSHEARZONE), STAT=ERROR)
call AllocateError("shearTEMP", ERROR)

if (allocated(modeldatactx%shearQpl)) then
    deallocate(modeldatactx%shearQpl)
endif
allocate(modeldatactx%shearQpl(NUMELSHEARZONE), STAT=ERROR)
call AllocateError("shearQpl", ERROR)

end subroutine

#endif

   
!********************************************************************
subroutine allocateSlipperyElementArrays()
!--------------------------------------------------------------------

use modelctx,        only: getrank
use modeltopology,   only: ndof, nen
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldefinition, only: nElemsWithSlipperyNodes, &
                           nLocalElemsWithSLipperyNodes
#endif
use meshdatamodule,  only: meshdatactx

implicit none

integer :: ERROR


if (allocated(modeldatactx%LMX))then
    deallocate(modeldatactx%LMX)
endif
#ifdef EXPERIMENTAL_ReducedSlipAlloc
!write(*,*) "a rank", getrank(), "allocates LMX", ndof, nen, nElemsWithSlipperyNodes
ALLOCATE(modeldatactx%LMX(NDOF,NEN,nElemsWithSLipperyNodes), STAT = ERROR)
#else
!write(*,*) "b rank", getrank(), "allocates LMX", meshdatactx%neglobal
ALLOCATE(modeldatactx%LMX(NDOF,NEN,meshdatactx%neglobal), STAT = ERROR)
#endif

if (ERROR .ne. 0) then
    write(0,*) "Could not allocate modeldata LMX; error: ", ERROR
    stop "Exiting..."
endif

if (allocated(modeldatactx%DXE))then
    deallocate(modeldatactx%DXE)
endif
#ifdef EXPERIMENTAL_ReducedSlipAlloc
ALLOCATE(modeldatactx%DXE(NDOF,NEN,nElemsWithSlipperyNodes), STAT = ERROR)
#else
ALLOCATE(modeldatactx%DXE(NDOF,NEN,meshdatactx%neglobal), STAT = ERROR)
#endif


if (ERROR .ne. 0) then
    write(0,*) "Could not allocate modeldata DXE; error: ", ERROR
    stop "Exiting..."
endif

if (allocated(modeldatactx%DELX))then
    deallocate(modeldatactx%DELX)
endif
#ifdef EXPERIMENTAL_ReducedSlipAlloc
ALLOCATE(modeldatactx%DELX(NDOF,NEN,nElemsWithSlipperyNodes), STAT = ERROR)
#else
ALLOCATE(modeldatactx%DELX(NDOF,NEN,meshdatactx%neglobal), STAT = ERROR)
#endif

if (ERROR .ne. 0) then
    write(0,*) "Could not allocate modeldata DELX; error: ", ERROR
    stop "Exiting..."
endif

if (allocated(modeldatactx%DELDXE))then
    deallocate(modeldatactx%DELDXE)
endif
#ifdef EXPERIMENTAL_ReducedSlipAlloc
ALLOCATE(modeldatactx%DELDXE(NDOF,NEN,nElemsWithSlipperyNodes), STAT = ERROR)
#else
ALLOCATE(modeldatactx%DELDXE(NDOF,NEN,meshdatactx%neglobal), STAT = ERROR)
#endif

if (ERROR .ne. 0) then
    write(0,*) "Could not allocate modeldata DELDXE; error: ", ERROR
    stop "Exiting..."
endif


!#ifdef EXPERIMENTAL_ReducedSlipAlloc
!if (allocated(modeldatactx%SlipElemLookup)) then
!    deallocate(modeldatactx%SlipElemLookup)
!endif
!allocate(modeldatactx%SlipElemLookup(meshdatactx%neglobal))
!#endif

!if (ERROR .ne. 0) then
!    write(0,*) "Could not allocate modeldata slipElemLookup; error: ", ERROR
!    stop "Exiting..."
!endif

modeldatactx%LMX = 0
modeldatactx%DXE = 0d0
modeldatactx%DELX = 0d0
modeldatactx%DELDXE = 0d0

!#ifdef EXPERIMENTAL_ReducedSlipAlloc
!modeldatactx%SlipElemLookup = 0
!#endif

end subroutine allocateSlipperyElementArrays
!********************************************************************

#ifdef EXPERIMENTAL_ReducedSlipAlloc
!********************************************************************
integer function elementSlipperyID(ElementID, callerID)

use meshdatamodule,  only: meshdatactx
use modelctx, only: getrank
!--------------------------------------------------------------------
! returns 0 if the element has no slippery nodes
! otherwise, returns the position of this element's slippery
! entries in de LMX,DXE,DELX and DELDXE array
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
integer :: ElementID
integer :: callerID
!--------------------------------------------------------------------

if (ElementID.le.0) then
    write(*,*) "caller", callerID,"request slippery sequence of a non-positive element number", ElementID
    stop "This should not happen; Please contact model support; Leaving GTecton"
else if (ElementID.gt.meshdatactx%neglobal) then
    write(*,*) "caller", callerID,"requests slippery sequence of element", ElementID, &
               "but the mesh has only ", meshdatactx%neglobal, "elements."
    stop "This should not happen; Please contact model support; Leaving GTecton"
else
    elementSlipperyID = modeldatactx%SlipElemLookup(ElementID)
endif

end function
!********************************************************************
#endif


!********************************************************************
subroutine deallocateslipperynodes(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%NSLIPglobal))then
   DEALLOCATE(modeldat%NSLIPglobal, STAT=ERROR)
endif

if (allocated(modeldat%NSLIP))then
   DEALLOCATE(modeldat%NSLIP, STAT=ERROR)
endif

if (allocated(modeldat%IDSLP))then
   DEALLOCATE(modeldat%IDSLP,STAT=ERROR)
endif

if (allocated(modeldat%IDSLE))then
   DEALLOCATE(modeldat%IDSLE,STAT=ERROR)
endif

if (allocated(modeldat%DIFORC))then
   DEALLOCATE(modeldat%DIFORC, STAT=ERROR)
endif

if (allocated(modeldat%IWINX))then
   DEALLOCATE(modeldat%IWINX,STAT=ERROR)
endif

if (allocated(modeldat%WINX))then
   DEALLOCATE(modeldat%WINX, STAT=ERROR)
endif

if (allocated(modeldat%DX))then
   DEALLOCATE(modeldat%DX,STAT=ERROR)
endif

if (allocated(modeldat%IDX))then
   DEALLOCATE(modeldat%IDX, STAT=ERROR)
endif

if (allocated(modeldat%IDXglobal))then
   DEALLOCATE(modeldat%IDXglobal,STAT=ERROR)
endif

if (allocated(modeldat%LMX))then
   DEALLOCATE(modeldat%LMX, STAT = ERROR)
endif

if (allocated(modeldat%DXE))then
   DEALLOCATE(modeldat%DXE, STAT = ERROR)
endif

if (allocated(modeldat%DELX))then
   DEALLOCATE(modeldat%DELX, STAT = ERROR)
endif

if (allocated(modeldat%DELDXE))then
   DEALLOCATE(modeldat%DELDXE, STAT = ERROR)
endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
if (allocated(modeldatactx%slipElemLookup)) then
    deallocate(modeldatactx%slipElemLookup)
endif
#endif

end subroutine deallocateslipperynodes
!********************************************************************

!**** fault parallel elements ***************************************
!********************************************************************
subroutine allocatefaultparallel(modeldat,nslskew)
use debugmodule, only: allocateError
use meshdatamodule, only: meshdatactx
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat

integer :: nslskew, ERROR
!--------------------------------------------------------------------
!ALLOCATE(modeldat%NSELSD(2,NSLSKEW), STAT=ERROR)
if (allocated(modeldat%NSELSD))then
    DEALLOCATE(modeldat%NSELSD,STAT=ERROR)
endif
ALLOCATE(modeldat%NSELSD(2,meshdatactx%nvglobal), STAT=ERROR)

call allocateError("modeldat%NSELSD", ERROR)

end subroutine allocatefaultparallel
!********************************************************************
   

!********************************************************************
subroutine deallocatefaultparallel(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%NSELSD))then
   DEALLOCATE(modeldat%NSELSD,STAT=ERROR)
endif

end subroutine deallocatefaultparallel
!--------------------------------------------------------------------

!**** thermal anomaly arrays ****************************************
!--------------------------------------------------------------------
subroutine allocatethermalanomaly(modeldat,ntanom)

use DebugModule, only: AllocateError
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ntanom, ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%ITANOM))then
    DEALLOCATE(modeldat%ITANOM,STAT=ERROR)
endif
ALLOCATE(modeldat%ITANOM(3,NTANOM),STAT=ERROR)

if (allocated(modeldat%TANOM))then
    DEALLOCATE(modeldat%TANOM,STAT=ERROR)
endif
ALLOCATE(modeldat%TANOM(NTANOM),STAT=ERROR)

call allocateError("Thermal anomaly", ERROR) 

modeldat%ITANOM = 0
modeldat%TANOM = 0d0

!--------------------------------------------------------------------
end subroutine allocatethermalanomaly
!********************************************************************
   

!********************************************************************
subroutine deallocatethermalanomaly(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%ITANOM))then
   DEALLOCATE(modeldat%ITANOM,STAT=ERROR)
endif

if (allocated(modeldat%TANOM))then
   DEALLOCATE(modeldat%TANOM,STAT=ERROR)
endif

end subroutine deallocatethermalanomaly
!--------------------------------------------------------------------


!**** thermal flux boundary conditions ******************************
subroutine allocatethermalfluxbcs(modeldat,nflx)

use DebugModule, only: AllocateError
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: nflx,ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%IFLX))then
    DEALLOCATE(modeldat%IFLX,STAT=ERROR)
endif
ALLOCATE(modeldat%IFLX(NFLX),STAT=ERROR)

if (allocated(modeldat%IFLS))then
    DEALLOCATE(modeldat%IFLS,STAT=ERROR)
endif
ALLOCATE(modeldat%IFLS(NFLX),STAT=ERROR)

if (allocated(modeldat%BFLX))then
    DEALLOCATE(modeldat%BFLX,STAT=ERROR)
endif
ALLOCATE(modeldat%BFLX(NFLX),STAT=ERROR)

call allocateError("Thermal flux", ERROR)

modeldat%IFLX = 0
modeldat%IFLS = 0
modeldat%BFLX = 0d0

end subroutine allocatethermalfluxbcs
!********************************************************************
   

!********************************************************************
subroutine deallocatethermalfluxbcs(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%IFLX))then
   DEALLOCATE(modeldat%IFLX,STAT=ERROR)
endif

if (allocated(modeldat%IFLS))then
   DEALLOCATE(modeldat%IFLS,STAT=ERROR)
endif

if (allocated(modeldat%BFLX))then
   DEALLOCATE(modeldat%BFLX,STAT=ERROR)
endif

end subroutine deallocatethermalfluxbcs
!********************************************************************


!**** thermal winkler fluxes ****************************************
!********************************************************************
subroutine allocatethermalwinkler(modeldat,ntwink)

use DebugModule, only: AllocateError
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ntwink,ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%ITWINK))then
    DEALLOCATE(modeldat%ITWINK,STAT=ERROR)
endif
ALLOCATE(modeldat%ITWINK(3,NTWINK), STAT=ERROR)

if (allocated(modeldat%TWINK))then
    DEALLOCATE(modeldat%TWINK,STAT=ERROR)
endif
ALLOCATE(modeldat%TWINK(NTWINK), STAT=ERROR)

call allocateError("Thermal Winkler", ERROR)

modeldat%ITWINK = 0
modeldat%TWINK = 0

end subroutine allocatethermalwinkler
!********************************************************************
   

!********************************************************************
subroutine deallocatethermalwinkler(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
integer :: ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%ITWINK))then
   DEALLOCATE(modeldat%ITWINK, STAT=ERROR)
endif

if (allocated(modeldat%TWINK))then
   DEALLOCATE(modeldat%TWINK, STAT=ERROR)
endif

end subroutine deallocatethermalwinkler
!********************************************************************


!**** thermalmodeldata **********************************************
!********************************************************************
subroutine allocatethermalmodeldata(modeldat, numnpGlobal, numel)

use DebugModule,    only: AllocateError
use meshdatamodule, only: meshdatactx
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
PetscInt :: numnpGlobal, numel, ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%IDT))then
    DEALLOCATE(modeldat%IDT,STAT=ERROR)
endif
!   ALLOCATE(modeldat%IDT(numnpGlobal),  STAT=ERROR)
ALLOCATE(modeldat%IDT(meshdatactx%nvlocal),  STAT=ERROR)

if (allocated(modeldat%TLST))then
    DEALLOCATE(modeldat%TLST,STAT=ERROR)
endif
ALLOCATE(modeldat%TLST(numnpGlobal),  STAT=ERROR)

if (allocated(modeldat%T))then
    DEALLOCATE(modeldat%T,STAT=ERROR)
endif
ALLOCATE(modeldat%T(numnpGlobal),  STAT=ERROR)

if (allocated(modeldat%T0))then
    DEALLOCATE(modeldat%T0,STAT=ERROR)
endif
ALLOCATE(modeldat%T0(numnpGlobal),  STAT=ERROR)


if (allocated(modeldat%ThermalFlux))then
    DEALLOCATE(modeldat%ThermalFlux,STAT=ERROR)
endif
ALLOCATE(modeldat%ThermalFlux(3,numel), STAT=ERROR)
call allocateError("Thermal modeldata", ERROR)

modeldat%IDT = 0
modeldat%T = 0d0
modeldat%ThermalFlux = 0d0

end subroutine allocatethermalmodeldata
!********************************************************************
   

!********************************************************************
subroutine deallocatethermalmodeldata(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
PetscInt :: ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%IDT))then
   DEALLOCATE(modeldat%IDT,  STAT=ERROR)
endif

if (allocated(modeldat%T))then
   DEALLOCATE(modeldat%T,  STAT=ERROR)
endif

if (allocated(modeldat%TLST))then
   DEALLOCATE(modeldat%TLST,  STAT=ERROR)
endif

if (allocated(modeldat%ThermalFlux))then
   DEALLOCATE(modeldat%ThermalFlux,  STAT=ERROR)
endif

end subroutine deallocatethermalmodeldata
!********************************************************************


!**** thermalbcdata *************************************************
!********************************************************************
subroutine allocatethermalbcdata(modeldat, numnp)

use DebugModule, only: AllocateError
!--------------------------------------------------------------------
!   use meshdatamodule, only: meshdatactx
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
PetscInt :: numnp, ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%ITBC))then
   DEALLOCATE(modeldat%ITBC,STAT=ERROR)
endif
ALLOCATE(modeldat%ITBC(1,numnp),  STAT=ERROR)


if (allocated(modeldat%ITBCGLOBAL))then
   DEALLOCATE(modeldat%ITBCGLOBAL,STAT=ERROR)
endif
ALLOCATE(modeldat%ITBCGLOBAL(numnp),  STAT=ERROR)

if (allocated(modeldat%TBC))then
   DEALLOCATE(modeldat%TBC,STAT=ERROR)
endif
ALLOCATE(modeldat%TBC(numnp),  STAT=ERROR)

call allocateError("Thermal boundary conditions", ERROR)

modeldat%ITBC = 0
modeldat%ITBCGLOBAL = 0
modeldat%TBC = 0d0

end subroutine allocatethermalbcdata
!********************************************************************
   

!********************************************************************
subroutine deallocatethermalbcdata(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
PetscInt :: ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%ITBC))then
   DEALLOCATE(modeldat%ITBC, STAT=ERROR)
endif

if (allocated(modeldat%ITBCGLOBAL))then
   DEALLOCATE(modeldat%ITBCGLOBAL, STAT=ERROR)
endif

if (allocated(modeldat%TBC))then
   DEALLOCATE(modeldat%TBC, STAT=ERROR)
endif

end subroutine deallocatethermalbcdata
!********************************************************************


!**** mechbcdata ****************************************************
!********************************************************************
subroutine allocatemechbcdata(modeldat, numnp, ndof, nlink)
use debugmodule,    only: allocateError
use modelctx,       only: getrank
!   use meshdatamodule, only: NvAllRelevant
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
PetscInt :: numnp, ndof, nlink, ERROR
!--------------------------------------------------------------------

if (allocated(modeldat%IBONDGLOBAL))then
    DEALLOCATE(modeldat%IBONDGLOBAL,STAT=ERROR)
endif
ALLOCATE(modeldat%IBONDGLOBAL(ndof, numnp),  STAT=ERROR)
call allocateError("IBONDGLOBAL", ERROR)

if (allocated(modeldat%IBOND))then
    DEALLOCATE(modeldat%IBOND,STAT=ERROR)
endif
ALLOCATE(modeldat%IBOND(ndof, numnp),  STAT=ERROR)
call allocateError("IBOND", ERROR)

if (allocated(modeldat%BOND))then
    DEALLOCATE(modeldat%BOND,STAT=ERROR)
endif
ALLOCATE(modeldat%BOND(ndof, numnp),  STAT=ERROR)
call allocateError("BOND", ERROR)

if (allocated(modeldat%LINK))then
    DEALLOCATE(modeldat%LINK,STAT=ERROR)
endif
ALLOCATE(modeldat%LINK(3, nlink),  STAT=ERROR)
call allocateError("LINK", ERROR)

if (allocated(modeldat%LINKID))then
    DEALLOCATE(modeldat%LINKID,STAT=ERROR)
endif
ALLOCATE(modeldat%LINKID(2-ndof, nlink),  STAT=ERROR)
call allocateError("LINKID", ERROR)

modeldat%IBONDGLOBAL = 0
modeldat%IBOND = 0
modeldat%BOND = 0d0
modeldat%LINK = 0
modeldat%LINKID = 0

end subroutine allocatemechbcdata
!********************************************************************
   

!********************************************************************
subroutine deallocatemechbcdata(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
PetscInt :: ERROR
!--------------------------------------------------------------------
!   DEALLOCATE(modeldat%IBONDGLOBAL, STAT=ERROR)
if (allocated(modeldat%IBOND))then
   DEALLOCATE(modeldat%IBOND, STAT=ERROR)
endif

if (allocated(modeldat%BOND))then
   DEALLOCATE(modeldat%BOND, STAT=ERROR)
endif

if (allocated(modeldat%LINK))then
   DEALLOCATE(modeldat%LINK, STAT=ERROR)
endif

if (allocated(modeldat%LINKID))then
   DEALLOCATE(modeldat%LINKID, STAT=ERROR)
endif

end subroutine deallocatemechbcdata
!********************************************************************


!**** mechwinklerbcdata *********************************************
!********************************************************************
subroutine allocatemechwinklerbcdata(modeldat, numnp, ndof)

use DebugModule, only: AllocateError
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
PetscInt :: numnp, ndof, ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%IWINK))then
   DEALLOCATE(modeldat%IWINK,STAT=ERROR)
endif
ALLOCATE(modeldat%IWINK(ndof, numnp),  STAT=ERROR)

if (allocated(modeldat%WINK))then
   DEALLOCATE(modeldat%WINK,STAT=ERROR)
endif
ALLOCATE(modeldat%WINK(ndof, numnp),  STAT=ERROR)

call allocateError("Winkler BC", ERROR)

modeldat%IWINK = 0
modeldat%WINK = 0d0

end subroutine allocatemechwinklerbcdata
!********************************************************************
   

!********************************************************************
subroutine deallocatemechwinklerbcdata(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
PetscInt :: ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%IWINK))then
   DEALLOCATE(modeldat%IWINK, STAT=ERROR)
endif

if (allocated(modeldat%WINK))then
   DEALLOCATE(modeldat%WINK, STAT=ERROR)
endif

end subroutine deallocatemechwinklerbcdata
!********************************************************************


!**** mechmaxwellbcdata *********************************************
!********************************************************************
subroutine allocatemechmaxwellbcdata(modeldat, numnp, ndof)

use DebugModule, only: AllocateError
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
PetscInt :: numnp, ndof, ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%FMAXW))then
   DEALLOCATE(modeldat%FMAXW,STAT=ERROR)
endif
ALLOCATE(modeldat%FMAXW(ndof, numnp),  STAT=ERROR)

if (allocated(modeldat%FTOT))then
   DEALLOCATE(modeldat%FTOT,STAT=ERROR)
endif
ALLOCATE(modeldat%FTOT(ndof, numnp),  STAT=ERROR)

call allocateError("Maxwell boundary conditions", ERROR)

modeldat%FMAXW = 0d0
modeldat%FTOT = 0d0

end subroutine allocatemechmaxwellbcdata
!********************************************************************
   

!********************************************************************
subroutine deallocatemechmaxwellbcdata(modeldat)
!--------------------------------------------------------------------
implicit none
!--------------------------------------------------------------------
type (modeldata) modeldat
PetscInt :: ERROR
!--------------------------------------------------------------------
if (allocated(modeldat%FMAXW))then
   DEALLOCATE(modeldat%FMAXW,  STAT=ERROR)
endif

if (allocated(modeldat%FTOT))then
   DEALLOCATE(modeldat%FTOT,  STAT=ERROR)
endif
end subroutine deallocatemechmaxwellbcdata
!********************************************************************
   

!********************************************************************

!-----routine to
integer function INLINK (NUMDOF,NODE)

USE MODELDEFINITION

implicit none

integer          :: NUMDOF, NODE, i, n, idof
!
! Checks if NDOF,NODE occurs in LINK
! LINK(1) = IDOF
! LINK(2) = number of node with displacement equal to that of
! LINK(3)
!

! returns ID of linked node if the node is linked.
! if not linked, returns -1

INLINK = -1
i = 0

100   i = i + 1

if (i.gt.NLINK) then
    return
endif

idof = modeldatactx%LINK(1,i)

if (idof.ne.NUMDOF) then
    goto 100
endif

n = modeldatactx%LINK(2,i)

if (n.ne.NODE) then
    goto 100
endif

INLINK = i
return

end function INLINK


!********************************************************************
logical function HasSlipperyNodes()
!--------------------------------------------------------------------
use modeldefinition, only: NUMSLP, &  ! number of slippery line entries in TECIN.DAT
                           NUMSN      ! number of physical slippery nodes in the partition

implicit none
! bit weird to have this as a separate function,
! but the advantages is that the selection criterium can be changed
! throughout the code by changing only this.
!--------------------------------------------------------------------
!   integer :: nslippery, nsplit
!--------------------------------------------------------------------
HasSlipperyNodes = .false.

if ((NUMSLP.gt.0) .or. &
    (NUMSN .gt.0)) then
    HasSlipperyNodes = .true.
endif

end function HasSlipperyNodes
!********************************************************************

END MODULE MODELDATAMODULE
