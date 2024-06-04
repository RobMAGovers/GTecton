MODULE AOMODULE
!***************************************************************
! The AO Module administrates the numbering of the grid:
! the mesh points and the elements. Each have global numbers
! and local numbers  (AO = Application Ordering)
!***************************************************************

use petscksp

   implicit none

#ifdef SPARSE

#define PETSC_AVOID_DECLARATIONS
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscao.h"
#include "petsc/finclude/petscis.h"
#undef PETSC_AVOID_DECLARATIONS

   ! application ordering 
   ! np: nodal points
   ! el: elements
   ! eq: equations
   ! eqt: equations thermal
   AO :: aonp, aoel, aoeq, aoeqt

   IS :: isglobal, islocal
   IS :: isglobald, islocald
   IS :: isglobalt, islocalt
   PetscInt, ALLOCATABLE :: vertices(:)
   PetscInt, ALLOCATABLE :: verticesmask(:)
   PetscInt              :: nvertices
   PetscInt, ALLOCATABLE :: equations(:)
   PetscInt, ALLOCATABLE :: equationsmask(:,:)
   PetscInt, ALLOCATABLE :: equationsmaskx(:,:)   ! e(X)tra; for slippery nodes
   PetscInt              :: nequations
   ISLocalToGlobalMapping :: isl2g

   ! thermal data
   PetscInt, ALLOCATABLE :: equationst(:)
   PetscInt, ALLOCATABLE :: equationsmaskt(:)
   PetscInt, ALLOCATABLE :: equationsmaskxt(:)    ! e(X)tra; for slippery nodes
   PetscInt              :: nequationst
   ISLocalToGlobalMapping :: isl2gt

   ! some datavectors for scattering purposes
   Vec :: ndofVec
   Vec :: nsdVec
   Vec :: oneVec
   Vec :: twoVec
   Vec :: ndofVecLoc
   Vec :: nsdVecLoc
   Vec :: oneVecLoc
   Vec :: twoVecLoc
   VecScatter :: ndofScatter
   VecScatter :: nsdScatter
   VecScatter :: oneScatter
   VecScatter :: twoScatter

   ! global variables for the variable parallel numbering routines
   integer :: rstart ! indicates start of memory where the various vectors are stored.
   PetscInt, ALLOCATABLE :: npordering(:)
   PetscInt, ALLOCATABLE :: elordering(:)


   logical,save :: ndofVecSet
   logical,save :: nsdVecSet
   logical,save :: oneVecSet
   logical,save :: twoVecSet



   PUBLIC  :: createEquationNumbering
   public  :: createMeshNumbering
!   public  :: createEquationNumbering


!      PUBLIC  :: createscattering
   PUBLIC  :: writeinputdata   
   PUBLIC  :: writemodeldata
   PUBLIC  :: scatterdata
   PUBLIC  :: scatterNDOFdata
   PUBLIC  :: scatterONEdata
   PUBLIC  :: scatterTWOdata
   PUBLIC  :: createlocaltoglobalmapping
   PUBLIC  :: getallvectorentries
   PUBLIC  :: writemeshdata
   PUBLIC  :: cleanAOdata

   PUBLIC  :: FillLocInd
   PUBLIC  :: FillLocElt
   PUBLIC  :: FillGlo2LocVertex
   PUBLIC  :: FillGlo2LocElement
   PUBLIC  :: FillFullGlo2LocElement


   PUBLIC  :: CreateVerticesMask

   PRIVATE :: CreateEquationsMask
   PRIVATE :: writeAOdatasummary
   PRIVATE :: scatterONEdataInit
   PRIVATE :: scatterTWOdataInit
   PRIVATE :: scatterNDOFdataInit
   PRIVATE :: scatterNSDdataInit
   PRIVATE :: scatterNSDdata
   PRIVATE :: writeAOmapping
   PRIVATE :: writeAOmeshdata
   PRIVATE :: renumbernodalpoints
   PRIVATE :: renumberequations

   contains


!***************************************************************
subroutine FillLocInd()
!***************************************************************
! GloInd has been read from the TECIN input data
! Local index has to be constructed
!***************************************************************
use  meshdatamodule, only:  meshdatactx
use modelctx
use debugmodule,     only: parallellog, &
                           iecho, &
                           debug, &
                           allocateError

!***************************************************************
implicit none
!***************************************************************
!#include "petsc/finclude/petscsysdef.h"
!***************************************************************
integer          :: i
integer          :: ierr, ERROR
!***************************************************************


!call MPI_BARRIER(0,ierr) ! 0 is mpi_comm_world

rstart = 0 ! from the module

!write(*,*) "rank", getrank(), "gloInd: ", meshdatactx%gloInd

! determine a continuing count over the processors
! example, let nvlocal over 4 procs be:  10   12   9   15
! rstart will be   10    22    31    46
call MPI_Scan(meshdatactx%Nvlocal, rstart, 1, &
MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)


rstart = rstart - meshdatactx%Nvlocal

! npordering directly from the module
ALLOCATE(npordering(meshdatactx%Nvlocal), STAT=ERROR)
call allocateError("npordering",ERROR)


! down from 1-based (Fortran) to 0-based (Petsc)
do i=1,meshdatactx%Nvlocal
    npordering(i) = rstart + i-1
enddo

do i=1,meshdatactx%Nvlocal
    meshdatactx%gloInd(i) = &
    meshdatactx%gloInd(i) - 1
enddo

! generate the aonp, which is essentially 
! a translation table between global and local
! indices.

call AOCreateBasic(MPI_COMM_WORLD, &
                   meshdatactx%Nvlocal, &
                   meshdatactx%gloInd, &
                   npordering, &
                   aonp, &
                   ierr)

CHKERRQ(ierr)

if (ierr.ne.0) then
    write(*,*) "AOmodule could not AOcreatebasic in subroutine FillLocInd, error: ", ierr
endif

DEALLOCATE(npordering, STAT=ERROR)

do i=1,meshdatactx%Nvlocal
    meshdatactx%locInd(i) = &
    meshdatactx%gloInd(i)
enddo

call AOApplicationToPetsc(aonp, &
                          meshdatactx%Nvlocal, &
                          meshdatactx%locInd, ierr)

do i=1,meshdatactx%Nvlocal
    meshdatactx%gloInd(i) = meshdatactx%gloInd(i) + 1
enddo

if (debug) then
!    write(*,*) 'rank ', getrank(), 'has locInd: ', meshdatactx%locInd 
!    write(*,*) 'rank ', getrank(), 'has gloInd: ', meshdatactx%gloInd
endif

!call MPI_BARRIER(0,ierr) ! 0 is mpi_comm_world


if (iecho.eq.8) then
    call parallellog("fillocind", "loc index completed")
endif

return

end subroutine
!***************************************************************


!***************************************************************
subroutine FillLocElt()
!
! GloInd has been read from the TECIN input data. Local index has to be constructed
!
use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use modelctx,        only: getrank, getsize
use debugmodule,     only: parallellog, &
                           debug, &
                           allocateError, &
                           delay

implicit none
#include "petsc/finclude/petsc.h" 
!-locl
integer          :: i, j
integer          :: ierr, ERROR

call MPI_Scan(meshdatactx%Nelocal, rstart, 1, &
  MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

rstart = rstart - meshdatactx%Nelocal

! elordering directly from the module
ALLOCATE(elordering(meshdatactx%Nelocal), STAT=ERROR)
call allocateError("elordering", ERROR)

do i=1,meshdatactx%Nelocal
    elordering(i) = rstart + i-1
enddo

do i=1,meshdatactx%Nelocal
    meshdatactx%gloElt(i) = meshdatactx%gloElt(i) - 1
enddo

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!do i=0,getsize()-1
!    call MPI_BARRIER(0, ierr)

!    call delay
!    if (getrank() .eq. i) then
!        write(*,*) "gloElt", meshdatactx%gloElt
!    endif
!    call delay
!enddo


!write(*,*) "meshdatactx%Nelocal",meshdatactx%Nelocal


!do i=1,meshdatactx%Nelocal!
!    write(*,*) "meshdatactx%gloElt",i,meshdatactx%gloElt(i)
!enddo

!write(*,*) "rank", getrank(),"nelocal: ", meshdatactx%Nelocal
!write(*,*) "rank", getrank(),"glo Elt", meshdatactx%gloElt


!write(*,*) "nelocal: ", meshdatactx%Nelocal
!write(*,*) "gloElt: ", meshdatactx%gloElt(1:10)
!write(*,*) "gloElt: ", &
!meshdatactx%gloElt(meshdatactx%Nelocal-10:meshdatactx%Nelocal)


if (meshdatactx%gloElt(meshdatactx%Nelocal).eq.-1 ) then
    write(*,*) "Error. gloElt array has not been filled completely"
    write(*,*) "Possible cause: element input is partitioned into more"
    write(*,*) "partitions than are used for the actual run."
    stop "Aborting run"
endif

call AOCreateBasic(MPI_COMM_WORLD, meshdatactx%Nelocal, &
                   meshdatactx%gloElt, elordering, aoel, ierr)

if (ierr.ne.0) then
    write(*,*) "AOmodule could not AOcreatebasic in subroutine FillLocElt, error: ", ierr
endif

CHKERRQ(ierr)


DEALLOCATE(elordering, STAT=ERROR)

do i=1,meshdatactx%Nelocal
    meshdatactx%locElt(i) = meshdatactx%gloElt(i)
enddo

call MPI_BARRIER(MPI_COMM_WORLD, ierr)


call AOApplicationToPetsc(aoel, meshdatactx%Nelocal, &
                          meshdatactx%locElt, ierr)

call AODestroy(aoel, ierr)

do i=1,meshdatactx%Nelocal
    meshdatactx%gloElt(i) = meshdatactx%gloElt(i) + 1
enddo

return

end subroutine FillLocElt
!***************************************************************



!***************************************************************
   subroutine FillGlo2LocVertex
!***************************************************************

   ! Now that GloInd and LocInd are known, we can build 
   ! a look up able to find the local number.

   ! meshdatactx%Glo2LocVertex(i) will hold the local
   ! index of the vertex 

!***************************************************************
   use meshdatamodule, only: meshdatactx
   use debugmodule,    only: iecho, &
                             parallellog, &
                             AllocateError
   use modelctx,       only: getrank, &
                             getsize
!***************************************************************
   implicit none
!***************************************************************

#include "mpif.h"

   integer :: ERROR
   integer :: i

   integer, allocatable :: offsets(:)
   integer, allocatable :: nvlocalOfAll(:)
   integer, allocatable :: Loc2Glo(:)

!***************************************************************

ERROR = 0

if (iecho.eq.8) then
    call parallelLog("FillGlo2LocVertex","allocating Loc2Glo for vertices")
endif

!write(*,*) ''  
! somehow this prevents the code from crashing on the following allocate.
! your guess is as good as mine... * confuzzled *

   allocate(Loc2Glo(meshdatactx%nvglobal), STAT=ERROR)
   call AllocateError("Loc2Glo", ERROR)

if (iecho.eq.8) then
    call parallelLog("FillGlo2LocVertex","allocated Loc2Glo for vertices")
endif



   ALLOCATE(meshdatactx%Glo2LocVertex(meshdatactx%nvglobal), STAT=ERROR)
   call AllocateError("meshdatactx%Glo2LocVertex", ERROR)

   ! each processor has its global indices known, and the local
   ! indices have been computed from this.
   ! use these to construct a translation table.

   ! we have to gather this from all the processors,
   ! but not all blocks have the same size, because not all
   ! processors have the same number of vertices.
   ! So we have to compute the offsets first, based on
   ! meshdatactx%nvlocal
   
   ! 1: make sure that all processors know the number of vertices 
   ! of each other processor
   allocate(nvlocalOfAll(getsize()))
   call mpi_allgather(meshdatactx%nvlocal, 1, MPI_int, nvlocalOfAll, 1, mpi_int, MPI_COMM_WORLD, ERROR)



   ! 2: compute the offset from that
   ! which is a cumulative sum of all sizes.
   allocate(offsets(getsize()))
   offsets(1) = 0
   do i=2,getsize()
       offsets(i) = offsets(i-1) + nvlocalOfAll(i-1)
   enddo


! 3: with these offsets, gather all global indices to the 
! translation array from local to global
do i=1,getsize()
    call mpi_gatherv(meshdatactx%gloInd, & ! sending array
                     meshdatactx%nvlocal,& ! number of elements sent
                     MPI_int,  &           ! type sent
                     Loc2Glo,&             ! destination array
                     nvlocalOfAll,&        ! number of elements received, from each proc
                     offsets,&             ! offset from start of receiving array
                     mpi_int,&             ! type sent
                     i-1, &                ! recipient processor
                     MPI_COMM_WORLD, &     ! communicator
                     error)                ! error

    ! barrier to prevent wayward traffic from becoming entangled
    call mpi_barrier(MPI_COMM_WORLD, error)
   enddo

   deallocate(offsets)
   deallocate(nvlocalOfAll)

   ! from local to global, construct the global to local array
   do i=1,meshdatactx%nvglobal
       meshdatactx%Glo2LocVertex(Loc2Glo(i)) = i
   enddo


!    write(*,*) "rank", getrank(), "glo2loc lookup: ", meshdatactx%Glo2LocVertex

   deallocate(Loc2Glo)

!   if (iecho.eq.8) then
!       call parallellog("fillglo2locvertex", "lookup table vertex done")
!   endif

call MPI_BARRIER(MPI_COMM_WORLD,error) ! 0 is mpi_comm_world


   return

   end subroutine FillGlo2LocVertex
!***************************************************************


! should be more something like this,
! for some reason that does not work.

!***************************************************************
   subroutine FillFullGlo2LocElement
!***************************************************************
   ! Now that GloElt and LocElt are known, we can build 
   ! a look up table to find the local number
!***************************************************************
   use meshdatamodule, only: meshdatactx
   use debugmodule,    only: AllocateError, &
                             iecho, &
                             parallelLog, &
                             debug
   use modelctx,       only: getrank, &
                             getsize
!***************************************************************
   implicit none
#include "mpif.h" 
!***************************************************************
   integer :: ERROR
   integer :: i
   integer, allocatable :: offsets(:)
   integer, allocatable :: nelocalOfAll(:)
   integer, allocatable :: Loc2Glo(:)

!***************************************************************

! The global-to-local lookup table for the vertices is different
! from the lookup for elements, because it has a different purpose.

! The vertex look-up is needed for reading the nodal boundary conditions
! outside the partition boundaries

! The element look-up had the goal of finding positions in IEN
! for a global element number for the boundary conditions.

! Because elements do not cross domain boundaries, we only need to set
! the lookup values for the locally owned elements.

! Except for the fault parallel elements, which is treated per 
! slippery point, and for every slippery point considers the
! elements boundaring on that points. Some or even all of those elements
! can be in another partition (although all is unlikely)

! Because the global values still range over the total number of the
! elements, the lookup table needs to have the size of the global 
! number of elements.


call MPI_BARRIER(MPI_COMM_WORLD,ERROR)


! this lookup table is built in a similar fashion 
! as the vertex table, right above. For detailed comment on the method,
! see there.

ALLOCATE(Loc2Glo(meshdatactx%neglobal), STAT=ERROR)
call AllocateError("Loc2Glo", ERROR)

ALLOCATE(meshdatactx%FullGlo2LocElement(meshdatactx%neglobal), STAT=ERROR)
call AllocateError("meshdatactx%Glo2LocElement", ERROR)


! 1: make sure that all processors know the number of vertices 
! of each other processor

allocate(nelocalOfAll(getsize()))
call mpi_allgather(meshdatactx%nelocal, 1, MPI_int, nelocalOfAll, 1, mpi_int, MPI_COMM_WORLD, ERROR)


! 2: compute the offset from that
! which is a cumulative sum of all sizes.

allocate(offsets(getsize()))
offsets(1) = 0
do i=2,getsize()
    offsets(i) = offsets(i-1) + nelocalOfAll(i-1)
enddo

! 3: with these offsets, gather all global indices to the 
! translation array from local to global

do i=1,getsize()
    call mpi_gatherv(meshdatactx%gloElt, & ! sending array
                     meshdatactx%nelocal,& ! number of elements sent
                     MPI_int,  &           ! type sent
                     Loc2Glo,&             ! destination array
                     nelocalOfAll,&        ! number of elements received, from each proc
                     offsets,&             ! offset from start of receiving array
                     mpi_int,&             ! type sent
                     i-1, &                ! recipient processor
                     MPI_COMM_WORLD, &     ! communicator
                     error)                ! error
    ! barrier to prevent wayward traffic from becoming entangled
    call mpi_barrier(MPI_COMM_WORLD, error)
enddo

deallocate(offsets)
deallocate(nelocalOfAll)

! from local to global, construct the global to local array
do i=1,meshdatactx%neglobal
    meshdatactx%FullGlo2LocElement(Loc2Glo(i)) = i
enddo

deallocate(Loc2Glo)

!write(*,*) "rank", getrank(), "bravely made elem lookup table: ", meshdatactx%FullGlo2LocElement


if (iecho.eq.8) then
    call parallellog("fillglo2locElement", "lookup table elements done")
endif

call MPI_BARRIER(MPI_COMM_WORLD,error)

return

end subroutine FillFullGlo2LocElement
!***************************************************************


!***************************************************************
subroutine FillGlo2LocElement
!***************************************************************
! Now that GloElt and LocElt are known, we can build 
! a look up able to find the local number
!***************************************************************
use meshdatamodule, only: meshdatactx
use modelctx, only: getrank
use debugmodule, only: allocateError
!***************************************************************
implicit none
!***************************************************************
integer :: ERROR
integer :: i
!***************************************************************

ALLOCATE(meshdatactx%Glo2LocElement(meshdatactx%neglobal), STAT=ERROR)
call allocateError("glo2loElement", ERROR)

meshdatactx%Glo2LocElement = -1

! write(*,*) "gloElt: ", meshdatactx%GloElt


do i=1,meshdatactx%nelocal
    meshdatactx%Glo2LocElement(meshdatactx%GloElt(i)) = i
enddo

return

end subroutine
!***************************************************************





!***************************************************************
subroutine countMechanicalEquations &
         (tmpeq, toteqs, totideqs, totidxeqs,ndof)

! count the mechanical equations for the normal nodes (modeldata ID)
! and those of the slippery nodes (modeldata IDX) and combine the
! two into the tmpeq, which will be used in AOmodule to distribute the
! equations over the processors.

use modeldefinition, only: NEQlocal
use modeldatamodule, only: modeldatactx, &
                           HasSlipperyNodes 
use meshdatamodule,  only: meshdatactx
use modelctx,        only: getrank
use debugmodule,     only: debug, &
                           allocateError

implicit none

integer, allocatable :: tmpeq(:)
integer :: next, toteqs, totideqs, totidxeqs,ndof
integer :: i,k
integer :: ERROR

ALLOCATE(tmpeq(NEQlocal), STAT=ERROR)
call allocateError("tmpeq", ERROR)

do i=1,NEQlocal
    tmpeq(i) = 0
enddo

next = 0 ! internal use
toteqs = 0 ! not used, remains 0
totideqs = 0 ! nNonSlipperyEquations
totidxeqs = 0 ! nSLipperyEquations

if (debug) then
!    write(*,*) "rank", getrank(),"has modeldat ID", modeldatactx%ID
!    write(*,*) "rank", getrank(),"has modeldat IDX", modeldatactx%IDX
endif


do i=1,meshdatactx%Nvlocal
    do k=1,NDOF

        if (HasSlipperyNodes()) then
            if (modeldatactx%IDX(k,i).ne.0) then
                totidxeqs = totidxeqs + 1
                next = next + 1
                tmpeq(next) = modeldatactx%IDX(k,i) - 1
            endif
        endif

        if (modeldatactx%ID(k,i).ne.0) then
            next = next + 1
            tmpeq(next) = modeldatactx%ID(k,i) - 1
        endif

    enddo
enddo

if (debug) then
!    write(*,*) "rank", getrank(),"used to construct tmpeq", tmpeq
endif

toteqs = next

end subroutine

!***************************************************************

!***************************************************************
subroutine countThermalEquations(tmpeqt, toteqst, totideqst, totidxeqst)

use modeldefinition, only: NTEQlocal
use modeldatamodule, only: modeldatactx
use modelctx,        only: hasdiffthermal, getrank
use meshdatamodule,  only: meshdatactx
use debugmodule,     only: allocateError

implicit none

integer, allocatable :: tmpeqt(:)
integer :: next, toteqst, totideqst, totidxeqst
integer :: i
integer :: ERROR

ALLOCATE(tmpeqt(NTEQlocal), STAT=ERROR)
call allocateError("tmpqeqt", ERROR)

tmpeqt = 0

next = 0
toteqst = 0
totideqst = 0
totidxeqst = 0

!write(*,*) "thread", getrank(), "CTE has diff thermal: ", hasdiffthermal
!write(*,*) "thread", getrank(), "CTE has modeldatactx%IDTX: ", modeldatactx%IDTX
!write(*,*) "thread", getrank(), "CTE has modeldatactx%IDT: ", modeldatactx%IDT

do i=1,meshdatactx%Nvlocal
    if (hasdiffthermal) then
        if (modeldatactx%IDTX(i).ne.0) then
            totidxeqst = totidxeqst + 1
            next = next + 1
            tmpeqt(next) = modeldatactx%IDTX(i) - 1
        endif
    endif

    if (modeldatactx%IDT(i).ne.0) then
        totideqst = totideqst + 1
        next = next + 1
        tmpeqt(next) = modeldatactx%IDT(i) - 1
    endif
enddo

toteqst = next

end subroutine


!***************************************************************

subroutine buildLM(NEN, NDOF)

    use meshdatamodule,  only: meshdatactx
    use modeldatamodule, only: modeldatactx
    use debugmodule,     only: AllocateError
    use modelctx,        only: getrank

    implicit none
    
    integer :: ERROR
    integer :: i,j,k
    integer :: nn
    integer :: NEN, NDOF

    ALLOCATE(modeldatactx%lm(NDOF,NEN,meshdatactx%Nelocal),STAT=ERROR)
    call AllocateError("Modeldata LM",ERROR)

    do i=1,meshdatactx%Nelocal
        do j=1,NEN
            nn = meshdatactx%ien(j,i)
            do k=1,NDOF
!                if (modeldatactx%IDglobal(k,nn).ne.0) then
!                    totelmeqs = totelmeqs + 1
!                endif
                modeldatactx%lm(k,j,i) = modeldatactx%IDglobal(k,nn)
            enddo
        enddo
    enddo 


    write(*,*) 'rank', getrank(), 'completed LM: ', modeldatactx%lm

end subroutine


!***************************************************************
subroutine CountEquations &
         (NEN, NDOF, totelmeqs, totelmeqsx, totelmeqst, totelmeqsxt)

! the LM arrays contain equation numbers per vertex per element
! So the total number of equations is not the total number of 
! euations.

! Note also that it counts over four quadrature points (NEN).
! So in 2D, the equation of the 3rd point, which is mentioned twice in IEN,
! is counted twice.

! However, the differential array LMX has a zero in the 
! 4th entry, even in 2D, so a slippery equation on the 3rd vertex is 
! only counted once for every elemet bordering on the node.

#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldefinition, only: nLocalElemsWithSLipperyNodes
use modeldatamodule, only: modeldatactx, &
                           HasSlipperyNodes, &
                           elementSlipperyID
#else
use modeldatamodule, only: modeldatactx, &
                           HasSlipperyNodes
#endif
use meshdatamodule,  only: meshdatactx
use modelctx,        only: hasthermal, &
                           hasdiffthermal, &
                           getrank
use debugmodule,     only: allocateError, &
                           debug

implicit none

integer :: i,j,k
integer :: nn
integer :: NEN, NDOF
integer :: totelmeqs, totelmeqsx, totelmeqst, totelmeqsxt
integer :: ERROR


integer :: globalElemID, slipperyElemID

totelmeqs = 0
totelmeqsx = 0
totelmeqst = 0
totelmeqsxt = 0




! hmmm... why is LM not yet allcated, but LMX is...
ALLOCATE(modeldatactx%lm(NDOF,NEN,meshdatactx%Nelocal),STAT=ERROR)

call allocateError("modeldat LM", ERROR)


!write(*,*) "rank", getrank(), "has LMX", modeldatactx%lmx

do i=1,meshdatactx%Nelocal
    do j=1,NEN
        nn = meshdatactx%ien(j,i)
        do k=1,NDOF
            if (modeldatactx%IDglobal(k,nn).ne.0) then
                totelmeqs = totelmeqs + 1
            endif
!           build LM array somewhere else
            modeldatactx%lm(k,j,i) = modeldatactx%IDglobal(k,nn)
        enddo
    enddo
enddo

!write(*,*) "rank", getrank(), "has regular eqs: ", totelmeqs
#ifdef EXPERIMENTAL_ReducedSlipAlloc
!write(*,*) "rank", getrank(), "has slip elems", nLocalElemsWithSLipperyNodes
#endif

if (HasSlipperyNodes()) then
    totelmeqsx = 0
    ! old way, walk through all local nodes of LMX.
    ! but LMX has been allocated globally...

#ifdef EXPERIMENTAL_ReducedSlipAlloc
    do i=1,meshdatactx%Nelocal

!        write(*,*) "rank",getrank(), "reads", i, "of", size(meshdatactx%gloElt,1)
        globalElemID = meshdatactx%gloElt(i)
!        write(*,*) "rank",getrank(), "glob elm ID", globalElemID
        slipperyElemID = elementSlipperyID(globalElemID,1)
!        write(*,*) "rank",getrank(), "grabs slip node", i, globalElemID, slipperyElemID
        
!        write(*,*) "size LMX 1", size(modeldatactx%lmx,1)
!        write(*,*) "size LMX 2", size(modeldatactx%lmx,2)
!        write(*,*) "size LMX 3", size(modeldatactx%lmx,3)

        if (slipperyElemID.gt.0) then
            do j=1,NEN
                do k=1,NDOF
                    ! lmx has been allocated in subroutine allocateslipperynodes
                    ! and filled in subroutine LOCALX in tecin.dat
                    ! but only carries data for the slippery nodes.
                    if (modeldatactx%lmx(k,j,slipperyElemID).ne.0) then
                        totelmeqsx = totelmeqsx + 1
                    endif
                enddo
            enddo
        else
            ! the element has no slippery nodes. Do nothing
        endif
    enddo



#else

    do i=1,meshdatactx%Nelocal
        do j=1,NEN
            do k=1,NDOF
                ! lmx has been allocated in subroutine allocateslipperynodes
                ! and filled in subroutine LOCALX in tecin.dat
                if (modeldatactx%lmx(k,j,i).ne.0) then
                    totelmeqsx = totelmeqsx + 1
                endif
            enddo
        enddo
    enddo

#endif

endif


!write(*,*) "rank", getrank(), "has slippery eqs: ", totelmeqsx



if (hasthermal) then
    totelmeqst = 0
    do i=1,meshdatactx%Nelocal
        do j=1,NEN
            ! lmt has been allocated in rdthrm
            if (modeldatactx%lmt(j,i).ne.0) then
                totelmeqst = totelmeqst + 1
            endif
        enddo
    enddo
endif

if (hasdiffthermal) then
    totelmeqsxt = 0
    do i=1,meshdatactx%Nelocal
        do j=1,NEN
            ! lmtx has been allocated in rdthrm
            if (modeldatactx%lmtx(j,i).ne.0) then
                totelmeqsxt = totelmeqsxt + 1
            endif
        enddo
    enddo
endif

end subroutine
!***************************************************************


!***************************************************************
subroutine buildAdjancyArrays &
         (NSD, NDOF, totadjeqs, totadjeqsx, totadjeqst, totadjeqsxt)

use modeldatamodule, only: modeldatactx, &
                           HasSlipperyNodes
use  meshdatamodule, only: meshdatactx
use modelctx,        only: hasthermal, &
                           hasdiffthermal, & 
                           getrank
use debugmodule,     only: allocateError, &
                           iecho, &
                           parallellog
use constants,       only: NNeighborsMax2D, &
                           NNeighborsMax3D

implicit none

integer :: i,j,k,nn
integer :: maxnbrs, nsd, ndof
integer :: ERROR
integer :: totadjeqs, totadjeqsx, totadjeqst, totadjeqsxt

totadjeqs = 0
totadjeqsx = 0
totadjeqst = 0
totadjeqsxt = 0

if (NSD.eq.2) then
    maxnbrs = NNeighborsMax2D
else
    maxnbrs = NNeighborsMax3D
endif

if (iecho.eq.8) then
    call ParallelLog("buildAdjancyArrays", "allocating arrays")
endif


ALLOCATE(modeldatactx%AdjID(NDOF, maxnbrs, meshdatactx%Nvlocal),STAT=ERROR)
call allocateError("AdjID", ERROR)

if (HasSlipperyNodes()) then
    ALLOCATE(modeldatactx%AdjIDX(NDOF, maxnbrs, meshdatactx%Nvlocal),STAT=ERROR)
    call allocateError("AdjIDX", ERROR)
endif

if (hasthermal) then
    ALLOCATE(modeldatactx%AdjIDT(maxnbrs,meshdatactx%Nvlocal),STAT=ERROR)
    call allocateError("AdjIDT", ERROR)

    if (hasdiffthermal) then
        ALLOCATE(modeldatactx%AdjIDTX(maxnbrs,meshdatactx%Nvlocal), STAT=ERROR)
        call allocateError("AdjIDTX", ERROR)
    endif
endif

if (iecho.eq.8) then
    call ParallelLog("buildAdjancyArrays", "allocated arrays")
endif


!     fil AdjID array (equation number per dof per vertex neighbor per
!     vertex)
modeldatactx%AdjID = 0
do i=1,meshdatactx%Nvlocal
  do j=1,meshdatactx%itot(i)
      nn = meshdatactx%AdjM(j,i)
      do k=1,NDOF
          if (modeldatactx%IDglobal(k,nn).ne.0) then
              totadjeqs = totadjeqs + 1
          endif
          modeldatactx%AdjID(k,j,i) = modeldatactx%IDglobal(k,nn)
      enddo
  enddo
enddo

if (iecho.eq.8) then
    call ParallelLog("buildAdjancyArrays", "finished constructing AdjID")
endif


!     like above for the splippery nodes
if (HasSlipperynodes()) then
  modeldatactx%AdjIDX = 0
  do i=1,meshdatactx%Nvlocal
      do j=1,meshdatactx%itot(i) ! loop over neighbors
          nn = meshdatactx%AdjM(j,i)
          do k=1,NDOF
              if (modeldatactx%IDXglobal(k,nn).ne.0) then
                  totadjeqsx = totadjeqsx + 1
              endif
              modeldatactx%AdjIDX(k,j,i) = modeldatactx%IDXglobal(k,nn)
          enddo
      enddo

  enddo
endif

if (iecho.eq.8) then
    call ParallelLog("buildAdjancyArrays", "finished constructing AdjIDX")
endif

if (hasthermal) then

    if (iecho.eq.8) then
        call ParallelLog("buildAdjancyArrays", "allocating therma adjancency arrays")
    endif


    modeldatactx%AdjIDT = 0
    do i=1,meshdatactx%Nvlocal
        do j=1,meshdatactx%itot(i)
            nn = meshdatactx%AdjM(j,i)
            if (modeldatactx%IDTglobal(nn).ne.0) then
                totadjeqst = totadjeqst + 1
            endif
            modeldatactx%AdjIDT(j,i) = modeldatactx%IDTglobal(nn)
        enddo
    enddo

    if (hasdiffthermal) then
        modeldatactx%AdjIDTX = 0
        do i=1,meshdatactx%Nvlocal
            do j=1,meshdatactx%itot(i)
                nn = meshdatactx%AdjM(j,i)
                if (modeldatactx%IDTXglobal(nn).ne.0) then
                    totadjeqsxt = totadjeqsxt + 1
                endif
                modeldatactx%AdjIDTX(j,i) = modeldatactx%IDTXglobal(nn)
            enddo
        enddo
    endif

    if (iecho.eq.8) then
        call ParallelLog("buildAdjancyArrays", "Finished thermal adjacency")
    endif


endif


end subroutine
!***************************************************************

subroutine createMeshNumbering(NEN, NDOF)
! after the vertices and elements have been read from TECIN.DAT,
! we can already number the mesh elements.
! This is necessary so that when the nodal boundary conditions
! are read, we know which points are the ghost points of each
! partition.

use meshdatamodule, only: meshdatactx
use modelctx,   only: getrank, getsize
use debugmodule,    only: AllocateError, &
                          iecho, &
                          ParallelLog, &
                          delay
use constants,  only: NNeighborsMax2D, &
                          NNeighborsMax3D
implicit none

integer :: NEN, NDOF
integer :: i, j
integer :: jstart
integer :: ERROR
integer :: nb

integer, allocatable :: tmpAdjM(:)
integer, allocatable :: tmpIen(:)


! we will copy the meshdata%ien and th meshdata$adjm into local
! arrays, where we will manipulate them, without influencing
! the main module arrays, which are needed later in the creation
! of the equation masking in subroutine createEquationNumbering.

! This should be nicely cleanly separated, but for now we are glad
! when it just works.
integer, allocatable :: localIen(:,:)
integer, allocatable :: localAdjm(:,:)
integer              :: maxnbrs

maxnbrs = meshdatactx%nMaxNeighbors

if (iecho.eq.8) then
    call ParallelLog("CreateMeshNum", "at beginning")
endif

allocate(localIen(NEN,meshdatactx%nelocal))
allocate(localAdjM(maxnbrs,meshdatactx%nvlocal))

localIen  = meshdatactx%ien
localAdjM = meshdatactx%AdjM

call FillLocInd()

! now that locind has been determined, we can build a gloInd2locInd vector.

call FillGlo2LocVertex()
call FillGlo2LocElement()
call FillFullGlo2LocElement()

!   element verieties are handle in subroutine ELNODE in tecin.dat
!    call FillLocElt()
!    call FillGlo2LocElement()

! we will copy the meshdata%ien and th meshdata$adjm into local
! arrays, where we will manipulate them, without influencing
! the main module arrays, which are needed later in the creation
! of the equation masking in subroutine createEquationNumbering.

! This should be nicely cleanly separated, but for now we are glad
! when it just works.

! alter adjancency
ALLOCATE(tmpAdjM(meshdatactx%Nvneighborstotal), STAT=ERROR)
call AllocateError("createMeshNumbering tmpAdjM", ERROR)


tmpAdjM = 0

jstart = 0
do i=1,meshdatactx%Nvlocal
    if (meshdatactx%itot(i).gt.maxnbrs) then
        write(*,*) "Something weird. Does partition.info contain the correct data?"
    endif
    do j=1,meshdatactx%itot(i)
        tmpAdjM(j+jstart) = localAdjM(j,i)-1
    enddo
    jstart = jstart  + meshdatactx%itot(i)
enddo

   ! aonp from module
call AOApplicationToPetsc(aonp, meshdatactx%Nvneighborstotal, &
tmpAdjM, ERROR)

jstart = 0
do i=1,meshdatactx%Nvlocal
    do j=1,meshdatactx%itot(i)
        localAdjM(j,i)= tmpAdjM(j+jstart)
    enddo
    jstart = jstart  + meshdatactx%itot(i)
enddo

deallocate(tmpAdjM)

!---------------------------------------------------------------------------

! alter element vertex numbers
ALLOCATE(tmpIen(meshdatactx%Nelocal*NEN), STAT=ERROR)
call AllocateError("createMeshNumbering tmpIen", ERROR)


do i=1,meshdatactx%Nelocal
    do j=1,NEN
        tmpIen((i-1)*NEN + j) = localIen(j,i)-1
    enddo
enddo


call AOApplicationToPetsc(aonp, meshdatactx%Nelocal*NEN, &
tmpIen, ERROR)

do i=1,meshdatactx%Nelocal
    do j=1,NEN
        localIen(j,i) = tmpIen((i-1)*NEN + j)
    enddo
enddo

deallocate(tmpIen)

!---------------------------------------------------------------
call CreateVerticesMask(meshdatactx,ndof, localIen, NEN, localAdjM, maxnbrs,1)

deallocate(localIen)
deallocate(localAdjM)

return

end subroutine



!***************************************************************
subroutine createEquationNumbering &
         (meshdat, modeldat, nen,ndof, nsd)
!***************************************************************
! Routine to create parallel numbering.
! Because this relates directly to the position 
! of entries in the stiffness matrix/load vector, this is 
! done in local numbering.
!***************************************************************
#ifdef EXPERIMENTAL_ReducedSlipAlloc
USE MODELDATAMODULE, only: modeldata, &
                           HasSlipperyNodes, &
                           elementSlipperyID
use meshdatamodule,  only: meshdatactx
#else
USE MODELDATAMODULE, only: modeldata, &
                           HasSlipperyNodes
#endif

USE MODELDEFINITION
USE MESHDATAMODULE,  only: meshdata
USE MODELCTX
use debugmodule,     only: iecho, &
                           parallelLog, &
                           delay, &
                           allocateError, &
                           debug
use constants,       only: NNeighborsMax2D, &
                           NNeighborsMax3D
!***************************************************************
   implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsysdef.h"
!#include "petsc/finclude/petscao.h"
!***************************************************************
   type (meshdata) :: meshdat
   type (modeldata) :: modeldat

   integer i, j, k,jstart, nb,n,nn,master,slave
   integer ERROR, nen,ndof, nsd, next,totelmeqs,totelmeqsx, &
 totadjeqs, totveqs, totadjeqsx
   integer totadjeqst, totadjeqsxt, totelmeqst, totelmeqsxt
   integer toteqs, totideqs, totidxeqs
   integer toteqst, totideqst, totidxeqst
   integer :: maxnbrs
!      logical :: DoStuffHere ! for debugging
PetscInt, ALLOCATABLE :: tmpeq(:),  &
                         tmpid(:),  &
                         tmpidx(:)
PetscInt, ALLOCATABLE :: tmpeqt(:),  &
                         tmpidt(:),  &
                         tmpidtx(:)
PetscInt, ALLOCATABLE :: eqordering(:), &
                         eqorderingt(:), &
                         tmp(:), &
                         tmp2(:), &
                         tmp3(:), &
                         tmp4(:), &
                         tmpfnnp(:) !, & ! which is what?
!                               npordering(:), &
!                               elordering(:)

integer :: globalElemID, slipperyElemID, iElem

PetscErrorCode ierr


integer :: idx, weight
integer :: slipSequenceNr, slipElemID, nodeSequenceNr
integer :: iSlipElemEntry, islipelem, globalNodeID, nSlipEntries


!***************************************************************


!write(*,*) "thread", getrank(), "entered createEquationNumbering"

call mpi_barrier(MPI_COMM_WORLD, ierr)

if (iecho.eq.8) then
    call ParallelLog("CreateParNum", "at beginning")
endif

   ! toteqs, tmpeq : relate to unique equation numbering on this processor
   ! totideqs, tmpid : related to equations numbers in ID array
   ! totidxeqs, tmpidx : related to equations numbers in IDX array

   ! here modeldat%ID had NO linked nodes yet (i.e. NEQlocal agrees
   ! with number of non-zeros on ID en IDX combined)
   ! the resulting array tmpeq contains only unique numbers


   ! Construct the local nunmbering of the vertices and elements.
!      call FillLocInd()
!      call FillLocElt()

!      call FillGlo2LocVertex()
!      call FillGlo2LocElement()

!write(*,*) "start rank", getrank(),"LMX: ", modeldat%lmx


! count number of equation of vertices belonging to this partition only.
! (no ghost points etc)
call countMechanicalEquations(tmpeq, toteqs, totideqs, totidxeqs, ndof)

!    write(*,*) "rank", getrank(),"totid(x)eqs", tmpeq, toteqs,totideqs, totidxeqs

if (iecho.eq.8) then
    call parallelLog("createEquationNumbering", "counted mechanical eqs")
endif


if (hasthermal) then
    call countThermalEquations(tmpeqt, toteqst, totideqst, totidxeqst)

!#ifdef SPARSE
!    do i=0,getsize()
!#endif
!        call delay()
!#ifdef SPARSE
!        call mpi_barrier(0,ierr)
!#endif

!        if (i.eq.getrank()) then
!            write(*,*) "thread", getrank(), "has counted thermal eqs: ", tmpeqt
!        endif

!#ifdef SPARSE
!    enddo
!#endif

endif



if (iecho.eq.8) then
    call parallelLog("createEquationNumbering", "counted thermal eqs")
endif

!      write(FILE_outputf,*) 'ID before all'
!      do i=1,meshdat%Nvlocal
!          write(FILE_outputf,*) (modeldat%ID(j,i),j=1,NDOF)
!      enddo
!      write(FILE_outputf,*) '------------------'

!     update (local) ID array to match IDglobal (corrects for
!     linked nodes) count TOTAL number of non-zeros at the same time

!if (debug) then
!    write(*,*) "rank", getrank(),"Nvlocal", meshdat%Nvlocal
!    write(*,*) "rank", getrank(),"gloInd", meshdat%gloInd
!    write(*,*) "rank", getrank(),"locInd", meshdat%locInd ! is 0-based here.
!    write(*,*) "rank", getrank(),"modeldat IDX global at begin CEN 1: ", modeldat%IDXglobal
!    write(*,*) "rank", getrank(),"modeldat ID global at begin CEN 1: ", modeldat%IDglobal
!endif

! count mechanical equations again, and use it

! to set the ID so that the global equation number 
! will be placed in the position of the global vertex indices
! of this partition.

! If partition x starts with vertice 8, 10, 103, ...
! and those have global eq nr 16,17, 20,21, 106,107 ...
! (maybe shifted a bit upward with repect to the 2n, 2n+1 scheme to make room for slippery nodes)
! Then ID will be:
! 16, 17, 20, 21, 106, 107...

do i=1,meshdat%Nvlocal
    do k=1,NDOF
        if (modeldat%IDglobal(k,meshdat%gloInd(i)).ne.0) then
            totideqs = totideqs + 1
        endif
        modeldat%ID(k,i) = modeldat%IDglobal(k,meshdat%gloInd(i))
    enddo
enddo


if (debug) then
!    write(*,*) 'rank ', getrank(), ' has base ID ', modeldat%ID
endif

if (iecho.eq.8) then
    call parallelLog("createEquationNumbering", "counting eqs")
endif


! Count again
! results of these counts are local, and as such differ per thread/partition
! Appears a bit double,
! What this subroutine mainly does, it is set the LM, LMX, LMT, LMXT arrays!
!
! Not the subtle notation 'totel'... ipv 'total'
! This is short for total over the elements.
! A vertex belonging to six elements, will have it equations counted six times as well.
! As such totelmeqs, totelmeqsx, totelmeqst, totelmeqsxt will always be suspicially high,
! although that is is by design and OK.
call CountEquations &
   (NEN, NDOF, totelmeqs, totelmeqsx, totelmeqst, totelmeqsxt)


if (iecho.eq.8) then
    call parallelLog("createEquationNumbering", "building adjacency array")
endif

call buildAdjancyArrays &
   (NSD, NDOF, totadjeqs, totadjeqsx, totadjeqst, totadjeqsxt)



!     numbering from 0 is required because for the AOCreateBasic routine
!     the complete(! set of points to be mapped, say N, need to run
!     from 0 up to an including N-1
do i=1,meshdat%Nvlocal
    meshdat%gloInd(i) = meshdat%gloInd(i) - 1
enddo

do i=1,meshdat%Nelocal
    meshdat%gloElt(i) = meshdat%gloElt(i) - 1
enddo

if (iecho.eq.8) then
    call parallelLog("createEquationNumbering", "Gathering element offsets")
endif



!     equations
! determine a continuing count over the equations of all threads
! example, let neqlocal over 4 procs be:  10   12   9   15
! rstart will be   10    22    31    46
call MPI_Scan(NEQlocal, rstart, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

! rstart indicates the start of the beginning of the next partition.
! (MPI_Scan returns 'total up until this partition')
! Substract number of own equations do indicate start of own partition in the global sequence
rstart = rstart - NEQlocal




!if (debug) then
!    write(*,*) "rank", getrank(), "found rstart: ", rstart
!endif

ALLOCATE(eqordering(NEQlocal), STAT=ERROR)
call allocateError("eqordering", ERROR)


! eqordering becomes a local array of the equation numbers over all processors.
! That inclused the slippery and
do i=1,NEQlocal
    eqordering(i) = rstart + i - 1
enddo

!if (debug) then
!    write(*,*) "rank", getrank(), "make eqordering: ", eqordering
!endif


if (hasthermal) then

    ! Thermal has it own rstart, because it has its own matrix/vector system.
    ! As such, it does interfere with the regular rstart
    call MPI_Scan(NTEQlocal, rstart, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)
    rstart = rstart - NTEQlocal
    ALLOCATE(eqorderingt(NTEQlocal), STAT=ERROR)
    call allocateError("eqorderingt", ERROR)
    eqorderingt = 0


!    write(*,*) "rank", getrank(), "has NTEQlocal", NTEQlocal

    do i=1,NTEQlocal
        eqorderingt(i) = rstart + i - 1
    enddo

endif

!     create application ordering from user defined numbering
!     to petsc defined numbering:
!     user numbering: 1 ... MAX random over all processors
!     petsc numbering: 1 ... n1 on processor 1,
!                     n1+1 ... n2 on processor 2

if (iecho.eq.8) then
    call parallelLog("createEquationNumbering", "Creating aoeq lookup table")
endif

! make aoeq, based on the eqordering determined above

call AOCreateBasic(MPI_COMM_WORLD, NEQlocal, tmpeq, eqordering, aoeq, ierr)

if (ierr.ne.0) then
    write(*,*) "AOmodule could not AOcreatebasic in subroutine createEquationNumbering, error: ", ierr
endif

CHKERRQ(ierr)



if (iecho.eq.8) then
    call parallelLog("createEquationNumbering", "create equation ordering")
endif


!write(*,*) "rank", getrank(), "has tmpeqt", tmpeqt
!write(*,*) "rank", getrank(), "has eqorderingt", eqorderingt



if (hasthermal) then
    call AOCreateBasic(MPI_COMM_WORLD, NTEQlocal, tmpeqt, eqorderingt, aoeqt,ierr)
endif

!call AOview(aoeqt, PETSC_VIEWER_STDOUT_WORLD, ierr )
!call delay()


! no longer required, as we proceed with aoeq (and possibly aoeqt), from here
DEALLOCATE(eqordering, STAT=ERROR)
DEALLOCATE(tmpeq, STAT=ERROR)
if (hasthermal) then
    DEALLOCATE(eqorderingt, STAT=ERROR)
    DEALLOCATE(tmpeqt, STAT=ERROR)
endif

ALLOCATE(tmp(meshdat%Nvneighborstotal), STAT=ERROR)
call allocateError("tmp", ERROR)

ALLOCATE(tmp2(meshdat%Nelocal*NEN), STAT=ERROR)
call allocateError("tmp2", ERROR)


!write(*,*) "rank", getrank(), "has AdjM", meshdat%AdjM

! tmp becomes an array with all the neighbor indices stored in it sequentially
jstart = 0
do i=1,meshdat%Nvlocal
    do j=1,meshdat%itot(i)
        tmp(j+jstart) = meshdat%AdjM(j,i)-1
    enddo
    jstart = jstart  + meshdat%itot(i)
enddo

!write(*,*) "rank", getrank(), "has jstart", jstart      
!write(*,*) "rank", getrank(), "has tmp", tmp


!     copy element vertex numbers
do i=1,meshdat%Nelocal
    do j=1,NEN
        tmp2((i-1)*NEN + j) = meshdat%ien(j,i)-1
    enddo
enddo

!write(*,*) "rank", getrank(), "has ien", meshdat%ien



!      call writeAOmeshdata(FILE_outputf, myrank, meshdat, tmp, nen, 0)
!      do i=1,meshdat%Nelocal
!          write(FILE_outputf,*) "Equation numbers of local element ",
!     .      meshdat%locElt(i), " are : ",
!     >      ((modeldat%lm(k,j,i), k=1,2), j=1,NEN)
!      enddo
!      write(FILE_outputf,'(8I5)') (tmp3(i),i=1,totelmeqs)

!     create Application to Petsc orderings
!      call AOApplicationToPetsc(aonp, meshdat%Nvlocal, &
!          meshdat%locInd, ierr)


if (iecho.eq.8) then
    call parallelLog("AOmod", "ordering App to PETSc")
endif


! aonp is the ordering made above.
! tmp(2) is being reordered according to this context,
! and returned in its reordered state.
call AOApplicationToPetsc(aonp, meshdat%Nvneighborstotal, &
   tmp, ierr)

call AOApplicationToPetsc(aonp, meshdat%Nelocal*NEN, &
   tmp2, ierr)




if (iecho.eq.8) then
    call parallelLog("AOmod", "ordered App to PETSc")
endif


! linked nodes, little less generic because of double array LINK
if (NLINK.gt.0) then
!          write(FILE_outputf,*) 'linked nodes before'
!          do i=1,NLINK
!              write(FILE_outputf,*) i, (modeldat%LINK(j,i),j=2,3)
!          enddo
    ALLOCATE(tmpfnnp(NLINK), STAT = ERROR)
    call allocateError("tmpfnnp", ERROR)

    do i=1,NLINK
        tmpfnnp(i) = modeldat%LINK(2,i) - 1
    enddo

    call AOApplicationToPetsc(aonp,NLINK, tmpfnnp,ierr)

    do i=1,NLINK
        modeldat%LINK(2,i) = tmpfnnp(i)
    enddo

    do i=1,NLINK
        tmpfnnp(i) = modeldat%LINK(3,i) - 1
    enddo

    call AOApplicationToPetsc(aonp,NLINK, tmpfnnp,ierr)

    do i=1,NLINK
        modeldat%LINK(3,i) = tmpfnnp(i)
    enddo

!          write(FILE_outputf,*) 'linked nodes after'
!          do i=1,NLINK
!              write(FILE_outputf,*) i, (modeldat%LINK(j,i),j=2,3)
!          enddo
    DEALLOCATE(tmpfnnp, STAT=ERROR)
endif


! slippery nodes, little less generic because of double array
! NSLIP
if (HasSlipperyNodes()) then
    if (iecho.eq.8) then
        call parallelLog("AOmod", "starting slip special")
    endif

!    write(*,*) "rank", getrank(), "has numslp", numslp


    ALLOCATE(tmpfnnp(NUMSLP), STAT = ERROR)
    call allocateError("tmpfnnp", ERROR)

    do i=1,NUMSLP
        tmpfnnp(i) = modeldat%NSLIP(2,i)-1
    enddo

!    write(*,*) "rank", getrank(), "has tmpfnnp", tmpfnnp

    call AOApplicationToPetsc(aonp, NUMSLP, tmpfnnp, ierr)

!    write(*,*) "rank", getrank(), "has reordered tmpfnnp", tmpfnnp       


    do i=1,NUMSLP
        modeldat%NSLIP(2,i) = tmpfnnp(i)
    enddo

    DEALLOCATE(tmpfnnp)

    call renumbernodalpoints(aonp, modeldat%IDSLP, NUMSN, ierr)

    if (iecho.eq.8) then
        call parallelLog("AOmod", "done slip special")
    endif

endif



if (NUMFN.gt.0) then
    ! faulted nodes; not faulted parallel elems.
    ALLOCATE(tmpfnnp(numfn), STAT = ERROR)
    call allocateError("tmpfnnp", ERROR)

    do i=1,NUMFN
        tmpfnnp(i) = modeldat%NFAULT(2,i)-1
    enddo
    call AOApplicationToPetsc(aonp, NUMFN, &
    tmpfnnp, ierr)
    do i=1,NUMFN
        modeldat%NFAULT(2,i) = tmpfnnp(i)
    enddo
    DEALLOCATE(tmpfnnp)
    if (iecho.eq.8) then
        call parallelLog("AOmod", "done fault special")
    endif

endif

! surface nodes
if (NSURF.gt.0) then
    call renumbernodalpoints(aonp, modeldat%ISURF, NSURF, ierr)
endif

! thermal anomaly data, little less generic because of double array
! ITANOM
if (NTANOM.gt.0) then
    ALLOCATE(tmpfnnp(NTANOM), STAT = ERROR)
    call allocateError("tmpfnnp", ERROR)

    do i=1,NTANOM
        tmpfnnp(i) = modeldat%ITANOM(1,i)-1
    enddo
    call AOApplicationToPetsc(aonp, NTANOM, &
    tmpfnnp, ierr)
    do i=1,NTANOM
        modeldat%ITANOM(1,i) = tmpfnnp(i)
    enddo
    DEALLOCATE(tmpfnnp)
endif



!    if (iecho.eq.8) then
!        call parallelLog("AOmod", "Destroying aonp")
!    endif

!    call mpi_barrier(MPI_COMM_WORLD, error)

!   call AODestroy(aonp, ierr)


!    if (iecho.eq.8) then
!        call parallelLog("AOmod", "Destroyed aonp")
!    endif

!    call mpi_barrier(MPI_COMM_WORLD, error)




!      call AOApplicationToPetsc(aoel, meshdat%Nelocal, &
!          meshdat%locElt, ierr)

!      call AODestroy(aoel, ierr)

!      if (NUMSLP.gt.0) then
!          write(FILE_outputf,*) totidxeqs, meshdat%Nvlocal, NUMSLP
!          do i=1,meshdat%Nvlocal
!              write(FILE_outputf,*) (modeldat%IDX(k,i),k=1,NDOF),' --- ', &
!                             (modeldat%ID(k,i),k=1,NDOF)
!          enddo
!      else
!          do i=1,meshdat%Nvlocal
!              write(FILE_outputf,*) (modeldat%ID(k,i),k=1,NDOF)
!          enddo
!      endif

if (iecho.eq.8) then
     call parallelLog("AOmod", "renumbering normal equations")
endif

call renumberequations(aoeq, modeldat%ID, size(modeldat%ID,2), &
                          totideqs, &
            meshdat%Nvlocal,NDOF, ierr,FILE_outputf,"normal")


if (iecho.eq.8) then
    call parallelLog("AOmod", "renumbered normal equations")
endif

if (HasSlipperyNodes()) then

    if (iecho.eq.8) then
        call parallelLog("AOmod", "slippery eqs renumbering")
    endif

    call renumberequations(aoeq, modeldat%IDX, size(modeldat%IDX,2), &
                           totidxeqs, &
            meshdat%Nvlocal,NDOF, ierr,FILE_outputf,"slippery")

    if (iecho.eq.8) then
        call parallelLog("AOmod", "done slippery eqs renumbering")
    endif
endif


if (hasthermal) then 
!    write(*,*) "calling renum eq with totideqst", totideqst, meshdat%Nvlocal
!    write(*,*) "size 1 of IDT: ", size(modeldat%IDT,1)
!    write(*,*) "size 2 of IDT: ", size(modeldat%IDT,2)

    call renumberThermalEquations(aoeqt, &
                           modeldat%IDT, size(modeldat%IDT,1), &
                           totideqst, &
      meshdat%Nvlocal,1, ierr,FILE_outputf,"thermal")
endif

if (hasdiffthermal) then
    call renumberThermalEquations(aoeqt, &
                           modeldat%IDTX, size(modeldat%IDTX,1), &
                           totidxeqst, &
                 meshdat%Nvlocal,1, ierr,FILE_outputf,"diffthermal")
endif


!      if (parout) then
!          write(FILE_outputf,*) 'summary of equation numbers (ID arrays)', &
!           ' after renumbering'
!          if (NUMSLP.ne.0) then
!              write(FILE_outputf,*) 'IDX and ID:'
!              do i=1,meshdat%Nvlocal
!                  write(FILE_outputf,*) i,(modeldat%IDX(j,i),j=1,NDOF), &
!              ' --- ', (modeldat%ID(j,i),j=1,NDOF)
!              enddo
!          else
!              write(FILE_outputf,*) 'ID'
!              do i=1,meshdat%Nvlocal
!                  write(FILE_outputf,*) i,(modeldat%ID(j,i),j=1,NDOF)
!              enddo
!          endif

!          if (hasthermal) then
!              if (hasdiffthermal) then
!                  write(FILE_outputf,*) 'IDTX and IDT:'
!                  do i=1,meshdat%Nvlocal
!                       write(FILE_outputf,*) modeldat%IDTX(i), ' --- ', &
!                          modeldat%IDT(i)
!                  enddo
!              else
!                  write(FILE_outputf,*) 'IDT:'
!                  do i=1,meshdat%Nvlocal
!                       write(FILE_outputf,*) modeldat%IDT(i)
!                  enddo
!              endif
!          endif
!      endif

if (iecho.eq.8) then
!     call parallelLog("AOmod", "building tmp4")
endif

!     fill array with equations numbers for neighbouring equations
ALLOCATE(tmp4(totadjeqs), STAT=ERROR)
call allocateError("tmp4", ERROR)
next = 0
do i=1,meshdat%Nvlocal
    do j=1,meshdat%itot(i)
        do k=1,NDOF
            if (modeldat%AdjID(k,j,i).ne.0) then
                next = next + 1
                tmp4(next) = modeldat%AdjID(k,j,i)-1
            endif
        enddo
    enddo
enddo

if (iecho.eq.8) then
!        call parallelLog("AOmod", "finished tmp4")
endif

call AOApplicationToPetsc(aoeq, totadjeqs, tmp4, ierr)

if (iecho.eq.8) then
        call parallelLog("AOmod", "building adjID")
endif

! alter equation numbers
next = 0
do i=1,meshdat%Nvlocal
    do j=1,meshdat%itot(i)
        do k=1,NDOF
            if (modeldat%AdjID(k,j,i).ne.0) then
                next = next + 1
                modeldat%AdjID(k,j,i) = tmp4(next) + 1
            endif
        enddo
    enddo
enddo

DEALLOCATE(tmp4, STAT=ERROR)


if (iecho.eq.8) then
    call parallelLog("AOmod", "renumbering slippery nodes")
endif

!write(*,*) "hasslip", HasSlipperyNodes()

if (HasSlipperyNodes()) then

    ALLOCATE(tmp4(totadjeqsx), STAT=ERROR)
    call allocateError("tmp4", ERROR)

!    write(*,*) "rank",getrank(),'doing loop 4'


    if (iecho.eq.12) then
        write(*,*) "AdjIDX before reforming: ", modeldat%AdjIDX
    endif

    next = 0
    do i=1,meshdat%Nvlocal
        do j=1,meshdat%itot(i)
            do k=1,NDOF
                if (modeldat%AdjIDX(k,j,i).ne.0) then
                    next = next + 1
                    tmp4(next) = modeldat%AdjIDX(k,j,i)-1
                endif
            enddo
        enddo
    enddo


!    write(*,*) "rank",getrank(),'done loop 4'

    if (iecho.eq.8) then
        call parallelLog("AOmod", "AOApp2PETSc running on tmp4")
    endif


    call AOApplicationToPetsc(aoeq, totadjeqsx, &
                              tmp4, ierr)


    if (iecho.eq.8) then
        call parallelLog("AOmod", "AOApp2PETSc ran on tmp4")
    endif


    ! alter equation numbers
    next = 0
    do i=1,meshdat%Nvlocal
        do j=1,meshdat%itot(i)
            do k=1,NDOF
                if (modeldat%AdjIDX(k,j,i).ne.0) then
                    next = next + 1
                    modeldat%AdjIDX(k,j,i) = tmp4(next) + 1
                endif
            enddo
        enddo
    enddo
    DEALLOCATE(tmp4, STAT=ERROR)

    if (iecho.eq.12) then
        write(*,*) "AdjIDX after reforming: ", modeldat%AdjIDX
    endif



    if (iecho.eq.8) then
        call parallelLog("AOmod", "slippery nodes renumbered")

    endif


else
!    write(*,*) 'no slippery nodes'

endif

if (iecho.eq.8) then
    call parallelLog("AOmod", "renumbered slippery nodes")
endif

!    write(*,*) "rank",getrank(),"WAAAARGH"


if (hasthermal) then
!     fill array with equations numbers for neighbouring equations
!          write(FILE_outputf,*) 'ADJIDT'
!          do i=1,meshdat%Nvlocal
!              write(FILE_outputf,*) i, (modeldat%AdjIDT(j,i),
!     .            j=1,meshdat%itot(i))
!          enddo
    ALLOCATE(tmp4(totadjeqst), STAT=ERROR)
    call allocateError("tmp4", ERROR)
    next = 0
    do i=1,meshdat%Nvlocal
        do j=1,meshdat%itot(i)
            if (modeldat%AdjIDT(j,i).ne.0) then
                next = next + 1
                tmp4(next) = modeldat%AdjIDT(j,i)-1
            endif
        enddo
    enddo
    call AOApplicationToPetsc(aoeqt, totadjeqst, &
                              tmp4, ierr)
    ! alter equation numbers
    next = 0
    do i=1,meshdat%Nvlocal
        do j=1,meshdat%itot(i)
            if (modeldat%AdjIDT(j,i).ne.0) then
                next = next + 1
                modeldat%AdjIDT(j,i) = tmp4(next) + 1
            endif
        enddo
    enddo
!          do i=1,meshdat%Nvlocal
!              write(FILE_outputf,*) i, (modeldat%AdjIDT(j,i),
!     .            j=1,meshdat%itot(i))
!          enddo
    DEALLOCATE(tmp4, STAT=ERROR)

endif


if (hasdiffthermal) then
!     fill array with equations numbers for neighbouring equations
    ALLOCATE(tmp4(totadjeqsxt), STAT=ERROR)
    call allocateError("tmp4", ERROR)
    next = 0
    do i=1,meshdat%Nvlocal
        do j=1,meshdat%itot(i)
            if (modeldat%AdjIDTX(j,i).ne.0) then
                next = next + 1
                tmp4(next) = modeldat%AdjIDTX(j,i)-1
            endif
        enddo
    enddo
    call AOApplicationToPetsc(aoeqt, totadjeqsxt, &
                              tmp4, ierr)
    ! alter equation numbers
    next = 0
    do i=1,meshdat%Nvlocal
        do j=1,meshdat%itot(i)
            if (modeldat%AdjIDTX(j,i).ne.0) then
                next = next + 1
                modeldat%AdjIDTX(j,i) = tmp4(next) + 1
            endif
        enddo
    enddo
    DEALLOCATE(tmp4, STAT=ERROR)
endif ! hasdiffthermal



!   copy element equation numbers
ALLOCATE(tmp3(totelmeqs), STAT=ERROR)
call allocateError("tmp3", ERROR)

next = 0
do i=1,meshdat%Nelocal
    do j=1,NEN
        do k=1,NDOF
            if (modeldat%lm(k,j,i).ne.0) then
                next = next + 1
                tmp3(next) = modeldat%lm(k,j,i)-1
            endif
        enddo
    enddo
enddo


call AOApplicationToPetsc(aoeq, totelmeqs,tmp3, ierr)


next = 0
do i=1,meshdat%Nelocal
    do j=1,NEN
        do k=1,NDOF
            if (modeldat%lm(k,j,i).ne.0) then
                next = next + 1
                modeldat%lm(k,j,i) = tmp3(next)+1

!                modeldat%lm(k,j,i) = equations(tmp3(next)+1)+1
            endif
        enddo
    enddo
enddo



 DEALLOCATE(tmp3, STAT=ERROR)

! HasSlipperyNodes() gives a global return
! Even partitions without slippery nodes will go into 
! this if block if other partitions do.
if (HasSlipperyNodes()) then




    if (iecho.eq.8) then
        call parallelLog("AOmod", "building LMX")
    endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
    ! we have the lookup table, generated from
    ! global numbering, while reading the slippery elements.
    ! Also, we have the LMX table in the sequence 
    ! in which they were read

    ! Now we have the local element numbering.
    ! For now, we will reorder the global LMX array
    ! to match the sequence of the local ordering.

    ! the equations array is not yet filled, here.
    ! after all, we are creating the equation numbers, and 
    ! the equations array can only be constructed once that
    ! has been completed.

    ! same goes for the equationsmask(x).


!write(*,*) "rank", getrank()," has locInd: ", meshdatactx%locInd
!write(*,*) "rank", getrank()," has gloInd: ", meshdatactx%gloInd

!write(*,*) "rank", getrank()," has prevamp LMX ", modeldat%lmx



!if (getrank().eq.2) then
!    write(*,*) "rank", getrank(),"gloelt: ", meshdatactx%gloelt
!    write(*,*) "rank", getrank(),"locelt: ", meshdatactx%locelt
!    write(*,*) "rank", getrank(),"nslipglobal: ", modeldat%nslipglobal
!    write(*,*) "rank", getrank(),"idxGlobal: ", modeldat%IDXglobal
!    write(*,*) "rank", getrank(),"global IEN: ", meshdatactx%globalIen
!    write(*,*) "rank", getrank(),"equationsmaskx :", equationsmaskx   ! is empty
!    do iElem =1,size(modeldat%lmx,3)
!        write(*,*) "rank",getrank(),"LMX", iElem, ":", modeldat%lmx(:,:,iElem)
!    enddo
!    write(*,*) "rank", getrank(),"LMX: ", modeldat%lmx
!    write(*,*) "rank", getrank(),"slimelepmLookup",modeldat%SlipElemLookup
!    write(*,*) "rank", getrank(),"NSLIP", modeldat%NSLIP


!    write(*,*) "sizes: ", size(modeldat%nslipglobal,1), size(modeldat%nslipglobal,2)

    nSlipEntries = size(modeldat%nslipglobal,2)
    do iSlipElemEntry =1,nSlipEntries ! number of slippery node entries in tecin
                       ! and as such entries in nslipglobal


        slipElemID = modeldat%nslipglobal(1,iSlipElemEntry)

!        write(*,*) "check slip elem entry", iSlipElemEntry, "which has elem", slipElemID

        if (modeldat%SlipElemLookup(slipElemID).gt.0) then
            slipSequenceNr = modeldat%SlipElemLookup(slipElemID)

!            write(*,*) "Element ", slipElemID,"has sequence nr", slipSequenceNr

            ! we have the global node nuer from nslipglobal.
            ! determine whether this is node 1,2,3 or 4 of the element
            ! using the shameful gloalIEN we have constructed in subroutine ELNODE
            globalNodeID = modeldat%nslipglobal(2,iSlipElemEntry)


!            write(*,*) "and global node ", globalNodeID
!            write(*,*) "looking up this node ID in ", meshdatactx%globalIEN(:,slipElemID)

            nodeSequenceNr = 0
            do i = 1,3
                if (globalNodeID.eq.meshdatactx%globalIEN(i,slipElemID)) then
                    nodeSequenceNr = i
                endif
            enddo
            if (NDOF.eq.3) then
                if (globalNodeID.eq.meshdatactx%globalIEN(4,slipElemID)) then
                    nodeSequenceNr = 4
                endif
            endif

!            write(*,*) "its node", globalNodeID, "is in position", nodeSequenceNr

            if (nodeSequenceNr.ge.1 .and. nodeSequenceNr.le.4) then
                if      (modeldat%nslipglobal(3,iSlipElemEntry).ne.0) then ! first weight

                    idx = modeldat%IDXglobal(1,globalNodeID)
                    weight = modeldat%NSLIPglobal(3,iSlipElemEntry)

                    modeldat%lmx(1,nodeSequenceNr,slipSequenceNr) = &
                        ISIGN(idx, weight)
                endif

                if (modeldat%nslipglobal(4,iSlipElemEntry).ne.0) then ! second weight

                    idx = modeldat%IDXglobal(2,globalNodeID)
                    weight = modeldat%NSLIPglobal(4,iSlipElemEntry)
                    
                    modeldat%lmx(2,nodeSequenceNr,slipSequenceNr) = &
                        ISIGN(idx, weight)
                endif
                
                if (NDOF.eq.3) then
                    if (modeldat%nslipglobal(5,iSlipElemEntry).ne.0) then ! third weight
    
                        idx = modeldat%IDXglobal(3,globalNodeID)
                        weight = modeldat%NSLIPglobal(5,iSlipElemEntry)
                    
                        modeldat%lmx(3,nodeSequenceNr,slipSequenceNr) = &
                            ISIGN(idx, weight)
                    endif
                endif

            else
                write(*,*) "ERROR, could nog find slippery node ",globalNodeID
                write(*,*) "in element", iSlipElem
                write(*,*) "nodeSequenceNr should be in [1,4] but is:", nodeSequenceNr
                stop "Leaving GTecton"    
            endif

        else
            write(*,*) "ERROR, nSlipGloal thinks that",slipElemID,"is a slippery element"
            write(*,*) "but SlipElemLookup disagree."
            write(*,*) "Please contact model support"
            stop "Leaving GTecton"
        endif

    enddo

    deallocate(meshdatactx%globalIEN)

!    do iElem =1,size(modeldat%lmx,3)
!        write(*,*) "rank",getrank(),"LMX", iElem, ":", modeldat%lmx(:,:,iElem)
!    enddo
!endif ! rank 2. Just for debugging


!    stop "So long"

!write(*,*) "rank", getrank()," has revamped LMX ", modeldat%lmx


#else

    ALLOCATE(tmp3(totelmeqsx), STAT=ERROR)

    ! noe thar
    call allocateError("tmp3", ERROR)

    next = 0

    do i=1,meshdat%Nelocal
        do j=1,NEN
            do k=1,NDOF
                if (modeldat%lmx(k,j,i).ne.0) then
                    next = next + 1
                    tmp3(next) = ABS(modeldat%lmx(k,j,i))-1
                endif
            enddo
        enddo
    enddo

    call AOApplicationToPetsc(aoeq,totelmeqsx,tmp3,ierr)
    next = 0

    do i=1,meshdat%Nelocal
        ! find global number from local number
        ! and find slippery ID from global number.
        do j=1,NEN
            do k=1,NDOF
                if (modeldat%lmx(k,j,i).ne.0) then
                    next = next + 1
                    modeldat%lmx(k,j,i) = &
                    ISIGN(tmp3(next)+1,modeldat%lmx(k,j,i))
                endif
            enddo
        enddo


    enddo

    DEALLOCATE(tmp3,STAT=ERROR)


    if (iecho.eq.8) then
        call parallelLog("AOmod", "done LMX building")
    endif


#endif







endif




!----- thermal ----------------

if (hasthermal) then
    ALLOCATE(tmp3(totelmeqst), STAT=ERROR)
    call allocateError("tmp3", ERROR)
    next = 0
    do i=1,meshdat%Nelocal
        do j=1,NEN
            if (modeldat%lmt(j,i).ne.0) then
                next = next + 1
                tmp3(next) = modeldat%lmt(j,i)-1
            endif
        enddo
    enddo
    call AOApplicationToPetsc(aoeqt, totelmeqst,tmp3, ierr)
    next = 0
    do i=1,meshdat%Nelocal
        do j=1,NEN
            if (modeldat%lmt(j,i).ne.0) then
                next = next + 1
                modeldat%lmt(j,i) = tmp3(next)+1
            endif
        enddo
    enddo
    DEALLOCATE(tmp3, STAT=ERROR)
endif


if (hasdiffthermal) then
    ALLOCATE(tmp3(totelmeqsxt), STAT=ERROR)
    call allocateError("tmp3", ERROR)
    next = 0
    do i=1,meshdat%Nelocal
        do j=1,NEN
            if (modeldat%lmtx(j,i).ne.0) then
                next = next + 1
                tmp3(next) = ABS(modeldat%lmtx(j,i))-1
            endif
        enddo
    enddo
    call AOApplicationToPetsc(aoeqt, totelmeqsxt,tmp3, ierr)
    next = 0
    do i=1,meshdat%Nelocal
        do j=1,NEN
            if (modeldat%lmtx(j,i).ne.0) then
                next = next + 1
                modeldat%lmtx(j,i) = &
      ISIGN(tmp3(next)+1, modeldat%lmtx(j,i))
            endif
        enddo
    enddo
    DEALLOCATE(tmp3, STAT=ERROR)
endif

jstart = 0

do i=1,meshdat%Nvlocal
    do j=1,meshdat%itot(i)
        meshdat%AdjM(j,i)= tmp(j+jstart)
    enddo
    jstart = jstart  + meshdat%itot(i)
enddo



! alter element vertex numbers
do i=1,meshdat%Nelocal
    do j=1,NEN
        meshdat%ien(j,i) = tmp2((i-1)*NEN + j)
    enddo
enddo


if (NSD.eq.2) then
    maxnbrs = NNeighborsMax2D
else
    maxnbrs = NNeighborsMax3D
endif


call CreateEquationsMask(meshdat, modeldat,NDOF)

!   if (iecho.eq.8) then
!       call parallelLog("AOmod", "equationsmask finished")
!   endif




   ! change everything back into base 1 numbering
do i=1,meshdat%Nvlocal
    meshdat%gloInd(i) = meshdat%gloInd(i)+1
    meshdat%locInd(i) = meshdat%locInd(i)+1
enddo

do i=1,meshdat%Nelocal
    meshdat%gloElt(i) = meshdat%gloElt(i)+1
    meshdat%locElt(i) = meshdat%locElt(i)+1
enddo



!      if (parout) then
!          write(FILE_outputf,*) 'Overview of vertices and equation masks:'
!          write(FILE_outputf,*) 'verticesmask'
!          do i=1,meshdat%Nvglobal
!              write(FILE_outputf,*) i, verticesmask(i)
!          enddo
!          write(FILE_outputf,*) 'vertices ',
!     .         (vertices(i),i=1,meshdat%Nvglobal)
!          write(FILE_outputf,*) 'equations ',
!     .         (equations(i),i=1,NDOF*meshdat%Nvglobal)

!          if (NUMSLP.gt.0) then
!              write(FILE_outputf,*) 'equationsmaskx and equationsmask:'
!              do i=1,meshdat%Nvglobal
!                  write(FILE_outputf,*) i, (equationsmaskx(j,i),j=1,NDOF), &
!                ' --- ', (equationsmask(j,i),j=1,NDOF)
!              enddo
!          else
!              write(FILE_outputf,*) 'equationsmask:'
!              do i=1,meshdat%Nvglobal
!                  write(FILE_outputf,*) i, (equationsmask(j,i),j=1,NDOF)
!              enddo
!          endif
!          if (hasthermal) then
!              if (hasdiffthermal) then
!                  write(FILE_outputf,*) 'equationsmaskxt and equationsmaskt:'
!                  do i=1,meshdat%Nvglobal
!                      write(FILE_outputf,*) i, equationsmaskxt(i), ' --- ', &
!                        equationsmaskt(i)
!                  enddo
!              else
!                  write(FILE_outputf,*) 'equationsmaskt:'
!                  do i=1,meshdat%Nvglobal
!                      write(FILE_outputf,*) i, equationsmaskt(i)
!                  enddo
!              endif
!              write(FILE_outputf,*) 'equationst ',
!     .         (equationst(i),i=1,meshdat%Nvglobal)
!          endif
!      endif


do i=1,meshdat%Nvlocal
    nb = meshdat%LocInd(i)
    do j=1,NDOF
        modeldat%ID(j,i) = equationsmask(j,nb)
    enddo
enddo



   if (HasSlipperyNodes()) then
  do i=1,meshdat%Nvlocal
      nb = meshdat%LocInd(i)
      do j=1,NDOF
          modeldat%IDX(j,i) = equationsmaskx(j,nb)
      enddo
  enddo
   endif


!write(*,*) "thread", getrank(), "has thermal: ", hasthermal

if (hasthermal) then
    do i=1,meshdat%Nvlocal
        nb = meshdat%LocInd(i)
        modeldat%IDT(i) = equationsmaskt(nb)
    enddo
!    write(*,*) "thread", getrank(), "has locind ", meshdat%LocInd
!    write(*,*) "thread", getrank(), "has equationsmaskt ", equationsmaskt
!    write(*,*) "thread", getrank(), "has modeldat%IDT", modeldat%IDT
endif

   if (hasdiffthermal) then
  do i=1,meshdat%Nvlocal
      nb = meshdat%LocInd(i)
      modeldat%IDTX(i) = equationsmaskxt(nb)
  enddo
   endif

   do i=1,meshdat%Nelocal
  do j=1,NEN
      nb = meshdat%ien(j,i)+1
      do k=1,NDOF
          if (modeldat%lm(k,j,i).ne.0) then
              modeldat%lm(k,j,i) = equationsmask(k,nb)
          endif
      enddo
  enddo
   enddo


!   if (iecho.eq.8) then
!       call parallelLog("AOmod", "building LMX")
!   endif


!---------- slippery ------------------
if (HasSlipperyNodes()) then

    do i=1,meshdat%Nelocal

#ifdef EXPERIMENTAL_ReducedSlipAlloc
!        write(*,*) "gloElt: ", meshdatactx%gloElt
        globalElemID = meshdatactx%gloElt(i) ! add this point gloelt has been rescaled to start at 1
                                             ! so there is no longer a need to +1
!        write(*,*) "gives global ID", globalElemID
        slipperyElemID = elementSlipperyID(globalElemID,2) ! funtion call, not an array
!        write(*,*) "gives slippery ID", slipperyElemID
#endif


        do j=1,NEN
!            write(*,*) 'checking ', j, i
            nb = meshdat%ien(j,i)+1
!            write(*,*) 'checked', nb
            do k=1,NDOF
!                write(*,*) 'LMX: ', k, j, i
!                write(*,*) 'is ', modeldat%lmx(k,j,i)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
                 ! lmx alays non-zero in new design

!                if (modeldat%lmx(k,j,i).ne.0) then
!                    write(*,*) 'equationmask ', k, nb, size(equationsmaskx,1),size(equationsmaskx,2)
!                    write(*,*) 'is ', equationsmaskx(k,nb)
!!!                    modeldat%lmx(k,j,slipperyElemID) = &
!!!                    ISIGN(equationsmaskx(k,nb),modeldat%lmx(k,j,slipperyElemID))
!                endif
#else
                if (modeldat%lmx(k,j,i).ne.0) then
!                    write(*,*) 'equationmask ', k, nb, size(equationsmaskx,1),size(equationsmaskx,2)
!                    write(*,*) 'is ', equationsmaskx(k,nb)
                    modeldat%lmx(k,j,i) = &
                    ISIGN(equationsmaskx(k,nb),modeldat%lmx(k,j,i))
                endif
#endif
            enddo
        enddo
!        write(*,*) 'rank', getrank(), 'elem', i, 'late has LMX', modeldat%lmx(:,:,i)
    enddo
endif


!write(*,*) "rank", getrank(),"LMX completed"


!if (iecho.eq.8) then
!     call parallelLog("AOmod", "LMX completed")
!endif


if (hasthermal) then
    do i=1,meshdat%Nelocal
        do j=1,NEN
            nb = meshdat%ien(j,i)+1
            if (modeldat%lmt(j,i).ne.0) then
                modeldat%lmt(j,i) = equationsmaskt(nb)
            endif
        enddo
    enddo
endif

   if (hasdiffthermal) then
  do i=1,meshdat%Nelocal
      do j=1,NEN
          nb = meshdat%ien(j,i)+1
          if (modeldat%lmtx(j,i).ne.0) then
              modeldat%lmtx(j,i) = &
      ISIGN(equationsmaskxt(nb),modeldat%lmtx(j,i))
          endif
      enddo
  enddo
   endif



!      do i=1,meshdat%Nelocal
!          write(FILE_outputf,*) "Equation numbers of local element ",
!     .      meshdat%locElt(i), " are : ",
!     >      ((modeldat%lm(k,j,i), k=1,2), j=1,NEN)
!      enddo

   ! here adapt nodal point numbers in AdjM and ien
   ! such that they are numbered locally



   do i=1,meshdat%Nvlocal
  do j=1,meshdat%itot(i)
      nb = meshdat%AdjM(j,i)
      meshdat%AdjM(j,i) = verticesmask(nb+1)! -1
  enddo
   enddo


!write(*,*) 'VM: ', verticesmask(1:100)
do i=1,meshdat%Nelocal
!    write(*,*) 'IEN before masking ', meshdat%ien(1:4,i)
    do j=1,NEN
        nb = meshdat%ien(j,i)
        meshdat%ien(j,i) = verticesmask(nb+1)! - 1
    enddo
!    write(*,*) 'IEN after masking ', meshdat%ien(1:4,i)

enddo



   !adapt linked nodes array
if (NLINK.gt.0) then
    if (iecho.eq.8) then
        call parallelLog("AOmod", "Linked node array")
    endif
!          write(FILE_outputf,*) 'linked nodes before petsc -> local'
!          do i=1,NLINK
!              write(FILE_outputf,*) i, (modeldat%LINK(j,i),j=2,3)
!          enddo
    do i=1,NLINK
        nb = modeldat%LINK(2,i)
        modeldat%LINK(2,i) = verticesmask(nb+1)
        nb = modeldat%LINK(3,i)
        modeldat%LINK(3,i) = verticesmask(nb+1)
    enddo
!          write(FILE_outputf,*) 'linked nodes after'
!          do i=1,NLINK
!              write(FILE_outputf,*) i, (modeldat%LINK(j,i),j=2,3)
!          enddo
endif

! adapt slippery nodes and elements
if (HasSlipperyNodes()) then
!    if (iecho.eq.8) then
!        call parallelLog("AOmod", "building NSLIP en IDSLP")
!    endif

    do i=1,NUMSLP

!        write(*,*) "rank", getrank(), "renumbery slippery entry i"

        nb = modeldat%NSLIP(2,i)

!        write(*,*) "rank", getrank(), "a"

        modeldat%NSLIP(2,i) = verticesmask(nb+1)

#ifdef SPARSE
        if (getsize().eq.1) then
            ! for single processor runs, local and global numbering are the same
            nb = modeldat%NSLIP(1,i)
        else
            nb = meshdat%Glo2LocElement(modeldat%NSLIP(1,i))
!            write(*,*) "rank", getrank(), "b", nb

        endif
#else
        nb = modeldat%NSLIP(1,i)
#endif

        !  this gets nb between 1 and meshdat%Nelocal
        if (nb.gt.0) then
            modeldat%NSLIP(1,i) = nb
        endif
!        write(*,*) "rank", getrank(), "d"

    enddo

!    write(*,*) "rank", getrank(), "finished slip loop"


    ! IDSLP contains the global vertex number.
    ! replace this by the local number,
    ! or by zero if the number does not occur in this partition
!    if (iecho.eq.8) then
!        call parallelLog("AOmod", "setting IDLP")
!    endif

!    write(*,*) "rank", getrank(), "has NUMSN", NUMSN

    do i=1,NUMSN
!        write(*,*) "rank", getrank(), "request IDSLP", i , " of", size(modeldat%IDSLP,1)
        nb = modeldat%IDSLP(i)
!        write(*,*) "rank", getrank(), "request vertmask ", i+1, "of", size(verticesmask,1)
        modeldat%IDSLP(i) = verticesmask(nb+1)
    enddo

!     write(*,*) "rank", getrank(), "starts renumbering elements"

    if (iecho.eq.8) then
        call parallelLog("AOmod", "renumbering elements for slippery nodes")
    endif
    call renumberelements(modeldat%IDSLE, NUMSE,"Slippery")
    if (iecho.eq.8) then
        call parallelLog("AOmod", "finished renumbering elements for slippery nodes")
    endif


!    write(*,*) "rank", getrank(), "Slipper has NSKSKEW: ", NSLSKEW

    ! with multiple processors, let numbers correspond to local numbers
    ! so that the references to equation numbers (which are on local nrs)
    ! can be applied later on.
    if (NSLSKEW.gt.0) then
        do i=1,NSLSKEW
#ifdef SPARSE
            if (getsize().eq.1) then
                ! for one processors, local and global numbering are the same
                nb = modeldat%NSELSD(1,i)
            else  
!                write(*,*) "rank", getrank(), "checks NSELSD ", i, "of", size(modeldat%NSELSD,2)
! PANIC, this enters -1's in NSELSD, which is incorrect.
!                nb = meshdat%Glo2LocElement(modeldat%NSELSD(1,i))
!                write(*,*) "rank", getrank(), "makes choice between", &
!                            meshdat%Glo2LocElement(modeldat%NSELSD(1,i)), "and" , &
!                            meshdat%FullGlo2LocElement(modeldat%NSELSD(1,i))
                 nb = meshdat%FullGlo2LocElement(modeldat%NSELSD(1,i))
            endif    
#else
            nb = modeldat%NSELSD(1,i)
#endif

!            write(*,*) "rank", getrank(), "found nb: ", nb

            modeldat%NSELSD(1,i) = nb
!            modeldat%NSELSD(1,i) = elementmask(nb+1)

        enddo

! dumb experiment
!        call renumberelements(modeldat%NSELSD, NSLSKEW,"FaultParElems")


    endif

!    if (iecho.eq.8) then
!        call parallelLog("AOmod", "completed building NSLIP en IDSLP")
!    endif


endif

!--------------------------------------------------------------------------------
#ifdef DEBUG
write(*,*) "rank",getrank(),"AOmodule made IDSLP: ", modeldat%IDSLP
write(*,*) 'doing faulted nodes? NUMFN=', NUMFN
#endif

! renumber faulted nodes and elements from global to local
if (NUMFN.gt.0) then
    if (iecho.eq.8) then
        call parallelLog("AOmod", "building NFAULT")
    endif
    do i=1,NUMFN
        nb = modeldat%NFAULT(2,i)                       ! global node number
        modeldat%NFAULT(2,i) = verticesmask(nb+1)       ! local node number

#ifdef DEBUG
        write(*,*) 'NFAULT(',i,') global node number ',nb,' set to local number ', verticesmask(nb+1)
#endif

#ifdef SPARSE
        ! get global element number
        if (getsize().eq.1) then
            nb = modeldat%NFAULT(1,i)
        else  
            nb = meshdat%Glo2LocElement(modeldat%NFAULT(1,i))
        endif
#else
        nb = modeldat%NFAULT(1,i)
#endif

#ifdef DEBUG
        write(*,*) 'NFAULT(',i,') global node number ',modeldat%NFAULT(1,i),&
         ' set to local number ',nb
#endif
        if (nb.gt.0) then
            modeldat%NFAULT(1,i) = nb
        endif
    enddo
endif

!--------------------------------------------------------------------------------
! adapt surface nodes
if (NSURF.gt.0) then
    if (iecho.eq.8) then
        call parallelLog("AOmod", "building surface nodes")
    endif

    do i=1,NSURF
        nb = modeldat%ISURF(i)
        modeldat%ISURF(i) = verticesmask(nb+1)
    enddo
endif


   ! renumber the element boundary contition from global numbering to
   ! local per processor numbering
   ! pre stresses
if (NPRE.gt.0) then
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumber elements for prestresses") 
    endif
    call renumberelements(modeldat%ISELM, NPRE,"Presstress")
!    call renumberelements(modeldat%ISELM, NPREglobal,"Presstress")
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumbered elements for prestresses")
    endif
endif

! pressure bc
if (NUMPR.gt.0) then
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumber elements for pressures")
    endif
    call renumberelements(modeldat%IELNO, NUMPR,"Pressure")
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumbered elements for pressures")
    endif
endif

! stress derived tractions
if (NUMSTR.gt.0) then
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumber elements for stresses")
    endif
    call renumberelements(modeldat%IELSTR, NUMSTR,"Stress")
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumbered elements for stresses")
    endif
endif

! winkler restoring pressures
if (NUMWNK.gt.0) then
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumber elements for Winkler pressures")
    endif
    call renumberelements(modeldat%IWELM, NUMWNK,"Winkler restoring")
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumbered elements for Winkler pressures")
    endif
endif

#ifdef SHEARZONE_TRACTIONS

if (NUMELSHEARZONE.gt.0) then
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumber elements for shear zone tractions")
    endif
    call renumberelements(modeldat%ELNOSHEARZONE, NUMELSHEARZONE,"shear zone tractions")
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumbered elements for shear zone tractions")
    endif
endif

#endif


! thermal anomaly
if (NTANOM.gt.0) then
    do i=1,NTANOM
        nb = modeldat%ITANOM(1,i)
        modeldat%ITANOM(1,i) = verticesmask(nb+1)
    enddo
endif

! heat flux
if (NFLX.gt.0) then
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumber elements for heat flux")
    endif
    call renumberelements(modeldat%IFLX, NFLX,"heat flux")
    if (iecho.eq.8) then
        call parallelLog("renumberEquations", "renumbered elements for heat flux")
    endif
endif

! thermal winkler flux
if (NTWINK.gt.0) then
    do i=1,NTWINK

#ifdef SPARSE
        if (getsize().eq.1) then
            nb = modeldat%ITWINK(1,i)
        else  
            nb = meshdat%Glo2LocElement(modeldat%ITWINK(1,i))
        endif
#else
        nb = modeldat%ITWINK(1,i)
#endif

        if (nb.gt.0) then
            modeldat%ITWINK(1,i) = nb
        endif
    enddo
endif


!      call writeAOmeshdata(FILE_outputf, myrank, meshdat, tmp, nen, 2)
!      if (parout) call writeAOdatasummary(FILE_outputf, getrank(), &
!           meshdat,modeldat,nen,ndof)

if (iecho.eq.8) then
    call parallelLog("AOmod", "Completed equation numbering")
endif

!if (allocated(verticesmask)) then
!    deallocate(verticesmask)
!endif




return
end subroutine createEquationNumbering

!***************************************************************


!***************************************************************
subroutine CreateVerticesMask(meshdat,ndof,ien,NEN,AdjM,maxnbrs, caller)
!***************************************************************
!  ien and adjm are not taken from meshdata structure
!  because the call to this sub has been isolated from
!  createParallelNumbering en modifications to 
!  meshdatactx%ien and ...%AdjM required them to be 
!  separate arguments for flexibility
!  This is admittedly not the most elegant solution...

USE MESHDATAMODULE
use modelctx,       only: getrank, &
                          getsize
use debugmodule,    only: parallelLog, &
                          iecho, &
                          allocateError, &
                          delay
!***************************************************************
implicit none
!***************************************************************
type (meshdata) :: meshdat
integer i,j,nb,k, ERROR, islink,ndof
integer :: NEN, maxnbrs
integer :: caller ! for debugging

integer :: ien(NEN, meshdat%Neglobal)
integer :: adjM(maxnbrs, meshdat%Nvlocal)

!***************************************************************
if (allocated(vertices)) then
!    deallocate(vertices)
    write(*,*) 'vertices already allocated, something odd'
    ! do nothing
else
    ALLOCATE(vertices(meshdat%Nvglobal), STAT=ERROR)
!    ALLOCATE(vertices(meshdat%Nvlocal), STAT=ERROR)
    call AllocateError("vertices", ERROR)
endif

if (allocated(verticesmask)) then
    write(*,*) 'verticesmask already allocated, something odd'

!    deallocate(verticesmask)
    ! do nothing
else
!    ALLOCATE(verticesmask(meshdat%Nvlocal), STAT=ERROR)
    ALLOCATE(verticesmask(meshdat%Nvglobal), STAT=ERROR)
    call AllocateError("verticesmask", ERROR)
endif

if (iecho.eq.8) then 
!    call ParallelLog("CreatevertMask", "at beginning")
endif

!   write(*,*) 'preVM rank', getrank(), 'at beginning of verticesmask'
!   write(*,*) 'preVM rank', getrank(), 'meshdat%locInd ', meshdat%locInd
!   write(*,*) 'preVM rank', getrank(), 'meshdat%itot ', meshdat%itot
!   write(*,*) 'preVM rank', getrank(), 'meshdat%AdjM ', AdjM
!   write(*,*) 'preVM rank', getrank(), 'meshdat%ien ', ien

do i=1,meshdat%Nvglobal
    verticesmask(i) = 0
    vertices(i) = 0
enddo
nvertices = 0

meshdat%NGhostPoints = 0
meshdat%NvRest = 0

   ! first add vertices that are owned
do i=1,meshdat%Nvlocal
!          write(*,*) 'rank ', getrank(), 'says: the ',i,'th point is: ' , meshdat%locInd(i)+1
    vertices(nvertices+1) = meshdat%locInd(i)+1 ! +1 to go from 0-based to 1-based
!    write(*,*) 'set ', vertices(nvertices+1), 'to', meshdat%locInd(i)+1
    nvertices = nvertices + 1
    verticesmask(meshdat%locInd(i)+1) = nvertices ! +1 samo samo
enddo

!write(*,*) 'rank ', getrank(), 'VM after own : ', verticesmask
!write(*,*) 'rank ', getrank(), "adjM says; itot: ", meshdat%itot
!write(*,*) 'rank ', getrank(), 'Adjancency: ', adjM
!write(*,*) 'rank ', getrank(), 'mesh Adjancency: ', meshdat%AdjM


   ! add neighbor vertices
do i=1,meshdat%Nvlocal

!    write(*,*) "adjM says; AdjM: ", adjM(:,i)


    do j=1,meshdat%itot(i) ! loop over the neighbors

!        nb = meshdat%AdjM(j,i) ! nb = *local* number of neighbor
        nb = adjM(j,i)

!        write(*,*) "read neighbor index: ", i, j, nb


!        nb = adjM(j,i)
        if (nb.lt.0) then
            write(*,*) "Rank", getrank(),&
                       "says: Neighbour",j,&
                       "of local point", i,&
                       "is 0, while the point should have", meshdat%itot(i), &
                       "neighbors."
            write(*,*) "Most likely cause: Mesh incorrect."
            stop "Leaving GTecton..."
        endif

        if (nb.ge.size(verticesmask,1)) then
            write(*,*) "Rank", getrank(),&
                       "says: Neighbour",j,&
                       "of local point", i,&
                       "is ",nb,", while the mesh only has ", meshdat%nvglobal, &
                       "points."
            write(*,*) "Most likely cause: Mesh incorrect."
            stop "Leaving GTecton..."
        endif


!              write(*,*) 'rank ', getrank(), 'says: neighbor ',j,' of ',i, ' is ' ,nb+1

        if (verticesmask(nb+1).eq.0) then

!                  write(*,*) 'rank ', getrank(), 'unknown node!'
            vertices(nvertices+1) = nb+1
!            write(*,*) 'set ', vertices(nvertices+1), 'to', nb+1
            nvertices = nvertices + 1
            meshdat%Nghostpoints = meshdat%NGhostPoints + 1
            verticesmask(nb+1) = nvertices
        endif
    enddo
enddo


!write(*,*) 'rank ', getrank(), 'VM after neighbors: ', verticesmask


   ! add vertices of elements
   ! this is required since an element may belong to processor i but
   ! its vertices to other processors, in which case these vertices
   ! are not necessarily neighbors of a vertex owned by processor i
   ! (although in practise, this case will not occur frequently)
do i=1,meshdat%Nelocal
    do j=1,4!NEN !fix this
!      nb = meshdat%ien(j,i)
        nb = ien(j,i)
!              write(*,*) 'rank ', getrank(), 'says: node ',j,' of elem',i, ' is ' ,nb+1
        if (verticesmask(nb+1).eq.0) then
!                  write(*,*) 'rank ', getrank(), 'unknown node!'
            vertices(nvertices+1) = nb+1
!            write(*,*) 'set ', vertices(nvertices+1), 'to', nb+1
            nvertices = nvertices + 1
            meshdat%NvRest = meshdat%NvRest + 1
            verticesmask(nb+1) = nvertices
        endif
    enddo
enddo


!write(*,*) 'rank ', getrank(), 'vertices after all: ', vertices


!write(*,*) 'rank ', getrank(), 'VM after all: ', verticesmask


meshdat%NvAllRelevant = meshdat%Nvlocal + meshdat%Nghostpoints + meshdat%NvRest
!   write(*,*) 'rank', getrank(), ' has', meshdat%NGhostPoints, ' ghost points'
!   write(*,*) 'rank', getrank(), ' has', meshdat%NvRest, ' rest points'
!   write(*,*) 'rank', getrank(), ' has', meshdat%NvAllRelevant, ' total relevant vertices.'


! do i=0,getsize()
!call delay()
!if (i.eq.getrank()) then
!   write(*,*) 'caller ', caller, 'rank', getrank(), 'vertices     ', vertices
!endif
!enddo

!do i=0,getsize()
!call delay()
!if (i.eq.getrank()) then
!   write(*,*) 'caller ', caller, 'rank', getrank(), 'verticesmask ', verticesmask
!endif
!enddo

end subroutine CreateVerticesMask
!*****************************************************************

!***************************************************************
subroutine CreateEquationsMask(meshdat, modeldat,ndof)
!***************************************************************
#ifdef EXPERIMENTAL_ReducedSlipAlloc
USE MODELDATAMODULE, only: modeldata, &
                           equationreference, &
                           HasSlipperyNodes, &
                           elementSlipperyID
USE MODELDEFINITION, only: NUMSLP, &
                           NEQ, &
                           NUMSLPglobal, &
                           nElemsWithSlipperyNodes
#else
USE MODELDATAMODULE, only: modeldata, &
                           equationreference, &
                           HasSlipperyNodes
USE MODELDEFINITION, only: NUMSLP, &
                           NEQ, &
                           NUMSLPglobal
#endif
USE MESHDATAMODULE,  only: meshdata
USE MODELCTX
use debugmodule,     only: parallellog, &
                           debug, &
                           iecho
use modeltopology,   only: nen
use algebra,         only: arraycontains
!***************************************************************
implicit none
!***************************************************************
!     this subroutine creates, for each processor, an array 
!     containing the indices corresponding to the equation numbers 
!     of all of the locally needed vertices.
!     
!     Those are not only the vertices that have been assigned to
!     the partition of this processor, but also:
!
!     1: vertices of other partitions but belonging to any element in the partition
!     2: neighboring vertices belonging to all vertices in the partition
!     3: neighboring vertices of the vertices in (1).  (not sure about this)
!   
!     There are 4 equationsmasks used.
!     Two for the thermal problem, and two for the mechanical.
!     - One of each of those for the normal nodal equations,
!     - The other is for slippery nodes. This gives extra degrees
!       of freedom which would be memory-inefficient to administrate 
!       in the 'normal' mask.
!
!     They are:   
!     - equationsmask     (normal mechanical equations)
!     - equationsmaskx    (eXtra slippery nodes equations)
!     - equationsmaskt    (Thermal equation)
!     - equationsmaskxt   (equations for eXtra differential Thermal thingies)



type (meshdata) :: meshdat
type (modeldata) :: modeldat

integer ndof,isdof,i,j,k,nb, ERROR, islink, lind, linkto, cand
integer next
logical doslip

integer :: sizeIDX1
integer :: sizeIDX2
integer :: sizeadjIDX1
integer :: sizeadjIDX2
integer :: sizeadjIDX3 
integer :: sizeeqmaskx1
integer :: sizeeqmaskx2
integer :: sizeLMX1
integer :: sizeLMX2
integer :: sizeLMX3
integer :: sizeequations

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer :: globalElemID
integer :: slipperySequenceNr
#endif

!***************************************************************
!     reuse IDglobal and IBONDglobal
!      DEALLOCATE( modeldat%IDglobal )

!write(*,*) 'about to allocate equations to: ', NDOF*meshdat%Nvglobal
!write(*,*) 'with NUMSLP(glob): ',NUMSLP, NUMSLPglobal

!TODO... not good yet, should make room for extra equations,
!not for every slippery node entry, because that number is higher.
!In TECIN, the nodes must be counted somehow... Is this worth it?

! the mask is split up between normal equations and slippery equations.
! However, the equations array holds both.

! The check for the slippery nodes is always executed before the 
! the administration of the normal nodes. As such, the slippery equations
! of a vertex are always *before* the normal equations.


ALLOCATE(equations(NDOF*meshdat%Nvglobal+NUMSLP), STAT=ERROR)
ALLOCATE(equationsmask(NDOF,meshdat%Nvglobal), STAT=ERROR)

equations = 0
equationsmask = 0
nequations = 0

doslip = HasSlipperyNodes()

if (doslip) then
    ALLOCATE(equationsmaskx(NDOF,meshdat%Nvglobal), STAT=ERROR)
    equationsmaskx = 0
endif

sizeequations = size(equations,1)

sizeIDX1    = size(modeldat%IDX,1)
sizeIDX2    = size(modeldat%IDX,2)  
sizeadjIDX1 = size(modeldat%adjIDX,1)
sizeadjIDX2 = size(modeldat%adjIDX,2)
sizeadjIDX3 = size(modeldat%adjIDX,3) 

sizeeqmaskx1 =  size(equationsmaskx,1)
sizeeqmaskx2 =  size(equationsmaskx,2)

sizeLMX1 = size(modeldat%lmx, 1)
sizeLMX2 = size(modeldat%lmx, 2)
sizeLMX3 = size(modeldat%lmx, 3)

!if (debug) then
!write(*,*) 'rank', getrank(), '#### equations: ', equations

!    do i=1,meshdat%Nelocal
!    write(*,*) 'rank', getrank(), ' #### ',i,'LM: ', modeldat%LM(:,:,i)
!    enddo
!endif

if (iecho.eq.8) then
    call parallellog("createequationsmask", "masking own equations")
endif


! equations vertices in this partition.
do i=1,meshdat%Nvlocal
    do j=1,NDOF

        if (doslip) then
            if (modeldat%IDX(j,i).ne.0) then
                ! add an extra degree of freedom for the slippery node.
                equations(nequations+1) = modeldat%IDX(j,i)
                nequations = nequations + 1
                equationsmaskx(j,meshdat%locInd(i)+1) = nequations
            endif
        endif

        if (modeldat%ID(j,i).ne.0) then
            call arraycontains(equations, nequations, modeldat%ID(j,i), next)
            if (next.le.nequations) then
                equationsmask(j,meshdat%locInd(i)+1) = next
            else
                equations(nequations+1) = modeldat%ID(j,i)
                nequations = nequations + 1
                equationsmask(j,meshdat%locInd(i)+1) =  nequations
            endif
        endif
    enddo
enddo

!if (debug) then
!write(*,*) 'rank', getrank(), 'own equations done: ' , equations 
!endif

if (iecho.eq.8) then
    call parallellog("createequationsmask", "masking neighboring equations")
endif


! equations of the neighbours of the vertices of this partition
! this is done in the sequence of the adjacency array, 
! which is not necessarily the same as the sequence of the vertex IDs
do i=1,meshdat%Nvlocal
!          write(*,*) 'checking neighbours of node ', i
    do j=1,meshdat%itot(i)         ! loop over neighbors
        nb = meshdat%AdjM(j,i)+1   ! neighbor index
!              write(*,*) 'checking neighbour', nb
        do k=1,NDOF

            if (doslip) then
                if (equationsmaskx(k,nb).eq.0) then
                    if (modeldat%AdjIDX(k,j,i).ne.0) then

                        call arraycontains(equations, nequations, modeldat%AdjIDX(k,j,i), next)

                        if (next.le.nequations) then
                            equationsmaskx(k,nb) = next
                        else
                            equations(nequations +1) = modeldat%AdjIDX(k,j,i)
                            nequations = nequations + 1
                            equationsmaskx(k,nb) = nequations
                        endif

                    endif
                endif
            endif

            if (equationsmask(k,nb).eq.0) then
                if (modeldat%AdjID(k,j,i).ne.0) then
                    call arraycontains(equations, nequations, modeldat%AdjID(k,j,i), next)
                    if (next.le.nequations) then
                        equationsmask(k,nb) = next
                    else
                        equations(nequations +1) = modeldat%AdjID(k,j,i)
                        nequations = nequations + 1
                        equationsmask(k,nb) = nequations
                    endif
                endif
            endif

        enddo
    enddo
enddo

!if (debug) then
!write(*,*) 'rank', getrank(), 'neighbor point equations done: ' , equations 
!endif

!      write(*,*) 'after adding neighbors of owned vertices'
!      write(*,*) 'equationsmask', equationsmask
!      write(*,*) 'equationsmaskx', equationsmaskx


!      do i=1,meshdat%Nvlocal
!      if (doslip) write(FILE_outputf,*) i,(modeldat%IDX(j,i),j=1,NDOF)
!      enddo

! add equations of elements
! this is required since an element may belong to processor is but
! its vertices to other processors, in which case these vertices
! are not necessarily neighbors of a vertex owned by processor i


if (iecho.eq.8) then
    call parallellog("createequationsmask", "masking remote equations")
endif


if (debug) then
!    do i=1,meshdat%Nelocal
!        write(*,*) 'rank', getrank(), ' ####',i,' LMX: ', modeldat%LMX(:,:,i)
!    enddo
endif


#ifdef EXPERIMENTAL_ReducedSlipAlloc
!do i=1,nElemsWithSlipperyNodes
do i=1,meshdat%Nelocal
#else
do i=1,meshdat%Nelocal
#endif

    do j=1,NEN
        nb = meshdat%ien(j,i)+1

        do k=1,NDOF


            if (doslip) then
                ! do the slippery equations first!
                if (equationsmaskx(k,nb).eq.0) then

#ifdef EXPERIMENTAL_ReducedSlipAlloc
!                     write(*,*) "rank", getrank(),"a gloelt", meshdat%gloElt


!                    globalElemID = meshdat%gloElt(i)+1 ! gloElt starts counting at 0, hence +1
!                    write(*,*) "checking global elem nr of ", i, "which is", globalElemID

                    slipperySequenceNr = elementSlipperyID(i,3)
!                    write(*,*) "linked to slip sequence nr", slipperySequenceNr

                    if (slipperySequenceNr.gt.0) then

                        if (modeldat%lmx(k,j,slipperySequenceNr).ne.0) then
                            call arraycontains(equations, &
                                             nequations, &
                                             ABS(modeldat%lmx(k,j,slipperySequenceNr)), &
                                             next)

                            if (next.le.nequations) then
                                ! ABS(modeldat%lmx(k,j,i)) occurs in the equations array
                                equationsmaskx(k,nb) = next
                            else

                                equations(nequations+1) = ABS(modeldat%lmx(k,j,slipperySequenceNr))
                                nequations = nequations + 1
                                equationsmaskx(k,nb) = nequations
                            endif

                        endif
                    endif

#else
          
                    if (modeldat%lmx(k,j,i).ne.0) then
                        call arraycontains(equations, &
                                         nequations, &
                                         ABS(modeldat%lmx(k,j,i)), &
                                         next)

                        if (next.le.nequations) then
                            ! ABS(modeldat%lmx(k,j,i)) occurs in the equations array
                            equationsmaskx(k,nb) = next
                        else

                            equations(nequations+1) = ABS(modeldat%lmx(k,j,i))

                            nequations = nequations + 1
                            equationsmaskx(k,nb) = nequations

                        endif                 
                    endif
#endif
                endif
            endif

!write(*,*) "rank", getrank(),"now do normal"


            ! and do the normal ones after second
            if (equationsmask(k,nb).eq.0) then

                if (modeldat%lm(k,j,i).ne.0) then
                    call arraycontains(equations, nequations, &
                                       modeldat%lm(k,j,i), next)
                    if (next.le.nequations) then
                        equationsmask(k,nb) = next
                    else
                        equations(nequations+1) = &
                        modeldat%lm(k,j,i)
                        nequations = nequations + 1
                        equationsmask(k,nb) = nequations
                    endif
                endif
            endif

!write(*,*) "rank", getrank(),"normal done"


        enddo
    enddo
enddo


!write(*,*) 'rank', getrank(), 'added equations of own elements: ' , equations


if (iecho.eq.8) then
    call parallellog("createequationsmask", "masking done")
endif


if (.not. hasthermal) then
    goto 111
endif

!*******************************************************************
! possibly do the thermal equation as well.

ALLOCATE(equationst(meshdat%Nvglobal), STAT=ERROR)
ALLOCATE(equationsmaskt(meshdat%Nvglobal), STAT=ERROR)
equationst = 0
equationsmaskt = 0
nequationst = 0

if (hasdiffthermal) then
    ALLOCATE(equationsmaskxt(meshdat%Nvglobal), STAT=ERROR)
    equationsmaskxt = 0
endif

do i=1,meshdat%Nvlocal
    if (hasdiffthermal) then
        if (modeldat%IDTX(i).ne.0) then
            equationst(nequationst+1) = modeldat%IDTX(i)
            nequationst = nequationst + 1
            equationsmaskxt(meshdat%locInd(i)+1) = &
                 nequationst
        endif
    endif
    if (modeldat%IDT(i).ne.0) then
        equationst(nequationst+1) = modeldat%IDT(i)
        nequationst = nequationst + 1
        equationsmaskt(meshdat%locInd(i)+1) = nequationst
    endif
enddo

do i=1,meshdat%Nvlocal
    do j=1,meshdat%itot(i)
        nb = meshdat%AdjM(j,i)+1
        if (hasdiffthermal) then
            if (equationsmaskxt(nb).eq.0) then
                if (modeldat%AdjIDTX(j,i).ne.0) then
                    equationst(nequationst +1) = &
                    modeldat%AdjIDTX(j,i)
                    nequationst = nequationst + 1
                    equationsmaskxt(nb) = nequationst
                endif
            endif
        endif
        if (equationsmaskt(nb).eq.0) then
            if (modeldat%AdjIDT(j,i).ne.0) then
                equationst(nequationst +1) = &
                modeldat%AdjIDT(j,i)
                nequationst = nequationst + 1
                equationsmaskt(nb) = nequationst
            endif
        endif
    enddo
enddo


! add equations of elements
! this is required since an element may belong to processor i but
! its vertices to other processors, in which case these vertices
! are not necessary neighbors of a vertex owned by processor i
do i=1,meshdat%Nelocal
    do j=1,4!NEN !fix this
        nb = meshdat%ien(j,i)+1
        if (hasdiffthermal) then
            if (equationsmaskxt(nb).eq.0) then
                if (modeldat%lmtx(j,i).ne.0) then
                    call arraycontains(equationst, nequationst, &
                                       ABS(modeldat%lmtx(j,i)), next)
                    if (next.le.nequationst) then
                        equationsmaskxt(nb) = next
                    else
                        equationst(nequationst+1) = &
                        ABS(modeldat%lmtx(j,i))
                        nequationst = nequationst + 1
                        equationsmaskxt(nb) = nequationst
                    endif
                endif
            endif
        endif
        if (equationsmaskt(nb).eq.0) then
            if (modeldat%lmt(j,i).ne.0) then
                call arraycontains(equationst, nequationst, &
                                   modeldat%lmt(j,i), next)
                if (next.le.nequationst) then
                    equationsmaskt(nb) = next
                else
                    equationst(nequationst+1) = &
                    modeldat%lmt(j,i)
                    nequationst = nequationst + 1
                    equationsmaskt(nb) = nequationst
                endif
            endif
        endif
    enddo
enddo


111 allocate(equationreference(NEQ))

do i=1,NEQ
    equationreference(i) = equations(i)
enddo


if (iecho.eq.8) then
    call parallellog("createequationsmask", "masking sub done")
endif

end subroutine CreateEquationsMask
!***************************************************************

!***************************************************************
subroutine checkreadwritearray(message,rank,max,idx)
! small subroutine to validate whether a read or write goes out of bounds
! Because it can happen without a proper error.
implicit none

integer :: rank, max, idx
character(len=1) :: message

if ((idx.gt.max) .or. (idx.lt.1)) then
    write(*,*) 'readwrite attempt FAIL ', message, '  idx:  ', idx, '  to max  ', max, '  on rank: ', rank
else
!          write(*,*) 'readwrite attempt OK   ', message, '  idx:  ', idx, '  to max  ', max, '  on rank: ', rank
endif


end subroutine



!***************************************************************




!***************************************************************
   subroutine writeAOdatasummary(FILE_outputfd, myrank, meshdat,modeldat, &
   nen,ndof)
!***************************************************************
   USE MODELDATAMODULE
   USE MODELDEFINITION
   USE MESHDATAMODULE
   USE MODELCTX
!***************************************************************
   implicit none
!***************************************************************
   type (meshdata) :: meshdat
   type (modeldata) :: modeldat

   integer i, j, k
   integer myrank, FILE_outputfd, nen,ndof,shift
   integer eqnsdof(ndof), eqnsth
!***************************************************************

   write(FILE_outputf, 31) meshdat%Nvlocal, NEQlocal, myrank
31 format(1x, 'Mapping of ', I5, ' vertices and ', &
  I5, ' mechanical equations owned by rank ', I5)

   if (HasSlipperyNodes()) then
  write(FILE_outputf,*) '   In case of slippery nodes for each vertex', &
     ' first row: slippery equations; ', &
     ' second row: normal equations '
   endif

   write(FILE_outputf,*) '     global ( eqn numb )', &
             '        petsc ( eqn numb )', &
             '        local ( eqn numb )'
   shift = 0

   if (HasSlipperyNodes()) then
  do i=1,meshdat%Nvlocal
      do k=1,NDOF
          eqnsdof(k) = 0
          if (equationsmaskx(k,vertices(i)).ne.0) then
              shift = shift + 1
              eqnsdof(k) = equations(shift)
          endif
      enddo
      write(FILE_outputf, 32) meshdat%gloInd(i), &
          (modeldat%IDXglobal(j,meshdat%gloInd(i)),j=1,NDOF), &
         vertices(i), (eqnsdof(j),j=1,NDOF), &
        i, (equationsmaskx(j,vertices(i)),j=1,NDOF)
      do k=1,NDOF
          eqnsdof(k) = 0
          if (equationsmask(k,vertices(i)).ne.0) then
              shift = shift + 1
              eqnsdof(k) = equations(shift)
          endif
      enddo
      write(FILE_outputf, 34) &
       (modeldat%IDglobal(j,meshdat%gloInd(i)),j=1,NDOF), &
!     .         (modeldat%ID(j,i),j=1,NDOF), &
       (eqnsdof(j),j=1,NDOF), &
       (equationsmask(j,vertices(i)),j=1,NDOF)
  enddo
   else
  do i=1,meshdat%Nvlocal
      do k=1,NDOF
          eqnsdof(k) = 0
          if (equationsmask(k,vertices(i)).ne.0) then
              shift = shift + 1
              eqnsdof(k) = equations(shift)
          endif
      enddo
      write(FILE_outputf, 32) meshdat%gloInd(i), &
       (modeldat%IDglobal(j,meshdat%gloInd(i)),j=1,NDOF),&
!     .         vertices(i), (modeldat%ID(j,i),j=1,NDOF), &
       vertices(i), (eqnsdof(j),j=1,NDOF), &
       i, (equationsmask(j,vertices(i)),j=1,NDOF)
  enddo
   endif

   write(FILE_outputf,*) '   vertices and mechanical equations ', &
                 'owned by neighbors'
   if (NLINK.gt.0) then
  write(FILE_outputf,*) '                        ', &
             '        petsc ( eqn numb*)', &
             '        local ( eqn numb )'
   else
  write(FILE_outputf,*) '                        ', &
             '        petsc ( eqn numb )', &
             '        local ( eqn numb )'
   endif

   if (nvertices.gt.meshdat%Nvlocal) then
  shift = NEQlocal
  if (HasSlipperyNodes()) then
      do i=meshdat%Nvlocal+1,nvertices
          do k=1,NDOF
              eqnsdof(k) = 0
              if (equationsmaskx(k,vertices(i)).ne.0) then
                  shift = shift + 1
                  eqnsdof(k) = equations(shift)
              endif
          enddo
          write(FILE_outputf, 33) &
           vertices(i), (eqnsdof(j),j=1,NDOF), &
           i, (equationsmaskx(j,vertices(i)),j=1,NDOF)
          do k=1,NDOF
              eqnsdof(k) = 0
              if (equationsmask(k,vertices(i)).ne.0) then
                  shift = shift + 1
                  eqnsdof(k) = equations(shift)
              endif
          enddo
          write(FILE_outputf, 35) &
           (eqnsdof(j),j=1,NDOF), &
           (equationsmask(j,vertices(i)),j=1,NDOF)

      enddo
  else
      do i=meshdat%Nvlocal+1,nvertices
          do k=1,NDOF
              eqnsdof(k) = 0
              if (equationsmask(k,vertices(i)).ne.0) then
                  shift = shift + 1
                  eqnsdof(k) = equations(shift)
              endif
          enddo
          write(FILE_outputf, 33) &
           vertices(i), (eqnsdof(j),j=1,NDOF), &
           i, (equationsmask(j,vertices(i)),j=1,NDOF)
      enddo
  endif
   endif
   if (NLINK.gt.0) then
  write(FILE_outputf,*) ' *) petsc equation numbers not correct ', &
      'in case of linked nodes'
   endif


   if (.not. hasthermal) goto 121

   write(FILE_outputf, 41) NTEQlocal, myrank
41 format(1x, 'and mapping of ', I5, &
  ' thermal equations owned by rank ', I5)
   if (hasdiffthermal) then
  write(FILE_outputf,*) '   In case of differential ', &
 ' temperatures switch first row: differential equation; ', &
 ' second row: normal equation '
   endif

   write(FILE_outputf,*) '     global ( eqn )', &
             '        petsc ( eqn )', &
             '        local ( eqn )'
   shift = 0
   if (hasdiffthermal) then
  do i=1,meshdat%Nvlocal
      if (equationsmaskxt(vertices(i)).ne.0) then
          shift = shift + 1
          eqnsth =  equationst(shift)
      endif
      write(FILE_outputf, 42) meshdat%gloInd(i), &
          modeldat%IDTXglobal(meshdat%gloInd(i)),&
!     .           vertices(i), modeldat%IDTX(i), &
         vertices(i), eqnsth, &
        i, equationsmaskxt(vertices(i))
      if (equationsmaskt(vertices(i)).ne.0) then
          shift = shift + 1
          eqnsth =  equationst(shift)
      endif
      write(FILE_outputf, 44) &
       modeldat%IDTglobal(meshdat%gloInd(i)),&
!     .         modeldat%IDT(i), &
       eqnsth, &
       equationsmaskt(vertices(i))
  enddo
   else
  do i=1,meshdat%Nvlocal
      if (equationsmaskt(vertices(i)).ne.0) then
          shift = shift + 1
          eqnsth =  equationst(shift)
      endif
      write(FILE_outputf, 42) meshdat%gloInd(i), &
       modeldat%IDTglobal(meshdat%gloInd(i)),&
!     .         vertices(i), modeldat%IDT(i), &
       vertices(i), eqnsth, &
       i, equationsmaskt(vertices(i))
  enddo
   endif

   write(FILE_outputf,*) '   vertices and thermal equations ', &
                 'owned by neighbors'
   write(FILE_outputf,*) '                   ', &
             '        petsc ( eqn )', &
             '        local ( eqn )'
   if (nvertices.gt.meshdat%Nvlocal) then
  shift = NTEQlocal
  if (hasdiffthermal) then
      do i=meshdat%Nvlocal+1,nvertices
          eqnsth = 0

          if (equationsmaskxt(vertices(i)).ne.0) then
              shift = shift + 1
              eqnsth = equationst(shift)
          endif

          write(FILE_outputf, 43) vertices(i), eqnsth, i, (equationsmaskx(j,vertices(i)),j=1,NDOF)
          eqnsth = 0

          if (equationsmaskt(vertices(i)).ne.0) then
              shift = shift + 1
              eqnsth = equationst(shift)
          endif

          write(FILE_outputf, 45) eqnsth, equationsmaskt(vertices(i))

      enddo
  else
      do i=meshdat%Nvlocal+1,nvertices
          eqnsth = 0

          if (equationsmaskt(vertices(i)).ne.0) then
              shift = shift + 1
              eqnsth = equationst(shift)
          endif

          write(FILE_outputf, 43) vertices(i), eqnsth, i, equationsmaskt(vertices(i))
      enddo
  endif
   endif



32 format(3x, I8, '  (',2I5,') || ', I8, '  (',2I5,') || ', &
      I8, '  (',2I5,')' )
34 format(11x, '  (',2I5,') || ', 8x, '  (',2I5,') || ', &
      8x, '  (',2I5,')' )
33 format(29x, I8, '  (',2I5,') || ', &
      I8, '  (',2I5,')' )
35 format(37x, '  (',2I5,') || ', &
      8x, '  (',2I5,')' )
42 format(3x, I8, '  (',I5,') || ', I8, '  (',I5,') || ', &
      I8, '  (',I5,')' )
44 format(11x, '  (',I5,') || ', 8x, '  (',I5,') || ', &
      8x, '  (',I5,')' )
43 format(24x, I8, '  (',I5,') || ', &
      I8, '  (',I5,')' )
45 format(32x, '  (',I5,') || ', &
      8x, '  (',I5,')' )

 121  write(FILE_outputf,* ) '  local neighbors of local vertices'
   do i=1,meshdat%Nvlocal
  write(FILE_outputf, 14) i, (meshdat%AdjM(j,i),j=1,meshdat%itot(i))
   enddo
14 format(3x, I5, ' | ', 15I8)

   write(FILE_outputf, 21) meshdat%Nelocal, myrank
21 format(1x,'Mapping of ',I5, ' elements owned by rank ', I5)
   write(FILE_outputf,*) '   global      petsc      '
   do i=1,meshdat%Nelocal
  write(FILE_outputf, 22) meshdat%gloElt(i), meshdat%locElt(i)
   enddo
22 format(3x, I8, ' || ', I8)
   write(FILE_outputf, *) 'Local vertices for local elements'
   do i=1,meshdat%Nelocal
  write(FILE_outputf, 23) i, (meshdat%ien(j,i),j=1,nen)
   enddo
23 format(3x, I5, ' | ', 10I12)
   return
   end subroutine writeAOdatasummary
!***************************************************************



!***************************************************************
subroutine scatterdata(meshdat, modeldat, myrank, FILE_outputf)
!***************************************************************
! this routine encapsulates the scattering of data to neighboring
! nodes, if needed. this is needed for data like coordinates of
! vertices, boundary conditions, slippery nodes, etc.
!***************************************************************
USE MESHDATAMODULE
USE MODELDATAMODULE
USE MODELTOPOLOGY
USE MODELDEFINITION
use modelctx,          only: getrank ! for debugging
use debugmodule,       only: iecho, &
                             parallelLog, &
                             delay, longdelay, &
                             verified_mpi_barrier, &
                             allocateError
#ifdef SPARSE
use sPETScmodule,       only: PETSC_error
#endif
!***************************************************************
implicit none
!***************************************************************
type (meshdata)  :: meshdat
type (modeldata) :: modeldat
integer          :: myrank, ERROR, hier
integer, intent(OUT) ::  FILE_outputf
!***************************************************************

call mpi_barrier(MPI_COMM_WORLD, error)

ndofVecSet   = .false.
nsdVecSet    = .false.
oneVecSet    = .false.
twoVecSet    = .false.


! do rotations
if (iecho.eq.8) then
    call parallellog("scatterdata", "scattering rotations")
endif

call scatterTWOdataInit(meshdat,modeldat,FILE_outputf)

!call verified_mpi_barrier("Finished scatter 2 init")

call scatterTWOdata(meshdat,modeldat, modeldat%SKEW,FILE_outputf,-1)
!call scatterNDOFdata(meshdat,modeldat, getrank(), modeldat%SKEW,FILE_outputf,-1,0)


! we do not want to start scattering when not all the processors
! have finished allocating modeldat%SKEW
!call mpi_barrier(MPI_COMM_WORLD, error)

!if (iecho.eq.8) then
!    call parallellog("scatterTWOdata", "done first call")
!endif

!call verified_mpi_barrier("Deallocing")


if (allocated(modeldat%SKEW)) then
    ! should be alloced, but just for robustness sake
    DEALLOCATE(modeldat%SKEW)
endif

call mpi_barrier(MPI_COMM_WORLD, error)


!call verified_mpi_barrier("Dealloced")



ALLOCATE(modeldat%SKEW(2,nvertices), stat =ERROR)
call allocateError("modeldata SKEW", ERROR)

!if (iecho.eq.8) then
!    call parallellog("scatterTWOdata", "realloced")
!endif


call mpi_barrier(MPI_COMM_WORLD, error)

modeldat%SKEW = 0d0

!        write(*,*) "rank", getrank(), " nulled"

! we do not want to start scattering when not all the processors
! have finished allocating modeldat%SKEW
!call mpi_barrier(MPI_COMM_WORLD, error)

!write(*,*) "rank", getrank(), " another two"


call scatterTWOdata(meshdat,modeldat, modeldat%SKEW,FILE_outputf,1)
!call scatterNDOFdata(meshdat,modeldat, getrank(), modeldat%SKEW,FILE_outputf,1, 0)




call scatterNSDdataInit(meshdat,modeldat,myrank,FILE_outputf)
call scatterNSDdata(meshdat,modeldat,myrank, meshdat%X,FILE_outputf,-1)

DEALLOCATE(meshdat%X,STAT=ERROR)
ALLOCATE(meshdat%X(nsd,nvertices),STAT=ERROR)
meshdat%X = 0d0

! we do not want to start scattering when not all the processors
! have finished allocating meshdat%X
call mpi_barrier(MPI_COMM_WORLD, error) 

call scatterNSDdata(meshdat,modeldat,myrank, meshdat%X,FILE_outputf,1)

! do boundary conditions
call scatterNDOFdataInit(meshdat,modeldat,myrank,FILE_outputf)
call scatterNDOFdata(meshdat,modeldat,myrank, modeldat%BOND,FILE_outputf,-1,9)

DEALLOCATE(modeldat%BOND,STAT=ERROR)
ALLOCATE(modeldat%BOND(ndof,nvertices),STAT=ERROR)
modeldat%BOND = 0d0

! we do not want to start scattering when not all the processors
! have finished allocating modeldat%BOND
call mpi_barrier(MPI_COMM_WORLD, error) 

call scatterNDOFdata(meshdat,modeldat,myrank, modeldat%BOND,FILE_outputf,1,10)


!write(*,*) "rank", getrank(), " slip?", NUMSLPglobal


if (NUMSLPglobal.ne.0) then


    if (allocated(modeldat%DIFORC)) then
!        write(*,*) "rank", getrank(), "has diforc allocated"
    else
!        write(*,*) "rank", getrank(), "has diforc not allocated"
    endif

!    write(*,*) "rank", getrank(), "a"

    call scatterNDOFdata(meshdat,modeldat,myrank, modeldat%DIFORC,FILE_outputf,-1,11)

!    write(*,*) "rank", getrank(), "b"

    DEALLOCATE(modeldat%DIFORC,STAT=ERROR)
!    ALLOCATE(modeldat%DIFORC(ndof,nvertices),STAT=ERROR)
    ALLOCATE(modeldat%DIFORC(ndof,meshdatactx%nvglobal),STAT=ERROR)

!    write(*,*) "rank", getrank(), "c"


    modeldat%DIFORC = 0d0

    ! we do not want to start scattering when not all the processors
    ! have finished allocating modeldat%DIFORC
    call mpi_barrier(MPI_COMM_WORLD, error) 

!    write(*,*) "rank", getrank(), "d"

    call scatterNDOFdata(meshdat,modeldat,myrank, modeldat%DIFORC,FILE_outputf,1,12)
!    write(*,*) "rank", getrank(), "z"

endif

!write(*,*) "rank", getrank(), " slip OK"


! do displacements

if (iecho.eq.8) then
    call parallellog("scatterdata", "scattering displacements")
endif 
call scatterONEdataInit(meshdat,modeldat,myrank,FILE_outputf)

DEALLOCATE(modeldat%D,STAT=ERROR)
ALLOCATE(modeldat%D(ndof,nvertices),STAT=ERROR)
modeldat%D = 0d0

if (NUMSLPglobal.ne.0) then

    DEALLOCATE(modeldat%DX,STAT=ERROR)
    ! this is using up extra memory for every vertex,
    ! while most likely only a tony minority of nodes
    ! is slippery. This should be optimised.
    ALLOCATE(modeldat%DX(ndof,nvertices),STAT=ERROR)
    modeldat%DX = 0d0
endif

! we do not want to start scattering when not all the processors
! have finished allocating modeldat%D and ...DX
call mpi_barrier(MPI_COMM_WORLD, error)

if (MODE.ge.3) then
    ! thermal round
    if (ITMODE.eq.1 .or. ITMODE.eq.3) then
        call scatterONEdata(meshdat,modeldat,myrank, modeldat%T,FILE_outputf,-1)
    endif

    DEALLOCATE(modeldat%T,STAT=ERROR)
    ALLOCATE(modeldat%T(nvertices),STAT=ERROR)
    modeldat%T = 0d0

    ! we do not want to start scattering when not all the processors
    ! have finished allocating modeldat%T
    call mpi_barrier(MPI_COMM_WORLD, error)

    if (ITMODE.eq.1 .or. ITMODE.eq.3) then
        call scatterONEdata(meshdat,modeldat,myrank, modeldat%T,FILE_outputf,1)
    endif

    if (NUMFNglobal+NUMSLPglobal.gt.0.and.IDIFT.eq.1) then
        DEALLOCATE(modeldat%TX,STAT=ERROR)
        ALLOCATE(modeldat%TX(nvertices),STAT=ERROR)
        modeldat%TX = 0d0
    endif

endif

if (iecho.eq.8) then
    call parallellog("sactterdata", "scattering done ")
endif 

return
end subroutine scatterdata


!***************************************************************
   subroutine scatterONEdataInit(meshdat,modeldat,myrank,FILE_outputf)
!***************************************************************
   USE MESHDATAMODULE
   USE MODELDATAMODULE
   USE MODELTOPOLOGY
!***************************************************************
   implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h"
!#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h"
!***************************************************************
   type (meshdata)  :: meshdat
   type (modeldata) :: modeldat
   PetscInt, ALLOCATABLE :: svertices(:)
   PetscErrorCode ierr
   integer i,j,myrank,FILE_outputf,ERROR
!***************************************************************
   ! create global vector for holding local data
   if (.not.oneVecSet) then

  call VecCreate(MPI_COMM_WORLD, oneVec, ierr)
  call VecSetSizes(oneVec, meshdat%Nvlocal, &
         meshdat%Nvglobal, ierr)
  call VecSetFromOptions(oneVec,ierr)
  call VecCreateSeq(MPI_COMM_SELF,nvertices, &
                oneVecLoc,ierr)
  call ISCreateStride(MPI_COMM_SELF, nvertices, &
  0, 1,islocald, ierr)

  ALLOCATE(svertices(nvertices),STAT=ERROR)
  do i=1,nvertices
      svertices(i) = vertices(i)-1
  enddo

  call ISCreateBlock(MPI_COMM_SELF,1, nvertices, &
  svertices,PETSC_COPY_VALUES,isglobald,ierr)

  call VecScatterCreate(oneVec, isglobald, oneVecLoc, &
  islocald, oneScatter, ierr)

  DEALLOCATE(svertices, STAT=ERROR)
  call ISDestroy(isglobald,ierr)
  call ISDestroy(islocald,ierr)
  oneVecSet = .true.
   endif
   end subroutine scatterONEdataInit
!***************************************************************

!***************************************************************
subroutine scatterTWOdataInit(meshdat, &
                              modeldat, &
                              FILE_outputf)
! builds the scattercontext 'twoScatter' which is used
! in subroutine scatterTWOdata to do the actual scattering.

!***************************************************************
USE MESHDATAMODULE   !,  only: meshdata
USE MODELDATAMODULE  !, only: modeldata
USE MODELTOPOLOGY
use debugmodule,        only: delay, &
                              verified_mpi_barrier
#ifdef SPARSE
use spetscmodule, only: PETSC_error
#else
use petscmodule, only: PETSC_error
#endif

use modelctx,           only: getrank
!***************************************************************
implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscsysdef.h"  ! for PetscInt and PetscErrorCode
!#include "petsc/finclude/petscvec.h"
!!#include "petsc/finclude/petscvecdef.h"

#include "petsc/finclude/petscvec.h" 
! for vecCreate and such


#include "petsc/finclude/petscis.h" 
! for the indexSet stuff

!***************************************************************
type (meshdata)  :: meshdat
type (modeldata) :: modeldat
PetscInt, ALLOCATABLE :: svertices(:)
PetscErrorCode   :: ierr
integer i,j,FILE_outputf,ERROR
!***************************************************************
! create global vector for holding local data

ierr=0


!write(*,*) 'rank', getrank(), "twoVecSet", twoVecSet
!write(*,*) 'rank', getrank(), "with nvertices", nvertices
!write(*,*) 'rank', getrank(), "nvloc, nvglob", meshdat%Nvlocal, meshdat%Nvglobal

if (.not.twoVecSet) then


!    call PETSc_VecCreate(twoVec)
    call VecCreate(MPI_COMM_WORLD, twoVec, ierr)
    call PETSC_error("scatterTWOdataInit","VecCreate",ierr)

    call VecSetSizes(twoVec,  &
                     meshdat%Nvlocal*2, &
                     meshdat%Nvglobal*2,  &
                     ierr)
    call PETSC_error("scatterTWOdataInit","VecSetSizes",ierr)



    call VecSetFromOptions(twoVec,ierr)
    call PETSC_error("scatterTWOdataInit","VecSetFromOptions",ierr)

    ! remember that nvertices is the total number of relevant vertices
    ! (including neighboring vertices and those of owned elements
    ! and on multiple processors, will be higher than meshdat%Nvlocal
    call VecCreateSeq(MPI_COMM_SELF, &
                      nvertices*2, &
                      twoVecLoc, &
                      ierr)
    call PETSC_error("scatterTWOdataInit","VecCreateSeq",ierr)

! from PETSc documentation:
!ISCreateStride(MPI_Comm comm,PetscInt n,PetscInt first,PetscInt step,IS *is)

!Input Parameters
!    comm   - the MPI communicator
!    n      - the length of the locally owned portion of the index set
!    first  - the first element of the locally owned portion of the index set
!    step   - the change to the next index

!Output Parameter
!   is     - the new index set 

! so IS islocald will contain running integers over the processors,
! a bit like the local ID's, but the size is different.

    call ISCreateStride(MPI_COMM_SELF,  &
                        nvertices*2, &
                        0,  &
                        1, &
                        islocald,  &
                        ierr)
    call PETSC_error("scatterTWOdataInit","ISCreateStride",ierr)

    ALLOCATE(svertices(2*nvertices),STAT=ERROR)

    do i=1,nvertices
        do j=1,2
            svertices(2*(i-1)+j) = 2*(vertices(i)-1) + j-1
        enddo
    enddo

!    write(*,*) "rank", getrank(), "scatter 2 init f", svertices


! create blocks of size 1? Is that useful at all?

!ISCreateBlock(MPI_Comm comm,PetscInt bs,PetscInt n,const PetscInt idx[],PetscCopyMode mode,IS *is)

!Input Parameters
!    comm     - the MPI communicator
!    bs      - number of elements in each block
!    n       - the length of the index set (the number of blocks)
!    idx     - the list of integers, one for each block and count of block not indices
!    mode     - see PetscCopyMode, only PETSC_COPY_VALUES and PETSC_OWN_POINTER are supported in this routine

!Output Parameter
!   is -the new index set 

    call ISCreateBlock(MPI_COMM_SELF, &
                       1,  &
                       2*nvertices, & 
                       svertices, &  
                       PETSC_COPY_VALUES, &
                       isglobald, &
                       ierr)
    call PETSC_error("scatterTWOdataInit","ISCreateBlock",ierr)



!    write(*,*) "rank", getrank(), "scatter 2 init g"


!    call delay()

!    call verified_mpi_barrier("About to vecScatterCreate with parameters twovec")

!    write(*,*) 'twovec'
!    call vecview(twoVec,PETSC_VIEWER_STDOUT_world, ierr)
!    call delay()
!    write(*,*) 'finished vecview with error', ierr

!    call verified_mpi_barrier("About to vecScatterCreate with parameters twovecloc")

!    write(*,*) 'twovecloc'
!    call vecview(twoVecLoc,PETSC_VIEWER_STDOUT_self, ierr)
!    call delay()
!    write(*,*) 'finished vecview with error', ierr

!    call verified_mpi_barrier("About to vecScatterCreate with parameters isglobald")

!    write(*,*) 'isglobald', getrank()
!    call ISview(isglobald,PETSC_VIEWER_STDOUT_self, ierr)
!    call delay()
!    write(*,*) 'finished isview with error', ierr

!    call verified_mpi_barrier("About to vecScatterCreate with parameters islocald")

!    write(*,*) 'islocald', getrank()
!    call ISview(islocald,PETSC_VIEWER_STDOUT_self, ierr) 
!    call delay()
!    write(*,*) 'finished isview with error', ierr

!    call delay()


!PetscErrorCode  VecScatterCreate(Vec xin,IS ix,Vec yin,IS iy,VecScatter *newctx)

!Input Parameters
!    xin     - a vector that defines the shape (parallel data layout of the vector) of vectors from which we scatter
!    yin     - a vector that defines the shape (parallel data layout of the vector) of vectors to which we scatter
!    ix     - the indices of xin to scatter (if NULL scatters all values)
!    iy     - the indices of yin to hold results (if NULL fills entire vector yin)

!Output Parameter
!   newctx -location to store the new scatter context 

! So this creates a scatter context that dictates how we have to scatter global
! data to local data.

! twovec and twovecloc are NOT of the same length!
! twovec is of the size of the global number of elements
! twovecloc is of the size of number of relevant vertices for this partition.

!    call verified_mpi_barrier("About to vecScatterCreate with parameters")


!    call verified_mpi_barrier("Here we go")


    call VecScatterCreate(twoVec, &
                          isglobald, &
                          twoVecLoc, &
                          islocald, &
                          twoScatter, &
                          ierr)
    call PETSC_error("scatterTWOdataInit","VecScatterCreate",ierr)


!    call verified_mpi_barrier("vecScatterCreate completed")


    DEALLOCATE(svertices, STAT=ERROR)


    call ISDestroy(isglobald,ierr)
    call ISDestroy(islocald,ierr)

    twoVecSet = .true.

endif

end subroutine scatterTWOdataInit
!***************************************************************



!***************************************************************
subroutine scatterNDOFdataInit(meshdat,modeldat,myrank,FILE_outputf)
!***************************************************************
USE MESHDATAMODULE
USE MODELDATAMODULE
USE MODELTOPOLOGY
use modelctx,       only: getrank
use debugmodule,    only: iecho, &
                          parallelLog
!***************************************************************
   implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h"
!#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h"
!***************************************************************
   type (meshdata)  :: meshdat
   type (modeldata) :: modeldat
   PetscInt, ALLOCATABLE :: svertices(:)
   PetscErrorCode ierr
   integer i,j,myrank,FILE_outputf,ERROR
!***************************************************************
   ! create global vector for holding local data


!write(*,*) "rank", getrank(), "ndofVecSet", ndofVecSet

if (.not.ndofVecSet) then

    call VecCreate(MPI_COMM_WORLD, ndofVec, ierr)

    call VecSetSizes(ndofVec, meshdat%Nvlocal*ndof, &
                     meshdat%Nvglobal*ndof, ierr)

    call VecSetFromOptions(ndofVec,ierr)

    call VecCreateSeq(MPI_COMM_SELF,ndof*nvertices, &
                      ndofVecLoc,ierr)

    call ISCreateStride(MPI_COMM_SELF, ndof*nvertices, &
                        0, 1,islocald, ierr)

    ALLOCATE(svertices(ndof*nvertices),STAT=ERROR)


    if (iecho.eq.8) then
        call parallelLog("scatterNDOFdata", "filling svertices")
    endif


    do i=1,nvertices
        do j=1,ndof
            svertices(ndof*(i-1)+j) = ndof*(vertices(i)-1) + j-1
        enddo
    enddo

    if (iecho.eq.8) then
        call parallelLog("scatterNDOFdata", "filling of svertices OK")
    endif


    call ISCreateBlock(MPI_COMM_SELF,1, ndof*nvertices, &
                       svertices,PETSC_COPY_VALUES,isglobald,ierr)

    call VecScatterCreate(ndofVec, isglobald, ndofVecLoc, &
                          islocald, ndofScatter, ierr)

    DEALLOCATE(svertices, STAT=ERROR)

    call ISDestroy(isglobald,ierr)
    call ISDestroy(islocald,ierr)

    ndofVecSet = .true.

endif

end subroutine scatterNDOFdataInit

!***************************************************************


!***************************************************************
subroutine scatterNDOFdata(meshdat,modeldat,myrank, &
    NDOFdata,FILE_outputf,IFLAG,callerID)
!***************************************************************
USE MESHDATAMODULE
USE MODELDATAMODULE
USE MODELTOPOLOGY
use modelctx,        only: getrank
!***************************************************************
   implicit none
!***************************************************************
   ! IFLAG = 0 run everything, this assumes that the data array has
   ! proper length (local size plus stride)
   ! IFLAG < 0 the data array still has local dimensions without
   ! extra for stride: perform only transfer of local data to global
   ! array
   ! IFLAG > 0 the data is already in global array and needs to be
   ! scattered back
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h"
!***************************************************************
type (meshdata)  :: meshdat
type (modeldata) :: modeldat

PetscScalar, pointer   :: tmpdata(:)
integer                :: callerID
PetscErrorCode ierr
integer i,j,myrank,ERROR,FILE_outputf,IFLAG
double precision :: NDOFdata(NDOF,*)
!***************************************************************
! transfer the local data in local array into global vector

!if (callerID.eq.7 .or. callerID.eq.8) then
!write(*,*) "rank", getrank(), "in scatterNDOF with IFLAG", IFLAG
!endif

if (IFLAG.le.0) then
    ! create a pointer tmpdata to the PETSc vector ndofVec
    call VecGetArrayF90(ndofVec, tmpdata,ierr)
!    write(*,*) "rank", getrank(), "geht loss"

    ! by editing tmpdata, we edit ndofVec implicitly
    do i=1,meshdat%Nvlocal
        do j=1,ndof
            tmpdata(ndof*(i-1)+j) = NDOFdata(j,i)
        enddo
    enddo
!    write(*,*) "rank", getrank(), "geht nicht mehr loss"

    ! let PETSc process the edits.
    call VecRestoreArrayF90(ndofVec, tmpdata,ierr)
endif

!if (callerID.eq.7 .or. callerID.eq.8) then
!write(*,*) "rank", getrank(), "second bit of scatterNDOF", iflag
!endif

if (IFLAG.ge.0) then
    call VecScatterBegin(ndofScatter, ndofVec, ndofVecLoc, &
                         INSERT_VALUES, SCATTER_FORWARD, ierr)
!    if (callerID.eq.7 .or. callerID.eq.8) then
!        write(*,*) "rank", getrank(), "vescatter begun", ierr
!    endif

    call VecScatterEnd(ndofScatter, ndofVec, ndofVecLoc, &
                       INSERT_VALUES, SCATTER_FORWARD, ierr)

!    if (callerID.eq.7 .or. callerID.eq.8) then
!        write(*,*) "rank", getrank(), "vecscatter done", ierr
!    endif

    call VecGetArrayF90(ndofVecLoc, tmpdata,ierr)
    do i=1,nvertices
        do j=1,ndof
            NDOFdata(j,i) = tmpdata(ndof*(i-1)+j)
        enddo
    enddo
    call VecRestoreArrayF90(ndofVecLoc, tmpdata,ierr)

!    if (callerID.eq.7 .or. callerID.eq.8) then
!    write(*,*) "rank", getrank(), "done"
!    endif

endif

nullify(tmpdata)

end subroutine scatterNDOFdata
!***************************************************************


!***************************************************************
   subroutine scatterNSDdataInit(meshdat,modeldat,myrank,FILE_outputf)
!***************************************************************
   USE MESHDATAMODULE
   USE MODELDATAMODULE
   USE MODELTOPOLOGY
!***************************************************************
   implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h"
!#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h"
!***************************************************************
   type (meshdata)  :: meshdat
   type (modeldata) :: modeldat
   PetscInt, ALLOCATABLE :: svertices(:)
   PetscErrorCode ierr
   integer i,j,myrank,ERROR,FILE_outputf
!***************************************************************


   if (.not.nsdVecSet) then


  call VecCreate(MPI_COMM_WORLD, nsdVec, ierr)

!          write(*,*) 'vec created'

  call VecSetSizes(nsdVec, meshdat%Nvlocal*nsd, &
         meshdat%Nvglobal*nsd, ierr)

  call VecSetFromOptions(nsdVec,ierr)

  call VecCreateSeq(MPI_COMM_SELF,nsd*nvertices, &
                nsdVecLoc,ierr)

  call ISCreateStride(MPI_COMM_SELF, nsd*nvertices, &
  0, 1,islocald, ierr)

  ALLOCATE(svertices(nsd*nvertices),STAT=ERROR)

  do i=1,nvertices
      do j=1,nsd
          svertices(nsd*(i-1)+j) = nsd*(vertices(i)-1) + j-1
      enddo
  enddo

  call ISCreateBlock(MPI_COMM_SELF, &
                     1,                  & ! block size &
                     nsd*nvertices,      & ! number of blocks &
                     svertices,          & ! list of integers &
                     PETSC_COPY_VALUES,  & ! values or pointers &
                     isglobald, &
                     ierr)

  call VecScatterCreate(nsdVec, isglobald, nsdVecLoc, &
  islocald, nsdScatter, ierr)

  call ISDestroy(isglobald,ierr)
  call ISDestroy(islocald,ierr)

  DEALLOCATE(svertices, STAT=ERROR)
  nsdVecSet = .true.
   endif

   end subroutine scatterNSDdataInit
!***************************************************************


!***************************************************************
   subroutine scatterNSDdata(meshdat,modeldat,myrank, &
    NSDdata,FILE_outputf,IFLAG)
!***************************************************************
   USE MESHDATAMODULE,  only: meshdata
   USE MODELDATAMODULE, only: modeldata
   USE MODELTOPOLOGY
!***************************************************************
   implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h"
!***************************************************************
   type (meshdata)  :: meshdat
   type (modeldata) :: modeldat
   PetscScalar, pointer   :: tmpdata(:)
   PetscErrorCode ierr
   integer i,j,myrank,FILE_outputf, IFLAG
   double precision :: NSDdata(NSD,*)
!***************************************************************
   if (IFLAG.le.0) then
  call VecGetArrayF90(nsdVec, tmpdata,ierr)
  do i=1,meshdat%Nvlocal
      do j=1,nsd
          tmpdata(nsd*(i-1)+j) = NSDdata(j,i)
      enddo
  enddo
  call VecRestoreArrayF90(nsdVec, tmpdata,ierr)
   endif

   if (IFLAG.ge.0) then
  call VecScatterBegin(nsdScatter, nsdVec, nsdVecLoc, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(nsdScatter, nsdVec,  nsdVecLoc, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecGetArrayF90(nsdVecLoc, tmpdata,ierr)
  do i=1,nvertices
      do j=1,nsd
          nsddata(j,i) = tmpdata(nsd*(i-1)+j)
      enddo
  enddo
  call VecRestoreArrayF90(nsdVecLoc, tmpdata,ierr)
   endif

   nullify(tmpdata)

   end subroutine scatterNSDdata
!***************************************************************


!***************************************************************
   subroutine scatterONEdata(meshdat,modeldat,myrank, &
    ONEdata,FILE_outputf,IFLAG)
!***************************************************************
   USE MESHDATAMODULE
   USE MODELDATAMODULE
   USE MODELTOPOLOGY
!***************************************************************
   implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h"
!***************************************************************
   type (meshdata)  :: meshdat
   type (modeldata) :: modeldat
   PetscScalar, pointer   :: tmpdata(:)
   PetscErrorCode ierr
   integer i,j,myrank,FILE_outputf, IFLAG
   double precision :: ONEdata(*)
!***************************************************************
   if (IFLAG.le.0) then
  call VecGetArrayF90(oneVec, tmpdata,ierr)
  do i=1,meshdat%Nvlocal
      tmpdata(i) = ONEdata(i)
  enddo
  call VecRestoreArrayF90(oneVec, tmpdata,ierr)
   endif

   if (IFLAG.ge.0) then
  call VecScatterBegin(oneScatter, oneVec, &
  oneVecLoc, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecScatterEnd(oneScatter, oneVec, &
  oneVecLoc, INSERT_VALUES, SCATTER_FORWARD, ierr)
  call VecGetArrayF90(oneVecLoc, tmpdata,ierr)
  do i=1,nvertices
      ONEdata(i) = tmpdata(i)
  enddo
  call VecRestoreArrayF90(oneVecLoc, tmpdata,ierr)
   endif

   nullify(tmpdata)

   end subroutine scatterONEdata
!***************************************************************


!***************************************************************
subroutine scatterTWOdata(meshdat,modeldat, &
    TWOdata,FILE_outputf,IFLAG)
!***************************************************************
USE MESHDATAMODULE
USE MODELDATAMODULE
USE MODELTOPOLOGY
use modelctx,       only: getrank
use debugmodule,    only: delay
!***************************************************************
   implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h"
!***************************************************************
type (meshdata)  :: meshdat
type (modeldata) :: modeldat

PetscScalar, pointer   :: tmpdata(:)
PetscErrorCode ierr

integer i,j,FILE_outputf, IFLAG
double precision :: TWOdata(2,*)
!***************************************************************

!write(*,*) 'scattertwodata 1', iflag

!call delay()

if (IFLAG.le.0) then
    call VecGetArrayF90(twoVec, tmpdata,ierr)
    do i=1,meshdat%Nvlocal
        do j=1,2
            tmpdata(2*(i-1)+j) = TWOdata(j,i)
        enddo
    enddo
    call VecRestoreArrayF90(twoVec, tmpdata,ierr)
endif

!write(*,*) 'scattertwodata 2'

! only start scattering when all threads have completed their twoVec
call mpi_barrier(MPI_COMM_WORLD, ierr)


if (IFLAG.ge.0) then
    call VecScatterBegin(twoScatter, &
                         twoVec, &
                         twoVecLoc, &
                         INSERT_VALUES, &
                         SCATTER_FORWARD, &
                         ierr)

    call VecScatterEnd(twoScatter, &
                       twoVec, &   
                       twoVecLoc, &
                       INSERT_VALUES, &  
                       SCATTER_FORWARD, &
                       ierr)

    call VecGetArrayF90(twoVecLoc, &
                        tmpdata, &
                        ierr)
    do i=1,nvertices
        do j=1,2
            twodata(j,i) = tmpdata(2*(i-1)+j)
        enddo
    enddo

    call VecRestoreArrayF90(twoVecLoc, &
                            tmpdata, &
                            ierr)

endif

nullify(tmpdata)

!write(*,*) "rank", getrank(), 'done scatterTWO'

end subroutine scatterTWOdata
!***************************************************************



!***************************************************************
   subroutine getallvectorentries(FILE_outputf,vec,IFLAG)
!***************************************************************
   USE MODELDEFINITION
!***************************************************************
   implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h"
!#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h"
#include "mpif.h"
!***************************************************************
   PetscScalar, pointer   :: tmp(:)
   Vec                    :: vec,vecloc
   VecScatter scatter
   PetscErrorCode ierr
   integer FILE_outputf,i, IFLAG
!***************************************************************

!      write(*,*) 'create vecloc'
   if (IFLAG.eq.0) then
  call VecCreateSeq(MPI_COMM_SELF,NEQglobal,vecloc,ierr)
   else
  call VecCreateSeq(MPI_COMM_SELF,NTEQglobal,vecloc,ierr)
   endif

   call VecScatterCreateToAll(vec,scatter, vecloc,ierr)

   call VecScatterBegin(scatter, vec, vecloc, INSERT_VALUES, SCATTER_FORWARD,ierr)
   call VecScatterEnd  (scatter, vec, vecloc, INSERT_VALUES, SCATTER_FORWARD,ierr)

   call VecGetArrayF90(vecloc,tmp,ierr)
   call VecRestoreArrayF90(vecloc,tmp,ierr)

!      write(*,*) 'destroy vecloc'
   call VecDestroy(vecloc,ierr)

   call VecScatterDestroy(scatter,ierr)

   nullify(tmp)

   end subroutine getallvectorentries
!***************************************************************


!***************************************************************
#ifdef OBSOLETE
   subroutine createscattering(meshdat,myrank,FILE_outputf,nen,nsd)
!***************************************************************
   USE MESHDATAMODULE
!***************************************************************
   implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h"
!#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h"
!#include "mpif.h"
!***************************************************************
   type (meshdata)       :: meshdat
   integer myrank, FILE_outputf, nen, nsd, ERROR,i,j
   PetscErrorCode        :: ierr
   PetscScalar, pointer  :: tmpcoords(:)
   PetscInt, pointer     :: ispoints(:)
   PetscInt, ALLOCATABLE :: svertices(:)
!***************************************************************
!          write(*,*) 'create localcoords'

   call VecCreateSeq(MPI_COMM_SELF,nsd*nvertices, &
  meshdat%localcoords, ierr)
   call ISCreateStride(MPI_COMM_SELF, nsd*nvertices, &
  0, 1,islocal, ierr)
   ALLOCATE(svertices(nsd*nvertices),STAT=ERROR)
   do i=1,nvertices
  do j=1,nsd
      svertices(nsd*(i-1)+j) = nsd*(vertices(i)-1) + j-1
  enddo
   enddo

   call ISCreateBlock(MPI_COMM_SELF,1, nsd*nvertices, &
  svertices,isglobal,ierr)


   DEALLOCATE(svertices, STAT=ERROR)

   call ISGetIndicesF90(islocal, ispoints, ierr)
   write(FILE_outputf,*) 'local', (ispoints(i),i=1,nsd*nvertices)
   call ISRestoreIndicesF90(islocal, ispoints, ierr)
   call ISGetIndicesF90(isglobal, ispoints, ierr)
   write(FILE_outputf,*) 'global', (ispoints(i),i=1,nsd*nvertices)
   call ISRestoreIndicesF90(isglobal, ispoints, ierr)

   call VecScatterCreate(coords, isglobal, meshdat%localcoords, &
  islocal, meshdat%scatter, ierr)
   call VecScatterBegin(meshdat%scatter, coords, &
  meshdat%localcoords, INSERT_VALUES, SCATTER_FORWARD, ierr)
   call VecScatterEnd(meshdat%scatter, coords, &
  meshdat%localcoords, INSERT_VALUES, SCATTER_FORWARD, ierr)

!      write(FILE_outputf,*) 'coords', meshdat%Nvlocal
!      call printpetsccoordinates(coords, nsd,
!     .   meshdat%Nvlocal, FILE_outputf)
!      call printpetsccoordinates(meshdat%localcoords,
!     .  nsd, nvertices, FILE_outputf)

   call ISDestroy(isglobal,ierr)
   call ISDestroy(islocal,ierr)

   DEALLOCATE(meshdat%X,STAT=ERROR)
!      ALLOCATE(meshdat%localX(nsd,nvertices),STAT=ERROR)
   ALLOCATE(meshdat%X(nsd,nvertices),STAT=ERROR)

   call VecGetArrayF90(meshdat%localcoords, tmpcoords,ierr)
   do i=1,nvertices
  do j=1,nsd
      meshdat%X(j,i) = tmpcoords(nsd*(i-1)+j)
  enddo
   enddo
   call VecRestoreArrayF90(meshdat%localcoords, tmpcoords,ierr)

   nullify(tmpcoords)
   nullify(ispoints)

   end subroutine createscattering
#endif 
! OBSOLETE
!***************************************************************


!-----routine to print the mapping from global (application) to local (petsc) numbering
!***************************************************************
subroutine writeAOmapping(FILE_outputf, myrank, meshdat)
!***************************************************************
USE MESHDATAMODULE, only: meshdata
!***************************************************************
type (meshdata) :: meshdat
integer FILE_outputf, myrank,  i
!***************************************************************

do i=1,meshdat%Nvlocal
    write(FILE_outputf,11) meshdat%gloInd(i), meshdat%locInd(i)
11 format(1x,'Global vertex ', I8, ' maps on local vertex ', I8)
enddo

do i=1,meshdat%Nelocal
    write(FILE_outputf,12) meshdat%gloElt(i), meshdat%locElt(i)
12 format(1x,'Global element ', I8, ' maps on local element ', I8)
enddo

end subroutine writeAOmapping
!***************************************************************



!-----routine to print mesh data prior and after application orderings
!***************************************************************
   subroutine writeAOmeshdata(FILE_outputf, myrank, meshdat, &
  nbrlist, nen, flag)
!***************************************************************
   USE MESHDATAMODULE, only: meshdata
!***************************************************************
   type (meshdata) :: meshdat
   integer FILE_outputf, myrank, flag, nen, i, j
   integer nbrlist(*)
!***************************************************************

   if (flag.eq.0) then
  write(FILE_outputf,*) 'Before AOApplicationToPetsc'
   else if (flag.eq.1) then
  write(FILE_outputf,*) ' ----- ----- ----- -----  '
  write(FILE_outputf,*) 'After AOApplicationToPetsc'
   else if (flag.eq.2) then
  write(FILE_outputf,*) ' ----- ----- ----- -----  '
  write(FILE_outputf,*) 'After contiguous ordering '
  goto 100
   endif

   write(FILE_outputf, 11)
  11  format(5x,'local indices are :')
   write(FILE_outputf, *) (meshdat%locInd(i), i=1,meshdat%Nvlocal)

   if (flag.ne.2) goto 150

100    write(FILE_outputf,*) 'The array vertices is : ', &
(vertices(i),i=1,nvertices)
   do i=1,meshdat%Nvglobal
  write(FILE_outputf,14) i-1,i-1,vertices(i),i-1,verticesmask(i)
   enddo
14 format(2x, 'vertex ', I8, '; vertices[',I8,']: ', I8, &
' and verticesmask[',I8,']: ', I8 )


150   do i=1,meshdat%Nvlocal
  write(FILE_outputf,*) "Neighbors of local vertex ", &
    meshdat%locInd(i), " are : ", &
    (meshdat%AdjM(j,i), j=1,meshdat%itot(i))
   enddo

   if (flag.eq.2) goto 200

   write(FILE_outputf,13) meshdat%Nvneighborstotal
  13  format(3x, I8, ' neighbors :')
   write(FILE_outputf,*) (nbrlist(i),i=1,meshdat%Nvneighborstotal)

   write(FILE_outputf, 12)
  12  format(5x,'local elements are :')
   write(FILE_outputf, *) (meshdat%locElt(i),i=1,meshdat%Nelocal)

200   do i=1,meshdat%Nelocal
  write(FILE_outputf,*) "Vertices of local element ", &
    meshdat%locElt(i), " are : ", &
    (meshdat%ien(j,i), j=1,NEN)
   enddo


   end subroutine writeAOmeshdata
!***************************************************************


!-----routien to write mesh data after reading input
!***************************************************************
   subroutine writemeshdata(FILE_outputf, myrank, meshdat)
!***************************************************************
   USE MESHDATAMODULE
!***************************************************************
   type (meshdata) :: meshdat
   integer FILE_outputf, myrank
   integer i,j
!***************************************************************
   write(FILE_outputf,*) 'MESH INFORMATION FOR PROCESSOR ', myrank
   write(FILE_outputf,*) '-----------------------------------'
   do i=1,meshdat%Nvlocal
  write(FILE_outputf,11) meshdat%gloInd(i),meshdat%X(1,i), &
                  meshdat%X(2,i)
  write(FILE_outputf,12) meshdat%itot(i), &
             (meshdat%AdjM(j,i),j=1,meshdat%itot(i))
   enddo
   do i=1,meshdat%Nelocal
  write(FILE_outputf,13) meshdat%gloElt(i), &
             (meshdat%ien(j,i),j=1,3), meshdat%mat(i)
   enddo
11 format(1x,'vertex ', i5, ' has coordinates ', 2G14.6)
12 format(4x,'with ', i5, ' neighbours: ', 10I7)
13 format(1x,'Element ', i5, ' has vertices ', 3I8, &
     ' and material ', I5)

   end subroutine writemeshdata
!***************************************************************


!***************************************************************
   subroutine createlocaltoglobalmapping(mat, vec, IFLAG)
!***************************************************************
!#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscmat.h"
!#include "petsc/finclude/petscvecdef.h"
!#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscis.h"
!#include "mpif.h"
!***************************************************************
!     IFLAG=0 mech
!     IFLAG=1 thermal

   Mat mat
   Vec vec
   PetscErrorCode ierr
   integer IFLAG,i
!***************************************************************
   if (IFLAG.eq.0) then
  equations = equations - 1

  call ISLocalToGlobalMappingCreate(MPI_COMM_SELF, 1, nequations, &
    equations, PETSC_COPY_VALUES,isl2g, ierr)

! use identical row mapping and column mapping.
  call MatSetLocalToGlobalMapping(mat, isl2g, isl2g, ierr)
  call VecSetLocalToGlobalMapping(vec, isl2g, ierr)
   else
  equationst = equationst - 1

  call ISLocalToGlobalMappingCreate(MPI_COMM_SELF, 1, nequationst, &
    equationst, PETSC_COPY_VALUES,isl2gt, ierr)

! use identical row mapping and column mapping.
  call MatSetLocalToGlobalMapping(mat, isl2gt, isl2gt, ierr)
  call VecSetLocalToGlobalMapping(vec, isl2gt, ierr)
   endif

   end subroutine createlocaltoglobalmapping
!***************************************************************


!------------------------------------------------------------------------------
!     renumbers elements from application to local petsc
!------------------------------------------------------------------------------
!***************************************************************
subroutine renumberelements(ellist, ellistSize, caller)
!***************************************************************
USE MESHDATAMODULE
use modelctx,      only: getsize, getrank
use debugmodule,   only: iecho, parallelLog
!***************************************************************
implicit none
!***************************************************************
PetscInt ellistSize
PetscInt ellist(*)
character(len=*) caller
integer i, elloc


if (ellistSize.eq.0) then
    ! this partition does not have any slippery nodes.
    ! It does not mean there are no slippery nodes at all in the domain...
    return
endif

!write(*,*) 'rank', getrank(), 'Renumbering elements for ', trim(caller)
!write(*,*) 'rank', getrank(), 'with ellistnum: ', ellistnum
!write(*,*) "rank", getrank(), "entered renum elems for elem", trim(caller)

do i=1,ellistSize

#ifdef SPARSE
    if (getsize().eq.1) then
        elloc = ellist(i)
    else  
!        write(*,*) "rank", getrank(), "reading ellist: ", i, ellist(i)
        elloc = meshdatactx%Glo2LocElement(ellist(i))
    endif
#else
    elloc = ellist(i)
#endif

    if (elloc.gt.0) then
        ellist(i) = elloc
    endif
enddo

!write(*,*) 'Renumbering completed'

end subroutine renumberelements
!***************************************************************



!------------------------------------------------------------------------------
!     renumbers nodal points from application to petsc
!------------------------------------------------------------------------------
!***************************************************************
subroutine renumbernodalpoints(aonp, nplist, nplistnum, ierr)
!***************************************************************
implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscao.h"
!***************************************************************
AO                    :: aonp
PetscInt              :: nplistnum
PetscInt              :: nplist(*)
PetscErrorCode        :: ierr
PetscInt, ALLOCATABLE :: tmparray(:)
integer               :: ERROR,i
!***************************************************************
ALLOCATE(tmparray(nplistnum), STAT=ERROR)
! assumes nplist is base 1, convert to base 0 for Petsc
do i=1,nplistnum
    tmparray(i) = nplist(i) - 1
enddo

call AOApplicationToPetsc(aonp, nplistnum,tmparray,ierr)

do i=1,nplistnum
    nplist(i) = tmparray(i)
enddo

DEALLOCATE(tmparray)

end subroutine renumbernodalpoints
!***************************************************************



!------------------------------------------------------------------------------
!     renumbers equations numbers from application to petsc
!------------------------------------------------------------------------------

!call renumberequations(aoeq, modeldat%ID, size(modeldat%ID,2), &
!                          totideqs, &
!            meshdat%Nvlocal,NDOF, ierr,FILE_outputf,"normal")


!***************************************************************
subroutine renumberequations(aovariant, &     ! renumbering table
                             eqarray, &       ! array to be renumbered
                             nPoints, &       ! size of this array   ** debugging only
                             eqnum, &         ! number of equations  ** debugging only
                             szarray, &       ! number of points in this partition
                             ndof, &
                             ierr, &
                             FILE_outputf, &
                             callerName)


use modelctx,    only: getrank
use debugmodule, only: debug

!***************************************************************
implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscao.h"
#include "petsc/finclude/petscsys.h"
!***************************************************************
AO                    :: aovariant
PetscInt              :: eqnum, ndof, szarray
integer               :: nPoints
PetscInt              :: eqarray(NDOF,nPoints)
PetscErrorCode        :: ierr
PetscInt, ALLOCATABLE :: tmpid(:)
integer               :: ERROR, i, next, k, FILE_outputf
character(len=*)      :: callerName
!***************************************************************

if (allocated(tmpid)) then
    deallocate(tmpid)
endif
ALLOCATE(tmpid(eqnum), STAT=ERROR)
if (error.ne.0) then
    stop "Could not allocate tmpid; leaving"
endif

next = 0

do i=1,szarray
    do k=1,NDOF
!        write(*,*) callerName, "trying to access", k,i, "of", size(eqarray,1), size(eqarray,2)
        if (eqarray(k,i).ne.0) then
            next = next + 1
!            write(*,*) "temping array loc", k, i,"of", &
!                         size(eqarray,1), nPoints
!            write(*,*) "tmpid ", next, "of", eqnum
            tmpid(next) = eqarray(k,i) - 1
        endif
    enddo
enddo


!if (debug) then
!write(*,*) 'rank ', getrank(), 'AOApp2PETSc ', trim(callerName), tmpid
!call AOview(aovariant, PETSC_VIEWER_STDOUT_WORLD, ierr)
!endif

call AOApplicationToPetsc(aovariant, eqnum, tmpid, ierr)

!if (debug) then
!write(*,*) 'rank ', getrank(), 'AOApp2PETSc done'
!endif

next = 0
do i=1,szarray
    do k=1,NDOF
        if (eqarray(k,i).ne.0) then
            next = next + 1
!            write(*,*) "writing array loc", k, i

            eqarray(k,i) = tmpid(next) + 1
        endif
    enddo
enddo


!write(*,*) "rank", getrank(), "has eqarray", eqarray

DEALLOCATE(tmpid)

!write(*,*) "leaving renumberequations"


end subroutine renumberequations
!***************************************************************



!***************************************************************
subroutine renumberThermalEquations(aovariant, &     ! renumbering table
                             eqarray, &       ! array to be renumbered
                             nPoints, &       ! size of this array   ** debugging only
                             eqnum, &         ! number of equations  ** debugging only
                             szarray, &       ! number of points in this partition
                             ndof, &
                             ierr, &
                             FILE_outputf, &
                             callerName)


use modelctx,    only: getrank
use debugmodule, only: debug

!***************************************************************
implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
!#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscao.h"
#include "petsc/finclude/petscsys.h"
!***************************************************************
AO                    :: aovariant
PetscInt              :: eqnum, ndof, szarray
integer               :: nPoints
PetscInt              :: eqarray(nPoints)
PetscErrorCode        :: ierr
PetscInt, ALLOCATABLE :: tmpid(:)
integer               :: ERROR, i, next, FILE_outputf
character(len=*)      :: callerName
!***************************************************************

if (allocated(tmpid)) then
    deallocate(tmpid)
endif
ALLOCATE(tmpid(eqnum), STAT=ERROR)
if (error.ne.0) then
    stop "Could not allocate tmpid; leaving"
endif

next = 0

do i=1,szarray
        if (eqarray(i).ne.0) then
            next = next + 1
            tmpid(next) = eqarray(i) - 1
        endif
enddo

call AOApplicationToPetsc(aovariant, eqnum, tmpid, ierr)

next = 0
do i=1,szarray
        if (eqarray(i).ne.0) then
            next = next + 1
            eqarray(i) = tmpid(next) + 1
        endif
enddo

DEALLOCATE(tmpid)

end subroutine renumberThermalEquations
!***************************************************************








!***************************************************************
   subroutine writeinputdata(meshdat, modeldat, myrank, FILE_outputf, &
  nen,ndof,nstr, flag)
!***************************************************************
   USE MODELDATAMODULE, only: modeldata, &
                         HasSlipperyNodes
   USE MODELDEFINITION
   USE MESHDATAMODULE,  only: meshdata
!***************************************************************
   implicit none
!***************************************************************
!#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscsys.h"
   type (meshdata) :: meshdat
   type (modeldata) :: modeldat
   integer myrank, FILE_outputf, nen, ndof, nstr,flag
   integer i,j
   logical doprint
!***************************************************************


   if (flag.eq.0) write(FILE_outputf, 11)
11 format(1x,'input data in global application numbering')
   if (flag.eq.1) write(FILE_outputf, 12)
12 format(1x,'input data in local petsc numbering')


   write(FILE_outputf,*) 'EQUATON NUMBERS'

   if (HasSlipperyNodes()) then
  do i=1,meshdat%Nvlocal
      write(FILE_outputf,33) i, (modeldat%IDX(j,i),j=1,NDOF), &
       (modeldat%ID(j,i),j=1,NDOF)
  enddo
   else
  do i=1,meshdat%Nvlocal
      write(FILE_outputf,32) i, (modeldat%ID(j,i),j=1,NDOF)
  enddo
   endif

   write(FILE_outputf,*) 'BOUNDARY CONDITIONS'

   do i=1,meshdat%Nvlocal
  write(FILE_outputf,41) (modeldat%IBOND(j,i), &
                          modeldat%BOND(j,i),j=1,NDOF)
   enddo

   if (NLINK.gt.0) then
  write(FILE_outputf,*) 'DEFINITION OF LINKED COORDINATES'
  do i=1,NLINK
      write(FILE_outputf,32) (modeldat%LINK(j,i),j=1,3)
  enddo
   endif

   if (NUMROT.gt.0) then
  write(FILE_outputf,*) 'DEFINITION OF LOCAL DOF ROTATIONS'
  do i=1,meshdat%Nvlocal
      doprint = modeldat%SKEW(1,i).ne.0d0 .or. &
                modeldat%SKEW(2,i).ne.0d0
      if (doprint) then
          write(FILE_outputf, 31) meshdat%gloInd(i), i, &
            (modeldat%SKEW(j,i),j=1,2)
      endif
  enddo
   endif

   if (NPRE.gt.0) then
  write(FILE_outputf,*) 'DEFINITION OF INITIAL STRESSES'
  do i=1,NPRE
      write(FILE_outputf,21) modeldat%ISELM(i),modeldat%ISTIME(i), &
         (modeldat%STN0(j,i),j=1,NSTR)
  enddo
   endif

   if (NUMPR.gt.0) then
  write(FILE_outputf,*) 'DEFINITION OF PRESSURE BOUNDARY CONDITIONS'
  do i=1,NUMPR
      write(FILE_outputf,22) modeldat%IELNO(i), modeldat%ISIDE(i), &
       modeldat%PRES(i)
  enddo
   endif

   if (NUMSTR.gt.0) then
  write(FILE_outputf,*) 'DEFINITION OF STRESS BOUNDARY CONDITIONS'
  do i=1,NUMSTR
      write(FILE_outputf,23) modeldat%IELSTR(i), modeldat%ISSIDE(i), &
       (modeldat%ISTR(j,i),j=1,2),(modeldat%STRS(j,i),j=1,NSTR)
  enddo
   endif

   if (NUMWNK.gt.0) then
  write(FILE_outputf,*) 'DEFINITION OF WINKLER RESTORING PRESSURES'
  do i=1,NUMWNK
      write(FILE_outputf,24) modeldat%IWELM(i), modeldat%IWSIDE(i), &
       modeldat%IWTIME(i), modeldat%WPRES(i)
  enddo
   endif

   if (NUMFN.gt.0) then
  write(FILE_outputf,*) 'DEFINITION OF FAULTED (SPLIT) NODES'
  do i=1,NUMFN
      write(FILE_outputf,25) (modeldat%NFAULT(j,i),j=1,3), &
         (modeldat%FAULT(j,i),j=1,NDOF)
  enddo
   endif

   if (HasSlipperyNodes()) then
  write(FILE_outputf,*) 'DEFINITION OF SLIPPERY NODES'
  do i=1,NUMSLP
      write(FILE_outputf,26) (modeldat%NSLIP(j,i),j=1,5), &
         (modeldat%DIFORC(j,i),j=1,NDOF)
  enddo
  write(FILE_outputf,*) '   SLIPPERY NODES'
  do i=1,NUMSN
      write(FILE_outputf,'(i6)') modeldat%IDSLP(i)
  enddo
  write(FILE_outputf,*) '   ELEMENTS WITH SLIPPERY NODES'
  do i=1,NUMSE
      write(FILE_outputf,'(i6)') modeldat%IDSLE(i)
  enddo

  if (NSLSKEW.gt.0) then
      write(FILE_outputf,*) '    FAULT PARALLEL SLIPPERY ELEMENT'
      do i=1,NSLSKEW
          write(FILE_outputf,'(2i6)') (modeldat%NSELSD(j,i),j=1,2)
      enddo
  endif

  if (NWINKX.gt.0) then
      write(FILE_outputf,*) '    DEFINITION OF DIFFERENTIAL WINKLER ', &
                          'FORCES'
      do i=1,meshdat%Nvlocal
          doprint = .false.
          do j=1,NDOF
              doprint = doprint .or. (modeldat%IWINX(j,i).ne.0)
          enddo
          if (doprint) then
              write(FILE_outputf,27) meshdat%gloInd(i),i, &
             (modeldat%IWINX(j,i),j=1,2), &
             (modeldat%WINX(j,i),j=1,NDOF)
          endif
      enddo
  endif
   endif

   if (NSURF.gt.0) then
  write(FILE_outputf,*) 'DEFINITION OF SURFACE NODES'
  do i=1,NSURF
      write(FILE_outputf,'(i6)') modeldat%ISURF(i)
  enddo
   endif

   if (NTANOM.gt.0) then
  write(FILE_outputf,*) 'DEFINITION OF THERMAL ANOMALIES'
  do i=1,NTANOM
      write(FILE_outputf, 24) (modeldat%ITANOM(j,i),j=1,3), &
         modeldat%TANOM(i)
  enddo
   endif

   if (NFLX.gt.0) then
  write(FILE_outputf,*) 'DEFINITION OF HEAT FLUX'
  do i=1,NFLX
      write(FILE_outputf,22) modeldat%IFLX(i), modeldat%IFLS(i), &
             modeldat%BFLX(i)
  enddo
   endif

   if (NTWINK.gt.0) then
  write(FILE_outputf,*) 'DEFINITION OF THERMAL WINKLER BC'
  do i=1,NTWINK
      write(FILE_outputf,24) (modeldat%ITWINK(j,i),j=1,3), &
           modeldat%TWINK(i)
  enddo
   endif


21 format(2i6,6g14.6)
22 format(2i6,1g14.6)
23 format(4i6,6g14.6)
24 format(3i6,1g14.6)
25 format(3i6,3g14.6)
26 format(5i6,5g14.6)
27 format(4i6,4g14.6)
31 format(2i6,3g14.6)
32 format(4i6)
33 format(7i6)
41 format(6(i6,G14.6))
   end subroutine writeinputdata
!***************************************************************


!***************************************************************
!-----routine to write model data needed by postprocessing mergefiles
subroutine writemodeldata(meshdat, modeldat,outputcontroldat, &
                          timestepdat, luin)
!***************************************************************
USE MESHDATAMODULE,  only: meshdata
USE MODELDATAMODULE, only: modeldata, modeldatactx
USE MODELDEFINITION
USE MODELTOPOLOGY
USE TIMESTEPMODULE
USE MODELCTX
USE FILEMODULE,      only: outdir
use debugmodule, only: iecho, &
                       parallelLog, delay
!***************************************************************
implicit none
!***************************************************************
type (modeldata) modeldat
type (meshdata) meshdat
type(outputcontroldata) :: outputcontroldat
type(timestepdata) :: timestepdat

character(len=255) modelpartname
character(len=40) frmt
integer luin
integer i,j,rnk, inneri, innerj, logsize, lnblk

external lnblk
!***************************************************************

if (iecho.eq.8) then
    call ParallelLog("writemodeldata", "just entered")
endif

if (getrank().eq.0) then
    call openf(luin,trim(outdir)//'/'//'modeldata.dat','unknown')

#ifdef EXPERIMENTAL_ReducedSlipAlloc
    write(luin, '(2i12)') meshdat%nvglobal, meshdat%neglobal
#endif

#ifdef EXPERIMENTALTRANSIENTDISCREEP
    write(luin, '(6i12)') NSD,NDOF,NSTR,NEN,NUMAT,NSTRTAY
#else
    write(luin, '(5i12)') NSD,NDOF,NSTR,NEN,NUMAT
#endif

    write(luin, '(8i12)') MODE, ICVIS, FaultVelo, NSED, &
      NSURFglobal, NUMFNglobal,NUMSLPglobal,NUMSNglobal

    if (mode.ge.3) then
        write(luin, '(2i12)') IADVEC, IDIFT
    endif

    write(luin, '(I12)') NMPRT

    if (NMPRT.gt.0) then
        write(luin, 3) (outputcontroldat%IMPRINT(i),i=1,NMPRT)
    endif

    if (mode.ge.3) then
        write(luin, '(I12)') NTPRT
        if (NTPRT.gt.0) then
            write(luin, 3) (outputcontroldat%ITPRINT(i),i=1,NTPRT)
        endif
    endif

  3       format (51i12)
    write(luin, '(i12)') NINTG
    write(luin, 3) timestepdat%MAXSTP(1:NINTG)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
    write(luin,*) nElemsWithSlipperyNodes
    write(luin, '(50I12)') modeldatactx%SlipElemLookup
#endif

    close(luin)

endif


if (iecho.eq.8) then
    call ParallelLog("writemodeldata", "written timestuff to luin")
endif

modelpartname = 'modeldata.dat.'

rnk = getrank()
logsize = log10(dble(getsize())) + 1

do i=1,logsize
    inneri = 10**(logsize-i)
    innerj = rnk/inneri
    modelpartname = modelpartname(1:lnblk(modelpartname))//char(48+innerj)
    rnk = rnk - innerj * inneri
enddo

modelpartname = trim(outdir)//'/'//trim(modelpartname)

call openf(luin,modelpartname,'unknown')

if (iecho.eq.8) then
    call ParallelLog("writemodeldata", "writes indices to luin")
endif

write(frmt,4) NEN
 4    format('(',I1,'I12)')

write(luin, '(I12)') meshdat%Nvlocal
write(luin, 3) (meshdat%gloInd(i),i=1,meshdat%Nvlocal)
write(luin, '(I12)') meshdat%Nelocal
write(luin, 3) (meshdat%gloElt(i),i=1,meshdat%Nelocal)

if (iecho.eq.8) then
    call ParallelLog("writemodeldata", "writes meshdat IEN to luin")
endif

do i=1,meshdat%Nelocal
    write(luin, frmt) (meshdat%IEN(j,i),j=1,NEN)
enddo

if (iecho.eq.8) then
    call ParallelLog("writemodeldata", "checking surface nodes")
endif


if (NSURFglobal.gt.0) then
    write(luin, '(I12)') NSURF
    write((luin),3) (modeldatactx%ISURFORDER(i),i=1,NSURF)
endif


#ifdef EXPERIMENTAL_ReducedSlipAlloc
write(luin, *) nLocalElemsWithSLipperyNodes
#endif

close(luin)

if (iecho.eq.8) then
    call ParallelLog("writemodeldata", "leaving")
endif

end subroutine writemodeldata
!***************************************************************


   subroutine cleanAOdata()
  integer :: ierr
  ! remove all vecs. They are global throughout
  ! the module.

  call VecDestroy(ndofVec, ierr)
  call VecDestroy(nsdVec, ierr)
  call VecDestroy(oneVec, ierr)
  call VecDestroy(twoVec, ierr)
  call VecDestroy(ndofVecLoc, ierr)
  call VecDestroy(nsdVecLoc, ierr)
  call VecDestroy(oneVecLoc, ierr)
  call VecDestroy(twoVecLoc, ierr)

  call VecScatterDestroy(ndofScatter, ierr)
  call VecScatterDestroy(nsdScatter, ierr)
  call VecScatterDestroy(oneScatter, ierr)
  call VecScatterDestroy(twoScatter, ierr)

  call ISLocalToGlobalMappingDestroy(isl2g, ierr)
  call ISLocalToGlobalMappingDestroy(isl2gt, ierr)


   end subroutine


#endif  
! endif sparse, over the entire module
END MODULE AOMODULE
