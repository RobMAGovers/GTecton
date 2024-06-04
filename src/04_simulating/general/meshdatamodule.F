
MODULE MESHDATAMODULE

USE FILEMODULE
use modeltopology, only: nSidesPerElem
use petscksp

implicit none

#define PETSC_AVOID_DECLARATIONS

!#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"

!#include "petsc/finclude/petscao.h"
!#include "petsc/finclude/petscis.h"

#undef PETSC_AVOID_DECLARATIONS
   TYPE meshdata

!         global quantity of vertices and elements
PetscInt :: Nvglobal!, NUMNP
PetscInt :: Neglobal!, NUMEL
!         local quantity of vertices and elements
PetscInt :: Nvlocal
PetscInt :: Nelocal



!integer, allocatable :: NvlocalOfAllPartitions(:)
!integer, allocatable :: NelocalOfAllPartitions(:)

!         offset per elements to refer to the local numbering 
!         vector. So, example, 
!         partition        Nelocal        offset
!         0                3            0
!          1                5            4
!          2                9            9
!          3                2            18

integer, allocatable :: elementOffset(:)

!         total quantity of neighbour vertices of this CPU
  PetscInt :: NGhostPoints

!         total quantity of neighbors of each point in the partition.
!         vertices can be counted multiple times in here.
  PetscInt :: Nvneighborstotal

!         total quantity of vertices that belong to elements
!         of this CPU, but are not neighbors of vertices of 
!         this CPU
  PetscInt :: NvRest

!         total number of relevant vertices for this CPU.
!         That includes:
!         1: its own vertices (nvlocal)
!         2: its neighbouring vertices (Nvneighborstotal)
!         3: other vertices of elements in this partition
  PetscInt :: NvAllRelevant

!         for each vertex the local indices of its neighbours
!         AdjM(i,j) is the local index of the i-th neighbour of vertex j.
!         Indices start at 0
  integer :: nMaxNeighbors
  PetscInt, ALLOCATABLE :: AdjM(:,:)


!         for each vertex the quantity of neighbours
  PetscInt, ALLOCATABLE :: itot(:)

!         for each vertex its ID, needed for plt3d
  PetscInt, ALLOCATABLE :: vertexID(:)


  PetscInt, ALLOCATABLE :: pointTags(:,:)
!         for each line, its start and end point
  PetscInt, ALLOCATABLE :: lineEndPointIDs(:,:)
!         for each line, its marker
  PetscInt, ALLOCATABLE :: lineTags(:,:)

! all actual line segments (not just labeled ones)
  PetscInt :: nLineSegments
  PetscInt, ALLOCATABLE :: lineSegments(:,:)
  PetscInt, ALLOCATABLE :: linesInElements(:,:)

  PetscInt              :: nTriangles
  PetscInt, ALLOCATABLE :: triangleConnectivity(:,:)
!         for each triangle, its marker
  PetscInt, ALLOCATABLE :: triangleTags(:,:)
!         for each triangle, its boundaring elements
  PetscInt, ALLOCATABLE :: triangleElements(:,:)
  PetscInt, ALLOCATABLE :: triangleElementSide  ! side 1, 2, 3 or 4

  PetscInt, ALLOCATABLE :: tetrahedronTags(:,:)

!         for each element its vertex indices (ien)
  PetscInt, ALLOCATABLE :: ien(:,:)
  PetscInt, ALLOCATABLE :: globalIen(:,:)

!         for each element its material number
  PetscInt, ALLOCATABLE :: mat(:)

!         rank for each vertex (global) (vertex 2 partition)
  PetscInt, ALLOCATABLE :: v2p(:)
!         rank for each element(global) (element 2 partition)
  PetscInt, ALLOCATABLE :: e2p(:)

!         global vertex numbers (for local vertices only)
  PetscInt, ALLOCATABLE :: gloInd(:)
!         global element numbers (for local elems only)
  PetscInt, ALLOCATABLE :: gloElt(:)
!         local vertex numbers (for local vertices only)
  PetscInt, ALLOCATABLE :: locInd(:)
!         local element numbers (for local elems only)
  PetscInt, ALLOCATABLE :: locElt(:)

!         markers of the nodal points
#ifndef SPARSE
  ! only needed in plt3d
  PetscInt, ALLOCATABLE :: globalMarker(:)
#endif

!         For each element its neighbors, 
!         used by the sample command in the postprocessing
  PetscInt, ALLOCATABLE :: ElementNeighbours(:,:)
!         maximum number of neighbours
!         determined during reading of the tecin.dat.partf.nps
  PetscInt              :: maxNeighbourCount

!         look-up-table to translate global numbers to
!         local numbers, so that we do not have to scan 
!         the gloInd array every time we need to do this.
!         Entry k, will contain the local index
!         of the vertex/element with global index k
!         If this element is not in the partition,
!         it will contain -1.
!         One for the Elements...
  PetscInt, ALLOCATABLE :: Glo2LocElement(:)
  PetscInt, ALLOCATABLE :: FullGlo2LocElement(:)

!         ... and one for the vertices
  PetscInt, ALLOCATABLE :: Glo2LocVertex(:)


!         coordinates
  PetscReal, ALLOCATABLE :: X(:,:) !coords(:,:)

!          PetscReal, ALLOCATABLE :: localX(:,:)
!         local coordinates
!          Vec :: localcoords
!          Vec :: localdisps
!          Vec :: localtemps
!          VecScatter scatter
!          VecScatter scatterd
!          VecScatter scattert


! hold the local elementID for each neighbor of each element
! neighboringElemID(i,j) will be the ID of the element bordering on side i of element j
integer, allocatable :: neighboringElemID(:,:)

! hold the partition numbers of each beighbor of each element
! neighboringElemPartition(i,j) will be the ID of the partition of 
! the element bordering on side i of element j
integer, allocatable :: neighboringElemPartition(:,:)

! store the number of ghost elements in each partition which are relevant for this partition
integer, allocatable :: hasNElemInPartition(:)

integer, allocatable :: NvlocalOfAllPartitions(:)
integer, allocatable :: NelocalOfAllPartitions(:)

integer :: lowestLocalElemID
integer :: lowestLocalVertexID


   END TYPE meshdata


#ifdef EXPERIMENTAL_reducedNbrsMemory
type neighborlist
    integer, allocatable(neighbors)
end type
#endif



   integer, parameter :: plotDOFperElem = 3 
   ! number of degrees of freedom of an element
   ! for constant functions, this is 1,
   ! for linear functions on triangles, it is 3
   ! this will be changeble in the future, 
   ! such as for higher order elements and differently shaped elements,
   ! but for now it is just 3.

   integer, parameter :: NPlotEntriesPerElem = 3
   ! for example, effective stress, effective strain 
   ! all have 3 data entries per element


   type (meshdata) :: meshdatactx


!#include "petsc/finclude/petscvec.h"
!      AO aonp, aoel, aoeq
!      IS isglobal, islocal
   Vec coords
   Vec disps
   Vec ttemps
!      PetscInt, ALLOCATABLE :: vertices(:),
!     .       verticesmask(:)
!      PetscInt nvertices
!      PRIVATE :: writeAOmeshdata, writeAOmapping


logical :: elementNeighborTableComplete

contains


!------------------------------------------------------------------------------

logical function vertexOccursInElement(vertexID, elemID, nNodesPerElement)

! helper function for initializeSample (right above)

implicit none

integer :: vertexID, elemID
integer :: i
integer :: nNodesPerElement

vertexOccursInElement = .false.

do i=1,nNodesPerElement
    if (meshdatactx%IEN(i, elemID) .eq. vertexID) then
         vertexOccursInElement = .true.
         return
    endif
enddo

return

end function

!------------------------------------------------------------------------------
subroutine setnumbers(meshdat,numnpglo,numnploc)

implicit none

type (meshdata) :: meshdat

PetscInt :: numnpglo,numnploc

meshdat%Nvglobal = numnpglo
meshdat%Nvlocal = numnploc

end subroutine setnumbers
!------------------------------------------------------------------------------

subroutine allocatevertexdata(meshdat, maxnbrs, nsd)

implicit none

type (meshdata) :: meshdat
PetscInt  ::  maxnbrs, nsd, nvglobal, nvlocal
PetscInt ::  ERROR

nvglobal = meshdat%Nvglobal
nvlocal = meshdat%Nvlocal

if (allocated(meshdat%AdjM)) then
    deallocate(meshdat%AdjM)
endif
ALLOCATE(meshdat%AdjM(maxnbrs,nvlocal), STAT=ERROR)

if (allocated(meshdat%itot)) then
    deallocate(meshdat%itot)
endif
ALLOCATE(meshdat%itot(nvlocal), STAT=ERROR)

if (allocated(meshdat%v2p)) then
    deallocate(meshdat%v2p)
endif
ALLOCATE(meshdat%v2p(nvglobal), STAT=ERROR) ! only during readfe

if (allocated(meshdat%gloInd)) then
    deallocate(meshdat%gloInd)
endif
ALLOCATE(meshdat%gloInd(nvlocal), STAT=ERROR)

#ifndef SPARSE
if (allocated(meshdat%globalMarker)) then
    deallocate(meshdat%globalMarker)
endif
ALLOCATE(meshdat%globalMarker(nvglobal), STAT=ERROR)
#endif

if (allocated(meshdat%locInd)) then
    deallocate(meshdat%locInd)
endif
ALLOCATE(meshdat%locInd(nvlocal), STAT=ERROR)

if (allocated(meshdat%X)) then
    deallocate(meshdat%X)
endif
ALLOCATE(meshdat%X(nsd,nvlocal), STAT=ERROR)

meshdat%AdjM   = 0
meshdat%itot   = 0
meshdat%v2p    = 0
meshdat%gloInd = 0
meshdat%locInd = 0
meshdat%X      = 0d0

if (ERROR.ne.0) then
    write(*,*) 'ERROR during allocation of vertex data'
endif

end subroutine allocatevertexdata

!------------------------------------------------------------------------------

subroutine allocateelementdata(meshdat, nen)

use modelctx,      only: getrank, getsize
use modeltopology, only: nSidesPerElem

implicit none

#ifdef SPARSE
integer, allocatable :: allNeLocal(:)
integer              :: i
#endif

type (meshdata) :: meshdat
PetscInt :: nen, neglobal, nelocal
PetscInt :: ERROR

neglobal = meshdat%Neglobal
nelocal = meshdat%Nelocal

if (allocated(meshdat%ien))then
    deallocate(meshdat%ien)
endif
ALLOCATE(meshdat%ien(nen,nelocal), STAT=ERROR)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
! TODO
!******* this is veeeery bad temporary thing made for the ****
!******* setting up op the LMX. ******************************
!******* Eventually, the required global indices should ******
!******* be obtained by communication from the other threads *
if (allocated(meshdat%globalIen))then
    deallocate(meshdat%globalIen)
endif
ALLOCATE(meshdat%globalIen(nen,neglobal), STAT=ERROR)  
!*************************************************************
#endif

if (allocated(meshdat%mat))then
    deallocate(meshdat%mat)
endif 
ALLOCATE(meshdat%mat(nelocal), STAT=ERROR)

if (allocated(meshdat%e2p))then
    deallocate(meshdat%e2p)
endif
ALLOCATE(meshdat%e2p(neglobal), STAT=ERROR) ! only during readfe

if (allocated(meshdat%gloElt))then
    deallocate(meshdat%gloElt)
endif
ALLOCATE(meshdat%gloElt(nelocal), STAT=ERROR)

if (allocated(meshdat%locElt))then
    deallocate(meshdat%locElt)
endif
ALLOCATE(meshdat%locElt(nelocal), STAT=ERROR)

! ghost element management
if (allocated(meshdat%neighboringElemID))then
    deallocate(meshdat%neighboringElemID)
endif
ALLOCATE(meshdat%neighboringElemID(nSidesPerElem,nelocal), STAT=ERROR)

if (allocated(meshdat%neighboringElemPartition))then
    deallocate(meshdat%neighboringElemPartition)
endif
ALLOCATE(meshdat%neighboringElemPartition(nSidesPerElem,nelocal), STAT=ERROR)

if (allocated(meshdat%hasNElemInPartition))then
    deallocate(meshdat%hasNElemInPartition)
endif

#ifdef SPARSE
ALLOCATE(meshdat%hasNElemInPartition(getsize()), STAT=ERROR)
allocate(meshdat%elementOffset(getsize()), STAT=ERROR)
allocate(allNeLocal(getsize()))

! gathering all the Nelocals
call mpi_allgather(meshdatactx%Nelocal, 1, MPI_int, allNeLocal, 1, mpi_int, MPI_COMM_WORLD, ERROR)

! sum them to get offsets of local numbering to this partition
meshdat%elementOffset(1) = 0
do i=2,getsize()
    meshdat%elementOffset(i) = meshdat%elementOffset(i-1) + allNeLocal(i-1)
enddo

deallocate(allNeLocal)
#endif


if (ERROR.ne.0) then
    stop "Error while allocating the element data"
endif

end subroutine allocateelementdata

!------------------------------------------------------------------------------

subroutine deallocatevertexdata(meshdat)

implicit none

type (meshdata) :: meshdat
PetscInt ERROR

DEALLOCATE(meshdat%AdjM, STAT=ERROR)
DEALLOCATE(meshdat%itot, STAT=ERROR)
DEALLOCATE(meshdat%v2p, STAT=ERROR)
DEALLOCATE(meshdat%gloInd, STAT=ERROR)
DEALLOCATE(meshdat%locInd, STAT=ERROR)
DEALLOCATE(meshdat%X, STAT=ERROR)

if (ERROR.ne.0) then
    stop "Error while deallocating the vertex data"
endif

end subroutine deallocatevertexdata

!------------------------------------------------------------------------------

subroutine deallocateelementdata(meshdat)

implicit none

type (meshdata) :: meshdat
PetscInt ERROR

DEALLOCATE(meshdat%ien, STAT=ERROR)
DEALLOCATE(meshdat%mat, STAT=ERROR)
DEALLOCATE(meshdat%e2p, STAT=ERROR)
DEALLOCATE(meshdat%gloElt, STAT=ERROR)
DEALLOCATE(meshdat%locElt, STAT=ERROR)

DEALLOCATE(meshdat%neighboringElemID, STAT=ERROR)
DEALLOCATE(meshdat%neighboringElemPartition, STAT=ERROR)
DEALLOCATE(meshdat%hasNElemInPartition, STAT=ERROR)

if (ERROR.ne.0) then
    stop "Error while deallocating the element data"
endif

end subroutine deallocateelementdata
!------------------------------------------------------------------------------

subroutine coordsMinMax (NSD,CRDMIN,CRDMAX)

use constants, only: dfmax

! Routine for determining bounds on coordinates data in array meshdatactx%X

implicit none
!-pass
integer          :: NSD
double precision :: CRDMIN(NSD),CRDMAX(NSD)
integer          :: i,j

do i=1,NSD
   CRDMIN(i) =  dfmax
   CRDMAX(i) = -dfmax
enddo

do j=1,meshdatactx%nvglobal
    do i=1,NSD
        if (meshdatactx%X(i,j).gt.CRDMAX(i)) then
            CRDMAX(i) = meshdatactx%X(i,j)
        endif
        if (meshdatactx%X(i,j).lt.CRDMIN(i)) then
            CRDMIN(i) = meshdatactx%X(i,j)
        endif
    enddo
enddo

return

end subroutine
!-------------------------------------------------------------------------------



! both print subroutines could be useful for debugging, but they serve no other purpose.
!----------------------------------------------------------------------------
!      subroutine printcoordinates(vector, nsd, length,FILE_outputf)
!
!
!      implicit none
!      integer nsd, length,i,j,FILE_outputf
!      double precision, dimension(nsd,*) :: vector
!
!      do i=1,length
!          write(FILE_outputf,*) i, (vector(j,i),j=1,NSD)
!      enddo
!
!      end subroutine printcoordinates

!#ifdef SPARSE
!      subroutine printpetsccoordinates(vector, nsd, length,FILE_outputf)
!      implicit none
!!#include "petsc/finclude/petsc.h"
!!#include "petsc/finclude/petscsysdef.h"
!!#include "petsc/finclude/petscvecdef.h"
!!#include "petsc/finclude/petscvec.h"
!      integer nsd, length,i,j,FILE_outputf
!      Vec :: vector
!      PetscScalar, pointer ::  tmpcoords(:)
!      PetscErrorCode ierr
!      PetscInt vecsize
!
!      call VecGetSize(vector, vecsize, ierr)
!      call VecGetArrayF90(vector, tmpcoords,ierr)
!      do i=1,length
!          write(FILE_outputf,*) i,(tmpcoords(nsd*(i-1)+j),j=1,nsd)
!      enddo
!      call VecRestoreArrayF90(vector, tmpcoords,ierr)
!
!      end subroutine printpetsccoordinates
!#endif



!-----------------------------------------------------------------

! below here are two routines that find the neighbors of each of the elements.
! They are computed globally, to not suffer the delays from the partition boundaries.
! At the time, the solver is not yet initialised, and we have lots of memory available.

! After the computation of the neighbors, they are divided over the elements,
! and the numbering schemes are switched from glocal to local.


! They are separate in 2D and 3D, because of the element conventions.

subroutine BuildElementNeighbourTable2D
! this subroutine should be called before 
! any sampling is done.

! it creates a look-up table so that for all 
! elements it neighbours are known.

implicit none

integer, allocatable :: ElementsPerVertex(:,:)
integer, allocatable :: ElementCount(:)

integer :: i,j,k
integer :: thisVertex, thisElem
integer :: BaseVertex
integer :: v1, v2, v3

integer, parameter   :: maxElementsPerVertex = 30



allocate(meshdatactx%ElementNeighbours(meshdatactx%neglobal,3))

meshdatactx%ElementNeighbours = 0

! for every side of each element, the adjacent element
! will be stored, of 0 if the surface is on the domain boundary

! first we will create a lookup table by scanning through the adjacency
! to, for every point, store the elements adjacent to it


allocate(ElementsPerVertex(meshdatactx%nvglobal, maxElementsPerVertex))
allocate(ElementCount(meshdatactx%nvglobal))

ElementCount = 0

do i=1,meshdatactx%neglobal
    do j=1,3 ! loop over vertices ! for triangle
        thisVertex = meshdatactx%ien(j,i)
        ! increase the count
        ElementCount(thisVertex) = ElementCount(thisVertex) + 1
        ! and add it to the lookup table
        if (ElementCount(thisVertex).gt.maxElementsPerVertex) then
            write(*,*) "Vertex", thisVertex, "has more than a ", maxElementsPerVertex, &
                       " elements adjacent to it."
            write(*,*) "This can be tolerated by increasing allocation in subroutine "
            write(*,*) "BuildElementNeighbourTable2D in file TempSrc/meshdatamodule.F,"
            write(*,*) "but it also raises question marks to the quality of mesh."
            STOP "Leaving GTecton"
        endif

        ElementsPerVertex(thisVertex, ElementCount(thisVertex)) = i
    enddo
enddo

! using that lookup table, build a neighbour table
do i=1,meshdatactx%neglobal
    ! for element i

    ! abbreviate its vertices 
    v1 = meshdatactx%IEN(1,i)
    v2 = meshdatactx%IEN(2,i)
    v3 = meshdatactx%IEN(3,i)

    ! ***** find the element bordering on side 1 (vertex 1 and 2) *******
    BaseVertex = meshdatactx%IEN(1,i)

    ! must have vertex IEN(1,i) in common.
    do k=1,ElementCount(BaseVertex)
        ! loop through elements that have this vertex
        ! until the element is found that has also 
        ! vertex 2 and 4, but *not* 3, because then it would
        ! be the same vertex

        thisElem = ElementsPerVertex(BaseVertex,k)

        if (vertexOccursInElement(v2, thisElem, 3) .and. &
      .not. vertexOccursInElement(v3, thisElem, 3)) then

            ! found neighboring element
            meshdatactx%ElementNeighbours(i,1) = thisElem
        endif

    enddo

    ! ***** side 2; vertices 2 and 3 ***********
    BaseVertex = meshdatactx%IEN(2,i)

    do k=1,ElementCount(BaseVertex)
        ! check which element in ElementsPerVertex(i,:)
        ! has vertex 1, 3 and 2nd the same vertices as element i

        thisElem = ElementsPerVertex(BaseVertex,k)

        if (vertexOccursInElement(v3, thisElem, 3) .and. &
      .not. vertexOccursInElement(v1, thisElem, 3)) then

            ! found neighboring element 
            meshdatactx%ElementNeighbours(i,2) = thisElem
        endif

    enddo

    ! ***** side 3; vertices 1 and 3 ***********
    BaseVertex = meshdatactx%IEN(1,i)

    do k=1,ElementCount(BaseVertex)
        ! check which element in ElementsPerVertex(i,:)
        ! has vertex 1, 4 and 3th same vertices as element i
        thisElem = ElementsPerVertex(BaseVertex,k)

        if (vertexOccursInElement(v3, thisElem, 3) .and. &
      .not. vertexOccursInElement(v2, thisElem, 3)) then
            ! found neighboring element 
            meshdatactx%ElementNeighbours(i,3) = thisElem
        endif

    enddo

enddo

! and remove the original look-up
deallocate(ElementsPerVertex)
deallocate(ElementCount)


elementNeighborTableComplete = .true.


end subroutine
!-----------------------------------------------------------------

subroutine BuildElementNeighbourTable3D

! this subroutine should be called before 
! any sampling is done.

! it creates a look-up table so that for all 
! elements it neighbours are known.


implicit none

integer, allocatable :: ElementsPerVertex(:,:)
integer, allocatable :: ElementCount(:)

integer :: i,j,k
integer :: thisVertex, thisElem
integer :: BaseVertex
integer :: v1, v2, v3, v4

integer, parameter   :: maxElementsPerVertex=150

allocate(meshdatactx%ElementNeighbours(meshdatactx%neglobal,4))

meshdatactx%ElementNeighbours = 0

! for every side of each element, the adjacent element
! will be stored, of 0 if the surface is on the domain boundary

! first we will create a lookup table by scanning through the adjacency
! to, for every point, store the elements adjacent to it


allocate(ElementsPerVertex(meshdatactx%nvglobal, maxElementsPerVertex))
allocate(ElementCount(meshdatactx%nvglobal))

ElementCount = 0

do i=1,meshdatactx%neglobal
    do j=1,4 ! loop over vertices ! for tetrahedrons. Adapt later to other elements 
        thisVertex = meshdatactx%IEN(j,i)
        ! increase the count
        ElementCount(thisVertex) = ElementCount(thisVertex) + 1
        ! and add it to the lookup table
        if (ElementCount(thisVertex).gt.maxElementsPerVertex) then
            write(*,*) "Vertex", thisVertex, "has more than a ",maxElementsPerVertex, &
                       " elements adjacent to it."
            write(*,*) "This can be tolerated by increasing allocation in subroutine "
            write(*,*) "BuildElementNeighbourTable3D in file TempSrc/meshdatamodule.F," 
            write(*,*) "but it also raises question marks to the quality of mesh."
            STOP "Leaving GTecton"
        endif
        ElementsPerVertex(thisVertex, ElementCount(thisVertex)) = i
    enddo
enddo

! using that lookup table, build a neighbour table
do i=1,meshdatactx%neglobal
    ! for element i

    ! abbreviate its vertices 
    v1 = meshdatactx%IEN(1,i)
    v2 = meshdatactx%IEN(2,i)
    v3 = meshdatactx%IEN(3,i)
    v4 = meshdatactx%IEN(4,i)

    ! ***** find the element bordering on side 1 *******
    BaseVertex = meshdatactx%IEN(1,i)

    ! must have vertex IEN(1,i) in common.


    do k=1,ElementCount(BaseVertex)
        ! loop through elements that have this vertex
        ! until the element is found that has also 
        ! vertex 2 and 4, but *not* 3, because then it would
        ! be the same vertex

        thisElem = ElementsPerVertex(BaseVertex,k)

        if (vertexOccursInElement(v2, thisElem, 4) .and. &
            vertexOccursInElement(v4, thisElem, 4) .and. &
      .not. vertexOccursInElement(v3, thisElem, 4)) then

            ! found neighboring element
            meshdatactx%ElementNeighbours(i,1) = thisElem
        endif

    enddo



    ! ***** side 2; vertices 1, 3 and 2 ***********
    BaseVertex = meshdatactx%IEN(1,i)

    do k=1,ElementCount(BaseVertex)
        ! check which element in ElementsPerVertex(i,:)
        ! has vertex 1, 3 and 2nd the same vertices as element i

        thisElem = ElementsPerVertex(BaseVertex,k)

        if (vertexOccursInElement(v3, thisElem, 4) .and. &
            vertexOccursInElement(v2, thisElem, 4) .and. &
      .not. vertexOccursInElement(v4, thisElem, 4)) then
            ! found neighboring element 
            meshdatactx%ElementNeighbours(i,2) = thisElem
       endif

    enddo

    ! ***** side 3; vertices 1, 4 and 3 ***********
    BaseVertex = meshdatactx%IEN(1,i)

    do k=1,ElementCount(BaseVertex)
        ! check which element in ElementsPerVertex(i,:)
        ! has vertex 1, 4 and 3th same vertices as element i
        thisElem = ElementsPerVertex(BaseVertex,k)

        if (vertexOccursInElement(v4, thisElem, 4) .and. &
            vertexOccursInElement(v3, thisElem, 4) .and. &
      .not. vertexOccursInElement(v2, thisElem, 4)) then
            ! found neighboring element 
            meshdatactx%ElementNeighbours(i,3) = thisElem
        endif

    enddo



    ! ***** side 4; vertices 2, 3 and 4 ***********
    BaseVertex = meshdatactx%IEN(2,i)

    do k=1,ElementCount(BaseVertex)
        ! check which element in ElementsPerVertex(i,:)
        ! has vertex 2, 3 and 4th same vertices as element i

        thisElem = ElementsPerVertex(BaseVertex,k)

        if (vertexOccursInElement(v3, thisElem, 4) .and. &
            vertexOccursInElement(v4, thisElem, 4) .and. &
      .not. vertexOccursInElement(v1, thisElem, 4)) then
            ! found neighboring element 
            meshdatactx%ElementNeighbours(i,4) = thisElem
        endif

    enddo

enddo


! and remove the original look-up
deallocate(ElementsPerVertex)
deallocate(ElementCount)


elementNeighborTableComplete = .true.


end subroutine
!-----------------------------------------------------------------


END MODULE MESHDATAMODULE
