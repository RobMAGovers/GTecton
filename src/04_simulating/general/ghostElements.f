!------------------------------------------------------------------------------

subroutine getGhostElements()

use modelctx,    only: getrank, &
                       getsize
!use debugmodule,   only: delay
use modeltopology, only: NDOF

#ifdef SPARSE
use AOmodule,    only: vertices, &
                       verticesmask
#endif

implicit none

! we have all the mesh information from reading the file.
! Using that, we can find the ghost elements.

! This is done in three steps.

! 1: list all the points in elements belonging to this partition.
!    Note that this list can differ from the points in
!    this partition. There are no points here that are not in elements
!    of this partition. And there can be points in it that are not in 
!    this partition.

! 2: for each of those points, list all the elements boundaring on it


! 2b: for each of those points, get the maximum number of neighbors from 
!     the partition that ows the point. itot, containing this number,
!     only present on the owning partition.

! 3: for each of those points, compare the number of matching elements
!    with the number of neighbouring points

! 5: build a list of points that have elements missing. Those are the
!    elements in other partitions, which we are looking for.

! 6: exchange this list with other processors, and all threads
!    will see this list.

! 7: All processors will assemble a list of element IDs they have 
!    to send to other processors

! 8: All processors will receive the lists from the others, and have
!    the IDs of the ghost elements, with a minimum of ghost memory used.



integer :: nPointsInElements, maxNeighbors
integer :: highestPoint, thisPoint, nbrs

integer :: iElem, iPoint, iThread, iCorner, jPoint,i ! loop indices

integer, allocatable :: nElemsPerVertex(:)
integer, allocatable :: ElemIDPervertex(:,:)

logical, allocatable :: elemIDdifferentNeighbors(:)
integer :: nMissingNeighbors

integer :: thisRank

logical, allocatable :: knowPoints(:)

integer, allocatable :: pointsInElements(:)

integer :: nForeignPoints
integer, allocatable :: foreignPoints(:)

integer :: smallestLocInd  ! of nodes belonging this partition
integer :: highestLocInd   ! 

integer, allocatable :: nForeignPointsOfAll(:)
integer              :: nAllForeignNodes
integer, allocatable :: foreignPointsOfAll(:)

integer, allocatable :: offsets(:)

integer, allocatable :: nNeighbors(:)

integer :: error
integer, allocatable :: nNeighborsFull(:)
integer, allocatable :: nNeighborsMine(:)

integer, allocatable :: nNeighborsOfPointsInMyElems(:)

integer              :: nIncomplete           ! number of points with not all neighbors known
integer, allocatable :: nIncompletePerProc(:) ! same, but for all threads
integer              :: nAllIncomplete        ! number of all incomplete nodes over all partitions
integer, allocatable :: incompleteNodes(:)    ! nodes that do not have all elements
                                              ! attached to it in this partition.
integer, allocatable :: allIncompletePoints(:)

integer, allocatable :: nAllNeighbors(:)

integer              :: nOwnPointsInOwnElems
integer              :: nodesPerElement

integer, allocatable :: pointsOfAllBorderElements(:)
integer              :: nPointsPerElement

integer :: smallestPointIDInMyElems  ! of nodes in IEN, as such relevant for this 
integer :: highestPointIDInMyElems   ! partition, even though the node itself may
                                     ! belong to another partition

integer, allocatable :: nElemsPerPoint(:)    ! for every ghost point, stores the number of 
                                             ! elements that this partition has, bordering on it.
                                             ! Only in this partition.

integer, allocatable :: nAllElemsPerPoint(:) ! for every ghost point, stores the number of 
                                             ! elements that this partition has, bordering on it.
                                             ! Summed over all processors.

integer, allocatable :: elemIDsToThisPoint(:,:) ! the matching element IDs

!integer, allocatable :: pointsOfAllBorderElements(:,:,:) ! all border elements


integer :: nMaxElemsPerPoint ! maximum to allocate enough room to store all elems.
integer :: MostNeighborsPerPoint ! over all partitions


integer, allocatable :: partitionContributions(:,:)
integer, allocatable :: otherPartitionsHaveCount(:,:)

integer, allocatable :: MyNeighborPartitionElemCount(:,:)

integer, allocatable :: receiveElements(:)
integer, allocatable :: receiveElemVertices(:,:)
integer :: nReceiveElems
integer, allocatable :: receiveOffsets(:)


integer :: nSendElements
integer :: nSendNodes
integer :: nSendElementsInThisProc
integer, allocatable :: nSendElementsToEachThread(:)
integer, allocatable :: sendElements(:)
integer, allocatable :: sendElementsOffsets(:)
integer, allocatable :: sendNodes(:)

integer :: nReceiveElements
integer, allocatable :: nReceiveElementsFromEachThread(:)

integer :: nSharedIncomplete ! incomplete points shared with an arbitrary other partition

integer, allocatable :: SharedIncompletePointWithPartIthread(:)

integer :: checkpoint
integer :: nSendNdes
logical :: pointfound
integer :: searchID

integer :: startIndex, endIndex
integer :: startOffsetIndex, endOffsetIndex
integer :: iDestination

integer, allocatable :: eachRankReceivesNElems(:)
integer, allocatable :: receiveoffsetofeachthread(:)

integer, allocatable :: vertex1Send(:), vertex2Send(:), vertex3Send(:), vertex4Send(:)
integer, allocatable :: vertex1Receive(:), vertex2Receive(:), vertex3Receive(:), vertex4Receive(:)

integer, allocatable :: elementNumberingOffset(:), neLocalOfEachThread(:)


thisRank = getrank()

!write(*,*) "Entered getGhostElements"

!if (thisRank.eq.3) then
!write(*,*) "rank",getrank()," has full glo 2 loc ", meshdatactx%FullGlo2LocElement
!write(*,*) "rank",getrank()," has glo elt: ", meshdatactx%gloElt


!write(*,*) "rank",getrank()," has loc elt: ", meshdatactx%locElt

!write(*,*) "rank",getrank()," has glo ind: ", meshdatactx%gloind
!write(*,*) "rank",getrank()," has loc ind: ", meshdatactx%locind
!endif

!write(*,*) "rank",getrank()," has ien: ", meshdatactx%ien


! 1 a determine the number of points we need to allocate space for,
!     The vertex indices in LocalIen all start at 1, so we simple need
!     the highest value.

highestPoint = 0
nPointsInElements = meshdatactx%nvglobal
do iElem = 1,meshdatactx%nelocal
    do iPoint = 1,3
        thisPoint = meshdatactx%ien(iPoint, iElem)
        if (thisPoint .gt. highestPoint) then
            highestPoint = thisPoint
        endif
    enddo
enddo

! 1 b  Note that IEN may not contain all the points of the vertices
!      array. There are neighboring points that are not needed here.
!      To determine the proper number, count the number of different
!      values in IEN.

allocate(knowPoints(highestPoint))
knowPoints = .false.

do iElem = 1,meshdatactx%nelocal
    do iPoint = 1,3
        thisPoint = meshdatactx%ien(iPoint, iElem)
        knowPoints(thisPoint) = .true.
    enddo
enddo

!if (thisRank.eq.3) then
!write(*,*) "rank",getrank()," has knowPoints: ", knowPoints
!#ifdef SPARSE
!write(*,*) "rank",getrank()," has vertices: ", vertices
!#endif
!endif

! count the trues
nPointsInElements = 0

do iPoint =1,highestPoint
    if (knowPoints(iPoint)) then
        nPointsInElements = nPointsInElements + 1
    endif
enddo


if (thisRank.eq.3) then
!write(*,*) "rank",getrank(),"Allocating ", nPointsInElements
endif

! and gather their IDs in array

allocate(pointsInElements(nPointsInElements))

nPointsInElements = 0

do iPoint = 1, highestPoint
!    write(*,*) "check", iPoint, "of", knowPoints(iPoint)
    if (knowPoints(iPoint)) then
        nPointsInElements = nPointsInElements + 1
!        if (thisRank.eq.3) then
!            write(*,*) "adding ", iPoint, " to pos ", nPointsInElements, "in pie"
!        endif
#ifdef SPARSE
	pointsInElements(nPointsInElements) = vertices(iPoint)
#else
        ! when running on core, ipoint and vertices(ipoint) would be the same
        ! but vertices in fact does not exist, then.
        pointsInElements(nPointsInElements) = iPoint
#endif
    endif
enddo

deallocate(knowPoints)

if (thisRank.eq.3) then
!write(*,*) "rank",getrank(),"has own points ", pointsInElements
endif

! determine the number of foreign points
! by checking whether a pointInElements entry
! also occurs in LocInd
smallestLocInd = meshdatactx%locInd(1)
highestLocInd = meshdatactx%locInd(meshdatactx%nvlocal)

nForeignPoints = 0

if (thisRank.eq.3) then
!write(*,*) "rank",getrank(),"has low, high ", smallestLocInd, highestLocInd
endif



do iPoint = 1,nPointsInElements
    thispoint = pointsInElements(iPoint)
    if ((thisPoint .gt. highestLocInd) .or. &
        (thisPoint .lt. smallestLocInd)) then
        nForeignPoints = nForeignPoints + 1
    endif
enddo

if (thisRank.eq.3) then
!write(*,*) "rank",getrank(),"has n foreign points ", nForeignPoints
endif



allocate(foreignPoints(nForeignPoints))

! fill it
nForeignPoints = 0



do iPoint = 1,nPointsInElements
    thispoint = pointsInElements(iPoint)

    if ((thisPoint .gt. highestLocInd) .or. &
        (thisPoint .lt. smallestLocInd)) then
        nForeignPoints = nForeignPoints + 1
        foreignPoints(nForeignPoints) = thisPoint
    endif
enddo

!if (thisRank.eq.3) then
!write(*,*) "rank",getrank(),"has foreign points ", foreignPoints
!endif

! Now that we know the number of foreign nodes that this processor
! has, we still have no idea to which processor it belongs.
! What we will do is construct a global list with all the foreign
! nodes.

! The cores will share the vector, and each core will select the
! points it owns, and sets the numbr of neighbors in it.


allocate(nForeignPointsOfAll(getsize()))

call mpi_allgather(nForeignPoints, 1, MPI_int, nForeignPointsOfAll, 1, mpi_int, MPI_COMM_WORLD, error)
if (error.ne.0) then
    stop "mpi_allgather failed with error code"
endif

if (thisRank.eq.3) then
!write(*,*) "rank",getrank(),"has foreign points of each: ", nForeignPointsOfAll
endif

! determine the sum to fond the total number of forein points,
! so that we can allocate room for all the foreign points,
! and it's matching number neighbors

nAllForeignNodes = sum(nForeignPointsOfAll)

if (thisRank.eq.3) then
!    write(*,*) "rank",getrank(),"has all foreign points: ", nAllForeignNodes
endif

allocate(foreignPointsOfAll(nAllForeignNodes))

! Now we can gather all those foreign point IDs from all the processors
! and combonulate them into ForeignPointsOfAll

! First, we need the offset for the points to be added
! so that each partition adds values to the proper place
allocate(offsets(getsize()))
offsets(1) = 0
do iThread=2,getsize()
     offsets(iThread) = offsets(iThread-1) + nForeignPointsOfAll(iThread-1)
enddo

if (thisRank.eq.3) then
!    write(*,*) "rank",getrank(),"has offsets: ", offsets
endif

! gather the foreign point IDs
do iThread=1,getsize()
    call mpi_gatherv(foreignPoints, &      ! sending array
                     nForeignPoints,&	   ! number of elements sent
                     MPI_int,  &           ! type sent
                     foreignPointsOfAll,&  ! destination array
                     nForeignPointsOfAll,& ! number of elements received, from each proc
                     offsets,&             ! offset from start of receiving array
                     mpi_int,&             ! type sent
                     iThread-1, &          ! recipient processor
                     MPI_COMM_WORLD, &     ! communicator
                     error)                ! error

    ! barrier to prevent wayward traffic from becoming entangled
    call mpi_barrier(MPI_COMM_WORLD, error)
enddo

if (thisRank.eq.3) then
!    write(*,*) "rank",getrank(),"has all foreigners: ", foreignPointsOfAll
endif

! Now out threads can go to each of the processors, 
! and add the total number of neighbors for each point

allocate(nNeighbors(nAllForeignNodes))
nNeighbors = 0

do iPoint = 1,nAllForeignNodes
    thisPoint = foreignPointsOfAll(iPoint)

    if ((thisPoint .le. highestLocInd) .and. &
        (thisPoint .ge. smallestLocInd)) then
        ! the points belongs to this partition.
        ! find the number of neighbors

!        write(*,*) "rank",getrank(),"loop ID", iPoint, "pointID = ", thisPoint
#ifdef SPARSE
!        write(*,*) "has verticesmask: ", verticesmask(thisPoint)
        nNeighbors(iPoint) = meshdatactx%itot(verticesmask(thisPoint))
#endif
    endif
enddo

!write(*,*) "rank",getrank(),"has nNeighbors", nNeighbors

! Now add all the nNeighbors vectors:
! The easiest method here is using mpi_reduce, summing it.
! This also sums all the zeroes, though, which will be 
! overly expensive when on run on many partitions, when there are
! many zeroes.

! TODO: optimize this. May not make much difference...

allocate(nNeighborsFull(nAllForeignNodes))

do iThread=1,getsize()
    call mpi_reduce(nNeighbors, &       ! send
                    nNeighborsFull, &   ! receive
                    nAllForeignNodes,&  ! number of variables in array
                    MPI_int, &          ! type
                    MPI_sum, &          ! operator
                    iThread-1, &        ! thread of the receiving array
                    MPI_COMM_WORLD, &
                    error)

    if (error.ne.0) then
        stop "mpi_reduce in meshdatamodule produced an error"
    endif

enddo

!if (thisRank.eq.3) then
!write(*,*) "rank",getrank(),"has offsets", offsets
!write(*,*) "rank",getrank(),"has nNeighborsFull", nNeighborsFull
!endif

! from this nNeighborsFull array, get the number of neighbors for
! the foreign nodes of our own partition only.
allocate(nNeighborsMine(nForeignPoints))

! We can use the offsets used earlier to retrieve the neighbor count:
do iPoint = 1,nForeignPoints
!    write(*,*) "rank", getrank(), "grabs position", iPoint, offsets(getrank())+iPoint
    nNeighborsMine(iPoint) = nNeighborsFull(offsets(thisRank)+iPoint)
enddo

if (thisRank.eq.3) then
!    write(*,*) "rank",getrank(),"has nNeighborcount", nNeighborsMine
endif

! Combine the number of neighbors of all points in elements in our partition
! his includes points in our own partition, and the list of neighbors
! from other partitions, which we created above.

allocate(nAllNeighbors(nPointsInElements))

! sadly, we cannot be sure that all the vertices in our own partition 
! and in elements in our own partition are in the same order as they are
! in ien. The only completely reliable array is pointsInElements.

nOwnPointsInOwnElems = nPointsInElements - nForeignPoints
!write(*,*) "rank",getrank(),"has nOwnPointsInOwnElems", nOwnPointsInOwnElems
!write(*,*) "rank",getrank(),"has locInd", meshdatactx%locind(1:nOwnPointsInOwnElems)
!write(*,*) "rank",getrank(),"has points", pointsInElements(1:nOwnPointsInOwnElems)


! Neighbor count itot refers to points in LocInd, in the same sequence.

! points foreignPoints have their neigbbors stored in nNeighborsMine, and they are 
! also in the same sequence

! add neighborCount of own elements
do iPoint = 1,nOwnPointsInOwnElems
    nAllNeighbors(iPoint) = meshdatactx%itot(iPoint)
    ! the number of neighboring vertices is the same as the number of
    ! adjacent elements... ... EXCEPT for those on the edge of the domain.
    ! There is however no clear way to determine points that are on the edge
    ! of the domain, and because the number of boundary nodes will be
    ! relatively small, especially for large meshes, we will take them along.
enddo
! and that of out foreign points
do iPoint = 1,nForeignPoints
    nAllNeighbors(iPoint + nOwnPointsInOwnElems) = nNeighborsMine(iPoint)
enddo


! Now we have the number of neighbors of all our nodes,
! whether they are our own (in itot) or the foreign nodes
! (in nNeighborsFull)

! We now count actual neighbors we have, and see if we have points
! with neighbors missing. Those neighbors will have to be in other
! partitions.

! The points are stored above in: 


! We browse through IEN and add the number of neighbors to an array
allocate(nNeighborsOfPointsInMyElems(nPointsInElements))

nNeighborsOfPointsInMyElems = 0


if (ndof.eq.2) then
    nodesPerElement = 3
else
    nodesPerElement = 4
endif

do iElem = 1, meshdatactx%nelocal
    do iCorner = 1, nodesPerElement

#ifdef SPARSE
        thisPoint = vertices(meshdatactx%ien(iCorner, iElem)) ! everything here in local numbering
#else
        thisPoint = meshdatactx%ien(iCorner, iElem) ! everything here in local numbering
#endif
!        if (thisRank.eq.3) then
!            write(*,*) "rank",getrank(),"checks point", iCorner, iElem, thisPoint 
!        endif

        ! check if this points is one of ours
        do iPoint = 1, nPointsInElements

!            if (thisRank.eq.3) then
!                write(*,*) "rank",getrank(),"compairing point", iPoint, pointsInElements(iPoint)
!            endif

            if (thisPoint .eq. pointsInElements(iPoint)) then
!                if (thisRank.eq.3) then
!                    write(*,*) "rank",getrank(),"GOT ONE!"
!                endif

                nNeighborsOfPointsInMyElems(iPoint) = nNeighborsOfPointsInMyElems(iPoint) + 1
            endif
        enddo

    enddo
enddo

!if (thisRank.eq.3) then
!    write(*,*) "rank",getrank()," has nNeighborsOfPointsInMyElems", nNeighborsOfPointsInMyElems ! OK
!    write(*,*) "rank",getrank(),"all own Neighborcount", nAllNeighbors ! OK
!    write(*,*) "rank",getrank(),"has points in elements", pointsInElements ! tested and OK
!endif

! Count the number of incomplete points
! (most of the points will be complete because they are inside the partition)

nIncomplete = 0

do iPoint = 1,nPointsInElements
    if (nNeighborsOfPointsInMyElems(iPoint) .ne. nAllNeighbors(iPoint)) then
        nIncomplete = nIncomplete + 1
    endif
enddo

allocate(incompleteNodes(nIncomplete))

nIncomplete = 0

do iPoint = 1,nPointsInElements
    if (nNeighborsOfPointsInMyElems(iPoint) .ne. nAllNeighbors(iPoint)) then
        nIncomplete = nIncomplete + 1
#ifdef SPARSE
	incompleteNodes(nIncomplete) = vertices(iPoint)
        incompleteNodes(nIncomplete) = pointsInElements(iPoint)
#else
        incompleteNodes(nIncomplete) = iPoint
#endif
    endif
enddo


!write(*,*) "rank",getrank(),"has incomplete points", incompleteNodes ! tested and OK



! Now we can exchange between the processors how many
! incomplete nodes there are globally.

allocate(nIncompletePerProc(getsize()))

call mpi_allgather(nIncomplete, 1, MPI_int, nIncompletePerProc, 1, mpi_int, MPI_COMM_WORLD, error)
if (error.ne.0) then
    stop "mpi_allgather failed with error code"
endif

nAllIncomplete = sum(nIncompletePerProc)

! Now we can allocate room for all the incomplete nodes.

allocate(allIncompletePoints(nAllIncomplete))

! And we can share the incomplete node IDs

!allocate(offsets(getsize()))
offsets(1) = 0
do iThread=2,getsize()
     offsets(iThread) = offsets(iThread-1) + nIncompletePerProc(iThread-1)
enddo

!if (thisRank.eq.3) then
!    write(*,*) "rank",getrank(),"has incompleteNodes", incompleteNodes
!endif

! gather the foreign incomplete point IDs
do iThread=1,getsize()
    call mpi_gatherv(incompleteNodes, &    ! sending array
                     nIncomplete,&         ! number of elements sent
                     MPI_int,  &           ! type sent
                     allIncompletePoints,& ! destination array
                     nIncompletePerProc,&  ! number of elements received, from each proc
                     offsets,&             ! offset from start of receiving array
                     mpi_int,&             ! type sent
                     iThread-1, &          ! recipient processor
                     MPI_COMM_WORLD, &     ! communicator
                     error)                ! error

    ! barrier to prevent wayward traffic from becoming entangled
    call mpi_barrier(MPI_COMM_WORLD, error)
enddo

!write


!if (thisRank.eq.3) then
!    write(*,*) "rank",getrank(),"has incompleteNodes of each", allIncompletePoints ! tested and OK
!endif


! TODO sort and uniquefy the array. That will half the transport of data.
! But will need a way to look up as well. I wonder which is fastest...

! Now that all processors know the points with missing neightbors all over the domain),
! each of the processors can walk through this list,
! and see if it has elements that border on his point.

! determine highest and smallest points.
! This will provide a quick check on whether 
! this partition has an element boundaring on it

! we could simply check verticesmask, but that is a global array 
! that we eventually want to get rid of, to improve scaling.

smallestPointIDInMyElems = meshdatactx%nvglobal + 1
highestPointIDInMyElems = 0

if (NDOF.eq.2) then
    nPointsPerElement = 3
else
    nPointsPerElement = 4
endif


do iElem = 1,meshdatactx%nelocal
    do iPoint = 1,nPointsPerElement
#ifdef SPARSE
        thisPoint = vertices(meshdatactx%ien(iPoint, iElem))
#else
        thisPoint = meshdatactx%ien(iPoint, iElem)
#endif
        if (thisPoint .gt. highestPointIDInMyElems) then
            highestPointIDInMyElems = thisPoint
        endif
        if (thisPoint .lt. smallestPointIDInMyElems) then
            smallestPointIDInMyElems = thisPoint
        endif
    enddo
enddo


!write(*,*) "rank",getrank(),"has minmax points ", smallestPointIDInMyElems, highestPointIDInMyElems


! no we can do a precheck of the points that could have elements in this partition:

#ifdef SPARSE
!write(*,*) "rank",getrank(),"has vertices: ", vertices

!write(*,*) "rank",getrank(),"has verticesmask: ", verticesmask
#endif

! allocate room to add element numbers.
allocate(nElemsPerPoint(nAllIncomplete)) ! for every ghost point, stores the number of 
                                          ! elements that this partition has, bordering on it.
nElemsPerPoint = 0


if (ndof.eq.2) then
    nMaxElemsPerPoint = 10
else if (ndof.eq.3) then
    nMaxElemsPerPoint = 25
else
    stop "ndof should be 2 or 3"
endif

allocate(elemIDsToThisPoint(nAllIncomplete, nMaxElemsPerPoint))

elemIDsToThisPoint = 0


! find the elements this partition has on the incomplete Points

do iPoint = 1,nAllIncomplete
    thisPoint = allIncompletePoints(iPoint)

    ! quick check, toprevent unnecessary looping withing the if statement
    if ((thisPoint .le. highestPointIDInMyElems) .and. &
        (thisPoint .ge. smallestPointIDInMyElems)) then

        ! Now we know that this points might belong to an element here.
        ! Know that not all points between smallestPointIDInMyElems
        ! and highestPointIDInMyElems have to belong to elements.
        ! There can be gaps.

        ! We can step through IEN and see if we can find a matching element
        do iElem = 1, meshdatactx%nelocal
            do jPoint = 1, nPointsPerElement
#ifdef SPARSE
                if (vertices(meshdatactx%ien(jPoint, iElem)) .eq. thisPoint) then
#else
                if (meshdatactx%ien(jPoint, iElem) .eq. thisPoint) then
#endif
                    ! YES, WE HAVE FOUND A MISSING GHOST ELEMENT
                    ! well, not necessarily. It can be so that 
                    ! this element belongs to an incomplete point
                    ! of itself, and as such was already known.
                    ! A filtering mechanism might be more
                    ! computationally expensive then just taking
                    ! them along. 
                    ! TODO: see if a filter mechanism is worth it.

                    ! And of those elements, we need the points as well, 
                    ! which we can obtain from the vertices array. 

                    ! Store these elements locally on this processor.
                    ! Next we can count them, and after that, we can allocate
                    ! a global array in which to put them.
                    nElemsPerPoint(iPoint) = nElemsPerPoint(iPoint) + 1

                    if (nElemsPerPoint(iPoint) .gt. nMaxElemsPerPoint) then
                        stop "Whaaaa, found a star in the mesh. Not ehough room allocated"
                    endif

                    elemIDsToThisPoint(iPoint, nElemsPerPoint(iPoint)) = meshdatactx%locElt(iElem)


                endif
            enddo
        enddo

    endif
enddo

!if (thisRank.eq.3) then
!    write(*,*) "rank",getrank(),"has n own elems to all points", nElemsPerPoint 
!endif


!call delay()
do iThread = 1,6
    if (iThread-1 .eq. getrank()) then
        do iPoint = 1, nAllIncomplete
!            write(*,*) "rank", getrank(), &
!                        "elements to node", allIncompletePoints(iPoint), &
!                        "has elems", elemIDsToThisPoint(iPoint, 1:nElemsPerPoint(iPoint))  ! tested and OK
        enddo
!        write(*,*) "rank",getrank(),"has n own elems to all points", nElemsPerPoint
    endif
!    call delay()
enddo
!call delay()





! sum the count of all those points, by adding them over the processors

allocate(nAllElemsPerPoint(nAllIncomplete))

do iThread=1,getsize()
    call mpi_reduce(nElemsPerPoint, &    ! send
                    nAllElemsPerPoint, & ! receive
                    nAllIncomplete,&     ! number of variables in array
                    MPI_int, &           ! type
                    MPI_sum, &           ! operator
                    iThread-1, &         ! thread of the receiving array
                    MPI_COMM_WORLD, &
                    error)

    if (error.ne.0) then
        stop "mpi_reduce in meshdatamodule produced an error"
    endif

enddo

! Now we know the number of neighboring elements of every incomplete point

!if (thisRank.eq.3) then
!    write(*,*) "rank",getrank(),"has n all elems to all points", nAllElemsPerPoint ! WOHOO, it works! Teste, 
! write(*,*) "n points", size(nAllElemsPerPoint,1)
!endif

! If we put all the elements in a data structure, we get a large amount of data
! That does not scale well, and will limit the resolution.
! for example, a rectangular 2D domain divided into 4.655.387 elems
! and split into 16 partitions, has 24408 nAllIncomplete
! allocation of nAllIncomplete * nPartitions * maxNbrs * (3 or 4 points per elem) * 4 bytes per int

! 48 MB. At higher resolution and in 3D, this will become a limiting factor.

! So, we will first let every partition know how many elements it has, 
! contributing to every incomplete point. based on that, we will set up
! the data structure to send/receive only the required elements to each partition.
! This will save memory and transmission time, at the price of computation time.
! There is no such thing as free beer.

allocate(partitionContributions(nAllIncomplete, getsize()))

do iPoint = 1,nAllIncomplete
    ! +1 because rank counts zero based and array dimensions 1-based
    partitionContributions(iPoint,thisrank+1) = nElemsPerPoint(iPoint)
enddo

if (thisRank.eq.3) then
!    write(*,*) "partitionContributions: ", partitionContributions ! checked and OK
endif

! now we can exchange the data, by adding all these arrays

allocate(otherPartitionsHaveCount(nAllIncomplete, getsize()))

do iThread=1,getsize()
    call mpi_reduce(partitionContributions,   & ! send
                    otherPartitionsHaveCount, & ! receive
                    nAllIncomplete*getsize() ,& ! number of variables in array
                    MPI_int, &           ! type
                    MPI_sum, &           ! operator
                    iThread-1, &         ! thread of the receiving array
                    MPI_COMM_WORLD, &
                    error)

    if (error.ne.0) then
        stop "mpi_reduce in meshdatamodule produced an error"
    endif

enddo

if (thisRank.eq.3) then
!    write(*,*) "otherPartitionsHaveCount ", otherPartitionsHaveCount ! checked and OK
endif



deallocate(partitionContributions)


!write(*,*) "rank",getrank(),"has n from part", otherPartitionsHaveCount

!write(*,*) "rank",getrank(),"has offsets: ", offsets

! make a list of incomplete nodes and the partitions that
! contain the missing elements, and the number of elements 
! for the missing partition.

allocate(MyNeighborPartitionElemCount(nIncomplete, getsize()))

do iPoint = 1,nIncomplete
    do iThread = 1,getsize()
        MyNeighborPartitionElemCount(iPoint, iThread) = &
            otherPartitionsHaveCount(offsets(getrank()+1) + iPoint, iThread)
    enddo
enddo

! print list 
!do iThread = 1,getsize()
!    write(*,*) "rank",getrank(),"has to rank",iThread-1,"search", MyNeighborPartitionElemCount(:,iThread)
    ! checked and OK
!enddo

! Now every partition knows how much has to retrieve from every partition.

! allocate the retrieval array
nReceiveElems = 0

do iThread = 1,getsize()
    if (iThread-1 .ne. getrank()) then
        do iPoint = 1, nIncomplete
            nReceiveElems = nReceiveElems + 1
        enddo
    endif
enddo


!allocate(receiveElements(nReceiveElems))
!receiveElements = 0
allocate(receiveElemVertices(nReceiveElems, nPointsPerElement))
receiveElemVertices = 0


! Now we know which partition has how many elements.
! elemIDsToThisPoint contains the element nubers to incomplete nodes.

! The incomplete nodes of this partition is either on the domain edge,
! or it is on a partition boundary and is also an incomplete node of
! another partition.



nSendNodes=0
nSendElements=0
nSendElementsInThisProc = 0

allocate(sendNodes(2*nIncomplete))
allocate(sendElements(10*nIncomplete))
allocate(sendElementsOffsets(getsize()))
allocate(nSendElementsToEachThread(getsize()))
allocate(nReceiveElementsFromEachThread(getsize()))

sendNodes = 0
sendElements = 0
sendElementsOffsets = 0
nSendElementsToEachThread = 0
nReceiveElementsFromEachThread = 0

do iThread = 1,getsize()


    nSendElementsInThisProc = 0

    if ((iThread-1).ne.getrank()) then

    ! We are preparing the ghost elements for thread iThread
    ! Every partition knows what element IDs he shares with any points.
    ! Because the IENs are known here, we can also assemble the points needed here.

    ! for every partition, determine whether it has incomplete points 
    ! in common. Assemble them

!    allocate(nIncomplete)
!    write(*,*) "rank",getrank(),"has incomplete", incompleteNodes
!    write(*,*) "rank",getrank(),"has allIncomplete", allIncompletePoints
!    write(*,*) "rank",getrank(),"has offsets", offsets
!    write(*,*) "rank",getrank(),"has nIncompletePerProc", nIncompletePerProc
    ! this thread will determine which incomplete points it has in common with 
    ! partition iThread

!        sendElementsOffsets(iThread) = nSendElements

        nSharedIncomplete = 0
!    allocate(SharedIncompletePointWithPartIthread(nIncomplete))    

        do iPoint = 1, nIncomplete
            ! for every incomplete point in this partition
            do jPoint = 1, nIncompletePerProc(iThread)

!                write(*,*) "i, j: ", iPoint, jPoint

                if (incompleteNodes(iPoint) .eq. allIncompletePoints(offsets(iThread)+jPoint)) then
                    ! We have a point of this partition whose elements have to be sent to
                    ! partition iThread
                    thisPoint = incompleteNodes(iPoint)

!                    write(*,*) "searching elements that belong to point", thisPoint

                    ! there is room for optimization here...
                    searchID = 0
                    pointFound = .false.
                    do while (.not. pointFound)
                        searchID = searchID + 1
                        checkPoint = allIncompletePoints(searchID)

!                        write(*,*) "comparing with point", checkPoint

                        if (checkPoint.eq.thisPoint) then
!                            write(*,*) "hebbes!"
                            do iElem = 1,nElemsPerPoint(searchID)
                                pointFound = .true.
                                nSendElements = nSendElements + 1
                                nSendElementsInThisProc = nSendElementsInThisProc + 1
                                if (nSendElements .le. size(sendElements,1)) then
                                    sendElements(nSendElements) = elemIDsToThisPoint(searchID,iElem)
                                else
                                    write(*,*) "sendElements array too small", nSendElements
                                    stop
                                endif
                            enddo
                        endif
                    enddo

                    nSendNodes = nSendNodes + 1
                    sendNodes(nSendNodes) = incompleteNodes(iPoint)
                endif
            enddo
        enddo

!    write(*,*) "rank",getrank(), &
!               "shares with partition", iThread, &
!               "points", SharedIncompletePointWithPartIthread(1:nSharedIncomplete)
! checked and OK!

    endif

    nSendElementsToEachThread(iThread) = nSendElementsInThisProc

enddo

do iThread=2,getsize()
        sendElementsOffsets(iThread) = sendElementsOffsets(iThread-1) + \
                                       nSendElementsToEachThread(iThread-1)
enddo


do iThread=1,getsize()
    if (getrank().eq.iThread-1) then
        write(*,*) "rank",getrank()," local nelems", meshdatactx%NelocalOfAllPartitions
    endif
!    call delay()
enddo


! Determine offset of the element numbering, from the 
! local element numbers that has been read in tecin.dat.
allocate(elementNumberingOffset(getsize()))
elementNumberingOffset = 0

do iThread = 2,getsize()
    elementNumberingOffset(iThread) = \
        meshdatactx%NelocalOfAllPartitions(iThread) + \
        elementNumberingOffset(iThread-1)
enddo

do iThread=1,getsize()
    if (getrank().eq.iThread-1) then
        write(*,*) "rank",getrank()," local elems offsets", elementNumberingOffset
    endif
!    call delay()
enddo

! Using that, and array sendElements, we can find the vertices that are to be exchanged.
do iElem =1,nSendElements
        write(*,*) "rank", getrank(), "checks"

    vertex1Send(iElem) = meshdatactx%ien(1,\
        sendElements(iElem) - elementNumberingOffset(getrank()+1))
    vertex2Send(iElem) = meshdatactx%ien(2,\
        sendElements(iElem) - elementNumberingOffset(getrank()+1))
    vertex3Send(iElem) = meshdatactx%ien(3,\
        sendElements(iElem) - elementNumberingOffset(getrank()+1))
    vertex4Send(iElem) = meshdatactx%ien(4,\
        sendElements(iElem) - elementNumberingOffset(getrank()+1))
enddo

do iThread=1,getsize()
    if (getrank().eq.iThread-1) then
        write(*,*) "rank",getrank()," sends vertices1", vertex1Send
    endif
!    call delay()
enddo

do iThread=1,getsize()
    if (getrank().eq.iThread-1) then
        write(*,*) "rank",getrank()," sends vertices2", vertex2Send
    endif
!    call delay()
enddo

do iThread=1,getsize()
    if (getrank().eq.iThread-1) then
        write(*,*) "rank",getrank()," sends vertices3", vertex3Send
    endif
!    call delay()
enddo

do iThread=1,getsize()
    if (getrank().eq.iThread-1) then
        write(*,*) "rank",getrank()," sends vertices4", vertex4Send
    endif
!    call delay()
enddo





do iThread=1,getsize()
    if (getrank().eq.iThread-1) then
        write(*,*) "rank",getrank(),"nSendElementsToEachThread", nSendElementsToEachThread
    endif
!    call delay()
enddo

!write(*,*) "rank",getrank(),"has n own elems to all points", nElemsPerPoint
!write(*,*) "rank",getrank(),"send nSendNodes", nSendNodes
!write(*,*) "rank",getrank(),"sendNodes", sendNodes(1:nSendNodes)

!write(*,*) "rank",getrank(),"send nSendElements", nSendElements


do iThread=1,getsize()
    if (getrank().eq.iThread-1) then
        write(*,*) "rank",getrank(),"sendElementsOffsets", sendElementsOffsets
    endif
!    call delay()
enddo


!write(*,*) "rank",getrank(),"sendElements", sendElements(1:nSendElements)


!   scatter the elements to their destination processors
! First scatter the number of elements going to each processor
! So that it can allocate room for the elements.




do iThread=1,getsize()
    ! Send parts of this partition to receive elems in all other partitions
    call mpi_scatter(nSendElementsToEachThread, &
                     1, &
                     MPI_int, &
                     nReceiveElementsFromEachThread(iThread), &
                     1, &
                     MPI_int, &
                     iThread-1, &
                     MPI_COMM_WORLD, &
                     error)
enddo


do iThread=1,getsize()
    if (getrank().eq.iThread-1) then
        ! checked and OK
        write(*,*) "rank",getrank()," nReceiveElements", nReceiveElementsFromEachThread
    endif
!    call delay()
enddo

nReceiveElems = sum(nReceiveElementsFromEachThread)

allocate(receiveOffsets(getsize()))
receiveOffsets = 0

do iThread=2,getsize()
    receiveOffsets(iThread) = nReceiveElementsFromEachThread(iThread-1) + receiveOffsets(iThread-1)
enddo

!do iThread=1,getsize()
!    if (getrank().eq.iThread-1) then
!        ! checked and OK
!        write(*,*) "rank",getrank(),"receive Offsets", receiveOffsets
!    endif
!    call delay()
!enddo

!do iThread=1,getsize()
!    if (getrank().eq.iThread-1) then
!        ! checked and OK
!        write(*,*) "rank",getrank(),"nReceive", nReceiveElems
!    endif
!    call delay()
!enddo





! Now we all know the offsets, and we all know what to send to what.
! each processor can gather the elements numbers of its neighboring partitions.

! The whole advantage is that we no longer need to share ALL border
! elemets, but ONLY those of the neighboring partitions.


!allocate(receiveElements(nReceiveElems))

!from openMPI documentation:

!MPI_IGATHERV(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNTS,
!        DISPLS, RECVTYPE, ROOT, COMM, REQUEST, IERROR)

!    <type>    SENDBUF(*), RECVBUF(*)
!    INTEGER    SENDCOUNT, SENDTYPE, RECVCOUNTS(*), DISPLS(*)
!    INTEGER    RECVTYPE, ROOT, COMM, REQUEST, IERROR


!    call mpi_gatherv(incompleteNodes, &    ! sending array
!                     nIncomplete,&         ! number of elements sent
!                     MPI_int,  &           ! type sent
!                     allIncompletePoints,& ! destination array
!                     nIncompletePerProc,&  ! number of elements received, from each proc
!                     offsets,&             ! offset from start of receiving array
!                     mpi_int,&             ! type sent
!                     iThread-1, &          ! recipient processor
!                     MPI_COMM_WORLD, &     ! communicator
!                     error)                ! error




!### to validate the element ID communications. Checked and OK :-)
!### Takes nproc^2 seconds to check, because of delay statements.
!### Those statement keep the order of data in standard out readable.
!write(*,*) "What all processors are going to send to processor 0"
!do iThread=1,getsize()
!    if (getrank().eq.iThread-1) then
!        do iDestination=1,getsize()
!        ! checked and OK
!        startIndex = sendElementsOffsets(iDestination) + 1
!        endIndex   = sendElementsOffsets(iDestination) + nSendElementsToEachThread(iDestination)
!        if (startIndex.gt.0 .and. endIndex.gt.0 .and. startIndex.lt.endIndex) then
!            write(*,*) "rank",getrank(), &
!                       "sending to",iDestination-1, &
!                       "from: ", startIndex, &
!                       "to",endIndex, &
!                       "integers: ", sendElements(startIndex:endIndex)
!        else
!            write(*,*) "rank",getrank(),"sends NOTHING to",iDestination-1,"! MUHAHAHAHAHA"
!        endif
!        call delay()
!        enddo
!    endif
!    call delay()
!enddo






! each partition must know have much every other partition is going to receive
allocate(eachRankReceivesNElems(getsize()))

call mpi_allgather(nReceiveElems, 1, MPI_int, eachRankReceivesNElems, 1, mpi_int, MPI_COMM_WORLD, error)

!write(*,*) "rank", getrank(), "has all receives", eachRankReceivesNElems

call mpi_barrier(0,error)

! each partition must know the receive offsets of every partition.

allocate(receiveOffsetOfEachThread(getsize()**2))

call mpi_allgather(receiveOffsets, &
                   getsize(), &
                   MPI_int, &
                   receiveOffsetOfEachThread, &
                   getsize(), &
                   mpi_int, &
                   MPI_COMM_WORLD, &
                   error)

!write(*,*) "rank", getrank(), "has all offsets", receiveOffsetOfEachThread

! each partition already knows it sendsize and send vector.




call mpi_barrier(0,error)

!call delay()


!call delay()
!stop "wohoo"

!receiveOffsets = receiveOffsets + 1 ! offset counts as 0 or 1 based?... blijkbaar niet

do iThread=1,getsize()


    ! every thread does its own gathering iteration.

    ! Determine the offsets    


    if (getrank().eq.iThread-1) then
        ! allocate and initialize receiving array only when 
        ! we are in the loop iteration of this partition
        allocate(receiveElements(eachRankReceivesNElems(iThread)))
        receiveElements = 0
    endif

    ! receive quantities are dependent

    if (getrank().eq.iThread-1) then
!        write(*,*) "### rank",getrank(), "will receive elems into ", receiveElements
    endif


    do i=1,getsize()
        iDestination = i-1
        if (getrank().eq.iDestination) then
!            write(*,*) "***************************************************************************"
!            write(*,*) "*** iteration for gather by ", iThread -1, "and sender",iDestination, "*****"
!            write(*,*) "***************************************************************************"
!            write(*,*) "select nElems ", iThread, "from ", nSendElementsToEachThread
!            write(*,*) "sending ", nSendElementsToEachThread(iThread), "elements"

            startIndex = sendElementsOffsets(iThread) + 1
            endIndex   = sendElementsOffsets(iThread) + nSendElementsToEachThread(iThread)

!            write(*,*) "take them from sendElements from",startIndex,"to",endIndex

!            write(*,*) "sending ",     sendElements(startIndex:endIndex)
 !           write(*,*) "there are ", endIndex - startIndex + 1, "of those"

            startOffsetIndex = (iThread-1)*getsize()+1
            endOffsetIndex   = iThread*getsize()

!            write(*,*) "Receiving ", nReceiveElementsFromEachThread, "elements"
!            write(*,*) "Offset index: ", startOffsetIndex, "to", endOffsetIndex
!            write(*,*) "to offsets", receiveOffsetOfEachThread(startOffsetIndex:endOffsetIndex)

!            write(*,*) "length of receiveElems: ", size(receiveElements,1)

        endif
!        call mpi_barrier(0,error)

!        call delay()


    enddo


!    write(*,*) "rank",getrank(),"gathering"
!    call mpi_barrier(0,error)


!    call mpi_gatherv(henkie, &             ! data that is transmitted
!                     nSendNumbers, &       ! n data that each processors sends
!                     MPI_int, &            ! send type
!                     globalhenk, &         ! receiving array
!                     nReceiveNumbers, &    ! n data received from each proc
!                     receiveOffsets, &     ! offset in final vector
!                     MPI_int, &            ! receive type
!                     i-1, &                  ! the receiving thread
!                     MPI_COMM_WORLD, & 
!                     error)


    call mpi_gatherv(sendElements(startIndex:endIndex), &
                        nSendElementsToEachThread(iThread), &
                     MPI_int, &
                        receiveElements, &
                        nReceiveElementsFromEachThread, &
                        receiveOffsetOfEachThread(startOffsetIndex:endOffsetIndex), &
                        MPI_int, &
                     iThread-1, &
                        MPI_COMM_WORLD, &
                     error)

    if (getrank().eq.iThread-1) then
        write(*,*) "### rank",getrank(), "Received elems", receiveElements
    endif

!    write(*,*) "rank",getrank(),"gathering complete"


    call mpi_barrier(0,error)

!    call delay()


enddo




call mpi_barrier(0,error)

!do iThread=1,getsize()
!    if (getrank().eq.iThread-1) then
        ! checked and OK
!        write(*,*) "### rank",getrank(), "Received elems", receiveElements
!    endif
!    call delay()
!enddo



!call delay() ! so that every output can be printed



!deallocate(pointsOfAllBorderElements)
!deallocate(elemIDdifferentNeighbors)

!deallocate(nElemsPerVertex)
!deallocate(ElemIDPervertex)
!deallocate(offsets)

stop "This is the eeeheheeend"


end subroutine

