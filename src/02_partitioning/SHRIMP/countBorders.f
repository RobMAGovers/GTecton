program countBorders

! This reads a mesh and counts the number of triangles that separate the partitions.
! This is a norm for the quality of the partitioning.

! Assumes mesh to be 3D.

implicit none


integer :: nPoints, nElements, nPartitions
integer :: iPoint,  iElement,  iPartition
integer :: jElement
integer :: iElem, jElem

double precision, allocatable :: coords(:,:)
integer,          allocatable :: connectivity(:,:)
integer,          allocatable :: elemPartition(:)
integer,          allocatable :: nNeighborsPerPoint(:)
integer,          allocatable :: nElemsPerPoint(:)
integer,          allocatable :: elemNeighbors(:,:)
integer,          allocatable :: nNeighborsFound(:)

integer :: dummy, np, ne

integer :: nBorder, neighborsFound
integer :: thisPoint

character(len=1024) :: record

integer, allocatable :: tempIntArray(:)

type listOfLists
    integer, allocatable :: list(:)
end type

type(listOfLists), allocatable :: neighborLists(:)
integer :: nNeighbors

type(listOfLists), allocatable :: elemsPerPointLookup(:)

logical, external :: twoElemsAreNeighbors

nPoints = 0
nElements = 0

!------------------------------------------------------------
! Read input data; partitioned mesh, GTECTON style
! needs: partition.info
!        tecin.dat.partf.nps
!        tecin.dat.partf.elm
!------------------------------------------------------------


open(unit=42, file="partition.info")
read(42,*) nPartitions
do iPartition = 1, nPartitions
    read(42,*) dummy, np, ne
    nPoints = nPoints + np
    nElements = nElements + ne

enddo
close(42)

allocate(coords(3, nPoints))
allocate(connectivity(4,nElements))
allocate(elemPartition(nElements))

allocate(neighborLists(nPoints))
allocate(nNeighborsPerPoint(nPoints))

!open(unit=43, file="tecin.dat.partf.nps")
!do iPoint = 1, nPoints
!    read(43,"(a)") record
!    read(record,*) dummy, dummy,dummy, coords(:, iPoint), nNeighborsPerPoint(iPoint)
!    allocate(neighborLists(iPoint)%list(nNeighbors))
!    read(record,*) dummy,dummy,dummy, coords(:, iPoint), nNeighborsPerPoint(iPoint), neighborLists(iPoint)%list
!enddo
!close(43)


open(unit=43, file="tecin.dat.partf.elm")
do iElement = 1, nElements
    read(43,*) elemPartition(iElement),dummy,dummy, connectivity(:, iElement)
enddo
close(43)

!---------------------------------------------------------------------------
! we need to build a list of the neighboring elements
! for every elements. 
! First a lookup table is built from the neighboring point IDs, which is known.

! This lookup table gains an entry for every element that has this node
!---------------------------------------------------------------------------

allocate(elemsPerPointLookup(nPoints))
allocate(nElemsPerPoint(nPoints))

do iPoint = 1, nPoints
    allocate(elemsPerPointLookup(iPoint)%list(80))
    elemsPerPointLookup(iPoint)%list = 0
    nElemsPerPoint(iPoint) = 0
enddo


do iElement = 1, nElements
    do iPoint = 1, 4
        thisPoint = connectivity(iPoint, iElement)

        nElemsPerPoint(thisPoint) = nElemsPerPoint(thisPoint) + 1
        if (nElemsPerPoint(thisPoint) .gt. 80) then
            stop "found know of more than 80 elements to 1 points. WHAAAA"
        endif
        elemsPerPointLookup(thisPoint)%list(nElemsPerPoint(thisPoint)) = iElement
    enddo
enddo

! Now we can step through the lookup table.
! For every point, we loop through the elements
! that are adjacent to this point

allocate(elemNeighbors(4, nElements))
elemNeighbors = 0

allocate(nNeighborsFound(nElements))
nNeighborsFound = 0

do iPoint = 1, nPoints
    do iElement = 1, nElemsPerPoint(iPoint) - 1
        do jElement    = iElement, nElemsPerPoint(iPoint)
            iElem = elemsPerPointLookup(iPoint)%list(iElement)
            jElem = elemsPerPointLookup(iPoint)%list(jElement)

            if (twoElemsAreNeighbors(connectivity(:, iElem), connectivity(:, jElem) )) then
                ! yay, is neighbor!
                ! check whether it is not already known as a neighbor
                if (iElem .ne. elemNeighbors(1, jElem) .and. &
                    iElem .ne. elemNeighbors(2, jElem) .and. &
                    iElem .ne. elemNeighbors(3, jElem) .and. &
                    iElem .ne. elemNeighbors(4, jElem)) then
                    ! true new neighbor!
                    nNeighborsFound(jElem) = &
                    nNeighborsFound(jElem) + 1
                    elemNeighbors(nNeighborsFound(jElem), jElem) = iElem
                endif
                ! and the other way around
                if (jElem .ne. elemNeighbors(1, iElem) .and. &
                    jElem .ne. elemNeighbors(2, iElem) .and. &
                    jElem .ne. elemNeighbors(3, iElem) .and. &
                    jElem .ne. elemNeighbors(4, iElem)) then
                    ! true new neighbor!
                    nNeighborsFound(iElem) = &
                    nNeighborsFound(iElem) + 1
                    elemNeighbors(nNeighborsFound(iElem), iElem) = jElem
                endif
            endif
        enddo
    enddo
enddo

! now all neighbors are known.
! Count the times when they are from a different partition.

nBorder = 0

do iElement = 1, nElements
    do jElement = 1,4
        jElem = elemNeighbors(jElement, iElement)

        ! element can be on domain boundary,
        ! then not all four neighbors are filled,
        ! and we can check for 0 entries.
        if (jElem.ne.0) then

            ! prevent double count
            if (jElem .lt. iElement) then
                if (elemPartition(iElement) .ne. &
                    elemPartition(jElem)) then
                    nBorder = nBorder + 1
                endif
            endif
        endif
    enddo
enddo

write(*,*) "Number of border surfaces = ", nBorder


end program




logical function twoElemsAreNeighbors(conn1, conn2)

implicit none

integer :: conn1(4), conn2(4)

integer :: nMatch, i, j

! If three points in conn1 are also in conn2, there is a match

nMatch = 0

do i = 1,4
    do j = 1,4
        if (conn1(i).eq.conn2(j)) then
            nMatch = nMatch + 1
        endif
    enddo
enddo

if (nMatch .eq. 3) then
    twoElemsAreNeighbors = .true.
else
    twoElemsAreNeighbors = .false.
endif

end function
