program partitionQuality

! this program computes the quality of the partition

! it reads a partitioned gtecton mesh and returnd the amound of triangles that
! have a tetrahedron beloning to a different partition

implicit none

integer nPoints, nElems
integer iPoint,  iElem, iSide


character (len=19), parameter :: tecindatpartfnps = "tecin.dat.partf.nps"
character (len=19), parameter :: tecindatpartfelm = "tecin.dat.partf.elm"




type brokenComb
    integer, allocatable :: list(:)
end type

type(brokenComb), allocatable :: ElementsPerNode(:)
type(brokenComb), allocatable :: neighborNodes(:)

integer, allocatable :: nElementsPerNode(:)
integer, allocatable :: nneighbours(:)
integer, allocatable :: elementPartition(:)
integer, allocatable :: elemNeighborList(:,:)
integer, allocatable :: connectivity(:,:)
integer              :: d1,d2,d3    ! dummies
double precision     :: x,y,z ! coords are not relevant, so not stored.
integer, allocatable :: tempIntBuffer(:)


integer              :: checkElem, checkElemID
integer              :: pos
integer              :: pointsOfTriangle(3)
integer              :: thisPoint, basePoint
integer              :: listLength

logical, external :: numberInList

integer :: boundaryCount



!***** read tecin.dat.part.elm ***
    call countlines(trim(tecindatpartfelm),nelems, 32)
    nelems = nelems - 1 ! remove the 'end' line

    if (nelems.eq.-1) then
        write(*,*) "Partitioned element file empty: ", trim(tecindatpartfelm)
    endif

    allocate(connectivity(4,nelems))
    allocate(elementPartition(nelems))

    open(unit = 233, file=trim(tecindatpartfelm))
    do iElem=1,nElems
        read (233,*) elementPartition(iElem), d2, d3, connectivity(:,iElem)
    enddo
    close(233)

!***** read tecin.dat.part.nps ***

    call countlines(trim(tecindatpartfnps),nPoints, 42)
    nPoints = nPoints - 1 ! -1 to remove the 'end' line

    if (nPoints .eq. -1) then
        write(*,*) "Partitioned nodal point file empty: ", trim(tecindatpartfnps)
    endif

    allocate(tempIntBuffer(200))
    tempIntBuffer = 0
    allocate(neighborNodes(nPoints))
    allocate(nneighbours(nPoints))

    open(unit = 244, file=trim(tecindatpartfnps))
    do iPoint=1,nPoints
        read (244,*) d1, d2, d3, x, y, z,  nneighbours(iPoint), tempIntBuffer(1:nneighbours(iPoint))
        allocate( neighborNodes(iPoint)%list( nneighbours(iPoint) ))
        neighborNodes(iPoint)%list = tempIntBuffer(1:nneighbours(iPoint))
        tempIntBuffer = 0
    enddo
    close(244)

    deallocate(tempIntBuffer)

!***** build a list of elements per node, from the neighbor table

    allocate(nElementsPerNode(nPoints))
    allocate(ElementsPerNode(nPoints))
    do iPoint = 1, nPoints
        nElementsPerNode(iPoint) = 0
        allocate(ElementsPerNode(iPoint)%list(20))
        ElementsPerNode(iPoint)%list = 0
    enddo

    do iElem = 1, nElems
        do iPoint = 1,4
            thisPoint = connectivity(iPoint, iElem)

            if (.not. numberInList(ElementsPerNode(thisPoint)%list, & ! list
                                   nElementsPerNode(thisPoint), &     ! length
                                   iElem) ) then                      ! search for this nr      

                listLength = size(ElementsPerNode(thisPoint)%list,1)

                ! we have a new element for the list
                if ( ElementsPerNode(thisPoint)%list(listLength) .ne. 0) then

                    ! list full, make bigger, temporarily store data in tempIntBuffer
                    allocate(tempIntBuffer(nElementsPerNode(thisPoint)))
                    tempIntBuffer = ElementsPerNode(thisPoint)%list
                    deallocate(ElementsPerNode(thisPoint)%list)
                    allocate(ElementsPerNode(thisPoint)%list(nElementsPerNode(thisPoint)+20))
                    ElementsPerNode(thisPoint)%list(1:nElementsPerNode(thisPoint)) = tempIntBuffer
                    deallocate(tempIntBuffer)

                endif

                nElementsPerNode(thisPoint) = nElementsPerNode(thisPoint) + 1
                ElementsPerNode(thisPoint)%list(nElementsPerNode(thisPoint)) = iElem

            endif
        enddo
    enddo


!***** from this list, build a list of neighboring elements for each element

    allocate(elemNeighborList(4,nElems))
    elemNeighborList = 0
    
    do iElem = 1, nElems
        do iSide = 1,4

            ! select the poins that belong to this side
            ! note that is nog the GTECTON convention
            ! side 1, points 2,3,4
            ! side 2, points 1,3,4
            ! side 3, points 1,2,4
            ! side 4, points 1,2,3
            pos = 1
            do iPoint = 1,4
                if (iPoint .ne. iSide) then
                    pointsOfTriangle(pos) =  connectivity(iPoint,iElem)
                    pos = pos + 1
                endif
            enddo

            ! now find and adjacent element that from the list adjacent to the node
            ! that has all these three nodes.
            ! We grab elements from the list of basePoint.
            if (iSide .eq. 1) then
                basePoint = connectivity(2,iElem)
            else
                basePoint = connectivity(1,iElem)
            endif

            do checkElem = 1, nElementsPerNode(basePoint)
                checkElemID = ElementsPerNode(basePoint)%list(checkElem)
            
                if ( numberInList(connectivity(:,checkElemID), 4, pointsOfTriangle(1)) .and. &
                     numberInList(connectivity(:,checkElemID), 4, pointsOfTriangle(2)) .and. &
                     numberInList(connectivity(:,checkElemID), 4, pointsOfTriangle(3)) ) then

                    elemNeighborList(iSide, iElem) = checkElemID

                endif

            enddo
        enddo
    enddo

!***** correlate this neighbor list to the partitions these elements have to find the size of the boundary

    boundaryCount = 0

    do iElem = 1, nElems
        do iSide = 1,4
            ! first check whether there actually is a neighbor,
            ! or whether the element is on the edge of the mesh.
            if (elemNeighborList(iSide, iElem) .ne. 0) then
                if (elementPartition(iElem) .ne. &
                    elementPartition(elemNeighborList(iSide, iElem))) then
                    boundaryCount = boundaryCount + 1
                endif
            endif
        enddo
    enddo

    write(*,*) "Mesh partitioning has a boundary of: ", boundaryCount, " triangles"

end program

subroutine countlines(filename, nlines, filehandle)
    implicit none
    character(len=*)   :: filename
    integer            :: nlines, filehandle, ierr

    nlines = 0
    open(unit=filehandle, file=trim(filename), iostat=ierr)

    if (ierr.ne.0) then
        write(0,*) "partitionQuality wanted to count lines in ", trim(filename)
        write(0,*) "but could not open it, error", ierr
    endif


    do while (.true.)
        read (filehandle,*, end=10)
        nlines = nlines + 1
    enddo
10  close(filehandle)

end subroutine

logical function numberInList(list, listLength, number)

implicit none

integer :: listLength, number
integer :: list(listLength)

integer :: iList

do iList = 1,listLength
    if (list(iList) .eq. number) then
        numberInList = .true.
        return
    endif
enddo

numberInList = .false.
return

end function


