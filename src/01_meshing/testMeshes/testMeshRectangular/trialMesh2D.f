program trialMesh2D

implicit none


! This subroutine creates the tecin.dat.partf.elm en tecin.dat.partf.nps
! to test the parallel numbering scheme, and the neighbor search.

! We force our own partition, because a METIS partitioning
! is freuently irregular in shape, and needs to be drawn out by hand, 
! which is a lot of work.

! The mesh is rectangular, and start from the origin in the positive quadrant
! Each rectangle is divided into two triangles.

! The partitions are rectangular subdomains of size
! nRectPerPartx X nRectPerParty

! The nodal points belong to the partition of the element to the lower left of them
! Points on the x and y axis, which have no elements on their lower left, belong to
! to the partition of the element on their upper right.
! The two poor nodes that have neither elements on their upper right, not lower left,
! will belong the only element in which they occur ;-)

! partition count begins near the origin and goes up with increasing x.
! when we reach the end of the row, we move up in y and number another row, etc.

! Elements are numbered incrementing by rectangle in the same manner as the points.
! We could choose them to increment over the partitions, but then the
! global and local numbering would be the same. We need to test the difference.

! parameters
double precision,  parameter :: dx = 0.1d0 ! side length of mesh rectangles
double precision,  parameter :: dy = 0.1d0 ! side length of mesh rectangles

integer,           parameter :: nRectPerPartx = 3 ! must be > 0
integer,           parameter :: nRectPerParty = 3 ! must be > 0

integer,           parameter :: nPartitionsx  = 3 ! must be > 0
integer,           parameter :: nPartitionsy  = 2 ! must be > 0

! unpartitioned files
character(len=19), parameter :: vertexFileU   = "tecin.dat.nps"
character(len=19), parameter :: elementFileU  = "tecin.dat.elm"

character(len=38), parameter :: vertexFormatU  = "(I13,1x,1e25.17,1x,1e25.17)"
character(len=9),  parameter :: elementFormatU = "(6I12)"

! pratitioned files
character(len=19), parameter :: vertexFile    = "tecin.dat.partf.nps"
character(len=19), parameter :: elementFile   = "tecin.dat.partf.elm"

character(len=38), parameter :: vertexFormat  = "(I5,I13,1x,1e25.17,1x,1e25.17,I6,6I13)"
character(len=9),  parameter :: elementFormat = "(I5,6I12)"

double precision,  parameter :: eps           = 0.000001


! derived parameters
integer,           parameter :: nRectX      = nRectPerPartx * nPartitionsx
integer,           parameter :: nRectY      = nRectPerParty * nPartitionsy
integer,           parameter :: nPointsX    = nRectX + 1
integer,           parameter :: nPointsY    = nRectY + 1
integer,           parameter :: nPoints     = nPointsX * nPointsY 
integer,           parameter :: nRects      = nRectX * nRectY
integer,           parameter :: nElements   = 2 * nRects
integer,           parameter :: nPartitions = nPartitionsx * nPartitionsy

integer :: iPartX,  iPartY   ! loop index
integer :: iPointX, iPointY  ! loop index
integer :: iPoint, iElem     ! loop index
integer :: iRectX, iRectY
integer :: iPart

! helpr variable to find the proper partitions
integer :: xthBlock, ythBlock, thisPart

! helper variables to find vertices of each element
integer :: topLeftPoint 
integer :: topRightPoint
integer :: bottomLeftPoint
integer :: bottomRightPoint

! helper variables for partition.info
integer, allocatable :: PointsPerPartition(:)
integer, allocatable :: ElementsPerPartition(:)


character(len=1)   :: numberString


type point
    double precision :: x, y
    integer :: partition
    integer :: nNeighbors
    integer :: neighbors(6) ! never more than 6
end type point

type element
    integer :: v1, v2, v3
    integer :: partition
end type element


type(point),   allocatable :: points(:)
type(element), allocatable :: elements(:)

!------------- Vertices -----------------

! allocate space for the mesh
allocate(points(nPoints))
allocate(PointsPerPartition(nPartitions))

do iPart=1,nPartitions
    PointsPerPartition(iPart) = 0
enddo



! set point coordinates and partition
iPoint = 0
do iPointY=1,nPointsY
    do iPointX=1,nPointsX
        iPoint = iPoint + 1
        ! coordinates
        points(iPoint)%x = dx * dble(iPointX-1)
        points(iPoint)%y = dy * dble(iPointY-1)
        ! partition
        xthBlock = ceiling(abs(dble(iPointX-1)/dble(nRectPerPartx) - eps))
        ythBlock = ceiling(abs(dble(iPointY-1)/dble(nRectPerParty) - eps))
        points(iPoint)%partition = nPartitionsx * (ythBlock-1) + xthBlock - 1
        PointsPerPartition(points(iPoint)%partition+1) = PointsPerPartition(points(iPoint)%partition+1) + 1
    enddo
enddo

! set neighbors. Do this per direction.
! Either way, we would need 6 different cases.
! Might as well use the most transparent one

! explicitly set nNeighbors to zero
do iPoint=1,nPoints
    points(iPoint)%nNeighbors = 0
enddo

do iPoint = 1, nPoints - nPointsX
    !add top neighbors
    points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
    points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint + nPointsX
    ! top right
    if (mod(iPoint, nPointsX).ne.0) then
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint + nPointsX + 1
    endif
enddo

! right
do iPoint = 1, nPoints
    if (mod(iPoint, nPointsX).ne.0) then
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint + 1
    endif
enddo

! bottom
do iPoint = nPointsX + 1, nPoints
    points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
    points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint - nPointsX
    ! bottom left
    if (mod(iPoint-1, nPointsX).ne.0) then
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint - nPointsX - 1
    endif
enddo

! left
do iPoint = 1, nPoints           
    if (mod(iPoint-1, nPointsX).ne.0) then
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint - 1
    endif
enddo


open(unit=10, file=vertexFile) 
do iPoint=1,nPoints
    ! construct the format based on the number of neighbors
    write(10, vertexFormat) points(iPoint)%partition, &
                            iPoint, &
                            points(iPoint)%x, &
                            points(iPoint)%y, &
                            points(iPoint)%nNeighbors, &
                            points(iPoint)%neighbors(1:6)
enddo
close(10)

open(unit=20, file=vertexFileU)
do iPoint=1,nPoints
    ! construct the format based on the number of neighbors
    write(20, vertexFormatU) iPoint, &
                             points(iPoint)%x, &
                             points(iPoint)%y
enddo
close(20)



deallocate(points)

!------------- Elements -----------------

allocate(elements(nElements))
allocate(elementsPerPartition(nPartitions))

do iPart=1,nPartitions
    elementsPerPartition(iPart) = 0
enddo


! do this in the order of the rectangles, and for everty
! rectangle do the two triangles separately.

iElem=0
do iRectY = 1,nRectY
    do iRectX = 1,nRectX

        bottomLeftPoint  = (iRectX-1) + (iRectY-1) * nPointsX + 1
        bottomRightPoint = bottomLeftPoint + 1
        topLeftPoint     = bottomLeftPoint + nPointsX
        topRightPoint    = topLeftPoint + 1

        xthBlock = floor(dble(iRectX-1)/dble(nRectPerPartx)+eps)
        ythBlock = floor(dble(iRectY-1)/dble(nRectPerParty)+eps)
        
        thisPart = nPartitionsx * ythBlock + xthBlock        

        elementsPerPartition(thisPart+1) = elementsPerPartition(thisPart+1) + 2

        ! top left triangle of the rectangle
        iElem = iElem + 1
        elements(iElem)%partition = thisPart
        elements(iElem)%v1 = bottomLeftPoint
        elements(iElem)%v2 = topRightPoint
        elements(iElem)%v3 = topLeftPoint

        ! bottom right triangle of the rectangle
        iElem = iElem + 1
        elements(iElem)%partition = thisPart
        elements(iElem)%v1 = bottomLeftPoint
        elements(iElem)%v2 = bottomRightPoint
        elements(iElem)%v3 = topRightPoint

    enddo
enddo


open(unit=11, file=elementFile)
do iElem=1,nElements
    ! construct the format based on the number of neighbors
    write(11, elementFormat) elements(iElem)%partition, &
                             iElem, &
                             1, &
                             elements(iElem)%v1, &
                             elements(iElem)%v2, &
                             elements(iElem)%v3, &
                             elements(iElem)%v3
enddo
close(11)

open(unit=21, file=elementFileU)
do iElem=1,nElements
    ! construct the format based on the number of neighbors
    write(21, elementFormatU) iElem, &
                             1, &
                             elements(iElem)%v1, &
                             elements(iElem)%v2, &
                             elements(iElem)%v3, &
                             elements(iElem)%v3
enddo
close(21)


deallocate(elements)

!----------- write a dummy TECIN file with the correct number of elements and vertices

open(unit=12, file="dummyTECIN.dat")


write(12,"(a)") "Dummy TECIN to let plnplt read in the data"
write(12,"(3i12)") nPoints, nElements, 1
write(12,"(a)") "    2    0    0    0    0    0    0    0    0    0"
write(12,"(a)") "    0"
write(12,"(a)") ".so tecin.dat.partf.nps"
write(12,"(a)") "end vertices"
write(12,"(a)") ".so tecin.dat.partf.elm"
write(12,"(a)") "end elements"
write(12,"(a)") "           1    1    1"
write(12,"(a)") "           2    1    1"
write(12,"(a)") "end"
write(12,"(a)") "end"
write(12,"(a)") "end"
write(12,"(a)") "end"
write(12,"(a)") "end"
write(12,"(a)") "    0"
write(12,"(a)") "    0"
write(12,"(a)") "  sec"
write(12,"(a)") "  0.5"
write(12,"(a)") "    0    1    1    2    0    0    0    0    0    0"
write(12,"(a)") "    0"
write(12,"(a)") "           0           0           0           0           0           0           0           0"
write(12,"(a)") "           1        5.0E10           0.3  1.000000e+30           1.0           1.0           1.0  material 1"
write(12,"(a)") "end material properties"
write(12,"(a)") "           0.0           0.0           0.0"

close(12)




open(unit=13, file="partition.info")

write(13,"(i5)") nPartitions
do iPart=1,nPartitions
    write(13,"(i5,2i12)") iPart-1, PointsPerPartition(iPart), elementsPerPartition(iPart)

enddo


close(13)

end program

