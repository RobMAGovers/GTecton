module types

implicit none

type point
    double precision :: x, y, z
    integer :: partition
    integer :: nNeighbors
    integer :: neighbors(14) ! never more than 14
end type point

type element
    integer :: v(4)
    integer :: partition
end type element

end module


program trialMesh3D

use types

implicit none


! This subroutine creates the tecin.dat.partf.elm en tecin.dat.partf.nps
! to test the parallel numbering scheme, and the neighbor search.

! We force our own partition, because a METIS partitioning
! is freuently irregular in shape, and needs to be drawn out by hand, 
! which is a lot of work.

! The mesh is rectangular, and start from the origin in the positive quadrant
! Each rectangle is divided into two triangles.

! The partitions are rectangular subdomains of size
! nCubesPerPartx X nCubesPerParty

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
double precision,  parameter :: dz = 0.1d0 ! side length of mesh rectangles

integer,           parameter :: nCubesPerPartx = 3 ! must be > 0
integer,           parameter :: nCubesPerParty = 3 ! must be > 0
integer,           parameter :: nCubesPerPartz = 3 ! must be > 0

integer,           parameter :: nPartitionsx  = 3 ! must be > 0
integer,           parameter :: nPartitionsy  = 2 ! must be > 0
integer,           parameter :: nPartitionsz  = 4 ! must be > 0

! unpartitioned files
character(len=19), parameter :: vertexFileU   = "tecin.dat.nps"
character(len=19), parameter :: elementFileU  = "tecin.dat.elm"

character(len=38), parameter :: vertexFormatU  = "(I13,1x,1e25.17,1x,1e25.17,1x,1e25.17)"
character(len=9),  parameter :: elementFormatU = "(6I12)"

! pratitioned files
character(len=19), parameter :: vertexFile    = "tecin.dat.partf.nps"
character(len=19), parameter :: elementFile   = "tecin.dat.partf.elm"

character(len=54), parameter :: vertexFormat  = "(I5,I13,I13,1x,1e25.17,1x,1e25.17,1x,1e25.17,I6,14I13)"
character(len=9),  parameter :: elementFormat = "(I5,7I12)"

double precision,  parameter :: eps           = 0.000001


! derived parameters
integer,           parameter :: nCubesX      = nCubesPerPartx * nPartitionsx
integer,           parameter :: nCubesY      = nCubesPerParty * nPartitionsy
integer,           parameter :: nCubesZ      = nCubesPerPartz * nPartitionsz

integer,           parameter :: nPointsX    = nCubesX + 1
integer,           parameter :: nPointsY    = nCubesY + 1
integer,           parameter :: nPointsZ    = nCubesZ + 1

integer,           parameter :: nPoints     = nPointsX * nPointsY * nPointsZ
integer,           parameter :: nCubes      = nCubesX * nCubesY * nCubesZ
integer,           parameter :: nElements   = 6 * nCubes
integer,           parameter :: nPartitions = nPartitionsx * nPartitionsy * nPartitionsz

integer :: iPartX,  iPartY, iPartZ   ! loop index
integer :: iPointX, iPointY, iPointZ  ! loop index
integer :: iPoint, iElem     ! loop index
integer :: iCubesX, iCubesY, iCubesZ
integer :: iPart

integer :: stepX, stepY, stepZ

! helpr variable to find the proper partitions
integer :: xthBlock, ythBlock, zthBlock, thisPart

! helper variables to find vertices of each element
integer :: topLeftPoint 
integer :: topRightPoint
integer :: bottomLeftPoint
integer :: bottomRightPoint

! helper variables for partition.info
integer, allocatable :: PointsPerPartition(:)
integer, allocatable :: ElementsPerPartition(:)

integer              :: lll, hll, lhl, hhl, llh, hlh, lhh, hhh
double precision     :: maxx, maxy, maxz

character(len=1)   :: numberString

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
do iPointZ=1,nPointsZ
    do iPointY=1,nPointsY
        do iPointX=1,nPointsX
            iPoint = iPoint + 1
            ! coordinates
            points(iPoint)%x = dx * dble(iPointX-1)
            points(iPoint)%y = dy * dble(iPointY-1)
            points(iPoint)%z = dz * dble(iPointZ-1)
            ! partition
            xthBlock = ceiling(abs(dble(iPointX-1)/dble(nCubesPerPartx) - eps))
            ythBlock = ceiling(abs(dble(iPointY-1)/dble(nCubesPerParty) - eps))
            zthBlock = ceiling(abs(dble(iPointZ-1)/dble(nCubesPerPartz) - eps))

            points(iPoint)%partition = nPartitionsx * nPartitionsy * (zthBlock-1) + nPartitionsx * (ythBlock-1) + xthBlock - 1
            PointsPerPartition(points(iPoint)%partition+1) = PointsPerPartition(points(iPoint)%partition+1) + 1
        enddo
    enddo
enddo

! set neighbors. Do this per direction.
! Either way, we would need 6 different cases.
! Might as well use the most transparent one

! explicitly set nNeighbors to zero
do iPoint=1,nPoints
    points(iPoint)%nNeighbors = 0
enddo


stepX = 1
stepY = nPointsX
stepZ = nPointsX * nPointsY

maxX = dx * nCubesX
maxY = dy * nCubesY
maxZ = dz * nCubesZ


do iPoint = 1, nPoints

! orthogonal neighbors
    if (points(iPoint)%x .gt. eps) then
        ! there is a neighbor in the negative X direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint - stepX
    endif

    if (points(iPoint)%x .lt. maxX - eps) then
        ! there is a neighbor in the positive X direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint + stepX
    endif

    if (points(iPoint)%y .gt. eps) then
        ! there is a neighbor in the negative Y direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint - stepY
    endif

    if (points(iPoint)%y .lt. maxY - eps) then
        ! there is a neighbor in the positive Y direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint + stepY
    endif

    if (points(iPoint)%z .gt. eps) then
        ! there is a neighbor in the negative Z direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint - stepZ
    endif

    if (points(iPoint)%z .lt. maxZ - eps) then
        ! there is a neighbor in the positive Z direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint + stepZ
    endif

! 2D diagonal neighbors

    if (points(iPoint)%y .gt. eps .and. points(iPoint)%z .gt. eps) then
        ! there is a neighbor in the negative Y,Z direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint - stepY - stepZ
    endif

    if (points(iPoint)%y .lt. maxY - eps .and. points(iPoint)%z .lt. maxZ - eps) then
        ! there is a neighbor in the positive Y,Z direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint + stepY + stepZ
    endif


    if (points(iPoint)%x .gt. eps .and. points(iPoint)%z .gt. eps) then
        ! there is a neighbor in the negative X,Z direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint - stepX - stepZ
    endif

    if (points(iPoint)%x .lt. maxX - eps .and. points(iPoint)%z .lt. maxZ - eps) then
        ! there is a neighbor in the positive X,Z direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint + stepX + stepZ
    endif


    if (points(iPoint)%x .gt. eps .and. points(iPoint)%y .gt. eps) then
        ! there is a neighbor in the negative X,Y direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint - stepX - stepY
    endif

    if (points(iPoint)%x .lt. maxX - eps .and. points(iPoint)%y .lt. maxY - eps) then
        ! there is a neighbor in the positive X,Y direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint + stepX + stepY
    endif

! 3D diagonal neighbors

    if (points(iPoint)%x .gt. eps .and. points(iPoint)%y .gt. eps .and. points(iPoint)%z .gt. eps) then
        ! there is a neighbor in the negative X,Y,Z direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint - stepX - stepY - stepZ
    endif

    if (points(iPoint)%x .lt. maxX - eps .and. points(iPoint)%y .lt. maxY - eps .and. points(iPoint)%z .lt. maxZ - eps) then
        ! there is a neighbor in the positive X,Y,Z direction
        points(iPoint)%nNeighbors = points(iPoint)%nNeighbors + 1
        points(iPoint)%neighbors(points(iPoint)%nNeighbors) = iPoint + stepX + stepY + stepZ
    endif

enddo

open(unit=110, file=vertexFile) 
do iPoint=1,nPoints
    ! construct the format based on the number of neighbors
    write(110, vertexFormat) points(iPoint)%partition, &
                            iPoint, &
                            1, &
                            points(iPoint)%x, &
                            points(iPoint)%y, &
                            points(iPoint)%z, &
                            points(iPoint)%nNeighbors, &
                            points(iPoint)%neighbors
enddo
write(110,'(a)') "end"
close(110)


open(unit=120, file=vertexFileU)
do iPoint=1,nPoints
    ! construct the format based on the number of neighbors
    write(120, vertexFormatU) iPoint, &
                             points(iPoint)%x, &
                             points(iPoint)%y, &
                             points(iPoint)%z
enddo
write(120,'(a)') "end"
close(120)


!------------- Elements -----------------

allocate(elements(nElements))
allocate(elementsPerPartition(nPartitions))

do iPart=1,nPartitions
    elementsPerPartition(iPart) = 0
enddo


! do this in the order of the rectangles, and for everty
! rectangle do the two triangles separately.

iElem=0
do iCubesZ = 1,nCubesZ
    do iCubesY = 1,nCubesY
        do iCubesX = 1,nCubesX

            ! the points are marked with three letter combinations, l(ow) or h(igh),
            ! the three letter symbolzis xyz

            lll = (iCubesX-1) + (iCubesY-1) * nPointsX + (iCubesZ-1) * nPointsX * nPointsY + 1
            hll = lll + stepX
            lhl = lll + stepY
            llh = lll + stepZ
            hhl = hll + stepY
            hlh = hll + stepZ
            lhh = lhl + stepZ
            hhh = lhh + stepX

            xthBlock = floor(dble(iCubesX-1)/dble(nCubesPerPartx)+eps)
            ythBlock = floor(dble(iCubesY-1)/dble(nCubesPerParty)+eps)
            zthBlock = floor(dble(iCubesZ-1)/dble(nCubesPerPartz)+eps)

            thisPart = nPartitionsx * nPartitionsy * zthBlock + nPartitionsx * ythBlock + xthBlock        

            elementsPerPartition(thisPart+1) = elementsPerPartition(thisPart+1) + 6

        
            ! add six elements of a cube

            iElem = iElem + 1
            elements(iElem)%partition = thisPart
            elements(iElem)%v(1) = lll
            elements(iElem)%v(2) = llh
            elements(iElem)%v(3) = hlh
            elements(iElem)%v(4) = hhh

            iElem = iElem + 1
            elements(iElem)%partition = thisPart
            elements(iElem)%v(1) = lll
            elements(iElem)%v(2) = llh
            elements(iElem)%v(3) = lhh
            elements(iElem)%v(4) = hhh

            iElem = iElem + 1
            elements(iElem)%partition = thisPart
            elements(iElem)%v(1) = lll
            elements(iElem)%v(2) = lhl
            elements(iElem)%v(3) = hhl
            elements(iElem)%v(4) = hhh

            iElem = iElem + 1
            elements(iElem)%partition = thisPart
            elements(iElem)%v(1) = lll
            elements(iElem)%v(2) = lhl
            elements(iElem)%v(3) = lhh
            elements(iElem)%v(4) = hhh

            iElem = iElem + 1
            elements(iElem)%partition = thisPart
            elements(iElem)%v(1) = lll
            elements(iElem)%v(2) = hll
            elements(iElem)%v(3) = hhl
            elements(iElem)%v(4) = hhh

            iElem = iElem + 1
            elements(iElem)%partition = thisPart
            elements(iElem)%v(1) = lll
            elements(iElem)%v(2) = hll
            elements(iElem)%v(3) = hlh
            elements(iElem)%v(4) = hhh

        enddo
    enddo
enddo


open(unit=111, file=elementFile)
do iElem=1,nElements
    ! construct the format based on the number of neighbors
    write(111, elementFormat) elements(iElem)%partition, &
                              iElem, &
                              1, &
                              elements(iElem)%v
enddo
write(111,'(a)') "end"
close(111)

open(unit=121, file=elementFileU)
do iElem=1,nElements
    ! construct the format based on the number of neighbors
    write(121, elementFormatU) iElem, &
                               1, &
                               elements(iElem)%v
enddo
write(121,'(a)') "end"
close(121)


call writeVTU(nPoints, points, nElements, elements)



deallocate(elements)

!----------- write a dummy TECIN file with the correct number of elements and vertices

open(unit=112, file="dummyTECIN.dat")


write(112,"(a)") "Dummy TECIN to let plnplt read in the data"
write(112,"(3i12)") nPoints, nElements, 1
write(112,"(a)") "    2    0    0    0    0    0    0    0    0    0"
write(112,"(a)") "    0"
write(112,"(a)") ".so tecin.dat.partf.nps"
write(112,"(a)") ".so tecin.dat.partf.elm"
write(112,"(a)") "           1    1    1"
write(112,"(a)") "           2    1    1"
write(112,"(a)") "end"
write(112,"(a)") "end"
write(112,"(a)") "end"
write(112,"(a)") "end"
write(112,"(a)") "end"
write(112,"(a)") "    0"
write(112,"(a)") "    0"
write(112,"(a)") "  sec"
write(112,"(a)") "  0.5"
write(112,"(a)") "    0    1    1    2    0    0    0    0    0    0"
write(112,"(a)") "    0"
write(112,"(a)") "           0           0           0           0           0           0           0           0"
write(112,"(a)") "           1        5.0E10           0.3  1.000000e+30           1.0           1.0           1.0  material 1"
write(112,"(a)") "end material properties"
write(112,"(a)") "           0.0           0.0           0.0"

close(112)




open(unit=113, file="partition.info")

write(113,"(i5)") nPartitions
do iPart=1,nPartitions
    write(113,"(i5,2i12)") iPart-1, PointsPerPartition(iPart), elementsPerPartition(iPart)

enddo


close(13)

end program

!-------------------------------------------------------------------------------------------

subroutine writeVTU(nPoints, points, nElements, elements)

use types

implicit none

integer, parameter :: VTKfile = 72

integer :: nPoints, iPoint, nElements, iElement, ierr

type(point)   :: points(nPoints)
type(element) :: elements(nElements)

character(len=12) nPoint_char, nElem_char
integer :: pointID

write (nPoint_char, *) 4 * nElements
write (nElem_char, *) nElements

open(unit=VTKfile, file="mesh3D.vtu")

write(VTKfile,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
write(VTKfile,*) '   <UnstructuredGrid>'
write(VTKfile,*) '     <Piece NumberOfPoints="'//nPoint_char//'" NumberOfCells="'//nElem_char//' ">'
!############# write the partitioning of the points
write(VTKfile,*) '         <PointData Scalars="scalars">'
write(VTKfile,*) '            <DataArray type="UInt32" NumberOfComponents="1" Name="vertex partitioning" Format="ascii">'

! print the nodes per element.
do iElement = 1, nElements
     write(VTKfile,*) points(elements(iElement)%v(1))%partition, &
                      points(elements(iElement)%v(2))%partition, &
                            points(elements(iElement)%v(3))%partition, &
                            points(elements(iElement)%v(4))%partition

enddo

write(VTKfile,*) '            </DataArray>'
write(VTKfile,*) '         </PointData>'
write(VTKfile,*) '         <CellData Scalars="scalars">'
write(VTKfile,*) '            <DataArray type="UInt32" NumberOfComponents="1" Name="Element partitioning" Format="ascii">'

! each thread writes its own actual element partitioning

do iElement = 1, nElements
    write(VTKfile,*) elements(iElement)%partition
enddo

write(VTKfile,*) '              </DataArray>'
write(VTKfile,*) '          </CellData>'
write(VTKfile,*) '          <Points>'
write(VTKfile,*) '             <DataArray type="Float32" NumberOfComponents="3" Name="Vertex coordinates" Format="ascii">'

! print the nodes per element.
do iElement = 1, nElements
    do iPoint = 1, 4
        pointID = elements(iElement)%v(iPoint)
        write(VTKfile,*) points(pointID)%x, &
                         points(pointID)%y, &
                         points(pointID)%z
    enddo
enddo

write(VTKfile,*) '              </DataArray>'
write(VTKfile,*) '           </Points>'
write(VTKfile,*) '           <Cells>'
write(VTKfile,*) '              <DataArray type="UInt64" Name="connectivity" Format="ascii">'

    ! connectivity is not the actual connectivity of the mesh,
    ! but the separate connectivity for each element,
    ! in order to plot a sharp boundary between the different elements.
do iElement = 1, nElements
    write(VTKfile,*) 4 * iElement - 4, &
                     4 * iElement - 3, &
                     4 * iElement - 2, &
                     4 * iElement - 1
enddo

write(VTKfile,*) '              </DataArray>'
write(VTKfile,*) '              <DataArray type="UInt32" Name="offsets" Format="ascii">'

do iElement = 1, nElements
    write(VTKfile,*) 4 * iElement
enddo

write(VTKfile,*) '             </DataArray>'
write(VTKfile,*) '             <DataArray type="Int32" Name="types" Format="ascii">'
! element type... 10 is tetrahedron; see https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf, page 9
do iElement = 1, nElements
    write(VTKfile,*) 10
enddo


write(VTKfile,*) '             </DataArray>'
write(VTKfile,*) '          </Cells>'
write(VTKfile,*) '       </Piece>'
write(VTKfile,*) '    </UnstructuredGrid>'
write(VTKfile,*) ' </VTKFile>'

close(VTKfile)


end subroutine


