! This program reads a *partitioned* gmsh file,
! and writes the tecin.dat.partf.nps
! and the tecin.dat.partf.elm
! in the format suitable for GTECTON


module mesh

type meshpoint
integer          :: partition
double precision :: x,y,z
integer          :: marker
integer          :: nNeighbors
integer, allocatable :: neighborList(:)
end type meshpoint

type element
integer          :: partition
integer          :: vertices(4)
integer          :: marker
integer          :: material
end type element

type(meshpoint), allocatable :: meshpoints(:)
type(element),   allocatable ::    elements(:)

end module

!-------------------------------------

program main


use mesh, only: meshpoints, elements

implicit none

character(len=256) :: gmshFileName

integer, parameter :: tectonPointFile = 42
integer, parameter :: tectonElementFile = 43
integer, parameter :: partitionFile = 44

integer            :: nElements, iElem
integer            :: nPoints, iPoint
integer            :: iNeighbor
integer            :: nPartitions, iPartition

logical            :: verbose

integer, allocatable :: tempArray(:)

integer, allocatable :: pointsPerPartition(:)
integer, allocatable ::  elemsPerPartition(:)


integer            :: pointID, neighborID
logical            :: neighborIsUnique
integer            :: checkNode

integer            :: nCurrent, nNext

call readCommandLineArguments(gmshFileName, verbose)

!   read GMSH mesh data
call readmesh (gmshFileName, nPoints, nElements)

! build node neighbor list
! walk through the elements,
! and for every element, its nodes are each other's
! neighbors. Using these nodes, build up a list,
! but make that ndoes are nodes registered double,
! because they are each other's neighbor for multiple
! elements.


! start out with room for 10 neighors per point
do iPoint = 1,nPoints
    allocate(meshpoints(iPoint)%neighborList(10))
    meshpoints(iPoint)%neighborList = 0
enddo

! walk through all the elements, 
! and monitor the neighbors
do iElem = 1,nElements
    do iPoint = 1,4
        pointID = elements(iElem)%vertices(iPoint)
        do iNeighbor = 1,4
            if (iPoint .ne. iNeighbor) then
                neighborID = elements(iElem)%vertices(iNeighbor)

                ! see whether this is the first neighbor node.
                ! if so, we must check for uniqueness.
                ! Otherwise we can add it without further thinking.

                if (meshpoints(pointID)%nNeighbors .gt. 0) then
                    ! check whether this neighbor has not yet been 
                    ! registered from a connection with an earlier element
                    neighborIsUnique = .true.
                    checkNode = 0
                    do while (neighborIsUnique .and. checkNode.lt.meshpoints(pointID)%nNeighbors)
                        checkNode = checkNode + 1
                        if (meshpoints(pointID)%neighborList(checkNode) .eq. neighborID) then
                             neighborIsUnique = .false.
                        endif
                    enddo

                    if (neighborIsUnique) then
                        ! check whether the neighbor array is full.
                        if (meshpoints(pointID)%neighborList(size( meshpoints(pointID)%neighborList)) .ne. 0 ) then

                             ! the last entry of the array is a non-zero
                             ! Make the array 10 longer

                             nCurrent = meshpoints(pointID)%nNeighbors
                             nNext = nCurrent + 10

                             allocate(tempArray(nCurrent))
                             tempArray = meshpoints(pointID)%neighborList

                             deallocate(meshpoints(pointID)%neighborList)
                             allocate(meshpoints(pointID)%neighborList(nNext))
                             meshpoints(pointID)%neighborList = 0
                             meshpoints(pointID)%neighborList(1:nCurrent) = tempArray
                             deallocate(tempArray)

                         endif    

                         ! the array still had room, or now has extra room.
                         ! add the new neighbor
                         meshpoints(pointID)%nNeighbors = &
                         meshpoints(pointID)%nNeighbors + 1
                         
                         meshpoints(pointID)%neighborList(meshpoints(pointID)%nNeighbors) = neighborID
                    endif

                else
                    ! this is the first neighbor
                    meshpoints(pointID)%nNeighbors = 1
                    meshpoints(pointID)%neighborList(1) = neighborID
                endif
            endif
        enddo
    enddo
enddo

! determine the partition number of the nodal points, 
! by selecting the partition of an element it is part of

do iElem = 1,nElements
    do iPoint = 1,4
        pointID = elements(iElem)%vertices(iPoint)
        if (meshpoints(pointID)%partition .eq. 0) then
            ! Encountered an unpartitioned node. 
            ! assign partition from adjacent element.
            meshpoints(pointID)%partition = &
            elements(iElem)%partition
        else
            ! the nodal point was already assigned the partition from another node
        endif
    enddo
enddo

!******** generate tecin.dat.partf.nps *************

open(unit=tectonPointFile,   file="tecin.dat.partf.nps")

do iPoint = 1, nPoints

    ! we compute partitions starting at 1, but gtecton wants them to start at 0, 
    ! hence the -1
    write(tectonPointFile, "(I5,1X,I12,1X,I12,1X,1E25.17,1X,1E25.17,1X,1E25.17,1X,I4)", advance='no') &
      meshpoints(iPoint)%partition-1, &
      iPoint, &
      1, & ! add dummy marker. as gmsh uses multiple markers per point and this can get confusing.
      meshpoints(iPoint)%x, &
      meshpoints(iPoint)%y, &
      meshpoints(iPoint)%z, &
      meshpoints(iPoint)%nNeighbors
    do iNeighbor = 1,meshpoints(iPoint)%nNeighbors-1
        write(tectonPointFile,"(1X,I12)",advance='no') &
        meshpoints(iPoint)%neighborList(iNeighbor)
    enddo
    write(tectonPointFile,"(1X,I12)") &
        meshpoints(iPoint)%neighborList(meshpoints(iPoint)%nNeighbors)
enddo

write(tectonPointFile,"(a)") "end partitioned nodal point file"


close(tectonPointFile)

!******** generate tecin.dat.partf.elm *************

open(unit=tectonElementFile, file="tecin.dat.partf.elm")

do iElem = 1, nElements
    ! we compute partitions starting at 1, but gtecton wants them to start at 0, 
    ! hence the -1
    write(tectonElementFile,"(i5,6i12)") &
      elements(iElem)%partition-1, &
      iElem, &
      elements(iElem)%material, &
      elements(iElem)%vertices(1:4)

enddo

write(tectonElementFile,"(a)") "end partitioned element file"

close(tectonElementFile)

!******** generate partition.info ****************

nPartitions = maxval(meshpoints(:)%partition)

allocate(pointsPerPartition(nPartitions))
allocate( elemsPerPartition(nPartitions))

pointsPerPartition = 0
elemsPerPartition = 0

do iPoint = 1, nPoints
    pointsPerPartition(meshpoints(iPoint)%partition) = &
    pointsPerPartition(meshpoints(iPoint)%partition) + 1
enddo

do iElem = 1, nElements
    elemsPerPartition(elements(iElem)%partition) = &
    elemsPerPartition(elements(iElem)%partition) + 1
enddo

open(unit=partitionFile, file="partition.info")

write(partitionFile, "(i5)") nPartitions
do iPartition = 1, nPartitions
    write(partitionFile, "(i5,2i12)") iPartition, &
                                      pointsPerPartition(iPartition), &
                                                   elemsPerPartition(iPartition)
enddo

close(partitionFile)


!***************************************************

write(*,*) "maximum neighboring nodes to another node:", maxval(meshpoints(:)%nNeighbors)


end program
!-------------------------------------------------------------------------------


subroutine readmesh(gmshFileName, nPoints, iTetrahedron)

use mesh

implicit none

integer :: gmshFileHandle
integer :: iTetrahedron, nPoints

character(len=256) :: string, record, gmshFileName

integer :: dummy, one
integer :: iPoint, iElement, nElements
integer :: elementType

open(unit=gmshFileHandle, file=gmshFileName)


! read first line, jusy says '$MeshFormat'
read(gmshFileHandle,'(a11)') string

if (string(1:11).ne.'$MeshFormat') then
     write(*,*) "header should start with $MeshFormat, but starts with", string(1:11)
     stop "does not appear to be a valid gmsh file"
endif

! read version string
read(gmshFileHandle,'(a7)') string

if (string(1:3) .ne. '2.2') then
    write(*,*)    "Expecting Gmsh format version 2.2. Got version ", string(1:3), " instead"
    stop "Leaving..."
endif

! read $EndMeshFormat
read(gmshFileHandle,'(a14)') string
!--------------------------------------------------------------------
! read $Nodes
read(gmshFileHandle,'(a6)') string

! read number of mesh points
read(gmshFileHandle,*) nPoints

allocate(meshpoints(nPoints))

! read the actual points.
! 'dummy' is an incremental counter
do iPoint=1,nPoints
    read(gmshFileHandle,*) dummy, meshpoints(iPoint)%x, &
                                  meshpoints(iPoint)%y, &
                                  meshpoints(iPoint)%z
enddo

! read $EndNodes
read(gmshFileHandle,'(a9)') string
!--------------------------------------------------------------------
! read $Elements
read(gmshFileHandle,'(a9)') string

! read number of Elements
! Note that gmsh also considers triangles to be 
! elements, while we are only interested in the 
! tetrahedrons. This is indicated by the second integers of the row:

! the Gmsh constants
! (taken from: http://gmsh.info/doc/texinfo/gmsh.html)
!integer,          parameter :: elemIsPoint          = 15
!integer,          parameter :: elemIsLine           = 1
!integer,          parameter :: elemIsTriangle       = 2
!integer,          parameter :: elemIsTetrahedron    = 4


read(gmshFileHandle,*) nElements

allocate(elements(nElements))


! read the actual elements. The number if integers changes 
! for different types of elements, so we cannot read it directly.
! In stead, read it into a string and interpret it.
iTetrahedron = 0



do iElement = 1,nElements

    read(gmshFileHandle,"(a)") record
    read(record,*) dummy, elementType

    if (elementType.eq.4) then

        ! yay, we have a tetrahedron
        iTetrahedron =    iTetrahedron + 1

        read(record,*) dummy, &
                       elementType, &
                       dummy, &  ! number of markers
                       elements(iTetrahedron)%material, & ! for some reason all elements have a 1
                       elements(iTetrahedron)%marker, &
                       one, & ! for some reason all elements have a 1
                       elements(iTetrahedron)%partition, &
                       elements(iTetrahedron)%vertices(1:4)
    else
        ! line contains a non-tetrahedron element. Skip it
    endif


enddo

! now a part of the elements array is still empty, but we just don't care.

end subroutine





!-------------------------------------------------------------------------------
    subroutine readCommandLineArguments &
       (gmshFileName, verbose)
!
! process command line arguments (see man page)
!
    implicit none
!-pass
logical :: verbose
character(len=256) :: gmshFileName
!-locl
integer :: i
integer :: narg
character(len=256) :: string
!-init
    verbose = .false.
!
    narg = COMMAND_ARGUMENT_COUNT()

    i = 1

    do while (i.le.narg) 
        call get_command_argument(i,string)


        if (string(1:3).eq.'-i ') then

            if (i.eq.narg) then
                stop "Usage: gmsh2tecton -i meshfile"                
            endif

            i = i + 1
            call get_command_argument(i,string)

            gmshFileName = string

        else if (string(1:3).eq.'-v ') then
            verbose = .true.

        else

            write(*,*) ' unknown option ',trim(string)
            write(*,*)
            stop "Usage: gmsh2tecton -i meshfile"
        endif

        i = i + 1

    enddo

    if (narg.eq.0) then
        stop "Usage: gmsh2tecton -i meshfile"
    endif

!
    return
    end
