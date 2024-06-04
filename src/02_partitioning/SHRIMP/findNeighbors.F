subroutine findNeighbors


use meshdata

use moduleMPI, only: myProcessorID, &
                     nProcessors, &
                     MPI_COMM_WORLD, &
                     MPI_INT, &
                     MPI_STATUS_IGNORE
                     


implicit none

integer            :: sendToThread
integer            :: receiveFromThread

integer            :: blockSize
integer            :: highestCheckedIndex
integer, allocatable :: neighborListLength(:)
logical            :: neighborisknown
integer            :: nAllocation
integer            :: nBlocksTaken
integer            :: neighborof


logical, external  :: numberInList

integer, allocatable :: tempNeighbors(:)

integer            :: iPoint, iElement, iPass, iErr, iNeighbor
integer            :: pointID, localPointID, neighborID


allocate(neighborCount(nPointsLocal))
neighborCount = 0

allocate(neighborlist(nPointsLocal))
do iPoint=1,nPointsLocal
    allocate(neighborlist(iPoint)%neighbors(neighborBase))
    neighborlist(iPoint)%neighbors = 0
enddo

allocate(neighborListLength(nPointsLocal))
neighborListLength = 0

! these are the limits from the thread
minNodeThisPartition = vertexOffset(myProcessorID+1)+1
maxNodeThisPartition = minNodeThisPartition + nPointsLocal-1

allocate(connectivityB(nPointsPerElement,nElementsLocal0))
connectivityB = 0

do iPass = 1, nProcessors

    call passData(nPointsPerElement * nElementsLocal0, &
                  connectivityA, &
                  connectivityB, &
                  iPass)


    ! the elements have been moved. Compare the connectivity with existing points.
    ! This partition has nPointsLocal points, and start count at vertexOffset(myProcessorID+1)
    ! The connectivity can have any number.

    ! The simple solution would be to use a 2D array, nPoints x maxNeighbors, and
    ! gradually reallocate it larger as more neighbors are found. However, experience taught
    ! that almost all nodes have a fairly moderate amounf of neighbors (typically 6
    ! for triangles in 2D, or 14-ish voor tetrahedrons) while there are on occasion a few nodes
    ! that stand out due to unfortunate meshing. Examples have been seen with 32 neighbors
    ! in 2D, or over a 100 in 3D. As this array becomes of considerable size, 
    ! it is worth to be creative to keep it's size limited, contrary to allocating more than
    ! a hundred neighbors for each and every node.

    ! It can in theory have an endless amount of neighors for a specific point, without
    ! accomodating space for other points as well.

    ! we have an array neighborList allocated, that has blocks of neighborBase +
    ! 10% expansion room. If this room is full, more space is allocated
    ! all these blocks are off 

    ! set a pointer to the array that is relevant for this stage.
    if (modulo(iPass,2) .eq. 1) then
        connectivity => connectivityB
    else
        connectivity => connectivityA
    endif


    do iElement = 1,nElementsLocal0

        do iPoint = 1, nPointsPerElement
            pointID = connectivity(iPoint, iElement)

            if ((pointID .ge. minNodeThisPartition) .and. &
                (pointID .le. maxNodeThisPartition)) then

                ! We have a node in our partition!
                ! All the other nodes in this element are neighbors.
                ! Later on, when we implement hexahedron elements,
                ! this is no longer the case, but that is tomorrows worry.
                localPointID = pointID - minNodeThisPartition + 1
                do iNeighbor = 1,nPointsPerElement
                    neighborID = connectivity(iNeighbor, iElement)
                    ! make sure that the point is not added as a neighbor of itself
                    if (neighborID .ne. pointID) then
                        ! check whether the node is already there
                        if (.not. numberInList(neighborlist(localPointID)%neighbors, &
                                               neighborListLength(localPointID), &
                                               neighborID)) then
                            ! add the node. But we must think about where to add it!

                            ! if the array is full, allocate more space
                            if (neighborCount(localPointID) .eq. neighborListLength(localPointID)) then

                                ! create a temporary array and place the known neighbors there.
                                allocate(tempNeighbors(neighborCount(localPointID)))
                                tempNeighbors = neighborlist(localPointID)%neighbors

                                ! remove the old array that is too small, and allocate with bigger size and set to 0
                                deallocate(neighborlist(localPointID)%neighbors)
                                  allocate(neighborlist(localPointID)%neighbors(neighborListLength(localPointID) + &
                                           neighborBase))
                                neighborlist(localPointID)%neighbors = 0

                                ! put the known nodes back
                                neighborlist(localPointID)%neighbors(1:neighborListLength(localPointID)) = tempNeighbors

                                ! increased the known size of the list
                                neighborListLength(localPointID) = neighborListLength(localPointID) + neighborBase

                                ! and finally take out the trash
                                deallocate(tempNeighbors)

                            endif

                            ! now we have space, we can add the node
                            neighborCount(localPointID) = neighborCount(localPointID) + 1
                            neighborlist(localPointID)%neighbors(neighborCount(localPointID)) = neighborID

                        endif


                    endif

                enddo

            endif

        enddo
    enddo
enddo

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

