subroutine matchElements(outputFormat)

use meshdata
use enumerates
use moduleMPI, only: MPI_INT, &
                     MPI_COMM_WORLD, &
                     MPI_STATUS_IGNORE, &
                     MPI_DOUBLE, &
                     nProcessors, &
                     myProcessorID

implicit none

integer            :: outputFormat

integer            :: sendToThread
integer            :: receiveFromThread

integer            :: nPointsMatched

integer            :: iPoint, iElement, iErr, iPass
integer            :: pointID, partitionID



allocate(partitionOfElements(nElementsLocal))
! initialize to -1 (instead of 0) for debugging purposes.
partitionOfElements = -1


allocate(partitionOfElementPoints(nPointsPerElement,nElementsLocal))

! VTK files require a connectivity that each point unique per element.
! So every point needs tbe there multiple times. Once per element it takes part in.
allocate(coordsForMyElems(nDimensions, nPointsPerElement, nElementsLocal))


! initialize to -1 (instead of 0) for debugging purposes.
! easier to see whether a partition has not been assigned,
! then the doubt whether a partition 0 exists or not.
partitionOfElementPoints = -1

! for pick out the ones they already have
minNodeThisPartition = vertexOffset(myProcessorID+1) + 1
maxNodeThisPartition = vertexOffset(myProcessorID+1) + nPointsLocal

pointCoords => pointCoordsA
connectivity => connectivityA
partitionOfPoints => partitionOfPointsA

do iElement = 1,nElementsLocal
    do iPoint = 1,nPointsPerElement

        pointID = connectivity(iPoint, iElement)

        if ((pointID .ge. minNodeThisPartition) .and. &
            (pointID .le. maxNodeThisPartition)) then

            if (outputFormat.eq.outputFormatGMT .or. &
                outputFormat.eq.outputFormatVTK) then

                coordsForMyElems(:, iPoint, iElement) = &
                pointCoords(:, pointID - minNodeThisPartition + 1)
            endif


            partitionOfElementPoints(iPoint, iElement) = &
              partitionOfPoints(pointID - minNodeThisPartition + 1)

        endif
    enddo
enddo

! and then look at all the points that the other processors have.

! The data is passed around the processors between the array and the array2.
! These are: partitionOfPoints and partitionOfPoints2
!            pointCoords       and pointCoords2

! This is done in two phases. First, all even processors send, and the odd ones receive.
! After that, the odd ones send and the even ones receive.
! To ensure that existing data is not overwritten:

! on odd  phases: the data is moved from arrayA to arrayB
! on even phases: the data is moved from arrayB to arrayA

!write(*,*) "1thread", myProcessorID, " coordA: ", pointCoordsA
!write(*,*) "1thread", myProcessorID, " coordB: ", pointCoordsB

allocate(partitionOfPointsB(nPointsLocal0))
partitionOfPointsB = 0

if (outputFormat.eq.outputFormatGMT .or. &
    outputFormat.eq.outputFormatVTK) then
!    allocate(pointCoords2(nDimensions, nPointsLocal0))
    pointCoordsB = 0d0
endif


do iPass = 1, nProcessors-1

	call passData(nPointsLocal0, &
                  partitionOfPointsA, &
                  partitionOfPointsB, &
                  iPass)

    if (modulo(iPass,2) .eq. 1) then
        partitionOfPoints => partitionOfPointsB
    else
        partitionOfPoints => partitionOfPointsA
    endif

    ! these are the limits from the thread from which the data was received
    minNodeThisPartition = vertexOffset(modulo(nProcessors + myProcessorID - iPass, nProcessors)+1)+1
    if (iPass .eq. myProcessorID) then
        ! this thread now has the zero-partition
        maxNodeThisPartition = minNodeThisPartition + nPointsLocal0-1
    else
        maxNodeThisPartition = minNodeThisPartition + nBase-1
    endif

    if (outputFormat.eq.outputFormatGMT .or. &
        outputFormat.eq.outputFormatVTK) then

        ! GMT requires the coordinates of the vertices to plot the elements.
        ! As such, we need to send them around as well

        ! For the staggered sending and receiving, we use the same method
        ! as described above for the point partitions

	    call passData(nPointsLocal0 * nDimensions, &
    	              pointCoordsA, &
        	          pointCoordsB, &
            	      iPass)

	    if (modulo(iPass,2) .eq. 1) then
    	    pointCoords => pointCoordsB
	    else
    	    pointCoords => pointCoordsA
	    endif

    endif
    ! endif outputformatGMT


!write(*,*) "3thread", myProcessorID, "pass", iPass, " coordA: ", pointCoordsA
!write(*,*) "3thread", myProcessorID, "pass", iPass, " coordB: ", pointCoordsB



    ! now that the data has been transported, we can look
    ! and compare points.
    ! We compare with array2 on odd passes,
    ! and with array on even passes

    nPointsMatched = 0

    ! odd pass, the data is now in array2, and not in array1


    ! now that all the data has been exchanched and established,
    ! we can check for matches.
    do iElement = 1,nElementsLocal
        do iPoint = 1,nPointsPerElement

            pointID = connectivity(iPoint, iElement)

            if ((pointID .ge. minNodeThisPartition) .and. &
                (pointID .le. maxNodeThisPartition)) then

                nPointsMatched = nPointsMatched + 1

                if (outputFormat.eq.outputFormatGMT .or. &
                    outputFormat.eq.outputFormatVTK) then

                    coordsForMyElems(:, iPoint, iElement) = &
                    pointCoords(:, pointID - minNodeThisPartition + 1)

                endif

                partitionOfElementPoints(iPoint, iElement) = &
                    partitionOfPoints(pointID - minNodeThisPartition + 1)

            endif
        enddo
    enddo




    ! reset receiving arrays, just for robustness
!    if (modulo(iPass,2) .eq. 1) then
!        partitionOfPointsA = 0
!        pointCoordsA       = 0d0
!    else
!        partitionOfPointsB = 0
!        pointCoordsB	   = 0d0
!	endif


enddo

! now partitionOfElementPoints has been filled, and for each element we take the lowesy value to number the elements

do iElement = 1,nElementsLocal

	partitionID = minval(partitionOfElementPoints(:, iElement))

	partitionOfElements(iElement) = partitionID

    ElementCountPerPartitionThread(partitionID) = &
    ElementCountPerPartitionThread(partitionID) + 1

enddo

! the point partitions need to brought back to their home partition, prior to being written to file

call passData(nPointsLocal0, &
              partitionOfPointsA, &
              partitionOfPointsB, &
              iPass)

end subroutine
