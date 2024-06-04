subroutine readInput(inputFormat)

use enumerates

implicit none


integer          :: inputFormat

if      (inputFormat .eq. inputFormatSimple) then
	call readInputSimple()
else if (inputFormat .eq. inputFormatGTECTON) then
    call readInputGTECTON()
else
	write(*,*) "inputFormat not recognized"
endif


end subroutine


! when adding additional input readers, make sure that the following variables are set:
! nPointsGlobal    number of points in total
! nPointsLocal     number of points on all partitions beside the 1st one
! nPointsLocal0    local number of points on the zero partition, which can be different from other partitions
! vertexOffset     offset where we start counting points on each partition
! nElementsGlobal  (see points)
! nElementslocal   (see points)
! nElementslocal0  (see points)
! elementOffset    (see points)
! pointCoordsA     actual coordinates
! connectivityA    actual connectivity

subroutine readInputSimple()

use files
use meshdata
use moduleMPI, only: myProcessorID, &
                     nProcessors, &
                     MPI_COMM_WORLD, &
                     MPI_INTEGER

implicit none

write(*,*) "simple input not yet implemented"

end subroutine



subroutine readInputGTECTON()

use files
use meshdata
use moduleMPI, only: myProcessorID, &
                     nProcessors, &
                     MPI_COMM_WORLD, &
                     MPI_INTEGER

implicit none

! local variables

integer :: ierr, iPoint, iElement, iProcessor
character(len=1) :: dummyChar ! to quickly read just the first space of a line,
                              ! to count the number of lines
integer          :: idx ! entries in tecin.dat.nps/elm
integer          :: elementcounterforthispartition
integer          ::   pointCounterForThisPartition

open(unit=40, file=pointFile)
open(unit=41, file=elementFile)

if (myProcessorID.eq.0) then

    nPointsGlobal = 0

    do while (.true.)
        read(40,*,end=10) dummyChar
        nPointsGlobal = nPointsGlobal + 1
    enddo
    10 rewind(40)
    nPointsGlobal = nPointsGlobal - 1 ! remove the line with the end statement from the count

    nElementsGlobal = 0

    do while (.true.)
        read(41,*,end=20) dummyChar
        nElementsGlobal = nElementsGlobal + 1
    enddo
    20 rewind(41)
    nElementsGlobal = nElementsGlobal - 1 ! remove the line with the end statement from the count

endif

! determine how much points each of the partitions will get

call MPI_Bcast( &
    nElementsGlobal, &  ! variable to be broadcast
    1, &                ! number of variables being broadcast
    MPI_INTEGER, &	! type
    0, &                ! sending thread
    MPI_COMM_WORLD, &   ! EVERYBODY!
    ierr)

if (ierr.ne.0) then
    write(*,*) "Could not broadcast number of elements, error: ", ierr
endif

call MPI_Bcast( &
    nPointsGlobal, &
    1, &
    MPI_INTEGER, &
    0, &    ! sending thread
    MPI_COMM_WORLD, &
    ierr)

if (ierr.ne.0) then
    write(*,*) "Could not broadcast number of points, error: ", ierr
endif

! determine how many elements will be handled by each of the threads

nBase = int(floor (dble(nElementsGlobal) / dble(nProcessors)))

if (myProcessorID.eq.0) then
    ! thread 0 gets the rest.
    nElementsLocal = nElementsGlobal - (nProcessors-1) * nBase
else
    ! all other threads get 
    nElementsLocal = nBase
endif

nElementsLocal0 = nElementsGlobal - (nProcessors-1) * nBase

! determine the offset for each element

allocate(elementOffset(nProcessors))
elementOffset(1) = 0
elementOffset(2) = nElementsGlobal - (nProcessors-1) * nBase
do iProcessor = 3,nProcessors
    elementOffset(iProcessor) = elementOffset(iProcessor-1) + nBase
enddo

! in a similar way, determine the offset for the vertices

nBase = int(floor (dble(nPointsGlobal) / dble(nProcessors)))

if (myProcessorID.eq.0) then
    ! thread 0 gets the rest.
    nPointsLocal = nPointsGlobal - (nProcessors-1) * nBase
else
    ! all other threads get 
    nPointsLocal = nBase
endif

! get the number of points of the zero partition,
! later, the, when the points are passed around to match
! the elements, even the block of the zero partition must
! fit every where, which would not be possible if it
! the array would be too small.
nPointsLocal0 = nPointsGlobal - (nProcessors-1) * nBase

!write(*,*) "thread", myProcessorID, "has nPointsLocal", nPointsLocal

allocate(vertexOffset(nProcessors))
vertexOffset(1) = 0
vertexOffset(2) = nPointsGlobal - (nProcessors-1) * nBase
do iProcessor = 3,nProcessors
    vertexOffset(iProcessor) = vertexOffset(iProcessor-1) + nBase
enddo

!*****************************************************************************************
!   read the element file, and only pick up the elements for this
!   processor.
!*****************************************************************************************

allocate(connectivityA(nPointsPerElement,nElementsLocal0))
!allocate(connectivityB(nPointsPerElement,nElementsLocal))

allocate(elementMarkers(nElementsLocal))

connectivityA = 0

elementCounterForThisPartition = 0

do iElement = 1, nElementsGlobal


    if (iElement .gt. elementOffset(myProcessorID + 1) .and. &  ! +1 because processor count starts at 0, but array index at 1
        elementCounterForThisPartition .lt. nElementsLocal) then

        elementCounterForThisPartition = elementCounterForThisPartition + 1

        ! read 3 points for 2 dimensions, and all 4 points for 3 dimenions: nDimensions+1
        read(41,*) idx, elementMarkers(elementCounterForThisPartition), &
        connectivityA(1:nDimensions+1, elementCounterForThisPartition)

    else
        ! the elements does not belong to this thread.
        ! skip it.
        read(41,*) idx
    endif

enddo

close(41)

!*****************************************************************************************
!   read the nodal point file, and only pick up the points for this
!   processor.
!*****************************************************************************************

allocate(pointCoordsA(nDimensions,nPointsLocal0))
allocate(pointCoordsB(nDimensions,nPointsLocal0))

allocate(pointMarkers(nPointsLocal))

pointCoordsA = 0d0
pointCoordsB = 0d0

pointMarkers = 0

pointCounterForThisPartition = 0

do iPoint = 1, nPointsGlobal

    if (iPoint .gt. vertexOffset(myProcessorID + 1) .and. &  ! +1 because processor count starts at 0, but array index at 1
       pointCounterForThisPartition .lt. nPointsLocal) then
        ! whoei, this point is for me!    

        pointCounterForThisPartition = pointCounterForThisPartition + 1

         read(40,*) idx, pointMarkers(pointCounterForThisPartition), &
            pointCoordsA(1:nDimensions, pointCounterForThisPartition)
    else
        ! vertex for another thread. Skip this one
         read(40,*) idx
    endif

enddo

close(40)


end subroutine
