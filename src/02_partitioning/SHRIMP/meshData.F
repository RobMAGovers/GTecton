! this model contains the data that is actually used 
! to store the unpartitioned mesh and to partition it.


module meshData

implicit none


! variables that describe the topology of the mesh

integer            :: neighborBase
integer            :: nPointsPerElement
integer            :: nDimensions


! variables determined by reading of the input

integer :: nPointsGlobal, nElementsGlobal ! in the entire mesh
integer :: nPointsLocal,  nElementsLocal  ! stored in this processor
integer :: nPointsLocal0, nElementsLocal0 ! value of the 0 thread, known everywhere,
                                          ! so that when we swing the data around the processors,
                                          ! we do not take into account the zero-thread has a different
                                          ! amount of points/elements. We waste only a few bytes, but gain
                                          ! simpler code in exchange for it.


integer, allocatable :: pointMarkers(:)
integer, allocatable :: elementMarkers(:)

integer, allocatable :: elementOffset(:)
integer, allocatable ::  vertexOffset(:)

! variables that contain the rotating data

double precision, pointer             :: pointCoords(:,:)
double precision, target, allocatable :: pointCoordsA(:,:)
double precision, target, allocatable :: pointCoordsB(:,:) ! a mirror, for passing the
!                                                   ! partitions of the points around,
!                                                   ! without overwriting the originals.

integer, pointer             :: connectivity(:,:)
integer, target, allocatable :: connectivityA(:,:)
integer, target, allocatable :: connectivityB(:,:)  ! a mirror, for passing the
                                            ! partitions of the points around,
                                            ! without overwriting the originals.

integer, pointer             :: partitionOfPoints(:)
integer, target, allocatable :: partitionOfPointsA(:)
integer, target, allocatable :: partitionOfPointsB(:) ! a mirror, for passing the
                                              ! partitions of the points around,
                                              ! without overwriting the originals.


! variables that contain the partitioned mesh

integer      :: nbase  ! the base number of nodes per partition

!partitionOfPoints
!partitionNumber

! extra data to administrate parallellism

integer :: minNodeThisPartition
integer :: maxNodeThisPartition

! variables for the neighbor registration

integer, allocatable :: neighborCount(:) 

type intAllocArray
    integer, allocatable :: neighbors(:)
end type

type(intAllocArray), allocatable :: neighborlist(:)


integer, allocatable :: partitionOfElements(:)
integer, allocatable :: partitionOfElementPoints(:,:) ! gather the partition for each node per element here


integer, allocatable :: PointCountPerPartitionThread(:)
integer, allocatable ::	ElementCountPerPartitionThread(:)
integer, allocatable :: PointCountPerPartitionGlobal(:)
integer, allocatable ::	ElementCountPerPartitionGlobal(:)


! needed for output writing.

double precision, allocatable :: coordsForMyElems(:,:,:) ! needed for GMT and VTK compatible output


interface passData
    module procedure passDoubles, passIntegers1D, passIntegers2D
end interface

interface onePass
    module procedure onePassDoubles, onePassIntegers1D, onePassIntegers2D
end interface


contains

! the data is cycled round so that all processors 
! can get a look at the complete mesh, 
! without the any thread needing the entire mesh.
subroutine passDoubles(arraySize, A, B, iPass)

use moduleMPI, only: MPI_COMM_WORLD, &
                     MPI_STATUS_IGNORE, &
                     MPI_DOUBLE, &
                     nProcessors, &
                     myProcessorID

implicit none

integer          :: arraySize, iPass
double precision :: A(:,:), B(:,:)

integer          :: receiveFromThread, sendToThread
integer          :: iErr

! this one sends to me
receiveFromThread  = modulo(nProcessors + myProcessorID - 1, nProcessors)
! I send to this one
sendToThread       = modulo(myProcessorID + 1, nProcessors)

if (modulo(iPass,2) .eq. 1) then

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call onePass(arraySize, A, B, modulo(myProcessorID,2) .eq. 0, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call onePass(arraySize, A, B, modulo(myProcessorID,2) .eq. 1, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

else

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call onePass(arraySize, B, A, modulo(myProcessorID,2) .eq. 0, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call onePass(arraySize, B, A, modulo(myProcessorID,2) .eq. 1, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

endif

end subroutine

subroutine passIntegers1D(arraySize, A, B, iPass)

use moduleMPI, only: MPI_COMM_WORLD, &
                     MPI_STATUS_IGNORE, &
                     MPI_INT, &
                     nProcessors, &
                     myProcessorID

implicit none

integer          :: arraySize, iPass
integer          :: A(:), B(:)

integer          :: receiveFromThread, sendToThread
integer          :: iErr

! this one sends to me
receiveFromThread  = modulo(nProcessors + myProcessorID - 1, nProcessors)

! I send to this one
sendToThread       = modulo(myProcessorID + 1, nProcessors)

if (modulo(iPass,2) .eq. 1) then

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call onePass(arraySize, A, B, modulo(myProcessorID,2) .eq. 0, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call onePass(arraySize, A, B, modulo(myProcessorID,2) .eq. 1, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

else

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call onePass(arraySize, B, A, modulo(myProcessorID,2) .eq. 0, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call onePass(arraySize, B, A, modulo(myProcessorID,2) .eq. 1, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

endif

end subroutine

subroutine passIntegers2D(arraySize, A, B, iPass)

use moduleMPI, only: MPI_COMM_WORLD, &
                     MPI_STATUS_IGNORE, &
                     MPI_INT, &
                     nProcessors, &
                     myProcessorID

implicit none

integer          :: arraySize, iPass
integer          :: A(:,:), B(:,:)

integer          :: receiveFromThread, sendToThread
integer          :: iErr

! this one sends to me
receiveFromThread  = modulo(nProcessors + myProcessorID - 1, nProcessors)

! I send to this one
sendToThread       = modulo(myProcessorID + 1, nProcessors)

if (modulo(iPass,2) .eq. 1) then

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call onePass(arraySize, A, B, modulo(myProcessorID,2) .eq. 0, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call onePass(arraySize, A, B, modulo(myProcessorID,2) .eq. 1, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

else

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call onePass(arraySize, B, A, modulo(myProcessorID,2) .eq. 0, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
	call onePass(arraySize, B, A, modulo(myProcessorID,2) .eq. 1, sendToThread, receiveFromThread)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

endif

end subroutine

subroutine onePassDoubles(arraySize, sendArray, receiveArray, sendCondition, &
                          sendToThread, receiveFromThread)

use moduleMPI, only: MPI_COMM_WORLD, &
                     MPI_STATUS_IGNORE, &
                     MPI_DOUBLE, &
                     nProcessors, &
                     myProcessorID
implicit none

integer :: arraySize
double precision :: sendArray(:,:), receiveArray(:,:)
logical :: sendCondition
integer :: sendToThread, receiveFromThread

integer :: ierr

if (sendCondition) then

    call MPI_Send(sendArray,&
                  arraySize,&
                  MPI_DOUBLE,&
                  sendToThread,&
                  0,&
                  MPI_COMM_WORLD,&
                  ierr)
else

    call MPI_Recv(receiveArray, &
                  arraySize, &
                  MPI_DOUBLE, &
                  receiveFromThread, &
                  0, &
                  MPI_COMM_WORLD, &
                  MPI_STATUS_IGNORE, &
                  ierr)
endif

end subroutine

subroutine onePassIntegers1D(arraySize, sendArray, receiveArray, sendCondition, &
                             sendToThread, receiveFromThread)

use moduleMPI, only: MPI_COMM_WORLD, &
                     MPI_STATUS_IGNORE, &
                     MPI_INT, &
                     nProcessors, &
                     myProcessorID
implicit none

integer :: arraySize
integer :: sendArray(:), receiveArray(:)
logical :: sendCondition
integer :: sendToThread, receiveFromThread

integer	:: ierr

if (sendCondition) then

    call MPI_Send(sendArray,&
                  arraySize,&
                  MPI_INT,&
                  sendToThread,&
                  0,&
                  MPI_COMM_WORLD,&
                  ierr)
else

    call MPI_Recv(receiveArray, &
                  arraySize, &
                  MPI_INT, &
                  receiveFromThread, &
                  0, &
                  MPI_COMM_WORLD, &
                  MPI_STATUS_IGNORE, &
                  ierr)
endif

end subroutine

subroutine onePassIntegers2D(arraySize, sendArray, receiveArray, sendCondition, &
                             sendToThread, receiveFromThread)

use moduleMPI, only: MPI_COMM_WORLD, &
                     MPI_STATUS_IGNORE, &
                     MPI_INT, &
                     nProcessors, &
                     myProcessorID
implicit none

integer :: arraySize
integer :: sendArray(:,:), receiveArray(:,:)
logical :: sendCondition
integer :: sendToThread, receiveFromThread

integer	:: ierr

if (sendCondition) then

    call MPI_Send(sendArray,&
                  arraySize,&
                  MPI_INT,&
                  sendToThread,&
                  0,&
                  MPI_COMM_WORLD,&
                  ierr)
else

    call MPI_Recv(receiveArray, &
                  arraySize, &
                  MPI_INT, &
                  receiveFromThread, &
                  0, &
                  MPI_COMM_WORLD, &
                  MPI_STATUS_IGNORE, &
                  ierr)
endif

end subroutine


end module
