module moduleMPI

implicit none

#include "mpif.h"

integer :: ierr
integer :: nProcessors, myProcessorID

contains

subroutine startMPI()

call MPI_Init(ierr)

if (ierr.ne.0) then
    write(*,*) "Could not initialize MPI, error: ", ierr
    stop "leaving SHRIMP..."
endif

call MPI_Comm_size(MPI_COMM_WORLD, nProcessors,   ierr)

if (ierr.ne.0) then
    write(*,*) "Could not determine number of processors, error: ", ierr
    stop "leaving SHRIMP..."
endif

call MPI_Comm_rank(MPI_COMM_WORLD, myProcessorID, ierr)

if (ierr.ne.0) then
    write(*,*) "Could not determine processor ID, error: ", ierr
    stop "leaving SHRIMP..."
endif

end subroutine

!-----------------------------------------------------------------------------

subroutine stopMPI()

implicit none

call MPI_Finalize(ierr)

if (ierr.ne.0) then
    write(*,*) "Could not close MPI, error: ", ierr
    stop "leaving SHRIMP..."
endif


end subroutine

!-----------------------------------------------------------------------------

subroutine parallelLog(message)

implicit none

character(len=*) :: message
integer          :: iThread, iError

call mpi_barrier(MPI_COMM_WORLD,ierror)

call sleep(1)

if (myProcessorID.eq.0) then
    write(*,*) '-----------------------------------------------------------------------'
endif


do iThread = 1,nProcessors
    call sleep(1)
    call mpi_barrier(MPI_COMM_WORLD,ierror)
    if (myProcessorID .eq. iThread-1) then
        write(*,*) "thread", myProcessorID, "says: ", message
    endif
enddo


call sleep(1)
if (myProcessorID.eq.0) then
    write(*,*) '-----------------------------------------------------------------------'
endif

call mpi_barrier(MPI_COMM_WORLD,ierror)


end subroutine

!-----------------------------------------------------------------------------


end module moduleMPI

