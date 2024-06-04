module debugmodule

use, intrinsic :: iso_c_binding, only: c_int  ! for sleep

use petscksp

! this module contains the debug vars and serves initially 
! as a replacement for the dbug common block.
! In a later stadium, the routines that print the debug info
! should also be moved to here.

integer, parameter :: byterPerDouble  = 8
integer, parameter :: byterPerInteger = 4


! itest = 0 no additional tests (default)
! itest = 1 perform additional tests on input data
! itest = 2 take a single timestep, output and exit

logical :: debug
integer :: itest, lustat
integer :: iecho
logical :: delayInParallelLog ! set to true when iecho 88. This is set to true, and iecho is set to 8, then.
logical :: usedPetscOptionsTxt

double precision :: petsctime
double precision :: owntime

integer     :: cr, cm, rate
integer     :: time_start, time_finish, clock_rate
integer     :: time_allStart
integer     :: iteration
integer     :: NicolaiNodeID, NicolaiFileID
public      :: Verified_MPI_Barrier
public      :: initialise_clock
public      :: memory_usage
public      :: initialiseDebug
public      :: ParallelLog
public      :: delay, longdelay
public      :: startClock
public      :: LogTimeStamp
public      :: xit
public      :: isFileOpen
logical, parameter   :: NewSample = .true.
integer, allocatable :: matrixEntriesPerRow(:)
double precision :: loadvec(20)
integer          :: debugFileID
integer :: iAmAtTimestep

#ifdef PLOTNONZEROES
logical, allocatable :: matrixEntryNonZero(:,:)
#endif

interface

      !  should be unsigned int ... not available in Fortran
      !  OK until highest bit gets set.
     function FortSleep (seconds)  bind ( C, name="sleep" )
          import
          integer (c_int) :: FortSleep
          integer (c_int), intent (in), VALUE :: seconds
     end function FortSleep

end interface


contains

!----------------------------------

#ifdef TEST

subroutine testIt()
! this sub has space for quick tests. Unit tests, etc.

use algebra, only: areaOfTriangleIn3D

implicit none

double precision :: xa, ya, za, xb, yb, zb, xc, yc, zc
double precision :: area

write(*,*) "Testing mode. As a user, this should not happen. Compile without -DTEST"


xa=0d0
ya=1d0
za=0d0

xb=1d0
yb=0d0
zb=0d0

xc=0.5d0
yc=0.5d0
zc=1.0d0


write(*,*) "calling area"

area = areaOfTriangleIn3D(xa, ya, za, xb, yb, zb, xc, yc, zc)

write(*,*)  "area: ", area


! room for test calls here.



end subroutine testIt

#endif

!----------------------------------

subroutine isFileOpen(fileHandle, tag)

implicit none

character(len=*) :: tag
integer          :: fileHandle
logical          :: isOpen         


inquire(unit=fileHandle, opened=isOpen)

if (isOpen) then
    write(*,*) "file ", fileHandle, "is open at tag", trim(tag)
else
    write(*,*) "file ", fileHandle, "is not open at tag", trim(tag)
endif

end subroutine

!----------------------------------

subroutine CheckArrayForNans(vector, name)

implicit none

character(len=*) :: name
double precision :: vector(:)
logical :: vectorHasNans
integer :: i


vectorHasNans = .false.
do i=1,size(vector,1)
    if (vector(i).ne.vector(i)) then
        vectorHasNans = .true.
        write(*,*) "vector ",trim(name)," has entry ", i, "NaN."
    endif
enddo

if (vectorHasNans) then
    write(*,*) "This is not supposed to happen. Contact model support"
    stop "Aborting execution"
else
    write(*,*) "vector ",trim(name), "is OK"
endif

end subroutine

!----------------------------------

!subroutine setStatus()

!implicit none

!call openf (lustat,statusFileName,'old')

!end subroutine
!----------------------------------

subroutine Verified_MPI_Barrier(message)
! This is an MPI barrier with every thread reporting its arrival
! So that it is possible to monitor if ertain threats never reach the barrier,
! in case of 'the program keeps running, but nothing happens'

#ifdef SPARSE
use modelctx, only: getrank
#endif

implicit none
character(len=*) :: message

#ifdef SPARSE

integer          :: ierr

    write(*,*) "Thread", getrank(), "reached barrier; ", trim(message)
    call mpi_barrier(0,ierr)
    if (ierr.ne.0) then
        write(*,*) "Thread", getrank(), "failed Verified_MPI_Barrier with error", ierr
    endif
    call mpi_barrier(0,ierr)

#endif
end subroutine


!----------------------------------

subroutine LogTimeStamp(caller)
! This sub prints the wall time used between this call and its previous call.
! It is assumed that the function
! InitializeClock ahs been called

#ifdef SPARSE
use modelctx, only: getrank, &
                    getsize
#endif

implicit none

character(len=*) :: caller

call system_clock(time_finish)

#ifdef SPARSE
write(*,*) caller, ' on rank', getrank(), &
          ' took: ', real(time_finish-time_start)/real(clock_rate), &
          ' seconds'
#else
write(*,*) caller, &
          ' took: ', real(time_finish-time_start)/real(clock_rate), &
          ' seconds'
#endif

time_start = time_finish




end subroutine

!----------------------------------

subroutine startClock()
! used to start the beginning of the time measuring block.
implicit none

call system_clock(time_start)

end subroutine

!----------------------------------

subroutine finalTime()

#ifdef SPARSE
use modelctx, only: getrank, &
                getsize 
#endif

implicit none

call system_clock(time_finish)

#ifdef SPARSE
 write(*,*) ' rank', getrank(), &
       ' took total time: ', real(time_finish-time_allStart)/real(clock_rate), &
       ' seconds'
#else
 write(*,*) ' took total time: ', real(time_finish-time_allStart)/real(clock_rate), &
            ' seconds'
#endif

end subroutine

!----------------------------------

subroutine AllocateError(thingy, errorID)
! message to announce and allocation went awry

    implicit none

    integer          :: errorID
    character(len=*) :: thingy

    if (errorID.ne.0) then
        write(*,*) 'ERROR could not allocate ', trim(thingy)
        write(*,*) 'ERROR', errorID, ' Exiting...'

        call xit(1," ")
    endif

end subroutine


!----------------------------------

subroutine ParallelLog(caller, message)
! this sub writes log messages in order of the threads
! To ensure the sequence is obeyed in the log files, there is a small
! delay caused by mpi_barriers, until all the threas past this point
! Using this logging could be costly (but als very useful) 
! with very many threads or very many times.

! mpi_barrier is not enough to accomplish this, because 
! the threads still buffer their standard out a bit,
! causing the sequence to be inconsistent.

! More over, to check whether threads have died or not, it is essential that
! this function can run even when some threads no longer call it,
! which would be prevented by using MPI barriers.

#ifdef SPARSE
use modelctx, only: getrank, &
                    getsize
#endif

implicit none


integer          :: thisrank, i
#ifdef SPARSE
integer          :: ierror
#endif
character(len=*) :: caller
character(len=*) :: message

integer,dimension(8) :: values

#ifdef SPARSE
thisrank = getrank()
#else
thisrank = 0
#endif

!call delay()
#ifdef SPARSE
!call delay()
    call mpi_barrier(0,ierror)
#endif

if (thisrank.eq.0) then
    write(*,*) '-----------------------------------------------------------------------'
endif 

!if (delayInParallelLog) then
!    call delay()
!endif

#ifdef SPARSE
    call mpi_barrier(0,ierror)
!    call delay()

#endif


#ifdef SPARSE
do i=0,getsize()
#else
i = 0
#endif


!    call delay()
#ifdef SPARSE
    call mpi_barrier(0,ierror)
#endif

    if (i.eq.thisrank) then
        call date_and_time(VALUES=values)

        write(*,"(a4,' ',i2,':',i2,':',i2,a8,' ',a,' ',a6,i4,a7,' ',a)") "time" ,values(5:7),' caller ', &
                   trim(caller), ' rank ', thisrank, ' says: ', trim(message)
    endif

!    call delay()

#ifdef SPARSE
enddo
#endif

#ifdef SPARSE
call mpi_barrier(0,ierror)
#endif

if (thisrank.eq.0) then
    write(*,*) '-----------------------------------------------------------------------'
endif

!call delay()
#ifdef SPARSE
call mpi_barrier(0,ierror)
#endif


end subroutine

!---------------------------------------------------------------

subroutine delay()
! Used for the delay of the logging.
! Time in nanoseconds.
! The sleep function of Fortran takes only integer seconds.
! which is often unnecessarily long and inconvenient
implicit none

integer :: ierr

!call smallSleep(1e8, ierr)

ierr = FortSleep ( 1 )


end subroutine
!---------------------------------------------------------------


subroutine longdelay()
! used for the delay of the logging
implicit none

integer :: how_long

how_long = FortSleep ( 5 )


end subroutine
!---------------------------------------------------------------

subroutine initialiseDebug()

#ifdef SPARSE

use modelctx, only: getrank, getsize

implicit none

integer, external :: nextlu
integer           :: thisrank, nRanks
character(len=1024) :: debugFileName

call initialise_clock()


if (debug) then
    thisRank = getrank()
    nRanks = getsize()

    if (nRanks.ge.1000) then ! just in case, although this many log files will be unwise.

        if (thisrank.ge.1000) then
            write (debugFileName, "(A20,I4,A4)") "GTectonLogFileThread", getrank(), ".log"
        else if (thisrank.ge.100) then
            write (debugFileName, "(A21,I3,A4)") "GTectonLogFileThread0", getrank(), ".log"
        else if (thisrank.ge.10) then
            write (debugFileName, "(A22,I2,A4)") "GTectonLogFileThread00", getrank(), ".log"
        else
            write (debugFileName, "(A23,I1,A4)") "GTectonLogFileThread000", getrank(), ".log"
        endif

    else if (nRanks.ge.100) then

        if (thisrank.ge.100) then
            write (debugFileName, "(A20,I3,A4)") "GTectonLogFileThread", getrank(), ".log"
        else if (thisrank.ge.10) then
            write (debugFileName, "(A21,I2,A4)") "GTectonLogFileThread0", getrank(), ".log"
        else
            write (debugFileName, "(A22,I1,A4)") "GTectonLogFileThread00", getrank(), ".log"
        endif

    else if (nRanks.ge.10) then

        if (thisrank.ge.10) then
            write (debugFileName, "(A20,I2,A4)") "GTectonLogFileThread", getrank(), ".log"
        else
            write (debugFileName, "(A21,I1,A4)") "GTectonLogFileThread0", getrank(), ".log"
        endif

    else

        write (debugFileName, "(A20,I1,A4)") "GTectonLogFileThread", getrank(), ".log"

    endif

    debugFileID = nextlu(10)
    open(unit=debugFileID, file=trim(debugFileName))
    write (debugFileID,*) "Logging GTecton events of thread ", getrank()
endif

#endif /* SPARSE */

end subroutine

!---------------------------------------------------------------

subroutine initialise_clock()

implicit none

! all variables directly from module

time_start = 0.0
time_finish = 0.0

call system_clock(time_allStart)
call system_clock(count_rate=clock_rate)
call system_clock(count_max=cm)

rate = REAL(cr)

end subroutine


!-------------------------------------------------------------------------

subroutine memory_usage()
 ! printd the process using the most memory, and the quantity, and the process.
 ! on a compute node, this is essentially always gtecton.
 !    call csystem('ps -e -orss=,args= | sort -b -k1,1n | tail -n 1')
end subroutine

!-------------------------------------------------------------------------

subroutine xit (errorCode, errorMessage)


USE MODELCTX
use filemodule,  only: do_soelim, &
                       statusFileName
#ifdef SPARSE
use sPETScModule, only: userFinalize
#endif
use    iomodule,  only: stderr, luin


implicit none


logical there
integer errorCode
character(len=*) :: errorMessage

inquire (file=statusFileName,exist=there)

if (iecho.eq.8) then
    call ParallelLog("xit", "Entered sub")
endif


if (there .and. lustat.gt.0) then
!    call openf (unit=lustat,file=statusFileName,staus="old",position="append",action="write")
    if (errorCode.eq.0) then
#ifdef SPARSE
        if (getrank().eq.0) then
#endif
            open (unit=lustat,file=statusFileName,status="old",position="append",action="write")
            write(lustat,10)
10              format('Done')
#ifdef SPARSE
        endif
#endif
    else
        call openf (lustat,statusFileName,'unknown')
        write(lustat,20) errorCode
20              format('Exit ',I2)
    endif
    call closef (lustat)
endif


if (errorCode.ne.0) then
#ifdef SPARSE
    write(stderr,30) getrank()
30          format('Leaving application from rank ',I5)
#else
    write(stderr,*) 'Leaving application'
#endif
!    most xit calls do not yet give an error message
     if (errorMessage(1:1).ne." ") then
         write(*,*) "Error message: ", trim(errorMessage)
     endif
endif

! in case gtecton has soelimt in the /tmp directory,
! delete this file
if (do_soelim) then
    close(unit=luin,status='delete')
else
    close(unit=luin)
endif 

if (debug) then
    close(debugFileID)
endif

#ifdef SPARSE
if (iecho.eq.8) then
    call ParallelLog("xit", "Calling userfinalize")
endif

call UserFinalize()
! userfinalize also closes mpi, so can no longer do parallellog
!if (iecho.eq.8) then
!    call ParallelLog("xit", "Finished userfinalize")
!endif
#endif

call finalTime()

call exitp(errorCode)  ! from tools

return
end subroutine

!-------------------------------------------------------------------------

! this sub should be move to a module that is not dependent on debugmodule,
! to prevent cyclic dependency with AO-module

!#ifdef PLOTNONZEROES

!subroutine nonZeroMatrixPlot()

!#ifdef SPARSE
!use aomodule, only: nequations
!use modelctx, only: getrank
!#endif

!implicit none

!character(len=50) :: imageFile, rankText
!integer :: imageID
!integer, external :: nextlu
!integer :: iMatEntry, jMatEntry
!integer :: matrixSize
!integer :: distanceFromDiagonal


!#ifdef SPARSE
!write(rankText,'(a1)') getrank()
!imageFile="nonZeroStructure" // rankText // ".bpm"
!#else
!imageFile="nonZeroStructure.bpm"
!#endif

!#ifdef SPARSE
!! AOmodule exists. Get nequations from there.
!matrixSize = nequations
!#else
!matrixSize = 5
!#endif

!distanceFromDiagonal = 0

!imageID = nextlu(12)

!open(unit=imageID, file=imageFile)

!write(imageID,1) "P1"
!1 format (a2)
!write(imageID,2) "# This is the non zero structure of the solution matrix"
!2 format (a55)

!write(imageID,*) matrixSize,matrixSize

!do iMatEntry = 1,matrixSize
!    do jMatEntry = 1,matrixSize-1
!        if (matrixEntryNonZero(iMatEntry, jMatEntry)) then
!            write(imageID,10,advance="no") 1
!            distanceFromDiagonal = distanceFromDiagonal + abs(iMatEntry-jMatEntry)
!        else
!            write(imageID,10,advance="no") 0
!        endif
!    enddo
!    if (matrixEntryNonZero(iMatEntry, matrixSize)) then
!        write(imageID,20) 1
!        distanceFromDiagonal = distanceFromDiagonal + abs(iMatEntry-jMatEntry)
!    else
!        write(imageID,20) 0
!    endif
!enddo

!10 format (i1,1x)
!20 format (i1)

!close(imageID)

!write(*,*) "distanceFromDiagonal: ", distanceFromDiagonal

!end subroutine

!#endif

!-------------------------------------------------------------------------

end module
