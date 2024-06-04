! This is SHRIMP, the Scalable Highly Recursive Irregular Mesh Partitioner

! It can partition any mesh, in 2D or 3D, without having to resort to graph theory.

! It workflow consist of the following parts

! 1) Read input files of the unpartitioned mesh
! 2) Determine partitioning of the mesh points, 
! 3) Determine element partitioning
! 4) when necessary, compute neighbors of the points,
!    which can be used for finding ghost elements/points
! 5) Write the output in a desired format


program shrimp

use enumerates
use meshData
use moduleMPI, only: startMPI, &
                     stopMPI, &
                     nProcessors, &
                     myProcessorID, &
                     parallelLog



implicit none

#include "mpif.h"


! input parameters

integer            :: nSteps

integer, parameter :: neighborBase2D = 10 ! will typically be less, for base allocation
integer, parameter :: neighborBase3D = 20

integer, parameter :: firstPartitionNumber = 1  ! we start to count partition numbers at 1

integer            :: inputFormat, outputFormat
!integer, parameter :: outputformat        = outputformatSimple ! for formats, see enumerates module
!integer, parameter :: outputformat        = outputformatVTK ! for formats, see enumerates module
!integer, parameter :: outputformat        = outputformatGTECTON ! for formats, see enumerates module


! indices and such

integer :: iError


! initiate pointers to the flip arrays.
pointCoords => pointCoordsA
connectivity =>    connectivityA
partitionOfPoints => partitionOfPointsA

! this will eventually set nDimensions and filenames and such
call parseCommandLineArgs(inputFormat, outputFormat, nSteps)

call startMPI()

if (nDimensions.eq.2) then
    nPointsPerElement = 3
    neighborBase = neighborBase2D
else
    nPointsPerElement = 4
    neighborBase = neighborBase3D
endif

allocate(PointCountPerPartitionThread(2**nSteps))
allocate(ElementCountPerPartitionThread(2**nSteps))
PointCountPerPartitionThread = 0
ElementCountPerPartitionThread = 0


call readInput(inputFormat)

allocate(partitionOfPointsA(nPointsLocal0))
partitionOfPointsA = 0

call splitPointCloud(nSteps, &          ! calling depth, decreasing for deeper recursion 
                     nPointsLocal, &    ! number of points on this call on this thread
                     nPointsGlobal, &   ! number of points on this call on all threads
                     pointCoordsA, &    ! coordinates of local points
                     partitionOfPointsA, &
                     firstPartitionNumber)


call matchElements(outputFormat)

call findNeighbors()

call writeOutput(outputFormat, nSteps)

call MPI_Barrier(MPI_COMM_WORLD, iError)
call stopMPI()

end program

!---------------------------------------------------------------------------------

subroutine parseCommandLineArgs(inputFormat, outputFormat, nSteps)

use files
use enumerates
use meshdata,  only: nDimensions
use moduleMPI, only: myProcessorID


implicit none

integer            :: nargs, iarg, nSteps
integer            :: inputFormat, outputFormat
character(len=1024) :: arg

logical            :: pset, eset, iset, oset, nset, dset

integer, external  :: stringToInteger

nargs  = command_argument_count()

pset = .false.
eset = .false.
iset = .false.
oset = .false.
nset = .false.
dset = .false.

iarg = 1

do while (iarg.le.nargs)

    call get_command_argument(iarg , arg)       ! tecin.dat.partf.nps

    ! strings (filenames)
    if (trim(arg).eq."-p") then                ! -p    point file name
        iarg = iarg + 1
        call get_command_argument(iarg , arg)
        pointFile = trim(arg)
        pset = .true.
    else if (trim(arg).eq."-e") then             ! -e    element file name
        iarg = iarg + 1
        call get_command_argument(iarg , arg)
        elementFile = trim(arg)
        eset = .true.
    else if (trim(arg).eq."-i") then             ! -i    input formant
        iarg = iarg + 1
        call get_command_argument(iarg , arg)
        if (trim(arg).eq.'gtecton' .or. &
            trim(arg).eq.'GTECTON') then
            inputFormat = inputFormatGTECTON
        else if (trim(arg).eq.'simple' .or. &
                 trim(arg).eq.'SIMPLE') then
               inputFormat = inputFormatSimple
           endif
        iset = .true.
    else if (trim(arg).eq."-o") then             ! -o    output format
        iarg = iarg + 1
        call get_command_argument(iarg , arg)
        if      (trim(arg).eq.'vtu' .or. &
                 trim(arg).eq.'VTU' .or. & 
                 trim(arg).eq.'vtk' .or. &
                 trim(arg).eq.'VTK') then
            outputFormat = outputFormatVTK
        else if (trim(arg).eq.'gtecton' .or. &
                 trim(arg).eq.'GTECTON') then
            outputFormat = outputFormatGTECTON
        else if (trim(arg).eq.'gmt' .or. &
                 trim(arg).eq.'GMT') then
            outputFormat = outputFormatGMT
        else if (trim(arg).eq.'simple' .or. &
                 trim(arg).eq.'SIMPLE') then
            outputFormat = outputFormatSimple
        endif
        oset = .true.

    ! integers
    else if (trim(arg).eq."-n") then             ! -n    split steps
        iarg = iarg + 1
        call get_command_argument(iarg , arg)
        nSteps = stringToInteger(trim(arg))
        nset = .true.
    else if (trim(arg).eq."-d") then             ! -d    nDimenions
        iarg = iarg + 1
        call get_command_argument(iarg , arg)
        nDimensions = stringToInteger(trim(arg))
        dset = .true.
    else
        if (myProcessorID.eq.0) then
!            write(*,*) "argument "\\trim(arg)\\" not recognized."
            call usage()
        endif        
    endif

    iarg = iarg + 1

enddo

! let us be pedantic with errors



if (myProcessorID.eq.0) then
    if (.not. pset) then
        write(*,*) "argument -p not set: point file"
    endif
    if (.not. pset) then
        write(*,*) "argument -e not set: element/connectivity file"
    endif
    if (.not. pset) then
        write(*,*) "argument -i not set: input format. Options: simple, gtecton"
    endif
    if (.not. pset) then
        write(*,*) "argument -o not set: output format. Options: simple, gtecton, gmt (only for 2D), vtu (only for 3D)"
    endif
    if (.not. pset) then
        write(*,*) "argument -n not set: nSteps"
    endif
    if (.not. pset) then
        write(*,*) "argument -d not set: nDimensions"
    endif


    if (.not. pset .or. &
        .not. eset .or. &
        .not. iset .or. &
        .not. oset .or. &
        .not. nset .or. &
        .not. dset) then
    
        write(*,*) "Not all arguments set"
        call usage()
        stop "Leaving SHRIMP"
    endif
endif



end subroutine

!---------------------------------------------------------------------------------

integer function stringToInteger(string)

implicit none

character(len=1) :: string

read(string(1:1),*) stringToInteger

end function

!---------------------------------------------------------------------------------

subroutine usage()

implicit none

write(*,*) " "
write(*,*) "usage: "
write(*,*) " "
write(*,*) "shrimp -p [pointfile] -e [elementfile] -i [inputFormat] -o [outputFormat] -n [nSteps] -d [nDimensions]"
write(*,*) " "
write(*,*) "available input formats:  simple, gtecton, currently only GTECTON input accepted"
write(*,*) "available output formats: simple, gtecton, gmt (only for 2D), vtu (only for 3D)"
write(*,*) "nSteps: domain will be split into 2^nSteps partitions."
write(*,*) "nDimensions: 2 or 3"
write(*,*) " "

end subroutine
