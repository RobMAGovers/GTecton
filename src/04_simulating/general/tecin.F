subroutine READFE(NSD,NDOF,NSTR,NEN,ierr)

! Reads partition file and TECIN.DAT file. All the .so files have already 
! been groffed into it, using soelim
! There are several routines (such as IGEN adn FGEN) that are also used
! by rdthrm.F, which reads the thermal input (typically TEMPIN.DAT)

USE LABELMODULE
USE MESHDATAMODULE,  only: meshdatactx, &
                           allocatevertexdata, &
                           allocateelementdata, &
                           BuildElementNeighbourTable2D
USE MODELDATAMODULE, only: modeldatactx, &
                           allocatemechmodeldata, &
                           allocatemechbcdata, &
                           allocatemechwinklerbcdata, &
                           deallocatemechwinklerbcdata, &
                           allocatemechmaxwellbcdata, &
                           IDIAG
USE MODELDEFINITION
USE TIMESTEPMODULE
USE MATERIALSMODULE, only: LMAT
USE MODELCTX
USE FILEMODULE       ! only, do_soelim.... maybe more, statusFileName. filter out implicitnoning
use formatsmodule,   only: FMT_tecin_dat_MeshCount, &
                           FMT_tecin_dat_MODEline
use debugmodule,     only: iecho, &
                           parallellog, &
                           Verified_MPI_Barrier, &
                           delay, &
                           xit, debug
use constants,       only: NNeighborsMax2D, &
                           NNeighborsMax3D, &
                           elementTypeP1P0, &
                           elementTypeP2P1
use algebra,         only: idig
use modeltopology,   only: elementType
use iomodule          

#ifdef SPARSE
USE AOMODULE
#endif

implicit none

!-pass
integer NSD,NDOF,NSTR,NEN,ierr

#ifdef SPARSE

#include "petsc/finclude/petscsys.h"

PetscErrorCode petscierr
PetscMPIInt tmp, irank
PetscInt size_in, rank_in, numnplocal_in, numellocal_in
#endif

!-locl
character(len=1024) FIN,tmpnm
character(len=1024) record

character quote
integer ERROR
integer, external :: soelim,ifopen,lnblk ! from tools
integer, external :: csystem
logical there

#ifdef EXPERIMENTAL_CYCLIC
integer :: useCyclic
#endif

! for a check to see whether the number of
! elements/vertices in TECIN matches the sum
! of all entries in partition.info
#ifdef SPARSE
integer :: nElementsSumPartitionInfo
integer :: nVerticesSumPartitionInfo

integer :: iProc
#endif
integer :: ios
integer :: i, j

integer :: thisRank, nRanks
integer :: maxNbrs

logical :: emptyPartitions

emptyPartitions = .false.

if (iecho.eq.8) then
    call parallelLog("readfe","Start reading TECIN file")
endif


!-init
#ifdef SPARSE
call openf (88,statusFileName,'unknown')
write(88,*) 'Reading input data'
close (88)
#endif

!allocate(A(4000000000)) ! 4 gigabytes

quote = CHAR(34)

thisRank = getrank()
nRanks = getsize()

! Let us assume we do not run gtecton on more than 10k cores. 
! If so, the wise user will soelim before the run, in stead
! of clogging the /tmp space with countless FEM files.

if (iecho.eq.8) then
    call parallelLog("readfe","Constructing FIN")
endif

if      (nRanks.ge.10000) then ! just in case, although this many log files will be unwise.
    write(stdout,*) "This will create more than 10k files in de /tmp directory."
    write(stdout,*) "Most likely this will not fit and the run will crash when "
    write(stdout,*) "reading the incomplete files."
    write(stdout,*) "Solution: Use soelim and make one input file to be read"
    write(stdout,*) "by all threads."
    write(stdout,*) "Example: Soelim -r TECIN.DAT > fulltecin.dat"
    write(stdout,*) "And run: pln/axi/f3d [other options] fein=fulltecin.dat noso"
    stop "Aborting run"
else if (nRanks.ge.1000) then ! just in case, although this many log files will be unwise.

    if (thisrank.ge.1000) then
        write (FIN, "(A3,I4)") "FEM", getrank()
    else if (thisrank.ge.100) then
        write (FIN, "(A4,I3)") "FEM0", getrank()
    else if (thisrank.ge.10) then
        write (FIN, "(A5,I2)") "FEM00", getrank()
    else
        write (FIN, "(A6,I1)") "FEM000", getrank()
    endif

else if (nRanks.ge.100) then

    if (thisrank.ge.100) then
        write (FIN, "(A3,I3)") "FEM", getrank()
    else if (thisrank.ge.10) then
        write (FIN, "(A4,I2)") "FEM0", getrank()
    else
        write (FIN, "(A5,I1)") "FEM00", getrank()
    endif

else if (nRanks.ge.10) then

    if (thisrank.ge.10) then
        write (FIN, "(A3,I2)") "FEM", getrank()
    else
        write (FIN, "(A4,I1)") "FEM0", getrank()
    endif

else

    write (FIN, "(A3,I1)") "FEM", getrank()

endif

if (iecho.eq.8) then
    call parallelLog("readfe","completed FIN")
endif


! When pln, axi, f3d are run, the workpath is set from the command line
! arguments.
! When plnplt etc are run, it is not. In that case, automatically set it to the
! current directory
!if (.not. WORKPATHIsSet) then
!    WORKPATH = "."
!    WORKPATHIsSet = .true.
!endif!


!if (iecho.eq.8) then
!    call parallelLog("readfe","reading path environment variables")
!endif

! getenv MYTMP and PWD. If MYTMP does not exist it is created
call readPathsFromEnvironment()

!if (iecho.eq.8) then
!    call parallelLog("readfe","finished reading path environment variables")
!endif

FEMdir = trim(myTmpDir)
FIN    = trim(FEMdir) // "/" //tmpnm(FIN)

#ifdef SPARSE

if (iecho.eq.8) then
    call parallelLog("readfe","determining MPI rank")
endif

call MPI_Comm_rank (PETSC_COMM_WORLD,irank,ierr)

if (iecho.eq.8) then
    call parallelLog("readfe","MPI rank determined")
endif


!************************* read partition.info ************************

inquire(file=PARTINFO, exist=there)

if (there) then

    if (iecho.eq.8) then
        call parallelLog("readfe","reading partition info")
    endif
    ierr = ifopen(luin,PARTINFO,'old')
    if (ierr.ne.0) then
        write (0,*) 'Failed to open partition file: ', PARTINFO 
        return
    else
        if (iecho.eq.2) then
            write(stdout,*) 'Opened partition file: ', PARTINFO
        endif
    endif
 

    read(luin,48) size_in
    if (size_in.ne.getsize()) then
        write (0,*) 'Size_in differs from getsize. FAIL!'

        write(stderr,52) size_in, PARTINFO, getsize()
        call xit(1, "Number of used processors does not match partitioning")
    endif


    if (iecho.eq.8) then
        call parallelLog("readfe","Allocating partition arrays")
    endif
    allocate(meshdatactx%NvlocalOfAllPartitions(getsize()))
    if (ierr.ne.0) then
        write(stdout,*) "Could not allocate eshdatactx%NvlocalOfAllPartitions"
        write(stdout,*) "Error: ", ierr
    endif


    allocate(meshdatactx%NelocalOfAllPartitions(getsize()), stat=ierr)
    if (ierr.ne.0) then
        write(stdout,*) "Could not allocate eshdatactx%NelocalOfAllPartitions"
        write(stdout,*) "Error: ", ierr
    endif
    if (iecho.eq.8) then
        call parallelLog("readfe","Succesfully allocated partition arrays")
    endif


    do iProc=1,getsize()

        read(luin,*) rank_in, numnplocal_in, numellocal_in

		if (numnplocal_in .eq. 0) then
			write(stderr,*) "Partition", iProc, "has 0 nodal points."
            emptyPartitions = .true.
		endif

        if (numellocal_in .eq. 0) then
            write(stderr,*) "Partition", iProc, "has 0 elements."
            emptyPartitions = .true.
        endif

        meshdatactx%NvlocalOfAllPartitions(iProc) = numnplocal_in
        meshdatactx%NelocalOfAllPartitions(iProc) = numellocal_in

        if (rank_in.eq.getrank()) then
            meshdatactx%Nvlocal = numnplocal_in
            meshdatactx%Nelocal = numellocal_in
        endif

    enddo

	if (emptyPartitions) then
         write(stderr,*) "The mesh has one or more empty partitions."
         write(stderr,*) "Possible cause: too many partitions for a coarse mesh"
         stop "Leaving GTECTON..."
	endif




    close (luin)

    ! from the nelocals, establish the first local number of the partition.
    ! This number is used later on to access the proper position
    ! in the DIFORC array

    meshdatactx%lowestLocalElemID = 0
    do iProc=1,getrank()
        meshdatactx%lowestLocalElemID = &
           meshdatactx%lowestLocalElemID + &
           meshdatactx%NelocalOfAllPartitions(iProc)
    enddo

    meshdatactx%lowestLocalVertexID = 0
    do iProc=1,getrank()
        meshdatactx%lowestLocalVertexID = &
           meshdatactx%lowestLocalVertexID + &
           meshdatactx%NvlocalOfAllPartitions(iProc)
    enddo

else
    write (0,*) 'no partition info file ', PARTINFO
    return
endif

if (iecho.eq.8) then
    call parallelLog("readfe","finished reading partition info")
endif

48   format(I10)
52   format(1x, I10, ' partitions in partition info file',1x,A, &
   1x,'whereas',i10,1x,'processors are used')

#else
    call setrank(0)
    call setsize(1)
#endif /* SPARSE */


if (do_soelim) then
    inquire(file=FIN,exist=there)
    if (.not.there) then

        if (iecho.eq.8) then
            call parallelLog("readfe","Running Soelim")
        endif
        ierr = soelim (FEIN,FIN)
        if (ierr.ne.0) then
            write(stderr,5) ierr
5           format(1x,'soelim: exit ',I1)
            return
        endif
        if (iecho.eq.8) then
            call parallelLog("readfe","Soelim completed")
        endif
    endif
else
    FIN = FEIN
endif


if (iecho.eq.2) then
    write(stdout,*) 'Reading input file: ', FIN
endif
ierr = ifopen(luin,FIN,'old')
if (ierr.ne.0) then
   return
endif
 

if (iecho.eq.2) then
    write(stdout,10)
10  format(1x,'*** MECHANICAL PROBLEM DEFINITION ***'/ &
        1x,'read TITLE (A1024)')
endif

if (iecho.eq.8) then
    call parallelLog("readfe","Actually reading the file")
endif


read(luin,'(a)',err=1000,end=1100,iostat=ios) TITLE
if (TITLE(1:7) .eq. "element") then
    ! we have deviating element type
    if (TITLE(8:11) .eq. "p1p0") then
        elementType = elementTypeP1P0
    else if (TITLE(8:11) .eq. "p2p1") then
        elementType = elementTypeP2P1
    else
        write(stdout,*) "requested a custom element type: ", TITLE(8:11)
        write(stdout,*) "element type not recognized"
    endif
    ! previous line was apparently not the title, but the element type. 
    ! Read the title from the second line.
    read(luin,'(a)',err=1000,end=1100,iostat=ios) TITLE
endif


if (iecho.eq.2) then
    write(stdout,*) 'rank ', getrank(), 'reads NUMNP, NUMEL, NUMAT, maxNbrs', FMT_tecin_dat_MeshCount
    write(stderr,15)
endif
15   format(1x,'read NUMNP',1x,'(I12)')


!############ read second line ############
maxNbrs = 0
read(luin,fmt=FMT_tecin_dat_MeshCount,err=1000,end=1100,iostat=ios) &
     NUMNP, NUMEL, NUMAT, maxNbrs

#ifdef SPARSE
! check whether NUMNP en NUMEL correspond to the amounts in partition.info.
! Because these files are generated, this will usually be the case,
! but manual messing around *frowns at certain PhD student who will remain unnamed* ;)
! can create differences and cause nasty errors in strange places.

! only for sparse, because plnplt does not read partition.info

if (getsize().eq.1) then
    if (NUMNP.ne.numnplocal_in) then
        write(stdout,*) "Number of nodal points in partition file does not match" // &
        " number of nodal points in TECIN file"
        write(stdout,*) trim(PARTINFO), " has", numnplocal_in, " nodal points" // &
        " while ", trim(FEIN), " has", NUMNP, "nodal points" 
        call xit(1,"NUMNP mismatch")
    endif
    if (NUMEL.ne.numellocal_in) then
        write(stdout,*) "Number of elements in partition file does not match" // &
        " number of elements in TECIN file"
        write(stdout,*) trim(PARTINFO), " has", numellocal_in, " elements" // &
        " while ", trim(FEIN), " has", NUMEL, "elements"
        call xit(1,"NUMEL mismatch")
    endif
else
    nVerticesSumPartitionInfo = 0
    nElementsSumPartitionInfo = 0    

    ! sum partition info entries over the processors.
    call MPI_barrier(MPI_COMM_WORLD, ierr)

    call MPI_Scan(meshdatactx%Nvlocal, &
                  nVerticesSumPartitionInfo, &
                  1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    call MPI_barrier(MPI_COMM_WORLD, ierr)

    call MPI_Scan(meshdatactx%Nelocal, &
                  nElementsSumPartitionInfo, &
                  1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    call MPI_barrier(MPI_COMM_WORLD, ierr)

    if (getrank()+1.eq.getsize()) then
        ! because when summed over the processors, only the highest processor has the total sum        
        ! +1 because rank starts counting at 0, but size at 1
        if (NUMNP.ne.nVerticesSumPartitionInfo) then
            write(stdout,*) "Number of nodal points in partition file does not match" // &
            " number of nodal points in TECIN file"
            write(stdout,*) trim(PARTINFO), " has", nVerticesSumPartitionInfo, " nodal points" // &
            " while ", trim(FEIN), " has", NUMNP, "nodal points"
            call xit(1,"NUMNP mismatch")
        endif
        if (NUMEL.ne.nElementsSumPartitionInfo) then
            write(stdout,*) "Number of elements in partition file does not match" // &
            " number of elements in TECIN file"
            write(stdout,*) trim(PARTINFO), " has", nElementsSumPartitionInfo, " elements" // &
            " while ", trim(FEIN), " has", NUMEL, "elements"
            call xit(1,"NUMEL mismatch")
        endif
    endif
endif

#endif

if (iecho.eq.2 .or. debug) then
#ifdef SPARSE
    write(stdout,*) 'rank ', getrank(), 'reads NUMNP, NUMEL, NUMAT', NUMNP, NUMEL, NUMAT
#else
    write(stdout,*) 'rank 0 reads NUMNP, NUMEL, NUMAT', NUMNP, NUMEL, NUMAT
#endif
endif

meshdatactx%Nvglobal = NUMNP
meshdatactx%Neglobal = NUMEL

#ifndef SPARSE
meshdatactx%nvlocal = NUMNP
meshdatactx%nelocal = NUMEL
#endif

if (maxNbrs.eq.0) then
        if (NDOF.eq.2) then
                meshdatactx%nMaxNeighbors = NNeighborsMax2D
        else
        meshdatactx%nMaxNeighbors = NNeighborsMax3D
        endif
else
        meshdatactx%nMaxNeighbors = maxNbrs
endif

! calling meshdatamodule.F
call allocatevertexdata(meshdatactx, meshdatactx%nMaxNeighbors, NSD) 

call allocateelementdata(meshdatactx, NEN)

!############ read third line ############

#ifdef EXPERIMENTAL_CYCLIC
if (iecho.eq.2) then
    write(stderr,25)
25   format(1x,'read MODE,NINTG,NMPRT,MAXIT,NLINK,ICVIS,ISTART,', &
       'NELAST,NPRMAT,FaultTVelo,UseCyclic (FMT_tecin_dat_MODEline)')
endif
#else
if (iecho.eq.2) then
    write(stderr,25)
25   format(1x,'read MODE,NINTG,NMPRT,MAXIT,NLINK,ICVIS,ISTART,', &
       'NELAST,NPRMAT,FaultVelo (FMT_tecin_dat_MODEline)')
endif
#endif


read(luin,'(a)',err=1000,end=1100,iostat=ios) record

if (record(5:5).eq." ") then
    write(stdout,*) "Input reader says: MODE flag incorrect. Position 5 empty"
endif
if (record(10:10).eq." ") then
    write(stdout,*) "Input reader says: Number of time step groups flag incorrect. Position 10 empty"
endif
if (record(15:15).eq." ") then
    write(stdout,*) "Input reader says: Number of output times flag incorrect. Position 15 empty"
endif
if (record(20:20).eq." ") then
    write(stdout,*) "Input reader says: Max It flag incorrect. Position 20 empty"
endif
if (record(25:25).eq." ") then
    write(stdout,*) "Input reader says: Number of linked nodes incorrect. Position 25 empty"
endif
if (record(30:30).eq." ") then
    write(stdout,*) "Input reader says: ICVIS flag incorrect. Position 30 empty"
endif
if (record(35:35).eq." ") then
    write(stdout,*) "Input reader says: ISTART flag incorrect. Position 35 empty"
endif
if (record(40:40).eq." ") then
    write(stdout,*) "Input reader says: NELAST flag incorrect. Position 40 empty"
endif
if (record(45:45).eq." ") then
    write(stdout,*) "Input reader says: NPRMAT flag incorrect. Position 45 empty"
endif
if (record(50:50).eq." ") then
    write(stdout,*) "Input reader says: FaultVelo flag incorrect. Position 50 empty"
endif
#ifdef EXPERIMENTAL_CYCLIC
if (record(55:55).eq." ") then
    write(stdout,*) "Input reader says: UseCyclicSlip flag incorrect. Position 55 empty"
endif
#endif

#ifdef EXPERIMENTAL_CYCLIC
read(record,fmt=FMT_tecin_dat_MODEline,err=1000,end=1100,iostat=ios) & 
     MODE,NINTG,NMPRT,MAXIT,NLINK,ICVIS,ISTART,NELAST,NPRMAT,FaultVelo,UseCyclic
#else
read(record,fmt=FMT_tecin_dat_MODEline,err=1000,end=1100,iostat=ios) &  
     MODE,NINTG,NMPRT,MAXIT,NLINK,ICVIS,ISTART,NELAST,NPRMAT,FaultVelo
#endif


if (iecho.eq.2 .or. debug) then
#ifdef EXPERIMENTAL_CYCLIC
    write(stdout,*) "read ", MODE,NINTG,NMPRT,MAXIT,NLINK,ICVIS,ISTART,NELAST,NPRMAT,FaultVelo,UseCyclic
#else
    write(stdout,*) "read ", MODE,NINTG,NMPRT,MAXIT,NLINK,ICVIS,ISTART,NELAST,NPRMAT,FaultVelo
#endif
endif

if (debug) then
    write(stdout,*) "read MODE: ",MODE
    write(stdout,*) "read number of time step groups: ",NINTG
    write(stdout,*) "read number of output times: ",NMPRT
    write(stdout,*) "read maximum number of iteration between stiffness matrix recalculation: ",MAXIT
    write(stdout,*) "read number of linked nodes: ",NLINK
    write(stdout,*) "read viscosity flags: ",ICVIS
    write(stdout,*) "read ISTART",ISTART
    write(stdout,*) "read NELAST",NELAST
    write(stdout,*) "read NPRMAT",NPRMAT
    write(stdout,*) "read output velocity? ",FaultVelo
#ifdef EXPERIMENTAL_CYCLIC
    write(stdout,*) "read UseCyclic: ", UseCyclic
#endif
endif

#ifdef EXPERIMENTAL_CYCLIC
if (UseCyclic.eq.1) then
    UseCyclicSlip = .true.
else
    UseCyclicSlip = .false.
endif
#endif

LMAT = (NINTG.gt.0 .and. ICVIS.ne.0 .and. NUMAT.ne.NUMEL)

if (IDIG(ICVIS,4,2).eq.0) then
    ICVIS = IDIG(ICVIS,4,1)*1000 + IDIG(ICVIS,4,3)*10
endif

! calling timestepmodule.F
call allocateoutputcontroldata(outputcontroldatactx,0)

if (iecho.eq.2) then 
   write(stderr,35)
35   format(1x,'*** read IMPRINT data ***'/ &
 1x,'nread = MAX(1,NMPRT)'/ &
 1x,'read (10000i5) (IMPRINT(i), i=1,nread)')
endif

!############ read output times (when to write fedsk files) ############
if (iecho.eq.8) then
    call ParallelLog("Readfe", "reading timesteps")
endif

call RDTIME(outputcontroldatactx%IMPRINT, NMPRT)

if (iecho.eq.8) then
    call ParallelLog("Readfe", "finished reading timesteps")
endif


if (OUTFIL(1)) then
    if (ISTART.lt.0) then
        call lueof(luout(1))
    endif
     call writst(luout(1),title)
    ! write FEOUT.DAT
    write(luout(1),40) NUMNP,NSD,NDOF,MODE,ISTART,NELAST
endif

40   FORMAT(1x/1x/1x/ &
   ' C O N T R O L   I N F O R M A T I O N                 '   //5X,&
   ' NUMBER OF NODAL POINTS . . . . . . . . . . . (NUMNP) =',I12//5X,&
   ' NUMBER OF SPACE DIMENSIONS . . . . . . . . . (NSD  ) =',I8//5X,&
   ' NUMBER OF DEGREES OF FREEDOM PER NODE  . . . (NDOF ) =',I8//5X,&
   ' SOLUTION CODE  . . . . . . . . . . . . . . . (MODE ) =',I8/ 5X,&
   '    EQ.0,  MECHANICAL INPUT DATA CHECK ONLY            ',  / 5X,&
   '    EQ.1,  MECHANICAL RANK CHECK ONLY                  ',  / 5X,&
   '    EQ.2,  MECHANICAL SOLUTION ONLY                    ',  / 5X,&
   '    EQ.3,  THERMAL INPUT DATA CHECK ONLY               ',  / 5X,&
   '    EQ.4,  THERMAL RANK CHECK ONLY                     ',  / 5X,&
   '    EQ.5,  THERMAL SOLUTION ONLY                       ',  / 5X,&
   '    EQ.6,  FULL THERMO-MECHANICAL SOLUTION             ',  //5X,&
   ' RESTART CODE . . . . . . . . . . . . . . . .(ISTART) =',I8/5X,&
   '    EQ.0,  NO RESTART                                  ',  / 5X,&
   '    GT.0,  NO RESTART, WRITE RESTART FILE EACH'/5X,&
   '           ISTART STEPS'/5X,&
   '    LT.0,  RESTART, WRITE RESTART FILE EACH'/5X,&
   '           ABS(ISTART) STEPS'/5X,&
   ' SKIP ELASTIC SOLUTION SWITCH . . . . . . . .(NELAST) =',I8/)

if (OUTFIL(1)) write(luout(1),42)
#ifdef SPARSE
42   format(6X,'SOLVER TYPE ',20(1X,'.'),' = SPARSE'/1X)
#else
42   format(6X,'SOLVER TYPE ',20(1X,'.'),' = DIRECT'/1X)
#endif
if (OUTFIL(1)) write(luout(1),45) ICVIS
45   format(6x, &
 'MATERIAL PARAMETER SWITCH  . . . . . . . . . (ICVIS) =',4X, &
 I0.4)
if (OUTFIL(1)) write(luout(1),46)
46   format(9x, &
 'FIRST DIGIT POWERLAW CREEP (0=OFF, 1=ON)',/,9X, &
 'SECOND DIGIT VISCOPLASTIC CREEP (0=OFF, 1=ON)',/,9X, &
 'THIRD DIGIT VISCOUS WEAKENING (0=OFF, 1=STRAIN, 2=RATE)',/,9X, &
 'FOURTH DIGIT BRITTLE WEAKENING (0=OFF, 1=STRAIN, 2=RATE)')
if(OUTFIL(1)) WRITE(luout(1),50) NINTG,MAXIT,NLINK

50   FORMAT(/5X,&
   ' NUMBER OF CHANGES IN TIME STEP . . . . . . . (NINTG) =',I8/ 5X,&
   '    EQ.0,  ELASTIC SOLUTION ONLY                       ', // 5X,&
   ' NUMBER OF STEPS BETWEEN RFF . . . . . . . . .(MAXIT) =',I8//5X,&
   ' NUMBER OF LINKED NODAL DOF''S  . . . . . . . .(NLINK) =',I8/)

IF(OUTFIL(1)) WRITE(luout(1),55) NMPRT,FaultVelo

55   FORMAT(5X, &
   ' NUMBER OF MECHANICAL SOLUTION OUTPUTS IS . . (NMPRT) =',I8/ 5X, &
   '    GT.0, OUTPUT AT SELECTED TIME STEPS                ',  / 5X, &
   '    EQ.0, NO OUTPUTS                                   ',  / 5X, &
   '    LT.0, OUTPUTS AT EVERY ABS(NMPRT) STEP             ',  / 5X, &
   ' INTERPRET FAULT INPUT AS VELOCITY: . . . (FaultVelo) =',I8/ 5X, &
   '    EQ.0, FAULT INPUT IS DISPLACEMENT',/5X, &
   '    EQ.1, FAULT INPUT IS VELOCITY')

if (OUTFIL(1).and.NMPRT.GT.0) write(luout(1),60)
60   format(/5x, &
   ' SELECTED OUTPUTS OCCUR AT THE FOLLOWING TIMESTEPS:'/)

if (NLINK.gt.0 .and. getrank().gt.0) then
    write(stderr,*) 'Linked nodes not implemented in parallel code'
    call exitp(1)
endif

if (iecho.eq.8) then
    call ParallelLog("Readfe", "puttig out")
endif

call PUTOUT (outputcontroldatactx%IMPRINT,NMPRT)
!
!IPMPRM = IPMPRT + MAX(1,NMPRT)

if (iecho.eq.2) write(stderr,65)

65   format(1x,'*** read IMATPR data ***'/ &
 1x,'if (NPRMAT > 0) {'/ &
 5x,'read (51i5) (IMATPR(i), i=1,NPRMAT)'/1x,'}')


if (iecho.eq.8) then
    call ParallelLog("Readfe", "doing rdmprt")
endif

!##### read when to output matrix diagonal output ######
call RDMPRT (outputcontroldatactx%IMATPR,NPRMAT)
!
if (iecho.eq.8) then
    call ParallelLog("readfe", "done rdmprt")
endif


if (OUTFIL(1)) then
    WRITE(luout(1),70) NPRMAT
endif

70   FORMAT(/5X, &
   ' NUMBER OF MATRIX DIAGONAL OUTPUTS IS . . .  (NPRMAT) =',I5/ 5X, &
   '    GT.0, OUTPUT AT SELECTED TIME STEPS                ',  / 5X, &
   '    EQ.0, NO OUTPUTS                                   ',  / 5X, &
   '    LT.0, OUTPUTS AT EVERY ABS(NPRMAT) STEP             '       )

if (OUTFIL(1).and.NPRMAT.GT.0) then
    write(luout(1),60)
endif

call PUTOUT (outputcontroldatactx%IMATPR,NPRMAT)
!
if (NUMNP.le.0) then
    write(stderr,'(1x,"no nodes")')
    call xit(1,"no nodes")
endif


if (NINTG.ne.0 .and. MAXIT.eq.0) then
    write(stderr,'(1x,"MAXIT=0 error")')
    call xit(1,"MAXIT=0 error")
endif

if (NELAST.ne.0 .and. ISTART.lt.0) then
    write(stderr, &
    '(1x,"Conflicting NELAST and ISTART settings")')
    call xit(1,"Conflicting NELAST and ISTART settings")
 endif



if (iecho.eq.8) then
    call ParallelLog("readfe", "Allocating mech model data")
endif


!write(stdout,*) "allocating mechinal modeldata with numnp: ", numnp

call allocatemechmodeldata()


!###########################################################################
!            R E A D    T H E      M E S H
!###########################################################################

!########### read coordinates (typically: tecin.dat.partf.nps) #############
if (iecho.eq.8) then
    call ParallelLog("readfe", "Reading coordinates")
endif

call COORD (NSD,meshdatactx%nMaxNeighbors)

if (iecho.eq.8) then
    call ParallelLog("readfe", "Finished coord")
endif


!########### read element connectivity (typically: tecin.dat.partf.elm) ######
call ELNODE (NEN)


if (iecho.eq.8) then
    call ParallelLog("readfe", "Finished elnode")
endif

!call BuildElementNeighbourTable2D

!###########################################################################

#ifdef SPARSE

! create vertex and element masks, fill LM arrays, all other mesh administration
! needed for correctly entering the nodal boundary conditions
if (iecho.eq.8) then
    call ParallelLog("readfe", "creating mesh numbering")
endif
call createMeshNumbering(NEN, NDOF)
if (iecho.eq.8) then
    call ParallelLog("readfe", "finished mesh numbering")
endif

#endif

call allocatemechbcdata(modeldatactx, NUMNP, NDOF, NLINK)
if (iecho.eq.8) then
    call ParallelLog("readfe", "finished allocating mech mod data")
endif

!########### read mechanical boundary conditions #############
if (debug) then
    write(stdout,*) "reading nodal boundary conditions"
endif
if (iecho.eq.8) then
    call ParallelLog("readfe", "reading nodal boundary cond.")
endif

call BC (meshdatactx, modeldatactx, NDOF)

if (iecho.eq.8) then
    call ParallelLog("readfe", "finished reading nodal boundary cond.")
endif

!---------------------------------------------------------------------

!call allocatemechwinklerbcdata(modeldatactx, NUMNP, NDOF)
call allocatemechwinklerbcdata(modeldatactx, meshdatactx%Nvglobal, NDOF)

!########### read winkler forces (typically in: tecin.dat.bcs) #############
!call RWINK(modeldatactx%IWINK,modeldatactx%WINK, meshdatactx%X, NSD,NDOF,NUMNP,NWINK,0)
call RWINK(meshdatactx%X, NSD,NDOF,NUMNP,NWINK,0)


if (NWINK.gt.0) then
    call CHKWNK (modeldatactx%IWINK,modeldatactx%IBOND, NDOF,NUMNP)
else
    call deallocatemechwinklerbcdata(modeldatactx)
endif

if (iecho.eq.8) then
    call ParallelLog("readfe", "finished Winkler stuff")
endif


!---------------------------------------------------------------------
!       Maxwell boundary condition data
NMAXW = 0


if (NWINK.gt.0) then
    call allocatemechmaxwellbcdata(modeldatactx, numnp, ndof)
!    call READMAXWELL (modeldatactx%IWINK, modeldatactx%WINK, modeldatactx%FMAXW, NDOF)
    call READMAXWELL (NDOF)

endif
!
if (iecho.eq.8) then
    call ParallelLog("bc", "finished reading Maxwell boundary conditions")
endif


!       Set pointers for skew boundary conditions and read them

!ALLOCATE(modeldatactx%SKEW(2,meshdatactx%Nvlocal),STAT=ERROR)
ALLOCATE(modeldatactx%SKEW(2,meshdatactx%Nvglobal),STAT=ERROR)
do i=1,size(modeldatactx%SKEW,1)
    do j=1,size(modeldatactx%SKEW,2)
        modeldatactx%SKEW(i,j) = 0.0d0
    enddo
enddo


!########### read Euler angles??? ########

if (iecho.eq.8) then
    call ParallelLog("readfe", "reading angles")
endif

call SKEWRD (NSD)


if (iecho.eq.8) then
    call ParallelLog("readfe", "finished angles")
endif


! notice that the timesteproutines are always called,
! even when there are no timesteps.
call allocatetimestepdata(timestepdatactx,0)

!call Verified_MPI_Barrier("yep allocated timestepdata in tecin")

if (debug) then
    write(stdout,*) "Reading timesteps"
endif

call CONTRL(timestepdatactx, outputcontroldatactx)

!call Verified_MPI_Barrier("after contrl")


!call delay()

if (iecho.eq.8) then
    call ParallelLog("readfe", "calling elmlib")
endif

!#########  Read information on element boundary conditions #####
call ELMLIB (NDOF,NSD,NSTR,NEN)

if (iecho.eq.8) then
    call ParallelLog("readfe", "finished elmlib")
endif

NEQ = NEQlocal

if (.not.feplt) then
    allocate(IDIAG(NEQ))
endif

if (do_soelim) then
    if (iecho.eq.2) then
        write(stdout,*) 'removing soelimed tecin file'
    endif
#ifdef SPARSE
    ! if gtecton has run soelim itself, throw away the /tmp/FEM* file
    call MPI_barrier(MPI_COMM_WORLD, ierr)
    if (getrank().eq.0) then
        close (unit=luin, status='delete')
    else
        close (unit=luin)
    endif
#else
    close (unit=luin, status='delete')
#endif
else
    if (iecho.eq.2) then
        write(stdout,*) 'keeping soelimed tecin file'
    endif
    ! if not, than the user has done this him/herself.
    ! and it might be needed for future uses. Do not delete
    close (unit=luin)
endif

if (iecho.eq.1) then
    write(stderr,75) quote,FIN(1:lnblk(FIN)),quote
    write(stderr,75) quote,FEIN(1:lnblk(FEIN)),quote
75 format(1x,'EOF ',80a)
endif

if (iecho.eq.8) then
    if (do_soelim) then
        call ParallelLog("READFE", "Finished reading " // trim(FIN))
    else
        call ParallelLog("READFE", "Finished reading " // trim(FEIN))
    endif
endif

! remove temporary lookup tables

#ifdef SPARSE
if (allocated(meshdatactx%e2p)) then
!    deallocate(meshdatactx%e2p)
endif
if (allocated(meshdatactx%v2p)) then
! noo, still needed when reading thermal input data
!    deallocate(meshdatactx%v2p)
endif
#endif

RETURN

1000    write(stderr,1001) ios,FEIN(1:lnblk(FEIN))
 1001   format(1x,'READFE: fatal input error ',i12,' on ',80a)
call prerr('READFE')
call xit(1,"fatal input error")

1100    write(stderr,1101) ios,FEIN(1:lnblk(FEIN))
 1101   format(1x,'READFE: EOF error ',i12,' on ',80a)
write(stderr,1102)
 1102   format(1x,'READFE: there may be a problem with the scratch', &
 ' file in /tmp or $MYTMP')
call xit(1,"Unexpected end of file")

end subroutine
!
!
!***************************************************************************
!
!       THE SECONDARY SUBROUTINES ARE LISTED BELOW IN ALPHABETICAL ORDER
!
!***************************************************************************

!        subroutine ADJID (ID,IDX,NSLIPglobal,IDSLP,IDSLE,NDOF)
subroutine ADJID (NDOF)

USE MODELDATAMODULE
USE MESHDATAMODULE
USE MODELDEFINITION
USE MODELCTX,        only : getrank, getsize
USE ALGEBRA
use debugmodule,     only : debug, &
                            debugFileID, &
                            allocateError, &
                            iecho, &
                            parallelLog, &
                            xit
#ifdef SPARSE
use aomodule,        only : verticesmask
#endif
use iomodule,       only: stderr
!
implicit none
!
!       Adjusts ID array and creates IDX array for additional degrees
!       of freedom associated with free slip interfaces. Also sets
!       NUMSN and NUMSE.
!
!-pass
integer NDOF
!-locl
logical constraint
integer i,ii,j,n,nn, iError
integer myrank, nexte, nextn

! Matt addition
integer, allocatable :: isThisNodeSlippery(:)
integer, allocatable :: isThisElementSlippery(:)


!-init
NUMSN = 0
NUMSE = 0

if (NUMSLPglobal.eq.0) then
   return
endif


myrank = getrank()
!        call CLEAR(IDX,NDOF*NUMNP)
modeldatactx%IDX = 0
modeldatactx%IDXglobal = 0
NUMSN = 0
NUMSE = 0
NUMSNglobal = 0
NUMSEglobal = 0
nexte = 0
nextn = 0


!if (getsize().gt.1) then
allocate(isThisNodeSlippery(meshdatactx%Nvglobal), stat = iError)
if (ierror.ne.0) then
    call allocateError("tecin isThisNodeSlippery", iError)
endif
!endif

if (iecho.eq.8) then
    call ParallelLog("adjid", "making slippery node look-up table")
endif

do j = 1,meshdatactx%Nvglobal
    isThisNodeSlippery(j) = 0
enddo

do j = 1,NUMSLPglobal
        isThisNodeSlippery(modeldatactx%NSLIPglobal(2,j)) = j
enddo

if (iecho.eq.8) then
    call ParallelLog("adjid", "finished slippery node look-up table")
endif



do n=1,meshdatactx%Nvglobal


!old way

!!   Check if node N occurs in NSLIPglobal
!!   which is equivalent to checking if there is a slippery node entry 
!!   of this node in the TECIN file.
!if (getsize().gt.1) then
!    do j=1,NUMSLPglobal  !modeldatactx%NUMSLP
!        if (n.eq.modeldatactx%NSLIPglobal(2,j)) then
!            goto 100
!        endif
!    enddo
!else
! new way

    if (isThisNodeSlippery(n).ne.0) then
!        write(stdout,*) "rank", getrank(),"looking up slip", n, "of", size(isThisNodeSlippery, 1)
        j = isThisNodeSlippery(n)
        goto 100
    endif
!endif


!           Node N does not occur in NSLIPglobal: next node
    cycle

!
!           Node N does occur in NSLIPglobal, add to appropriate owner of N!!!!
100         numsnglobal = numsnglobal + 1


! TODO IS NOT GOOD... Slippery nodes must be known outside partition
!#ifdef SPARSE
!    if (verticesmask(meshdatactx%Glo2LocVertex(n)).ne.0) then
!#else
!    if (meshdatactx%v2p(n).eq.myrank) then
!#endif
        nextn = nextn + 1
        NUMSN = nextn

        if (debug) then
             write(debugFileID,*) "SETTING IDSLP", nextn, "to ", n
        endif

        modeldatactx%IDSLP(nextn) = n
!    endif

    constraint = .false.

    do i=1,NDOF
!       Constraint on degree of freedom i?
        if ((modeldatactx%NSLIPglobal(2+i,j).EQ.0).OR. &
            (modeldatactx%IDglobal(i,n)     .EQ.0)) then
             cycle
        endif
        
        constraint = .true.

#ifdef ZISFAULTNORMAL
        if (i.eq.NDOF) then
            write(stderr,*) "NSLIP entry number ",j," has nonzero idof=",i
            write(stderr,*) "This is inconsistent with the app because it", &
             "was compiled with -DZISFAULTNORMAL"
            call xit(1,"Exiting GTecton ...")
        endif
#endif

!LWITODO perhaps not do this for every DOF... check this
        ! copy IDglobal into locally owned IDX!

!!!!!! HIERRRRRRRRRRRRRRRRRRRRRRRRRRRRRR !!!!!!!
!        if (meshdatactx%v2p(n).eq.myrank) then
!#ifdef SPARSE
!        if (verticesmask(meshdatactx%Glo2LocVertex(n)).eq.0) then
!#else
        if (meshdatactx%v2p(n).eq.myrank) then
!#endif
#ifdef SPARSE
            if (getsize().eq.1) then
                nn = n
            else
!                nn = meshdatactx%Glo2LocVertex(n)
                nn = verticesmask(meshdatactx%Glo2LocVertex(n))
            endif
#else
            nn=n
#endif
            
!            write(stdout,*) 'rank ', myrank, 'from IDglobal: ', modeldatactx%IDglobal
!            write(stdout,*) 'rank ', myrank, 'add vertex ', n, nn, 'NDOF', i, modeldatactx%IDglobal(i,n)

            modeldatactx%IDX(i,nn) = modeldatactx%IDglobal(i,n)!ID(i,n)
            NEQlocal = NEQlocal + 1
        endif

!       write(stdout,*) 'rank ', myrank, 'resulting in IDX: ', modeldatactx%IDX

        ! copy IDglobal into global IDX, needed for LMX later
        modeldatactx%IDXglobal(i,n) = modeldatactx%IDglobal(i,n)

!       Modify equation numbers ID for this particular node
!       to make room for the slippery node equation.
        do ii=i,NDOF
            if (modeldatactx%IDglobal(ii,n).ne.0) then
                modeldatactx%IDglobal(ii,n)= modeldatactx%IDglobal(ii,n) + 1 
            endif
        enddo

        ! increase all equation numbers of higher vertices
        ! to make room for the equation of the slippery node.
        do nn=n+1,meshdatactx%Nvglobal
            do ii=1,NDOF
                if (modeldatactx%IDglobal(ii,nn).ne.0) then
                    modeldatactx%IDglobal(ii,nn) = modeldatactx%IDglobal(ii,nn) + 1
                endif
            enddo
        enddo
       
        NEQ = NEQ + 1
        NEQglobal = NEQglobal + 1

    enddo

    if (.not.constraint) then
!        write(stderr,1) j,modeldatactx%NSLIPglobal(2,j)
! 1           format(1x,'NSLIP(',I3,') (node=',i8,') FATAL problem'/ &
!         1x,'Combining a displacement or velocity bc with ', &
!         'a slippery dof in the same'/1x,'direction is ', &
!         'currently not possible')
        write(stderr,*) "node", n, "is both slippery and has a velocity BC"
        write(stderr,*) "That is currently not possible."
        write(stderr,*) "(Coordinates: ", meshdatactx%X(1,n), meshdatactx%X(2,n), ")"
        call xit(1,"Exiting GTecton...")
    endif

enddo


if (allocated(isThisNodeSlippery)) then
    deallocate(isThisNodeSlippery)
endif


allocate(isThisElementSlippery(meshdatactx%Neglobal), stat = iError)
if (ierror.ne.0) then
    call allocateError("tecin isThisElementSlippery", iError)
endif

if (iecho.eq.8) then
    call ParallelLog("adjid", "making slippery element look-up table")
endif

do j = 1,meshdatactx%Neglobal
    isThisElementSlippery(j) = 0
enddo

do j = 1,NUMSLPglobal
    isThisElementSlippery(modeldatactx%NSLIPglobal(1,j)) = j
enddo

if (iecho.eq.8) then
    call ParallelLog("adjid", "finished slippery element look-up table")
endif


do n=1,meshdatactx%Neglobal!NUMEL

! old way
!if (getsize().gt.1) then
!    do j=1,NUMSLPglobal!modeldatactx%NUMSLP
!        if (n.eq.modeldatactx%NSLIPglobal(1,j)) goto 400
!    enddo
!else
! NEW WAY
    if (isThisElementSlippery(n).ne.0) then
        goto 400
    endif
!endif

    cycle

400         NUMSEglobal = NUMSEglobal + 1
    if (meshdatactx%e2p(n).eq.myrank) then
        nexte = nexte + 1
        NUMSE = nexte
        modeldatactx%IDSLE(nexte) = n
    endif
enddo
!

if (allocated(isThisElementSlippery)) then
    deallocate(isThisElementSlippery)
endif


do i=1,meshdatactx%nvlocal
    do j=1,NDOF
     if (modeldatactx%IDX(j,i) .ne. modeldatactx%IDXglobal(j,i)) then
!                write(stdout,*) 'argh, IDX IDXglobal wrong', i, j, & 
!                            modeldatactx%IDX(j,i), modeldatactx%IDXglobal(j,i)
     endif
     if (modeldatactx%ID(j,i) .ne. modeldatactx%IDglobal(j,i)) then
!                write(stdout,*) 'argh, ID IDglobal wrong', i, j, & 
!                            modeldatactx%ID(j,i), modeldatactx%IDglobal(j,i)
     endif

    enddo
enddo


if (iecho.eq.8) then
    call ParallelLog("adjid", "finshed sub")
endif


return
end
!-------------------------------------------------------------------------------
subroutine BC (meshdat, modeldat, NDOF)

USE LABELMODULE
USE MESHDATAMODULE
USE MODELDATAMODULE
USE MODELDEFINITION
USE MODELCTX
USE ALGEBRA
USE FORMATSMODULE, only: FMT_tecin_dat_bcs_ibc_2D,  &
                         FMT_tecin_dat_bcs_ibc_3D,  &
                         FMT_tecin_dat_bcs_bc_2D,   &
                         FMT_tecin_dat_bcs_bc_3D,   &
                         FMT_tecin_dat_bcs_iwink,   &
                         FMT_tecin_dat_bcs_wink,    &
                         FMT_tecin_dat_bcs_euler,   &
                         FMT_tecin_dat_bcs_LinkedNodes
use debugmodule,   only: debug, &
                         iecho, &
                         ParallelLog, &
                         Verified_MPI_Barrier, xit
use constants,     only: isOpn, isPln, isAxi, isF3d
use iomodule          
!
! Program for reading and printing boundary condition data and
! establishing equation numbers. 
! It reads the tecin.dat.bcs file, as groffed into tecin.dat
!

implicit none
integer, parameter :: MAXBC=6
!-pass
type (meshdata)  :: meshdat
type (modeldata) :: modeldat
integer NDOF


!-locl
logical            CONSTRAINT
character(len=4)   LTYPE(MAXBC)
character(len=255) LEGEND(ndof)
! nc counts the number of bcs in each of the coordinate directions.
! It is used to detect whether any direction has no bcs.
! The model could slip entirely in that direction.
integer            nc(3),isum, nloc
integer ( kind = 4 ) h,m,s,values(8) ! for timestamps in output

integer :: nn, n, i, j, k ! loop counters
integer :: itype, idum, nread
!-init
data LTYPE/'FREE','DISP','VELO','FORC','INIT','RATE'/
data nc/0,0,0/

call CLEAR(modeldat%BOND,NDOF*meshdat%Nvlocal,"modeldat%BOND")
nn = 51
!
if (iecho.eq.8) then
    call ParallelLog("bc", "entered bc")
endif

if (iecho.eq.2) then
    write(stderr,10) FMT_tecin_dat_bcs_ibc_2D
10   format(1x,'*** Read boundary condition codes ***'/ &
    1x,'n=1; while (n != 0) {'/ &
    5x,'read ',A9,' n,(IBOND(i,n),i=1,NDOF)'/ &
    1x,'}')
endif
!
call CLEAR(modeldat%IBOND,       NDOF*meshdat%Nvlocal,"modeldat%IBOND")
call CLEAR(modeldat%IBONDglobal, NDOF*meshdat%Nvglobal,"modeldat%IBONDglobal")

!###########  read boundary condition codes  ###############

if (iecho.eq.8) then
!    call ParallelLog("bc", "reading codes")
endif

if (NDOF.eq.2) then
    ! modeldat%IBOND is being filled according to local
    ! numbering of the vertices.
    call IGEN(modeldat%IBOND, &
              NDOF, &
              meshdat%nvglobal, &
              FMT_tecin_dat_bcs_ibc_2D, &
              nread, &
              1, &
              "Boundary condition indices; 2D")
else if (NDOF.eq.3) then
    call IGEN(modeldat%IBOND, &
              NDOF, &
              meshdat%nvglobal, &
              FMT_tecin_dat_bcs_ibc_3D, &
              nread, &
              1, & 
              "Boundary condition indices; 3D")
else
    write(stdout,*) 'NDOF should be either 2 or 3; contact model support.'
endif

#ifdef Lukas
!write(stdout,*) 'rank ', getrank(), 'after calling IGEN, IBOND: ', modeldat%IBOND
#endif

if (iecho.eq.2) then
    write(stderr,20) FMT_tecin_dat_bcs_bc_2D
20   format(1x,'*** Read boundary condition data ***'/ &
    1x,'n=1; while (n != 0) {'/ &
    5x,'read ',A12,' n,(BOND(i,n),i=1,NDOF)'/ &
    1x,'}')
endif
!

!######## read boundary condition values ###########

if (iecho.eq.8) then
!    call ParallelLog("bc", "reading magnitudes")
endif


if (NDOF.eq.2) then
    ! modeldat%BOND is being filled according to local
    ! numbering of the vertices.
    call FGEN(modeldat%BOND,NDOF,meshdat%nvglobal,FMT_tecin_dat_bcs_bc_2D,idum, &
              "nodal boundary conditions, magnitudes, 2D")
else
    call FGEN(modeldat%BOND,NDOF,meshdat%nvglobal,FMT_tecin_dat_bcs_bc_3D,idum, &
              "nodal boundary conditions, magnitudes, 3D")
endif


#ifdef Lukas
write(stdout,*) 'rank ', getrank(), 'after calling FGEN, BOND: ', modeldat%BOND
#endif


if (iecho.eq.8) then
    call ParallelLog("bc", "done magnitudes")
endif

if (NLINK.gt.0) then

    if (iecho.eq.2) then
        write(stderr,30) FMT_tecin_dat_bcs_LinkedNodes
30      format(1x,'*** Read linked node data ***'/ &
         1x,'if (NLINK > 0) {'/ &
         5x,'idof=1; while (idof != 0 & idof != end) {'/ &
         9x,'read ',A8,' idof,LinkNode,LinksToNode,LinkNodeIncrement'/ &
         5x,'}'/1x,'}')
    endif

!   read extra boundary condition codes
    call LNKGEN (modeldat%LINK,NLINK,NDOF,NUMNP)

endif

if (iecho.eq.8) then
    call ParallelLog("bc", "cleared linkgen")
endif

if (feplt) then
    return
endif

! establish global equation numbers for temporary use

NEQglobal = 0
NEQlocal = 0

do n=1,meshdat%nvlocal
    do i=1,NDOF
        modeldat%ID(i,n) = 0
    enddo
enddo

call CLEAR(modeldat%IDglobal, NDOF*meshdat%Nvglobal,"modeldat%IDglobal")

if (iecho.eq.8) then
    call ParallelLog("bc", "modifying eqns numbers")
endif

if (iecho.eq.2) then
    write(stdout,*) '*** modifying equation numbers due to nodal BCS'
endif

!    build equation numbers, based on nodal boundary conditions
do n=1,meshdat%Nvglobal

    if (iecho.eq.2 .and. mod(n,99999) .eq. n/99999) then
        ! give a progress report every 100k vertices.
        ! otherwise some uses may report the software 'hangs' for 
        ! very large meshes, because this can take a while.
        call date_and_time ( values = values )
        h = values(5)
        m = values(6)
        s = values(7)
        write(stdout,111) n, meshdat%nvglobal, h,m,s
111     format(1x,'progress: ',i12,' of ',i12, ' on time ', i2, ':',i2,':',i2)
    endif

!    if (debug) then
!        write(stdout,*) 'vertex ',n, ' has ',(modeldatactx%IBONDglobal(i,n),i=1,NDOF)
!    endif

    do i=1,NDOF

!        if (debug) then
!            write(stdout,*) 'IBONDglobal: ',i,n, modeldatactx%IBONDglobal(i,n), &
!             modeldatactx%IBOND(i,n), modeldat%IBOND(i,n)
!        endif

        if (modeldat%IBONDglobal(i,n).eq.1 .or.  modeldat%IBONDglobal(i,n).eq.2 .or. &
         modeldat%IBONDglobal(i,n).eq.5) then 
            ! these boundary conditions give no degree of freedom, and as such no extra equation number
            modeldat%IDglobal(i,n) = 0
            nc(i) = nc(i) + 1
        else
            ! other boundary conditions (force (3) or free (0)) do not take a degree
            ! of freedom away. Hence NEQ[global/local] up

!           Check if IDOF,NODE occurs as slave entry in LINK array
            k = INLINK (i,n)
            if (k.le.0) then
                ! node is not a slave linked to a master node, so it gets its own equations.
                NEQglobal = NEQglobal + 1
                if (meshdat%v2p(n).eq.getrank()) then
                    NEQlocal = NEQlocal + 1
!                    if (debug) write(stdout,*) 'size of procs: ', getsize()

#ifdef SPARSE
                   if (getsize().eq.1) then
                       nloc = n
                   else
                       nloc = meshdat%Glo2LocVertex(n)
                   endif

!                   if (debug) write(stdout,*) "rank", getrank(), "thinks that node", n , "has local number", nloc

#else /* sparse */
                   ! this is a one processor run.
                   ! global and local numbering are the same.
                   nloc = n
#endif /* sparse */ 

                   modeldat%ID(i,nloc) = NEQglobal
                endif

                modeldat%IDglobal(i,n) =  NEQglobal
            else
                if (iecho.eq.5) then
                    write(stdout,*) 'Adding no equation for node ', n, &
                           'd.o.f.', i, &
                           'because node is slave node in a link with master',k
                endif

            endif
        endif

    enddo ! dofs loop

enddo ! nodal points loop

!call Verified_MPI_Barrier("renumbered eq")

if (debug) then
    if (parout) then
        write(FILE_outputf,*) 'after reading boundary conditions and allocation of equation numbers'
        write(FILE_outputf,*) '    global array w/ equation numbers:'
        do i=1,meshdat%Nvglobal
            write(FILE_outputf,*) i, (modeldat%IDglobal(j,i),j=1,NDOF)
        enddo
        write(FILE_outputf,*) '    local array w/ equation numbers (vertices still have global numbers):'
        do i=1,meshdat%Nvlocal
            write(FILE_outputf,*) meshdat%gloInd(i), (modeldat%ID(j,i),j=1,NDOF)
        enddo
    endif
endif

if (MODE.le.2 .or. MODE.eq.6) then

! check for each of the coordinate directions by how many positional BCS it is limited.

    ! TODO in axisymmetric code, with a hole around the symmetry axis, the model 
    ! cannot be fixed there, but that should be fine...  is that so?
    ! Radial component of hoop stress constrains the motion!
    if (simVariety .eq. isPln) then
        do i=1,NDOF
            if (nc(i).lt.1) then
                write(stderr,35) i
35              format(1x,"Warning: No displacement or velocity boundary condition on degree of freedom ",I1)
                write(stderr,*) "Only Winkler forces and pressures can restrict the motion..."
!            call xit(1,"No displacement or velocity boundary condition on one DOF")
            endif
        enddo
    else if (simVariety .eq. isOpn) then
        ! OPN, we have ndof 3, but no constraint in the third
        do i=1,2  ! only check the first 2 degrees of freedom
            if (nc(i).lt.1) then
                write(stderr,35) i
                write(stderr,*) "Only Winkler forces and pressures can restrict the motion..."
            endif
        enddo
    else if (simVariety .eq. isF3d) then
        do i=1,3  ! only check the first 2 degrees of freedom
            if (nc(i).lt.1) then
                write(stderr,35) i
                write(stderr,*) "Only Winkler forces and pressures can restrict the motion..."
            endif
        enddo
    else if (simVariety .eq. isAxi) then
        if (nc(1).lt.1) then
            write(stderr,35) 1
            write(stderr,*) "Only hoop stress, Winkler forces and pressures can restrict the motion..."
        endif
    endif

    isum = 0

    do i=1,NDOF
        if (nc(i).ge.2) isum = isum + 1
    enddo

    if (isum.lt.NDOF-1) then
        write(stderr,36)
36      format(1x,'Displacement/velocity boundary conditions ', &
         'don''t constrain rigid body rotation of the entire ', &
         'model')
        call xit(1,"Entire model allowed to rotate; not enough positional BCs")
    endif

endif
!
if (nread.le.0) then
    if (iecho.eq.2) then
        write(*,*) 'rank', getrank(), 'has no nodal bcs. Returning'
    endif
    return
endif
!
!       checking
if (iecho.eq.2) then
    write(*,*) '*** checking boundary condition integrity ***'
!            write(*,*) 'rank', getrank(),'looping over ', meshdat%Nvlocal
endif

do n=1,meshdat%Nvlocal

    ! progress-o-meter, because it can take a long time
    if (iecho.eq.2 .and. mod(n,99999) .eq. n/99999) then
        call date_and_time ( values = values )
        h = values(5)
        m = values(6)
        s = values(7)
        write(*,222) n, meshdat%nvglobal, h,m,s
222             format(1x,'a progress: ',i12,' of ',i12, ' on time ', i2, ':',i2,':',i2)
    endif

    do i=1,NDOF
        if (modeldat%IBOND(i,n) .lt. 0) then
            write(stderr,40) i,n
40               format(1x,'BC: IBOND(',i1,',',I5,') < 0')
            call xit(1,"IBOND")
        else if (modeldat%IBOND(i,n) .gt. MAXBC) then
            write(stderr,50) i,n,MAXBC-1
50               format(1x,'BC: IBOND(',i1,',',I5,') > ',I1)
            call xit(1,"IBOND")
        endif
    enddo

enddo


if (iecho.eq.2) then
   write(*,*) 'rank', getrank(), &
              '*** finished checking boundary condition integrity ***'
endif

!if (debug) then
!     write(*,*) "rank", getrank(), " has finished BC with modeldata ID: ", modeldat%ID
!endif

!
if (.not.OUTFIL(1)) then
    return
endif
!
do n=1,meshdat%Nvglobal
    if (nn.gt.50) then
        nn=1
        write(luout(1),60)
60           format(1x///2x,'N O D A L',3x,'F O R C E S',3x,'A N D', &
         3x,'D I S P L A C E M E N T S'// &
         6x,'KEY TO BOUNDARY CONDITION CODES:'// &
         10x,'FREE = UNCONSTRAINED DEGREE OF FREEDOM;'/ &
         17x,'THE NUMERICAL VALUE TO THE LEFT IS MEANINGLESS.'/ &
         10x,'DISP = FIXED DISPLACEMENT'/ &
         10x,'VELO = CONSTANT VELOCITY'/ &
         10x,'INIT = INITIAL DISPLACEMENT'/ &
         10x,'FORC = APPLIED FORCE'/ &
         10x,'RATE = CONSTANT STRAIN RATE'//)
        write(luout(1),70) (LABELD(I),I=1,NDOF)
70                format(4X,'NODE # ',6(4X,A4,3X,'CONDITION  ')/)
        write(luout(1),80)
80                format(/)
    endif
    CONSTRAINT = .false.

    do i=1,NDOF
        if (modeldat%IBOND(i,n) .ne. 0) then
            CONSTRAINT = .true.
        endif
    enddo

    if (CONSTRAINT) then
        do i=1,NDOF
            itype = modeldat%IBOND(i,n)+1
            LEGEND(i) = LTYPE(itype)
        enddo
        write(luout(1),90) n,(modeldat%BOND(i,n), &
                         LEGEND(i),i=1,NDOF)
90           format(1X,I8,2X,6(1PE12.5,2X,A4,4X))
        nn=nn+1
    endif
enddo
!
if (iecho.eq.8) then
    call ParallelLog("bc", "finished BC")
endif

return
end subroutine BC

!-----------------------------------------------------------------------

subroutine CHKWNK (IWINK,IBOND,NDOF,NUMNP)

use debugmodule, only: iecho
use iomodule          
!
!       input checks on nodal Winkler bcs
!
implicit none
!-pass
integer   :: NDOF,NUMNP
integer   :: IWINK, IBOND
dimension :: IWINK(NDOF,*),IBOND(NDOF,*)
!-locl
integer n,i
!
do n=1,NUMNP
    do i=1,NDOF
        if ((IBOND(i,n).eq.1 .or. IBOND(i,n).eq.2) &
         .and. IWINK(i,n).ne.0 .and. iecho.gt.0) then
             if (iecho.eq.2) then
                 write(stderr,10) i,n
10               format(1x,'Warning: nodal Winkler condition on dof', &
             1x,I1,' of node ',I8,' ignored, since there also ', &
             'is a displacement/velocity bc on this node')
             endif
        endif
    enddo
enddo
!
return
end subroutine
!-----------------------------------------------------------------------
subroutine CHKSLP (IEN,NSLIP,NUMSLP,NEN,NUMEL)

use modelctx,        only: getrank, &
                           getsize
USE MESHDATAMODULE
USE MODELDATAMODULE
use debugmodule,     only: xit
use iomodule,        only: stderr
!
!       Performs a couple of simple checks on slippery node definitions
!

implicit none
!-pass
integer NUMSLP,NEN,NUMEL
integer :: NSLIP, IEN
dimension NSLIP(5,*),IEN(NEN,NUMEL)
!-locl
integer  :: i, ii, in, nn, j, k, ielm, inod
integer, parameter :: MAXN=150000
logical InNodes, error
integer :: nodes(MAXN)
integer neloc, nelm, lnod



if (.not. HasSlipperyNodes()) then
    return
endif
!
error = .false.
!
!       Check if slippery nodes actually occur in specified elements
!       loop through the TECIN.DAT entries

do i=1,NUMSLP

    nelm = NSLIP(1,i) ! this is the global number.

    lnod = 1

#ifdef SPARSE
    if (getsize().eq.1) then
        neloc = nelm
    else
!        write(*,*) 'rank', getrank(), i, 'checking array ', nelm
        write(*,*) 'rank', getrank(), 'global ', nelm, ' is local: ', meshdatactx%Glo2LocElement(nelm)
        neloc = meshdatactx%Glo2LocElement(nelm)
!        write(*,*) 'rank', getrank(), i, 'Glo2LocElement ', neloc
    endif
#else
    neloc = nelm
#endif


    !----------------------
    ! check whether the element node combination exists in the mesh
   
    ! IEN is stored locally, but slip is in global element ID,
    ! so this goes wrong. Commented out for now but should be fixed; TODO

100  if (lnod.gt.NEN) then
        ! all vertices of this element have been checked, and the one
        ! for this entry was not among them. Input data is messed up, kind sir.
        write(*,*) "ERROR: Slippery entry", i, ". Node", NSLIP(2,i), "does not occur in element", nelm

!        write(stderr,2) i,NSLIP(2,i),nelm
!        write(*,2) i,NSLIP(2,i),nelm
! 2           format(1x,'RSLIP entry ',I12,': node ',I12,' does not', &
!         ' occur in element ',I12)
        call xit(1,"Slippery configuration incorrect")
    endif

    !   NSLIP(1,i) = element;  NSLIP(2,i) = vertex
    if (NSLIP(2,i).ne.IEN(lnod,neloc)) then
        lnod = lnod + 1
        goto 100
    endif
!----------------------

enddo


!       Check if each specified slippery node is defined as slippery
!       in all elements in which it occurs.
!       Set up list of nodes occuring in NSLIP.
nodes(1) = NSLIP(2,1)
nn = 1


do i=2,NUMSLP

    in = 1
    InNodes = .false.

300 if (in.le.nn) then
!       Search nodes-array to find out if NSLIP-node already
!       is in it.
        if (nodes(in).eq.NSLIP(2,i)) then
            InNodes = .true.
            goto 400
        endif
        in = in + 1
        goto 300
    endif

400 if (.not.InNodes) then
!       Not in nodes-array: add it.
        nn = nn + 1
        if (nn.gt.MAXN) then
            write(stderr,3) MAXN
 3               format(1x,'CHKSLP: scratch memory size = ',I3)
            call xit(1,"Too many slippery nodes")
        endif
        nodes(nn) = NSLIP(2,i)
    endif

enddo



!       Loop nodes
do in=1,nn
    do 600 j=1,meshdatactx%Neglobal !NUMEL

!               Find element in which node nodes(in) occurs
        do ii=1,NEN
            i = ii
            if (IEN(i,j).eq.nodes(in)) then
                goto 500
            endif
        enddo

        goto 600

!       Yep, it is in element J. Now check that element J
!       and node nodes(in) also appear as an entry in NSLIP
500     do k=1,NUMSLP
            if (NSLIP(1,k).eq.meshdatactx%gloElt(j).and. &
                NSLIP(2,k).eq.IEN(i,j)) then
                goto 600
            endif
        enddo

        write(stderr,4) getrank(),IEN(i,j),meshdatactx%gloElt(j)
 4           format(1x,'rank ', I5,' says: CHKSLP: node ',I12,' defined as slippery,',1x, &
         'but element ',I12,', which is'/1x,'not in the',1x, &
         'slippery nodes list, also contains this node')

        error = .true.

600 continue
enddo


!       Verify that double element/node entries do not occur
do i=1,NUMSLP-1

   ielm = NSLIP(1,i)
   inod = NSLIP(2,i)

   do j=i+1,NUMSLP
        if (NSLIP(1,j).eq.ielm .and. NSLIP(2,j).eq.inod) then
            write(stderr,5) i,j
 5               format(1x,'CHKSLP: NSLIP entries ',I12,' and ',I12, &
             ' refer to the same element/node combination')
!                    call xit(1," ")
        endif
   enddo

enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine COLHT (IDIAG,ID,IEN,LM,LMX,NDOF,NEN,NUMEL)
!
! Computes column heights
!
!USE MODELDATAMODULE, only: IDIAG
!USE MESHDATAMODULE
USE MODELDEFINITION, only: NUMSN, NEQ
use constants,       only: intmax
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: elementSlipperyID
#endif

implicit none

!-pass
integer NDOF,NEN, NUMEL
integer IDIAG(NEQ),LM(NDOF,NEN,*),LMX(NDOF,NEN,*),IEN(NEN,*), ID(NDOF,*)
!-locl
integer k,minimum,i,j,m,numeq
!


do k=1,NUMEL
    minimum = intmax
!   find smallest equation number in current element
    do j=1,NEN
        do i=1,NDOF

            numeq = LM(i,j,k)
            if (numeq.ne.0) then
                minimum = MIN0(minimum,numeq)
            endif
#ifdef EXPERIMENTAL_ReducedSlipAlloc
            if(elementSlipperyID(k,100).gt.0) then
                if (NUMSN*LMX(i,j,elementSlipperyID(k,100)).ne.0) then
                    minimum = MIN0(minimum,IABS(LMX(i,j,elementSlipperyID(k,100))))
                endif
            endif
#else
            if (NUMSN*LMX(i,j,k).ne.0) then
                minimum = MIN0(minimum,IABS(LMX(i,j,k)))
            endif
#endif

        enddo
    enddo

!   Set column heights
    do j=1,NEN
        do i=1,NDOF
            numeq = LM(i,j,k)
            if (numeq.gt.0) then
                m = numeq-minimum
                if (m.gt.IDIAG(numeq)) then
                    IDIAG(numeq) = m
                endif
            endif
        enddo
    enddo

    if (NUMSN.gt.0) then
        do j=1,NEN
            do i=1,NDOF
#ifdef EXPERIMENTAL_ReducedSlipAlloc
                numeq = IABS(LMX(i,j,elementSlipperyID(k,101)))
#else
                numeq = IABS(LMX(i,j,k))
#endif
                if (numeq.gt.0) then
                    m = numeq-minimum
                    if (m.gt.IDIAG(numeq)) then
                        IDIAG(numeq) = m
                    endif
                endif
            enddo
        enddo
    endif

enddo
!
return
end

!-------------------------------------------------------------------------------
subroutine COORD (NSD,nMaxNeighbors)

USE LABELMODULE
USE MESHDATAMODULE, only: meshdatactx
USE MODELDEFINITION
USE MODELCTX
USE FORMATSMODULE,  only: FMT_tecin_partf_nps_neighbour, &
                          FMT_tecin_partf_nps_1,         &
                          FMT_2D_tecin_partf_nps_2,      &
                          FMT_3D_tecin_partf_nps_2
use debugmodule,    only: iecho, &
                          ParallelLog, &
                          xit, debug
use constants,      only: NNeighborsMax2D, &
                          NNeighborsMax3D
use iomodule          
!
! Reads and prints coordinate data.
! The problem with reading the partitioned nodel points file (typically: tecin.dat.partf.nps)
! is that is not known in advance how many neighbours a certain vertex has.
! Hence first a part of the line must be read. In this part, the neighbourcount can be read.
! USing this count, the whole line can be read.
!
! This routine is both called by Gtecton (who reads the partitioned data,
! and grabs out only the part used by the processor it works on)
! as well as by the plot programs plnpln/plt3d/etc, who need every point,
! independent of partition.

implicit none
!-pass
integer :: NSD
integer :: nMaxNeighbors
!-locl
logical end
character(len=2064) record
character(len=69) fmt
double precision :: f(3)
integer tmp, nread, nextread, nbrs, nbrlist(100), maxNeighbourCount
integer ( kind = 4 ) h,m,s,values(8) ! for timestamps in output
integer :: ios
integer :: i, l, lg, n, nn, idummy
integer :: marker

integer, external :: lnblk ! from tools
integer :: maxNeighbors
integer :: nPartitions

!-init

l=0
lg=0
nread=0
nextread = 0
maxNeighbourCount=0

nPartitions = getsize()

!       typically tecin.dat.partf.nps
! partition, ID, label, x, y, [z], nNeighbors, neighborIDs
if (NSD.eq.3) then
    fmt = '(I5,1X,I12,1X,I12,1X,1E25.17,1X,1E25.17,1X,1E25.17,1X,I4,150(1X,I12))'
else if (NSD.eq.2) then
    fmt = '(I5,1X,I12,1X,I12,1X,1E25.17,1X,1E25.17,1X,I4,150(1X,I12))'
else
    write(*,*) 'NSD should be 2 or 3; contact model support'
endif

if (iecho.eq.8) then
    call ParallelLog("coord", "reading vertex coordinates")
else if (iecho.eq.2) then
    write(stderr,1) fmt
1   format(1x,'*** read nodal coordinate data ***'/ &
    1x,'n=1; while (n != 0){'/ &
    5x,'read ',A69,' partition_number,nodal_point_number,marker,(X(i,n),i=1,NSD)'/ &
    5x,'number_of_neighbor_points neighbor_point_numbers'/ &
    1x,'}')
endif
 
meshdatactx%Nvneighborstotal = 0

!############ read nodal point coordinates line ############
! Every line is read from the file once,
! But it is interpreted multiple times.
100     call readrc (luin,record,'COORD',end)

if (debug) then
    write(*,*) "read coord record: ", trim(record)
endif

if (end) then
    goto 200
endif

if (record(5:5).eq." ") then
    write(*,*) "Input reader says: Something wrong with the partition number of vertex", nread+1
    write(*,*) "Position 5 empty in line:"
    write(*,"(a)") trim(record)
endif

if (record(18:18).eq." ") then 
    write(*,*) "Input reader says: Something wrong with the vertex ID of vertex", nread+1
    write(*,*) "Position 18 empty in line:"
    write(*,"(a)") trim(record)
endif

if (record(31:31).eq." ") then
    write(*,*) "Input reader says: Something wrong with the vertex marker of vertex", nread+1
    write(*,*) "Position 31 empty in line:"
    write(*,"(a)") trim(record)
endif

if (record(57:57).eq." ") then
    write(*,*) "Input reader says: Something wrong with the X coordinate of vertex", nread+1
    write(*,*) "Position 57 empty in line:"
    write(*,"(a)") trim(record)
endif

if (record(83:83).eq." ") then
    write(*,*) "Input reader says: Something wrong with the Y coordinate of vertex", nread+1
    write(*,*) "Position 83 empty in line:"
    write(*,"(a)") trim(record)
endif

if (NSD.eq.3) then
    if (record(109:109).eq." ") then
        write(*,*) "Input reader says: Something wrong with the Z coordinate of vertex", nread+1
        write(*,*) "Position 109 empty in line:"
        write(*,"(a)") trim(record)
    endif
endif

if (debug) then
    write(*,*) "interpreting record with format", trim(fmt)
endif

read (record,fmt,err=300,end=400,iostat=ios) tmp,n,marker,(f(i),i=1,NSD)

if (debug) then
        write(*,*) "interpreted a: ", tmp,n,marker,(f(i),i=1,NSD)
endif

nread = nread + 1

! progress-O-meter
if (iecho.eq.2 .and. mod(nread,99999) .eq. nread/99999) then
     call date_and_time ( values = values )
     h = values(5)
     m = values(6)
     s = values(7)
     write(*,111) nread, meshdatactx%nvglobal, h,m,s
111          format(1x,'progress: ',i12,' of ',i12, ' on time ', i2, ':',i2,':',i2)
endif

! read the rank of this line
read (record,'(I5)',err=300,end=400,iostat=ios) tmp

#ifdef SPARSE
if (tmp .ge. nPartitions) then
    if (getrank() .eq. 0) then
    write(stderr,*) "The simulation is run over", nPartitions, "Partitions."
    write(stderr,*) "Numbering of partitions starts at 0."
    write(stderr,*) "Example: Four partitions are numbered 0, 1, 2 and 3."
    write(stderr,*) "Maximum partition number for this model is", nPartitions-1
    write(stderr,*) "but nodal point ", nread, "has partition", tmp, "which is too high"
    STOP "Leaving GTECTON..."
    endif
endif
#endif

meshdatactx%v2p(nread) = tmp

#ifndef SPARSE
! this code is executed when a plot program is run.
! and makes sure that it reads all the data points, 
! regardless of partition.
tmp = getrank()
#endif

! only read the details of the line when the vertex is in the partition
! of this processor
if (tmp.eq.getrank()) then
    nextread = nextread + 1

    if (debug) then
        write(*,*) "----------------------------"
        write(*,*) "coord says: read record ", trim(record)
    endif

    ! split a record into tmp, n, NSD function values and neighbours
    read (record,fmt,err=300,end=400,iostat=ios) tmp,n,marker,(f(i),i=1,NSD),nbrs

    if (debug) then
        write(*,*) "interpreted b: ", tmp,n,marker,(f(i),i=1,NSD),nbrs
    endif

    if (nbrs.eq.0) then
         write(*,*) 'Error: vertex ', n, ' appears to have no neighbors.'
         if (nextread.eq.1) then
             write(*,*) "perhaps trying to read 3D tecin file with a 2D version of GTecton?"
         endif
         goto 200 ! is error, a vertex always has neighbors.
    endif

    if (nbrs.gt.meshdatactx%maxNeighbourCount) then
        meshdatactx%maxNeighbourCount = nbrs
    endif

    read (record,fmt,err=300,end=400,iostat=ios) tmp,n,marker,(f(i),i=1,NSD),idummy,(nbrlist(i),i=1,nbrs)

    if (nextread.gt.meshdatactx%Nvglobal) then
        write(stderr,13) meshdatactx%Nvglobal, getrank()
13      format(1x,'COORD: more than ',I12,' entries on rank ', I4)
        call xit(1,"Too many vertices")
    endif

    if (nread.gt.meshdatactx%Nvglobal) then
        write(stderr,3) meshdatactx%Nvglobal
 3      format(1x,'COORD: more than ',I12,' entries')
        call xit(1,"Too many vertices")
    endif

    if (n.lt.1 .or. n.gt.meshdatactx%Nvglobal) then
        write(stderr,4) n, meshdatactx%Nvglobal
 4      format(1x,'COORD: entry',I12,' outside range [1,',I12,']')
        call xit(1,"Vertex index incorrect")
    endif

    if (nextread.gt.meshdatactx%nvlocal) then
        write(stderr,*) "thread", getrank() ," read local nodal point nr ", nextread, &
                        "but expected only ", meshdatactx%nvlocal, &
                        "of the global ", meshdatactx%nvglobal
        write(stderr,*) "Perhaps partition.info incorrect?"
        stop "Please correct input files. Leaving GTECTON..."
    endif

    do i=1,NSD
        meshdatactx%X(i,nextread) = f(i)
    enddo

#ifndef	SPARSE
    meshdatactx%globalMarker(nextread) = marker
#endif
    meshdatactx%gloInd(nextread) = n
    meshdatactx%Nvneighborstotal = meshdatactx%Nvneighborstotal + nbrs
    meshdatactx%itot(nextread)   = nbrs

    do i=1,nbrs
        ! this is the global numbering.
        ! in AOmodule, this will be translated to the local numbering.
         if (i.gt.nMaxNeighbors) then
            write(*,*) "A mesh point had too many neighbors."
            write(*,*) "neighbors more than",i-1, "while maximum is", maxNeighbors ,"."
            write(*,*) "This can be stretched, but a mesh with very elongated elements"
            write(*,*) "will give inaccurate results. Please review the mesh."
            stop "Leaving GTecton"
        else
            meshdatactx%AdjM(i,nextread) = nbrlist(i)
        endif
    enddo

! index of 3rd neighbor of tenth vertex of this partition:
! meshdatactx%AdjM(3,10)   

endif

goto 100

200 continue
if (nread.lt.NUMNP) then
    write(stderr,6) nread,NUMNP
 6  format(1x,'COORD: read ',I12,' nodal point coordinates, expected: ',I12)
    call xit(1,"Read less vertices than indicated in TECIN input.")
endif

if (nread.gt.NUMNP) then
    write(stderr,6) nread,NUMNP
    call xit(1,"Read more vertices than indicated in TECIN input.")
endif

if (.not.OUTFIL(1)) then
    if (iecho.eq.8) then
        call ParallelLog("coord", "done")
    endif
    return
endif

nn = 51 
do n=1,meshdatactx%Nvlocal !NUMNP
    nn = nn+1
    if (nn.gt.50) then
        nn = 1
        write(luout(1),7) (LABELC(I),I=1,NSD)
 7           format(///' N O D A L   C O O R D I N A T E   D A T A' &
         ///5X,'NODE NUMBER ',13X,3(A4,18X)/)
        write(luout(1),'(/)')
    endif
    write(luout(1),8) N,(meshdatactx%X(I,N),I=1,NSD)
 8       format(3X,I12,10X,3(1PE20.8))
enddo

return

300     write(stderr,9) ios
 9   format(1x,'COORD: read error ',g14.6)
write(stderr,10) n,(f(i),i=1,NSD)
10   format(1x,'COORD: previous record:',i5,7g14.6)
call prerr('COORD')
call xit(1,"Error when reading vertices")

400     write(stderr,11) ios
11   format(1x,'COORD: premature EOF error ',g14.6)
write(stderr,9) n,(f(i),i=1,NSD)
call prerr('COORD')
call xit(1,"Unexpected EOF when reading vertices")

end
!-------------------------------------------------------------------------------

subroutine EDATA (NDOF)

USE LABELMODULE
USE MATERIALSMODULE
USE MODELDEFINITION, only: icvis, nsed !, radius, ISHELL, numat, igrav
USE MODELCTX
use timestepmodule,  only: NINTG
use debugmodule,     only: iecho, debug
use algebra,         only: idig
!        use constants, only: dfmin,dlogmin,dfmax,dlogmax

!
! Subroutine to read and check data on material properties
!

implicit none
!-pass
integer NDOF
!-local

!       Read and check elasticity and constant viscosity data
call RDEDAT()

!       The three subs are read based on the variable ICVIS
!       which is a binary contraction of the four booleans
!       that appear as checkboxes on the tecinbuilder.
!       if it is 0000 (and appears as a single 0 in the tecinbuilder), all of them return data
!       if it is 1111, all the subs read data.

if (debug) then
    write(*,*) "EDATA says, ICVIS: ", icvis
endif 


!       Read and check powerlaw creep data
if (IDIG(ICVIS,4,1).eq.1) then
        if (iecho.eq.2) then
                write(*,*) "reading material powerlaw creep parameters"
        endif
        call RDCDAT()
endif

#ifdef EXPERIMENTALTRANSIENTDISCREEP
!       Read and check transient dislocation creep data
if (IDIG(ICVIS,4,1).eq.2) then
        if (iecho.eq.2) then
                write(*,*) "reading material transient dislocation creep parameters"
        endif
        call RDDISCDAT()
!       Read and check diffusion creep data
        if (iecho.eq.2) then
                write(*,*) "reading material diffusion creep parameters"
        endif
        call RDDIFFCDAT()
endif
#endif


!       Read and check plasticity data
if (IDIG(ICVIS,4,2).eq.1) then
    if (iecho.eq.2) then
        write(*,*) "reading material plasticity parameters"
    endif
        call RDPDAT()
endif

!       Read and check weakening data, fill PRPDIF
if (IDIG(ICVIS,4,3).eq.1  .or.  &
    IDIG(ICVIS,4,4).eq.1 ) then
    if (iecho.eq.2) then
        write(*,*) "reading material weakening parameters"
    endif
        call RDDDAT()
endif

!do i=1,numat
!    write(*,*) "in edata: PRPDIF", matdatactx%PRPDIF(:,i)
!enddo

!       Read en check on sediment transport parameters
!       only if NSED != 0.
call RDSEDI(PRPSED, NSED, NINTG)

!        if (iecho.eq.2) write(stderr,1)
!    1   format(1x,'if (NINTG != 0 & NSED != 0)',1x,
!     >   'read RHOSED,SEDDIF (2G14.0)')
!        if (NSED.ne.0 .and. NINTG.ne.0) then
!            read(luin,2,err=600,end=600,iostat=ios) (PRPSED(I),I=1,2)
!    2       format(2G14.0)
!            if (OUTFIL(1)) then
!            write(luout(1),3) PRPSED(1),PRPSED(2)
!    3           format(///1X,
!     >           'S E D I M E N T A R Y',3X,'T R A N S P O R T',3X,
!     >           'P A R A M E T E R S'//5X,
!     >           'DENSITY',5X,'EROSION'/17X,'CONSTANT'//5X,
!     >           F7.2,5X,1PE8.1)
!            endif
!        endif
!



!       Read and check on radius; only if ISHELL == 1
call RDRADIUS()
!        if (ISHELL.ne.0) then
!            if (iecho.eq.2) write(stderr,4)
!    4       format(1x,'read RADIUS (G14.0)')
!            read(luin,5,err=600,end=600,iostat=ios) RADIUS
!            read(luin,5,err=600,end=600,iostat=ios) matdatactx%RADIUS
!    5       format(3G14.0)
!            if (RADIUS.le.1.0) then
!                write(stderr,15)
!   15           format(1x, 'ISHELL.ne.0 AND RADIUS.le.1.0')
!                call prerr('EDATA')
!                call xit(1," ")
!            endif
!            if (OUTFIL(1)) then
!                if (ISHELL.ne.0) write(luout(1),7) RADIUS
!    7           format(///1X,'SPHERICAL SHELL RADIUS = ',1PG12.2/)
!            endif
!        endif
!
!       Read and check on gravity

!if (igrav.eq.1) then
        call RDGRAV(GRAV, NDOF)
!else if (igrav.eq.0) then
        ! no gravity, do nothing
!else
!        write(*,*) "sub edata says: Igrav not 0 or 1. This should not happen. "
!    write(*,*) "Please contact model support"
!endif
        

!        if (iecho.eq.2) write(stderr,6)
!    6   format(1x,'read GRAV(i), i=1,NDOF (3G14.0)')
!        read(luin,5,err=600,end=600,iostat=ios) (GRAV(I),I=1,NDOF)
!*
!        if (OUTFIL(1)) then
!            write(luout(1),8)
!    8       format(///1X,
!     >       'G R A V I T Y   L O A D   M U L T I P L I E R S'//)
!            write(luout(1),9) (LABELD(i),GRAV(i),i=1,NDOF)
!    9       format(5X,A4,
!     >       '-DIRECTION  . . . . . . . . . . . . . . . . =',1PG12.2/)
!        endif
!
return


!600     write(stderr,14) ios,luin
!14   format(1x,'EDATA: fatal input error ',i8,' on unit # ',i5)
!call prerr('EDATA')
!call xit(1,"")
end
!-------------------------------------------------------------------------------
subroutine RDSEDI(PRPSED,NSED,NINTG)

use debugmodule, only: iecho, xit, debug
use iomodule          

implicit none


double precision PRPSED(2)
integer NSED, NINTG, i, ios

if (iecho.eq.2) then
    write(stderr,1)
 1   format(1x,'if (NINTG != 0 & NSED != 0)',1x, &
 'read RHOSED,SEDDIF (2G14.0)')
endif

if (NSED.ne.0 .and. NINTG.ne.0) then
    read(luin,2,err=600,end=600,iostat=ios) (PRPSED(I),I=1,2)
 2       format(2G14.0)

    if (debug) then
        write(*,*) "read sediment data", PRPSED(1:2)
    endif



    if (OUTFIL(1)) then
        write(luout(1),3) PRPSED(1),PRPSED(2)
 3           format(///1X, &
         'S E D I M E N T A R Y',3X,'T R A N S P O R T',3X, &
         'P A R A M E T E R S'//5X, &
         'DENSITY',5X,'EROSION'/17X,'CONSTANT'//5X, &
         F7.2,5X,1PE8.1)
    endif
else
    if (debug) then
        write(*,*) "NSED = 0 .and. NINTG = 0; do not read sediment data"
    endif
endif

return

600     write(stderr,14) ios,luin
14   format(1x,'RDSEDI: fatal input error ',i8,' on unit # ',i5)
call prerr('RDSEDI')
call xit(1,"Input error in sediments")

end
!-------------------------------------------------------------------------------
subroutine RDRADIUS()

use modeldefinition, only: radius, ishell
use debugmodule,     only: iecho, xit, debug
use iomodule          

implicit none

integer ios
character(len=1024) record
logical :: formatError
 

if (ISHELL.ne.0) then

    formatError = .false.

    if (iecho.eq.2) then
        write(stderr,4)
 4      format(1x,'read RADIUS (G14.0)')
    endif

    read(luin,"(a)",err=600,end=500,iostat=ios) record
    
    if (len(trim(record)).gt.14) then
        write(*,*) "Input reader says: Something wrong with the radius entry"
        write(*,*) "Line contains too much data", trim(record)
        formatError = .true.
    endif

    if (record(14:14).eq." ") then
        write(*,*) "Input reader says: Something wrong with the radius entry"
        write(*,*) "character 14 is empty"
        formatError = .true.
    endif

    read(record,5,err=600,end=500,iostat=ios) RADIUS
 5       format(1G14.0)

    if (iecho.eq.2  .or. debug .or. formatError) then
        write(*,*) "Read radius", radius
    endif

    if (RADIUS.le.1.0) then
        write(stderr,15)
15           format(1x, 'ISHELL.ne.0 AND RADIUS.le.1.0')
        call prerr('RDRADIUS')
        call xit(1,"Radius should be larger then 1")
    endif

    if (OUTFIL(1)) then
        write(luout(1),7) RADIUS
 7           format(///1X,'SPHERICAL SHELL RADIUS = ',1PG12.2/)
    endif

else
    if (debug) then
        write(*,*) "ishell 0, do not read radius."
    endif
endif


return

500     write(stderr,13) ios,luin
13   format(1x,'RDRADIUS: reached while reading radius ',i8,' on unit # ',i5)
call prerr('RDRADIUS')
call xit(1,"Error when reading radius")

600     write(stderr,14) ios,luin
14   format(1x,'RDRADIUS: fatal input error ',i8,' on unit # ',i5)
call prerr('RDRADIUS')
call xit(1,"Unexpected EOF when reading radius")
end

!-------------------------------------------------------------------------------
subroutine RDGRAV(GRAV,NDOF)

USE LABELMODULE
use debugmodule, only: iecho, parallellog, xit
use iomodule          


implicit none

integer NDOF,i,ios
double precision, dimension(NDOF) :: GRAV

if (iecho.eq.8) then
    call parallelLog("rdgrav","Entered rdgrav")
endif

if (iecho.eq.2) then
    write(stderr,6)
 6   format(1x,'read GRAV(i), i=1,NDOF (3G14.0)')
endif

read(luin,5,err=600,end=600,iostat=ios) (GRAV(I),I=1,NDOF)
 5       format(3G14.0)
!
if (OUTFIL(1)) then
    write(luout(1),8)
 8       format(///1X, &
     'G R A V I T Y   L O A D   M U L T I P L I E R S'//)
    write(luout(1),9) (LABELD(i),GRAV(i),i=1,NDOF)
 9       format(5X,A4, &
     '-DIRECTION  . . . . . . . . . . . . . . . . =',1PG12.2/)
endif

return

600     write(stderr,14) ios,luin
14   format(1x,'RDGRAV: fatal input error ',i8,' on unit # ',i5)
call prerr('RDGRAV')
call xit(1,"Error when reading gravity")

end
!-------------------------------------------------------------------------------
subroutine RDEDAT ()

USE MATERIALSMODULE, only : materialgen, &
                            matdatactx
USE MODELDEFINITION
USE MODELCTX
use constants,       only : dfmin,dlogmin,dfmax,dlogmax
use debugmodule, only: iecho, xit
use timestepmodule,  only : nintg
use iomodule          

!
! reads and checks elasticity and constant viscosity data
!

implicit none

integer :: nn, i
double precision :: E, Pois, thick, dens, emhu, anpwr

!
!       Announce reading operation
if (iecho.eq.2) then
    write(stderr,1)
 1   format(1x,'*** mechanical material properties ***'/ &
 1x,'n=1; while (n != 0 & n != ''end'') {'/ &
 5x,'read (i12,6g14.0)',1x, &
 'n,E(n),POIS(n),EMHU(n),ANPWR(n),WT(n),TH(n)'/ &
 1x,'}')
endif
!
!       Read data
!        call MATGEN (PRPMAT,6,6,NUMAT,0)

!        write(*,*) 'RDEDAT calls materialgen'
call materialgen(matdatactx%PRPMAT,6,6,NUMAT,0)
!
!       Print data
if (OUTFIL(1)) then
    nn = 51
    do i=1,NUMAT
        if (nn.gt.50) then
            write(luout(1),2)
 2               format(//3x,'SET',6x,'YOUNGS',4x,'POISSON',4x, &
             'DENSITY',5x,'THICKNESS',4x,'VISCOSITY',4x,'POWER'/ &
             2x,'NUMBER',4x,'MODULUS',5x,'RATIO'/)
            nn = 0
        endif
        write(luout(1),3) i,matdatactx%PRPMAT(1,i), &
                            matdatactx%PRPMAT(2,i), &
                            matdatactx%PRPMAT(5,i), &
                            matdatactx%PRPMAT(6,i), &
                            matdatactx%PRPMAT(3,i), &
                            matdatactx%PRPMAT(4,i)

 3           format(1x,i5,3x,1Pg10.2,6x,0Pf4.2,2x,1Pg11.3,2x,1Pg10.2, &
         3x,1Pe10.2,5x,0Pf4.2)

        nn = nn + 1
    enddo
endif
!
!       Check data

do  i=1,NUMAT
    E = matdatactx%PRPMAT(1,i)

    if (E.lt.0.0) then
        write(stderr,4) i
 4           format(1x,'Young''s modulus < 0 in material set ',I2)
        call xit(1,"Young''s modulus should be higer than 0")
    endif

    POIS = matdatactx%PRPMAT(2,i)
    if (POIS+1d0.le.dfmin .or. 0.5d0-POIS.le.dfmin) then
        write(stderr,5) i
 5           format(1x,'Poisson''s ratio outside range <-1,0.5> ', &
         'in material set ',I2)
        call xit(1,"Poisson ration outside range <-1,0.5>")
    endif

    THICK = matdatactx%PRPMAT(6,i)
    if (NTYPE.lt.5 .and. THICK.le.dfmin) then
        write(stderr,6) THICK, i
 6           format(1x,"Element thickness", G14.6, "<= 0 in material set ",I2)
        call xit(1,"Element thickness <= 0 in material set")
    endif

    DENS = matdatactx%PRPMAT(5,i)
    if (IGRAV.ne.0 .and. DENS.le.0.0) then
        write(stderr,7) i
 7           format(1x,'Warning: mass density <= 0 in',1x, &
         'material set ',I2)
    endif

    if (MODE.ge.3 .and. DENS.le.0.0) then
        write(stderr,8)
 8           format(10x,'THERMAL DIFFUSION WILL NOT WORK'/)
    endif

    if (NINTG.gt.0) then

        EMHU = matdatactx%PRPMAT(3,i)
        if (EMHU.le.dfmin) then
            write(stderr,9) i
 9                format(1x,'Viscosity <= 0 in material set ',I2)
            call xit(1,"Viscosity <= 0 in material set")
        endif

        ANPWR = matdatactx%PRPMAT(4,i)
        if (ANPWR.le.dfmin) then
            write(stderr,10) i
10           format(1x,'Viscosity power <= 0 in material set ',I2)
            call xit(1,"Viscosity power <= 0 in material set")
        endif

    endif
enddo

return
end
!-------------------------------------------------------------------------------
integer function IDIG (number,length,i)

implicit none

integer   :: length
integer   :: number, i

integer   :: n
integer   :: j, jdig

!
! returns digit number I of integer NUMBER with LENGTH positions
!
n = IABS(number)
if (length.gt.i) then
    do j=1,length-i
        jdig = MOD(n,10)
        n = (n-jdig)/10
    enddo
    IDIG = MOD(n,10)
else if (length.eq.i) then
    IDIG = MOD(n,10)
else
    IDIG = 0
endif

return
end
!-------------------------------------------------------------------------------
subroutine RDCDAT ()

USE MATERIALSMODULE,  only: materialgen, &
                            matdatactx
USE MODELCTX
USE TIMESTEPMODULE,   only: NINTG
USE MODELDEFINITION,
use constants,        only: dfmin,dlogmin,dfmax,dlogmax,TK0,R
use debugmodule,      only: iecho, xit, debug
use algebra,          only: idig
use iomodule          

!
!       Reads and checks powerlaw creep data for variable viscosity
!

implicit none

integer :: i, n, nn
double precision :: qpl, qmax, apl, anpl, anplmax, v

!
!       Announce reading operation
if (iecho.eq.2) then
    write(stderr,1)
 1   format(1x,'if (NINTG > 0 & ICVIS) {'/ &
 5x,'*** powerlaw creep properties ***'/ &
 5x,'n=1; while (n != 0 & n != ''end'') {'/ &
 9x,'read (i12,5g14.0)',1x, &
 'n,QPL(n),APL(n),NPL(n),VISMIN(n),VISMAX(n)'/ &
 5x,'}'/1x,'}')
endif

if (IDIG(ICVIS,4,1).eq.0 .or. NINTG.eq.0) then
    return
endif

if (iecho.eq.2 .or. debug) then
    write(*,*) "Reading powerlaw creep data"
endif

!       Read data
call materialgen (matdatactx%PRPCRP,6,6,NUMAT,0)
!        call MATGEN (PRPCRP,6,6,NUMAT,0)
do i=1,NUMAT
    if (matdatactx%PRPCRP(4,i).lt.dfmin) then
        matdatactx%PRPCRP(4,i) = dfmin
    endif
    if (matdatactx%PRPCRP(5,i).lt.dfmin) then
        matdatactx%PRPCRP(5,i) = dfmax
    endif
enddo
!
!       Print data
if (OUTFIL(1)) then
    nn = 51
    do i=1,NUMAT
        if (nn.gt.50) then
            write(luout(1),2)
 2               format(//1X,'P O W E R L A W',3X,'C R E E P',3X, &
             'P A R A M E T E R S'// &
             3X,'SET',5X,'ACTIVATION',2X,'PRE-EXPONENT', &
             2X,'POWER',9X,'VISCOSITY'/2X,'NUMBER', &
             5X,'ENERGY',6X,'CONSTANT',14X,'MINIMUM   MAXIMUM'/)
            nn = 1
        endif
        write(luout(1),3) i,(matdatactx%PRPCRP(n,i),n=1,6)
 3                format(1x,i5,5x,1Pg9.1,3x,1Pg10.2,5x,0PF4.2,4x, &
         1PE9.2,1x,1PE9.2,1x,1PE9.2)
        nn = nn + 1
    enddo
endif

!       Check data
do i=1,NUMAT
    QPL = matdatactx%PRPCRP(1,i)
    QMAX = R * (TK0+1200) * LOG(dfmax)

    if (QPL.lt.0. .or. QPL.ge.QMAX) then
        write(stderr,5) QMAX,i
 5           format(1x,'Powerlaw activation energy out of range [0,', &
         1PE9.1,'> in material set ',I2)
        call xit(1,"Powerlaw activation energy out of range")
    endif

    APL = matdatactx%PRPCRP(2,i)
    if (APL.le.dfmin .or. APL.ge.dfmax) then
        write(stderr,6) dfmax,i
 6           format(1x,'Powerlaw pre-exponent out of range <0,', &
         1PE9.1,'> in material set ',I2)
        call xit(1,"Powerlaw pre-exponent out of range")
    endif

    ANPL = matdatactx%PRPCRP(3,i)
    ANPLMAX = dlogmax/6d0
    if (ANPL.lt.1.0 .or. ANPL.ge.ANPLMAX) then
        write(stderr,*) "Material", i, "has powerlaw power", ANPL
           write(stderr,*)    "power should be in range [1, ", ANPLMAX, ">"
!        write(stderr,7) INT(ANPLMAX),i
! 7           format(1x,'Powerlaw power out of range [1,', &
!         I2,'> in material set ',I2)
        call xit(1,"Powerlaw power out of range")
    endif

    V = matdatactx%PRPCRP(6,i) ! activation volume
    if (V.lt.0d0 .or. V.gt.1.0D-3) then
        write(stderr,*) "Material", i, "has activation volume", V
        write(stderr,*) "Activation volume should be in range [0, 0.001>"

!        write(stderr,8) 1.0D-3,i
! 8           format(1x,'Powerlaw activation volume out of range [0,', &
!         1PE9.1,'> in material set ',I2)
        call xit(1,"Powerlaw activation volume out of range")
    endif

enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine RDDDAT ()

USE MATERIALSMODULE, only : materialgen, &
                            matdatactx
USE MODELCTX
!USE TIMESTEPMODULE,  only : NINTG
USE MODELDEFINITION
use constants,       only: dfmin,dlogmin,dfmax,dlogmax
use debugmodule,     only: iecho, xit, debug
use algebra,         only: idig
use iomodule          
!
!       reads and checks strain (rate) weakening data
!

implicit none
!-pass
!-locl
character(len=11) vtype,ctype
integer :: i, n, nn
!
!       Announce reading operation
if (iecho.eq.2) then
    write(stderr,1)
 1       format(1x,'if (NINTG>0 & ICVIS) {'/ &
     5x,'*** strain (rate) weakening properties ***'/ &
     5x,'n=1; while (n != 0 & n != ''end'') {'/ &
     9x,'read (i12,4g14.0)',1x, &
     'n,VWEAK(n),VSIGMA(n),CWEAK(n),CSIGMA(n)'/ &
     5x,'}'/1x,'}')
endif
!
!if ((IDIG(ICVIS,4,3).eq.0  .and. & 
!     IDIG(ICVIS,4,4).eq.0) .or.  &
!               NINTG.eq.0)         then
!    return
!endif

if (IDIG(ICVIS,4,3).eq.0  .and. &
    IDIG(ICVIS,4,4).eq.0) then
    return
endif


if (iecho.eq.2 .or. debug) then
    write(*,*) "Reading strain (rate) weakening data"
endif

!       Read creep weakening data
!        call MATGEN (PRPDIF,9,7,NUMAT,0)
!        write(*,*) 'RDDDAT calls materialgen'
call materialgen (matdatactx%PRPDIF,9,7,NUMAT,0)

vtype = 'STRAIN     '
ctype = 'STRAIN     '

if (IDIG(ICVIS,4,3).eq.0) then
    do i=1,NUMAT
        matdatactx%PRPDIF(1,i) = 1d0
        matdatactx%PRPDIF(2,i) = 0d0
    enddo
else if (IDIG(ICVIS,4,3).eq.2) then
    vtype = 'STRAIN RATE'
endif

if (IDIG(ICVIS,4,4).eq.0) then
    do i=1,NUMAT
        matdatactx%PRPDIF(3,i) = 1d0
        matdatactx%PRPDIF(4,i) = 0d0
    enddo
else if (IDIG(ICVIS,4,4).eq.2) then
    ctype = 'STRAIN RATE'
endif

!       Print data

if (OUTFIL(1)) then
    nn = 51
    do i=1,NUMAT
        if (nn.gt.50) then
            write(luout(1),2)
 2               format(//1x,'S T R A I N',3X,'( R A T E )',3X, &
             'W E A K E N I N G   P A R A M E T E R S'/)
            write(luout(1),3) vtype,ctype
 3               format(3X,'SET',4X,'VISCOSITY ',A11,8X, &
            'YIELD STRESS ', &
             A11/10X,'WEAKENING',4X,'INTERVAL',8X,'WEAKENING', &
             4X,'INTERVAL'/)
            nn = 1
        endif
        write(luout(1),4) i,(matdatactx%PRPDIF(n,i),n=1,4)
 4           format(I8,2X,1PE9.2,3X,1PE9.2,8X,1PE9.2,3X,1PE9.2)
        nn = nn + 1
    enddo
endif

!       Check data

do i=1,NUMAT

    if (IDIG(ICVIS,4,3).ne.0) then
        if (matdatactx%PRPDIF(1,i).le.-dfmin) then
            write(stderr,10) i
10               format(1x,'Negative viscous weakening parameter', &
             ' in material data set ',I2)
            call xit(1,"Negative viscous weakening parameter")
        endif
        if (matdatactx%PRPDIF(1,i).ne.1d0 .and. &
            matdatactx%PRPDIF(2,i).le.dfmin) then
            write(stderr,11) i
11               format(1x,'Viscous weakening interval too small', &
             ' in material data set ',I2)
            call xit(1,"Viscous weakening interval too small")
        endif
    endif

    if (IDIG(ICVIS,4,4).ne.0) then
        if (matdatactx%PRPDIF(3,i).le.-dfmin) then
            write(stderr,12) i
12               format(1x,'Negative yield stress weakening', &
             ' parameter in material data set ',I2)
            call xit(1,"Negative yield stress weakening")
        endif
        if (matdatactx%PRPDIF(3,i).ne.1d0 .and. &
            matdatactx%PRPDIF(4,i).le.dfmin) then
            write(stderr,13) i
13               format(1x,'Yield stress weakening interval too ', &
             ' small in material data set ',I2)
            call xit(1,"Yield stress weakening interval too small")
        endif
    endif

enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine RDPDAT ()

USE MATERIALSMODULE, only : materialgen, &
                            matdatactx
USE MODELCTX
!USE TIMESTEPMODULE, only : NINTG
USE MODELDEFINITION
use constants, only: dfmin,dlogmin,dfmax,dlogmax
use debugmodule, only: iecho, xit, debug
use algebra,         only: idig
use iomodule          

!
! reads and checks plasticity data
!

implicit none
!-locl
character(len=14) PLTYP
character(len=5) law
dimension PLTYP(0:4)
external FLOAT
integer :: i, n, nn
double precision :: PSTRS, fluidy, YLDPWR, frangl
integer          :: NPTYPE
!-init
data PLTYP /'NO PLASTICITY ','    TRESCA    ','  VON MISES   ', &
 ' MOHR-COULOMB ','DRUCKER-PRAGER'/


!       Announce reading operation
if (iecho.eq.2) then
    write(stderr,1)
 1       format(1x,'if (NINTG>0 & ICVIS) {'/ &
     5x,'*** mechanical plastic properties ***'/ &
     5x,'n=1; while (n != 0 & n != ''end'') {'/ &
     9x,'read (4i5,5g14.0)',1x, &
     'n,NPTYPE(n),NFLOW(n),FRANGL(n),YIELD(n),FLUIDY(n),', &
     'YLDPWR(n),STNHRD(n)'/5x,'}'/1x,'}')
endif
!
!if (IDIG(ICVIS,4,2).eq.0 .or. NINTG.eq.0) then
!    return
!endif

if (IDIG(ICVIS,4,2).eq.0) then
    return
endif

if (iecho.eq.2 .or. debug) then
    write(*,*) "Reading plasticity data"
endif

!
!       Read data

call materialgen (matdatactx%PRPLAS,9,7,NUMAT,1)
!
if (iecho.eq.2 .or. debug) then
    write(*,*) "Finished reading plasticity data"
endif



if (OUTFIL(1)) then
    nn = 51
    do i=1,NUMAT

        if (nn.gt.50) then
            write(luout(1),2)
 2               format(// &
!                    10x,'PLASTICITY',3x,'FRICTION',2x,'YIELD', &
!                    3x,'FLUIDITY',3x,'YIELD',3x,'YIELD',5x,'STRAIN'/ &
!                    13x,'TYPE',8x,'ANGLE',4x,'YIELD',16X,'POWER',4x, &
!                    'LAW',4x,'HARDENING'/34X,'STRESS'/) &
             10x,'PLASTICITY',3x,'FRICTION',3x,'YIELD', &
             4x,'FLUIDITY',4x,'YIELD',4x,'YIELD',5x,'STRAIN'/ &
             13x,'TYPE',8x,'ANGLE',4x,'STRESS',15X,'POWER',4x, &
             'LAW',4x,'HARDENING'/)
            nn = 1
        endif

        NPTYPE = INT(matdatactx%PRPLAS(1,i))

        if (NPTYPE.lt.0 .or. NPTYPE.gt.4) then
            write(stderr,3) i
 3               format(1x,'Plasticity type out of range [0,4]', &
             ' in material set ',I2)
            call xit(1,"Plasticity type out of range [0,4]")
        endif

        if (INT(matdatactx%PRPLAS(6,i)).eq.0) then
            law = ' EXP '
        else
            law = 'POWER'
        endif

        write(luout(1),4) i,PLTYP(NPTYPE),(matdatactx%PRPLAS(n,i),n=2,5), &
         law,matdatactx%PRPLAS(7,i)
 4                format(1x,i5,2x,a14,4x,f4.1,3x,1PE8.2,3x,1PE8.2,4x, &
         0PF4.1,3x,a5,4x,1PE8.2)
        nn = nn + 1
    enddo
endif

!       Check data

do i=1,NUMAT

    NPTYPE = INT(matdatactx%PRPLAS(1,i))

    if (NPTYPE.ne.0) then

        PSTRS = matdatactx%PRPLAS(3,i)

        if (PSTRS.le.0.) then
            write(stderr,5) i
 5           format(1x,'Uniaxial yield stress <= 0 in', &
         ' material set',I2)
            call xit(1,"Uniaxial yield stress <= 0 in")
        endif

        FLUIDY = matdatactx%PRPLAS(4,i)
        if (FLUIDY.le.dfmin .or. FLUIDY.ge.dfmax) then
            write(stderr,6) dfmax,i
 6           format(1x,'Fluidity out of range <0,',1PE9.1,'> ', &
         'in material set ',I2)
            call xit(1,"Fluidity out of range")
        endif

        YLDPWR = matdatactx%PRPLAS(5,i)
        if (YLDPWR.lt.0d0 .or. YLDPWR.ge.dlogmax) then
            write(stderr,7) INT(dlogmax),i
 7           format(1x,'Plastic yield power out of range [0,', &
         I3,'> in material set ',I2)
            call xit(1,"Plastic yield power out of range")
        endif

        if (NPTYPE.eq.3 .or. NPTYPE.eq.4) then
            FRANGL = matdatactx%PRPLAS(2,i)
            if (ABS(FRANGL).gt.90d0) then
                write(stderr,8) i
 8               format(1x,'Friction angle out of range ', &
             '[-90,90] in material set ',I2)
!                    call xit(1," ")
            endif
        endif

        if (NPTYPE.eq.1 .or. NPTYPE.eq.3) then
            write(stderr,9)
 9           format(1x,'Tresca and Mohr-Coulomb criteria ', &
             'incompletely implemented')
            call xit(1,"Tresca and Mohr-Coulomb criteria incompletely implemented")
        endif

    endif
enddo
!
!       Initialize
do i=1,NUMAT
    matdatactx%PRPLAS(8,i) = matdatactx%PRPLAS(3,i)
    matdatactx%PRPLAS(9,i) = matdatactx%PRPLAS(4,i)
enddo
!
return
end

!-------------------------------------------------------------------------------
#ifdef EXPERIMENTALTRANSIENTDISCREEP

subroutine RDDISCDAT()

USE MATERIALSMODULE,  only: materialgen, &
                            matdatactx
USE MODELCTX
USE TIMESTEPMODULE,   only: NINTG
USE MODELDEFINITION,
use constants,        only: dfmin,dlogmin,dfmax,dlogmax,TK0,R
use debugmodule,      only: iecho, xit, debug
use algebra,          only: idig
use iomodule          

!
!       Reads and checks dislocation creep data for stress/state variable viscosity
! State is dependent on internal stress (i.e. Taylor stress) related to variable
! dislocation density
! based on Breithaupt et al. 2023 https://doi.org/10.1073/pnas.2203448120

implicit none

integer :: i, n, nn
!double precision :: qpl, qmax, apl, anpl, anplmax, v
double precision :: qpl,qmax,apl,RatioElast
!
!       Announce reading operation
if (iecho.eq.2) then
    write(stderr,1)
 1   format(1x,'if (NINTG > 0 & ICVIS) {'/ &
 5x,'*** transient/steady state dislocation creep properties ***'/ &
 5x,'n=1; while (n != 0 & n != ''end'') {'/ &
 9x,'read (i12,10g14.0)',1x, &
 'n,QDC(n),ADC(n),SigmaP(n),b(n),BetaThres(n),MTaylor(n),R0gb(n),R0pipe(n),VDC(n),RatioRhoInit(n)'/ &
 5x,'}'/1x,'}')
endif

if (IDIG(ICVIS,4,1).ne.2 .or. NINTG.eq.0) then
    return
endif

if (iecho.eq.2 .or. debug) then
    write(*,*) "Reading transient-steadystate dislocation creep data"
endif

!       Read data

call materialgen (matdatactx%PRPDISLC,10,10,NUMAT,0)

! for now, do not apply checks to input
!do i=1,NUMAT
!    if (matdatactx%PRPCRP(4,i).lt.dfmin) then
!        matdatactx%PRPCRP(4,i) = dfmin
!    endif
!    if (matdatactx%PRPCRP(5,i).lt.dfmin) then
!        matdatactx%PRPCRP(5,i) = dfmax
!    endif
!enddo
!
!       Print data
if (OUTFIL(1)) then
    nn = 51
    do i=1,NUMAT
        if (nn.gt.50) then
            write(luout(1),2)
 2               format(//1X,'T R A N S I E N T - S T E A D Y S T A T E D I S L O C A T I O N',3X,'C R E E P',3X, &
             'P A R A M E T E R S'// &
             3X,'SET',4X,'ACTIVATION',1X,'PRE-EXPONENT', &
             2X,'PEIERLS',4X,'BURGERS',14X,'BETA', &
             4X,'M, TAYLOR',5X,'RATE COEFF.', &
             7X,'RATE COEFF.',&  
             2x,'ACTIVATION' , &
             3x,'INITIAL RATIO'/,&  
             
             1X,'NUMBER',7X,'ENERGY',2X,'COEFFICIENT', &
             3X,'STRESS',5X,'LENGTH',2X,'THRESHOLD STRESS', &
             2X,'STRESS RATE',2X,'GRAIN BOUNDARY', &
             4X,'PIPE DIFFUSION',&  
             6X,'VOLUME' , &
             3X,'TAYLOR STRESS TO'/, &
             
             
             60X,'COEFFICIENT', &
             2X,'COEFFICIENT',7X,'DIFFUSION', &
             33X,'INITIAL STRESS'/)
            nn = 1
        endif
        write(luout(1),3) i,(matdatactx%PRPDISLC(n,i),n=1,10)
 3                format(1x,i5,5x,1Pg9.2,2x,1Pg9.2,2x,1Pg9.3,2x, &
         1Pg9.3,2x,1Pg9.3,11x,1Pg9.3,7x,1Pg9.3,9x,1Pg9.3,7x,&
         1Pg9.3,7x,1Pg9.3)
        nn = nn + 1
    enddo
endif

!       Check data, still incomplete
do i=1,NUMAT
    QPL = matdatactx%PRPDISLC(1,i)
    QMAX = R * (TK0+1200) * LOG(dfmax)

    if (QPL.lt.0. .or. QPL.ge.QMAX) then
        write(stderr,5) QMAX,i
 5           format(1x,'Dislocation creep activation energy out of range [0,', &
         1PE9.1,'> in material set ',I2)
        call xit(1,"Dislocation creep activation energy out of range")
    endif

    APL = matdatactx%PRPDISLC(2,i)
    if (APL.le.0 .or. APL.ge.dfmax) then
        write(stderr,6) dfmax,i
 6           format(1x,'Dislocation creep pre-exponent out of range <0,', &
         1PE9.1,'> in material set ',I2)
        call xit(1,"Dislocation creep pre-exponent out of range")
    endif
    
    RatioElast = matdatactx%PRPDISLC(10,i)

    if (RatioElast.le.0 .or. RatioElast.gt.1) then
    write(stderr,7) 1d0,i
 7           format(1x,'Ratio initial Taylor stress out of range <0,', &
         1PE9.1,'> in material set ',I2)
        call xit(1,"Ratio initial Taylor stress out of range")
    endif
    
! add checks for other parameters


!    ANPL = matdatactx%PRPCRP(3,i)
!    ANPLMAX = dlogmax/6d0
!    if (ANPL.lt.1.0 .or. ANPL.ge.ANPLMAX) then
!        write(stderr,*) "Material", i, "has powerlaw power", ANPL
!        write(stderr,*) "power should be in range [1, ", ANPLMAX, ">"
!        write(stderr,7) INT(ANPLMAX),i
! 7           format(1x,'Powerlaw power out of range [1,', &
!         I2,'> in material set ',I2)
!        call xit(1,"Powerlaw power out of range")
!    endif
!
!    V = matdatactx%PRPCRP(6,i) ! activation volume
!    if (V.lt.0d0 .or. V.gt.1.0D-3) then
!        write(stderr,*) "Material", i, "has activation volume", V
!        write(stderr,*) "Activation volume should be in range [0, 0.001>"

!        write(stderr,8) 1.0D-3,i
! 8           format(1x,'Powerlaw activation volume out of range [0,', &
!         1PE9.1,'> in material set ',I2)
!        call xit(1,"Powerlaw activation volume out of range")
!    endif

enddo
!
return
end

#endif
!-------------------------------------------------------------------------------

#ifdef EXPERIMENTALTRANSIENTDISCREEP
subroutine RDDIFFCDAT()

USE MATERIALSMODULE,  only: materialgen, &
                            matdatactx
USE MODELCTX
USE TIMESTEPMODULE,   only: NINTG
USE MODELDEFINITION,
use constants,        only: dfmin,dlogmin,dfmax,dlogmax,TK0,R
use debugmodule,      only: iecho, xit, debug
use algebra,          only: idig
use iomodule          

!
!       Reads and checks diffusion creep parameters

implicit none
!-local
integer :: i, n, nn
double precision :: qdiffc,adiffc,pExp,d,vdiffc,qmax

!
!       Announce reading operation
if (iecho.eq.2) then
    write(stderr,1)
 1   format(1x,'if (NINTG > 0 & ICVIS) {'/ &
 5x,'*** diffusion creep properties ***'/ &
 5x,'n=1; while (n != 0 & n != ''end'') {'/ &
 9x,'read (i12,5g14.0)',1x, &
 'n,Qdiffc(n),Adiffc(n),PExp(n),Vdiffc(n),d(n)'/ &
 5x,'}'/1x,'}')
endif

if (IDIG(ICVIS,4,1).ne.2 .or. NINTG.eq.0) then
    return
endif

if (iecho.eq.2 .or. debug) then
    write(*,*) "Reading diffusion creep data"
endif

!       Read data
call materialgen (matdatactx%PRPDIFFC,5,5,NUMAT,0)

! for now, do not apply checks to input
!do i=1,NUMAT
!    if (matdatactx%PRPCRP(4,i).lt.dfmin) then
!        matdatactx%PRPCRP(4,i) = dfmin
!    endif
!    if (matdatactx%PRPCRP(5,i).lt.dfmin) then
!        matdatactx%PRPCRP(5,i) = dfmax
!    endif
!enddo
!
!       Print data
if (OUTFIL(1)) then
    nn = 51
    do i=1,NUMAT
        if (nn.gt.50) then
            write(luout(1),2)
 2               format(//1X,'DIFFUSION',3X,'CREEP',3X, &
             'P A R A M E T E R S'// &
             3X,'SET',4X,'ACTIVATION',2X,'PRE-EXPONENT', &
             2X,'GRAIN SIZE',2X,'ACTIVATION',6X,'GRAIN' /, &
             6X,8X,'ENERGY',14X, &
             4X,'EXPONENT',6X,'VOLUME',7X,'SIZE' /)
            nn = 1
        endif
        write(luout(1),3) i,(matdatactx%PRPDIFFC(n,i),n=1,5)
 3                format(1x,i5,5x,1Pg9.2,5x,1Pg9.2,4x,1Pg9.2,2x,1Pg9.2,2x,1Pg9.2)
        nn = nn + 1
    enddo
endif

!       Check data, still incomplete
do i=1,NUMAT
    Qdiffc = matdatactx%PRPDIFFC(1,i)
    QMAX = R * (TK0+1200) * LOG(dfmax)

    if (Qdiffc.lt.0. .or. Qdiffc.ge.QMAX) then
        write(stderr,5) QMAX,i
 5           format(1x,'Dislocation creep activation energy out of range [0,', &
         1PE9.1,'> in material set ',I2)
        call xit(1,"Diffusion creep activation energy out of range")
    endif

    ADiffc = matdatactx%PRPDIFFC(2,i)
    if (ADiffc.le.dfmin .or. ADiffc.ge.dfmax) then
        write(stderr,6) dfmax,i
 6           format(1x,'Diffusion creep pre-exponent out of range <0,', &
         1PE9.1,'> in material set ',I2)
        call xit(1,"Dislocation creep pre-exponent out of range")
    endif
! add checks for other parameters

!    ANPL = matdatactx%PRPCRP(3,i)
!    ANPLMAX = dlogmax/6d0
!    if (ANPL.lt.1.0 .or. ANPL.ge.ANPLMAX) then
!        write(stderr,*) "Material", i, "has powerlaw power", ANPL
!        write(stderr,*) "power should be in range [1, ", ANPLMAX, ">"
!        write(stderr,7) INT(ANPLMAX),i
! 7           format(1x,'Powerlaw power out of range [1,', &
!         I2,'> in material set ',I2)
!        call xit(1,"Powerlaw power out of range")
!    endif
!
!    V = matdatactx%PRPCRP(6,i) ! activation volume
!    if (V.lt.0d0 .or. V.gt.1.0D-3) then
!        write(stderr,*) "Material", i, "has activation volume", V
!        write(stderr,*) "Activation volume should be in range [0, 0.001>"

!        write(stderr,8) 1.0D-3,i
! 8           format(1x,'Powerlaw activation volume out of range [0,', &
!         1PE9.1,'> in material set ',I2)
!        call xit(1,"Powerlaw activation volume out of range")
!    endif

enddo
!
return
end
#endif

!-----------------------------------------------------------------------
subroutine ELGEN (NEN)

!
! Reads and generates element nodal data and material set numbers
!

USE MODELCTX
USE MODELDEFINITION
USE MESHDATAMODULE
USE FORMATSMODULE, only: FMT_tecin_partf_elm
use debugmodule,   only: iecho, debug, &
                         parallellog, xit
use iomodule          

implicit none
!-pass
integer NEN
!-locl
logical end
character(len=86) record  ! LWITODO increase length
character(len=13) fmt
integer :: ia(8)
integer :: ID, material
integer tmp, nread, nextread
integer ( kind = 4 ) h,m,s,values(8) ! for timestamps in output
integer :: l, lg, i, ios, n
integer :: mat1
integer :: nPartitions

!-init
L     = 0
LG    = 0
nread = 0
nextread = 0
fmt = '(I5,I12,9I12)'
nPartitions = getsize()

if (iecho.eq.8) then
    call ParallelLog("ELGEN", "Reading elements")
endif

if (iecho.eq.2) then
    write(stderr,1) fmt
 1   format(1x,'*** reading element definitions ***'/ &
 1x,'elno=1; while (elno != 0) {'/ &
 5x,'read ',A11,' rank,elno,MAT(elno),(IEN(i,elno),i=1,NEN)'/ & !,increm'/ &
 1x,'}')
endif

end = .false.
loop_read_element_input: do while (.not.end) 

    call readrc (luin,record,'ELGEN',end)
    if (debug) then
        write(*,*) "ELGEN read element record: ", trim(record)
    endif
    if (end) then
        if (debug) then
            write(*,*) "detected an ''end'' statement"
        endif
        exit loop_read_element_input
    endif

    nread = nread + 1
    if (nread.gt.meshdatactx%Neglobal) then
        write(stderr,3) meshdatactx%Neglobal
 3      format(1x,'ELGEN: more than the expected ',I12,' element ', &
           'entries in input list')
        call xit(1,"Too many elements")
    endif

    if (iecho.eq.2 .and. mod(nread,99999) .eq. nread/99999) then
        call date_and_time ( values = values )
        h = values(5)
        m = values(6)
        s = values(7)
        write(*,10) nread, meshdatactx%neglobal, h,m,s
10      format(1x,'progress: ',i12,' of ',i12, ' on time ', i2, ':',i2,':',i2)
    endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc 
    read(record,'(I5,6I13)',err=300,end=400,iostat=ios) tmp, ID, material, &
        meshdatactx%globalIen(1:4, nread)
#else
    read(record,'(I5,6I13)',err=300,end=400,iostat=ios) tmp, ID, material 
#endif

#ifdef SPARSE
if (tmp .ge. nPartitions) then
    if (getrank() .eq. 0) then
    write(stderr,*) "The simulation is run over", nPartitions, "Partitions."
    write(stderr,*) "Numbering of partitions starts at 0."
    write(stderr,*) "Example: Four partitions are numbered 0, 1, 2 and 3."
    write(stderr,*) "Maximum partition number for this model is", nPartitions-1
    write(stderr,*) "but element ", nread, "has partition", tmp, "which is too high"
    STOP "Leaving GTECTON..."
    endif
endif
#endif

    meshdatactx%E2P(nread) = tmp    ! element partition number


#ifndef SPARSE
    tmp = getrank()
#endif

    if (tmp.eq.getrank()) then

        read(record,*,err=301,end=401,iostat=ios) tmp, n, mat1, (ia(i), i=1,NEN)

        if (n.eq.0) then        ! element number=0 also signals end of element list
            exit loop_read_element_input
        endif

        nextread = nextread + 1
        if (nextread.gt.meshdatactx%Nelocal) then
            write(stderr,20) meshdatactx%Nelocal, getrank()
20          format(1x,'ELGEN: more than the expected ',I12, &
                ' element entries in input list for rank ',I12)
            call xit(1,"Too many elements")
        endif

        if (mat1.lt.1 .or. mat1.gt.NUMAT) then
            write(stderr,30) mat1, n,NUMAT
30          format(1x,'ELGEN: material entry ',I12,' for element ',I12,' refers to ', &
               'material property group outside the range [1,',I12,']')
            call xit(1,"entry for material outside range")
        endif
 
        meshdatactx%mat(nextread) = mat1
        meshdatactx%gloElt(nextread) = n

        do i=1,NEN
            if (ia(i).lt.1 .or. ia(i).gt.meshdatactx%Nvglobal) then
                write(stderr,40) n,meshdatactx%Nvglobal
40              format(1x,'ELGEN: entry for element ',I12,' refers to ', &
                    'node numbers outside the range [1,',I12,']')
                call xit(1,"Element index out of range")
            endif
            meshdatactx%ien(i,nextread) = ia(i)
        enddo

        if (debug) then
            write(*,*) 'IEN after reading',n, meshdatactx%ien(1:4,nextread)
        endif

    else
        ! element does not belong to this partition.
        ! check whether it is a ghost element, meaning it is 
        ! in another partition, but borders on an element in this partition.
        ! It does not even have to border on a node in this partition. 
        ! nor does it have to have a boundary in common with an element in 
        ! this partition.
    
        ! imagine:
        ! 
        !  _______a_______a________
        !  \      /\      /\      /
        !   \    /  \ e  /  \    /      <-- elements of partition A
        !    \  /    \  /    \  /
        !     \/______\/______\/      . . . . . . . . . partition boundary  
        !     a\     a/\     a/
        !       \    /  \    /         <-- elements of partition B
        !        \  /    \  /
        !         \/______\/
        !
        !
        ! points a belongs to partition A, element e belongs to partition A,
        ! but is is ghost element of partition B.
        ! How can B detect this element?
        ! We loop over the elements in partition B,
        ! And over the nodes in this element.
        ! If the element in question has even a single point in 
        ! common, in global numbering, with an element in B, it is 
        ! a ghost element of B. It is irrelevant then to which partition the 
        ! point belongs.
    
        ! This is a quadratic search, which is rather unpleasant.
        ! to speed things up, we must illiminate our search quite 
        ! drastically. The only way we can do this, is so keep track of
        ! whether some elements have all of their surrounding elements
        ! (not just the neighboring, but also the poitns sharing
        ! elements) known. We can than flag that element off the search list.
        
        ! Oh, shit, does not work, because all elements of our own 
        ! partition are not yet known.
    
        ! must first read all the elements. and let each partition fill.
        ! Then we must read again to determine the ghost elements.
        ! That or a global array. 
    endif

enddo loop_read_element_input


if (nread.ne.meshdatactx%Neglobal) then
    write(stderr,50) nread,meshdatactx%Neglobal
50  format(1x,'ELGEN: read ',I12,' element entries, expected ',I12)
    call xit(1,"")
endif

if (nextread.ne.meshdatactx%Nelocal) then
    write(stderr,60) getrank(),meshdatactx%Nelocal,nextread
60  format(1x,'ELGEN: rank ',I3,' expected ',I12,' element entries, but read ',I12)
    call xit(1,"")
endif

if (debug) then
    write(*,*) "rank=",getrank(),"processed all elements, assimilated ",meshdatactx%Nelocal, &
     " elements for this partition"
endif

if (iecho.eq.8) then
    call ParallelLog("ELGEN", "Finished reading elements")
endif
 
return
 


300 write(stderr,*) "Error when reading first bit of element line"
    write(stderr,*) "rank", getrank(), "tmp, ID, material", tmp, ID, material
    write(stderr,*) "from record", trim(record)
    write(stderr,*) "Error: ", ios
call xit(1," ")

400 write(stderr,*) "Unexpected end of file when reading first bit of element line"
    write(stderr,*) "rank", getrank(), "tmp, ID, material", tmp, ID, material
    write(stderr,*) "from record", trim(record)
    write(stderr,*) "Error: ", ios
call xit(1," ")

301 write(stderr,*) "Error when reading full element line"
    write(stderr,*) "rank", getrank(), "tmp, ID, material", tmp, ID, material, ia
    write(stderr,*) "from record", trim(record)
    write(stderr,*) "Error: ", ios
call xit(1," ")

401 write(stderr,*) "Unexpected end of file when reading full element line"
    write(stderr,*) "rank", getrank(), "tmp, ID, material", tmp, ID, material, ia
    write(stderr,*) "from record", trim(record)
    write(stderr,*) "Error: ", ios
call xit(1," ")

end subroutine
!-------------------------------------------------------------------------------
subroutine ELMAT ()

USE MESHDATAMODULE,  ONLY: meshdatactx
USE MATERIALSMODULE, ONLY: matdatactx, LMAT
USE MODELDEFINITION
USE MODELCTX
use algebra,         only: idig
!
! Routine to assign material numbers to elements
!

implicit none

integer, parameter :: NMAT=6,NTEM=6,NCRP=6,NDIF=9,NPLS=9
#ifdef EXPERIMENTALTRANSIENTDISCREEP
integer, parameter :: Ndiffc=5,Ndislc=10
#endif
!-locl

integer :: i,j, imat

double precision, allocatable, dimension(:,:)  :: tmpmat
double precision, allocatable, dimension(:,:)  :: tmpcrp
double precision, allocatable, dimension(:,:)  :: tmptem
double precision, allocatable, dimension(:,:)  :: tmplas
double precision, allocatable, dimension(:,:)  :: tmpdif
#ifdef EXPERIMENTALTRANSIENTDISCREEP
double precision, allocatable, dimension(:,:)  :: tmpdislc
double precision, allocatable, dimension(:,:)  :: tmpdiffc
#endif

if (.not.LMAT) then
    return
endif

! do the properties one at a time, so that the tmp arrays do not
! have to be allocated all five at once.
!--------------------------

allocate(tmpmat(NMAT,NUMAT))

do j=1,NUMAT
    do i=1,NMAT
        tmpmat(i,j) = matdatactx%PRPMAT(i,j)
    enddo
enddo

do j=1,meshdatactx%Nelocal
    imat = meshdatactx%mat(j)!MAT(j)
    do i=1,NMAT  ! properties per material (not number of material)
        matdatactx%PRPMAT(i,j) = tmpmat(i,imat)
    enddo
enddo

deallocate(tmpmat)

!--------------------------

allocate(tmplas(NPLS,NUMAT))

do j=1,NUMAT
    if (IDIG(ICVIS,4,2).ne.0) then
        do i=1,NPLS
            tmplas(i,j) = matdatactx%PRPLAS(i,j)
        enddo
    endif
enddo

do j=1,meshdatactx%Nelocal
    imat = meshdatactx%mat(j)
    if (IDIG(ICVIS,4,2).ne.0) then
        do i=1,NPLS
            matdatactx%PRPLAS(i,j) = tmplas(i,imat)
        enddo
    endif
enddo

deallocate(tmplas)

!--------------------------

allocate(tmptem(NTEM,NUMAT))

do j=1,NUMAT
    if (MODE.gt.2) then
        do i=1,NTEM
            tmptem(i,j) = matdatactx%PRPTEM(i,j)
        enddo
    endif
enddo

do j=1,meshdatactx%Nelocal
    imat = meshdatactx%mat(j)
    if (MODE.gt.2) then
        do i=1,NTEM
            matdatactx%PRPTEM(i,j) = tmptem(i,imat)
        enddo
    endif
enddo

deallocate(tmptem)

!--------------------------


allocate(tmpcrp(NCRP,NUMAT))

do j=1,NUMAT
    if (IDIG(ICVIS,4,1).ne.0) then
        do i=1,NCRP
            tmpcrp(i,j) = matdatactx%PRPCRP(i,j)
        enddo
    endif
enddo

do j=1,meshdatactx%Nelocal
    imat = meshdatactx%mat(j)
    if (IDIG(ICVIS,4,1).ne.0) then
        do i=1,NCRP
            matdatactx%PRPCRP(i,j) = tmpcrp(i,imat)
        enddo
    endif
enddo

deallocate(tmpcrp)

!-------------------------

allocate(tmpdif(NDIF,NUMAT))

do j=1,NUMAT
    if (IDIG(ICVIS,4,3).ne.0) then
        do i=1,NDIF
            tmpdif(i,j) = matdatactx%PRPDIF(i,j)
        enddo
    endif
enddo

do j=1,meshdatactx%Nelocal
    imat = meshdatactx%mat(j)
    if (IDIG(ICVIS,4,3).ne.0) then
        do i=1,NDIF
            matdatactx%PRPDIF(i,j) = tmpdif(i,imat)
        enddo
    endif
enddo

deallocate(tmpdif)

!-----------------------

#ifdef EXPERIMENTALTRANSIENTDISCREEP
if (IDIG(ICVIS,4,1).eq.2) then
! dislocation creep parameters
allocate(tmpdislc(Ndislc,NUMAT))

do j=1,NUMAT
        do i=1,NDISLC
            tmpdislc(i,j) = matdatactx%PRPDISLC(i,j)
        enddo
enddo

do j=1,meshdatactx%Nelocal
        imat = meshdatactx%mat(j)
        do i=1,NDISLC
            matdatactx%PRPDISLC(i,j) = tmpdislc(i,imat)
        enddo
enddo

deallocate(tmpdislc)

! diffusion creep parameters
allocate(tmpdiffc(Ndiffc,NUMAT))

do j=1,NUMAT
        do i=1,NDIFFC
            tmpdiffc(i,j) = matdatactx%PRPDIFFC(i,j)
        enddo
enddo

do j=1,meshdatactx%Nelocal
        imat = meshdatactx%mat(j)
        do i=1,NDIffC
            matdatactx%PRPDIFFC(i,j) = tmpdiffc(i,imat)
        enddo
enddo

deallocate(tmpdiffc)
endif
#endif


return
end
!--------------------------------------------------------------------------------

subroutine ELMLIB (NDOF,NSD,NSTR,NEN)

USE MESHDATAMODULE,  only: meshdatactx

#ifdef SHEARZONE_TRACTIONS
USE MODELDATAMODULE, only: modeldatactx, HasSlipperyNodes, &
                           allocateinitstress, &
                           allocatepressurebcs, &
                           allocatestressbcs, &
                           allocatewinklerrestoring, &
                           allocatefaultednodes, &
                           allocateslipperynodes, &
                           allocateSlipperyElementArrays, &
                           allocatefaultparallel, &
                           allocateShearTraxArrays, &
                           showfbc, MEQ
#else
USE MODELDATAMODULE, only: modeldatactx, HasSlipperyNodes, &
                           allocateinitstress, &
                           allocatepressurebcs, &
                           allocatestressbcs, &
                           allocatewinklerrestoring, &
                           allocatefaultednodes, &
                           allocateslipperynodes, &
                           allocateSlipperyElementArrays, &
                           allocatefaultparallel, &
                                    showfbc, MEQ

#endif

#ifdef EXPERIMENTALTRANSIENTDISCREEP
use MODELTOPOLOGY, only:   NSTRTAY
#endif

#ifdef EXPERIMENTAL_FRICTION
USE MODELDEFINITION, only: ICVIS,IOPT,IPOINT,LGDEF, IRESDU, &
                           IGRAV,IVLIM,INCOMP,LinPL,NSED,ISHELL, &
                           NPREglobal, NSURFglobal, NSLSKEWglobal, &
                           NUMPRglobal, NUMSLPglobal, NUMSTRglobal, NUMWNKglobal, &
                           NUMFNglobal, NTYPE, iword, &
                           NPRE, NSLSKEW, NUMFN, NUMPR, NUMSLP, NUMSTR, NUMWNK, &
                           NWINKX, MTOT, NSURF, NUMEL, NUMNP, NUMAT, &
                           nFaultSegmentEntriesGlobal, feplt
#elif SHEARZONE_TRACTIONS
USE MODELDEFINITION, only: ICVIS,IOPT,IPOINT,LGDEF, IRESDU,&
                           IGRAV,IVLIM,INCOMP,LinPL,NSED,ISHELL, &
                           NPREglobal, NSURFglobal, NSLSKEWglobal, &
                           NUMPRglobal, NUMSLPglobal, NUMSTRglobal, NUMWNKglobal, &
                           NUMFNglobal, NTYPE, &
                           NPRE, NSLSKEW, NUMFN, NUMPR, NUMSLP, NUMSTR, NUMWNK, &
                           NWINKX, NSURF, NUMEL, NUMNP, NUMAT, &
                           NUMELSHEARZONE, feplt
#else
USE MODELDEFINITION, only: ICVIS,IOPT,IPOINT,LGDEF, IRESDU, &
                           IGRAV,IVLIM,INCOMP,LinPL,NSED,ISHELL, &
                           NPREglobal, NSURFglobal, NSLSKEWglobal, &
                           NUMPRglobal, NUMSLPglobal, NUMSTRglobal, NUMWNKglobal, &
                           NUMFNglobal, NTYPE, iword, &
                           NPRE, NSLSKEW, NUMFN, NUMPR, NUMSLP, NUMSTR, NUMWNK, &
                           NWINKX, MTOT, NSURF, NUMEL, NUMNP, NUMAT, feplt
#endif
USE MATERIALSMODULE, only: LMAT, &
                           allocatematdatamech
USE MODELCTX,        only: getrank
USE ALGEBRA,         only: CLEAR, idig
use formatsmodule,   only: FMT_tecin_dat_Nsurf   ,    &
                           FMT_tecin_dat_ELMswitches, &
                           FMT_tecin_dat_ELMcounts,   &
                           FMT_tecin_partf_elm,       &
                           FMT_tecin_nSurface_nodes
use debugmodule,     only: debug, &
                           iecho, &
                           parallelLog, &
                           allocateError, &
                           Verified_MPI_Barrier, &
                           xit,itest
use timestepmodule,  only: nintg
use iomodule          

#ifdef SPARSE
use AOmodule,        only     :  nequations
#else
use modeldefinition, only     :  NEQlocal
#endif

!
! Assembles data on element properties, adjusts ID arrays, and constructs LM
! and IDIAG arrays.
!

implicit none
!-pass
integer NDOF,NSTR,NEN
!-locl
logical docrp,dodif,doplas
character(len=35) ELMLBL
dimension ELMLBL(6)
integer ERROR
character(len=100) record
integer :: ios
integer :: NSD
integer :: i

! for side flipping
!integer :: side1, side2, side3, neighbors(3), sideID
!integer, allocatable :: slipSides(:)

!integer :: iPoint
!-init
data ELMLBL/ &
 'PLANE STRESS-STRAIN','OUT-OF-PLANE STRESS-STRAIN', &
 'AXISYMMETRIC','SPHERICAL','3D BRICKS','3D TETRAHEDRA'/
 
if (iecho.eq.2 .or. debug) then 
    write(stdout,'(100a)') 'read IOPT,IPOINT,LGDEF,IRESDU,IGRAV,IVLIM,INCOMP,LinPL,NSED,ISHELL (10i5)'
    write(stdout,*) "with format", FMT_tecin_dat_ELMswitches
endif 
    
!call Verified_MPI_Barrier("in elmlib")

! these NUM*** variables reside in the modeldefinition module
read(luin, FMT_tecin_dat_ELMswitches, err=200,end=210,iostat=ios) & 
   IOPT,IPOINT,LGDEF, IRESDU,IGRAV,IVLIM,INCOMP,LinPL,NSED,ISHELL

if (iecho.eq.2 .or. debug) then
    write(stdout,*) 'read: IOPT,IPOINT,LGDEF,IRESDU,IGRAV,IVLIM,INCOMP,LinPL,NSED,ISHELL (10i5)'
    write(stdout,*) 'IOPT                  :', IOPT
    write(stdout,*) 'IPOINT                :', IPOINT
    write(stdout,*) 'Large deformation     :', LGDEF
    write(stdout,*) 'Residual force update :', IRESDU
    write(stdout,*) 'Gravity present       :', IGRAV
    write(stdout,*) 'IVLIM                 :', IVLIM
    write(stdout,*) 'Incompressibily       :', INCOMP
    write(stdout,*) 'LinPL                 :', LinPL
    write(stdout,*) 'Sediments             :', NSED
    write(stdout,*) 'Domain on the shell   :', ISHELL
endif


if (ishell.eq.1) then
    ! the coordinates must not be too extravagant.
    ! Check whether this is OK
    do i=1,meshdatactx%nvlocal
        if ((meshdatactx%X(1,i).lt.-180) .or.  (meshdatactx%X(1,i).gt. 180)) then
            write(stderr,*) "model on the shell, but X coordinate of point", i, &
               "is", meshdatactx%X(1,i) 
            write(stderr,*) "Should be inside the [-180, 180] interval"
            call xit(1,' ')
        endif
        if ((meshdatactx%X(2,i).lt.-89) .or.  (meshdatactx%X(2,i).gt. 89)) then
            write(stderr,*) "model on the shell, but Y coordinate of point", i, &
               "is", meshdatactx%X(2,i)
            write(stderr,*) "Should be inside the [-89, 89] interval"
            call xit(1,' ')
        endif
        enddo
endif

if (iecho.eq.2 .or. debug) then
    write(stdout,*) 'read NSURF ', FMT_tecin_nSurface_nodes
endif

!########## read number of surface elements #####################
read(luin,FMT_tecin_dat_Nsurf,err=100,end=110,iostat=ios) NSURF


!############# read element config line ##########

read(luin,'(a)',err=300,end=310,iostat=ios) record

if (record(12:12).eq." ") then
    write(stderr,*) "Input reader says: Something wrong with the number of pressure BCS; position 12 empty"
endif
if (record(24:24).eq." ") then
    write(stderr,*) "Input reader says: Something wrong with the number of stress BCS; position 24 empty"
endif
if (record(36:36).eq." ") then
    write(stderr,*) "Input reader says: Something wrong with the number of Winkler pressure BCS; position 36 empty"
endif
if (record(48:48).eq." ") then
    write(stderr,*) "Input reader says: Something wrong with the number of slippery nodes; position 48 empty"
endif
if (record(60:60).eq." ") then
    write(stderr,*) "Input reader says: Something wrong with the number of fault parallel surfaces; position 60 empty"
endif
if (record(72:72).eq." ") then
    write(stderr,*) "Input reader says: Something wrong with the number of faulted nodes; position 72 empty"
endif
if (record(84:84).eq." ") then
    write(stderr,*) "Input reader says: Something wrong with the number of prestress BCS; position 84 empty"
endif

#ifdef EXPERIMENTAL_FRICTION
read(record,FMT_tecin_dat_ELMcounts,err=300,end=310,iostat=ios) &
   NUMPR,NUMSTR, NUMWNK,NUMSLP,NSLSKEW,NUMFN,NPRE,nFaultSegmentEntriesGlobal
#elif SHEARZONE_TRACTIONS
read(record,FMT_tecin_dat_ELMcounts,err=300,end=310,iostat=ios) &
   NUMPR,NUMSTR, NUMWNK,NUMSLP,NSLSKEW,NUMFN,NPRE, NUMELSHEARZONE
#else
read(record,FMT_tecin_dat_ELMcounts,err=300,end=310,iostat=ios) &
   NUMPR,NUMSTR, NUMWNK,NUMSLP,NSLSKEW,NUMFN,NPRE
#endif

if (iecho.eq.2) then
#ifdef EXPERIMENTAL_FRICTION
    write(stdout,*) 'read: ', &
      NUMPR,NUMSTR, NUMWNK,NUMSLP,NSLSKEW,NUMFN,NPRE,nFaultSegmentEntriesGlobal
#elif SHEARZONE_TRACTIONS
    write(stdout,*) 'read: ', &
      NUMPR,NUMSTR, NUMWNK,NUMSLP,NSLSKEW,NUMFN,NPRE, NUMELSHEARZONE
#else
    write(stdout,*) 'read: ', &
      NUMPR,NUMSTR, NUMWNK,NUMSLP,NSLSKEW,NUMFN,NPRE                       
#endif
endif

if (iecho.eq.2 .or. debug) then
    write(stdout,*) "rank", getrank(),"read: Number of pressure BCS: ", NUMPR
    write(stdout,*) "rank", getrank(),"read: Number of stress BCS: ", NUMSTR
    write(stdout,*) "rank", getrank(),"read: Number of Winkler pressures", NUMWNK
    write(stdout,*) "rank", getrank(),"read: Number of slippery node entries: ", NUMSLP
    write(stdout,*) "rank", getrank(),"read: Number of fault parallel entries: ", NSLSKEW
    write(stdout,*) "rank", getrank(),"read: Number of faulted node entries ", NUMFN
    write(stdout,*) "rank", getrank(),"read: Number of prestress BCS ", NPRE
#ifdef EXPERIMENTAL_FRICTION
    write(stdout,*) "rank", getrank(),"read: Number of fault segments",nFaultSegmentEntriesGlobal
#elif SHEARZONE_TRACTIONS
    write(stdout,*) "rank", getrank(),"read: Number of shear zone entries", NUMELSHEARZONE
#endif
endif

#ifdef EXPERIMENTAL_FRICTION
if (nFaultSegmentEntriesGlobal.ne.0 .and. NUMSLP.eq.0) then
    write(stderr,*) "friction entries are supplied, but no slippery nodes"
    write(stderr,*) "This inconsistent"
    call xit(1,' ')
endif
#elif SHEARZONE_TRACTIONS
if (NUMELSHEARZONE.ne.0 .and. NUMSLP.eq.0) then
    write(stderr,*) "Shear zone traction entries given, but no slippery nodes"
    write(stderr,*) "This inconsistent"
    call xit(1,' ')
endif
#endif

if (OUTFIL(1)) then
    write(luout(1),3)
 3       format(1X//1X,'S W I T C H E S'/)
    if (NTYPE.eq.1) write(luout(1),4) IOPT
 4       format(6X, &
     'ANALYSIS OPTION . . . . . . . . . . . . . . (IOPT ) =', &
     I5/8X,'EQ.0, PLANE STRAIN'/8X,'EQ.1, PLANE STRESS'/)
    if (NTYPE.eq.1) write(luout(1),5) ISHELL
 5       format(6X, &
     'SPHERICAL SHELL OPTION  . . . . . . . . .  (ISHELL) =', &
     I5/8X,'EQ.0, FLAT PLATE'/8X,'EQ.1, SPHERICAL PLATE'/)
    write(luout(1),6) IPOINT,NEN
 6       format(6X, &
     'R H S  INTEGRATION CODE . . . . . . . . .  (IPOINT) =', &
     I5/8X, &
     'EQ.1, MID-POINT GAUSSIAN QUADRATURE FOR FORCES',/8X, &
     'EQ.',I1,', GAUSSIAN QUADRATURE FOR FORCES'/)
    write(luout(1),7) LGDEF,IRESDU
 7       format(6X, &
     'LARGE DEFORMATION STRAIN EMPLOYED . . . . . (LGDEF) =', &
     I5/8X, &
     'EQ.0, LINEAR STRAIN-DISPLACEMENT MATRIX'/8X, &
     'EQ.1, UPDATED LARGE DEFORMATION FORMULATION'//6X, &
     'LOAD VECTOR UPDATED BY UNBALANCED FORCES   (IRESDU) =', &
     I5/8X, &
     'EQ.0, LOAD NOT UPDATED'/8X, &
     'EQ.1, LOAD NOT UPDATED, ECHO UNBALANCE'/8X, &
     'EQ.2, LOAD UPDATED (NORMAL FOR LGDEF=1)'/)
    write(luout(1),8) IGRAV
 8       format(6X, &
     'GRAVITY OPTION  . . . . . . . . . . . . . . (IGRAV) =', &
     I5/8X, &
     'EQ.0, GRAVITY OFF',/8X, &
     'EQ.1, GRAVITY ON'/)
    write(luout(1),9) IVLIM
 9       format(6X, &
     'LOWER BOUND VISCOSITY SWITCH  . . . . . . . (IVLIM) =', &
     I5/8X,'EQ.0, NO VISCOSITY LIMITS',/8X, &
     I5/8X,'EQ.1, STATIC VISCOSITY LIMITS',/8X, &
     'EQ.2, MINIMUM AND MAXIMUM VISCOSITY LIMITED BY'/14X, &
     'TIME STEP SIZE AND STATIC LIMITS'/)
     if (NTYPE.eq.1) write(luout(1),10) INCOMP
10        format(6X, &
      'INCOMPATIBLE MODES  . . . . . . . . . . .  (INCOMP) =', &
      I5/8X,'EQ.0, NO INCOMPATIBLE MODES'/8X, &
      'EQ.1, TWO ADDITIONAL DISPLACEMENT MODES'/)
     write(luout(1),11) LinPL
11        format(6X, &
      'Linearise the power law . . . . . . . . . . (LinPL) =', &
      I5/8X,'EQ.1, Linearise'/8X,'EQ.0, do not linearise'/)
    write(luout(1),12) NSED
12       format(6X, &
     'EROSION/SEDIMENTATION OPTION . . . . . . . . (NSED) =', &
     I5/8X, &
     'EQ.0. NO SEDIMENT TRANSPORT LOADS'/8X, &
     'EQ.1. SEDIMENT LOADS INCLUDED'//)
    write(luout(1),13) ELMLBL(NTYPE),NUMEL
13       format(1X///1X, &
     'E L E M E N T    S Y S T E M   D A T A'//6X, &
     'ELEMENT TYPE:  ',A35//6X, &
     'NUMBER OF ELEMENTS  . . . . . . . . . . . . (NUMEL) =', &
     I8/)
    write(luout(1),14) NSURF
14       format(6X, &
     'NUMBER OF UPPER SURFACE NODE ENTRIES   . . .(NSURF) =', &
     I8//)

    write(luout(1),15)
15       format(1X//1X,'E L E M E N T',3X,'B O U N D A R Y',3X, &
     'C O N D I T I O N S'/)

#ifdef EXPERIMENTAL_FRICTION
    write(luout(1),16) NUMPR,NUMSTR,NUMWNK,NUMSLP,NSLSKEW,NUMFN, &
     NPRE,nFaultSegmentEntriesGlobal
#elif SHEARZONE_TRACTIONS
    write(luout(1),16) NUMPR,NUMSTR,NUMWNK,NUMSLP,NSLSKEW,NUMFN, &
     NPRE, NUMELSHEARZONE
#else
    write(luout(1),16) NUMPR,NUMSTR,NUMWNK,NUMSLP,NSLSKEW,NUMFN, &
     NPRE
#endif

16 format(6X, &
     'NUMBER OF ELEMENT PRESSURE ENTRIES . . . .  (NUMPR) =', &
     I8//6X, &
     'NUMBER OF ELEMENT STRESS ENTRIES . . . . . (NUMSTR) =', &
     I8//6X, &
     'NUMBER OF WINKLER PRESSURE ENTRIES . . . . (NUMWNK) =', &
     I8//6X, &
     'NUMBER OF SLIPPERY NODE ENTRIES  . . . . . (NUMSLP) =', &
     I8//6X, &
     'NUMBER OF FAULT PARALLEL ELEMENT ENTRIES  (NSLSKEW) =', &
     I8//6X, &
     'NUMBER OF FAULTED NODE ENTRIES . . . . . .  (NUMFN) =', &
     I8//6X, &
#ifdef EXPERIMENTAL_FRICTION
     'NUMBER OF PRE-STRESS ENTRIES  . . . . . . .  (NPRE) =', &
     I8//6X, &
     'NUMBER OF FRICTION ENTRIES (nFaultSegmentEntriesGlobal) =', &
     I8/)
#elif SHEARZONE_TRACTIONS
     'NUMBER OF PRE-STRESS ENTRIES  . . . . . . .  (NPRE) =', &
     I8//6X, &
     'NUMBER OF SHEARZONE ENTRIES  . . . (NUMELSHEARZONE) =', &
     I8//6X)
#else
     'NUMBER OF PRE-STRESS ENTRIES  . . . . . . .  (NPRE) =', &
     I8//6X)
#endif

    write(luout(1),19) NUMAT
19  format(//1X, &
     'M A T E R I A L   S E T   D A T A'//6X, &
     'NUMBER OF MATERIAL SETS  . . . . . . . . .  (NUMAT) =', &
     I8//)
endif

call chknpar(NEN)

! copy the values read to their global counterparts
NUMPRglobal   = NUMPR
NUMSTRglobal  = NUMSTR
NUMWNKglobal  = NUMWNK
NUMSLPglobal  = NUMSLP
NSLSKEWglobal = NSLSKEW
NUMFNglobal   = NUMFN
NPREglobal    = NPRE
NSURFglobal   = NSURF


if (IDIG(ICVIS,4,4).ne.0) then
    if (allocated(modeldatactx%VPSTR)) then
        deallocate(modeldatactx%VPSTR)
    endif

    ALLOCATE(modeldatactx%VPSTR(NSTR,meshdatactx%nelocal), STAT=ERROR)
endif


docrp  = (IDIG(ICVIS,4,1).ne.0)
dodif  = (IDIG(ICVIS,4,3).ne.0 .or. IDIG(ICVIS,4,4).ne.0)
doplas = (IDIG(ICVIS,4,2).ne.0)


if (LMAT) then
    call allocatematdatamech(NUMEL, NDOF, docrp,dodif,doplas)
else
    call allocatematdatamech(NUMAT, NDOF, docrp,dodif,doplas)
endif

if (allocated(modeldatactx%STN)) then
    deallocate(modeldatactx%STN)
endif
ALLOCATE(modeldatactx%STN(NSTR,meshdatactx%nelocal), STAT=ERROR)
call AllocateError("Modeldata STN", ERROR)

if (allocated(modeldatactx%STR)) then
    deallocate(modeldatactx%STR)
endif
ALLOCATE(modeldatactx%STR(NSTR,meshdatactx%nelocal), STAT=ERROR)
call AllocateError("Modeldata STR", ERROR)

#ifdef EXPERIMENTALTRANSIENTDISCREEP
! taylor stress
  ALLOCATE(modeldatactx%STNTAY(NSTRTAY,meshdatactx%nelocal), STAT=ERROR)
  call AllocateError("Modeldata STNTAY", ERROR)
! taylor stress time derivative
  ALLOCATE(modeldatactx%STNTAYRATE(NSTRTAY,meshdatactx%nelocal), STAT=ERROR)
  call AllocateError("Modeldata STNTAYRATE", ERROR)
! threshold stress, scalar
  ALLOCATE(modeldatactx%STNTHRES(meshdatactx%nelocal), STAT=ERROR)
  call AllocateError("Modeldata STNTHRES", ERROR)
! dislocation creep effective strain rate, scalar
  ALLOCATE(modeldatactx%RATEDISLEFF(meshdatactx%nelocal), STAT=ERROR)
  call AllocateError("Modeldata RATEDISLEFF", ERROR)
! dislocation creep effective viscosity
    ALLOCATE(modeldatactx%INVEFFVISCDISL(meshdatactx%nelocal), STAT=ERROR)
  call AllocateError("Modeldata EFFVISCDISL", ERROR)
! diffusion creep effective viscosity  
    ALLOCATE(modeldatactx%INVEFFVISCDIFF(meshdatactx%nelocal), STAT=ERROR)
  call AllocateError("Modeldata EFFVISCDIFF", ERROR)
! composite (dislocation + diffusion) creep effective viscosity  
!    ALLOCATE(modeldatactx%EFFVISC(meshdatactx%nelocal), STAT=ERROR)
!  call AllocateError("Modeldata EFFVISC", ERROR)
#endif 


#ifndef EXPERIMENTAL_NORATESTORAGE
if (NINTG.gt.0) then
    if (allocated(modeldatactx%rate)) then
        deallocate(modeldatactx%rate)
    endif
    ALLOCATE(modeldatactx%RATE(NSTR,NUMEL), STAT=ERROR)
    call AllocateError("Modeldata strainrate", ERROR)
endif
#endif

 
if (showfbc) then
    ! generate alternative equation numbers
    ! MEQ                   number of nodal displacement bcs

    allocate(modeldatactx%forShowfbcIAeqnr(NDOF,NUMNP))
    allocate(modeldatactx%forShowfbcIAXeqnr(NDOF,NUMNP))

    modeldatactx%forShowfbcIAeqnr = 0
    modeldatactx%forShowfbcIAXeqnr = 0

    call SETNEQBC (modeldatactx%IBOND,modeldatactx%forShowfbcIAeqnr,NDOF,NUMNP,MEQ)

endif


!##########  read material data ##########
call EDATA  (NDOF)

!call Verified_MPI_Barrier("almost through elmlib")


!if (iecho.eq.8) then
!    call ParallelLog("readfe", "passed barrier")
!endif



call INCCHK (meshdatactx%IEN,NEN)
if (NPRE.gt.0) then
    call allocateinitstress(modeldatactx,NSTR,NPRE)
endif

npre = 0

   call RDPRE  (NSTR)

if (iecho.eq.8) then
    call ParallelLog("readfe", "read prestress")
endif


!     read element pressure bc
if (NUMPR.ne.0) then
    call allocatepressurebcs(modeldatactx, IABS(NUMPR))
    call ELPRES (modeldatactx%IELNO, modeldatactx%ISIDE, modeldatactx%PRES)
    call FLXCHK (modeldatactx%IELNO, modeldatactx%ISIDE, &
                 IABS(NUMPR),NEN,'pressure data')
endif

!     read stress boundary bc
if (NUMSTR.ne.0) then
    call allocatestressbcs(modeldatactx,nstr,numstr)
    call ELSTRS (modeldatactx%IELSTR,modeldatactx%ISSIDE, &
                 modeldatactx%STRS,modeldatactx%ISTR,NSTR)
    call FLXCHK (modeldatactx%IELSTR,modeldatactx%ISSIDE, &
                 NUMSTR,NEN,'stress data')
endif

!     read winkler restoring pressures
if (NUMWNK.ne.0) then
    call allocatewinklerrestoring(modeldatactx,NUMWNK)
    call ELWINK ()
    call FLXCHK (modeldatactx%IWELM,modeldatactx%IWSIDE, &
                 NUMWNK,NEN,'Winkler data')
endif

! read faulted nodes
if (NUMFNglobal.gt.0) then
    call allocatefaultednodes(modeldatactx, &
      ndof,nen,numfn,meshdatactx%Nelocal)

!   this allocates a bit too much actually
    call RFAULT (modeldatactx%NFAULT, modeldatactx%FAULT, &
      meshdatactx%IEN,NEN,NDOF)
endif


!     read slippery nodes
!     first put slippery node data (element+nodal point) in global array
!     since this is needed for later renumering of equations in adjid
!     routine --> this is the NSLIPglobal array; the NSLIP array
!     contains the per processor data, based on element
!     allocate arrays related to slippery nodes, this includes arrays
!     like IDX, LMX etc

if (HasSlipperyNodes()) then

    if (iecho.eq.8) then
        call ParallelLog("readfe", "reading slippery nodes")
    endif

    call allocateslipperynodes(modeldatactx,NDOF,NEN,NUMSLP, &
      meshdatactx%Nvlocal,meshdatactx%Nvglobal,meshdatactx%Nelocal)

    if (iecho.eq.8) then
        call ParallelLog("readfe", "calling rslip")
    endif

    call RSLIP (modeldatactx%DIFORC,meshdatactx%IEN,NEN,NDOF,NSD)

#ifdef Lukas
    ! correct slippery nodes for corners in faults.
    ! Note: this is scary:
    allocate(slipSides(meshdatactx%neglobal))
    slipSides = 0
    if (NUMSLPglobal .gt. 0) then
    
        ! step through the slippery array
        do iSlip = 1, NUMSLPglobal
                elemID = modeldatactx%NSLIPglobal(1,iSlip)
                sideID = modeldatactx%NSLIPglobal(3,iSlip)

                if (slipSides(elemID).eq.0) then
                        ! this element has not yet another slippery node to it.
                        slipSides(elemID) = sideID
                elseif (slipSides(elemID).eq.-sideID) then
                        ! there is a conflict in this elements.
                        slipSides(elemID) = 2
                elseif (slipSides(elemID).eq.2) then
                        ! this node has already got a prior conflict recorded. Do nothing.
                elseif (slipSides(elemID).eq.sideID) then
                        ! do nothing. This element already has the good side
!               else
!                       write(*,*) "Something is fishy with the weight of the slippery element"
!                       write(*,*) "Should be 1, 0, or -1, but is", slipSides(elemID)
                endif
        enddo

!        do iSlip = 1, NUMSLPglobal

!            elemID = modeldatactx%NSLIPglobal(1,iSlip)
!            sideID = modeldatactx%NSLIPglobal(3,iSlip)

!            if (slipSides(elemID).eq.2) then
                ! there is a conflicting element.

! GRRRRRRRR fuckit
!                neighbors = meshdatactx%ElementNeighbours(elemID,:)

!                side1 = slipSides(neighbors(1))
!                side2 = slipSides(neighbors(2))
!                side3 = slipSides(neighbors(3))

!                if ((abs(side1) + abs(side2) + abs(side3)).eq.3) then

!                    if (( (side1 + side2 + side3) * sideID) .eq. -1) then

!                        modeldatactx%NSLIPglobal(3,iSlip) = -sideID
!                    endif
!                endif

!            endif
!        enddo


endif


!---- finished flipping slips
#endif

    if (iecho.eq.8) then
        call ParallelLog("readfe", "done rslip")
    endif

if (iecho.eq.2) write(*,*) 'allocate slippery element arrays'
    call allocateSlipperyElementArrays()

    if (iecho.eq.8) then
        call ParallelLog("readfe", "read and stored slippery nodes")
    endif

endif

if (itest.eq.1) then
    if (iecho.eq.2) write(*,*) 'check coincident slippery and split nodes'
    ! verify whether elements/nodes appear in both split and slippery arrays
    call CHKFAULTSLIP(modeldatactx%NFAULT,modeldatactx%NSLIP,NUMFN,NUMSLP)
    if (iecho.eq.1) write(*,*) 'No coincident slippery and split nodes'
endif

! read angles
if (NSLSKEW.gt.0) then

    if (iecho.eq.2) write(*,*) 'allocate and read fault paralle interfaces'

    call allocatefaultparallel(modeldatactx,NSLSKEW)

    if (iecho.eq.8) then
        call ParallelLog("readfe", "calling slipperySkew")
    endif

! reading fault parallel interfaces
    call RSlipperySKEW (modeldatactx%NSELSD, &
                       modeldatactx%NSLIP, &
                       modeldatactx%SKEW, &
                        meshdatactx%IEN, NSD,NDOF,NEN)

    if (iecho.eq.8) then
        call ParallelLog("readfe", "finished slipperySkew")
    endif

endif


! read differential winkler forces
if (NUMSLPglobal.gt.0) then

    if (iecho.eq.8) then
        call ParallelLog("readfe", "starting rwink")
    endif

    if (iecho.eq.2) write(*,*) 'read differential Winkler data'
    call RWINK (meshdatactx%X, &
                NSD,NDOF, &
                meshdatactx%Nvlocal, &
                NWINKX,1)

    if (iecho.eq.8) then
        call ParallelLog("readfe", "finished rwink")
    endif
endif

#ifdef SHEARZONE_TRACTIONS
! add shear zone tractions
if (NUMELSHEARZONE.gt.0) then
    call allocateShearTraxArrays()
    call readShearTrax()
    if (allocated(modeldatactx%shearF)) then
        deallocate(modeldatactx%shearF)
    endif
#ifdef SPARSE
    allocate(modeldatactx%shearF(nequations), STAT=ERROR)
#else
    allocate(modeldatactx%shearF(NEQlocal), STAT=ERROR)
#endif /* SPARSE */
    call AllocateError("shearF", ERROR)
endif
#endif /* SHEARZONE_TRACTIONS */


#ifdef EXPERIMENTAL_FRICTION
if (nFaultSegmentEntriesGlobal.gt.0) then
    call readFriction()
endif
#endif

! read surface nodes
if (NSURF.gt.0) then
    if (allocated(modeldatactx%ISURF)) then
        deallocate(modeldatactx%ISURF)
    endif
    ALLOCATE(modeldatactx%ISURF(NSURF), STAT=ERROR)
    call AllocateError("Modeldata ISurface", ERROR)
    
    if (allocated(modeldatactx%ISURFORDER)) then
        deallocate(modeldatactx%ISURFORDER)
    endif
    ALLOCATE(modeldatactx%ISURFORDER(NSURF), STAT=ERROR)
    call AllocateError("Modeldata Isurf order", ERROR)
    
    if (allocated(modeldatactx%SURF)) then
        deallocate(modeldatactx%SURF)
    endif
    ALLOCATE(modeldatactx%SURF(NSD,NSURF), STAT=ERROR)
    call AllocateError("Modeldata Surface", ERROR)
    if (allocated(modeldatactx%SURF0)) then
        deallocate(modeldatactx%SURF0)
    endif
    ALLOCATE(modeldatactx%SURF0(NSD,NSURF), STAT=ERROR)
    call AllocateError("Modeldata Surface 0", ERROR)

    call RDSURF (modeldatactx%ISURF,NSURF)
    if (iecho.eq.8) then
        call ParallelLog("readfe", "finished rdsurf")
    endif
endif

if (feplt) then
    ! calling program is plt3d, plnplt, axiplt or opnplt: skip equations stuff
else
    call ADJID (NDOF)

    ! updates the local ID array for potential slippery nodes
    if (iecho.eq.2) write(*,*) 'update equation numbers for slippery nodes'
    call updatelocalequationid(meshdatactx, modeldatactx,ndof)
    if (iecho.eq.8) then
        call ParallelLog("readfe", "updated local equations")
    endif

    ! do linked nodes (adds no new equations): updates only the IDglobal array
    if (iecho.eq.2) write(*,*) 'update equation numbers for linked nodes'
    call addlinkednodes(meshdatactx, modeldatactx)
    if (iecho.eq.8) then
        call ParallelLog("readfe", "Added linked nodes")
    endif
endif

! add faulted nodes, add no new equations
if (iecho.eq.2) write(*,*) 'create index list for faulted nodes'
call LOCALF (modeldatactx%NFAULT,meshdatactx%IEN,NDOF,NUMFN,NUMEL,NEN)
if (iecho.eq.8) then
    call ParallelLog("readfe", "Localised nfault array")
endif

if (.not.feplt) then
    if (iecho.eq.2) write(*,*) 'localizing slippery node equation numbers'
    call LOCALX (modeldatactx%IDXglobal, modeldatactx%NSLIP, NUMSLP)
    if (iecho.eq.8) then
        call ParallelLog("readfe", "finished localx")
    endif
endif   

if (iecho.eq.8) then
    call ParallelLog("readfe", "read all tecin")
endif

return

100 write(stderr,101) ios
101 format(1x,'ELMLIB: read error ',i4,' while reading NUMEL,NUMAT,NSURF')
call prerr('ELMLIB')
call xit(1,"Error while reading element data")

110   write(stderr,111) ios
111   format(1x,'ELMLIB: premature EOF error ',i4,' while reading ', &
 'NUMEL,NUMAT,NSURF')
call prerr('ELMLIB')
call xit(1,"Unexpected EOF while reading element data")

200   write(stderr,201) ios
201   format(1x,'ELMLIB: read error ',i4,' while reading IOPT,', &
 'IPOINT,LGDEF,IRESDU,IGRAV,IVLIM,INCOMP,LinPL,NSED')
write(stderr,202)

202   format(1x,'ELMLIB: possibly G-TECTON is trying to read an ', &
 'older input file format?')
call prerr('ELMLIB')
call xit(1,"Error while reading element data; File format incorrect?")

210   write(stderr,211) ios
211   format(1x,'ELMLIB: premature EOF error ',i4,' while reading ', &
 'IOPT,IPOINT,LGDEF,IRESDU,IGRAV,IVLIM,INCOMP,LinPL,NSED')
call prerr('ELMLIB')
call xit(1,"Element switches")

300   write(stderr,301) ios
301   format(1x,'ELMLIB: read error ',i4,' while reading ', &
 'NUMPR,NUMSTR,NUMWNK,NUMSLP,NSLSKEW,NUMFN,NPRE')
write(stderr,202)
call prerr('ELMLIB')
call xit(1,"Element boundary condition quantities")

310   write(stderr,311) ios
311   format(1x,'ELMLIB: premature EOF error ',i4,' while reading ', &
 'NUMPR,NUMSTR,NUMWNK,NUMSLP,NSLSKEW,NUMFN,NPRE')
call prerr('ELMLIB')
call xit(1,"Unexpected EOF while reading element boundary condition quantities")

end subroutine
!-------------------------------------------------------------------------------
subroutine ELNODE (NEN)

!
! Subroutine reads and prints out the IEN element-node array
!

USE MESHDATAMODULE, only: meshdatactx
USE MODELCTX,       only: getrank
USE MODELDEFINITION
#ifdef SPARSE
use debugmodule,    only: parallellog, &
                          iecho, &
                          debug
#endif

use petscksp
use iomodule          

#ifdef SPARSE
use AOmodule
#endif

implicit none
!-pass
integer NEN
!-locl
integer          :: i,n,nn
character(len=4) :: HEAD(8)
!integer          :: nodelist(meshdatactx%Nvlocal)
!integer          :: j, ntot  ! for connectivity test, see below
!-init
data HEAD/8*'NODE'/

call ELGEN (NEN)

#ifdef SPARSE
! construct a lookup table to find local element vertices,
! for parallel runs

if (iecho.eq.8) then
    call ParallelLog("Elnode", "building lookup tables")
endif

if (debug) then
    write(*,*)  "rank", getrank(),"filling locelt"
endif

call FillLocElt()       

if (debug) then
    write(*,*)  "rank", getrank(),"done filling locelt"
endif

if (iecho.eq.8) then
    call ParallelLog("Elnode", "Finished loookup tables")
endif

#endif

!
!       test node connectivity
!
!do i=1,meshdatactx%Nelocal
!    do j=1,NEN
!        if (nodelist(meshdatactx%IEN(j,i)).ne.1) then
!            nodelist(meshdatactx%IEN(j,i)) = 1
!            ntot = ntot + 1
!        endif
!    enddo
!enddo
!if (ntot.lt.meshdatactx%Nvlocal) then
!    do i=1,NUMNP
!        if (nodelist(i).eq.0) then
!            write(stderr,4) i
! 4               format(' Node ',I12,' is not in any element')
!            call exitp(1)
!        endif
!    enddo
!endif
!!!!!!!!!!!!!!!!!!!!!

!        test connectivity old style; above is much faster
!        do 100 i=1,NUMNP
!            do n=1,NUMEL
!               do j=1,NEN
!                    if (IEN(j,n).eq.i) goto 100
!                enddo
!            enddo
!            write(stderr,4) i
!    4       format(' Node ',I8,' is not in any element')
!            call exitp(1)
!100     continue
!

if (.not.OUTFIL(1)) return
 
nn=0
do n=1,meshdatactx%Nelocal
    if (MOD(nn,50).eq.0) then
        write(luout(1),1) (HEAD(i),i,i=1,NEN)
 1           format(/// &
         ' E L E M E N T  D A T A ',//5X, &
         ' ELEMENT    MATERIAL  ',8(A4,I2,4X))
        write(luout(1),2)
 2           format(5X,'  NUMBER '/)
    endif
    nn = nn + 1
    write(luout(1),3) n,meshdatactx%MAT(n), &
           (meshdatactx%IEN(i,n),i=1,NEN)
 3       format(3X,I8,9(2X,I8))
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine INCCHK (IEN,NEN)

USE MODELDEFINITION, only: INCOMP, ntype
use meshdatamodule,  only: meshdatactx
use debugmodule,     only: xit
use iomodule          

!
!       Subroutine checks if only quadrilateral/hexahedron elements
!       are present in case of incompatible modes
!

implicit none
!-pass
integer NEN
integer IEN
dimension IEN(NEN,*)
!-local
integer :: n

if (INCOMP.eq.1) then
    do n=1,meshdatactx%nelocal
        if (IEN(3,n).eq.IEN(4,n)) then
            if (NTYPE.eq.1) then
                write(stderr,1) n
            else if (NTYPE.eq.5) then
                write(stderr,5) n
            else
                write(stderr,6)
6                       format(1x,'INCCHK: NTYPE error')
!                       call xit(1," ")
            endif
1                   format(1x,'Incompatible modes error in element ',I12, &
             ': only quadrilateral elements allowed')
5                   format(1x,'Incompatible modes error in element ',I12, &
             ': only brick elements allowed')
            call xit(1,"Incompatible modes in elements; only bricks allowed")
        endif
    enddo
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine ELPRES (IELNO,ISIDE,PRES)

USE MODELDEFINITION, only : NUMPR, NUMPRglobal
use debugmodule,     only : iecho, delay
use modelctx,        only : getrank
use formatsmodule,   only : FMT_tecin_pressure
use iomodule          
!
! Subroutine to read and print out pressure load data
!     reads per element per processor
!

implicit none
!-pass
integer :: IELNO, ISIDE
double precision :: PRES
dimension IELNO(*),ISIDE(*),PRES(*)
!-locl
integer :: n, nn

!
if (iecho.eq.2) then
    write(stderr,1) FMT_tecin_pressure
 1   format(1x,'if (NUMPR <> 0){'/ &
 5X,'*** read element pressure data ***'// &
 5X,'read ',A14,' IELNO(n),ISIDE(n),PRES(n)'/ &
 5X,'if (IELNO(n)==0) return'/1x,'}')
endif
!
if (NUMPRglobal.eq.0) then
    return
endif
!
call FLXGEN (IELNO,ISIDE,PRES,IABS(NUMPRglobal),NUMPR,1,FMT_tecin_pressure,'ElPres')

!
if (.not. OUTFIL(1)) then
    return
endif
!
nn = 0
do n=1,IABS(NUMPR)
    if (MOD(nn,50).eq.0) then
        write(luout(1),2)
 2           format(///1X, &
         'P R E S S U R E   L O A D   D A T A',//6X, &
         'ELEMENT    SIDE',10X,'PRESSURE',/7X, &
         'NUMBER   NUMBER'/)
    endif
    write(luout(1),3) IELNO(n),ISIDE(n),PRES(n)
 3       format(3X,I12,4X,I5,1PE20.8)
    nn = nn + 1
enddo
!
return
end
!-------------------------------------------------------------------------------

subroutine FGEN (A,M,MAX,FMT,nread,ReadingWhat)

USE MODELCTX,       only : getrank, getsize
#ifdef SPARSE
USE MESHDATAMODULE, only : meshdatactx
#endif
USE MODELDEFINITION
use debugmodule,    only : iecho, &
                           parallelLog, &
                           delay, xit, debug
#ifdef SPARSE
use AOmodule,       only : vertices, &
                           verticesmask
#endif
use iomodule          

!
!    reads and generates floating point node data
!

implicit none
!-pass
integer M,MAX  ! m is number of variables per point. 
               ! Usually this is NDOF, but for initial temp it is 1,
               ! independent of the number of physical dimensions.

character FMT*(*)
double precision :: A
dimension A(M,MAX)

!-locl
logical end
character(len=68) record
double precision :: f(3)
integer tmpn, n, nread
character(len=*) ReadingWhat ! for a clearer error message
integer :: i, ios
integer, external :: lnblk

!-init
nread=0


if (iecho.eq.8) then
    call ParallelLog("FGEN", "Entered FGEN")
endif

if (M.gt.3) then
    stop 'FGEN overflow'
endif

if (debug) then
#ifdef SPARSE
    call parallelLog("fgen","start reading " // trim(ReadingWhat))
#else
    write(*,*) "rank", getrank(), " fgen READING ", trim(ReadingWhat)
#endif
endif

100     call readrc (luin,record,'FGEN',end)

#ifdef SPARSE
!call mpi_barrier(MPI_COMM_WORLD,petscierr)
#endif


if (debug) then
    write(*,*) "rank", getrank(), 'FGEN ' ,trim(ReadingWhat), ' read record: ', trim(record)
endif

if (end) then
    if (debug) then
        write(*,*) "rank", getrank(), "FGEN ' ,trim(ReadingWhat), '  detected an ''end'' statement!"
    endif
    goto 200
endif

if (debug) then
!    write(*,*) '--------------------'
!    write(*,*) 'rank', getrank(), 'FGEN read record: ', trim(record)
    write(*,*) 'rank', getrank(), 'FGEN  ' ,trim(ReadingWhat), &
               ' wants to split it with FMT', FMT, ' and M ', M
endif

read (record,FMT,err=300,end=400,iostat=ios) n,(f(i),i=1,M)

if (debug) then
    write(*,*) 'rank', getrank(), 'FGEN  ' ,trim(ReadingWhat), &
               ' resulted in ', n,(f(i),i=1,M)
!    write(*,*) '--------------------'
endif

! check whether the node is in this parition.
! If it is not, skip it.
#ifdef SPARSE
if (n .gt. meshdatactx%nvglobal) then
	write(*,*) "rank ", getrank(), "detects error! While reading ",trim(ReadingWhat), &
               " , FGEN want to check the vertix2partition array on index n= ", n, &
               "While the mesh only has ", meshdatactx%nvglobal, "points!"
    write(*,*) "based on record: ", trim(record)
    write(*,*) "interpreted as: ", n,f(1:M)
endif



if (meshdatactx%v2p(n) .ne. getrank()) then
    goto 100
endif
!if (meshdatactx%v2p(i) .ne. getrank()) then
!    goto 200
!endif


#endif





if (n.eq.0) then
!    write(*,*) 'rank', getrank(), 'node id ', n, 'wrong!'
    goto 200
endif

nread = nread + 1
tmpn = n

#ifdef SPARSE
if (getsize().eq.1) then
    n = tmpn
else
    n = verticesmask(meshdatactx%Glo2LocVertex(tmpn))
endif
#else
n = tmpn
#endif


! experiment part 2
if (n.eq.-1) then
    goto 100
endif


!        n = nread

!if (nread.gt.MAX) then
!    write(stderr,2) MAX
! 2       format(1x,'FGEN: more than ',I12,' entries')
!    write(stderr,*) 'something wrong with ', callerMessage
!    call xit(1," ")
!endif

if (n.lt.1 .or. n.gt.MAX) then
    write(stderr,*) 'FGEN: entry ',n,'outside range [1,',MAX,']'
    write(stderr,*) 'something wrong with ', ReadingWhat
    call xit(1,"FGEN; entry out of range")
endif

!if (iecho.eq.8) then
!    call ParallelLog("FGEN", "writing IA")
!endif


do i=1,m
!    if (i.gt.size(A,1) .or. n.gt.size(A,2)) then
!        write(*,*) "FGEN writing A", i, n, &
!                   "of ", size(A,1), size(A,2), &
!                   "while reading", trim(ReadThis)
!        call delay()
!       stop "Aborting"
!    else
        A(i,n) = f(i)
!    endif

!    if (i.gt.size(f,1)) then    
!        write(*,*) "FGEN reading f", i, &
!                   "of ", size(f,1), &
!                   "while reading", trim(ReadThis)
!        stop "Aborting"
!    endif


!    A(i,n) = f(i)
enddo

goto 100


200 if (iecho.eq.8) then
    call ParallelLog("FGEN", "Leaving FGEN")
endif

#ifdef SPARSE
if (debug) then
    call parallelLog("fgen","finished reading " // trim(ReadingWhat))
endif
#endif

flush(stdout)


return

300     write(stderr,5) ios
 5   format(1x,'FGEN: read error ',i4)
write(stderr,*) 'something wrong with ', ReadingWhat

write(stderr,6) n,(f(i),i=1,m)
 6   format(1x,'FGEN: previous record:'/i12,7g14.6)

write(stderr,7) '"',record(1:lnblk(record)),'"'
 7   format(1x,'FGEN: current record: ',110a)

write(stderr,*) 'tried to interpret with format: ', fmt

call prerr('FGEN')
call xit(1,"Error in FGEN")

400     write(stderr,8) ios
 8   format(1x,'FGEN: premature EOF error ',I4)
write(stderr,5) n,(f(i),i=1,m)
call prerr('FGEN')
call xit(1,"Unexpected EOF in FGEN")

end
!-------------------------------------------------------------------------------
subroutine IGEN (IA,M,MAX,FMT,NREAD,flag,ReadingWhat)

! sub to read nodal boundary condition codes.
! The selection whether or not a BC is relevant
! for this partition, is made here as well.

USE MODELCTX,        only : getrank, getsize
USE MODELDATAMODULE, only : modeldatactx
USE MESHDATAMODULE,  only : meshdatactx
USE MODELDEFINITION
use debugmodule,     only : xit, debug, &
                            parallellog
#ifdef SPARSE
use aomodule,        only : vertices, &
                            verticesmask
#endif
use iomodule          
!
!       FLAG=1: put ibc into global array for temporary use

implicit none
!-pass
integer M,MAX,NREAD, flag  ! M = number of integers per line, after the index
                           ! MAX = maximum nuber of lines 
                           ! nread = read number of lines
                           ! flag = 1 for mechanical nodal; 
                           !        2 for thermal nodal
                           !        0 for rest such as (differential) Winkler forces
character :: FMT*(*)
!dimension IA(M,*)  ! IA = integer array
integer   :: IA(M,MAX)
character(len=*) :: ReadingWhat


!-locl
logical end
integer :: i, ib(3), ios
character(len=1024)record
!        integer getrank, tmpn
integer tmpn, n, nThermalBC
!-init
nread = 0
ib = 0
nThermalBC = 0

!write(*,*) "rank",getrank(),"IGEN reading ", trim(ReadingWhat), &
!           " with length", M, MAX, flag

if (M.gt.3) then
    stop 'IGEN overflow; M larger than three.'
endif

if (debug) then
#ifdef SPARSE
    call parallelLog("igen","start reading " // trim(ReadingWhat))
#else
    write(*,*) "rank", getrank(), "IGEN READING ", ReadingWhat
#endif
endif


100     call readrc (luin,record,'IGEN',end)

if (debug) then
    write(*,*) 'rank', getrank(), 'IGEN ' ,trim(ReadingWhat), ' read record: ', trim(record)
endif

if (end) then
    if (debug) then
        write(*,*) "rank", getrank(), "IGEN " ,trim(ReadingWhat), "  detected an ''end'' statement!"
    endif
    goto 200
endif

read(record,fmt,err=300,end=310,iostat=ios) n,(ib(i),i=1,M)

! n is the nodal point, ib are the values
if (n.gt.meshdatactx%nvglobal) then
    write(*,*) "Encountered nodal boundary condition of node", n
    write(*,*) "But the mesh only has", meshdatactx%nvglobal, "nodal points"
    stop "Exiting..."
endif

! some global arrays needed for establishing equation numbers.
! Nodal boundary condition can restrict the number of degrees of freedom
! And currently that needs to en known globally.

if (flag.eq.1) then
    ! mechanical boundary conditions
    do i=1,M
        modeldatactx%IBONDGLOBAL(i,n) = ib(i)
    enddo
else if (flag.eq.2) then
    ! thermal boundary conditions;
    do i=1,M  ! usually, M=1 (are there cases with M>1 ?)
        modeldatactx%ITBCGLOBAL(n) = ib(i)

        if (ib(i) .eq. 1) then
            nThermalBC = nThermalBC + 1
        endif
    enddo


    if (nThermalBC .eq. meshdatactx%nvglobal .and. &
              mode .eq. 5) then
        ! mode 5 = only thermal
        write(stderr, '(a)') "ERROR: The temperature of every nodal point has been"
        write(stderr, '(a)') "fixated by a thermal boundary condition."
        write(stderr, '(a)') "Did you perhaps mean to set initial temperatures to each"
        write(stderr, '(a)') "nodal point in stead of fixating it?" 
        write(stderr, '(a)') "Then an 'end' statement between the thermal boundary"
        write(stderr, '(a)') "conditions and the initial temperatures may be missing."
        write(stderr, '(a)') "See 'man 7 pln' for details on the input format."
        stop "Leaving GTECTON..."
    endif

else
    ! flag.eq.0  = (differential) Winkler forces
endif
! the arrays below here are only filled locally, so
! a partition selection is made

if (meshdatactx%v2p(n).ne.getrank()) then
!    write(*,*) 'rank ', getrank(),' read BCS value ', n,(ib(i),i=1,M), 'value ignored'
!    write(*,*) 'but verticesmask: ', verticesmask(n), vertices(n)
else
!    write(*,*) 'rank ', getrank(),' read BCS value ', n,(ib(i),i=1,M), 'value added'
!    write(*,*) 'but verticesmask: ', verticesmask(n), vertices(n)
endif


! experiment. Do not read only bcs from this partition, but of all partitions!
!write(*,*) " Glo2loc: ", meshdatactx%Glo2LocVertex(n), " but v2p: ", meshdatactx%v2p(n)

#ifdef SPARSE
if (verticesmask(meshdatactx%Glo2LocVertex(n)).eq.0) then
    ! vertices mask contains a nonzero for every 
    ! relevant vertex for this partition.
    ! If it is zero, skip storing this entry.
    goto 100
endif
#endif

if (n.eq.0) then
    goto 200
endif
!
nread = nread + 1
tmpn = n

#ifdef SPARSE
if (getsize().eq.1) then
    n = tmpn
else
    n = verticesmask(meshdatactx%Glo2LocVertex(tmpn))
!    write(*,*) "rank",getrank(),"global vertex",tmpn,"is local",n
!    n = n - meshdatactx%lowestLocalVertexID
!    write(*,*) "rank",getrank(),"scales down to", n
endif
#else
n = tmpn
#endif


if (n.lt.1 .or. n.gt.MAX) then
    write(stderr,3) n, MAX
 3       format(1x,'IGEN:',I12,' entry outside range [1,',I12,']')
    write(stderr,*) "while reading: ", trim(ReadingWhat)
    call xit(1,"IGEN outside range")
endif 


do i=1,M ! M is number of variables per point
       IA(i,n) = ib(i)
enddo



!if (flag.eq.10) then
!    write(*,*) "Normal Winkler Forces"
!    do i=1,M ! M is number of variables per point
!        IA(i,n) = ib(i)
!        modeldatactx%iwink(i,n) = ib(i)
!    enddo
!else if (flag.eq.11) then
!    write(*,*) "Differential Winkler Forces"
!    do i=1,M ! M is number of variables per point
!        IA(i,n) = ib(i)
!        write(*,*) "rank",getrank(),"writing iwinx", &
!                   size(modeldatactx%iwinx,1), &
!                   size(modeldatactx%iwinx,2), i, n
!        if (n .gt. size(modeldatactx%iwinx,2)) then
!            write(*,*) "rank",getrank(),"fuckupski ", n, size(modeldatactx%iwinx,2)
!        endif
!        modeldatactx%iwinx(i,n) = ib(i)
!        write(*,*) "rank",getrank(),"writing iwinx done"
!    enddo
!endif




goto 100

200     i=1
#ifdef SPARSE
if (debug) then
    call parallelLog("igen","finished reading " // trim(ReadingWhat))
endif
#endif

flush(stdout)

return


300 write(stderr,*) 'rank', getrank(), 'IGEN ' ,trim(ReadingWhat), " read error", ios

!     write(stderr,5) ios
!write(stderr,*) 'rank', getrank(), 'IGEN ' ,trim(ReadingWhat), " read error", ios

! 5   format(1x,'IGEN: read error ',i5)
!write(stderr,6) n,(ib(i),i=1,M)
! 6   format(1x,'IGEN: previous record: ',i12,7g14.6)
!call prerr('IGEN')
!call xit(1,"IGEN; read error")

310 write(stderr,*) 'rank', getrank(), 'IGEN ' ,trim(ReadingWhat), " premature end reached", ios

!     write(stderr,15) ios
!write(stderr,*) 'rank', getrank(), 'IGEN ' ,trim(ReadingWhat), " premature end reached", ios

!15   format(1x,'IGEN: read error; end reached ',i5)
!write(stderr,16) n,(ib(i),i=1,M)    
!16   format(1x,'IGEN: previous record: ',i12,7g14.6)
!call prerr('IGEN')
!call xit(1,"IGEN; read error")



end
!-------------------------------------------------------------------------------


subroutine LOCAL (ID,IEN,LM,LINK,NDOF,NEN,NUMEL,NUMNP,NLINK)

use debugmodule,     only: xit
use iomodule,        only: stderr
! replaced by code in createparallel../
!
! Subroutine to localize ID array, i.e. create LM from ID
!
implicit none

!-pass
integer :: NDOF,NUMNP,NEN,NUMEL,NLINK
integer :: ID(NDOF,NUMNP),IEN(NEN,NUMEL),LM(NDOF,NEN,NUMEL)
integer :: LINK(3,*)
!-locl
integer :: i,j,k,n,nn,master
!
do n=1,NLINK
    i = LINK(1,n)
    master = LINK(3,n)
    k = ID(i,master)

    if (k.gt.0) then
        ID(i,LINK(2,n)) = k
    else
        write(stderr,*) 'BC: linked node master ',master, &
         'generates no equation'
        call xit(1,"Linked node master generates no equation")
    endif

enddo
!
do k=1,NUMEL
    do j=1,NEN
        nn = IEN(j,k)
        do i=1,NDOF
               LM(i,j,k) = ID(i,nn)
        enddo
    enddo
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine LOCALF (NFAULT,IEN,NDOF,NUMFN,NUMEL,NEN)
 
USE MESHDATAMODULE,  only: meshdatactx
USE ALGEBRA,         only: clear
use modeldatamodule, only: modeldatactx
use modelctx,        only: getsize
use debugmodule,     only: xit
use iomodule,        only: stderr
!       this routine goes per processor
!
!       Localize NFAULT array for efficient computation

implicit none

!-pass
integer :: NDOF,NUMFN,NUMEL,NEN
integer :: NFAULT(3,*),IEN(NEN,*)
!-locl
integer :: i,j,k,nelm,nloc,node
logical :: NoMatch

if (NUMFN.le.0) then
    return
endif
 
call CLEAR (modeldatactx%LMF,NDOF*NEN*meshdatactx%Nelocal,"modeldatactx%LMF") 
 
iloop: do i=1,NUMFN

!   nelm and node are still numbered globally
    nelm = IABS(modeldatactx%NFAULT(1,i))

    ! find local element number... (what happens on partition boundaries?)
#ifdef SPARSE
    if (getsize().eq.1) then
        nloc = nelm
    else
        nloc = meshdatactx%Glo2LocElement(nelm)
    endif
#else
    nloc = nelm
#endif

    if (nloc.le.0) then
        cycle iloop
    endif

    node = modeldatactx%NFAULT(2,i)
    NoMatch = .true.
    do j=1,NEN
        if (IEN(j,nloc).eq.node) then
            NoMatch = .false.
            do k=1,NDOF
                modeldatactx%LMF(k,j,nloc) = i  ! entry of NFAULT
            enddo
        endif
    enddo

    if (NoMatch) then ! should never occur, we verified this in RFAULT
        write(stderr,1) node, nelm 
 1      format(1x,'LOCALF: node ',I12,' doesn''t occur in',1x, &
         'element ',I12)
        call xit(1,"LOCALF")
    endif

enddo iloop

return
end subroutine
!-------------------------------------------------------------------------------
subroutine LOCALX (IDXglobal, NSLIP, NUMSLP)

USE MESHDATAMODULE
USE ALGEBRA,         only: clear
use modeltopology,   only: ndof, nen
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: HasSlipperyNodes, &
                           modeldatactx, &
                           elementSlipperyID
#else
use modeldatamodule, only: HasSlipperyNodes, &
                           modeldatactx
#endif
use meshdatamodule,  only: meshdatactx
use modelctx,        only: getsize, &
                           getrank
use iomodule          
use debugmodule,     only: debug, xit
!
!       Localize IDXglobal and transfer sign information to LMX

implicit none

!-pass
integer :: NUMSLP
integer :: IDXglobal, nslip
dimension IDXglobal(NDOF,*),NSLIP(5,NUMSLP)
!-locl
integer i,j,nel,node,idof
integer neloc
!-init
if (.not. HasSlipperyNodes()) then
    return
endif

if (debug) then
    write(*,*) "rank",getrank(),"in localx, IEN:", meshdatactx%IEN 
    write(*,*) "rank",getrank(),"in localx, numslip",NUMSLP," and NSLIP", NSLIP
endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
modeldatactx%LMX = 0
#else
call CLEAR(modeldatactx%LMX,NDOF*NEN*meshdatactx%Nelocal,"modeldatactx%LMX")
#endif

! for every slippery node entry,
! check whether the node-element combinations are valid.

! NSLIP has already been stored by partition, while reading them.
! so that all NSLIP entries are relevant for this partition.
do i=1,NUMSLP

    nel  = NSLIP(1,i) ! global element nr.
    node = NSLIP(2,i) ! global node number


#ifdef SPARSE
    
    if (getsize().eq.1) then
        neloc = nel
    else
!?      neloc = meshdatactx%Glo2LocElement(nel) - meshdatactx%LocElt(1)
!?      neloc = nel
        neloc = meshdatactx%Glo2LocElement(nel)
    endif
#else
    neloc = nel
#endif

    if (debug) then
        write(*,*) 'tecin/localx rank', getrank(), 'global element# ', nel, ' is local# ',neloc
    endif

    do j=1,NEN
        ! Loop through the nodes and check whether it is in the element in question.
        ! IEN contains global vertex numbers but numbering of the elements according
        ! to what it found for its partitioin in tecin.dat.partf.elm.
        ! But slippery node entries contain global on global
        if (node.EQ.meshdatactx%IEN(j,neloc)) then
            goto 100
        endif
    enddo
    write(stderr,1) i,node,nel
 1  format(1x,'LOCALX: NSLIP(',I12,'): node=',I12,' does not occur in element ',I12)
    call xit(1,"Slippery nodes incorrect")


100 continue

!   Global node number NODE is in the same partition as the element with
!   global element number NEL and partition-local element number NELOC
    if (debug) write(*,*) 'rank', getrank(), "is convinced that node", node, "is in elem", nel

    do idof=1,NDOF
        if (debug) write(*,*) 'tecin/localx: IDX(idof=',idof,',node=',node,')=',IDXglobal(idof,node)

#ifdef EXPERIMENTAL_ReducedSlipAlloc

! test this!
        if (debug) write(*,*) "test this, accessing LMX ", elementSlipperyID(nel,102), &
         "of", meshdatactx%Nelocal

        modeldatactx%LMX(idof,j,elementSlipperyID(nel,102)) = ISIGN(IDXglobal(idof,node),NSLIP(2+idof,i))

        if (debug) write(*,*) "write LMX", idof, j, elementSlipperyID(nel,102), "to be", &
          idof, node, IDXglobal(idof,node), i, NSLIP(2+idof,i)

#else
        modeldatactx%LMX(idof,j,neloc) = ISIGN(IDXglobal(idof,node),NSLIP(2+idof,i))
        if (debug) then
            write(*,*) 'tecin/localx: IDX(idof=',idof,',node=',node,')=',IDXglobal(idof,node)
            if (IDXglobal(idof,node).ne.0) then
                write(*,*) 'tecin/localx: IEN(',j,',neloc=',neloc,')=',meshdatactx%IEN(j,neloc)
                write(*,*) 'tecin/localx: LMX(',idof,',',j,', neloc=',neloc,')=',modeldatactx%LMX(idof,j,neloc)
            endif
        endif
#endif

    enddo
enddo

return
end
!-------------------------------------------------------------------------------
subroutine SKEWRD (NSD)

USE MESHDATAMODULE,  only: meshdatactx
USE MODELDATAMODULE, only: modeldatactx
USE MODELDEFINITION, only: NUMROT
!USE ALGEBRA
USE CONSTANTS,       only: zero, pi, deg2rad
use debugmodule,     only: iecho, parallellog, Verified_MPI_Barrier
use modelctx,        only: getrank
use iomodule          
!
! Program to read and print the array SKEW that specifies local
! coordinate rotations at each node. Used for applying boundary
! conditions in directions other than the global coordinate system.
!
implicit none
!-pass
integer :: NSD
!-locl
character(len=13) fmt
integer :: nn, i, j, idum
!-init

!call Verified_MPI_Barrier("starting to read euler angles")

fmt = '(I12,3G14.0)'
!
if (iecho.eq.2) then
  write(stderr,1) fmt
 1   format(1x,'*** read rotated boundary condition data ***'/ &
 1x,'n=1; while (n != 0){'/ &
 5x,'read ',A12,' n,(SKEW(i,n),i=1,2)'/ &
 1x,'}')
endif


! skew is being read per processor
call FGEN(modeldatactx%SKEW,2,meshdatactx%nvglobal,fmt,idum,"Euler angles")

!call Verified_MPI_Barrier("done FGEN")


!NUMROT = 0

!do i=1,meshdatactx%Nvlocal
!    if ((modeldatactx%SKEW(1,i).ne.ZERO).or. &
!        (modeldatactx%SKEW(2,i).ne.ZERO)) then
!!        write(*,*) "rank", getrank(), "found SKEW", modeldatactx%SKEW(1,i), modeldatactx%SKEW(2,i)
!        NUMROT = NUMROT + 1
!    endif
!enddo

NUMROT = idum

if (NUMROT.eq.0) then
    return
endif
!
if (OUTFIL(1)) then
    nn = 0
    do i=1,meshdatactx%Nvlocal!NUMNP
        if ((modeldatactx%SKEW(1,i).ne.ZERO).or. &
            (modeldatactx%SKEW(2,i).ne.ZERO)) then

            if (MOD(nn,50).eq.0) then
                if (NSD.eq.2) then
                    write(luout(1),2)
 2                   format(1X///' L O C A L',4x, &
                 'C O O R D I N A T E',4x,'R O T A T I O N S'// &
                 5X,'ALPHA = CC ROTATION (DEGREES) FROM X AXIS', &
                 ' IN XY PLANE'//3X,'THIS ROTATES THE SKEW', &
                 ' COORDINATES INTO THE GLOBAL SYSTEM'//8X, &
                 'NODE #',10X,'ALPHA'/)
                else
                write(luout(1),3)
 3                   format(1X///' L O C A L',4x, &
                 'C O O R D I N A T E',4x,'R O T A T I O N S'// &
                 5X,'ALPHA = CC ROTATION (DEGREES) FROM X AXIS', &
                 ' IN XY PLANE'/6X,'BETA = CC ROTATION',1x, &
                 '(DEGREES) FROM X AXIS IN XZ PLANE'//3X, &
                 'ALPHA ROTATION IS APPLIED FIRST, THEN THE',1x, &
                 'BETA ROTATION'/3X,'THIS ROTATES THE SKEW',1x, &
                 'COORDINATES INTO THE GLOBAL SYSTEM'//8X, &
                 'NODE #',10X,'ALPHA',11X,'BETA'/)
                endif
            endif
            nn = nn + 1

            write(luout(1),4) i,(modeldatactx%SKEW(j,i),j=1,NSD-1)
 4                format(2X,I12,10X,2(F15.10))

        endif
    enddo
endif

!do i=1,meshdatactx%Nvlocal!,NUMNP
do i=1,meshdatactx%Nvglobal!,NUMNP
    if ((modeldatactx%SKEW(1,i).ne.ZERO).or.&
        (modeldatactx%SKEW(2,i).ne.ZERO)) then
         modeldatactx%SKEW(1,i) = deg2rad * modeldatactx%SKEW(1,i)
         modeldatactx%SKEW(2,i) = deg2rad * modeldatactx%SKEW(2,i)
    endif
enddo


return
end
!-------------------------------------------------------------------------------
subroutine CHKFAULTSLIP(NFAULT,NSLIP,NUMFN,NUMSLP)

! verifies that the same element/node combination does not appear in faulted and
! slippery node arrays

use debugmodule,     only : xit
use iomodule          

implicit none

!-pass
integer :: NFAULT,NSLIP,NUMFN,NUMSLP
dimension NFAULT(3,*),NSLIP(5,*)
!-locl
logical :: multiple
integer :: if,felem,fnode,is,selem,snode

if (NUMFN*NUMSLP.eq.0) then
    return
endif

multiple = .false.
do if=1,NUMFN
    felem = NFAULT(1,if)
    fnode = NFAULT(2,if)
    do is=1,NUMSLP
        selem = NSLIP(1,is)
        snode = NSLIP(2,is)
        if (selem.eq.felem .and. snode.eq.fnode) then
            multiple = .true.
            write(stderr,10) fnode,felem
10          format('*** Node ',I12,' in element ',I12,' is listed as both ', &
             'faulted and slippery ***')
        endif
    enddo
enddo
if (multiple) then
    call xit(1,'***** fatal dual node/elements ****')
endif

return

end
!-------------------------------------------------------------------------------
subroutine RFAULT(NFAULT,FAULT,IEN,NEN,NDOF)

USE MODELCTX,        only : getrank, getsize
USE MESHDATAMODULE
USE MODELDATAMODULE, only : modeldatactx
USE MODELDEFINITION
use formatsmodule,   only : FMT_tecin_faulted_nodes
use debugmodule,     only : iecho, xit, debug
use timestepmodule,  only : faultvelo
use iomodule          

! Reads, verifies and prints data on faulted nodes
!
!    logical VELBC is used to change interpretation of the data in FAULT
!    VELBC = .true.: FAULT means constant velocity (NFAULT(3,.) is not used)
!
!    VELBC = .false.: FAULT means differential displacement
!                     NFAULT(3,.) >= 0 differential displacement applied @ time step NFAULT(3,.)
!                     NFAULT(3,.)  < 0 periodic displacement applied every ABS(NFAULT(3,.) time step 

implicit none

!-pass
integer NDOF,NEN
integer   :: IEN
double precision :: fault
integer :: nFault
dimension NFAULT(3,*),FAULT(NDOF,*),IEN(NEN,*)
!-locl
logical end,VELBC
character(len=71) record ! 12 + 12 + 5 + 3 * 14
integer :: ia(3)
double precision :: f(3)
integer, external :: lnblk
integer :: neloc, nread, nn, i, j
integer :: ios
integer :: lnod
integer :: ielm, inod, nelm, itim

#ifdef DEBUG
write(*,*) "entered rfault"
#endif

if (iecho.eq.2) then
    write(stderr,1) FMT_tecin_faulted_nodes
 1  format(1x,'if (NUMFN > 0){'/ &
     5x,'*** read faulted nodes data ***'/ &
     5x,'while (NFAULT(1,i)!=0 & NFAULT(1,i)!=end) {'/ &
     9x,'read ',A17,1x, &
     'NFAULT(1,i),NFAULT(2,i),NFAULT(3,i),(FAULT(j,i),',1x, &
     'j=1,NDOF)'/5x,'}'/1x,'}')
endif
if (NUMFNglobal.le.0) then
    return
endif
 
VELBC = (faultvelo.ne.0)
nread = 0
NUMFN = 0

100 call readrc (luin,record,'RFAULT',end)

    if (debug) then
        write(*,*) 'read fault record: ', record
    endif

    if (end) then
#ifdef DEBUG
        write(*,*) "end found"
#endif
        goto 200
    endif

    read(record,FMT_tecin_faulted_nodes,err=500,end=600,iostat=ios) (ia(i),i=1,3), (f(i),i=1,NDOF)
    if (debug) then
        write(*,*) 'interpreted: ', (ia(i),i=1,3), (f(i),i=1,NDOF)
    endif

    if (ia(1).eq.0) then
#ifdef DEBUG
        write(*,*) "found end"
#endif
        goto 200
    endif

    nread = nread + 1
    if (nread.gt.NUMFNglobal) then
        write(stderr,3) NUMFNglobal
 3      format(1x,'RFAULT: more than the expected ',I12, 'entries in input list')
        call xit(1,"Too many fault entries")
    endif
    if (ia(1).lt.1 .or. ia(1).gt.meshdatactx%neglobal) then
        write(stderr,4) meshdatactx%neglobal
 4      format(1x,'RFAULT: element outside range [1,',I12,']')
        call xit(1,"Fault element ID outside range")
    endif
    if (ia(2).lt.1 .or. ia(2).gt.meshdatactx%nvglobal) then
        write(stderr,5) meshdatactx%nvglobal
 5      format(1x,'RFAULT: node outside range [1,',I12,']')
        call xit(1,"Fault vertex ID outside range")
    endif
    modeldatactx%NFAULTglobal(1,nread) = ia(1)  ! element in which split element occurs
    modeldatactx%NFAULTglobal(2,nread) = ia(2)  ! split node number
    if (VELBC) then
        modeldatactx%NFAULTglobal(3,nread) = 0  ! ignore input
    else
        modeldatactx%NFAULTglobal(3,nread) = ia(3)  ! onset time
    endif

#ifdef SPARSE
    if (meshdatactx%e2p(ia(1)).ne.getrank()) then
#ifdef DEBUG
        write(*,*) "tried to match rank: ", meshdatactx%e2p(ia(1)), "to", getrank()
        write(*,*) "rank not matching, skipping"
#endif
        ! this faulted node does not belong to this partition
        goto 100
    endif
#endif

    NUMFN = NUMFN + 1   ! increment faulted node entries for this partition
    NFAULT(1,NUMFN) = ia(1)
    NFAULT(2,NUMFN) = ia(2)
    if (VELBC) then
        NFAULT(3,NUMFN) = 0
    else
        NFAULT(3,NUMFN) = ia(3)
    endif

    do i=1,NDOF
        FAULT(i,NUMFN) = FAULT(i,NUMFN) + f(i)
    enddo

goto 100


200 if (nread.ne.NUMFNglobal) then
        write(stderr,6) nread,NUMFNglobal
6       format(1x,'RFAULT: read ',I12,' data, expected ',I12)
        call xit(1,"Not enough fault entries")
endif

nn=0
if (OUTFIL(1)) then
    do i=1,NUMFN
        if (MOD(nn,50).eq.0) then
            if (VELBC) then
                write(luout(1),7)
7               format(1X///2X,'F A U L T E D   N O D E   D A T A'// &
                 4X,'NODE',9X,'ELEM',4X,'DUMMY',8X,'SPLIT',12X,'SPLIT',12X,'SPLIT'/ &
                36X,'VELO-DOF1',8X,'VELO-DOF2',8X,'VELO-DOF3'//)
            else
                write(luout(1),8)
8              format(1X///2X,'F A U L T E D   N O D E   D A T A'// &
                6X,'KEY TO APPLICATION MODE CODES:'// &
                10X,'MODE.GE.0  = TIME OF DISPLACEMENT'/ &
                10X,'MODE.LT.0  = DISPLACEMENT EVERY -MODE TIME STEPS'// &
                4X,' NODE',3X,'ELEM',5X,'MODE',9X,'SPLIT',12X,'SPLIT',12X,'SPLIT'/ &
                32X,'DISP-DOF1',8X,'DISP-DOF2',8X,'DISP-DOF3'//)
            endif
        endif
        nn = nn + 1
        write(luout(1),9) NFAULT(2,i),NFAULT(1,i),NFAULT(3,i),(FAULT(j,i),j=1,NDOF)
9       format(I8,5X,I8,5X,I4,3(5X,1PG12.5))
    enddo
endif

!write(*,*) 'numfn: ', numfn, numfnglobal

! Check if faulted nodes actually occur in specified elements
nelm = 0
outer: do i=1,NUMFN
    nelm = NFAULT(1,i)
!    write(*,*) 'nelm: ',i, nelm

    if (nelm.lt.1 .or. nelm.gt.meshdatactx%Neglobal) then
        write(stderr,12) i
12      format(1x,'NFAULT entry ',I12,': element number out',1x, &
         'of range [1,NUMEL]')
        call xit(1,"Faulted nodes, element ID out of range ")
    endif

    lnod = 1
    do while (lnod.le.NEN) 

#ifdef SPARSE
        if (getsize().eq.1) then
            neloc = nelm
        else
            neloc = meshdatactx%Glo2LocElement(nelm)
        endif
#else
        neloc = nelm
#endif

        if (NFAULT(2,i).eq.IEN(lnod,neloc)) then
            cycle outer
        endif
        lnod = lnod + 1

    enddo

enddo outer

if (numfn.gt.0) then
    if (lnod.gt.NEN) then
        write(stderr,13) i,NFAULT(2,i),nelm
13  format(1x,'NFAULT entry ',I12,': node ',I12,' does not', &
         ' occur in element ',I12)
        call xit(1,"Faulted nodes, vertex ID out or range")
    endif
endif

! Verify that multiple element/node entries do not occur and that timing is
! correct
do i=1,NUMFN-1
   ielm = NFAULT(1,i)
   inod = NFAULT(2,i)
   itim = NFAULT(3,i)
   do j=i+1,NUMFN
        if (NFAULT(1,j).eq.ielm .and. NFAULT(2,j).eq.inod) then
            write(stderr,*) "Every element/node combination may occur only once"
            write(stderr,*) "in the list of faulted nodes."
            write(stderr,*) "Line", i, "and line", j, "appear to have the same combination."

            write(stderr,14) i,j
14          format(1x,'RFAULT: NFAULT entries ',I12,' and ',I12, &
             ' refer to the same element/node combination')
            call xit(1," ")
        endif
        if (NFAULT(2,j).eq.inod .and. NFAULT(3,j).ne.itim) then
            write(stderr,15) inod,i,itim,j,NFAULT(3,j)
15          format(1x,'Faulted node timing discrepancy for node ',I12,/ &
             'NFAULT(3,',I6,') = ',I9,' and NFAULT(3,',I6,') = ',I9)
            call xit(1,"RFAULT")
        endif
   enddo
enddo

return
!
500 write(stderr,501) ios
501 format(1x,'RFAULT: input error ',i4)
call prerr('RFAULT')
call xit(1,"Error when reading faults")
600 write(stderr,601) ios
601 format(1x,'RFAULT: premature EOF error ',i4)
call prerr('RFAULT')
call xit(1,"Unexpected EOF when reading faults")
end
!-------------------------------------------------------------------------------
subroutine RSLIP (DIFORC,IEN,NEN,NDOF,NSD)
!
! Reads and prints data on free slip interfaces. Every processor read all because of renumbering needed later
! (if each processor would read only his part, the global renumbering would fail).
! NSLIPglobal is a global array needed for later renumbering of equations, NSLIP is the processor-local array.
!
use modelctx,        only : getrank, getsize
use meshdatamodule,  only : meshdatactx
use modeldatamodule, only : modeldatactx
use modeldefinition
use formatsmodule,   only : fmt_tecin_slippery_nodes
use debugmodule,     only : iecho,debug,xit,itest
use algebra,         only : idig
use iomodule


implicit none
!-pass
integer   :: NDOF,NEN
integer   :: IEN
dimension IEN(NEN,*)
integer :: NSD
#ifdef SPARSE
!double precision :: DIFORC(NDOF,meshdatactx%nvlocal)
double precision :: DIFORC(NDOF,meshdatactx%nvglobal)
#else
double precision :: DIFORc
dimension DIFORC(NDOF,*)
#endif

!-locl
logical end
character(len=256) record
double precision :: f(3)
integer :: ia(5)
integer :: nploc, nread
integer :: i, j, k, nn
integer :: ios
integer :: localSequence
integer :: iSlipEntry
integer :: nodeID
integer :: weightx,weighty
logical :: nodesAreConsistent, firstInconsistentNode

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer :: iElem
logical, allocatable :: elemHasSlipNode(:)
allocate(elemHasSlipNode(meshdatactx%neglobal))
elemHasSlipNode = .false.
#endif

nread=0

if (iecho.eq.2) then
    write(stderr,1) FMT_tecin_slippery_nodes
1   format(1x,'if (NUMSLP > 0) {'/ &
           5x,'*** read slippery nodes data ***'/ &
           5x,'while (NSLIP(1,i)!=0 & NSLIP(1,i)!=end) {'/ &
           9x,'read ',A18,' (NSLIP(j,i),j=1,5),', &
           '(DIFORC(j,NSLIP(2,i)),j=1,NDOF)'/5x,'}'/1x,'}')
endif

if (NUMSLPglobal.le.0) then
    if (iecho.eq.2) then
        write(*,*) "No slippery nodes configured."
    endif
    return
endif

NUMSLP  = 0
NUMSE   = 0
NUMSN   = 0

#ifdef EXPERIMENTAL_ReducedSlipAlloc
nLocalElemsWithSLipperyNodes = 0
nElemsWithSlipperyNodes = 0
#endif

read_and_process: do

    call readrc (luin,record,'RSLIP',end)

    if (debug) then
        write(*,*) "rank ", getrank(), "rslip read record: ", trim(record)
    endif

    if (end) then
        if (debug) then
            write(*,*) "rank ", getrank(), "rslip found the end of slippery node list"
        endif
        exit read_and_process
    endif

    if (debug) then
        write(*,*) "rank ", getrank(), "interpreting with format ", FMT_tecin_slippery_nodes
    endif

    ! ia goes from 1 to 5, (elemID, vertexID, 3 weights... even if only 2 dimension, still 3 weights.
    read(record,FMT_tecin_slippery_nodes,err=300,end=400,iostat=ios) (ia(i),i=1,5), (f(i),i=1,NDOF) 

    if (debug) then
        write(*,*) "rank ", getrank(), "interpreted as ", ia(1:5), f(1:NDOF)
    endif

    if (ia(1).eq.0) then
        exit read_and_process
    endif

    nread = nread + 1

    if (nread.gt.NUMSLPglobal) then
        write(stderr,3) NUMSLPglobal
3       format(1x,'RSLIP: more than the expected ',I8, ' entries in input list')
        call xit(1,"Too many slippery nodes")
    endif

    if (ia(1).lt.1) then
        write(*,*) "Read slippery node element", ia(1), &
         "but element numbering starts at 1."
    endif

    if (ia(1).gt.meshdatactx%neglobal) then
        write(*,*) "Read slippery node element", ia(1), & 
         "but mesh contains only", meshdatactx%neglobal, "elements"
    endif

    if (ia(1).lt.1 .or. ia(1).gt.meshdatactx%neglobal) then
        write(stderr,4) meshdatactx%neglobal
4       format(1x,'RSLIP: element outside range [1,',I12,']')
        call xit(1,"Element number out of range")
    endif

    if (ia(2).lt.1) then
        write(*,*) "Read slippery node ", ia(2), & 
         "but point numbering starts at 1."
    endif

    if (ia(2).gt.meshdatactx%nvglobal) then
        write(*,*) "Read slippery node ", ia(2), &
         "but mesh contains only", meshdatactx%nvglobal, "nodal points"
    endif

    if (ia(2).lt.1 .or. ia(2).gt.meshdatactx%nvglobal) then
        write(stderr,*) "RSLIP: slippery node", ia(2), &
         "outside range [1", meshdatactx%nvglobal,"]" 
        call xit(1,"Vertex number out of range")
    endif

    ! store all slippery nodes in a global array
    ! TODO: Is this necessary? Seems bad for scaling.
    modeldatactx%NSLIPglobal(1,nread) = ia(1) ! element
    modeldatactx%NSLIPglobal(2,nread) = ia(2) ! node
    modeldatactx%NSLIPglobal(3,nread) = ia(3) ! weight 1
    modeldatactx%NSLIPglobal(4,nread) = ia(4) ! weight 2

    if (NDOF.eq.3) then
        modeldatactx%NSLIPglobal(5,nread) = ia(5) ! weight 3
    else
        modeldatactx%NSLIPglobal(5,nread) = 0
    endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
    !write(*,*) "accessing elemhasSlipnode", ia(1), "of", size(elemHasSlipNode,1)
    elemHasSlipNode(ia(1)) = .true.
#endif

    if (debug) then
        write(*,*) 'rank ', getrank(), 'has ia set to: ', ia
    endif

    if (debug) then
        if (meshdatactx%e2p(ia(1)).eq.getrank()) then
            ! write(*,*) 'elem; rank', getrank(), 'will pick up combination ', ia(1), ia(2)
        else
            ! write(*,*) 'elem; rank', getrank(), 'will ignore  combination ', ia(1), ia(2)
        endif
        if (meshdatactx%v2p(ia(2)).eq.getrank()) then
            ! write(*,*) 'node; rank', getrank(), 'should pick up  combination', ia(1), ia(2)
        else
            ! write(*,*) 'node; rank', getrank(), 'should ignore combination', ia(1), ia(2)
        endif
    endif

    ! fill locally owned NSLIP, based on element number
    ! *** e2p *** so, only of elements belonging to this partition
    ! entries are selected for this partition.

#ifdef SPARSE

    if (meshdatactx%e2p(ia(1)).eq.getrank()) then

#else
    ! if not SPARSE, then plnplt/plt3d is calling rslip, which is single processors,
    ! so that further selection is not needed and all points are automatically included.
#endif

        NUMSLPlocal = NUMSLPlocal + 1
        NUMSLP = NUMSLP + 1

        if (debug) then
            write(*,*) 'rank ', getrank(), 'NUMSLP', NUMSLP, 'takes up element ', ia(1)
        endif

        ! modeldatactx%NUMSE =  modeldatactx%NUMSE+1
        modeldatactx%NSLIP(1,NUMSLP) = ia(1)
        modeldatactx%NSLIP(2,NUMSLP) = ia(2)
        modeldatactx%NSLIP(3,NUMSLP) = ia(3)
        modeldatactx%NSLIP(4,NUMSLP) = ia(4)

        if (NDOF.eq.3) then
            modeldatactx%NSLIP(5,NUMSLP) = ia(5)
        else
            modeldatactx%NSLIP(5,NUMSLP) = 0
        endif

        if (debug) then
            write(*,*) 'rank ', getrank(), 'NUMSLP', NUMSLP, 'NSLIP: ', modeldatactx%NSLIP(1:5, NUMSLP)
        endif

        ! check whether the slip has a weight (direction) assigned to it.
        if (NDOF.eq.2) then
            if ((modeldatactx%NSLIP(3,NUMSLP).eq.0) .and. &
                (modeldatactx%NSLIP(4,NUMSLP).eq.0)) then
                write(stderr,*) 'Slippery entry ', nread, ' has x and y weights both 0.'
                write(stderr,*) 'At least one should be +/- 1'
                call xit(1,"Slippery node entry has no weight")
            endif
        else
            if ((modeldatactx%NSLIP(3,NUMSLP).eq.0) .and. &
                (modeldatactx%NSLIP(4,NUMSLP).eq.0) .and. &
                (modeldatactx%NSLIP(5,NUMSLP).eq.0)) then
                write(stderr,*) 'Slippery entry ', nread, ' has x, y and z weights all 0'
                write(stderr,*) 'At least one should be +/- 1'
                call xit(1,"Slippery node entry has no weight")
            endif
        endif

#ifdef SPARSE

    endif

#else
    ! no if-block, so no endif needed!
#endif

    if (meshdatactx%v2p(ia(2)).eq.getrank()) then

#ifdef SPARSE

        if (getsize().eq.1) then
            nploc = ia(2)
        else
            nploc = meshdatactx%Glo2LocVertex(ia(2))
        endif

#else

        nploc = ia(2)

#endif

        if (nploc.ge.0) then
            localSequence = nploc - meshdatactx%lowestLocalVertexID
            do i=1,NDOF

#ifdef DIFORC_SET

                if (f(i).ne.0d0) then 
                    DIFORC(i,localSequence) = f(i)
                endif

#else

                DIFORC(i,localSequence) = DIFORC(i,localSequence) + f(i)

#endif
            enddo
        endif
    endif

enddo read_and_process


if (nread.ne.NUMSLPglobal) then
    write(stderr,6) nread,NUMSLPglobal
6   format(1x,'RSLIP: read ',I12,' data, expected ',I12)
    call xit(1," ")
endif

if (getrank().eq.0) then
    if (OUTFIL(1)) then
        nn=0
        do j=1,NUMSLP
            if (MOD(nn,50).eq.0) then
               ! new header line
               write(luout(1),7)
 7           format(1X/// &
             15X,'S L I P P E R Y   N O D E   D A T A'//5X, &
              'ELEM',10X,'NODE',5X,'DOF WEIGHT',17X, &
              'DIFFERENTIAL FORCE'//5X,'NUM',11X,'NUM',4X, &
              '1    2    3',11X,'DOF1',11X,'DOF2',11X,'DOF3'//)
            endif

            nn = nn + 1
            k = modeldatactx%NSLIP(2,j)
            write(luout(1),8) (modeldatactx%NSLIP(i,j),i=1,5), &
             (DIFORC(i,k),i=1,NDOF)
 8           format(I8,7X,I8,3I5,3X,3(3X,G12.5))
        enddo
    endif
endif

! finished reading slippery nodes.
! validate that all nodes have both negative and positive weights contributed to them
! otherwise the thing PETSc will get an inconsistent matrix and will get really ill.
! validate this in a global array, as the data exists in global numbering

nodesAreConsistent = .true.
firstInconsistentNode = .true.

#ifdef EXPERIMENTAL_FRICTION
if (itest.eq.1 .or.  (nFaultSegmentEntriesGlobal.gt.0 .and.(.not.feplt)) ) then
#else
if (itest.eq.1) then
#endif
    do iSlipEntry=1,NUMSLPglobal
        nodeID = modeldatactx%NSLIPglobal(2,iSlipEntry)
        weightx = modeldatactx%NSLIPglobal(3,iSlipEntry)
        weighty = modeldatactx%NSLIPglobal(4,iSlipEntry)

        if (weightx * weighty .eq. -1) then
            write(stderr,*) "The x and y weights of slippery ", iSlipEntry
            write(stderr,*) "Are opposing. They should be the same"
            stop "Leaving GTecton..."
        endif

        ! checkSignDuality add digits:
        ! 1:  x has positive weights
        ! 2:  x has negative weights
        ! 3:  y has positive weights
        ! 4:  y has negative weights

        if      (weightx.eq. 1) then
            if (idig(modeldatactx%checkSignDuality(nodeID),4,1).eq.0) then
                     modeldatactx%checkSignDuality(nodeID) = &
                     modeldatactx%checkSignDuality(nodeID) + 1000
            endif
        else if (weightx.eq.-1) then
            if (idig(modeldatactx%checkSignDuality(nodeID),4,2).eq.0) then
                     modeldatactx%checkSignDuality(nodeID) = &
                     modeldatactx%checkSignDuality(nodeID) + 100
            endif
        endif
    
        if      (weighty.eq. 1) then
            if (idig(modeldatactx%checkSignDuality(nodeID),4,3).eq.0) then
                     modeldatactx%checkSignDuality(nodeID) = &
                     modeldatactx%checkSignDuality(nodeID) + 10
            endif 
        else if (weighty.eq.-1) then
            if (idig(modeldatactx%checkSignDuality(nodeID),4,4).eq.0) then
                     modeldatactx%checkSignDuality(nodeID) = & 
                     modeldatactx%checkSignDuality(nodeID) + 1
            endif
        endif
    enddo


    ! allowed entries for checkSignDuality
    !   0000
    !   1100
    !   0011
    !   1111

    nodesAreConsistent = .true.

    do nodeID=1,meshdatactx%nvglobal

        if (idig(modeldatactx%checkSignDuality(nodeID),4,1) .eq. 1 .and. &
         idig(modeldatactx%checkSignDuality(nodeID),4,2) .eq. 0) then
            write(stderr,*) "Error in the slippery node configuration."
            write(stderr,*) "Every node should have both elements with positive and "
            write(stderr,*) "elements with negative weights assigned to it, "
            write(stderr,*) "but node", nodeID, "only has positive x weights."
            nodesAreConsistent = .false.
        else if (idig(modeldatactx%checkSignDuality(nodeID),4,1) .eq. 0 .and. &
             idig(modeldatactx%checkSignDuality(nodeID),4,2) .eq. 1) then
            write(stderr,*) "Error in the slippery node configuration."
            write(stderr,*) "Every node should have both elements with positive and "
            write(stderr,*) "elements with negative weights assigned to it, "
            write(stderr,*) "but node", nodeID, "only has negative x weights."
            nodesAreConsistent = .false.
        endif

        if (idig(modeldatactx%checkSignDuality(nodeID),4,3) .eq. 1 .and. &
            idig(modeldatactx%checkSignDuality(nodeID),4,4) .eq. 0) then
            write(stderr,*) "Error in the slippery node configuration."
            write(stderr,*) "Every node should have both elements with positive and "
            write(stderr,*) "elements with negative weights assigned to it, "
            write(stderr,*) "but node", nodeID, "only has positive y weights."
            nodesAreConsistent = .false.
        else if (idig(modeldatactx%checkSignDuality(nodeID),4,3) .eq. 0 .and. &
                 idig(modeldatactx%checkSignDuality(nodeID),4,4) .eq. 1) then
            write(stderr,*) "Error in the slippery node configuration."
            write(stderr,*) "Every node should have both elements with positive and "
            write(stderr,*) "elements with negative weights assigned to it, "
            write(stderr,*) "but node", nodeID, "only has negative y weights."
            nodesAreConsistent = .false.
        endif

    enddo
endif

if (.not. nodesAreConsistent) then
    stop "Leaving GTecton..."
endif

if (debug) then
  do nodeID=1,meshdatactx%nvglobal
    i=modeldatactx%checkSignDuality(nodeID)
    if (i.ne.0) write(*,*) 'Rank=',getrank(),'readfriction: signduality(node=', &
                nodeID,')=',i
  enddo
endif

if (itest.eq.1 .and. iecho.gt.0) then
    write(*,*) "Rslip passed sign duality check"
endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc

allocate(modeldatactx%SlipElemLookup(meshdatactx%neglobal))
modeldatactx%SlipElemLookup = 0

nElemsWithSlipperyNodes = 0

do iElem=1,meshdatactx%neglobal
    if (elemHasSlipNode(iElem)) then
        nElemsWithSlipperyNodes = nElemsWithSlipperyNodes + 1

        if (meshdatactx%e2p(iElem).eq.getrank()) then
            nLocalElemsWithSLipperyNodes = nLocalElemsWithSLipperyNodes + 1
        endif


!        write(*,*) "accessing SlipElemLookup", nElemsWithSlipperyNodes, &
!                   "of", size(modeldatactx%SlipElemLookup,1)

!        modeldatactx%SlipElemLookup(nElemsWithSlipperyNodes) = iElem
        modeldatactx%SlipElemLookup(iElem) = nElemsWithSlipperyNodes

    endif
enddo

!    write(*,*) "rank", getrank(), "lookup:", modeldatactx%SlipElemLookup

write(*,*) "tecin says: nElemsWithSlipperyNodes =", nElemsWithSlipperyNodes


deallocate(elemHasSlipNode)
#endif

return

300 write(stderr,301) ios
301 format(1x,'RSLIP: input error ',i8)
    write(stderr,*) "read slippery entry: ", trim(record)
    write(stderr,*) "tried to interpret with format: ", FMT_tecin_slippery_nodes
    write(stderr,*) "this failed with error code: ", ios
    call xit(1,"Read error in slippery nodes")

400 write(stderr,401) ios
401 format(1x,'RSLIP: premature EOF error ',i8)
    call xit(1,"Unexpected EOF in slippery nodes")

end subroutine RSLIP
!-------------------------------------------------------------------------------
subroutine RSLipperySKEW (NSELSD,NSLIP,SKEW,IEN,NSD,NDOF,NEN)

USE MODELDEFINITION      ! for NSLSKEW among others
USE ALGEBRA
use meshdatamodule,  only: meshdatactx
USE CONSTANTS,       only: zero, & 
                           pi, &
                           halfpi
use debugmodule,     only: debug, &
                           iecho, xit
use formatsmodule,   only: FMT_tecin_FaultParallelSurfaces
use modelctx,        only: getrank
use iomodule          
!
! reads and prints data related to maintaining differential slippery
! displacements fault parallel.
!
implicit none
!-pass
integer NSD,NEN,NDOF
integer NSLIP(5,*),NSELSD(2,NSLSKEW),IEN(NEN,*)
double precision SKEW(2,meshdatactx%nvglobal)
!-locl
logical end,quad,modified
character(len=17)record
integer nread,i,j,k,ielm,iside,node,ios,ia(2),lien(4)
data modified /.false./

!       consider using other way of finding random numers
!       initialize random number generator (see NR p 274-275)
integer :: jran,iia, iic, iim
integer :: kk, kold, ns
real rran
double precision random


! random number required to deviate values of 0 skew very slightly from zero.
! so that they are included into the computations
iim = 6075
iia = 1366
iic = 1283
! set ran to random_number
call random_number(rran)
rran = rran * mod(NUMNP,iim-1)
jran = int(rran)


!write(*,*) "rank",getrank(),"At beginning of RSlipperySkew: ", SKEW
!write(*,*) "rank",getrank(),"At beginning of RSlipperySkew: ", numnp, meshdatactx%nvglobal


if (iecho.eq.2) then
    write(stderr,1) FMT_tecin_FaultParallelSurfaces
 1   format(1x,'if (NSLSKEW > 0) {'/ &
 5x,'*** read fault parallel skew data: elements,sides  ***'/ &
 5x,'while (NSELSD(1,i)!=0 & NSELSD(1,i)!=end) {'/ &
 9x,'read ',A7,' (NSELSD(1,i),NSELSD(2,i)'/5x,'}'/1x,'}')
endif 


if (NSLSKEW.eq.0) then
    if (debug) then
        write(*,*) "Zero fault parallel surfaces defined, skip reading them"
    endif
    ! NSLSKEW is from modeldefinition
    return
else
    if (debug) then
        write(*,*) "Reading fault parallel surfaces; expecting to find: ", NSLSKEW
    endif
endif

nread = 0

100     call readrc (luin,record,'RSLSKEW',end)

if (debug) then
    write(*,*) "Read record",nread,"faultpar surface ", trim(record)
endif

if (end) then
    goto 200
endif

read(record,FMT_tecin_FaultParallelSurfaces,err=800,end=900,iostat=ios) ia(1),ia(2)

if (debug) then
    write(*,*) 'RSLSKEW read fault parallel surface', ia(1), ia(2) 
endif

if (ia(1).eq.0) then
    ! zero used as alternative to indicate 'end' of list.
    goto 200
endif
!
if (ia(1).lt.1 .or. ia(1).gt.meshdatactx%neglobal) then
    write(stderr,3) meshdatactx%neglobal
    write(*,3) meshdatactx%neglobal
 3       format(1x,'Fault parallel surfaces (sub RSLSKEW): '//&
                'element number entry out of range [1,', &
                I5,']')
    call xit(1,"Element ID of fault entry greater than the number of elements in the mesh.")
endif
!
if (NSD.eq.2) then
    if (ia(2).lt.1 .or. ia(2).gt.4) then
        write(stderr,4) 4
        write(*,4) 4

 4           format(1x,'Fault parallel surfaces (sub RSLSKEW): '//&
                    'element side number out of range ', &
                    '[1,',I1,']')
        call xit(1,"Edge ID outside interval [1,4]")
    endif
else if (NSD.eq.3) then
    if (ia(2).lt.1 .or. ia(2).gt.6) then
        ! 6 for a cuboid element. 
        ! This does not validate correctly for a purely tetraoidal mesh.
        write(stderr,4) 6
        write(*,4) 6

        call xit(1,"Edge ID outside interval [1,6]")
    endif
else
    call xit(1,"Number of spatial dimensions should be 2 or 3")
endif
!
nread = nread + 1


if (nread.gt.NSLSKEW) then
    write(stderr,5) NSLSKEW
    write(*,5) NSLSKEW

 5       format(1x,'Fault parallel surfaces (sub RSLSKEW): '//&
               ' more than the expected ',I12,' entries ', &
               'in input list')
    call xit(1,"Too many fault parallel surfaces listed")
endif


!write(*,*) 'rank', getrank(), 'writing entry ', nread, 'in NSELSD', size(NSELSD,2)

NSELSD(1,nread) = ia(1)
NSELSD(2,nread) = ia(2)

goto 100
!
200     if (nread.ne.NSLSKEW) then ! NSLSKEW as read from TECIN.DAT
    write(stderr,6) nread,NSLSKEW
 6       format(1x,'Fault parallel surfaces (sub RSLSKEW):'// &
                ' read ',i12,' entries, expected ',i12)
    call xit(1," ")
endif
!---------- finished reading --------------

if (debug) then
    write(*,*) 'rank', getrank(), "finished reading NSESLD: ", NSELSD
endif


!---------- checking integrity and comparing with slippery node entries -------
call CLEAR (lien,4,"lien")
k = 0
ielm = 0

do i=1,NSLSKEW

    ielm  = NSELSD(1,i) 
    iside = NSELSD(2,i)
    quad  = (IEN(3,ielm).ne.IEN(4,ielm))

!   get local nodal point numbers along specified element side
    call SIDENP (iside,quad,lien)

    ! loop though the vertices of the element of the element in 
    ! the fault entry    
    do j=1,4
        if (lien(j).eq.0) then
            ! local node number zero. This is not supposed to happen.
            cycle
        endif

!       set global node number
        node = IEN(lien(j),ielm)

!                if (SKEW(1,node).eq.ZERO.and.SKEW(2,node).eq.ZERO) then
!                    write(stderr,7) node,ielm,iside
!    7               format(1x,'RSLSKEW: node ',I5,' is not skewed. ',
!     >              'Element = ',I5,' side = ',I1,'.')
!                    call xit(1," ")
!                endif
!                verify that this node is slippery
        kold = k

        do kk=1,NUMSLP
!            write(*,*) 'checking slippery entry: ', kk, 'with values', NSLIP(1,k), NSLIP(2,k)
            ! search through the slippery node entries, until we find the 
            ! one matching the fault parallel surface
            k = kk
            if (ielm.eq.NSLIP(1,k) .and. &
                node.eq.NSLIP(2,k)) then
!                write(*,*) 'match with slippery entry: ', NSLIP(1,k), NSLIP(2,k)
                goto 300
            endif
        enddo

!                write(stderr,8) node,ielm
!    8           format(1x,'RSLSKEW: node ',I5,' in element ',I5,
!     >          'is not defined as slippery.')
!                call xit(1," ")

300     if (j.eq.1) then
            cycle
        endif

!       Verify that slippery dofs are uniform along element sides
        ns = 0

        do kk=1,NDOF

            if (NSLIP(2+kk,k).ne.0) then
                ns = ns + 1
            endif


!            if (NSLIP(2+kk,kold).ne.NSLIP(2+kk,k)) then
            if (abs(NSLIP(2+kk,kold)) .ne. &
                abs(NSLIP(2+kk,k)))         then

                write(stderr,9) ielm,kold,k,NSLIP(2,kold),NSLIP(2,k)
                write(*,9) ielm,kold,k,NSLIP(2,kold),NSLIP(2,k)

 9                   format(1x,'RSLSKEW: non-uniform slippery ', &
                 'degree of freedom along element ',I12,':'/10x, &
                 'NSLIP entries ',I12,' and ',I12,' for nodes ', &
                 I12,' and ',I12,' have different dofs.')
                call xit(1,"Non-uniform slippery degrees of freedom")
            endif

        enddo

        if (ns.eq.NDOF) then
            write(stderr,10) iside,ielm,NDOF
            write(*,10) iside,ielm,NDOF

10               format(1x,'RSLSKEW: number of degrees of freedom ', &
             'along side ',I1,' of element ',I12,' is ',I1/10x, &
             'inconsistent with slip along a discontinuity.')
!            call xit(1,"Number of degrees of freedom incorrect")
        endif
    enddo
enddo



! Modify NSLIP so that local x-dofs are parallel to fault
do i=1,NSLSKEW
    ielm  = NSELSD(1,i)
    iside = NSELSD(2,i)
    quad  = (IEN(3,ielm).ne.IEN(4,ielm))
    call SIDENP (iside,quad,lien)
    do j=1,4
        if (lien(j).eq.0) then
            cycle
        endif
        node = IEN(lien(j),ielm)
        do kk=1,NUMSLP
            k = kk
            if (node.eq.NSLIP(2,k)) then

                if (NTYPE.eq.5) then
                    if (NSLIP(4,k).ne.0) then
                        write(stderr,12) ielm,node,k
12                           format(1x,'RSLSKEW: please specify ', &
                         'NSLIP data for element ',I12,', node ', &
                         I12,' (entry ',I12,' so that the local ', &
                         'y-axis is perpendicular to the ', &
                         'fault plane')
                        call xit(1," ")
                    endif
                else
                    if (NSLIP(3,k).eq.0) then
!                               Label nodes when modified
                        NSLIP(2,k) = -NSLIP(2,k)

                        call IntegerSWAP(NSLIP(3,k),NSLIP(4,k))

                        write(stderr,11) NSLIP(1,k),node
11                           format(1x,'RSLSKEW warning: modified ', &
                         'd.o.f. for element ',I12,', node ',I12, &
                         ','/18x,'so that the local x-axis is ', &
                         'parallel to the fault')
                        modified = .true.
                    endif
                endif
            endif
        enddo
    enddo
enddo


if (modified) then
!           Search labeled nodes and modify associated skew angle
    do 700 node=1,NUMNP
        do k=1,NUMSLP
            if (node.eq.-NSLIP(2,k)) then
                SKEW(1,node) = SKEW(1,node) - HALFPI
                if (SKEW(1,node).eq.ZERO) then
                    ! if original skew is zero, the angle 
                    ! is not updated in further timesteps.
                    ! to prevent this, the angle is deviated slightly from 0.
                    jran = mod(jran*iia+iic,iim)
                    random = DBLE(jran)/DBLE(iim)
                    SKEW(1,node) = -1d-11 + 2d-11*random!DBLE(RAND())
                endif
                goto 700
            endif
        enddo
    700 continue
!
!           Remove labels
    do i=1,NUMSLP
        if (NSLIP(2,i).lt.0) then
            NSLIP(2,i) = -NSLIP(2,i)
        endif
    enddo

    if (DEBUG) then
        write(stderr,13)
13           format(1x,'RSLSKEW: MODIFIED NSLIP ARRAY:')
        do k=1,NUMSLP
            write(stderr,14) k,(NSLIP(i,k),i=1,2+NDOF)
14               format(' NSLIP(*,',I8,') = ',I8,1X,I8,1X,3(1X,I3))
        enddo
        write(stderr,15)
15           format(/1x,'RSLSKEW: MODIFIED SKEW ARRAY:')
        do k=1,NUMNP
            if (SKEW(1,k).ne.ZERO .or. SKEW(2,k).ne.ZERO) &
             write(stderr,16) k,SKEW(1,k)*1.8d2/PI, &
              SKEW(2,k)*1.8d2/PI
16                format(' SKEW(*,',I8,') = ',1PG12.6,1X,1PG12.6)
        enddo
    endif

endif

!
if (OUTFIL(1)) then
    i=0
    do j=1,NSLSKEW
        if (MOD(i,50).eq.0) then
            write(luout(1),17)
17           format(1X///15X,'S L I P P E R Y',4X,'P A R A L L E L', &
                    4X,'S K E W',4X,'D A T A'//5X,'ELEM',7X,'SIDE'/ &
                   /6X,'NUM'//)
        endif
        i = i + 1
        write(luout(1),18) NSELSD(1,j),NSELSD(2,j)
18      format(I8,7X,I5)
    enddo
endif


!write(*,*) "rank",getrank(),"SKEW At end of RSlipperySkew: ", SKEW
!write(*,*) "rank",getrank(),"NSELSD At end of RSlipperySkew: ", NSELSD
!   NSELSD is just like in TECIN...

!
return
!
800     write(stderr,19) ios
19   format(1x,'RSLSKEW: input error ',i8)
call prerr('RSLSKEW')
call xit(1,"Input error when reading fault parallel surfaces")

900     write(stderr,20) ios
20   format(1x,'RSLSKEW: premature EOF error ',i8)
call prerr('RSLSKEW')
call xit(1,"Unexpected EOF when reading fault parallel surfaces")
!
end
!-------------------------------------------------------------------------------
subroutine IntegerSwap (n,m)
implicit none
integer :: i,n,m
i = n
n = m
m = i
return
end
!-------------------------------------------------------------------------------
!subroutine RWINK (IWINK,WINK,X,NSD,NDOF,NUMNP,NWINK,IFLAG)
subroutine RWINK (X,NSD,NDOF,NUMNP,NWINK,IFLAG)


#ifdef EXPERIMENTAL_CYCLIC
use modelDefinition, only: UseCyclicSlip
#endif
USE LABELMODULE
USE MODELDATAMODULE
USE MESHDATAMODULE
USE ALGEBRA
use debugmodule,     only: iecho, parallellog, Verified_MPI_Barrier
use modelctx,        only: getrank
use iomodule          
!
! Routine for reading and printing data on Winkler restoring forces.
!       WINK(NDOF,NUMNP)        = values of winkler restoring spring constant.
!       IWINK(NDOF,NUMNP)       = application mode:
!                                 IWINK =  0, no Winkler forces,
!                                 IWINK =  N, applied at time step N-1
!                                 IWINK = -N, applied up until step N,
!                                             removed thereafter
!

implicit none
!-pass
integer NSD,NDOF,NUMNP,iflag
integer :: nWink
double precision :: x
dimension X(NSD,*)
!-locl
logical           :: warn_about_phase
character(len=13) :: fmt
integer           :: i, n, nn
integer           :: idum
warn_about_phase = .true.
fmt = '(I12,3I5)  '

if (iecho.eq.2) then
    if (iflag.eq.0) then
        write(stderr,1) fmt
 1       format(1x,'*** Winkler condition codes ***'/ &
            1x,'n=1; while (n != 0) {'/ &
            5x,'read ',A8,' n,(modeldatactx%IWINK(i,n),i=1,NDOF)'/ &
            1x,'}')
    else
        write(stderr,101) fmt
 101   format(1x,'*** Differential Winkler condition codes ***'/ &
            1x,'n=1; while (n != 0) {'/ &
            5x,'read ',A8,' n,(modeldatactx%IWINX(i,n),i=1,NDOF)'/ &
            1x,'}')


    endif
endif

if (iflag.eq.0) then
    ! normal Winkler forces
    modeldatactx%IWINK(:,:) = 0
else
    ! differential Winkler forces
    modeldatactx%IWINX(:,:) = 0
endif

if (iflag.eq.0) then
    ! Winkler forces
    call IGEN(modeldatactx%IWINK,NDOF,meshdatactx%nvglobal,&
              fmt,NWINK,10,"Winkler Forces indices")
else
    ! differential Winkler forces
    call IGEN(modeldatactx%IWINX,NDOF,meshdatactx%nvglobal,&
              fmt,NWINK,11,"differential Winkler Forces")
endif

#ifdef EXPERIMENTAL_CYCLIC

if (UseCyclicSlip) then  ! from modeldefinition
    if (IFLAG.eq.1) then  ! only differential winkler forces -> winx
        do n=1,NUMNP
#ifdef ZISFAULTNORMAL
            do i=1,NDOF-1
#else
            do i=1,NDOF
#endif
                if (modeldatactx%IWINX(i,n).ne.0 .and. &
                ABS(modeldatactx%IWINX(i,n)).le.10) then
                    write(*,*) 'Error with Winkler forces'
                    write(*,*) 'vertex: ', n
                    write(*,*) 'd.o.f.: ', i
                    write(*,*) 'IWINK : ', modeldatactx%IWINX(i,n)
                    stop 'periodic fault unlocking error'
                endif
            enddo
#ifdef ZISFAULTNORMAL
            if (modeldatactx%IWINX(NDOF,n).ne.0 .and. warn_about_phase) then
                write(*,*) '*********** WARNING *************'
                write(*,*) 'IWINK(',NDOF,',',n,')= ', modeldatactx%IWINX(NDOF,n)
                write(*,*) 'IT IS ASSUMED THAT THIS ENTRY REFERS TO A PHASE DELAY'
                write(*,*) ''
                write(*,*) 'THIS IS A FEATURE OF COMPILING WITH -DCYCLIC ICW COMPILING WITH -DZISFAULTNORMAL'
                write(*,*) 'SIMILAR WARNINGS FOR OTHER NODES WILL NOT BE SHOWN'
                write(*,*) '******* END WARNING *************'
                warn_about_phase = .false.
            endif
#endif
        enddo
    endif
endif
#endif

fmt = '(I12,3G14.0)'

if (iecho.eq.2) then
    if (iflag.eq.0) then
        write(stderr,2) fmt
 2   format(1x,'*** Read Winkler spring constants ***'/ &
 1x,'n=1; while (n != 0) {'/ &
 5x,'read ',A12,' n,(WINK(i,n),i=1,NDOF)'/ &
 1x,'}')
    else
        write(stderr,102) fmt
 102   format(1x,'*** Read differential Winkler spring constants ***'/ &
 1x,'n=1; while (n != 0) {'/ &
 5x,'read ',A12,' n,(WINX(i,n),i=1,NDOF)'/ &
 1x,'}')

    endif
endif


if (iflag.eq.0) then
    modeldatactx%WINK = 0d0
    call FGEN(modeldatactx%WINK,NDOF,meshdatactx%nvglobal,fmt,&
              idum,"Winkler forces, magnitures")
else
   ! Nicolai
    modeldatactx%WINX = 0d0
    call FGEN(modeldatactx%WINX,NDOF,meshdatactx%nvglobal,fmt,&
              idum,"differential Winkler forces, magnitudes")
endif


!
if (.not.OUTFIL(1) .or. NWINK.le.0) then
    return
endif
!
nn = 51
do n=1,NUMNP
    if (nn.gt.50) then
        if (IFLAG.EQ.0) write(luout(1),3)
 3           format(1X///2x, &
         'W I N K L E R   R E S T O R I N G   F O R C E S'//)
        if (IFLAG.EQ.1) write(luout(1),4)
 4           format(1X///2x, &
         'D I F F E R E N T I A L',3x,'W I N K L E R',3x, &
         'R E S T O R I N G',3x,'F O R C E S'//)
        if (IFLAG.EQ.2) write(luout(1),5) 1
 5           format(1X///2x, &
         'D I F F E R E N T I A L',3x,'W I N K L E R',3x, &
         'F O R C E S',3x,'O N'/2x,'T R I P L E',3x, &
         'J U N C T I O N',3x,'D O F ',I1//)
        if (IFLAG.EQ.3) write(luout(1),5) 2
        write(luout(1),7)(LABELD(I),I=1,NDOF)
 7           format(6x,'KEY TO APPLICATION MODE CODES:'// &
         10x,'MODE.EQ.0, NO WINKLER RESTORING FORCES;',1x, &
         'THE NUMERICAL'/ &
         21x,'VALUE TO THE LEFT IS MEANINGLESS.'/ &
         10x,'MODE.GT.0, FORCES APPLIED IN STEP MODE-1'/ &
         10x,'MODE.LT.0, FORCES REMOVED IN STEP IABS(MODE)'/// &
         4X,'NODE # ',A4,'-WINKLR     MODE',5(2X,A4, &
         '-WINKLR     MODE')/)
        write(luout(1),8)
 8           format(/)
        nn = 1
    endif
    do i=1,NDOF
!                IF (IWINK(i,n)) 200,100,200
        if (iflag.eq.0) then
            if (modeldatactx%IWINK(i,n).ne.0) goto 200
        else
            if (modeldatactx%IWINX(i,n).ne.0) goto 200
        endif
    enddo

    cycle

    200 if (iflag.eq.0) then
        write(luout(1),9) n,(modeldatactx%WINK(i,n),modeldatactx%IWINK(i,n),i=1,NDOF)
    else
        write(luout(1),9) n,(modeldatactx%WINX(i,n),modeldatactx%IWINX(i,n),i=1,NDOF)
    endif
 9       format(1X,I8,2X,6(1PE12.5,2X,I5,3X))



    nn = nn + 1

enddo
!
return
end

!-----------------------------------------------------------------------

#ifdef SHEARZONE_TRACTIONS

subroutine readShearTrax()

use iomodule,       only    :   luin, stderr
use MODELDEFINITION,only    :   NUMELSHEARZONE, NUMNPSHEARZONE, NUMSLP, NTYPE
use MODELDATAMODULE,only    :   modeldatactx
USE MESHDATAMODULE, only    :   meshdatactx
use MODELTOPOLOGY,  only    :   NDOF
use debugmodule,    only    :   iecho, xit, debug
use MODELCTX,       only    :   getrank
use formatsmodule,  only    :   FMT_tecin_ShearTrax, FMT_tecin_ShearTraxQT
use algebra,        only    :   arrayContains

implicit none

integer            ::  nread=0,elno,n,nodenum,i,num_slp_el_nodes,n_nodes_per_side,weight
integer            ::  NodeSnp(3),ielm,side,iostat
integer,allocatable::  ShearZoneNodes(:)
double precision   ::  eta,thick,temperature,activQ
character(len=256) ::  record
logical            ::  end

NUMNPSHEARZONE = 0
if (NUMELSHEARZONE.eq.0) return ! If there are no shear zone entries in the model
weight = 0

call readrc(luin,record,'readShearTrax',end)
read(record,'(I5)') modeldatactx%shearnpl
if ((iecho.eq.2 .or. debug) .and. getrank().eq.0) write(*,*) &
       'Read exponent for power-law creep shear zone viscosity:', modeldatactx%shearnpl

if ((iecho.eq.2 .or. debug) .and. getrank().eq.0) then
    write(*,1) NUMELSHEARZONE, FMT_tecin_ShearTrax
1   format(1x,'if (NUMPR <> 0){'/ &
      5X,'*** Attempting to read ',I6,' shear zone traction entries ***'// &
      5X,'read ',A12,' ELNOSHEARZONE(n),shearETA(n),shearTHICK(n)'/ &
      5X,'if (IELNO(n)==0) return'/1x,'}')
endif

if (NUMSLP.le.0) then
    if (getrank().eq.0) write(stderr,10) NUMSLP,NUMELSHEARZONE
10  format(1x,"Shear traction input error: NUMSLP = ",I8," while NUMELSHEARZONE = ",I8)
    call exitp(1)
endif

end=.false. ! Necessary because the readrc above will set end as positive if shearnpl is 0
loop_read_and_verify: do while (nread.le.NUMELSHEARZONE .and. .not.end)
    call readrc(luin,record,'readShearTrax',end)
    if (end) then
        if (debug) write(*,*) "RANK:", getrank(), ", readShearTrax found the end of the shear zone " // &
         "tractions list after reading", nread, "of", NUMELSHEARZONE, "expected entries"
        if (nread.ne.NUMELSHEARZONE) then
            write(stderr,*) "RANK: ", getrank(), ", readShearTrax expected ", NUMELSHEARZONE, &
                            "shear zone traction entries, but found only ", NREAD
            call xit(1,"Too few shear zone entries")
        endif
    else
        if (modeldatactx%shearnpl.gt.0) then
            read(record,FMT_tecin_ShearTraxQT) elno, eta, thick, temperature, activQ ! Read (working) variables from record
        else
            read(record,FMT_tecin_ShearTrax) elno, eta, thick ! Read (working) variables from record
            activQ      = 0d0
            temperature = 0d0 ! Initialise unused array to avoid compiler warnings
        endif
        if (elno.eq.0) then ! Take 0 as a signal of the end of the input data section (as an alternative to "end")
            end = .true.
        else
            nread = nread+1 ! Count number of entries read

            ! Check number of nodes on the current element that are slippery
            if (debug) write(*,*) "Rank", getrank(), "checking whether any of the", NUMSLP, &
              " slippery node entries is in element ", elno

            num_slp_el_nodes = 0 ! Reset number of slippery nodes in element
            loop_slippery_node_entries: do n=1,NUMSLP
                if (modeldatactx%NSLIP(1,n).eq.elno) then
                    nodenum = modeldatactx%NSLIP(2,n)
                    if (debug) write(*,*) "Rank", getrank(), "node ", nodenum, "is slippery and in the element"

                    ! verify consistency of weights
                    if (weight.eq.0) then
                        weight = modeldatactx%NSLIP(3,n)
                    else
                        if (weight.ne. modeldatactx%NSLIP(3,n)) then
                            if (getrank().eq.0) write(stderr,*) "readShearTrax: element ",elno, &
                              " has different weights than previous shear zone entries"
                            call exitp(1)
                        endif
                    endif
                    
                    do i=1,NDOF
                        if (modeldatactx%DIFORC(i,nodenum) .ne. 0.0) then
                            write(stderr,*) "readShearTrax, rank ", getrank(), ": the slippery node number ", &
                            nodenum, ", in element number", elno, &
                            ", to which you want to apply shear tractions, already has differential forces applied"
                            call xit(1,"Differential forces and shear zone tractions applied to same node")
                        endif
                    enddo 
                    num_slp_el_nodes = num_slp_el_nodes + 1
                endif
            enddo loop_slippery_node_entries
            if (debug) write(*,*) "Number of slippery nodes in element: ", num_slp_el_nodes

            modeldatactx%ELNOSHEARZONE(nread)   =   elno
            modeldatactx%shearETA(nread)        =   eta
            modeldatactx%shearTHICK(nread)      =   thick
            modeldatactx%shearQpl(nread)        =   activQ
            modeldatactx%shearTEMP(nread)       =   temperature
        endif
    endif
enddo loop_read_and_verify

if (.not.end) then
    if (getrank().eq.0) then
        write(stderr,*) "RANK: ", getrank(), ", readShearTrax expected ", NUMELSHEARZONE, &
          "shear zone traction entries, but did not find 'end' or '0' after entry no. ", nread
    endif
    call exitp(1)
endif

if ((iecho.eq.2 .or. debug) .and. getrank().eq.0) write(*,*) &
    "Finished reading all expected shear zone traction entries"

! identify side numbers for the given elements and count the partition-local number of nodes in the shear zone
if (NTYPE.eq.6) then
    n_nodes_per_side = 3
else
    n_nodes_per_side = 2
endif
allocate(ShearZoneNodes(NUMELSHEARZONE*3),stat=iostat) ! Upper limit from 3D: three nodes per element face
if (iostat.ne.0) then
    write(stderr,80) getrank(),iostat,NUMELSHEARZONE*3
80  format("Rank ",I3," encountered error ",I5," while allocating ShearZoneNodes array with dimension ",I8)
    call exitp(1)
endif
ShearZoneNodes = 0
loop_shear_zone_elements: do ielm=1,NUMELSHEARZONE
    elno = modeldatactx%ELNOSHEARZONE(ielm)
    num_slp_el_nodes = 0 
    do n=1,NUMSLP
        if (modeldatactx%NSLIP(1,n).eq.elno) then
            num_slp_el_nodes = num_slp_el_nodes + 1
            if (num_slp_el_nodes.gt.n_nodes_per_side) then
                write(stderr,*) 'readShearTrax fatal: too many nodes'
                call exitp(1)
            endif
            NodeSnp(num_slp_el_nodes) = modeldatactx%NSLIP(2,n)
        endif
    enddo
    if (num_slp_el_nodes .lt. n_nodes_per_side) then
        write(stderr,*) "readShearTrax, rank ", getrank(), ": element number ", elno, "should have ", &
           n_nodes_per_side, "slippery nodes, but only has ", num_slp_el_nodes
        call xit(1,"Wrong shear traction input (not an interface element")
    endif
    call npside(NodeSnp,meshdatactx%IEN(1,elno),side)
    modeldatactx%SIDESHEARZONE(ielm) = side
    if (debug) then
        write(*,70) getrank(),elno,modeldatactx%SIDESHEARZONE(ielm)
70      format(" RANK ",I3,": readShearTrax element ",I8," has side ",I1, &
            " bounding the shear zone")
    endif
    increment_number_shear_zone_nodes: do i=1,n_nodes_per_side
        if (NUMNPSHEARZONE.gt.0) then
            call arrayContains(ShearZoneNodes,NUMNPSHEARZONE,NodeSnp(i),n)
        else
            n = 1
        endif
        if (n.gt.NUMNPSHEARZONE) then
            NUMNPSHEARZONE = NUMNPSHEARZONE + 1 
            ShearZoneNodes(NUMNPSHEARZONE) = NodeSnp(i)
        endif
    enddo increment_number_shear_zone_nodes
enddo loop_shear_zone_elements

if (debug) then
    write(*,72) getrank(),NUMNPSHEARZONE
72  format(" RANK ",I3," readShearTrax: ",I8," nodes along the shear zone")
endif

if (NUMNPSHEARZONE.gt.0) then
    ! copy shear zone node numbers to (more compact) NPNOSHEARZONE array
    if (allocated(modeldatactx%NPNOSHEARZONE)) then
        deallocate(modeldatactx%NPNOSHEARZONE)
    endif
    allocate(modeldatactx%NPNOSHEARZONE(NUMNPSHEARZONE), STAT=iostat)
    if (iostat.ne.0) then
        write(stderr,73) getrank(),iostat,NUMNPSHEARZONE
73      format("Rank ",I3," encountered error ",I5," while allocating NPNOSHEARZONE array with dimension ",I8)
        call exitp(1)
    endif
    do i=1,NUMNPSHEARZONE
        modeldatactx%NPNOSHEARZONE(i) = ShearZoneNodes(i)
    enddo

    ! Allocate total differential force array
    if (allocated(modeldatactx%DFMAXSHEARZONE)) then
        deallocate(modeldatactx%DFMAXSHEARZONE)
    endif
    allocate(modeldatactx%DFMAXSHEARZONE(NUMNPSHEARZONE), STAT=iostat)
    if (iostat.ne.0) then
        write(stderr,74) getrank(),iostat,NUMNPSHEARZONE
74      format("Rank ",I3," encountered error ",I5," while allocating DFMAXSHEARZONE array with dimension ",I8)
        call exitp(1)
    endif
endif

deallocate(ShearZoneNodes)

return
end subroutine

#endif /* SHEARZONE_TRACTIONS */
!-----------------------------------------------------------------------


#ifdef EXPERIMENTAL_FRICTION
subroutine readFriction

use modeldefinition, only: nFaultSegmentEntriesGlobal, &
                           radius, &
                           ishell, &
                           NUMSLPglobal, &
                           NWINKX, &
                           NTYPE
use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use formatsmodule,   only: FMT_tecin_dat_friction_2D, &
                           FMT_tecin_dat_friction_3D
use iomodule
use debugmodule,     only: iecho, parallelLog, debug
use modeltopology,   only: NDOF
use modelctx,        only: getrank ! for debugging
use algebra,         only: arrayContains, &
                           distanceBetweenTwoPoints, &
                           distanceBetweenTwoPointsOnASphere, &
                           areaOfTriangleIn3D
use materialsmodule, only: matdatactx, &
                           MaterialIDOfElement
use petscksp        ! for MPI stuff
use frictionmodule,  only: friction

implicit none
integer            :: iostat, iEntry, iError
integer            :: entriesRead,i
logical            :: end
character(len=255) :: record
double precision, allocatable :: exchangeCoordinates(:)
double precision, allocatable :: allCoordinates1(:)
double precision, allocatable :: allCoordinates2(:)
double precision, allocatable :: allCoordinates3(:)
double precision   :: thickness
double precision   :: nodalFriction
double precision   :: nodalStrength
integer            :: pos, z
integer            :: iDOF, iPoint, iElem, iFrictionEntry
integer, allocatable :: frictionElements(:)
integer, allocatable :: frictionAllElements(:)
logical              :: point1inElem
logical              :: point2inElem
logical              :: point3inElem
double precision     :: sideLength
integer              :: materialID
double precision, allocatable :: area(:)

NWINKX=0

if (iecho.eq.2) then
    write(*,*) '*****************************'
    if (NDOF.eq.2) then
        write(*,10) FMT_tecin_dat_friction_2D
10      format(1x,"Reading friction entries with format ",a16/ &
         5x,"Fault number, two node numbers of this fault segment, ", &
         "static friction value, differential Winkler constant")
    else
        write(*,11) FMT_tecin_dat_friction_3D
11      format(1x,"Reading friction entries with format ",a16/ &
         5x,"Fault number, three node numbers of this fault segment, ", &
         "static friction value, differential Winkler constant")
    endif
endif 

if (debug) print *,'tecin/readfriction: nFaultSegmentEntriesGlobal=',nFaultSegmentEntriesGlobal

allocate(friction(nFaultSegmentEntriesGlobal),stat=iostat)
if (iostat.ne.0) then
    write(stderr,*) "tecin/readfriction: Could not allocate friction array of size", &
                     nFaultSegmentEntriesGlobal
    write(stderr,*) "** FATAL ** tecin/readfriction: error=", iostat
    call exitp(1)
endif

! ********** IT WOULD BE MORE LOGICAL TO HAVE ELM#/SIDE# COMBINATIONS ******
do iEntry=1,nFaultSegmentEntriesGlobal+1
    call readrc (luin,record,"readFriction",end)
    if (end) then
        if (iEntry .eq. nFaultSegmentEntriesGlobal+1) then
            goto 50
        endif
    endif
    if (debug) write(*,*) 'tecin/readfriction: ',record
    if (NDOF.eq.2) then
        read(record, FMT_tecin_dat_friction_2D, err=100,end=110) &
                                                      friction(iEntry)%faultNumber,& 
                                                      friction(iEntry)%point(1),&
                                                      friction(iEntry)%point(2),&
                                                      friction(iEntry)%staticFriction,&
                                                      friction(iEntry)%Winkler
    else if (NDOF.eq.3) then
        read(record, FMT_tecin_dat_friction_3D, err=100,end=110) &
                                                      friction(iEntry)%faultNumber,&
                                                      friction(iEntry)%point(1),&
                                                      friction(iEntry)%point(2),&
                                                      friction(iEntry)%point(3),&
                                                      friction(iEntry)%staticFriction,&
                                                      friction(iEntry)%Winkler
    endif
    NWINKX=NWINKX+1 
enddo

50 if (iecho.eq.8) then
    call parallelLog("readFriction","finished reading friction entries")
endif

! number of nodes with diff Winklers on each fault is one more than number of
! segments. works only if last fault number is the highest:
NWINKX=NWINKX+friction(nFaultSegmentEntriesGlobal)%faultNumber
if (debug) then
    write(*,*) "readFriction: (**SUSPECT **) NWINKX = ",NWINKX
endif

! verify input entries. 
! ********** NONE OF THE TESTS BELOW WOULD BE NEEDED WHEN ELM/SIDE WOULD BE
! INPUT ****

!***********************
! TEST 1: verify that node entries appear in the NSLIP array 
!***********************
do iEntry=1,nFaultSegmentEntriesGlobal
    call arrayContains(modeldatactx%NSLIPglobal(2,:), &
                       numslpglobal, &
                       friction(iEntry)%point(1), &
                       pos)
    if (pos.gt.numslpglobal) then
        write(stderr,*) "readFriction: node ", friction(iEntry)%point(1), &
                   " is in friction entry ", iEntry, &
                   " but is not a slippery node"
        call exitp(1)
    endif

    call arrayContains(modeldatactx%NSLIPglobal(2,:), &
                       numslpglobal, &
                       friction(iEntry)%point(2), &
                       pos)
    if (pos.gt.numslpglobal) then
        write(stderr,*) "readFriction: node ", friction(iEntry)%point(2), &
                   " is in friction entry ", iEntry, &
                   " but is not a slippery node"
        call exitp(1)
    endif

    if (NDOF.eq.3) then
        call arrayContains(modeldatactx%NSLIPglobal(2,:), &
                           numslpglobal, &
                           friction(iEntry)%point(3), &
                           pos)
        if (pos.gt.numslpglobal) then
            write(stderr,*) " readFriction: node", friction(iEntry)%point(3), &
                       " is in friction entry ", iEntry, &
                       " but is not a slippery node"
            call exitp(1)
        endif
    endif

enddo

if (iecho.eq.8) then
    call parallelLog("readFriction","passed friction test 1")
endif

!***********************
! TEST 2: make sure that every combination of friction nodes occurs in two elements
! Because the total IEN cannot fit on one processor, we will let each processor find its 
! own elements
!***********************
allocate(frictionElements(nFaultSegmentEntriesGlobal),stat=iostat)
if (iostat.ne.0) then
    write(stderr,*) "tecin/readfriction: Could not allocate frictionElements array of size", &
                     nFaultSegmentEntriesGlobal
    write(stderr,*) "** FATAL ** tecin/readfriction: error=", iostat
    call exitp(1)
endif

frictionElements = 0
! each partition checks whether its own elements fit, and if so, we administrate this
! with an increment in frictionElements. In the end, we sum this over the processors,
! and it should be a 2 for every entry.

do iEntry=1,nFaultSegmentEntriesGlobal
    do iElem=1,meshdatactx%nelocal
        point1inElem = .false.
        point2inElem = .false.
        point3inElem = .false.

        do iPoint=1,4
            if (meshdatactx%IEN(iPoint,iElem).eq.friction(iEntry)%point(1)) then
                point1inElem = .true.
            endif
            if (meshdatactx%IEN(iPoint,iElem).eq.friction(iEntry)%point(2)) then
                point2inElem = .true.
            endif
            if (NDOF.eq.3) then
                if (meshdatactx%IEN(iPoint,iElem).eq.friction(iEntry)%point(3)) then
                    point3inElem = .true.
                endif
            endif
        enddo

        if (NDOF.eq.2) then
            if (point1inElem .and. point2inElem) then
                frictionElements(iEntry) = frictionElements(iEntry) + 1
            endif
        else
            if (point1inElem .and. point2inElem .and. point3inElem) then
                frictionElements(iEntry) = frictionElements(iEntry) + 1
            endif
        endif
    enddo
enddo

allocate(frictionAllElements(nFaultSegmentEntriesGlobal),stat=iostat)
if (iostat.ne.0) then
    write(stderr,*) "tecin/readfriction: Could not allocate frictionAllElements array of size", &
                     nFaultSegmentEntriesGlobal
    write(stderr,*) "** FATAL ** tecin/readfriction: error=", iostat
    call exitp(1)
endif


#ifdef SPARSE

call mpi_allreduce(frictionElements,       &  ! send
                   frictionAllElements,    &  ! receive
                   nFaultSegmentEntriesGlobal, &  ! number of variables in array
                   MPI_int,                &  ! type
                   MPI_sum,                &  ! operator
                   MPI_COMM_WORLD,         &
                   iError)
#else
frictionAllElements = frictionElements
#endif

deallocate(frictionElements)

do iEntry=1,nFaultSegmentEntriesGlobal
    if (frictionAllElements(iEntry).eq.2) then
        ! everything OK, the section has two adjacent elements
    else
        write(stderr,*) "Friction data incorrect"
        write(stderr,*) "Friction entry ", iEntry
        if (NDOF.eq.2) then
            write(stderr,*) "has points", friction(iEntry)%point(1), &
                              "and", friction(iEntry)%point(2)
        else
            write(stderr,*) "has points", friction(iEntry)%point(1), &
                              "and", friction(iEntry)%point(2), &
                              "and", friction(iEntry)%point(3)
        endif
        write(stderr,*) "But there is/are only", frictionAllElements(iEntry), &
                   "elements having both this combination of points."
        call exitp(1)
    endif

enddo

deallocate(frictionAllElements)

if (iecho.eq.8) then
    call parallelLog("readFriction","passed friction test 2")
endif

!***********************
! END TEST 2
!***********************


if (outfil(1)) then
  do iFrictionEntry=1,nFaultSegmentEntriesGlobal
    if (mod(iFrictionEntry-1,50).eq.0) then
      write(luout(1),*) " "
      write(luout(1),"(a)") " F R I C T I O N "
      write(luout(1),*) " "
      if (NDOF.eq.2) then
        write(luout(1),"(a)") &
         "       FAULT       START        END          FRICTION       WINKLER"
        write(luout(1),"(a)") &
         "                   POINT       POINT"
      else
        write(luout(1),"(a)") &
         "       FAULT      POINT 1     POINT 2     POINT 3       FRICTION      WINKLER"
      endif
      write(luout(1),*) " "
    endif
    if (NTYPE.ne.6) then
      write(luout(1),"(I12,I12,I12,5X,E14.6,E14.6)") &
       friction(iFrictionEntry)%faultNumber, &
       friction(iFrictionEntry)%point(1), &
       friction(iFrictionEntry)%point(2), &
       friction(iFrictionEntry)%staticFriction, &
       friction(iFrictionEntry)%Winkler
    else
      write(luout(1),"(I12,I12,I12,I12,5X,E14.6,E14.6)") &
       friction(iFrictionEntry)%faultNumber, &
       friction(iFrictionEntry)%point(1), &
       friction(iFrictionEntry)%point(2), &
       friction(iFrictionEntry)%point(3), &
       friction(iFrictionEntry)%staticFriction, &
       friction(iFrictionEntry)%Winkler
    endif
  enddo
endif


!****
! Get coordinates to set fault segment area
!****
allocate(exchangeCoordinates(NDOF*nFaultSegmentEntriesGlobal))
exchangeCoordinates = 0d0

   ! exchange coordinates of first point
do iEntry=1,nFaultSegmentEntriesGlobal
    ! (array, length, value, pos)
    call arrayContains(meshdatactx%gloInd, &
                       meshdatactx%nvLocal, &
                       friction(iEntry)%point(1), &
                       pos)

    if (pos.le.meshdatactx%nvLocal) then
        do iDOF=1,NDOF
            exchangeCoordinates(NDOF*(iEntry-1)+iDOF) = &
                 meshdatactx%X(iDOF,pos)
        enddo
    else
        ! this point does not belong to this partition
    endif

enddo

allocate(allCoordinates1(NDOF*nFaultSegmentEntriesGlobal))

#ifdef SPARSE
call mpi_allreduce(exchangeCoordinates,    &  ! send
                   allCoordinates1,        &  ! receive
                   NDOF*nFaultSegmentEntriesGlobal, &  ! number of variables in array
                   MPI_double,             &  ! type
                   MPI_sum,                &  ! operator
                   MPI_COMM_WORLD,         &
                   iError)

#else
allCoordinates1 = exchangeCoordinates
#endif

! and exchange coordinates of second point
exchangeCoordinates = 0d0

do iEntry=1,nFaultSegmentEntriesGlobal
    ! (array, length, value, pos)
    call arrayContains(meshdatactx%gloInd, &
                       meshdatactx%nvLocal, &
                       friction(iEntry)%point(2), &
                       pos)


    if (pos.le.meshdatactx%nvLocal) then
        do iDOF=1,NDOF
            exchangeCoordinates(NDOF*(iEntry-1)+iDOF) = &
                 meshdatactx%X(iDOF,pos)
        enddo
    else
        ! this point does not belong to this partition
    endif

enddo

allocate(allCoordinates2(NDOF*nFaultSegmentEntriesGlobal))

#ifdef SPARSE
call mpi_allreduce(exchangeCoordinates,    &  ! send
                   allCoordinates2,        &  ! receive
                   NDOF*nFaultSegmentEntriesGlobal, &  ! number of variables in array
                   MPI_double,             &  ! type
                   MPI_sum,                &  ! operator
                   MPI_COMM_WORLD,         &
                   iError)
#else
allCoordinates2 = exchangeCoordinates
#endif

if (ndof.eq.3) then
    ! and exchange third coordinate when needed
    exchangeCoordinates = 0d0

    do iEntry=1,nFaultSegmentEntriesGlobal
        ! (array, length, value, pos)
        call arrayContains(meshdatactx%gloInd, &
                           meshdatactx%nvLocal, &
                           friction(iEntry)%point(2), &
                           pos)


        if (pos.le.meshdatactx%nvLocal) then
            do iDOF=1,NDOF
                exchangeCoordinates(NDOF*(iEntry-1)+iDOF) = &
                     meshdatactx%X(iDOF,pos)
            enddo
        else
            ! this point does not belong to this partition
        endif

    enddo

    allocate(allCoordinates3(NDOF*nFaultSegmentEntriesGlobal))

#ifdef SPARSE
    call mpi_allreduce(exchangeCoordinates,    &  ! send
                       allCoordinates3,        &  ! receive
                       NDOF*nFaultSegmentEntriesGlobal, &  ! number of variables in array
                       MPI_double,             &  ! type
                       MPI_sum,                &  ! operator
                       MPI_COMM_WORLD,         &
                       iError)
#else
allCoordinates3 = exchangeCoordinates
#endif

endif


!****
! now compute segment area
!****
thickness = matdatactx%PRPMAT(6,1)  ! sixth entry of first material

do iEntry=1,nFaultSegmentEntriesGlobal
    if (NDOF.eq.2) then
        if (ISHELL.ne.0) then
            ! distance on the shell
            friction(iEntry)%sideArea = thickness * &
                distanceBetweenTwoPointsOnASphere(allCoordinates1(NDOF*(iEntry-1)+1), &
                                                  allCoordinates1(NDOF*(iEntry-1)+2), &
                                                  allCoordinates2(NDOF*(iEntry-1)+1), &
                                                  allCoordinates2(NDOF*(iEntry-1)+2), &
                                                  radius)
             if (debug) then
                write(*,*) "distance from point",allCoordinates1(NDOF*(iEntry-1)+1), &
                  allCoordinates1(NDOF*(iEntry-1)+2)," to ",&
                  allCoordinates2(NDOF*(iEntry-1)+1),allCoordinates2(NDOF*(iEntry-1)+2)
                write(*,*) "is", friction(iEntry)%sideArea / thickness
                write(*,*) "set sideArea (shell) to", friction(iEntry)%sideArea
            endif
        else
            ! cartesian distance in 2D
            friction(iEntry)%sideArea = thickness * &
                   distanceBetweenTwoPoints(allCoordinates1(NDOF*(iEntry-1)+1), &
                                            allCoordinates1(NDOF*(iEntry-1)+2), &
                                            0d0, &  ! z
                                            allCoordinates2(NDOF*(iEntry-1)+1), &
                                            allCoordinates2(NDOF*(iEntry-1)+2), &
                                            0d0)  ! adding two zeroes for z coords
             if (debug) then
                write(*,*) "distance from point",allCoordinates1(NDOF*(iEntry-1)+1), &
                  allCoordinates1(NDOF*(iEntry-1)+2)," to ",&
                  allCoordinates2(NDOF*(iEntry-1)+1),allCoordinates2(NDOF*(iEntry-1)+2)
                write(*,*) "is", friction(iEntry)%sideArea / thickness
                write(*,*) "set sideArea (plane) to", friction(iEntry)%sideArea
            endif
        endif
    else
        ! cartesian distance in 3D
        friction(iEntry)%sideArea = &
               areaOfTriangleIn3D(allCoordinates1(NDOF*(iEntry-1)+1), &
                                  allCoordinates1(NDOF*(iEntry-1)+2), &
                                  allCoordinates1(NDOF*(iEntry-1)+3), &
                                  allCoordinates2(NDOF*(iEntry-1)+1), &
                                  allCoordinates2(NDOF*(iEntry-1)+2), &
                                  allCoordinates2(NDOF*(iEntry-1)+3), &
                                  allCoordinates3(NDOF*(iEntry-1)+1), &
                                  allCoordinates3(NDOF*(iEntry-1)+2), &
                                  allCoordinates3(NDOF*(iEntry-1)+3))
    endif
enddo

deallocate(allCoordinates1)
deallocate(allCoordinates2)
if (NDOF.eq.3) then
    deallocate(allCoordinates3)
endif

! use the areas to translate the static Frictions and Winkler constants (in
! the friction struct) into strengths and differential winkler forces that act 
! on the (slippery) nodal points

do iEntry=1,nFaultSegmentEntriesGlobal
        nodalFriction = friction(iEntry)%Winkler * friction(iEntry)%sideArea
        nodalStrength = friction(iEntry)%staticFriction * friction(iEntry)%sideArea
        if (NDOF.eq.2) then        
        friction(iEntry)%NodalFrictionContribution(1) = nodalFriction *2d0 / 2d0
        friction(iEntry)%NodalFrictionContribution(2) = nodalFriction *2d0 / 2d0
        friction(iEntry)%NodalStrengthContribution(1) = nodalStrength / 2d0
        friction(iEntry)%NodalStrengthContribution(2) = nodalStrength / 2d0
        else
        friction(iEntry)%NodalFrictionContribution(1) = nodalFriction *2d0 / 3d0
        friction(iEntry)%NodalFrictionContribution(2) = nodalFriction *2d0 / 3d0
        friction(iEntry)%NodalFrictionContribution(3) = nodalFriction *2d0 / 3d0
        friction(iEntry)%NodalStrengthContribution(1) = nodalStrength / 3d0
        friction(iEntry)%NodalStrengthContribution(2) = nodalStrength / 3d0
        friction(iEntry)%NodalStrengthContribution(3) = nodalStrength / 3d0
        endif
enddo

! These frictions must be added to the WINX, IWINX and DIFFSTRENGTH vectors
! on the suitable place.

! In 3D, this does require a local euler angle rotation of the axes, so that you
! keep only two slippery degrees of freedom on your faults. It has to be like this:
! you cannot have opening or 'intrusive' faults in 3D (or in plane strain for
! that matter), and slip is only along the fault plane --> only two slippery DOF

! Not forget to multiply with a factor 2, as the differential Winkler constants
! as input should be twice the value imposed physically

do iEntry=1,nFaultSegmentEntriesGlobal

  do iPoint=1,2
    do iDOF=1,2

      if (modeldatactx%checkSignDuality(friction(iEntry)%point(iPoint)).eq.1111) then
       !!!! two degrees of freedom on slippery nodes
                   
           modeldatactx%IWINX(iDOF,friction(iEntry)%point(iPoint)) = 1 

           modeldatactx%WINX(iDOF,friction(iEntry)%point(iPoint)) = &
             modeldatactx%WINX(iDOF,friction(iEntry)%point(iPoint)) + &
             friction(iEntry)%NodalFrictionContribution(iPoint) *2.

           modeldatactx%DIFFSTRENGTH(friction(iEntry)%point(iPoint)) = &
             modeldatactx%DIFFSTRENGTH(friction(iEntry)%point(iPoint)) + &
             friction(iEntry)%NodalStrengthContribution(iPoint)

      else if (modeldatactx%checkSignDuality(friction(iEntry)%point(iPoint)).eq.1100) then 
       !!!! degree of freedom on x-slippery node

           modeldatactx%IWINX(1,friction(iEntry)%point(iPoint)) = 1
           modeldatactx%IWINX(2,friction(iEntry)%point(iPoint)) = 0 
           modeldatactx%WINX(1,friction(iEntry)%point(iPoint)) = &
             modeldatactx%WINX(1,friction(iEntry)%point(iPoint)) + &
             friction(iEntry)%NodalFrictionContribution(iPoint)
           modeldatactx%WINX(2,friction(iEntry)%point(iPoint)) = 0

           modeldatactx%DIFFSTRENGTH(friction(iEntry)%point(iPoint)) = &
             modeldatactx%DIFFSTRENGTH(friction(iEntry)%point(iPoint)) + &
             friction(iEntry)%NodalStrengthContribution(iPoint)

      else if (modeldatactx%checkSignDuality(friction(iEntry)%point(iPoint)).eq.11) then
       !!!! degree of freedom on y-slippery node

          modeldatactx%IWINX(1,friction(iEntry)%point(iPoint)) = 0
          modeldatactx%IWINX(2,friction(iEntry)%point(iPoint)) = 1 

          modeldatactx%WINX(1,friction(iEntry)%point(iPoint)) = 0
           modeldatactx%WINX(2,friction(iEntry)%point(iPoint)) = &
             modeldatactx%WINX(2,friction(iEntry)%point(iPoint)) + &
             friction(iEntry)%NodalFrictionContribution(iPoint)

           modeldatactx%DIFFSTRENGTH(friction(iEntry)%point(iPoint)) = &
             modeldatactx%DIFFSTRENGTH(friction(iEntry)%point(iPoint)) + &
             friction(iEntry)%NodalStrengthContribution(iPoint)
      endif
    enddo          
  enddo
enddo

if (debug) then
  do iEntry=1,nFaultSegmentEntriesGlobal
    do iPoint=1,2
      do iDOF=1,2
        i=modeldatactx%IWINX(iDOF,friction(iEntry)%point(iPoint))
        write(*,*) 'IWINX(idof=',idof,',node=',friction(iEntry)%point(iPoint), &
          ')=',i
        if (i.ne.0) write(*,*) 'WINX(idof=',idof,',node=',friction(iEntry)%point(iPoint), &
          ')=',modeldatactx%WINX(iDOF,friction(iEntry)%point(iPoint))
      enddo          
      write(*,*) 'DIFFSTRENGTH(node=',friction(iEntry)%point(iPoint), &
        ')=',modeldatactx%DIFFSTRENGTH(friction(iEntry)%point(iPoint))
    enddo
  enddo
endif

!deallocate(modeldatactx%checkSignDuality)

return

100 write(*,*) "Could not interpret record", trim(record)
    write(*,*) "Error: ", iostat
    stop "Leaving Gtecton"

110 write(*,*) "Reached premature end when interpreting record", trim(record)
    stop "Leaving Gtecton"

end subroutine

#endif
!-----------------------------------------------------------------------
!subroutine READMAXWELL (IWINK,WINK,FMAXW,NDOF)
subroutine READMAXWELL (NDOF)

use modeldatamodule, only: modeldatactx
USE LABELMODULE
USE MODELDEFINITION
use modelctx,        only: getrank
use debugmodule,     only: iecho, xit
use iomodule          
!
! Routine for reading and printing data on Maxwell boundary conditions
!

implicit none
!-pass
integer NDOF
!integer :: iWink
!double precision :: wink, fMaxw
!dimension IWINK(NDOF,*),WINK(NDOF,*),FMAXW(NDOF,*)
!-locl
logical constraint
character(len=13)fmt
integer :: m, mg,i, n, nn, max
!-init

!        write(0,*) 'Reading Maxwell boundary conditions'

m  = 0
mg = 0

fmt = '(2I12,3G14.0)'

if (iecho.eq.2) write(stderr,10) fmt
10   format(1x,'*** Read Maxwell bc magnitudes ***'/ &
 1x,'n=1; while (n != 0) {'/ &
 5x,'read ',A12,' node,increment,(FMAXW(i,n),i=1,NDOF)'/ &
 1x,'}')
!
max = NUMNP
!max = meshdatactx%nvglobal
call FGEN (modeldatactx%FMAXW,NDOF,max,fmt,NMAXW,"MaxwellBCS")

if (iecho.eq.2) then
    write(stderr,*) 'reading ', NMAXW, ' Maxwell boundary conditions'
endif


if (NMAXW.eq.0) then
    return
endif

write(*,*) 'Rank ', getrank(), ' reads Maxwell boundary conditions'


do n=1,NUMNP
    do i=1,NDOF
        if (modeldatactx%FMAXW(i,n).gt.0d0) then
            if (modeldatactx%IWINK(i,n).eq.0   .or. &
                modeldatactx%WINK(i,n) .le.0d0) then
                write(stderr,20) i,n
20                   format(1x,'RMAXW: Maxwell constraint on FMAXW(', &
                 I1,',',I5,'), not on corresponding Winkler ', &
                 'entry')
                call xit(1," ")
            endif
        endif
        if (modeldatactx%FMAXW(i,n).lt.0d0) then
            write(stderr,30) i,n
30               format(1x,'RMAXW: negative Maxwell constraint on ', &
             'FMAXW(',I1,',',I5,')')
            call xit(1," ")
        endif
    enddo
enddo

if (.not.OUTFIL(1)) return
!
nn = 51
do n=1,NUMNP
    if (nn.gt.50) then
        write(luout(1),60)
60           format(1X///2x, &
         'M A X W E L L   R E S T O R I N G   F O R C E S'//)
        write(luout(1),70)(LABELD(I),I=1,NDOF)
70           format(6x,'KEY TO APPLICATION MODE CODES:'// &
         10x,'MODE.EQ.0, NO MAXWELL RESTORING FORCES;',1x, &
         'THE NUMERICAL'/ &
         21x,'VALUE TO THE LEFT IS MEANINGLESS.'/ &
         10x,'MODE.GT.0, FORCES APPLIED IN STEP MODE-1'/ &
         10x,'MODE.LT.0, ELASTIC FORCES REMOVED IN STEP ',1X, &
         'IABS(MODE)'/// &
         4X,'NODE # ',A4,'-MAXWELL    MODE',5(2X,A4, &
         '-MAXWELL    MODE')/)
        write(luout(1),80)
80           format(/)
        nn = 1
    endif
    constraint = (modeldatactx%FMAXW(1,n).gt.0d0)
    do i=2,NDOF
        constraint = (constraint.or.(modeldatactx%FMAXW(i,n).gt.0d0))
    enddo
    if (constraint) then
        write(luout(1),90) n,(modeldatactx%FMAXW(i,n),modeldatactx%IWINK(i,n),i=1,NDOF)
90           format(1X,I8,2X,6(1PE12.5,2X,I5,3X))
        nn = nn + 1
    endif
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine writst(lunit,string)
!
! writes a string to logical unit
!
implicit none
!-pass
character         :: string*(*)
integer, external :: lnblk
integer           :: i, l
integer           :: lunit
!,adjustl
!
call adjustl(string)
l=lnblk(string)

if (l.gt.0) then
    write(lunit,1) (string(i:i),i=1,l)
 1       format(1x,80a)
    write(lunit,'(1x)')
endif

return
end
!-------------------------------------------------------------------------------
subroutine chknpar(NEN)

USE MODELCTX
USE MODELDEFINITION
use timestepmodule,  only: NINTG
USE MESHDATAMODULE,  only: meshdatactx
use algebra,         only: idig
use debugmodule,     only: xit
use iomodule          

!
implicit none
!
! Checks contents of NPAR parameters.
! Todo: which are those?
!
integer NEN

!-locl
logical :: srfrd
integer :: ierr
!
ierr = 0

if (meshdatactx%nvglobal.lt.1) then
    write(stderr,1)
 1       format(1x,'NUMEL < 1')
    ierr = 1
endif

if (NUMAT.lt.1) then
    write(stderr,2)
 2       format(1x,'NUMAT < 1')
    ierr = 1
endif

if (IPOINT.ne.1 .and. IPOINT.ne.NEN) then
    write(stderr,3) NEN
 3       format(1x,'IPOINT not equal 1 or ',I1)
    ierr = 1
endif

!       establish when surface nodes are required input
srfrd = (NINTG.gt.0 .and. (IDIG(ICVIS,4,3).ne.0 .or. NSED.ne.0))
srfrd = (srfrd .and. ((NTYPE.eq.1 .and. IOPT.eq.0) .or. NTYPE.eq.2) )
srfrd = (srfrd .and. NELAST.EQ.0 .and. IGRAV.EQ.0)

if (NSURF.le.0 .and. srfrd) then
    if (NSED.ne.0) then
        write(stderr,4)
 4           format(1x,'Sediment tracking requires NSURF>0')
    endif
    if (IDIG(ICVIS,4,3).ne.0) then
        write(stderr,21)
21           format(1x,'Density stripping requires NSURF>0')
    endif
    ierr = 1
endif

if (NUMSTR.lt.0) then
    write(stderr,5)
 5       format(1x,'NUMSTR < 0')
    ierr = 1
endif

if (NUMPR.lt.0) then
    if (NSURF.eq.0) write(stderr,6)
 6       format(1x,'NUMPR<0 & NSURF=0: no isostatic corrections')
    if (LGDEF.eq.0) write(stderr,*) 'NUMPR<0 and LGDEF=0: ', &
     'no isostatic corrections'
endif

if (NUMSLP.lt.0) then
    write(stderr,7)
 7       format(1x,'NUMSLP < 0')
    ierr = 1
endif

if (NUMFN.lt.0) then
    write(stderr,8)
 8       format(1x,'NUMFN < 0')
    ierr = 1
endif

if (INCOMP.lt.0) then
   write(stderr,9)
 9      format(1x,'INCOMP < 0')
   ierr=1
endif

if (LGDEF.ne.0 .and. IRESDU.lt.2) then
   write(stderr,10)
10      format(1x,'Warning: LGDEF "on" and IRESDU "off"')
endif

if (NSLSKEW.lt.0) then
    write(stderr,11)
11       format(1x,'NSLSKEW < 0')
    ierr=1
endif

if (NSLSKEW.gt.0) then
!    if (NUMROTglobal.eq.0) then
!        write(stderr,*) "Read",NSLSKEW, "Skewed angles, but number of rotations is zero"
!    endif

!    if (NUMSLPglobal.eq.0) then
!        write(stderr,*) "Read", NSLSKEW, "Skewed angles, but number of slippery node entries is zero"
!    endif
endif

if (NUMWNK.lt.0) then
    write(stderr,13)
13       format(1x,'NUMWNK < 0')
    ierr=1
endif

if (ierr.eq.1) then
    call xit(1,"BC parameters ")
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine RDPRE (NSTR)

USE MODELDATAMODULE
USE MESHDATAMODULE,  only: meshdatactx
USE MODELCTX,        only: getrank
USE MODELDEFINITION
USE LABELMODULE
use formatsmodule,   only: FMT_tecin_initial_stresses
use debugmodule,     only: iecho, xit, debug
use iomodule          
!
! program to read and print initial stress data
!     pre-stresses are read per element, per processor
!

implicit none
!-pass
integer NSTR
!-locl
integer ios
logical end,constraint
character (len=106)record !12 + 6 + 6*14
double precision :: f(6)
integer :: NumRd, LstRd, mg
integer :: i, it, n, nn
integer :: nel
!-init

!        write(*,*) 'rank', getrank(), &
!                   'Reading initial stresses', NPREglobal

mg    = 0
LstRd = 0
NumRd = 0

if (iecho.eq.2) then
    write(stderr,2) FMT_tecin_initial_stresses
 2   format(1x,'if (NPRE <> 0){'/ &
     5X,'*** read pre-stress data ***'/ &
     5X,'read ',A15,1x, &
      'modeldatactx%ISELM(n)'// &
      'modeldatactx%ISTIME(n)'// &
      '(modeldatactx%STN0(i,n,i=1,NSTR)' / &
     5X,'if (modeldatactx%ISELM(n)==0) return'/1x,'}')
endif


if (NPREglobal.eq.0) then
!            write(*,*) 'NPREglobal is 0, returning', NSTR
    return
endif
!


! read stress lines until one belonging to this rank is found...
100     call readrc (luin,record,'RDPRE',end)

if (debug) then
    write(*,*) "rdpre read record", record
endif

if (end) then
    if (debug) then
        write(*,*) 'rank',getrank(),'end found of prestress'
    endif
    goto 200
endif

if (debug) then
           write(*,*) 'rdpre read record: ', trim(record)
           write(*,*) 'interpreting with format: ', FMT_tecin_initial_stresses
endif

read (record,FMT_tecin_initial_stresses,iostat=ios,err=300,end=300) &
      nel,it, (f(i),i=1,NSTR)


if (debug) then
    write(*,*) "rdpre intepreted", nel,it, (f(i),i=1,NSTR)
endif

if (nel.gt.meshdatactx%neglobal) then
    write(*,*) 'Good heavens! Found prestress of element ', nel, &
               'But the model contains only ', meshdatactx%neglobal, &
               'elements.'
endif



if (nel.eq.0) then
    goto 200
endif

#ifdef SPARSE
!     check for owner of element nel; e2p is a global array (e2p = element 2 partition)
!        write(*,*) 'checking partition', nel, meshdatactx%e2p(nel), getrank()
if (meshdatactx%e2p(nel).ne.getrank()) then
!            write(*,*) 'element partition does not match'
    goto 100
endif
#endif


! ... yes, the stress line belongs tot his thread, do not go to 100, but continue.

NumRd = NumRd + 1

NPRE = NumRd

if (NumRd.gt.NPREglobal) then
    write(stderr,12) NPREglobal
12       format(1x,'rdpre: number of initial stresses exceeds ',I12)
    write(stderr,*) "Perhaps ''end'' statement missing?"
    call xit(1," ")
endif
!
do i=1,NSTR
!            write(*,*) 'accessing STN0 ', i,     'of', size(modeldatactx%STN0,1)
!            write(*,*) 'accessing STN0 ', Numrd, 'of', size(modeldatactx%STN0,2)
    modeldatactx%STN0(i,NumRd) = f(i)
enddo


!        write(*,*) 'accessing ISELM ', NumRd, &
!                   'of', size(modeldatactx%ISELM,1)        
modeldatactx%ISELM(NumRd)  = nel

!        write(*,*) 'accessing ISTIME ', NumRd, &
!                   'of', size(modeldatactx%ISTIME,1)
modeldatactx%ISTIME(NumRd) = it
!
LstRd = NumRd
!        write(*,*) 'finished a round of RDPRE'
goto 100   ! .. read next line


! read all the stress lines. Continued only from 'end' or from zero elems.

!        write(*,*) '*** rank', getrank(), '***', NPREglobal, NPRE, NumRd


! 200     ! write(*,*) 'rank ', getrank(), 'read all stresses: ', NumRd
!        if (NumRd.ne.NPRE) then
!            write(stderr,13) NumRd,NPRE
!   13       format(1x,'RDPRE: read ',I8,' entries, expected ',I8)
!            write(*,*) 'rank', getrank() , 'RDPRE: read ',NumRd,' entries, expected ',NPRE, 'exiting.'
!            call xit(1," ")
!        endif


! when necessary, write output to file.
200     if (.not.OUTFIL(1)) then
!            write(*,*) 'rank', getrank(),'returning from rdpre'
    return
endif
!
nn = 0
do n=1,NPRE

    constraint = .false.

    do i=1,NSTR
        constraint = (constraint .or. modeldatactx%STN0(i,n).ne.0.0)
    enddo

    if (constraint) then
        if (MOD(nn,50).eq.0) then

            write(luout(1),14) (LABELS(I),I=1,NSTR)
14               format(///1X,'I N I T I A L',3x,'S T R E S S',3x, &
             'D A T A'/// &
             6x,'KEY TO APPLICATION MODE CODES:'//10x, &
             'MODE.LE.0, NO PRE-STRESS; THE NUMERICAL'/ &
             21x,'VALUE TO THE LEFT IS MEANINGLESS.'/ &
             10x,'MODE.GT.0, PRE-STRESS APPLIED IN STEP MODE-1' &
             ///1x,'ELEM',5x,'MODE',9x,6(A3,13x))

            write(luout(1),16)
16               format(/)

        endif

        write(luout(1),18) modeldatactx%ISELM(n), &
                           modeldatactx%ISTIME(n), &
                          (modeldatactx%STN0(i,n),i=1,NSTR)
18           format(1x,i8,4x,i5,6(4x,1PE12.5))

        nn=nn+1
    endif

enddo
!

!write(*,*) 'rank ', getrank(), 'returning from rdpre 2'
return

300     write(*,*) 'rank ', getrank(), 'did a 300; error'
call prerr('rdpre')
write(stderr,301) ios,record

301     format(1x,'rdpre: IO error ',I5/1x,'last record=',a106)

!write(*,*) 'rank ', getrank(), 'returning from rdpre 3'

call xit(1," ")

end
!-------------------------------------------------------------------------------
subroutine FLXGEN (IELM,ISIDE,FLX,NFLUX,NFLXREAD,NDIM,FMT, callerName)

USE MODELCTX,        only: getrank
USE MESHDATAMODULE
USE MODELDEFINITION
use debugmodule,     only: xit, debug
use iomodule          
!
! Routine to read and generate floating point flux data
! (Also used by the thermal flux)
!
!       NFLUX            = number of elements with flux boundary condition
!       NDIM             = number of flux components on each element side
!       IELM (NFLUX)     = element numbers array
!       ISIDE (NFLUX)    = element side numbers array
!       FLX (NDIM,NFLUX) = element flux array
!

implicit none
!-pass
character FMT*(*)
integer NFLUX,NDIM
integer :: iElm, iSide
double precision :: flx
dimension IELM(*),ISIDE(*),FLX(NDIM,*)
!-locl
logical end
character(len=105)record
double precision  :: f(6)
integer, external :: lnblk
integer           :: NFLXREAD, numflx, lstflx, mg
character(len=*)  :: callerName
integer           :: i
integer           :: nel, is, ios
!-init
mg     = 0
numflx = 0
lstflx = 0
NFLXREAD = 0
!
if (NFLUX.le.0) then
    return
endif
!
100     call readrc (luin,record,'FLXGEN',end)


if (debug) then
    write(*,*) callerName, "has let flxgen read", trim(record)
endif

if (end) then
    if (debug) then
        write(*,*) "found 'end' statement"
    endif
    goto 200
endif

read (record,FMT,err=300,end=300,iostat=ios) nel,is, (f(i),i=1,NDIM)



if (debug) then
    write(*,*) callerName, "has let flxgen read", nel, is, (f(i),i=1,NDIM)
endif


if (nel.eq.0) then
    goto 200
endif

if (nel.gt.NUMEL) then
     write(stderr,9) nel
 9        format(1x,'flxgen: non-existing element ',I12)
     call exitp(1)
endif

!     only add element when element number nel matches this processor,
!     e2p is a global array, no need to look up local number of nel

if (meshdatactx%e2p(nel).ne.getrank()) then
    goto 100
endif


numflx = numflx + 1

if (numflx.gt.NFLUX) then
    write(stderr,6) NFLUX
 6       format(1x,'flxgen: number of flux bc''s exceeds ',I12)
    write(stderr,*) "Perhaps ''end'' statement missing?"

    call xit(1," ")
endif
!
do i=1,NDIM
    FLX(i,numflx) = f(i)
enddo

!        write(*,*) 'in between: ', numflx, NFLXREAD, NFLUX

IELM(numflx)  = nel
ISIDE(numflx) = is
!
lstflx = numflx
NFLXREAD = numflx


goto 100

!
200     if (numflx.ne.NFLXREAD) then
    write(stderr,7) numflx,NFLXREAD
 7       format(1x,'FLXGEN: read ',I12,' entries, expected ',I12)
    call xit(1," ")
endif



return
!
300     write(stderr,8) ios
 8   format(1x,'flxgen: fatal input error ',i8)
call prerr('flxgen')
call xit(1," ")

end
!-------------------------------------------------------------------------------
subroutine ELSTRS (IELSTR,ISSIDE,STRS,ISTR,NSTR)

USE MODELDATAMODULE
USE MODELDEFINITION, only : NUMSTR, NUMSTRglobal
USE MESHDATAMODULE,  only : meshdatactx
USE MODELCTX,        only : getrank
USE LABELMODULE
use debugmodule,     only : iecho, xit, debug
use iomodule          
!
! Subroutine to read and print out stress bc data
!
! traction lines are:
! element, sideID[1,2,3], Tstart, Tend, xx, yy, xy, [3 more z terms]
! stress increases from 0 at Tstart-1 to the full stress value on Tend,
! and will stay on that value for the rest of the simulation

implicit none
!-pass
integer NSTR
double precision :: STRS
integer          :: IELSTR,ISSIDE,ISTR
dimension STRS(6,*),IELSTR(*),ISSIDE(*),ISTR(2,*)
!-locl
logical lend
character(len=117)record
character(len=20) fmt
double precision :: f(6)

integer, external :: lnblk

integer :: side123
integer :: mg, numflx, lstflx, nread
integer :: i, n, nn
integer :: it0, it1, nel
integer :: ios
!-init
mg     = 0
numflx = 0
lstflx = 0

fmt = '(I12,I5,2I6,6G14.0)'
 
if (iecho.eq.2) then
    write(stderr,2) fmt
2   format(1x,'if (NUMSTR > 0){'/ &
     5X,'*** read element stress data ***'/ &
     5X,'read ',A19,1x, &
     'IELSTR(n),ISSIDE(n),ISTR(1,n),ISTR(2,n),', &
     '(STRS(i,n),i=1,6)'/5X,'if (IELSTR(n)==0) return'/1x,'}')
endif
 
if (NUMSTRglobal.le.0) then
    return
endif

nread = 0
NUMSTR = 0
 
100 call readrc (luin,record,'ELSTRS',lend)

    if (debug) then
        write(*,*) "Elstrs read record: ", trim(record)
        write(*,*) "with end: ", lend
    endif

    if (lend) then
        goto 200
    endif

    read (record,fmt,err=300,end=300,iostat=ios) nel,side123,it0,it1, &
     (f(i),i=1,6)

    if (debug) then
        write(*,*) "interpreted as: ", nel,side123,it0,it1, (f(i),i=1,6)
    endif

    if (nel.gt.meshdatactx%neglobal) then
        write(stderr,*) "Found stress of element ", nel, &
               "But the model contains only ", meshdatactx%neglobal, "elements."
    endif


    if (nel.eq.0) then
        goto 200
    endif

    if (it1.lt.it0) then
        it1 = it0
    endif

    !     check for owner of element nel, if owner is not this processor,
    !     skip otherwise
    nread = nread + 1
    if (meshdatactx%e2p(nel).ne.getrank()) then
        goto 100
    endif
 
    NUMSTR = NUMSTR+1
    numflx = numflx + 1

    if (numflx.gt.NUMSTRglobal) then
        write(stderr,14) NUMSTRglobal
14      format(1x,'ELSTRS: number of stress bc''s exceeds ',I8)
        write(stderr,*) "Perhaps ''end'' statement missing?"
        call xit(1," ")
    endif
 
    do i=1,6
        STRS(i,numflx) = f(i)
    enddo

    IELSTR(numflx) = nel
    ISSIDE(numflx) = side123
    ISTR(1,numflx) = it0
    ISTR(2,numflx) = it1
 
    lstflx = numflx
goto 100
 
200 continue

if (nread.ne.NUMSTRglobal) then
    write(stderr,15) nread,NUMSTRglobal
15  format(1x,'ELSTRS: read ',I8,' element entries, expected ', &
     I8)
    call xit(1," ")
endif
 
if (.not.OUTFIL(1)) return
 
nn = 0
do n=1,NUMSTR
    if (MOD(nn,50).eq.0) then
        write(luout(1),16) 'SXX','SYY','SZZ','SXY','SXZ','SYZ'
16      format(///1X,'S T R E S S',4X,'B O U N D A R Y',4X, &
         'C O N D I T I O N',4X,'D A T A'/// &
         6x,'KEY TO APPLICATION TIMES:'//10x, &
         'TIME.LT.0, NO STRESS BOUNDARY CONDITION;',1x, &
         'THE NUMERICAL'/ &
         21x,'VALUE TO THE LEFT IS MEANINGLESS.'/ &
         10x,'TIME.GE.0, STRESSES APPLIED IN STEP TIME'/// &
         4x,'ELEMENT',3x,'SIDE',5x,'TIME',7x,6(12x,A3))
        write(luout(1),18)
18      format(/)
    endif

    it0 = ISTR(1,n)
    it1 = ISTR(2,n)

    if (it0.eq.it1) then
        write(luout(1),19) IELSTR(n),ISSIDE(n),it0, STRS(1:6,n)
19      format(2x,i8,3x,i5,3x,i6,7x,6(1x,1PG14.4))
    else
        write(luout(1),20) IELSTR(n),ISSIDE(n),it0,it1, STRS(1:6,n)
20      format(2x,i8,3x,i5,3x,i6,'-',i6,6(1x,1PG14.4))
    endif

    nn=nn+1
enddo
 
return

300  write(stderr,22) ios
22   format(1x,'elstrs: fatal input error ',i8)
call prerr('elstrs')
call xit(1," ")

end subroutine
!-------------------------------------------------------------------------------
subroutine ELWINK ()

USE MODELDEFINITION, only : NUMWNK, NUMWNKglobal
USE MESHDATAMODULE,  only : meshdatactx
use modeldatamodule, only : modeldatactx
USE MODELCTX,        only : getrank
USE LABELMODULE
use debugmodule,     only : iecho, xit, debug
use formatsmodule,   only : FMT_tecin_WinklerRestoringPress
use iomodule          
!
! Subroutine to read and print out restoring pressure data
!
!       NUMWINK = number of elements with flux boundary condition
!       IWELM   = element numbers array
!       IWSIDE  = element side numbers array
!       IWTIME  = Winkler mode data
!       WPRES   = Winkler pressure values
!       ng      = element generation increment
!

implicit none

!-locl
logical end
character(len=49) record
external lnblk
!        integer getrank, nread
integer nread
double precision :: f
integer :: mg, numflx, lstflx
integer :: n, nn
integer :: nel, is, it
integer :: ios
!-init
mg     = 0
numflx = 0
lstflx = 0
nread = 0

!
if (iecho.eq.2) write(stderr,2) FMT_tecin_WinklerRestoringPress
 2   format(1x,'if (NUMWNK < > 0){'/ &
 5X,'*** read restoring pressure data ***'/ &
 5X,'read ',A17,1x, &
 'IWELM(n),IWSIDE(n),IWTIME(n),WPRES(n)'/ &
 5X,'if (IWELM(n)==0) return'/1x,'}')
!

if (NUMWNKglobal.eq.0) then
    return
endif
!

100     call readrc (luin,record,'ELWINK',end)

if (end) then
    goto 200
endif

read (record,FMT_tecin_WinklerRestoringPress,err=300,end=400,iostat=ios) nel,is,it,f
! element ID, side nr, onset time, magnitude

if (nel.gt.meshdatactx%neglobal) then
    write(*,*) 'Good heavens! Found Winkler pressure of element ', nel, &
               'but the model contains only ', meshdatactx%neglobal, 'elements.'
    goto 200
endif


if (nel.eq.0) then
    ! alternative 'end'
    goto 200
endif

nread = nread + 1 ! global counter, should match NUMWNKglobal

if (meshdatactx%e2p(nel).ne.getrank()) then
   goto 100
endif

numflx = numflx + 1
if (numflx.gt.IABS(NUMWNKglobal)) then
    write(stderr,14) IABS(NUMWNKglobal)
14       format(1x,'elwink: number of Winkler bc''s exceeds ',I12)
    write(stderr,*) "Perhaps ''end'' statement missing?"

    call xit(1," ")
endif
!
if (debug) then
    write(*,*) "getrank", getrank(), "read winkler restoring pressure: ", nel, is, it, f
endif

modeldatactx%WPRES(numflx)  = f
modeldatactx%IWELM(numflx)  = nel
modeldatactx%IWSIDE(numflx) = is
modeldatactx%IWTIME(numflx) = it
!
lstflx = numflx

goto 100
!
200     NUMWNK = numflx

if (nread.ne.IABS(NUMWNKglobal)) then
    write(stderr,15) nread,IABS(NUMWNKglobal)
15       format(1x,'elwink: read ',I12,' entries, expected ',I12)
    call xit(1," ")
endif




if (.not.OUTFIL(1)) then
    return
endif
!
nn = 0
do n=1,IABS(NUMWNK)
    if (MOD(nn,50).eq.0) then
        write(luout(1),16)
16           format(///1x,'R E S T O R I N G',3X,'P R E S S U R E', &
         3x,'D A T A'//6x,'KEY TO APPLICATION MODE CODES:'// &
         10x,'MODE.EQ.0, NO WINKLER RESTORING PRESSURES;',1x, &
         'THE NUMERICAL'/ &
         21x,'VALUE TO THE LEFT IS MEANINGLESS.'/ &
         10x,'MODE.EQ.-1, PRESSURES ALWAYS APPLIED'/ &
         10x,'MODE.GT.0 , PRESSURES APPLIED FROM STEP MODE-1 ', &
         'ONWARD'/ &
         10x,'MODE.LT.0 , PRESSURES APPLIED UNTIL STEP ', &
         'IABS(MODE)'/ &
         ///2x,'ELEMENT #',3x,'SIDE #',6x,'MODE',9x,'VALUE'/)
    endif
    write(luout(1),20) modeldatactx%IWELM(n),modeldatactx%IWSIDE(n),modeldatactx%IWTIME(n), &
     modeldatactx%WPRES(n)
20       format(3x,i8,4x,i5,4x,i6,4x,1PE12.5)
    nn=nn+1
enddo
!
return
!
300     write(stderr,22) ios
22   format(1x,'elwink: fatal input error ',i12)
call prerr('elwink')
call xit(1," ")

400     write(stderr,23) ios
23   format(1x,'elwink: premature end of list ',i12)
call prerr('elwink')
call xit(1," ")

end
!-------------------------------------------------------------------------------
subroutine FLXCHK (IELM,ISIDE,NFLUX,NEN,DATA)

USE MODELDEFINITION
use debugmodule,     only: xit
use iomodule          

!
! Routine to do some checking on flux data
!

implicit none
!-pass
integer NFLUX,NEN
character(len=*) DATA
integer   :: iElm, iSide
dimension IELM(*),ISIDE(*)
integer, external :: lnblk
!-local
integer :: i

if (NFLUX.le.0) then
    return
endif
!
do i=1,NFLUX

    if (IELM(i).le.0 .or. IABS(IELM(i)).gt.NUMEL) then
        write(*,*) "FLXCHK found flux entry ",i 
        write(*,*) "with element nr", IABS(IELM(i))
        write(*,*) "while it should be in interval < 0, nElems (", NUMEL,") ]"

        write(stderr,1) DATA(1:lnblk(DATA))
 1           format(1x,'flxchk: element number out of range in ',60a)
        call xit(1," ")
    endif

    if (ISIDE(i).le.0 .or. IABS(ISIDE(i)).gt.6) then
                write(*,*) "FLXCHK found flux entry ",i
        write(*,*) "with side nr", IABS(ISIDE(i))
        write(*,*) "while it should be in interval < 0,6 ] "

        write(stderr,2) DATA(1:lnblk(DATA))
 2           format(1x,'flxchk: element side number out of range in', &
         1x,60a)
        call xit(1," ")
    endif

enddo
!
return
end
!
!-------------------------------------------------------------------------------

subroutine shared (IEN,IJSUR,NUMEL,NEN,NSUR)
!
! sets up IJSUR array

use debugmodule, only: iecho
use iomodule,    only: stdout
!
implicit none
!-pass
integer :: NEN,NSUR,NUMEL
integer :: IEN(NEN,NUMEL),IJSUR(NSUR,NUMEL)
!local
integer :: n, nn
!
if (iecho.eq.1) then
    write(stdout,1,advance='no')
 1       format(1x,'Setting up element neighbours array ... ')
    call flushf(stdout)
endif
!
do nn=1,NUMEL
    n = nn
    call SHAREL (IEN,IJSUR,NUMEL,n)
enddo
!
!        if (iecho.eq.1) write(stderr,2)
!    2   format('done')
!
return
end
!-------------------------------------------------------------------------------
subroutine RDSURF (ISURF,NSURF)

USE MODELDATAMODULE
USE MESHDATAMODULE, only : meshdatactx
USE MODELCTX,       only : getrank
USE LABELMODULE
use debugmodule,     only: iecho, xit
use modeldefinition, only: NSURFglobal
use iomodule
!
! Reads surface nodal points data
!
implicit none
!-pass
integer :: iSurf
dimension ISURF(*)
integer nSurf
!-locl
logical end
character(len=24)record
character(len=6) fmt
integer, external :: lnblk
integer :: i, L, nread
integer :: ios
integer :: n, ng
!        integer getrank
!-init
L     = 0
nread = 0
NSURF = 0
fmt = '(2I12)'
!
if (iecho.eq.2) write(stderr,10) fmt
10   format('1x,if (NSURF > 0){'/ &
 5X,'*** read surface data ***'/ &
 5X,'read ',A5,1x,'ISURF(n),ng'/ &
 5X,'if (ISURF(n)==0 || ISURF(n)==end) return'/1x,'}')
!
if (NSURF.le.0) return
!
100     call readrc (luin,record,'RDSURF',end)
if (end) goto 400
read(record,fmt,err=600,end=700,iostat=ios) n,ng
if (n.eq.0) goto 400
!
nread = nread + 1
if (meshdatactx%v2p(n).ne.getrank()) goto 100
!
if (nread.gt.NSURFglobal) then
    write(stderr,4) NSURFglobal
 4           format(1x,'rdsurf: number input > ',I8)
    call xit(1," ")
endif

if (N.lt.L) then
    write(stderr,14) N, L
14           format(1x,'rdsurf: current node number ', I8, &
    ' < previous node number ',I8,'. Node numbers should be ', &
    'in ascending order.')
    call xit(1," ")
endif
L  = N
! increase local surface nodes by one
NSURF = NSURF + 1
! add global number to array (needed for postprocessing
! purposes)
modeldatactx%ISURFORDER(nsurf) = nread
! add node number to ISURF array
ISURF(nsurf) = n

goto 100
!
400     if (nread.ne.NSURFglobal) then
    write(stderr,5) nread,NSURFglobal
 5       format(1x,'rdsurf: generated ',I8,' data, expected ',I8)
    call xit(1," ")
endif
!
if (.not.OUTFIL(1)) return
!
n = 0
do i=1,NSURF
    if (MOD(n,50).eq.0) then
        write(luout(1),6)
 6           format(//1x, &
         'S U R F A C E',3x,'N O D E S'//5x,'NODE'/)
    endif
    n = n + 1
    write(luout(1),7) ISURF(i)
 7       format(1x,I8)
enddo
!
return
!
600     write(stderr,8) ios,NSURF,nread
 8   format(1x,'rdsurf: fatal read error ',i8,'.'/1x, &
 'Expected ',I8,' surface node entries, succesfully read ',i8)
call prerr('rdsurf')
call xit(1," ")
700     write(stderr,9) NSURF,nread
 9   format(1x,'rdsurf: fatal EOF error.'/1x, &
 'Expected ',I8,' surface node entries, succesfully read ',i8)
call xit(1," ")
end
!-------------------------------------------------------------------------------
subroutine LNKGEN (LINK,NLINK,NDOF,NUMNP)

USE LABELMODULE
USE FORMATSMODULE, only: FMT_tecin_dat_bcs_LinkedNodes
use debugmodule,   only: xit, debug
use iomodule          

implicit none
!-pass
integer NLINK,NDOF,NUMNP
integer :: LINK
dimension LINK(3,*)

!-locl
logical end
character(len=37) record  ! 3 * 12 +1 or the newline
character(len=1) cdof
external lnblk
integer :: LstNod, LstRd, NowRd, NowNod
integer :: MLink, idof
integer :: i, nn, nnn, mm
integer :: node1, node2
integer :: ios
!-init
LstNod = 0
LstRd  = 0
NowRd  = 0

if (NLINK.le.0) return

100 call readrc (luin,record,'LNKGEN',end)

if (end) then
    goto 400
endif

if (debug) then
    write(*,*) "linked nodes read record ", trim(record)
endif

read(record,FMT_tecin_dat_bcs_LinkedNodes,err=700,end=710,iostat=ios)  &
   idof,NowNod,MLink

if (debug) then
    write(*,*) "interpreted as ", idof,NowNod,MLink
endif

if (idof.eq.0) then
    goto 400
endif

NowRd = NowRd + 1

if (NowRd.gt.NLINK) then
    write(stderr,*) 'reading',NowRd,'linked node entries.'
    write(stderr,*) 'more than NLINK: ', NLINK
    call xit(1," ")
endif

LINK(1,NowRd) = idof
LINK(2,NowRd) = NowNod
LINK(3,NowRd) = MLink

LstNod = NowNod
LstRd  = NowRd
goto 100

400     if (NowRd.ne.NLINK) then
    write(stderr,70) NowRd,NLINK
70       format(1x,'LNKGEN: generated ',I12,' links, expected ',I12)
    call xit(1," ")
endif
!
! Some more checking
do i=1,NLINK
    idof  = LINK(1,i)
    node1 = LINK(2,i)
    node2 = LINK(3,i)

    if (node1.eq.node2) then
        write(stderr,80) i
80      format(1x,'LNKGEN: LINK(',I12,') links to itself ')
        call xit(1," ")
    else if (node1.lt.1 .or. node1.gt.NUMNP) then
        write(stderr,90) i
90           format(1x,'LNKGEN: LINK(',I12,') node number error')
        call xit(1," ")
    else if (node2.lt.1 .or. node2.gt.NUMNP) then
        write(stderr,90) i
        call xit(1," ")
    else if (idof.lt.1 .or. idof.gt.NDOF) then
        write(stderr,110) i
  110           format(1x,'LNKGEN: LINK(',I12,') DOF error')
        call xit(1," ")
    endif

enddo
!
if (.not.OUTFIL(1)) then
    return
endif

nnn = 51

do i=1,NLINK
!
! DOF     NODE             NODE
!
!  X     11234  EQUAL TO  11345
!
    if (nnn.gt.50) then
        write(luout(1),120)
  120           format(/1X/1X/1X,'L I N K E D',3X,'N O D A L',3X, &
         'D I S P L A C E M E N T S'/1X/1X,'DOF',5X,'NODE',13X, &
         'NODE'/)
        nnn = 0
    endif
    nnn = nnn + 1
    cdof = LABELC(LINK(1,i))(3:3)
    nn   = LINK(2,i)
    mm   = LINK(3,i)
    write(luout(1),130) cdof,nn,mm
  130       format(2X,A1,5X,I5,2X,'EQUAL TO',2X,I5)
enddo

return
!
700     write(stderr,4) ios
 4   format(1x,'LNKGEN: fatal input error ',i8)
call prerr('LNKGEN')
call xit(1," ")

710     write(stderr,5) ios
 5   format(1x,'LNKGEN: reached premature end of boundary conditions ',i8)
call prerr('LNKGEN')
call xit(1," ")


end
!-----------------------------------------------------------------------
subroutine readrc (lu,record,routine,end)
!
use debugmodule,     only: xit, debug
use iomodule,        only: stderr
use modelctx,        only: getrank

implicit none

!-pass
logical           :: end
integer           :: lu
character(len=*)  :: record,routine
!-locl
integer           :: f,l
integer, external :: fnblnk,lnblk
integer           :: i, ios

!-init
end = .false.
!
read(lu,1,err=100,end=200,iostat=ios) record
1   format(a)

if (debug) then
    write(*,*) "rank", getrank(), "readrc read record###", trim(record)
endif

f = fnblnk(record)
l = lnblk(record)

if (l.le.f) then
    end = .true.
else
    if (record(f:f+2).eq.'end') then
        do i=1,len(record)
            if (record(i:i).eq."0" .or. &
                record(i:i).eq."1" .or. &
                record(i:i).eq."2" .or. &
                record(i:i).eq."3" .or. &
                record(i:i).eq."4" .or. &
                record(i:i).eq."5" .or. &
                record(i:i).eq."6" .or. &
                record(i:i).eq."7" .or. &
                record(i:i).eq."8" .or. &
                record(i:i).eq."9") then
                write(*,*) "read record: ", trim(record)
                write(*,*) "end statement contains numeric entries"
                write(*,*) "Warning: perhaps linebreak after end statement missing?"
            endif
        enddo

        end = .true.
        record = ' '
    endif
endif

if (debug) then
    write(*,*) "rank", getrank(), " reached end? ", end
endif

return
 
100  write(stderr,2) ios,routine(1:lnblk(routine))
2    format(1x,'read error ',I4,' in routine ',10a)
call prerr(routine)
call xit(1," ")

200  write(stderr,3) ios,routine(1:lnblk(routine)),record
3    format(1x,'premature EOF error ',I4,' in routine ',10a, &
      'record', 72a)
call prerr(routine)
call xit(1," ")

end
!----------------------------------------------------------------------------
subroutine SETNEQBC (IBOND,IA,NDOF,NUMNP,MEQ)
!
! identifies equation numbers that are not actually generated (bcs)

implicit none

!-pass
integer :: NDOF,NUMNP,MEQ
integer :: IA(NDOF,NUMNP)
integer :: IBOND(NDOF,*)
!-locl
integer :: n,i
!
MEQ = 0
do n=1,NUMNP
     do i=1,NDOF
        if (IBOND(i,n).eq.1 .or. &
            IBOND(i,n).eq.2 .or. &
            IBOND(i,n).eq.5) then

            MEQ = MEQ + 1
            IA(i,n) = MEQ
        endif
    enddo
enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
!subroutine xit (i,message)
!!subroutine xit (i)

!USE MODELCTX
!use filemodule,  only: statusFileName, &
!                       do_soelim
!use debugmodule, only: lustat, &
!                       finalTime, &
!                       UserFinalize


!logical          :: there
!integer          :: i
!character(len=*) :: message

!inquire (file=statusFileName,exist=there)

!if (there .and. lustat.gt.0) then
!    call openf (lustat,statusFileName,'old')
!    if (i.eq.0) then
!                if (getrank().eq.0) then
!                write(lustat,10)
!10              format('Done')
!                endif
!    else
!        write(lustat,20) i
!20              format('Exit ',I2)
!    endif
!    call closef (lustat)
!endif


!if (i.ne.0) then
!    write(stderr,*) "Leaving application from rank ", getrank()
!    if (trim(message).ne." ") then
!        write(stderr,*) "With message: ", trim(message)
!    endif
!!    write(*,*) "Leaving application from rank ", getrank()
!!    if (trim(message).ne." ") then
!!        write(*,*) "With message: ", trim(message)
!!    endif
!endif

!! in case gtecton has soelimt in the /tmp directory,
!! delete this file
!if (do_soelim) then
!    close(unit=luin,status='delete')
!else
!    close(unit=luin)
!endif

!#ifdef SPARSE
!call UserFinalize()
!#endif

!call finalTime()

!call exitp(i)
!return
!end
!-------------------------------------------------------------------------------
subroutine updatelocalequationid(meshdat, modeldat,ndof)

USE MODELDATAMODULE, only: modeldata
USE MESHDATAMODULE,  only: meshdata
use modelctx,        only: getrank

implicit none

type (meshdata) :: meshdat
type (modeldata) :: modeldat

integer i,k,ndof

!     update local ID array to comply with IDglobal that might hold
!     slippery nodes

do i=1,meshdat%Nvlocal
    do k=1,NDOF
        modeldat%ID(k,i) = &
        modeldat%IDglobal(k,meshdat%gloInd(i))
    enddo
enddo

return
end

!-----------------------------------------
subroutine addlinkednodes(meshdat, modeldat)

USE MODELDATAMODULE
uSE MODELDEFINITION
USE MESHDATAMODULE
use debugmodule,     only: xit
use iomodule,        only: stderr
implicit none


type (meshdata) :: meshdat
type (modeldata) :: modeldat

integer n,i,k,master

!     linked nodes update in IDglobal array
!
if (NLINK.gt.0) then
    do n=1,NLINK
        i = modeldat%LINK(1,n)
        master = modeldat%LINK(3,n)
        k = modeldat%IDglobal(i,master)
        if (k.gt.0) then
            modeldat%IDglobal(i,modeldat%LINK(2,n)) = k
        else
            write(stderr,*) 'BC: linked node master ',master, &
           'generates no equation'
            call xit(1," ")
       endif
   enddo
endif

return
end
