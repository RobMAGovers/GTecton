MODULE FILEMODULE

character(len=255) :: PARTINFO
character(len=255) :: WORKPATH
logical            :: WORKPATHIsSet
character(len=255) :: OUTDIR
character(len=255) :: BaseFedsk
character(len=255) :: statusFileName
character(len=255) :: GTectonRCFileName

character(len=255) :: FullInputMechanical  ! combobulated TECIN file. Formerly known as FIN
character(len=255) :: FullInputThermal     ! combobulated TEMPIN file. Formerly also known as FIN

character(len=255) :: myTmpDir
character(len=255) :: FemDir               ! myTmpDir/FEM

integer            :: meshFileID           ! gmesh file, read by plt3d

integer, parameter :: NumbersInMB = 1048576/4   ! divided by four for a 4 bit number

! whether or not the tecin.dat contains .so inclusion lines.
! this requires more hard disk space but obsoletes a separate soelim phase before the run.
logical :: do_soelim

public :: index_filename
public :: readPathsFromEnvironment

contains


function index_filename (basename, ntimesteps, thistimestep, parallel)

#ifdef SPARSE
USE modelctx, only: getsize, getrank
#endif

implicit none


! this function build a fedsk filename based on a basename provided by the caller.
! It will either create a filename with both the timestep and the processor idea, when parallel == .true.
! or it will create a filename with only the timestep, when parallel == .false.
! The first is used by the pln and f3d routines to write their output.
! The second is used by plnplt and plt3d to read the merged data.

integer          :: ntimesteps, thistimestep
character(len=*) :: basename
character(len=100) :: index_filename
logical          :: parallel

integer          :: local_thistimestep
integer          :: inneri, innerj, logsize, i
#ifdef SPARSE
integer          :: nprocs, j, rnk, thisproc
#endif

local_thistimestep = thistimestep

index_filename = basename  ! set them equal, then append numbers to complete name
index_filename = trim(index_filename)//'.'

if (.not. parallel) then

    if (ntimesteps.ge.1) then
        logsize = log10(dble(ntimesteps)) + 1
    else
        logsize = 1
    endif

    do i=1,logsize
        inneri = 10**(logsize-i)
        innerj = local_thistimestep/inneri
        index_filename = trim(index_filename)//char(48+innerj)
        local_thistimestep = local_thistimestep - innerj * inneri
    enddo
    
else ! parallel .eq. .true.

#ifdef SPARSE
    rnk = getrank()

    nprocs = getsize()
    thisproc = getrank()

    do i=1,nprocs

        rnk = i-1
        logsize = log10(real(nprocs-1)) + 1

        do j=1,logsize
            inneri = 10**(logsize-j)
            innerj = rnk/inneri
            index_filename = trim(index_filename)//char(48+innerj)
            rnk = rnk - innerj * inneri
        enddo
    enddo
#endif

endif


end function
    
!--------------------------------------------------------------

subroutine readPathsFromEnvironment()

! both pln and plnplt require the knowledge of workpath and the
! temporary path, possible set by the user.

!use debugmodule, only: iecho

implicit none
logical :: there
integer :: length,status,csystem
character (len=1000) :: command
external csystem

!if (iecho.ne.0) then
!write(*,*) "In readpaths: WORKPATHIsSet: ", WORKPATHIsSet
!write(*,*) "on entering readpaths: workpath:", trim(workpath)
!endif

if (.not. WORKPATHIsSet) then
    call Get_Environment_Variable("PWD",WORKPATH,length,status)
else
    status = 0
endif

!if (iecho.ne.0) then
!    write(*,*) "filemodule, sub readPathsFromEnvironment "
!    write(*,*) "reads workpath with status", status
!    write(*,*) "and is now: ", trim(workpath)
!endif

if (status.eq.1) then
! PWD variable does not exist
  WORKPATH = '.'
else if (status.eq.0) then
! PWD variable succesfully loaded into WORKPATH variable
  continue
else if (status.eq.-1) then
  write(0,1)
1 format('Program error: getenv(PWD) does not fit into variable')
  WORKPATH = '.'
else
  write(0,2)
2 format('Program error: getenv(PWD) failed')
  WORKPATH = '.'
endif

call Get_Environment_Variable("MYTMP",myTmpDir,length,status)

if (status.eq.1) then
! MYTMP variable does not exist
  myTmpDir = WORKPATH//"/FEM"
else if (status.eq.0) then
! MYTMP variable succesfully loaded into myTmpDir variable
  continue
else if (status.eq.-1) then
  ! this is weird and unlikely to happen
  write(0,3)
3 format('Program error: getenv(MYTMP) does not fit into variable')
  myTmpDir = WORKPATH//'/FEM'
else
  write(0,4)
4 format('Program error: getenv(MYTMP) failed')
  myTmpDir = WORKPATH//'/FEM'
endif

inquire(file=trim(myTmpDir),exist=there)


if (.not.there) then
!  if (iecho.ne.0) then
!    write(*,*) "tmpdir does not exist", trim(myTmpDir), ". Making it..."
!  endif
  command = 'mkdir -p '// trim(myTmpDir)
  status = csystem (command)
  if (status.ne.0) then
    write(0,*) 'readPathsFromEnvironment: could not mkdir ' // trim(myTmpDir)
    call exitp(1)
  endif
else
  ! tmp directory already exists. Do nothing.
endif

end subroutine

!--------------------------------------------------------------

END MODULE
