subroutine setio()

USE MODELCTX
USE MODELDEFINITION
USE FILEMODULE         ! among which:  workpath, statusFileName, myTmpDir
use debugmodule,     only: lustat, &
                              itest, &
                              debug, &
                              iecho, &
                              delayInParallelLog, &
                              usedPetscOptionsTxt, xit, &
                              NicolaiNodeID, NicolaiFileID
use gitmodule, only: getGITtag
use formatsmodule,   only: tecindatnps_long
use petscksp
use test,         only: runTests
use iomodule


implicit none

#include "coord.i"
#ifdef SPARSE
#include "petsc/finclude/petscsys.h"
#endif


!-locl
character(len=180) arg
character(len=255) FECHO
character(len=24)  fdate
character(len=255) dirname
character(len=255)  string    
character(len=7)  gittag
logical there,loginfo
integer fnblnk, lnblk, csystem, ierr
character(len=1)   quote
external fnblnk, lnblk,dirname,csystem
logical outdir_specified

integer, external :: iflu, nextlu
integer :: narg, iarg
integer :: logsize
integer :: lutmp
integer :: m, n
integer :: rnk


#ifdef SPARSE
 PetscMPIInt irank

 call MPI_Comm_rank (PETSC_COMM_WORLD,irank,pierr)
#endif

!-init
stderr    = iflu('stderr')
stdin     = iflu('stdin')
stdout    = iflu('stdout')

luin      = nextlu(10)
luout(1)  = nextlu(luin+1)
luout(2)  = nextlu(luout(1)+1)
lustat    = nextlu(luout(2)+1)


if (debug) then
    write(*,*) 'set stderr   to ', stderr
    write(*,*) 'set stdin    to ', stdin
    write(*,*) 'set stdout   to ', stdout
    write(*,*) 'set luin     to ', luin
    write(*,*) 'set luout(1) to ', luout(1)
    write(*,*) 'set luout(2) to ', luout(2)
    write(*,*) 'set lustat   to ', lustat
endif



if (debug) then
    FILE_outputf     = nextlu(lustat+1)!+getrank()
    part_name = 'outputf.'
    rnk = getrank()
    logsize = log10(dble(getsize())) + 1

!        do i=1,logsize
!            inneri = 10**(logsize-i)
!            innerj = rnk/inneri
!            part_name = part_name(1:lnblk(part_name))//char(48+innerj)
!            rnk = rnk - innerj * inneri
!        enddo

!       write(*,*) 'part_name: ', part_name

endif

! Set the default for the various input parameters
! depending on command line arguments, those can still be 
! modified below. In case no arguments concerning that
! term have been provided, the defaults survive.

FEIN      = 'TECIN.DAT'
FEOUT     = 'FEOUT.DAT'
TEMPIN    = 'TEMPIN.DAT'
FEDSK     = 'FEDSK.DAT'
PARTINFO  = 'partition.info'
WORKPATHIsSet = .false.

! these variables are read only for pln, etc.
! plnplt needs them read 

FECHO     = ' '
OUTFIL(1) = .false.
OUTFIL(2) = .true.

iecho     = 0
debug     = .false.
itest     = 0
ascout    = .false.
quote     = CHAR(34)
loginfo   = .false.
optijd    = .false.
do_soelim = .true.
outdir_specified = .false.
tecindatnps_long = .false.
NicolaiNodeID = 0

narg  = COMMAND_ARGUMENT_COUNT()

if (narg.ge.1) then

  do iarg=1,narg
    ! keep in mind that the reading of the arguments
    ! is case sensitive
    call get_command_argument(iarg,arg)

    if (arg(1:4).eq.'echo') then
        n = INDEX(arg,'=')
        if (n.gt.0) then
            m = lnblk(arg)
            n = n+1
            if (arg(n:n).eq.'0') then
                iecho = 0
            else if (arg(n:n).eq.'1') then
                iecho = 1
            else if (arg(n:n).eq.'2') then
                iecho = 2      ! gives information about the reading of the mechanical input.
            else if (arg(n:n).eq.'3') then
                iecho = 3      ! echoes info about writing to FEDSK.DAT
            else if (arg(n:n).eq.'4') then
                iecho = 4      ! info thermal input file
            else if (arg(n:n).eq.'5') then
                iecho = 5 ! for overview of petsc interaction
                if (arg(n:n+1).eq.'55') then
                    iecho = 55 ! same as echo 5, put shows alllllllll petsc stuff.
                               ! this can be much. Every value inserted into matrix and vector, etc...
                               ! (note that at the time of writing this is not nearly done)
                endif
            else if (arg(n:n).eq.'6') then
                iecho = 6 ! for timing/profiling stuff
            else if (arg(n:n).eq.'7') then
                iecho = 7 ! for memory usage
            else if (arg(n:n).eq.'8') then
                if (arg(n:n+1).eq.'88') then
                    ! 88 inserts time breaks in the parallellog, getting a very easily readable
                    ! logfile, at the price of runtime.
                    delayInParallelLog = .true.
                else
                    delayInParallelLog = .false.
                endif

                iecho = 8 ! for thread monitoring (this will take a bit longer but logging will be detailed,
                          ! parallelLog in the debug module. See there for more documentation.
                
            else if (arg(n:n).eq.'9') then
                iecho = 9
                debug = .true.
            else
                write(stderr,10) 0,8
10              format(1x,'echo argument should be in',1x, &
                    'range ',I1,'-',I1)
                call xit(1," ")
            endif
        else
            iecho = 1
        endif

!temp fix for optijd flag, for setting residual to zero after restart

    elseif (arg(1:6).eq.'optijd') then
        optijd = .true.
    elseif (arg(1:4).eq.'noec') then
        iecho = 0
    elseif (arg(1:2).eq.'bi') then
        ascout = .false.
    elseif (arg(1:2).eq.'as') then
        ascout = .true.
    elseif (arg(1:4).eq.'noas') then
        ascout = .false.
    elseif (arg(1:4).eq.'nobi') then
        ascout = .true.
    elseif (arg(1:4).eq.'noso') then
        do_soelim = .false.
    elseif (arg(1:4).eq.'doso') then
        do_soelim = .true.
    elseif (arg(1:5).eq.'debug') then
        debug = .true.
    elseif (arg(1:4).eq.'test') then
        n = INDEX(arg,'=')
        if (n.gt.0) then
            m = lnblk(arg)
            n = n+1
            if (arg(n:n).eq.'0') then
                itest = 0
            else if (arg(n:n).eq.'1') then
                itest = 1
            else if (arg(n:n).eq.'2') then
                itest = 2
            else
                write(stderr,*) 'test=',arg(n:n),' not recognized'
                call xit(1," ")
            endif
        else
            itest = 2 ! run single time step and exit
        endif
    elseif (arg(1:5).eq.'feout' .or. arg(1:5).eq.'FEOUT') then
        if (getsize().gt.1) then
            OUTFIL(1) = .false.
            write(stderr, 14)
14                format(1x, 'feout not possible with multiple cpus')
        else
            OUTFIL(1) = .true.
            n = INDEX(arg,'=')
            if (n.gt.0) then
                m = lnblk(arg)
                n = n+1
                if (m.gt.n) then
                    FEOUT=arg(n:m)
                else
                    write(stderr,15)
15                        format(1x,'feout argument error')
                    call xit(1," ")
                endif
            endif
        endif
    elseif (arg(1:5).eq.'long') then
        tecindatnps_long = .true.
    elseif (arg.eq.'parout') then
        parout = .true.
    elseif (arg.eq.'nofeout') then
        OUTFIL(1) = .false.
    elseif (arg(1:8).eq.'workpath') then
        WORKPATHIsSet = .true.
        n = INDEX(arg,'=')
        if (n.gt.0) then
        m = lnblk(arg)
        n = n+1
        if (n.gt.0 .and. m.gt.n) then
            WORKPATH=arg(n:m)
        else
            write(stderr,18)
18                           format(1x,'workpath argument error')
            call xit(1," ")
        endif
        else
            write(stderr,18)
            call xit(1," ")
        endif
    elseif (arg(2:7).eq.'icolai') then ! n chopped off to avoid case sensitivity mistake
        n = INDEX(arg,'=')
        if (n.gt.0) then
            m = lnblk(arg)
            n = n+1
            if (n.gt.0 .and. m.ge.n) then
                read(arg(n:m),*) NicolaiNodeID
                NicolaiFileID=nextlu(42)
                open(unit=NicolaiFileID, file="Nicolai.log")
            else
                write(stderr,16)
16                           format(1x,'Nicolai monitor nodes')
                call xit(1," ")
            endif
        else
            write(stderr,16)
            call xit(1," ")
        endif

    elseif (arg(1:4).eq.'fein') then
        n = INDEX(arg,'=')
        if (n.gt.0) then
            m = lnblk(arg)
            n = n+1
            if (n.gt.0 .and. m.gt.n) then
                FEIN=arg(n:m)
                if (iecho.eq.2) then
                    write (*,*) "From arg: FEIN: ", trim(FEIN)
                endif
            else
                write(stderr,20)
20                               format(1x,'fein argument error')
                call xit(1," ")
            endif
        else
            write(stderr,20)
            call xit(1," ")
        endif
    elseif (arg(1:8).eq.'partinfo') then
        n = INDEX(arg,'=')
        if (n.gt.0) then
            m = lnblk(arg)
            n = n+1
            if (n.gt.0 .and. m.gt.n) then
                PARTINFO=arg(n:m)
!                        write (0,*) 'From arg: PARTINFO: ', PARTINFO
            else
                write(stderr,25)
22                       format(1x,'partinfo argument error')
                call xit(1," ")
            endif
        else
            write(stderr,25)
            call xit(1," ")
        endif
    elseif (arg(1:6).eq.'tempin') then
        n = INDEX(arg,'=')
        if (n.gt.0) then
            m = lnblk(arg)
            n = n+1
            if (n.gt.0 .and. m.gt.n) then
                TEMPIN=arg(n:m)
!                        write (0,*) 'From arg: TEMPIN: ', TEMPIN
            else
                write(stderr,25)
25                       format(1x,'tempin argument error')
                call xit(1," ")
            endif
        else
            write(stderr,25)
            call xit(1," ")
        endif
    elseif (arg(1:6).eq.'outdir') then
        n = INDEX(arg,'=')
        if (n.gt.0) then
            m = lnblk(arg)
            n = n+1
            if (n.gt.0 .and. m.gt.n) then
                OUTDIR=arg(n:m)
                outdir_specified = .true.
            else
                write(stderr,27)
27                       format(1x,'outdir argument error')
            endif
        else
            write(stderr,27)
            call xit(1," ")
        endif
    elseif (arg(1:5).eq.'fedsk') then
        n = INDEX(arg,'=')
        if (n.gt.0) then
            m = lnblk(arg)
            n = n+1
            if (n.gt.0 .and. m.gt.n) then
                FEDSK=arg(n:m)
!       WIENAND begin
                FEDSKRT = lnblk(FEDSK)
                if(FEDSK(FEDSKRT:FEDSKRT).eq.'/') then
                    fedskdir=fedsk
                    fedsk = fedsk(1:fedskrt) // 'FEDSK.DAT'
                    fedskrt = lnblk(fedsk)
                    string = 'mkdir -p '//fedskdir(1:lnblk(fedskdir))
                    ierr = csystem(string)
                else
                    fedskdir=dirname(fedsk)
                    if(fedskdir.ne.'.') then
                        string = 'mkdir -p ' // fedskdir(1:lnblk(fedskdir))
                        ierr = csystem(string)
                    endif
                endif
!       WIENAND end
            else
                write(stderr,30)
30              format(1x,'fedsk argument error')
                call xit(1," ")
            endif
        else
            write(stderr,30)
            call xit(1," ")
        endif

    elseif (arg(1:4).eq.'out=') then
        m = lnblk(arg)
        n = 5
        if (m.le.n) then
            write(stderr,32)
32                   format(1x,'out= argument error')
            call xit(1," ")
        endif
        FECHO = arg(n:m)
        stderr = nextlu(lustat+1)
        stdout = stderr

#ifdef SPARSE
        else if (arg.eq.'-log_info') then
            loginfo = .true.
#else
        else
            write(stderr,35) 'Unknown flag: ',quote, &
             arg(fnblnk(arg):lnblk(arg)),quote
35               format(1x,80A)
            write(stderr,35) 'Use ',quote,'man pln',quote, &
             ' for a list of options'
            call xit(1," ")
#endif

        endif
    enddo
endif

if (iecho.eq.1) then
    call getGITtag(gitTag)
    write(stdout,37)
37  format(' GTECTON version ',a7)
endif

call readPathsFromEnvironment()

! check whether the directory already exists and generate it when it is not.
if (getrank().eq.0) then
    inquire(file=myTmpDir, exist=there)

    if (.not. there) then
        ! -p = do not whine about extra parent directories you have to make
        ierr = csystem('mkdir -p '//trim(myTmpDir))
        if (ierr.ne.0) then
            write(stderr,*) "Could not generate temporary directory directory: ", myTmpDir
            write(stderr,*) "Leaving application."
            call xit(1," ")
        endif
    end if

endif


!*****************************************************************
! All command line arguments have been translated to variables.
! All the file variables that remain unset, will be set here


if(outdir_specified) then

!    write(*,*) 'outdir specified'

    OUTDIR = trim(WORKPATH)//'/'//trim(OUTDIR)

    ! check whether the outdit exists. If it does not, make it.
    ! This has to be done only once, so only the ero-thread does this

    if (getrank().eq.0) then
        inquire(file=OUTDIR, exist=there)

        if (.not. there) then
            ierr = csystem('mkdir -p '//trim(OUTDIR))
        end if

        if (ierr.ne.0) then
            write(stderr,*) 'Could not generate output directory. Leaving application.'
            call xit(1," ")
        endif
    endif
else
    OUTDIR = WORKPATH
endif

if (len(trim(workpath)) < 255) then  
    ! workpath has been set, prepend it.
    FEIN = trim(workpath) // '/' // trim(FEIN)
    TEMPIN = trim(workpath)//'/'//trim(TEMPIN)
    FEOUT = trim(outdir)//'/'//trim(FEOUT)
    FEDSK = trim(outdir)//'/'//trim(FEDSK)
    PARTINFO = trim(workpath)//'/'//trim(PARTINFO)
endif

if (debug) then
    part_name = trim(outdir)//'/'//trim(part_name)
endif

if (iecho.eq.2) then
    write (0,*) 'at the end of setting IO:'
    write (0,*) 'FEIN     : ', trim(FEIN)
    write (0,*) 'TEMPIN   : ', trim(TEMPIN)
    write (0,*) 'FEOUT    : ', trim(FEOUT)
    write (0,*) 'FEDSK    : ', trim(FEDSK)
    write (0,*) 'PARTINFO : ', trim(PARTINFO)
endif

FEDSKRT   = lnblk(FEDSK)

statusFileName    = 'STATUS'
GTectonRCFileName = 'GTecton.rc'

if (lnblk(FECHO).gt.0) then
call openf(stdout,trim(outdir)//'/'//FECHO,'unknown')
write(stdout,'(a24)') 'date not implemented'
endif

 if (.not.ascout .and. iecho.eq.3) then
!        write(stderr,40)
!   40       format(1x,'Results written to binary file')
iecho = 1
 else
!        write(stderr,41)
!   41       format(1x,'Results written to ascii file')
 endif

if (iecho.eq.2) then
    write(stderr,50)
50       format(1x,'*** input file format analysis ***')
elseif (iecho.eq.3) then
    write(stderr,55) '*** file ',quote,FEDSK(1:lnblk(FEDSK)), &
     quote,' contents analysis ***'
55       format(1x,80A)
elseif (iecho.eq.4) then
    write(stderr,55) '*** file ',quote, &
       TEMPIN(1:lnblk(TEMPIN)),quote,' contents analysis ***'
endif

! assure that only thread 0 writes to FEOUT.DAT
OUTFIL(1) = (OUTFIL(1).and.(getrank().eq.0))
if (OUTFIL(1)) call openf(luout(1),FEOUT,'unknown')

inquire(file=trim(workpath)//"/"//"BTOT.DAT ",exist=there)
if (there) then
!    write(*,*) 'Still a BTOT.DAT left from a previous run...'
endif

inquire(file=FEIN,exist=there)

if (.not.there) then
    if (iecho.eq.2) then
        write(*,*) "Could not find FEIN ", trim(FEIN)
    endif
    write(stderr,'(1x,''file not found '',80A)') &
     quote,FEIN(1:lnblk(FEIN)),quote
    call xit(1,"could not find input file")
else
    if (iecho.eq.2) then
        write(*,*) "Found FEIN: ", trim(FEIN)
    endif
endif

if (getrank().eq.0) then
     call openf (lustat,statusFileName,'unknown')
     write(lustat,*) 'Initializing'
60   format('Busy: 0')
     call closef (lustat)
endif

if (iecho.ne.0 .and. (.not.loginfo)) then
    write(stderr,80)
80  format(1x,'consider using the "-log_info" flag to see how', &
     1X,'PETSc performs')
endif

return
1000 write(*,*) "Could not write to file ",trim(GTectonRCFileName)
call xit(1," ")
 end
!-------------------------------------------------------------------------------
subroutine writeGTectonrc()

use modeldefinition, only: mode
use gitmodule,       only: getGIThash, getGITtag, getGITdate
use debugmodule,     only: usedPetscOptionsTxt
use fileModule,      only: GTectonRCFileName



#ifdef SPARSE
use spetscmodule,    only: M_SLES, T_SLES
#endif

implicit none

integer           :: GtectonRCFileHandle
integer           :: ierr

integer, external :: nextlu

character(len=40) :: gitHash
character(len=7)  :: gitTag
character(len=20) :: gitDate


double precision  :: relativeTol, absoluteTol, divergenceTol
integer           :: maxIterations


! for IX, IY and IZ
#include "coord.i"

GtectonRCFileHandle = nextlu(20)

open(unit = GtectonRCFileHandle, &
     file = trim(GTectonRCFileName), &
     iostat = ierr, action = 'write')

if (ierr.ne.0) then
    write(*,*) "Could not open file GTecton.rc to store build configuration data"
    stop "Leaving GTECTON..."
endif

! call generated ../src/TempSrc/getGITdata.F in its own shared object to return the hash.
call getGIThash(gitHash)
call getGITtag(gitTag)
call getGITdate(gitDate)

!--------------

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "These binaries were compiled with the following version of the source code:"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "Version : " // gitTag
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "GIT hash: " // gitHash
write(GtectonRCFileHandle, "(a)") "from"
write(GtectonRCFileHandle, "(a)") "GIT date: " // gitDate
write(GtectonRCFileHandle, "(a)") ""

!--------------

write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "EXPERIMENTAL_CYCLIC allows periodic and non-periodic earthquake cycles"
write(GtectonRCFileHandle, "(a)") ""

write(GtectonRCFileHandle, "(a)", advance="no") "EXPERIMENTAL_CYCLIC: "
#ifdef EXPERIMENTAL_CYCLIC
 write(GtectonRCFileHandle,*) '.true. '
#else
 write(GtectonRCFileHandle,*) '.false.'
#endif
write(GtectonRCFileHandle, *) "    "

!--------------

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "ZISFAULTNORMAL is relevant when EXPERIMENTAL_CYCLIC=.true."
write(GtectonRCFileHandle, "(a)") ""

#ifdef EXPERIMENTAL_CYCLIC
write(GtectonRCFileHandle, "(a)", advance="no") "ZISFAULTNORMAL: "
#ifdef ZISFAULTNORMAL
 write(GtectonRCFileHandle,*) '.true.: Euler angles need to align fault slip parallel to x&y directions, and'
 write(GtectonRCFileHandle,*) '        the y(2d) or z(3d) entry of IWINX is used to specify the cycle phase'
#else
 write(GtectonRCFileHandle,*) '.false.'
#endif
#endif
write(GtectonRCFileHandle, *) "    "

!--------------

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "LINPL for stress dependence in powerlaw viscosity calculation"
write(GtectonRCFileHandle, "(a)") ""

write(GtectonRCFileHandle, "(a)", advance="no") "LINPL: "
#ifdef LINPL
 write(GtectonRCFileHandle,*) '.true.: linearized powerlaw stress dependence'
#else
 write(GtectonRCFileHandle,*) '.false.: powerlaw stress dependence'
#endif
write(GtectonRCFileHandle, *) "    "

!--------------

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "EXPERIMENTAL_ReducedSlipAlloc is used to store the slippery nodes more efficiently."
write(GtectonRCFileHandle, "(a)") ""

write(GtectonRCFileHandle, "(a)", advance="no") "EXPERIMENTAL_ReducedSlipAlloc: "
#ifdef EXPERIMENTAL_ReducedSlipAlloc
 write(GtectonRCFileHandle,*) '.true. '
#else
 write(GtectonRCFileHandle,*) '.false.'
#endif
write(GtectonRCFileHandle, *) "    "

!---------------


#ifdef SHEARZONE_TRACTIONS
#ifdef EXPERIMENTAL_FRICTION
write(stderr,*) 'This version of GTECTON was compiled with -DSHEARZONE_TRACTIONS ',&
 'and -DEXPERIMENTAL_FRICTION'
write(stderr,*) 'These compile flags are mutually incompatible'
call exitp(0)
#endif
#endif

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "SHEARZONE_TRACTIONS is used for velocity dependent tractions on slippery node interfaces."
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)", advance="no") "SHEARZONE_TRACTIONS: "
#ifdef SHEARZONE_TRACTIONS
 write(GtectonRCFileHandle,*) '.true. '
#else
 write(GtectonRCFileHandle,*) '.false.'
#endif
write(GtectonRCFileHandle, *) " "

!--------------------------------

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "EXPERIMENTAL_FRICTION is used for stick-slip for a constant kinematic traction."
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)", advance="no") "EXPERIMENTAL_FRICTION: "
#ifdef EXPERIMENTAL_FRICTION
 write(GtectonRCFileHandle,*) '.true. '
#else
 write(GtectonRCFileHandle,*) '.false.'
#endif
write(GtectonRCFileHandle, *) " "



!-------------

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "LongData: mantissa size of ascii floating point data in FEDSK files"
write(GtectonRCFileHandle, "(a)") ""

write(GtectonRCFileHandle, "(a)", advance="no") "LongData was set to: "
#ifdef LongData
 write(GtectonRCFileHandle,*) '.true.: E25.16'
#else
 write(GtectonRCFileHandle,*) '.false.: E15.6'
#endif

write(GtectonRCFileHandle, "(a)") ""

!--------------

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "DIFORC_SET hack for convenience of some users"
write(GtectonRCFileHandle, "(a)") ""

write(GtectonRCFileHandle, "(a)", advance="no") "DIFORC_SET = "
#ifdef DIFORC_SET
 write(GtectonRCFileHandle,*) '.true.: sets the differential force while reading input'
#else
 write(GtectonRCFileHandle,*) '.false.: adds differential force while reading input'
#endif

write(GtectonRCFileHandle, *) " "

!-------------

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "MUMPS option to use external MUMPS solver in stead of PETSc"
write(GtectonRCFileHandle, "(a)") ""

write(GtectonRCFileHandle, "(a)", advance="no") "MUMPS: "
#ifdef MUMPS
 write(GtectonRCFileHandle,*) '.true.: uses external MUMPS. Its location has to be included in CMakeLists.txt'
#else
 write(GtectonRCFileHandle,*) '.false.: uses default PETSc'
#endif

write(GtectonRCFileHandle, *) " "

!-------------

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "TIMING returns profiling information at the cost of run time."
write(GtectonRCFileHandle, "(a)") ""

write(GtectonRCFileHandle, "(a)", advance="no") "TIMING: "
#ifdef TIMING
 write(GtectonRCFileHandle,*) '.true. '
#else
 write(GtectonRCFileHandle,*) '.false. '
#endif

write(GtectonRCFileHandle, *) " "

!-------------

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, *)     "Coordinate directions : ix=",IX," iy=",IY," iz=",IZ
write(GtectonRCFileHandle, "(a)") ""

!--------------

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
#ifdef SPARSE 
!otherwise there us no PETSC stuff

write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)") "Extra PETSc parameters can be read from file using the argument 'usedPetscOptionsTxt'"
write(GtectonRCFileHandle, "(a)") ""
write(GtectonRCFileHandle, "(a)", advance="no") "usedPetscOptionsTxt: "

if (usedPetscOptionsTxt) then
 write(GtectonRCFileHandle,*) '.true. '
else
 write(GtectonRCFileHandle,*) '.false. '
endif
write(GtectonRCFileHandle, "(a)") ""

if (mode.eq.2 .or. mode.eq.6) then
    !there is a mechanical component
    call KSPGetTolerances(M_SLES, relativeTol, absoluteTol, divergenceTol, maxIterations, ierr)

    if (ierr.ne.0) then
        write(0,*) "Warning: Setio could not retrieve mechanical solver tolerances, error: ", ierr
    else
        write(GtectonRCFileHandle, "(a)") "PETSc parameters used for the mechanical solve:"
        write(GtectonRCFileHandle, "(a)") ""
        write(GtectonRCFileHandle, *) "Relative   tolerance: ", relativeTol
        write(GtectonRCFileHandle, *) "Absolute   tolerance: ", absoluteTol
        write(GtectonRCFileHandle, *) "Divergence tolerance: ", divergenceTol
        write(GtectonRCFileHandle, *) "Maximum number of iterations: ", maxIterations
        write(GtectonRCFileHandle, "(a)") ""
    endif
endif

if (mode.eq.5 .or. mode.eq.6) then
    !there is a thermal component
    call KSPGetTolerances(T_SLES, relativeTol, absoluteTol, divergenceTol, maxIterations, ierr)

    if (ierr.ne.0) then
        write(0,*) "Warning: Setio could not retrieve thermal solver tolerances, error: ", ierr
    else
        write(GtectonRCFileHandle, "(a)") "PETSc parameters used for the thermal solve:"
        write(GtectonRCFileHandle, "(a)") ""
        write(GtectonRCFileHandle, *) "Relative   tolerance: ", relativeTol 
        write(GtectonRCFileHandle, *) "Absolute   tolerance: ", absoluteTol 
        write(GtectonRCFileHandle, *) "Divergence tolerance: ", divergenceTol 
        write(GtectonRCFileHandle, *) "Maximum number of iterations: ", maxIterations 
        write(GtectonRCFileHandle, "(a)") ""
    endif

endif

write(GtectonRCFileHandle, "(a)") "--------------------------------------------------------------------------"
#endif


!--------------

close(GtectonRCFileHandle)

end subroutine

!-------------------------------------------------------------------------------


subroutine OPNDSK (openlu,io)

USE TIMESTEPMODULE,  only: NINTG,NMPRT,NTPRT,MAXITIME, &
                           nstep, ITWORK,IMWORK,IMTWRK,refact

USE MODELDEFINITION
USE MODELCTX
use modeldatamodule, only: HasSlipperyNodes
use modeltopology,   only: nen, nsd, ndof, nstr
use algebra,         only: idig
use debugmodule,     only: xit
use iomodule,        only: ascout, luout, outfil, stderr, fedskrt, fedsk

!
implicit none
!-pass
 logical openlu
integer :: io ! to pass the return value back.
!    dimension RATE(NSTR,*),SURF(NSD,*),FLUX(NSD,*)

!-local
character(len=1)  :: a
logical           :: there
integer           :: nmech,ntemp,mxint
integer           :: rnk, logsize, inneri, innerj
integer, external :: INMAX,NINT
integer, external :: lnblk
integer           :: i, ierr
integer           :: ncntr, ndsk
integer           :: nmstp, ntstp
double precision  :: xdsk
!-init
 mxint = INMAX()
!
 if (.not. OUTFIL(2)) return
!
!     in case of a restart, no need to open fedsk output files
!     since output is written in one file per time step

!w    if (ISTART.lt.0) then
!w        inquire(file=FEDSK,exist=there)
!w        if (.not.there) then
!w        write(stderr,1) '"',FEDSK(1:lnblk(FEDSK)),'"'
!w   1           format(1x,'RESTART- file not found: ',258a)
!w        call xit(1," ")
!w        endif
!w    endif
!
#ifndef __ia64
!    estimate prospective FEDSK disk usage
 ndsk = 0

if (NINTG.gt.0.and. &
   (MODE.eq.2.or.MODE.eq.6)) then
    if (NMPRT.ge.0) then
        nmech = NMPRT
    else
        nmech = MAXITIME/IABS(NMPRT)
    endif
else
    nmech = 0
endif

if (NINTG.gt.0.and. (MODE.eq.5.or.MODE.eq.6)) then
    if (NTPRT.ge.0) then
        ntemp = NTPRT
    else
        ntemp = MAXITIME/IABS(NTPRT)
    endif
else
    ntemp = 0
endif

if (MODE.eq.2.or.MODE.eq.6) then
    !  Elastic data block
    call INTINCR (ndsk,NUMNP*NDOF,IERR) ! displacements
    if (IERR.ne.0) goto 800
    if (HasSlipperyNodes()) then
        call INTINCR (ndsk,NUMEL*NEN*NDOF,IERR) ! slippery
        if (IERR.ne.0) goto 800
    endif
    call INTINCR (ndsk,NUMEL*NSTR,IERR) ! stresses
    if (IERR.ne.0) goto 800
    call INTINCR (ndsk,NUMEL*NSTR,IERR) ! strains
    if (IERR.ne.0) goto 800
endif

if (MODE.eq.5.or.MODE.eq.6) then
    !        Temperature initialization data block
    call INTINCR (ndsk,NUMNP,IERR) ! temperatures
    if (IERR.ne.0) goto 800
    if (NUMFN+NUMSLP.gt.0 .and. IDIFT.eq.1) then
        call INTINCR (ndsk,NUMNP,IERR) ! differential temps
        if (IERR.ne.0) goto 800
    endif
    call INTINCR (ndsk,NUMEL*NSD,IERR) ! fluxes
    if (IERR.ne.0) goto 800
endif

nmstp = 0 ! disk space for mechanical timesteps
ntstp = 0 ! disk space for thermal timesteps

if (MODE.eq.6.and.IDIG(ICVIS,4,1).gt.0) then
    call INTINCR (nmstp,2*NUMEL,IERR) ! viscosities
    if (IERR.ne.0) goto 800
endif

if (MODE.eq.6.and.IDIG(ICVIS,4,2).gt.0) then
    call INTINCR (nmstp,2*NUMEL,IERR) ! plastic
    if (IERR.ne.0) goto 800
endif

call INTINCR (nmstp,NUMNP*NDOF,IERR) ! displs
 
if (IERR.ne.0) goto 800
 
!if (IVELOUT.ne.0) then
call INTINCR (nmstp,NUMNP*NDOF,IERR) ! velos
if (IERR.ne.0) goto 800
!endif

if (HasSlipperyNodes()) then
    call INTINCR (nmstp,NUMEL*NEN*NDOF,IERR)
    if (IERR.ne.0) goto 800
!if (IVELOUT.ne.0) then
    call INTINCR (nmstp,NUMEL*NEN*NDOF,IERR)
    if (IERR.ne.0) goto 800
!endif
endif

call INTINCR (nmstp,NUMEL*NSTR,IERR) ! stresses
if (IERR.ne.0) goto 800

call INTINCR (nmstp,NUMEL*NSTR,IERR) ! strains
if (IERR.ne.0) goto 800

call INTINCR (nmstp,NUMEL*NSTR,IERR) ! rates
if (IERR.ne.0) goto 800

if (NSURF.gt.0.and.NSED.ne.0) then
    call INTINCR (nmstp,NSURF*NSD,IERR) ! surfs
    if (IERR.ne.0) goto 800
endif

if (IADVEC.ne.0) then
    call INTINCR (ntstp,NUMNP*NDOF,IERR) ! velos
    if (IERR.ne.0) goto 800
endif

call INTINCR (ntstp,NUMNP,IERR) ! temps

if (IERR.ne.0) goto 800

if (NUMFN+NUMSLP.gt.0 .and. IDIFT.eq.1) then
    call INTINCR (ntstp,NUMNP,IERR) ! diftemps
    if (IERR.ne.0) goto 800
endif

call INTINCR (ntstp,NUMEL*NSD,IERR) ! fluxes
if (IERR.ne.0) goto 800


!if (DBLE(nmech)*DBLE(nmstp).gt.DBLE(mxint)) goto 900

!call INTINCR (ndsk,nmech*nmstp,IERR)
!if (IERR.ne.0) goto 900

!if (DBLE(ntemp)*DBLE(ntstp).gt.DBLE(mxint)) goto 900
!call INTINCR (ndsk,ntemp*ntstp,IERR)

!if (IERR.ne.0) goto 900

!if (ndsk.ge.mxint/(2*IWORD)) goto 900 ! see notes



!
!    Check here for availability of disk space
!    write (0,*) 'Established disk requirement: ', ndsk, ' MB.'
xdsk = DBLE(ndsk)/DBLE(2**20)
ndsk = NINT(xdsk)
!    if (iecho.eq.1) write(stdout,2) ndsk
!    2   format(1x,'"FEDSK.DAT" requires ',I8,' Mbytes of disk space')
#endif
!
if (ISTART.ge.0) then
!       WIENANDf begin
!       now adapt the FEDSK filename such that is included some counterto indicate the time stepping.
!       For consequence, the length of this counter is determined to be the same throughtout the simulation,
!       to prevent things like:

!         fedsk.par.0
!         fedsk.par.100
!         fedsk.par.200

!       but instead get:

!         fedsk.par.000
!         fedsk.par.100
!         fedsk.par.200

    if (MAXITIME.ge.1) then
        ncntr = int(log10(dble(MAXITIME))) + 1
    else
        ncntr = 1
    endif

    FEDSK = FEDSK(1:FEDSKRT)//'.'

    do i=1,ncntr
        FEDSK = FEDSK(1:FEDSKRT+i)//char(48)
      enddo

  ! next add a counter for the rank

    FEDSK = FEDSK(1:FEDSKRT+ncntr+1)//'.'
    rnk = getrank()

    if (getsize().gt.1) then
        logsize = log10(dble(getsize())) + 1
    else ! getsize = 0
        logsize = 1 ! a zero takes one character
    endif

    do i=1,logsize
        inneri = 10**(logsize-i)
        innerj = rnk/inneri
        FEDSK = FEDSK(1:FEDSKRT+ncntr+1+i)//char(48+innerj)
        ! remove last digit to proceed to one more significant digit
        rnk = rnk - innerj * inneri
    enddo


!          if (getsize().gt.1) then
!              FEDSK = FEDSK(1:FEDSKRT+ncntr+1)//'.'
!              rnk = getrank()
!              logsize = log10(dble(getsize())) + 1
!              do i=1,logsize
!                  inneri = 10**(logsize-i)
!                  innerj = rnk/inneri
!                  FEDSK = FEDSK(1:FEDSKRT+ncntr+1+i)//char(48+innerj)
!                  rnk = rnk - innerj * inneri
!              enddo
!          endif

!        write(*,*) 'Thread ', getrank(),' of ',getsize(),' says:'
!        write(*,*) 'FEDSK after counter: ', trim(FEDSK)

!       WIENAND end

    if (openlu) then
        if (ascout) then
            call openf (luout(2),FEDSK,'unknown')
        else
            call fbopen (luout(2),FEDSK,'unknown')
        endif
    endif


!*************************************************************************
!*************************************************************************
!          Here we write the header in the fedsk output files!
!*************************************************************************
!*************************************************************************

!          if (ascout) then
!              write(luout(2),3,err=100,iostat=io) VERSION
!    3           format(1x,'G-TECTON ',A15)
!          else
!              write(luout(2),err=200,iostat=io) VERSION
!          endif




!      else
!        if (ascout) then
!        call lueof(luout(2))
!        else
!10        read(luout(2),end=20,err=300,iostat=io) a
!        goto 10
!20        backspace(luout(2),err=400,iostat=io)
!        if (MODE.eq.2 .or. MODE.eq.6) then
!            if (NSURF.le.0) then
!            read(luout(2),err=500,iostat=io)
!     >             ((RATE(i,n),i=1,NSTR),n=1,NUMEL)
!            else
!            read(luout(2),err=600,iostat=io)
!     >             ((SURF(i,n),i=1,NSD),n=1,NSURF)
!            endif
!        else if (MODE.eq.5) then
!            read(luout(2),err=700,iostat=io)
!     >             ((FLUX(i,j),i=1,NSD),j=1,NUMEL)
!        endif
!        endif
endif
!
return
!

100    write(stderr,101) io,'"',FEDSK(1:lnblk(FEDSK)),'"'
  101   format(1x,'formatted write error ',I5,' on ',258a)
 if (.not.openlu) return
 call xit(1," ")

200    write(stderr,201) io,'"',FEDSK(1:lnblk(FEDSK)),'"'
  201   format(1x,'unformatted write error ',I5,' on ',258a)
 if (.not.openlu) return
 call xit(1," ")

300    write(stderr,301) io,'"',FEDSK(1:lnblk(FEDSK)),'"'
  301   format(1x,'unformatted read error ',I5,' on ',258a)
 if (.not.openlu) return
 call xit(1," ")

400    write(stderr,401) io,'"',FEDSK(1:lnblk(FEDSK)),'"'
  401   format(1x,'unformatted backspace error ',I5,' on ',258a)
 if (.not.openlu) return
 call xit(1," ")

500    write(stderr,501) io,'"',FEDSK(1:lnblk(FEDSK)),'"'
  501   format(1x,'unformatted read RATE error ',I5,' on ',258a)
 if (.not.openlu) return
 call xit(1," ")

600    write(stderr,601) io,'"',FEDSK(1:lnblk(FEDSK)),'"'
  601   format(1x,'unformatted read SURF error ',I5,' on ',258a)
 if (.not.openlu) return
 call xit(1," ")

700    write(stderr,701) io,'"',FEDSK(1:lnblk(FEDSK)),'"'
  701   format(1x,'unformatted read FLUX error ',I5,' on ',258a)
 if (.not.openlu) return
 call xit(1," ")

800    write(stderr,801) NINT(DBLE(mxint)/DBLE(2**29))
  801   format(1x,'Anticipated size of "FEDSK.DAT" larger than ', &
 I2,' Gbytes.')
 call xit(1," ")

900    write(stderr,901) NINT(DBLE(mxint)/DBLE(2**29)), &
   nmech,nmstp*4*IWORD/2**20,ntemp,ntstp*4*IWORD/2**20
  901   format(1x,'Anticipated size of "FEDSK.DAT" larger than ', &
 I2,' Gbytes.'/1x,I5,' mechanical outputs of ',I4, &
 ' MBytes each'/1x,I5,' thermal outputs of ',I4,' MBytes each')
 call xit(1," ")

 end
!-------------------------------------------------------------------------------
 subroutine INTINCR (INTGR,INCREMENT,IERR)
!
!    increments an integer when there is a danger of overflow
!
!-pass
 integer INTGR,INCREMENT,IERR
!-locl
 integer imx,INMAX
 save imx
 external INMAX
 data imx /0/
!
 if (imx.eq.0) then
imx = INMAX()
 endif

 if (DBLE(INTGR)+DBLE(INCREMENT).gt.DBLE(imx)) then
INTGR = imx
!IERR = 1   ! error if we get much output. Well, screw that, we want it anyway!
IERR=0
 else
INTGR = INTGR + INCREMENT
IERR = 0
 endif
!
 return
 end
