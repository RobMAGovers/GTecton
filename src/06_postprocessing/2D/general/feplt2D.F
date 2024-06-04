program main

! Program for creating GMT plot files from GTECTON output. The code is
! closely linked with the TECTON code itself. Many of the variable and
! common block names are copied from GTECTON.
!
!
! NOTE: WE HAVE BEGUN REPLACING STATIC BY DYNAMIC MEMORY ALLOCATION. 
! THE OTHER MAJOR CHANGE THAT IS UNDERWAY IS REPLACING COMMON BLOCKS
! BY MODULES. WHILE MAKING THESE CHANGES, THE CLASSICAL FORTRAN77-TYPE
! FORMATTING IS UPDATED.
! CONSEQUENCE IS THAT MUCH OF THE INFORMATION BELOW IS PARTLY INCORRECT:
! VARIABLE AND ARRAY NAMES ARE RETAINED AS MUCH AS POSSIBLE, BUT PARTICULARLY
! THE ARRAY DIMENSIONS AND THEREFORE THE WAY THAT ARRAYS ARE BEING PASSED TO
! SUBROUTINES AND FUNCTIONS HAS, OR IS BEING CHANGED.

!-------------------------------------------------------------------------------
!
!       XPT(NPLDOF,NPLT,NUMTR)          plot coordinates (NUMTR <= NSUR*NUMEL)
!       IJSUR(NSUR,NUMEL)               list of elements surrounding an element
!
!-------------------------------------------------------------------------------
!       Pointer variables in FePlt
!-------------------------------------------------------------------------------
!
!       POINTER         ARRAY           LENGTH
!
!       IPXPT           XPT             NPLDOF*NPLT*NSUR*NUMEL*IWORD
!
!-------------------------------------------------------------------------------
!       Scalars used in FePlt
!-------------------------------------------------------------------------------
!       NSUR = NEN      Number of elements bordering an element
!       NPLT            Number of points in an elementary contouring body
!                       (triangle in 2D)
!       NPLDOF          Number of degrees of freedom associated with an
!                       elementary contouring body. Typically NPLDOF is
!                       determined from NDOF + 1 function value, minus the
!                       effect of projection.
!       SETRIA          Logical denoting whether contour body corners,
!                       that are set up from element centers need to be
!                       calculated.
!       SETSHARE        logical denoting whether a list needs to be made of
!                       element neighbours.
!       READSY          logical: bc-symbols need to be read or not.
!       NUMTR           number of contour triangles generated from element mesh
!       coormn          coordinates minimum
!       coormx          coordinates maximum
!       dispmn(2)       displacement minima
!       dispmx(2)       displacement maxima
!       strsmn(2)       stress minima
!       strsmx(2)       stress maxima
!       ratemn(2)       strain rate minima
!       ratemx(2)       strain rate maxima
!       tempmn          temperatures minimum
!       tempmx          temperatures maximum
!       fluxmn(2)       heat flux minima
!       fluxmx(2)       heat flux maxima
!       deformed        logical: deformed or undeformed geometry.
!       GMT             logical: controls output format (GMT or P)
!       grnum           logical: elements numbered or not
!       npnum           logical: nodal points numbered or not.
!       viscnu          logical: print viscosities as numbers or not.
!       ascin           logical: input file ASCII or binary
!       ascout          logical: "P" output files ASCII or binary
!       sbhor           logical: horizontal scalebar or not.
!       vislog          logical: logarithmic viscosities or not.
!       ratlog          logical: logarithmic strain rates or not.
!       BritDt          logical: time step checking brittle range or not
!       sdev            logical: deviatoric strain or not
!       Sscale          Stress scaling factor
!       Escale          Strain scaling factor
!       Rscale          Strain rate scaling factor
!       Dscale          Displacments scaling factor
!       Uscale          Uplift scaling factor
!       Vscale          Scaling factor for calculating linearized viscosities
!       deffac          Amplification factor of deformation
!       difstn          logical: true differential stress or not
!       Ssb             Stress scale bar size
!       Esb             Strain scale bar size
!       Rsb             Rate scale bar size
!       Dsb             Displacement scale bar size
!       hatangle        Arrow head half top angle
!       hatsize         Arrow head top size
!       sysi            Symbol size
!       sysc            logical: scale bc symbols or not
!       mattsym         use Matt's new symbols designed to operate with GMT 5 macros
!       itime           Integer time
!       ielbas          Element base number
!       ielstp          Element increment number
!       inpbas          Node base number
!       inpstp          Node increment number
!       iarcont         Integer: filled arrows (1) or not (0)
!       iDisp           1st digit 1=displacement
!                       1st digit 2=velocity
!                       2nd digit 0=total
!                       2nd digit 1=average
!                       2nd digit 2=slippery
!                       2nd digit 3=faulted
!                       3rd digit 0=effective
!                       3rd digit 1=x
!                       3rd digit 2=y
!                       3rd digit 3=z
!       npmon           logical: monitor nodal parameters or not
!       elmon           logical: monitor element parameters or not
!       nnel            Number of monitored elements
!       nnp             Number of monitored nodes
!       triang          logical: triangle contouring or not
!       irate           Strain rate type: 0=viscous, 1=total
!       viscli          logical: linearized viscosities or not
!       iDtfl           Time step check mode:
!                       0 = stability ratio
!                       1 = oscillation limit ratio
!                       2 = relative modal accuracy
!       ncol            Column used for force monitor
!       nstn            Stress component that is integrated for force monitor
!       append          logical: append data to file or not
!       Mout            logical: mechanical output at requested time or not
!       Tout            logical: thermal output at requested time or not
!-------------------------------------------------------------------------------

USE MODELDEFINITION, only: IWORD, ntype,feplt, nsur
use modeltopology,   only: NSD
USE CONSTANTS,       only: MB,eps,dfmin, dfmax,dlogmin,dlogmax
use timestepmodule,  only: maxitime
use debugmodule,     only: iecho, &
                           initialiseDebug, &
                           LogTimeStamp, &
                           startClock
use filemodule,      only: do_soelim, &
                           WORKPATHIsSet, &
                           readPathsFromEnvironment
use iomodule,        only: stderr, stdin, bufedit, ilu, lucmd, fcmd
use meshdatamodule,  only: elementNeighborTableComplete, &
                           meshdatactx
use modeldatamodule, only: allocatePlotData
use materialsmodule, only: lmat
use errors,          only: exec, ok, prgerr, stop, syntax, eof

implicit none

#include "lenchar.i" 
! For the lenchar parameter


!-locl
integer, parameter    :: MAXPAR=10000 ! maximum number of command line arguments
character(len=LENCHAR):: cmdln,prgname
character(len=LENCHAR), allocatable :: cmd(:)

integer               :: parseCommandLineArguments,csystem
logical               :: savecmd,rewcf
integer, external     :: lnblk, ifopen, igetcmd, iprinfo ! from tools
integer               :: nrec
integer               :: ierr
double precision      :: TIME
integer               :: lnb
integer               :: ncmd ! number of command line arguments
!-init
data cmdln /' '/
feplt = .false.
IWORD = 2

elementNeighborTableComplete = .false.

allocate(cmd(MAXPAR))

! get temp path
call readPathsFromEnvironment()

call initialiseDebug()

MAXITIME = 0
TIME = 0d0
LMAT = .false.

iecho = 0
do_soelim = .true.
WORKPATHIsSet = .false.

#ifdef SPARSE
call MPI_INIT(ierr)
#endif 

! set geometrical factors
call pltcon()

! establish I/O unit numbers. 
call setio()

! read shell arguments
call shell(prgname,ierr)

! set defaults
call pltdef

if (lnblk(fcmd).le.0) then
    call sayhello()
endif

! 'restart'-jump point
100     call pltdef

call bufopn()   ! open command-line history buffer
rewcf = .true.
nrec = 0
if (lnblk(fcmd).gt.0) then
    ierr = ifopen(lucmd(ilu),fcmd,'old')
    if (ierr.ne.0) then
        ierr = EXEC
        lucmd(ilu) = stdin
    endif
else
    lucmd(ilu) = stdin
endif

! read a command line
200     ierr = igetcmd(prgname,cmdln,rewcf)
nrec = nrec + 1
rewcf = .false.

if (ierr.eq.SYNTAX) then
    write(stderr,10) nrec
10       format(1x,'Syntax error or unknown command in line ',I3/1x, &
     'type "?" for a list of available commands')
    goto 200
endif


if (ierr.eq.EOF) then
    call xeceof (ierr)
    if (ierr.eq.STOP .or. ierr.eq.PRGERR) then
        write(stderr,*) "Encountered unexpected End Of File"
        goto 1000
    endif
    goto 200
endif

savecmd = (ilu.eq.0)

!       parse command line
ierr = parseCommandLineArguments(cmdln,cmd,ncmd,MAXPAR)

if (ncmd.lt.1) then
    if (lucmd(ilu).eq.stdin) then
        savecmd = .false.
    endif
    ! we have had an open line. Skip this
!    write(stderr,*) "Encountered less than one command argument"
!   not necessary to report every open line in a plnplt script as an error,
!   unless somebody disagrees...
    goto 1000
endif

if (ierr.ne.OK) then
    write(stderr,10) nrec
    goto 200
endif

lnb=lnblk(cmdln)

if (iecho.eq.1) then
    if (lnb.gt.0) then
        write(stderr,20) cmdln(1:lnb)
    else
        write(stderr,20)
    endif
endif

20   format(1x,1024a)
!
!       classify and take action
if ( cmdln(1:1).eq.'?' ) then
    if (lnb.gt.1) then
        cmd(1) = cmdln(2:lnb)
        cmdln = cmd(1)
        ierr = parseCommandLineArguments(cmdln,cmd(2),ncmd,MAXPAR)
        if (ierr.ne.OK) goto 200
        ncmd = ncmd + 1
    endif
    ierr = iprinfo(cmd,ncmd)
    savecmd = .false.
else if ( cmd(1)(1:1).eq.'!' ) then
    if (lnb.ge.2) then
        ierr = csystem(cmdln(2:lnb))
        if (ierr.ne.OK) ierr = EXEC
    endif
else if ( cmd(1)(1:4).eq.'femi' ) then
    if (iecho.eq.6) call startClock()
    call xfemi (cmd,ncmd,ierr)
	! now that femi has been executed, we now the number of elements and such,
	! and we can allocate the work array for the plotting routines, XPT
    call allocatePlotdata(meshdatactx%nvglobal, &
                          meshdatactx%neglobal, & 
                          nsur, nsd)
    if (iecho.eq.6) call LogTimeStamp("xfemi")
else if ( cmd(1)(1:5).eq.'sides' ) then
    if (iecho.eq.6) call startClock()
    call xsides (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xsides")
elseif ( cmd(1)(1:2).eq.'wf'.and.NTYPE.eq.2 ) then
    if (iecho.eq.6) call startClock()
    call xwfemi (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xwfemi")
elseif ( cmd(1)(1:4).eq.'femo' ) then
    if (iecho.eq.6) call startClock()
    call xfemo (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xfemo")
elseif ( cmd(1)(1:4).eq.'core' ) then
    if (iecho.eq.6) call startClock()
    call xcore (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xcore")
elseif ( cmd(1)(1:4).eq.'cons' ) then
    if (iecho.eq.6) call startClock()
    call xcons (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xcons")
elseif ( cmd(1)(1:3).eq.'np' ) then
    if (iecho.eq.6) call startClock()
    call xnp (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xnp")
elseif ( cmd(1)(1:3).eq.'gr '.or. cmd(1)(1:3).eq.'gri' ) then
    if (iecho.eq.6) call startClock()
    call xgrid (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xgrid")
elseif ( cmd(1)(1:4).eq.'elmc' ) then
    if (iecho.eq.6) call startClock()
    call xelmc (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xelmc")
elseif ( cmd(1)(1:4).eq.'stre' ) then
    if (iecho.eq.6) call startClock()
    call xstress (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xstress")
elseif ( cmd(1)(1:4).eq.'stra' ) then
    if (iecho.eq.6) call startClock()
    call xstrain (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xstrain")
elseif ( cmd(1)(1:4).eq.'save' ) then
    if (iecho.eq.6) call startClock()
    call xsave (ierr)
    if (iecho.eq.6) call LogTimeStamp("xsave")
elseif ( cmd(1)(1:4).eq.'sedi' ) then
    if (iecho.eq.6) call startClock()
    call xsedi (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xsedi")
elseif ( cmd(1)(1:2).eq.'ra' ) then
    if (iecho.eq.6) call startClock()
    call xrate (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xrate")
elseif ( cmd(1)(1:2).eq.'te' ) then
    if (iecho.eq.6) call startClock()
    call xte (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xte")
elseif ( cmd(1)(1:2).eq.'di' ) then
    if (iecho.eq.6) call startClock()
    call xdispl (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xdispl")
elseif ( cmd(1)(1:2).eq.'ve' ) then
    if (iecho.eq.6) call startClock()
    call xvelo (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xvelo")
elseif ( cmd(1)(1:2).eq.'fl' ) then
    if (iecho.eq.6) call startClock()
    call xflux (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xflux")
elseif ( cmd(1)(1:9).eq.'invariant' ) then
    if (iecho.eq.6) call startClock()
    call xInvariant (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xInvariant")
elseif ( cmd(1)(1:2).eq.'pr' ) then
    if (iecho.eq.6) call startClock()
    call xpr (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xpr")
elseif ( cmd(1)(1:3).eq.'vol' ) then
    if (iecho.eq.6) call startClock()
    call xVolumetricStrain (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xVolumetricStrain")
elseif ( cmd(1)(1:5).eq.'error' ) then
    if (iecho.eq.6) call startClock()
    call xError (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xError")
elseif ( cmd(1)(1:5).eq.'estre' ) then
    if (iecho.eq.6) call startClock()
    call xestress (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xestress")
elseif ( cmd(1)(1:5).eq.'estra' ) then
    if (iecho.eq.6) call startClock()
    call xestrain (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xestrain")
elseif ( cmd(1)(1:2).eq.'es' .or. cmd(1)(1:11).eq.'elementsize') then
    if (iecho.eq.6) call startClock()
    call xElementSize (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xes")

elseif ( cmd(1)(1:6) .eq.'eratio' .or. &
         cmd(1)(1:12).eq.'elementratio') then
    if (iecho.eq.6) call startClock()
    call xElementRatio (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xer")

elseif ( cmd(1)(1:5) .eq.'slips') then
    if (iecho.eq.6) call startClock()
    call xSlipSides (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xSlipSides")

elseif ( cmd(1)(1:3).eq.'era' ) then
    if (iecho.eq.6) call startClock()
    call xeffr (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xeffr")

elseif ( cmd(1)(1:2) .eq.'er') then
    if (iecho.eq.6) call startClock()
    call xElementRatio (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xer")

elseif ( cmd(1)(1:3).eq.'efl' ) then
    if (iecho.eq.6) call startClock()
    call xeflux (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xeflux")
elseif ( cmd(1)(1:3).eq.'edi' ) then
    if (iecho.eq.6) call startClock()
    call xedisp (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xedisp")
elseif ( cmd(1)(1:3).eq.'eve' ) then
    if (iecho.eq.6) call startClock()
    call xevel (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xevel")
elseif ( cmd(1)(1:4).eq.'gref' ) then
    if (iecho.eq.6) call startClock()
    call xgref (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xgref")
elseif ( cmd(1)(1:3).eq.'gra' ) then
    if (iecho.eq.6) call startClock()
    call xgrav (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xgrav")
elseif ( cmd(1)(1:4).eq.'heat' ) then
    if (iecho.eq.6) call startClock()
    call xheat (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xheat")
elseif ( cmd(1)(1:4).eq.'surf' ) then
    if (iecho.eq.6) call startClock()
    write(*,*) "Command ''surf'' is no longer supported."
!    call xsrflx (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xsrflx")
elseif ( cmd(1).eq.'read' ) then
    if (iecho.eq.6) call startClock()
    call xread (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xread")
elseif ( cmd(1)(1:4).eq.'samp' ) then
    if (iecho.eq.6) call startClock()
    call xsample (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xsampl")
elseif ( cmd(1)(1:2).eq.'yi' ) then
    if (iecho.eq.6) call startClock()
    call xyield (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xyield")
elseif ( cmd(1)(1:2).eq.'ma' ) then
    if (iecho.eq.6) call startClock()
    call xma (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xma")
elseif ( cmd(1)(1:4).eq.'mesh' ) then
    if (iecho.eq.6) call startClock()
    call xmesh (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xma")
elseif ( cmd(1)     .eq.'cp' ) then
    if (iecho.eq.6) call startClock()
    call xcp (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xcp")
elseif ( cmd(1)(1:2).eq.'du' ) then
    if (iecho.eq.6) call startClock()
    write(*,*) "Command ''du'' is no longer supported."
!    call xdu (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xdu")
elseif ( cmd(1)(1:2).eq.'vi' ) then
    if (iecho.eq.6) call startClock()
    call xvi (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xvi")
elseif ( cmd(1)(1:2).eq.'fo' ) then
    write(*,*) "Sorry, Sahib, command ''fo'' no longer used."
!    if (iecho.eq.6) call startClock()
!     call xecfo (cmd,ncmd,ierr)
!    if (iecho.eq.6) call LogTimeStamp("xecfo")
elseif ( cmd(1)(1:2).eq.'be' ) then
    if (iecho.eq.6) call startClock()
    write(*,*) "Command ''be'' is no longer supported."
!    call xbeta (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xbeta")
elseif ( cmd(1)     .eq.'dt' ) then
    if (iecho.eq.6) call startClock()
    call xdt (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xdt")
elseif ( cmd(1)     .eq.'nbc' ) then
    if (iecho.eq.6) call startClock()
    call xnbc (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xnbc")
elseif ( cmd(1)     .eq.'ebc' ) then
    if (iecho.eq.6) call startClock()
    call xebc (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xebc")
elseif ( cmd(1)     .eq.'tc' ) then
    if (iecho.eq.6) call startClock()
    call xtc (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xtc")
elseif ( cmd(1)(1:2).eq.'sc' ) then
    if (iecho.eq.6) call startClock()
    write(*,*) "Command ''sc'' is no longer supported."
!    call xsc (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xsc")
elseif ( cmd(1)(1:2).eq.'ec' ) then
    if (iecho.eq.6) call startClock()
    call echo(cmdln,cmd,ncmd)
    if (iecho.eq.6) call LogTimeStamp("echo")
elseif ( cmd(1)(1:4).eq.'shsl' ) then
    if (iecho.eq.6) call startClock()
    call xshsl (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xshsl")
elseif ( cmd(1)(1:5).eq.'bound' ) then
    if (iecho.eq.6) call startClock()
    call xBoundaries (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xBoundaries")
elseif ( cmd(1)(1:2).eq.'sh' ) then  ! show
    if (iecho.eq.6) call startClock()
    if (lucmd(ilu).eq.stdin) then
        savecmd = .false.
    endif
    call shcmd (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("shcmd")
elseif ( cmd(1)(1:2).eq.'q ' .or. &
         cmd(1)(1:2).eq.'qu' .or. &
         cmd(1).eq.'q!'.or. cmd(1).eq.'quit!' .or. &
         cmd(1)(1:4).eq.'exit') then
    if (lucmd(ilu).eq.stdin) savecmd = .false.
    call xecq (cmd,ncmd,ierr)
elseif ( cmd(1)(1:1).eq.'#' ) then
    continue
elseif ( cmd(1)     .eq.'wait' ) then
    call wait (ierr)
elseif ( cmd(1)(1:2).eq.'se' ) then
    call setcmd (cmd,ncmd,ierr)
elseif ( cmd(1)(1:2).eq.'mo' ) then
    if (iecho.eq.6) call startClock()
    call xmo (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xmo")
elseif ( cmd(1)(1:4).eq.'nomo' ) then
    if (iecho.eq.6) call startClock()
    call xnomo (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xnomo")
elseif ( cmd(1)(1:2).eq.'u ' ) then
    if (iecho.eq.6) call startClock()
    call xecu (cmd,ncmd,ierr)
    savecmd = .false.
    if (iecho.eq.6) call LogTimeStamp("xecu")
elseif ( cmd(1)(1:2).eq.'h '.or. cmd(1)(1:4).eq.'hist' ) then
    if (iecho.eq.6) call startClock()
    call xech (cmd,ncmd,ierr)
    savecmd = .false.
    if (iecho.eq.6) call LogTimeStamp("xech")
elseif ( cmd(1)(1:3).eq.'cf ' ) then
    if (iecho.eq.6) call startClock()
    call xcf (cmd,ncmd,ierr)
    if (ierr.eq.OK) rewcf = .true.
    if (iecho.eq.6) call LogTimeStamp("xcf")
elseif ( cmd(1)(1:2).eq.'e ' .or. cmd(1)(1:4).eq.'edit' .or. &
 cmd(1).eq.'e!') then
    savecmd = .false.
    if (iecho.eq.6) call startClock()
    call xedit (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xedit")
elseif ( cmd(1)(1:2).eq.'w ' .or. cmd(1).eq.'w!') then
    savecmd = .false.
    if (iecho.eq.6) call startClock()
    call xecw (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xecw")
elseif ( cmd(1)(1:3).eq.'we ') then
    savecmd = .false.
    if (iecho.eq.6) call startClock()
    call xecw (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xecw")
    if (ierr.eq.OK) call xedit (cmd,ncmd,ierr)
elseif ( cmd(1)(1:3).eq.'wq ') then
    savecmd = .false.
    if (iecho.eq.6) call startClock()
    call xecw (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xecw")
    if (ierr.eq.OK) call xecq (cmd,ncmd,ierr)
elseif ( cmd(1)(1:2).eq.'r ' .or. cmd(1)(1:3).eq.'res' .or. &
 cmd(1).eq.'r!') then
    savecmd = .false.
    if (iecho.eq.6) call startClock()
    call xecr (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xecr")
!    if (ierr.ne.OK) then
!        write(stderr,*) "Something with the timestamping went awry."
!        goto 1000
!    endif
    rewcf = .true.
    goto 100
elseif (cmd(1)(1:4).eq.'rot '.or.cmd(1)(1:4).eq.'rota') then
    if (iecho.eq.6) call startClock()
    call xrot (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xrot")
elseif (cmd(1)(1:4).eq.'rotr') then
    if (iecho.eq.6) call startClock()
    call xrotr (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xrotr")
elseif (cmd(1)(1:5).eq.'good') then
    if (iecho.eq.6) call startClock()
    call xgood (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xgood")
elseif (cmd(1)(1:6).eq.'regrid') then
    if (iecho.eq.6) call startClock()
    write(*,*) "Regridding temporarily disabled"
!    call xregrid (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xregrid")
elseif (cmd(1)(1:4).eq.'conn') then
    if (iecho.eq.6) call startClock()
    call xconn (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xconn")
elseif (cmd(1)(1:5).eq.'slip ' .or. cmd(1)(1:5).eq.'slipp') then
    if (iecho.eq.6) call startClock()
    call xslip (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xslip")
elseif (cmd(1)(1:4).eq.'dens') then
    if (iecho.eq.6) call startClock()
    call xdens (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xdens")
elseif (cmd(1)(1:4).eq.'part') then
    if (iecho.eq.6) call startClock()
    call xpartition (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xpart")
#ifdef EXPERIMENTAL_FRICTION
elseif (cmd(1)(1:4).eq.'fric') then
    if (iecho.eq.6) call startClock()
    call xfriction (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xpart")
#endif
else
    write(stderr,*) "Command ", trim(cmd(1)), " not recognized"
    write(stderr,10) nrec
    ierr = SYNTAX
endif

1000  if (lucmd(ilu).eq.stdin) then
    if (ierr.eq.OK .and. savecmd) then
        call bufput (cmdln,ierr)
        bufedit = .true.
    endif
else
    if ((ierr.eq.OK .or. ierr.eq.EXEC) .and. savecmd) then
        call bufput (cmdln,ierr)
    endif

    if (ierr.eq.SYNTAX .or. ierr.eq.EXEC) then

        if (lucmd(ilu).ne.stdin) then
              write(stderr,30) nrec
30           format(1x,'input error in line ',I3)
             write(*,*) "possible cause: trying to write output to a directory that does not exist."
        endif

        ierr = OK
        call wait(ierr)

        if (ierr.eq.EOF) then
            ilu = 0
            lucmd(ilu)  = stdin
        endif
    endif

endif

if (ierr.eq.STOP) then
    call bufcls()
    call exitp(0)
endif

if (ierr.eq.PRGERR) then
    call exitp(1)
endif

ierr = OK
goto 200
end
!-----------------------------------------------------------------------
subroutine sayhello()

use iomodule,        only: stderr

write(stderr,*) "This software is protected under the Creative Commons Attribution-NonCommercial-ShareAlike license,"
write(stderr,*)	"version 4.0, international: https://creativecommons.org/licenses/by-nc-sa/4.0/"
write(stderr,*) " "
write(stderr,*) "For more information, contact:"
write(stderr,*) "Rob Govers,"
write(stderr,*) "Department of Earth Sciences,"
write(stderr,*) "Utrecht University, Netherlands."
write(stderr,*) "E-mail: r.govers@uu.nl"
return
end
