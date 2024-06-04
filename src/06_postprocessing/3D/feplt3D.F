program main

! **************************************************************************************************
! *                                                                                                *   
! * Plt3d is a research tool aimed at providing ways to                                            *
! *     1) output files for direct visualization of FEM in- and output data                        *     
! *     2) monitor/probe the FEM output for nodes or elements with specific tags                   *                
! *     3) monitor/probe FEM output for specific node or element numbers                           *        
! *                                                                                                *
! * Most of the output is written as Visualization Toolkit (VTK) files. In the current version, we *
! * use the VTK variant for SERIAL unstructured data (VTU). The data format is limited to ASCII,   *
! * i.e., there is no BINARY support.                                                              *
! *                                                                                                *
! * The functionality of the current version of the code is very limited. The most important       *
! * limitation is that the change in output data structure (npbased to elmbased) that comes with   *
! * the introduction of faulted or slippery nodes is hardly recognized and implemented.            *
! * Monitoring (aims 2 and 3 above) is highly incompletely supported.                              *
! *                                                                                                *   
! **************************************************************************************************
!
! The code is closely linked with GTECTON and GMSH. Many of the variable and
! common block names are copied from GTECTON.
!
! NOTE: WE HAVE BEGUN REPLACING STATIC BY DYNAMIC MEMORY ALLOCATION. 
! THE OTHER MAJOR CHANGE THAT IS UNDERWAY IS REPLACING COMMON BLOCKS
! BY MODULES. WHILE MAKING THESE CHANGES, THE CLASSICAL FORTRAN77-TYPE
! FORMATTING IS UPDATED.
! CONSEQUENCE IS THAT MUCH OF THE INFORMATION BELOW IS PARTLY INCORRECT:
! VARIABLE AND ARRAY NAMES ARE RETAINED AS MUCH AS POSSIBLE, BUT PARTICULARLY
! THE ARRAY DIMENSIONS AND THEREFORE THE WAY THAT ARRAYS ARE BEING PASSED TO
! SUBROUTINES AND FUNCTIONS HAS, OR IS BEING CHANGED.
!
!-------------------------------------------------------------------------------
!    Arrays used in FePlt
!-------------------------------------------------------------------------------
!
!    XPT(NPLDOF,NPLT,NUMTR)        plot coordinates (NUMTR <= NEN*NUMEL)
!    IJSUR(NEN,NUMEL)        list of elements surrounding an element
!
!-------------------------------------------------------------------------------
!    Pointer variables in FePlt
!-------------------------------------------------------------------------------
!
!    POINTER        ARRAY        LENGTH
!
!    IPXPT        XPT        NPLDOF*NPLT*NEN*NUMEL*IWORD
!
!-------------------------------------------------------------------------------
!    Scalars used in FePlt
!-------------------------------------------------------------------------------
!    NPLT        Number of points in an elementary contouring body
!            (triangle in 2D)
!    NPLDOF        Number of degrees of freedom associated with an
!            elementary contouring body. Typically NPLDOF is
!            determined from NDOF + 1 function value, minus the
!            effect of projection.
!    SETRIA        Logical denoting whether contour body corners,
!            that are set up from element centers need to be
!            calculated.
!    SETSHARE    logical denoting whether a list needs to be made of
!            element neighbours.
!    READSY        logical: bc-symbols need to be read or not.
!    NUMTR        number of contour triangles generated from element mesh
!    coormn        coordinates minimum
!    coormx        coordinates maximum
!    dispmn(2)    displacement minima
!    dispmx(2)    displacement maxima
!    strsmn(2)    stress minima
!    strsmx(2)    stress maxima
!    ratemn(2)    strain rate minima
!    ratemx(2)    strain rate maxima
!    tempmn       temperatures minimum
!    tempmx       temperatures maximum
!    fluxmn(2)    heat flux minima
!    fluxmx(2)    heat flux maxima
!    deformed    logical: deformed or undeformed geometry.
!    grnum        logical: elements numbered or not
!    npnum        logical: nodal points numbered or not.
!    viscnu        logical: print viscosities as numbers or not.
!    ascin        logical: input file ASCII or binary
!    ascout        logical: "P" output files ASCII or binary
!    sbhor        logical: horizontal scalebar or not.
!    vislog        logical: logarithmic viscosities or not.
!    ratlog        logical: logarithmic strain rates or not.
!    verbose        logical: extensive output or not.
!    BritDt        logical: time step checking brittle range or not
!    stndev        logical: deviatoric strain or not
!    Sscale        Stress scaling factor
!    Escale        Strain scaling factor
!    Rscale        Strain rate scaling factor
!    Dscale        Displacments scaling factor
!    Uscale        Uplift scaling factor
!    Vscale        Scaling factor for calculating linearized viscosities
!    deffac        Amplification factor of deformation
!    difstn        logical: true differential stress or not
!    Ssb        Stress scale bar size
!    Esb        Strain scale bar size
!    Rsb        Rate scale bar size
!    Dsb        Displacement scale bar size
!    hatangle    Arrow head half top angle
!    hatsize        Arrow head top size
!    sysi        Symbol size
!    sysc        logical: scale bc symbols or not
!    itime        Integer time
!    ielbas        Element base number
!    ielstp        Element increment number
!    inpbas        Node base number
!    inpstp        Node increment number
!    iarcont        Integer: filled arrows (1) or not (0)
!    iEffD        0 = total effective displacements
!            1 = normal effective displacements
!            2 = slippery node effective displacements
!            3 = faulted node effective displacements
!    iDisp        1 = normal displacement arrows
!            2 = slippery node displacement arrows
!            3 = faulted node displacement arrows
!    npmon        logical: monitor nodal parameters or not
!    elmon        logical: monitor element parameters or not
!    nnel        Number of monitored elements
!    nnp        Number of monitored nodes
!    triang        logical: triangle contouring or not
!    irate        Strain rate type: 0=viscous, 1=total
!    viscli        logical: linearized viscosities or not
!    iDtfl        Time step check mode:
!            0 = stability ratio
!            1 = oscillation limit ratio
!            2 = relative modal accuracy
!    ncol        Column used for force monitor
!    nstn        Stress component that is integrated for force monitor
!    append        logical: append data to file or not
!    Mout        logical: mechanical output at requested time or not
!    Tout        logical: thermal output at requested time or not
!-------------------------------------------------------------------------------
!

use filemodule,      only: do_soelim
USE MODELDEFINITION
USE TIMESTEPMODULE
use meshdatamodule,  only: elementNeighborTableComplete, &
                           BuildElementNeighbourTable3D, &
                           meshdatactx
use modeltopology,   only: NSD
USE MODELCTX
use constants,       only: eps,dfmin,dlogmin,dfmax,dlogmax
use debugmodule,     only: lustat, debug, itest, NewSample, iecho, &
                           startClock, &
                           logTimeStamp
use errors,          only: exec, OK, prgerr, stop, syntax, eof
use VTKmodule,       only: VTKInititialisePlotoptions, nOutputTags, connDone, gmshRead, &
                           outputTags, meshAdjusted
use iomodule,        only: icoord, bufedit, fcmd, ilu
use modeltopology,   only: NEN
use iomodule,        only: stderr, stdin, crdfil,elmfil,ctrfil, lucmd
use modeldatamodule, only: allocatePlotData
use materialsmodule, only: lmat

implicit none

#include "lenchar.i" /* For the lenchar parameter */




!-locl
integer, parameter :: MAXPAR=10000
!#include "lenchar.i"
character(len=LENCHAR) cmdln,prgname

character(len=LENCHAR), allocatable :: cmd(:)

!dimension cmd(MAXPAR)
integer parseCommandLineArguments,csystem
logical savecmd,rewcf
integer, external :: lnblk, igetcmd, ifopen, iprinfo ! from tools

real :: start, finish

integer :: ierr
integer :: lnb
integer :: ncmd

!-init
data cmdln /' '/
feplt       = .false.
IWORD       = 2  ! 2 4-byte blocks per getal. Double precision, dus 
LMAT        = .false.
iecho       = 0
do_soelim   = .true.
connDone    = .false.
gmshRead    = .false.
outputTags  = 0
nOutputTags = 0

meshAdjusted = .false.  ! from vtkmodule

allocate(cmd(MAXPAR))

NEN = 4

! initialize sample flag
! To make sure that lookup table 
! for ssaple command is ony created once.
elementNeighborTableComplete = .false.


! initialize the VTK plot parameters.
! By default, they (except the displacement)
! are all set to false.
! by setting more to true, more is plotted:
call VTKInititialisePlotoptions()

!    set geometrical factors
call pltcon()

!       establish I/O unit numbers
!        note thet plt3d has its own setio. It does not have the same 
!        setio has used by gtecton. LWItodo: this might be fixed.
call setio()

!       read shell arguments
call shell(prgname,ierr)

!    set defaults
call pltdef

!    say hello
if (lnblk(fcmd).le.0) then
    call sayhello()
endif

!    `restart'-jump point
100    call pltdef

call bufopn()

rewcf = .true.
if (lnblk(fcmd).gt.0) then
    ierr = ifopen(lucmd(ilu),fcmd,'old')
    if (ierr.ne.0) then
       ierr = EXEC
       lucmd(ilu) = stdin
    endif
else
    lucmd(ilu) = stdin
endif


!    read a command line
200     ierr = igetcmd(prgname,cmdln,rewcf)


rewcf = .false.
if (ierr.eq.SYNTAX) then
    write(stderr,10)
10       format(1x,'Syntax error or unknown command'/1x, &
     'type "?" for a list of available commands')
    goto 200
endif

if (ierr.eq.EOF) then
    call xeceof (cmd,ncmd,ierr)
    if (ierr.eq.STOP .or. ierr.eq.PRGERR) then
        goto 1000
    endif
    goto 200
endif
savecmd = (ilu.eq.0)

if (iecho.eq.8) then
    write(*,*) "Parsing command line arguments" 
endif


!    parse command line
ierr = parseCommandLineArguments(cmdln,cmd,ncmd,MAXPAR)

if (ncmd.lt.1) then
    if (lucmd(ilu).eq.stdin) then
        savecmd = .false.
    endif
    goto 1000
endif

if (ierr.ne.OK) then
    write(stderr,10)
    goto 200
endif

! write content of command line to standard error, on echo
lnb=lnblk(cmdln)
if (iecho.eq.1) then
    if (lnb.gt.0) then
        write(stderr,20) cmdln(1:lnb)
    else
        write(stderr,20)
    endif
endif
20    format(1x,1024a)



!    classify and take action
if ( cmd(1)(1:1).eq.'?' ) then
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
    ierr = csystem(cmdln(2:lnblk(cmdln)))
    if (ierr.ne.OK) ierr = EXEC
endif


else if ( cmd(1)(1:8).eq.'readmesh' ) then
    if (iecho.eq.6) call startClock()
    call xReadMesh (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("readmesh")
else if ( cmd(1)(1:4).eq.'femi' ) then
    if (iecho.eq.6) call startClock()
    call xfemi (cmd,ncmd,ierr)
    call allocatePlotdata(meshdatactx%nvglobal, &
                          meshdatactx%neglobal, &
                          nsur, nsd)
    if (iecho.eq.6) call LogTimeStamp("femi")
elseif ( cmd(1)(1:4).eq.'femo' ) then
    if (iecho.eq.6) call startClock()
    call xfemo (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("femo")
elseif ((cmd(1)(1:3).eq.'vtu').or. &
        (cmd(1)(1:5).eq.'write').or. &
        (cmd(1)(1:8).eq.'vtuwrite')) then

    if (iecho.eq.6) call startClock()
    call xvtu (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("writing VTU")

elseif ((cmd(1)(1:4).eq.'time').or. &
        (cmd(1)(1:7).eq.'vtutime')) then
    if (iecho.eq.6) call startClock()
    call xvtutime (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("time")
elseif ( cmd(1)(1:4).eq.'conn' ) then
    if (iecho.eq.6) call startClock()
    call xconn (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("conn")
elseif ( cmd(1)(1:4).eq.'samp' ) then
    if (.not. elementNeighborTableComplete .and. Newsample) then
        if (iecho.eq.6) call startClock() 
        call BuildElementNeighbourTable3D() ! in meshdatamodule
        if (iecho.eq.6) call LogTimeStamp("Sample lookoup table")
    endif
    if (iecho.eq.6) call startClock()
    call xsample (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("sampling")

!        elseif ( cmd(1)(1:4).eq.'cons' ) then
!            call xcons (cmd,ncmd,ierr)
elseif ( cmd(1)(1:3).eq.'np ' ) then
    if (iecho.eq.6) call startClock() 
    call xnp (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("np")
elseif ( cmd(1)(1:3).eq.'gr '.or. &
         cmd(1)(1:3).eq.'gri' ) then
    if (iecho.eq.6) call startClock() 
    call xgrid (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("grid")
elseif ( cmd(1)(1:8).eq.'traction' ) then
    if (iecho.eq.6) call startClock()
    call xtraction (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xtraction")
elseif ( cmd(1)(1:4).eq.'elmc' ) then
    if (iecho.eq.6) call startClock() 
    call xelmc (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("elmcenter")
elseif ( cmd(1)(1:4).eq.'stre' ) then
    if (iecho.eq.6) call startClock() 
    call xstress (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("stress")
#ifdef EXPERIMENTALTRANSIENTDISCREEP
elseif ( cmd(1)(1:12).eq.'taylorstress' ) then
    if (iecho.eq.6) call startClock() 
    call xtaylorstress (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("taylorstress")
!elseif ( cmd(1)(1:12).eq.'thresholdstress' ) then
!    if (iecho.eq.6) call startClock() 
!    call xthresholdstress (cmd,ncmd,ierr)
!    if (iecho.eq.6) call LogTimeStamp("thresholdstress")    
#endif    
elseif ( cmd(1)(1:4).eq.'stra' ) then
    if (iecho.eq.6) call startClock() 
    call xstrain (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("strain")
!        elseif ( cmd(1)(1:4).eq.'save' ) then
!            call xsave (cmd,ncmd,ierr)
!        elseif ( cmd(1)(1:4).eq.'sedi' ) then
!            call xsedi (cmd,ncmd,ierr)
elseif ( cmd(1)(1:2).eq.'ra' ) then
    if (iecho.eq.6) call startClock() 
    call xrate (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("strain rate")
elseif ( cmd(1)(1:2).eq.'te' ) then
    if (iecho.eq.6) call startClock() 
    call xte (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("te")
elseif ( cmd(1)(1:2).eq.'di' ) then
    if (iecho.eq.6) call startClock() 
    call xdispl (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("displacement")
elseif ( cmd(1)(1:2).eq.'fl' ) then
    if (iecho.eq.6) call startClock() 
    call xflux (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("flux")
!        elseif ( cmd(1).eq.'dr' ) then
!            call xdr (cmd,ncmd,ierr)
!        elseif ( cmd(1)(1:2).eq.'pr' ) then
!            call xpr (cmd,ncmd,ierr)
elseif ( cmd(1)(1:5).eq.'estre' ) then
    if (iecho.eq.6) call startClock() 
    call xestress (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("effective stress")
elseif ( cmd(1)(1:5).eq.'estra' ) then
    if (iecho.eq.6) call startClock() 
    call xestrain (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("effective strain")
elseif ( cmd(1)(1:3).eq.'era' ) then
    if (iecho.eq.6) call startClock() 
    call xeffr (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("effective strain rate")
!        elseif ( cmd(1)(1:3).eq.'efl' ) then
!            call xeflux (cmd,ncmd,ierr)
elseif ( cmd(1)(1:3).eq.'edi' ) then
    if (iecho.eq.6) call startClock() 
    call xedisp (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("effective displacement")
!        elseif ( cmd(1)(1:4).eq.'gref' ) then
!            call xgref (cmd,ncmd,ierr)
!        elseif ( cmd(1)(1:3).eq.'gra' ) then
!            call xgrav (cmd,ncmd,ierr)
!        elseif ( cmd(1)(1:4).eq.'heat' ) then
!            call xheat (cmd,ncmd,ierr)
!        elseif ( cmd(1)(1:4).eq.'surf' ) then
!            call xsrflx (cmd,ncmd,ierr)
!        elseif ( cmd(1).eq.'read' ) then
!            call xread (cmd,ncmd,ierr)
!        elseif ( cmd(1)(1:4).eq.'resa' ) then
!            call xrsmpl (cmd,ncmd,ierr)
!        elseif ( cmd(1)(1:2).eq.'yi' ) then
!            call xyield (cmd,ncmd,ierr)
elseif ( cmd(1)(1:7).eq.'fullgra' ) then
    ! compute gravity of comlete domain
    call xfullgrav (cmd,ncmd,ierr)
elseif ( cmd(1)(1:2).eq.'ma' ) then
    if (iecho.eq.6) call startClock() 
    call xmat (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("mat")
elseif ( cmd(1)     .eq.'cp' ) then
    if (iecho.eq.6) call startClock() 
    call xcp (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("copy")
elseif ( cmd(1)(1:2).eq.'vi' ) then
    if (iecho.eq.6) call startClock() 
    call xvi (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("vi")
elseif ( cmd(1)(1:2).eq.'ve' ) then
    if (iecho.eq.6) call startClock() 
    call xvelo (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("velocity")
!        elseif ( cmd(1)(1:2).eq.'fo' ) then
!            call xecfo (cmd,ncmd,ierr)
!        elseif ( cmd(1)(1:2).eq.'be' ) then
!            call xbeta (cmd,ncmd,ierr)
!        elseif ( cmd(1)     .eq.'dt' ) then
!            call xdt (cmd,ncmd,ierr)
elseif ( cmd(1)     .eq.'nbc' ) then
    if (iecho.eq.6) call startClock() 
    call xnbc (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("nbc")
elseif ( cmd(1)(1:3).eq.'deb' ) then
    if (iecho.eq.6) call startClock() 
    call xdebrah (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("debrah")
elseif ( cmd(1)(1:3).eq.'cou' ) then
    if (iecho.eq.6) call startClock() 
    call xcoulomb (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("coulomb")
elseif ( cmd(1)     .eq.'ebc' ) then
    if (iecho.eq.6) call startClock() 
    call xebc (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("ebc")
elseif ( cmd(1)     .eq.'tc' ) then
    if (iecho.eq.6) call startClock() 
    call xtc (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("tc")
elseif ( cmd(1)     .eq.'rt' ) then
    if (iecho.eq.6) call startClock() 
    call xreslv (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("reslv")
elseif ( cmd(1)(1:3).eq.'rot' ) then
    if (iecho.eq.6) call startClock() 
    call xrot (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("rot")
!        elseif ( cmd(1)(1:2).eq.'sc' ) then
!            call xsc (cmd,ncmd,ierr)
elseif ( cmd(1)(1:2).eq.'ec' ) then
    if (iecho.eq.6) call startClock() 
    call echo(cmdln,cmd,ncmd)
    if (iecho.eq.6) call LogTimeStamp("echo")
elseif ( cmd(1)(1:2).eq.'sh' ) then
    if (lucmd(ilu).eq.stdin) then
        savecmd = .false.
    endif
        call shcmd (cmd,ncmd,ierr)
elseif ( cmd(1)(1:2).eq.'q '.or. &
         cmd(1)(1:2).eq.'qu' .or. &
         cmd(1).eq.'q!'.or. &
         cmd(1).eq.'quit!'  .or. &
         cmd(1)(1:4).eq.'exit') then
        if (lucmd(ilu).eq.stdin) then
             savecmd = .false.
        endif
        call xecq (cmd,ncmd,ierr)
elseif ( cmd(1)(1:1).eq.'#' ) then
    continue
elseif ( cmd(1)     .eq.'wait' ) then
    call wait (ierr)
elseif ( cmd(1)(1:2).eq.'se' ) then
    if (iecho.eq.6) call startClock() 
    call setcmd (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("set")
elseif ( cmd(1)(1:2).eq.'mo' ) then
    if (iecho.eq.6) call startClock() 
    call xmo (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("monitor")
elseif ( cmd(1)(1:4).eq.'nomo' ) then
    if (iecho.eq.6) call startClock() 
    call xnomo (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("no monitor")
elseif ( cmd(1)(1:2).eq.'u ' ) then
    if (iecho.eq.6) call startClock() 
    call xecu (cmd,ncmd,ierr)
    savecmd = .false.
    if (iecho.eq.6) call LogTimeStamp("ecu")
elseif ( cmd(1)(1:2).eq.'h '.or. &
         cmd(1)(1:4).eq.'hist' ) then
    if (iecho.eq.6) call startClock() 
    call xech (cmd,ncmd,ierr)
    savecmd = .false.
    if (iecho.eq.6) call LogTimeStamp("hist")
elseif ( cmd(1)(1:3).eq.'cf ' ) then
    if (iecho.eq.6) call startClock() 
    call xcf (cmd,ncmd,ierr)
    if (ierr.eq.OK) then
        rewcf = .true.
    endif
    if (iecho.eq.6) call LogTimeStamp("cf")
elseif ( cmd(1)(1:2).eq.'e ' .or. &
         cmd(1)(1:4).eq.'edit' .or. &
         cmd(1).eq.'e!') then
    if (iecho.eq.6) call startClock()
    savecmd = .false.
    call xedit (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("edit")
elseif ( cmd(1)(1:2).eq.'w ' .or. cmd(1).eq.'w!') then
    if (iecho.eq.6) call startClock()
    savecmd = .false.
    call xecw (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("ecw")
elseif ( cmd(1)(1:3).eq.'we ') then
    if (iecho.eq.6) call startClock()
    savecmd = .false.
    call xecw (cmd,ncmd,ierr)
    if (ierr.eq.OK) then
        call xedit (cmd,ncmd,ierr)
    endif
    if (iecho.eq.6) call LogTimeStamp("we")
elseif ( cmd(1)(1:3).eq.'wq ') then
    if (iecho.eq.6) call startClock()
    savecmd = .false.
    call xecw (cmd,ncmd,ierr)
    if (ierr.eq.OK) then
         call xecq (cmd,ncmd,ierr)
    endif
    if (iecho.eq.6) call LogTimeStamp("ecq")
elseif ( cmd(1)(1:2).eq.'r ' .or. &
         cmd(1)(1:3).eq.'res' .or. &
         cmd(1).eq.'r!') then
    if (iecho.eq.6) call startClock()
    savecmd = .false.
    call xecr (cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("ecr")
    if (ierr.ne.OK) then
        goto 1000
    endif
    rewcf = .true.
    goto 100
elseif (cmd(1)(1:4).eq.'xsec') then
        !call BuildElementNeighbourTable3D() ! in meshdatamodule
    call xcrosssec(cmd,ncmd,ierr)
    if (iecho.eq.6) call LogTimeStamp("xsec")
else
    write(stderr,*) "Command ", trim(cmd(1)), " not recognized"
    write(stderr,10)
    ierr = SYNTAX
endif

1000    if (lucmd(ilu).eq.stdin) then
if (ierr.eq.OK .and. savecmd) then
call bufput (cmdln,ierr)
bufedit = .true.
endif
 else
if ((ierr.eq.OK.or.ierr.eq.EXEC) .and. savecmd) &
       call bufput (cmdln,ierr)
if (ierr.eq.SYNTAX .or. ierr.eq.EXEC) then
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

 if (ierr.eq.PRGERR) call exitp(1)

 ierr = OK
 goto 200
 end

!-----------------------------------------------------------------------

subroutine sayhello()

use iomodule,           only: stderr

write(stderr,*)	"This software is protected under the Creative Commons Attribution-NonCommercial-ShareAlike license,"
write(stderr,*) "version 4.0, international: https://creativecommons.org/licenses/by-nc-sa/4.0/"
write(stderr,*) " "
write(stderr,*) "For more information, contact:"
write(stderr,*) "Rob Govers," 
write(stderr,*) "Department of Earth Sciences,"  
write(stderr,*) "Utrecht University, Netherlands." 
write(stderr,*) "E-mail: r.govers@uu.nl"  

return
 end
