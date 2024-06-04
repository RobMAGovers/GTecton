subroutine shcmd (cmd,ncmd,ierr)

use modeldatamodule, only: modeldatactx
use modeldefinition, only: ntanom, numstr, numwnk
use timestepmodule,  only: timestepdatactx, nstep
use iomodule,        only: stderr, stdout, luout
use monitorVars,     only: elmon, nnel, nnp, npmon, nsid, sidemon, monel, monsid, monnp
use errors,          only: exec, ok, syntax

 implicit none
!-pass
 integer ncmd,ierr
 character cmd*(*)
 dimension cmd(*)




!-locl
 integer csystem
 character(len=512) tmpnm,tmpfile
 character(len=600) syscmd
 logical numeric
 integer, external :: lnblk, nextlu, ifopen, indx
integer :: lu, ish, m, n, i
real, external :: chreal

!# sh[ow] parameter
 lu = stdout
 2    format(1x,'Usage: sh[ow] parameter [> file]')
 if (ncmd.lt.2) then
ierr = SYNTAX
write(stderr,2)
goto 1000
 endif
 if (ncmd.ge.3) then
if (cmd(3).eq.'>') then
if (ncmd.lt.4) then
    ierr = SYNTAX
    write(stderr,2)
    goto 1000
endif
lu = luout(1)
ierr = ifopen(lu,cmd(4),'unknown')
if (ierr.ne.0) then
    lu = stdout
    ierr = OK
endif
endif
 endif

 if (cmd(2)(1:2).eq.'mo') then
if (npmon) then
write(lu,4) nnp
 4        format(1x,I3,' monitored nodal points:')
write(lu,6) (monnp(i),i=1,nnp)
 6        format(10(1X,I5,1X)/)
endif
if (elmon) then
write(lu,8) nnel
 8        format(1x,I3,' monitored elements:')
write(lu,6) (monel(i),i=1,nnel)
endif
if (sidemon) then
write(lu,10) nsid
10        format(1x,I3,' monitored element/sides:')
write(lu,12) (monsid(1,i),monsid(2,i),i=1,nsid)
12        format(9(1X,I5,'/',I1)/)
endif
goto 1000
 endif

 if (cmd(2)(1:2).eq.'fe') then
if (lu.eq.stdout) then
lu = nextlu(0)
tmpfile = tmpnm('less')
ierr = ifopen(lu,tmpfile,'unknown')
if (ierr.ne.0) then
    lu = stdout
    ierr = OK
    tmpfile = ' '
endif
else
tmpfile = ' '
endif
call shfe (lu)
if (lnblk(tmpfile).gt.0) then
call closef(lu)
syscmd = 'more '// tmpfile(1:lnblk(tmpfile))
ierr = csystem(syscmd)
if (ierr.ne.0) then
    write(stderr,13) ierr,'(',syscmd(1:lnblk(syscmd)), &
           ')'
13            format(1x,'exit ',I3,': csystem ',80A)
    ierr = EXEC
    goto 1000
endif
syscmd = 'rm -f '//tmpfile(1:lnblk(tmpfile))
ierr = csystem(syscmd)
if (ierr.ne.0) then
    write(stderr,13) ierr,'(',syscmd(1:lnblk(syscmd)), &
           ')'
    ierr = EXEC
    goto 1000
endif
lu = stdout
goto 1000
endif
 endif

 if (cmd(2)(1:2).eq.'ti') then
14        format(1x,'Usage: echo time[(itime)] [> file]')
n = indx(cmd(2),'(')
if (n.le.0) then
call wrealt (-NSTEP,lu)
goto 1000
endif
m = indx(cmd(2),')')
if (m.le.0 .or. m.le.n) then
write(stderr,14)
ierr = SYNTAX
goto 1000
endif
if (cmd(2)(n:m).eq.'(itime)') then
ish = NSTEP
else
if (numeric(cmd(2)(n+1:m-1))) then
    ish = INT(chreal(cmd(2)(n+1:m-1)))
else
    write(stderr,14)
    ierr = SYNTAX
    goto 1000
endif
endif
call wrealt (ish,lu)

goto 1000
 endif

if (cmd(2)(1:2).eq.'ca') then
    if (lu.eq.stdout) then
        lu = nextlu(0)
        tmpfile = tmpnm('less')
        ierr = ifopen(lu,tmpfile,'unknown')
        if (ierr.ne.0) then
            lu = stdout
            ierr = OK
            tmpfile = ' '
        endif
    else
        tmpfile = ' '
    endif
    call wcalct (timestepdatactx%icalc,lu)
    if (lnblk(tmpfile).gt.0) then
        call closef(lu)
        syscmd = 'more '// tmpfile(1:lnblk(tmpfile))
        ierr = csystem(syscmd)
        if (ierr.ne.0) then
            write(stderr,13) ierr,'(',syscmd(1:lnblk(syscmd)),')'
            ierr = EXEC
            goto 1000
        endif
        syscmd = 'rm -f '//tmpfile(1:lnblk(tmpfile))
        ierr = csystem(syscmd)
        if (ierr.ne.0) then
            write(stderr,13) ierr,'(',syscmd(1:lnblk(syscmd)),')'
            ierr = EXEC
            goto 1000
        endif
        lu = stdout
        goto 1000
    endif
endif

 if (cmd(2)(1:2).eq.'pl') then
if (lu.eq.stdout) then
lu = nextlu(0)
tmpfile = tmpnm('less')
ierr = ifopen(lu,tmpfile,'unknown')
if (ierr.ne.0) then
    ierr = OK
    tmpfile = ' '
    lu = stdout
endif
else
tmpfile = ' '
endif
call shplt (lu)
if (lnblk(tmpfile).gt.0) then
call closef(lu)
syscmd = 'more '// tmpfile(1:lnblk(tmpfile))
ierr = csystem(syscmd)
if (ierr.ne.0) then
    write(stderr,13) ierr,'(',syscmd(1:lnblk(syscmd)), &
           ')'
    ierr = EXEC
    goto 1000
endif
syscmd = 'rm -f '//tmpfile(1:lnblk(tmpfile))
ierr = csystem(syscmd)
if (ierr.ne.0) then
    write(stderr,13) ierr,'(',syscmd(1:lnblk(syscmd)), &
           ')'
    ierr = EXEC
    goto 1000
endif
lu = stdout
endif
goto 1000
 endif

 if (cmd(2).eq.'bc') then
if (lu.eq.stdout) then
lu = nextlu(0)
tmpfile = tmpnm('less')
ierr = ifopen(lu,tmpfile,'unknown')
if (ierr.ne.0) then
    lu = stdout
    ierr = OK
    tmpfile = ' '
endif
else
tmpfile = ' '
endif

call shstrs (modeldatactx%STRS,modeldatactx%IELSTR, &
       modeldatactx%ISSIDE,modeldatactx%istr,NUMSTR,lu)

call shwink (modeldatactx%IWELM, modeldatactx%IWSIDE,modeldatactx%IWTIME,NUMWNK,lu)

call shanom (modeldatactx%itanom,lu)

if (lnblk(tmpfile).gt.0) then
call closef(lu)
syscmd = 'more '// tmpfile(1:lnblk(tmpfile))
ierr = csystem(syscmd)
if (ierr.ne.0) then
    write(stderr,13) ierr,'(',syscmd(1:lnblk(syscmd)), &
           ')'
    ierr = EXEC
    goto 1000
endif
syscmd = 'rm -f '//tmpfile(1:lnblk(tmpfile))
ierr = csystem(syscmd)
if (ierr.ne.0) then
    write(stderr,13) ierr,'(',syscmd(1:lnblk(syscmd)), &
           ')'
    ierr = EXEC
    goto 1000
endif
lu = stdout
endif
goto 1000
 endif

 if (cmd(2).eq.'dt') then
if (lu.eq.stdout) then
lu = nextlu(0)
tmpfile = tmpnm('less')
ierr = ifopen(lu,tmpfile,'unknown')
if (ierr.ne.0) then
    lu = stdout
    ierr = OK
    tmpfile = ' '
endif
else
tmpfile = ' '
endif

call shgrps (lu)

if (lnblk(tmpfile).gt.0) then
call closef(lu)
syscmd = 'more '// tmpfile(1:lnblk(tmpfile))
ierr = csystem(syscmd)
if (ierr.ne.0) then
    write(stderr,13) ierr,'(',syscmd(1:lnblk(syscmd)), &
           ')'
    ierr = EXEC
    goto 1000
endif
syscmd = 'rm -f '//tmpfile(1:lnblk(tmpfile))
ierr = csystem(syscmd)
if (ierr.ne.0) then
    write(stderr,13) ierr,'(',syscmd(1:lnblk(syscmd)), &
           ')'
    ierr = EXEC
    goto 1000
endif
lu = stdout
endif
goto 1000
 endif

 write(stderr,995) cmd(2)(1:lnblk(cmd(2)))
  995    format(1x,'show: unknown option "',80A)
 write(stderr,996)
  996    format(1x,'type "? show" for a list of options')

1000    if (lu.ne.stdout) call closef(lu)
 return
 end
!-------------------------------------------------------------------------------
 subroutine shfe (lu)

use meshdatamodule,  only: meshdatactx
use modeldatamodule, only: modeldatactx
use modeldefinition
use timestepmodule,  only: NMPRT, maxit, maxtit, NINTG, &
                           NPRMAT, NTCALC, NTPRT, FaultVelo
use plotControl,     only: jShell

implicit none
!-pass
 integer lu


 write(lu,10) meshdatactx%nvglobal, &
              meshdatactx%neglobal, &
              JSHELL,iopt
10   format(1x, &
 'Number of nodal points                 (NUMNP) = ',I12/1x, &
 'Number of elements                     (NUMEL) = ',I12/1x, &
 'Spherical shell element switch        (ISHELL) = ',I12/5x, &
 '= 0, plane elements'/5x, &
 '= 1, spherical elements'/1x, &
 'Two-dimensional stress state            (IOPT) = ',I12/5x, &
 '= 0, plane strain'/5x, &
 '= 1, plane stress')
if (JSHELL.eq.1) write(lu,15) RADIUS
15   format(1x,&
 'Radius                                         = ',1PG14.6)

write(lu,20) MODE,ISTART
20   format(1x, &
 'Solver switch                           (MODE) = ',I12/5x, &
 '= 0, mechanical data check mode'/5x, &
 '= 1, mechanical rank check mode',/5x, &
 '= 2, mechanical solution only'/5x, &
 '= 3, thermal data check mode'/5x, &
 '= 4, thermal rank check mode'/5x, &
 '= 5, thermal solution only'/5x, &
 '= 6, full thermo-mechanical solution mode'/1x, &
 'Restart switch                        (ISTART) = ',I12/5x, &
 '= 0, no restart'/5x, &
 '> 0, no restart, write restart file each ISTART steps'/5x, &
 '< 0, restart, write restart file each -ISTART steps')

write(lu,22) NELAST
22   format(1x, &
 'Skip elastic solution switch          (NELAST) = ',I12/5x, &
 '= 0, compute elastic solution'/5x, &
 '= 1, skip elastic solution')

write(lu,30) NUMAT
if (MODE.eq.2 .or. MODE.eq.6) then
    write(lu,32) NOCOMPR
    if (NINTG.gt.0) write(lu,40) ICVIS,IVLIM
endif

30   format(1x, &
 'Number of material properties          (NUMAT) = ',I12)
32   format(1x, &
 'Incompressible elasticity switch     (NOCOMPR) = ',I12/5x, &
 '= 0, compressible'/5x, &
 '= 1, incompressible')
40   format(1x, &
 'Material property update switch        (ICVIS) = ',I12/5x, &
 '= 0, constant viscosity'/5x, &
 'first digit; powerlaw viscosity'/5x, &
 'second digit; visco-plasticity'/5x, &
 'third digit; viscosity weakening'/1x, &
 'fourth digit; yield stress weakening'/1x, &
 'Viscosity limit switch                 (IVLIM) = ',I12/5x, &
 '= 0, no limit'/5x, &
 '= 1, static viscosity limits'/5x, &
 '= 2, time step controlled viscosity limts')

write(lu,50) IPOINT
50   format(1x, &
 'Numerical integration accuracy        (IPOINT) = ',I12/5x, &
 '= 1, low order'/5x, &
 '= 4, accurate')
if (MODE.le.2.or.MODE.eq.6) write(lu,52) INCOMP
52   format(1x, &
 'Incompatible mode for quad elements   (INCOMP) = ',I12/5x, &
 '= 0, off'/5x, &
 '= 1, on')
if (MODE.ge.3) write(lu,54) ITMODE,IDIFT
54   format(1x, &
 'Temperatures initialization switch    (ITMODE) = ',I12/5x, &
 '= 0, steady state'/5x, &
 '= 1, read from file'/5x, &
 '= 2,3 compute with function "temp0.F"'/1x, &
 'Differential temperatures switch       (IDIFT) = ',I12/5x, &
 '= 0, no temperature diffs accross faults'/5x, &
 '= 1, allow differential temperatures')
if (MODE.ge.3.and.MODE.le.5) write(lu,56) IADVEC
56   format(1x, &
 'Advecting grid switch                 (IADVEC) = ',I12/5x, &
 '= 0, no advection'/5x, &
 '= 1, advection')
write(lu,60) NINTG
60   format(1x, &
 'Number of time step groups             (NINTG) = ',I12)
if (NINTG.gt.0) then
    if (MODE.eq.2.or.MODE.eq.6) then
#ifdef EXPERIMENTAL_CYCLIC
        if (UseCyclicSlip) then
            write(lu,72) 1
        else
            write(lu,72) 0
        endif
#endif
    endif
    if (MODE.eq.5.or.MODE.eq.6) then
        write(lu,80) NTCALC,NTPRT,MAXTIT 
    endif
    if (MODE.eq.2.or.MODE.eq.6) then
        write(lu,90) LGDEF,IRESDU
    endif
endif

write(lu,70) NMPRT,MAXIT,NPRMAT,FAULTVELO
70   format(1x, &
 'Mechanical solution output switch      (NMPRT) = ',I12/5x, &
 '>= 0, number of outputs'/5x, &
 '<0  , output every -NMPRT time steps'/1x, &
 'Maximum number of time steps without'/1x, &
 'reforming mechanical stiffness matrix  (MAXIT) = ',I12/1x, &
 'Mechanical matrix diagonal outputs    (NPRMAT) = ',I12/5x, &
 '>= 0, number of outputs'/5x, &
 '<0  , output every -NPRMAT time steps'/1x, &
 'Velocity output switch             (FAULTVELO) = ',I12/5x, &
 ' = 0, Fault input read as displacement'/5x, &
 ' = 1, Fault input read as velocity')
72  format(5x, &
 'Cyclic locking/unlocking              (CYCLIC) = ',I12)
80   format(1x, &
 'Thermal calculation switch            (NTCALC) = ',I12/5x, &
 '>= 0, number of calculations'/5x, &
 '<0  , update every -NTCALC time steps'/1x, &
 'Thermal solution output switch         (NTPRT) = ',I12/5x, &
 '>= 0, number of outputs'/5x, &
 '<0  , output every -NTPRT time steps'/1x, &
 'Maximum number of time steps without'/1x, &
 'reforming thermal stiffness matrix    (MAXTIT) = ',I12)
90   format(1x, &
 'Large deformation formulation switch   (LGDEF) = ',I12/5x, &
 '= 0, off'/5x, &
 '= 1, on'/1x, &
 'Residual force update switch          (IRESDU) = ',I12/5x, &
 '= 0,1, off'/5x, &
 '= 2, on')
if (MODE.le.2 .or. MODE.eq.6) write(lu,100) NPRE,IGRAV
  100   format(1x, &
 'Number of pre-stress loads              (NPRE) = ',I12/1x, &
 'Gravity switch                         (IGRAV) = ',I12/5x, &
 '= 0, off'/5x, &
 '= 1, on')
if ((MODE.le.2 .or. MODE.eq.6).or. &
 (IADVEC.ne.0.and.MODE.ge.3.and.MODE.le.5)) write(lu,110) NUMFN
  110   format(1x, &
 'Number of faulted nodes                (NUMFN) = ',I12)
if (MODE.le.2 .or. MODE.eq.6) write(lu,120) NUMSN,NUMSLP
  120   format(1x, &
 'Number of slippery nodes               (NUMSN) = ',I12/1x, &
 'Number of slippery node entries       (NUMSLP) = ',I12)
if (MODE.le.2 .or. MODE.eq.6) write(lu,130) NLINK,NUMSTR,NUMPR, &
 NUMWNK
  130   format(1x, &
 'Number of linked nodal dofs            (NLINK) = ',I12/1x, &
 'Number of stress boundary conditions  (NUMSTR) = ',I12/1x, &
 'Number of pressure boundary conditions (NUMPR) = ',I12/1x, &
 'Number of Winkler pressures           (NUMWNK) = ',I12)
if (MODE.le.2 .or. MODE.eq.6) write(lu,140) NSURF,NSED
  140   format(1x, &
 'Number of surface nodes                (NSURF) = ',I12/1x, &
 'Sediment transport loads switch         (NSED) = ',I12/5x, &
 '= 0, off'/5x, &
 '= 1, on')
if (MODE.ge.3) write(lu,150) NFLX,NTWINK,NTANOM
  150   format(1x, &
 'Number of heat flow boundary conditions (NFLX) = ',I12/1x, &
 'Number of thermal Winkler bc''s        (NTWINK) = ',I12)
if (MODE.eq.6.and.NINTG.gt.0) write(lu,155) ITPRES
  155   format(1x, &
 'Temperature induced body force switch (ITPRES) = ',I12/5x, &
 '= 0, on'/5x, &
 '= 1, off')
if (MODE.ge.3.and.NINTG.gt.0) write(lu,160) NTANOM,ITECT
  160   format(1x, &
 'Temperature anomaly switch            (NTANOM) = ',I12/5x, &
 '>0, number of anomaly inputs'/5x, &
 '<0, add anomaly in "addan.F" at time -NTANOM'/1x, &
 'Anomaly induced body force switch      (ITECT) = ',I12/5x, &
 '= 0, on if ITPRES=0'/5x, &
 '= 1, off')

return
end
!-------------------------------------------------------------------------------
double precision function TIMEQV (MAXSTP,DELT,NINTG,ITIME)

!    Finds real time equivalent of integer time ITIME

 implicit none
!-pass
 integer NINTG,ITIME
integer          :: maxstp
double precision :: delt
 dimension MAXSTP(*),DELT(*)
!-local
integer :: it1, intg, it0
double precision :: t1, t0


 TIMEQV = -1d0
 if (ITIME.le.0) TIMEQV = 0d0
 if (NINTG.le.0) return

 IT1  = 0
 T1   = 0d0
 INTG = 1
100    if (INTG.gt.NINTG) goto 300
 IT0 = IT1
 T0  = T1
 IT1 = IT0 + MAXSTP(INTG)
 T1  = T0 + DBLE(MAXSTP(INTG))*DELT(INTG)
 if ( (ITIME-IT0)*(ITIME-IT1) .le. 0 ) goto 200
 INTG = INTG + 1
 goto 100

200    TIMEQV = T0 + DBLE(ITIME-IT0) * (T1-T0)/DBLE(IT1-IT0)
300    return
 end
!-------------------------------------------------------------------------------
 subroutine shplt (lu)

use timestepmodule,  only: DELTP,ALFAP, nstep
use modeltopology,   only: NSD
use meshdatamodule,  only: meshdatactx
use iomodule,        only: ifein, ifeout, ascout, mout, tout, append
use extrema,         only: stressPositiveSmallest,       stressPositiveBiggest, &
                           stressNegativeSmallest,       stressNegativeBiggest, &
                           strainPositiveSmallest,       strainPositiveBiggest, &
                           strainNegativeSmallest,       strainNegativeBiggest, &
                           strainRatePositiveSmallest,   strainratePositiveBiggest, &
                           strainRateNegativeSmallest,   strainRateNegativeBiggest, &
                           displacementPositiveSmallest, displacementPositiveBiggest, &
                           displacementNegativeSmallest, displacementNegativeBiggest
use postProcessData, only: GMT, deffac, deformed
use monitorVars,     only: elmon, npmon
use plotControl        !only: a lot


! show plot-program related parameters

 implicit none
!-pass
 integer lu



!-locl
double precision :: fac
logical :: lfein,lfeout,arcont
integer :: i
!
lfein  = (ifein .ne.0)
lfeout = (ifeout.ne.0)
arcont= (iarcont.ne.0)

write(lu,10) lfein,lfeout,ascout,append,GMT,NSTEP
10   format(1x,/'INPUT/OUTPUT:'/1x, &
 'Finite element input  file has been read = ',L1/1x, &
 'Finite element output file has been read = ',L1/1x, &
 'Ascii output files                       = ',L1/1x, &
 'Append output data to file               = ',L1/1x, &
 '"GMT" output format                      = ',L1/1x, &
 'Time                               itime = ',I5 )

write(lu,20) deformed,grnum,npnum,viscnu
20   format(1x,/'MESH:'/1x, &
 'Deformed                = ',L1/1x, &
 'Grid numbers            = ',L1/1x, &
 'Node numbers            = ',L1/1x, &
 'Viscosity numbers       = ',L1/)

write(lu,30) elmon,npmon
30   format(1x,'MONITORING:'/1x, &
 'monitoring elements     = ',L1/1x, &
 'monitoring nodes        = ',L1/)

fac = 0d0
if (deformed) then
    fac = deffac
endif

!obsolete
! write(lu,40) Sscale,Ssb,Escale,Esb,Rscale,Rsb,Dscale,Dsb,Fscale,Fsb,fac
40   format(1x,'SCALING:'/1x, &
 'Stress:       scale factor  = ',1PG12.4/1x, &
 '              scale bar     = ',1PG12.4/1x, &
 'Strain:       scale factor  = ',1PG12.4/1x, &
 '              scale bar     = ',1PG12.4/1x, &
 'Strain rate:  scale factor  = ',1PG12.4/1x, &
 '              scale bar     = ',1PG12.4/1x, &
 'Displacement: scale factor  = ',1PG12.4/1x, &
 '              scale bar     = ',1PG12.4/1x, &
 'Heat flux:    scale factor  = ',1PG12.4/1x, &
 '              scale bar     = ',1PG12.4/1x, &
 'Deformation:  amplification = ',1PG12.4/)

if (ifein.ne.0) then
    write(lu,50) (crdmin(i),crdmax(i),i=1,NSD)
50   format(1x,'RANGES:'/1x, &
 'Coordinates:  ',1Pg9.1,' <= X <= ',1Pg9.1/1x, &
 14X,1Pg9.1,' <= Y <= ',1Pg9.1/1x, &
 14X,1Pg9.1,' <= Z <= ',1Pg9.1)
endif

if (Mout) then
    write(lu,60) displacementPositiveSmallest, displacementPositiveBiggest, &
                 stressNegativeBiggest,        stressNegativeSmallest, &
                 stressPositiveSmallest,       stressPositiveBiggest, &
                 strainNegativeBiggest,        strainNegativeSmallest, &
                 strainPositiveSmallest,       strainPositiveBiggest, &
                 strainRateNegativeBiggest,    strainRateNegativeSmallest, &
                 strainRatePositiveSmallest,   strainratePositiveBiggest

60   format(1x, &
 'Displacement:                              ' , &
                   (1Pg9.1,' <= D <= ',1Pg9.1,2x)/1x, &
 'Stress:       ',2(1Pg9.1,' <= S <= ',1Pg9.1,2x)/1x, &
 'Strain:       ',2(1Pg9.1,' <= S <= ',1Pg9.1,2x)/1x, &
 'Strain rate:  ',2(1Pg9.1,' <= R <= ',1Pg9.1,2x))
endif

if (Tout) then
    write(lu,70) tempmn,tempmx, fluxmn(2),fluxmx(2),fluxmn(1),fluxmx(1)
70   format(1x, &
 'Temperature:  ',1Pg9.1,' <= T <= ',1Pg9.1/1x, &
 'Heat flux:    ',2(1Pg9.1,' <= F <= ',1Pg9.1,2x))
endif

if (Mout.or.Tout) then
    write(lu,'(1X)')
endif

write(lu,80) vislog,viscli,ratlog
80   format(1x,'TRANSFORMATIONS AND SETTINGS:'/1x, &
 'Logarithmic viscosities    = ',L1/1x, &
 'Linearized  viscosities    = ',L1/1x, &
 'Logarithmic strain rates   = ',L1)

write(lu,120) DELTP,ALFAP
  120   format(1x,'TIME STEP:'/1x, &
 'Time step size = ',1PG12.4/1x, &
 'Alpha          = ',1PG12.4)

write(lu,130) hatsize,hatangle,arcont,sysi,sysc,sbhor
  130   format(1x,'LAYOUT:'/1x, &
 'Arrow: hatsize       = ',1PG12.4/1x, &
 '       top angle     = ',1PG12.4,' degrees'/1x, &
 '       contoured     = ',L1/1x, &
 'Symbol size          = ',1PG12.4/1x, &
 'Scaling of bc symbol = ',L1/1x, &
 'Horizontal scale bar = ',L1/)

!write(lu,140) ielbas,meshdatactx%neglobal,ielstp,inpbas,meshdatactx%nvglobal,inpstp
!  140   format(1x,'INCREMENTING:'/1x, &
! 'element number ElNo:'/1x, &
! '    for (ElNo=',I5,'; ElNo<=',I5,'; ElNo=ElNo+',I5,')'/1x, &
! 'nodal point number NpNo:'/1x, &
! '    for (NpNo=',I5,'; NpNo<=',I5,'; NpNo=NpNo+',I5,')'/)

return
end
!-------------------------------------------------------------------------------
double precision function RealTime (itime)

USE CONSTANTS,      only: year
use timestepmodule, only: timestepdatactx, nintg

implicit none

!-pass
integer itime

!-local
integer :: n, igrp, i
double precision :: x

!-init
n    = 0
igrp = 1
i    = 1
x    = 0d0
RealTime=0d0

if (NINTG.le.0) then
    RealTime=-1.D0
    goto 300
endif

100    if (igrp.gt.NINTG) then
    RealTime=-2.D0
    goto  300
endif

200    if (i.gt.timestepdatactx%MAXSTP(igrp)) then
    i    = 1
    igrp = igrp + 1
    goto 100
endif

n = n + 1
x = x + timestepdatactx%DELT(igrp)

if (n.eq.itime) then
    RealTime = x/YEAR
    goto 300
endif

i = i + 1
goto 200

300    return
end
!-------------------------------------------------------------------------------
subroutine shstrs (STRS,IELSTR,ISSIDE,ISTR, &
                   NUMSTR,lu)

USE ALGEBRA
use timestepmodule, only: timestepdatactx, nintg
use modeltopology,  only: ndof, nstr
use label,          only: labels

implicit none

!-pass
integer          :: NUMSTR, ielstr, isside, istr, lu
double precision :: strs
dimension        :: IELSTR(*),ISSIDE(*),ISTR(2,*), &
                    STRS(6,*)
!-locl
logical          :: updateStress
double precision :: SUM(6)
integer          :: i, iel, jel, jside, n, nn, nprt
integer          :: iside
double precision :: time
double precision, external :: realTime
!-init
nprt=-1

if (NUMSTR.le.0) return

do n=1,NUMSTR
    if (nprt.ge.51 .or. nprt.lt.0) then
        if (nprt.gt.0) then
            write(lu,2)
        endif
        write(lu,1) (LABELS(i),i=1,6)
        write(lu,2)
        nprt=0
    endif

    iel = IELSTR(n)
    iside = ISSIDE(n)

    call clear (SUM,6,"SUM")
    call UPDATE (SUM,STRS(1,n),6)

    if (n.gt.1) then
        do nn=1,n-1
            jel = IELSTR(nn)
            jside = ISSIDE(nn)
            if (iel.eq.jel .and. iside.eq.jside) then
                updateStress = .false.
            else
                updateStress = .true.
            endif
        enddo
    endif

    if (updateStress) then

        time = RealTime (ISTR(2,n))

        if (time.lt.0.0) then
            write(lu,4) iel,iside
        else
            write(lu,3) iel,iside,time,(SUM(i),i=1,6)
        endif

        nprt=nprt+1

        if (n.lt.NUMSTR) then
            do nn=n+1,NUMSTR

                jel = IELSTR(nn)
                jside = ISSIDE(nn)

                if (iel.eq.jel .and. iside.eq.jside) then

                    call UPDATE (SUM,STRS(1,nn),6)
                    time = RealTime (ISTR(2,nn))

                    if (time.lt.0.0) then
                        write(lu,4) iel,iside
                    else
                        write(lu,3) iel,iside,time,(SUM(i),i=1,6)
                    endif
                    nprt=nprt+1
                endif
            enddo
        endif

    endif
enddo


return

 1    format(1X,'ELEM#',1X,'SIDE',2X,'AFTER YEAR',6X,6(A4,8X))
 2    format(/)
 3    format(1x,I5,2X,I1,4X,1PG10.2,6(2X,1PG10.2))
 4    format(1x,I5,2x,I1,4X,'NO STRESS')
! ELEM# SIDE  AFTER TIME      SXX         SYY          SXY
! IIIII  8    FFFFFFFFFF  XXXXXXXXXX  XXXXXXXXXX  XXXXXXXXXXX
end subroutine
!-------------------------------------------------------------------------------

! todo test this
subroutine shwink (IWELM,IWSIDE,IWTIME, NUMWNK,lu)

use timestepmodule, only: timestepdatactx, nintg

implicit none

!-pass
integer :: NUMWNK
integer :: iwelm, iwside, iwtime
dimension IWELM(*),IWSIDE(*),IWTIME(*)

!-locl
integer :: lu, i, iel, iside, itime, m, n
double precision :: time
double precision, external :: realTime

if (NUMWNK.eq.0) then
    return
endif

n = 1
!    find an element that has not been processed
100    if (n.le.ABS(NUMWNK)) then

    if (IWELM(n).gt.0) then
        goto 200
    endif

    n = n + 1
    goto 100

    200    iel = IWELM(n)
    iside = IABS(IWSIDE(n))
    m = n
300 if (m.le.IABS(NUMWNK)) then
        if (IWELM(m).eq.iel .and. IABS(IWSIDE(m)).eq.iside) then

            if (IWTIME(m).gt.0) then
                itime = IWTIME(m)-1
            else
                itime = -IWTIME(m)
            endif

            if (itime.gt.0) then
                time = RealTime (itime)
            else
                time = 0d0
            endif

            if (time.lt.0) then
                if (time.lt.-1.1d0) then
                    if (IWTIME(m).gt.0) then
                        write(lu,1) iel,iside
 1                       format(1x,'Winkler pressure on element ',I5, &
                     ' (side ',I1,') is never applied')
                    else
                        write(lu,2) iel,iside
 2                       format(1x,'Winkler pressure on element ',I5, &
                         ' (side ',I1,') is always applied')
                    endif
                endif
            else
                if (IWTIME(m).gt.0) then
                    write(lu,3) iel,iside,time
 3                   format(1x,'Winkler pressure on element ',I5, &
                 ' (side ',I1,') applied at ',1PG14.6,' year')
                else
                    write(lu,4) iel,iside,time
 4                   format(1x,'Winkler pressure on element ',I5, &
                 ' (side ',I1,') removed at ',1PG14.6,' year')
                endif
            endif
            IWELM(m) = -IWELM(m)
        endif
        m = m + 1
        goto 300
    endif
    n = n + 1
    goto 100

endif

do i=1,ABS(NUMWNK)
       IWELM(i) = IABS(IWELM(i))
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine shanom (ITANOM,lu)

use modeldefinition, only: ntanom
use timestepmodule,  only: timestepdatactx, &
                           nintg

implicit none

!-pass
integer :: itanom(3,ntanom)
integer :: lu


!-locl
integer :: i, itime, mode, node
double precision :: time

double precision, external :: realTime

if (NTANOM.le.0) then
    return
endif

do i=1,NTANOM

    node  = ITANOM(1,i)
    mode  = ITANOM(2,i)
    itime = ITANOM(3,i)

    time  = RealTime (itime)

    if (time.lt.0d0) then
        if (mode.eq.0) then
            write(lu,1) node
 1               format(1x,'Temperature of node ',I5, &
             ' is never replaced')
        else
            write(lu,2) node
 2               format(1x,'Temperature of node ',I5, &
             ' is never increased')
        endif
    else
        if (mode.eq.0) then
            write(lu,3) node,time
 3               format(1x,'Temperature of node ',I5, &
             ' is replaced at ',1PG14.6,' year')
        else
            write(lu,4) node,time
 4               format(1x,'Temperature of node ',I5, &
             ' is increased at ',1PG14.6,' year')
        endif
    endif
enddo

return
end
!-------------------------------------------------------------------------------
subroutine shgrps (lu)

USE CONSTANTS, only: year
use timestepmodule, only: timestepdatactx, nintg


implicit none
!-pass
integer :: lu
!-local
double precision :: dt, t0, t1
integer :: igrp

logical :: doMechanic, doThermal



if (NINTG.le.0) then
    return
endif

t0=0d0

do igrp=1,NINTG

    dt = timestepdatactx%DELT(igrp)/YEAR

    t1 = t0 + dt*timestepdatactx%MAXSTP(igrp)

    ! decide whether to show thermal or mechanical timesteps
    ! for now, just show mechanical only.   

    doMechanic = allocated(timestepdatactx%ALFA)
    doThermal  = allocated(timestepdatactx%TALF)


    if (doMechanic .and. doThermal) then
        write(lu,1) t0, t1, dt, timestepdatactx%ALFA(igrp), &
                                timestepdatactx%TALF(igrp)

 1       format(1x,'from ',1PG10.2,' to ',1PG10.2,' year: dt= ', &
     1PG14.6,' year (AlphaM=',0PF4.2,' AlphaT=',F4.2,')')
    else if (doMechanic) then
        ! and no mechanic
                write(lu,2) t0, t1, dt, timestepdatactx%ALFA(igrp)

 2       format(1x,'from ',1PG10.2,' to ',1PG10.2,' year: dt= ', &
     1PG14.6,' year (AlphaM=',0PF4.2,')')
    else
        ! only doThermal
                write(lu,3) t0, t1, dt, timestepdatactx%TALF(igrp)

 3       format(1x,'from ',1PG10.2,' to ',1PG10.2,' year: dt= ', &
     1PG14.6,' year (AlphaT=',0PF4.2,')')
    endif

    t0=t1
enddo


return
end

