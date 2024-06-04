subroutine shcmd (cmd,ncmd,ierr)

use modeltopology,   only: ndof, nstr
use modeldefinition, only: ntanom, numwnk, numstr
use timestepmodule,  only: timestepdatactx, outputcontroldatactx, nstep, nintg
use modeldatamodule, only: modeldatactx
use iomodule,        only: stderr, stdout, luout
use monitorVars,     only: elmon, nnel, npmon, nnp, monnp, monel
use errors,          only: ok, exec, syntax

implicit none
!-pass
integer ncmd,ierr
character cmd*(*)
dimension cmd(*)




!-locl
integer csystem
character(len=512) syscmd,tmpnm,tmpfile
logical numeric
integer, external :: lnblk, nextlu, ifopen, indx
double precision,    external :: chreal
integer :: i, m, n, lu, ish

!# sh[ow] parameter
 lu = stdout
 1    format(1x,'Usage: sh[ow] parameter [> file]')

if (ncmd.lt.2) then
    ierr = SYNTAX
    write(stderr,1)
    goto 1000
endif

if (ncmd.ge.3) then
    if (cmd(3).eq.'>') then
        if (ncmd.lt.4) then
            ierr = SYNTAX
            write(stderr,1)
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
        write(lu,2) nnp
 2        format(1x,I3,' monitored nodal points:')
        write(lu,3) (monnp(i),i=1,nnp)
 3        format(10(1X,I8))
    endif
    if (elmon) then
        write(lu,4) nnel
 4        format(1x,I3,' monitored elements:')
        write(lu,3) (monel(i),i=1,nnel)
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
        call wrealt (outputcontroldatactx%IMPRINT, &
                     outputcontroldatactx%ITPRINT, &
                     timestepdatactx%MAXSTP, &
       timestepdatactx%DELT,NINTG,-NSTEP,lu)
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
    call wrealt (outputcontroldatactx%IMPRINT, &
                     outputcontroldatactx%ITPRINT, &
                           timestepdatactx%MAXSTP, &
                timestepdatactx%DELT, &
       NINTG,ish,lu)
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

    call shstrs (modeldatactx%STRS, &
                 modeldatactx%ielstr, &
                 modeldatactx%ISSIDE, &
                 modeldatactx%istr,NUMSTR,lu)

    call shwink (modeldatactx%IWELM, &
                 modeldatactx%IWSIDE, &
                 modeldatactx%IWTIME,NUMWNK,lu)

    call shanom (modeldatactx%ITANOM, lu)

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

    call showTimestepGroups (lu)

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
use modeltopology,   only: nen

implicit none
!-pass
integer lu

!

write(lu,10) NEN, meshdatactx%nvglobal, meshdatactx%neglobal
10 format(1x, &
 'Number of nodes per element              (NEN) = ',I8/1x, &
 'Number of nodal points                 (NUMNP) = ',I8/1x, &
 'Number of elements                     (NUMEL) = ',I8)

write(lu,20) MODE,ISTART
20 format(1x, &
 'Solver switch                           (MODE) = ',I8/5x, &
 '= 0, mechanical data check mode'/5x, &
 '= 1, mechanical rank check mode',/5x, &
 '= 2, mechanical solution only'/5x, &
 '= 3, thermal data check mode'/5x, &
 '= 4, thermal rank check mode'/5x, &
 '= 5, thermal solution only'/5x, &
 '= 6, full thermo-mechanical solution mode'/1x, &
 'Restart switch                        (ISTART) = ',I8/5x, &
 '= 0, no restart'/5x, &
 '> 0, no restart, write restart file each ISTART steps'/5x, &
 '< 0, restart, write restart file each -ISTART steps')

write(lu,22) NELAST
22 format(1x, &
 'Skip elastic solution switch          (NELAST) = ',I8/5x, &
 '= 0, compute elastic solution'/5x, &
 '= 1, skip elastic solution'/5x, &
 '= 2, initialize stresses with pre-stresses')

write(lu,30) NUMAT
30 format(1x, &
 'Number of material properties          (NUMAT) = ',I8)

if (MODE.eq.2 .or. MODE.eq.6) then
    write(lu,32) NOCOMPR
32  format(1x, &
     'Incompressible elasticity switch     (NOCOMPR) = ',I8/5x, &
     '= 0, compressible'/5x, &
     '= 1, incompressible')
    if (NINTG.gt.0) then
        write(lu,40) ICVIS,IVLIM
40      format(1x, &
         'Material property update switch        (ICVIS) = ',I8/5x, &
         '= 0, constant viscosity'/5x, &
         'first digit powerlaw viscosity'/5x, &
         'second digit visco-plasticity'/5x, &
         'third digit viscous weakening'/ &
         1x, &
         'Viscosity limit switch                 (IVLIM) = ',I8/5x, &
         '= 0, no limit'/5x, &
         '= 1, static viscosity limits'/5x, &
         '= 2, time step controlled viscosity limts')
    endif
endif

write(lu,50) IPOINT
50   format(1x, &
 'Numerical integration accuracy        (IPOINT) = ',I8/5x, &
 '= 1, low order'/5x, &
 '> 1, accurate')

if (MODE.ge.3) then
    write(lu,54) ITMODE,IDIFT
54  format(1x, &
     'Temperatures initialization switch    (ITMODE) = ',I8/5x, &
     '= 0, steady state'/5x, &
     '= 1, read from file'/5x, &
     '= 2,3 compute with function "temp0.F"'/1x, &
     'Differential temperatures switch       (IDIFT) = ',I8/5x, &
     '= 0, no temperature diffs accross faults'/5x, &
     '= 1, allow differential temperatures')
endif

if (MODE.ge.3.and.MODE.le.5) then
    write(lu,56) IADVEC
endif

56   format(1x, &
 'Advecting grid switch                 (IADVEC) = ',I8/5x, &
 '= 0, no advection'/5x, &
 '= 1, advection')

write(lu,60) NINTG

60   format(1x, &
 'Number of time step groups             (NINTG) = ',I8)

write(lu,70) NMPRT,MAXIT,NPRMAT,FAULTVELO
70   format(1x, &
 'Mechanical solution output switch      (NMPRT) = ',I8/5x, &
 '>= 0, number of outputs'/5x, &
 '<0  , output every -NMPRT time steps'/1x, &
 'Maximum number of time steps without'/1x, &
 'reforming mechanical stiffness matrix  (MAXIT) = ',I8/1x, &
 'Mechanical matrix diagonal outputs    (NPRMAT) = ',I8/5x, &
 '>= 0, number of outputs'/5x, &
 '<0  , output every -NPRMAT time steps'/1x, &
 'Fault velocity switch              (FaultVelo) = ',I8/5x, &
 ' = 0, Fault boundary condition = displacement'/5x, &
 ' = 1, Fault boundary condition = velocity')

#ifdef EXPERIMENTAL_CYCLIC
if (NINTG.gt.0 .and.  (MODE.eq.2.or.MODE.eq.6)) then
    if (UseCyclicSlip) then
        write(lu,72) 1
    else
        write(lu,72) 0
    endif
72  format(1x, &
     'Cyclic locking/unlocking              (CYCLIC) = ',I12/5x, &
     ' = 0, off'/5x, &
     ' = 1, on')
endif
#endif

if (NINTG.gt.0) then
    if (MODE.eq.5.or.MODE.eq.6) then
        write(lu,80) NTCALC,NTPRT,MAXTIT
    endif
    if (MODE.eq.2.or.MODE.eq.6) then
        write(lu,90) LGDEF,IRESDU
    endif
endif  

80   format(1x, &
 'Thermal calculation switch            (NTCALC) = ',I8/5x, &
 '>= 0, number of calculations'/5x, &
 '<0  , update every -NTCALC time steps'/1x, &
 'Thermal solution output switch         (NTPRT) = ',I8/5x, &
 '>= 0, number of outputs'/5x, &
 '<0  , output every -NTPRT time steps'/1x, &
 'Maximum number of time steps without'/1x, &
 'reforming thermal stiffness matrix    (MAXTIT) = ',I8)

90   format(1x, &
 'Large deformation formulation switch   (LGDEF) = ',I8/5x, &
 '= 0, off'/5x, &
 '= 1, on'/1x, &
 'Residual force update switch          (IRESDU) = ',I8/5x, &
 '= 0,1, off'/5x, &
 '= 2, on')

if (MODE.le.2 .or. MODE.eq.6) write(lu,100) NPRE,IGRAV

  100   format(1x, &
 'Number of pre-stress loads              (NPRE) = ',I8/1x, &
 'Gravity switch                         (IGRAV) = ',I8/5x, &
 '= 0, off'/5x, &
 '= 1, on')

if ((MODE.le.2 .or. MODE.eq.6).or. &
 (IADVEC.ne.0.and.MODE.ge.3.and.MODE.le.5)) write(lu,110) NUMFN

  110   format(1x, &
 'Number of faulted nodes                (NUMFN) = ',I8)

if (MODE.le.2 .or. MODE.eq.6) write(lu,120) NUMSN,NUMSLP

  120   format(1x, &
 'Number of slippery nodes               (NUMSN) = ',I8/1x, &
 'Number of slippery node entries       (NUMSLP) = ',I8)

if (MODE.le.2 .or. MODE.eq.6) write(lu,130) NLINK,NUMSTR,NUMPR, &
 NUMWNK

  130   format(1x, &
 'Number of linked nodal dofs            (NLINK) = ',I8/1x, &
 'Number of stress boundary conditions  (NUMSTR) = ',I8/1x, &
 'Number of pressure boundary conditions (NUMPR) = ',I8/1x, &
 'Number of Winkler pressures           (NUMWNK) = ',I8)

if (MODE.le.2 .or. MODE.eq.6) write(lu,140) NSURF,NSED

  140   format(1x, &
 'Number of surface nodes                (NSURF) = ',I8/1x, &
 'Sediment transport loads switch         (NSED) = ',I8/5x, &
 '= 0, off'/5x, &
 '= 1, on')

if (MODE.ge.3) write(lu,150) NFLX,NTWINK,NTANOM

  150   format(1x, &
 'Number of heat flow boundary conditions (NFLX) = ',I8/1x, &
 'Number of thermal Winkler bc''s        (NTWINK) = ',I8)

if (MODE.eq.6.and.NINTG.gt.0) write(lu,155) ITPRES

  155   format(1x, &
 'Temperature induced body force switch (ITPRES) = ',I8/5x, &
 '= 0, on'/5x, &
 '= 1, off')

if (MODE.ge.3.and.NINTG.gt.0) write(lu,160) NTANOM,ITECT

  160   format(1x, &
 'Temperature anomaly switch            (NTANOM) = ',I8/5x, &
 '>0, number of anomaly inputs'/5x, &
 '<0, add anomaly in "addan.F" at time -NTANOM'/1x, &
 'Anomaly induced body force switch      (ITECT) = ',I8/5x, &
 '= 0, on if ITPRES=0'/5x, &
 '= 1, off')
!
 return
 end
!-------------------------------------------------------------------------------
function TIMEQV (ITIME)

!    Finds real time equivalent of integer time ITIME

use timestepmodule, only: timestepdatactx, nintg

implicit none

!-pass
integer ITIME

!-local
integer          :: IT0, IT1, IntG
double precision :: timeQV, T0, T1

TIMEQV = -1.d0

if (ITIME.le.0) then
    TIMEQV = 0.d0
endif

if (NINTG.le.0) then
    return
endif

IT1  = 0
T1   = 0.d0
INTG = 1

100    if (INTG.gt.NINTG) goto 300

IT0 = IT1
T0  = T1
IT1 = IT0 + timestepdatactx%MAXSTP(INTG)
T1  = T0 + DBLE(timestepdatactx%MAXSTP(INTG))*timestepdatactx%DELT(INTG)

if ( (ITIME-IT0)*(ITIME-IT1) .le. 0 ) goto 200

INTG = INTG + 1

goto 100

200    TIMEQV = T0 + DBLE(ITIME-IT0) * (T1-T0)/DBLE(IT1-IT0)
300    return

end
!-------------------------------------------------------------------------------
 subroutine shplt (lu)

! show plot-program related parameters

use timestepmodule,  only: DELTP,ALFAP, nstep

use extrema,         only: stressPositiveSmallest,       stressPositiveBiggest, &
                           stressNegativeSmallest,       stressNegativeBiggest, &
                           strainPositiveSmallest,       strainPositiveBiggest, &
                           strainNegativeSmallest,       strainNegativeBiggest, &
                           strainRatePositiveSmallest,   strainratePositiveBiggest, &
                           strainRateNegativeSmallest,   strainRateNegativeBiggest, &
                           displacementPositiveSmallest, displacementPositiveBiggest, &
                           displacementNegativeSmallest, displacementNegativeBiggest

use modeltopology,   only: NSD
use iomodule,        only: ifein,ifeout,ascout,append, mout,tout
use postProcessData, only: deffac,deformed
use monitorVars,     only: elmon, npmon
use plotControl,     only: Sscale, Ssb, &
                           Escale, Esb, &
                           Rscale, Rsb, &
                           Dscale, Dsb, &
                           Fscale, Fsb, &
                           grnum, &
                           hatsize, hatangle, iarcont, icompens, &
                           idisp, idtfl, iEffD, iRate, npnum, ratlog, &
                           sbhor, sysI, sysC, &
                           viscli, viscnu, vislog, &
                           Ygref, &
                           tempmn,tempmx, fluxmn,fluxmx,fluxmn,fluxmx, &
                           crdmin, crdmax

implicit none
!-pass
integer lu



!-locl
#include "lenchar.i"
character(len=1)      lf
character(len=LENCHAR)text
logical lfein,lfeout,arcont
double precision :: fac
integer          :: i

integer, external :: lnblk
! needed for the extrema
integer :: iDOF, iVertex, iElement
double precision :: effectiveSum

!-init
lf = CHAR(10)

lfein  = (ifein .ne.0)
lfeout = (ifeout.ne.0)
arcont= (iarcont.ne.0)

 write(lu,10) lfein,lfeout,ascout,append,NSTEP
10   format(1x,/'INPUT/OUTPUT:'/1x, &
 'Finite element input  file has been read = ',L1/1x, &
 'Finite element output file has been read = ',L1/1x, &
 'Ascii output files                       = ',L1/1x, &
 'Append output data to file               = ',L1/1x, &
 'Time                               itime = ',I8 )

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
fac = 0.d0

if (deformed) then
    fac = deffac
endif

write(lu,40) Sscale,Ssb,Escale,Esb,Rscale,Rsb,Dscale,Dsb,Fscale, &
 Fsb,fac
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

!if (Mout) then
!    write(lu,60) strsmn(2),strsmx(2),strsmn(1),strsmx(1), &
! ratemn(2),ratemx(2),ratemn(1),ratemx(1), &
! dispmn(2),dispmx(2),dispmn(1),dispmx(1)
!60   format(1x, &
! 'Stress:       ',2(1Pg9.1,' <= S <= ',1Pg9.1,2x)/1x, &
! 'Strain rate:  ',2(1Pg9.1,' <= R <= ',1Pg9.1,2x)/1x, &
! 'Displacement: ',2(1Pg9.1,' <= D <= ',1Pg9.1,2x))
!endif

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
    write(lu,70) tempmn,tempmx, &
 fluxmn(2),fluxmx(2),fluxmn(1),fluxmx(1)
70   format(1x, &
 'Temperature:  ',1Pg9.1,' <= T <= ',1Pg9.1/1x, &
 'Heat flux:    ',2(1Pg9.1,' <= F <= ',1Pg9.1,2x))
endif

if (Mout.or.Tout) then
    write(lu,'(1X)')
endif

if (irate.eq.0) then
    text='total'
else if (irate.eq.1) then
    text='viscous'
else if (irate.eq.2) then
    text='plastic'
endif

write(lu,80) vislog,viscli,ratlog,text(1:lnblk(text))
80   format(1x,'TRANSFORMATIONS AND SETTINGS:'/1x, &
 'Logarithmic viscosities    = ',L1/1x, &
 'Linearized  viscosities    = ',L1/1x, &
 'Logarithmic strain rates   = ',L1/1x, &
 'Strain rate                = ',7A)

if (iEffD.eq.0) then
    text = 'total'
elseif (iEffD.eq.1) then
    text = 'normal'
elseif (iEffD.eq.2) then
    text = 'slippery node'
elseif (iEffD.eq.3) then
    text = 'split node'
else
    text = 'BULLSHIT'
endif

write(lu,90) text(1:lnblk(text))
90   format(1x, &
 'Effective displacement type= ',13A)
if (icompens.eq.0) then
    text = 'gravity: raw (erosion surface)'
else if (icompens.eq.1) then
    write(text,91) Ygref
91       format('gravity: free air on Y = ',1PG12.1)
else if (icompens.eq.2) then
    write(text,92) Ygref
92       format('gravity: Bouguer on Y = ',1PG12.1)
else
    text = 'gravity: BULLSHIT'
endif

write(lu,100) text(1:lnblk(text))
  100   format(1x,80a)
if (iDisp.eq.1) then
    text = 'normal'
elseif (iDisp.eq.2) then
    text = 'slippery node'
elseif (iDisp.eq.3) then
    text = 'faulted node'
else
    text = 'BULLSHIT'
endif

write(lu,110) text(1:lnblk(text)),lf
  110   format(1x, &
 'displacements type         = ',14A)
if (iDtfl.eq.0) then
    text='stability ratio'
else if (iDtfl.eq.1) then
    text='oscillation limit ratio'
else
    text='relative modal accuracy'
endif

write(lu,120) DELTP,ALFAP,text(1:lnblk(text)),lf
  120   format(1x,'TIME STEP:'/1x, &
 'Time step size = ',1PG12.4/1x, &
 'Alpha          = ',1PG12.4/1x, &
 '"dt" command checks ',23A)

write(lu,130) hatsize,hatangle,arcont,sysi,sysc,sbhor
  130   format(1x,'LAYOUT:'/1x, &
 'Arrow: hatsize       = ',1PG12.4/1x, &
 '       top angle     = ',1PG12.4,' degrees'/1x, &
 '       contoured     = ',L1/1x, &
 'Symbol size          = ',1PG12.4/1x, &
 'Scaling of bc symbol = ',L1/1x, &
 'Horizontal scale bar = ',L1/)

!write(lu,140) ielbas,NUMEL,ielstp,inpbas,NUMNP,inpstp
!  140   format(1x,'INCREMENTING:'/1x, &
! 'element number ElNo:'/1x, &
! '    for (ElNo=',I8,'; ElNo<=',I8,'; ElNo=ElNo+',I8,')'/1x, &
! 'nodal point number NpNo:'/1x, &
! '    for (NpNo=',I8,'; NpNo<=',I8,'; NpNo=NpNo+',I8,')'/)

return
end
!-------------------------------------------------------------------------------
double precision function RealTime (itime)

USE CONSTANTS,      only: YEAR
use timestepmodule, only: timestepdatactx, &
                          nintg

implicit none
!-pass
integer itime
!-local
integer          :: n, igrp, i
double precision :: x


!-init
n    = 0
igrp = 1
i    = 1
x    = 0.d0
RealTime=0.d0

if (NINTG.le.0) then
    RealTime=-1.D0
    goto 300
endif

100    if (igrp.gt.NINTG) then
    RealTime=-2.D0
    goto 300
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
use modeltopology,  only: nstr
use label,          only: labels

implicit none
!-pass
integer          ::NUMSTR
integer          :: ielstr, isside, istr
double precision :: strs
dimension        :: IELSTR(*),ISSIDE(*),ISTR(2,*), &
                    STRS(NSTR,*)
integer          :: lu
!-locl
logical          :: updateStress
double precision :: SUM(6)
integer          :: i, iel, jel, iside, jside, n, nn, nprt
double precision :: time
double precision, external :: RealTime
!-init
nprt=-1

if (NUMSTR.le.0) then
    return
endif

do n=1,NUMSTR
    if (nprt.ge.51 .or. nprt.lt.0) then
        if (nprt.gt.0) then
            write(lu,2)
        endif
        write(lu,1) (LABELS(i),i=1,NSTR)
        write(lu,2)
        nprt=0
    endif

    iel = IELSTR(n)
    iside = ISSIDE(n)

    call clear (SUM,6,"SUM")
    call UPDATE (SUM,STRS(1,n),NSTR)

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
            write(lu,3) iel,iside,time,(SUM(i),i=1,NSTR)
        endif

        nprt=nprt+1

        if (n.lt.NUMSTR) then
            do nn=n+1,NUMSTR

                jel = IELSTR(nn)
                jside = ISSIDE(nn)

                if (iel.eq.jel .and. iside.eq.jside) then

                    call UPDATE (SUM,STRS(1,nn),NSTR)
                    time = RealTime (ISTR(2,nn))

                    if (time.lt.0.0) then
                        write(lu,4) iel,iside
                    else
                        write(lu,3) iel,iside,time,(SUM(i),i=1,NSTR)
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
 3    format(1x,I8,2X,I1,4X,1PG10.2,6(2X,1PG10.2))
 4    format(1x,I8,2x,I1,4X,'NO STRESS')
! ELEM# SIDE  AFTER TIME      SXX         SYY          SXY
! IIIII  8    FFFFFFFFFF  XXXXXXXXXX  XXXXXXXXXX  XXXXXXXXXXX
 end
!-------------------------------------------------------------------------------
subroutine shwink (IWELM,IWSIDE,IWTIME, &
                   NUMWNK,lu)

use modeltopology,  only: ndof
use timestepmodule, only: timestepdatactx, nintg

implicit none
!-pass
 integer NUMWNK
integer :: iwelm, iwside, iwtime
 dimension IWELM(*),IWSIDE(*),IWTIME(*)
integer :: lu

!-locl
integer :: i, n, iel, m, iside
integer :: itime
double precision :: time
double precision, external :: RealTime


if (NUMWNK.le.0) then
    return
endif

n = 1
!    find an element that has not been processed
100    if (n.gt.NUMWNK) then
    goto 400
endif

if (IWELM(n).gt.0) then
    goto 200
endif

n = n + 1
goto 100

200    iel = IWELM(n)
iside = IWSIDE(n)
m = n

300    if (m.le.NUMWNK) then
    if (IWELM(m).eq.iel .and. IWSIDE(m).eq.iside) then

        if (IWTIME(m).gt.0) then
            itime = IWTIME(m)-1
        else
            itime = -IWTIME(m)
        endif

        if (itime.gt.0) then
            time = RealTime (itime)
        else
            time = 0.d0
        endif

        if (time.lt.0) then
            if (time.lt.-1.1) then
                if (IWTIME(m).gt.0) then
                    write(lu,1) iel,iside
 1                       format(1x,'Winkler pressure on element ',I8, &
                     ' (side ',I1,') is never applied')
                else
                    write(lu,2) iel,iside
 2                       format(1x,'Winkler pressure on element ',I8, &
                     ' (side ',I1,') is always applied')
                endif
            endif
        else
            if (IWTIME(m).gt.0) then
                write(lu,3) iel,iside,time
 3                   format(1x,'Winkler pressure on element ',I8, &
                 ' (side ',I1,') applied at ',1PG14.6,' year')
            else
                write(lu,4) iel,iside,time
 4                   format(1x,'Winkler pressure on element ',I8, &
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

400 continue
do i=1,NUMWNK
     IWELM(i) = IABS(IWELM(i))
enddo

return
end
!-------------------------------------------------------------------------------
subroutine shanom (ITANOM,lu)

use modeldefinition, only: ntanom
use timestepmodule, only: timestepdatactx, &
                          nintg

 implicit none
!-pass
integer itanom
dimension ITANOM(3,*)
integer :: lu

!-locl
integer :: mode, node, itime
integer :: i
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
    if (time.lt.0.) then
        if (mode.eq.0) then
            write(lu,1) node
 1               format(1x,'Temperature of node ',I8, &
             ' is never replaced')
        else
            write(lu,2) node
 2               format(1x,'Temperature of node ',I8, &
             ' is never increased')
        endif
    else
        if (mode.eq.0) then
            write(lu,3) node,time
 3               format(1x,'Temperature of node ',I8, &
             ' is replaced at ',1PG14.6,' year')
        else
            write(lu,4) node,time
 4               format(1x,'Temperature of node ',I8, &
             ' is increased at ',1PG14.6,' year')
        endif
    endif
enddo

return
end
!-------------------------------------------------------------------------------

subroutine showTimestepGroups (lu)

USE CONSTANTS,       only: year
use timestepmodule,  only: timestepdatactx, nintg

implicit none
!-pass
integer          :: lu

!-local
integer          :: igrp
double precision :: t0, t1, dt

logical          :: doMechanic, doThermal



if (NINTG.le.0) then
    return
endif

t0=0.d0

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
