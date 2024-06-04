subroutine RDTHRM (T0compute,ierr)
! This subroutine reads the thermal input,
! containing initial values and boundary condition.

USE MATERIALSMODULE, only: allocatematdatatemp, &
                           matdatactx
USE MESHDATAMODULE,  only: meshdatactx
USE MODELDATAMODULE, only: modeldatactx, &
                           allocatethermalmodeldata, &
                           allocatethermalbcdata, &
                           allocatethermalanomaly, &
                           allocatethermalfluxbcs, &
                           allocatethermalwinkler, &
                           itdiag
USE MODELDEFINITION, only: iadvec, IDIFT, &
                           itect, itmode, &
                           itpres, mode, nflx, nflxglobal, &
                           ntanom, ntanomglobal, &
                           nteq, nteqglobal, nteqlocal, &
                           ntwink, ntwinkglobal, &
                           numat, numel, numfn, numnp, numslp, &
                           allTemperaturesFixed
USE MODELTOPOLOGY,   only: nen, nsd, ndof
USE TIMESTEPMODULE,  only: outputcontroldatactx, &
                           timestepdatactx, &
                           allocateoutputcontroldata, &
                           allocatetimestepdata, &
                           rdtime, &
                           putout, &
                           chkout, &
                           rdtstp, &
                           NINTG, &
                           NTCALC, &
                           MAXTIT, &
                           NTPRT
USE MODELCTX
USE ALGEBRA,         only: CLEAR
use filemodule,      only: FullInputThermal, &
                           workpath, &
                           MyTmpDir
use debugmodule,     only: iecho, &
                           parallelLog, &
                           xit, &
                           debug
use petscksp
use iomodule          

!-------------------------------------------------------------------------------
! Routine to read thermal input data.
! Secondary routines residing in other libraries are:
!    LOCAL     - "tecin.f"
!    WRITST    - "tecin.f"
!    FGEN      - "tecin.f"
!    IGEN      - "tecin.f"
!-------------------------------------------------------------------------------

implicit none

!-pass
logical :: T0compute
integer :: ierr

#ifdef SPARSE
#include "petsc/finclude/petscsys.h"
PetscErrorCode petscierr
#endif

!-locl
character(len=255) :: tmpnm
character          :: quote
integer, external  :: soelim,ifopen, lnblk, csystem
#ifdef SPARSE
logical            :: there
#endif
integer            :: ERROR

integer            :: i, j, nn
integer            :: ios

integer            :: ndoftemp
!-init

if (iecho.eq.8) then
    call parallelLog("rdthrm","just entered")
endif

numnp = meshdatactx%nvglobal

quote = CHAR(34)

if (MyTmpDir(1:1).eq." ") then
    ! tmpdir not set, use default
     FullInputThermal   = tmpnm(trim(WORKPATH)//"/FET") ! FET + date/time stamp
else
    ! tmpdir set, use it.
    FullInputThermal   = tmpnm(trim(MyTmpDir)//"/FET")
endif

hasthermal = .false.

if (MODE.lt.3) then
    return
endif

hasthermal = .true.

if (iecho.eq.8) then
    call parallelLog("rdthrm","possible soelim")
endif

#ifdef SPARSE
if (iecho.eq.8) then
    call parallelLog("rdthrm","checking existence of "//FullInputThermal)
endif

inquire(file=FullInputThermal,exist=there)

if (.not.there) then

!    write(*,*) "firing system call: " // &
!               "soelim -r -I "//trim(WORKPATH)//" "//trim(TEMPIN)//" >| "//trim(FullInputThermal)

    ! include workpath to find possible files that are to be included, there.
    ierr = csystem("soelim -r -I "//trim(WORKPATH)//" "//trim(TEMPIN)//" >| "//trim(FullInputThermal))

    if (ierr.ne.0) then
        write(stderr,4) ierr
        return
    endif
else
    ! it already exists, do nothing
endif



#else

! for plnplt, and such
!write(*,*) "workpath : ", trim(WORKPATH)
!write(*,*) "tempin   : ", trim(TEMPIN)
!write(*,*) "fullinput: ", trim(FullInputThermal)

ierr = csystem("soelim -r -I "//trim(WORKPATH)//" "//trim(TEMPIN)//" >| "//trim(FullInputThermal))


if (ierr.ne.0) then
     write(stderr,4) ierr
     write(stderr,*) "ERROR: failed to execute soelim command:"
     write(stderr,*) "soelim -r -I "//trim(WORKPATH)//" "//trim(TEMPIN)//" >| "//trim(FullInputThermal)
     write(stderr,*) "returned error code:  ", ierr
     STOP "Leaving GTECTON..."
endif

#endif

if (iecho.eq.8) then
    call parallelLog("rdthrm","finished possible soelim")
endif

4  format(1x,'soelim: exit ', I1)
ierr = ifopen(luin,FullInputThermal,'old')

if (ierr.ne.0) then
    return
endif

if (iecho.eq.4) then
    write(stderr,5)
 5   format(1x,79('*')/1x,'*** THERMAL PROBLEM DEFINITION ***'/ &
 1x,'read (a80) TITLE')
endif

read(luin,'(A)',err=100,end=110,iostat=ios) TTITLE

if (iecho.eq.4) then
    write(stderr,10)
10    format(1x,'read (6i5) ITMODE,NTCALC,NTPRT,MAXTIT,IADVEC,IDIFT')
endif

read(luin,15,err=100,end=110,iostat=ios) &
   ITMODE,NTCALC,NTPRT,MAXTIT,IADVEC,IDIFT
15    format(15I5)

if (iecho.eq.8) then
    call parallelLog("rdthrm","read title and options")
endif


call allocateoutputcontroldata(outputcontroldatactx,1)
call allocatetimestepdata(timestepdatactx,1)


if (iecho.eq.4) then
    write(stderr,20)
20   format(1x,'*** read ICALC data ***'/1x, &
 'nread = MAX(NTCALC,1)'/1x, &
 'read (51i5) (ICALC(i), i=1,nread)')
endif

if (iecho.eq.8) then
    call parallelLog("rdthrm","reading time")
endif

call RDTIME(timestepdatactx%ICALC, NTCALC)

if (iecho.eq.8) then
    call parallelLog("rdthrm","finished reading time")
endif


if (iecho.eq.4) then
    write(stderr,21)
21   format(1x,'*** read ITPRINT data ***'/1x, &
 'nread = MAX(NTPRT,1)'/1x, &
 'read (51i5) (ITPRINT(i), i=1,nread)')
endif

call RDTIME(outputcontroldatactx%ITPRINT, NTPRT)

if (OUTFIL(1)) then
    write(luout(1),23)
23       format(1x/1x/1x/1x,79('*')/1x/1x,'H E A T',3x, &
     'E Q U A T I O N S',3x,'I N P U T'/1x/1x)
    call writst(luout(1),TITLE)
endif

if (OUTFIL(1)) then
    write(luout(1),25) ITMODE
25   format(///1x, &
 'T E M P E R A T U R E   C O N T R O L   D A T A      ',  //5X, &
 'INITIAL TEMPERATURE MODE . . . . . . . . . (ITMODE) =',I5/ 5X, &
 '  EQ.0,  CALCULATE STEADY STATE DISTRIBUTION         ',  / 5X, &
 '  EQ.1,  READ INITIAL DISTRIBUTION FROM FILE         ',  / 5X, &
 '  EQ.2,  CALCULATE INITIAL DISTRIUBTION              ',  / 5X, &
 '  EQ.3,  READ PARAMETERS AND CALCULATE               ',  /   )
endif

if (ITMODE.eq.2) then
    call PRFUNC (0)
endif

if (OUTFIL(1)) then
    write(luout(1),26) IADVEC,IDIFT,MAXTIT
26   format(//5X, &
 'ADVECTING GRID SWITCH  . . . . . . . . . . (IADVEC) =',I5/ 5X, &
 '  EQ.0,  NORMAL CALCULATION                          ',  / 5X, &
 '  EQ.1,  ADVECTING GRID CALCULATION                  ',  //5X, &
 'DIFFERENTIAL TEMPERATURES SWITCH  . . . . . (IDIFT) =',I5/ 5X, &
 '  EQ.0,  NO DIFFERENTIAL TEMPERATURES                ',  / 5X, &
 '  EQ.1,  DIFFERENTIAL TEMPERATURES                   ',  //5X, &
 'NO. OF STEPS BETWEEN STIFFNESS REFORMS . . (MAXTIT) =',I5/   )
endif

if (IADVEC.ne.0) then
    if (MODE.ne.3 .and. MODE.ne.4 .and. MODE.ne.5) then
        write(stderr,27)
27           format(1x,'Advecting grid calculation only for',1x, &
         'MODE=3, 4 or 5')
        call xit(1," ")
    endif
endif

if (OUTFIL(1)) then
    write(luout(1),30) NTCALC
30   format(5X, &
 'THE NUMBER OF TEMPERATURE CALCS IS . . . . (NTCALC) =',I5/ 5X, &
 '  GT.0, CALCULATION AT SELECTED TIME STEPS           ',  / 5X, &
 '  EQ.0, NO THERMAL CALCULATION                       ',  / 5X, &
 '  LT.0, CALCULATION EVERY ABS(NTCALC) STEP           '       )
endif

if (OUTFIL(1).and.NTCALC.GT.0) then
    write(luout(1),31)
31    format(5x,'TIME STEPS SELECTED FOR CALCULATION:'/)
   call PUTOUT (timestepdatactx%ICALC, &
           NTCALC)
endif

if (OUTFIL(1)) then
    WRITE(luout(1),32) NTPRT
32   FORMAT(//5X, &
 'THE NUMBER OF TEMPERATURE OUTPUTS IS . . .  (NTPRT) =',I5/ 5X, &
 '  GT.0, OUTPUT AT SELECTED TIME STEPS                ',  / 5X, &
 '  EQ.0, NO OUTPUTS                                   ',  / 5X, &
 '  LT.0, OUTPUTS AT EVERY ABS(NTPRT) STEP             '       )
endif

if (OUTFIL(1).and.NTPRT.GT.0) then
    write(luout(1),33)
33   format(5x, &
 'SELECTED OUTPUTS OCCUR AT THE FOLLOWING TIME STEPS:'/)
endif

call PUTOUT (outputcontroldatactx%ITPRINT, &
             NTPRT)

!    check if requested output times are actually calculated
call CHKOUT(timestepdatactx, outputcontroldatactx)

if (iecho.eq.8) then
    call parallelLog("rdthrm","allocating thermal modeldata")
endif

call allocatethermalmodeldata(modeldatactx, numnp, numel)
call allocatethermalbcdata(modeldatactx, numnp)

if (iecho.eq.8) then
    call parallelLog("rdthrm","allocated thermal modeldata")
endif

if (NINTG.gt.0) then
    ALLOCATE(modeldatactx%TLST(meshdatactx%Nvlocal), STAT=ERROR)
endif

!    read initial distribution if requested
call CLEAR (modeldatactx%T,NUMNP,"modeldatactx%T")

call RDTEMP (modeldatactx%T,meshdatactx%Nvglobal,ITMODE)

!     calculate initial temperatures if requested
if (ITMODE.eq.2 .or. (ITMODE.eq.3.and.T0compute)) then
    call NITTEM (modeldatactx%T,meshdatactx%X,NSD,NUMNP)
endif

if (iecho.eq.8) then
    call parallelLog("rdthrm","allocating thermal boundary condition data")
endif

call allocatethermalbcdata(modeldatactx, numnp)


if (iecho.eq.8) then
    call parallelLog("rdthrm","reading thermal bc")
endif



CALL RDTBC(meshdatactx%nvglobal,ITMODE)

if (iecho.eq.8) then
    call parallelLog("rdthrm","finished reading thermal bc")
endif

!    allocate space for differential temperature arrays

   hasdiffthermal = ((NUMSLP+NUMFN).gt.0.and.IDIFT.eq.1)


!    Eliminate degrees of freedom associated with boundary conditions
!       and compute NTEQ. Set equation numbers.
ALLOCATE(modeldatactx%IDTglobal(meshdatactx%Nvglobal),STAT=ERROR)
if (ERROR.ne.0) then
    write(*,*) 'could not allocate modeldata IDTglobal. Error ', ERROR
    call xit(1," ")
endif

ALLOCATE(modeldatactx%IDTXglobal(meshdatactx%Nvglobal),STAT=ERROR)
if (ERROR.ne.0) then
    write(*,*) 'could not allocate modeldata IDTXglobal. Error ', ERROR
    call xit(1," ")
endif

ALLOCATE(modeldatactx%IDTX(meshdatactx%Nvlocal),STAT=ERROR)
if (ERROR.ne.0) then
    write(*,*) 'could not allocate modeldata IDTX. Error ', ERROR
    call xit(1," ")
endif

modeldatactx%IDTglobal = 0
modeldatactx%IDTXglobal = 0
modeldatactx%IDT = 0
modeldatactx%IDTX =  0

call SETTID (modeldatactx%IDTglobal, &
             modeldatactx%IDTXglobal, &
             modeldatactx%ITBCGLOBAL)

!    establish implicit params array
call RDTSTP(timestepdatactx%TALF,NINTG)

!    allocate space for thermal properties
call allocatematdatatemp(NUMEL)


!    read thermal material properties

if (iecho.eq.8) then
    call parallelLog("rdthrm","reading thermal material properties")
endif

matdatactx%PRPTEM = 0d0
call RDTMAT(NSD,NUMAT)

if (iecho.eq.8) then
    call parallelLog("rdthrm","finished reading thermal material properties")
endif

!    if (LARGELM) then
!        if (iecho.eq.4) write(stderr,40) '(5I8)'
!   40        format(1x,'read NFLX,NTANOM,ITPRES,ITECT,NTWINK ',A5)
!        read(luin,'(5i8)',err=100,end=110,iostat=ios) NFLX,NTANOM, &
!               ITPRES,ITECT,NTWINK
!    else
!        if (iecho.eq.4) write(stderr,40) '(5I5)'
!        read(luin,'(5I5)',err=100,end=110,iostat=ios) NFLX,NTANOM, &
!               ITPRES,ITECT,NTWINK
!    endif

if (iecho.eq.4) then 
   write(stderr,40) '(5I12)'
40        format(1x,'read NFLX,NTANOM,ITPRES,ITECT,NTWINK ',A5)
endif


if (debug) then
    write(*,*) "reading NFLX,NTANOM,ITPRES,ITECT,NTWINK; 5i12"
endif


read(luin,'(5i12)',err=100,end=110,iostat=ios) NFLX,NTANOM,ITPRES,ITECT,NTWINK

if (debug) then
    write(*,*) "read ", NFLX,NTANOM,ITPRES,ITECT,NTWINK
endif

!
if (OUTFIL(1)) then
    write(luout(1),50) ITPRES,ITECT,NFLX,NTWINK,NTANOM

50   format(///1x, &
 'T E M P E R A T U R E   S W I T C H E S',//5X, &
 'THERMAL PRESSURE SWITCH  . . . . . . . . . (ITPRES) =',I8/ 5X, &
 '  EQ.0, ADD THERMAL CONTRIBUTIONS TO LOADS           ',  / 5X, &
 '  EQ.1, DO NOT ADD THERMAL CONTRIBUTIONS             ',  //5X, &
 'TECTONIC CONTRIBUTION FROM THERMAL ANOMALY  (ITECT) =',I8/ 5X, &
 '  HAS MEANING ONLY WHEN ITPRES = 0                   ',  / 5X, &
 '  EQ.0, ADD THERMAL EXPANSION CONTRIBUTIONS          ',  / 5X, &
 '  EQ.1, NEGLECT THERMAL EXPANSION CONTRIBUTIONS      ',  //5X, &
 'NO. OF HEAT FLUX BOUNDARY CONDITIONS . . . . (NFLX) =',I8//5X, &
 'NO. OF THERMAL WINKLER FLUXES  . . . . . . (NTWINK) =',I8//5X, &
 'NO. OF THERMAL ANOMALY DATA  . . . . . . . (NTANOM) =',I8/ 5X, &
 '  GT.0, READ NTANOM ANOMALY DATA FROM INPUT          ',  / 5X, &
 '  LT.0, CALCULATE ANOMALY AFTER ITIME=IABS(NTANOM)   ',  /   )

endif

if (MODE.ne.6) then
    if (ITPRES.ne.0 .or. ITECT.ne.0) then
        write(stderr,52)
    endif
52  format('rdthrm: ITPRES!=0 or ITECT!=0 makes no sense')
endif

if (NTANOM.lt.0) then
    call PRFUNC (1)
endif

NTANOMglobal  = NTANOM
NFLXglobal    = NFLX
NTWINKglobal  = NTWINK

if (NTANOM.gt.0) then
!   allocate space for thermal anomaly arrays
    call allocatethermalanomaly(modeldatactx,ntanom)
    call RDTAN(modeldatactx%ITANOM,modeldatactx%TANOM)
else
    call CHKANM(timestepdatactx%ICALC,NTCALC,NTANOM)
endif


 
if (NFLX.gt.0) then
    call allocatethermalfluxbcs(modeldatactx,NFLX)
    modeldatactx%IFLX = 0
    modeldatactx%IFLS = 0
    modeldatactx%BFLX = 0d0
    call RDFLUX (modeldatactx%IFLS, &
                 modeldatactx%IFLX, &
                 modeldatactx%BFLX)
endif

if (iecho.eq.8) then
    call parallelLog("rdthrm","finished reading flux")
endif


if (NTWINK.gt.0) then

    if (iecho.eq.8) then
        call parallelLog("rdthrm","allocating thermal winkler data")
    endif

    call allocatethermalwinkler(modeldatactx, NTWINK)

    if (iecho.eq.8) then
        call parallelLog("rdthrm","allocated thermal winkler data")
    endif

    call RDWFLX (modeldatactx%ITWINK, modeldatactx%TWINK)

    if (iecho.eq.8) then
        call parallelLog("rdthrm","finished reading flux")
    endif

endif

#ifdef SPARSE
!   if (PetSc_Init) then
!   call MPI_Barrier(PETSC_COMM_WORLD,petscierr)
!   endif

if (iecho.eq.1) then
    write(stderr,60) quote,FullInputThermal(1:lnblk(FullInputThermal)),quote
    write(stderr,60) quote,TEMPIN(1:lnblk(TEMPIN)),quote
endif

#else

if (iecho.eq.1) then
    write(stderr,60) quote,FullInputThermal(1:lnblk(FullInputThermal)),quote
    write(stderr,60) quote,TEMPIN(1:lnblk(TEMPIN)),quote
endif

#endif
60    format(1x,'EOF ',80a)

!    fill LMT array and LMTX array
!       this setns only the LMT array
ALLOCATE(modeldatactx%LMT(NEN,meshdatactx%Nelocal), STAT=ERROR)
if (ERROR.ne.0) then
    write(*,*) 'could not allocate modeldata LMT. Error ', ERROR
    call xit(1," ")
endif
modeldatactx%LMT = 0d0


ALLOCATE(modeldatactx%LMTX(NEN,meshdatactx%Nelocal), STAT=ERROR)
if (ERROR.ne.0) then
    write(*,*) 'could not allocate modeldata LMTX. Error ', ERROR
    call xit(1," ")
endif
modeldatactx%LMTX = 0d0


do i=1,meshdatactx%Nelocal
    do j=1,NEN
        nn = meshdatactx%IEN(j,i)
        modeldatactx%lmt(j,i) = modeldatactx%IDTglobal(nn)
    enddo
enddo

call LOCLTX (modeldatactx%IDTXglobal,&
              meshdatactx%IEN,  &
             modeldatactx%LMTX, &
             modeldatactx%NSLIP, &
             modeldatactx%NFAULT, &
             modeldatactx%FAULT)

NTEQ = NTEQlocal

!allocate(itdiag(nteq))
! Lukas vd Wiel; 22 Okt 2021
! temporarily set to global, to pass the parallel thermal benchmarks
!write(*,*) "allocating itdiag NTEQglobal",NTEQglobal, allTemperaturesFixed
allocate(itdiag(NTEQglobal), STAT=ERROR)
if (ERROR.ne.0) then
    write(*,*) 'could not allocate itdiag. Error ', ERROR
    call xit(1," ")
endif

itdiag = 0

if (.not. allTemperaturesFixed) then
	call TCOLHT(modeldatactx%LMT,modeldatactx%LMTX,NEN)
endif

if (iecho.eq.8) then
    call parallelLog("rdthrm","finished reading thermal input")
endif

return

100    write(stderr,101) ios,quote,FullInputThermal(1:lnblk(FullInputThermal)),quote
  101    format(1x,'RDTHRM: fatal input error ',i8,' on ',80a)

call prerr('RDTHRM')
call xit(1," ")

110    write(stderr,111) quote,FullInputThermal(1:lnblk(FullInputThermal)),quote
  111    format(1x,'RDTHRM: premature EOF error on ',80a)

call xit(1," ")

end subroutine

!-----------------------------------------------------------------------
!        secondary routines follow in alphabetic order
!-----------------------------------------------------------------------
!    subroutine CHKOUT(ICALC,ITPRINT,NINTG,NTCALC,NTPRT)
!
!c Routine checking if requested output times match calculation times and
!c if times are in ascending order.
!
!c SNGL    implicit real (a-h,o-z)
!    implicit double precision (a-h,o-z)
!c-pass
!    integer NINTG,NTCALC,NTPRT
!    dimension ICALC(*),ITPRINT(*)
!c-comm
!
!    if (NINTG.le.0) return
!
!    if (NTCALC.gt.1) then
!c        check ascending order
!        do 100 i=2,NTCALC
!        if (ICALC(i-1).ge.ICALC(i)) then
!            write(stderr,10)
!   10            format(1x,'ICALC not ascending')
!            call xit(1," ")
!        endif
!100        continue
!    endif
!    if (NTPRT.gt.1) then
!c        check ascending order
!        do 200 i=2,NTPRT
!        if (ITPRINT(i-1).ge.ITPRINT(i)) then
!            write(stderr,20)
!   20            format(1x,'ITPRINT not in ascending order')
!            call xit(1," ")
!        endif
!200        continue
!    endif
!    if (NTCALC.eq.0 .and. NTPRT.ne.0) then
!            write(stderr,25)
!   25       format(1x,'Number of thermal calculations = 0, so',1x,
!     >       'NTPRT should be = 0 too')
!            call xit(1," ")
!        endif
!
!    if (NTCALC.gt.0 .and. NTPRT.gt.0) then
!c        check if all ITPRINT numbers occur in ICALC also
!        j = 1
!        do 400 i=1,NTPRT
!300        if (j.gt.NTCALC) then
!            write(stderr,30) i,ITPRINT(i)
!   30               format(1x,'ITPRINT(',I2,')=',I3,
!     >               ' does not occur in ICALC')
!            call xit(1," ")
!        endif
!        ISTEP = ITPRINT(i)
!        if (ISTEP.eq.ICALC(j)) goto 400
!        j = j + 1
!        goto 300
!400        continue
!    else if (NTCALC.ge.0 .and. NTPRT.lt.0) then
!        MTPRT = INT(REAL(ICALC(NTCALC))/REAL(-NTPRT))
!        if (MTPRT.lt.1) then
!                write(stderr,35) -NTPRT
!   35           format(1x,'No temperature calculation every ',I4,1x,
!     >           'time steps')
!                call xit(1," ")
!            endif
!        j = 1
!        do 600 i=1,MTPRT
!500        if (j.gt.NTCALC) then
!            write(stderr,40) i
!   40               format(1x,'ITPRINT(',I2,')',
!     >               ' does not occur in ICALC')
!            call xit(1," ")
!        endif
!        ISTEP = i * (-NTPRT)
!        if (ISTEP.eq.ICALC(j)) goto 600
!        j = j + 1
!        goto 500
!600        continue
!    else if (NTCALC.lt.0 .and. NTPRT.gt.0) then
!        MTCALC = INT(REAL(ITPRINT(NTPRT))/REAL(-NTCALC))
!        j = 1
!        do 800 i=1,NTPRT
!700        if (j.gt.MTCALC) then
!            write(stderr,30) i,ITPRINT(i)
!            call xit(1," ")
!        endif
!        ISTEP = j * (-NTCALC)
!        if (ISTEP.eq.ITPRINT(i)) goto 800
!        j = j + 1
!        goto 700
!800        continue
!    else if (NTCALC.lt.0 .and. NTPRT.lt.0) then
!        if (MOD(-NTPRT,-NTCALC).ne.0) then
!        write(stderr,50)
!   50        format(1x,'NTPRT and NTCALC do not match')
!        call xit(1," ")
!        endif
!    endif
!
!    return
!    end
!-----------------------------------------------------------------------
subroutine CHKANM(ICALC,NTCALC,NTANOM)

use timestepmodule, only: nintg
use debugmodule,     only: xit
use iomodule          

! Routine checking if requested anomaly addition time occurs in ICALC

implicit none
!-pass
integer   :: icalc,NTCALC,NTANOM
dimension :: ICALC(*)
!-local
integer   :: i

!-init
i=1

if (NINTG.le.0 .or. NTANOM.ge.0) then
    return
endif

!    check if all -NTANOM occurs in ICALC
if (NTCALC.eq.0) then
    write(stderr,1)
 1        format(1x,'NTCALC=0: no thermal anomaly added')
    call xit(1," ")
else if (NTCALC.gt.0) then

100  if (i.gt.NTCALC) then
        write(stderr,2) -NTANOM
 2           format(1x,'ITIME = ',I6,' doesn''t occur in ICALC:'/1x, &
         'Thermal anomaly not added')
        call xit(1," ")
    endif

    if (ICALC(i).eq.-NTANOM) then
        goto 1000
    endif

    i = i + 1
    goto 100
else
    if (MOD(-NTANOM,-NTCALC).ne.0) then
        write(stderr,3) -NTCALC,-NTANOM
 3           format(1x,'Thermal calculation every ',I4,' time steps'/ &
         1x,'Anomaly added at step     ',I4/ &
         1x,'Obviously, no thermal anomaly is added at all')
!        call xit(1," ")
    endif
endif

1000    return
end subroutine

!-----------------------------------------------------------------------

subroutine LOCLTX (IDTXglobal,IEN,LMTX,NSLIP,NFAULT, &
                    FAULT)

#ifdef SPARSE
USE MESHDATAMODULE,  only: meshdatactx
#endif
!USE MODELDATAMODULE
USE MODELDEFINITION, only: numfn, numslp, idift, numel
USE ALGEBRA,         only: clear
use modelctx,        only: getsize
use modeltopology,   only: ndof, nen
use debugmodule,     only: xit
use iomodule          
!
!    Localize IDTX and transfer sign information to LMTX
!
implicit none
!-pass

integer          :: idtxglobal, ien, nslip, nfault
double precision :: fault
dimension        :: IDTXglobal(*),IEN(NEN,*),NFAULT(3,*), &
                    FAULT(NDOF,*), NSLIP(5,*)
integer          :: lmtx(nen, numel)

!-locl
integer          :: i,j,jj,nel,node, neloc, referenceJ
double precision :: fltmax
!
if (NUMSLP+NUMFN.eq.0.or.IDIFT.eq.0) then
    return
endif
!
call CLEAR (LMTX,NEN*NUMEL,"LMTX")
!
if (NUMSLP.gt.0) then
    do i=1,NUMSLP
        nel  = NSLIP(1,i)
        node = NSLIP(2,i)

#ifdef SPARSE
        if (getsize().eq.1) then
            neloc = nel
        else
            neloc = meshdatactx%Glo2LocElement(nel)
        endif
#else
        neloc = nel
#endif


        do jj=1,NEN
            j = jj
            if (node.eq.IEN(j,neloc)) then
                referenceJ = j
                goto 100
            endif
        enddo

        write(stderr,1) i,node,nel
 1          format(1x,'LOCLTX: NSLIP(',I3,'): node=',I12, &
         ' does not occur in element ',I12)
        call xit(1," ")

100     LMTX(referenceJ,neloc) = ISIGN(IDTXglobal(node),NSLIP(3,i))
    enddo
endif
!
if (NUMFN.gt.0) then
   do i=1,NUMFN
      nel  = NFAULT(1,i)
      node = NFAULT(2,i)

      do j=1,NEN
          if (node.EQ.IEN(j,nel)) goto 200
      enddo

      write(stderr,2) i,node,nel
 2         format(1x,'LOCLTX: NFAULT(',I3,'): node=',I12, &
         ' does not occur in element ',I12)
      call xit(1," ")

200   fltmax = FAULT(1,i)

      do j=2,NDOF
          if (FAULT(j,i).gt.fltmax .or. FAULT(j,i).lt.-fltmax) then
              fltmax = FAULT(j,i)
          endif
      enddo

      if (fltmax.gt.0.) then
          LMTX(j,nel) =  IDTXglobal(node)
      else if (fltmax.lt.0.) then
          LMTX(j,nel) = -IDTXglobal(node)
      else
          write(stderr,3) node,nel,i
 3             format(1x,'LOCLTX: no slip on faulted node ',I12, &
          ' of element ',I5,' (entry ',I12,')')
          call xit(1," ")
      endif

   enddo
endif
!
return
end subroutine

!-----------------------------------------------------------------------

subroutine RDFLUX (IFLS,IFLX,BFLX)

USE MODELDEFINITION, only: NFLX, NFLXglobal
use formatsmodule,   only: FMT_tempin_flux
use debugmodule,     only: iecho
use iomodule          

implicit none


integer, parameter :: MAXP=50

!-pass
integer            :: iflx, ifls
double precision   :: bflx
dimension          :: IFLX(*),IFLS(*),BFLX(*)
!-local
integer            :: n, nn


!-init
nn = MAXP + 1

if (iecho.eq.4) then
    write(stderr,1) FMT_tempin_flux
 1   format(1x,'if (NFLX > 0){'/ &
 5X,'*** read element normal heat flux data ***'/ &
 5X,'read ',A14,' IFLX(n),IFLS(n),BFLX(n)'/ &
 5X,'if (IFLX(n)==0) return'/1x,'}')
endif

if (NFLXglobal.le.0) then
    return
endif

call FLXGEN(IFLX,IFLS,BFLX,NFLXglobal,NFLX,1,FMT_tempin_flux,"RDFLUX")

if (.not.OUTFIL(1)) then
    return
endif

do n=1,NFLX
    if (nn.gt.MAXP) then
        nn = 0
        write(luout(1),2)
 2           format(///1x, &
         'N O R M A L    S U R F A C E    H E A T    F L U X' &
         //5X,'ELEM #',3X,'SIDE',11X,'FLUX'/)
    endif

    write(luout(1),3) IFLX(n),IFLS(n),BFLX(n)
 3        format(3X,I8,4X,I1,2X,5X,1PG14.6)
    nn = nn + 1
enddo

return
end subroutine

!-------------------------------------------------------------------------------

subroutine RDWFLX (ITWINK,TWINK)

USE MODELDATAMODULE
USE MESHDATAMODULE,  only: meshdatactx
USE MODELCTX,        only: getrank
USE MODELDEFINITION
USE LABELMODULE
use debugmodule,     only: iecho, xit
use iomodule          

implicit none
!-pass

integer            :: ITWINK
double precision   :: TWINK
dimension          :: ITWINK(3,*),TWINK(*)
!-locl
logical            :: end
character(len=48)  :: record
character(len=16)  :: fmt
integer            :: nreadglobal
integer            :: lstflx, numflx
integer            :: i, it, is, ios, nel, n, nn
double precision   :: f

integer,external   :: lnblk
!      integer getrank, nreadglobal

!-init
numflx = 0
lstflx = 0
nreadglobal = 0
!
if (NTWINK.le.0) then
    return
endif

fmt = '(2I12,2I5,G14.0)'
!
if (iecho.eq.4) then
    write(stderr,1) fmt
 1   format(1x,'*** read thermal Winkler bc data ***'/1x, &
 'elm=1; while (elm != 0) {'/1x, &
 '    read ',A15,' elm,increment,iside,itime,Winkler'/1x, &
 '}')
endif

100    call readrc (luin,record,'RDWFLX',end)

if (end) then
    goto 200
endif

read (record,2,err=300,end=300,iostat=ios) nel,is, &
                                           it,f
 2    format(4I5,G14.0)

if (nel.eq.0) then
    goto 200
endif

nreadglobal = nreadglobal + 1

if (meshdatactx%e2p(nel).ne.getrank()) then
    goto 100
endif


numflx = numflx + 1
NTWINK = numflx

if (numflx.gt.NTWINKglobal) then
    write(stderr,6) NTWINKglobal
 6     format(1x,'RDWFLX: number of flux bc''s exceeds ',I8)
    call xit(1," ")
endif

TWINK(numflx) = f
ITWINK(1,numflx) = nel
ITWINK(2,numflx) = is
ITWINK(3,numflx) = it

lstflx = numflx

goto 100

200 if (nreadglobal.ne.NTWINK) then
    write(stderr,7) nreadglobal,NTWINK
 7       format(1x,'RDWFLX: read ',I8,' entries, expected ',I8)
    call xit(1," ")
endif

if (.not.OUTFIL(1)) then
    return
endif

nn = 0
do n=1,NTWINK
    if (MOD(nn,50).eq.0) then
        write(luout(1),8)
 8           format(///1x, &
         'W I N K L E R   S U R F A C E   F L U X E S' &
         //5X,'ELEM #',3X,'SIDE',3X,'TIME',11X,'FLUX'/)
    endif
    write(luout(1),9) (ITWINK(i,n),i=1,3),TWINK(n)
 9       format(3X,I8,4X,I1,2X,I7,5X,1PG14.6)
    nn = nn + 1
enddo

return

300    write(stderr,10) ios
10    format(1x,'RDWFLX: fatal input error ',i8)
call prerr('RDWFLX')
call xit(1," ")

end subroutine

!-------------------------------------------------------------------------------

subroutine RDTAN (ITANOM,TANOM)

USE MODELDEFINITION, only: NTANOM, NTANOMglobal
USE MESHDATAMODULE,  only: meshdatactx
USE MODELCTX,        only: getrank
use debugmodule,     only: iecho, xit
use iomodule          

! Subroutine to read and print out thermal anomaly data

implicit none
!-pass
integer            :: ITANOM
double precision   :: TANOM
dimension          :: ITANOM(3,*),TANOM(*)
!-locl
logical            ::  end
character(len=49)  :: record
character(len=7)   :: type
character(len=18)  :: fmt
integer, external  :: lnblk

integer            :: nreadglobal

integer            :: im, it, inp, ios
integer            :: mg, ng, n, nn, np, nnp, num, lst
integer            :: nfirst, nlast
double precision   :: df, ff, f

!-init
mg  = 0
num = 0
lst = 0
nreadglobal = 0

fmt = '(2I12,I5,I6,G14.0)'

if (iecho.ge.2) then
    write(stderr,2) fmt
 2   format('1x,if (NTANOM > 0){'/ &
 5X,'*** read thermal anomaly data ***'/ &
 5X,'read ',A17,1x, &
 'ITANOM(1,n),ITANOM(2,n),ITANOM(3,n),TANOM(n)'/ &
 5X,'if (ITANOM(1,n)==0 || ITANOM(1,n)=end) return'/1x,'}')
endif

if (NTANOMglobal.le.0) then
    return
endif

100    call readrc (luin,record,'RDTAN',end)
if (end) then
    goto 300
endif

read (record,fmt,err=500,end=500,iostat=ios) np,im,it,f

if (np.eq.0) then
    goto 300
endif

nreadglobal = nreadglobal + 1

if (meshdatactx%v2p(np).ne.getrank()) then
    goto 100
endif

if (mg.ne.0) then
!        generate data

    if (im.ne.ITANOM(2,lst)) then
        write(stderr,7)
 7        format(1x,'rdtan: mode increment not implemented')
        call xit(1," ")
    endif

    if (it.ne.ITANOM(3,lst)) then
        write(stderr,8)
 8        format(1x,'rdtan: time increment not implemented')
        call xit(1," ")
    endif

    nnp = (np-ITANOM(1,lst))/mg
    df = (f-TANOM(lst))/nnp
    ff = TANOM(lst)

    nfirst = ITANOM(1,num) + mg
    nlast  = np - mg

    do inp=nfirst,nlast,mg
        num = num + 1

        if (num.gt.NTANOMglobal) then
            write(stderr,10) NTANOMglobal
10            format(1x,'rdtan: number of inputs exceeds ',I8)
            call xit(1," ")
        endif

        ITANOM(1,num) = inp
        ITANOM(2,num) = ITANOM(2,lst)
        ITANOM(3,num) = ITANOM(3,lst)

        ff = ff + df
        TANOM(num) = ff
    enddo
endif

num = num + 1
NTANOM = num
if (num.gt.NTANOMglobal) then
    write(stderr,12) NTANOMglobal
12        format(1x,'rdtan: number of nodal points exceeds ',I8)
    call xit(1," ")
endif

ITANOM(1,num) = np
ITANOM(2,num) = im
ITANOM(3,num) = it
TANOM(num)    = f

lst = num
mg  = ng

goto 100

300    if (nreadglobal.lt.NTANOMglobal) then
    write(stderr,14) NTANOMglobal
14       format(1x,'rdtan: actual number of anomaly inputs is less', &
     1x,'than ',I8)
    call xit(1," ")
endif

if (.not.OUTFIL(1)) then
    return
endif

nn = 51

do n=1,NTANOM
    if (nn.gt.50) then
        nn=0
        write(luout(1),15)
15           format(///1X,'T H E R M A L',4X,'A N O M A L Y',4x, &
         'D A T A'//11X,'TIME.GE.0 = APPLY AFTER THIS TIME'/11X, &
         'TIME.LT.0 = APPLY FROM ABS(TIME) ONWARD'///5x, &
         'NODE',5x,'TIME',5x,'TEMPERATURE CHANGE'/)
    endif

    if (ITANOM(2,n).eq.0) then
        type = 'REPLACE'
    else
        type = 'ADD    '
    endif

    write(luout(1),16) ITANOM(1,n),ITANOM(3,n),TANOM(n),type
16        format(1x,i8,3x,i6,9x,f7.1,3x,A7)

    nn=nn+1
enddo

return

500    write(stderr,18) ios
18    format(1x,'rdtan: fatal input error ',i8)
call prerr('rdtan')
call xit(1," ")

end
!-------------------------------------------------------------------------------
subroutine RDTBC (NUMNP,ITMODE)
! reads thermal boundry conditions.

use modeldefinition, only: allTemperaturesFixed
USE MODELDATAMODULE
USE MESHDATAMODULE
USE MODELCTX
use debugmodule, only: iecho, xit
use iomodule          

implicit none

! overwrite NDOF from modeltopology.
integer, parameter :: NDOF=1, MAXP=50


!-pass
integer            :: NUMNP,ITMODE
!-locl
integer            :: idum, n, nn, ntbc

character(len=11)fmt ! index format only 8, but float format is longer

fmt = '(I12,I5)   '

modeldatactx%ITBC = 0
modeldatactx%TBC = 0d0
modeldatactx%ITBCglobal = 0d0

if (iecho.eq.4) then
    write(stderr,1) fmt


1   format(1x,'*** read thermal boundary condition codes ***'/ &
    1x,'node=1; while (node != 0) {'/ &
    5x,'read ',A8,' node,increment,ITBC(1,node)'/ &
    1x,'}')

endif

! second entry is 1 in stead of NDO, because heat is not specific on
! direction.
call IGEN (modeldatactx%ITBC,1,NUMNP,fmt,NTBC,2,"Thermal BCS")


do n=1,NUMNP
    if (modeldatactx%ITBC(1,n).lt.0 .or. &
        modeldatactx%ITBC(1,n).gt.1) then

          write(stderr,*) 'RDTBC: incorrect boundary condition code for node ',n
          write(stderr,*) 'RDTBC: code ', modeldatactx%ITBC(1,n) ,'found. Should be 0 or 1'
          call xit(1," ")
      endif
enddo

fmt = '(I12,G14.0)'

if (iecho.eq.4) then
    write(stderr,3) fmt
 3   format(1x,'*** read thermal boundary condition data ***'/1x, &
 'node=1; while (node != 0) {'/1x, &
 '    read ',A11,' node,increment,TBC(node)'/1x, &
 '}')
endif

call FGEN (modeldatactx%TBC,NDOF,NUMNP,fmt,idum,"thermal boundary conditions")

allTemperaturesFixed = .true.

if (ITMODE.ne.0) then
    do n=1,meshdatactx%Nvlocal!NUMNP
        if (modeldatactx%ITBC(1,n).eq.1) then
            modeldatactx%T(n) = modeldatactx%TBC(n)
		else
			allTemperaturesFixed = .false.
        endif
    enddo
endif

if (OUTFIL(1)) then
    nn=MAXP+1
    do n=1,NUMNP

        if (nn.gt.MAXP) then
            nn = 1
            write(luout(1),6)
        endif

        if (modeldatactx%ITBC(1,n).eq.1) then
            write(luout(1),4) n,'TEMP',modeldatactx%TBC(n)
 4            format(6X,I8,4X,A4,8X,1PG14.6)
            nn = nn + 1
        else if (modeldatactx%ITBC(1,n).eq.3) then
            write(luout(1),4) n,'LOAD',modeldatactx%TBC(n)
            nn = nn + 1
        else if (modeldatactx%ITBC(1,n).eq.0) then
            continue
        else
            write(stderr,5) n,modeldatactx%ITBC(1,n)
 5               format(1x,'RDTBC: unknown nodal boundary condition', &
             1x,'type ',I3)
            call exitp(1)
        endif

    enddo
endif

return

 6   format(1x/1x/1x/1x, &
 'N O D A L   B O U N D A R Y   C O N D I T I O N S'/1x/ &
 5X,'NODE #',17X,'CONDITION'/)
end
!-------------------------------------------------------------------------------
subroutine RDTEMP (T,NUMNP,ITMODE)
! Read initial conditions for temperatures at nodes.

use modelctx,         only: getrank
use meshdatamodule,   only: meshdatactx
use debugmodule,      only: iecho
use iomodule          

implicit none
!-pass
integer          :: NUMNP,ITMODE
double precision :: T(NUMNP)
!-locl
integer :: nread
character(len=12)fmt
fmt = '(I12,G14.0)'

if (iecho.eq.4) then
    write(stderr,1) fmt
 1   format(1x,'if (ITMODE == 1 || ITMODE == 3) {'/ &
 5x,'*** read nodal temperatures ***'/ &
 5x,'node=1; while (node != 0) {'/ &
 9x,'read ',A11,' node,T(node)'/ &
 5x,'}'/1x,'}')
endif

if (ITMODE.ne.1 .and. ITMODE.ne.3) then
    return
endif

call FGEN(T,1,NUMNP,fmt,nread,"initial temperatures")

if (nread .ne. meshdatactx%nvlocal) then
	write(stderr,*) "Thread", getrank(), &
                    "has ", meshdatactx%nvlocal, &
                    "local points, but has read only", nread, &
                    "initial temperatures"

!    write(stderr,2) NUMNP,nread
! 2       format(1x,'RDTEMP: expected ',I12,' temperature inputs,',1x, &
!     'read ',I12)
endif

return
end subroutine

!-------------------------------------------------------------------------------
subroutine RDTMAT (NSD,NUMAT)

USE MATERIALSMODULE, only : materialgen, &
                            matdatactx
use debugmodule, only: iecho, xit
use iomodule          

implicit none


integer, parameter :: MAXP=50

!-pass
integer            :: NSD,NUMAT
!-locl
logical            :: lerr
integer            :: i, j, nn

if (iecho.eq.4) then
    write(stderr,1)
 1   format(1x,'*** thermal material properties ***'/ &
 'n=1; while (n != 0 & n != ''end'') {'/ &
 4x,'read (2i5,6g14.0)',1X, &
 'n,HEAT(n),CP(n),TEXP(n),(COND(i,n),i=1,NSD)'/ &
 '}')
endif

!    call MATGEN (PRPTEM,6,3+NSD,NUMAT,0)
!      conductivities are dependent on direction, hence NSD
call materialgen(matdatactx%PRPTEM,6,3+NSD,NUMAT,0)

!    Check conductivities
lerr = .false.
do i=1,NUMAT
    do j=1,NSD
        if (matdatactx%PRPTEM(3+j,i).le.1e-6) then
            write(stderr,2) j,i
 2               format(1x,'RDTMAT: Conductivity(',I1,') of',1x, &
             'material number ',I6,' <= 1E-6')
            lerr = .true.
        endif
    enddo
enddo

if (lerr) then
    call xit(1," ")
endif


if (OUTFIL(1)) then
    nn = MAXP + 1
    do i=1,NUMAT
        if (nn.gt.MAXP) then
            nn = 1
                if (NSD.eq.1) write(luout(1),4)
                if (NSD.eq.2) write(luout(1),5)
                if (NSD.eq.3) write(luout(1),6)
        endif
        write(luout(1),3) i,(matdatactx%PRPTEM(j,i),j=1,3+NSD)
 3           format(i5,2x,1Pe8.2,5x,0Pf6.1,5x,1Pe7.1,6x,0Pf5.2,3x, &
         0Pf5.2,3x,0Pf5.2)
        nn = nn + 1
    enddo
endif

return

 4   format(///1x, &
 'T H E R M A L   M A T E R I A L   P R O P E R T I E S'// &
 2x,'SET',4x,'HEAT',6x,'SPECIFIC',4x,'THERMAL',5x, &
 'CONDUCTIVITY'/6x,'PRODUCTION',5x,'HEAT',4x,'EXPANSIVITY',8x, &
 'X'/)
 5   format(///1x, &
 ' T H E R M A L   M A T E R I A L   P R O P E R T I E S'// &
 2x,'SET',4x,'HEAT',6x,'SPECIFIC',4x,'THERMAL',8x, &
 'CONDUCTIVITY'/6x,'PRODUCTION',5x,'HEAT',4x,'EXPANSIVITY',8x, &
 'X',6x,'Y'/)
 6   format(///1x, &
 ' T H E R M A L   M A T E R I A L   P R O P E R T I E S'// &
 2x,'SET',4x,'HEAT',6x,'SPECIFIC',4x,'THERMAL',13x, &
 'CONDUCTIVITY'/6x,'PRODUCTION',5x,'HEAT',4x,'EXPANSIVITY',8x, &
 'X',6x,'Y',6x,'Z'/)

end subroutine

!-------------------------------------------------------------------------------

subroutine SETTID(IDTglobal,IDTXglobal,ITBCGLOBAL)

USE MODELDATAMODULE
USE MESHDATAMODULE
USE MODELDEFINITION
USE MODELCTX
USE ALGEBRA,         only: clear

implicit none

!
!   because of global numbering, the arrays need to be global as well.
!   At a later stage the per processor equation numbers are determined.
!   nfault and nslip are global arrays
! Sets equation numbers in IDT and IDTX arrays
!
!-pass
integer IDTglobal(*),IDTXglobal(*),ITBCGLOBAL(*)!,NSLIP(5,*),NFAULT(3,*)
!-locl
integer n,j
integer nloc
!-init
NTEQglobal = 0
NTEQlocal = 0

call CLEAR (IDTglobal,NUMNP,"IDTglobal")

if (NUMSNglobal+NUMFNglobal.gt.0.and.IDIFT.eq.1) then
    call CLEAR (IDTXglobal,NUMNP,"IDTXglobal")
endif



!write(*,*) "thread", getrank(), "v2p", meshdatactx%v2p

do n=1,NUMNP  ! loop through the global elements...

    ! ITBCGLOBAL is 1 for every point for which a thermal boundary condition
    ! has been set. It is zero for all other points.
    if (ITBCGLOBAL(n).ne.1) then
        ! this point has no thermal boundary condition, so it needs a thermal
        ! equation number.

#ifdef SPARSE
        if (getsize().eq.1) then
            nloc = n
        else
            nloc = meshdatactx%Glo2LocVertex(n)
        endif
#else
        nloc = n
#endif


        ! check whether there is differential trouble
        if (NUMSLP.gt.0.and.IDIFT.eq.1) then
            do j=1,NUMSLP
                if (n.eq.modeldatactx%NSLIPglobal(2,j)) goto 100
            enddo
        endif

        if (NUMFN.gt.0.and.IDIFT.eq.1) then
            do j=1,NUMFN
                if (n.eq.modeldatactx%NFAULTglobal(2,j)) goto 100
            enddo
        endif

        ! no differential trouble, skip the next bit
        goto 200

        ! yes differential trouble, add an extra differential thermal equation 
100     NTEQglobal = NTEQglobal + 1
        IDTXglobal(n) = NTEQglobal

        if (meshdatactx%v2p(n).eq.getrank() .and. nloc.ne.-1) then
            NTEQlocal = NTEQlocal + 1
            modeldatactx%IDTX(nloc) = NTEQglobal
        endif
!
200     NTEQglobal = NTEQglobal + 1


        IDTglobal(n)  = NTEQglobal

        if (meshdatactx%v2p(n).eq.getrank() .and. nloc.ne.-1) then
            NTEQlocal = NTEQlocal + 1

!            write(*,*) "thread", getrank(), &
!                       "writing IDT", nloc-meshdatactx%lowestLocalVertexID, &
!                       "of", size(modeldatactx%IDT,1)

!            modeldatactx%IDT(nloc) = NTEQglobal
            modeldatactx%IDT(nloc-meshdatactx%lowestLocalVertexID) = NTEQglobal

        endif

    endif
enddo

return
end subroutine

!-----------------------------------------------------------------------

subroutine storerr (N,string)

use debugmodule,     only: xit
use iomodule          

implicit none
!-pass
integer    :: N
character  :: string*(*)
!-locl
integer    :: l
integer, external :: lnblk

l = LEN(string)

if (l.gt.0) then
    l = lnblk(string)

    if (l.gt.0) then
        write(stderr,1) N,string(1:lnblk(string))
 1           format(//1x,'Insufficient storage by ',I12, &
         ' single precision words'/ &
         1x,'Error while allocating space for ',80A)
    else
        write(stderr,2) N
 2           format(//1x,'Insufficient storage by ',I12, &
         ' single precision words')
    endif
else
    write(stderr,2) N
endif

call xit(1," ")
return
end subroutine

!-------------------------------------------------------------------------------

subroutine TCOLHT (LMT,LMTX,NEN)

USE MODELDEFINITION, only: idift, numfn, numslp, nteq
use modeldatamodule, only: itdiag
USE ALGEBRA,         only: clear
use constants,       only: intmax
use modelctx,        only: getrank
use meshdatamodule,  only: meshdatactx

!
! calculates column heights
!
implicit none
!-pass
integer   :: NEN
integer   :: lmt, lmtx
dimension :: LMT(NEN,*),LMTX(NEN,*)
!-locl
integer   :: k,minimum,j, m, numeq


ITDIAG = 0

!write(*,*) "tcolht rank", getrank(), "numel:", numel

!do k=1,NUMEL
do k=1,meshdatactx%nelocal


    minimum = intmax
!   Find smallest equation number in current element
    do j=1,NEN
        numeq = LMT(j,k)
        if (numeq.ne.0) then
            minimum = MIN0(minimum,numeq)
        endif

        if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
            if (LMTX(j,k).ne.0) then
                minimum = MIN0(minimum,IABS(LMTX(j,k)))
             endif
        endif

    enddo

!   Set column heights
    do j=1,NEN
        numeq = LMT(j,k)
        if (numeq.ne.0) then
            m = numeq - minimum
!            write(*,*) "checking Thermal Eq ", numeq, minimum, m
            if (m.gt.ITDIAG(numeq)) then
                ITDIAG(numeq) = m
            endif
        endif
    enddo

    if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
        do j=1,NEN
            numeq = IABS(LMTX(j,k))
            if (numeq.ne.0) then
                m = numeq - minimum
                if (m.gt.ITDIAG(numeq)) then
                    ITDIAG(numeq) = m
                endif
            endif
        enddo
    endif
enddo
!
return
end subroutine

!-----------------------------------------------------------------------

subroutine PRFUNC (MODE)

! subroutine to print function to calculate initial temperatures
use debugmodule,     only: xit
use iomodule          

implicit none
!-pass
integer            :: mode
!-locl
character(len=100) :: file,record,getfil
integer            :: lu

! todo really? not nextlu? What if 44 is already assigned?
lu = 44

if (.not.OUTFIL(1)) then
    return
endif

if (MODE.eq.0) then
    file = getfil('temp0.F')
else
    file = getfil('addan.F')
endif

if (MODE.eq.0) then
    write(luout(1),15)
15       format(1x,'PROCEDURE FOR CALCULATING INITIAL TEMPERATURES:' &
     //5X,'+',74('-'),'+')
else
    write(luout(1),17)
17       format(1x,'PROCEDURE FOR CALCULATING TEMPERATURE ANOMALIES:' &
     //5X,'+',74('-'),'+')
endif

call openf (lu,file,'old')

100    read(lu,20,err=300,end=200) record
20    format(A72)
write(luout(1),30) record
30    format(5X,'| ',A72,' |')
goto 100

200    write(luout(1),35)
35    format(5X,'+',74('-'),'+')
call closef(lu)
return

300    call prerr(file)
call xit(1," ")

end subroutine

!-----------------------------------------------------------------------

subroutine nittem (T,X,NSD,NUMNP)

! Routine for calculating initial temperature distribution

implicit none
!-pass
integer                    :: NSD,NUMNP

double precision           :: T(NUMNP),X(NSD,NUMNP)
!-locl
double precision, external :: temp0
integer                    :: n

do n=1,NUMNP
    T(n) = temp0 (X(1,n),T(n),NSD)
enddo

return
end
