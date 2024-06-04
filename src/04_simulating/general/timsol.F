subroutine TIMSOL ()

!-------------------------------------------------------------------------------
!       Subroutine to controlling the time stepping of mechanical and thermal
!       problem.
!-------------------------------------------------------------------------------

USE MATERIALSMODULE
USE MESHDATAMODULE
USE MODELDATAMODULE
USE MODELDEFINITION
USE TIMESTEPMODULE
USE MODELCTX
USE ALGEBRA,          only: IDFORM, idig
USE CONSTANTS,        only: year, &
                            writeRestart
use filemodule,       only: statusFileName
use debugmodule,      only: lustat, &
                            itest, &
                            debug, &
                            petsctime, &
                            owntime, &
                            time_start, &
                            iecho, &
                            LogTimeStamp, &
                            startClock, &
                            delay, &
                            parallelLog, &
                            CheckArrayForNans, &
                            iAmAtTimestep, isfileopen
use modeltopology,    only: NSD,NDOF,NSTR,NEN
use iomodule

#ifdef SPARSE
USE AOMODULE
#endif

implicit none


!-locl
integer :: IDMAT(NSTR)
!        logical there,output
logical there
character(len=50) string
integer csystem,ierr
integer rnk, logsize, inneri, innerj
external stoper,csystem
logical DoOutM, DoOutT

integer :: i, ii, j, inner, igrp, jlocal, jnnerj, ncntr

integer :: nfirst, naxstp
double precision  :: tdt

!-init
data NFIRST,NAXSTP,TDT/0,0,0d0/
ITPRT  = 1
IMPRT  = 1
IPRMAT = 1
ITCALC = 1
IMTWRK = 0

! Construct IDMAT, addresses of square array diagonals
call IDFORM(IDMAT,NSTR)

! initialize mechanical problem-related variables, pointers and arrays
call MINIT (NSTR)

! initialize thermal problem-related variables, pointers and arrays
!call TINIT ()   ! but it doesn't...


if (IADVEC.ne.0) then
    if (iecho.eq.8) then
        call ParallelLog("timsol", "calling load")
    endif

    call LOAD (modeldatactx%ID, &
    modeldatactx%IBOND,modeldatactx%BOND,modeldatactx%D, &
    modeldatactx%DELD, modeldatactx%B,meshdatactx%X, &
    NUMNP,NEQ,0)
endif

! initialize time stepping variables
call INITIM (timestepdatactx, outputcontroldatactx, NSTEP,NFIRST,NAXSTP,IGRP)

if (debug) write(*,*) '**** done INITIM ', NSTEP,NFIRST,NAXSTP,IGRP


! initialize timers
petsctime = 0.0
owntime = 0.0
TDT = TIME

! Loop over time step *groups*
loop_time_step_groups: do ii=IGRP,NINTG
    i = ii

    if (debug) write(*,*) '**** looping over timestep group ', i, 'with stepsize: ', DELTP

!   Define constants to control stepping in current group
    if (i.ne.IGRP) then
        call CONST(timestepdatactx,i,NAXSTP,NFIRST)
    endif

!   Loop over timesteps in current timestep group
    loop_time_steps_in_group: do j=NFIRST,NAXSTP

        if (iecho.eq.8) then
            call ParallelLog("timsol", "starting a new timestep")
        endif

        NSTEP = NSTEP + 1
        TIME  = TIME  + DELTP
        TDT   = TDT   + DELTP

        if (iecho.eq.8) then
            call ParallelLog("timsol", "determining whether output is necessary")
        endif

        ! DoOutM: output for mechanical timestep
        if (NMPRT.lt.0) then
            DoOutM = (mod(NSTEP,-NMPRT).eq.0)
        else
            DoOutM = output(outputcontroldatactx%IMPRINT, NMPRT,IMPRT,NSTEP)        
        endif

        ! DoOutT: output for thermal timestep
        DoOutT = output(outputcontroldatactx%ITPRINT, NTPRT,ITPRT,NSTEP)

        ! Construct filename for FEM output and open file
        if (DoOutM.or.DoOutT) then
            ncntr = log10(dble(MAXITIME)) + 1
            FEDSK = FEDSK(1:FEDSKRT)//'.'

            jlocal = j
            do inner=1,ncntr
                inneri = 10**(ncntr-inner)
                jnnerj = jlocal/inneri
                FEDSK = FEDSK(1:FEDSKRT+inner)//char(48 + jnnerj)
                jlocal = jlocal - jnnerj * inneri
            enddo

            ! add counter for the processor ID
            FEDSK = FEDSK(1:FEDSKRT+ncntr+1)//'.'
            rnk = getrank()

            if (getsize().gt.1) then
                logsize = log10(dble(getsize())) + 1
            else
                logsize = 1
            endif

            do i=1,logsize
                inneri = 10**(logsize-i)
                innerj = rnk/inneri
                FEDSK = FEDSK(1:FEDSKRT+ncntr+1+i)//char(48+innerj)
                rnk = rnk - innerj * inneri
            enddo

            if (iecho.eq.3) then
                write(*,*) 'opening FEDSK file: ', trim(FEDSK)
            endif

            if (ascout) then
                call openf (luout(2),FEDSK,'unknown')
            else
                call fbopen (luout(2),FEDSK,'unknown')
            endif

        endif
 
        call followtime (0)

        if (OUTFIL(1)) then
            write(luout(1),10) TIME,TIME/YEAR, TIME/YEAR*1d-6
10                 format(//1x,80('-')/1x, &
              'OUTPUT WHICH FOLLOWS IS AT TIME =', &
              1PG8.2,' (',1PG8.2,' YEAR, ',0PF7.1,' MA)')
        endif

        call buildyearstring (NSTEP,TIME,string)

        if (iecho.eq.1) then
            write(stdout,*) trim(string)
        endif
        if (debug) write(*,*) "calling tmarch with TDT: ", TDT
        if (iecho.eq.8) then
            call ParallelLog("timsol", "calling tmarch")
        endif

!       Calculate temperature distribution
        call TMARCH (TDT,NSD,NDOF,NEN,NSTEP)

        if (iecho.eq.8) then
            call ParallelLog("timsol", "finished tmarch")
        endif

!       determine whether implicit code needs to be invoked
        IMWORK = 0
        if (REFACT.or.(ALFAP.ne.0.0.and.MOD(j,MAXIT).eq.0)) then
            IMWORK = 1
        endif
        REFACT = .false.


#ifdef EXPERIMENTALTRANSIENTDISCREEP
! the jacobian matrix for dislocation creep has not been implemented yet, so switch off implicit option
        IMWORK = 0
        if (iecho.eq.1) write(stdout,35)
35            format(1x,'setting time marching to explicit by default, no jacobian defined')
        if (MODE.eq.6 .and. IDIG(ICVIS,4,1).eq.2) then

           if (NSTEP.eq.1) then
              if (iecho.eq.8) then
                 call ParallelLog("timsol", "calling SetThresholdStress")
              endif
           ! only initialise internal stresses for the first time step
              if (iecho.eq.1) write(stdout,31)
31            format(1x,'initialising internal stresses: threshold stress')
              ! define threshold                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
              call SetThresholdStress(modeldatactx%STNTHRES)
              if (iecho.eq.8) then
                 call ParallelLog("timsol", "finished SetThresholdStress")
              endif
              if (iecho.eq. 1) write(stdout,32)
32            format(1x,'initialising internal stresses: Taylor stress')
              if (iecho.eq. 1) write(stdout,33) 
33            format(1x,'set initial Taylor stress as fraction of initial elastic stress')
              if (iecho.eq.8) then
                 call ParallelLog("timsol", "calling InitTaylorStess")
              endif
              ! initialise Taylor stress, as function of elastic stresses
              call InitTaylorStess(modeldatactx%STNTAY)
                  
              if (iecho.eq.8) then
                 call ParallelLog("timsol", "finished InitTaylorStess")
              endif
           endif
        endif
#endif


!       Calculate material properties
        if (ICVIS.ne.0 .and.(MODE.eq.2.or.MODE.eq.6) ) then
            if (iecho.eq.1) write(stdout,30)
30          format(1x,'Calculating material properties')
            call SETPRP (outputcontroldatactx%IMPRINT, &
             NSD,NDOF,NEN,NSTR)
        endif

!       Write calculated viscosities and powers to file
        if (IDIG(ICVIS,4,1).ne.0 .and. MODE.eq.6 .and. DoOutM) then
            if (LMAT) then
                call WrtCrp (matdatactx%PRPMAT,NUMEL)
            else
                call WrtCrp (matdatactx%PRPMAT,NUMAT)
            endif
        endif

#ifdef EXPERIMENTALTRANSIENTDISCREEP
! Write non-linear viscosity parameters to file
!           if (iecho.eq.1) write(stdout,34)
!34          format(1x,'still to do: output of non-linear viscosity parameters to file')
#endif

!       Write viscoplastic parameters to file
        if (IDIG(ICVIS,4,2).ne.0 .and. (MODE.eq.2.or.MODE.eq.6) .and. DoOutM) then
            if (LMAT) then
                call WrtPls (NUMEL)
            else
                call WrtPls (NUMAT)
            endif
        endif

        if (iecho.eq.8) then
            call ParallelLog("timsol", "Calling Mmarch")
        endif

!       Perform mechanical calculation
        call MMARCH (IDMAT,NSD,NDOF,NSTR,NEN,NSTEP)

        if (iecho.eq.6) then
            call LogTimeStamp("mmarch step")
        endif

        if (debug) then
            call ParallelLog("timsol", "survived an mmarch round")
        endif

        if (itest.eq.2) then
            return
        endif
 
!       Write restart file?
        if (ISTART.ne.0) then
            if (MOD(NSTEP,IABS(ISTART)).eq.0) then
                if (iecho.eq.8) then
                    call ParallelLog("timsol", "resting the io")
                endif
                call RESTIO (modeldatactx%D,modeldatactx%DXE, &
                             modeldatactx%STN,modeldatactx%STR,modeldatactx%BTOT, &
                             modeldatactx%T,modeldatactx%FTOT,modeldatactx%VPSTR, &
                             NDOF,NEN,NSTR,writeRestart)
            endif
        endif

        if (iecho.gt.0) then
            call openf (lustat,statusFileName,'unknown')
            write(lustat,40) NSTEP
40          format('Busy: completed ',I6,' steps')
            call closef (lustat)
        else
            if (MOD(NSTEP,10).eq.0) then
!               update every ten timesteps.
                call openf (lustat,statusFileName,'unknown')
                write(lustat,40) NSTEP
                call closef (lustat)
            endif
        endif

        tdt = 0d0

        if (iecho.eq.8) then
            call ParallelLog("timsol", "Finished timestep.")
        endif

    enddo loop_time_steps_in_group

enddo loop_time_step_groups

if (iecho.eq.8) then
    call ParallelLog("timsol", "Finished all timesteps.")
endif

! force write of restart file at end of run
!if (ISTART.ne.0) then
!    call RESTIO (modeldatactx%D,modeldatactx%DXE, &
!    modeldatactx%STN,modeldatactx%STR,modeldatactx%BTOT, &
!    modeldatactx%T,modeldatactx%FTOT,modeldatactx%VPSTR, &
!    NDOF,NEN,NSTR,0)
!endif

call followtime (MAX(NSTEP,1))

return
end subroutine
!
!-------------------------------------------------------------------------------
!
subroutine buildyearstring (NSTEP,TIME,string)

USE CONSTANTS, only: year
use modelctx,  only: getrank
use iomodule, only: stderr

!
implicit none

character(len=*) string
integer nstep
double precision :: time




!write(*,*) "buildyearstring has time/year", TIME, year

 if (getrank().eq.0) then

 if (NSTEP .lt. 10) then
write(string,1,err=100) NSTEP,TIME/YEAR
 1       format('Working on time(',I1,')=',1PE13.4,' year')
 else if (NSTEP .lt. 100) then
write(string,2,err=100) NSTEP,TIME/YEAR
 2       format('Working on time(',I2,')=',1PE13.4,' year')
 else if (NSTEP .lt. 1000) then
write(string,3,err=100) NSTEP,TIME/YEAR
 3       format('Working on time(',I3,')=',1PE13.4,' year')
 else if (NSTEP .lt. 10000) then
write(string,4,err=100) NSTEP,TIME/YEAR
 4       format('Working on time(',I4,')=',1PE13.4,' year')
 else if (NSTEP .lt. 100000) then
write(string,5,err=100) NSTEP,TIME/YEAR
 5       format('Working on time(',I5,')=',1PE13.4,' year')
 else if (NSTEP .lt. 1000000) then
write(string,6,err=100) NSTEP,TIME/YEAR
 6       format('Working on time(',I6,')=',1PE13.4,' year')
 else if (NSTEP .lt. 10000000) then
write(string,7,err=100) NSTEP,TIME/YEAR
 7       format('Working on time(',I7,')=',1PE13.4,' year')
 else if (NSTEP .lt. 100000000) then
write(string,8,err=100) NSTEP,TIME/YEAR
 8       format('Working on time(',I8,')=',1PE13.4,' year')
 else if (NSTEP .lt. 1000000000) then
write(string,9,err=100) NSTEP,TIME/YEAR
 9       format('Working on time(',I9,')=',1PE13.4,' year')
 else
write(string,10,err=100) NSTEP,TIME/YEAR
10       format('Working on time(',I12,')=',1PE13.4,' year')
 endif

 endif


!
 return

100    write(stderr,101)
  101   format(1x,'buildyearstring: internal write error')

 return
 end
!-------------------------------------------------------------------------------

subroutine INITIM (timestepdat, outputcontroldat, NSTEP,NFIRST,NLAST,IGRP)

use modeldefinition, only: mode, istart
USE TIMESTEPMODULE,  only: timestepdata, outputcontroldata, &
                           deltp, alfap, &
                           NMPRT, IMPRT, NTPRT, ITPRT, ntcalc, ITCALC, nintg, TALFP
USE MODELCTX
use debugmodule,     only: xit
use iomodule,        only: stderr, stdout

!
! Initializes time integration variables. In case of restart, NSTEP specifies
! output that was succesfully read from the restart file. IGRP, NFIRST, IMPRT,
! ITPRT, ITCALC and TIME are set to next time step quantities.

implicit none

!-pass
   type(timestepdata)      :: timestepdat
   type(outputcontroldata) :: outputcontroldat
integer          :: NSTEP,NFIRST,NLAST,IGRP
!-local
double precision :: time
integer          :: n, i, mstep

! for stderr

TIME = 0d0

if (ISTART.eq.0 .or. NSTEP.eq.0) then
    IGRP   = 1
    NFIRST = 1
    NLAST  = timestepdat%MAXSTP(1)
else
    mstep = 0
    NFIRST = 0
 
    do n=1,NINTG
        IGRP = n
        do i=1,timestepdat%MAXSTP(n)
            NFIRST = NFIRST + 1
            TIME  = TIME + timestepdat%DELT(IGRP)
            mstep = mstep + 1
            if (mstep.eq.NSTEP) goto 100
        enddo
    enddo

    goto 200

100 NFIRST = NFIRST + 1
    NLAST = NFIRST + timestepdat%MAXSTP(IGRP) - i - 1

    if (NFIRST.gt.NLAST) then
        IGRP = IGRP + 1
        if (IGRP.gt.NINTG) goto 200
        NLAST = NFIRST + timestepdat%MAXSTP(IGRP) - 1
    endif

! mechanical solution
if (MODE.eq.2.or.MODE.eq.6) then
    if (NSTEP.ge. &
       outputcontroldat%IMPRINT(NMPRT)) then
!           11-05-2011
!           this may seem an error, allowing calculations beyond
!           the last output moment. But this may be required in
!           case of restarts with a lot of calculations but
!           no need for many outputs.
        write(stderr,2) NSTEP, &
         outputcontroldat%IMPRINT(NMPRT)
!                goto 300
    endif
    do while (outputcontroldat%IMPRINT(IMPRT).le.NSTEP .and. IMPRT.lt.NMPRT)
        IMPRT = IMPRT + 1
    enddo
endif

! temperature solution
if (MODE.eq.5.or.MODE.eq.6) then
    if (NSTEP.ge. &
        outputcontroldat%ITPRINT(NTPRT)) then
!               this might be an error,
!               however, using an initial temperature field
!               without further thermal calculations is
!               perfectly ok
        write(stdout,3) NSTEP, &
          outputcontroldat%ITPRINT(NTPRT)
!                goto 400
    endif
    do while (outputcontroldat%ITPRINT(ITPRT).le.NSTEP .and. ITPRT.lt.NTPRT)
        ITPRT = ITPRT + 1
    enddo
    do while (timestepdat%ICALC(ITCALC).le.NSTEP .and. ITCALC.lt.NTCALC)
        ITCALC = ITCALC + 1
    enddo
endif
 endif
 
DELTP = timestepdat%DELT(IGRP)
ALFAP = timestepdat%ALFA(IGRP)

if (MODE.ge.3) then
    TALFP = timestepdat%TALF(IGRP)
endif
!
return
!
200     write(stderr,1) NSTEP+1
 1   format(1x,'RESTART: NSTEP = ',I7,' not in time step groups')
call xit(1," ")

300     write(stderr,2) NSTEP, &
   outputcontroldat%IMPRINT(NMPRT)
 2   format(1x,'RESTART Warning: NSTEP (',I7, &
 ') >= IMPRINT(NMPRT) (',I7,'):'// &
 ' no mechanical outputs will be written which makes little sense')
call xit(1," ")

400     write(stdout,3) NSTEP, &
          outputcontroldat%ITPRINT(NTPRT)
 3   format(1x,'RESTART Warning: NSTEP (',I7, &
 ') >= ITPRINT(NTPRT) (',I7,'):'// &
 ' no thermal outputs will be written which makes little sense')
call xit(1," ")

!
end
!-------------------------------------------------------------------------------
subroutine CONST (timestepdat,ITG,LAST_TIME_STEP_NUMBER,FIRST_TIME_STEP_NUMBER)

! Defines time step parameters:
!    ITG                        timestep group number
!    ALFAP                      implicitness for mechanical problem
!    TALFP                      implicitness for thermal problem
!    DELTP                      time step size
!    FIRST_TIME_STEP_NUMBER     overall number of the first timestep in the current timestep group
!    LAST_TIME_STEP_NUMBER      overall number of the first timestep in the current timestep group

USE TIMESTEPMODULE,  only: DELTP,  &
                           ALFAP,  &
                           TALFP,  &
                           timestepdata
USE MODELDEFINITION, only: MODE
 
 
implicit none
!-pass
type(timestepdata) :: timestepdat
integer            :: ITG, LAST_TIME_STEP_NUMBER, FIRST_TIME_STEP_NUMBER

DELTP  = timestepdat%DELT(ITG)
ALFAP  = timestepdat%ALFA(ITG)

if (MODE.ge.3) then
   TALFP  = timestepdat%TALF(ITG)
endif

FIRST_TIME_STEP_NUMBER = LAST_TIME_STEP_NUMBER + 1
LAST_TIME_STEP_NUMBER  = FIRST_TIME_STEP_NUMBER + timestepdat%MAXSTP(ITG) - 1

   return
   end
!-------------------------------------------------------------------------------
subroutine MINIT (NSTR)

!USE MODELCTX
USE MODELDEFINITION, only: mode
USE ALGEBRA,         only: clear
use modeldatamodule, only: modeldatactx

! Initialize variables, pointers and arrays of the mechanical problem

implicit none
!-pass
integer :: NSTR


if (MODE.eq.2 .or. MODE.eq.6) then

!        Clear total strain rate array
#ifndef EXPERIMENTAL_NORATESTORAGE
    modeldatactx%rate = 0d0
#endif
!    call CLEAR(modeldatactx%rate,NUMEL*NSTR)

endif

return
end
!-------------------------------------------------------------------------------
