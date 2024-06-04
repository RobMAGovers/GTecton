MODULE TIMESTEPMODULE

!     this module handles both time step parameters, information
!     concerning output of data, and parameters controlling other
!     actions related to time stepping

implicit none

public :: allocateoutputcontroldata
public :: deallocateoutputcontroldata
public :: allocatetimestepdata
public :: deallocatetimestepdata
public :: CONTRL
public :: RDTIME ! interprets timestep parameters as read from TECIN.DAT
public :: PUTOUT
public :: RDTSTP    
public :: RDMPRT
public :: chktim
public :: chkout
public :: output
public :: SetTimeStepSize
public :: TimDif

double precision TIME
integer :: MAXITIME ! total number of iterations that the simulation runs
double precision :: DELTP,ALFAP,TALFP
integer IMPRT,ITPRT,IPRMAT,FaultVelo
integer ITCALC
integer MAXIT, MAXTIT

integer NMPRT, NPRMAT
integer NTPRT

! todo, replace common block stepct with this
integer :: nstep, ITWORK,IMWORK,IMTWRK
logical :: refact


TYPE outputcontroldata
    integer, ALLOCATABLE :: IMPRINT(:)
    integer, ALLOCATABLE :: ITPRINT(:)
    integer, ALLOCATABLE :: IMATPR(:)
END TYPE outputcontroldata

integer NINTG, NTCALC

TYPE timestepdata
    integer, ALLOCATABLE :: MAXSTP(:)
    double precision, ALLOCATABLE :: DELT(:)
    double precision, ALLOCATABLE :: ALFA(:)
    double precision, ALLOCATABLE :: TALF(:)
    integer, ALLOCATABLE :: ICALC(:)
END TYPE timestepdata

type(outputcontroldata) :: outputcontroldatactx
type(timestepdata)      ::      timestepdatactx

contains

!-----------------------------------------------------

subroutine allocateoutputcontroldata(outputcontroldat,flag)

use debugmodule, only: xit

implicit none
type(outputcontroldata) :: outputcontroldat
integer flag,ERROR,tmp

if (flag.eq.0) then
    tmp = MAX(1,NMPRT)

    if (allocated(outputcontroldat%IMPRINT))then
           deallocate(outputcontroldat%IMPRINT)
    endif
    ALLOCATE(outputcontroldat%IMPRINT(tmp), STAT=ERROR)

    tmp = MAX(1,NPRMAT)

    if (allocated(outputcontroldat%IMATPR))then
        deallocate(outputcontroldat%IMATPR)
    endif
    ALLOCATE(outputcontroldat%IMATPR(tmp),  STAT=ERROR)
else
    tmp = MAX(1,NTPRT)

    if (allocated(outputcontroldat%ITPRINT))then
        deallocate(outputcontroldat%ITPRINT)
    endif
    ALLOCATE(outputcontroldat%ITPRINT(tmp), STAT=ERROR)
endif

if (ERROR.ne.0) then
    write(*,*) 'could not allocate outputcontroldata. Error ', ERROR
    call xit(1," ")
endif


end subroutine allocateoutputcontroldata

!-----------------------------------------------------

subroutine deallocateoutputcontroldata(outputcontroldat,flag)

implicit none

type(outputcontroldata) :: outputcontroldat
integer flag,ERROR

if (flag.eq.0) then
    if (allocated(outputcontroldat%IMPRINT))then
        deallocate(outputcontroldat%IMPRINT)
    endif
    DEALLOCATE(outputcontroldat%IMPRINT, STAT=ERROR)
    if (allocated(outputcontroldat%IMATPR))then
        deallocate(outputcontroldat%IMATPR)
    endif
    DEALLOCATE(outputcontroldat%IMATPR, STAT=ERROR)
else
    if (allocated(outputcontroldat%ITPRINT))then
        deallocate(outputcontroldat%ITPRINT)
    endif
        DEALLOCATE(outputcontroldat%ITPRINT, STAT=ERROR)
endif

end subroutine deallocateoutputcontroldata

!------timestepdata---------------------------------------------------
subroutine allocatetimestepdata(timestepdat,flag)

use modelctx,    only: getrank
use debugmodule, only: parallellog, &
                       iecho, &
                       Verified_MPI_Barrier, &
                       allocateError

implicit none

type(timestepdata) :: timestepdat
integer            :: flag, ERROR,tmp

if (iecho.eq.8) then
    call ParallelLog("allocateTimeStepdata", "started")
endif


if (flag.eq.0) then
    tmp = MAX(1,NINTG)
    if (allocated(timestepdat%MAXSTP)) then
        deallocate(timestepdat%MAXSTP)
    endif
    ALLOCATE(timestepdat%MAXSTP(tmp), STAT=ERROR)

    if (allocated(timestepdat%DELT)) then
        deallocate(timestepdat%DELT)
    endif
    ALLOCATE(timestepdat%DELT(tmp), STAT=ERROR)

    if (allocated(timestepdat%ALFA)) then
        deallocate(timestepdat%ALFA)
    endif
    ALLOCATE(timestepdat%ALFA(tmp), STAT=ERROR)
else
    tmp = MAX(1,NINTG)
    if (allocated(timestepdat%TALF)) then
        deallocate(timestepdat%TALF)
    endif
    ALLOCATE(timestepdat%TALF(tmp), STAT=ERROR)

    tmp = MAX(1,NTCALC)
    if (allocated(timestepdat%ICALC)) then
        deallocate(timestepdat%ICALC)
    endif
    ALLOCATE(timestepdat%ICALC(tmp), STAT=ERROR)
endif
call allocateError("timestep data", ERROR)


end subroutine allocatetimestepdata

!-----------------------------------------------------

subroutine deallocatetimestepdata(timestepdat,flag)

implicit none

type(timestepdata) :: timestepdat
integer flag, ERROR

if (flag.eq.0) then
    if (allocated(timestepdat%MAXSTP)) then
        deallocate(timestepdat%MAXSTP)
    endif
    DEALLOCATE(timestepdat%MAXSTP, STAT=ERROR)
    if (allocated(timestepdat%DELT)) then
        deallocate(timestepdat%DELT)
    endif
    DEALLOCATE(timestepdat%DELT, STAT=ERROR)
    if (allocated(timestepdat%ALFA)) then
        deallocate(timestepdat%ALFA)
    endif
    DEALLOCATE(timestepdat%ALFA, STAT=ERROR)
else
    if (allocated(timestepdat%TALF)) then
        deallocate(timestepdat%TALF)
    endif
    DEALLOCATE(timestepdat%TALF, STAT=ERROR)
    if (allocated(timestepdat%ICALC)) then
        deallocate(timestepdat%ICALC)
    endif
    DEALLOCATE(timestepdat%ICALC, STAT=ERROR)
endif

end subroutine deallocatetimestepdata

!-----routine to read time step parameters; called from READFE, in TECIN.DAT.

subroutine CONTRL (timestepdat, outputcontroldat)

USE CONSTANTS, only: year
use modelctx,  only: getrank
use debugmodule, only: iecho, Verified_MPI_Barrier, parallellog, xit, debug
use iomodule,       only: luin, luout, stderr, outfil
!
! program to read in number of time step, time step and alpha info
!
!-pass
type(timestepdata)      :: timestepdat
type(outputcontroldata) :: outputcontroldat

double precision :: factor

!-locl
integer fnblnk,lnblk,i, ios, nread,m,n
character(len=5) UNIT
character(len=2) tmp
dimension UNIT(5000)
external fnblnk,lnblk
character(len=25000) :: timesteprecord
!-init
MAXITIME = 0

if (iecho.eq.8) then
    call parallellog("Contrl","Entered contrl")
endif

if (NINTG.gt.5000) then
    stop 'CONTRL: overflow on "unit" array'
endif
!
if (iecho.eq.2) then
    write(stderr,1)
endif

!################## read timestep data from TECIN.DAT ###################
 1 format(1x,'*** read time step data ***'/ &
 1x,'nread = MAX(1,NINTG)'/ &
 1x,'read (5000i5  ) (timestepdat%MAXSTP(i), i=1,nread)'/ &
 1x,'read (5000f5.0) (timestepdat%DELT(i), i=1,nread)'/ &
 1x,'read (5000a5  ) (UNIT(i), i=1,nread)'/ &
 1x,'read (5000f5.0) (timestepdat%ALFA(i), i=1,nread)')
!
nread = MAX(NINTG,1)

!-----------------------------------------------------------
if (debug) then
    write(*,*) "Trying to read", nread, "timestep groups"
endif
!-----------------------------------------------------------

read(luin,'(a)') timesteprecord
if (debug) then
    write(*,*) "Read time step group record: ", trim(timesteprecord)
endif

read(timesteprecord,3,err=600,end=700,iostat=ios) &
 (timestepdat%MAXSTP(i),i=1,nread)
 3 format(5000i5)

if (debug) then
    write(*,*) "Extracted time step groups: ", timestepdat%MAXSTP(1:nread)
endif
!-----------------------------------------------------------

read(luin,'(a)') timesteprecord
if (debug) then
    write(*,*) "Read time duration record: ", trim(timesteprecord)
endif

read(timesteprecord,4,err=610,end=710,iostat=ios) &
 (timestepdat%DELT(i),i=1,nread)
 4 format(5000f5.0)

if (debug) then
    write(*,*) "Extracted time step durations: ", timestepdat%DELT(1:nread)
endif


!-----------------------------------------------------------


! units are translated to multiples of seconds below, hence
! they are not part of the timestepdat structure.

read(luin,'(a)') timesteprecord
if (debug) then
    write(*,*) "Read time step unit record: ", trim(timesteprecord)
endif

read(timesteprecord,5,err=620,end=720,iostat=ios) &
 (UNIT(i),i=1,nread)
 5 format(5000A5)

if (debug) then
    write(*,*) "Extracted time step units: ", UNIT(1:nread)
endif

!-----------------------------------------------------------

read(luin,'(a)') timesteprecord
if (debug) then
    write(*,*) "Read time step alpha record: ", trim(timesteprecord)
endif

read(timesteprecord,6,err=630,end=730,iostat=ios) &
 (timestepdat%ALFA(i),i=1,nread)
 6 format(5000F5.0)

if (debug) then
    write(*,*) "Extracted time steps alpha: ", timestepdat%ALFA(1:nread)
endif


!-----------------------------------------------------------


if (iecho.eq.2) then
    if (nread.gt.1) then
        write(stderr,*) 'read ', nread, 'time step groups'
    else
        write(stderr,*) 'read ', nread, 'time step group'
    endif
    write(stderr,*) 'step count  : ', timestepdat%MAXSTP(1:nread)
    write(stderr,*) 'stepsizes   : ', timestepdat%DELT(1:nread)
    do i=1,nread
        write(stderr,*) 'unit string : ', UNIT(i)
    enddo
    write(stderr,*) 'alfa        : ', timestepdat%ALFA(1:nread)
endif

! check if there are no empty timestep groups.
! if so, maybe there are too many time step groups given
do i=2,nread
    if (timestepdat%MAXSTP(i) .eq. 0d0) then
         write(*,*) '********* WARNING *************'
         write(*,*) 'timestepgroup ', i, 'is empty'
         write(*,*) 'perhaps more timestepgroups given than actually set?'
    endif
enddo

do i=1,NINTG
    MAXITIME = MAXITIME + timestepdat%MAXSTP(i)
enddo

n = NMPRT

do i=n,1,-1
    if (outputcontroldat%IMPRINT(i).gt.MAXITIME) then
       NMPRT = NMPRT - 1
    endif
enddo

if (n.ne.NMPRT) then
    write(stderr,21) MAXITIME, n, NMPRT
  21        format(1x,'CONTRL warning: output is requested at time ', & 
 'steps larger than the total number of time steps ', I7, '.', &
 ' The number of outputs NMPRT has been reduced from ', I7, &
 ' to ',I7,'.')
endif
!
if (OUTFIL(1)) then
    write(luout(1),7)
 7     format(///1x,'T I M E   S T E P   I N F O R M A T I O N'// &
4x,'TIME STEP',6x,'# OF STEPS',7x,'STEP SIZE',13x,'ALFA'/ &
5x,'GROUP #',8x,'IN GROUP'//)

    do i=1,nread
        write(luout(1),8) i,timestepdat%MAXSTP(i), &
     timestepdat%DELT(i),UNIT(i), &
     timestepdat%ALFA(i)
 8              format(7x,i5,12x,i5,6x,1PG12.3,2x,a5,5x,0Pf4.2)
    enddo

endif
!

if (iecho.eq.8) then
    call parallellog("Contrl","halfway contrl")
endif

do i=1,NINTG

    n = fnblnk(UNIT(i))
    m = lnblk(UNIT(i))
    m = MIN0(m+1,5)

    if (n.lt.m) then
        tmp = UNIT(i)(n:m)
        if (tmp(1:2).eq.'se') then
            factor = 1d0
        else if (tmp.eq.'ms') then
            factor = 1D-3
        else if ((tmp.eq.'ye') .or. (tmp.eq.'year')) then
            factor = YEAR
        else if (tmp.eq.'Ma'.or.tmp.eq.'My') then
            factor = YEAR * 1D6
        else
            write(stderr,13) UNIT(i)
13             format(1x,'Unit of time not recognized: ',A5)
           factor = 1d0
        endif
    else
        factor = 1d0
    endif

    timestepdat%DELT(i) = timestepdat%DELT(i)*factor
enddo


if (iecho.eq.8) then
    call parallellog("Contrl","Leaving contrl")
endif

return


600   write(stderr,*) 'CONTRL: fatal input error when reading nr of timesteps in group. Error: ', ios
goto 1000
700   write(stderr,*) 'CONTRL: premature end when reading nr of timesteps in group. Error: ', ios
goto 1000

610   write(stderr,*) 'CONTRL: fatal input error when reading duration of timesteps. Error: ', ios
goto 1000
710   write(stderr,*) 'CONTRL: premature end when reading duration of timesteps. Error: ', ios
goto 1000

620   write(stderr,*) 'CONTRL: fatal input error when reading time unit of timesteps. Error: ', ios
goto 1000
720   write(stderr,*) 'CONTRL: premature end when reading time unit of timesteps. Error: ', ios
goto 1000

630   write(stderr,*) 'CONTRL: fatal input error when reading alpha value of timesteps. Error: ', ios
goto 1000
730   write(stderr,*) 'CONTRL: premature end when reading alpha value of timesteps. Error: ', ios
goto 1000



1000  call triostat(ios)
call prerr('CONTRL')
call xit(1," ")

end subroutine CONTRL

!-----routine to read time from input------------------------

subroutine RDTIME (ITIME,N)

use debugmodule, only: xit
use iomodule,           only: luin, stderr

implicit none

!-pass
integer N,NN,i,nread,ios
integer ITIME(MAX(N,1))
!-local
character(len=50000) record

nread = MAX(N,1)

read(luin,"(a)",err=100,end=110,iostat=ios) record

! if the number of time outputs in TECIN.DAT is positive,
! we have to print on the timesteps in the record.
! Otherwise we need to every every -NMPRT timesteps

if (NMPRT.gt.0) then

    read(record,3,err=100,end=110,iostat=ios) (ITIME(i),i=1,nread)
3 format(10000i5)

    if (MAXITIME.gt.0) then ! maxitime comes directly from the module
        NN = N
        do i=NN,1,-1
            if (ITIME(i).gt.MAXITIME) then
                N = N - 1
            endif
        enddo
    endif

    ! check whether output timestep are in the proper order:
    if (N.gt.1) then
        do i=1,N-1
            if (ITIME(i+1).lt.ITIME(i)) then
                write(*,*) "Output times should be in strictly ascending order"
                write(*,*) "Output time ",i,ITIME(i)," is higher than ",i+1,ITIME(i+1)
                stop "Exiting... "
            endif
        enddo
    endif

else
    ITIME(1) = 0
endif

return

100   write(stderr,101) ios
  101 format(1x,' rdtime: input error ',I8)
call prerr('RDTIME')
call xit(1," ")

110   write(stderr,111) 
  111 format(1x,' rdtime: premature EOF error')
call xit(1," ")

end subroutine RDTIME

!----some additional routines related to reading and plotting of time
subroutine PUTOUT (IPRINT,NPRINT)

use iomodule, only: luout, outfil
!-pass
integer NPRINT,i
integer IPRINT(*)

if ((.not.OUTFIL(1)).or.(NPRINT.le.0)) then
    return
endif

write(luout(1),2) (IPRINT(i),i=1,NPRINT)
 2 format(6x,10I6)

return
end subroutine PUTOUT
!-------------------------------------------------------------------------------
subroutine RDTSTP (TALF,NINTG)

use debugmodule, only: iecho, xit
use iomodule,           only: luin, luout, stderr, outfil

!-pass
integer NINTG, ios, nread,i
double precision TALF(*)
!
if (iecho.eq.4) then
    write(stderr,1)
 1     format(1x,'*** RDTSTP: read implicit parameters: ***'/ &
 1x,'nread = MAX(1,NINTG)'/ &
 1x,'read (TALF(i),i=1,nread) (50f5.0)')
endif
!
nread = MAX(NINTG,1)
read(luin,3,err=100,end=110,iostat=ios) (TALF(i),i=1,nread)
 3 format(50f5.0)

if (NINTG.gt.0 .and. OUTFIL(1)) then
     write(luout(1),4)
 4     format(/// &
  ' T H E R M A L   I M P L I C I T N E S S   P A R A M S'// &
  5X,'TIME STEP      ALFA'/ &
  5X,'  GROUP #'//)
    do i=1,NINTG
        write(luout(1),5) i,TALF(i)
 5         format(9X,I5,4X,1PG10.4)
    enddo
endif

return

100   write(stderr,101) ios
  101 format(1x,'RDTSTP: fatal input error ',i8)
call prerr('RDTSTP')
call xit(1," ")

110   write(stderr,111)
  111 format(1x,'RDTSTP: premature EOF error')
call xit(1," ")

end subroutine RDTSTP
!-------------------------------------------------------------------------------
subroutine RDMPRT (IMATPR,NPRMAT)

use debugmodule, only: xit
use iomodule,           only: luin
 ! read time step groups from tecin.dat

implicit none
!-pass
integer NPRMAT,i
integer IMATPR(*)


if (NPRMAT.le.0) then
    return
endif

read(luin,1,err=100,end=100) (IMATPR(i),i=1,NPRMAT)
 1 format(51I5)

return

100   call prerr('rdmprt')
call xit(1," ")

end subroutine RDMPRT


!------ TimDif ---------------------


subroutine timdif (TMECH,TTHERM)

use constants, only: zero

!
! returns the (real) time difference between current and previous thermal and
! mechanical output times, TTHERM and TMECH
!
implicit none

!-pass
double precision :: TTHERM
double precision :: TMECH
!-locl
logical :: Mout,Tout
integer :: itime, igrp
double precision :: time
integer :: i

!-init
TTHERM = ZERO
TMECH  = ZERO

time   = ZERO
itime  = 0

if (NSTEP.le.0) then
    goto 100
endif

do igrp=1,NINTG
    do i=1,timestepdatactx%MAXSTP(igrp)

        itime = itime + 1
!            NSTEP is time of current output
        time  = time  + timestepdatactx%DELT(igrp)

        if (itime.eq.NSTEP) then
            goto 100
        endif

        call chktim(itime,Mout,Tout)

        if (Mout) then
            TMECH  = time
        endif

        if (Tout) then
            TTHERM = time
        endif

    enddo
enddo


100 TMECH  = time - TMECH
TTHERM = time - TTHERM

return

end subroutine




!----- set timestep size------------------------------------

subroutine SetTimeStepSize ()

use modeldefinition, only: mode
use iomodule,        only: stderr

! Routine for setting time step size at current time

implicit none

!-local
integer :: igrp, LastStep

!-init

igrp   = 1
LastStep = 0

DELTP  = 0d0
ALFAP  = 0d0
TALFP  = 0d0

if (NSTEP.le.0) then
    return
endif

100 if (igrp.le.NINTG) then
    LastStep = LastStep + timestepdatactx%MAXSTP(igrp)
    if (NSTEP.le.LastStep) then
        goto 200
    endif
    igrp = igrp + 1
    goto 100
else
    write(stderr,1) NSTEP
1   format(1x,'setstp: NSTEP not found: ',I5)
    call exitp(1)
endif

200 DELTP = timestepdatactx%DELT(igrp)

if (MODE.eq.2 .or. MODE.eq.6) then
    ALFAP = timestepdatactx%ALFA(igrp)
endif

if (MODE.eq.5 .or. MODE.eq.6) then
    TALFP = timestepdatactx%TALF(igrp)
endif

return
end subroutine

!-----check whether this timestep requires output (fedsk) files------------------------------------
subroutine chktim (IREQ,Mout,Tout)

USE MODELCTX
USE MODELDEFINITION, only: MODE

! Check if requested time is in output.

!-pass
integer :: IREQ
logical :: Mout,Tout
!-locl
logical :: Mech,Temp
integer :: LastStep
integer :: im, it, i

!-init
im   = 0
it   = 0

Mech = .false.
Temp = .false.

 ! elastic output is always provided.
if (IREQ.eq.0) then
    if (MODE.eq.2 .or. MODE.eq.6) then
        Mech = .true.
    endif
    if (MODE.eq.5 .or. MODE.eq.6) then
        Temp = .true.
    endif
    goto 600
endif

 ! NINTG = number of time step groups
if (NINTG.le.0) then
    goto 600
endif

LastStep = 0

do i=1,NINTG
   LastStep = LastStep + timestepdatactx%MAXSTP(i)
enddo

if (IREQ.gt.LastStep) then
   goto 600
endif

if (MODE.ne.2 .and. MODE.ne.6) then
    goto 300
endif

if (NMPRT.gt.0) then
    i = 1
100 if (i.le.NMPRT) then
        if (outputcontroldatactx%IMPRINT(i).eq.IREQ) then
            im = i
            goto 200
        endif
        i = i + 1
        goto 100
    endif
    goto 300
200        Mech = .true.
else if (NMPRT.lt.0) then
    if (MOD(IREQ,IABS(NMPRT)).eq.0) then
        Mech = .true.
    endif
endif

300 i = 1

if (MODE.ne.5 .and. MODE.ne.6) then
    goto 600
endif

400 if (NTPRT.gt.0) then
    if (i.le.NTPRT) then
         if (outputcontroldatactx%ITPRINT(i).eq.IREQ) then
             it = i
             goto 500
         endif
         i = i + 1
         goto 400
    endif
    goto 600
500     Temp = .true.
else if (NTPRT.lt.0) then
    if (MOD(IREQ,IABS(NTPRT)).eq.0) then
        Temp = .true.
    endif
endif

600 if (Temp .or. Mech) then
    IMPRT = im
    ITPRT = it
endif

Mout  = Mech
Tout  = Temp

return
end subroutine



!-----check on thermal output timesteps----------------------------------------------
subroutine CHKOUT(timestepdat,outputcontroldat)

use debugmodule, only: xit
use iomodule,    only: stderr

! Routine checking if requested output times match calculation times and
! if times are in ascending order.

implicit none
!-pass
type(timestepdata) :: timestepdat
type(outputcontroldata) :: outputcontroldat

integer i,j,ISTEP,MTPRT,MTCALC

if (NINTG.le.0) then
    return
endif

if (NTCALC.gt.1) then
!        check ascending order
    do i=2,NTCALC
           if (timestepdat%ICALC(i-1).ge. &
               timestepdat%ICALC(i)) then
           write(stderr,10)
10             format(1x,'ICALC not ascending')
           call xit(1," ")
       endif
    enddo
endif

if (NTPRT.gt.1) then
!        check ascending order
    do i=2,NTPRT
        if (outputcontroldat%ITPRINT(i-1).ge. &
            outputcontroldat%ITPRINT(i)) then
            write(stderr,20)
20             format(1x,'ITPRINT not in ascending order')
            call xit(1," ")
        endif
    enddo
endif

if (NTCALC.eq.0 .and. NTPRT.ne.0) then
    write(stderr,25)
25     format(1x,'Number of thermal calculations = 0, so',1x, &
       'NTPRT should be = 0 too')
    call xit(1," ")
endif

if (NTCALC.gt.0 .and. NTPRT.gt.0) then
!        check if all ITPRINT numbers occur in ICALC also
    j = 1
    do i=1,NTPRT
300     if (j.gt.NTCALC) then
            write(stderr,30) i,outputcontroldat%ITPRINT(i)
30             format(1x,'ITPRINT(',I2,')=',I3, &
            ' does not occur in ICALC')
           call xit(1," ")
        endif
        ISTEP = outputcontroldat%ITPRINT(i)
        if (ISTEP.eq.timestepdat%ICALC(j)) then
            cycle
        endif
        j = j + 1
        goto 300
    enddo

else if (NTCALC.ge.0 .and. NTPRT.lt.0) then

    MTPRT = INT(REAL(timestepdat%ICALC(NTCALC))/ &
                REAL(-NTPRT))
    if (MTPRT.lt.1) then
        write(stderr,35) -NTPRT
35         format(1x,'No temperature calculation every ',I4,1x, &
      'time steps')
       call xit(1," ")
    endif
    j = 1
    do i=1,MTPRT
500     if (j.gt.NTCALC) then
            write(stderr,40) i
40             format(1x,'ITPRINT(',I2,')', &
         ' does not occur in ICALC')
           call xit(1," ")
        endif
        ISTEP = i * (-NTPRT)
        if (ISTEP.eq.timestepdat%ICALC(j)) then
            cycle
        endif
        j = j + 1
        goto 500
    enddo

else if (NTCALC.lt.0 .and. NTPRT.gt.0) then

    MTCALC = INT(REAL(outputcontroldat%ITPRINT( &
                 NTPRT))/ &
            REAL(-NTCALC))
    j = 1
    do i=1,NTPRT
700     if (j.gt.MTCALC) then
            write(stderr,30) i,outputcontroldat%ITPRINT(i)
            call xit(1," ")
        endif
        ISTEP = j * (-NTCALC)
        if (ISTEP.eq.outputcontroldat%ITPRINT(i)) then
            cycle
        endif
        j = j + 1
        goto 700
    enddo

else if (NTCALC.lt.0 .and. NTPRT.lt.0) then
    if (MOD(-NTPRT, &
       -NTCALC).ne.0) then
        write(stderr,50)
50         format(1x,'NTPRT and NTCALC do not match')
        call xit(1," ")
    endif
endif

return
end subroutine CHKOUT
!-------------------------------------



logical function output (IARRAY,NPRINT,IPRINT,iSTEP)

! determines where there should be output at this level.
! By determined whether istep occurs in iarray

implicit none

integer NPRINT,IPRINT,iSTEP
integer IARRAY(*)
integer i

output = .false.

if (NPRINT.eq.0 .or. IPRINT.eq.0) then
    return
endif

do i=1,NPRINT
    if (IARRAY(i).eq.iSTEP) then
        output = .true.
        return
    endif
enddo

!check for negative output, which should give a true every NPRINT timesteps
if (NPRINT.lt.0) then
    if (mod(iStep, -nPrint).eq.0) then
        output = .true.
        return
    endif
endif


return
end function output


END MODULE TIMESTEPMODULE
