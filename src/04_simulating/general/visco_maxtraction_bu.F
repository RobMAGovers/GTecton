subroutine MMARCH (IDMAT,NSD,NDOF,NSTR,NEN,NSTEP)

USE MESHDATAMODULE,  only:  meshdatactx
USE MODELDATAMODULE, only: modeldatactx, &
                           IDIAG
USE MODELDEFINITION
USE MATERIALSMODULE, only:   matdatactx, PRPSED
USE TIMESTEPMODULE
USE MODELCTX
USE ALGEBRA
use constants,       only: eps
use filemodule,      only: workpath, &
                           outdir
use debugmodule,     only: lustat, &
                           itest, &
                           debug, &
                           petsctime, &
                           owntime, &
                           time_start, &
                           time_finish, &
                           rate, &
                           iecho, &
                           LogTimeStamp, &
                           parallelLog, &
                           checkArrayForNans, &
                           delay, &
                           iAmAtTimestep, & 
                           isFileOpen, &
                           NicolaiNodeID, &
                           NicolaiFileID
use petscksp
use iomodule          

#ifdef SPARSE
USE AOMODULE
#endif

!-------------------------------------------------------------------------------
! Program to solve the time dependent mechanical problem.
!-------------------------------------------------------------------------------

implicit none
!***** pass
integer NSD,NDOF,NSTR,NEN
integer :: IDMAT(NSTR)

#ifdef SPARSE
#include "petsc/finclude/petscsys.h"
  PetscErrorCode ierr
#else
  integer :: ierr
#endif

!**** local variables
logical DoOut,PrtMat
#ifdef checkLoadvecForNans
logical :: loadVecHasNans
#endif
external lnblk
integer :: inneri, innerj, logsize, rnk, i, lnblk, iVecEntry
double precision :: ACCURCY, grav(ndof), psumsq, factor, YFORCE
save YFORCE

logical :: iterate
integer :: ITERNUMBER
double precision :: tempDELTP
#ifdef EXPERIMENTAL_FRICTION
  integer :: k,j
  double precision :: CONVERGENCE_VALUE
#endif
character(len=255) residufile
character(len=5) chartime

!! init
#ifdef EXPERIMENTAL_FRICTION
 data YFORCE/0d0/ 
 CONVERGENCE_VALUE=1.0D-1
#endif

if (MODE.ne.2 .and. MODE.ne.6) then
    return
endif

!***** whether or not to write a fedsk file
if (NMPRT.lt.0) then
    DoOut = (mod(NSTEP,-NMPRT).eq.0)
else
    DoOut = output(outputcontroldatactx%IMPRINT, NMPRT,IMPRT,NSTEP)
endif

! whether or not to write matrix data to output
PrtMat = output (outputcontroldatactx%IMATPR,  NPRMAT,IPRMAT,NSTEP)

if (itest.eq.2) then
    DoOut = .true.
    PrtMat = .true.
endif

!**** verify time step size against stability criteria
if (DoOut) then
    call StpChk (modeldatactx%STN,meshdatactx%mat, matdatactx%PRPMAT,NSTEP)
endif

! Check whether Winkler forces are being imposed or removed
if (NWINK .gt.0) then
    call CKWINK(modeldatactx%IWINK,NSTEP,NDOF*NUMNP,IMWORK)
endif

if (NWINKX.gt.0) then
    call CKWINX(modeldatactx%IWINX,NSTEP,NDOF,NUMNP,IMWORK)
endif

if (NUMWNK.ne.0) then
    call CKWINK(modeldatactx%IWTIME,NSTEP,IABS(NUMWNK),IMWORK)
endif

if (iecho.eq.6) then
    call LogTimeStamp("Checked Winkler")
endif

! Update skew angles for slippery displacements
if (NSLSKEW.gt.0) then
    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "Updating angles")
    endif
    call EQUATE (modeldatactx%OSKEW,modeldatactx%SKEW,2*NUMNP)

    call SlipperySkew(modeldatactx%IDSLP,modeldatactx%SKEW, &
                      modeldatactx%NSELSD,meshdatactx%IEN,meshdatactx%X, &
                      modeldatactx%D,modeldatactx%TFAULT,modeldatactx%DXE, &
                      modeldatactx%LMF)
    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "Finished updating angles")
    endif
    if (iecho.eq.6) then
        call LogTimeStamp("Done fault parallel surfaces")
    endif
endif

if (IMWORK.ne.0) then
    if (showfbc) then
        call CLEAR(modeldatactx%S,MEQ*NDOF*NUMNP)
    endif

    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "Impliciting")
    endif

    ! Lagrangian (or requested) update of stiffness matrix
    call IMPLCT(modeldatactx%A,IDIAG,modeldatactx%B,modeldatactx%ID, &
                modeldatactx%IDX,meshdatactx%X,modeldatactx%D, &
                modeldatactx%SKEW,modeldatactx%IWINK, &
                modeldatactx%WINK,modeldatactx%IWINX,modeldatactx%WINX, &
                meshdatactx%IEN,modeldatactx%LM,modeldatactx%LMX, &
                meshdatactx%mat,matdatactx%PRPMAT, &
                modeldatactx%DMAT,modeldatactx%STN,IDMAT,modeldatactx%LMF, &
                modeldatactx%TFAULT,matdatactx%prplas,modeldatactx%IWELM, &
                modeldatactx%IWSIDE,modeldatactx%IWTIME,modeldatactx%WPRES, &
                modeldatactx%DXE,modeldatactx%NSLIP,modeldatactx%S, &
                modeldatactx%forShowfbcIAeqnr, &
                modeldatactx%forShowfbcIAXeqnr, &
                modeldatactx%BETB,NSD,NDOF,NSTR,NEN,NSTEP)
    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "Done impliciting")
    endif
    if (iecho.eq.6) then
        call LogTimeStamp("Done implicit")
    endif

endif

!****  Print matrix diagonal if required
if (PrtMat) then
    ! printp only works for the direct solver.
    call PRINTP (modeldatactx%A,IDIAG,modeldatactx%ID, &
    modeldatactx%IDX,NDOF,0)
endif

if (IRESDU.gt.0 .or. DoOut) then
!*****    Calculate deviations from global force balance in previous step
!*****    Overwrites B and BRES

    call CLEAR (modeldatactx%B,NEQlocal,"modeldatactx%B")
    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "calling residu")
    endif

    !***** compute load vector B equivalent to stresses and ...
    !***** .... residu, BRES=BTOT-B. 
    call RESIDU (PSUMSQ,ACCURCY)

    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "Finished residu")
    else if (iecho.eq.6) then
        call LogTimeStamp("Done Residual force update")
    endif

    if (optijd) then
        !**** optijd=.true. immediately after restart (first MMARCH call). The equivalent load vector of 
        !**** stresses is set here: BTOT = B. Implicitly, this means that the force residu before 
        !**** restart was zero. Therefore, set the residual load vector BRES=0
        call EQUATE(modeldatactx%BTOT,modeldatactx%B,NEQlocal)
        call CLEAR(modeldatactx%BRES,NEQlocal,"modeldatactx%BRES")
        optijd=.false.
    endif

    if (DoOut) then

        write(chartime,'(I5.5)') NSTEP
        residufile = 'residu.dat.'//chartime//'.'

        rnk = getrank()
        logsize = log10(dble(getsize())) + 1
        do i=1,logsize
            inneri = 10**(logsize-i)
            innerj = rnk/inneri
            residufile = residufile(1:lnblk(residufile))//char(48+innerj)
            rnk = rnk - innerj * inneri
        enddo

        residufile = trim(outdir)//'/'//trim(residufile)

        call WRBRES (modeldatactx%BRES,meshdatactx%X, &
                     modeldatactx%D,modeldatactx%SKEW, &
                     meshdatactx%IEN,modeldatactx%LM,modeldatactx%LMX, &
                     modeldatactx%LMF,modeldatactx%TFAULT, &
                     modeldatactx%DXE,NSD,NDOF,NEN,residufile)!, 'residux.dat', 'residuy.dat')

        if (iecho.eq.6) then
            call LogTimeStamp("Did out")
        endif
    endif

    if (DoOut) then
        write(stderr,4) NSTEP-1,ACCURCY*1e2,NSTEP
4       format(1x,'Residual force from time step ',I5,' is ',F7.1,' %. ',&
         'You requested output after the upcoming time step (',I5,')')
        if (ACCURCY.gt.5d-2) then
            if (IRESDU.ne.2) then
                write(stderr,6) IRESDU
6               format('*** WARNING ***'/'*** The solution will not be balanced/accurate because ',&
                 'the residual force update is off (IRESDU=',I1,')')
            else
                write(stderr,8) 
8               format('*** WARNING ***'/'*** The residual force will be updated, but given that ',&
                 'it is greater than 5%,'/'*** it is likely that the solution will not be balanced/accurate.'/&
                 '*** The residual force file "residu.dat..." has more information')
            endif
        endif
    else
        if (iecho.eq.1) then
            if (IRESDU.eq.1) then
                write(stdout,10) NSTEP-1,ACCURCY*1d2
10              format(1x,'Residual force from time step ',I5,' = ',F7.1, ' %')
            else if (IRESDU.eq.2) then
                write(stdout,11) NSTEP-1,ACCURCY*1d2
11              format(1x,'Residual force from time step ',I5,' = ',F7.1, ' % (updated in next step)')
            endif
        endif
    endif

endif ! of calculation of residual force 

!*****    Clear the global load increment array B=0
!#ifdef SPARSE
     modeldatactx%B = 0d0
!#else
!         call CLEAR (modeldatactx%B,NEQlocal,"modeldatactx%B")
!#endif /* SPARSE */

#ifdef SHEARZONE_TRACTIONS
! Add velocity-dependent shear zone tractions
if (NUMELSHEARZONE.gt.0) then ! If there are any shear zone entries applying to local elements
    ! if (DELTP.gt.0.0) then          ! Apply slip rate--dependent tractions only if the timestep size is not 0
    !                                 ! (which would yield an infinite slip rate)
        if (NTYPE.eq.1) then        ! If using pln
#ifdef EXPERIMENTAL_SHEARTRAX_DAMPING
            call SHEARTRAX(modeldatactx%B,modeldatactx%shearF,modeldatactx%BTOT,meshdatactx%X,modeldatactx%D,modeldatactx%IDX, &
             meshdatactx%IEN, modeldatactx%TFAULT,modeldatactx%NSLIP,modeldatactx%DX,DELTP,modeldatactx%shearETA, &
             modeldatactx%shearTHICK,modeldatactx%ELNOSHEARZONE,matdatactx%PRPMAT,modeldatactx%shearnpl, & 
             modeldatactx%shearQpl,modeldatactx%shearTEMP)
#else
            call SHEARTRAX(modeldatactx%B,modeldatactx%BTOT,meshdatactx%X,modeldatactx%D,modeldatactx%IDX, &
             meshdatactx%IEN, modeldatactx%TFAULT,modeldatactx%NSLIP,modeldatactx%DX,DELTP,modeldatactx%shearETA, &
             modeldatactx%shearTHICK,modeldatactx%ELNOSHEARZONE,matdatactx%PRPMAT,modeldatactx%shearnpl, & 
             modeldatactx%shearQpl,modeldatactx%shearTEMP)
#endif
        else if (NTYPE.eq.6) then   ! If using f3d
            call SHEARTRAX(modeldatactx%B,modeldatactx%BTOT,meshdatactx%X,modeldatactx%D,modeldatactx%IDX,meshdatactx%IEN,modeldatactx%TFAULT,&      
             modeldatactx%NSLIP,modeldatactx%DX,DELTP,modeldatactx%shearETA,modeldatactx%shearTHICK,&
             modeldatactx%ELNOSHEARZONE,modeldatactx%shearnpl,modeldatactx%shearQpl,modeldatactx%shearTEMP)
        endif                  
        modeldatactx%shearF = modeldatactx%B ! Save B vector to later add to BTOT
    ! endif
endif
#endif /* SHEARZONE_TRACTIONS */


!**** Clear displacement increment: DELD=0
modeldatactx%DELD = 0d0

!**** Clear faulted node displacement increment: DFAULT=0
if (NUMFN.gt.0) then
    modeldatactx%DFAULT = 0d0
endif

!**** Clear slippery node displacement increment: DX=0
if (NUMSN.gt.0) then
    modeldatactx%DX = 0d0
endif

!***** add pseudo-forces from nodal boundary conditions
if (iecho.eq.8) then
    call ParallelLog("visco", "calling load")
endif

call LOAD (modeldatactx%ID,modeldatactx%IBOND,modeldatactx%BOND, &
           modeldatactx%D, modeldatactx%DELD, modeldatactx%B, &
           meshdatactx%X,NUMNP,NEQlocal,1)


if (iecho.eq.8) then
    call ParallelLog("visco", "left load")
else if (iecho.eq.6) then
    call LogTimeStamp("Done load")
endif

if (DoOut .and. debug) then
    write(luout(1),20)
20  format(//1x,'LOAD VECTOR NODAL BC''S'/)
    call PRINTB (modeldatactx%B,modeldatactx%ID, &
                 modeldatactx%IDX,NDOF)
endif

!*****    Update faulted node displacements
if (iecho.eq.8) then
    call ParallelLog("visco", "entering loadf")
endif

call LOADF (meshdatactx%IEN,modeldatactx%NFAULT, modeldatactx%FAULT, &
            modeldatactx%DFAULT,NUMEL)


if (iecho.eq.8) then
    call ParallelLog("visco", "left loadf")
elseif (iecho.eq.6) then
    call LogTimeStamp("Done loadf ")
endif

!*****    If thermo-mechanical calculation, update global load vector with
!*****    body force contributions from thermal expansion.
if (IMTWRK.ne.0 .and. ITPRES.eq.0) then
    if (iecho.eq.1) then
        write(stdout,30)
30      format(1x,'Loads due to thermal expansion are being added')
    endif

    call TBDYLD(meshdatactx%mat,meshdatactx%IEN,meshdatactx%X, &
                modeldatactx%D,modeldatactx%LMF,modeldatactx%TFAULT, &
                modeldatactx%LMX, &
                modeldatactx%TLST,modeldatactx%T,modeldatactx%TX, &
                modeldatactx%TLSTX,matdatactx%PRPMAT,matdatactx%PRPTEM,GRAV, &
                modeldatactx%B,modeldatactx%LM,modeldatactx%SKEW, &
                modeldatactx%LMTX,modeldatactx%DXE)

    if (iecho.eq.6) then
        call LogTimeStamp("Done tbdyld ")
    endif

    if (DoOut .and. debug) then
        write(luout(1),40)
40         format(//1x,'LOAD VECTOR AFTER THERMAL LOADS'/)
        call PRINTB (modeldatactx%B,modeldatactx%ID, &
                   modeldatactx%IDX,NDOF)
    endif
endif

!******************************************************************
!***** OBSOLETE CALL TO MPI_BARRIER ?????
!******************************************************************
#ifdef SPARSE
call MPI_barrier(MPI_COMM_WORLD, ierr)
#endif

! MPI_BARRIER because not all threads may have differential winkler nodes.
! In that case, the threads without do enter the if statement below,
! and will be too fast ad adding the stress below, while
! the stress caused by the diff. Winkler forces is not yet known.

!write(*,*) "rank", getrank(), "has nwinkx", nwinkx
if (NWINKX.gt.0) then
    ! that is local nwinx. So some partitions may go through here, 
    ! while others skip it

! write(*,*) "BTOT @ pre winkx unlock", modeldatactx%BTOT
    call UNLOCK()
! write(*,*) "BTOT @ post winkx unlock", modeldatactx%BTOT
endif

!******************************************************************
!***** OBSOLETE CALL TO MPI_BARRIER ?????
!******************************************************************
#ifdef SPARSE
call MPI_barrier(MPI_COMM_WORLD, ierr)
#endif

!write(*,*) "rank", getrank(), "has finished unlocking loop"


!*****    Add time-dependent stress forces, if present, to force vector
if (NUMSTR.ne.0) then
      call ADDSTR (modeldatactx%B, &
                meshdatactx%X, &
               modeldatactx%D, &
                meshdatactx%IEN, &
               modeldatactx%LM, &
                meshdatactx%mat, &
                 matdatactx%PRPMAT, &
               modeldatactx%IELSTR, &
               modeldatactx%ISSIDE, &
               modeldatactx%STRS, &   ! stress boundary condition, not stress during the run, which is STN
               modeldatactx%ISTR, &
               modeldatactx%LMF, &
               modeldatactx%TFAULT, &
               modeldatactx%SKEW, &
               modeldatactx%DXE)
    if (iecho.eq.6) then
        call LogTimeStamp("Done addstr")
    endif

    if (DoOut .and. debug) then
        write(luout(1),50)
50         format(//1x,'LOAD VECTOR AFTER STRESS LOADS'/)
        call PRINTB (modeldatactx%B,modeldatactx%ID, &
                   modeldatactx%IDX,NDOF)
    endif
endif

!*****    Add time-dependent prestress forces to global load vector
if (NPRE.gt.0) then
    call PRESTN (modeldatactx%B,modeldatactx%STN0, &
                 modeldatactx%ISELM,modeldatactx%ISTIME,meshdatactx%X, &
                 modeldatactx%D,modeldatactx%SKEW,meshdatactx%IEN, &
                 modeldatactx%LM,modeldatactx%LMX,meshdatactx%mat, &
                 matdatactx%PRPMAT,modeldatactx%LMF,modeldatactx%TFAULT, &
                 modeldatactx%DXE,NSD,NSTR,NEN,NDOF)
    if (iecho.eq.6) then
        call LogTimeStamp("Done prestress")
    endif
endif

!*****    Add loads from sediment redistribution
call sdimnt (modeldatactx%ISURF,modeldatactx%SURF,modeldatactx%B, &
             modeldatactx%ID,GRAV,PRPSED,DELTP)

!*****    Add Maxwell boundary conditions to global load
if (NMAXW.gt.0) then
    call MAXWEL (modeldatactx%FMAXW,modeldatactx%FTOT, &
                 modeldatactx%B,modeldatactx%ID,NDOF)
endif

!    Update the applied forces (excluding pseudoforces due to displacement
!    boundary conditions, split nodes, viscous flow or plastic yield).
#ifdef SPARSE

!******************************************************************
!***** Some partitions may have prestress and stress,
!***** while other partitions may not.
!***** Make sure that the threads of those that do not can catch up
!***** with the ones that don't
!******************************************************************
call MPI_barrier(MPI_COMM_WORLD, ierr)
call SyncLoadVector(0)

!      write(FILE_outputf,*) 'B after sync'
!      do i=1,NEQext
!           write(FILE_outputf,*) i, equations(i), modeldatactx%B(i)
!      enddo
#endif

!***** add B to BTOT
call UPDATE (modeldatactx%BTOT,modeldatactx%B,NEQlocal)
if (iecho.eq.6) then
    call LogTimeStamp("Done update")
endif

!*****    Update global force vector with residuals if required
if (IRESDU.gt.1) then
    call UPDATE (modeldatactx%B, modeldatactx%BRES,NEQlocal)
endif

#ifdef SPARSE
    call CLEAR(modeldatactx%BRES,nequations,"modeldatactx%BRES")
#else
    call CLEAR (modeldatactx%BRES,NEQlocal,"modeldatactx%BRES")
#endif 

if (DoOut .and. debug) then
    write(luout(1),55)
55  format(//1x,'LOAD VECTOR AFTER RESIDUAL LOADS'/)
    call PRINTB (modeldatactx%B,modeldatactx%ID, &
                 modeldatactx%IDX,NDOF)
endif

! **** Update incremental load vector B with contribution from visco-plastic pseudoforces
! **** and from velocity boundary conditions
if (iecho.eq.8) then
    call ParallelLog("visco", "calling formef")
endif
call FORMEF (modeldatactx%B,meshdatactx%X,modeldatactx%D, &
             modeldatactx%DELD,modeldatactx%SKEW,meshdatactx%IEN, &
             modeldatactx%LM,modeldatactx%LMX, &
             meshdatactx%mat,modeldatactx%BETA,modeldatactx%BETB, &
             matdatactx%PRPMAT,matdatactx%prplas,modeldatactx%DMAT, &
             modeldatactx%STN,modeldatactx%LMF,modeldatactx%TFAULT, &
             modeldatactx%DXE,NSD,NDOF,NSTR,NEN)

if (iecho.eq.8) then
    call ParallelLog("visco", "left formef")
else if (iecho.eq.6) then
    call LogTimeStamp("Done formef")
endif

if (DoOut.and. debug) then
    write(luout(1),60)
60  format(//1x,'LOAD VECTOR AFTER ADDING LOADS FROM VISCO-ELASTIC AND/OR ', &
     'VISCO-PLASTIC RELAXATION, AND FROM VELOCITY BCS'/)
    call PRINTB (modeldatactx%B,modeldatactx%ID, &
                 modeldatactx%IDX,NDOF)
endif

!*****    Add pseudoforces from time dependent split node displacements.
!*****    DFAULT has split node displacement increments during current time step
if (NUMFN.gt.0) then
    if (iecho.eq.1) then
        write(stdout,65)
65      format(1x,'Split node loads are being added'/)
    endif

    call FORMF (modeldatactx%B,meshdatactx%X,modeldatactx%SKEW, &
        meshdatactx%IEN,modeldatactx%LM,modeldatactx%D, &
        modeldatactx%LMX,meshdatactx%mat,matdatactx%PRPMAT, &
        modeldatactx%DMAT,modeldatactx%LMF,modeldatactx%NFAULT, &
        modeldatactx%DFAULT,modeldatactx%TFAULT, &
        modeldatactx%DXE,NSD,NDOF,NSTR,NEN)
endif

#ifdef SPARSE
if (IMWORK.ne.0) then
  call BuildMat(0)
endif
#endif

ITERNUMBER=0
iterate = .true. ! single step is also called "iteration"
tempDELTP=DELTP  ! store time step size 

#ifdef EXPERIMENTAL_FRICTION

if (nFaultSegmentEntriesGlobal.gt.0) then
    !***** in the loop below DELD and DX are the iterative displacement updates
    if (debug) write(*,*) 'mmarch: zeroing cumulative Disp and DX arrays'
    modeldatactx%cumulativeDisp = 0d0  ! equivalent to CLEAR subroutine (also works in parallel)
    modeldatactx%cumulativeDX = 0d0
endif

#endif /* EXPERIMENTAL_FRICTION */


!*******************************************************************
!*****  start friction iterations block
!
! Note that pln also enters this loop if there is no friction!
! Then it just does a single iteration because the boolean
! 'iterate' is set to false at the appropriate place below:
!
!*******************************************************************


do while (iterate .and. ITERNUMBER.lt.1000)

    if (iecho.eq.8) then
        call parallellog("mmarch","starting friction loop")
    endif

    ITERNUMBER = ITERNUMBER + 1
    if (debug) write(*,*) 'mmarch: starting iteration ',ITERNUMBER
    if (ITERNUMBER.eq.1000) then
        write(0,*) 'NSTEP=',NSTEP,': number of iterations = 1000'
        call stoper()
    endif

    !***** Solve linear system of equations and return result in B. BRES used for temp storage
    call SOLVE (modeldatactx%A,modeldatactx%B,modeldatactx%BRES,IDIAG, &
                modeldatactx%ID,modeldatactx%IDX, &
                NUMNP,NEQlocal,NUMSN,NDOF,.FALSE.,.TRUE.)

    if (iecho.eq.6) then
        call LogTimeStamp("Done solve")
    endif

    !***** Copy displacement increment from solution vector B to DELD
    call DISP (modeldatactx%B,modeldatactx%DELD,modeldatactx%ID, &
               modeldatactx%LINK,NDOF,NUMNP,NLINK)


    !***** Copy slippery displacement increment from solution vector B to DX
    if (NUMSLP.NE.0) then
        call DISP (modeldatactx%B, modeldatactx%DX,modeldatactx%IDX, &
                   modeldatactx%LINK,NDOF,NUMNP,0)

        if (iecho.eq.6) then
            call LogTimeStamp("Done disp")
        endif
    endif

#ifdef SPARSE
    call scatterNDOFdata(meshdatactx,modeldatactx, &
                         getrank(),modeldatactx%DELD,FILE_outputf,0,1)

    if (NUMSLPglobal.ne.0) then
        call scatterNDOFdata(meshdatactx,modeldatactx, &
                             getrank(),modeldatactx%DX,FILE_outputf,0,2)
    endif

    if (iecho.eq.6) then
        call LogTimeStamp("Done scattering")
    endif
#endif /* SPARSE */

    if (NWINK.gt.0) then
        call WINKLF(modeldatactx%BTOT, &
                    modeldatactx%DELD,modeldatactx%ID, &
                    modeldatactx%IWINK,modeldatactx%WINK,NDOF,NUMNP,0)
        if (iecho.eq.6) then
            call LogTimeStamp("Done winklf")
        endif
    endif

    if (NWINKX.gt.0) then
        call WINKLF(modeldatactx%BTOT, &
                    modeldatactx%DX,modeldatactx%IDX, & 
                    modeldatactx%IWINX,modeldatactx%WINX,NDOF,NUMNP,1)
        if (iecho.eq.6) then
            call LogTimeStamp("Done nwinkx")
        endif
    endif

    if (NUMWNK.ne.0) then
        call WINKLP (modeldatactx%BTOT,modeldatactx%DELD, &
                     modeldatactx%IWELM,modeldatactx%IWSIDE,modeldatactx%IWTIME, &
                     modeldatactx%WPRES,meshdatactx%IEN,meshdatactx%mat, &
                     meshdatactx%X,modeldatactx%D,matdatactx%PRPMAT, &
                     modeldatactx%LM,modeldatactx%LMF,modeldatactx%TFAULT, &
                     modeldatactx%SKEW,modeldatactx%DXE,modeldatactx%LMX, &
                     modeldatactx%NSLIP,modeldatactx%DX)
        if (iecho.eq.6) then
            call LogTimeStamp("Done NUMWNK")
        endif
    endif

#ifdef SPARSE
    call SyncLoadVector(-1)    ! BTOT scatter
    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "Finished syncloadvec")
    endif
#endif

    !*****    Update Maxwell load vector
    if (NMAXW.gt.0) then
        if (iecho.eq.8) then
            call ParallelLog("Mmarch", "calling MXWLD")
        endif
    
        call MXWLD (modeldatactx%IWINK,modeldatactx%WINK, &
                    modeldatactx%FTOT,modeldatactx%FMAXW, &
                    modeldatactx%DELD,NDOF)
        if (iecho.eq.8) then
            call ParallelLog("Mmarch", "Finished MXWLD")
        endif
    endif

    if (showfbc) then
        call getneqbc (modeldatactx%IBOND,modeldatactx%S, &
                       modeldatactx%forShowfbcIAeqnr, &
                       modeldatactx%DELD, &
                       modeldatactx%BOND, &
                       NDOF,NUMNP,MEQ)
        if (MOD(NSTEP,1).eq.0) then
            call prtneqbc (modeldatactx%IBOND, &
                           modeldatactx%forShowfbcIAeqnr, &
                           modeldatactx%BOND,NDOF,NUMNP,NSTEP)
        endif
    endif

#ifdef EXPERIMENTAL_FRICTION

    if (nFaultSegmentEntriesGlobal.gt.0) then
        !***** Update displacement increments during this time step with the iteration displacement increment
        call UPDATE (modeldatactx%cumulativeDisp,modeldatactx%DELD,NDOF*meshdatactx%nvlocal)
        call UPDATE (modeldatactx%cumulativeDX,modeldatactx%DX,NDOF*meshdatactx%nvlocal)

        !***** clear incremental load vector B
        call CLEAR (modeldatactx%B,NEQlocal,"modeldatactx%B")

        !***** Update DIFORC from -DX*WINX and ...
        !***** ... compare differential Winkler forces with strength and ...
        !***** ... add overshoot to B-vector.
        call MAXTRACTION (modeldatactx%DX,modeldatactx%IDX,modeldatactx%IWINX,&
                          modeldatactx%WINX,modeldatactx%DIFFSTRENGTH,modeldatactx%B,&
                          modeldatactx%DIFORC,NDOF,NUMNP,NEQlocal,ITERNUMBER,&
                          YFORCE)

!***********************************************
!! THIS CANNOT WORK IN PARALLEL WITHOUT SCATTER OF B
!! see function globalSumSquare in algebra module that
!! should be able to handle this.
!***********************************************
        !***** Convergence? 
        PSUMSQ = DOT(modeldatactx%B,modeldatactx%B,NEQlocal)
        if (PSUMSQ.GT.0d0) then
            PSUMSQ = SQRT(PSUMSQ)
            if (iecho.eq.1) then
                if (YFORCE.gt.1d-16) write(*,*) 'NSTEP=',NSTEP,' traction overshoot = ',&
                                     PSUMSQ/YFORCE*100.,'%'
            endif
        endif
        iterate = (PSUMSQ.gt.CONVERGENCE_VALUE)
!******************
iterate = .false.
!***********************
    else
        iterate = .false.
    endif
#else
    iterate = .false.   ! no iteration actually ...
#endif /* EXPERIMENTAL_FRICTION */

    !**** Rotate displacement increments into global coordinate system
    call RDISP(modeldatactx%DELD,modeldatactx%SKEW,NDOF,NUMNP)

    if (NUMSLP.gt.0) then
        !**** Rotate differential displacement increments into global coordinate system
        call RDISP(modeldatactx%DX,modeldatactx%SKEW,NDOF,NUMNP)
    endif

    !**** Rotate split displacement increments into global coordinate system
    call ROTFLT(modeldatactx%NFAULT,modeldatactx%DFAULT, &
                modeldatactx%SKEW,NDOF,NUMFN)

#ifdef SPARSE
!    if (iecho.eq.8) then
!        call ParallelLog("Mmarch", "Scattering DELD and DX")
!    endif
     call scatterNDOFdata(meshdatactx,modeldatactx, &
                          getrank(),modeldatactx%DELD,FILE_outputf,0,3)

     if (NUMSLPglobal.ne.0) then
           call scatterNDOFdata(meshdatactx,modeldatactx, &
                                getrank(),modeldatactx%DX,FILE_outputf,0,4)
     endif

#endif /* SPARSE */

     !***** Update stresses
     call STRESN (meshdatactx%X,modeldatactx%D,modeldatactx%DELD, &
                  modeldatactx%DMAT,modeldatactx%BETA,modeldatactx%STN, &
                  meshdatactx%IEN,modeldatactx%LMF,modeldatactx%TFAULT, &
                  modeldatactx%DFAULT,modeldatactx%LMX, &
                  modeldatactx%DX,modeldatactx%SKEW, &
                  modeldatactx%DXE)
     if (iecho.eq.8) then
         call ParallelLog("Mmarch", "Updated stress")
     endif

     !*****    Update strains
     call UPSTRA (meshdatactx%X,modeldatactx%D,modeldatactx%DELD, &
                  modeldatactx%STR,meshdatactx%IEN,modeldatactx%LMF, &
                  modeldatactx%TFAULT,modeldatactx%DFAULT,modeldatactx%LMX, &
                  modeldatactx%DX, modeldatactx%SKEW,modeldatactx%DXE)
     if (iecho.eq.8) then
         call ParallelLog("Mmarch", "Updated strains")
     endif

    !***  update strainrates
    call TRATE (meshdatactx%X,modeldatactx%D,modeldatactx%DELD, &
                modeldatactx%RATE,meshdatactx%IEN,modeldatactx%LMF, &
                modeldatactx%TFAULT,modeldatactx%DFAULT,modeldatactx%LMX, &
                modeldatactx%DX, modeldatactx%SKEW,modeldatactx%DXE, &
                NSD,NDOF,NSTR,NEN)
     if (iecho.eq.8) then
         call ParallelLog("Mmarch", "Updated strain rates")
     endif

     !*****  Update visco-plastic strains
     if (IDIG(ICVIS,4,4).ne.0) &
      call UPVPSTR (modeldatactx%BETB,modeldatactx%VPSTR,NSTR)

     !***** Increment total displacements D with increments during this time step DELD
     !***** D = D + DELD
     call UPDATE (modeldatactx%D,modeldatactx%DELD,NDOF*meshdatactx%nvlocal)

     if (iecho.eq.8) then
         call ParallelLog("Mmarch", "updated D")
     endif
 
     !***** update total faulted node displacements
     if (NUMFN.gt.0) then
         !**** add DFAULT to TFAULT
         call UPDATE (modeldatactx%TFAULT,modeldatactx%DFAULT,NDOF*NUMFN)
         call CLEAR (modeldatactx%DFAULT,NDOF*NUMFN)
     endif

     !***** Copy differential displacement during time step DX to equivalent element-based DELX
     if (NUMSLP.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc
        call CLEAR (modeldatactx%DELX,NEN*NDOF*nElemsWithSlipperyNodes)
#else
        call CLEAR (modeldatactx%DELX,NEN*NDOF*NUMEL)
#endif
        call UPDXE (modeldatactx%DX,meshdatactx%IEN, &
                    modeldatactx%LMX,modeldatactx%SKEW,modeldatactx%IDSLE, &
                    modeldatactx%DELX,NDOF,NEN)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
        call UPDATE (modeldatactx%DXE,modeldatactx%DELX, &
                     nElemsWithSlipperyNodes*NEN*NDOF)
#else
        call UPDATE (modeldatactx%DXE,modeldatactx%DELX, &
                     NUMEL*NEN*NDOF)
#endif
        if (iecho.eq.8) then
            call ParallelLog("Mmarch", "updated DXE")
        endif
     endif

     !*****    Increment surface nodal coordinates by changes DELD
     call ADDSRF (modeldatactx%ISURF,modeldatactx%SURF,meshdatactx%X, &
                  modeldatactx%DELD,NSD,NDOF,NSURF,1)


#ifdef EXPERIMENTAL_FRICTION

     if (iterate) then
!        write(*,*) 'NSTEP,IT,BNORM=',NSTEP,ITERNUMBER,SQRT(PSUMSQ)
         DELTP=0d0  ! avoids viscoelastic stress relaxation (STRESN, UPVPSTR) during iterations
         if (IRESDU.ge.1) then
             !***** correct BTOT vector for WINKLF update
             call UPDATE (modeldatactx%BTOT,modeldatactx%B,NEQlocal)
         endif
         !**** clear displacement increments to avoid loading by displacement or velocity bcs during iterations
         modeldatactx%DELD = 0d0
 
         if (NUMSN.gt.0) then
             modeldatactx%DX = 0d0
          endif
         if (NUMFN.gt.0) then
             modeldatactx%DFAULT = 0d0
         endif
     else
         if (nFaultSegmentEntriesGlobal.gt.0) then
           if (iecho.gt.0) then
          write(*,6109) NSTEP,ITERNUMBER,SQRT(PSUMSQ)
6109      format('**** NSTEP ',I3,' required ',I3,' iterations, BNORM=',1PG12.2)
           endif
         endif
         DELTP=tempDELTP

         !**** compute displacements during time step in global coordinates
         if (nFaultSegmentEntriesGlobal.gt.0) then

             modeldatactx%DELD = modeldatactx%cumulativeDisp    

             ! rotate displacement to matck Euler angles
             ! becasue cumulative Disp was not rotated
             call RDISP(modeldatactx%DELD,modeldatactx%SKEW,NDOF,NUMNP)

         endif

!!!!!!!
!        call RDISP(modeldatactx%DELD,modeldatactx%SKEW,NDOF,NUMNP)


         if (NUMSLPglobal.ne.0) then
             !*** compute element-based differential displacements during time step in global coordinates DELX
             if (nFaultSegmentEntriesGlobal.gt.0) then
                 modeldatactx%DX = modeldatactx%cumulativeDX
                 ! DX has been rotated, but it is overwritten by
                 ! cumulativeDX, which was not yet rotated,
                 ! so we need to rotate again.
                 
                 call RDISP(modeldatactx%DX,modeldatactx%SKEW,NDOF,NUMNP)

             endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
             call CLEAR (modeldatactx%DELX,NEN*NDOF*nElemsWithSlipperyNodes)
#else
             call CLEAR (modeldatactx%DELX,NEN*NDOF*NUMEL)
#endif
             ! set DELX, based on values in DX
             call UPDXE (modeldatactx%DX,meshdatactx%IEN, &
                         modeldatactx%LMX,modeldatactx%SKEW,modeldatactx%IDSLE, &
                         modeldatactx%DELX,NDOF,NEN)
         endif
     endif
#endif /* EXPERIMENTAL_FRICTION */

    if (iterate) then
        if (iecho.eq.8) then
            write(*,*) "friction iter number: ", iternumber
            call parallellog("mmarch","starting another friction iteration")
        endif
    endif


enddo ! iterate
!*******************************************************************
!*****  end friction iteration block
!*******************************************************************



!***** write data to FEOUT and FEDSK files
if (DoOut) then 
    call PRINTD (modeldatactx%D,NUMNP,0)
    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "written displacement")
    endif

#ifdef EXPERIMENTAL_FRICTION

    if (NUMSLPglobal.gt.0) then
        call PRINTD (modeldatactx%DIFORC,NUMNP,2)
    else
        ! modeldatactx%DIFORC not allocated
    endif

    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "written differential forces")
    endif
#endif /* EXPERIMENTAL_FRICTION */

    !***** compute and write velocities
    if (DELTP.gt.EPS) then
        factor = 1d0/DELTP
        call BMULT (modeldatactx%DELD,NDOF*NUMNP,factor)
    endif
    call PRINTD (modeldatactx%DELD,NUMNP,1)
    if (DELTP.gt.EPS) then
        ! put it back where we found it.
        factor = DELTP
        call BMULT (modeldatactx%DELD,NDOF*NUMNP,factor)
    endif
    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "written velocities")
    endif

    !***** write faulted node displacements
    call PSPLIT (modeldatactx%TFAULT,modeldatactx%NFAULT,NDOF,NUMFN)
    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "printed faulted node displacements")
    endif

    if (NUMSLPglobal.ne.0) then
        !***** write differential displacements
        call PRINTDXE (modeldatactx%DXE, &
                       meshdatactx%IEN, &
                       modeldatactx%IDSLE,0)

        !***** compute and write differential velocities
        factor = 1d0/DELTP
        if (DELTP.lt.EPS) then
            factor = 1d0
        endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
        call BMULT (modeldatactx%DELX,NDOF*NEN*nElemsWithSlipperyNodes,factor)
#else
        call BMULT (modeldatactx%DELX,NDOF*NEN*NUMEL,factor)
#endif

        call PRINTDXE (modeldatactx%DELX, &
                          meshdatactx%IEN, &
                          modeldatactx%IDSLE,1)

        factor = DELTP
        if (DELTP.lt.EPS) then
            factor = 1d0
        endif
#ifdef EXPERIMENTAL_ReducedSlipAlloc
        call BMULT (modeldatactx%DELX,NDOF*NEN*nElemsWithSlipperyNodes,factor)
#else
        call BMULT (modeldatactx%DELX,NDOF*NEN*NUMEL,factor)
#endif
        if (iecho.eq.8) then
            call ParallelLog("Mmarch", "printed slippery stuff")
        endif
    endif

    !**** print stress
    call PRINTS (modeldatactx%STN,NSTR,NUMEL)
    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "printed stress")
    endif

    !**** print strain
    call PRINTE (modeldatactx%STR,NSTR,NUMEL)
    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "printed strains")
    endif

#ifndef EXPERIMENTAL_NORATESTORAGE
    !*****  Compute and write total strain rate

!**********************************************************************************
!************* ROUTINE IS NOT CORRECT IN THAT STRAIN IS COMPUTED FROM D+DX
!************* IN AN ELEMENT GEOMETRY THAT HAS ALREADY BEEN UPDATED WITH D+DX
!**********************************************************************************

!    call TRATE (meshdatactx%X,modeldatactx%D,modeldatactx%DELD, &
!                modeldatactx%RATE,meshdatactx%IEN,modeldatactx%LMF, &
!                modeldatactx%TFAULT,modeldatactx%DFAULT,modeldatactx%LMX, &
!                modeldatactx%DX, modeldatactx%SKEW,modeldatactx%DXE, &
!                NSD,NDOF,NSTR,NEN)

!    if (iecho.eq.8) then
!        call ParallelLog("Mmarch", "done strain rate")
!    elseif (iecho.eq.6) then
!        call LogTimeStamp("Done trate")
!    endif

    call PRINTR (modeldatactx%RATE,NSTR,NUMEL)

    if (iecho.eq.6) then
        call LogTimeStamp("printed some")
    elseif (iecho.eq.8) then
        call ParallelLog("Mmarch", "written output")
    endif
#endif

    if (NSED.ne.0) then
        call PRTSED (modeldatactx%SURF, modeldatactx%ISURF,NSD,NSURF)
    endif
    if (iecho.eq.6) then
        call LogTimeStamp("printed sediments")
    endif

   if (OUTFIL(1)) then
!        call flushf(luout(1))
    endif

    if (OUTFIL(2)) then
!        call flushf(luout(2))
    endif
    if (iecho.eq.6) then
        call LogTimeStamp("flushed")
    endif

    call DUMPBTOT (modeldatactx%BTOT,NEQ)
    if (iecho.eq.6) then
        call LogTimeStamp("printed dumpbtot")
    endif
endif ! DoOut (writing data block)

if (PrtMat) then
    IPRMAT=IPRMAT+1
endif
if (IPRMAT.GT.NPRMAT) then
    IPRMAT=NPRMAT
endif
if (DoOut) then
    IMPRT=IMPRT+1
endif
if (IMPRT.GT.NMPRT) then
    IMPRT=NMPRT
endif

if (iecho.eq.6) then
    call LogTimeStamp("Done mmarch call")
elseif (iecho.eq.8) then
    call ParallelLog("Mmarch", "Leaving mmarch")
endif

return
end subroutine
!-------------------------------------------------------------------------------
!        Secondary subroutines are listed below in alphabetical order
!-------------------------------------------------------------------------------

#ifdef EXPERIMENTAL_FRICTION

subroutine MAXTRACTION (DX,IDX,IWINX,WINX,DIFFSTRENGTH, &
                        B,DIFORC,NDOF,NUMNP,NEQ,ITERNUMBER, &
                        YFORCE)

use modelctx,         only: getrank ! for debuggign
use timestepmodule
use algebra
use meshdatamodule,   only: meshdatactx
use debugmodule,      only: NicolaiNodeID, &
                            NicolaiFileID, &
                            debug
use iomodule          

implicit none
! loads the B-vector in case the shear traction
! on a fault surpasses its shear strength

!-pass
integer NEQ,NUMNP,NDOF
double precision :: DX(NDOF,meshdatactx%nvlocal) ! local confirmed
double precision :: WINX(NDOF,meshdatactx%nvglobal) ! global confirmed
double precision :: DIFFSTRENGTH(meshdatactx%nvglobal) ! global confirmed
double precision :: B(NEQ)
double precision :: DIFORC(NDOF,meshdatactx%nvlocal) ! local confirmed
double precision :: YFORCE
integer          :: IDX(NDOF,meshdatactx%nvglobal) ! global confirmed
integer          :: IWINX(NDOF,meshdatactx%nvglobal) ! global confirmed
integer          :: ITERNUMBER
!-local
logical firstentry
integer idof,j,k
double precision :: DMAG,DF_NORM,DF_component,pow,frac

integer , allocatable :: loc2gloVertex(:)
integer :: globalID, iVertex
save firstentry
data firstentry/.true./

!make an inverse lookup table
if (.not. allocated(loc2gloVertex)) then
    allocate(loc2gloVertex(meshdatactx%nvglobal))
    do iVertex =1,meshdatactx%nvglobal
        loc2gloVertex(meshdatactx%Glo2LocVertex(ivertex)) = iVertex
    enddo
endif

!***** update the total differential force from the springs
do j=1,NUMNP
    globalID = loc2gloVertex(j+meshdatactx%lowestLocalVertexID)
    if (debug) write(*,*) 'Rank ',getrank(),' MAXTRACTION: local node#=',j, &
               ', global node#=',globalID
#ifdef ZISFAULTNORMAL
    do idof=1,NDOF-1
#else
    do idof=1,NDOF
#endif
        if(IWINX(idof,globalID).ne.0) then
            DIFORC(idof,j)=DIFORC(idof,j)-DX(idof,j)*WINX(idof,globalID)
            if (debug) &
             write(*,*) 'Rank ',getrank(),' MAXTRACTION: DIFORC(idof=',idof,',node=',j,')=', &
                        DIFORC(idof,j) 
        else
            if (debug) &
             write(*,*) 'Rank ',getrank(),' MAXTRACTION: IWINX(idof=',idof,', node=',j,')=', &
                        IWINX(idof,globalID)
        endif
    end do
end do

! compute and impose overshoot
do j=1,NUMNP
  globalID = loc2gloVertex(j+meshdatactx%lowestLocalVertexID)

#ifdef ZISFAULTNORMAL
! RG 2018-11-23:
  do idof=1,NDOF-1
#else
  do idof=1,NDOF
#endif

    if(IWINX(idof,globalID).ne.0) then

      !***** determine the norm of the total differential force
      DF_NORM = SQRT(DOT(DIFORC(1,j),DIFORC(1,j),NDOF))

      !**** calulate whether the slippery nodes are under- or overshooting
      DMAG = DF_NORM - DIFFSTRENGTH(globalID)

      if (debug) write(*,1) DF_NORM,globalID,DIFFSTRENGTH(globalID),DMAG
1     format('MAXTRACTION: norm of DIFORC=',1PE16.2,', DIFFSTRENGTH(node=',I7,')=', &
       1PE16.2,'), overshoot=',1PE16.2)

      if (DMAG.gt.1d0) then
        !**** shear traction surpasses the shear strength, update B-vector

        if (IDX(idof,j).le.0) goto 300

        if (j.gt.size(IDX,2)) then
          write(stderr,*) "Error on thread",getrank(),": MAXTRACTION tries to read the modeldat IDX entry", j
          write(stderr,*) "while there are only: ", size(IDX,2), "nodes in the partition"
          write(stderr,*) "This should not have happened. Contact model support"
          stop "Exiting"
        endif

        if (firstentry) YFORCE = YFORCE + DIFFSTRENGTH(globalID)
        !***** project overshoot on the -DIFORC vector
        DF_component = -DMAG*DIFORC(idof,j)/DF_NORM

        if (debug) write(*,*) 'MAXTRACTION overshoot: DF_component=',DF_component

        !***** subtract overshoot (component) from the load vector
        B(IDX(idof,j)) = B(IDX(idof,j)) + DF_component
        if (debug) write(*,*) 'MAXTRACTION: updated B(ieq=',IDX(idof,j),')=',B(IDX(idof,j))

        if (getrank().eq.0) then ! prevent comfusing output in parallel runs
          if (NicolaiNodeID.gt.0) then
            write(NicolaiFileID, "(a,i12,a,2e25.16,a,2e25.16)") "iter  ", ITERNUMBER, &
               "  has diforc  ", DIFORC(1,NicolaiNodeID), &
               DIFORC(2,NicolaiNodeID), &
               "  and DX  ",     DX(1,NicolaiNodeID), &
               DX(2,NicolaiNodeID)
          endif
        endif

        DIFORC(idof,j) = DIFORC(idof,j) + DF_component
        if (debug) write(*,*) 'MAXTRACTION: updated DIFORC(idof=',idof,',node=',j,')=', &
                   DIFORC(idof,j) 
      end if
300   continue
    end if
  end do
end do

IF (firstentry) then
    if (YFORCE.lt.0d0) then
        write(stderr,*) 'maxtraction: negative differential force error'
        call exitp(1)
    endif
    firstentry = .false.
endif
if (YFORCE.lt.1d-16) write(stderr,*) 'warning: yield strength < 1e-16'

! do j=1,NUMNP
!   do i=1,NDOF
!     if(IWINX(i,j).ne.0) then
!         write(*,20) NSTEP,i,j,DIFORC(i,j),WINX(i,j)
! 20      format('NSTEP=',I3,' DIFORC(',I1,',',I3,')=',1PG12.3,' WINX=',1PG12.3)
!      end if
!    end do
! end do

end subroutine

#endif

!-------------------------------------------------------------------------------
subroutine CKWINK (IWINK,NSTEP,LNGTH,IMWORK)

! Routine for checking if implicit code needs to be invoked because
! Winkler forces are being changed.

! Routine is called to check iwink, and iwtime, so that
! we cannot read the winkler variable directly from modeldatactx.

 integer NSTEP,LNGTH,IMWORK
 dimension IWINK(LNGTH)

!    refactored anyway ?
 if (IMWORK.ne.0) goto 200

do 100 k=1,LNGTH
    LL = IWINK(k)
    if (LL.lt.0) then
        if (NSTEP.ge.-LL) goto 100
    else if (LL.gt.0) then
        if (NSTEP.lt.(LL-1)) goto 100
    else
        goto 100
    endif
    IMWORK = 1
100 continue

200 return
end
!-------------------------------------------------------------------------------
subroutine CKWINX (IWINX,NSTEP,NDOF,LNGTH,IMWORK)

#ifdef EXPERIMENTAL_CYCLIC
use modeldefinition, only: UseCyclicSlip
#endif
use iomodule          

! Routine for checking if implicit code needs to be invoked because
! Differential Winkler forces are being changed.

 integer NSTEP,NDOF,LNGTH,IMWORK
 dimension IWINX(NDOF,LNGTH)
!locl
 integer i,k,itshift
integer :: nstep_minus_itshift

!    refactored anyway ?
if (IMWORK.ne.0) goto 200

itshift = 0

do 100 k=1,LNGTH
#ifdef EXPERIMENTAL_CYCLIC
#ifdef ZISFAULTNORMAL
    do i=1,NDOF-1
         itshift = IWINX(NDOF,k)
#else
    do i=1,NDOF
#endif
         LL = IWINX(i,k)
         if (UseCyclicSlip) then
            nstep_minus_itshift = NSTEP - itshift
            if (nstep_minus_itshift.lt.0) then
                nstep_minus_itshift = nstep_minus_itshift + IABS(LL)
            endif
!             Interpretation of IWINX as periodic
              if (LL.lt.0) then
!                 typically locked, except when MOD(NSTEP,ABS(LL)) < 10
                  if (MOD(nstep_minus_itshift,IABS(LL)).ne.0) goto 100
!                 if (MOD(NSTEP,IABS(LL)).lt.10) goto 100
              else if (LL.gt.0) then
!                 typically unlocked, except when MOD(NSTEP,LL) < 10
                  if (MOD(nstep_minus_itshift,LL).ne.10) goto 100
!                 if (MOD(NSTEP,LL).ge.10) goto 100
              else
                  goto 100
              endif
        else
              if (LL.lt.0) then ! LL is IWINX entry
                  if (NSTEP.ge.-LL) goto 100
              else if (LL.gt.0) then
                  if (NSTEP.lt.(LL-1)) goto 100
              else
                  goto 100
              endif
        endif
    enddo
#else
    do i=1,NDOF
         LL = IWINX(i,k)
         if (LL.lt.0) then ! LL is IWINX entry
             if (NSTEP.ge.-LL) goto 100
         else if (LL.gt.0) then
             if (NSTEP.lt.(LL-1)) goto 100
         else
             goto 100
         endif
    enddo
#endif

IMWORK = 1

100    continue

200    return
 end
!-------------------------------------------------------------------------------
 subroutine DUMPBTOT (BTOT,NEQ)

 use filemodule, only: outdir
use debugmodule,     only: xit
use iomodule          

!
 implicit none
!-pass
 integer NEQ
double precision :: BTOT(NEQ)

!-local
integer :: lutmp
integer, external :: nextlu
integer :: i, ios
!
lutmp = nextlu(0)

call fbopen (lutmp,trim(outdir)//'/'//'BTOT.DAT','unknown')

write(lutmp,err=100,iostat=ios) (BTOT(i),i=1,NEQ)

call closef(lutmp)
!
 return
100    write(stderr,1) ios
 1   format(1x,'DUMPBTOT: write error ',I5)
 call xit(1," ")
 end
!-------------------------------------------------------------------------------
subroutine FORMDM (DMAT,DELTP,ALFAP,STN,IDMAT,PRPLAS,PRPMAT, BETB,NSTR)

!
! Augments the inverse of the material matrix:
! (D)-1 + ALFAP*DELTP*(BETA)-1
! and inverts it (Zienckiewics, p.476).
!
implicit none
!-pass

integer NSTR
double precision :: ALFAP, DeltP

double precision :: STN, DMAT, PRPLAS, PRPMAT, BETB
integer          :: IDMAT
dimension STN(NSTR),IDMAT(NSTR),DMAT(NSTR,NSTR),PRPLAS(9), PRPMAT(6),BETB(NSTR)

!-locl
double precision :: BETAJ(6,6)
integer          :: i, j

! create Jacobian of visco-plastic strain rate (stress) tensor
! Makeing BETAJ
call FORMJB (BETAJ,STN,PRPMAT,PRPLAS,BETB)

do i=1,NSTR
    do j=1,NSTR
        DMAT(i,j) = DMAT(i,j) + ALFAP*DELTP * BETAJ(i,j)
    enddo
enddo

!       invert updated material matrix
call INVMAT (DMAT,IDMAT,NSTR)

!
return
end
!-------------------------------------------------------------------------------
subroutine FORMEF (B,X,D,DELD,SKEW,IEN,LM,LMX,MAT, &
   BETA,BETB,PRPMAT,PRPLAS,DMAT,STN,LMF,TFAULT,DXE, &
   NSD,NDOF,NSTR,NEN)

! Updates the load vector B with stress changes from viscoelastic and/or
! viscoplastic flow, and with forces equivalent to velocity boundary conditions
 
use constants,       only: zero, &
                           useDXE
use meshdatamodule,  only: meshdatactx
USE TIMESTEPMODULE
USE MATERIALSMODULE
USE MODELDEFINITION
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: modeldatactx, &
                           elementSlipperyID
#else
use modeldatamodule, only: modeldatactx
#endif
USE ALGEBRA
use modelctx,        only: getrank ! debug only
use iomodule          
use debugmodule
 
implicit none

!-pass
integer NDOF,NSD,NEN,NSTR

integer          :: mat
double precision :: B, X, D, DELD, SKEW, PRPMAT, PRPLAS
double precision :: BETA, BETB, DMAT, STN, TFAULT

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer :: LMX(NDOF,NEN,nElemsWithSlipperyNodes)
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
integer :: LMX(NDOF,NEN,meshdatactx%neglobal)
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif

integer :: IEN, LM, LMF

dimension B(*),X(NSD,*),D(NDOF,*),DELD(NDOF,*),SKEW(2,*), &
   IEN(NEN,*),LM(NDOF,NEN,*),MAT(*),PRPMAT(6,*), &
   PRPLAS(9,*),BETA(NSTR,NUMEL),BETB(NSTR,*),DMAT(NSTR,NSTR,*), &
   STN(NSTR,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*)
!-locl
logical QUAD,NONZERO

double precision :: S(576),STEMP(576),PL(24),DL(24),XL(24),EVP(6),EE(6)

integer :: n, nn, m, ii
integer :: ierr

double precision :: thick

do nn=1,NUMEL

    n = nn
    m = MaterialIDOfElement(n)

    ! clear local load vector
    call CLEAR (PL,NEE,"PL")

    ! localize coordinates
    call LCOORD (X,XL,IEN(1,n),nn)
    ! large deformation update
    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,1)        
        call ADDSNE (DL,n,NDOF,NEN, useDXE)
        call REZONE (XL,DL,1d0)
    endif

    ! calculate local stress changes EVP from strains
    call FORMVE (STN(1,n),BETA(1,n),BETB(1,n),DMAT(1,1,n),EVP, &
                 EE,PRPMAT(1,m),PRPLAS(1,m),DELTP,NSTR)

    ! if (debug) write(*,*) 'Rank=',getrank(),' FORMEF: norm of EVP(elem#=',n,')=', &
    !  SQRT(DOT_PRODUCT(EVP,EVP))

    ! calculate effective force at each node
    QUAD = (IEN(3,n).ne.IEN(4,n))
    THICK = PRPMAT(6,m)
    call EFORCE (EVP,PL,XL,THICK,QUAD,ierr)

    if (ierr.ne.0) then
        write(stderr,1) n
 1      format(1x,'formef: error in element ',I12)
        call stoper()
    endif

    if (NUMROT.gt.0) then
        call RPFORC (PL,SKEW,IEN(1,n),NDOF,NEN)
    endif

    ! Localize velocity boundary conditions
    call LDISBC (DL,DELD,meshdatactx%IEN(1,n),n,NDOF,NEN,meshdatactx%gloElt(nn))

    NONZERO = .false.
    do ii=1,NDOF*NEN
        NONZERO = (NONZERO .or. (DL(ii).ne.ZERO))
    enddo

    if (NONZERO) then
        ! clear local stiffness matrix
        call CLEAR (S,NEE*NEE)
     
        ! calculate upper triangle of local stiffness matrix
        call STIFF (DMAT(1,1,n),XL,S,THICK,QUAD,ierr)
        if (ierr.ne.0) then
            write(stderr,1) n
            call stoper()
        endif
     
        ! make stiffness matrix symmetric
        call SYMMET (S,NEE)
     
        ! if rotated bcs rotate stiffness matrix
        if (NUMROT.gt.0) then
            CALL RSTIFF (S,STEMP,IEN,SKEW,NDOF,NEN,NEE,n,n)
        endif

        ! adjust load vector PL for forces equivalent to velocity bcs
        ! PL(i) = PL(i) - S(i,j)*DL(j)
        call DISBC (DL,PL,S,NEE)
    endif

    ! update global load vector.. B = B + PL * LM(1,1,n)
    call ADDFOR (B,PL,LM(1,1,n),NEE,nn,10)
    ! if (debug) write(*,*) 'Rank=',getrank(),' FORMEF: PL(elem#=',n,')=',PL(1:NEE)

    if (NUMSN.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc
        if(elementSlipperyID(n,110).gt.0) then
            call ADDFOR (B,PL,LMX(1,1,elementSlipperyID(n,110)),NEE,nn,11)
        else
            ! this element has no slippery nodes. Skip it.
        endif
#else
        call ADDFOR (B,PL,LMX(1,1,n),NEE,nn,11)
        ! if (debug) write(*,*) 'Rank=',getrank(),' FORMEF: PL(elem#=',n,')=',PL(1:NEE)
#endif
    endif

enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine FORMVE (STN,BETA,BETB,DMAT,EVP,EE,PRPMAT,PRPLAS, DELTP,NSTR)
!
! Subroutine to compute stress change EVP due to viscous and
! viscoplastic strain (EE). EE is computed from strain rate BETA.
!

use modeldatamodule,  only: modeldatactx
use modeldefinition
use algebra,          only: vmprd
use debugmodule

implicit none
!-pass
integer          :: NSTR
double precision :: STN(NSTR),BETA(NSTR),EVP(NSTR),EE(NSTR),BETB(NSTR), &
                    DMAT(NSTR,NSTR),PRPMAT(6),PRPLAS(9)
double precision :: deltp
!-local
integer :: i

! compute strain rates from stress field -> beta, beta*stress
call FORMBT (STN,BETA,BETB,PRPMAT,PRPLAS)

do i=1,NSTR
    EE(i) = DELTP * BETA(i)
enddo

! compute EVP(I) = DMAT(I,J)*EE(J) (DMAT = Material matrix, EE = local strains)
call VMPRD (DMAT,EE,EVP,NSTR,NSTR)

return
end
!-------------------------------------------------------------------------------
subroutine IMPLCT (A,IDIAG,B,ID,IDX,X,D,SKEW,IWINK, &
 WINK,IWINX,WINX,IEN,LM,LMX,MAT,PRPMAT,DMAT,STN, &
 IDMAT,LMF,TFAULT,PRPLAS,IWELM,IWSIDE,IWTIME,WPRES, &
 DXE,NSLIP,S,IA,IAX,BETB,NSD,NDOF,NSTR,NEN,NSTEP)

use iomodule          
   use meshdatamodule,  only: meshdatactx
#ifdef EXPERIMENTAL_ReducedSlipAlloc
   use modeldatamodule, only: elementSlipperyID, modeldatactx
#else
   use modeldatamodule, only: modeldatactx
#endif
   USE TIMESTEPMODULE
   USE MATERIALSMODULE
   USE MODELDEFINITION
   USE MODELCTX,       only: getrank
   USE ALGEBRA,        only: clear, symmet
   use debugmodule, only: debug, iecho, parallelLog
use constants,      only: useDXE
#ifdef SPARSE
use spetscmodule, only:  M_A, userFinalize
#else
!use petscmodule, only:  M_A
#endif
!
! Initiates implicit solution mode
!
implicit none
!-pass
integer NDOF,NSD,NEN,NSTR

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer :: LMX(NDOF,NEN,nElemsWithSlipperyNodes)
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
integer :: LMX(NDOF,NEN,meshdatactx%neglobal)
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif



integer :: id, idiag, idx, iwink, iwinx, ien, lm, idmat(nstr), lmf
integer :: iwelm, iwside, iwtime, nslip
integer :: ia(ndof, numnp), iax(ndof, numnp), nstep, mat

double precision :: a, b, x, d, skew, wink, winx, prpmat
double precision :: dmat, stn, tfault, prplas, wpres, s, betb

dimension A(*),IDIAG(NEQ),B(*),ID(NDOF,*),IDX(NDOF,*),X(NSD,*), &
 D(NDOF,*),SKEW(2,*),IWINK(NDOF,*),WINK(NDOF,*), &
 IWINX(NDOF,*),WINX(NDOF,*),IEN(NEN,*),LM(NDOF,NEN,*), &
 MAT(*),PRPMAT(6,*),DMAT(NSTR,NSTR,*), &
 STN(NSTR,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
 PRPLAS(9,*),IWELM(*),IWSIDE(*),IWTIME(*),WPRES(*), &
 NSLIP(5,*),S(MEQ,*), BETB(NSTR,*)

#ifdef SPARSE

!#include "petsc/finclude/petscsysdef.h"
!#include "petsc/finclude/petscvecdef.h"
!#include "petsc/finclude/petscmatdef.h"
!#include "petsc/finclude/petsckspdef.h"

#endif
!-locl
double precision :: SL(576),STEMP(576),P(24),XL(24),DL(24)
logical QUAD
integer :: i,ii,j,jj,k,m,n,nn
double precision :: E, pois, thick
integer :: ierr

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer, parameter :: dummyLMX(8) = (/0,0,0,0,0,0,0,0/)
integer :: slipSequenceNr
#endif

if (getrank().eq.0) then 
    if (iecho.eq.1) then
!        write(*,*) 'Reforming mechanical stiffness matrix'
    endif
endif

!       clear stiffness and material matrices
#ifdef SPARSE
call MatZeroEntries(M_A,ierr)

if (ierr.ne.0) then
    write(stderr,2) ierr
 2       format(1x,'IMPLCT: MatZeroEntries error ',I4)
    call UserFinalize()
endif
#else
call CLEAR(A,NSIZEA,"A")
#endif
!
    if (iecho.eq.8) then
        call ParallelLog("Mmarch", "Implicit loops over elements")
    endif


do nn=1,NUMEL

    n = nn
    m = MaterialIDOfElement(n)
    E    = matdatactx%PRPMAT(1,m)
    POIS = matdatactx%PRPMAT(2,m)
    THICK= matdatactx%PRPMAT(6,m)

    !**** clear element load and stiffness matrix
    call CLEAR (P,NEE,"P")
    call CLEAR (SL,NEE*NEE)
 
    !**** set up local material matrix
    call MATERL (DMAT(1,1,n),E,POIS)

    !***** invert local material matrix and put result back into DMAT
    call INVMAT (DMAT(1,1,n),IDMAT,NSTR)

    !***** create [D^-1 + DELTP*ALFAP*S]^-1 (strain rate Jacobian S)
    call FORMDM (DMAT(1,1,n),DELTP,ALFAP,STN(1,n),IDMAT, &
     PRPLAS(1,m),matdatactx%PRPMAT(1,m),BETB(1,n),NSTR)

    !***** localize coordinate data
    call LCOORD (X,XL,IEN(1,n),nn)

    !***** large deformation update
    if (LGDEF.ne.0) then
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,2)
        call ADDSNE (DL,n,NDOF,NEN, useDXE)
        call REZONE (XL,DL,1d0)
    endif

    !***** calculate upper triangle of local stiffness matrix
    QUAD = (IEN(3,n).ne.IEN(4,n))
    call STIFF (DMAT(1,1,n),XL,SL,THICK,QUAD,ierr)

    if (ierr.ne.0) then
        write(stderr,3) n
3       format(1x,'implct: error in element ',I5)
        call stoper()
    endif

    !***** symmetrize local stiffness matrix
    call SYMMET(SL,NEE)

    !***** rotate local stiffness matrix to global coordinate system
    if (NUMROT.gt.0) then
        call RSTIFF(SL,STEMP,IEN,SKEW,NDOF,NEN,NEE, n,n)
    endif

    !***** add local stiffness to global arrays
#ifdef EXPERIMENTAL_ReducedSlipAlloc

    slipSequenceNr = elementSlipperyID(n,111)

    if (slipSequenceNr.gt.0) then
        call ADDSTF (A, &
                 SL, &   
                 IDIAG, &
                 modeldatactx%LM(1,1,n), &
                 modeldatactx%LMX(1,1,slipSequenceNr), &
                 NEE, &
                 NUMSN, &
                 meshdatactx%gloElt(n))
    else
        call ADDSTF (A, &
                 SL, &
                 IDIAG, &
                 modeldatactx%LM(1,1,n), &
                 dummyLMX, &
                 NEE, &
                 NUMSN, &
                 meshdatactx%gloElt(n))
    endif


!    call ADDSTF (A,SL,IDIAG,LM(1,1,n),elementSlipperyID(n,111),NEE,NUMSN,n)
    !note that elementSlipperyID(n) can be 0
#else
    call ADDSTF (A,SL,IDIAG,LM(1,1,n),LMX(1,1,n),NEE,NUMSN,n)
#endif

    if (showfbc) then
        call addneqbc (SL,IEN(1,n),S,IA,IAX,NEN,NEE, NDOF,MEQ)
    endif
enddo

if (iecho.eq.8) then
    call ParallelLog("Mmarch", "Implicit does Winkler stuff")
endif

! Nicolai needs this
! add regular Winkler forces to global stiffness
if (NWINK .gt.0) then
    call WINKLR (ID ,IWINK,WINK,NDOF,NUMNP,0)
endif

! and the differential Winkler forces
if (NWINKX.gt.0) then
    call WINKLR (IDX,IWINX,WINX,NDOF,NUMNP,1)
endif

! Winkler pressures.
if (NUMWNK.ne.0) then
    call ADDWNK (IWELM,IWSIDE,IWTIME,WPRES,IEN,MAT, &
 X,D,PRPMAT,A,LM,IDIAG,LMF,TFAULT,SKEW, &
 DXE,LMX,NSLIP)
endif

#ifndef SPARSE
!       factorize global stiffness
call SOLVE (A,B,B,IDIAG,ID,IDX,NUMNP,NEQ,NUMSN, &
 NDOF,.TRUE.,.FALSE.)
#endif
 
return
end

!----------------------------------------------------------------------------

   subroutine printlocaldata(ien, xl, dl,n)

   USE MODELCTX
   USE MESHDATAMODULE
   USE MODELTOPOLOGY

   integer IEN(*),i,j,n
   double precision XL(NSD,*),DL(NDOF,*)

!      write(FILE_outputf,*) 'some data for local element ',n, &
!        ' which is global element ', meshdatactx%gloElt(n)

!      write(FILE_outputf,*) (IEN(i),i=1,NEN)

!      do i=1,NEN
!          write(FILE_outputf,*) (XL(j,i),j=1,NSD)
!      enddo

!      do i=1,NEN
!          write(FILE_outputf,*) (DL(j,i),j=1,NDOF)
!      enddo

   return
   end
!-------------------------------------------------------------------------------
 subroutine PRINTR (RATE,NSTR,NUMEL)

 USE LABELMODULE
 use debugmodule, only: iecho, xit
use modeltopology, only: ndof
use formatsmodule, only: FMT_fedsk_stressStrainRate_nstr3, &
                         FMT_fedsk_stressStrainRate_nstr4, &
                         FMT_fedsk_stressStrainRate_nstr5, &
                         FMT_fedsk_stressStrainRate_nstr6
use iomodule

! Routine to print total strain rates

 implicit none
!-pass
 integer NSTR,NUMEL
double precision :: RATE
 dimension RATE(NSTR,*)
!-locl
! character(len=19) fmt
integer :: m, nn, i, n

!-init
 nn = 51

if (iecho.eq.3) then
    write(stdout,1)
 1   format(1x,'if (NUMEL <> 0) {'/ &
 5x,'for (n=1;n<=NUMEL;n++) {'/ &
 9x,'write n,(RATE(i,n),i=1,NSTR) (i5,x,6e16.7)'/ &
 5x,'}'/1x,'}')
endif

if (NUMEL.eq.0) then
    return
endif

if (OUTFIL(2)) then
    if (ascout) then

        do n=1,NUMEL
            if      (NSTR.eq.6) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr6) n,(rate(i,n),i=1,NSTR)
            else if (NSTR.eq.5) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr5) n,(rate(i,n),i=1,NSTR)
            else if (NSTR.eq.4) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr4) n,(rate(i,n),i=1,NSTR)
            else if (NSTR.eq.3) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr3) n,(rate(i,n),i=1,NSTR)
            endif

        enddo
    else
        write(luout(2),err=100) ((RATE(i,n),i=1,NSTR),n=1,NUMEL)
    endif
endif

if (.not. OUTFIL(1)) then
    return
endif

do n=1,NUMEL
    m = n

    if (nn.gt.50) then
        nn = 0
        write(luout(1),3) (LABELR(i),i=1,NSTR)
 3           format(///1X,'T O T A L   S T R A I N   R A T E'/// &
         3X,' ELEM #',5x,6(a4,10x))
        write(luout(1),4)
 4           format(1X)
    endif

    write(luout(1),5) m,(RATE(i,m),i=1,NSTR)
 5        format(1x,i8,1x,6(1pe14.5))
    nn=nn+1
enddo

return

100    call prerr('printr')
 call xit(1," ")
 end
!-------------------------------------------------------------------------------
 subroutine PRINTE (STR,NSTR,NUMEL)

USE LABELMODULE
use debugmodule,   only: iecho, xit
use modeltopology, only: ndof
use formatsmodule, only: FMT_fedsk_stressStrainRate_nstr3, &
                         FMT_fedsk_stressStrainRate_nstr4, &
                         FMT_fedsk_stressStrainRate_nstr5, &
                         FMT_fedsk_stressStrainRate_nstr6
use debugmodule,   only: isFileOpen
use iomodule

! Routine to print finite strains

implicit none
!-pass
 integer NSTR,NUMEL
double precision :: STR
 dimension STR(NSTR,*)
!-locl
double precision :: factor
! character(len=19) fmt
integer :: nn,n,i,m
!-init
 nn = 51

if (iecho.eq.3) then
    write(stdout,1)
 1   format(1x,'if (NUMEL <> 0) {'/ &
 5x,'for (n=1;n<=NUMEL;n++) {'/ &
 9x,'write n,(STR(i,n),i=1,NSTR) (i5,x,6e16.7)'/ &
 5x,'}'/1x,'}')
endif

if (NUMEL.eq.0) then
    return
endif



if (OUTFIL(2)) then
    if (ascout) then
        do n=1,NUMEL
            if (NSTR.eq.6) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr6) n,(STR(i,n),i=1,NSTR)
            else if (NSTR.eq.5) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr5) n,(STR(i,n),i=1,NSTR)
            else if (NSTR.eq.4) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr4) n,(STR(i,n),i=1,NSTR)
            else if (NSTR.eq.3) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr3) n,(STR(i,n),i=1,NSTR)
            endif
        enddo
    else
        write(luout(2),err=100) ((STR(i,n),i=1,NSTR),n=1,NUMEL)
    endif
endif


if (.not. OUTFIL(1)) then
    return
endif

do n=1,NUMEL
    m = n

    if (nn.gt.50) then
        nn = 0
        write(luout(1),3) (LABELR(i),i=1,NSTR)
 3           format(///1X,'F I N I T E   S T R A I N'/// &
         3X,' ELEM #',5x,6(a4,10x))
        write(luout(1),4)
 4           format(1X)
    endif

    write(luout(1),5) m,(STR(i,m),i=1,NSTR)
 5        format(1x,i8,1x,6(1pe14.5))
    nn=nn+1
enddo

return

100    call prerr('printe')
 call xit(1," ")
 end
!-------------------------------------------------------------------------------
subroutine PRTSED (SURF,ISURF,NSD,NSURF)

USE LABELMODULE
use debugmodule, only: iecho, xit
use timestepmodule, only: NINTG
use iomodule

! Routine to print actual surface

 implicit none
!-pass
 integer NSD,NSURF
double precision :: SURF
integer :: isurf
 dimension SURF(NSD,*),ISURF(*)
!-local
integer :: n,nn,m,i

!-init
 nn = 51

if (iecho.eq.3) write(stdout,1)
 1   format(1x,'if (NSURF <> 0) {'/ &
 5x,'for (n=1;n<=NSURF;n++) {'/ &
 9x,'write n,(SURF(i,n),i=1,NSD) (i5,x,3e16.7)'/ &
 5x,'}'/1x,'}')

 if (NSURF.le.0.or.NINTG.le.0) return

 if (OUTFIL(2)) then
!         write (*,*) 'Going to write actual surface'
if (ascout) then
do n=1,NSURF
            write(luout(2),2) n,(SURF(i,n),i=1,NSD)
 2               format(I8,1x,3(1pg15.6))
enddo
else
write(luout(2),err=100) &
       ((SURF(i,n),i=1,NSD),n=1,NSURF)
endif
 endif

 if (.not. OUTFIL(1)) return

 do n=1,NSURF
m = ISURF(n)
if (nn.gt.50) then
    nn = 0
    write(luout(1),3) (LABELD(i),i=1,NSD)
 3           format(///1X,'S U R F A C E',3X,'U P L I F T'/// &
         3X,' NODE #',5x,3(a4,10x))
write(luout(1),4)
 4           format(1X)
endif
write(luout(1),5) m,(SURF(i,n),i=1,NSD)
 5       format(1x,i8,1x,3(1pe14.5))
nn=nn+1
 enddo
 return

100    call prerr('PRTSED')
 call xit(1," ")
 end

!-------------------------------------------------------------------------------

subroutine STRESN (X,D,DELD,DMAT,BETA,STN,IEN,LMF,TFAULT,DFAULT, &
   LMX,DX,SKEW,DXE)


#ifdef EXPERIMENTAL_ReducedSlipAlloc
   use modeldatamodule, only: modeldatactx, elementSlipperyID
   USE MODELDEFINITION, only: LGDEF, numel, nElemsWithSlipperyNodes
#else
use modeldatamodule, only: modeldatactx
USE MODELDEFINITION, only: LGDEF, numel
use meshdatamodule,  only: meshdatactx
#endif
   USE TIMESTEPMODULE
   USE ALGEBRA
use modelctx, only: getrank ! for debugging only 
use constants, only: useDXE
use iomodule          
use modeltopology, only: ndof, nen, nsd, nstr
!
! Program to compute changes in stress and update STN
!
 implicit none

!-pass

integer :: ien,lmf
double precision :: x,d,deld,dmat,beta,stn,tfault,dfault,dx,skew

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer          :: LMX(NDOF,NEN,nElemsWithSlipperyNodes)
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
integer          :: LMX(NDOF,NEN,meshdatactx%neglobal)
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif

dimension X(NSD,*),D(NDOF,*),DELD(NDOF,*),DMAT(NSTR,NSTR,*), &
   BETA(NSTR,*),STN(NSTR,*),IEN(NEN,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),DFAULT(NDOF,*),DX(NDOF,*), &
   SKEW(2,*)
!-locl
double precision :: ee(6),delstr(6),db(6),xl(24),dl(24)
integer :: n,nn,i,ierr

do nn=1,NUMEL
    n = nn
    call LCOORD (X,xl,IEN(1,n),n)

    if (LGDEF.ne.0) then
        !**** Update element geometry with displacements at end of previous time step
        call LDISP (dl,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,3)
!#ifdef EXPERIMENTAL_ReducedSlipAlloc
        call ADDSNE (dl,n,NDOF,NEN,useDXE)
!#else
!       call ADDSNE (dl,DXE(1,1,n),NDOF,NEN)
!#endif
        call REZONE (xl,dl,1d0)  ! xl = xl + dl * 1d0

    endif


    !***** Compute infinitesimal strain increment (calling with DELD (delta D), not D itself)
#ifdef EXPERIMENTAL_ReducedSlipAlloc
    if(elementSlipperyID(n,112).gt.0) then
        call STRAIN (xl,IEN(1,n),DELD,DX,DFAULT, &
                     LMX(1,1,elementSlipperyID(n,112)),n,SKEW, &
                     ee,NDOF,NSD,NEN,NSTR,nn,ierr)
    endif
#else
    call STRAIN (xl,IEN(1,n),DELD,DX,DFAULT, &
       LMX(1,1,n),n,SKEW, &
       ee,NDOF,NSD,NEN,NSTR,nn,ierr)

#endif

    if (ierr.ne.0) then
        write(stderr,1) n
1       format(1x,'stresn: error in element ',I12)
        call stoper()
    endif
 
    !*****        Calculate elastic strain change
    do i=1,NSTR
        delstr(i) = ee(i) - DELTP * modeldatactx%BETA(i,n)
    enddo
 
    !***** Compute product of DMAT with column vector DELSTR into DB
    call VMPRD (DMAT(1,1,n),delstr,db,NSTR,NSTR)

    !*****        Update stresses
    call UPDATE (STN(1,n),db,NSTR)

enddo

return
end
!-------------------------------------------------------------------------------
subroutine UPSTRA (X,D,DELD,STR,IEN,LMF,TFAULT,DFAULT, &
                   LMX,DX,SKEW,DXE)

USE MODELDEFINITION
USE ALGEBRA,         only: UPDATE
use modeltopology,   only: NSD, NDOF,NSTR,NEN
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: elementSlipperyID
#else
use meshdatamodule,  only: meshdatactx
#endif
use constants,       only: useDXE
use iomodule          
!
! Program to update strains
!
implicit none

!-pass

integer :: ien, lmf
double precision :: x,d,dx,deld,str,tfault,dfault,skew

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer          :: LMX(NDOF,NEN,nElemsWithSlipperyNodes)
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
integer          :: LMX(NDOF,NEN,meshdatactx%neglobal)
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif



dimension X(NSD,*),D(NDOF,*),DELD(NDOF,*),STR(NSTR,*), &
 IEN(NEN,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),DFAULT(NDOF,*), &
 DX(NDOF,*), SKEW(2,*)

!-locl
integer nn,n, ierr
double precision :: ee(6),xl(24),dl(24)
!
do nn=1,NUMEL
    n = nn
!           Compute element geometry at end of previous time step
    call LCOORD (X,xl,IEN(1,n),nn)

    if (LGDEF.ne.0) then
        call LDISP (dl,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,4)
        call ADDSNE (DL,n,NDOF,NEN, useDXE)
        call REZONE (xl,dl,1d0)
    endif

!           Compute infinitesimal strain increment

#ifdef EXPERIMENTAL_ReducedSlipAlloc
    if(elementSlipperyID(n,113).gt.0) then
        call STRAIN (xl,IEN(1,n),DELD,DX,DFAULT, &
                     LMX(1,1,elementSlipperyID(n,113)),n,SKEW, &
                     ee,NDOF,NSD,NEN,NSTR,nn,ierr)
    endif
#else
    call STRAIN (xl,IEN(1,n),DELD,DX,DFAULT, &
       LMX(1,1,n),n,SKEW, &
       ee,NDOF,NSD,NEN,NSTR,nn,ierr)

#endif



    if (ierr.ne.0) then
        write(stderr,1) n
 1           format(1x,'upstra: error in element ',I6)
        call stoper()
    endif

    call UPDATE (STR(1,n),ee,NSTR)
enddo
!
return
end
!-------------------------------------------------------------------------------


subroutine UPVPSTR (BETB,VPSTR,NSTR)

USE TIMESTEPMODULE
USE MODELDEFINITION
!
! updates visco-plastic strains
!
implicit none

!-pass
integer NSTR
double precision :: betb, vpstr
dimension BETB(NSTR,*),VPSTR(NSTR,*)
!-locl
integer n,i
!
do n=1,NUMEL
    do i=1,NSTR
        VPSTR(i,n) = VPSTR(i,n) + DELTP*BETB(i,n)
    enddo
enddo
!
return
end
!
!-------------------------------------------------------------------------------
subroutine TBDYLD (MAT,IEN,X,D,LMF,TFAULT,LMX, &
                   TLST,T,TX,TLSTX,PRPMAT,PRPTEM,GRAV,B,LM,SKEW,LMTX,DXE)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
USE MATERIALSMODULE, only: LMAT, &
                           matdatactx, &
                           MaterialIDOfElement
use modeldatamodule, only: elementSlipperyID
USE MODELDEFINITION, only: numel,numsn,numrot,lgdef,nee, &
                           nElemsWithSlipperyNodes, &
                           igrav
#else
USE MATERIALSMODULE, only: LMAT, &
                           matdatactx, &
                           MaterialIDOfElement
!use modeldatamodule, only: elementSlipperyID
USE MODELDEFINITION, only: numel,numsn,numrot,lgdef,nee, igrav
use meshdatamodule,  only: meshdatactx
#endif


USE ALGEBRA
use constants, only: one, useDXE
use modeltopology, only: NEN,NSD,NDOF,NSTR
use iomodule          
!
! add body forces due to thermal expansion
!
 implicit none
!-pass

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer          :: LMX(NDOF,NEN,nElemsWithSlipperyNodes)
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
integer          :: LMX(NDOF,NEN,meshdatactx%neglobal)
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif


integer :: ien,lmf,lm,lmtx, mat
double precision :: x,d,tfault,tlst,t,prpmat,grav,b,prptem,skew,tx,tlstx

 dimension MAT(*),IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),TLST(*),T(*),PRPMAT(6,*),GRAV(NDOF), &
   B(*),LM(NDOF,NEN,*),PRPTEM(6,*),SKEW(2,*), &
   TX(*),LMTX(NEN,*),TLSTX(*)

!-locl
 logical quad
double precision :: pl(24),xl(24),dl(24),dmatl(36),st(6),ee(6)
double precision :: e,pois,densit,texp,thick
integer          :: ierr
double precision :: t0,t1,drho
integer          :: m,n
double precision, external :: Tcentr

!
 do n=1,NUMEL
m = MaterialIDOfElement(n)

e      = matdatactx%PRPMAT(1,m)
pois   = matdatactx%PRPMAT(2,m)
densit = matdatactx%PRPMAT(5,m)
texp   = matdatactx%PRPTEM(3,m)
thick  = matdatactx%PRPMAT(6,m)
ierr   = 0
!
quad = (IEN(3,n).ne.IEN(4,n))
call CLEAR (pl,NEE,"pl")
!
call LCOORD (X,xl,IEN(1,n),n)

if (LGDEF.ne.0) then
    call LDISP  (dl,D,IEN(1,n),NDOF,NEN)
    call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,5)
    call ADDSNE (DL,n,NDOF,NEN, useDXE)
    call REZONE (xl,dl,1d0)
endif
!
t0 = Tcentr (TLST,TLSTX,xl,IEN(1,n),LMTX(1,n),ierr)

if (ierr.ne.0) then
    write(stderr,1) n
 1           format(1x,'tbdyld: error in element ',I5)
    call stoper()
endif

t1 = Tcentr (T,TX,xl,IEN(1,n),LMTX(1,n),ierr)

if (ierr.ne.0) then
    write(stderr,1) n
    call stoper()
endif
!
!        add density change contribution to local load
drho = densit * texp * (t0-t1)

! generate pl from GRAV and the rest
! pl is a *change* in gravitaty induced force
if (igrav.gt.0) then
    call GRAVLD (pl,xl,GRAV,thick,drho,quad,ierr)
endif

if (ierr.ne.0) then
    write(stderr,1) n
    call stoper()
endif

!        calculate elastic material matrix
!            call MATERL (dmatl,e,pois)
!        calculate strain increment due to thermal expansion
!        call TSTRN (ee,t1-t0,texp)
!        calculate elastic stress increment
!        call VMPRD (dmatl,ee,st,NSTR,NSTR)
!        add thermal pressure contribution to local load
!        call EFORCE (st,pl,xl,thick,quad,ierr)
!        if (ierr.ne.0) then
!        write(stderr,1) n
!        call stoper()
!        endif
!
if (NUMROT.gt.0) then
    call RPFORC (pl,SKEW,IEN(1,n),NDOF,NEN)
endif

call ADDFOR (B,pl,LM(1,1,n),NEE,n,12)

if (NUMSN.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc
    if(elementSlipperyID(n,114).gt.0) then
        call ADDFOR (B,PL,LMX(1,1,elementSlipperyID(n,114)),NEE,n,13)
    else
        ! this element has no slippery nodes. Skip it.
    endif
#else
    call ADDFOR (B,PL,LMX(1,1,n),NEE,n,13)
#endif
endif


 enddo
!
 return
 end

#ifndef EXPERIMENTAL_NORATESTORAGE
!-------------------------------------------------------------------------------
subroutine TRATE (X,D,DELD,RATE,IEN,LMF,TFAULT,DFAULT,LMX, &
   DX,SKEW,DXE,NSD,NDOF,NSTR,NEN)


use modeldatamodule, only: modeldatactx
USE TIMESTEPMODULE

#ifdef EXPERIMENTAL_ReducedSlipAlloc
USE MODELDEFINITION, only: numel,lgdef, nElemsWithSlipperyNodes
use modeldatamodule, only: elementSlipperyID
#else
USE MODELDEFINITION, only: numel,lgdef
use meshdatamodule,  only: meshdatactx
#endif
use iomodule          
use constants, only: eps, useDXE
!use meshdatamodule, only: meshdatactx

!
! Computes total strain rate
! of each element
!
 implicit none

!-pass
integer NSD,NDOF,NSTR,NEN

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer          :: LMX(NDOF,NEN,nElemsWithSlipperyNodes)
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
integer          :: LMX(NDOF,NEN,meshdatactx%neglobal)
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif

integer :: ien,lmf
double precision :: x,d,deld,rate,skew,tfault,dfault,dx

dimension X(NSD,meshdatactx%nvlocal),D(NDOF,meshdatactx%nvlocal),&
          DELD(NDOF,*),RATE(NSTR,*), &
   IEN(NEN,*),SKEW(2,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   DFAULT(NDOF,*),DX(NDOF,meshdatactx%nvlocal)
!-locl
double precision :: ee(6),xl(24),dl(24)
double precision :: factor
integer          :: i,n,nn,ierr

if (DELTP.gt.EPS) then
    factor = 1d0/DELTP
else
    factor = 1d0
endif
 
do nn=1,NUMEL
    n = nn
!   Compute element geometry at end of previous time step
    call LCOORD (X,xl,IEN(1,n),nn)

    if (LGDEF.ne.0) then
        call LDISP (dl,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,6)
        call ADDSNE (DL,n,NDOF,NEN, useDXE)
        call REZONE (xl,dl,1d0)
    endif

!        Compute infinitesimal strain increment
#ifdef EXPERIMENTAL_ReducedSlipAlloc
    if(elementSlipperyID(n,115).gt.0) then
        call STRAIN (xl,IEN(1,n),DELD,DX,DFAULT, &
                     LMX(1,1,elementSlipperyID(n,115)),n,SKEW, &
                     ee,NDOF,NSD,NEN,NSTR,nn,ierr)
    endif
#else
    call STRAIN (xl,IEN(1,n),DELD,DX,DFAULT, &
       LMX(1,1,n),n,SKEW, &
       ee,NDOF,NSD,NEN,NSTR,nn,ierr)
#endif

    if (ierr.ne.0) then
        write(stderr,1) n
 1      format(1x,'trate: error in element ',I6)
        call stoper()
    endif

    do i=1,NSTR
        RATE(i,n) = ee(i)*factor
    enddo

enddo
!
return
end

#else

!-------------------------------------------------------------------------------

subroutine StrainrateOfElem(X,D,DELD,strainRate,IEN,LMF,TFAULT,DFAULT,LMX, &
   DX,SKEW,DXE, elemID)

USE TIMESTEPMODULE

#ifdef EXPERIMENTAL_ReducedSlipAlloc
USE MODELDEFINITION, only: numel,lgdef, nElemsWithSlipperyNodes
use modeldatamodule, only: elementSlipperyID
#else
USE MODELDEFINITION, only: numel,lgdef
use meshdatamodule,  only: meshdatactx
#endif

use modeltopology, only: NSD,NDOF,NSTR,NEN
use constants, only: eps, useDXE 
!
! Computes total strain rate
!
implicit none

!-pass

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer          :: LMX(NDOF,NEN,nElemsWithSlipperyNodes)
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
integer          :: LMX(NDOF,NEN,meshdatactx%neglobal)
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif

integer :: ien,lmf
double precision :: x,d,deld,skew,tfault,dfault,dx

dimension X(NSD,*),D(NDOF,*),DELD(NDOF,*), &
   IEN(NEN,*),SKEW(2,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   DFAULT(NDOF,*),DX(NDOF,*)

double precision :: strainRate(nstr)
!-locl 
double precision :: xl(24),dl(24)
double precision :: factor
integer          :: i,n,nn,ierr

integer          :: elemID

if (DELTP.gt.EPS) then
    factor = 1d0/DELTP
else
    factor = 1d0
endif  

! Compute element geometry at end of previous time step
call LCOORD (X,xl,IEN(1,elemID),elemID)

if (LGDEF.ne.0) then
    call LDISP (dl,D,IEN(1,elemID),NDOF,NEN)
    call AddFaultDisplacement (dl,elemID,TFAULT,NDOF,NEN,6)
    call ADDSNE (DL,elemID,NDOF,NEN, useDXE)
    call REZONE (xl,dl,1d0)
endif

! Compute infinitesimal strain increment
#ifdef EXPERIMENTAL_ReducedSlipAlloc
if(elementSlipperyID(n,115).gt.0) then
    call STRAIN (xl,IEN(1,elemID),DELD,DX,DFAULT, &
                 LMX(1,1,elementSlipperyID(elemID,115)),elemID,SKEW, &
                 strainRate,NDOF,NSD,NEN,NSTR,elemID,ierr)
endif
#else  
call STRAIN (xl,IEN(1,elemID),DELD,DX,DFAULT, &
             LMX(1,1,elemID),n,SKEW, &
             strainRate,NDOF,NSD,NEN,NSTR,elemID,ierr)
#endif

if (ierr.ne.0) then
     write(stderr,1) n
1    format(1x,'trate: error in element ',I6)
     call stoper()
endif

do i=1,NSTR
    strainRate(i) = strainRate(i)*factor
enddo

return 
end    

#endif

!-------------------------------------------------------------------------------
subroutine MAXWEL (FMAXW,FTOT,B,ID,NDOF)

USE TIMESTEPMODULE
USE MODELDEFINITION, only: NUMNP
!
! Adds loads from Maxwell boundary conditions to load vector
!
 implicit none
!-pass
integer          :: id,ndof
double precision :: b,fmaxw,ftot
 dimension B(*),ID(NDOF,*),FMAXW(NDOF,*),FTOT(NDOF,*)
!-local
integer :: n,i,k

 do n=1,NUMNP
do i=1,NDOF
    k = ID(i,n)
    if (k.gt.0) then
        B(k) = B(k) - (FMAXW(i,n)*DELTP)*FTOT(i,n)
    endif
enddo
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine MXWLD (IWINK,WINK,FTOT,FMAXW,DELD,NDOF)
   USE TIMESTEPMODULE
 USE MODELDEFINITION, only: NUMNP
 use constants, only: dfmin,dlogmin,dfmax,dlogmax

!
! Sets or updates load on nodes with Maxwell boundary conditions.
!
 implicit none
!-pass
integer :: ndof
integer :: iwink
double precision :: wink, ftot, deld, fmaxw
 dimension IWINK(NDOF,*),WINK(NDOF,*),FTOT(NDOF,*), &
   DELD(NDOF,*),FMAXW(NDOF,*)
!-local
integer :: i,n,ll

 do n=1,NUMNP
do i=1,NDOF
LL = IWINK(i,n)

if (LL.gt.0) then

    if (NSTEP.lt.(LL-1)) then
        FTOT(i,n) = 0d0
    else if (NSTEP.eq.(LL-1)) then
        FTOT(i,n) = -WINK(i,n)*DELD(i,n)
    else
        FTOT(i,n) = FTOT(i,n) -WINK(i,n)*DELD(i,n) &
           - (FMAXW(i,n)*DELTP)*FTOT(i,n)
    endif

    if (ABS(FTOT(i,n)).le.dfmin) then
        FTOT(i,n) = 0d0
    endif

else if (LL.lt.0) then

    if (NSTEP.eq.0) then
        FTOT(i,n) = -WINK(i,n)*DELD(i,n)
    else if (NSTEP.lt.IABS(LL)) then
        FTOT(i,n) = FTOT(i,n) -WINK(i,n)*DELD(i,n) &
           - (FMAXW(i,n)*DELTP)*FTOT(i,n)
    else if (NSTEP.ge.IABS(LL)) then
        FTOT(i,n) = FTOT(i,n) &
           - (FMAXW(i,n)*DELTP)*FTOT(i,n)
    endif

    if (ABS(FTOT(i,n)).le.dfmin) then 
        FTOT(i,n) = 0d0
    endif

else

    FTOT(i,n) = 0d0

endif

enddo
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
subroutine UNLOCK ()

#ifdef EXPERIMENTAL_CYCLIC
USE MODELDEFINITION, only: NUMNP, UseCyclicSlip, numsn
#else
USE MODELDEFINITION, only: NUMNP, numsn
#endif

use modeltopology,   only: ndof

 use meshdatamodule, only: meshdatactx
use modeldatamodule, only: modeldatactx

#ifdef SPARSE
use modelctx,        only: getrank, getsize
#endif
use debugmodule,     only: delay
use timestepmodule,  only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

! Remove forces locking a fault, except for initial applied forces. Useful
! for studying dynamic fault slip. In case of complete detachment, both delta F
! and F across the fault are removed.
! In the current version, unlocking is periodic

! Essentially, this unlocks the slippery nodes, and as such is the opposite  
! of the winklf routine above

implicit none

!-local
integer :: i, j, k, l, diffWinkTiming, iProc
integer :: myRank
integer :: itshift
integer :: nstep_minus_itshift


if (NUMSN.le.0) then
    return
endif

itshift = 0

! for every vertex, check various reasons why it
! should *not* be unlocked/locked

! ToDO, maybe nvlocal???
!do j=1,meshdatactx%nvglobal
do j=1,meshdatactx%nvlocal
#ifdef EXPERIMENTAL_CYCLIC
#ifdef ZISFAULTNORMAL
!RG 2018-11-23: time shift in periodic locking/unlocking
    do i=1,NDOF-1
        itshift = modeldatactx%IWINX(NDOF,j)
#else
    do i=1,NDOF
#endif
#else
    do i=1,NDOF
#endif
#ifdef EXPERIMENTAL_CYCLIC
        diffWinkTiming = modeldatactx%IWINX(i,j)
        if (UseCyclicSlip) then
            ! Interpretation of IWINK as periodic

            nstep_minus_itshift = NSTEP - itshift
            if (nstep_minus_itshift.lt.0) then
                nstep_minus_itshift = nstep_minus_itshift + IABS(diffWinkTiming)
            endif

            if (diffWinkTiming.lt.0) then
                ! typically locked, except when MOD(NSTEP,ABS(diffWinkTiming)) < 10
                if (MOD(nstep_minus_itshift,IABS(diffWinkTiming)).ne.0) then
!                    write(*,*) "skipping; reason 1",j,i
                    goto 100
                endif
            else if (diffWinkTiming.gt.0) then
                ! typically unlocked, except when MOD(NSTEP,diffWinkTiming) < 10
                if (MOD(nstep_minus_itshift,diffWinkTiming).ne.10) then
!                    write(*,*) "skipping; reason 2",j,i
                    goto 100
                endif
            else
                ! diffWinkTiming == 0; no differential Winkler issues
!                write(*,*) "skipping; reason 3",j,i
                goto 100
            endif
        else
            if (diffWinkTiming.gt.0) then
                ! locked from diffWinkTiming onward, so add diforc before this
                if (NSTEP.ge.diffWinkTiming) goto 100
            else if (diffWinkTiming.lt.0) then
                ! locked until and including -diffWinkTiming-1, so add diforc after this
                if (NSTEP.le.(-diffWinkTiming-1)) goto 100
            endif
        endif
#else
        if (diffWinkTiming.gt.0) then
            ! locked from diffWinkTiming onward, so add diforc before this
            if (NSTEP.ge.diffWinkTiming) goto 100
        else if (diffWinkTiming.lt.0) then
            ! locked until and including -diffWinkTiming-1, so add diforc after this
            if (NSTEP.le.(-diffWinkTiming-1)) goto 100
        endif
#endif

! The node has passed the test, ie, it is unlocked 

        k = modeldatactx%IDX(i,j)

        if (k.gt.0) then
            modeldatactx%B(k) = -modeldatactx%BTOT(k) + &
                                 modeldatactx%DIFORC(i,j)

!           write(stdout,1) j,i,-BTOT(k),DIFORC(i,j)
!    1      format(1X,'unlock: node,dof,dF,-diforce= ', &
!            2I5,2(1PG18.6))

        endif
100     continue

    enddo
enddo
!
return
end
!-------------------------------------------------------------------------------
 subroutine WINKLF (BTOT,DELD,ID,IWINK,WINK,NDOF,NUMNP,winkflag)

#ifdef EXPERIMENTAL_CYCLIC
use modelDefinition, only: UseCyclicSlip
#endif
use meshdatamodule
use modeldatamodule, only: modeldatactx
use modelctx,        only: getrank ! for debuggign
use timestepmodule,  only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

! Routine to compute Winkler restoring forces from the displacements and
! add them to the global force vector BTOT.

! Essentially, this locks the slippery nodes, and as such is the opposite
! of the unlock routine above


! winkFLAG = 0; Nodal Winkler forces
! winkFLAG = 1; Differential Winkler forces

implicit none

!-pass
integer NDOF,NUMNP,winkflag, id
double precision :: btot, deld, wink
integer :: iwink
 dimension BTOT(*),DELD(NDOF,*),ID(NDOF,*),IWINK(NDOF,*), &
   WINK(NDOF,*)
!-local
integer :: i, j, k, l, LL, iend, itshift
integer :: nstep_minus_itshift

! RG 2018-11-23
iend = NDOF
#ifdef EXPERIMENTAL_CYCLIC
#ifdef ZISFAULTNORMAL
if (winkFLAG.ne.0) iend = NDOF-1
#endif
#endif

itshift = 0

!do 400 j=1,NUMNP
!do 400 j=1,meshdatactx%nvlocal
do j=1,meshdatactx%nvglobal
    do i=1,iend
        LL = IWINK(i,j)

        ! select negative conditions, which prevent the
        ! differential Winkler forces from taking effect.

        if (winkFLAG.eq.0) then
            ! normal winkler forces
            if (LL.lt.0) then
                ! Winkler forces applied up to time LL
                if (NSTEP.ge.-LL) goto 300
            else if (LL.gt.0) then
                ! Winkler forces applied starting at time LL-1
                if (NSTEP.lt.(LL-1)) goto 300
            else
                goto 300
            endif
        else
            ! differential winkler forces
#ifdef EXPERIMENTAL_CYCLIC
            if (UseCyclicSlip) then
#ifdef ZISFAULTNORMAL
                itshift = IWINK(NDOF,j)
#endif
            nstep_minus_itshift = NSTEP - itshift
            if (nstep_minus_itshift.lt.0) then
                nstep_minus_itshift = nstep_minus_itshift + IABS(LL)
            endif
!               Interpretation of IWINK as periodic
                if (LL.lt.0) then
!                   typically locked, except when MOD(NSTEP,ABS(LL)) < 10
                    if (MOD(nstep_minus_itshift,IABS(LL)).lt.10) goto 300
                else if (LL.gt.0) then
!                   typically unlocked, except when MOD(NSTEP,LL) < 10
                    if (MOD(nstep_minus_itshift,LL).ge.10) goto 300
                else
                    goto 300
                endif
            else
!               Non-periodic IWINK
                if (LL.lt.0) then 
                    ! Winkler forces applied up to time LL   
                    if (NSTEP.ge.-LL) goto 300
                else if (LL.gt.0) then
                    ! Winkler forces applied starting at time LL-1
                    if (NSTEP.lt.(LL-1)) goto 300
                else
                    goto 300
                endif
            endif
#else
            if (LL.lt.0) then
                ! Winkler forces applied until NSTEP -LL   
                if (NSTEP.ge.-LL) goto 300
            else if (LL.gt.0) then
                ! Winkler forces applied starting at time LL-1
                if (NSTEP.lt.(LL-1)) goto 300
            else
                goto 300
            endif
#endif
        endif

        ! all the reasons that the (differential) Winkler force
        ! would not be applied had been treated. If the flow reaches here,
        ! apply the force

        ! ID both for regular and for differential Winkler forces.
        ! The argument to Winklf is either ID or IDX.

        if (winkflag.eq.0) then
            if (j.gt.size(modeldatactx%ID,2)) then
                write(*,*) "Error on thread",getrank(),": WINKLF tries to read the modeldat ID entry", j
                write(*,*) "while there are only: ", size(modeldatactx%ID,2), "nodes in the partition"
                write(*,*) "This should not have happened. Contact model support"
                stop "Exiting"
            endif
            k =  modeldatactx%ID(i,j)           
        else
            if (j.gt.size(modeldatactx%IDX,2)) then
                write(*,*) "Error on thread",getrank(),": WINKLF tries to read the modeldat IDX entry", j
                write(*,*) "while there are only: ", size(modeldatactx%IDX,2), "nodes in the partition"
                write(*,*) "This should not have happened. Contact model support"
                stop "Exiting"
            endif
            k =  modeldatactx%IDX(i,j)
        endif

        if (k.le.0 .or. k.gt.size(modeldatactx%BTOT,1)) goto 300

!!  ONDUIDELIJK WAAROM NIET DE WINK UIT HET ARGUMENT VD SUBROUTINE WORDT !GEBRUIKT
        if (winkflag.eq.0) then
            ! regular Winkler forces
               modeldatactx%BTOT(k)   = modeldatactx%BTOT(k) - &
                         modeldatactx%WINK(i,j)*modeldatactx%DELD(i,j)
        else
            ! differential Winkler forces
               modeldatactx%BTOT(k)   = modeldatactx%BTOT(k) - &
                         modeldatactx%WINX(i,j)*modeldatactx%DX(i,j)

        endif
300 continue

  enddo
enddo

return
end
!-------------------------------------------------------------------------------
subroutine StpChk (STN,MAT,PRPMAT,NSTEP)

USE TIMESTEPMODULE
USE MATERIALSMODULE
USE MODELDEFINITION
use constants,       only: dfmin,dlogmin,dfmax,dlogmax
use modeltopology,   only: nstr
use iomodule          

! calculates and outputs stability, oscillation and accuracy characteristics

implicit none
!-pass
integer :: nstep, mat
double precision :: STN,PRPMAT
dimension STN(NSTR,*),MAT(*),PRPMAT(6,*)

!-local
double precision :: E,pois,emhu,anpwr
integer          :: n, nn,m
double precision :: Emax
double precision :: Ac,Oc,St
double precision, external :: Accur, Oscill, Stabil, eigMax
double precision :: eigVal

!-init
Emax = -dfmax

if (.not.OUTFIL(1)) then
    return
endif

do nn=1,NUMEL

    n = nn
    m = MaterialIDOfElement(n)
    E     = matdatactx%PRPMAT(1,m)
    POIS  = matdatactx%PRPMAT(2,m)
    EMHU  = matdatactx%PRPMAT(3,m)
    ANPWR = matdatactx%PRPMAT(4,m)
    Eigval = EigMax (STN(1,n),E,POIS,EMHU,ANPWR,ALFAP)

    if (Eigval.gt.Emax) then
        Emax = Eigval
    endif
enddo

Ac = Accur  (ALFAP,Emax,DELTP)
Oc = Oscill (ALFAP,Emax,DELTP)
St = Stabil (ALFAP,Emax,DELTP)

write(luout(1),1) Ac,Oc,St,NSTEP,DELTP
 1   format(///1X, &
 'T I M E',3X,'S T E P',3X,'C H A R A C T E R I S T I C S'// &
 3x,'RELATIVE MODAL ERROR = ',1PG12.2// &
 3x,'OSCILLATION LIMIT  . = ',1PG12.2// &
 3x,'STABILITY RATIO  . . = ',1PG12.2// &
 3x,'NSTEP  . . . . . . . = ',1PG12.2// &
 3x,'DELTP  . . . . . . . = ',1PG25.15/)

return
end
!-------------------------------------------------------------------------------
! subroutine RESIDU (B,BRES,BTOT,STN,X,D,SKEW,IEN,LM, &
!   LMX,LMF,TFAULT,MAT,PRPMAT,DXE,OSKEW,ID,IDX,NSTR, &
!   NSD,NDOF,NEN,PSUMSQ,ACCURCY)

subroutine residu(PSUMSQ,ACCURCY)

USE MATERIALSMODULE, only: MaterialIDOfElement, &
                           matdatactx
USE MODELDEFINITION, only: lgdef, &
                           NSLSKEW, &
                           ntype, &
                           numrot, &
                           numsn, &
                           numel, numnp, &
                           NEQlocal, NEQglobal
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: modeldatactx, &
                           elementSlipperyID
#else
use modeldatamodule, only: modeldatactx
#endif
USE MESHDATAMODULE,  only: meshdatactx
USE MODELCTX,        only: getrank  ! for debugging
USE ALGEBRA,         only: DOT, &
                           clear, & 
                           equate, & 
                           vmprd, &
                           transp, &
                           formrt

use constants,       only: dfmin,dlogmin,dfmax,dlogmax,eps, zero, one,useDXE
use debugmodule,     only: lustat, itest, debug, iecho, &
                           CheckArrayForNans
use modeltopology,   only: NSTR,NSD,NDOF,NEN
use iomodule          
use timestepmodule,  only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

!
! Computes the deviation from overall force balance
!
implicit none

!-pass
!integer NSTR,NSD,NDOF,NEN
!integer IEN(NEN,*)

!double precision :: B, BRES, BTOT, STN, X, D, SKEW, tfault
!double precision :: prpmat, dxe, oskew
!integer          :: LM, LMX, LMF, id, idx, mat
!dimension B(*),BRES(*),BTOT(*),STN(NSTR,*),X(NSD,*),D(NDOF,*), &
!   SKEW(2,*),LM(NDOF,NEN,*),LMX(NDOF,NEN,*), &
!   TFAULT(NDOF,*),LMF(NDOF,NEN,*),PRPMAT(6,*),MAT(*), &
! DXE(NDOF,NEN,*),OSKEW(2,*),ID(NDOF,*),IDX(NDOF,*)

double precision :: accurcy

!-locl
integer          :: NEE
logical quad,firstentry
double precision :: thick
integer i,j,l,lu,m,n
double precision :: pl(24),xl(24),dl(24),bl(3),rot(3,3),stnl(6)
double precision :: psumsq, bsumsq
integer :: ierr
save firstentry
data firstentry/.true./
!-init

psumsq = ZERO
bsumsq = ZERO
NEE    = NDOF * NEN

pl = 0d0
xl = 0d0
dl = 0d0
bl = 0d0
rot = 0d0
stnl = 0d0

modeldatactx%B = 0d0

! only do contributions to elements within this partition
do n=1,meshdatactx%nelocal

    if (iecho.eq.9) then
    write(*,*) "------------ element: ", n, "----------------"
    endif

    m = MaterialIDOfElement(n)

    thick = matdatactx%PRPMAT(6,m)
    quad = (meshdatactx%IEN(3,n).ne.meshdatactx%IEN(4,n))

    call CLEAR (pl,NEE,"pl")
    call LCOORD (meshdatactx%X,xl,meshdatactx%IEN(1,n),n)

    if (LGDEF.ne.0) then
        call LDISP  (dl,modeldatactx%D,meshdatactx%IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (dl,n,modeldatactx%TFAULT,NDOF,NEN,7)
        call ADDSNE (dl,n,NDOF,NEN, useDXE)
        call REZONE (xl,dl,ONE)
    endif

    call EQUATE (stnl,modeldatactx%STN(1,n),NSTR)

    call EFORCE (stnl,pl,xl,thick,quad,ierr)

    if (ierr.ne.0) then
        write(stderr,1) n
 1         format(1x,'residu: error in element ',I5)
        call stoper()
    endif

    if (NUMROT.ne.0) then
        if (NSLSKEW.gt.0) then
            CALL RPFORC (pl,modeldatactx%OSKEW,meshdatactx%IEN(1,n),NDOF,NEN)
        else
            CALL RPFORC (pl,modeldatactx%SKEW,meshdatactx%IEN(1,n),NDOF,NEN)
        endif
    endif

!*******************************************************************************
!***** PSUMSQ DOES NOT WORK IN PARALLEL:
!*******************************************************************************
    do i=1,NEE
        psumsq = psumsq + pl(i)*pl(i)
    enddo

    call ADDFOR (modeldatactx%B,pl,modeldatactx%LM(1,1,n),NEE,n,17)

    if (NUMSN.gt.0) then

        ! call ADDFOR (B,     pl,LMX(1,1,n),NEE,n,18)
        ! abs, because both pl and LMX have the correct sign. They would cancel
        ! each other out. The side of LMX is considered leading, because it is
        ! directly from the modeldata structure.
        ! call ADDFOR (modeldatactx%B,abs(pl),modeldatactx%LMX(1,1,n),NEE,n,18)
        ! write(*,*) "rank", getrank(),"B before slippery addfor", modeldatactx%B

#ifdef EXPERIMENTAL_ReducedSlipAlloc
        if(elementSlipperyID(n,116).gt.0) then
            call ADDFOR (modeldatactx%B,PL,modeldatactx%LMX(1,1,elementSlipperyID(n,116)),NEE,n,18)
        else
            ! this element has no slippery nodes. Skip it.
        endif
#else
        call ADDFOR (modeldatactx%B,PL,modeldatactx%LMX(1,1,n),NEE,n,18)
#endif

    endif
enddo ! loop over the local vertices

if (debug) then
    if (OUTFIL(1)) then
        lu = luout(1)
    else
        lu = stdout
    endif
endif

#ifdef SPARSE
!   write(*,*) "calling syncloadvec 4" 
!******************************************************************'
!***** OBSOLETE???? REASON?
!******************************************************************'
call SyncLoadVector(0)

#endif

do i=1,NEQlocal
    modeldatactx%BRES(i) = modeldatactx%BTOT(i) - modeldatactx%B(i)
    if (debug) then
        write(lu,2) i,modeldatactx%BTOT(i),modeldatactx%B(i),modeldatactx%BRES(i)
 2       format(1x,'NEQ=',I5,': BTOT=',1PE20.12,'  BSTN=',1PE20.12, &
     '  BRES=',1PE20.12)
    endif
enddo

! if (firstentry) then
  ! write(*,*) '******************************************************************'
  ! write(*,*) '!***** RESIDU: MISSING SCATTER OF BRES????'
  ! write(*,*)' !******************************************************************'
  ! firstentry=.false.
! endif
bsumsq = DOT(modeldatactx%BRES, modeldatactx%BRES, NEQlocal)

!***************************************************
!***** CALCULATION OF ACCURCY WILL FAIL IN PARALLEL
!***************************************************
ACCURCY = ZERO

if (psumsq.gt.dfmin) then
    accurcy = SQRT(bsumsq/psumsq)
endif

if (NUMROT.eq.0 .or. NSLSKEW.eq.0 .or. NSTEP.eq.0) then
    goto 200
endif

!*****    rotate BRES to updated coordinate system
do 100 i=1,NUMNP

    if (modeldatactx%OSKEW(1,i).eq.ZERO .and. &
        modeldatactx%OSKEW(2,i).eq.ZERO) then
        goto 100
    endif

    call CLEAR (bl,NDOF,"bl")

    do j=1,NDOF
        l = modeldatactx%ID(j,i)

        if (l.gt.0) then
            bl(j) = modeldatactx%BRES(l)
        else
            l = modeldatactx%IDX(j,i)
            if (l.gt.0) then
                bl(j) = modeldatactx%BRES(l)
            endif
        endif
    enddo

    call FORMRT (modeldatactx%OSKEW(:,i),rot,NDOF)
    call VMPRD  (rot,bl,pl,NDOF,3)
    call FORMRT (modeldatactx%SKEW(:,i),rot,NDOF)
    call TRANSP (rot,3)
    call VMPRD  (rot,pl,bl,NDOF,3)

    do j=1,NDOF
        l = modeldatactx%ID(j,i)
        if (l.gt.0) then
            modeldatactx%BRES(l) = bl(j)
        else
            l = modeldatactx%IDX(j,i)
            if (l.gt.0) then
                modeldatactx%BRES(l) = bl(j)
            endif
        endif
    enddo

100 continue

200 continue
if (OUTFIL(1)) then
    if (accurcy.lt.dfmin) then
        write(luout(1),3)
 3      format(/1X,'RESIDUAL FORCE AT EQUILIBRIUM'/)
    else
        write(luout(1),4) accurcy
 4      format(/1X,'RESIDUAL FORCE DEVIATION FROM EQUILIBRIUM ', &
         '= ',1PG12.3/)
   endif
endif

psumsq = bsumsq

return
end

!-----------------------------------------------------------------------
subroutine WRBRES (BRES,X,D,SKEW,IEN,LM,LMX,LMF, &
                   TFAULT,DXE,NSD,NDOF,NEN,file)!,filex,filey)

USE MODELDEFINITION
USE ALGEBRA
use meshdatamodule, only: meshdatactx
use constants,      only: useDXE

#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: elementSlipperyID
#endif
!
! Outputs force vector
!
implicit none
!-pass
character(len=*) file!,filex,filey
integer NSD,NDOF,NEN

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer          :: LMX(NDOF,NEN,nElemsWithSlipperyNodes)
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
integer          :: LMX(NDOF,NEN,meshdatactx%neglobal)
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif

integer :: IEN,LM,LMF
double precision :: BRES,X,D,TFAULT,SKEW
dimension BRES(*),X(NSD,*),D(NDOF,*),IEN(NEN,*),LM(NDOF,NEN,*), &
         TFAULT(NDOF,*),LMF(NDOF,NEN,*), SKEW(2,*)

!-locl
double precision :: xl(24),dl(24),bl(24)
integer :: lien(8),llien(8)
integer :: i,j,k,n,nn,nbl,keq,idof,jnode,nelem
integer :: lu
integer, external :: nextlu

! integer luxx, luyy

if (NTYPE.ge.5) then ! f3d
    lu = nextlu(0)
    call openf(lu,file,'unknown')
    write(lu,1) NDOF,NEN*NUMEL
 1       format('object 1 class array type float rank 1 shape ',I1, &
     ' items ',I5,' ASCII data follows')

    do n=1,NUMEL
        call LCOORD (X,XL,IEN(1,n),n)

        if (LGDEF.ne.0) then
            call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,8)
            call ADDSNE (DL,n,NDOF,NEN,useDXE)
            call REZONE (XL,DL,1d0)
        endif
        do j=1,NEN
            write(lu,2) (xl(i+(j-1)*NDOF),i=1,NDOF)
 2                format(1PG14.6,1X,1PG14.6,1X,1PG14.6)
        enddo
    enddo

    write(lu,3) NEN,NUMEL
 3       format(/'object 2 class array type int rank 1 shape ',I1, &
     ' items ',I5,1X,'ASCII data follows')

    do n=1,NUMEL
        do j=1,NEN
            lien(j) = (n-1)*NEN+j
        enddo

        llien(1)=lien(1)-1
        llien(2)=lien(2)-1
        llien(3)=lien(5)-1
        llien(4)=lien(6)-1
        llien(5)=lien(4)-1
        llien(6)=lien(3)-1
        llien(7)=lien(8)-1
        llien(8)=lien(7)-1
    
        write(lu,4) (llien(j),j=1,8)
 4           format(8(1X,I6))
    enddo

    write(lu,5) NDOF,NEN*NUMEL
 5       format('attribute "element type" string "cubes"'/ &
     'attribute "ref" string "positions"'//'object 3 class ', &
     'array type float rank 1 shape ',I1,' items ',I5, &
     ' ASCII data follows')

    do n=1,NUMEL
        call CLEAR (bl,NEN*NDOF)
        do j=1,NEN
            do i=1,NDOF
                k = LM(i,j,n)
                if (k.ne.0) then
                    bl(i+(j-1)*NDOF) = BRES(k)
                else
                    if (NUMSLP.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc
                        if(elementSlipperyID(n,117).gt.0) then
                            k = IABS(LMX(i,j,elementSlipperyID(n,117)))
                        else
                            k = 0
                        endif
#else
                        k = IABS(LMX(i,j,n))
#endif
                        if (k.ne.0) then
                            bl(i+(j-1)*NDOF) = BRES(k)
                        endif
                    endif
                endif
            enddo
        enddo
        write(lu,6) ((bl(i+(j-1)*NDOF),i=1,NDOF),j=1,NEN)
 6           format(24(1X,1PG14.6))
    enddo

    write(lu,7)
 7       format('attribute "dep" string "positions"'// &
     'object "irregular positions irregular',1x, &
     'connections" class field'/'component "positions"',1x, &
     'value 1'/'component "connections" value 2'/ &
     'component "data" value 3'/'end')
    call closef(lu)
!
elseif (NTYPE.eq.1) then
!   only triangles
    do nelem=1,NUMEL
        if (IEN(3,nelem).ne.IEN(4,nelem)) then
            return
        endif
    enddo

    lu = nextlu(0)
    call openf(lu,file,'unknown')

!   output per element: node coordinates and magnitude of residual force per node
    do nelem=1,NUMEL
        call LCOORD (X,XL,IEN(1,nelem),nelem)
        if (LGDEF.ne.0) then
            call LDISP  (DL,D,IEN(1,nelem),NDOF,NEN)
            call AddFaultDisplacement (DL,nelem,TFAULT,NDOF,NEN,9)
            call ADDSNE (DL,nelem,NDOF,NEN,useDXE)
            call REZONE (XL,DL,1d0)
        endif

        call CLEAR (bl,NEN,"bl")
        do jnode=1,3
            do idof=1,NDOF
                keq = LM(idof,jnode,nelem)
                if (keq.ne.0) then
                    bl(jnode) = bl(jnode) + BRES(keq)**2
                else
                    if (NUMSLP.gt.0) then

#ifdef EXPERIMENTAL_ReducedSlipAlloc
                        if(elementSlipperyID(nelem,118).gt.0) then
                            keq = IABS(LMX(idof,jnode,elementSlipperyID(nelem,118)))
                        else
                            keq = 0
                        endif
#else
                        keq = IABS(LMX(idof,jnode,nelem))
#endif
                        if (keq.ne.0) then
                            bl(jnode) = bl(jnode) + BRES(keq)**2
                        endif
                    endif
                endif
                bl(jnode) = sqrt(bl(jnode))
            enddo
            write(lu,*) xl(1+(jnode-1)*NSD),xl(2+(jnode-1)*NSD),bl(jnode)
        enddo
    enddo
    call closef(lu)

endif

return
end
!-------------------------------------------------------------------------------

 subroutine RESTIO (D,DXE,STN,STR,BTOT,T,FTOT,VPSTR, &
 NDOF,NEN,NSTR,IFLAG)
   USE MESHDATAMODULE
   USE MODELDATAMODULE
   USE MODELDEFINITION
   USE MODELCTX
USE FILEMODULE,      only: workpath
use debugmodule,     only: xit, iecho
use iomodule          
use timestepmodule,  only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

!
! reads or writes restart file
!
 implicit none
!-pass




integer :: ndof, nstr, nen, iflag
double precision :: d, dxe, stn, btot, t, str, ftot, vpstr
 dimension D(NDOF,*),DXE(NDOF,NEN,*),STN(NSTR,*),BTOT(*),T(*), &
   STR(NSTR,*),FTOT(NDOF,*),VPSTR(NSTR,*)
!-locl
   integer lu,ios,n,i, j,ERROR
   character(len=255) restartfile
integer, external :: nextlu, idig

integer :: nranks, thisRank
!
thisRank = getrank()
nRanks = getsize()


if (iecho.ne.0) then
    write(*,*) "restio in visco.f says: workpath=", trim(workpath)
endif

if      (nRanks.ge.10000) then ! just in case, although this many log files will be unwise.
    write(*,*) "You are about to create over 10k restart files"
    write(*,*) "This not very wise"
    write(*,*) "You should talk to a professional about this..."
    write(*,*) "Or GTecton is broken. Then please talk to model support first."
    stop "Aborting run"
else if (nRanks.ge.1000) then ! just in case, although this many log files will be unwise.

    if (thisrank.ge.1000) then
        write (restartfile, "(A,I4)") trim(workpath)//"/RESTART.DAT.",getrank()
    else if (thisrank.ge.100) then
        write (restartfile, "(A,I3)") trim(workpath)//"/RESTART.DAT.0",getrank()
    else if (thisrank.ge.10) then
        write (restartfile, "(A,I2)") trim(workpath)//"/RESTART.DAT.00",getrank()
    else
        write (restartfile, "(A,I1)") trim(workpath)//"/RESTART.DAT.000",getrank()
    endif

else if (nRanks.ge.100) then

    if (thisrank.ge.100) then
        write (restartfile, "(A,I3)") trim(workpath)//"/RESTART.DAT.",getrank()
    else if (thisrank.ge.10) then
        write (restartfile, "(A,I2)") trim(workpath)//"/RESTART.DAT.0",getrank()
    else
        write (restartfile, "(A,I1)") trim(workpath)//"/RESTART.DAT.00",getrank()
    endif

else if (nRanks.ge.10) then

    if (thisrank.ge.10) then
        write (restartfile, "(A,I2)") trim(workpath)//"/RESTART.DAT.",getrank()
    else
        write (restartfile, "(A,I1)") trim(workpath)//"/RESTART.DAT.0",getrank()
    endif

else

    write (restartfile, "(A,I1)") trim(workpath)//"/RESTART.DAT.",getrank()

endif


!   write(restartfile,21) trim(workpath)//"/RESTART.DAT"//getrank()
! 21   format(A,I1)


 lu = nextlu(luout(2)+10)
 if (.true.) then
if (IFLAG.eq.0) then
    call openf(lu,trim(restartfile),'unknown')
    write(lu,1,err=100,iostat=ios) NSTEP
1            format(I10)
    if (MODE.eq.5 .or. MODE.eq.6) then
        write(lu,2,err=100,iostat=ios) (T(n),n=1,NUMNP)
    endif
    if (MODE.eq.2 .or. MODE.eq.6) then
    write(lu,2,err=100,iostat=ios) &
           ((D(i,n),i=1,NDOF),n=1,NUMNP)


    if (NUMSE.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc
        write(lu,2,err=100,iostat=ios) &
           (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,nElemsWithSlipperyNodes)
#else
        write(lu,2,err=100,iostat=ios) &
           (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
#endif
    endif


    write(lu,2,err=100,iostat=ios) &
           ((STN(i,n),i=1,NSTR),n=1,NUMEL)
    write(lu,2,err=100,iostat=ios) &
           ((STR(i,n),i=1,NSTR),n=1,NUMEL)
            if (IDIG(ICVIS,4,4).ne.0) &
             write(lu,2,err=100,iostat=ios) &
             ((VPSTR(i,n),i=1,NSTR),n=1,NUMEL)
    write(lu,2,err=100,iostat=ios) (BTOT(i),i=1,NEQ)
    if (NMAXW.gt.0) write(lu,2,err=100,iostat=ios) &
           ((FTOT(i,n),i=1,NDOF),n=1,NUMNP)
    endif
2            format(5E28.18)
else
    call openf(lu,trim(restartfile),'old')
    read(lu,1,err=200,iostat=ios) NSTEP
    if (MODE.eq.5 .or. MODE.eq.6) then
       read(lu,2,err=200,iostat=ios) (T(n),n=1,NUMNP)
   endif
    if (MODE.eq.2 .or. MODE.eq.6) then
    read(lu,2,err=200,iostat=ios) &
           ((D(i,n),i=1,NDOF),n=1,NUMNP)

    if (NUMSE.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc
        read(lu,2,err=200,iostat=ios) &
           (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,nElemsWithSlipperyNodes)
#else 
        read(lu,2,err=200,iostat=ios) &
           (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
#endif
    endif

    read(lu,2,err=200,iostat=ios) &
           ((STN(i,n),i=1,NSTR),n=1,NUMEL)
    read(lu,2,err=200,iostat=ios) &
           ((STR(i,n),i=1,NSTR),n=1,NUMEL)
            if (IDIG(ICVIS,4,4).ne.0) &
             read(lu,2,err=200,iostat=ios) &
             ((VPSTR(i,n),i=1,NSTR),n=1,NUMEL)
    read(lu,2,err=200,iostat=ios) (BTOT(i),i=1,NEQ)
    if (NMAXW.gt.0) read(lu,2,err=200,iostat=ios) &
           ((FTOT(i,n),i=1,NDOF),n=1,NUMNP)
    endif
    write(stdout,3) NSTEP
 3            format(1x,'EOF trim(workpath)//"/RESTART.DAT": NSTEP = ',I9)
endif
   else if (getrank().gt.1) then
  if (IFLAG.eq.0) then
!              call mpi_writerestart(meshdatactx, modeldatactx,
!     .          mode, nstep)
  else
!              call mpi_readrestart(meshdatactx, modeldatactx,
!     .          mode, nstep)
  endif
 else
if (IFLAG.eq.0) then
    call fbopen(lu,trim(workpath)//"/RESTART.DAT",'unknown')
    write(lu,err=100,iostat=ios) NSTEP
    if (MODE.eq.5 .or. MODE.eq.6) &
       write(lu,err=100,iostat=ios) (T(n),n=1,NUMNP)
    if (MODE.eq.2 .or. MODE.eq.6) then
    write(lu,err=100,iostat=ios) &
           ((D(i,n),i=1,NDOF),n=1,NUMNP)

    if (NUMSE.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc
        write(lu,err=100,iostat=ios) &
           (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,nElemsWithSlipperyNodes)
#else
        write(lu,err=100,iostat=ios) &
           (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
#endif
    endif


    write(lu,err=100,iostat=ios) &
           ((STN(i,n),i=1,NSTR),n=1,NUMEL)
    write(lu,err=100,iostat=ios) &
           ((STR(i,n),i=1,NSTR),n=1,NUMEL)
            if (IDIG(ICVIS,4,4).ne.0) &
             write(lu,err=100,iostat=ios) &
             ((VPSTR(i,n),i=1,NSTR),n=1,NUMEL)
    write(lu,err=100,iostat=ios) (BTOT(i),i=1,NEQ)
    if (NMAXW.gt.0) write(lu,err=100,iostat=ios) &
           ((FTOT(i,n),i=1,NDOF),n=1,NUMNP)
    endif
else
    call fbopen(lu,'RESTART.DAT','old')
    read(lu,err=200,iostat=ios) NSTEP
    if (MODE.eq.5 .or. MODE.eq.6) &
       read(lu,err=200,iostat=ios) (T(n),n=1,NUMNP)
    if (MODE.eq.2 .or. MODE.eq.6) then
    read(lu,err=200,iostat=ios) &
           ((D(i,n),i=1,NDOF),n=1,NUMNP)


    if (NUMSE.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc
        read(lu,err=200,iostat=ios) &
           (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,nElemsWithSlipperyNodes)
#else
        read(lu,err=200,iostat=ios) &
           (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
#endif
    endif

    read(lu,err=200,iostat=ios) &
           ((STN(i,n),i=1,NSTR),n=1,NUMEL)
    read(lu,err=200,iostat=ios) &
           ((STR(i,n),i=1,NSTR),n=1,NUMEL)
            if (IDIG(ICVIS,4,4).ne.0) &
             read(lu,err=200,iostat=ios) &
             ((VPSTR(i,n),i=1,NSTR),n=1,NUMEL)
    read(lu,err=200,iostat=ios) (BTOT(i),i=1,NEQ)
    if (NMAXW.gt.0) read(lu,err=200,iostat=ios) &
           ((FTOT(i,n),i=1,NDOF),n=1,NUMNP)
    endif
    write(stdout,3) NSTEP
endif
 endif
 call closef(lu)
!
 return
!
100    write(stderr,101) ios
  101   format(1x,'RESTIO: write error ',I5)
 call xit(1," ")
200    write(stderr,201) ios
  201   format(1x,'RESTIO: read error ',I5)
 call xit(1," ")
 end
