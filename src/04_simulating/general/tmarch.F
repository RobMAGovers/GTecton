subroutine tmarch (TDT,NSD,NDOF,NEN,NSTEP)

   USE MESHDATAMODULE
   USE MODELDATAMODULE
   USE MODELDEFINITION
   USE MATERIALSMODULE
   USE TIMESTEPMODULE,  only: itcalc, timestepdatactx, outputcontroldatactx, &
                              itprt, maxtit, ntcalc, ntprt, talfp, output, &
                              ITWORK,IMWORK,IMTWRK,refact

   USE MODELCTX
   USE ALGEBRA
   use debugmodule,     only: lustat, itest, debug, iecho, ParallelLog
use iomodule,           only: stdout, stderr, outfil, luout

#ifdef SPARSE
USE AOMODULE
#endif

!  TDT = time belonging to this timestep, in seconds.



!-------------------------------------------------------------------------------
!     Subroutine to perform "thermal time stepping"
!-------------------------------------------------------------------------------

implicit none
!-pass
   integer NSD,NDOF,NEN
!    dimension ICALC(*),ITPRINT(*)

!-locl

double precision, save :: STEPLAST,ALFLST
integer,          save :: NENTRY
!      logical DoOut, output
   logical DoOut

double precision       :: tdt
integer                :: nstep, its, nd
integer                :: ierr

!-init



!        write(*,*) 'rank', getrank(),'entered timestepping with', TDT


   data STEPLAST,ALFLST,NENTRY /0d0,0d0,0/
   IMTWRK = 0
!
!     Is MODE set to thermal time stepping?
if (MODE.lt.5) then
    go to 1000
endif

!write(*,*) "check for output with icalc: ", timestepdatactx%ICALC,NTCALC,ITCALC,NSTEP

!     Calculation at current time step required?
if (.not.(output(timestepdatactx%ICALC, NTCALC,ITCALC,NSTEP))) then
    ! there is no caluclation here, so we cannot write output either.
!    write(*,*) "skipping temp calc here"
    goto 1000
else
!    write(*,*) "yay, we have time calc here"
endif
!
!     Output at current time step required?

DoOut  = output (outputcontroldatactx%ITPRINT, NTPRT,ITPRT,NSTEP)

!write(*,*) "computed doout: ", doout


if (itest.eq.2) then
    DoOut = .true.
endif
!
   if (iecho.eq.1) write(stdout,10)
10 format(1x,'Thermal calculation')
!
   call EQUATE (modeldatactx%TLST,modeldatactx%T,NUMNP)

 if (NUMSLP+NUMFN.gt.0.and.IDIFT.eq.1) then
call EQUATE (modeldatactx%TLSTX,modeldatactx%TX,NUMNP)
 endif
!
!     Add thermal anomalies

   call TEMADD (meshdatactx%X,modeldatactx%D, modeldatactx%T, &
modeldatactx%TLST,modeldatactx%ITANOM, modeldatactx%TANOM, &
NSTEP,NUMNP,NSD,NDOF)
!
!     set flag for thermo-mechanical coupling
   if (MODE.eq.6) IMTWRK = 1
!
!     determine whether A needs to be (re)factored
   ITWORK = 0
   if (ABS(STEPLAST-TDT).gt.1E-6.or.ABS(ALFLST-TALFP).gt.1E-6 .or. &
 (LGDEF.ne.0.and.MODE.eq.6) .or. MOD(NENTRY,MAXTIT).eq.0 .or. &
 IADVEC.ne.0) ITWORK = 1
!
if (IADVEC.ne.0) then
!     geometry changes from advection

if (iecho.eq.8) then
    call ParallelLog("tmarch", "calling load")
endif


  call LOAD (modeldatactx%ID,modeldatactx%IBOND, &
   modeldatactx%BOND,modeldatactx%D,modeldatactx%DELD, &
   modeldatactx%B,meshdatactx%X,NUMNP,NEQ,1)
  call UPDATE (modeldatactx%D,modeldatactx%DELD,NDOF*NUMNP)

!        write(*,*) 'timestep writing: ', DoOut
  if (DoOut) then
        call PRINTD (modeldatactx%D,NUMNP,0)
  endif
   endif
!
#ifdef SPARSE
   call CLEAR (modeldatactx%AA,NTEQlocal,"modeldatactx%AA")
#endif
!     Set stiffness matrix
   call FORMTA (modeldatactx%AA,itdiag,meshdatactx%X, &
meshdatactx%IEN,modeldatactx%LMT,modeldatactx%LMTX, &
meshdatactx%mat,matdatactx%PRPMAT, matdatactx%PRPTEM, &
modeldatactx%BT,modeldatactx%T,modeldatactx%TX, &
modeldatactx%D,modeldatactx%LMF,modeldatactx%TFAULT, &
modeldatactx%SKEW,modeldatactx%DXE,TDT,NSD,NEN,NDOF,ITWORK)

!     Add normal boundary flux contributions to load vector

   if (NFLX.gt.0) then
  call ADDFLX (modeldatactx%BT,meshdatactx%IEN,meshdatactx%X, &
   modeldatactx%D,meshdatactx%mat,matdatactx%PRPMAT, &
   modeldatactx%IFLS,modeldatactx%IFLX,modeldatactx%BFLX, &
   modeldatactx%LMT,modeldatactx%LMTX,modeldatactx%LMF, &
   modeldatactx%TFAULT,modeldatactx%SKEW,modeldatactx%DXE)
   endif

#ifndef SPARSE

 if (ITWORK.ne.0) then
if (iecho.eq.1) write(stdout,20)
20        format(1x,'Factoring thermal stiffness matrix')
call LUSOLVE (modeldatactx%AA,modeldatactx%BT,itdiag,NTEQ,.TRUE., &
       .FALSE.)
nd = 0
if (IDIFT.eq.1) nd = NUMSLP+NUMFN

call PIVOTS (modeldatactx%AA,itdiag,modeldatactx%IDT, &
     modeldatactx%IDTX, 1,NTEQ,nd,0,NUMNP,MODE)
call CONDNR (modeldatactx%AA,itdiag,NTEQ,ierr)
 endif
#endif
!

#ifdef SPARSE
   if (ITWORK.ne.0) then
  if (iecho.eq.1) write(stdout,30)
30     format(1x,'Inverting thermal stiffness matrix')
  call SyncLoadVector(1)
  call BuildMat(1) ! Assemble matrix and prepare for solution.
  call InitSolve(1)

!********************************************************************
!                       solving 
!write(*,*) "solving a temperature timestep"
!********************************************************************
  call ITSOLVE (modeldatactx%BT,modeldatactx%AA,its,1)
!write(*,*) "finished solving a temperature timestep"

  call EQUATE (modeldatactx%BT,modeldatactx%AA,NTEQ)
!********************************************************************
!********************************************************************
  if (iecho.gt.0) write(stdout,40) its
40       format(1X,'Sparse solution of conduction equation required', &
     1X,I3,' iterations')
   endif
#else
!     Backsubstitute inverted stiffness matrix
   if (iecho.eq.1) write(stdout,30)
30 format(1x,'Computing temperatures')
   call LUSOLVE (modeldatactx%AA,modeldatactx%BT,itdiag,NTEQ,.FALSE.,.TRUE.)
#endif
!
!     localize temperatures in T and TX
    call TEMPS(modeldatactx%BT,modeldatactx%T,modeldatactx%TX, &
          modeldatactx%IDT,modeldatactx%IDTX)

#ifdef SPARSE

!      call scattertemperatures(meshdatactx,modeldatactx,
!     .                         getrank(),modeldatactx%T,FILE_outputf,0)
   call scatterONEdata(meshdatactx,modeldatactx, &
           getrank(),modeldatactx%T,FILE_outputf,0)
   if (NUMFNglobal+NUMSLPglobal.gt.0.and.IDIFT.eq.1) &
   call scatterONEdata(meshdatactx,modeldatactx, &
           getrank(),modeldatactx%TX,FILE_outputf,0)
#endif



!     output requested?
if (DoOut) then

    !     print temperatures
    call PRINTT(modeldatactx%T,NUMNP)

! differential temperatures perhaps useful later.
!     print differential temperatures
!    call PRINTTX (modeldatactx%TX,modeldatactx%IDSLP, &
!                  modeldatactx%NFAULT)

!     compute the fluxes in each element and print them
    call EFLUX (meshdatactx%X,meshdatactx%IEN,meshdatactx%mat, &
                 modeldatactx%T,modeldatactx%TX,matdatactx%PRPTEM, &
                modeldatactx%FLUX,modeldatactx%D,modeldatactx%LMF, &
                modeldatactx%TFAULT, &
                 modeldatactx%LMTX,modeldatactx%SKEW,modeldatactx%DXE, &
                 NSD,NEN,NDOF)

    ! and print the fluxes
    call PRINTF(modeldatactx%FLUX,NSD,NUMEL)

    ITPRT = MIN0(ITPRT+1,NTPRT)

    if (OUTFIL(1)) call flushf(luout(1))
    if (OUTFIL(2)) call flushf(luout(2))
endif

   ITCALC   = MIN0(ITCALC+1,NTCALC)
   STEPLAST = TDT
   ALFLST   = TALFP
   TDT      = 0d0
   NENTRY   = NENTRY + 1
1000  return
   end
!-----------------------------------------------------------------------
subroutine ADDLAB (SC,SK,AL,BL,TDT,NEN)

use constants, only: one
use timestepmodule, only: talfp

! Routine creating local AA and AB arrays

implicit none
!-pass
integer NEN
double precision :: SC(NEN,NEN),SK(NEN,NEN),AL(NEN,NEN),BL(NEN,NEN)
double precision :: TDT
!-local
integer :: i,j
!-init

!write(*,*) "addlab has", tdt, talfp
!write(*,*) "addlab has SC", SC

do j=1,NEN
    do i=1,j
!        write(*,*) "Addlab adding to A", i, j, SC(i,j), TDT
!        write(*,*) "Addlab adding to A", i, j, SK(i,j)*TALFP, SC(i,j)/TDT
!        write(*,*) "Addlab adding to B", i, j, SC(i,j)/TDT, SK(i,j)*(one-TALFP)

        AL(i,j) = SK(i,j)*TALFP + SC(i,j)/TDT
        BL(i,j) = SC(i,j)/TDT - SK(i,j)*(one-TALFP)
    enddo
enddo

return
end
!-----------------------------------------------------------------------
subroutine FORMTA (AA,ITDIAG,X,IEN,LMT,LMTX,MAT,PRPMAT,PRPTEM, &
                   BT,T,TX,D,LMF,TFAULT,SKEW,DXE,TDT, &
                   NSD,NEN,NDOF,ITWORK)

!USE TIMESTEPMODULE
USE MATERIALSMODULE
USE MODELDEFINITION
USE ALGEBRA
use debugmodule,     only: iecho
use constants,       only: one, useDXE
use meshdatamodule,  only: meshdatactx
#ifdef SPARSE
use spetscmodule,    only: T_A
#endif
use iomodule,        only: stderr, stdout

!use modeltopology, only: ndof, nen

! Reform the thermal stiffness matrix


implicit none
!-pass


integer :: NSD,NEN,ITWORK,NDOF

#ifdef EXPERIMENTAL_ReducedSlipAlloc
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif

double precision :: tdt

integer          :: itdiag, ien, lmt, lmtx, mat, lmf
double precision :: AA, x, prpmat, prptem, bt, t, tx, d, tfault, skew

 dimension AA(*),ITDIAG(*),X(NSD,*),IEN(NEN,*), &
 LMT(NEN,*),LMTX(NEN,*),MAT(*),PRPMAT(6,*),PRPTEM(6,*), &
 BT(*),T(*),TX(*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*)

#ifdef SPARSE

!#include "petsc/finclude/petscvecdef.h"
!#include "petsc/finclude/petscmatdef.h"
!#include "petsc/finclude/petsckspdef.h"

#endif
!-locl
 logical QUAD
 integer nn,n,m,i
double precision :: SK(64),SC(64),TL(8),XL(24),AL(64),BL(64),PL(8),DL(24), COND(3)
double precision :: heat, thick, texp, rhocp
integer          :: ierr
double precision :: T1
double precision, external :: Tcentr
!-init
 ierr = 0

SK = 0d0
SC = 0d0
TL = 0d0
XL = 0d0
AL = 0d0
BL = 0d0
PL = 0d0
DL = 0d0
COND = 0d0

!
if (ITWORK.ne.0) then
    if (iecho.eq.1) then
        write(stdout,10)
10        format(1x,'Reforming thermal stiffness matrix')
    endif
#ifdef SPARSE
call MatZeroEntries(T_A,ierr)
#else
call CLEAR(AA,NAA,"AA")
#endif
endif

call CLEAR(BT,NTEQext,"BT")
!
do nn=1,NUMEL

!    write(*,*) "******** adding stuff for element ", nn, "************"

    n    = nn

    m = MaterialIDOfElement(n)

    do i=1,NSD
        COND(i) = PRPTEM(3+i,m)
    enddo

THICK = PRPMAT(6,m)
HEAT  = PRPTEM(1,m)
TEXP  = PRPTEM(3,m)
QUAD = (IEN(3,n).ne.IEN(4,n))
!
!        clear local stiffness matrices
call CLEAR(SK,NEN*NEN,"SK")
call CLEAR(SC,NEN*NEN,"SC")
call CLEAR(AL,NEN*NEN,"AL")
call CLEAR(BL,NEN*NEN,"BL")
!
!        localize (deformed) element coordinates
call LCOORD (X,XL,IEN(1,n),nn)
if (LGDEF.ne.0) then
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
    call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,10)
    call ADDSNE (DL,n,NDOF,NEN,useDXE)
    call REZONE (XL,DL,1d0)
endif
!
!        Set Rho*Cp
T1 = Tcentr (T,TX,XL,IEN(1,n),LMTX(1,n),ierr)
    if (ierr.ne.0) then
        write(stderr,1) n
        call stoper()
    endif

! write(*,*) "computing rhocp ", PRPMAT(5,m), PRPTEM(2,m), one, t1, texp


    RHOCP = ABS(PRPMAT(5,m)* PRPTEM(2,m)*(ONE-T1*TEXP))

! write(*,*) "found rhocp", rhocp

!
!        construct the local conductivity stiffness matrix SK
call TSTIFF (XL,QUAD,SK,COND,THICK,ierr)

!write(*,*) "Made conduct matrix ", SK(1:11)


if (ierr.ne.0) then
    write(stderr,1) n
 1        format(1x,'formta: error in element ',I5)
    call stoper()
endif
!
!        construct the local heat capacity stiffness matrix SC
call CPSTIFF (XL,QUAD,SC,RHOCP,THICK,ierr)
!write(*,*) "Made capacity matrix ", SC(1:11)


if (ierr.ne.0) then
    write(stderr,1) n
    call stoper()
endif


!write(*,*) "calling addlab with TDT: ", TDT, talfp
!        add SC and SK into AL and BL
call ADDLAB (SC,SK,AL,BL,TDT,NEN)

!        write(*,*) "al before symmet: ", AL(1:11)
!        write(*,*) "bl before symmet: ", bL(1:11)


!
!        make AL and BL symmetric
call SYMMET(AL,NEN)
call SYMMET(BL,NEN)

!        write(*,*) "al after symmet: ", AL(1:11)
!        write(*,*) "bl after symmet: ", bL(1:11)

!
!        localize temperatures into TL
call LTEMP (IEN(1,n),T,TL,NEN)
!
!        Adjust for differential temperatures
if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
       call LTEMPX (LMTX(1,n),IEN(1,n),TX,TL,NEN)
endif
!
!        calculate product of BL and TL into local load vector PL
call VMPRD (BL,TL,PL,NEN,NEN)
!
!        localize prescribed nodal temperature conditions
call LTEMPBC (TL,T,LMT(1,n),LMTX(1,n),IEN(1,n),NEN)
!
!        adjust the load vector for these conditions
call TEMPBC(TL,PL,AL,NEN)
!
!        add heat production contribution to load vector (teken?)
call LHEAT (XL,QUAD,PL,HEAT,THICK,ierr)

if (ierr.ne.0) then
    write(stderr,1) n
    call stoper()
endif
!
!        add local stiffness and load to global arrays
    call ADDTLD (BT,PL,LMT(1,n),NEN)
    if (NUMSLP+NUMFN.gt.0.and.IDIFT.eq.1) &
       call ADDTLD (BT,PL,LMTX(1,n),NEN)



 !       write(*,*) "calling adtstf with al: ", AL(1:11)


!   write(*,*) '**** itwork; ', itwork, ' and elem: ', nn, ' *****'
if (ITWORK.ne.0) then
     call ADTSTF (AA,AL,ITDIAG,LMT(1,n), LMTX(1,n),NEN)

! call from steady
!     call ADTSTF (AA,sl,ITDIAG,LMT(1,n),LMTX(1,n),NEN)

endif
 enddo
!
 return
 end
