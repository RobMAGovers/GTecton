!********************************************************************
! This file computes the steady state solution.
! Also the initial solution for a time evolution solution
!
! Before the call to STEADY has been made, 
!  The model must have been read
!  and put into the proper structure: [model/mesh/mat]datactx.
!********************************************************************

subroutine STEADY ()!NSD,NDOF,NEN)

USE MESHDATAMODULE
USE MODELDATAMODULE
USE MODELDEFINITION
USE MATERIALSMODULE
USE MODELTOPOLOGY
use iomodule
USE MODELCTX
USE ALGEBRA
use debugmodule, only: lustat, debug, itest, iecho, parallelLog
#ifdef SPARSE
USE AOMODULE
#endif 
!
!-----------------------------------------------------------------------
!    Subroutine to set up and solve steady-state heat equation
!    Secondary routines residing in other libraries are:
!     STORERR - "rdthrm.f"
!     CLEAR, LUSOLVE, LCOORD, SYMMET - "elastc.f"
!     ADDFLX, TSTIFF, LHEAT, THGRAD - thermal topology libs
!-----------------------------------------------------------------------
!
implicit none
!-pass
integer ERROR, i, ierr
!-local
integer :: its, nd

!-init
if (iecho.eq.8) then
    call parallelLog("steady","welcome")
endif

if (MODE.lt.4) then
    return
endif
!
!    Allocate space for FLUX array
ALLOCATE(modeldatactx%FLUX(NSD,meshdatactx%Nelocal),STAT=ERROR)

!
!    Steady-state calculation required?
if (ISTART.lt.0) goto 200 ! restart: goto flux calculation
if (ITMODE.ne.0) goto 100 ! no steady state T calculation
!
if (iecho.eq.1) write(stdout,2)
 2    format(1x,'Steady-state thermal calculation has started')

#ifndef SPARSE
if (IADVEC.ne.0) then
    write(stderr,28)
28       format(1x,'Steady state conduction/advection solution',1x, &
     'only with sparse solvers')
    call exitp(1)
endif
#endif

if (.not. allTemperaturesFixed) then
call CLEAR (modeldatactx%BT,NTEQext,"modeldatactx%BT")
#ifdef SPARSE
call CLEAR (modeldatactx%AA,NTEQ,"modeldatactx%AA (1)")
#else
call CLEAR(modeldatactx%AA,NAA,"modeldatactx%AA (2)")
#endif
endif
!
!    Transfer temperature boundary conditions to temperatures vector

if (iecho.eq.8) then
    call parallelLog("steady","calling tload")
endif


call TLOAD (modeldatactx%ITBC, &
            modeldatactx%TBC, &
            modeldatactx%T, &
            modeldatactx%BT, &
            modeldatactx%IDT, &
            NUMNP)

if (iecho.eq.8) then
    call parallelLog("steady","done tload")
endif

!
#ifdef SPARSE
! this is crucial
call scatterONEdata(meshdatactx,modeldatactx, &
                    getrank(),modeldatactx%T,FILE_outputf,0)
#endif 
!    Construct the global thermal stiffness matrix and add prescribed
!    temperature conditions and internal heating contributions to global
!    load vector.
 
!   if (iecho.eq.1) write(stdout,4)
!    4    format(1x,'Constructing thermal stiffness matrix')

if (iecho.eq.8) then
    call parallelLog("steady","calling mktft")
endif



call MKSTFT (modeldatactx%AA, &
             itdiag,  &
             modeldatactx%BT, &
              meshdatactx%X, &
             modeldatactx%T, &
              meshdatactx%IEN, &
             modeldatactx%LMT, &
             modeldatactx%LMTX, &
              meshdatactx%mat, &
               matdatactx%PRPMAT, &
               matdatactx%PRPTEM, &
             modeldatactx%BOND, &
             modeldatactx%NFAULT, &
             modeldatactx%FAULT, &
             modeldatactx%LMF, &
             NSD, &
             NDOF, &
             NEN)

if (iecho.eq.8) then
    call parallelLog("steady","calling temwnk")
endif


call TEMWNK (modeldatactx%AA, &
             modeldatactx%LMT, &
             itdiag, &
             modeldatactx%ITWINK,  &
             modeldatactx%WINK, &
              meshdatactx%IEN, &
              meshdatactx%mat, &
               matdatactx%PRPMAT,  &
              meshdatactx%X, &
             modeldatactx%D, &
             modeldatactx%LMF, &
             modeldatactx%FAULT, &
             modeldatactx%SKEW, &
             modeldatactx%DXE)
!
!     Add normal boundary flux contributions to load vector

if (NFLX.gt.0) then

    if (iecho.eq.8) then
        call parallelLog("steady","adding fluxes")
    endif

    call ADDFLX (modeldatactx%BT, &
                  meshdatactx%IEN, &
                  meshdatactx%X, &
                 modeldatactx%D, &
                  meshdatactx%mat, &
                   matdatactx%PRPMAT, &
                 modeldatactx%IFLS, &
                 modeldatactx%IFLX, &
                 modeldatactx%BFLX, &
                 modeldatactx%LMT, &
                 modeldatactx%LMTX, &
                 modeldatactx%LMF, &
                 modeldatactx%TFAULT, &
                 modeldatactx%SKEW, &
                 modeldatactx%DXE)

endif
!
!    Invert the stiffness matrix

!    if (iecho.eq.1) write(stdout,8)
!    8    format(1x,'Inverting thermal stiffness matrix')

!
#ifdef SPARSE
call SyncLoadVector(1)
call BuildMat(1) ! Assemble matrix and prepare for solution.

if (debug) then
    call ViewMatrix(1)
endif

call InitSolve(1)

   


!***********************************************************************************
!            solve the theraml equation system
!***********************************************************************************
!      write(*,*) 'before thermal solve: ', NTEQglobal, NTEQlocal

   
do i=1,NTEQlocal
!   write(*,*) 'load: ', i, ' value: ', modeldatactx%BT(i)
enddo

if (iecho.eq.8) then
    call parallelLog("steady","calling solver")
endif

call ITSOLVE (modeldatactx%BT,modeldatactx%AA,its,1) ! AA=solution vector

!      call ITSOLVE (modeldatactx%BT,modeldatactx%AA,its,NTEQglobal,1) ! AA=solution vector
!***********************************************************************************
!***********************************************************************************

if (iecho.gt.0) then
    write(stdout,9) its
 9   format(1X,'Sparse solution of conduction equation required ',I3, &
 ' iterations')
endif

call EQUATE (modeldatactx%BT,modeldatactx%AA,NTEQlocal) ! put solution in load vector

#else

call LUSOLVE (modeldatactx%AA,modeldatactx%BT,itdiag,NTEQ,.TRUE., &
              .FALSE.)
nd = 0
if (IDIFT.eq.1) then
    nd = NUMSLP+NUMFN
endif

call PIVOTS (modeldatactx%AA,itdiag,modeldatactx%IDT, &
             modeldatactx%IDTX, & 
             1,NTEQ,nd,0,NUMNP,MODE)

call CONDNR (modeldatactx%AA,itdiag,NTEQ,ierr)

!    thermal rank check only
if (MODE.eq.4) then
    return
endif
!    backsubstitute the global temperatures array into BT
if (iecho.eq.1) then
    write(stdout,10)
10    format(1x,'Calculating temperatures')
endif

call LUSOLVE (modeldatactx%AA,modeldatactx%BT,itdiag,NTEQ,.FALSE., &
              .TRUE.)
#endif
!
!    localize temperatures in T and TX

if (iecho.eq.8) then
    call parallelLog("steady","calling temps")
endif


call TEMPS (modeldatactx%BT, &
            modeldatactx%T,  &
            modeldatactx%TX, &
            modeldatactx%IDT,&
            modeldatactx%IDTX)

#ifdef SPARSE
!      call scattertemperatures(meshdatactx,modeldatactx,
!     .                         getrank(),modeldatactx%T,FILE_outputf,0)

if (iecho.eq.8) then
    call parallelLog("steady","doing scatterstuff")
endif


call scatterONEdata(meshdatactx,modeldatactx, &
                    getrank(),modeldatactx%T,FILE_outputf,0)

if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
    call scatterONEdata(meshdatactx,modeldatactx, &
                        getrank(),modeldatactx%TX,FILE_outputf,0)
endif

#endif
!
100    continue
!
!    Add thermal anomalies to initial temperatures
if (iecho.eq.8) then
    call parallelLog("steady","calling temadd")
endif

if (NTANOM.gt.0) then
    call TEMADD (meshdatactx%X, &
                modeldatactx%D, &
                modeldatactx%T, &
                modeldatactx%TLST, &
                modeldatactx%ITANOM, &
                modeldatactx%TANOM, &
                0, &
                NUMNP, &
                NSD, &
                NDOF)
endif
!
!    print temperatures

if (iecho.eq.8) then
    call parallelLog("steady","printing output; calling printt")
endif


call PRINTT(modeldatactx%T,NUMNP)
!    print temperatures at discontinuities

if (iecho.eq.8) then
    call parallelLog("steady","calling printx")
endif

! differential temperatures perhaps useful later.
!   call PRINTTX (modeldatactx%TX, &
!            modeldatactx%IDSLP,modeldatactx%NFAULT)

!    compute the fluxes in each element

200   if (iecho.eq.8) then
    call parallelLog("steady","computing flux per element")
endif

if (NFLX.gt.0) then
    call EFLUX (meshdatactx%X, &
                meshdatactx%IEN, &
                meshdatactx%mat, &
               modeldatactx%T, &
               modeldatactx%TX, &
                 matdatactx%PRPTEM, &
               modeldatactx%FLUX, &
               modeldatactx%D, &
               modeldatactx%LMF, &
               modeldatactx%TFAULT, &
               modeldatactx%LMTX, &
               modeldatactx%SKEW, &
               modeldatactx%DXE, &
               NSD, &
               NEN, &
               NDOF)
endif

if (iecho.eq.8) then
    call parallelLog("steady","finished computing fluxed")
endif

!    print fluxes
if (ISTART.ge.0) then
    call PRINTF(modeldatactx%FLUX,NSD,NUMEL)
endif

return
end subroutine
!***********************************************************************

!-------------------------------------------------------------------------------
!    secondary routines are listed below in alphabetic order
!-------------------------------------------------------------------------------
subroutine ADDTLD (BT,PL,LMT,NEN)

USE MODELDEFINITION
USE MODELCTX
!
! Routine to add local load contributions to global vector
!
implicit none
!-pass
integer NEN

integer :: lmt
double precision :: bt, pl
dimension BT(*),PL(NEN),LMT(NEN)

!local
integer :: i, k

do i=1,NEN
    k = IABS(LMT(i))
      if (k.ne.0) then
        BT(k) = BT(k) + PL(i)*ISIGN(1,LMT(i))
    endif
enddo
!
return
end


!***********************************************************************
subroutine ADTSTF (AA,SL,ITDIAG,LMT,LMTX,NEN)
!***********************************************************************
USE MODELDEFINITION
USE MODELCTX
!***********************************************************************
! Routine to add local thermal element stiffness to global array
!***********************************************************************
implicit none
!-pass
 integer NEN


dimension AA(*),SL(NEN,*),ITDIAG(*)
integer :: ITDIAG
double precision :: AA, SL

integer :: LMT(NEN),LMTX(NEN)


#ifdef SPARSE

#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscmatdef.h"
!#include "petsc/finclude/petscvecdef.h"
!#include "petsc/finclude/petsckspdef.h"

   PetscMPIInt irank
#endif
!-locl
 logical diftemp
 integer i,j,k,l,m, rows,cols, isgn
 double precision v
!***********************************************************************
!-init
 diftemp = (NUMSLP+NUMFN.gt.0.and.IDIFT.eq.1)

#ifdef SPARSE
do j=1,NEN
    m = LMT(j)
    if (m.eq.0) then
        cycle
    endif
    do i=1,NEN
        k = LMT(i)
        if (k.eq.0) then
            goto 100
        endif
        v = SL(i,j)
        call value2Tmatrix(k-1,m-1,v,1)
!              call MatSetValuesLocal(T_A,1,k-1,1,m-1,v,ADD_VALUES,ierr)
!              if (ierr.ne.0) then
!                  write(stderr,1) ierr
!    1             format(1x,'ADTSTF: MatSetValuesLocal error ',I4)
!                  call stoper()
!              endif
100     if (.not.diftemp) then
            cycle
        endif
        k = IABS(LMTX(i))
        if (k.eq.0) then
            cycle
        endif
        v = SL(i,j)*ISIGN(1,LMTX(i))
        call value2Tmatrix(k-1,m-1,v,2)
!              call MatSetValuesLocal(T_A,1,k-1,1,m-1,v,ADD_VALUES,ierr)
!              if (ierr.ne.0) then
!                  write(stderr,1) ierr
!                  call stoper()
!              endif
    enddo
enddo

!
if (diftemp) then
    do j=1,NEN
        m = IABS(LMTX(j))
        if (m.eq.0) cycle

        do i=1,NEN
            k = LMT(i)
            if (k.eq.0) cycle
            v = SL(i,j)*ISIGN(1,LMTX(j))
            call value2Tmatrix(k-1,m-1,v,3)
!            call MatSetValuesLocal(T_A,1,k-1,1,m-1,v,ADD_VALUES,ierr)
!            if (ierr.ne.0) then
!                write(stderr,1) ierr
!                call stoper()
!            endif
            k = IABS(LMTX(i))
            if (k.eq.0) cycle
            v = SL(i,j)*ISIGN(1,LMTX(i))*ISIGN(1,LMTX(j))
            call value2Tmatrix(k-1,m-1,v,4)
!            call MatSetValuesLocal(T_A,1,k-1,1,m-1,v,ADD_VALUES,ierr)
!            if (ierr.ne.0) then
!                write(stderr,1) ierr
!                call stoper()
!            endif

        enddo
    enddo
endif

#else

do j=1,NEN
    m = LMT(j)
    if (m.gt.0) then
        L = ITDIAG(m) - m
        do i=1,NEN
            k = LMT(i)
            if (k.ne.0 .and. k.le.m) then
                k = L + k
                AA(k) = AA(k) + SL(i,j)
            endif
            if (diftemp) then
                k = IABS(LMTX(i))
                if (k.ne.0 .and. k.le.m) then
                    k = L + k
                    AA(k) = AA(k) + SL(i,j)*ISIGN(1,LMTX(i))
                endif
            endif
        enddo
    endif
enddo
!
if (diftemp) then
    do j=1,NEN
        m = IABS(LMTX(j))
        if (m.gt.0) then
            L = ITDIAG(m) - m
            isgn = ISIGN(1,LMTX(j))
            do i=1,NEN
                k = LMT(i)
                if (k.gt.0 .and. k.le.m) then
                    k = L + k
                    AA(k) = AA(k) + isgn*SL(i,j)
                endif
                k = IABS(LMTX(i))
                if (k.gt.0 .and. k.le.m) then
                    k = L + k
                    AA(k) = AA(k)+isgn*SL(i,j)*ISIGN(1,LMTX(i))
                endif
            enddo
        endif
    enddo
endif
#endif
!
return
end
!***********************************************************************

#ifdef SPARSE

subroutine value2Tmatrix(k,m,v,callerID)

use spetscmodule, only: T_A
use petscksp,     only: ADD_VALUES
use iomodule,     only: stderr

! This function is an interface between PETSC and and GTecton.
! This inserts values into the stiffness matrix

implicit none


!#include "petsc/finclude/petscsysdef.h"
!#include "petsc/finclude/petsckspdef.h"

! for ADD_VALUES
!#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petscvec.h"
! for MatSetValuesLocal
#include "petsc/finclude/petscmat.h"
! for M_A

! stderr  

integer          :: k,m, ierr, callerID
double precision :: v

!    write(*,*) 'T_A; ',callerID,' adding value ', v, ' to pos: ', k,m
! write(*,*) "adtstf Adding ", v, "to thermal matrix: ", k, m

call MatSetValuesLocal(T_A,1,k,1,m,v,ADD_VALUES,ierr)
!    write(*,*) 'done, with result: ', ierr

if (ierr.ne.0) then
    write(stderr,*) 'ADDSTF: MatSetValuesLocal error ', ierr
    write(stderr,*) 'writing value ',v,' to pos ', k,m, ' of mat T_A'
    call stoper()
endif

end subroutine

#endif


!***********************************************************************
subroutine EFLUX (X,IEN,MAT,T,TX,PRPTEM,FLUX,D,LMF,TFAULT, &
                  LMTX,SKEW,DXE,NSD,NEN,NDOF)
!***********************************************************************
USE MATERIALSMODULE
USE MODELDEFINITION, only: lgdef, numfn, numslp, idift, numel
USE ALGEBRA, only: CLEAR
use constants, only: useDXE
use iomodule, only: stderr
!***********************************************************************
! Routine to calculate heat fluxes in all elements

implicit none
!-pass
integer NSD,NEN,NDOF

integer :: ien, mat, lmf, lmtx
double precision :: x, t, tx, prptem, flux, d, tfault, skew, dxe
   dimension X(NSD,*),IEN(NEN,*),MAT(*),T(*),TX(*),PRPTEM(6,*), &
 FLUX(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
    LMTX(NEN,*),SKEW(2,*),DXE(NDOF,NEN,*)

!***********************************************************************
!***********************************************************************
!-locl
logical          :: faulted
double precision :: XL(24),TL(8),DL(24)
integer          :: m, n, nn, ierr
!***********************************************************************
do nn=1,NUMEL
    n = nn
    call CLEAR(XL,24,"XL")
!         localize (deformed) element coordinatesc
    call LCOORD (X,XL,IEN(1,n),nn)
    if (LGDEF.NE.0) then
        call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,11)
        call ADDSNE (DL,n,NDOF,NEN,useDXE)
        call REZONE (XL,DL,1d0)
    endif
!         localize temperatures
    call LTEMP (IEN(1,n),T,TL,NEN)
!         Modify for discontinuities
    if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
        call LTEMPX (LMTX(1,n),IEN(1,n),TX,TL,NEN)
    endif

    m = MaterialIDOfElement(n)

!         calculate center heat flow

    call LFLUX (XL,TL,IEN(1,n),PRPTEM(1,m),FLUX(1,n),ierr)
    if (ierr.ne.0) then
        write(stderr,1) n
 1         format(1x,'eflux: error in element ',I5)
        call stoper()
    endif
enddo
!
return
end subroutine
!***********************************************************************


!***********************************************************************
subroutine LTEMP (IEN,T,TL,NEN)
!***********************************************************************
! subroutine to localize temperature data
!***********************************************************************
implicit none
!-pass
integer NEN
integer :: ien
double precision :: t, tl
dimension IEN(NEN),T(*),TL(NEN)
!-local
integer :: j,k
!***********************************************************************
do j=1,NEN
    k     = IEN(j)
    TL(j) = T(k)
enddo
!***********************************************************************
return
end subroutine
!***********************************************************************



!***********************************************************************
subroutine LTEMPX (LMTX,IEN,TX,TL,NEN)
!***********************************************************************
! subroutine to adjust for differential temperatures
!***********************************************************************
implicit none
!-pass
integer NEN
integer :: ien, lmtx
double precision :: tx, tl
dimension IEN(NEN),TX(*),TL(NEN),LMTX(NEN)
!-local
integer :: j,k
!***********************************************************************
do j=1,NEN
    if (LMTX(j).ne.0) then
        k = IEN(j)
        TL(j) = TL(j) + ISIGN(1,LMTX(j))*TX(k)
    endif
enddo

return
end subroutine
!***********************************************************************



!***********************************************************************
subroutine LTEMPBC (TL,T,LMT,LMTX,IEN,NEN)
!***********************************************************************
USE MODELDEFINITION, only: numfn, numslp, idift
USE ALGEBRA, only: CLEAR
!***********************************************************************
! Routine to localize prescribed temperature boundary conditions.
! TL(NEN) = element prescribed displacement boundary condition array.
!***********************************************************************
implicit none
!-pass
integer NEN
integer :: lmt, lmtx, ien
double precision :: tl, t
dimension TL(NEN),T(*),LMT(NEN),LMTX(NEN),IEN(NEN)
!-local
integer :: i,k
!***********************************************************************
call CLEAR (TL,NEN,"TL")

do i=1,NEN
    if (LMT(i).gt.0) then
        cycle
    endif
    if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
        if (LMTX(i).gt.0) then
            cycle
        endif
    endif
    k = IEN(i)
    TL(i) = T(k)
enddo

return
end subroutine
!***********************************************************************


!***********************************************************************
subroutine MKSTFT (AA,ITDIAG,BT,X,T,IEN,LMT,LMTX,MAT,PRPMAT, &
                   PRPTEM,BOND,NFAULT,TFAULT,LMF,NSD,NDOF,NEN)
!***********************************************************************
USE MATERIALSMODULE, only: materialidofelement
!   USE MODELDATAMODULE, only: numfn, numslp
!   USE MESHDATAMODULE,  only: numel
USE MODELDEFINITION, only: numel, iadvec, idift, numfn, numslp
!   USE MODELCTX
USE ALGEBRA, only: clear, symmet
use iomodule, only: stderr
!***********************************************************************
! Subroutine to form stiffness matrix and adjust load vector for prescribed
! temperature boundary conditions.

implicit none
!-pass
integer NDOF,NSD,NEN

integer :: itdiag, ien, lmt, lmtx, mat, lmf, nfault
double precision :: aa, bt, t, prpmat, prptem, bond, tfault,x

dimension AA(*),ITDIAG(*),BT(*),X(NSD,*),T(*),IEN(NEN,*), &
 LMT(NEN,*),LMTX(NEN,*),MAT(*),PRPMAT(6,*),PRPTEM(6,*), &
   BOND(NDOF,*),NFAULT(3,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*)
!***********************************************************************

!-locl
double precision :: sl(64),pl(8),tl(8),xl(24),cond(3),vl(24)
logical quad
integer hier, n, m, nn, ierr,i
double precision :: rhocp, heat, thick


do nn=1,NUMEL
    n    = nn
    m = MaterialIDOfElement(n)
    do i=1,NSD
        cond(i) = PRPTEM(3+i,m)
    enddo
    heat = PRPTEM(1,m)
    thick = PRPMAT(6,m)
    quad = (IEN(3,n).ne.IEN(4,n))

!   localize element coordinates
    call LCOORD(X,xl,IEN(1,n),nn)

!   clear local stiffness and load vector
    call CLEAR (sl,NEN*NEN,"sl")
    call CLEAR (pl,NEN,"pl")

!   construct the local stiffness matrix
    call TSTIFF (xl,quad,sl,cond,thick,ierr)
    if (ierr.ne.0) then
        write(stderr,1) n
 1         format(1x,'mkstft: error in element ',I5)
        call stoper()
    endif

!         make it symmetric
    call SYMMET(sl,NEN)

    if (IADVEC.ne.0) then
        rhocp = ABS(PRPMAT(5,m))* PRPTEM(2,m)
!             localize velocities
        call LCOORD (BOND,vl,IEN(1,n),nn)

        call AddFaultDisplacement (vl,n,TFAULT,NDOF,NEN,12)
        call ADVSTF (sl,xl,vl,pl,cond,rhocp,thick,heat,quad, &
                     ierr)
        if (ierr.ne.0) then
            write(stderr,1) n
            call stoper()
        endif
    endif

!   localize prescribed nodal temperature conditions
    call LTEMPBC (tl,T,LMT(1,n),LMTX(1,n),IEN(1,n),NEN)

!   adjust the load vector for these conditions
    call TEMPBC(tl,pl,sl,NEN)

!   add heat production contribution to load vector
    call LHEAT (xl,quad,pl,heat,thick,ierr)
    if (ierr.ne.0) then
        write(stderr,1) n
        call stoper()
    endif

!         add local load to global arrays
    call ADDTLD (BT,PL,LMT(1,n),NEN)
    if (NUMSLP+NUMFN.gt.0.and.IDIFT.eq.1) then
       call ADDTLD (BT,PL,LMTX(1,n),NEN)
    endif

!         add local stiffness to global stiffness matrix
!          write(*,*) '****   element ', n, nn, '********'
    call ADTSTF (AA,sl,ITDIAG,LMT(1,n),LMTX(1,n),NEN)
enddo

return
end subroutine
!***********************************************************************




!***********************************************************************
subroutine printf (FLUX,NSD,NUMEL)
!***********************************************************************
! Routine to print element thermal fluxes
!***********************************************************************

use debugmodule, only: iecho, xit
use iomodule,       only: ascout, outfil, luout, stdout

implicit none

integer, parameter :: MAXP=50
!-pass
integer NSD,NUMEL
double precision :: FLUX(NSD,NUMEL)
!-locl
integer :: map,nn,j,i
character(len=19)fmt
!***********************************************************************
if (iecho.eq.3) then
    write(stdout,1)
 1   format(1x,'for (n=1;n<=NUMEL;n++) {'/ &
 5x,'write n,(FLUX(i,n),i=1,NSD) (I12,X,3PG15.7)'/1x,'}')
endif


if (OUTFIL(2)) then
    if (ascout) then
        fmt = '(I12,1X,3(1PG15.6))'
        ! ascii files get a separate line for each element
        do j=1,NUMEL
            write(luout(2),fmt) j,(FLUX(i,j),i=1,NSD)
        enddo
    else
        ! ... but binary files get all the elements in one chunk
        write(luout(2),err=200) ((FLUX(i,j),i=1,NSD),j=1,NUMEL)
    endif
else
!        write(*,*) 'Not going to write heat stuff'
endif

if (OUTFIL(1)) then
    nn = MAXP + 1
    do j=1,NUMEL
        if (nn.gt.MAXP) then
            nn = 1
            if (NSD.eq.1) then
                write(luout(1),3)
 3               format(///1x,'H E A T   F L U X E S'// &
             1x,'ELEM #',14X,'QX'/)
            else if (NSD.eq.2) then
                write(luout(1),4)
 4               format(///1x,'H E A T   F L U X E S'// &
             1x,'ELEM #',14X,'QX',17X,'QY'/)
            else if (NSD.eq.3) then
                write(luout(1),5)
 5               format(///1x,'H E A T   F L U X E S'// &
             1x,'ELEM #',14X,'QX',17X,'QY',17X,'QZ'/)
            endif
        endif
        write(luout(1),6) j,(FLUX(i,j),i=1,NSD)
 6        format(1x,I6,8X,3(1PG15.6,4X))
        nn = nn + 1
    enddo
endif

return

200    call prerr('PRINTF')
call xit(1," ")
end subroutine
!***********************************************************************


!***********************************************************************
subroutine PRINTT (T,NUMNP)
!***********************************************************************
! Routine to print temperatures
!***********************************************************************
use debugmodule, only: iecho, parallelLog, xit
use iomodule,    only: ascout, luout, outfil, stdout

implicit none

integer, parameter :: MAXP=50

!-pass
integer NUMNP
double precision :: T(NUMNP)

! for luout(2), pointing to the fedsk file, and ascout

!-locl
character(len=16)fmt
integer :: iPoint
integer :: n,nn
!***********************************************************************
if (iecho.eq.3) then
    write(stdout,1)
 1   format(1x,'for (n=1;n<=NUMNP;n++) {'/ &
 5x,'write n,T(n) (I5,X,1PG15.6)'/1x,'}')
endif

if (iecho.eq.8) then
    call parallelLog("printt","welcome")
endif


if (OUTFIL(2)) then
    if (ascout) then
        fmt = '(I12,1X,1PG15.6)'
        ! text files get the index of the nodal point 'n'...
        do iPoint=1,numnp
            write(luout(2),fmt) iPoint,T(iPoint)
        enddo
    else
        ! ... but binary files do not. It saves space.
        write(luout(2),err=200) (T(N),N=1,NUMNP)
    endif
endif

if (OUTFIL(1)) then
    nn = MAXP + 1
    do n=1,NUMNP
        if (nn.gt.MAXP) then
            nn = 1
            write(luout(1),3)
 3               format(///1x,'N O D A L   T E M P E R A T U R E S'// &
             5X,'NODE #',18X,'TEMPERATURE'/)
        endif
        write(luout(1),4) n,T(n)
 4            format(5X,I6,17X,1PG15.6)
        nn = nn + 1
    enddo
endif

return

200    call prerr('PRINTT')
call xit(1," ")
end subroutine
!***********************************************************************



!***********************************************************************
subroutine PRINTTX (TX,IDSLP,NFAULT)
! currently not in use
! But do not remove! May be useful later.

!***********************************************************************
USE MODELDEFINITION, only: numnp, numfn, numsn, idift
use iomodule,        only: ascout, outfil, luout
!***********************************************************************
!
! Routine to print differential temperatures
!
!***********************************************************************
implicit none
!-pass
integer :: idslp, nfault
double precision :: tx
dimension TX(*),IDSLP(*),NFAULT(3,*)
!***********************************************************************
!***********************************************************************
!-locl
character(len=15)fmt
integer :: j, m, n, nout
!
if ((.not.OUTFIL(1)).and.(.not.OUTFIL(2))) then
    return
endif

if (NUMFN+NUMSN.le.0.or.IDIFT.eq.0) then
    return
endif
!
if (OUTFIL(2)) then
!        write (*,*) 'Print differential temperatures'
    if (ascout) then
        fmt = '(I8,1X,1PG15.6)'
        ! ascii files get the index of the node
        write(luout(2),fmt) (n,TX(n),n=1,NUMNP)
    else
        ! binary files do not get the index
        write(luout(2)) (TX(n),n=1,NUMNP)
    endif
endif
!
if (OUTFIL(1)) then
    nout = 0
    if (NUMSN.gt.0) then
        do n=1,NUMSN
            if (MOD(nout,50).eq.0) then
                write(luout(1),2)
 2               format(///1X,'D I F F E R E N T I A L',3X, &
             'T E M P E R A T U R E S'///7X, &
             'NODE',5X,'TEMPERATURE'/)
            endif
            m = IDSLP(n)
            write(luout(1),3) m,TX(m)
 3          format(3X,I8,5X,F11.5)
            nout = nout + 1
        enddo
    endif

    if (NUMFN.gt.0) then
        do n=1,NUMNP
            do j=1,NUMFN
                if (NFAULT(2,j).eq.n) then
                    goto 100
                endif
            enddo
            cycle
100         if (MOD(nout,50).eq.0) then
                write(luout(1),2)
            endif
            write(luout(1),3) n,TX(n)
            nout = nout + 1
        enddo
    endif
endif
!
return
end subroutine
!***********************************************************************


!***********************************************************************
subroutine TEMADD (X,D,T,TLST,ITANOM,TANOM,ITIME,NUMNP,NSD,NDOF)
!***********************************************************************
! Subroutine to add thermal anomaly data to global temperatures
!***********************************************************************
USE MODELDEFINITION, only: NFLX,NTANOM,ITPRES,ITECT,NTWINK
use debugmodule, only: iecho, xit
use iomodule,     only: stdout, stderr

implicit none
!-pass
!    integer NTANOM,NUMNP,ITIME,NSD,NDOF
integer NUMNP,ITIME,NSD,NDOF

integer :: itanom
double precision :: t, tlst, tanom, x, d
dimension T(*),TLST(*),ITANOM(3,*),TANOM(*),X(NSD,*),D(NDOF,*)
!-local
logical :: done
integer :: n, np
!***********************************************************************
!***********************************************************************
!-external
double precision, external :: AddAn ! add anomaly
!***********************************************************************
DONE = .false.
if (NTANOM.eq.0) then
    return
endif

if (NTANOM.gt.0) then

    do n=1,NTANOM

        if (ITANOM(3,n).lt.0) then
!               Start applying at ITIME
            if (ITIME.lt.-ITANOM(3,n)) then
                cycle
            endif
        else
!               Apply once at ITIME
            if (ITANOM(3,n).ne.ITIME) then
                cycle
            endif
        endif

        DONE = .true.
        np = ITANOM(1,n)

        if (np.ge.1 .and. np.le.NUMNP) then
            if (ITANOM(2,n).eq.0) then
!                   Replace
                T(np) = TANOM(n)
!                    No tectonic contribution (ITECT<>0) if T(i)=TLST(i)
                if (ITECT.ne.0) TLST(np) = TANOM(n)
            else
!                   Add
                T(np) = T(np) + TANOM(n)
                if (ITECT.ne.0) TLST(np) = TLST(np) + TANOM(n)
            endif
        else
            write(stderr,1) np
 1               format(1x,'temadd: nodal point out of range: np=', &
             I6)
        call xit(1," ")
        endif
    enddo

else

    if (IABS(NTANOM).ne.ITIME) then
        return
    endif

    DONE = .true.
    do n=1,NUMNP
!        ADD analytical anomaly
        T(n) = T(n) + AddAn (X(1,n),D(1,n),NDOF,NSD)
        if (ITECT.ne.0) then
            TLST(n) = T(n)
        endif
    enddo
endif

if (DONE .and. iecho.eq.1) then
    write(stdout,2)
 2        format(1x,'Thermal anomaly added')
endif
return
end subroutine
!***********************************************************************



!***********************************************************************
subroutine TEMPBC (TL,P,S,NEN)
!***********************************************************************
! routine to adjust the load vector for essential boundary conditions
!***********************************************************************
implicit none
!-pass
integer NEN
double precision :: TL(NEN),P(NEN),S(NEN,NEN)
!-local
integer          :: i,j
!***********************************************************************
do j=1,NEN
    do i=1,NEN
        P(i) = P(i) - S(i,j)*TL(j)
    enddo
enddo
return
end
!***********************************************************************


!***********************************************************************
subroutine TEMPS (BT,T,TX,IDT,IDTX)
!***********************************************************************
USE MODELDEFINITION, only: numnp, numfn, numslp, idift
!***********************************************************************
!     subroutine to transfer computed temperatures into T or TX
!***********************************************************************
implicit none
!-pass

double precision :: BT,T,TX
integer          :: IDT,IDTX
dimension BT(*),T(*),TX(*),IDT(*),IDTX(*)
!-local
integer :: i,k

!***********************************************************************
do i=1,NUMNP
    k = IDT(i)
    if (k.gt.0) then
        T(i) = BT(k)
    endif

    if (NUMFN+NUMSLP.gt.0.and.IDIFT.eq.1) then
        k = IDTX(i)
        if (k.gt.0) then
            TX(i) = BT(k)
        endif
    endif
enddo

return
end
!***********************************************************************

!***********************************************************************
subroutine TLOAD (ITBC,TBC,T,BT,IDT,NUMNP)
!***********************************************************************
! subroutine to transfer nodal boundary conditions into appropriate vectors
!***********************************************************************
implicit none
!-pass
integer NUMNP
integer          :: itbc, idt,k,n
double precision :: tbc, t, bt
dimension ITBC(1,NUMNP),TBC(NUMNP),T(NUMNP),BT(*),IDT(*)
!***********************************************************************

do n=1,NUMNP
    if (ITBC(1,n).eq.1) then
        T(n) = TBC(n)
    else if (ITBC(1,n).eq.3) then
        k = IDT(n)
        if (k.gt.0) then
            write(*,*) 'tload adds ', TBC(n), ' to load vec entry: ', k
            BT(k) = BT(k) + TBC(n)
        endif
    endif
enddo

return
end
!***********************************************************************
