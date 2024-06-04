subroutine INITAL ()

USE MESHDATAMODULE,  only: meshdatactx
USE MODELDATAMODULE, only: modeldatactx, &
                           HasSlipperyNodes, &
                           idiag, meq, showfbc
USE MODELDEFINITION
USE MODELTOPOLOGY
USE MATERIALSMODULE, only: matdatactx, grav
USE TIMESTEPMODULE
USE MODELCTX
USE ALGEBRA
USE CONSTANTS,       only: MB

#ifdef PLOTNONZEROES
use debugmodule,     only: debug, &
                           parallellog, &
                           startClock, &
                           logTimeStamp, &
                           iecho, &
                           debug, &
                           allocateError, &
                           nonZeroMatrixPlot, &
                           CheckArrayForNans, &
                           xit
#else
use debugmodule,     only: debug, &
                           parallellog, &
                           startClock, &
                           logTimeStamp, &
                           iecho, &
                           debug, &
                           allocateError, &
                           CheckArrayForNans, &
                           delay, &
                           xit
#endif

#ifdef SPARSE
USE AOMODULE
#endif
use iomodule          

!
! Subroutine to construct and find solution to the elastic mechanical problem
!

implicit none

!-locl
!        logical output,PrtMat
logical PrtMat
integer ERROR
double precision, external :: FLOAT
double precision :: accuracy, psumsq

!-init
DELTP = 0d0

#ifdef DEBUG
!write(*,*) "rank", getrank(), "at start inital, skew", modeldatactx%skew
#endif

if (iecho.eq.8) then
    call ParallelLog("elastc", "Computing elastic solution")
endif


!#ifdef SPARSE
!    allocate(modeldatactx%B(nequations))
!#else
!    allocate(modeldatactx%B(NEQlocal))
!#endif




!       start timer
call followtime (0)
!
ALLOCATE(modeldatactx%DMAT(NSTR,NSTR,NUMEL),STAT=ERROR)
call allocateError("modeldat DMAT", ERROR)


if (NINTG.gt.0) then
    ALLOCATE(modeldatactx%BETA(NSTR,NUMEL), STAT=ERROR)
    call allocateError("modeldatactx%BETA", ERROR)

    ALLOCATE(modeldatactx%BETB(NSTR,NUMEL), STAT=ERROR)
    call allocateError("modeldatactx%BETB", ERROR)

    if (NSLSKEW.gt.0) then
        ALLOCATE(modeldatactx%OSKEW(2,NUMNP), STAT=ERROR)
        call allocateError("modeldatactx%OSKEW", ERROR)
    endif
endif

if (MODE.ne.1 .and. MODE.ne.2 .and. MODE.ne.6) then
    if (iecho.eq.2 .or. iecho.eq.8) then
        write(*,*) "mode not 1, 2 or 6"
        write(*,*) "no mechanical solution required"
        write(*,*) "Returning."
    endif
    return
endif

!
#ifndef SPARSE
call CLEAR(modeldatactx%A,NSIZEA,"modeldatactx%A")
#endif

if (ISTART.ge.0) then
!   No restart: clear arrays
    call CLEAR (modeldatactx%D,NDOF*meshdatactx%nvlocal,"modeldatactx%D")
    call CLEAR (modeldatactx%STN,NSTR*NUMEL,"modeldatactx%STN")
    call CLEAR (modeldatactx%STR,NSTR*NUMEL,"modeldatactx%STR")

    if (IDIG(ICVIS,4,4).ne.0) then
        call CLEAR (modeldatactx%VPSTR,NSTR*NUMEL,"modeldatactx%VPSTR")
    endif

    if (NUMSLP.gt.0) then
#ifdef EXPERIMENTAL_ReducedSlipAlloc
        call CLEAR (modeldatactx%DXE,NDOF*NEN*nElemsWithSlipperyNodes,"modeldatactx%DXE")
#else
        call CLEAR (modeldatactx%DXE,NDOF*NEN*NUMEL,"modeldatactx%DXE")
#endif
    endif

endif


#ifdef DEBUG
!write(*,*) "a update size1 ", size(modeldatactx%D,1), size(modeldatactx%D,2)
write(*,*) "rank", getrank(),"has NSLSKEW", NSLSKEW
#endif

if (NSLSKEW.gt.0) then
    ! NSLSKEW is the global number of fault parallel element entries

#ifdef DEBUG
!    write(*,*) "rank", getrank(),"call slipperyskew with IDSLP",modeldatactx%IDSLP
    write(*,*) "rank", getrank(),"call slipperyskew with SKEW: ", modeldatactx%SKEW
    write(*,*) "rank", getrank(),"call slipperyskew with IEN: ", meshdatactx%IEN
#endif

    ! modeldatactx%IDSLP contains the local indices
    ! of the angles nodes on this partition.
    call SlipperySkew(modeldatactx%IDSLP,modeldatactx%SKEW, &
       modeldatactx%NSELSD,meshdatactx%IEN, &
        meshdatactx%X,modeldatactx%D,modeldatactx%TFAULT, &
        modeldatactx%DXE,modeldatactx%LMF)


#ifdef DEBUG
    write(*,*) "rank", getrank(), "done formk SKEW ", modeldatactx%SKEW
#endif

    if (iecho.eq.8) then
        call ParallelLog("elastc", "Set angles when required")
    endif
endif

#ifdef DEBUG
!write(*,*) "rank", getrank(), "after slipperyskew, skew", modeldatactx%skew
#endif


if (ISTART.lt.0) then
!    write(*,*) 'istart'
    goto 100
endif


PrtMat = output (outputcontroldatactx%IMATPR, &
                 NPRMAT,IPRMAT,NSTEP)

if (NELAST.ne.0) then
 
!           update STN array with pre-stresses
    call PRTRF (modeldatactx%STN0, modeldatactx%ISELM, &
        modeldatactx%ISTIME, modeldatactx%STN)

!           compute equivalent load vector B from STN


#ifdef SPARSE
    call CLEAR(modeldatactx%B,nequations,"modeldatactx%B")
#else
    call CLEAR (modeldatactx%B,NEQlocal,"modeldatactx%B")
#endif
! SPARSE

    call RESIDU (PSUMSQ,ACCURaCY)

!   set BTOT equal to B, or read it from "BTOT.DAT" if it exists
    call SETBTOT (modeldatactx%BTOT,modeldatactx%B,NEQlocal)

endif

#ifdef SPARSE
call CLEAR(modeldatactx%B,nequations,"modeldatactx%B")
#else
call CLEAR (modeldatactx%B,NEQlocal,"modeldatactx%B")
#endif
! SPARSE
!


if (iecho.eq.8) then
    call ParallelLog("elastc", "calling load")
endif

!       Transfer boundary conditions into global load displacement vector
!     from various boundary condition sources, they are combined
! into   modeldatactx%D for the displacement and velocity BC, 
! and in modeldatactx%B for the forces
call LOAD (modeldatactx%ID,modeldatactx%IBOND,modeldatactx%BOND, &
           modeldatactx%D, modeldatactx%DELD, modeldatactx%B, &
           meshdatactx%X,NUMNP,NEQlocal,0)

if (iecho.eq.8) then
    call ParallelLog("elastc", "Added boundary condition to load vector")
endif

! Set faulted node displacements
call LOADF (meshdatactx%IEN,modeldatactx%NFAULT, modeldatactx%FAULT, &
            modeldatactx%DFAULT,NUMEL)

if (iecho.eq.8) then
    call ParallelLog("elastc", "Added fault contribution to load vector")
endif

!       Add differential forces across internal free interfaces

if (NUMSLP.gt.0) then
    if (iecho.eq.1) then
         write(stdout,25)
    endif
25       format(1x,'Slippery node forces are being added')
    call LOADX (modeldatactx%B, modeldatactx%IDX, &
                modeldatactx%DIFORC, &
                NUMSLP,NUMNP)

#ifdef EXPERIMENTAL_FRICTION

! in case friction is employed in concert with a differential force
! across the internal interface, the DIFORC array needs to be cleared
! after loading into the load vector, so that the array can be used later
! on in the time stepping algorithms

call CLEAR (modeldatactx%DIFORC,NDOF*NUMNP)

#endif

endif

if (iecho.eq.8) then
    call ParallelLog("elastc", "Added slippery contribution to load vector")
endif

!       Add forces from initial stresses to global load vector
if (NELAST.eq.0 .and. &
    NPRE  .gt.0) then

    call PRESTN (modeldatactx%B, modeldatactx%STN0, &
     modeldatactx%ISELM, modeldatactx%ISTIME,meshdatactx%X, &
     modeldatactx%D, modeldatactx%SKEW,meshdatactx%IEN, &
     modeldatactx%LM, modeldatactx%LMX,meshdatactx%mat, &
     matdatactx%PRPMAT,modeldatactx%LMF, modeldatactx%TFAULT, &
     modeldatactx%DXE)

    if (iecho.eq.8) then
        call ParallelLog("elastc", "Added prestress to load vector")
    endif

endif

!       Add pressure forces, if present, to global load vector
if (NUMPR.ne.0) then
    if (iecho.eq.1) then
        write(stdout,35)
35       format(1x,'Pressure loads are being added')
    endif

    call ADDPR(modeldatactx%B, &
                meshdatactx%X, &
               modeldatactx%D, &
                meshdatactx%IEN, &
               modeldatactx%LM, &
                meshdatactx%mat, &
                 matdatactx%PRPMAT, &
               modeldatactx%IELNO, &
               modeldatactx%ISIDE, &
               modeldatactx%PRES, &
               modeldatactx%LMF,  &
               modeldatactx%TFAULT, &
               modeldatactx%SKEW,  &
               modeldatactx%DXE)
    if (iecho.eq.8) then
        call ParallelLog("elastc", "Added pressure to load vector")
    endif
endif


! Add stress forces, if present, to global load vector
if (NUMSTR.ne.0) then
!	write(*,*) "calling addstr", size(modeldatactx%STRS,1), size(modeldatactx%STRS,2)
    call ADDSTR (modeldatactx%B, &
                  meshdatactx%X, &
                 modeldatactx%D, &
                  meshdatactx%IEN, &
                 modeldatactx%LM, &
                  meshdatactx%mat, &
                   matdatactx%PRPMAT, &
                 modeldatactx%IELSTR, &
                 modeldatactx%ISSIDE, &
                 modeldatactx%ISTR, &
                 modeldatactx%LMF, &
                 modeldatactx%TFAULT, &
                 modeldatactx%SKEW, &
                 modeldatactx%DXE)

    if (iecho.eq.8) then
        call ParallelLog("elastc", "Added stress to load vector")
    endif

endif



if (IGRAV.ne.0) then
    if (iecho.eq.1) write(stdout,45)
45       format(1x,'Gravity body forces are being added')

    call GLOAD (modeldatactx%B,meshdatactx%X,modeldatactx%SKEW, &
      GRAV,meshdatactx%IEN,modeldatactx%LM,modeldatactx%LMX, &
      meshdatactx%mat,matdatactx%PRPMAT,modeldatactx%T, &
      modeldatactx%TX,modeldatactx%LMTX,matdatactx%PRPTEM, &
      NSD,NDOF,NEN, igrav)

endif

if (iecho.eq.8) then
    call ParallelLog("elastc", "Added gravity to load vector")
endif


!
!       Save global load for later use (excluding pseudoforces due to
!       displacement bcs
!        call UPDATE (A(IPBTOT),A(IPB),NEQ)

! take care of stride
#ifdef SPARSE

! TODO: document what is going on here.
!write(*,*) "calling syncloadvec 1"
!********** HERE WE NEED A CLEAR COMMENT THAT EXPLAINS THE FUNCTIONALITY OF SYNCLOADVECTOR
call SyncLoadVector(0)
if (iecho.eq.8) then
    call ParallelLog("elastc", "Synced load vector")
endif

#endif

!call CheckArrayForNans(modeldatactx%B,"b before updating Btot")
!call CheckArrayForNans(modeldatactx%BTOT,"btot before update with B")


call UPDATE(modeldatactx%BTOT,modeldatactx%B,NEQ)

!#ifdef Lukas
!write(*,*) "BB rank", getrank()," B after update", modeldatactx%B(1:20)
!#endif



!call CheckArrayForNans(modeldatactx%BTOT,"btot after update with B")


if (iecho.eq.8) then
    call ParallelLog("elastc", "Updates total load vector")
endif


!write(*,*) "post update size1 ", size(modeldatactx%D,1), size(modeldatactx%D,2), modeldatactx%D(1:2, 1:5)

!
!       Construct global stiffness matrix and add displacement bcs
if (iecho.eq.1) then
    write(stdout,50)
50       format(1x,'Stiffness matrix under construction')
endif
    
if (iecho.eq.6) then
    call startClock()
endif


if (iecho.eq.8) then
    call ParallelLog("elastc", "Starting formk")
endif


#ifdef DEBUG
!write(*,*) "rank", getrank(), "pre formk SKEW ", modeldatactx%SKEW
#endif

#ifdef Lukas
!write(*,*) "CC rank", getrank()," B", modeldatactx%B
#endif





call FORMK( &
 modeldatactx%A,IDIAG,modeldatactx%B,meshdatactx%X,modeldatactx%D, &
 modeldatactx%SKEW,meshdatactx%IEN,modeldatactx%lm, &
 modeldatactx%lmx,meshdatactx%mat, &
 modeldatactx%DMAT,modeldatactx%S, &
 modeldatactx%forShowfbcIAeqnr, &
 modeldatactx%forShowfbcIAXeqnr, & ! A(IPIA),A(IPIAX), &
 NSD,NDOF,NSTR,NEN)



if (iecho.eq.8) then
    call ParallelLog("elastc", "Completed formk")
endif


if (iecho.eq.6) then
    call logTimeStamp("Combobulated stiffness matrix (formk)")
endif


!       Modify stiffness matrix for winkler restoring forces, if present
if (NWINK.gt.0) then
   call WINKLR(modeldatactx%ID, &
     modeldatactx%IWINK,modeldatactx%WINK,NDOF,NUMNP,0)
endif


! diff winkler
if (NWINKX.gt.0) then
   call WINKLR(modeldatactx%IDX, &
     modeldatactx%IWINX,modeldatactx%WINX,NDOF,NUMNP,1)
endif
!
!       Modify stiffness matrix for winkler restoring pressures, if present
#ifdef SPARSE
if (NUMWNK.ne.0) then
    call ADDWNK (modeldatactx%IWELM, &
                 modeldatactx%IWSIDE, &
                 modeldatactx%IWTIME, &
                 modeldatactx%WPRES, &
                 meshdatactx%IEN, &
                 meshdatactx%mat, &
                 meshdatactx%X, &
                 modeldatactx%D, &
                 matdatactx%PRPMAT, &
                 modeldatactx%A, &
                 modeldatactx%LM, &
                 IDIAG, &
                 modeldatactx%LMF, &
                 modeldatactx%TFAULT, &
                 modeldatactx%SKEW, &
                 modeldatactx%DXE, &
                 modeldatactx%LMX, &
                 modeldatactx%NSLIP)

    if (iecho.eq.8) then
        call ParallelLog("elastc", "Added winkler pressures to stiffness matrix")
    endif

endif
#endif

!       Modify stiffness matrix for initial displacement bcs

! do stuff in the A vector
call INITD (IDIAG,modeldatactx%A,NDOF,NUMNP,NEQlocal,0)

! modify the load vector
call INITD (IDIAG,modeldatactx%A,NDOF,NUMNP,NEQlocal,1)

!       Add faulted node pseudoforces to global load vector
if (NUMFN.gt.0) then

    if (iecho.eq.8) then
        call ParallelLog("elastc", "Adding faulted entries to load vector")
    endif


    if (iecho.eq.1) write(stdout,55)
55       format(1x,'Faulted node loads are being added')


    call FORMF (modeldatactx%B, &
                 meshdatactx%X, &
                modeldatactx%SKEW, &
                 meshdatactx%IEN, &
                modeldatactx%LM, &
                modeldatactx%D, &
                modeldatactx%LMX, &
                 meshdatactx%mat, &
                  matdatactx%PRPMAT, &
                modeldatactx%DMAT, &
                modeldatactx%LMF, &
                modeldatactx%NFAULT, &
                modeldatactx%DFAULT, &
                modeldatactx%TFAULT, &
                modeldatactx%DXE, &
                NSD,NDOF,NSTR,NEN)

    if (iecho.eq.8) then
        call ParallelLog("elastc", "Done FormF")
    endif

endif

!
#ifdef SPARSE
call BuildMat(0) ! Assemble matrix and prepare for solution.
if (iecho.eq.8) then
    call ParallelLog("elastc", "Finished matrix assembly")
endif
#endif
! SPARSE
!
if (debug) then
    if (iecho.eq.8) then
        call ParallelLog("elastc", "Printing debug output")
    endif
!   Print matrix
#ifdef SPARSE
    call ViewMatrix(0)
#else
    call PRINTA (modeldatactx%A,IDIAG,NEQ)
#endif
!   Print load vector
    call PRINTB (modeldatactx%B, modeldatactx%ID, modeldatactx%IDX, NDOF)
endif


if (iecho.eq.8) then
    call ParallelLog("elastc", "printed matrix")
endif


!       Print matrix diagonal if required
if (PrtMat) then

#ifdef SPARSE
    if (.not.debug) then
        call ViewMatrix(0)
    endif
#else
    call PRINTP (modeldatactx%A,IDIAG,modeldatactx%ID,modeldatactx%IDX,NDOF,0)
#endif

endif



#ifdef PLOTNONZEROES
call nonZeroMatrixPlot() ! from debug module
#endif




#ifdef checkLoadvecForNans
!call CheckArrayForNans(modeldatactx%B, "load pre elastic solve")
#endif 



#ifdef Lukas
write(*,*) "EE rank", getrank(),"B", modeldatactx%B
#endif



#ifndef SPARSE
!       this version of gtecton is SPARSE, omit this
!       Factorize the global stiffness matrix

! we practically never do this...

if (iecho.eq.8) then
    call ParallelLog("elastc", "Starting a dense solve")
endif
call SOLVE (modeldatactx%A,&
            modeldatactx%B, &
            modeldatactx%BRES, &
            IDIAG,modeldatactx%ID, &
            modeldatactx%IDX, NUMNP,NEQ,NUMSN,NDOF,.TRUE.,.FALSE.)
#endif

!       if mechanical rank check only, stop here
if (MODE.eq.1) then
    return
endif

if (iecho.eq.8) then
    call ParallelLog("elastc", "Starting a sparse solve")
endif

call SOLVE (modeldatactx%A, &
            modeldatactx%B, &
            modeldatactx%BRES, &
            IDIAG, &
            modeldatactx%ID, &
            modeldatactx%IDX, &
            NUMNP, &
            NEQlocal, &
            NUMSN, &
            NDOF, &
            .FALSE., &
            .TRUE.)



if (iecho.eq.8) then
     call parallelLog("inital", "finished solve")
endif

!       Localize the displacements in D(NDOF,NUMNP) and DX
call DISP (modeldatactx%B, &
           modeldatactx%D, &
           modeldatactx%ID, &
           modeldatactx%LINK, &
           NDOF, &
           NUMNP, &
           NLINK)

if (iecho.eq.8) then
     call parallelLog("inital", "localised displacement")
endif


if (NUMSLP.ne.0) then
    call DISP (modeldatactx%B, &
               modeldatactx%DX, &
               modeldatactx%IDX, &
               modeldatactx%LINK, &
               NDOF, &
               NUMNP, &
               0)
endif


#ifdef SPARSE
!      call scatterdisplacements(meshdatactx,modeldatactx,
!     .                         getrank(),modeldatactx%D,FILE_outputf,0)
   call scatterNDOFdata(meshdatactx,modeldatactx, &
             getrank(),modeldatactx%D,FILE_outputf,0,5)


if (iecho.eq.8) then
     call parallelLog("inital", "scattered NDOF data")
endif



   if (NUMSLPglobal.ne.0) then
  call scatterNDOFdata(meshdatactx,modeldatactx, &
             getrank(),modeldatactx%DX,FILE_outputf,0,6)
   endif
#endif
! SPARSE



if (NWINK.gt.0) then 
    if (iecho.eq.8) then
         call parallelLog("inital", "Doing WinklF")
    endif
    call WINKLF(modeldatactx%BTOT, &
                modeldatactx%D, &
                modeldatactx%ID, &
                modeldatactx%IWINK, &
                modeldatactx%WINK, &
                NDOF, &
                 meshdatactx%Nvlocal, &
                0)  ! 0 to indicate regular Winkler forces
    if (iecho.eq.8) then
         call parallelLog("inital", "Done WinkLF")
    endif
endif


if (NWINKX.gt.0) then 
    if (iecho.eq.8) then
         call parallelLog("inital", "Doing differential Winkler forces")
    endif

    call WINKLF(modeldatactx%BTOT, &
                modeldatactx%DX, &
                modeldatactx%IDX, &
                modeldatactx%IWINX, &
                modeldatactx%WINX, &
                NDOF, &
                 meshdatactx%Nvlocal, &
                1) ! 1 to indicate differential Winkler forces

    if (iecho.eq.8) then
         call parallelLog("inital", "Done differential Winkler forces")
    endif
endif

!       Add forces due to Winkler pressures to BTOT
if (NUMWNK.ne.0) then
    if (iecho.eq.8) then
         call parallelLog("inital", "Doing Winkler pressures")
    endif
    call WINKLP (modeldatactx%BTOT,modeldatactx%D, &
     modeldatactx%IWELM,modeldatactx%IWSIDE,modeldatactx%IWTIME, &
     modeldatactx%WPRES,meshdatactx%IEN,meshdatactx%mat, &
     meshdatactx%X,modeldatactx%D,matdatactx%PRPMAT, &
     modeldatactx%LM,modeldatactx%LMF,modeldatactx%TFAULT, &
     modeldatactx%SKEW,modeldatactx%DXE,modeldatactx%LMX, &
     modeldatactx%NSLIP,modeldatactx%DX)
    if (iecho.eq.8) then
         call parallelLog("inital", "Done Winkler pressures")
    endif
endif
!

!       Set Maxwell boundary conditions load vector
if (NMAXW.gt.0) then
    call MXWLD (modeldatactx%IWINK,modeldatactx%WINK, &
     modeldatactx%FTOT,modeldatactx%FMAXW,modeldatactx%D,NDOF)
endif
!
if (showfbc) then ! Show Force Bounday Conditions
    call GETNEQBC (modeldatactx%IBOND, &
                   modeldatactx%S, &
                   modeldatactx%forShowfbcIAeqnr, &
                   modeldatactx%D, &
                   modeldatactx%BOND, &
                   NDOF,NUMNP,MEQ)

    call PRTNEQBC (modeldatactx%IBOND, &
                   modeldatactx%forShowfbcIAeqnr,&
                   modeldatactx%BOND,&
                   NDOF,NUMNP,NSTEP)
endif

!
!       Rotate skew coordinates to global system
!if (NUMROT.ne.0) then
! note that not every partition may have 
! nodal points with rotated coordinates. 
! So some threads may take this loop, while others may not.

    call RDISP (modeldatactx%D,modeldatactx%SKEW,NDOF,NUMNP)

    if (HasSlipperyNodes()) then
        call RDISP(modeldatactx%DX,modeldatactx%SKEW,NDOF,NUMNP)
    endif

    call ROTFLT(modeldatactx%NFAULT,modeldatactx%DFAULT, &
       modeldatactx%SKEW,NDOF,NUMFN)


#ifdef SPARSE

    call scatterNDOFdata(meshdatactx,modeldatactx, &
             getrank(),modeldatactx%D,FILE_outputf,0,7)

    if (iecho.eq.8) then
         call parallelLog("inital", "scattered more")
    endif


    if (NUMSLPglobal.ne.0) then
        call scatterNDOFdata(meshdatactx,modeldatactx, &
             getrank(),modeldatactx%DX,FILE_outputf,0,8)
    endif
#endif

    if (iecho.eq.8) then
        call parallelLog("inital", "Finished Euler angle stuff")
    endif




! MPI barrier so that the threads with Euler angles have opportunity to catch up
! with the ones that did not.

if (NUMSLP.gt.0) then
    call UPDXE (modeldatactx%DX,meshdatactx%IEN,modeldatactx%LMX, &
     modeldatactx%SKEW,modeldatactx%IDSLE,modeldatactx%DXE, &
     NDOF,NEN)
endif

if (NUMFN.gt.0) then
    ! Copy DFAULT into TFAULT   (maybe: add dfault to tfault, 
    !                           for faulted displacement in multiple times)
    call EQUATE (modeldatactx%TFAULT,&
                 modeldatactx%DFAULT, &
                 NDOF*NUMFN)
endif

!       Store surface nodal point coordinates
100     call ADDSRF (modeldatactx%ISURF, modeldatactx%SURF, &
             meshdatactx%X, modeldatactx%D,NSD,NDOF,NSURF,0)

if (NUMPR.lt.0 .and.  &
    NSURF.gt.0 .and. &
    LGDEF.ne.0) then
    call EQUATE (modeldatactx%SURF0, &
                 modeldatactx%SURF,&
                 NSD*NSURF)
endif

if (ISTART.lt.0) then
    goto 300
endif

if (iecho.eq.8) then
     call parallelLog("inital", "writing displacements to fedsk files")
endif

! write the displacement to file.
call PRINTD (modeldatactx%D,NUMNP,0)

! write differential forces when needed
#ifdef EXPERIMENTAL_FRICTION
    if (iecho.eq.8) then
         call parallelLog("inital", "writing diff forces to fedsk files")
    endif

    if (NUMSLPglobal.gt.0) then
        call PRINTD(modeldatactx%DIFORC,NUMNP,2)
    else
        ! modeldatactx%DIFORC niet gealloceerd
    endif
#endif

if (iecho.eq.8) then
     call parallelLog("inital", "writing split nodes (if they are there) to fedsk files")
endif

! print faulted nodes
call PSPLIT(modeldatactx%TFAULT,modeldatactx%NFAULT, &
            NDOF,NUMFN)

if (NUMSLPglobal.ne.0) then
    call PRINTDXE(modeldatactx%DXE,meshdatactx%IEN, &
                  modeldatactx%IDSLE,0)
endif


! Compute stresses and strains in each element


if (iecho.eq.8) then
     call parallelLog("inital", "computing stress")
endif


call STRESS (meshdatactx%X, &
             meshdatactx%IEN, &
            modeldatactx%D, &
            modeldatactx%STN, &  ! stress goes here &
            modeldatactx%LMF, &
            modeldatactx%TFAULT, &
            modeldatactx%LMX, &
            modeldatactx%DX, &
            modeldatactx%DMAT, &
            modeldatactx%SKEW)

if (iecho.eq.8) then
     call parallelLog("inital", "computing strain")
endif


call STRAIS (meshdatactx%X, &
             meshdatactx%IEN, &
            modeldatactx%D, &
            modeldatactx%STR, & ! strain goes here &
            modeldatactx%LMF, &
            modeldatactx%TFAULT, &
            modeldatactx%LMX, &
            modeldatactx%DX, &
            modeldatactx%SKEW)

if (iecho.eq.8) then
     call parallelLog("inital", "printing stress and strain")
endif


! Print stresses 
call PRINTS (modeldatactx%STN,NSTR,NUMEL)

! Output strains
call PRINTE (modeldatactx%STR,NSTR,NUMEL)

300 continue

! Force refactoring in first viscous step
REFACT = .true.

if (iecho.eq.5) then
    write(*,*) 'rank', getrank(), ' has done solvecaller'
endif

if (iecho.eq.8) then
     call parallelLog("inital", "finished writing stuff to fedsk files")
endif

return
end subroutine
!-------------------------------------------------------------------------------

subroutine ADDFOR (B,P,LM,NEE, elemID, callerID)

use constants,      only: one
use debugmodule,    only: iecho
use modelctx,       only: getrank ! debug
use meshdatamodule, only: meshdatactx ! debug
#ifdef SPARSE
use AOmodule,       only: equations ! debug
#endif

! P = Local mechanical load vector.
! LM = LM or LMX, depending on caller, so cannot LM from modeldata

implicit none
!-pass
integer NEE
!    integer elemID
integer :: LM
double precision :: B, P(24)
dimension B(*),LM(*)

integer :: elemID, callerID ! have been added for debugging purposes.

double precision, external :: FLOAT   ! from tools
integer  :: j,k, rank
double precision :: temp    

 
 rank = getrank()

#ifdef Lukas
!write(*,*) "addfor says: LM: ", LM(1:NEE)
!write(*,*) "addfor says: P: ", P(1:NEE)
#endif

do j=1,NEE
    k = IABS(LM(j))

    if (k.ne.0) then
        ! This sign is relevant with split nodes, 
        ! when there, LM(j) can become negative.
        ! Otherwise, LM is always positive.
        !      write(*,*) 'addfor says: B(k) ', B(k)

        if (callerID.eq.1 .or. &
            callerID.eq.3 .or. &
            callerID.eq.5 .or. &
            callerID.eq.7.or. &
            callerID.eq.11 .or. &
            callerID.eq.13 .or. &
            callerID.eq.18) then
            ! contributions from slippery nodes
            ! LMX is given as argument, passed as LM, hence LM in equation

!#ifdef Lukas
!            write(*,*) "addfor slip set B(",k,") from ", B(k), &
!                       "by adding ", P(j), "times", SIGN(one,FLOAT(LM(j))), &
!                       "to",B(k) + P(j) * SIGN(one,FLOAT(LM(j)))
!#endif

            B(k) = B(k) + P(j) * SIGN(one,FLOAT(LM(j)))

!            B(k) = B(k) + P(j) 
        else
            ! contributions from regular nodes
!#ifdef Lukas
!            if (callerID.eq.4) then
!            write(*,*) "rank", getrank(), &
!                       "addfor regular set B(",k,") from ", B(k), &
!                       "by adding ", P(j), &
!                       "to",B(k) + P(j) * SIGN(one,FLOAT(LM(j)))
!            endif
!#endif


            B(k) = B(k) + P(j) * SIGN(one,FLOAT(LM(j)))
        endif

!        B(k) = B(k) + P(j)

    endif
enddo


return
end

!-------------------------------------------------------------------------------
subroutine AddFaultDisplacement (DL,elemID,DFAULT,NDOF,NEN, callerID)

USE MODELDEFINITION, only: NUMFN
use modeldatamodule, only: modeldatactx

! Adds displacements at faulted nodes to the local displacement vector

implicit none
!-pass
integer          :: NDOF,NEN
integer          :: elemID
integer          :: callerID ! every calling function has a different integer,
                             ! to trace when things go wrong.
double precision :: DL(NDOF,NEN)
double precision :: DFAULT(NDOF,*)
!-locl
integer :: i,j,k

if (NUMFN.le.0) then
    return
endif

do i=1,NEN
    k = modeldatactx%LMF(1,i,elemID) ! retrieve index of the FAULT array
    if (k.gt.0) then
        do j=1,NDOF
            DL(j,i) = DL(j,i) + DFAULT(j,k)
        enddo
    endif
enddo

return
end
!-------------------------------------------------------------------------------
#ifdef EXPERIMENTAL_ReducedSlipAlloc
subroutine ADDSN (DL,DX,IEN,LMX,SKEW,NDOF,NEN, elemID)
#else
subroutine ADDSN (DL,DX,IEN,SKEW,NDOF,NEN, elemID)
#endif



USE MODELDEFINITION
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: modeldatactx, &
                           elementSlipperyID
#else
use modeldatamodule, only: modeldatactx
#endif
USE MODELCTX
USE ALGEBRA,         only: formrt, clear, transp, equate
use constants,       only: zero, one
use debugmodule,     only: xit
use iomodule

!
! adds displacements due to slip across internal interfaces to the local
! displacement vectors
!

implicit none
!-pass
integer NDOF,NEN
integer :: ien
#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer :: lmx
dimension :: LMX(NEE)
#endif
double precision :: dl, dx, skew
dimension DL(NDOF,*),DX(NDOF,*),IEN(NEN),SKEW(2,*)

!-locl
integer i,k,j, m
double precision :: ROT(3,3),DXL(3),TEMP(3)
external FLOAT
integer, intrinsic :: SIGN
double precision :: sum
integer :: elemID
integer :: sgn

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer :: slipperySequenceNr
#endif

!
!#ifdef Lukas
!write(*,*) "###### addsn starts with DL: ", DL(1:3,1:3)
!#endif
!write(*,*) "###### addsn starts with IEN: ", IEN
!write(*,*) "###### addsn starts with LMX: ", LMX

!write(*,*) "##### computing alternaitive elemID", elemID
!write(*,*) "##### computing alternaitive LMX 1",size(modeldatactx%LMX,1)
!write(*,*) "##### computing alternaitive LMX 2",size(modeldatactx%LMX,2)
!write(*,*) "##### computing alternaitive LMX 3",size(modeldatactx%LMX,3)
!write(*,*) "###### addsn has alternative LMX ", modeldatactx%LMX(:,:,elemID)


! DL = localised displacement
! NEN = number of vertices per element

do j=1,NEN
    k = IEN(j)

!    write(*,*) "addsn has NUMROT: ", NUMROT


    if (NUMROT.gt.0) then
!        back-rotate diff displacement vector to local skewed system
        call CLEAR (DXL,NDOF,"DXL")

        if (modeldatactx%SKEW(1,k).ne.0.0 .or. &
            modeldatactx%SKEW(2,k).ne.0.0) then

            call FORMRT (modeldatactx%SKEW(:,k),ROT,ndof)
            call TRANSP (ROT,3)
            call EQUATE (TEMP,DX(1:NDOF,k),NDOF)

            do i=1,NDOF
                sum = ZERO
                do m=1,NDOF
                    sum = sum + ROT(i,m)*TEMP(m)
                enddo
                DXL(i) = sum
            enddo

        else
            call EQUATE (DXL,DX(1:NDOF,k),NDOF)
        endif


!    write(*,*) "addsn says: A DXL: ", DXL

        do i=1,NDOF

#ifdef EXPERIMENTAL_ReducedSlipAlloc
            slipperySequenceNr = elementSlipperyID(elemID,20)
            if(slipperySequenceNr.gt.0) then
                if (modeldatactx%LMX(i, j, slipperySequenceNr).ne.0) then
                    DXL(i) = DXL(i) * &
                     SIGN(ONE,dble(modeldatactx%LMX(i, j, slipperySequenceNr)))
                endif
            endif
#else
            if (modeldatactx%LMX(i, j, elemID).ne.0) then
                DXL(i) = DXL(i) * SIGN(ONE,dble(modeldatactx%LMX(i, j, elemID)))
            endif
#endif
        enddo



!    write(*,*) "addsn says: AA DXL: ", DXL


        if (SKEW(1,k).ne.0.0 .or. SKEW(2,k).ne.0.0) then
            call TRANSP (ROT,3)
            call EQUATE (TEMP,DXL,NDOF)

            do i=1,NDOF
                sum = ZERO
                do m=1,NDOF
                    sum = sum + ROT(i,m)*TEMP(m)
                enddo
                DXL(i) = sum
            enddo

        endif

!    write(*,*) "addsn says: B DXL: ", DXL


        do i=1,NDOF
!        write(*,*) "**** addsne  ", i,j, DL(i,j), " + " , i, DXL(i)
            DL(i,j) = DL(i,j) + DXL(i)
        enddo
!!
    else 
    ! numrot = 0
       ! add slippery displacement to regular displacement
        do i=1,NDOF
            !l = i+(j-1)*NDOF
#ifdef EXPERIMENTAL_ReducedSlipAlloc
            slipperySequenceNr = elementSlipperyID(elemID,21)
            if(slipperySequenceNr.gt.0) then
                if (modeldatactx%LMX(i, j, slipperySequenceNr).ne.0) then
                    sgn = SIGN(ONE,dble(modeldatactx%LMX(i, j, slipperySequenceNr)))
                    DL(i,j) = DL(i,j) + sgn*DX(i,k)
                endif
            endif
#else
            if (modeldatactx%LMX(i, j, elemID).ne.0) then
                sgn = SIGN(ONE,dble(modeldatactx%LMX(i, j, elemID)))
                DL(i,j) = DL(i,j) + sgn*DX(i,k)
            endif
#endif
        enddo
!!
      endif
 
enddo

!        write(*,*) "addsn says DL: ", DL(1:3,1:8)


#ifdef Lukas
!write(*,*) "###2### addsn starts with DL: ", DL(1:2,1:4)
!write(*,*) "inbetweenies: ", NTYPE, IEN(1:8)
#endif

if (NTYPE.eq.1) then
    if (IEN(3).eq.IEN(4)) then
        DL(1,4) = DL(1,3)
        DL(2,4) = DL(2,3)
    endif
else if (NTYPE.eq.2) then
    if (IEN(3).eq.IEN(4)) then
        DL(1,4) = DL(1,3)
        DL(2,4) = DL(2,3)
        DL(3,4) = DL(3,3)
    endif
else if (NTYPE.eq.5) then
    if (IEN(3).eq.IEN(4) .and. & 
        IEN(7).eq.IEN(8)) then
        DL(1,4) = DL(1,3)
        DL(2,4) = DL(2,3)
        DL(3,4) = DL(3,3)
        DL(1,8) = DL(1,7)
        DL(2,8) = DL(2,7)
        DL(3,8) = DL(3,7)
    endif
else if (NTYPE.eq.6) then
    continue
else
    write(stderr,3) NTYPE
 3        format(1x,'ADDSN: NTYPE = ',I1,' not supported')
    call xit(1," ")
endif

#ifdef Lukas
!write(*,*) "###3### addsn starts with DL: ", DL(1:2,1:4)
#endif

!
return
end
!-------------------------------------------------------------------------------
subroutine ADDSNE (DL,elemID,NDOF,NEN,whichOne)

! Adds slippery node displacements (in DXE) to local displacement vector (in DL)
! or slippery node velocity (in DELDXE) to local velocity vector (also DL here)

use constants,       only: useDXE
USE MODELDEFINITION, only: NUMSLP
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: modeldatactx, &
                           elementSlipperyID
#else
use modeldatamodule, only: modeldatactx
#endif

implicit none

!-pass
integer          :: whichOne   !   0 = DXE; 1 = DELDXE
integer          :: NDOF,NEN,elemID
double precision :: DL(NDOF,NEN)
!-locl
integer          :: i,j
integer          :: slipperySequenceNr

if (NUMSLP.eq.0) then
    ! the values in DXE will all be zero. No need to add.
    return
endif

!write(*,*) "calling addsne with elemID", elemID
!write(*,*) "called for whichone ", whichOne


#ifdef EXPERIMENTAL_ReducedSlipAlloc
slipperySequenceNr = elementSlipperyID(elemID,22)
#endif

if (whichOne.eq.useDXE) then
! displacement
#ifdef EXPERIMENTAL_ReducedSlipAlloc
    
    if(slipperySequenceNr.gt.0) then
        do j=1,NEN
            do i=1,NDOF
                DL(i,j) = DL(i,j) + modeldatactx%DXE(i,j,slipperySequenceNr)
            enddo
        enddo 
    endif
#else
    do j=1,NEN
        do i=1,NDOF
            DL(i,j) = DL(i,j) + modeldatactx%DXE(i,j,elemID)
        enddo
    enddo
#endif
else

#ifdef EXPERIMENTAL_ReducedSlipAlloc
    if(slipperySequenceNr.gt.0) then
        do j=1,NEN
            do i=1,NDOF
                DL(i,j) = DL(i,j) + modeldatactx%DELDXE(i,j,slipperySequenceNr)
            enddo
        enddo
    endif
#else
    do j=1,NEN
        do i=1,NDOF
            DL(i,j) = DL(i,j) + modeldatactx%DELDXE(i,j,elemID)
        enddo
    enddo
#endif



endif



return
end
!-------------------------------------------------------------------------------
 subroutine ADDSRF (ISURF,SURF,X,D,NSD,NDOF,NSURF,IFLAG)
!
! Note that faulted or slippery displacements of surface nodes will not
! be passed into the SURF arrays, because faulted and slippery displacements
! are element based.
!
implicit none
!-pass
integer          :: NSD,NDOF,NSURF,IFLAG, ISURF
double precision :: surf, x, D
dimension        :: ISURF(*),SURF(NSD,*),X(NSD,*),D(NDOF,*)
!-local
integer          :: i, j, n

if (NSURF.le.0) then
    return
endif

do j=1,NSURF
    n = ISURF(j)
    do i=1,NSD
        if (IFLAG.eq.0) then
            SURF(i,j) = X(i,n) + D(i,n)
        else
            SURF(i,j) = SURF(i,j) + D(i,n)
        endif
    enddo
enddo

return
end
!-------------------------------------------------------------------------------
subroutine ADDSTF (A,SL,IDIAG,LM,LMX,NEE,NUMSN,globalElemID)

USE MODELCTX
use modeldefinition, only: neq
use modeldatamodule, only: HasSlipperyNodes, &
                           modeldatactx
use debugmodule,     only: iecho, &
                           parallellog, &
                           debug, delay
!use meshdatamodule,  only: meshdatactx
use constants,       only: one
!use petscksp

#ifdef SPARSE
use AOmodule,        only: equations, &
                            equationsmaskx
use spetscmodule,     only: PETSc_MatSetValuesLocal, M_A
#endif

! Program to add element stiffness to global stiffness matrix

! LM and LMX of only a single element is passed on to this
! subroutine. If the element has slippery nodes, LMX has the 
! proper references to the equationsx array. If not, it is full
! of zeros. Effectively this is the same as it was prior
! to slippery element allocation refactoring


 implicit none
!-pass
! NEE = NEN * NDOF
integer   :: NEE,NUMSN,LM(NEE)

integer :: LMX(NEE)

integer       :: IDIAG(NEQ)
double precision A, SL(NEE,NEE)
dimension :: A(*)
integer :: globalElemID

#ifdef SPARSE

#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvecdef.h"
!#include "petsc/finclude/petscmatdef.h"
!#include "petsc/finclude/petsckspdef.h"

   PetscMPIInt irank
#endif
!-local
integer i,j,k,l,m
external FLOAT
#ifdef SPARSE
double precision :: v ! value to be set into stiffness matrix
#endif
double precision :: sgn
integer          :: ierr

#ifdef SPARSE

!write(*,*) "ADDSTF called for elem", globalElemID, "with LMX", LMX


! step through the local stiffness matrix, and if
! positions in tht matrix correspond to a position
! in the global matrix (usually the case; only for position/
! velocity boundary conditions, lines and columns of the local
! stiffness matrix are ommited. (see the checks k==0 and m==0)

do j=1,NEE

    k = LM(j)

    if (k.eq.0) then
        cycle ! skip in case no d.o.f. (boundary cond...)
    endif

    do i=1,NEE

        m = LM(i)

        if (m.eq.0) then
            goto 100  ! skip in case no d.o.f. (boundary cond...)
        endif
        v = SL(i,j)

        call PETSc_MatSetValuesLocal(M_A,k-1,m-1,v,1,globalElemID,i,j,ierr)

100     if (.not. HasSlipperyNodes()) then
            cycle
        endif

        m = IABS(LMX(i))

        if (m.eq.0) then
            cycle
        endif

!if (getrank().eq.0) then
! USEFUL FOR DEBUGGING. PLEASE LEAVE IN COMMENT.
!        write(*,*) "addstf 2; LMX: ", LMX

!        write(*,*) "rank", getrank(), &
!                   "element ", globalElemID

!        write(*,*) "rank", getrank(), globalElemID, &
!                   "has local ID ", meshdatactx%Glo2LocElement(globalElemID) 

!        write(*,*) "rank", getrank(), &
!                   "has nodes", meshdatactx%ien(:,meshdatactx%Glo2LocElement(globalElemID))

!        write(*,*) "rank", getrank(), &
!                   "element ", globalElemID, &
!                   "has pos in IEN", meshdatactx%Glo2LocElement(globalElemID) - &
!                                     meshdatactx%elementOffset(getrank()+1)

!endif


        v = SL(i,j)*SIGN(one,dble(LMX(i)))
!        write(*,*) "Found v", v, "from ", i, j, SL(i,j), "and LMX", LMX(i)
        call PETSc_MatSetValuesLocal(M_A,k-1,m-1,v,2,globalElemID,i,j,ierr)
    enddo
enddo

#else

do j=1,NEE

    k = LM(j)

    if (k.eq.0) then
        cycle
    endif
    L = IDIAG(k) - k
    do i=1,NEE

        m = LM(i)

        if (m.eq.0 .or. m.gt.k) then
            goto 100
        endif
        m = L + m
        A(m) = A(m) + SL(i,j)
100     if (.not. HasSlipperyNodes()) then
            cycle
        endif

        m = IABS(LMX(i))

        if (m.eq.0 .or. m.gt.k) then
            cycle
        endif
        m = L + m
        A(m) = A(m) + SL(i,j)*SIGN(one,dble(LMX(i)))
    enddo
enddo
#endif

!********* slippery nodes ****************

if (.not. HasSlipperyNodes()) then
    goto 2200
endif

#ifdef SPARSE

    do j=1,NEE

        k = IABS(LMX(j))

        if (k.eq.0) then
            cycle
        endif

        sgn = SIGN(one,dble(LMX(j)))

        do i=1,NEE

            m = LM(i)

            if (m.eq.0) then
                goto 600
            endif

            v = sgn*SL(i,j)
            call PETSc_MatSetValuesLocal(M_A,k-1,m-1,      v,5,globalElemID,i,j,ierr)
!            call PETSc_MatSetValuesLocal(M_A,k  ,m  ,      v,5,globalElemID,i,j)

600         m = IABS(LMX(i))

            if (m.eq.0) then
                cycle
            endif


        v = sgn*SL(i,j)*SIGN(one,dble(LMX(i)))


!            write(*,*) "prepping 6; ", sgn, dble(LMX(i)), SIGN(one,dble(LMX(i))), SL(i,j), "v", v


        call PETSc_MatSetValuesLocal(M_A,k-1,m-1,v,6,globalElemID,i,j,ierr)


!            v = sgn*SL(i,j)*SIGN(one,dble(LMX(i)))

!            call PETSc_MatSetValuesLocal(M_A,k-1,m-1,      v,6,globalelemID,i,j,ierr)
!            call PETSc_MatSetValuesLocal(M_A,k  ,m ,       v,6,globalelemID,i,j)

    enddo
enddo

#else

    do j=1,NEE

        k = IABS(LMX(j))

        if (k.eq.0) then
            cycle
        endif
        L = IDIAG(k) - k
        sgn = SIGN(one,dble(LMX(j)))

        do i=1,NEE

            m = LM(i)
            if (m.eq.0 .or. m.gt.k) then
                goto 600
            endif
            m = L + m
            A(m) = A(m) + sgn*SL(i,j)

600         m = IABS(LMX(i))
            if (m.eq.0 .or. m.gt.k) then
                cycle
            endif
            m = L + m
            A(m) = A(m) + sgn*SL(i,j)*SIGN(one,dble(LMX(i)))

    enddo
enddo

#endif

2200 return

end

!-------------------------------------------------------------------------------

subroutine CONDNR (A,IDIAG,NEQ,ierr)

use constants, only: dfmin,dfmax,eps
use iomodule          

!    Routine to determine the formal conditioning number of the
!    array D of factorization of A = U(TRANSPOSE) * D * U
!    return status ierr,    0=OK
!                1=error

implicit none
!-pass
integer NEQ
double precision :: A
integer :: IDIAG(NEQ)
dimension A(*)
integer :: ierr
!-local
double precision :: wmin, wmax, CondNo
double precision :: stfns
integer          :: i, k

!-init
wmin =  dfmax
wmax = -dfmax
ierr = 0

do i=1,NEQ
    k = IDIAG(i)
    stfns = A(k)
    if (stfns.lt.wmin) then
        wmin = stfns
    endif
    if (stfns.gt.wmax) then
        wmax = stfns
    endif
enddo

if (ABS(wmax).lt.dfmin) then
    CondNo = dfmax
else
    CondNo = wmin/wmax
endif

if (CondNo.lt.EPS) then
    write(stdout,1) CondNo
 1        format(///1x,'*** ILL-CONDITIONED STIFFNESS MATRIX ***'/ &
     5X,'RECIPROCAL FORMAL CONDITION NUMBER = ',1PE12.2//)
    if (OUTFIL(1)) then
        write(luout(1),1) CondNo
    endif
    ierr = 1
endif

return
end

!-------------------------------------------------------------------------------
subroutine DISBC(DL,P,S,NEE)

use constants, only: dfmin

!.... PROGRAM TO ADJUST LOAD VECTOR FOR PRESCRIBED DISPLACEMENT
!     BOUNDARY CONDITIONS

implicit none
!-pass

integer          :: NEE, i, j
double precision :: DL(NEE),P(NEE),S(NEE,NEE)


do j=1,NEE
    if (ABS(DL(j)).gt.dfmin) then
        do i=1,NEE
            P(i) = P(i) - S(i,j)*DL(j)
        enddo
    endif
enddo


return
end

!-------------------------------------------------------------------------------
subroutine DISP (B,D,ID,LINK,NDOF,NUMNP,NLINK)

USE MODELDATAMODULE, only: modeldatactx, &
                           inlink
USE MESHDATAMODULE,  only: meshdatactx
USE MODELDEFINITION, only: NEQglobal
USE MODELCTX,        only: getrank

!***** Transfers nodal displacement increments into from the solution vector B into D 

implicit none
!-pass
integer NDOF,NUMNP,NLINK
double precision :: B, D
integer          :: ID, LINK
dimension B(*),D(NDOF,*),ID(NDOF,*),LINK(3,*)
!-local
integer :: i, j, k, m


! is this check of the node is linked necessary? isnot it so that upon
! reading of input, in subroutine local, that the equation number of a
! linked node j and dof i is set to that of its master (linked-to) node
! at LINK(3,n) where j=LINK(2,n)?

!      do j=1,meshdatactx%Nvlocal
!          do i=1,NDOF
!     check if ndof i, node j is linked
!              m = INLINK(LINK,NLINK,i,
!     .            meshdatactx%gloInd(i))
!              if (m.le.0) then
!                  k = modeldatactx%ID(i,j)
!              else
!                  k = modeldatactx%ID(i,unknown
!              endif
!              if (k.gt.0) D(i,j) = B(k)
!          enddo
!      enddo

!write(*,*) 'rank', getrank(),' DISP has numnp: ', numnp, ' and ID ', ID(1:2, 1:meshdatactx%nvlocal)



do j=1,meshdatactx%nvlocal
    do i=1,NDOF

!       is this check necessary? isnot ID(i,j) automatically correct in
!       case of a linked node? see routine tecin.addlinkednodes

! possible need to get motion from another partition, when a linked node
! is located in another partition?...

        m = INLINK (i,j)     !LINK,NLINK,i,j)

        if (m.le.0) then
            k = ID(i,j)
        else
            k = ID(i,LINK(3,m))
        endif
    
        if (k.gt.0) then 
!            write(*,*) 'disp putting ', j, i, k, B(k)
            D(i,j) = B(k)
        endif

    enddo
enddo

return
end

!-----------------------------------------------------------------------

subroutine FORMF (B,X,SKEW,IEN,LM,D,LMX,MAT,PRPMAT,DMAT,LMF, &
    NFAULT,DFAULT,TFAULT,DXE,NSD,NDOF,NSTR,NEN)

USE MATERIALSMODULE, only: LMAT, MaterialIDOfElement
USE MESHDATAMODULE

#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: elementSlipperyID
#endif

USE MODELDEFINITION
use modeldatamodule, only: HasSlipperyNodes
USE MODELCTX
USE ALGEBRA
use constants,       only: useDXE
use timestepmodule,  only: nstep, faultvelo

use iomodule          
!
! Generates forces due to faulted node displacements (in DFAULT)
!
!    logical VELBC is used to change interpretation of the data in FAULT
!    VELBC = .true.: FAULT means constant velocity (NFAULT(3,.) is not used)
!
!    VELBC = .false.: FAULT means differential displacement
!                     NFAULT(3,.) >= 0 differential displacement applied @ time step NFAULT(3,.)
!                     NFAULT(3,.)  < 0 periodic displacement applied every ABS(NFAULT(3,.) time step 
!

implicit none
!-pass
integer NSD,NDOF,NSTR,NEN
integer :: IEN, LM, LMX, LMF, NFAULT, MAT
double precision :: X, B, PRPMAT, SKEW, DFAULT, DMAT, D, TFAULT, DXE
dimension X(NSD,*),B(*),PRPMAT(6,*),SKEW(2,*),MAT(*),IEN(NEN,*), &
    LM(NDOF,NEN,*),NFAULT(3,*),DFAULT(NDOF,*),DMAT(NSTR,NSTR,*), &
    LMX(NDOF,NEN,*),D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
    DXE(NDOF,NEN,*)
!-locl
double precision :: S(576),STEMP(576),P(24),XL(24),DL(24)
logical doit,QUAD,VELBC
double precision :: thick
integer          :: i, ii, j, k, n, m
integer          :: ierr
integer          :: vertexID,nodeID

!-init
vertexID = 0
nodeID = 0
VELBC = (faultvelo.ne.0)

do ii=1,NUMFN   ! NUMFN is partition-local number of slippery elements and nodes (AOmodule unclear on ghost nodes and ghost elements)
    i = ii
    ! is a force update needed in the current timestep (doit)?
    if (VELBC) then
        doit = .true.
    else
        if (NFAULT(3,i).ge.0) then
            doit = (NSTEP.eq.NFAULT(3,i))
        else
            doit = (MOD(NSTEP,-NFAULT(3,i)).eq.0)
        endif
    endif

    if (doit) then
        n = IABS(NFAULT(1,i))           ! partition-local element number
        nodeID = NFAULT(2,i)            ! partition-local node number
        m = MaterialIDOfElement(n)
        THICK = PRPMAT(6,m)
        QUAD=(IEN(3,n).ne.IEN(4,n))

        call LCOORD (meshdatactx%X, XL, meshdatactx%IEN(1,N),0)
        if (LGDEF.ne.0) then
            call LDISP (DL,D,IEN(1,n),NDOF,NEN)
            call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,49)
            call ADDSNE (DL,n,NDOF,NEN,useDXE)
            call REZONE (XL,DL,1d0)
        endif

!       calculate local stiffness matrix
        call CLEAR(S,NEE*NEE,"S")   ! NEE=NDOF*NEN
        call STIFF (DMAT(1,1,n),XL,S,THICK,QUAD,ierr, nodeID)
        if (ierr.ne.0) then
            write(stderr,1) n
1           format(1x,'formf: error in element ',I5)
            call stoper()
        endif
        call SYMMET(S,NEE)

!       perform rotation of local axes to global axes directions
        if (NUMROT.ne.0) then
            call RSTIFF(S,STEMP,IEN,SKEW,NDOF,NEN, NEE,n,n)
        endif

#ifdef DEBUG
        write(*,*) 'formf: element stiffness matrix S:'
        do i=1,NEE
            write(*,2) ((S((i-1)*NEE+j),j=1,NEE)
2           format(576(1PE15.2,1X))
        enddo
#endif

!       now find which vertex (1,2,..,NEN) this is of the current element
        vertexID = 0
        do j=1,NEN
            if (IEN(j,n).eq.nodeID) then
                vertexID = j
            endif
        enddo
        if (vertexID.eq.0) then
        endif

!       load fault displacements into DL
        call CLEAR(DL,NDOF*NEN,"DL")
!       get DFAULT index of this particular node and element
        k = LMF(1,vertexID,n)      
        call LFLTEQ (DL,NFAULT(1,ii),DFAULT(1,k),IEN(1,n),NEN,NDOF)

#ifdef DEBUG
        write(*,*) 'formf: split node displacements loaded into DL: ', DL(1:NEE)
#endif

!       compute load contribution from this split node (P = -S.DL)
        call CLEAR(P,NEE,"P")
        call VMPRD(S,DL,P,NEE,NEE)
        do k=1,NEE
            P(k) = -P(k)
        enddo

#ifdef DEBUG
        write(*,*) 'formf: load contribution by ', P(1:NEE)
#endif

!       add the element-local load to normal dofs of the (partition-)global load vector
        call ADDFOR (B,P,LM(1,1,n),NEE,ii,0)

        if (HasSlipperyNodes()) then 
#ifdef EXPERIMENTAL_ReducedSlipAlloc
            if(elementSlipperyID(n,23).gt.0) then
                call ADDFOR (B,P,LMX(1,1,elementSlipperyID(n,23)),NEE,ii,1)
            endif
#else
!           add the element-local load to slippery dofs of the (partition-)global load vector
            call ADDFOR (B,P,LMX(1,1,n),NEE,ii,1)
#endif
        endif

    endif ! doit

enddo ! loop over NFAULT entries

return
end
!-------------------------------------------------------------------------------
subroutine LFLTEQ(DL,NFAULT,DFAULT,IEN,NEN,NDOF)
!
!    subroutine to localize faulted node displacements
!
implicit none
!-pass
integer   :: NEN, NDOF

double precision :: DL, DFAULT
integer   :: nfault, ien
dimension :: DL(ndof,*),IEN(*),DFAULT(*),NFAULT(*)
!-locl
integer i,j
!
do i=1,NEN
    if (NFAULT(2).eq.IEN(i)) then
        do j=1,NDOF
            DL(j,i) = DFAULT(j)
        enddo
    endif
enddo
!
return
end

!-------------------------------------------------------------------------------
subroutine FORMK (A,IDIAG,B,X,D,SKEW,IEN,LM,LMX,MAT, &
                     DMAT,S,IA,IAX,NSD,NDOF,NSTR,NEN)

USE MATERIALSMODULE, only : LMAT, &
                          matdatactx, &
                          MaterialIDOfElement
USE MODELDEFINITION

#ifdef EXPERIMENTAL_ReducedSlipAlloc
   use modeldatamodule, only: hasSlipperyNodes, &
                              modeldatactx, &
                              elementSlipperyID, &
                              meq, showfbc
#else
   use modeldatamodule, only: hasSlipperyNodes, &
                              modeldatactx, &
                                       meq, showfbc
#endif

USE MODELCTX
USE ALGEBRA
USE MESHDATAMODULE
use debugmodule,     only: lustat, &
                           itest, &
                           iecho, &
                           debug, &
                           parallelLog, &
                           debugFileID, &
                           delay

#ifdef SPARSE
use AOmodule, only: equationsmask, equationsmaskx
#endif
use iomodule          

! Routine to form the stiffness matrix K and add forces due to
! prescribed displacement boundary conditions do the load vector.

! this routine will call the function addstf, which inserts the local stiffness contributions
! to the global stiffness matrix in PETSc


implicit none
!-pass
integer NSD,NDOF,NEN,NSTR


integer :: IDIAG, IEN, LM, LMX, IA(ndof,numnp), IAX(ndof,numnp), mat
double precision :: A, B, X, D, SKEW, S
dimension :: A(*),IDIAG(NEQ),B(*),X(NSD,*),D(NDOF,*),SKEW(2,*), &
             IEN(NEN,*),LM(NDOF,NEN,*),LMX(NDOF,NEN,*),MAT(*), &
             S(MEQ,*)
double precision :: DMAT(NSTR,NSTR,NUMEL)

!-locl
character bs
logical QUAD
! this arrays have rather large size, to make them compatible
! with cubic meshes. With tetrahedrons, the higher indices will
! not be used.
! SL = Stiffness Local
! DL = Displacement local
! local = only for one particular element
double precision :: SL(576),STEMP(576),P(24),DL(24),XL(24)
integer :: i, j, n, nn, m, ierr
double precision :: E, Pois, thick
integer, parameter :: dummyLMX(8) = (/0,0,0,0,0,0,0,0/)

#ifdef EXPERIMENTAL_ReducedSlipAlloc

integer :: slipSequenceNr
integer :: globalElmID

#endif

!-init
bs = CHAR(8)
call CLEAR (DMAT,NSTR*NSTR*NUMEL,"DMAT")

if (debug) then
    !write(*,*) "rank",getrank(),"LM used in formk: ", modeldatactx%LM
endif

if (iecho.eq.8) then
    call ParallelLog("elastc", "Entering formk element loop")
endif


if (getrank().eq.0) then
!   write(*,*) "priori incantatum", modeldatactx%LMX
#ifdef SPARSE
!   write(*,*) "eqmask: ", equationsmask
!   write(*,*) "eqmaskx: ", equationsmaskx
#endif
endif

do nn=1,meshdatactx%nelocal


    if (HasSlipperyNodes()) then

#ifdef EXPERIMENTAL_ReducedSlipAlloc

        ! n goes from 1 to the local element count.
        ! Note that this is not the same as the local element number,
        ! which count continuously through the partitions.
    
        ! We have the LMX array sorted by global element ID.
        globalElmID = meshdatactx%gloElt(nn)
        !write(*,*) "gloelt: ", meshdatactx%gloelt
        !write(*,*) "gives globalelm", nn, globalElmID
        slipSequenceNr = elementSlipperyID(globalElmID,24)

#endif
    endif


    if (debug) then
    !   write(debugFileID,*) "----------------------------------------------------"
    !   write(debugFileID,*) "rank",getrank(),"Adding stiffness for element ", nn
    !   write(debugFileID,*) "----------------------------------------------------"
    endif

    n = nn
    call CLEAR(SL,NEE*NEE,"SL")
    call CLEAR(P,NEE,"P")

!   localize element coordinates
    call LCOORD (X,XL,IEN(1,n),nn)

!   construct the local material matrix and insert it
!   in the global array DMAT(NSTR,NSTR,NUMEL)
    m = MaterialIDOfElement(n)
    E     = matdatactx%PRPMAT(1,m)
    POIS  = matdatactx%PRPMAT(2,m)
    THICK = matdatactx%PRPMAT(6,m)
    call MATERL (DMAT(1,1,n),E,POIS)

!   construct the upper local stiffness matrix
    QUAD = (IEN(3,n).ne.IEN(4,n))
    call STIFF (DMAT(1,1,n),XL,SL,THICK,QUAD,ierr, n)
    if (ierr.ne.0) then
        write(stderr,1) n
 1           format(1x,'formk: error in element ',I5)
        call stoper()
    endif

!   symmetrize it
    call SYMMET (SL,NEE)

!   rotate stiffness for skew bcs
    if (NUMROT.gt.0) then
        if (debug) write(*,*) "b rank", getrank(), "NUMROT: ", NUMROT

        ! n, element ID, will be used to get vertex from IEN array, and as such must use local numbering
        call RSTIFF (SL,STEMP,IEN,SKEW,NDOF,NEN,NEE,n,meshdatactx%locElt(nn))
    endif
    
    call LDISBC (DL, modeldatactx%D, meshdatactx%IEN(1,n), n, NDOF, NEN, meshdatactx%gloElt(nn))

    ! build P from DL and SL

    if (debug) then
        write(*,*) "rank", getrank(), "for elem",nn,"calls DISBC with SL ", SL(1:64)
        write(*,*) "rank", getrank(), "for elem",nn,"calls DISBC with DL ", DL(1:8) 
    endif

    call DISBC (DL,P,SL,NEE)

    if (debug) then
        write(stdout,30) n
30      format(/1x/1x,'Local stiffness matrix of element ',I5)
        call putmat (SL,NEE,NEE,stdout)
        write(stdout,40) n
40           format(/1x,'Local load vector of element ',I5)
        call putvec (P,NEE,stdout)
    endif

!   add the local arrays to the global stiffness and load
    call ADDFOR (B,P,LM(1,1,N),NEE,meshdatactx%gloElt(nn),2)   ! from disp bcs

    if (HasSlipperyNodes()) then

!       write(*,*) "rank", getrank(), "enter slippery nodes, calling addfor"

#ifdef EXPERIMENTAL_ReducedSlipAlloc

        if (slipSequenceNr.gt.0) then
            ! TODO fix call and tine
            call ADDFOR (B,P,modeldatactx%LMX(1,1,slipSequenceNr),NEE,nn,3)
        else
            ! the element has no slippery nodes, no need to add slippery stuff
        endif
#else
        call ADDFOR (B,P,modeldatactx%LMX(1,1,nn),NEE,nn,3)
#endif

        if (debug) write(*,*) "rank", getrank(), "leaves slippery nodes, after addfor"

    endif

    ! send global element nr as last argument for debugging.
    ! the enries in the SL matrix are inserted into the stiffness matrix


#ifdef EXPERIMENTAL_ReducedSlipAlloc

!    write(*,*) "rank", getrank(),"----------------------------------------------------------"
!    write(*,*) "rank", getrank(),"calling addstf for element", n, slipSequenceNr
!    write(*,*) "rank", getrank(),"with LM :", modeldatactx%LM(:,:,n)
    if (slipSequenceNr.gt.0) then
!        write(*,*) "rank", getrank(),"with LMX:", modeldatactx%LMX(:,:,slipSequenceNr)
    else
!        write(*,*) "rank", getrank(),"with LMX:",dummyLMX
    endif

    if (slipSequenceNr.gt.0) then
!        write(*,*) "addstf called wtih slipseq", slipSequenceNr
        call ADDSTF (A, &
                 SL, &
                 IDIAG, &
                 modeldatactx%LM(1,1,n), &
                 modeldatactx%LMX(1,1,slipSequenceNr), & 
                 NEE, &
                 NUMSN, &
                 meshdatactx%gloElt(n))
    else
!        write(*,*) "addstf called wtih dummy"
        call ADDSTF (A, &
                 SL, &   
                 IDIAG, &
                 modeldatactx%LM(1,1,n), &
                 dummyLMX, &
                 NEE, &
                 NUMSN, &
                 meshdatactx%gloElt(n))
    endif

!    write(*,*) "rank", getrank(),"----------------------------------------------------------"

#else

    if (hasSlipperyNodes()) then
        call ADDSTF (A,SL,IDIAG,modeldatactx%LM(1,1,n),modeldatactx%LMX(1,1,n),NEE, &
           NUMSN, meshdatactx%gloElt(n))
    else
        call ADDSTF (A,SL,IDIAG,modeldatactx%LM(1,1,n),dummyLMX,NEE, &
           NUMSN, meshdatactx%gloElt(n))
    endif


#endif /* EXPERIMENTAL_ReducedSlipAlloc */


    if (showfbc) then
        call ADDNEQBC (SL,meshdatactx%IEN(1,n),S,IA,IAX,NEN,NEE,NDOF,MEQ)

    endif
enddo ! loop over elements

if (iecho.eq.8) then
    call ParallelLog("elastc", "Completed formk element loop")
endif

if (iecho.eq.1) then
    write(stdout,'(1X)')
endif

if (debug) then
    do i=1,MEQ
        write(*,80) i,(S(i,j),j=1,NDOF*NUMNP)
80      format(1X,I1,1X,8(1PG12.4))
    enddo
endif
 
return
end
!-----------------------------------------------------------------------

double precision function getaval (A,IDIAG,NEQ,I,J)
!
! Routine to get element I,J from matrix A.
!
use constants, only: zero
use iomodule          
!use modeldefinition, only: neq

implicit none
!-pass
double precision :: A
integer :: NEQ
integer          :: IDiag(neq)
dimension A(*)
integer :: i, j
!-locl
integer :: ii, jj, jh, ij
!
if (I.lt.1.or.I.gt.NEQ) then
    write(stderr,1) I,NEQ
 1        format(1x,'getaval: attempt to access element i=',I9, &
     ' (should be in range [1,',I9,']')
    call stoper()
endif

if (J.lt.1.or.J.gt.NEQ) then
    write(stderr,2) J,NEQ
 2        format(1x,'getaval: attempt to access element j=',I9, &
     ' (should be in range [1,',I9,']')
    call stoper()
endif
!
if (J.ge.I) then
    ii = I
    jj = J
else
    ii = J
    jj = I
endif

 if (ii.eq.1 .and. jj.eq.1) then
getaval = A(IDIAG(1))
 else
jh = IDIAG(jj) - IDIAG(jj-1)
if (ii.ge.jj-jh+1 .and. ii.le.jj) then
    ij = IDIAG(jj) - (jj-ii)
    getaval = A(ij)
else
    getaval = ZERO
endif
 endif
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine GLOAD (B,X,SKEW,GRAV,IEN,LM,LMX,MAT, &
   PRPMAT,T,TX,LMTX,PRPTEM,NSD,NDOF,NEN, igrav)

use constants,       only: one
USE MATERIALSMODULE, only : LMAT, &
                            MaterialIDOfElement
use modeldatamodule, only: modeldatactx
USE MODELDEFINITION, only: mode, numel, numrot, nee
!   USE MODELCTX
USE ALGEBRA,         only: clear
!   USE MESHDATAMODULE
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: HasSlipperyNodes, &
                           elementSlipperyID
#else
use modeldatamodule, only: HasSlipperyNodes
#endif
use iomodule,        only: stderr

!
! add body forces due to gravity.
! called once, handles all elements at once.

 implicit none
!-pass
 integer NSD,NDOF,NEN
double precision :: B, X, SKEW, GRAV, T, TX, PRPTEM, PRPMAT
integer          :: IEN, LM, LMX, LMTX, MAT
 dimension B(*),X(NSD,*),SKEW(2,*),GRAV(NDOF),IEN(NEN,*), &
   LM(NDOF,NEN,*),LMX(NDOF,NEN,*),PRPMAT(6,*),MAT(*), &
   T(*),TX(*),LMTX(NEN,*),PRPTEM(6,*)
!-locl
logical quad
double precision ::  pl(24),xl(24), dl(24)
double precision :: t0, rho
integer          :: n, nn, m
double precision :: densit, thick, texp
integer          :: ierr, igrav
double precision, external :: Tcentr ! from plnlib/f3dlib/etc.
!-init
 ierr = 0
!

do nn=1,NUMEL

    n = nn
    m = MaterialIDOfElement(n)

    densit = PRPMAT(5,m)
    thick  = PRPMAT(6,m)
    call LCOORD (X,xl,IEN(1,n),nn)
    call LCOORD (modeldatactx%D,dl,IEN(1,n),nn)


    if (MODE.ge.3) then
        texp   = PRPTEM(3,m)
        t0 = Tcentr (T,TX,xl,IEN(1,n),LMTX(1,n),ierr)
        if (ierr.ne.0) then
            write(stderr,1) n
 1           format(1x,'GLOAD: error in element ',I5)
            call stoper()
        endif
        rho = densit*(ONE-texp*t0)
    else
        rho = densit
    endif

    quad = (IEN(3,n).ne.IEN(4,n))
    pl = 0d0

    ! set pl (contribution of gravity on points)
    ! using gravity direction/magnitude in GRAV, and the element properties
    call GRAVLD (pl,(XL+dl),GRAV,THICK,RHO,quad,ierr)

    if (ierr.ne.0) then
        write(stderr,1) n
        call stoper()
    endif

    if (NUMROT.gt.0) then
        call RPFORC (pl,SKEW,IEN(1,n),NDOF,NEN)
    endif

    call ADDFOR (B,pl,LM(1,1,n),NEE,nn,4)   ! gravity

    if (HasSlipperyNodes()) then
!   CALL ADDFOR (B,pl,LMX(1,1,n),NEE,nn,5)  ! gravity
#ifdef EXPERIMENTAL_ReducedSlipAlloc
        if(elementSlipperyID(nn,25).gt.0) then
            call ADDFOR (B,PL,LMX(1,1,elementSlipperyID(nn,25)),NEE,nn,5)
        else
            ! this element has no slippery nodes. Skip it.
        endif
#else
        call ADDFOR (B,PL,LMX(1,1,nn),NEE,nn,5)
#endif

    endif

enddo
!
return
end

!-------------------------------------------------------------------------------
 subroutine INITD (IDIAG,A,NDOF,NUMNP,NEQ,iflag)

use constants,       only: zero,dfmin,dlogmin,dfmax,dlogmax
use modelctx,        only: getrank ! for debugging only
use modeldatamodule, only: modeldatactx
use debugmodule,     only: iecho
use timestepmodule,  only: refact

#ifdef SPARSE
use spetscmodule, only: M_A
#endif
use iomodule

!
! Routine for adjusting stiffness matrix for initial displacement bc
!    iflag = 0: modify A
!          = 1: modify B
!
 implicit none
!-pass
integer NDOF,NUMNP,NEQ
integer :: idiag
double precision :: A
dimension IDIAG(NEQ), A(*)
integer :: iflag
!-locl
integer :: do_initd
double precision :: diagmax
save do_initd,diagmax
integer :: indx
integer :: i, ii, ieqn, ierr, idof, n, j
double precision stfns
integer, external :: IJSKY ! from elsewhere in elastc.F
!-init


data do_initd/1/
REFACT = .false.
!
if (do_initd.eq.1) then
    do n=1,NUMNP
        do i=1,NDOF
            if (modeldatactx%IBOND(i,n).eq.4) goto 100
        enddo
    enddo
    do_initd = 0
    return
100        continue
endif

if (do_initd.eq.0) then
   return
endif


!
#ifdef SPARSE
 write(stderr,*) 'INITD: not prepared for sparse solution'
 call stoper()
#endif

!    Determine maximum diagonal value to use as scaling factor for stability
if (do_initd.eq.1) then
    diagmax = -dfmax
    do i=1,NEQ

#ifdef SPARSE
        call MatGetValues (M_A,1,i-1,1,i-1,stfns,ierr)
        if (ierr.ne.0) then
            write(stderr,1) ierr
 1                format(1x,'INITD: MatGetValues error ',I4)
            call stoper()
        endif
#else
        indx = IJSKY(IDIAG,NEQ,i,i)

        if (indx.eq.0) then
            write(stderr,1) i,i
 1               format(1x,'INITD: A(',I4,',',I4,') = 0')
            call stoper()
        endif

        stfns = A(indx)
#endif

        if (diagmax.lt.stfns) then
            diagmax = stfns
        endif

    enddo
    do_initd = 2
endif

!   Modify equations, retain matrix topology
 do n=1,NUMNP

!        write(*,*) '------- next point --------'

do idof=1,NDOF

if (modeldatactx%IBOND(idof,n).ne.4) then
    cycle
endif


REFACT = .true.
ieqn = modeldatactx%ID(idof,n)

!       replace b.c. by a r.h.s. load in each equation
do ii=1,NEQ
    ! fetch index for skyline storing method.
    indx = IJSKY (IDIAG,NEQ,ii,ieqn)
    stfns = ZERO

    if (indx.gt.0) then
        stfns = A(indx)
    endif


    if (iflag.eq.1) then
        write(*,*) 'rank',getrank(),'vtx',n,'pos',ii,'stfns BOND ', stfns, modeldatactx%BOND(idof,n)
        modeldatactx%B(ii)=modeldatactx%B(ii)-stfns*modeldatactx%BOND(idof,n)
    endif

enddo

!        Nullify column elements
if (ieqn.gt.1) then
    do i=1,ieqn-1
        indx = IJSKY (IDIAG,NEQ,i,ieqn)
        if (indx.gt.0.and.iflag.eq.0) then
            A(indx) = ZERO
        endif
    enddo
endif

!        Nullify row elements
if (ieqn.lt.NEQ) then
    do j=ieqn+1,NEQ
        indx = IJSKY (IDIAG,NEQ,ieqn,j)
        if (indx.gt.0.and.iflag.eq.0) then
            A(indx) = ZERO
        endif
    enddo
endif

indx = IJSKY (IDIAG,NEQ,ieqn,ieqn)
if (indx.le.0) then
    write(stderr,2)
 2            format(1x,'INITD fatal: indexing error')
    call stoper()
endif

if (iflag.eq.0) then
    A(indx) = diagmax
endif

if (iflag.eq.1) then
    modeldatactx%B(ieqn) = modeldatactx%BOND(idof,n)*diagmax
endif

    enddo
enddo
!
 if (iecho.eq.1 .and. REFACT ) then
write(stdout,3)
3    format(1x,'Initial displacement loads are being added')
 endif


 return
 end
!-------------------------------------------------------------------------------
integer function IJSKY (IDIAG,NEQ,ii,jj)

use iomodule          

! function returning index ij for skyline storage model
implicit none

!-pass
integer :: ii,jj,NEQ
integer :: IDIAG(NEQ)
!-local
integer :: i,j
integer :: jh, imin

if (ii.lt.1 .or. ii.gt.NEQ) then
    write(stderr,1) NEQ
 1        format(1x,'IJSKY fatal: 1 <= i <= ',I12)
    call stoper()
endif

if (jj.lt.1 .or. jj.gt.NEQ) then
    write(stderr,2) NEQ
 2        format(1x,'IJSKY fatal: 1 <= j <= ',I12)
    call stoper()
endif

if (ii.gt.jj) then
    i = jj
    j = ii
else
    i = ii
    j = jj
endif

if (j.eq.1) then
    IJSKY = 1
else
    jh = IDIAG(j)-IDIAG(j-1)
    imin = j - jh + 1
    if (i.lt.imin) then
        IJSKY = 0
    else
        IJSKY = IDIAG(j) - (j-i)
    endif
endif

return
end function
!-------------------------------------------------------------------------------
subroutine LDISBC (DL,DISPL,IEN,elemID,NDOF,NEN,globElemID)

! Program to localize prescribed displacement boundary conditions
! on one particular element.

! It is called both by elastic and by viscous,
! DISPL is either modeldatactx%D, or modeldatactx%DELD in the caller

use algebra, only: clear
use modelctx, only: getrank
use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use debugmodule
use iomodule

implicit none
!-pass
integer          :: NDOF,NEN, elemID
double precision :: DL(NDOF,NEN),DISPL(NDOF,meshdatactx%nvglobal)
integer          :: IEN(NEN)
integer          :: globElemID ! for debugging
!-locl
integer :: i,j,np,numeq

call CLEAR (DL,NDOF*NEN,"DL")

do j=1,NEN
    np = IEN(j)
    if (np.le.0) then
        write(stderr,*) "LDISBC **fatal**: IEN(",j,", partition local elem#=", &
         elemID," (global elem#=",globElemID,") )=",np
        call exitp(1)
    endif
    do i=1,NDOF
        ! if (debug) write(*,*) "Rank",getrank(),"LDISBC: LM(idof=",i,", ",j, &
        !            ", elm#=",elemID,")=",modeldatactx%LM(i,j,elemID)
        numeq = modeldatactx%LM(i,j,elemID)
        if (numeq.le.0) then
        !     if (debug) write(*,*) "Rank",getrank(),"LDISBC: set DL(idof=",i, &
        !                ",",j,")=",DISPL(i,np)
            DL(i,j) = DISPL(i,np)
        endif
    enddo
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine LDISP (DL,D,IEN,NDOF,NEN)

! Program to localize displacements

USE ALGEBRA, only: clear

implicit none
!-pass
integer          :: NEN,NDOF
integer          :: IEN(NEN)
double precision :: DL(NDOF,NEN),D(NDOF,*)

integer :: i,j,k

do i=1,NDOF
    do j=1,NEN
        DL(i,j) = 0.0d0
    enddo
enddo

do j=1,NEN
    k = IEN(j)
    ! ldisp is called with D or DELD as D, so must not
    ! use modeldata D, but a local variable
    DL(:,j) = D(:,k)
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine LSTRES (SL,STN,NSTR)
!
! routine to localize stress, strain or strain rate
!

implicit none
!-pass
integer NSTR, i
double precision :: STN(NSTR),SL(NSTR)

do i=1,NSTR
    SL(i) = STN(i)
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine LOAD (ID,IBOND,BOND,D,DELD,B,X,NUMNP, NEQ,IFLAG)

USE TIMESTEPMODULE
USE ALGEBRA
use debugmodule,     only: iecho, parallellog
use modelctx,        only: getrank ! for debugging
use meshdatamodule,  only: meshdatactx
use modeltopology,   only: NSD, NDOF
use modeldefinition, only: NEQGlobal

! program to transfer nodal boundary conditions into the appropriate vector:

!    B(NEQGlobal)      = global right-hand side vector
!    D(NDOF,NUMNP)     = global displacements
!    X(NSD,NUMNP)      = global coordinates
!    DELD(NDOF,NUMNP)  = changes in displacement (constant velocity)

!    IFLAG = 0    first entry (elastic solution)
!          > 0    viscous bcs


 implicit none
!-pass
 integer NUMNP,NEQ,IFLAG

integer          :: ID(NDOF,meshdatactx%nvglobal),IBOND(NDOF,meshdatactx%nvglobal)
double precision :: BOND(NDOF,meshdatactx%nvglobal), &
                    B(NEQGlobal), &
                    D(NDOF,meshdatactx%nvglobal), &
                    DELD(NDOF,meshdatactx%nvglobal), &
                    X(NSD,meshdatactx%nvglobal)



!-locl
integer FILE_outputf
double precision :: BONDL(3)
integer :: i, j, k, n

if (iecho.eq.8) then
    call ParallelLog("load", "entering loop")
endif



! This loop should ideally be over the local numbers only, but then
! the boundary conditions of points are not transmitted well across partition boundaries.
do j=1,meshdatactx%nvglobal
!do j=1,meshdatactx%nvlocal

do i=1,NDOF

#ifdef Lukas
        write(*,*) "rank",getrank(),'IBOND: ', i, j, IBOND(i,j)
#endif

if (IBOND(i,j).eq.1) then
!           displacement BC
    if (IFLAG.eq.0) then
        D(i,j) = BOND(i,j)
!                write(*,*) 'set displacement BC to: vertex: ', j, BOND(i,j)
    endif

else if (IBOND(i,j).eq.2) then
!           velocity BC
    DELD(i,j) = DELTP*BOND(i,j)
    if (IFLAG.eq.0) then
        D(i,j) = 0d0
    endif

else if (IBOND(i,j).eq.3) then
!           nodal force bc
!            write(*,*) 'NUMNP',j,'DOF',i,'IFLAG',IFLAG 
    if (IFLAG.eq.0.or.IFLAG.eq.3) then
        call CLEAR (BONDL,NDOF,"BONDL")
        BONDL(i) = BOND(i,j)
        do k=1,NDOF
            n = ID(k,j)
!                    write(*,*) 'NUMNP',j,'DOF',i,'DOF',k, 'n: ', n
            if (n.gt.0) then
                B(n) = B(n) + BONDL(k)
!                   write(*,*) "rank",getrank(),"LOAD added ", BONDL(k), "to B", i, B(i)
            endif
        enddo
    elseif (IFLAG.eq.5) then
        call CLEAR (BONDL,NDOF,"BONDL")
        BONDL(i) = BOND(i,j)
        do k=1,NDOF
            n = ID(k,j)
            if (n.gt.0) then
                B(n) = B(n) - BONDL(k)
!                   write(*,*) "rank",getrank(),"LOAD substracted ", BONDL(k), "from B", i, B(i)
            endif
        enddo
    endif

else if (IBOND(i,j).eq.5) then
!           strain rate BC (only works when x=0 is fixed boundary)

    if (IFLAG.eq.0) then
        D(i,j) = 0d0
    endif
    DELD(i,j) = DELTP*BOND(i,j)*(X(i,j)+D(i,j))
endif


 enddo ! loop over NDOF
  enddo ! loop over NUMNP

if (iecho.eq.8) then
    call ParallelLog("load", "leaving sub")
endif


 return
 end
!-------------------------------------------------------------------------------
subroutine LOADX (B,IDX,DIFORC,NUMSLP,NUMNPlocal)

use modelctx,      only: getrank
use modeltopology, only: NDOF
!
!    Transfer differential forces into global load vector
!
implicit none
!-pass
integer :: NUMSLP,NUMNPlocal
integer :: idx
double precision :: b, diforc
dimension B(*),IDX(NDOF,NUMNPlocal),DIFORC(NDOF,NUMNPlocal)
!-locl
 integer i,j,k

if (NUMSLP.le.0) then
    return
endif

! what happens if this partition has no slippery nodes? Then NUMNPlocal = 0, no loop
do j=1,NUMNPlocal
    do i=1,NDOF
        k = IDX(i,j)
        if (IDX(i,j).gt.0) then
            ! seems OK
            B(k) = DIFORC(i,j) ! set, not add
        endif
    enddo
enddo

!
return
end
!-------------------------------------------------------------------------------
subroutine LOADF (IEN,NFAULT,FAULT,DFAULT,NUMEL)

use modeldefinition, only: numfn
USE TIMESTEPMODULE
USE ALGEBRA
use modeltopology,   only: ndof, nen
use modeldatamodule, only: modeldatactx
use debugmodule,     only: xit
use timestepmodule,  only: nstep, ITWORK,IMWORK,IMTWRK,refact, faultvelo

use iomodule

!    set faulted node displacements DFAULT increment
!
!    logical VELBC is used to change interpretation of the data in FAULT
!    VELBC = .true.: FAULT means constant velocity (NFAULT(3,.) is not used)
!
!    VELBC = .false.: FAULT means differential displacement
!                     NFAULT(3,.) >= 0 differential displacement applied @ time step NFAULT(3,.)
!                     NFAULT(3,.)  < 0 periodic displacement applied every ABS(NFAULT(3,.) time step 
!
implicit none

!-pass
integer :: NFAULT
double precision :: FAULT, DFAULT
integer :: IEN(NEN,*)
integer :: NUMEL
dimension NFAULT(3,*),FAULT(NDOF,*),DFAULT(NDOF,NUMFN)

!-locl
integer j,k,mode
logical doit,VELBC,NoMatch
double precision :: weight
integer i,nodeID,elemID,vertexID

if (NUMFN.le.0) then
    return
endif

VELBC = (faultvelo.ne.0)
call CLEAR(modeldatactx%DFAULT,NDOF*NUMFN,"modeldatactx%DFAULT")

do k=1,NUMFN   ! NUMFN is partition-local number of slippery elements and nodes (AOmodule unclear on ghost nodes and ghost elements)
    if (VELBC) then
        doit = .true.
        weight = DELTP
    else
        if (NFAULT(3,k).ge.0) then
            doit = (NSTEP.eq.NFAULT(3,k))
        else
            doit = (MOD(NSTEP,-NFAULT(3,k)).eq.0)
        endif
        weight = 1d0
    endif

    if (doit) then
        nodeID = NFAULT(2,k) ! partition-local node number
        elemID = NFAULT(1,k) ! partition-local element number
        NoMatch = .true.
        do j=1,NEN
            if (IEN(j,elemID).eq.nodeID) then
                vertexID = j
                NoMatch = .false.
            endif
        enddo 

        if (NoMatch) then
            write(stderr,*) 'LOADF: node ',nodeID,' doesn''t occur in element ',elemID
            call xit(1,"LOADF")
        endif
        i = modeldatactx%LMF(1,vertexID,elemID) ! returns the DFAULT index for this node and element
        if (i.lt.1) then
            write(stderr,*) 'LOADF: LMF index < 1 of node ',nodeID,' and element ',elemID
            call xit(1,"LOADF")
        endif
        do j=1,NDOF
            modeldatactx%DFAULT(j,i) = modeldatactx%DFAULT(j,i) + weight*FAULT(j,k)
        enddo
    endif

enddo

return
end
!-------------------------------------------------------------------------------
 subroutine PIVOTS (A,IDIAG,ID,IDX,NDOF,NEQ,NUMSN, &
   NUMNP,MODE)

 use constants, only: dfmin
 use modeldatamodule, only: HasSlipperyNodes
use iomodule 
 
!
!    Routine to determine the number of zero and negative terms in the
!    array D of factorization of A = U(TRANSPOSE) * D * U
!

 implicit none
!-pass
 integer MODE,NDOF,NEQ,NUMSN,NUMNP
integer :: idiag, id, idx
double precision :: a
 dimension A(*),IDIAG(NEQ),ID(NDOF,NUMNP),IDX(NDOF,NUMNP)

!-local
integer :: iz, in, i, k
integer :: idof
integer :: node

!-init
 iz = 0
 in = 0
!
 do i=1,NEQ
k = IDIAG(i)
if (ABS(A(k)).le.dfmin .or. A(k).lt.-dfmin) then
    call EQNODE (ID,NDOF,NUMNP,i,node,idof)

    if (node.le.0 .and. HasSlipperyNodes()) then
        call EQNODE (IDX,NDOF, NUMNP,i,node,idof)
    endif

    if (node.le.0) then
        write(stderr,1) i
 1            format(1x,'PIVOTS: cannot identify node',1x, &
             'for equation number ',I5)
        call stoper()
    endif

    if (ABS(A(k)).le.dfmin) then
        write(stderr,2) i,node,idof
 2               format(1x,'zero pivot: NEQ=',I5,4x,'NODE=',I5,4x, &
             'IDOF=',I5)
        write(stderr,*) 'k=',k
        iz=iz+1
    else
        write(stderr,3) i,node,idof
 3            format(1x,'negative pivot: NEQ=',I5,4x,'NODE=',I5, &
             4x,'IDOF=',I5)
        in=in+1
        write(stderr,*) 'k=',k
    endif

endif
 enddo
!
 if (iz.ne.0 .or. in.ne.0) then
write(stderr,4) iz,in
 4        format(///1x,'*** ILL-CONDITIONED STIFFNESS MATRIX ***'/ &
     1x,I5,' ZERO PIVOTS FOUND'/ &
     1x,I5,' NEGATIVE PIVOTS FOUND'//)
if (MODE.ne.1 .and. MODE.ne.4) call stoper()
 endif
!
 return
 end
!-------------------------------------------------------------------------------
subroutine EQNODE (ID,NDOF,NUMNP,IEQ,NODE,IDOF)

! Finds the index NDX and dof number IDOF belonging to equation number IEQ
! When called with the ID(NDOF,NUM=NUMNP) array, NDX is the node number
! When called with the IDX(NDOF,NUM=NUMSN) array, NDX refers to the entry number in the IDSLP array (ie, IDSLP(NDX) is the node number


implicit none
!-pass
integer :: NDOF,NUMNP,IEQ,NODE,IDOF
integer :: ID(NDOF,NUMNP)
!-local
integer :: i,j
!-init
NODE=-1
IDOF=-1

do j=1,NUMNP
    do i=1,NDOF
        if (ID(i,j).eq.IEQ) then
            NODE = j
            IDOF = i
            return
        endif
    enddo
enddo

end subroutine

!-------------------------------------------------------------------------------
 subroutine PRINTA (A,IDIAG,NEQ)

use iomodule          
!
! Routine to print column compacted matrix A
!     A(NA)        matrix in column compacted form
!    IDIAG(NEQ)    diagonal adresses of A(NA)
!    NEQ        number of equations
!

 implicit none
!-pass
integer :: NEQ, idiag
double precision :: a
dimension A(*),IDIAG(NEQ)
!-locl
 character(len=20) f
integer :: i, j
double precision, external :: getaval
!-init
 j=0
!
 if (.not.OUTFIL(1)) return
!
 write(f,1) NEQ
 1   format('(1x,',I3,'(1PG12.4,1X))')
 do i=1,NEQ
if (MOD(i-1,50).eq.0) write(luout(1),2) NEQ,NEQ
 2        format(///' M A T R I X : ',I3,' X ',I3/)
write(luout(1),fmt=f) &
       (getaval(A,IDIAG,NEQ,i,j),j=1,NEQ)
 enddo
!
 return
 end subroutine
!-------------------------------------------------------------------------------
subroutine PRINTB (B,ID,IDX,NDOF)

USE MODELDEFINITION
use constants, only: dfmin,dlogmin,dfmax,dlogmax
use iomodule

! Prints load vector

implicit none
!-pass
integer NDOF
integer :: id, idx
double precision :: B
dimension B(*),ID(NDOF,*),IDX(NDOF,*)
!-local
integer :: nn, i, j
integer :: ieq

nn = 51

do j=1,NUMNP
    do i=1,NDOF
        ieq = ID(i,j)
        if (ieq.le.0) goto 100
        if (nn.gt.50) then
            write(luout(1),5)
5           format(//1x,'NODE',7x,'DOF',7x,'NEQ',9x,'B(NEQ)'/)
            nn = 0
        endif
        write(luout(1),6) j,i,IEQ,B(IEQ)
6       format(1x,i8,5x,i5,5x,i8,5x,1PG20.8)
        nn = nn + 1
100     if (NUMSN.le.0) cycle
        ieq = IDX(i,j)
        if (ieq.le.0) cycle
        if (nn.gt.50) then
             write(luout(1),5)
             nn = 0
        endif
        write(luout(1),6) j,i,IEQ,B(IEQ)
        nn = nn + 1
    enddo
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine PRINTD(D,NUMNP,IFLAG)

USE LABELMODULE
USE MODELCTX,        only: getrank
USE MESHDATAMODULE,  only: meshdatactx
USE MODELDATAMODULE, only: modeldatactx,  &
                           equationreference
USE MODELDEFINITION, only: NEQ
use debugmodule,     only: debug, iecho, xit
use formatsmodule,   only: FMT_fedsk2D_displacement, &
                           FMT_fedsk3D_displacement
use modeltopology,   only: NDOF
use iomodule


!    Routine to print:
!     (differential) displacements   (IFLAG=0) or 
!     (differential) velocities      (IFLAG=1) or
!     dfferential strengths          (IFLAG=2)

!    Note that this functions is called for both the displacement
!    modeldatactx%D and the differential displacement (modeldatactx%DELD)
!    hence we must use D here as the input variable, in stead of 
!    a member of the modeldatactx struct.

implicit none
!-pass
integer NUMNP
double precision :: D(NDOF,NUMNP)
integer :: iflag
!-locl
integer :: vertexid
integer :: i, n, nn, m
!-init

nn = 50
!
if (iecho.eq.3) then
    write(stdout,1)
 1   format(1x,'if (NUMNP <> 0) {'/ &
 5x,'for (n=1;n<=NUMNP;n++) {'/ &
 9x,'write n,(D(i,n),i=1,NDOF) (I5,X,6E15.6)'/ &
 5x,'}'/1x,'}')
endif
!
if (NUMNP.le.0) then
    return
endif


!***********************************************
!  Adding displacement to fedsk files
!***********************************************

if (OUTFIL(2)) then
    if (ascout) then

        do n=1,NUMNP
            if (NDOF.eq.2) then
                write(luout(2),FMT_fedsk2D_displacement) n, D(1,n), D(2,n)
            else if (NDOF.eq.3) then
                write(luout(2),FMT_fedsk3D_displacement) n, D(1,n), D(2,n), D(3,n)
            else
                write(*,*) 'NDOF unknown: ', NDOF, 'contact model support'
            endif

        enddo

    else
        write(luout(2),err=200) ((D(I,N),I=1,NDOF),N=1,NUMNP)
    endif
endif
!
if (.not.OUTFIL(1)) then
    return
endif
!
do n=1,NUMNP
    m = n
    nn = nn + 1

    if (nn.gt.50) then
        nn = 1

        if (IFLAG.eq.0) then
            write(luout(1),3) (LABELD(i),i=1,NDOF)
        else if (IFLAG.eq.1) then
            write(luout(1),7) (LABELD(i),i=1,NDOF)
        else if (IFLAG.eq.2) then
            write(luout(1),8) (LABELD(i),i=1,NDOF)
        else
            write(*,*) "Printd called with unknown flag: ", iflag
            write(*,*) "This is not supposed to happen. Please contact model support."
        endif
 3           format(///1X,' D I S P L A C E M E N T S '/// &
         4X,'NODE NUMBER',15X,6(A4,14X))
 7           format(///1X,' V E L O C I T I E S '/// &
         4X,'NODE NUMBER',15X,6(A4,14X))
 8           format(///1X,' D I F F E R E N T I A L     F O R C E S '/// &
         4X,'NODE NUMBER',15X,6(A4,14X))
        write(luout(1),4)
 4           format(/)
    endif


    if (.not.debug) then
        write(luout(1),5)  m,(D(i,n),i=1,NDOF)
 5            format(3X,I8,10X,6(1PE18.5))
    else
        write(luout(1),6)  m,(D(i,n),i=1,NDOF)
 6            format(3X,I8,10X,6(1PG20.8))
    endif

enddo
!

return

200    call prerr('printd')
call xit(1," ")
end subroutine
!-------------------------------------------------------------------------------
subroutine PRINTDXE (DXE,IEN,IDSLE,IFLAG)

USE LABELMODULE

#ifdef EXPERIMENTAL_ReducedSlipAlloc
USE MODELDEFINITION,only: nElemsWithSlipperyNodes, &
                          NUMSEglobal
#else
USE MODELDEFINITION,only: NUMSEglobal
#endif

use formatsmodule,  only: FMT_fedsk2D_diff_displacement, &
                          FMT_fedsk3D_diff_displacement
use meshdatamodule, only: meshdatactx
use debugmodule,    only: iecho, xit

use modeldatamodule,only: modeldatactx
use modeltopology,  only: ndof, nen 
use modelctx,       only: getsize
use iomodule

!
! Prints element-based differential displacements across free interfaces
! Displacements in all elements are printed as this is the only
! practical way of dealing with regridding (where NUMSE is modified).
!

! This sub is caled for for differential displacement and for
! differential velocities. Hence DXE cannot be replaced by the DXE
! from the modeldata structure

implicit none
!-pass
integer :: iflag
integer :: ien, idsle
double precision :: dxe
#ifdef EXPERIMENTAL_ReducedSlipAlloc
dimension :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
dimension :: DXE(NDOF,NEN,*)
#endif

dimension IEN(NEN,*),IDSLE(*)
!-locl
integer :: n, nn, i, j, nout
character(len=23) fmt

integer :: globalElemID
integer :: slipperySequenceNr

!
if (iecho.eq.3) then
write(stdout,1)
 1   format(1x,'differential displacements across', &
 ' free interfaces'/ &
 1x,'for (n=1;n<=NUMEL;n++) {'/ &
 5x,'write n,j,(DXE(i,j,n),i=1,NDOF,j=1,NEN)',1x, &
 '(I5,X,3E15.6)'/1x,'}')
endif
!

if (NUMSEglobal.le.0) then
!    write(*,*) "NUMSE 0"
    ! NUMSE ? slippery elements?
    return
endif

if (.not.OUTFIL(2)) then
    goto 100
endif

if (NDOF.eq.2) then
    fmt = FMT_fedsk2D_diff_displacement
else if (NDOF.eq.3) then
    fmt = FMT_fedsk3D_diff_displacement
else
    write(*,*) 'NDOF not 2 or 3; this should not happen.'
endif


!write(*,*) "DXE: ", DXE

#ifdef EXPERIMENTAL_ReducedSlipAlloc
! test this. IEN probably not correct
if (ascout) then
    do nn=1,meshdatactx%nelocal
        globalElemID = meshdatactx%gloelt(nn)
        
        ! for each element...
        slipperySequenceNr = modeldatactx%SlipElemLookup(globalElemID)

        if (slipperySequenceNr.gt.0) then
            do j=1,NEN
                write(luout(2),fmt,err=200) globalElemID,IEN(j,nn), DXE(1:NDOF,j,slipperySequenceNr)
            enddo
        else
            ! This element has no slippery node; nothing to write
        endif

    enddo
else  
    write(luout(2),err=200) (((DXE(i,j,n),i=1,NDOF),j=1,NEN), n=1,nElemsWithSlipperyNodes)
endif
#else
if (ascout) then
    do nn=1,meshdatactx%nelocal
        ! for each element...
        n = nn
        do j=1,NEN
!        do j=1,3      
            ! ... and for each vertex of that element, write diff. displ.
            write(luout(2),fmt,err=200) n,IEN(j,n), DXE(1:NDOF,j,n)
        enddo
!        write(luout(2),fmt,err=200) n,IEN(4,n), 0, 0

    enddo
else
    write(luout(2),err=200) (((DXE(i,j,n),i=1,NDOF),j=1,NEN), n=1,meshdatactx%nelocal)
endif
#endif

100 if (.not.OUTFIL(1)) then
    return
endif


if (getsize().eq.1) then

nout = 51
do nn=1,NUMSEglobal

    if (nn.gt.0 .and. nn.le.size(modeldatactx%IDSLE,1)) then
        n = modeldatactx%IDSLE(nn)

        if (n.gt.0 .and. n.le.size(meshdatactx%IEN,2)) then
        

            if (nout.gt.50) then

                if (IFLAG.eq.0) then
                   write(luout(1),3,err=200) (LABELD(i),i=1,NDOF)
                else
                    write(luout(1),7,err=200) (LABELD(i),i=1,NDOF)
                endif

 3       format(///1X, &
             'D I F F E R E N T I A L',3x, &
             'D I S P L A C E M E N T S'///1X,'ELEMENT',5X, &
             'NODE',6X,6(A4,11X))
 7       format(///1X, &
             'D I F F E R E N T I A L',3x, &
             'V E L O C I T I E S'///1X,'ELEMENT',5X, &
             'NODE',6X,6(A4,11X))
            write(luout(1),4,err=200)
 4       format(/)
                nout = 0
            endif

            do j=1,NEN
                nout = nout + 1
#ifdef EXPERIMENTAL_ReducedSlipAlloc
!        slipperySequenceNr = modeldatactx%SlipElemLookup(globalElemID)
!        write(luout(1),5,err=200) n,IEN(j,n), &
!         (DXE(i,j,slipperySequenceNr),i=1,NDOF)
                write(luout(1),*) "DXE in FEOUT not yet supported in reduced memory mode."
#else
                write(luout(1),5,err=200) n,meshdatactx%IEN(j,n), &
                                          (DXE(i,j,n),i=1,NDOF)
#endif
 5              format(I8,1X,I8,3(1PG15.6))
            enddo



        else
            write(*,*) "PRINTDXE has found illegal n", n
            write(*,*) "wants to read IEN entry ", n
            write(*,*) "while there are only ", size(meshdatactx%IEN,2)
        endif


    else
        write(*,*) "PRINTDXE has found illegal nn", nn
        write(*,*) "Wants to check IDSLE (slippery element) entry ", nn
        write(*,*) "while there are only", size(modeldatactx%IDSLE,1), "entries"
    endif


enddo

endif

return

200    call prerr('printdxe')
call xit(1," ")

end

!-------------------------------------------------------------------------------
subroutine PRINTP(A,IDIAG,ID,IDX,NDOF,IMODE)

USE MODELDEFINITION
use iomodule, only: stdout, luout, outfil

!
! Program to print array D after factorzation A=(U)T*D*U
!
!    IMODE    ACTION
!      0    requested output of displacements
!      1    requested output of temperatures


implicit none
!-pass
integer NDOF,IMODE
double precision :: a
double precision :: d
integer :: idiag(neq), id, idx

dimension A(*),ID(NDOF,*),IDX(NDOF,*)

#ifndef SPARSE

!-locl
integer nn,nq,nsn,n,j,i,idof,node
!-init
nn = 0


if (IMODE.eq.0) then
    nq  = NEQ
    nsn = NUMSN
else if (IMODE.eq.1) then
    nq  = NTEQ
    nsn = NUMSN+NUMFN
endif
!
do n=1,nq
    do j=1,NUMNP
        node = j
        do i=1,NDOF
            idof = i
            if (ID(i,j).eq.n) then
                goto 100
            endif
            if (nsn.gt.0) then
                if (IDX(i,j).eq.n) then
                    goto 100
                endif
            endif
        enddo
    enddo
enddo
!
100        if (MOD(nn,50).eq.0) then

if (OUTFIL(1)) then
    write(luout(1),1)
else
    write(stdout,1)
endif
 1           format(///1x,'ARRAY',3x,'D(NEQ)',3x,'OF',3x, &
         'FACTORIZATION',3x,'A=(U)T*D*U'/// &
         4x,'NEQ',1x,'NODE',1x,'DOF',13x,'D'//)
endif

nn = nn + 1
i  = IDIAG(n)
d  = A(i)

if (OUTFIL(1)) then
    write(luout(1),2) n,node,idof,d
else
    write(stdout,2) n,node,idof,d
endif

 2        format(1X,I8,I8,I5,4X,1PE20.8)
200    continue
!
#endif
 return
 end

!-------------------------------------------------------------------------------

subroutine PRINTS (STN,NSTR,NUMEL)

USE LABELMODULE
use debugmodule,     only: iecho, parallelLog, xit
use formatsmodule,   only: FMT_fedsk_stressStrainRate_nstr3, &
                           FMT_fedsk_stressStrainRate_nstr4, &
                           FMT_fedsk_stressStrainRate_nstr5, &
                           FMT_fedsk_stressStrainRate_nstr6
use iomodule

!     Print stress


 implicit none
!-pass
 integer NSTR,NUMEL
double precision :: stn
 dimension STN(NSTR,*)

! for outfil

!-locl
! character(len=19) fmt
integer :: i, n, nn
!-init
 nn = 51
!
if (iecho.eq.3) then
    write(stdout,1)
 1   format(1x,'if (NUMEL <> 0) {'/ &
 5x,'for (n=1;n<=NUMEL;N++) {'/ &
 9x,'write n,(STN(i,n),i=1,NSTR) (I5,X,6E15.6)'/ &
 5x,'}'/1x,'}')
endif

if (iecho.eq.8) then
    call ParallelLog("PRINTS", "Hello")
endif
!
if (NUMEL.eq.0) then
    write(*,*) "Error: PRINTS is called to write zero elements. Returning"
    return
endif

!***********************************************************
! Write Stress to fedsk file
!***********************************************************

if (OUTFIL(2)) then

    if (iecho.eq.8) then
        call ParallelLog("PRINTS", "Writing stress to fedsk")
    endif

    if (ascout) then
        do n=1,NUMEL

!            write(*,*) "writin stress", n,(STN(i,n),i=1,NSTR)

            if (NSTR.eq.6) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr6) n,(STN(i,n),i=1,NSTR)
            else if (NSTR.eq.5) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr5) n,(STN(i,n),i=1,NSTR)
            else if (NSTR.eq.4) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr4) n,(STN(i,n),i=1,NSTR)
            else if (NSTR.eq.3) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr3) n,(STN(i,n),i=1,NSTR)
            endif

        enddo
    else
        write(luout(2),err=200) ((STN(I,N),I=1,NSTR),N=1,NUMEL)
    endif
endif

if (.not.OUTFIL(1)) then
    if (iecho.eq.8) then
        call ParallelLog("PRINTS", "Not writing FEOUT, leaving sub")
    endif
    return
endif

do n=1,NUMEL
    nn = nn + 1
    if (nn.gt.50) then
        nn = 1
        write(luout(1),3) (LABELS(I),I=1,NSTR)
3           format(///2X,'S T R E S S'///3X,' ELEM #',5X,6(A4,10X),A4)
        write(luout(1),4)
 4           format(/)
    endif
    write(luout(1),5) n,(STN(i,n),i=1,NSTR)
 5       format(1X,I8,1X,6(1PE14.5))
enddo

return

200 if (iecho.eq.8) then
    call ParallelLog("PRINTS", "Failed writing stress")
endif

call prerr('prints')
write(*,*) "An error occured while writing the stress to file. Leaving GTecton."
call xit(1," ")
end subroutine

!-------------------------------------------------------------------------------

#ifdef EXPERIMENTALTRANSIENTDISCREEP

subroutine PRINTTAYSTR (STNTAY,NSTRTAY,NUMEL)

USE LABELMODULE
use debugmodule,     only: iecho, parallelLog, xit
use formatsmodule,   only: FMT_fedsk_stressStrainRate_nstr1, &
                           FMT_fedsk_stressStrainRate_nstr3, &
                           FMT_fedsk_stressStrainRate_nstr4, &
                           FMT_fedsk_stressStrainRate_nstr5, &
                           FMT_fedsk_stressStrainRate_nstr6
use iomodule

!     Print Taylor stress


 implicit none
!-pass
 integer NSTRTAY,NUMEL
double precision :: stntay
 dimension STNTAY(NSTRTAY,*)

! for outfil

!-locl
! character(len=19) fmt
integer :: i, n, nn
!-init
 nn = 51
!
if (iecho.eq.3) then
    write(stdout,1)
 1   format(1x,'if (NUMEL <> 0) {'/ &
 5x,'for (n=1;n<=NUMEL;N++) {'/ &
 9x,'write n,(STNTAY(i,n),i=1,NSTRTAY) (I5,X,6E15.6)'/ &
 5x,'}'/1x,'}')
endif

if (iecho.eq.8) then
    call ParallelLog("PRINTS", "Hello")
endif
!
if (NUMEL.eq.0) then
    write(*,*) "Error: PRINTS is called to write zero elements. Returning"
    return
endif

!***********************************************************
! Write Taylor Stress to fedsk file
!***********************************************************

if (OUTFIL(2)) then

    if (iecho.eq.8) then
        call ParallelLog("PRINTTAYSTR", "Writing Taylor stress to fedsk")
    endif

    if (ascout) then
        do n=1,NUMEL

            if (NSTRTAY.eq.6) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr6) n,(STNTAY(i,n),i=1,NSTRTAY)
            else if (NSTRTAY.eq.5) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr5) n,(STNTAY(i,n),i=1,NSTRTAY)
            else if (NSTRTAY.eq.4) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr4) n,(STNTAY(i,n),i=1,NSTRTAY)
            else if (NSTRTAY.eq.3) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr3) n,(STNTAY(i,n),i=1,NSTRTAY)
            else if (NSTRTAY.eq.1) then
                write(luout(2),FMT_fedsk_stressStrainRate_nstr1) n,(STNTAY(i,n),i=1,NSTRTAY)
            endif

        enddo
    else
        write(luout(2),err=200) ((STNTAY(I,N),I=1,NSTRTAY),N=1,NUMEL)
    endif
endif

if (.not.OUTFIL(1)) then
    if (iecho.eq.8) then
        call ParallelLog("PRINTTAYSTR", "Not writing FEOUT, leaving sub")
    endif
    return
endif

do n=1,NUMEL
    nn = nn + 1
    if (nn.gt.50) then
        nn = 1
        write(luout(1),3) (LABELS(I),I=1,NSTRTAY)
3           format(///2X,' T A Y L O R S T R E S S'///3X,' ELEM #',5X,6(A4,10X),A4)
        write(luout(1),4)
 4           format(/)
    endif
    write(luout(1),5) n,(STNTAY(i,n),i=1,NSTRTAY)
 5       format(1X,I8,1X,6(1PE14.5))
enddo

return

200 if (iecho.eq.8) then
    call ParallelLog("PRINTTAYSTR", "Failed writing Taylor stress")
endif

call prerr('printtaystr')
write(*,*) "An error occured while writing the Taylor stress to file. Leaving GTecton."
call xit(1," ")
end subroutine

#endif

!-------------------------------------------------------------------------------

subroutine PSPLIT (TFAULT,NFAULT,NDOF,NUMFN)

USE LABELMODULE
use iomodule          

!    Prints out faulted node displacements

implicit none
!-pass
integer :: nfault, ndof, numfn
double precision :: tfault
dimension TFAULT(NDOF,*),NFAULT(3,*)

!-local
integer :: nn, i, j

!-init
nn = 0

if (NUMFN.eq.0) then
    return
endif

if (OUTFIL(1)) then
    do j=1,NUMFN

        if (MOD(nn,50).eq.0) then
            write(luout(1),3) (LABELD(i),i=1,NDOF)
 3           format(///1x,'S P L I T',3x,'N O D E',3x, &
         'D I S P L A C E M E N T S'//4X,'ELEMENT',7X,'NODE', &
         15X,3(A4,14X))
        endif

        if (MOD(nn,50).eq.0) then
            write(luout(1),'(/)')
            nn = nn + 1
            write(luout(1),4) IABS(NFAULT(1,j)),NFAULT(2,j), (TFAULT(i,j),i=1,NDOF)
 4           format(1X,I8,4X,I8,6X,6(1PG18.5))
        endif
    enddo
endif
!
return
!
end subroutine
!-------------------------------------------------------------------------------
subroutine PRINTX (NDOF,NUMSN,IFLAG)

USE LABELMODULE
use debugmodule,     only: iecho, xit
use formatsmodule,   only: FMT_fedsk2D_diff_displacement, &
                           FMT_fedsk3D_diff_displacement

use modeldatamodule, only: modeldatactx
use iomodule
! Prints out differential displacements (IFLAG=0) or velocities (IFLAG=1)
! across free interfaces


implicit none
!-pass
integer          :: NDOF,NUMSN
!double precision :: DX
integer :: iflag
!-locl
character(len=23) fmt
integer :: i, n, m
!

if (iecho.eq.3) then
    write(stdout,1)
 1   format(1x,'differential displacements across free interfaces'/ &
            1x,'for (n=1;n<=NUMSN;n++) {'/ &
            5x,'write n,(DX(i,n),i=1,NDOF) (I5,X,3E15.6)'/ &
            1x,'}')
endif

if (NUMSN.le.0) then
    return
endif

if (OUTFIL(2)) then

    if (ascout) then

        if (NDOF.eq.2) then
            fmt = FMT_fedsk2D_diff_displacement
        else
            fmt = FMT_fedsk3D_diff_displacement
        endif

!#ifdef Lukas
!        fmt = '(I12,1X,3(1PG25.16))'
!#else
!        fmt = '(I12,1X,3(1PG15.6))'
!#endif

        write(*,*) "modeldat DX:", modeldatactx%DX
        write(*,*) "IDslip:", modeldatactx%IDSLP


        do n=1,NUMSN
            m=modeldatactx%IDSLP(n)
            write(luout(2),fmt) m,(modeldatactx%DX(i,m),i=1,NDOF)
        enddo
    else
        write(luout(2),err=200) &
           ((modeldatactx%DX(i,modeldatactx%IDSLP(n)),i=1,NDOF), n=1,NUMSN)
    endif
endif

if (.not.OUTFIL(1)) then
    return
endif

if (IFLAG.eq.0) then
    write(luout(1),3) (LABELD(i),I=1,NDOF)
else
    write(luout(1),7) (LABELD(i),I=1,NDOF)
endif

 3   format(///1X, &
 'D I F F E R E N T I A L   D I S P L A C E M E N T S'/// &
 4X,'NODE NUMBER',15X,6(A4,14X))

 7   format(///1X, &
 'D I F F E R E N T I A L   V E L O C I T I E S'/// &
 4X,'NODE NUMBER',15X,6(A4,14X))
write(luout(1),4)

 4   format(/)


! this is a bit iffy. idslp array is only filled for local entries. The rest is 0
do n=1,NUMSN
!    write(*,*) "a", n
    m = modeldatactx%IDSLP(N)
    if (m.gt.0) then
!        write(*,*) "idslp", modeldatactx%IDSLP
!       write(*,*) "modeldat DX:", modeldatactx%DX
!        write(*,*) "modeldat DX:", size(modeldatactx%DX,1), size(modeldatactx%DX,2)

        write(luout(1),5) m,(modeldatactx%DX(i,m),i=1,NDOF)
 5        format(7x,I8,3(3x,1PG15.6))
    endif
enddo


return
200    call prerr('printx')
 call xit(1," ")
 end
!-------------------------------------------------------------------------------
subroutine RDISP (D,SKEW,NDOF,NUMNP)
!
! Rotates displacement vector to global coordinate directions
!
USE ALGEBRA,         only: formrt, equate, vmprd
use constants, only: eps
use modelctx, only: getrank


implicit none
!-pass
integer NDOF,NUMNP
double precision :: d, skew
dimension D(NDOF,*),SKEW(2,*)
!-locl
double precision :: ROT(3,3),DL(3)
integer          :: n

do n=1,NUMNP
    if (ABS(SKEW(1,n)).gt.EPS .or. &
        ABS(SKEW(2,n)).gt.EPS) then
        call FORMRT (SKEW(:,n),ROT,NDOF)
        call EQUATE (DL,D(1:NDOF,n),NDOF)
        call VMPRD (ROT,DL,D(1,n),NDOF,3)
    else
        ! the node is not rotated.
    endif
enddo
!
return
end
!-------------------------------------------------------------------------------
 subroutine ROTFLT (NFAULT,DFAULT,SKEW,NDOF,NUMFN)
!
!    Rotate faulted node displacements to global coordinate directions
!
USE ALGEBRA,         only: formrt, equate, vmprd
 use constants, only: eps


 implicit none
!-pass
 integer NDOF,NUMFN
integer :: nfault
double precision :: dfault, skew
 dimension DFAULT(NDOF,*),NFAULT(3,*),SKEW(2,*)
!-locl
 integer i,n
double precision rot(3,3),dl(3)
!
 if (NUMFN.le.0) then
    return
endif
!
 do i=1,NUMFN
n = NFAULT(2,i)
if (ABS(SKEW(1,n)).gt.EPS .or.  &
    ABS(SKEW(2,n)).gt.EPS) then
    call FORMRT(SKEW(:,n),rot,NDOF)
    call EQUATE(dl,DFAULT(1:NDOF,i),NDOF)
    call VMPRD (rot,dl,DFAULT(1,i),NDOF,3)
endif
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine RPFORC (P,SKEW,IEN,NDOF,NEN)
!
! Rotates element force vector into skewed coordinate system
!
 USE ALGEBRA
 use constants, only: eps


 implicit none
!-pass
 integer NDOF,NEN
integer :: ien
double precision :: p, skew
 dimension P(NDOF,*),SKEW(2,*),IEN(*)
!-locl
double precision :: rot(3,3),pl(3)
integer          :: i, k
!
 do i=1,NEN
k = IEN(i)
if (ABS(SKEW(1,k)).gt.EPS .or. ABS(SKEW(2,k)).gt.EPS) then
    call FORMRT (SKEW(:,k),rot,NDOF)
    call TRANSP (rot,3)
    call EQUATE (pl,P(1:NDOF,i),NDOF)
    call VMPRD  (rot,pl,P(1,i),NDOF,3)
endif
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine RSTIFF(S,STEMP,IEN,SKEW,NDOF,NEN,NEE,localElemID, globalElemID)

! Rotates the local stiffness matrix-s coordinate system with respect to the
! global coordinates at specified nodes for imposition of skew boundary
! conditions.  Kt= (R-1) K R.

USE ALGEBRA,         only: formrt, equate, vmprd, transp
use meshdatamodule, only: meshdatactx
use modelctx, only: getrank ! only for debugging
use debugmodule, only: debug, debugFileID
use constants, only: zero, pi, halfpi

implicit none
!-pass
integer NDOF,NEE,NEN
integer :: localElemID
integer :: globalElemID
integer :: IEN(NEN,meshdatactx%nelocal)
double precision :: S(NEE,NEE),STEMP(NEE,NEE)

double precision :: SKEW
dimension :: SKEW(2,*)
!-locl
double precision :: ROT(3,3)
 integer i,j,k,l,ll
logical :: quad
double precision :: sum
!-init

do i=1,NEN
    ! loop over the vertices of this element
    ! get vertex id in partition-local numering
    k = IEN(i,localElemID)
    if (SKEW(1,k).ne.zero .or. SKEW(2,k).ne.zero) then
        call FORMRT(SKEW(:,k),ROT,NDOF)
        call EQUATE(STEMP,S,NEE*NEE)
        ll = NDOF*(i-1)
        do j=1,NEE
            do k=1,NDOF
                sum = zero
                do L=1,NDOF
                    sum = sum + STEMP(j,ll+L)*ROT(L,k)
                enddo
                S(j,ll+k) = sum
            enddo
        enddo
        call TRANSP(ROT,3)
        call EQUATE(STEMP,S,NEE*NEE)
        do j=1,NDOF
            do k=1,NEE
                sum = zero
                do L=1,NDOF
                    sum = sum + ROT(j,L)*STEMP(ll+L,k)
                enddo
                S(ll+j,k) = sum
            enddo
        enddo
    endif
enddo

return
end
!-------------------------------------------------------------------------------
 subroutine SOLVE (A,B,X,IDIAG,ID,IDX,NUMNP,NEQ, &
   NUMSN,NDOF,FACT,BACK)

USE MODELCTX
USE ALGEBRA,         only: EQUATE
USE MESHDATAMODULE,  only: meshdatactx
USE MODELDATAMODULE, only: modeldatactx, &
                           equationreference
USE MODELDEFINITION, only: NEQglobal, mode
use debugmodule,     only: parallelLog, &
                           iecho
use iomodule
!    USE AOmodule,        only: equations

implicit none
!-pass
integer NUMNP,NEQ,NUMSN,NDOF
logical FACT,BACK
double precision :: a, b, x
integer          :: id, idx, IDIAG(NEQ)
dimension A(*),B(*),X(*),ID(NDOF,*),IDX(NDOF,*)
integer          :: ierr

!local
integer :: its

if (iecho.eq.8) then
    call ParallelLog("solve", "Printing stuff")
endif

if (getrank()==0) then
    if (iecho.eq.1) then
#ifdef SPARSE
        write(stdout,1,advance='no')
 1       format(1x,'Solving linear system ...')
#else
        if (FACT.and.BACK) then
            write(stdout,1,advance='no')
 1           format(1x,'Solving linear system ...')
        else if (FACT) then
            write(stdout,2,advance='no')
 2           format(1x,'Factoring matrix ...')
        else if (BACK) then
            write(stdout,3,advance='no')
 3           format(1x,'Back-substituting ...')
        endif
#endif
        call flushf(stdout)
    endif
endif

if (iecho.eq.8) then
    call parallelLog("Solve", "calling Initsolve")
endif

#ifdef SPARSE

call InitSolve (0)

if (iecho.eq.8) then
    call parallelLog("Solve", "calling ItSolve")
endif

!write(*,*) "rank", getrank(),"B before kspsolve", B(1:10)
call ITSOLVE (B,X,its,0)
!write(*,*) "rank", getrank(),"X after kspsolve", X(1:10)

! set B equal to X.
! B will be newly set for every next iteration
! so is overwriten by X to return.
call EQUATE (B,X,NEQ)

#else
call lusolve (A,B,IDIAG,NEQ,FACT,BACK)
#endif

#ifdef SPARSE
!    if (iecho.gt.0) write(stdout,5) its
!    5   format(1X,'Sparse solution of mechanical equilibrium equation',
!     >   1x,'required ',I3,' iterations')
#else
if (.not.BACK) then
    call PIVOTS (A,IDIAG,ID,IDX,NDOF,NEQ,NUMSN, NUMNP,MODE)
    call CONDNR (A,IDIAG,NEQ,ierr)
endif
#endif

return
end
!-------------------------------------------------------------------------------
SUBROUTINE LUSOLVE (A,B,IDIAG,NEQ,FACT,BACK)
USE ALGEBRA
use constants, only: zero

implicit none
!
!.... PROGRAM TO PERFORM  A =(U)T*D*U FACTORIZATION AND/OR
!        BACKSUBSTITION OF SYMMETRIC POSITIVE DEFINITE SYSTEM OF EQUATIO
!        A(NA) COEFFICIENT MATRIX STORED IN COMPACTED COLUMN FORM
!        (AFTER FACTORIZATION CONTAINS D AND U)
!        B(NEQ)     = RIGHT SIDE VECTOR(AFTER BACKSUBSTITUTION
!                        COLUMN SOLUTION)
!        IDIAG(NEQ) = ADDRESSES OF DIAGONAL TERMS IN A(NA)
!        NEQ        = NUMBER OF EQUATIONS
!        FACT       = .TRUE. , FACTOR A(NA)
!                     .FALSE. , DO NOT FACTOR A(NA)
!        BACK       = .TRUE. , FORWARD REDUCE B(NEQ) AND BACKSUBSTITUTE
!                     .FALSE. , DO NOT FORWARD REDUCE B(NEQ) OR
!                        BACKSUBSTITUTE
!
LOGICAL          :: FACT,BACK

double precision :: a, b
integer          :: idiag, neq
DIMENSION        :: A(*),B(*),IDIAG(NEQ)

!-local
integer          :: j, jr, jd, jh, iss, ie, k, i, ir, id, ih
double precision :: d

!.... FACTOR A, REDUCE B
!
 ! LWI: replaced IS by ISS, because IS is a reserved word.
JR = 0

DO J = 1,NEQ

    JD = IDIAG(J)
    JH = JD - JR
    ISS = J - JH + 2
!      IF(JH-2) 600,300,100
    if (JH.lt.2) then
        goto 600
    elseif (JH.eq.2) then
        goto 300
    endif

    IF(.NOT.FACT) then
        GO TO 500
    endif

    IE = J - 1
    K  = JR + 2
    ID = IDIAG(ISS - 1)
!
!.... REDUCE ALL EQUATIONS EXCEPT DIAGONAL
!
    DO I = ISS,IE
        IR = ID
        ID = IDIAG(I)
        IH = MIN0(ID-IR-1,I-ISS+1)
        IF(IH.GT.0) then
            A(K) = A(K) - DOT(A(K-IH),A(ID-IH),IH)
        endif
        K = K + 1
    enddo
!
!.... REDUCE DIAGONAL TERM
!
300 IF(.NOT.FACT) then
      GO TO 500
    endif

    IR = JR + 1
    IE = JD - 1
    K = J - JD

    DO I = IR,IE

        ID = IDIAG(K+I)
        IF(A(ID).EQ.zero) then
            cycle
        endif

        D = A(I)
        A(I) = A(I)/A(ID)
        A(JD) = A(JD) - D*A(I)
    enddo

!.... REDUCE RHS
500 IF(BACK) then
        B(J) = B(J) - DOT(A(JR+1),B(ISS-1),JH-1)
    endif

600 JR = JD
enddo

    IF(.NOT.BACK) then
        RETURN
    endif
!
!.... DIVIDE BY DIAGONAL PIVOTS
!
    DO I = 1,NEQ
        ID = IDIAG(I)
        IF(A(ID).NE.zero) then
            B(I) = B(I)/A(ID)
        endif
    enddo
!
!.... BACKSUBSTITUTE
!
    J = NEQ
    JD = IDIAG(J)
800 D = B(J)
    J = J - 1

    IF(J.LE.0) then
        RETURN
    endif

    JR = IDIAG(J)

    IF(JD-JR.LE.1) then
        GO TO 1000
    endif
    ISS = J - JD + JR + 2
    K = JR - ISS + 1
    DO  I = ISS,J
        B(I) = B(I) - A(I+K)*D
    enddo

1000 JD = JR
   GO TO 800
   END
!
!-------------------------------------------------------------------------------
#ifdef SPARSE
subroutine ITSOLVE (B,X,ITS,IFLAG)
! IFLAG = 0 for mechanical solve
! IFLAG = 1 for thermal solve

! NEQ = *local* number of equations

use modeldatamodule, only: modeldatactx ! for debugging only. Can be removed
use aomodule       , only: equations
USE MODELCTX
USE MODELDEFINITION, only: NEQlocal, &
                           NEQglobal, &
                           NTEQlocal  !, &
!                                   NTEQglobal
use constants,       only: eps, dfmax
use debugmodule,     only: lustat, &
                           itest, &
                           debug, &
                           petsctime, &
                           owntime, &
                           time_start, &
                           time_finish, &
                           rate, &
                           memory_usage, &
                           iecho, &
                           LogTimeStamp, &
                           startClock, &
                           parallelLog, &
                           delay, longdelay
use meshdatamodule,  only: meshdatactx
use spetscmodule

!
implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"



!-pass
!    Not accounted for possibilty that integer differes from PetscInt
 integer ITS,IFLAG
PetscScalar B(*),X(*)
!-locl
logical Mentry,Tentry
VecScatter vscat
PetscErrorCode ierr
KSPConvergedReason reason
PetscMPIInt irank
PetscLogDouble time_sol0,time_sol1
PetscScalar, pointer   :: tmpvec(:)
integer i,j
save Mentry,Tentry
data Mentry,Tentry/.false.,.false./
! debug variables below
!PetscReal :: rnorm
double precision :: rnorm


if (IFLAG.eq.0) then
!******************************************************************
!       Solving Mechanical system
!******************************************************************

    if (iecho.eq.8) then
        call parallelLog("ITsolve", "preparing KSPSolve")
    endif

    call VecZeroEntries(M_B, ierr)
    if (ierr.ne.0) then
        write(*,*) "VecZeroEntries returned an error: ", ierr
        stop "Leaving GTecton..."
    endif

    call VecAssemblyBegin(M_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyBegin 1 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif

    call VecAssemblyEnd(M_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyEnd 1 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif


! The B vector that has been given as an argument, contains 
! load vector contributions to the nodes of this partition, but,
! on multiple partition runs, also possibly information of neighboring nodes.
! The matrix partitioning does need to match the domain partitioning exactly,
! although in practise, due to the METIS partitioning, they will likely be 
! very similar but different

! We will let PETSc handle getting this information to the proper processors.

! VecSetValue gets a global index, independent of whether the equation with that
! index belongs to this partition or not.

! VecAssembly[Begin/End] orders PETSc to figure it all out.

    i=1
    do while (1.eq.1)
!       equations array is filled with -1 entries, except for the only real equations entries.
        if (equations(i).eq.-1) then
            ! only real eqs numbers are filled up with positive values.
            exit
        else
            ! more elegant to do this with VacSetvalue*s*, TODO
            ! or with a scatter context. (see below for inspiration


            if (b(i).ne.b(i)) then
                write(*,*) "GTECTON tried to insert a NaN in position", equations(i), &
                           "of the load vector. Leaving GTecton..."
                stop "Leaving GTecton..."
            else if (b(i).gt.dfmax) then
                write(*,*) "GTECTON tried to insert a +infinity in position", equations(i), &
                           "of the load vector. Leaving GTecton..."
                stop "Leaving GTecton..."
            else if (b(i).lt.-dfmax) then
                write(*,*) "GTECTON tried to insert a -infinity in position", equations(i), &
                           "of the load vector. Leaving GTecton..."
                stop "Leaving GTecton..."
            endif



            ! adding B(i) to position equations(i) of PETSc vector M_B
            call VecSetValue(M_B, equations(i), B(i), ADD_VALUES, ierr)

            if (ierr.ne.0) then
                write(*,*) 'VecSetValue returned error: ', ierr
            endif

        endif           
        i=i+1
    enddo


    call VecAssemblyBegin(M_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyBegin 2 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif

    call VecAssemblyEnd(M_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyEnd 2 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif


    if (iecho.eq.8) then
        call parallelLog("ITsolve", "prepared load vector B")
    endif



      call VecGetArrayF90(M_X,tmpvec,ierr)
      do i=1,NEQlocal
          tmpvec(i) = X(i)
     enddo
      call VecRestoreArrayF90(M_X,tmpvec,ierr)


!*** record time stuff
!       if (getrank().eq.0) then
!           call system_clock(time_finish)
!           write(*,*) 'entering solve'
!           write(*,*) 'measured: ', time_finish, time_start, rate
!           write(*,*) 'rank', getrank(),'Gtecton time: ', dble(time_finish-time_start)/dble(rate)
!       endif
!********************

!    if (iecho.eq.5) then
!            write(*,*) '++++++++++++++++'
!            call KSPview (M_SLES,PETSC_VIEWER_STDOUT_SELF, ierr)
!            write(*,*) '++++++++++++++++'
!    endif


!#ifdef Lukas
!    call delay()
!            write(*,*) 'stiffness matrix just before solve'
!             call matview (M_A,PETSC_VIEWER_STDOUT_world, ierr)
!            write(*,*) 'load vector just before solve'
!             call vecview (M_B,PETSC_VIEWER_STDOUT_world, ierr)
!            write(*,*) 'X just before solve'
!            call vecview (M_X,PETSC_VIEWER_STDOUT_world, ierr)
!    call delay()
!#endif

!write(*,*) "rank", getrank(), "vector B", modeldatactx%B(1:4)


    if (iecho.eq.6) then
        call LogTimeStamp("Just before KSPsolve")
    endif

    if (iecho.eq.8) then
        call parallelLog("ITsolve", "calling KSPSolve")
    endif


    !*****************************************************************
    !*****************************************************************
    !       C A L L I N G    T H E     P E T S C   S O L V E R
    !*****************************************************************
    !*****************************************************************
    !*****************************************************************
    call KSPSolve(M_SLES,M_B,M_X,ierr)
    !*****************************************************************
    !*****************************************************************



    if (iecho.eq.8) then
        call parallelLog("ITsolve", "finished KSPSolve")
    endif

    if (iecho.eq.6) then
        call LogTimeStamp("Having PETSc solve a linear system") 
    endif

    if (iecho.eq.9) then
         write(*,*) 'finished KSPsolve'
    endif
    

    call KSPGetConvergedReason(M_SLES,reason,ierr)

    call HANDLE_PETSC_CONVERGENCE(reason)

!    write(*,*) "converged reason: ", reason


    !**********************************************************************

    call KSPGetIterationNumber(M_SLES,ITS,ierr)

    
    if (iecho.eq.5) then
        write(*,*) 'Solving system took ', ITS, ' iterations', ierr
    endif


    call VecGetArrayF90(M_X,tmpvec,ierr)

    do i=1,NEQlocal
        X(i) = tmpvec(i)

        if (X(i).ne.X(i)) then
            write(*,*) "PETSc has returned a NaN on position", i, &
                       "of the meachnical solution"
!            stop "Leaving GTECTON."
        else if (X(i).gt.dfmax) then
            write(*,*) "PETSc has returned an infty on position", i, &
                       "of the mechanical solution, ref val:", dfmax
!            stop "Leaving GTECTON."
        else if (X(i).lt.-dfmax) then
            write(*,*) "PETSc has returned an -infty on position", i, &
                       "of the mechanical solution, ref val:", -dfmax
!            stop "Leaving GTECTON."
        endif

    enddo
    call VecRestoreArrayF90(M_X,tmpvec,ierr)

else ! IFLAG 1
    ! Do the thermal computation
    call VecGetArrayF90(T_X,tmpvec,ierr)
    do i=1,NTEQlocal
        tmpvec(i) = X(i)
    enddo
    call VecRestoreArrayF90(T_X,tmpvec,ierr)


    call VecGetArrayF90(T_B,tmpvec,ierr)
    do i=1,NTEQlocal
        tmpvec(i) = B(i)
    enddo
    call VecRestoreArrayF90(T_B,tmpvec,ierr)

    if (iecho.eq.5) then
        write(*,*) 'rank', getrank(), 'enters thermal'
    endif


    !******************************************************************
    !       Solving Thermal system
    !******************************************************************


!#ifdef Lukas
!    call delay()
!            write(*,*) 'stiffness matrix just before solve'
!             call matview (T_A,PETSC_VIEWER_STDOUT_world, ierr)
!            write(*,*) 'load vector just before solve'
!             call vecview (T_B,PETSC_VIEWER_STDOUT_world, ierr)
!            write(*,*) 'X just before solve'
!            call vecview (M_X,PETSC_VIEWER_STDOUT_world, ierr)
!    call delay()
!#endif

!    write(*,*) "calling for a thermal solve"

    if (iecho.eq.8) then
        call parallelLog("ITsolve", "calling KSPSolve for thermal")
    endif

    call KSPSolve(T_SLES,T_B,T_X,ierr)

    if (iecho.eq.8) then
        call parallelLog("ITsolve", "finished KSPSolve for thermal")
    endif


!#ifdef Lukas
!    call delay()
!     call vecview (T_X,PETSC_VIEWER_STDOUT_world, ierr)
!    call delay()
!    call delay()!
!    call MPI_BARRIER(MPI_COMM_WORLD,iERR)

!#endif



    call KSPGetConvergedReason(T_SLES,reason,ierr)

    call HANDLE_PETSC_CONVERGENCE(reason)

    call KSPGetIterationNumber(T_SLES,ITS,ierr)

    if (iecho.eq.5) then
        write(*,*) 'convergence after: ', ITS, 'iterations'
    endif

    call VecGetArrayF90(T_X,tmpvec,ierr)

    do i=1,NTEQlocal
    !              write(FILE_outputf,*) i, X(i)
        X(i) = tmpvec(i)
        if (X(i).ne.X(i)) then
            write(*,*) "PETSc has returned a NaN on position", i, &
                       "of the thermal solution"
            stop "Leaving GTECTON."
        else if (X(i).gt.dfmax) then
            write(*,*) "PETSc has returned an infty on position", i, &
                       "of the thermal solution"
            stop "Leaving GTECTON."
        else if (X(i).lt.-dfmax) then
            write(*,*) "PETSc has returned an -infty on position", i, &
                       "of the thermal solution"
            stop "Leaving GTECTON."
        endif

    enddo

    call VecRestoreArrayF90(T_X,tmpvec,ierr)

endif



#ifdef TIMING
   if (irank.eq.0) then
!          call PetscGetTime(time_sol1, ierr)
       time_soltot = time_soltot + (time_sol1 - time_sol0)
   endif
#endif
! TIMING


if (iecho.eq.5) then
    write(*,*) 'rank', getrank(), 'done thermal'
endif

if (iecho.eq.8) then
    call parallelLog("ITsolve", "Leaving subroutine")
endif


return
end subroutine

#endif
! SPARSE
!----------------------------------------------------------------------

subroutine HANDLE_PETSC_CONVERGENCE(reason)

use debugmodule, only: iecho
use petscksp
use iomodule

implicit none

!#include "/net/home/gtecton/openmpi/x86_64/include/mpif.h"
! for stderr

!#include "petsc/finclude/petsckspdef.h"

#include "petsc/finclude/petscksp.h"  
! for the convergence codes.

!#include "petsc/finclude/petsc.h"
!
!
!

KSPConvergedReason reason

!-----------  converged reason; show only with iecho 5 (petsc log) --------------------------

if (iecho.eq.5) then
  if (reason.eq.KSP_CONVERGED_RTOL) then
      write(*,2)
      write(stderr,2)
 2         format(1X,'ITSOLVE: Reached relative convergence')
  else if (reason.eq.KSP_CONVERGED_ATOL) then
      write(*,4)
      write(stderr,4)
 4         format(1X,'ITSOLVE: Reached absolute convergence')
  else if (reason.eq.KSP_CONVERGED_ITS) then
      write(*,5)
      write(stderr,5)
 5         format(1X,'ITSOLVE: Reached CG new curve')
  else if (reason.eq.KSP_CONVERGED_CG_NEG_CURVE) then
      write(*,6)
      write(stderr,6)
 6         format(1X,'ITSOLVE: Reached CG negative curve.')
  else if (reason.eq.KSP_CONVERGED_CG_CONSTRAINED) then
      write(*,7)
      write(stderr,7)
 7         format(1X,'ITSOLVE: Reached CG constrained.')
  else if (reason.eq.KSP_CONVERGED_STEP_LENGTH) then
      write(*,8)
      write(stderr,8)
 8         format(1X,'ITSOLVE: Reached step length convergence')
  elseif (reason.eq.KSP_CONVERGED_HAPPY_BREAKDOWN) then
      write(*,9)
      write(stderr,9)
 9         format(1X,'ITSOLVE: Reached happy breakdown')
  endif
endif
!-----------  diverged reason; is an error; show always -------------------------------------------------
if (reason.eq.KSP_DIVERGED_INDEFINITE_PC) then
    write(*,12)
    write(stderr,12)
 12       format(1X,'ITSOLVE: Divergence because of indefinite',1X, &
         'preconditioner')
    call stoper()

else if (reason.eq.KSP_DIVERGED_ITS) then
    write(*,14)
    write(stderr,14)
 14       format(1X,'PETSc solver: Ran out of iterations before any',1x, &
         'convergence criteria was reached ')
    call stoper()

else if (reason.eq.KSP_DIVERGED_DTOL) then
    write(*,15)
    write(stderr,15)
15  format(1X,'PETSc solver returned KSP_DIVERGED_DTOL'/&
           1X,'norm of the residual >= dtol * (norm of the load vector)')
    write(stderr,*) "Possible cause: Matrix too stiff"
    write(stderr,*) "Suggestion: Reduce viscosity or increase magnitude of boundary conditions."
    call stoper()

else if (reason.eq.KSP_DIVERGED_BREAKDOWN) then
    write(*,16)
    write(stderr,16)
 16       format(1X,'PETSc solver: A breakdown in the Krylov method',1x, &
      'was detected so the method could not continue to',1x, &
      'enlarge the Krylov space')
    call stoper()

else if (reason.eq.KSP_DIVERGED_BREAKDOWN_BICG) then
    write(*,17)
    write(stderr,17)
 17       format(1X,'PETSc solver: A breakdown in the KSPBICG method', &
      ' was detected so the method could not continue to',1x, &
      'enlarge the Krylov space')
    call stoper()

else if (reason.eq.KSP_DIVERGED_NONSYMMETRIC) then
    write(*,18)
    write(stderr,18)
 18       format(1X,'PETSc solver: It appears the operator or',1X, &
      'preconditioner is not symmetric and this Krylov ', &
      'method (KSPCG, KSPMINRES, KSPCR) requires symmetry')
    call stoper()
endif

 end subroutine



!-------------------------------------------------------------------------------
subroutine STRAIN (XL,IEN,D,DX,TFAULT,LMX,FaultElemID, &
         SKEW,EE,NDOF,NSD,NEN,NSTR,elemID,IERR)

use modeldatamodule, only: modeldatactx
USE MODELDEFINITION
USE MODELCTX
!
!    Computes infinitesimal strain tensor
!

implicit none
!-pass
integer :: NSTR,NDOF,NSD,NEN
integer :: elemID, FaultElemID
integer ::          ien, lmx
double precision :: xl, d, dx, tfault, ee, skew
dimension :: XL(NSD,NEN),IEN(NEN),D(NDOF,*),DX(NDOF,*), &
    TFAULT(NDOF,*), &
    LMX(NDOF,NEN),EE(NSTR),SKEW(2,*)
integer :: ierr
!-locl
double precision :: dl(24)

dl = 0d0
! localize displacement
call LDISP (dl,D,IEN,NDOF,NEN)

! add faulted node displacement
call AddFaultDisplacement (dl,FaultElemID,TFAULT,NDOF,NEN,50)

if (NUMSLP.gt.0) then

    ! construct. Adds slip from slippery nodes to dl
#ifdef EXPERIMENTAL_ReducedSlipAlloc
    call ADDSN  (dl,DX,IEN,LMX,modeldatactx%SKEW,NDOF,NEN, elemID)
#else
    call ADDSN  (dl,DX,IEN,modeldatactx%SKEW,NDOF,NEN, elemID)
#endif

endif

call BDELD (XL,DL,EE,IEN(3).ne.IEN(4),IERR)

return
end subroutine


!-------------------------------------------------------------------------------
subroutine STRAINRate (XL,IEN,D,DX,TFAULT,LMX,FaultElemID, &
         SKEW,EE,NDOF,NSD,NEN,NSTR,elemID,IERR)

use modeldatamodule, only: modeldatactx
USE MODELDEFINITION
USE MODELCTX
use meshdatamodule,  only: meshdatactx
!use debugmodule,     only: iecho
!
!    Computes infinitesimal strain tensor
!
! output = EE


implicit none
!-pass
integer :: NSTR,NDOF,NSD,NEN
integer :: elemID, FaultElemID
integer ::          ien, lmx
double precision :: xl, d, dx, tfault, ee, skew
dimension :: XL(NSD,NEN),IEN(NEN),D(NDOF,*),DX(NDOF,*), &
    TFAULT(NDOF,*), &
    LMX(NDOF,NEN),EE(NSTR),SKEW(2,*)
integer :: ierr
!-locl
double precision :: dl(24)

!write(*,*) "STRAIN called"


dl = 0d0

! todo clean up this bloody mess

! localize displacement
call LDISP (dl,D,IEN,NDOF,NEN)

!write(*,*) "STRAIN SAYS: after ldisp; dl: ", elemID,dl(1:8)

call AddFaultDisplacement (dl,FaultElemID,TFAULT,NDOF,NEN,50)


!if (iecho.eq.9) then
!write(*,*) "STRAIN SAYS: before addsn; dl: ", dl(1:6)
!endif

if (NUMSLP.gt.0) then
#ifdef Lukas
!    write(*,*) "STRAIN says: calling ADDSN with LMX", LMX
!    write(*,*) "STRAIN says: calling ADDSN with MD LMX", modeldatactx%LMX(1:NDOF,1:NEN,FaultElemID)
!    write(*,*) "STRAIN says: calling ADDSN with IEN", IEN
#endif
    ! construct. Adds effect of DX to dl
#ifdef EXPERIMENTAL_ReducedSlipAlloc
!    call ADDSN  (dl,DX,IEN,LMX,modeldatactx%SKEW,NDOF,NEN, elemID)
#else
!    call ADDSN  (dl,DX,IEN,modeldatactx%SKEW,NDOF,NEN, elemID)
#endif
!    call ADDSN(dl,DX,IEN,modeldatactx%LMX(1:NDOF,1:NEN,elemID),modeldatactx%SKEW,NDOF,NEN, elemID)
endif

!write(*,*) "STRAIN SAYS: after addsn; dl: ", elemID,dl(1:8)

call BDELD (XL,DL,EE,IEN(3).ne.IEN(4),IERR)

!write(*,*) "STRAIN SAYS: after bdeld; dl: ", elemID,dl(1:8)
!write(*,*) "strain says; EE:  ", elemID,EE(1:3)


return
end subroutine



!-------------------------------------------------------------------------------
subroutine STRESS (X,IEN,D,STN,LMF,TFAULT,LMX,DX, &
                   DMAT,SKEW)

USE MODELDEFINITION
use iomodule          
USE MODELCTX
USE ALGEBRA
use meshdatamodule, only: meshdatactx
use modeltopology, only: NDOF,NSD,NSTR,NEN
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: elementSlipperyID
#endif
!
!    Computes stress tensor
!

implicit none
!-pass
!integer :: NSTR,NDOF,NSD,NEN
integer          :: LMX, ien, lmf
double precision :: STN,TFAULT,X,D,DX,DMAT,SKEW
dimension :: STN(NSTR,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
LMX(NDOF,NEN,*),X(NSD,*),IEN(NEN,*),D(NDOF,*),DX(NDOF,*), &
DMAT(NSTR,NSTR,*),SKEW(2,*)
!-locl
integer :: i, ierr
integer :: n, nn
double precision :: ee(6),xl(24),stnl(6)

integer :: thisElem

#ifdef Lukas
!write(*,*) "COMPUTING STRESS"
#endif

do nn=1,meshdatactx%nelocal

    thisElem = meshdatactx%locElt(nn)
 
    n = nn

    ! make a small array with 
    call LCOORD (X,xl,meshdatactx%IEN(1,n),n)

#ifdef EXPERIMENTAL_ReducedSlipAlloc

!TODO, dit gaat nog keihard mis

    if (elementSlipperyID(n,26).gt.0) then
        call STRAIN (xl,meshdatactx%IEN(1,n),D,DX,TFAULT,LMX(1,1,elementSlipperyID(n,26)), &
                     n,SKEW,ee,NDOF,NSD,NEN, &
                     NSTR,n,ierr)
    endif

#else
    call STRAIN (xl,meshdatactx%IEN(1,n),D,DX,TFAULT,LMX(1,1,n), &
                 n,SKEW,ee,NDOF,NSD,NEN, &
                 NSTR,n,ierr)
#endif

    if (ierr.ne.0) then
        write(stderr,1) n
 1         format(1x,'stress: error in element ',I12)
        call stoper()
    endif

#ifdef Lukas
!    write(*,*) "stress after strain ; dmat:  ", dmat(1:NSTR, 1:NSTR, nn)
#endif


!   Calculate product of matrix DMAT and vector EE
    call VMPRD (DMAT(1,1,n),ee,stnl,NSTR,NSTR)

#ifdef Lukas
!    write(*,*) 'stress says: stnl: ', stnl
#endif

    do i=1,NSTR
        STN(i,n) = STN(i,n) + stnl(i)
!        STN(i,n) = stnl(i)
    enddo


enddo
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine STRAIS (X,IEN,D,STN,LMF,TFAULT,LMX,DX,SKEW)

USE MODELDEFINITION
use modeltopology, only: NSTR,NDOF,NSD,NEN
use meshdatamodule, only: meshdatactx
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: elementSlipperyID
#endif
use iomodule
!
!    Computes strain tensor in all elements
!

implicit none
!-pass

integer          :: lmf, lmx, ien
double precision :: stn, tfault, x, d, dx, skew

 dimension STN(NSTR,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   LMX(NDOF,NEN,*),X(NSD,*),IEN(NEN,*),D(NDOF,*),DX(NDOF,*), &
   SKEW(2,*)
!-locl
double precision :: xl(24)
integer :: n, nn, ierr


do nn=1,NUMEL
    n = nn
    call LCOORD (X,xl,IEN(1,n),nn)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
    if (elementSlipperyID(n,28).gt.0) then
        call STRAIN (xl,meshdatactx%IEN(1,n),D,DX,TFAULT,LMX(1,1,elementSlipperyID(n,28)), &
                     n,SKEW,STN(1,n),NDOF,NSD,NEN, &
                     NSTR,nn,ierr)
    endif
#else 
    call STRAIN (xl,meshdatactx%IEN(1,n),D,DX,TFAULT,LMX(1,1,n), &
                 n,SKEW,STN(1,n),NDOF,NSD,NEN, &
                 NSTR,nn,ierr)
#endif


    if (ierr.ne.0) then
        write(stderr,1) n
 1           format(1x,'strais: error in element ',I6)
        call stoper()
    endif
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine WINKLR(ID,IWINK,WINK,NDOF,NUMNP,winkflag)

use timestepmodule,  only: nstep, ITWORK, refact

#ifdef SPARSE
use spetscmodule,    only: PETSc_MatSetValuesLocal, M_A
#endif

#ifdef EXPERIMENTAL_CYCLIC
use modeldefinition, only: UseCyclicSlip
#endif
!
! Routine to implement Winkler restoring forces on nodes designated by
! IWINK. This program adds a constant given by WINK to the corresponding
! diagonal elements of the global stiffness matrix.
!
! winkFLAG = 0: Nodal Winkler forces
! winkFLAG = 1: Differential Winkler forces



 implicit none
!-pass
 integer NDOF,NUMNP, winkflag
integer :: iwink, id
double precision :: wink
! dimension A(*),IDIAG(NEQ),
dimension :: ID(NDOF,*),IWINK(NDOF,*),WINK(NDOF,*)
!local
integer :: i,j,k,LL,itshift
integer :: ierr
integer :: ipdiag
integer :: nstep_minus_itshift

#ifdef SPARSE


#include "petsc/finclude/petscvec.h"

   PetscMPIInt irank
#endif

do j=1,NUMNP
#ifdef EXPERIMENTAL_CYCLIC
#ifdef ZISFAULTNORMAL
! RG 23/11/2018: allow for time delay in cyclic locking/unlocking
    do i=1,NDOF-1
        itshift = IWINK(NDOF,j)
#else
    do i=1,NDOF
        itshift = 0
#endif 
#else
    do i=1,NDOF
#endif
        LL = IWINK(i,j)

        if (winkflag.eq.0) then
            ! ordinary winkler forces
            if (LL.lt.0) then
                if (NSTEP.ge.-LL) cycle
            else if (LL.gt.0) then
                if (NSTEP.lt.(LL-1)) cycle
            else
                cycle
            endif
        else
            ! differential Winkler forces
#ifdef EXPERIMENTAL_CYCLIC
            if (UseCyclicSlip) then
                nstep_minus_itshift = NSTEP - itshift
                if (nstep_minus_itshift.lt.0) then
                    nstep_minus_itshift = nstep_minus_itshift + IABS(LL)
                endif
                ! Interpretation of IWINK as periodic
                if (LL.lt.0) then
                    ! typically locked, except when MOD(NSTEP,ABS(LL)) < 10
                    if (MOD(nstep_minus_itshift,IABS(LL)).lt.10) cycle
                else if (LL.gt.0) then
                    ! typically unlocked, except when MOD(NSTEP,LL) < 10
                    if (MOD(nstep_minus_itshift,LL).ge.10) cycle
                else
                    cycle
                endif
            else
                if (LL.lt.0) then 
                    ! no Winkler force when NSTEP >= -IWINX
                    if (NSTEP.ge.-LL) cycle
                else if (LL.gt.0) then
                    ! no Winkler force when NSTEP < (IWINX-1)
                    if (NSTEP.lt.(LL-1)) cycle
                else
                    cycle
                endif
            endif
#else
            if (LL.lt.0) then 
                ! no Winkler force when NSTEP >= -IWINX
                if (NSTEP.ge.-LL) cycle
            else if (LL.gt.0) then
                ! no Winkler force when NSTEP < (IWINX-1)
                if (NSTEP.lt.(LL-1)) cycle
            else
                cycle
            endif
#endif
        endif ! added RG2021-02-26

        k = ID(i,j)
        if (k.le.0) cycle ! no equation associated with this node

#ifdef SPARSE
        if (winkflag.eq.0) then
            ! normal winler force
            call PETSc_MatSetValuesLocal(M_A,k-1,k-1,WINK(i,j),77,0,0,0,ierr)
        else
            ! differential winkler forces
            call PETSc_MatSetValuesLocal(M_A,k-1,k-1,2d0 * WINK(i,j),17,0,0,0,ierr)
        endif
#endif
    enddo
enddo
 
return
end subroutine
!-------------------------------------------------------------------------------
subroutine followtime (mode)

use debugmodule, only: iecho
use iomodule


implicit none
!-pass
integer mode
!-locl
real tarray(2)
character(len=24) date,fdate
integer :: iflag
integer :: tTotal
save Ttotal,iflag
!external dtime
data iflag/0/

if (iecho.eq.1) then

! LWItodo, set fdate to th 2008 std.
!        date = 'date not implemented'
    if (mode.eq.0) then
        if (iflag.ne.-12345) then
            iflag =-12345
            Ttotal=0
!            call dtime(tarray)
        else
!            call dtime(tarray)
!            user=DBLE(tarray(1))
!            sys =DBLE(tarray(2))
!            Ttotal=Ttotal+user+sys

!                write(stdout,2) date(12:19),user,sys,user+sys
!    2               format(/1x,a8,' - previous time step: usr=',F7.1, &
!                     ' sys=',F7.1,' elapsed=',F7.1)
        endif
    else
!        call dtime(tarray)
!        user=DBLE(tarray(1))
!        sys =DBLE(tarray(2))
!        Ttotal=Ttotal+user+sys
        write(stdout,3) date(12:19),dble(Ttotal),dble(Ttotal/mode)

 3           format(/1x,a8,' - FE calculation ready:'/ &
         1x,'total elapsed time=',F9.1,' real time average=', &
         F7.1,' per step')

    endif
endif

return
end subroutine
!-------------------------------------------------------------------------------
subroutine stoper()

USE MATERIALSMODULE
USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR
USE MODELDATAMODULE
USE MODELCTX
use debugmodule,     only: xit
use iomodule          
use timestepmodule,  only: nstep

!
! stops program on error, and prevents the use from finding out anything
! about the error, which is kind of annoying.
!

implicit none
!
call flushf(stdout)

call INTCON ()

if (OUTFIL(1)) then
    write(*,*) "Closing feout file"
    call closef(luout(1))
endif
call closef(luout(2))

call xit(1," ")

!!!!!!!!!!!! NO ATTEMPT TO WRITE FECORE FILE --- IN PARALLEL USELESS
luout(1) = 11
call openf (luout(1),'fecore','unknown')

OUTFIL(1) = .true.
OUTFIL(2) = .false.
!
write(luout(1),1) NSTEP
 1   format(/1x,'STOPER: ERROR FORCED EXIT'//5X,'NSTEP = ',I6//)

call PRINTID ()
!
if (MODE.eq.5 .or. MODE.eq.6) then
    call PRINTT (modeldatactx%T,NUMNP)
    call PRINTF (modeldatactx%FLUX,NSD,NUMEL)
endif

if (LMAT) then
    call WRTCRP (matdatactx%prpmat,NUMEL)
    call WRTPLS (NUMEL)
else
    call WRTCRP (matdatactx%prpmat,NUMAT)
    call WRTPLS (NUMAT)
endif

if (MODE.eq.2 .or. MODE.eq.6) then
    !write displacement
    call PRINTD (modeldatactx%D,NUMNP,0)
    ! write differential displacement
!    call PRINTX (modeldatactx%DX,modeldatactx%IDSLP, &
!        NDOF,NUMSN,0)
    call PRINTX (NDOF,NUMSN,0)
    call PRINTS (modeldatactx%STN,NSTR,NUMEL)
endif

call closef(luout(1))
!
call xit(1," ")

end subroutine
!-------------------------------------------------------------------------------
subroutine SETBTOT (BTOT,B,NEQ)

USE ALGEBRA,     only: equate
use filemodule,  only: outdir
use debugmodule, only: iecho, xit
use iomodule

!

implicit none
!-pass
integer NEQ
double precision :: B(NEQ),BTOT(NEQ)
!-locl
logical :: there
integer :: i
integer :: lutmp
integer :: ios
integer, external :: nextlu
!
inquire(file=trim(outdir)//'/'//'BTOT.DAT',exist=there)

if (there) then
    lutmp = nextlu(0)

    call fbopen (lutmp,trim(outdir)//'/'//'BTOT.DAT','old')

    read(lutmp,err=100,end=200,iostat=ios) (BTOT(i),i=1,NEQ)

    call closef(lutmp)

    if (iecho.ne.0) then
        write(stderr,1)
    endif
 1       format(1x,'reading "BTOT.DAT"')

else
    write(stderr,2)
 2       format(1x,'Assuming that initial model stresses',1x, &
     'represent equilibrium')
    call EQUATE (BTOT,B,NEQ)
endif
!
return

100     write(stderr,101) ios
  101   format(1x,'SETBTOT: read error ',I5)
call xit(1," ")

200     write(stderr,201) ios
  201   format(1x,'SETBTOT: EOF error ',I5)
call xit(1," ")

end subroutine
!-------------------------------------------------------------------------------
subroutine PRESTN (B,STN0,ISELM,ISTIME,X,D,SKEW,IEN,LM,LMX,MAT, &
                   PRPMAT,LMF,TFAULT,DXE)

USE MATERIALSMODULE, only: MaterialIDOfElement
USE MODELDEFINITION
USE ALGEBRA
use modeltopology,   only: NSD,NSTR,NEN,NDOF

#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule, only: HasSlipperyNodes, &
                           elementSlipperyID
#else
use modeldatamodule, only: HasSlipperyNodes
#endif

use debugmodule,     only: iecho
use constants,       only: useDXE
use iomodule          
use timestepmodule,  only: nstep, refact

!
! Adds forces from initial stresses (boundary conditions) to global load vector.
!

 implicit none
!-pass

integer ::          ien, lm, lmx, iselm, istime, lmf, mat
double precision :: b, x, skew, stn0, prpmat, d, tfault, dxe

 dimension B(*),X(NSD,*),SKEW(2,*),STN0(NSTR,*),IEN(NEN,*), &
   LM(NDOF,NEN,*),LMX(NDOF,NEN,*),MAT(*),PRPMAT(6,*),D(NDOF,*), &
   ISELM(*),ISTIME(*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
 DXE(NDOF,NEN,*)
!-locl
 logical quad
double precision :: pl(24),xl(24),dl(24)
integer :: it, m, n, nn, numadd
double precision :: thick
integer   :: ierr
!-init
 numadd = 0
!

if (NPRE.eq.0) then
    return
endif

do nn=1,NPRE
    n = ISELM(nn)
    it  = ISTIME(nn)

    ! check whether we apply the stress at the proper timestep.
    if ((it-1).ne.NSTEP .or. it.eq.0) then
        cycle
    endif

numadd = numadd + 1
quad = (IEN(3,n).ne.IEN(4,n))

call CLEAR  (pl,NEE,"pl")
call LCOORD (X,xl,IEN(1,n),0)

if (LGDEF.ne.0) then
    call LDISP  (dl,D,IEN(1,n),NDOF,NEN)
    call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,51)
    call ADDSNE (dl,n,NDOF,NEN,useDXE)
    call REZONE (xl,dl,1d0)
endif

m = MaterialIDOfElement(n)

thick = PRPMAT(6,m)

call EFORCE (STN0(1,nn),pl,xl,thick,quad,ierr)

if (ierr.ne.0) then
    write(stderr,1) n
 1           format(1x,'prestn: error in element ',I5)
    call stoper()
endif

if (NUMROT.ne.0) call RPFORC (pl,SKEW,IEN(1,n),NDOF,NEN)

call ADDFOR (B,pl,LM(1,1,n),NEE,nn,6)  ! prestress



if (HasSlipperyNodes()) then

#ifdef EXPERIMENTAL_ReducedSlipAlloc
    if(elementSlipperyID(n,29).gt.0) then
        call ADDFOR (B,PL,LMX(1,1,elementSlipperyID(n,29)),NEE,nn,7)
    else
        ! this element has no slippery nodes. Skip it.
    endif
#else
    call ADDFOR (B,PL,LMX(1,1,n),NEE,nn,7)
#endif

endif

enddo
!
 if (iecho.eq.1 .and. numadd.gt.0) write(stdout,2)
 2   format(1x,'Pre-stress loads are being added')
!
 return
 end
!-------------------------------------------------------------------------------
subroutine UPDXL (DL,DX,IEN,LMX,SKEW,NDOF,NEN, elemID)

USE MODELDEFINITION
USE ALGEBRA,         only: formrt, transp, vmprd
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldatamodule,   only: modeldatactx, &
                             elementSlipperyID
#endif
use constants,         only: zero, one
!
!       update element based differential displacements DL

!   This routine is called for every element that has slippery nodes.


implicit none
!-pass
integer NDOF,NEN
integer :: ien, lmx
double precision :: dl, dx, skew
dimension DL(NDOF,NEN),DX(NDOF,*),IEN(NEN),LMX(NDOF,NEN), &
   SKEW(2,*)
!-locl
double precision :: rot(3,3),dxl(3),dxm(3)
external FLOAT
integer :: elemID
integer :: i,j,k
integer :: sgn
intrinsic SIGN

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer :: slipperySequenceNr

slipperySequenceNr = elementSlipperyID(elemID,30)
if (slipperySequenceNr.eq.0) then
    ! this element has no slippery nodes
    return
endif

#endif

do j=1,NEN

    dxl(1) = 0.0d0
    dxl(2) = 0.0d0

    ! loop through the nodes of this element.
    k = IEN(j)

    if (NUMROT.eq.0.or.&
       ! not rotated
       (SKEW(1,k).eq.ZERO.and. &
        SKEW(2,k).eq.ZERO)) then

        do i=1,NDOF
#ifdef EXPERIMENTAL_ReducedSlipAlloc
            if(slipperySequenceNr.gt.0) then
                sgn = SIGN(ONE,dble(modeldatactx%LMX(i,j,slipperySequenceNr)))
                DL(i,j) = DL(i,j) + sgn*DX(i,k)
            else
                ! element has no slippery node, no need to add slippery displacement.
                ! more over, it no longer has a place in LMX
            endif
#else
            if (LMX(i,j).ne.0) then
                ! this element/vertex/DOF is slippery
                sgn = SIGN(ONE,dble(LMX(i,j)))
                DL(i,j) = DL(i,j) + sgn*DX(i,k)
            endif
#endif
        enddo
    else
        ! rotated

!       rotate differential displacements to local (skewed) system
        call FORMRT (SKEW(:,k),rot,NDOF)
        call TRANSP (rot,3)
        call VMPRD  (rot,DX(1,k),dxl,NDOF,3)

!       add (local!) sign information to differential displacements
        do i=1,NDOF
#ifdef EXPERIMENTAL_ReducedSlipAlloc
!            slipperySequenceNr = elementSlipperyID(elemID,31)
            if(slipperySequenceNr.gt.0) then
                if (modeldatactx%LMX(i,j,slipperySequenceNr).ne.0) then
                    dxl(i) = dxl(i)*SIGN(ONE,dble(modeldatactx%LMX(i,j,slipperySequenceNr)))
                else
                    dxl(i) = 0.0d0
                endif
            else
                dxl(i) = 0.0d0
            endif
#else
            if (LMX(i,j).ne.0) then
                dxl(i) = dxl(i)*SIGN(ONE,dble(LMX(i,j)))
            else
                dxl(i) = 0.0d0
            endif
#endif
        enddo

!       back-rotate differential displacement to global system
        call TRANSP (rot,3)
        call VMPRD  (rot,dxl,dxm,NDOF,3)
        do i=1,NDOF
            DL(i,j) = DL(i,j) + dxm(i)
        enddo

    endif
    

enddo

!
return
end
!-------------------------------------------------------------------------------
subroutine UPDXE (DX,IEN,LMX,SKEW,IDSLE,DXE,NDOF,NEN)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
USE MODELDEFINITION, only: nElemsWithSlipperyNodes, &
                           nLocalElemsWithSLipperyNodes
use modeldatamodule, only: modeldatactx, &
                           elementSlipperyID
#else
USE MODELDEFINITION, only: numse
use modeldatamodule, only: modeldatactx
#endif

use meshdatamodule,  only: meshdatactx


#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modelctx,        only: getrank
#endif
!
!       update differential displacements DXE
!
! DXE is either modeldatactx%DXE (when called from elastc)
! or modeldatactx%DELX, (when called from visco)
! so do not replace with modeldatactx%something.

! However, both are of exactly the same size,
! so we can set the size explicitly

implicit none
!-pass
integer NDOF,NEN

integer :: ien(nen, meshdatactx%nelocal)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer          :: LMX(NDOF,NEN,nElemsWithSlipperyNodes)
double precision :: DXE(NDOF,NEN,nElemsWithSLipperyNodes)
integer          :: idsle(nLocalElemsWithSLipperyNodes)
#else
integer          :: LMX(NDOF,NEN,meshdatactx%neglobal)
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
integer          :: idsle
dimension        :: IDSLE(*)
#endif

double precision :: dx, skew
dimension        :: DX(NDOF,*),SKEW(2,*)

!-locl
integer nn,n

#ifdef EXPERIMENTAL_ReducedSlipAlloc
integer :: iSlipElem, globalElemID, localElemID, positionInIEN, elemID
integer :: slipperySequenceNr

!if (getrank().eq.0) then
!write(*,*) "rank", getrank(), "UPDXE says; 1 ",nElemsWithSlipperyNodes
!write(*,*) "rank", getrank(), "UPDXE says; nsliploc ",nLocalElemsWithSLipperyNodes
!write(*,*) "rank", getrank(), "UPDXE says; IDSLE ",IDSLE(1:nLocalElemsWithSLipperyNodes)
!write(*,*) "rank", getrank(), "UPDXE says; lookup: ",modeldatactx%SlipElemLookup
!write(*,*) "rank", getrank(), "UPDXE says; gloelt: ", meshdatactx%gloelt
!write(*,*) "rank", getrank(), "UPDXE says; locelt: ", meshdatactx%locelt
!write(*,*) "rank", getrank(), "UPDXE says--------------------------------------------------"
!endif

do iSlipElem = 1,nLocalElemsWithSLipperyNodes

!    write(*,*) "rank", getrank(), "UPDXE says; IDSLE:", IDSLE(1:nLocalElemsWithSLipperyNodes)


    elemID = IDSLE(iSlipElem) ! elemID = local counter.
                              ! Note that this is not the same as the local
                              ! element number! The local element number counts
                              ! through the various partitions.
                              ! This one starts at 1 for every partition!


    globalElemID = meshdatactx%gloelt(elemID)
    slipperySequenceNr = elementSlipperyID(globalElemID, 170)

!if (getrank().eq.0) then
!    write(*,*) "rank ", getrank(), "in UPDXE step", iSlipElem, globalElemID, slipperySequenceNr
!    write(*,*) "rank ", getrank(), "with LMX ", LMX(:,:,slipperySequenceNr)
!    write(*,*) "rank ", getrank(), "with DXE ", DXE(:,:,slipperySequenceNr)
!    write(*,*) "rank ", getrank(), "with IEN ", IEN(:,elemID)
!endif

    call UPDXL (DXE(1,1,slipperySequenceNr), & ! to be tested.
                DX, &
                IEN(1,elemID), &  ! tested, OK
                LMX(1,1,slipperySequenceNr), & ! 
                modeldatactx%SKEW, &
                NDOF, &
                NEN, &
                globalElemID)
enddo

#else

do nn=1,NUMSE
    n = IDSLE(nn)
    ! set DXE based on values in DX
    call UPDXL (DXE(1,1,n),DX,IEN(1,n),LMX(1,1,n),modeldatactx%SKEW,NDOF, NEN, n)
     
enddo

#endif
!
return
end
!-------------------------------------------------------------------------------
double precision function Accur (alpha,lambda,DeltaT)

USE ALGEBRA
use constants, only: one, two, three, six, dlogmax

! calculates relative modal error


implicit none

double precision, parameter :: small=1.D-30
double precision, parameter :: great=1.D+06
!-pass
! SNGL    real alpha,lambda,DeltaT
double precision :: alpha,lambda,DeltaT
!-locl
! SNGL    real logLam,logDel
double precision :: logLam,logDel
logical :: ov

logLam = TENLOG(lambda)
logDel = TENLOG(DeltaT)

! overflow check
ov = .false.
ov = (three*logLam              .gt. dlogmax .or. ov)
ov = (three*logDel              .gt. dlogmax .or. ov)
ov = (three*logDel+three*logDel .gt. dlogmax .or. ov)

if (ov) then
    Accur = great
else
    Accur = ABS( one - two * alpha ) / two &
            * lambda**two * DeltaT &
          + ABS( one + three * alpha * (alpha-one) ) / six &
        * lambda**three * DeltaT**two
endif

if (Accur.gt.great) then
    Accur = great
endif

if (Accur.lt.small) then
    Accur = small
endif

return
end
!-------------------------------------------------------------------------------
double precision function Oscill (alpha,lambda,DeltaT)

use constants, only: zero, one, dfmin

implicit none

double precision, parameter :: small=1.D-30
double precision, parameter :: great=1.D+06
!-pass
double precision :: alpha,lambda,DeltaT
!-locl
double precision :: dT

if ( (one-alpha) .gt. dfmin ) then
    ! conditional oscillation behavior
    if (lambda.gt.dfmin) then
        ! calculate time step at oscillation limit
        dT = one / ( (one - alpha) * lambda )
        Oscill = DeltaT/dT
        if (Oscill.gt.great) then
            Oscill = great
        endif
    else
        Oscill = great
    endif
else
    ! unconditional oscillation stability
    Oscill = zero
endif

if (Oscill.lt.small) then
    Oscill = small
endif

return
end
!-------------------------------------------------------------------------------
double precision function Stabil (alpha,lambda,DeltaT)

use constants, only: zero, one, two, dfmin

implicit none

double precision, parameter :: small=1.D-30
double precision, parameter :: great=1.D+06
!-pass
! SNGL    real alpha,lambda,DeltaT
double precision :: alpha,lambda,DeltaT
!-locl
! SNGL    real dT
double precision :: dT

if (alpha.lt.0.5) then
    ! conditional stability
    if (lambda.gt.dfmin) then
        ! calculate time step that is required for stability
        dT = two / ( (one - two*alpha) * lambda )
        Stabil = DeltaT/dT
        if (Stabil.gt.great) then
            Stabil = great
        endif
    else
        Stabil = great
    endif
else
    ! unconditional stability
    Stabil = zero
endif

if (Stabil.lt.small) then
    Stabil = small
endif

return
end
!-------------------------------------------------------------------------------

integer function ieq(i,n,NDOF)

implicit none

integer :: i,n,NDOF

ieq = (n -1) * NDOF + i

end function

!-------------------------------------------------------------------------------

subroutine ADDNEQBC (SL,IEN,S,IA,IAX,NEN,NEE,NDOF,MEQ)

USE MODELDEFINITION, only: NUMFN, NUMNP
!
! adds elements of local stiffness matrix to the system of equations which
! do not actually become part of the stiffness matrix.
!
! todo... does whaaaaat?
!
implicit none
!-pass
integer   :: NEN,NEE,NDOF,MEQ
integer   :: ien, ia(NDOF,NUMNP), iax(NDOF,NUMNP)
double precision :: sl, s
dimension :: SL(NEE,NEE),IEN(NEN),S(MEQ,*)
!-locl
integer   :: nn,n,i,leq,geq,mm,m,j,keq,heq
integer, external :: ieq

do nn=1,NEN
    n = IEN(nn)
    do i=1,NDOF
        geq = IA(i,n)
        if (geq.ne.0) then
            leq = ieq(i,nn,NDOF)
            do mm=1,NEN
                m = IEN(mm)
                do j=1,NDOF
                    heq = ieq(j,m,NDOF)
                    keq = ieq(j,mm,NDOF)
                    S(geq,heq) = S(geq,heq) + SL(leq,keq)
                enddo
            enddo
        endif

        geq = 0

        if (NUMFN.gt.0) then
            geq = IAX(i,n)
        endif

        if (geq.ne.0) then
            leq = ieq(i,nn,NDOF)
            do mm=1,NEN
                m = IEN(mm)
                do j=1,NDOF
                    heq = ieq(j,m,NDOF)
                    keq = ieq(j,mm,NDOF)
                    S(geq,heq) = S(geq,heq) + SL(leq,keq)
                enddo
            enddo
        endif
    enddo
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine GETNEQBC (IBOND,S,IA,DELD,FBC,NDOF,NUMNP,MEQ)

implicit none
!-pass
integer          :: ndof
integer          :: meq, numnp
double precision :: s, deld, fbc
integer          :: ibond, ia(ndof, numnp)
dimension S(MEQ,*),DELD(NDOF,*),IBOND(NDOF,*),FBC(MEQ)
!-locl
integer n,i,m,j,geq,heq
integer, external :: ieq

do n=1,NUMNP
    do i=1,NDOF
        if (IBOND(i,n).eq.1 .or. &
            IBOND(i,n).eq.2 .or. &
            IBOND(i,n).eq.5) then

            geq = IA(i,n)
            FBC(geq) = 0D0
            do m=1,NUMNP
                do j=1,NDOF
                    heq = ieq(j,m,NDOF)
                    FBC(geq) = FBC(geq) + S(geq,heq)*DELD(j,m)
                enddo
            enddo
        endif
    enddo
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine PRTNEQBC (IBOND,IA,FBC,NDOF,NUMNP,NSTEP)
!
implicit none
!-pass
integer :: ibond, ndof, numnp
double precision :: fbc
integer :: ia(ndof, numnp)
dimension IBOND(NDOF,*),FBC(*)
integer :: nstep
!-locl
integer :: lutmp
save lutmp
integer n,i,geq
integer, external :: ieq, nextlu

!-init
data lutmp/-1/
!
if (lutmp.lt.0) then
    lutmp = nextlu(0)
    call openf(lutmp,'prtneqbc.tmp','unknown')
endif
!
do n=1,NUMNP

    do i=1,NDOF
        if (IBOND(i,n).eq.1 .or. &   ! displacement
            IBOND(i,n).eq.2 .or. &   ! velocity
            IBOND(i,n).eq.5) then    ! ?

            geq = IA(i,n)

            write(lutmp,1) NSTEP,i,n,FBC(geq)
 1            format(1x,I5,1x,'B(',I1,',',I5,') = ',1PG15.6)

        endif
    enddo
enddo
!
return
end
!-------------------------------------------------------------------------------
#ifdef SPARSE
 subroutine ViewMatrix (IFLAG)
!
! Opens a local X-Display and to display the nonzero matrix topology.
! Color coding remains a bit vague. The sequence {ViewerDrawGetDraw,
! DrawGetPopup, DrawScalePopup(popup,0D0,1D0)} was not successful to
! create a contour scale window.

use iomodule          
use petscksp
!

implicit none

#ifdef PETSC232
#include "petsc/finclude/petscviewer.h"
#else
#include "petsc/finclude/petscviewer.h"
#endif

!-pass
integer :: iflag
!-locl
integer :: ierr

PetscViewer XDisplay
!
call PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER, &
   'Nonzero matrix structure',0,0,300,300,XDisplay,ierr)

if (ierr.ne.0) then
    write(stderr,1) ierr
 1        format(1X,'ViewMatrix: ViewerDrawOpen error ',I4)
    call stoper()
endif

if (IFLAG.eq.0) then
!    call MatView(M_A,XDisplay,ierr)
else
!    call MatView(T_A,XDisplay,ierr)
endif

if (ierr.ne.0) then
    write(stderr,2) ierr
 2        format(1X,'ViewMatrix: MatView error ',I4)
    call stoper()
endif

! write(stdout,3,advance='no')
! 3   format(1x,'Continue > ')
! read(stdin,'(a1)') c

call PetscViewerDestroy(XDisplay,ierr)

if (ierr.ne.0) then
    write(stderr,4) ierr
 4        format(1X,'ViewMatrix: ViewerDisplay error ',I4)
    call stoper()
endif

return
end subroutine
#endif
! SPARSE

!-------------------------------------------------------------------------------
#ifdef SPARSE
   subroutine PrintVector(IFLAG)

   USE MODELDATAMODULE
   USE MODELDEFINITION
   USE AOMODULE, only: getallvectorentries, &
                        VecRestoreArrayF90, &
                        VecGetArrayF90
   USE MODELCTX
use spetscmodule, only: M_B, T_B

   implicit none




#include "petsc/finclude/petscvec.h"

   PetscErrorCode ierr
   PetscScalar, pointer   :: tmp(:)
   integer IFLAG,i

   if (IFLAG.eq.0) then

  call VecGetArrayF90(M_B,tmp,ierr)
  do i=1,NEQ
      modeldatactx%B(i) = tmp(i)
  enddo
  call VecRestoreArrayF90(M_B,tmp,ierr)

  call getallvectorentries(FILE_outputf,M_B,0)
   else

  call VecGetArrayF90(T_B,tmp,ierr)
  do i=1,NTEQ
      modeldatactx%BT(i) = tmp(i)
  enddo
  call VecRestoreArrayF90(T_B,tmp,ierr)

  call getallvectorentries(FILE_outputf,T_B,1)

   endif

   return
   end
#endif
! SPARSE

!-------------------------------------------------------------------------------
#ifdef SPARSE

subroutine FillLoadVector(IFLAG,ISTRIDE)

USE MODELDATAMODULE, only: modeldatactx
USE MODELDEFINITION
USE AOMODULE,       only: nequationst, &
                           nequations, &
                           ADD_VALUES
USE MODELCTX
use debugmodule, only: iecho, parallelLog
use spetscmodule, only: M_B, T_B

implicit none
PetscErrorCode ierr
integer IFLAG,i, ERROR, neqns
logical ISTRIDE
PetscInt, ALLOCATABLE :: indexarray(:)

if (IFLAG.eq.0 .or. IFLAG.eq.-1) then

    neqns = neq
    if (ISTRIDE) then
        neqns = nequations
    endif

    if (iecho.eq.8) then
        call ParallelLog("FillLoadVector", "setting M_B to 0")
    endif

    call VecZeroEntries(M_B,ierr)

    if (iecho.eq.8) then
        call ParallelLog("FillLoadVector", "Finished setting M_B to 0")
    endif

    ALLOCATE(indexarray(nequations),STAT=ERROR)
    do i=1,neqns
        indexarray(i) = i-1
    enddo

    if (IFLAG.eq.0) then
        call VecSetValuesLocal(M_B,neqns,indexarray, &
                               modeldatactx%B,ADD_VALUES,ierr)
    else
        call VecSetValuesLocal(M_B,neqns,indexarray, &
                               modeldatactx%BTOT,ADD_VALUES,ierr)
    endif

    if (iecho.eq.8) then
        call ParallelLog("FillLoadVector", "added modeldata variables")
    endif

    DEALLOCATE(indexarray)
else
    neqns = nteq

    if (ISTRIDE) then
        neqns = nequationst
    endif

    call VecZeroEntries(T_B,ierr)

    ALLOCATE(indexarray(nequationst),STAT=ERROR)
    do i=1,neqns
        indexarray(i) = i-1
    enddo

    call VecSetValuesLocal(T_B,neqns,indexarray, &
                           modeldatactx%BT,ADD_VALUES,ierr)
    DEALLOCATE(indexarray)
endif

return
end subroutine

#endif /* SPARSE */
!-------------------------------------------------------------------------------
#ifdef SPARSE

subroutine SyncLoadVector(IFLAG)

USE MODELDATAMODULE, only: modeldatactx
uSE MODELDEFINITION
USE AOMODULE
USE MODELCTX
use debugmodule, only: iecho, parallellog
use spetscmodule, only: M_B, T_B

implicit none

#include "petsc/finclude/petscvec.h"

PetscErrorCode         :: ierr
integer                :: IFLAG,i, ERROR
PetscScalar, pointer   :: tmp(:)

if (iecho.eq.8) then
    call ParallelLog("SyncLoadVector", "Entered subroutine")
endif

call FillLoadVector(IFLAG,.true.)


! The M_B vector is a PETSC vector distributed over the threads.
! It is created carrying only local space for local equations,
! see subroutine UserInitializeWD in tecdat.F.

! However, the B vector (used to be) global.
! Entries put in the B vector of equations from other threads
! were still only present in the processor that put it there, 
! and not for the other thread where it actually is about.

! This subroutine makes sure that the M_B entries of the PETSc
! vector M_B are distributed over the relevant threads.


if (IFLAG.eq.0 .or. IFLAG.eq.-1) then
    ! There is a mechanical part.

    if (iecho.eq.8) then
        call ParallelLog("SyncLoadVector", " Assembling PETSc load vector ")
    endif

    ! step 1: make sure that all local contributions to M_B are available 
    !         for all other processes. Each thread can throw data in M_B, using 

    call VecAssemblyBegin(M_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyBegin 4 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif

    call VecAssemblyEnd(M_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyEnd 4 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif

    if (iecho.eq.8) then
        call ParallelLog("SyncLoadVector", "Finished Assembling PETSc load vector")
    endif

    ! step 2: The M_B array cannot be acessed directly, because it is not just an array, 
    !         but a PETSc object.
    !         Put the M_B object data into a local tmp array, which CAN be accessed.
    !         And it has been 'get', PETSc can no longer do any operations on it, 
    !         until it has been restored (see below).

    call VecGetArrayF90(M_B,tmp,ierr)

    ! step 3: Update B and BTOT in the modeldata structure.
    !         This step puts data in it that other threads contributed to the load vector,
    !         but that concern the domain section of this thread.

    if (IFLAG.eq.0) then
        modeldatactx%B = 0d0
        do i=1,NEQ
            modeldatactx%B(i) = tmp(i)
        enddo
    else
        ! iflag -1
        modeldatactx%BTOT = 0d0
        do i=1,NEQ
            modeldatactx%BTOT(i) = tmp(i)
        enddo
    endif

    ! step 4: Close up the M_B PETSc object again, 
    !         So that PETSc can do operations on it once more.

    call VecRestoreArrayF90(M_B,tmp,ierr)


    if (iecho.eq.8) then
        call ParallelLog("SyncLoadVector", "Stuffed PETSc vector into PETSc")
    endif


else
    ! There is a thermal part.
    call VecAssemblyBegin(T_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyBegin 5 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif

    call VecAssemblyEnd(T_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyBegin 5 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif


    modeldatactx%BT = 0d0
    call VecGetArrayF90(T_B,tmp,ierr)
    do i=1,NTEQ
        modeldatactx%BT(i) = tmp(i)
    enddo
    call VecRestoreArrayF90(T_B,tmp,ierr)
endif


return
end
#endif /* SPARSE */
!-------------------------------------------------------------------------------
#ifdef SPARSE

subroutine BuildMat (IFLAG)

use debugmodule, only: iecho, &
                       parallelLog, xit
use spetscmodule, only: M_A, M_B, T_A, T_B

use petscksp,    only: MAT_FINAL_ASSEMBLY, &
                       MAT_SYMMETRIC, &
                       PETSC_TRUE
use iomodule          

implicit none
#include "petsc/finclude/petscmat.h"  /* For matrix types */
#include "petsc/finclude/petscsys.h"  /* for PETSC_TRUE */
!-locl
logical Mfirst,Tfirst
integer :: iflag, ierr
save Mfirst,Tfirst
!-init
data Mfirst,Tfirst /.true.,.true./

if (IFLAG.eq.0) then

    if (iecho.eq.8) then
        call ParallelLog("BuildMat", "Assembling mechanical matrix")
    endif

    ! we have a mechanical initialization
    call MatAssemblyBegin(M_A,MAT_FINAL_ASSEMBLY,ierr)
    if (ierr.ne.0) then
        write(stderr,10) ierr
10         format(1x,'BuildMat: MatAssemblyBegin error ',I5)
        call xit(1,"MatAssemblyBegin failed")

    endif

    call MatAssemblyEnd(M_A,MAT_FINAL_ASSEMBLY,ierr)
    if (ierr.ne.0) then
        write(stderr,20) ierr
20           format(1x,'BuildMat: MatAssemblyEnd error ',I5)
        call xit(1,"MatAssemblyEnd failed")
    endif

    if (iecho.eq.8) then
        call ParallelLog("BuildMat", "Assembling load vector")
    endif

    call FillLoadVector(IFLAG, .true.)
    
    call VecAssemblyBegin(M_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyBegin 6 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif

    call VecAssemblyEnd(M_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyEnd 6 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif

    if (iecho.eq.8) then
        call ParallelLog("BuildMat", "Assembling load vector completed")
    endif


    if (Mfirst) then
        call MatSetOption(M_A,MAT_SYMMETRIC,PETSC_TRUE,ierr)
        if (ierr.ne.0) then
            write(stderr,30) ierr
30          format(1x,'BuildMat: MatSetOption error ',I4)
            call xit(2,"MatSetOption failed")
        endif
        Mfirst = .false.
    endif

else

    call MatAssemblyBegin(T_A,MAT_FINAL_ASSEMBLY,ierr)
    if (ierr.ne.0) then
        write(stderr,10) ierr
        call stoper()
    endif

    call MatAssemblyEnd(T_A,MAT_FINAL_ASSEMBLY,ierr)
    if (ierr.ne.0) then
        write(stderr,20) ierr
        call stoper()
    endif

    call FillLoadVector(IFLAG, .true.)
    call VecAssemblyBegin(T_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyBegin 7 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif

    call VecAssemblyEnd(T_B,ierr)
    if (ierr.ne.0) then
        write(*,*) "VecAssemblyEnd 7 returned an error: ", ierr
        stop "Leaving GTecton..."
    endif

    call PrintVector(IFLAG)

    Tfirst = .false.
endif
 
return
end subroutine

#endif /* SPARSE */
!-------------------------------------------------------------------------------
subroutine PRINTID ()

USE MODELDEFINITION, only: numnp, numsn
use modeltopology,   only: NDOF
use modeldatamodule, only: modeldatactx
use iomodule,        only: luout

implicit none

integer :: i, n, nn
!
!-pass
!
nn = 0
do n=1,NUMNP
    do i=1,NDOF
        if (MOD(nn,50).eq.0) then
           write(luout(1),1)
 1           format(//'    NODE   DOF       NEQ')
        endif
        write(luout(1),2) n,i,modeldatactx%ID(i,n)
 2      format(I8,4X,I1,3X,I8)
        nn = nn + 1
    enddo
enddo

if (NUMSN.gt.0) then
    nn = 1
    do n=1,NUMNP
        do i=1,NDOF
            if (MOD(nn,50).eq.0 .or. (n.eq.1.and.i.eq.1)) then
                write(luout(1),3)
 3               format(//' SL NODE   DOF       NEQ')
            endif
            if (modeldatactx%IDX(i,n).ne.0) then
                write(luout(1),2) n,i,modeldatactx%IDX(i,n)
                nn = nn + 1
            endif
        enddo
    enddo
endif
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine PRTRF (STN0,ISELM,ISTIME,STN)

USE MODELDEFINITION, only: npre
use modeltopology,   only: nstr
use timestepmodule,  only: nstep

implicit none

!-pass
integer :: ISELM,ISTIME
double precision :: STN0,STN
dimension STN0(NSTR,*),ISELM(*),ISTIME(*),STN(NSTR,*)
!-local
integer :: i, j, n

if (NPRE.eq.0) then
    return
endif
!
do i=1,NPRE
    n = ISELM(i)
    if (NSTEP.eq.ISTIME(i)-1) then
        do j=1,NSTR
            STN(j,n) = STN(j,n) + STN0(j,i)
        enddo
    endif
enddo

return
end subroutine

!-------------------------------------------------------------------------------
subroutine INVMAT (A,IDL,NDIM)
!

use algebra, only: clear

! Subroutine to invert a symmetric matrix
!

implicit none
!-pass
! TODO verify, is IDL an int?
integer          :: NDIM,i,j, num
double precision :: A(NDIM,NDIM)
integer          :: IDL(NDIM)
!-locl
double precision :: AI(9),ALIN(45)
!
num = 1
!   copy A to ALIN
do j=1,NDIM
    do i=1,j
        ALIN(num) = A(i,j)
        num = num + 1
    enddo
enddo
!
! factorize ALIN - result remains in ALIN
call CLEAR (AI,NDIM,"AI")
call LUSOLVE (ALIN,AI,IDL,NDIM,.TRUE.,.FALSE.)

do i=1,NDIM
    call CLEAR (AI,NDIM,"AI")
    AI(i) = 1d0
!   invert matrix and put result into AI
    call LUSOLVE (ALIN,AI,IDL,NDIM,.FALSE.,.TRUE.)
    do j=1,NDIM
        A(j,i) = AI(j)
    enddo
enddo
!
return
end subroutine

!------------------------------------------------------------------------------
