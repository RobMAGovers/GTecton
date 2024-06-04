   program main
!************************************************************************
   USE TIMESTEPMODULE
   USE MATERIALSMODULE
   USE MODELTOPOLOGY
   USE MODELDEFINITION
   USE MERGERMODULE
   USE IOMODULE
   USE FILEMODULE
   USE MODELDATAMODULE
use debugmodule, only: xit
!************************************************************************
   implicit none
!************************************************************************
  
   type (modeldata) :: read_in, write_out
   type (matdata)   ::  mat_in,   mat_out

   integer ERROR
   integer nvmax, nemax
   integer t_nvlocal, t_nelocal, t_rank
   integer i, ireq,j, ierr,ifopen,k
   integer nsteps, maxtime
   integer, allocatable :: outputs(:)
   logical there, readout
!************************************************************************
   integer nextlu
   external nextlu

   logical :: MergeOnlyOneTimestep
   integer :: MergeOnlyTimestep

   logical :: DeleteOldFiles

! test stuff:
   integer :: timeID, procID


   luin = nextlu(33)

   ! process command line arguments
   call parseCommandLineArgs(MergeOnlyOneTimestep, MergeOnlyTimestep, DeleteOldFiles)

   nvglobal = 0
   neglobal = 0

!****** read partition.info ***

   inquire(file=PARTINFO, exist=there)

   if (there) then
  ierr = ifopen(luin,PARTINFO,'old')
  read(luin,*) nprocs

  ALLOCATE(nvlocal(nprocs), STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate nvlocal. Error ', ERROR
      call xit(1," ")
  endif

  ALLOCATE(nelocal(nprocs), STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate nelocal. Error ', ERROR
      call xit(1," ")
  endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
  ALLOCATE(nLocalElemsWithSLipperyNodesOfEachPartition(nprocs), STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate nLocalElemsWithSLipperyNodesOfEachPartition. Error ', ERROR
      call xit(1," ")
  endif
#endif


  do i=1,nprocs
      read(luin,*) t_rank, t_nvlocal, t_nelocal
      nvglobal = nvglobal + t_nvlocal
      neglobal = neglobal + t_nelocal
      nvlocal(t_rank+1) = t_nvlocal
      nelocal(t_rank+1) = t_nelocal
  enddo

  close(luin)
   else
!          write(stderr,101) partitioninfo
  write(stderr,101) PARTINFO
  call cleanup()
  stop
   endif

 101  format(1x,'Error: No file',1x,A,1x,'found, leaving application')

!*** finished reading partinfo

!     write(*,*) "nvlocal", nvlocal
!     write(*,*)	"nelocal", nelocal


!      nvmax = nvlocal(1)
!      nemax = nelocal(1)
!      do i=2,nprocs
!          nvmax = MAX(nvmax,nvlocal(i))
!          nemax = MAX(nemax,nelocal(i))
!      enddo

   nvmax = maxval(nvlocal)
   nemax = maxval(nelocal)

!      write(*,*) 'nvmax: ', nvmax
!      write(*,*) 'nemax: ', nemax


   ALLOCATE(gloInd(nvmax,nprocs), STAT=ERROR)
   if (ERROR.ne.0) then
  write(*,*) 'Could not allocate gloInd. Error ', ERROR
  call xit(1," ")
   endif

   ALLOCATE(gloElt(nemax,nprocs), STAT=ERROR)
   if (ERROR.ne.0) then
  write(*,*) 'Could not allocate gloElt. Error ', ERROR
  call xit(1," ")
   endif


   call readmodel() !from tecin.F

!      call printdata()

   maxtime = 0

   do i=1,NINTG
  maxtime = maxtime + timestepdatactx%MAXSTP(i)
   enddo

   mechrd = readout(outputcontroldatactx%IMPRINT,0,NMPRT, mode.eq.2.or.mode.eq.6)
   temprd = readout(outputcontroldatactx%ITPRINT,0,NTPRT, mode.eq.5.or.mode.eq.6)

   
!      write(*,*) 'first call'
   call initmodeldata(write_out, mat_out, nvglobal, neglobal,numat)
!      write(*,*) 'second call'
   call initmodeldata(read_in,mat_in,nvmax,nemax,numat)
!      write(*,*) 'done second'

   

!**** If there are multiple timesteps, merge all the time steps.

   if (MergeOnlyOneTimestep) then


  mechrd = readout(outputcontroldatactx%IMPRINT,MergeOnlyTimestep,NMPRT, mode.eq.2.or.mode.eq.6)
  temprd = readout(outputcontroldatactx%ITPRINT,MergeOnlyTimestep,NTPRT, mode.eq.5.or.mode.eq.6)

  if (mechrd.or.temprd) then

      write(*,*) 'Only merging timestep', MergeOnlyTimestep

      do j=1,nprocs

          call openfile(0, MergeOnlyTimestep, j-1, nprocs, maxtime)

          call readfemout(read_in, mat_in, write_out, mat_out, MergeOnlyTimestep,j)

          if (DeleteOldFiles) then
              write(*,*) "deleting partitioned outputfile"
              close(luin, status="delete")
          else
              close(luin)
          endif
      enddo

      call openfile(1, MergeOnlyTimestep, -1, 0, maxtime)
      call writefemout(write_out,mat_out,MergeOnlyTimestep)
      close(luout(2))
  else
      write(*,*) 'There is no data of timestep ', ireq
  endif


   else

  write(*,*) 'Merging timestep 0'


 if (mechrd.or.temprd) then

  do j=1,nprocs
      ! when reading binaries, the nodal point data is being
      ! in array blocks. The size of these arrays must be set to match
      ! the number of nodal points processed by this specific
      ! processor (all procs can (and likely will) have different
      ! amounts). nv = number of vertices; ne = number of elements
!              call initmodeldata(read_in, &
!                                 mat_in, &
!                                 nvlocal(j), &
!                                 nelocal(j), &
!                                 numat)

      luin = nextlu(44)

      call openfile(0, 0, j-1, nprocs, maxtime)
      ! j-1 because j is from 1 to nprocs, while file numbering
      ! is from 0 to nprocs-1. Inside, the pointer luin
      ! is computed

      call readfemout(read_in,  &
                      mat_in,  &
                      write_out,  &
                      mat_out, &
                      0, &
                      j)


      if (DeleteOldFiles) then
          close(luin, status="delete")
      else
          close(luin)
      endif

  enddo

  call openfile(1, 0, -1, 0, maxtime)
  call writefemout(write_out,mat_out,0)
  close(luout(2))
   endif







do ireq=1,maxtime

  mechrd = readout(outputcontroldatactx%IMPRINT,ireq,NMPRT, mode.eq.2.or.mode.eq.6)
  temprd = readout(outputcontroldatactx%ITPRINT,ireq,NTPRT, mode.eq.5.or.mode.eq.6)

  if (mechrd.or.temprd) then

      write(*,*) 'Merging timestep', ireq

      do j=1,nprocs

          call openfile(0, ireq, j-1, nprocs, maxtime)

          call readfemout(read_in, mat_in, write_out, mat_out, ireq,j)

          close(luin)

      enddo

      call openfile(1, ireq, -1, 0, maxtime)
      call writefemout(write_out,mat_out,ireq)
      close(luout(2))
  endif

enddo

   endif


   call freemodeldata(write_out, mat_out)
   call freemodeldata(read_in,mat_in)

   call cleanup()
   stop
   end

!********************************************************************


   logical function readout(printlist, record, nprt,lmode)
   implicit none
   integer :: printlist(*), record, nprt
   logical :: lmode
   integer :: iprt

   readout = .false.

   if (.not. lmode) goto 300

   if (record.eq.0) then
  readout = .true.
  goto 300
   endif

   if (nprt.eq.0) goto 300

   if (nprt.lt.0) then
  if (MOD(record,IABS(nprt)).eq.0) readout = .true.
  goto 300
   endif

   do iprt = 1,nprt
  if (printlist(iprt).eq.record) then
      readout = .true.
      goto 300
  endif
   enddo

 300  return

   end function readout

!********************************************************************

subroutine readmodel()

USE TIMESTEPMODULE
USE MATERIALSMODULE
USE MODELTOPOLOGY
USE MODELDEFINITION
USE MERGERMODULE
USE IOMODULE
use debugmodule,    only: xit, iecho
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use meshdatamodule,    only: meshdatactx
use modeldatamodule,   only: modeldatactx
#endif

implicit none

character(len=255) modelpartname
character(len=40) frmt
integer i,j,k,xx, ERROR
logical there

integer logsize, inneri, innerj, rnk


call openf (luin,'modeldata.dat','unknown')

#ifdef EXPERIMENTAL_ReducedSlipAlloc
    read(luin, '(2i12)') meshdatactx%nvglobal, meshdatactx%neglobal

if (iecho.eq.2) then
    write(*,*) "read number of vertices: ", meshdatactx%nvglobal
    write(*,*) "read number of elements: ", meshdatactx%neglobal
endif

#endif

#ifdef EXPERIMENTALTRANSIENTDISCREEP
read(luin, '(6i12)') NSD,NDOF,NSTR,NEN,NUMAT,NSTRTAY
#else
read(luin, '(5i12)') NSD,NDOF,NSTR,NEN, NUMAT
#endif

if (iecho.eq.2) then
    write(*,*) "read NSD:  ", NSD
    write(*,*) "read NDOF: ", NDOF
    write(*,*) "read NSTR: ", NSTR
    write(*,*) "read NEN:  ", NEN
    write(*,*) "read NUMAT:", NUMAT
#ifdef EXPERIMENTALTRANSIENTDISCREEP
    write(*,*) "read NSTRTAY: ", NSTRTAY
#endif
endif

read(luin, '(8i12)') MODE, ICVIS, FaultVelo, NSED, &
                     NSURFglobal, NUMFNglobal,NUMSLPglobal,NUMSNglobal

if (iecho.eq.2) then
    write(*,*) "read MODE:  ", MODE
    write(*,*) "read ICVIS: ", ICVIS
    write(*,*) "read FaultVelo: ", FaultVelo
    write(*,*) "read NSED:  ", NSED
    write(*,*) "read NSURFglobal:", NSURFglobal
    write(*,*) "read NUMFNglobal:", NUMFNglobal
    write(*,*) "read NUMSLPglobal:", NUMSLPglobal
    write(*,*) "read NUMSNglobal:", NUMSNglobal
endif


if (mode.ge.3) then
      read(luin, '(2i12)') IADVEC, IDIFT

    if (iecho.eq.2) then
        write(*,*) "read IADVEC: ", IADVEC
        write(*,*) "read IDIFT: ", IDIFT
    endif

endif

read(luin, '(I12)') NMPRT

if (iecho.eq.2) then
    write(*,*) "read NMPRT: ", NMPRT
endif


if (NMPRT.gt.0) then
      ALLOCATE(outputcontroldatactx%IMPRINT(MAX(1,NMPRT)), &
             STAT=ERROR)
      if (ERROR.ne.0) then
        write(*,*) 'Could not allocate IMPRINT. Error ', ERROR
           call xit(1," ")
      endif

     read(luin, 3) (outputcontroldatactx%IMPRINT(i),i=1,NMPRT)

    if (iecho.eq.2) then
        write(*,*) "read IMPRINT", outputcontroldatactx%IMPRINT
    endif

endif

if (mode.ge.3) then
    read(luin, '(I12)') NTPRT
    if (NTPRT.gt.0) then
          ALLOCATE(outputcontroldatactx%ITPRINT(MAX(1,NTPRT)), &
                 STAT=ERROR)
        if (ERROR.ne.0) then
            write(*,*) 'Could not allocate ITPRINT. Error ', ERROR
            call xit(1," ")
        endif

        read(luin, 3) (outputcontroldatactx%ITPRINT(i),i=1,NTPRT)
    endif
endif

  3   format (51i12)

read(luin, '(i12)') NINTG ! number of timestep groups

if (iecho.eq.2) then
    write(*,*) "read: ", NINTG, "time step groups"
endif


ALLOCATE(timestepdatactx%MAXSTP(MAX(1,NINTG)),STAT=ERROR)
if (ERROR.ne.0) then
    write(*,*) 'Could not allocate MAXSTP. Error ', ERROR
    call xit(1," ")
endif

! number of time steps in each group
read(luin, 3) (timestepdatactx%MAXSTP(i),i=1,NINTG)

if (iecho.eq.2) then
    write(*,*) "read: ", timestepdatactx%MAXSTP, "time steps in those groups"
endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
read(luin,*) nElemsWithSlipperyNodes
allocate(modeldatactx%SlipElemLookup(meshdatactx%neglobal))
read(luin,'(50i12)') modeldatactx%SlipElemLookup(1:meshdatactx%neglobal)

!write(*,*) "read slip lookup", modeldatactx%SlipElemLookup

#endif

LMAT = (NINTG.gt.0 .and. ICVIS.ne.0)
close(luin)



!      write(*,*) "completed reading"
!      write(stdout, '(5i5)') NSD,NDOF,NSTR,NEN, NUMAT
!      write(stdout, '(10i5)') MODE, ICVIS, IVELOUT, NSED,NSURF,
!     .        NSURFglobal, NUMFNglobal,NUMSLPglobal,NUMSNglobal
!
!      if (mode.ge.3) then
!          write(stdout, '(2i5)') IADVEC, IDIFT
!      endif
!      write(stdout, '(I5)') NMPRT
!      if (NMPRT.gt.0)
!     >  write(stdout, 3) (outputcontroldatactx%IMPRINT(i),i=1,NMPRT)
!      if (mode.ge.3) then
!          write(stdout, '(I5)') NTPRT
!          if (NTPRT.gt.0)
!     >      write(stdout, 3) (outputcontroldatactx%ITPRINT(i),i=1,NTPRT)
!      endif
!      write(stdout, '(i5)') NINTG
!      write(stdout, 3) (timestepdatactx%MAXSTP(i),i=1,NINTG)


   ALLOCATE(IEN(NEN,neglobal), STAT=ERROR)
   if (ERROR.ne.0) then
  write(*,*) 'Could not allocate IEN. Error ', ERROR
  call xit(1," ")
   endif


   if (NSURFglobal.gt.0) then
  ALLOCATE(ISURFORDER(NSURFGLOBAL,nprocs), STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate NSURFGLOBAL. Error ', ERROR
      call xit(1," ")
  endif

  ALLOCATE(NSURFS(nprocs), STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate NSURFS. Error ', ERROR
      call xit(1," ")
  endif
   endif

   write(frmt,4) NEN
 4    format('(',I1,'I12)')

   do i=1,nprocs


  modelpartname = 'modeldata.dat.'

  rnk = i-1

  if (nprocs.eq.1) then
    logsize = 1
  else
    logsize = log10(real(nprocs-1)) + 1
  endif

  ! when merging from only 1 processors, logsize = -infty
  ! fix this:
  if (logsize < 1) then
      logsize = 1
  endif

  do j=1,logsize
      inneri = 10**(logsize-j)
      innerj = rnk/inneri
      modelpartname = trim(modelpartname)//char(48+innerj)
      rnk = rnk - innerj * inneri
  enddo

!          modelpartname = trim(modelpartname)


  inquire(file=modelpartname, exist=there)

  ! set maximum memory MTOT from defined value MAXA included from size.i
  if (.not.there) then
      write(*,*) 'file: ', trim(modelpartname), ' not found.; exiting...'
      call xit(1," ")
  else 
!              write(*,*) 'file: ', trim(modelpartname), ' found.'
  endif

  call openf(luin,modelpartname,'old')

  read(luin, '(I12)')  xx
  read(luin, 3) (gloInd(j,i),j=1,xx)
  read(luin, '(I12)')  xx
  read(luin, 3) (gloElt(j,i),j=1,xx)

  do j=1,xx
      read(luin,frmt) (IEN(k,gloElt(j,i)),k=1,NEN)
  enddo

  if (NSURFglobal.gt.0) then
      read(luin, '(I12)') NSURFS(i)
      read(luin, 3) (ISURFORDER(j,i),j=1,NSURFS(i))
  endif


#ifdef EXPERIMENTAL_ReducedSlipAlloc
    read(luin, *) nLocalElemsWithSLipperyNodesOfEachPartition(i)
! nLocalElemsWithSLipperyNodes
#endif


  close(luin)
   enddo

   return
   end

!********************************************************************

   subroutine initmodeldata(modeldat, matdat,numnp, numel, idummy)

   USE MERGERMODULE
   USE MODELDATAMODULE
   USE MODELTOPOLOGY
   USE MATERIALSMODULE
   USE MODELDEFINITION, only : numat, nsurfglobal
use debugmodule, only: xit
#ifdef EXPERIMENTAL_ReducedSlipAlloc
use modeldefinition, only: nElemsWithSlipperyNodes
#endif

   implicit none

   type (modeldata) :: modeldat
   type (matdata) :: matdat

    ! idummy is numat in the argument, which also comes from materials module
   integer ERROR, numnp, numel, idummy

!      write(*,*) 'initmodeldata says: numnp: ', numnp
!      write(*,*) 'initmodeldata says: numel: ', numel


   ! temp
   if (temprd) then

! displacements

  if (ALLOCATED(modeldat%D)) then
      DEALLOCATE(modeldat%D,    STAT=ERROR)
  endif

  ALLOCATE(modeldat%D(ndof,numnp),STAT=ERROR)

  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata D. Error ', ERROR
      call xit(1," ")
  endif

  modeldat%D = 0d0

! temperatures

  if (ALLOCATED(modeldat%T)) then
      DEALLOCATE(modeldat%T,    STAT=ERROR)
  endif

  ALLOCATE(modeldat%T(numnp),STAT=ERROR)

  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata T. Error ', ERROR
      call xit(1," ")
  endif

! something else

  if (ALLOCATED(modeldat%TX)) then
      DEALLOCATE(modeldat%TX,   STAT=ERROR)
  endif

  ALLOCATE(modeldat%TX(numnp),STAT=ERROR)

  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata TX. Error ', ERROR
      call xit(1," ")
  endif

! heat flux

  if (ALLOCATED(modeldat%FLUX)) then
      DEALLOCATE(modeldat%FLUX, STAT=ERROR)
  endif

  ALLOCATE(modeldat%FLUX(nsd,numel),STAT=ERROR)

  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata FLUX Error ', ERROR
      call xit(1," ")
  endif

   endif

   ! mech
   if (mechrd) then

!--- material data

  if (ALLOCATED(matdat%PRPMAT)) then
      DEALLOCATE(matdat%PRPMAT,STAT=ERROR)
  endif

  if (LMAT) then
      ALLOCATE(matdat%PRPMAT(6,numel), STAT=ERROR)
  else
      ALLOCATE(matdat%PRPMAT(6,numat), STAT=ERROR)
  endif

  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate matdata PRPMAT. Error ', ERROR
      call xit(1," ")
  endif


  if (ALLOCATED(matdat%PRPLAS)) then
      DEALLOCATE(matdat%PRPLAS,STAT=ERROR)
  endif

  if (LMAT) then
      ALLOCATE(matdat%PRPLAS(9,numel), STAT=ERROR)
  else
      ALLOCATE(matdat%PRPLAS(9,numat), STAT=ERROR)
  endif

  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate matdata PRPLAS. Error ', ERROR
      call xit(1," ")
  endif

!--- D (dispalcement)

  if (ALLOCATED(modeldat%D)) then
      DEALLOCATE(modeldat%D, STAT=ERROR)
  endif

  ALLOCATE(modeldat%D(ndof,numnp),STAT=ERROR)

  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata D. Error ', ERROR
      call xit(1," ")
  endif

  modeldat%D = 0d0
!--- DELD  (velocity)

  if (ALLOCATED(modeldat%DELD)) then
      DEALLOCATE(modeldat%DELD, STAT=ERROR)
  endif

  ALLOCATE(modeldat%DELD(ndof,numnp),STAT=ERROR)

  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata DELD. Error ', ERROR
      call xit(1," ")
  endif

  modeldat%DELD = 0d0

!--- DXE (differential displacement)

  if (ALLOCATED(modeldat%DXE)) then
      DEALLOCATE(modeldat%DXE, STAT=ERROR)
  endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
  ALLOCATE(modeldat%DXE(ndof,nen,nElemsWithSlipperyNodes),STAT=ERROR)
#else 
  ALLOCATE(modeldat%DXE(ndof,nen,numel),STAT=ERROR)
#endif

  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata DXE. Error ', ERROR
      call xit(1," ")
  endif

  modeldat%DXE = 0d0

!--- DELX (differential velocity)

  if (ALLOCATED(modeldat%DELX)) then
      DEALLOCATE(modeldat%DELX, STAT=ERROR)
  endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
  ALLOCATE(modeldat%DELX(ndof,nen,nElemsWithSlipperyNodes),STAT=ERROR)
#else
  ALLOCATE(modeldat%DELX(ndof,nen,numel),STAT=ERROR)
#endif

  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata DELX. Error ', ERROR
      call xit(1," ")
  endif
  modeldat%DELX = 0d0

!--- stress

  if (ALLOCATED(modeldat%STN)) then
      DEALLOCATE(modeldat%STN, STAT=ERROR)
  endif
  ALLOCATE(modeldat%STN(nstr,numel),STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata STN. Error ', ERROR
      call xit(1," ")
  endif
  modeldat%STN = 0d0


#ifdef EXPERIMENTALTRANSIENTDISCREEP
!--- taylor stress

  if (ALLOCATED(modeldat%STNTAY)) then
      DEALLOCATE(modeldat%STNTAY, STAT=ERROR)
  endif
  ALLOCATE(modeldat%STNTAY(nstrtay,numel),STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata STNTAY. Error ', ERROR
      call xit(1," ")
  endif
  modeldat%STNTAY = 0d0
#endif

!--- strain

  if (ALLOCATED(modeldat%STR)) then
      DEALLOCATE(modeldat%STR, STAT=ERROR)
  endif
  ALLOCATE(modeldat%STR(nstr,numel),STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata STR. Error ', ERROR
      call xit(1," ")
  endif
  modeldat%STR = 0d0

!--- strain rate

#ifndef EXPERIMENTAL_NORATESTORAGE
  if (ALLOCATED(modeldat%RATE)) then
      DEALLOCATE(modeldat%RATE, STAT=ERROR)
  endif
  ALLOCATE(modeldat%RATE(nstr,numel),STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata RATE. Error ', ERROR
      call xit(1," ")
  endif
#endif

!--- surface nodes

  if (ALLOCATED(modeldat%SURF)) then
      DEALLOCATE(modeldat%SURF, STAT=ERROR)
  endif
  ALLOCATE(modeldat%SURF(nsd,nsurfglobal),STAT=ERROR)
  if (ERROR.ne.0) then
      write(*,*) 'Could not allocate modeldata SURF. Error ', ERROR
      call xit(1," ")
  endif

   endif

   return
   end

!********************************************************************

   subroutine freemodeldata(modeldat, matdat)

   USE MERGERMODULE
   USE MODELDATAMODULE, only: modeldata
   USE MATERIALSMODULE, only: matdata

   implicit none

   type (modeldata) :: modeldat
   type (matdata) :: matdat

   integer ERROR

   ! temp
   if (temprd) then
  if (ALLOCATED(modeldat%D))      DEALLOCATE(modeldat%D,      STAT=ERROR)
  if (ALLOCATED(modeldat%T))      DEALLOCATE(modeldat%T,      STAT=ERROR)
  if (ALLOCATED(modeldat%TX))     DEALLOCATE(modeldat%TX,     STAT=ERROR)
  if (ALLOCATED(modeldat%FLUX))   DEALLOCATE(modeldat%FLUX,   STAT=ERROR)
   endif

   ! mech
   if (mechrd) then
  if (ALLOCATED(  matdat%PRPMAT)) DEALLOCATE(  matdat%PRPMAT, STAT=ERROR)
  if (ALLOCATED(  matdat%PRPLAS)) DEALLOCATE(  matdat%PRPLAS, STAT=ERROR)
  if (ALLOCATED(modeldat%D))      DEALLOCATE(modeldat%D,      STAT=ERROR)
  if (ALLOCATED(modeldat%DELD))   DEALLOCATE(modeldat%DELD,   STAT=ERROR)
  if (ALLOCATED(modeldat%DXE))    DEALLOCATE(modeldat%DXE,    STAT=ERROR)
  if (ALLOCATED(modeldat%DELX))   DEALLOCATE(modeldat%DELX,   STAT=ERROR)
  if (ALLOCATED(modeldat%STN))    DEALLOCATE(modeldat%STN,    STAT=ERROR)
  if (ALLOCATED(modeldat%STR))    DEALLOCATE(modeldat%STR,    STAT=ERROR)
#ifndef EXPERIMENTAL_NORATESTORAGE
  if (ALLOCATED(modeldat%RATE))   DEALLOCATE(modeldat%RATE,   STAT=ERROR)
#endif
#ifdef EXPERIMENTALTRANSIENTDISCREEP
  if (ALLOCATED(modeldat%STNTAY)) DEALLOCATE(modeldat%STNTAY,   STAT=ERROR)
#endif
  if (ALLOCATED(modeldat%SURF))   DEALLOCATE(modeldat%SURF,   STAT=ERROR)
   endif

   return
   end

!********************************************************************

   subroutine readandwrite(read_in, mat_in, write_out, mat_out, nrank, nsize)

   USE IOMODULE
   USE MODELDATAMODULE, only: modeldata
   USE MATERIALSMODULE, only: matdata
   USE MERGERMODULE

   implicit none

   type (modeldata) :: read_in, write_out
   type (matdata)   ::  mat_in,   mat_out

   integer nrank, nsize

   if (nrank.lt.0 .or. nrank.gt.nsize) then
  write(stdout,31) nrank, nsize-1
 31       format(1x,'rank',1x,I5,1x,'not in range [0,', &
   I5)
  return
   endif



   return
   end

!********************************************************************

   subroutine openfile(outfile, nstep, nrank, nsize, maxitime)

   USE IOMODULE

!     generic routine to open files for reading, writing, with or
!     without processor ids
!     outfile: flag for target file
!            0=file is input, i.e. parallel data
!            1=file is output
!     nstep: time step
!     nrank: rank of input file (only with outfile=0)
!     nsize: number if separate files per time step
!     maxitime: maximum time

   implicit none

   integer  :: outfile
   integer  :: nstep, nsize, nrank, maxitime

   integer  :: i,ncntr,jlocal,inner,jnnerj, inneri,innerj
   integer  :: rnk, logsize, ierr, lnblk, ifopen, ifbopen
   logical  :: isopen, isthere
   external :: lnblk

   integer  :: nextlu
   external :: nextlu

   if (MAXITIME.ge.1) then
  ncntr = log10(dble(MAXITIME)) + 1
   else
  ncntr = 1
   endif

   FEDSK = FEDSK(1:FEDSKRT)//'.'
   jlocal = nstep

   do inner=1,ncntr
  inneri = 10**(ncntr-inner)
  jnnerj = jlocal/inneri
  FEDSK = FEDSK(1:FEDSKRT+inner)//char(48 + jnnerj)
  jlocal = jlocal - jnnerj * inneri
   enddo

!      write(*,*) 'aaa', outfile, nrank, nsize, FEDSK

   if (outfile.eq.1) then

!          write(*,*) 'opening output: ', FEDSK

  inquire(unit=luout(2),opened=isopen)

  if (isopen) close(luout(2))

  if (ascin) then
!              write (*,*) 'Opening FEDSK output file as text file: ', FEDSK
      call openf (luout(2),FEDSK(1:lnblk(FEDSK)),'unknown')
  else
!              write (*,*) 'Opening FEDSK output file as binary file: ', FEDSK
      call fbopen (luout(2),FEDSK(1:lnblk(FEDSK)),'unknown')
  endif

  write(stdout, 31) fedsk(1:lnblk(fedsk))
 31       format(1x,'Output file',1x,A,1x,'opened for writing')

   else if (outfile.eq.0 .and. nrank.ge.0 .and. nsize.gt.0) then

  ! open parallel input files

  FEDSK = FEDSK(1:FEDSKRT+ncntr+1)//'.'
  rnk = nrank

  logsize = log10(dble(nsize)) + 1

  do i=1,logsize
      inneri = 10**(logsize-i)
      innerj = rnk/inneri
      FEDSK = FEDSK(1:FEDSKRT+ncntr+1+i)//char(48+innerj)
      rnk = rnk - innerj * inneri
  enddo

  inquire(unit=luin,opened=isopen)

  if (isopen) then
!              write (*,*) 'file pointer is occupied', luin
!              write (*,*) 'closing luin: ', luin
      close(luin)
!              write (*,*) 'closed luin'
  else
!              write (*,*) 'file pointer is free', luin
  endif

!          write(*,*) 'looking for file: ', FEDSK

  inquire(file=FEDSK(1:lnblk(FEDSK)), exist=isthere)

  if (.not. isthere) then

      write(stderr, 41) FEDSK(1:lnblk(FEDSK))
 41           format(1x, 'Input file with name',1x,A,1x,'does not',1x, &
        'exist. Leaving application')

      stop
  else
!              write(*,*) 'YAY! file: ', trim(FEDSK), 'found!'
  endif

!          write(*,*) 'opening file: ', FEDSK(1:lnblk(FEDSK))

  if (ascin) then
      ierr = ifopen(luin,FEDSK(1:lnblk(FEDSK)),'old')
      if (ierr.ne.0) return
  else
      ierr = ifbopen(luin,FEDSK(1:lnblk(FEDSK)),'old')
      if (ierr.ne.0) return
  endif

!          write(stdout, 32) fedsk(1:lnblk(fedsk))
! 32       format(5x,'Input file',1x,A,1x,'opened for reading')

   else
  write(*,*) 'this should not happen; contact model support.'
   endif


   return
   end


!********************************************************************

   subroutine parseCommandLineArgs(MergeOnlyOneTimestep, MergeOnlyTimestep, DeleteOldFiles)

! interpret command line arguments.

   USE IOMODULE
   USE FILEMODULE
   use debugmodule, only: iecho

   implicit none

   character(len=255) arg
   integer narg,iarg,n,m
   integer nextlu, lnblk,iflu

   logical :: MergeOnlyOneTimestep
   integer :: MergeOnlyTimestep
   logical :: DeleteOldFiles

   ! set defaults
   narg = command_argument_count()
   FEIN      = 'TECIN.DAT'
   TEMPIN    = 'TEMPIN.DAT'
   FEDSK     = 'FEDSK.DAT'
   FEDSKRT   = lnblk(FEDSK)
   PARTINFO  = 'partition.info'
   stderr    = iflu('stderr')
   stdin     = iflu('stdin')
   stdout    = iflu('stdout')
   luin      = nextlu(10)
   luout(2)  = nextlu(luin+1)
   iecho     = 0
   ascin     = .false.

   MergeOnlyOneTimestep = .false.
   MergeOnlyTimestep    = -1
   DeleteOldFiles       = .false.

   m = 0

   ! change those defaults based on command line arguments

   if (narg.ge.1) then
  do iarg=1,narg
  call get_command_argument(iarg,arg)

  if (arg(1:5).eq.'fedsk') then
      n = INDEX(arg,'=')
      if (n.gt.0) then
          m = lnblk(arg)
          n=n+1
          if (n.gt.0 .and. m.gt.n) then
              FEDSK=arg(n:m)
              FEDSKRT = lnblk(FEDSK)
              if(FEDSK(FEDSKRT:FEDSKRT).eq.'/') then
                  fedsk = fedsk(1:fedskrt)//'FEDSK.DAT'
                  fedskrt = lnblk(fedsk)
              endif
          else
              write(stderr, 30)
  30                  format(1x,'fedsk argument error')
          endif
      else
          write(stderr,30)
      endif

  elseif (arg(1:8).eq.'partinfo') then
      n = INDEX(arg,'=')
      if (n.gt.0) then
          m = lnblk(arg)
          n = n+1
          if (n.gt.0 .and. m.gt.n) then
              PARTINFO=arg(n:m)
!                      write (0,*) 'partinfo from arg: ', PARTINFO 
          else
              write(stderr,26)
26                 format(1x,'partinfo argument error')
              return
          endif
      else
          write(stderr,26)
          return !call xit(1," ")
      endif

  elseif (arg(1:4).eq.'echo') then
      n = INDEX(arg,'=')
      if (n.gt.0) then
          if (n.gt.0 .and. m.gt.n) then
              read( arg(6:6), '(i1)' )  iecho
              if (iecho.gt.0) then
                  write(*,*) 'echo set to ', iecho
              endif
          else
              write(stderr,27)
27                 format(1x,'echo argument error')
              return
          endif

      endif

  elseif (arg(1:8).eq.'workpath') then
      n = INDEX(arg,'=')
      if (n.gt.0) then
          m = lnblk(arg)
          n = n+1
          if (n.gt.0 .and. m.gt.n) then
              WORKPATH=arg(n:m)
!                      write (0,*) 'workpath from arg: ', WORKPATH
          else
              write(stderr,28)
28                 format(1x,'workpath argument error')
              return
          endif
      else
          write(stderr,28)
          return!call xit(1," ")
      endif

  elseif (arg(1:6).eq.'tempin') then
      n = INDEX(arg,'=')
      if (n.gt.0) then
          m = lnblk(arg)
          n = n+1
          if (n.gt.0 .and. m.gt.n) then
              TEMPIN=arg(n:m)
          else
              write(stderr,25)
25                 format(1x,'tempin argument error')
              return !call xit(1," ")
          endif
      else
          write(stderr,25)
          return !call xit(1," ")
      endif

  elseif (arg(1:8).eq.'timestep') then
      n = INDEX(arg,'=')
      if (n.gt.0) then
          m = lnblk(arg)
          n = n+1 ! +1 to skip the '=' sign
          if (n.gt.0 .and. m.ge.n) then
              read(arg(n:m), '(i10)' ) MergeOnlyTimestep
              MergeOnlyOneTimestep = .true.
          else
              write(stderr,23)
23                 format(1x,'tempin argument error')
              return !call xit(1," ")
          endif
      else 
          write(stderr,23)
          return !call xit(1," ")
      endif

  elseif (arg(1:2).eq.'bi') then
      ascin = .false.
  elseif (arg(1:2).eq.'as') then
      ascin = .true.
  elseif (arg(1:4).eq.'noas') then
      ascin = .false.
  elseif (arg(1:4).eq.'nobi') then
      ascin = .true.
  elseif (arg(1:5).eq.'clean') then
      DeleteOldFiles = .true.

  else
      write(stderr,31)
  31          format(1x,'unknown argument')
  endif

  enddo
   else
  write(stdout,32) fedsk
  ! assume FEDSK.DAT.XXX.YYY is the name of outputs
   endif

   write(stdout,32) trim(fedsk)
  32  format(1x,'Output collected from files with base "',A,'"')

   return
   end

!********************************************************************

   subroutine cleanup()

   USE IOMODULE

   logical isopen

   inquire(unit=luout(2),opened=isopen)
   if (isopen) close(luout(2))

   inquire(unit=luin,opened=isopen)
   if (isopen) close(luin)

   return
   end
