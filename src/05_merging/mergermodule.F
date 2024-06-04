MODULE MERGERMODULE

integer nprocs
integer, allocatable :: nvlocal(:)
integer, allocatable :: nelocal(:)
integer, allocatable :: gloInd(:,:)
integer, allocatable :: gloElt(:,:)
integer, allocatable :: IEN(:,:)
integer, allocatable :: ISURFORDER(:,:)
! do not confuse this with the NSURF from MODELDEFINITION
integer, allocatable :: NSURFS(:)

integer :: nvglobal, neglobal
logical :: temprd, mechrd
!****** functions ************************************
public  :: writefemout
public  :: readfemout
public  :: printdata

contains

!*****************************************************
subroutine printdata()
!*****************************************************
implicit none
!*****************************************************
integer i,j
!*****************************************************
!      write(*,*) 'n vertices local: ', (nvlocal(i),i=1,nprocs)
!      write(*,*) 'n elems    local: ', (nelocal(i),i=1,nprocs)
!      write(*,*) 'n vertices global:', nvglobal
!      write(*,*) 'n elems    global:', neglobal
!*****************************************************
end subroutine printdata


!*****************************************************
subroutine readfemout(read_in,  &
                      mat_in,  &
                      write_out,  &
                      mat_out, &
                      ireq, &
                      nrank)
!*****************************************************
! Read the fedsk files. 
! The fedsk files can contains various amounts of data,
! and they are read type by type
!
! This subroutine is called for every partition
!*****************************************************
USE MODELDATAMODULE
use meshdatamodule,  only: meshdatactx
USE MODELDEFINITION
USE MATERIALSMODULE
USE TIMESTEPMODULE
USE MODELTOPOLOGY
USE IOMODULE
use debugmodule,     only: iecho
use algebra,         only: idig  
use formatsmodule,   only: FMT_fedsk2D_diff_displacement, &
                           FMT_fedsk3D_diff_displacement, &
                           FMT_fedsk2D_displacement, &
                           FMT_fedsk3D_displacement, &
                           FMT_fedsk_stressStrainRate_nstr3, &  
                           FMT_fedsk_stressStrainRate_nstr4, &
                           FMT_fedsk_stressStrainRate_nstr5, &
                           FMT_fedsk_stressStrainRate_nstr6
#ifdef EXPERIMENTALTRANSIENTDISCREEP
use formatsmodule,   only:  FMT_fedsk_stressStrainRate_nstr1
#endif
use modelctx,        only: getrank


!*****************************************************
implicit none
!*****************************************************
type (modeldata) :: read_in, write_out
type (matdata)   ::  mat_in,   mat_out

double precision :: read_num
character(len=1024) :: frmt
character(len=80)   :: datatype ! to pass on to potential error
character(len=80)   :: readLineNrString
integer          :: io, nglob, nrank
integer          :: IREQ, n, m, i, j, jj, k
integer                :: nGlobalSequence

numnp = nvlocal(nrank)
numel = nelocal(nrank)
!*****************************************************

if ((.not.temprd).and.(.not.mechrd)) then
    write (*,*) 'Neither mechanical nor temperature ', &
                'data found. Exiting.'
    return
endif

write(*,101) nrank, ireq
101  format(6x,'Reading output: processor',1x,I5,1x, &
                 'and time step',1x,I10)

!**** read temperatures ******************************


if (temprd) then
!          write (*,*) 'reading main temp block now', ireq, iadvec

    if (iecho.eq.2) then 
        write(*,*) '***** Reading thermal data *****'
    endif


    !----------- only for higher timesteps than 1; no idea what this is ----------
  
    if (ireq.gt.0 .and. IADVEC.ne.0) then

        datatype = "advection something"

        read_in%D = 0d0

        if (ascin) then

!                  write (*,*) 'reading temp D from txt file'                  
!                  write(frmt,80) NDOF

            do n=1,NUMNP
                read(luin,fmt=frmt,err=1000,end=2000,iostat=io) &
                     j,(read_in%D(i,n),i=1,NDOF)
                if (j.ne.n) then
                    write(stderr,85) n
                    return
                endif
            enddo
        else
!                  write (*,*) 'reading temp from binary file', NUMNP
            read(luin,err=1000,end=2000,iostat=io) &
                ((read_in%D(i,n),i=1,NDOF),n=1,NUMNP)

!            if (io.ne.0) then
!                write(*,*) "Error while reading D ", io
!                stop "Exiting mergefiles..."
!            endif

        endif

        do n=1,NUMNP
            nglob = gloInd(n,nrank)
            do j=1,NDOF
                write_out%D(j,nglob) = read_in%D(j,n)
            enddo
        enddo

        if (iecho.eq.2) then
            write(stderr,90)
        endif

    endif

   !----------- for all timesteps timesteps, including 0 ----------


      read_in%T = 0d0

  !-------------------- temperatures -------------------
    datatype = "temperatures"

    if (ascin) then
!              write (*,*) 'reading temp from txt file'
        frmt='(I12,1X,E15.0)'

        do n=1,NUMNP
            read(luin,frmt,err=1000,end=2000,iostat=io) &
                 j,read_in%T(n)

            if (j.ne.n) then
                write(stderr,10) n
10                   format(1x,'femout: synchronization error in', &
               1x,'temperatures data block for element ',I12)
               return
           endif
       enddo
    else
!              write (*,*) 'reading temp from binary file'
        read(luin,err=1000,end=2000,iostat=io) &
            (read_in%T(n),n=1,NUMNP)
    endif

    do n=1,NUMNP
        nglob = gloInd(n,nrank)
        write_out%T(nglob) = read_in%T(n)
    enddo

    if (iecho.eq.2) then
        write(stderr,15)
15     format(1x,'temperatures have been read')
    endif

!-------------------- differential temperatures -------------------
! currently disabled.
! Do not remove. 
! Might be required later.

!datatype = "differential temperatures"
!          write (*,*) 'reading diff heat'


!  if (NUMFNglobal+NUMSNglobal.gt.0.and.IDIFT.eq.1) then


!              read_in%TX = 0d0!call clear (TX,NUMNP,"TX")
!      frmt='(I12,1X,F11.0)'

!      if (ascin) then
!                  write (*,*) 'reading diff temp from txt file'
!          read(luin,frmt,err=1000,end=2000,iostat=io) &
!           (j,read_in%TX(j),m=1,NUMNP)

!        write(*,*) "read diff temp", n, j, read_in%TX(:)


!      else
!                  write (*,*) 'reading diff temp from bin file'
!          read(luin,err=1000,end=2000,iostat=io) &
!           (read_in%TX(j),j=1,NUMNP)
!      endif

!      do n=1,NUMNP
!          nglob = gloInd(n,nrank)
!          write_out%TX(nglob) = read_in%TX(n)
!      enddo

!      if (iecho.eq.2) write(stderr,30)

!30         format(1x,'differential temperatures have been read')
!  endif


!-------------------- heat flux -------------------
    read_in%FLUX = 0d0

    datatype = "heat flux"

    if (ascin) then
!              write (*,*) 'reading heat flux from txt file'
        write(frmt,35) NSD
35      format('(I12,1X,',I1,'E15.0)')

        if (iecho.eq.2) then
            write(0,*) 'format for flux: ' , frmt
        endif

        do n=1,NUMEL
            read(luin,fmt=frmt,err=1000,end=2000,iostat=io) &
                  j,(read_in%FLUX(i,n),i=1,NSD)


            if (j.ne.n) then
                write(stderr,40) n
40              format(1x,'femout: synchronization error in', &
                    1x,'heat flux data block for element ',I12)
                return
            endif

        enddo
    else
        read(luin,err=1000,end=2000,iostat=io) &
             ((read_in%FLUX(i,n),i=1,NSD),n=1,NUMEL)
    endif



    do n=1,NUMEL
        nglob = gloElt(n,nrank)
        do j=1,NSD
            write_out%FLUX(j,nglob) = read_in%FLUX(j,n)
        enddo
    enddo

!-------------------- done reading thermal data -------------------



    if (iecho.eq.2) then
        write(stderr,45)
45      format(1x,'heat fluxes have been read')
    endif

endif



!****** Reading mechanical data *******************************
if (mechrd) then

!          write (0,*) 'reading main mech block, with IREQ', IREQ

    if (iecho.eq.2) then 
        write(*,*) '***** Reading mechanical data *****'
    endif

  

    if (IREQ.gt.0) then

        if (LMAT) then
            m = NUMEL
        else
            m = NUMAT
        endif

!****** Viscosities ********************

        datatype = "viscosities"

!              write (*,*) 'perhaps reading viscosities?'
!              write (*,*) 'IDIG 1: ', IDIG(ICVIS,4,1)

        if (IDIG(ICVIS,4,1).ne.0 .and. MODE.eq.6) then
            if (ascin) then
!                      write (*,*) 'reading viscosities from txt file'                      
                do n=1,m
                    read(luin,fmt=50,err=1000,end=2000,iostat=io) &
                         j,mat_in%PRPMAT(3,n),mat_in%PRPMAT(4,n)
50                  format(I12,1X,E15.0,1X,F4.0)

!                if (io.ne.0) then
!                    write(*,*) "Error ",io,"while reading Viscosities ", n
!                    stop "Exiting mergefiles..."
!                 endif

                    if (j.ne.n) then
                        write(stderr,55) n
55                      format(1x,'femout: synchronization error', &
                               1x,'in viscosity data block for material ', &
                               I12)
                        return
                    endif
                enddo
            else
                read(luin,err=1000,end=2000,iostat=io) &
                  (mat_in%PRPMAT(3,n),mat_in%PRPMAT(4,n),n=1,m)
            endif

            if (iecho.eq.2) then
                write(stderr,60)
60             format(1x,'viscosities have been read')
            endif

          ! in case .not.LMAT, all outputs have same data for
          ! PRPMAT, take only 0th

            if (LMAT) then
                do n=1,m
                    nglob=gloElt(n,nrank)
                    mat_out%PRPMAT(3,nglob) = mat_in%PRPMAT(3,n)
                    mat_out%PRPMAT(4,nglob) = mat_in%PRPMAT(4,n)
                enddo
            else
                if (nrank.eq.0) then
                    do n=1,m
                        mat_out%PRPMAT(3,n) = mat_in%PRPMAT(3,n)
                        mat_out%PRPMAT(4,n) = mat_in%PRPMAT(4,n)
                    enddo
                endif
            endif
        endif

!************** plasticity params ************************

        datatype = "plasticity parameters"

!              write (*,*) 'IDIG 2', IDIG(ICVIS,4,2)

        if (IDIG(ICVIS,4,2).ne.0) then

            if (ascin) then
!                      write (*,*) 'reading plastic params from txt file'
                do n=1,m
                    read(luin,fmt=65,err=1000,end=2000,iostat=io) &
                         j,mat_in%PRPLAS(3,n),mat_in%PRPLAS(4,n)
65                  format(I12,1X,E15.0,1X,E15.0)
                    if (j.ne.n) then
                        write(stderr,70) n
70                      format(1x,'femout: synchronization error', &
                               1x,'in plastic data block for material ', &
                               I12)
                        return
                    endif
                enddo
            else
!               write (*,*) 'reading plastic params from bin file'
                read(luin,err=1000,end=2000,iostat=io) &
                     (mat_in%PRPLAS(3,n),mat_in%PRPLAS(4,n),n=1,m)
            endif

            if (iecho.eq.2) then
                write(stderr,75)
75              format(1x,'plastic parameters have been read')
            endif

            if (LMAT) then
                do n=1,m
                    nglob=gloElt(n,nrank)
                    mat_out%PRPLAS(3,nglob) = mat_in%PRPLAS(3,n)
                    mat_out%PRPLAS(4,nglob) = mat_in%PRPLAS(4,n)
                enddo
            else
                if (nrank.eq.0) then
                    do n=1,m
                        mat_out%PRPLAS(3,n) = mat_in%PRPLAS(3,n)
                        mat_out%PRPLAS(4,n) = mat_in%PRPLAS(4,n)
                    enddo
                endif
            endif

        endif
    endif  ! (IREQ.gt.0)



!          if (IVELOUT.eq.0) then
!              call EQUATE (read_in%DELD,read_in%D,NDOF*NUMNP)
!          else
!              call CLEAR (read_in%DELD,NDOF*NUMNP)
!          endif

!          call CLEAR (read_in%D,NDOF*NUMNP)
!          read_in%D = 0d0

! *************** displacements *********************
    datatype = "displacements"
    if (ascin) then

        if (iecho.eq.2) then
            write (*,*) 'reading displacements from txt file', NUMNP, luin
        endif

        if (NDOF.eq.2) then
            frmt = FMT_fedsk2D_displacement
        endif

        if (NDOF.eq.3) then
            frmt = FMT_fedsk3D_displacement
        endif

        if (iecho.eq.2) then
            write(*,*) '***** Reading displacements *****'
            write(*,*) 'format: ', frmt
        endif

        do n=1,NUMNP

            read(luin,frmt,err=1000,end=2000,iostat=io) j,(read_in%D(i,n),i=1,NDOF)

            if (j.ne.n) then
                write(stderr,85) n
85              format(1x,'femout: synchronization error in', &
                       1x,'displacements data block for node ',I8)
                return
            endif
        enddo

    else
        read(luin,err=1000,end=2000,iostat=io) &
            ((read_in%D(i,n),i=1,NDOF),n=1,NUMNP)
    endif

    do n=1,NUMNP
        nglob = gloInd(n,nrank)
        do j=1,NDOF
            write_out%D(j,nglob) = read_in%D(j,n)
        enddo
    enddo


    if (iecho.eq.2) then
        write(stderr,90)
90      format(1x,'displacements have been read')
    endif

#ifdef EXPERIMENTAL_FRICTION

! *************** differential forces *********************

    if (NUMSNglobal.gt.0) then

        if (.not. allocated(read_in%DIFORC)) then
            ALLOCATE(read_in%DIFORC(ndof,sum(nvlocal)))
        endif
        if (.not. allocated(write_out%DIFORC)) then
            ALLOCATE(write_out%DIFORC(ndof,sum(nvlocal)))
        endif

        datatype = "differential forces"
        if (ascin) then

            if (iecho.eq.2) then
                write (*,*) 'reading diff forces from txt file', NUMNP, luin
            endif

            ! diff forces have the same format
            if (NDOF.eq.2) then
                frmt = FMT_fedsk2D_displacement
            endif

            if (NDOF.eq.3) then
                frmt = FMT_fedsk3D_displacement
            endif

            if (iecho.eq.2) then
                write(*,*) '***** Reading diff forces *****'
                write(*,*) 'format: ', frmt
            endif

            do n=1,NUMNP

                read(luin,frmt,err=1000,end=2000,iostat=io) &
                    j,(read_in%Diforc(i,n),i=1,NDOF)

                if (j.ne.n) then
                    write(stderr,86) n
86                  format(1x,'femout: synchronization error in', &
                           1x,'diff force data block for node ',I8)
                    return
                endif
            enddo

        else
            read(luin,err=1000,end=2000,iostat=io) &
                ((read_in%Diforc(i,n),i=1,NDOF),n=1,NUMNP)
        endif

        do n=1,NUMNP
            nglob = gloInd(n,nrank)
            do j=1,NDOF
                write_out%Diforc(j,nglob) = read_in%Diforc(j,n)
            enddo
        enddo

        if (iecho.eq.2) then
            write(stderr,91)
91          format(1x,'differential forces have been read')
        endif
    endif  ! NUMSNglobal.gt.0

#endif // EXPERIMENTAL_FRICTION


!**************** Velocities ****************************
    if (IREQ.gt.0) then
        datatype = "velocities"

        if (iecho.eq.2) then
            write(*,*) '***** Reading velocities *****'
        endif


        if (ascin) then
!                  write (*,*) 'reading velocities from txt file'
            do n=1,NUMNP
                read(luin,fmt=frmt,err=1000,end=2000,iostat=io) &
                     j,(read_in%DELD(i,n),i=1,NDOF)

!           if (io.ne.0) then
!                write(*,*) "Error ",io,"while reading velocities of point ", n
!                stop "Exiting mergefiles..."
!           endif



                if (j.ne.n) then
                    write(stderr,95) n
95                  format(1x,'femout: synchronization error', &
                          'in velocities data block for node ',I8)
                    return
                endif
            enddo
        else
!           write (*,*) 'reading velocities from bin file'
            read(luin,err=1000,end=2000,iostat=io) &
                ((read_in%DELD(i,n),i=1,NDOF),n=1,NUMNP)
        endif

        do n=1,NUMNP
            nglob = gloInd(n,nrank)
            do j=1,NDOF
                write_out%DELD(j,nglob) = read_in%DELD(j,n)
            enddo
        enddo

        if (iecho.eq.2) then
            write(stderr,96)
96         format(1x,'velocities have been read')
        endif

    endif

!*************** Differential displacement ****************

    if (NUMSNglobal.gt.0) then
        datatype = "differential displacement"

!              if (IVELOUT.eq.0) then
!                   call EQUATE(read_in%DELX,
!     .                 read_in%DXE,NDOF*NEN*NUMEL)
!              else
!                  call CLEAR (read_in%DELX,NDOF*NEN*NUMEL)
!              endif
!              call CLEAR (read_in%DXE,NDOF*NEN*NUMEL)

        if (iecho.eq.2) then
            write(*,*) '***** Differential displacements *****'
        endif

        if (ascin) then

            if (NDOF.eq.2) then
                frmt = FMT_fedsk2D_diff_displacement
            else
                frmt = FMT_fedsk3D_diff_displacement
            endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc

        ! we have to determine the number of elements with slippery nodes.
        ! We do not read in TECIN.DAT, which contains this information.

!         write(*,*) "diff displ", nLocalElemsWithSlipperyNodesOfEachPartition(nrank)

            do n=1,nLocalElemsWithSlipperyNodesOfEachPartition(nrank)
                do k=1,NEN
                    read(luin,fmt=frmt,err=1000,end=2000, iostat=io) &
                         j,jj,(read_in%DXE(i,k,n),i=1,NDOF)
!                write(*,*) "read diff displ", j,jj,(read_in%DXE(i,k,n),i=1,NDOF)
                enddo
            enddo

!        write(*,*) "filled DXE: ", read_in%DXE

        else

!          write (*,*) 'Reading diff displacement from bin'

            read(luin,err=1000,end=2000,iostat=io) &
                 (((read_in%DXE(i,j,n),i=1,NDOF),j=1,NEN), &
                n=1,nLocalElemsWithSlipperyNodesOfEachPartition(nrank))

#else

            do n=1,NUMEL
                do k=1,NEN
                    read(luin,fmt=frmt,err=1000,end=2000, iostat=io) &
                          j,jj,(read_in%DXE(i,k,n),i=1,NDOF)


!                  read(luin,*,err=1000,end=2000, iostat=io) &
!                            j,jj,(read_in%DXE(i,k,n),i=1,NDOF)
                    if (j.ne.n) then
                        write(stderr,100) n
  100                   format(1x,'femout: synchronization',1x, &
                               'error in slippery nodes data block', &
                               ' for node ',I8)
                        return
                    endif
                enddo
            enddo
        else

!          write (*,*) 'Reading diff displacement from bin'

            read(luin,err=1000,end=2000,iostat=io) &
                 (((read_in%DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
#endif

        endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc
!********* This is fragile. Test this well!
!        write(*,*) "input DXE", read_in%DXE
        datatype = "differential delta displacement; DXE"


        write(*,*) "read n slip local", nLocalElemsWithSlipperyNodesOfEachPartition


        do n=1,nLocalElemsWithSlipperyNodesOfEachPartition(nrank)
                ! walk through the entrie
!              nglob = gloElt(elementSlipperyID(n,50),nrank)
!              nGlobalSequence = elementSlipperyID(nglob,51)
            do k=1,NEN
                do i=1,NDOF
!                      write_out%DXE(i,k,nGlobalSequence) = &
!                       read_in%DXE(i,k,elementSlipperyID(n,52))
                      write_out%DXE(i,k,n) = &
                       read_in%DXE(i,k,n)
                enddo
            enddo
        enddo
#else
        do n=1,NUMEL
            nglob = gloElt(n,nrank)
            do k=1,NEN
                do i=1,NDOF
                    write_out%DXE(i,k,nglob) = &
                      read_in%DXE(i,k,n)
                enddo
            enddo
        enddo
#endif

        if (iecho.eq.2) then
            write(stderr,105)
  105       format(1x,'differential displacements have been read')
        endif

!************ Differential velocities *********************
        if (IREQ.gt.0) then
            datatype = "differential velocities"

            if (iecho.eq.2) then
                write(*,*) '***** Differential velocities *****'
            endif

           if (ascin) then

#ifdef EXPERIMENTAL_ReducedSlipAlloc

                do n=1,nLocalElemsWithSlipperyNodesOfEachPartition(nrank)
                    do k=1,NEN
                        read(luin,fmt=frmt,err=1000,end=2000, &
                             iostat=io) j,jj, &    
                              (read_in%DELX(i,k,n),i=1,NDOF)

                    enddo
                enddo
#else

!                      write (*,*) 'reading diff. vel. from txt file'
                do n=1,NUMEL
                    do k=1,NEN
                        read(luin,fmt=frmt,err=1000,end=2000, &
                             iostat=io) j,jj, &
                              (read_in%DELX(i,k,n),i=1,NDOF)

!           if (io.ne.0) then
!                write(*,*) "Error ",io,"while reading differential velocities of point ", k
!                write(*,*) "of element", n
!                stop "Exiting mergefiles..."
!           endif



                        if (j.ne.n) then
                            write(stderr,107) n
  107                               format(1x,'femout: synchronization', &
                             ' error in slippery velocities', &
                             ' for node ',I8)
                            return
                        endif
                    enddo
                enddo

#endif

            else ! if ascin
!                      write (*,*) 'reading diff. vel. from txt file'

#ifdef EXPERIMENTAL_ReducedSlipAlloc
                read(luin,err=1000,end=2000,iostat=io) &
                     (((read_in%DELX(i,j,n),i=1,NDOF),j=1,NEN) &
                     ,n=1,nLocalElemsWithSlipperyNodesOfEachPartition(nrank))
#else
                read(luin,err=1000,end=2000,iostat=io) &
                     (((read_in%DELX(i,j,n),i=1,NDOF),j=1,NEN) &
                      ,n=1,NUMEL)
#endif

!           if (io.ne.0) then
!                write(*,*) "Error ",io,"while reading differential displacements"
!                stop "Exiting mergefiles..."
!           endif

            endif

! insert the locally read data in its proper place in the global array.
! This is difficult in reduced cllocation mode, because:
! - The read LMX data is not directly dorrelated to the node sequence
! - same for the written data.
! So what we have to do is:
! 1: read from the input file
! 2: determine the local node sequence number from he slippery sequence
! 3: determine the global point number from the local number.
! 4: write the data

#ifdef EXPERIMENTAL_ReducedSlipAlloc
!********* This is fragile. Test this well!
          ! read all the element that contain slippery nodes,
          ! as that is the storage convention in the fedsk files.
            do n=1,nLocalElemsWithSLipperyNodesOfEachPartition(nrank)

              ! of this element, determine the global element 


              ! contain the slippery sequence of the global element.
!              write(*,*) "rank", getrank(), "has slipLookup: ", modeldatactx%SlipElemLookup
!                write(*,*) "rank", getrank(), "has gloElt ", gloElt


!              nglob = gloElt(elementSlipperyID(n,53),nrank)
!              write(*,*) "rank", getrank(), "found slippery element", n, "which has glob", nglob
!              nGlobalSequence = elementSlipperyID(nglob,54)
!              write(*,*) "rank", getrank(), "found global sequence", nGlobalSequence

              ! found proper index, transfer the element
!              do k=1,NEN
!                  do i=1,NDOF
!                      write_out%DELX(i,k,nGlobalSequence) = &
!                       read_in%DELX(i,k,elementSlipperyID(n,55))
!                  enddo
!              enddo

            enddo
#else
            do n=1,NUMEL
                nglob = gloElt(n,nrank)
                do k=1,NEN
                    do i=1,NDOF
                        write_out%DELX(i,k,nglob) = &
                          read_in%DELX(i,k,n)
                    enddo
                enddo
            enddo
#endif

            if (iecho.eq.2) then
                write(stderr,108)
  108           format(1x,'differential velocities have been read')
            endif
        endif
    endif  ! (NUMSNglobal.gt.0)

!************* Stress *************************
    datatype = "stress"
!          call clear (read_in%STN,NSTR*NUMEL)

    if (iecho.eq.2) then
        write(*,*) '***** Reading stresses *****'
        write(*,*) 'per line ', NSTR,' entries', NUMEL, meshdatactx%nelocal
    endif

    if (ascin) then
!              write (0,*) 'Read stress from txt'
!      write(frmt,130) NSTR
!  130         format('(I12,1X,',I1,'E15.0)')
        do n=1,NUMEL

    ! format adapts to long or short data format
            if (nstr.eq.3) then
                 read(luin,fmt=FMT_fedsk_stressStrainRate_nstr3,err=1000,end=2000,iostat=io) &
                     j,(read_in%STN(i,n),i=1,NSTR)

            else if (nstr.eq.4) then
                 read(luin,fmt=FMT_fedsk_stressStrainRate_nstr4,err=1000,end=2000,iostat=io) &
                     j,(read_in%STN(i,n),i=1,NSTR)

            else if (nstr.eq.5) then
!        write(*,*) "read nstr 5 entry", n
                 read(luin,fmt=FMT_fedsk_stressStrainRate_nstr5,err=1000,end=2000,iostat=io) &
                     j,(read_in%STN(i,n),i=1,NSTR)
!        write(*,*) "read", j, (read_in%STN(i,n),i=1,NSTR)

            else if (nstr.eq.6) then
                 read(luin,fmt=FMT_fedsk_stressStrainRate_nstr6,err=1000,end=2000,iostat=io) &
                     j,(read_in%STN(i,n),i=1,NSTR)
            else
                write(0,*) "NSTR", nstr, "not known. Should be 3, 4, 5 or 6. Please contact model support"
            endif

!            write(*,*) "rank", nrank, "read stress", j,(read_in%STN(i,n),i=1,NSTR)

            if (io.ne.0) then
                write(*,*) "Error ",io,"while reading stress of element ", n
                stop "Exiting mergefiles..."
            endif



!          if (iecho.eq.2) then
!              write(*,*) 'read stress ', j,(read_in%STN(i,n),i=1,NSTR)
!          endif

            if (j.ne.n) then
                write(stderr,135) n
  135           format(1x,'femout: synchronization error in', &
                      1x,'stresses data block for element ',I12)
                return
            endif
        enddo
    else
!              write (*,*) 'Read stress from bin'
        read(luin,err=1000,end=2000,iostat=io) &
             ((read_in%STN(i,n),i=1,NSTR),n=1,NUMEL)

!           if (io.ne.0) then
!                write(*,*) "Error ",io,"while reading stresses"
!                stop "Exiting mergefiles..."
!           endif

    endif

    do n=1,NUMEL
        nglob = gloElt(n,nrank)
        do i=1,NSTR
            write_out%STN(i,nglob) = read_in%STN(i,n)
            if (iecho.eq.2) then
                write(*,*) 'set output stress ', i, nglob, write_out%STN(i,nglob)
            endif
        enddo
    enddo

    if (iecho.eq.2) then
        write(stderr,140)
  140     format(1x,'stresses have been read')
    endif

#ifdef EXPERIMENTALTRANSIENTDISCREEP
!************* Taylor Stress *************************
    if (IREQ.gt.0) then ! for time steps larger than 0
    datatype = "taylor stress"
!          call clear (read_in%STNTAY,NSTRTAY*NUMEL)
        


    if (iecho.eq.2) then
        write(*,*) '***** Reading taylor stresses *****'
        write(*,*) 'per line ', NSTRTAY,' entries', NUMEL, meshdatactx%nelocal
    endif

    if (ascin) then
    ! acsii output
        do n=1,NUMEL

    ! format adapts to long or short data format
            if (nstrtay.eq.1) then
                 read(luin,fmt=FMT_fedsk_stressStrainRate_nstr1,err=1000,end=2000,iostat=io) &
                     j,(read_in%STNTAY(i,n),i=1,NSTRTAY)
            else if (nstrtay.eq.3) then
                 read(luin,fmt=FMT_fedsk_stressStrainRate_nstr3,err=1000,end=2000,iostat=io) &
                     j,(read_in%STNTAY(i,n),i=1,NSTRTAY)
            else if (nstrtay.eq.4) then
                 read(luin,fmt=FMT_fedsk_stressStrainRate_nstr4,err=1000,end=2000,iostat=io) &
                     j,(read_in%STNTAY(i,n),i=1,NSTRTAY)
            else if (nstrtay.eq.5) then
                 read(luin,fmt=FMT_fedsk_stressStrainRate_nstr5,err=1000,end=2000,iostat=io) &
                     j,(read_in%STNTAY(i,n),i=1,NSTRTAY)
            else if (nstrtay.eq.6) then
                 read(luin,fmt=FMT_fedsk_stressStrainRate_nstr6,err=1000,end=2000,iostat=io) &
                     j,(read_in%STNTAY(i,n),i=1,NSTRTAY)
            else
                write(0,*) "NSTRTAY", nstrtay, "not known. Should be 1, 3, 4, 5 or 6. Please contact model support"
            endif


            if (io.ne.0) then
                write(*,*) "Error ",io,"while reading taylor stress of element ", n
                stop "Exiting mergefiles..."
            endif





            if (j.ne.n) then
                write(stderr,136) n
  136           format(1x,'femout: synchronization error in', &
                      1x,'taylor stresses data block for element ',I12)
                return
            endif
        enddo
    else
    ! binary output
!              write (*,*) 'Read taylor stress from bin'
        read(luin,err=1000,end=2000,iostat=io) &
             ((read_in%STNTAY(i,n),i=1,NSTRTAY),n=1,NUMEL)



    endif

    do n=1,NUMEL
        nglob = gloElt(n,nrank)
        do i=1,NSTRTAY
            write_out%STNTAY(i,nglob) = read_in%STNTAY(i,n)
            if (iecho.eq.2) then
                write(*,*) 'set output taylor stress ', i, nglob, write_out%STNTAY(i,nglob)
            endif
        enddo
    enddo

    if (iecho.eq.2) then
        write(stderr,141)
  141     format(1x,'taylor stresses have been read')
    endif
    endif
#endif

!***************** Strain **********************

    datatype = "strain"

    if (iecho.eq.2) then
        write(*,*) '***** Reading strains *****'
        write(*,*) 'per line ', NSTR,' entries'
    endif


    if (ascin) then
!              write (0,*) 'Read strain from txt'
        do n=1,NUMEL


            ! format adapts to long or short data format                                  
            if (nstr.eq.3) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr3,err=1000,end=2000,iostat=io) &
                     j,(read_in%STR(i,n),i=1,NSTR)
            else if (nstr.eq.4) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr4,err=1000,end=2000,iostat=io) &
                     j,(read_in%STR(i,n),i=1,NSTR)
            else if (nstr.eq.5) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr5,err=1000,end=2000,iostat=io) &
                     j,(read_in%STR(i,n),i=1,NSTR)
            else if (nstr.eq.6) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr6,err=1000,end=2000,iostat=io) &
                     j,(read_in%STR(i,n),i=1,NSTR)
            else
                write(0,*) "NSTR", nstr, "not known. Should be 3, 4, 5 or 6. Please contact model support"
            endif

!          read(luin,fmt=frmt,err=1000,end=2000,iostat=io) j,(read_in%STR(i,n),i=1,NSTR)

!           if (io.ne.0) then
!                write(*,*) "Error ",io,"while reading strain of element ", n
!                stop "Exiting mergefiles..."
!           endif


            if (j.ne.n) then
                write(stderr,142) n
  142           format(1x,'femout: synchronization error in', &
                       1x,'strains data block for element ',I12)
                return
            endif
        enddo
    else
!              write (*,*) 'Read strains from bin'
        read(luin,err=1000,end=2000,iostat=io) &
            ((read_in%STR(i,n),i=1,NSTR),n=1,NUMEL)

!           if (io.ne.0) then
!                write(*,*) "Error ",io,"while reading strains of elements "
!                stop "Exiting mergefiles..."
!           endif


    endif

!          write (0,*) 'done reading strain'

!          don't scale the strain as this is only a number pushing
!          routine
!          do n=1,NUMEL
!              read_in%STR(3,n) = read_in%STR(3,n)*5D-1
!          enddo

    do n=1,NUMEL
        nglob = gloElt(n,nrank)
        do i=1,NSTR
            write_out%STR(i,nglob) = read_in%STR(i,n)
        enddo
    enddo


    if (iecho.eq.2) then
        write(stderr,143)
  143   format(1x,'strains have been read')
    endif



#ifndef EXPERIMENTAL_NORATESTORAGE

!write(*,*) "determining whether to read rates", ireq

!***************** Strain rates **********************
    if (IREQ.gt.0) then



        datatype = "strain rates"

        if (iecho.eq.2) then
            write(*,*) '***** Reading strain rates *****'
            write(*,*) 'per line ', NSTR,' entries'
        endif


        if (ascin) then

            write(frmt,145) NSTR
  145       format('(I12,1X,',I1,'E15.0)')

            do n=1,NUMEL


                ! format adapts to long or short data format                                  
                if (nstr.eq.3) then
                     read(luin,fmt=FMT_fedsk_stressStrainRate_nstr3,err=1000,end=2000,iostat=io) &
                         j,(read_in%RATE(i,n),i=1,NSTR)
                else if (nstr.eq.4) then
                     read(luin,fmt=FMT_fedsk_stressStrainRate_nstr4,err=1000,end=2000,iostat=io) &
                         j,(read_in%RATE(i,n),i=1,NSTR)
                else if (nstr.eq.5) then
                     read(luin,fmt=FMT_fedsk_stressStrainRate_nstr5,err=1000,end=2000,iostat=io) &
                         j,(read_in%RATE(i,n),i=1,NSTR)
                else if (nstr.eq.6) then
                     read(luin,fmt=FMT_fedsk_stressStrainRate_nstr6,err=1000,end=2000,iostat=io) &
                         j,(read_in%RATE(i,n),i=1,NSTR)
                else
                    write(0,*) "NSTR", nstr, "not known. Should be 3, 4, 5 or 6. Please contact model support"
                endif

!              read(luin,fmt=frmt,err=1000,end=2000,iostat=io)   j,(read_in%RATE(i,n),i=1,NSTR)

!           if (io.ne.0) then
!                write(*,*) "Error ",io,"while reading strain rate of element ", n, "."
!                stop "Exiting mergefiles..."
!           endif


!                      write(*,*) 'strain rate: ',j,(read_in%RATE(i,n),i=1,NSTR)

                if (j.ne.n) then
                    write(stderr,150) n
  150               format(1x,'femout: synchronization error', 1x,'in rates data block for element ',I8)
                    return
                endif
            enddo
        else
!                  write (*,*) 'Read strain rates from bin'
            read(luin,err=1000,end=2000,iostat=io)  ((read_in%RATE(i,n),i=1,NSTR),n=1,NUMEL)

!           if (io.ne.0) then
!                write(*,*) "Error ",io,"while reading strain rates of elements."
!                stop "Exiting mergefiles..."
!           endif


        endif

!               write (*,*) 'finished reading strain rates'

!              don't scale as this is only a number pushing routine
!              do n=1,NUMEL
!                  read_in%RATE(3,n) = read_in%RATE(3,n)*5D-1
!              enddo

        do n=1,NUMEL
            nglob = gloElt(n,nrank)
            do i=1,NSTR
                write_out%RATE(i,nglob) = read_in%RATE(i,n)
            enddo
        enddo

        if (iecho.eq.2) then
            write(stderr,155)
  155         format(1x,'total strain rates have been read')
        endif

!******************** Surface uplifts **********************

!!              if (NSURF.gt.0.and.(NSED.gt.0.or.IFILEVERSION.lt.2002))
!              ! fileversion always after 2002 for this code ...
!              if (NSURFglobal.gt.0.and.NSED.gt.0) then
!
!                  if (ascin) then
!                      write (*,*) 'Read surface uplifts from txt'
!                      write(frmt,160) NDOF
!  160                 format('(I8,1X,',I1,'E15.0)')
!                      do n=1,NSURFS(nrank)
!                          read(luin,fmt=frmt,err=1000,end=2000,
!     >                     iostat=io) j,(read_in%SURF(i,n),i=1,NSD)
!                          if (j.ne.n) then
!                              write(stderr,165) n
!  165                         format(1x,'femout: synchronization err',
!     >                          1x,'in block for surface node ',I8)
!                              return
!                          endif
!                      enddo
!                  else
!                      write (*,*) 'Read surface uplifts from bin'
!                      read(luin,err=1000,end=2000,iostat=io)
!     >      ((read_in%SURF(i,n),i=1,NSD),n=1,NSURFS(nrank))
!                  endif
!
!                  do n=1,NSURFS(nrank)
!                      nglob = ISURFORDER(n,nrank)
!                      do i=1,NSD
!                          write_out%SURF(i,nglob) = read_in%SURF(i,n)
!                      enddo
!                  enddo
!
!                  if (iecho.eq.2) write(stderr,170)
!  170             format(1x,'surface uplifts have been read')
!
!              endif ! NSURFglobal.gt.0.and.NSED.gt.0
    endif  ! IREQ > 0

#endif // EXPERIMENTAL_NORATESTORAGE

endif   ! mechrd == true



666   return


!.... ....................... R E A D   E R R O R S ............................

1000    write(stderr,*) "femout: Error occurred while trying to read ", trim(datatype), &
                        " of rank ", nrank, &
                        "of timestep ", IREQ, &
                        "error: ", io
if (trim(datatype).eq."displacements" .and. nrank.eq.1 .and. IREQ.eq.0 .and. io.eq.5010) then
    write(stderr,*) "------------------------------------------------------------------"
    write(stderr,*) "Possible cause: merging ascii output in binary mode or the other way around."
    write(stderr,*) "Mergefiles reads fedsk files default in binary mode."
    write(stderr,*) "Add option 'as' to overwrite and read in ascii mode."
    write(stderr,*) "------------------------------------------------------------------"
endif


        stop "Leaving mergefiles..."
return

2000    write(stderr,*) "femout: Unexpected EOF while trying to read ", trim(datatype), &
                        " of rank ", nrank, &
                        "of timestep ", IREQ, &
                        "error: ", io
        write(stderr,*) "------------------------------------------------------------------"
        write(stderr,*) "Possible reason: trying to read an ascii fedsk file in binary mode."
        write(stderr,*) "Mergefiles reads fedsk files default in binary mode."
        write(stderr,*) "Add option 'as' to overwrite and read in ascii mode."
        write(stderr,*) "------------------------------------------------------------------"

        stop "Leaving mergefiles..."
return


!*****************************************************
   end subroutine readfemout


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine writefemout(write_out, mat_out,ireq)

USE MODELDATAMODULE
USE MODELDEFINITION
USE MATERIALSMODULE
USE TIMESTEPMODULE
USE MODELTOPOLOGY
USE IOMODULE
use debugmodule, only: iecho
use algebra,        only: idig
use formatsmodule,   only: FMT_fedsk2D_diff_displacement, &
                           FMT_fedsk3D_diff_displacement, &
                           FMT_fedsk2D_displacement, &
                           FMT_fedsk3D_displacement, &
                           FMT_fedsk_stressStrainRate_nstr3, &
                           FMT_fedsk_stressStrainRate_nstr4, &
                           FMT_fedsk_stressStrainRate_nstr5, &
                           FMT_fedsk_stressStrainRate_nstr6
#ifdef EXPERIMENTALTRANSIENTDISCREEP
use formatsmodule,   only: FMT_fedsk_stressStrainRate_nstr1
#endif 

implicit none

type (modeldata) :: write_out
type (matdata) :: mat_out

character(len=40) frmt
integer IREQ, n, m, i, j, jj, k

integer :: slipperySequenceNr
integer :: iElem, iSlipElem, slipElemID


numnp = nvglobal
numel = neglobal

if ((.not.temprd).and.(.not.mechrd)) then
    return
endif

   write(*,101)  ireq
 101  format(2x,'WRITING OUTPUT: time step',1x,I10)
!      if (ireq.eq.0) then
!         if (ascin) then
!              write(luout(2),'(A50)',err=1000) fileversion
!          else
!              write(luout(2),err=1000) fileversion
!          endif
!      endif

!     write temperatures
if (temprd) then

    if (ireq.gt.0 .and. IADVEC.ne.0) then
        if (ascin) then
            write(frmt,80) NDOF
            do n=1,NUMNP
                write(luout(2),fmt=frmt,err=1000) &
                  n,(write_out%D(i,n),i=1,NDOF)
            enddo
        else
            write(luout(2),err=1000) &
              ((write_out%D(i,n),i=1,NDOF),n=1,NUMNP)
        endif
    endif

  if (ascin) then
      frmt='(I12,1X,1PG15.6)'
      do n=1,NUMNP
          write(luout(2),frmt,err=1000) &
           n,write_out%T(n)
      enddo
  else
      write(luout(2),err=1000) &
        (write_out%T(n),n=1,NUMNP)
  endif

!-------------------- differential temperatures -------------------
! currently disabled.
! Do not remove. 
! Might be required later.
!  if (NUMFNglobal+NUMSNglobal.gt.0.and.IDIFT.eq.1) then
!      frmt='(I12,1X,1PG15.6)'
!      if (ascin) then
!          write(luout(2),frmt,err=1000) &
!           (n,write_out%TX(n),n=1,NUMNP)
!      else
!          write(luout(2),err=1000) &
!           (write_out%TX(n),n=1,NUMNP)
!      endif
!  endif

  if (ascin) then
      write(frmt,35) NSD
35         format('(I12,1X,',I1,'(1PG15.6))')
      do n=1,NUMEL
          write(luout(2),fmt=frmt,err=1000) &
            n,(write_out%FLUX(i,n),i=1,NSD)
      enddo
  else
      write(luout(2),err=1000) &
           ((write_out%FLUX(i,n),i=1,NSD),n=1,NUMEL)
  endif
   endif

   if (mechrd) then

!          write(*,*) 'writing mechanical output'

  ! timesteps above 0
  if (IREQ.gt.0) then

      if (LMAT) then
          m = NUMEL
      else
          m = NUMAT
      endif

      ! powerlaw creep viscosity
      if (IDIG(ICVIS,4,1).ne.0 .and. MODE.eq.6) then
!                  write (*,*) 'powerlaw creep'
          if (ascin) then
              do n=1,m
                  write(luout(2),fmt=50,err=1000) &
          n,mat_out%PRPMAT(3,n),mat_out%PRPMAT(4,n)
50                  format(I12,1X,1PG15.6,1X,0PF4.2)
              enddo
          else
              write(luout(2),err=1000) &
             (mat_out%PRPMAT(3,n),mat_out%PRPMAT(4,n),n=1,m)
          endif
      endif

      ! visco plastic update switch
      if (IDIG(ICVIS,4,2).ne.0) then
!                  write (*,*) 'write  visco plastic update'
          if (ascin) then
              do n=1,m
                  write(luout(2),fmt=65,err=1000) &
            n,mat_out%PRPLAS(3,n),mat_out%PRPLAS(4,n)
65                     format(I12,1X,1PG15.6,1X,1PG15.6)
              enddo
          else
              write(luout(2),err=1000) &
              (mat_out%PRPLAS(3,n),mat_out%PRPLAS(4,n),n=1,m)
          endif
      endif
  endif

  !displacement
  if (ascin) then
!              write (*,*) 'write displacement'
!      write(frmt,80) NDOF
80         format('(I12,1X,',I1,'(1PG15.6))')
        if (NDOF.eq.2) then
            frmt = FMT_fedsk2D_displacement
        else
            frmt = FMT_fedsk3D_displacement
        endif

      do n=1,NUMNP
          write(luout(2),fmt=frmt,err=1000) &
           n,(write_out%D(i,n),i=1,NDOF)
      enddo
  else
      write(luout(2),err=1000) &
         ((write_out%D(i,n),i=1,NDOF),n=1,NUMNP)
  endif

#ifdef EXPERIMENTAL_FRICTION
if (NUMSNglobal.gt.0) then


  !differential forces
  if (ascin) then
!              write (*,*) 'write displacement'
!      write(frmt,80) NDOF
!80         format('(I12,1X,',I1,'(1PG15.6))')
        if (NDOF.eq.2) then
            frmt = FMT_fedsk2D_displacement
        else
            frmt = FMT_fedsk3D_displacement
        endif

      do n=1,NUMNP
          write(luout(2),fmt=frmt,err=1000) &
           n,(write_out%Diforc(i,n),i=1,NDOF)
      enddo
  else
      write(luout(2),err=1000) &
         ((write_out%Diforc(i,n),i=1,NDOF),n=1,NUMNP)
  endif
endif ! NUMSNglobal.gt.0
#endif


  ! write velocities
!          write(*,*) 'determine whether to read DELD ', ireq, ivelout
  if (IREQ.gt.0) then
!              write (*,*) 'write diff dispalcement DEL-D'
      if (ascin) then
            ! format same as for displacements.

          do n=1,NUMNP
              write(luout(2),fmt=frmt,err=1000) &
       n,(write_out%DELD(i,n),i=1,NDOF)
          enddo
      else
          write(luout(2),err=1000) &
          ((write_out%DELD(i,n),i=1,NDOF),n=1,NUMNP)
      endif
  endif

! write differential displacements
  if (NUMSNglobal.gt.0) then
      if (iecho.eq.2) then
          write (*,*) 'write slippery DXE'
      endif

      if (ascin) then


        if (NDOF.eq.2) then
            frmt = FMT_fedsk2D_diff_displacement
        else
            frmt = FMT_fedsk3D_diff_displacement
        endif

#ifdef EXPERIMENTAL_ReducedSlipAlloc

!          write(frmt,97) NDOF
97             format('(2(I12,1X),',I1,'(1PG15.6))')

!        write (*,*) "going to write DXE", write_out%DXE


!        write (*,*) "lookuptable: ", modeldatactx%SlipElemLookup

!        write (*,*) "NUMEL: ", NUMEL
        
        ! walk through the lookup table, and if we encounter 
        ! an element, we take it from the DXE
        iSlipElem = 0
        do iElem=1,NUMEL
            slipperySequenceNr = modeldatactx%SlipElemLookup(iElem)
            if (slipperySequenceNr.gt.0) then
                ! yay, we got one!
                iSlipElem = iSlipElem + 1
                do k=1,NEN
                    write(luout(2),fmt=frmt,err=1000) &
                     iElem, &
                     IEN(k,iElem), &
                     (write_out%DXE(i,k,iSlipElem),i=1,NDOF)
                enddo
            endif
        enddo

      else
          write(luout(2),err=1000) &
  (((write_out%DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,nElemsWithSlipperyNodes)
      endif


#else

!          write(frmt,97) NDOF
97             format('(2(I12,1X),',I1,'(1PG15.6))')

          do n=1,NUMEL
              do k=1,NEN
                  write(luout(2),fmt=frmt,err=1000) &
           n,IEN(k,n),(write_out%DXE(i,k,n),i=1,NDOF)
              enddo
          enddo
      else
          write(luout(2),err=1000) &
  (((write_out%DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
      endif

#endif

      if (IREQ.gt.0) then
          if (iecho.eq.2) then
              write (*,*) 'write slippery DEL-X'
          endif

            ! same format as for diff. displacements

#ifdef EXPERIMENTAL_ReducedSlipAlloc
          if (ascin) then

        iSlipElem = 0
        do iElem=1,NUMEL
            slipperySequenceNr = modeldatactx%SlipElemLookup(iElem)
            if (slipperySequenceNr.gt.0) then
                ! yay, we got one!
                iSlipElem = iSlipElem + 1

                  do k=1,NEN

                      write(luout(2),fmt=frmt,err=1000) &
                           iElem, &
                           IEN(k,iElem), &
                           (write_out%DELX(i,k,iSlipElem),i=1,NDOF)

                  enddo
               endif
              enddo
          else
              write(luout(2),err=1000) &
                 (((write_out%DELX(i,j,n),i=1,NDOF),j=1,NEN) &
                  ,n=1,nElemsWithSlipperyNodes)
          endif
#else    
          if (ascin) then
              do n=1,NUMEL
                  do k=1,NEN
                      write(luout(2),fmt=frmt,err=1000) &
                 n,IEN(k,n),(write_out%DELX(i,k,n),i=1,NDOF)
                  enddo
              enddo
          else
              write(luout(2),err=1000) &
                 (((write_out%DELX(i,j,n),i=1,NDOF),j=1,NEN) &
                  ,n=1,NUMEL)
          endif
#endif
      endif

  endif

! write stress
  if (ascin) then
!              write (*,*) 'write stress'
      write(frmt,130) NSTR
  130         format('(I12,1X,',I1,'(1PG15.6))')
      do n=1,NUMEL
!          write(*,*) n,(write_out%STN(i,n),i=1,NSTR)



         if (nstr.eq.3) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr3) &
             n,(write_out%STN(i,n),i=1,NSTR)
         else if (nstr.eq.4) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr4) &
             n,(write_out%STN(i,n),i=1,NSTR)
         else if (nstr.eq.5) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr5) &
             n,(write_out%STN(i,n),i=1,NSTR)
         else if (nstr.eq.6) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr6) &
             n,(write_out%STN(i,n),i=1,NSTR)
         else
           write(0,*) "NSTR", nstr, "not known. Should be 3, 4, 5 or 6. Please contact model support"

          endif

      enddo
  else
      write(luout(2),err=1000) &
         ((write_out%STN(i,n),i=1,NSTR),n=1,NUMEL)
  endif


#ifdef EXPERIMENTALTRANSIENTDISCREEP
! write taylor stress 
  if (IREQ.gt.0) then
       if (ascin) then
!              write (*,*) 'write taylor stress'
      write(frmt,131) NSTRTAY
  131         format('(I12,1X,',I1,'(1PG15.6))')
           do n=1,NUMEL




         if (nstrtay.eq.1) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr1) &
             n,(write_out%STNTAY(i,n),i=1,NSTRTAY)
         else if (nstrtay.eq.3) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr3) &
             n,(write_out%STNTAY(i,n),i=1,NSTRTAY)
         else if (nstrtay.eq.4) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr4) &
             n,(write_out%STNTAY(i,n),i=1,NSTRTAY)
         else if (nstrtay.eq.5) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr5) &
             n,(write_out%STNTAY(i,n),i=1,NSTRTAY)
         else if (nstrtay.eq.6) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr6) &
             n,(write_out%STNTAY(i,n),i=1,NSTRTAY)
         else
           write(0,*) "NSTRTAY", nstrTAY, "not known. Should be 1, 3, 4, 5 or 6. Please contact model support"

          endif

        enddo
      else
      write(luout(2),err=1000) &
         ((write_out%STNTAY(i,n),i=1,NSTRTAY),n=1,NUMEL)
  endif
  endif
#endif



! write strain
  if (ascin) then
!              write (*,*) 'write strain'
      do n=1,NUMEL
!          write(luout(2),fmt=frmt,err=1000) &
!             n,(write_out%STR(i,n),i=1,NSTR)


         if (nstr.eq.3) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr3) &
             n,(write_out%STR(i,n),i=1,NSTR)
         else if (nstr.eq.4) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr4) &
             n,(write_out%STR(i,n),i=1,NSTR)
         else if (nstr.eq.5) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr5) &
             n,(write_out%STR(i,n),i=1,NSTR)
         else if (nstr.eq.6) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr6) &
             n,(write_out%STR(i,n),i=1,NSTR)
         else
           write(0,*) "NSTR", nstr, "not known. Should be 3, 4, 5 or 6. Please contact model support"

         endif


      enddo
  else
      write(luout(2),err=1000) &
         ((write_out%STR(i,n),i=1,NSTR),n=1,NUMEL)
  endif


  if (IREQ.gt.0) then

#ifndef EXPERIMENTAL_NORATESTORAGE
      if (ascin) then
!                  write (*,*) 'write rates'

       do n=1,NUMEL


         if (nstr.eq.3) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr3) &
             n,(write_out%RATE(i,n),i=1,NSTR)
         else if (nstr.eq.4) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr4) &
             n,(write_out%RATE(i,n),i=1,NSTR)
         else if (nstr.eq.5) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr5) &
             n,(write_out%RATE(i,n),i=1,NSTR)
         else if (nstr.eq.6) then
          write(luout(2),fmt=FMT_fedsk_stressStrainRate_nstr6) &
             n,(write_out%RATE(i,n),i=1,NSTR)
         else
           write(0,*) "NSTR", nstr, "not known. Should be 3, 4, 5 or 6. Please contact model support"
         endif



        enddo

!          write(frmt,145) NSTR
!  145             format('(I12,1X,',I1,'(1PG15.6))')
!          do n=1,NUMEL
!              write(luout(2),fmt=frmt,err=1000) &
!        n,(write_out%RATE(i,n),i=1,NSTR)
!          enddo
      else
          write(luout(2),err=1000) &
      ((write_out%RATE(i,n),i=1,NSTR),n=1,NUMEL)
      endif
#endif

!              if (NSURF.gt.0.and.(NSED.gt.0.or.IFILEVERSION.lt.2002))
      ! fileversion always after 2002 for this code ...

      if (NSURFglobal.gt.0.and.NSED.gt.0) then
          if (ascin) then
!                      write (*,*) 'write surface nodes'
              write(frmt,160) NDOF
  160                 format('(I12,1X,',I1,'(1PG15.6))')
              do n=1,NSURFglobal
                  write(luout(2),fmt=frmt,err=1000) &
                 n,(write_out%SURF(i,n),i=1,NSD)
              enddo
          else
              write(luout(2),err=1000) &
    ((write_out%SURF(i,n),i=1,NSD),n=1,NSURFglobal)
          endif
      endif

  endif
   endif

   return

1000  write(stderr,1001) IREQ
 1001 format(1x,'femout: write error '/ &
 1x,'femout: Error occurred while trying to write data',1x, &
 'block ',I5)
stop "Leaving mergefiles"

   return

   end subroutine writefemout


   END MODULE MERGERMODULE
