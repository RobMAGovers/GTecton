!-------------------------------------------------------------------------------
subroutine extrme ()

use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use modeltopology,   only: NDOF, NSTR, NSD
use iomodule,        only: mout, tout
use timestepmodule,  only: nstep
use plotcontrol,     only: dispmn, dispmx, &
                           strsmn, strsmx, &
                           ratemn, ratemx, &
                           tempmn, tempmx, &
                           fluxmn, fluxmx

!  determines extrema in arrays

implicit none

if (Mout) then
    call vecxtr (modeldatactx%D,NDOF,meshdatactx%nvglobal,dispmn,dispmx,0)
    call vecxtr (modeldatactx%STN,NSTR,meshdatactx%neglobal,strsmn,strsmx,1)
    if (NSTEP.gt.0) then
#ifndef EXPERIMENTAL_NORATESTORAGE
        call vecxtr (modeldatactx%rate,NSTR,meshdatactx%neglobal,ratemn,ratemx,0)
#endif
    endif
endif

if (Tout) then
    call sclxtr (modeldatactx%T,meshdatactx%nvglobal,tempmn,tempmx)
    call vecxtr (modeldatactx%Thermalflux, & 
                 NSD,meshdatactx%neglobal,fluxmn,fluxmx,0)
endif

return
end
!-------------------------------------------------------------------------------
 subroutine fcp (Ireq,ierr)

use iomodule,        only: lufo
use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use materialsmodule, only: matdatactx
use modeldefinition, only: mode
use timestepmodule,  only: NMPRT, NTPRT, outputcontroldatactx
use errors,          only: ok
!
! Routine to copy a single data block
!
implicit none
!-pass
 integer ierr,Ireq



!-locl
 logical temprd,mechrd,RDTOO
!-init
 ierr = OK
!
!    write(*,*) 'fcp determines mechtemp with MODE', MODE


!    determine if mechanical and thermal blocks need to be read
 mechrd = RDTOO (outputcontroldatactx%imprint,IREQ,NMPRT, &
                 MODE.eq.2.or.MODE.eq.6)

 temprd = RDTOO (outputcontroldatactx%itprint,IREQ,NTPRT, &
                 MODE.eq.5.or.MODE.eq.6)

!    write(*,*) 'fcp says: mech/temp', mechrd, temprd

!
!       read file
#ifdef EXPERIMENTAL_NORATESTORAGE
call femout (modeldatactx%d, &
     modeldatactx%STN, &
     modeldatactx%STR, &
    modeldatactx%DXE,matdatactx%prpmat, &
    matdatactx%PRPLAS,&
     modeldatactx%T,modeldatactx%TX, modeldatactx%Thermalflux, &
    modeldatactx%SURF,modeldatactx%idsle,&
     modeldatactx%DELD, &
    modeldatactx%DX,IREQ,mechrd,temprd,ierr)
#else
!    read the output file
call femout (modeldatactx%d, &
     modeldatactx%STN, &
     modeldatactx%STR, &
    modeldatactx%DXE,matdatactx%prpmat, &
    matdatactx%PRPLAS,&
     modeldatactx%RATE,modeldatactx%T,modeldatactx%TX, modeldatactx%Thermalflux, &
    modeldatactx%SURF,modeldatactx%idsle,&
     modeldatactx%DELD, &
    modeldatactx%DX,IREQ,mechrd,temprd,ierr)
#endif
 if (ierr.ne.OK) return
!
!    write file
#ifndef EXPERIMENTAL_NORATESTORAGE
call fewrite (modeldatactx%D,&
              modeldatactx%STN,&
              modeldatactx%STR,modeldatactx%DXE,modeldatactx%idsle, &
   matdatactx%PRPMAT,meshdatactx%IEN,matdatactx%PRPLAS, &
   modeldatactx%rate,modeldatactx%T, modeldatactx%Thermalflux, &
   Ireq,mechrd,temprd,lufo,ierr)
#else
call fewrite (modeldatactx%D,&
              modeldatactx%STN,&
              modeldatactx%STR,modeldatactx%DXE,modeldatactx%idsle, &
   matdatactx%PRPMAT,meshdatactx%IEN,matdatactx%PRPLAS, &
   modeldatactx%T, modeldatactx%Thermalflux, &
   Ireq,mechrd,temprd,lufo,ierr) 
#endif


!
 return
 end
!-----------------------------------------------------------------------

#ifdef EXPERIMENTAL_NORATESTORAGE

subroutine femout(D,STN,STR,DXE,PRPMAT,PRPLAS,T,TX, &
    FLUX,SURF,IDSLE,DELD,DELDXE,IREQ,mechrd,temprd, &
    ierr)

#elseif EXPERIMENTALTRANSIENTDISCREEP

subroutine femout(D,STN,STNTAY,STR,DXE,PRPMAT,PRPLAS,T,TX, &
    FLUX,SURF,IDSLE,DELD,DELDXE,IREQ,mechrd,temprd, &
    ierr)
    

#else

subroutine femout(D,STN,STR,DXE,PRPMAT,PRPLAS,RATE,T,TX, &
    FLUX,SURF,IDSLE,DELD,DELDXE,IREQ,mechrd,temprd, &
    ierr)

#endif

use constants,       only: dfmax
USE ALGEBRA
USE MODELDEFINITION
uSE MESHDATAMODULE
USE MODELDATAMODULE
use formatsmodule,   only: FMT_fedsk2D_displacement, &  
                           FMT_fedsk2D_velocities, &
                           FMT_fedsk2D_diff_displacement, &
                           FMT_fedsk3D_displacement, &
                           FMT_fedsk3D_velocities, &
                           FMT_fedsk_stressStrainRate_nstr3, &
                           FMT_fedsk_stressStrainRate_nstr4, &
                           FMT_fedsk_stressStrainRate_nstr5, &
                           FMT_fedsk_stressStrainRate_nstr6, &
                           FMT_fedsk3D_diff_displacement
use debugmodule,     only: iecho
use errors,          only: exec, ok
use modeltopology,   only: NSD,NDOF,NSTR,NEN

use extrema,         only: stressPositiveSmallest,       stressPositiveBiggest, &
                           stressNegativeSmallest,       stressNegativeBiggest, &
                           strainPositiveSmallest,       strainPositiveBiggest, &
                           strainNegativeSmallest,       strainNegativeBiggest, &
                           strainRatePositiveSmallest,   strainratePositiveBiggest, &
                           strainRateNegativeSmallest,   strainRateNegativeBiggest, &
                           displacementPositiveSmallest, displacementPositiveBiggest, &
                           displacementNegativeSmallest, displacementNegativeBiggest
use timestepmodule,  only: NINTG, MAXIT, DELTP
use materialsmodule, only: LMAT
use iomodule,        only: ascin, stderr, luin
#ifdef EXPERIMENTALTRANSIENTDISCREEP
use formatsmodule,   only: FMT_fedsk_stressStrainRate_nstr1
use modeltopology,   only: NSTRTAY
#endif  

!
! routine to read FEDSK files
!

implicit none
!-pass
logical mechrd,temprd
integer IREQ

integer          :: idsle
double precision :: prpmat, t, prplas, surf, deld, tx, deldxe
#ifndef EXPERIMENTAL_NORATESTORAGE
double precision :: rate
#endif

double precision :: d(NDOF,NUMNP), FLUX(NSD,NUMEL)
double precision :: DXE(NDOF,NEN,NUMEL)
double precision :: STN(NSTR,NUMEL)
double precision :: STR(NSTR,NUMEL)
#ifdef EXPERIMENTALTRANSIENTDISCREEP
double precision :: STNTAY(NSTRTAY,NUMEL)
#endif

dimension PRPMAT(6,*), &
  T(*),PRPLAS(9,*),SURF(NSD,*), DELD(NDOF,*),&
    IDSLE(*),TX(*), &
    DELDXE(NDOF,NEN,*)

#ifndef EXPERIMENTAL_NORATESTORAGE
dimension   RATE(NSTR,*)
#endif

!-locl
character(len=255) frmt
character(len=255) datatype ! to pass the datatype on the a possible error
character(len=1024) dataRecord
double precision   :: effectiveSum,maxv
integer            :: iVertex, iElement, iDOF

integer            :: i, iElem, j, jj, k, m, n
integer            :: io
integer            :: ierr

!-init
ierr  = OK
!
stressPositiveSmallest = dfmax
stressNegativeSmallest = -dfmax
stressPositiveBiggest = 0d0
stressNegativeBiggest = 0d0

strainPositiveSmallest = dfmax
strainNegativeSmallest = -dfmax
strainPositiveBiggest = 0d0 
strainNegativeBiggest = 0d0

strainRatePositiveSmallest = dfmax
strainRateNegativeSmallest = -dfmax
strainRatePositiveBiggest = 0d0
strainRateNegativeBiggest = 0d0

displacementPositiveSmallest = dfmax
displacementNegativeSmallest = -dfmax
displacementPositiveBiggest = 0d0 
displacementNegativeBiggest = 0d0


if ((.not.temprd).and.(.not.mechrd)) then
    return
endif

!*****************************************************
!    read temperatures data block
!*****************************************************

if (iecho.eq.2) then
     write(stderr,*) 'Reading output data'
endif

if (temprd) then
    if (IREQ.gt.0 .and. IADVEC.ne.0) then
        call clear(D,NDOF*NUMNP,"D")
        if (ascin) then
            write(frmt,80) NDOF
80          format('(I12,1X,',I1,'E15.0)')
            do n=1,NUMNP
                read(luin,fmt=frmt,err=1000,end=2000,iostat=io) &
                 j,(D(i,n),i=1,NDOF)
                if (j.ne.n) then
                    write(stderr,85) n
                    ierr = EXEC
                    return
                endif
            enddo
        else  ! read as binary
            read(luin,err=1000,end=2000,iostat=io) &
               ((D(i,n),i=1,NDOF),n=1,NUMNP)
        endif
        if (iecho.eq.2) then
            write(stderr,90)
        endif
    endif
    call clear (T,NUMNP,"T")
    if (ascin) then
        frmt='(I12,1X,E15.0)'
        do n=1,NUMNP
            read(luin,frmt,err=1000,end=2000,iostat=io) &
           j,T(n)
            if (j.ne.n) then
                write(stderr,10) n
10              format(1x,'femout [a]: synchronization error in', &
                1x,'temperatures data block for element ',I12)
                ierr = EXEC
                return
            endif
        enddo
    else
        read(luin,err=1000,end=2000,iostat=io) (T(n),n=1,NUMNP)
    endif

    if (iecho.eq.2) then 
        write(stderr,15)
15        format(1x,'temperatures have been read')
    endif

    if (NUMFN+NUMSN.gt.0.and.IDIFT.eq.1) then

        call clear (TX,NUMNP,"TX")
        frmt='(I12,1X,F11.0)'
        if (ascin) then
            read(luin,frmt,err=1000,end=2000,iostat=io) &
               (j,TX(j),m=1,NUMNP)
        else
            read(luin,err=1000,end=2000,iostat=io) (TX(j),j=1, &
               NUMNP)
        endif

        if (iecho.eq.2) then 
            write(stderr,30)
30            format(1x,'differential temperatures have been read')
        endif
    endif

!************ read the thermal fluxes

    call clear (FLUX,NSD*NUMEL,"FLUX")

    if (ascin) then
        write(frmt,35) NSD
35      format('(I12,1X,',I1,'E15.0)')
        do n=1,NUMEL
            read(luin,fmt=frmt,err=1000,end=2000,iostat=io) j,(FLUX(i,n),i=1,NSD)
            if (iecho.eq.22) then
                write(*,*) 'read fluxes: ', n, (FLUX(1:NSD,n))
            endif
            if (j.ne.n) then
                write(stderr,40) n
40                   format(1x,'femout [b]: synchronization error in', &
                 1x,'heat flux data block for element ',I12)
                 ierr = EXEC
                 return
            endif
        enddo
    else
        read(luin,err=1000,end=2000,iostat=io) &
             ((FLUX(i,n),i=1,NSD),n=1,NUMEL)
    endif

    if (iecho.eq.2) then
        write(stderr,45)
45        format(1x,'heat fluxes have been read')
    endif

endif

!*****************************************************
!  Read mechanical output
!*****************************************************

if (mechrd) then

    !************ PRPMAT **********************
    if (IREQ.gt.0) then
        if (LMAT) then
            m = NUMEL
        else
            m = NUMAT
        endif

        if (IDIG(ICVIS,4,1).ne.0 .and. MODE.eq.6) then
            if (ascin) then
                do n=1,m
                    read(luin,fmt=50,err=1000,end=2000,iostat=io) &
                     j,PRPMAT(3,n),PRPMAT(4,n)
50                  format(I12,1X,E15.0,1X,F4.0)
                    if (j.ne.n) then
                        write(stderr,55) n
55                      format(1x,'femout [c]: synchronization error', &
                         1x,'in viscosity data block for material ', &
                         I12)
                        ierr = EXEC
                        return
                    endif
                enddo
            else
                read(luin,err=1000,end=2000,iostat=io) &
                (PRPMAT(3,n),PRPMAT(4,n),n=1,m)
            endif
            if (iecho.eq.2) write(stderr,60)
60          format(1x,'viscosities have been read')
        endif

        !************ PRPLAS **********************
        if (IDIG(ICVIS,4,2).ne.0) then
            if (ascin) then
                do n=1,m
                    read(luin,fmt=65,err=1000,end=2000,iostat=io) &
                     j,PRPLAS(3,n),PRPLAS(4,n)
65                  format(I12,1X,E15.0,1X,E15.0)
                    if (j.ne.n) then
                        write(stderr,70) n
70                      format(1x,'femout [d]: synchronization error', &
                         1x,'in plastic data block for material ', &
                         I12)
                        ierr = EXEC
                        return
                    endif
                enddo
            else
                read(luin,err=1000,end=2000,iostat=io) &
                (PRPLAS(3,n),PRPLAS(4,n),n=1,m)
            endif
            if (iecho.eq.2) write(stderr,75)
75          format(1x,'plastic parameters have been read')
        endif
    endif


    !****************** displacement ********************
    datatype = "displacements"
    call clear_matrix(modeldatactx%D, NDOF, NUMNP)
    if (ascin) then
        if (NDOF.eq.2) then
            frmt = FMT_fedsk2D_displacement
        else
            frmt = FMT_fedsk3D_displacement
        endif
        do n=1,NUMNP
            read(luin,fmt=frmt,err=1000,end=2000,iostat=io) &
              j,(modeldatactx%D(i,n),i=1,NDOF)
            if (j.ne.n) then
                write(stderr,85) n
85              format(1x,'femout [e]: synchronization error in', &
                 1x,'displacements data block for node ',I12)
                ierr = EXEC
                return
            endif
        enddo
    else
        read(luin,err=1000,end=2000,iostat=io) &
        ((modeldatactx%D(i,n),i=1,NDOF),n=1,NUMNP)
    endif
    if (iecho.eq.2) write(stderr,90)
90  format(1x,'displacements have been read')

    ! determine displacement extrema
    do iVertex=1,NUMNP
        effectiveSum = 0d0
        do iDOF=1,NDOF
            effectiveSum = effectiveSum + modeldatactx%D(iDOF,iVertex)**2
        enddo
        effectiveSum = sqrt(effectiveSum)
        if (effectiveSum.ge.0d0) then
            if (effectiveSum.gt.displacementPositiveBiggest) then
                displacementPositiveBiggest = effectiveSum
            endif
            if (effectiveSum.lt.displacementPositiveSmallest) then
                displacementPositiveSmallest = effectiveSum
            endif
        else if (effectiveSum.le.0d0) then
            if (effectiveSum.lt.displacementNegativeBiggest) then
                displacementNegativeBiggest = effectiveSum
            endif
            if (effectiveSum.gt.displacementNegativeSmallest) then
                displacementNegativeSmallest = effectiveSum
            endif
        endif
    enddo

#ifdef EXPERIMENTAL_FRICTION

    !******************** Differential forces  *******************
    if (numslpglobal.gt.0) then
        datatype = "diff forces"
        call clear_matrix(modeldatactx%Diforc, NDOF, NUMNP)
        if (ascin) then
            if (NDOF.eq.2) then
                frmt = FMT_fedsk2D_displacement
            else
                frmt = FMT_fedsk3D_displacement
            endif
        
            do n=1,NUMNP
                read(luin,fmt=frmt,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%Diforc(i,n),i=1,NDOF)
                if (j.ne.n) then
                    write(stderr,86) n
86                  format(1x,'femout [e]: synchronization error in', &
                     1x,'differential forces data block for node ',I12)
                    ierr = EXEC
                    return
                endif
            enddo
        else
            read(luin,err=1000,end=2000,iostat=io) &
             ((modeldatactx%Diforc(i,n),i=1,NDOF),n=1,NUMNP)
        endif
    endif
    if (iecho.eq.2) write(stderr,91)
91  format(1x,'differential forces have been read')

#endif

    call clear_matrix(modeldatactx%DELD, NDOF, NUMNP)
    if (IREQ.gt.0) then
        !******************** velocity  *******************
        datatype = "velocities"
        if (ascin) then
            do n=1,NUMNP
                if (NDOF.eq.2) then
                    read(luin,fmt=FMT_fedsk2D_velocities,err=1000,end=2000,iostat=io) &
                     j,(modeldatactx%DELD(i,n),i=1,NDOF)
                else
                    read(luin,fmt=FMT_fedsk3D_velocities,err=1000,end=2000,iostat=io) &
                     j,(modeldatactx%DELD(i,n),i=1,NDOF)
                endif
                if (j.ne.n) then
                    write(stderr,95) n
95                  format(1x,'femout [e]: synchronization error', &
                     'in velocities data block for node ',I12)
                    ierr = EXEC
                    return
                endif
            enddo
        else
            read(luin,err=1000,end=2000,iostat=io) &
             ((modeldatactx%DELD(i,n),i=1,NDOF),n=1,NUMNP)
        endif
        if (DELTP.eq.0d0) then
            ! determine maximum displacement increment
            maxv = -1d0
            do n=1,NUMNP
                maxv = MAX(SQRT(DELD(1,n)**2+DELD(2,n)**2+DELD(3,n)**2),maxv)
            enddo
            if (maxv.gt.1D-2) then
                write(stderr,92) IREQ,maxv
92              format('femout **WARNING**: displacement increment during zero-size ',&
                 'timestep ',I5,' more than 1cm (',1PE12.6,' m).'/ &
                 'Zeroing velocities')
            endif
            call clear_matrix(modeldatactx%DELD, NDOF, NUMNP)
        endif
        if (iecho.eq.2) write(stderr,96)
96      format(1x,'velocities have been read')
    endif 

    !******************** DXE (differential displacement)  *******************
    datatype = "differential displacements"
    if (NUMSN.gt.0) then

#ifdef EXPERIMENTAL_ReducedSlipAlloc
        call CLEAR (DXE,NDOF*NEN*nElemsWithSlipperyNodes)
        if (ascin) then
            if (NDOF.eq.2) then
                frmt = FMT_fedsk2D_diff_displacement
            else if (NDOF.eq.3) then
                frmt = FMT_fedsk3D_diff_displacement
            else
                write(stderr,*) "diff disp says: NDOF not 2 or 3; contact model support"
            endif
    
            if (iecho.eq.2) then
                write(stderr,*) 'Going to read diff displ with format: ', frmt
            endif
    
            do n=1,nElemsWithSlipperyNodes
                do k=1,NEN  ! quick and dirty fix
                    read(luin,"(a)",err=1000,end=2000, iostat=io) dataRecord
    
                    read(dataRecord,fmt=frmt,err=1000,end=2000, iostat=io) j,jj,(DXE(i,k,n),i=1,NDOF)
    
                enddo
            enddo
        else
            ! read as binaries
            read(luin,err=1000,end=2000,iostat=io) &
            (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
        endif
#else
        call CLEAR (DXE,NDOF*NEN*NUMEL,"DXE")
        if (ascin) then
            if (NDOF.eq.2) then
                frmt = FMT_fedsk2D_diff_displacement
            else if (NDOF.eq.3) then
                frmt = FMT_fedsk3D_diff_displacement
            endif
            do n=1,NUMEL
                do k=1,NEN
                    read(luin,fmt=frmt,err=1000,end=2000, iostat=io) j,jj,(DXE(i,k,n),i=1,NDOF)
                    if (j.ne.n) then
  100                   format(1x,'femout [f]: synchronization',1x, &
                        'error in slippery nodes data block', &
                        ' for node ',I12)
                        ierr = EXEC
                        return
                    endif
                enddo
            enddo
        else
            ! read as binaries
            read(luin,err=1000,end=2000,iostat=io) &
            (((DXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
        endif
#endif

        if (iecho.eq.2) write(stderr,105)
  105   format(1x,'slippery displacements have been read')

!******************** DELDXE (differential velocities)  *******************
        datatype = "differential velocities"

        if (IREQ.gt.0) then

#ifdef EXPERIMENTAL_ReducedSlipAlloc
            call CLEAR (modeldatactx%DELDXE,NDOF*NEN*nElemsWithSlipperyNodes,"modeldatactx%DELDXE")
#else
            call CLEAR (modeldatactx%DELDXE,NDOF*NEN*NUMEL,"modeldatactx%DELDXE")
#endif
            modeldatactx%DELDXE = 0d0

#ifdef EXPERIMENTAL_ReducedSlipAlloc
            if (ascin) then 
                do n=1,nElemsWithSlipperyNodes
                    do k=1,NEN

                        read(luin,*) dataRecord

                        write(*,*) "read record     ***"// dataRecord(1:50)//"***"


                        read(dataRecord,fmt=frmt,err=1000,end=2000, iostat=io) &
                         j,jj, (modeldatactx%DELDXE(i,k,n),i=1,NDOF)

                        write(*,*) "read diff. velo ", j, jj, modeldatactx%DELDXE(1:NDOF,k,n)

                    enddo
                enddo
            else
                read(luin,err=1000,end=2000,iostat=io) &
                (((modeldatactx%DELDXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,nElemsWithSlipperyNodes)
            endif
            if (DELTP.eq.0d0) then
                ! determine maximum slippery displacement increment
                maxv = -1d0
                do n=1,nElemsWithSlipperyNodes
                    do j=1,NEN
                        maxv = MAX(SQRT(modeldatactx%DELDXE(1,j,n)**2+modeldatactx%DELDXE(2,j,n)**2+modeldatactx%DELDXE(3,j,n)**2),maxv)
                    enddo
                enddo
                if (maxv.gt.1D-2) then
                    write(stderr,106) IREQ,maxv
106                  format('femout **WARNING**: slippery displacement increment during zero-size ',&
                     'timestep ',I5,' more than 1cm (',1PE12.6,' m).'/ &
                     'Zeroing slippery velocities')
                endif
                call CLEAR (modeldatactx%DELDXE,NDOF*NEN*nElemsWithSlipperyNodes)
            endif
#else
            if (ascin) then
                do n=1,NUMEL
                    do k=1,NEN
                        read(luin,fmt=frmt,err=1000,end=2000, &
                         iostat=io) j,jj, &
                         (modeldatactx%DELDXE(i,k,n),i=1,NDOF)
                        if (j.ne.n) then
                            write(stderr,107) n
  107                       format(1x,'femout [g]: synchronization', &
                             ' error in slippery velocities', &
                             ' for node ',I12)
                            ierr = EXEC
                            return
                        endif
                    enddo
                enddo
            else
                read(luin,err=1000,end=2000,iostat=io) &
                (((modeldatactx%DELDXE(i,j,n),i=1,NDOF),j=1,NEN),n=1,NUMEL)
            endif
            if (DELTP.eq.0d0) then
                ! determine maximum slippery displacement increment
                maxv = -1d0
                do n=1,NUMEL
                    do j=1,NEN
                        maxv = MAX(SQRT(modeldatactx%DELDXE(1,j,n)**2+modeldatactx%DELDXE(2,j,n)**2+modeldatactx%DELDXE(3,j,n)**2),maxv)
                    enddo
                enddo
                if (maxv.gt.1D-2) then
                    write(stderr,106) IREQ,maxv
106                 format('femout **WARNING**: slippery displacement increment during zero-size ',&
                     'timestep ',I5,' more than 1cm (',1PE12.6,' m).'/ &
                     'Zeroing slippery velocities')
                endif
                call CLEAR (modeldatactx%DELDXE,NDOF*NEN*NUMEL,"modeldatactx%DELDXE")
            endif
#endif

            if (iecho.eq.2) write(stderr,108)
  108       format(1x,'differential velocities have been read')
        endif

    endif

    !********************** stress ************
    datatype = "stress"
    call clear (STN,NSTR*NUMEL,"STN")
    if (ascin) then
        write(frmt,130) NSTR
  130   format('(I12,1X,',I1,'E15.0)')
        do n=1,NUMEL
            if (nstr.eq.3) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr3,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STN(i,n),i=1,NSTR)
            else if (nstr.eq.4) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr4,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STN(i,n),i=1,NSTR)
            else if (nstr.eq.5) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr5,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STN(i,n),i=1,NSTR)
            else if (nstr.eq.6) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr6,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STN(i,n),i=1,NSTR)
            endif
    
            if (j.ne.n) then
                write(stderr,135) n
  135           format(1x,'femout [j]: synchronization error in', &
                1x,'stresses data block for element ',I12)
                ierr = EXEC
                return
            endif
        enddo
    else
        read(luin,err=1000,end=2000,iostat=io) &
        ((modeldatactx%STN(i,n),i=1,NSTR),n=1,NUMEL)
    endif
    if (iecho.eq.2) write(stderr,140)
140 format(1x,'stresses have been read')

    ! determine stress extrema
    do iElem=1,NUMEL
        effectiveSum = 0d0
        do iDOF=1,NDOF
            effectiveSum = effectiveSum + modeldatactx%STN(iDOF,iElem)
        enddo
        if (effectiveSum.ge.0d0) then
            if (effectiveSum.gt.stressPositiveBiggest) stressPositiveBiggest = effectiveSum
            if (effectiveSum.lt.stressPositiveSmallest) stressPositiveSmallest = effectiveSum
        else if (effectiveSum.le.0d0) then
            if (effectiveSum.lt.stressNegativeBiggest) stressNegativeBiggest = effectiveSum
            if (effectiveSum.gt.stressNegativeSmallest) stressNegativeSmallest = effectiveSum
        endif
    enddo

#ifdef EXPERIMENTALTRANSIENTDISCREEP

    !********************** taylor stress ************
        
    if (IREQ.gt.0) then
    datatype = "taylor stress"
    call clear (STNTAY,NSTRTAY*NUMEL,"STNTAY")
    if (ascin) then
        write(frmt,131) NSTRTAY
  131   format('(I12,1X,',I1,'E15.0)')
        do n=1,NUMEL
            if (nstrtay.eq.1) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr1,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STNTAY(i,n),i=1,NSTRTAY)
            elseif (nstrtay.eq.3) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr3,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STNTAY(i,n),i=1,NSTRTAY)
            else if (nstrtay.eq.4) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr4,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STNTAY(i,n),i=1,NSTRTAY)
            else if (nstrtay.eq.5) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr5,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STNTAY(i,n),i=1,NSTRTAY)
            else if (nstrtay.eq.6) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr6,err=1000,end=2000,iostat=io) &
                 j,(modeldatactx%STNTAY(i,n),i=1,NSTRTAY)
            endif
    
            if (j.ne.n) then
                write(stderr,136) n
  136           format(1x,'femout [j]: synchronization error in', &
                1x,'stresses data block for element ',I12)
                ierr = EXEC
                return
            endif
        enddo
    else
        read(luin,err=1000,end=2000,iostat=io) &
        ((modeldatactx%STNTAY(i,n),i=1,NSTRTAY),n=1,NUMEL)
    endif
    if (iecho.eq.2) write(stderr,141)
141 format(1x,'taylor stresses have been read')

!    ! determine stress extrema
!    do iElem=1,NUMEL
!        effectiveSum = 0d0
!        do iDOF=1,NDOF
!            effectiveSum = effectiveSum + modeldatactx%STN(iDOF,iElem)
!        enddo
!        if (effectiveSum.ge.0d0) then
!            if (effectiveSum.gt.stressPositiveBiggest) stressPositiveBiggest = effectiveSum
!            if (effectiveSum.lt.stressPositiveSmallest) stressPositiveSmallest = effectiveSum
!        else if (effectiveSum.le.0d0) then
!            if (effectiveSum.lt.stressNegativeBiggest) stressNegativeBiggest = effectiveSum
!            if (effectiveSum.gt.stressNegativeSmallest) stressNegativeSmallest = effectiveSum
!        endif
!    enddo
    endif
#endif    




    !****************** strains ****************************
    datatype = "strain"
    call clear (STR,NSTR*NUMEL,"STR")
    if (ascin) then
        do n=1,NUMEL
            ! format adapts to long or short data format
            if (nstr.eq.3) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr3,err=1000,end=2000,iostat=io) &
                j,(modeldatactx%STR(i,n),i=1,NSTR)
            else if (nstr.eq.4) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr4,err=1000,end=2000,iostat=io) &
                j,(modeldatactx%STR(i,n),i=1,NSTR)
            else if (nstr.eq.5) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr5,err=1000,end=2000,iostat=io) &
                j,(modeldatactx%STR(i,n),i=1,NSTR)
            else if (nstr.eq.6) then
                read(luin,fmt=FMT_fedsk_stressStrainRate_nstr6,err=1000,end=2000,iostat=io) &
                j,(modeldatactx%STR(i,n),i=1,NSTR)
            endif
            if (j.ne.n) then
                write(stderr,142) n
  142           format(1x,'femout [k]: synchronization error in', &
                 1x,'strains data block for element ',I12)
                ierr = EXEC
                return
            endif
        enddo
     else
        read(luin,err=1000,end=2000,iostat=io) &
         ((modeldatactx%STR(i,n),i=1,NSTR),n=1,NUMEL)
     endif
     if (iecho.eq.2) write(stderr,143)
143  format(1x,'strains have been read')

     ! Scaling off-diagonal strain components (engineering strain to infinitesimal)
     if (NTYPE.eq.1) then
         ! pln
         do n=1,NUMEL
             modeldatactx%STR(3,n) = modeldatactx%STR(3,n)*5D-1
         enddo
     elseif (NTYPE.eq.2) then
         ! opn
         do n=1,NUMEL
             modeldatactx%STR(3,n) = modeldatactx%STR(3,n)*5D-1
             modeldatactx%STR(4,n) = modeldatactx%STR(4,n)*5D-1
             modeldatactx%STR(5,n) = modeldatactx%STR(5,n)*5D-1
         enddo
     elseif (NTYPE.eq.3) then
         ! axi
         do n=1,NUMEL
             modeldatactx%STR(4,n) = modeldatactx%STR(4,n)*5D-1
         enddo
     elseif (NTYPE.eq.6) then
         ! f3d
         do n=1,NUMEL
             modeldatactx%STR(4,n) = modeldatactx%STR(4,n)*5D-1
             modeldatactx%STR(5,n) = modeldatactx%STR(5,n)*5D-1
             modeldatactx%STR(6,n) = modeldatactx%STR(6,n)*5D-1
         enddo
     endif

     do iElem=1,NUMEL
         effectiveSum = 0d0
         do iDOF=1,NDOF
             effectiveSum = effectiveSum + modeldatactx%STR(iDOF,iElem)
         enddo
         if (effectiveSum.ge.0d0) then
             if (effectiveSum.gt.strainPositiveBiggest) then
                 strainPositiveBiggest = effectiveSum
             endif
             if (effectiveSum.lt.strainPositiveSmallest) then
                 strainPositiveSmallest = effectiveSum
             endif
         else if (effectiveSum.le.0d0) then
             if (effectiveSum.lt.strainNegativeBiggest) then
                 strainNegativeBiggest = effectiveSum
             endif
             if (effectiveSum.gt.strainNegativeSmallest) then
                 strainNegativeSmallest = effectiveSum
             endif
         endif
     enddo

#ifndef EXPERIMENTAL_NORATESTORAGE

     !****************** strain rates ****************************
     datatype = "strain rates"
     if (IREQ.gt.0) then
        if (ascin) then
            write(frmt,145) NSTR
  145       format('(I12,1X,',I1,'E15.0)')
            do n=1,NUMEL
                if (nstr.eq.3) then
                    read(luin,fmt=FMT_fedsk_stressStrainRate_nstr3,err=1000,end=2000,iostat=io) &
                        j,(modeldatactx%RATE(i,n),i=1,NSTR)
                else if (nstr.eq.4) then
                    read(luin,fmt=FMT_fedsk_stressStrainRate_nstr4,err=1000,end=2000,iostat=io) &
                        j,(modeldatactx%RATE(i,n),i=1,NSTR)
                else if (nstr.eq.5) then
                    read(luin,fmt=FMT_fedsk_stressStrainRate_nstr5,err=1000,end=2000,iostat=io) &
                        j,(modeldatactx%RATE(i,n),i=1,NSTR)
                else if (nstr.eq.6) then
                    read(luin,fmt=FMT_fedsk_stressStrainRate_nstr6,err=1000,end=2000,iostat=io) &
                        j,(modeldatactx%RATE(i,n),i=1,NSTR)
                endif
                if (j.ne.n) then
                    write(stderr,150) n
  150               format(1x,'femout [l]: synchronization error', &
                     1x,'in rates data block for element ',I12)
                    ierr = EXEC
                    return
                endif
            enddo
        else
            read(luin,err=1000,end=2000,iostat=io) &
             ((modeldatactx%RATE(i,n),i=1,NSTR),n=1,NUMEL)
        endif

        if (iecho.eq.2) write(stderr,155)
  155   format(1x,'total strain rates have been read')

        ! Scaling off-diagonal strain rate components 
        if (NTYPE.eq.1) then
            ! pln
            do n=1,NUMEL
                modeldatactx%RATE(3,n) = modeldatactx%STR(3,n)*5D-1
            enddo
        elseif (NTYPE.eq.2) then
            ! opn
            do n=1,NUMEL
                modeldatactx%RATE(3,n) = modeldatactx%STR(3,n)*5D-1
                modeldatactx%RATE(4,n) = modeldatactx%STR(4,n)*5D-1
                modeldatactx%RATE(5,n) = modeldatactx%STR(5,n)*5D-1
            enddo
        elseif (NTYPE.eq.3) then
            ! axi
            do n=1,NUMEL
                modeldatactx%RATE(4,n) = modeldatactx%STR(4,n)*5D-1
            enddo
        elseif (NTYPE.eq.6) then
            ! f3d
            do n=1,NUMEL
                modeldatactx%RATE(4,n) = modeldatactx%STR(4,n)*5D-1
                modeldatactx%RATE(5,n) = modeldatactx%STR(5,n)*5D-1
                modeldatactx%RATE(6,n) = modeldatactx%STR(6,n)*5D-1
            enddo
        endif
    
        ! determine strain rate extrema
        do iElem=1,NUMEL
            effectiveSum = 0d0
            do iDOF=1,NDOF
                effectiveSum = effectiveSum + modeldatactx%RATE(iDOF,iElem)
            enddo
            if (effectiveSum.ge.0d0) then
                if (effectiveSum.gt.strainratePositiveBiggest) then
                    strainratePositiveBiggest = effectiveSum
                endif
                if (effectiveSum.lt.strainratePositiveSmallest) then
                    strainratePositiveSmallest = effectiveSum
                endif
            else if (effectiveSum.le.0d0) then
                if (effectiveSum.lt.strainrateNegativeBiggest) then
                    strainrateNegativeBiggest = effectiveSum
                endif
                if (effectiveSum.gt.strainrateNegativeSmallest) then
                    strainrateNegativeSmallest = effectiveSum
                endif
            endif
        enddo
    endif
#endif

    !****************** surface uplifts ****************************

    datatype = "surface uplifts"
    if (NSURF.gt.0.and.NSED.gt.0) then
        if (ascin) then
            write(frmt,160) NDOF
  160       format('(I12,1X,',I1,'E15.0)')
            do n=1,NSURF
                read(luin,fmt=frmt,err=1000,end=2000, &
                iostat=io) j,(SURF(i,n),i=1,NSD)
                if (j.ne.n) then
                    write(stderr,165) n
  165               format(1x,'femout [m]: synchronization err', &
                   1x,'in block for surface node ',I12)
                    ierr = EXEC
                    return
                endif
            enddo
        else
            read(luin,err=1000,end=2000,iostat=io) &
            ((SURF(i,n),i=1,NSD),n=1,NSURF)
        endif
        if (iecho.eq.2) write(stderr,170)
  170   format(1x,'surface uplifts have been read')
    endif

endif ! read mechanical results

return 

!.... ....................... R E A D   E R R O R S ............................

1000    write(stderr,*) 'femout: Error occurred while trying to read ', trim(datatype), &
                        ' of timestep ', IREQ, &
                        ' error: ', io
ierr = EXEC
return

2000    write(stderr,*) 'femout: Unexpected EOF while trying to read ', trim(datatype), &
                        ' of timestep ', IREQ, &
                        ' error: ', io

ierr = EXEC
return

end subroutine

!-------------------------------------------------------------------------------

! subroutine fewrite (D,STN,STR,DXE,IDSLE,IDNT,PRPMAT, &
!   IEN,PRPLAS,RATE,T,FLUX,IREQ,MECHRD,TEMPRD,LUFO,IERR)

subroutine fewrite (D,STN,STR,DXE,IDSLE,PRPMAT, &
                    IEN,PRPLAS,RATE,T,FLUX,IREQ,MECHRD,TEMPRD,LUFO,IERR)

USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use timestepmodule,  only: NINTG,MAXIT
use materialsmodule, only: lmat
use debugmodule,     only: iecho
use errors,          only: exec, ok
#ifndef EXPERIMENTAL_ReducedSlipAlloc
use meshdatamodule,  only: meshdatactx
#endif
use iomodule,        only: ascout, stderr, luin
!
! routine to write FE output file
!
implicit none
!-pass
#ifdef EXPERIMENTAL_ReducedSlipAlloc
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif

logical MECHRD,TEMPRD
integer IREQ, lufo, ierr

integer          :: idsle, ien !,idnt
double precision :: d,stn,prpmat,rate,t,flux,prplas,str

dimension D(NDOF,*),STN(NSTR,*),PRPMAT(6,*), &
 RATE(NSTR,*),T(*),FLUX(NSD,*),PRPLAS(9,*),IDSLE(*), &
   IEN(NEN,*),STR(NSTR,*) ! IDNT(*)
!-locl
character(len=40)frmt
integer :: i, j, m, n, nn, io
!-init
ierr  = OK


if ((.not.temprd).and.(.not.mechrd)) then
    return
endif

!    write temperatures data block
if (TEMPRD) then

    if (ascout) then
        write(lufo,1,err=1000,iostat=io) (n,T(n),n=1,NUMNP)
 1       format(I12,1X,1PG15.6)
    else
        write(lufo,err=1000,iostat=io) (T(n),n=1,NUMNP)
    endif

    if (iecho.eq.2) then
        write(stderr,2)
 2        format(1x,'temperatures have been written')
    endif

    if (ascout) then
        write(frmt,3) NSD
 3       format('(I12,1X,',I1,'(1PG15.6))')
        do n=1,NUMEL
            write(lufo,fmt=frmt,err=1000,iostat=io) &
               n,(FLUX(i,n),i=1,NSD)
        enddo
    else
        write(lufo,err=1000,iostat=io) ((FLUX(i,n),i=1,NSD),n=1,NUMEL)
    endif

    if (iecho.eq.2) then
        write(stderr,4)
 4        format(1x,'heat fluxes have been written')
    endif

endif

if (mechrd) then
    if (IREQ.gt.0) then

        if (LMAT) then
            m = NUMEL
        else
            m = NUMAT
        endif

        if (ascout) then
            write(lufo,fmt=5,err=1000,iostat=io) &
             (n,PRPMAT(3,n),PRPMAT(4,n),n=1,m)
 5             format(I12,1X,1PG15.6,1X,0PF4.2)
        else
            write(lufo,err=1000,iostat=io) &
             (PRPMAT(3,n),PRPMAT(4,n),n=1,m)
        endif

        if (iecho.eq.2) then
            write(stderr,6)
 6            format(1x,'viscosities have been written')
        endif

        if (ascout) then
             write(lufo,fmt=7,err=1000,iostat=io) &
           (n,PRPLAS(3,n),PRPLAS(4,n),n=1,m)
 7            format(I12,1X,1PG15.6,1X,1PG15.6)
        else
            write(lufo,err=1000,iostat=io) &
           (PRPLAS(3,n),PRPLAS(4,n),n=1,m)
        endif

        if (iecho.eq.2) then
            write(stderr,8)
 8            format(1x,'plastic parameters have been written')
        endif
    endif

    if (ascout) then
        write(frmt,9) NDOF
 9       format('(I12,1X,',I1,'(1PG15.6))')
        do n=1,NUMNP
            write(lufo,fmt=frmt,err=1000,iostat=io) &
            n,(D(i,n),i=1,NDOF)
        enddo
    else
        write(lufo,err=1000,iostat=io) &
             ((D(i,n),i=1,NDOF),n=1,NUMNP)
    endif

    if (iecho.eq.2) then
        write(stderr,10)
10        format(1x,'displacements have been written')
    endif

    if (NUMSN.gt.0) then
        if (ascout) then
            write(frmt,11) NDOF
11            format('(I12,1X,I5,1X,',I1,'(1PG15.6))')
            do nn=1,NUMSE
                n = IDSLE(nn)
                do j=1,NEN
                    write(lufo,fmt=frmt,err=1000,iostat=io) &
                    n,IEN(j,n),(DXE(i,j,n),i=1,NDOF)
                enddo
            enddo
        else
             write(lufo,err=1000,iostat=io) (((DXE(i,j,IDSLE(n)), &
               i=1,NDOF),j=1,NEN),n=1,NUMSE)
        endif

        if (iecho.eq.2) then
            write(stderr,12)
12            format(1x,'differential displacements written')
        endif
    endif


    if (ascout) then
        write(frmt,13) NSTR
13        format('(I12,1X,',I1,'(1PG15.6))')
        do n=1,NUMEL
            write(lufo,fmt=frmt,err=1000,iostat=io) &
           n,(STN(i,n),i=1,NSTR)
        enddo
    else
        write(lufo,err=1000,iostat=io) &
             ((STN(i,n),i=1,NSTR),n=1,NUMEL)
    endif

    if (iecho.eq.2) then
        write(stderr,14)
14        format(1x,'stresses have been written')
    endif

    do n=1,NUMEL
        STR(4,n) = STR(4,n)*2d0
        STR(5,n) = STR(5,n)*2d0
        STR(6,n) = STR(6,n)*2d0
    enddo

    if (ascout) then

        do n=1,NUMEL
            write(lufo,fmt=frmt,err=1000,iostat=io) &
           n,(STR(i,n),i=1,NSTR)
        enddo

    else
        write(lufo,err=1000,iostat=io) &
             ((STR(i,n),i=1,NSTR),n=1,NUMEL)
    endif

    do n=1,NUMEL
        STR(4,n) = STR(4,n)*5D-1
        STR(5,n) = STR(5,n)*5D-1
        STR(6,n) = STR(6,n)*5D-1
    enddo

    if (iecho.eq.2) then
        write(stderr,15)
15        format(1x,'strains have been written')
    endif

    if (IREQ.gt.0) then

        do n=1,NUMEL
            RATE(4,n) = RATE(4,n)*2d0
            RATE(5,n) = RATE(5,n)*2d0
            RATE(6,n) = RATE(6,n)*2d0
        enddo

        if (ascout) then
            write(frmt,16) NSTR
16            format('(I12,1X,',I1,'(1PG15.6))')
            do n=1,NUMEL
                write(lufo,fmt=frmt,err=1000,iostat=io) &
                  n,(RATE(i,n),i=1,NSTR)
            enddo
        else
            write(lufo,err=1000,iostat=io) &
               ((RATE(i,n),i=1,NSTR),n=1,NUMEL)
        endif

        do n=1,NUMEL
            RATE(4,n) = RATE(4,n)*5D-1
            RATE(5,n) = RATE(5,n)*5D-1
            RATE(6,n) = RATE(6,n)*5D-1
        enddo

        if (iecho.eq.2) then
            write(stderr,17)
17            format(1x,'total strain rates have been written')
        endif
    endif
endif

900    return 
!.... ................... W R I T E   E R R O R S .........................
1000    write(stderr,901) io
 901    format(1x,'fewrite: write error ',I3)
call prerr('fewrite')
ierr = EXEC
goto 900
end subroutine
!-------------------------------------------------------------------------------
subroutine FOUTCP (IREQ,NREQ,IERR)

USE MODELDEFINITION, only: MODE
use iomodule,        only: stderr, luin, lufo
use errors,          only: ok, exec

! Routine to copy FE output at requested time IREQ.

implicit none
!-pass
integer IREQ,NREQ,IERR




if (MODE.ne.2 .and. MODE.ne.5 .and. MODE.ne.6) then
    write(stderr,1) MODE
 1        format(1x,'FoutCp: MODE = ',I1,' no data on output file')
    ierr = EXEC
    return
endif
!
if (IREQ.eq.0) then
    rewind(luin)
    rewind(lufo)
    if (ierr.ne.OK) return
    call OPNDSK(.false.,ierr)
    if (ierr.ne.0) ierr = EXEC
    if (ierr.ne.OK) return
endif
!
!    copy data block
call fcp (IREQ,IERR)

if (IERR.ne.OK) return
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine FOUTRD (IREQ,IERR)

USE ALGEBRA
USE MODELDEFINITION, only: MODE, numsn, nsurf
USE MODELDATAMODULE
USE MESHDATAMODULE
USE filemodule,      only: index_filename,BaseFedsk
use timestepmodule,  only: chktim,SetTimeStepSize,timdif,maxitime, &
                           timestepdatactx, nstep
use constants,       only: one
use debugmodule,     only: iecho
use errors,          only: exec, ok
use modeltopology,   only: nsd, ndof
use iomodule,        only: fedsk, mout, tout, lufo, luin, ascin, ifeout, stderr
use postProcessData, only: deffac,setria,deformed

! Routine to read FE output at requested time IREQ.
! ifeout = 0, (re-)opened output file - initialize pointers &
!     > 0, re-entering the routine.
! NSTEP = step that has been read in last reading operation

implicit none
!-pass
integer IREQ,IERR



!-locl
integer NSKIP,ITIME
integer :: i,j
logical there
double precision :: minone
integer :: icoord
integer, external :: ifopen, ifbopen, lnblk
!-init
data ITIME/0/
minone = -one

if (ifeout.eq.0) then
    rewind(luin)
    NSTEP  = -1
endif

if (MODE.ne.2 .and. MODE.ne.5 .and.  MODE.ne.6) then
    write(stderr,1) MODE
1   format(1x,'FoutRd: MODE = ',I1,' no data on output file')
    ierr = EXEC
    goto 100
endif

NSKIP = IREQ - NSTEP
if (NSKIP.eq.0   .and.  ifeout.ne.0) then
    return
endif

! Check if requested time is in output
call chktim (IREQ,MOUT,TOUT)

if ((.not.MOUT).and.(.not.TOUT)) then
    ierr = EXEC
    goto 100
endif

FEDSK = index_filename(baseFEDSK, MAXITIME, IREQ, .false.)
inquire(file=FEDSK, exist=there)
if (.not.there) then
    write(stderr, 4) '"', FEDSK(1:lnblk(FEDSK)), '"'
4   format(1x, 'FE-output file NOT found: ', 258a)
    goto 100
endif

! close old file...
close(luin)

! ...and open new one
if (ascin) then
    ierr = ifopen (luin, trim(FEDSK), 'old')
    if (ierr.ne.0) then
        ierr = EXEC
        goto 100
    endif
else
    ierr = ifbopen(luin, FEDSK, 'old')
    if (ierr.ne.0) then
        ierr = EXEC
        goto 100
    endif
endif

if (IREQ.eq.0) then
    modeldatactx%DELD = 0d0
    if (NUMSN.gt.0) then
        modeldatactx%DX = 0d0
    endif
    if (ierr.ne.OK) then
        return
    endif
endif

! set timestep data to current value
NSTEP = IREQ
call SetTimeStepSize()

!**** actually read the fedsk file *******
call fread (ireq, ierr)
if (ierr.ne.OK) then
    ifeout = 0
    modeldatactx%d = 0d0
    goto 100
endif

if (iecho.eq.1) then
    write(stderr,2) NSTEP
2   format(1x,'FE-output file has been read successfully',1x, &
     '(itime=',i5,')')
endif

ifeout = 1

! set faulted node displacements and velocities
call setfaulteddispvelo (NSTEP,modeldatactx%NFAULT,modeldatactx%FAULT, &
  modeldatactx%TFAULT,modeldatactx%DFAULT, &
  timestepdatactx%MAXSTP,timestepdatactx%DELT, &
  modeldatactx%SKEW)

icoord = 0
call extrme ()

if (deformed) then
    setria=.true.
endif

if (NSTEP.eq.0) then
   call ADDSRF (modeldatactx%ISURF,modeldatactx%SURF, meshdatactx%X, modeldatactx%D, NSD,NDOF,NSURF,0)
endif

100    return
end
!-------------------------------------------------------------------------------
subroutine fread (ireq, ierr)

USE MODELDEFINITION, only: MODE, numfn, numslp, numrot, numel
USE ALGEBRA,         only: equate, clear, update
use modeldatamodule, only: modeldatactx
use meshdatamodule,  only: meshdatactx
use timestepmodule,  only: chktim, nstep
use materialsmodule, only: matdatactx
use modeltopology,   only: ndof, nen
use plotControl,     only: irelt
use errors,          only: ok

! Routine to read a single data block

implicit none




!-locl
logical temprd,mechrd,RDTOO
integer ierr
integer :: numnp
integer :: ireq

numnp = meshdatactx%nvglobal

!-init
ierr = OK
IREQ = NSTEP

! determine if mechanical block needs to be read
mechrd = (MODE.eq.2.or.MODE.eq.6)
temprd = (MODE.eq.5.or.MODE.eq.6)
call chktim (IREQ,Mechrd,Temprd)

if (temprd) then
    call EQUATE (modeldatactx%TLST, modeldatactx%T,NUMNP)
    if (NUMFN+NUMSLP.gt.0) then
         call EQUATE (modeldatactx%TLSTX,modeldatactx%TX, NUMNP)
    endif
endif

! read the output file
#ifdef EXPERIMENTAL_NORATESTORAGE
call femout (modeldatactx%d, &
      modeldatactx%STN, &
      modeldatactx%STR, &
      modeldatactx%DXE,matdatactx%prpmat, &
      matdatactx%PRPLAS,&
      modeldatactx%T,modeldatactx%TX, modeldatactx%Thermalflux, &
      modeldatactx%SURF,modeldatactx%idsle,&
      modeldatactx%DELD, &
      modeldatactx%DX,IREQ,mechrd,temprd,ierr)
#elseif EXPERIMENTALTRANSIENTDISCREEP
call femout (modeldatactx%d, &
      modeldatactx%STN, &
      modeldatactx%STNTAY, &
      modeldatactx%STR, &
      modeldatactx%DXE,matdatactx%prpmat, &
      matdatactx%PRPLAS,&
      modeldatactx%T,modeldatactx%TX, modeldatactx%Thermalflux, &
      modeldatactx%SURF,modeldatactx%idsle,&
      modeldatactx%DELD, &
      modeldatactx%DX,IREQ,mechrd,temprd,ierr)
#else
!    read the output file
call femout (modeldatactx%d, &
      modeldatactx%STN, &
      modeldatactx%STR, &
      modeldatactx%DXE,matdatactx%prpmat, &
      matdatactx%PRPLAS,&
      modeldatactx%RATE,modeldatactx%T,modeldatactx%TX, modeldatactx%Thermalflux, &
      modeldatactx%SURF,modeldatactx%idsle,&
      modeldatactx%DELD, &
      modeldatactx%DX,IREQ,mechrd,temprd,ierr)
#endif

if (ierr.ne.OK) goto 100

NSTEP = IREQ

if (NSTEP.eq.irelt .and. temprd) then
   call EQUATE(modeldatactx%T0,modeldatactx%T,NUMNP)
endif

100    return
end
!-------------------------------------------------------------------------------
subroutine faultequalize (NFAULT,FAULT,NUMFN,NDOF)
!
! recomputes faulted displacement/velocity such that it is equally distributed
! on both sides of a fault. This is relevant for plotting displacements, which are
! relative to the normal/average node position
!
use iomodule,    only: stderr
use debugmodule, only: iecho

implicit none
!-pass
integer :: NFAULT(3,*),NUMFN,NDOF
double precision :: FAULT(NDOF,*)
!-locl
logical :: doecho
integer :: nf,node,idof,j
double precision :: minslp(3),maxslp(3),shift,slip,ffe

do nf=1,NUMFN

    ! already processed nodes are/will be labeled with a negative node number entry
    if (NFAULT(2,nf).gt.0) then

        ! find maximum and minimum faulted displacement for this node
        node = NFAULT(2,nf)
        do idof=1,NDOF
            minslp(idof) = FAULT(idof,nf)
            maxslp(idof) = FAULT(idof,nf)
        enddo
        NFAULT(2,nf) = -NFAULT(2,nf)
        do j=nf+1,NUMFN
            if (NFAULT(2,j).eq.node) then
                do idof=1,NDOF
                   if (FAULT(idof,j).lt.minslp(idof)) minslp(idof)=FAULT(idof,j)
                   if (FAULT(idof,j).gt.maxslp(idof)) maxslp(idof)=FAULT(idof,j)
                enddo
                NFAULT(2,j) = -NFAULT(2,j)
            endif
        enddo
        
!       test
        do idof=1,NDOF
            if (minslp(idof).ne.0d0 .and. maxslp(idof).ne.0d0) then
                if (minslp(idof).eq.maxslp(idof)) then
                    write(stderr,10) nf,node,idof
10                  format('Input error for NFAULT entry ',I6,':', &
                     'no differential slip for node ',I12,' dof=',I1,/ & 
                     'Probably, entries for one fault side are missing')
                    call exitp(1)
                endif
            endif
        enddo

        ! recompute all differential displacement entries for this node
        do j=nf,NUMFN
            if (IABS(NFAULT(2,j)).eq.node) then
                doecho = .false.
                do idof=1,NDOF
                    shift = ( maxslp(idof) + minslp(idof) ) * 0.5d0
                    slip  = ( maxslp(idof) - minslp(idof) ) * 0.5d0
                    doecho = (doecho .or. (shift.gt.0d0 .and. iecho.eq.2) )
                    ffe = ABS(FAULT(idof,j)-shift) ! = magnitude of expected new value
                    if (ABS(ffe-slip).gt.1d-25) then
                        write(stderr,20) j,node,idof
20                      format('Input error for NFAULT entry ',I6,':', &
                         'variable differential slip in FAULT for node ',I12,' dof=',I1)
                        call exitp(1)
                    endif
                    FAULT(idof,j) = SIGN(slip,FAULT(idof,j)-shift)
                enddo
                if (doecho) then
                    write(*,15) j,NFAULT(1,j),IABS(NFAULT(2,j)), &
                     (FAULT(idof,j),idof=1,NDOF)
15                  format('FAULT entry ',I6,' for element ',I12,' and node ',&
                     I12,': split redefined to',3(1x,1PE12.4))
                endif
            endif
        enddo

    endif

enddo

! make node numbers positive again
do nf=1,NUMFN
    NFAULT(2,nf) = IABS(NFAULT(2,nf))
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine setfaulteddispvelo (IREQ,NFAULT,FAULT,TFAULT,DFAULT,MAXSTP,DELT,SKEW)

! loads faulted node displacements and velocities at requested time step number IREQ
! TFAULT holds the total differential displacement of split nodes
! DFAULT hold the differential velocity between IREQ-1 and IREQ

use ALGEBRA,         only: EQUATE,VMPRD,FORMRT
use CONSTANTS,       only: ZERO,EPS
use MODELTOPOLOGY,   only: NDOF
use TIMESTEPMODULE,  only: NINTG,faultvelo
use MODELDEFINITION, only: NUMROT,NUMFN,MODE
implicit none
!-pass
integer :: IREQ,NFAULT,MAXSTP
double precision :: FAULT,TFAULT,DFAULT,DELT,SKEW
dimension NFAULT(3,*),FAULT(NDOF,*),TFAULT(NDOF,*),DFAULT(NDOF,*),MAXSTP(*),DELT(*),SKEW(2,*)
!-locl
logical :: VELBC,doit
integer :: iflt,idof,node,itime,igrp,grpindx
double precision :: rot(3,3),dl(3)

VELBC = (faultvelo.ne.0)

if ((MODE.ne.2 .and. MODE.ne.6) .or. NUMFN.eq.0) return

! initialize TFAULT with zero-time displacements, and  DFAULT=0 (large array,
! so don't do a wholesale zeroing)
do iflt=1,NUMFN
    do idof=1,NDOF
        TFAULT(idof,iflt) = ZERO
        DFAULT(idof,iflt) = ZERO
    enddo
    if (.not.VELBC) then
        if (NFAULT(3,iflt).le.0) then
            do idof=1,NDOF
                TFAULT(idof,iflt) = FAULT(idof,iflt)
            enddo
        endif
    endif
enddo

if (NINTG.eq.0 .or. IREQ.eq.0) then
    if (NUMROT.gt.0) then
        ! Euler rotations on TFAULT before we go
        do iflt=1,NUMFN
            node = NFAULT(2,iflt)
            if (ABS(SKEW(1,node)).gt.EPS .or. ABS(SKEW(2,node)).gt.EPS) then
                call FORMRT(SKEW(:,node),rot,NDOF)
                call EQUATE(dl,TFAULT(1:NDOF,iflt),NDOF)
                call VMPRD (rot,dl,TFAULT(1,iflt),NDOF,3)
            endif
        enddo
    endif
    return
endif

! step through timestep groups and increment TFAULT to the value at itime=IREQ
do iflt=1,NUMFN
    itime = 0
    igrp = 1
    loop_ireq: do while (itime.lt.IREQ)
        grpindx = 1 
        loop_grps: do while (grpindx.le.MAXSTP(igrp))
            itime = itime + 1
            if (VELBC) then
                do idof=1,NDOF
                    TFAULT(idof,iflt) = TFAULT(idof,iflt) + DELT(igrp)*FAULT(idof,iflt)
                enddo
            else
                if (NFAULT(3,iflt).ge.0) then
                    doit = (itime.eq.NFAULT(3,iflt))
                else
                    doit = (MOD(itime,-NFAULT(3,iflt)).eq.0)
                endif
                if (doit) then
                    do idof=1,NDOF
                        TFAULT(idof,iflt) = TFAULT(idof,iflt) + FAULT(idof,iflt)
                    enddo
                endif
            endif
            if (itime.eq.IREQ) then
                exit loop_ireq
            endif
            grpindx = grpindx + 1
        enddo loop_grps
        igrp = igrp + 1
    enddo loop_ireq
    if (VELBC) then
        do idof=1,NDOF
            DFAULT(idof,iflt) = FAULT(idof,iflt)
        enddo
    else
        ! there is no realistic velocity that can be attributed to incidental displacement
        ! therefore, only attribute a velocity for NFAULT(3,:) =  -1
        do idof=1,NDOF
            if (NFAULT(3,iflt).eq.-1 .and. DELT(igrp).gt.EPS) then
                DFAULT(idof,iflt) = FAULT(idof,iflt)/DELT(igrp)
            else
                DFAULT(idof,iflt) = ZERO
            endif
        enddo
    endif
    ! apply Euler rotations
    node = NFAULT(2,iflt)
    if (ABS(SKEW(1,node)).gt.EPS .or. ABS(SKEW(2,node)).gt.EPS) then
        call FORMRT(SKEW(:,node),rot,NDOF)
        call EQUATE(dl,TFAULT(1:NDOF,iflt),NDOF)
        call VMPRD (rot,dl,TFAULT(1,iflt),NDOF,3)
        call EQUATE(dl,DFAULT(1:NDOF,iflt),NDOF)
        call VMPRD (rot,dl,DFAULT(1,iflt),NDOF,3)
    endif
enddo

return
end
!-------------------------------------------------------------------------------
subroutine fskip (ireq,NSKIP,ierr)

use errors,     only: ok

! Routine to skip NSKIP data blocks

implicit none

!-pass
integer ireq,NSKIP,ierr
!-locl
 integer i

if (NSKIP.le.0) then
    return
endif

do i=1,NSKIP
    call fread (ireq,ierr)
    if (ierr.ne.OK) then
        return
    endif
enddo

return

end subroutine
!-------------------------------------------------------------------------------
integer function iSoFil (FIN,FOUT,NF,MAXF)

use iomodule,        only: stderr
use errors,          only: ok, exec, prgerr


! retrieves all files nested in FIN (by ".so")

implicit none
!-pass
integer MAXF,NF,ierr
character FIN*(*),FOUT*(*)
dimension FOUT(MAXF)
!-locl
integer, parameter :: maxnest=10
#include "lenchar.i"
logical there
integer fnblnk,lui(maxnest)
character(len=LENCHAR) tmp,file,dirname,home
character(len=255)record
integer, external :: flnblk, lnblk,nextlu
integer :: m, n
integer :: ilevel, l
!-init
 lui(1) = nextlu(0)
 NF     = 0
 ierr   = OK
 ilevel = 1

 l = LEN(FIN)

 if (l.le.0) then
write(stderr,1)
 1        format(1x,'iSoFil: length of FIN error')
ierr = PRGERR
goto 1000
 endif

 m = LEN(FOUT(1))

 if (m.le.l) then
write(stderr,2)
 2        format(1x,'iSoFil: length of FOUT error')
ierr = PRGERR
goto 1000
 endif

 n = fnblnk(FIN)
 m = lnblk(FIN)

 if (n.ge.m) then
write(stderr,3)
 3        format(1x,'iSoFil: empty input file error')
ierr = EXEC
goto 1000
 endif

 file = FIN(n:m)

100    inquire (file=file,exist=there)

 if (.not.there) then
write(stderr,4) '"',file(n:m),'"'
 4        format(1x,'iSoFil: file does''nt exist: ',80A)
ierr = EXEC
goto 1000
 endif

 open(unit=lui(ilevel),status='old',file=file,err=2000)
 rewind(lui(ilevel),err=2000)
 NF = NF + 1

 if (NF.gt.MAXF) then
write(stderr,5) MAXF
 5       format(1x,'iSoFil: number of include files exceeds',1x,', maximum = ',i3)
ierr = PRGERR
goto 1000
 endif

 FOUT(NF) = file(n:m)

200    if (ilevel.gt.0) then
300        read(lui(ilevel),6,err=2000,end=400) record
 6        format(A80)
if (record(1:4).ne.'.so ') goto 300
! ".so "
!        increment nesting level
ilevel = ilevel + 1
if (ilevel.gt.maxnest) then
write(stderr,7) maxnest
 7        format(1x,'iSoFil: maximum ".so" nesting depth = ',i3)
ierr = EXEC
goto 1000
endif
!        read file name
l = lnblk(record)
if (l.le.4) then
write(stderr,8) '"',file(n:m),'"'
 8        format(1x,'iSoFil: ".so "-error in ',80A)
ierr = EXEC
goto 1000
endif
lui(ilevel) = nextlu(0)
tmp = record(5:l)
n   = fnblnk(tmp)
m   = l - 4
if (tmp(n:n).eq.'/') then
file = tmp(n:m)
else
home = dirname(FIN)
file = home(1:lnblk(home)) // '/' // tmp(n:m)
n    = fnblnk(file)
m    = lnblk(file)
endif
goto 100
! "EOF"
400        close(unit=lui(ilevel),err=2000)
ilevel = ilevel - 1
goto 200
 endif

1000    iSoFil = ierr
 return
2000    write(stderr,9) '"',file(n:m),'"'
 9    format(1x,'iSoFil: I/O error on file ',80a)
 call prerr('iSoFil')
 ierr = EXEC
 goto 1000
 end

!-------------------------------------------------------------------------------

 subroutine itst1 (A,N,istep)

use iomodule, only: stderr

implicit none

integer :: N

 integer A
 dimension A(N)

integer istep,ii,jj

write(stderr,'(''A('',I3,'')'')') N
jj = MAX0(istep,1)

do ii=1,N
    if (MOD(ii,jj).eq.0) then
        write(stderr,'(I4,I12)') ii,A(ii)
    endif
enddo

return
end
!-------------------------------------------------------------------------------
subroutine itst2 (A,N,M,istep)

use iomodule, only: stderr

implicit none

integer :: n, m

integer :: A
 dimension A(N,M)
!-local
integer :: ii, jj, istep, nn

write(stderr,'(''A('',I3,'','',I3,'')'')') N,M
jj = MAX0(istep,1)

do ii=1,M
    if (MOD(ii,jj).eq.0) then
        write(stderr,'(I4,8I12)') ii,(A(nn,ii),nn=1,N)
    endif
enddo

return
end
!-------------------------------------------------------------------------------
subroutine lstout(IMPRINT,ITPRINT,MAXSTP)

USE MODELDEFINITION, only: MODE,ICVIS,ISTART,NELAST
use iomodule,        only: stderr
use timestepmodule,  only: nmprt, ntprt, nintg

implicit none

! Check if there are calculations after last result

!-pass
integer         :: IMPRINT(*),ITPRINT(*),MAXSTP(*)
!-locl
integer         :: LstMch,LstTmp,LstStp,LstPrt
integer         :: i

!-init
LstMch = 0
LstTmp = 0
LstStp = 0

if (NINTG.eq.0) then
    return
endif

do i=1,NINTG
    LstStp = LstStp + MAXSTP(i)
enddo

if (LstStp.eq.0) then
    write(stderr,10) NINTG
10       format(1x,'Inconsistency in time step data: NINTG = ',I3,1x, &
     'and last step number = 0')
    return
endif

if (MODE.eq.2 .or. MODE.eq.6) then
    if (NMPRT.gt.0) then
        LstMch = IMPRINT(NMPRT)
    else if (NMPRT.lt.0) then
        LstMch = INT(LstStp/NMPRT) * NMPRT
    endif
endif

if (MODE.eq.5 .or. MODE.eq.6) then
    if (NTPRT.gt.0) then
        LstTmp = ITPRINT(NTPRT)
    else if (NTPRT.lt.0) then
        LstTmp = INT(LstStp/NTPRT) * NTPRT
    endif
endif

LstPrt = MAX(LstMch,LstTmp)

if (LstPrt.lt.LstStp) then
    if (MODE.eq.2) then
        write(stderr,20) LstStp,LstMch
20           format(1x,'Too much work: last calculation at time',1x, &
         'step ',I6/1x, &
         '               last mechanical output at     ',I6)
    else if (MODE.eq.5) then
        write(stderr,30) LstStp,LstTmp
30           format(1x,'Too much work: last calculation at time',1x, &
         'step ',I6/1x, &
         '               last thermal output at        ',I6)
    else
        write(stderr,40) LstStp,LstMch,LstTmp
40           format(1x,'Too much work: last calculation at time',1x, &
         'step ',I6/1x, &
         '               last mechanical output at     ', &
         I6/1x,'               last thermal output at        ', &
         I6)
    endif
endif

return
end subroutine
!-------------------------------------------------------------------------------
integer function openform (lu,f,ascii)

use errors,     only: exec

implicit none

#include "lenchar.i"
!-pass
integer   lu
character f*(*)
logical   ascii
!-locl
integer           :: ierr
integer, external :: ifopen, ifbopen
logical           :: AsciF

ascii = AsciF(f)

if (ascii) then
    ierr = ifopen (lu,f,'old')
else
    ierr = ifbopen(lu,f,'old')
endif

if (ierr.ne.0) then
    ierr = EXEC
endif

50    openform = ierr
return
end function
!-------------------------------------------------------------------------------
subroutine rdsymb (numsy,SYMB,NSYM,MAXSYM)

USE ALGEBRA
use debugmodule, only: iecho
use iomodule,    only: stderr

implicit none

!-pass
integer numsy,NSYM,MAXSYM
double precision SYMB(2,MAXSYM)

!-locl
character(len=2)charnu
character(len=255) fname
character(len=78) PU

integer :: luloc
integer, external :: lnblk,fnblnk, nextlu
integer :: i, j, l, n
integer :: ios, itype

NSYM = 0
call get_environment_variable('PU',PU)
l = lnblk(PU)

if (l.le.0) then
    write(stderr,1)
 1        format(1x,'rdsymb: P-util directory $PU undefined')
    call exitp(1)
endif

write(charnu,2,err=1000,iostat=ios) numsy
 2    format(I2)
fname = PU(1:l) // '/symbols/' // &
charnu(fnblnk(charnu):lnblk(charnu))
luloc = nextlu(0)
call openf(luloc,fname,'old')
call clear(SYMB,2*MAXSYM,"symb")
l = 0
100    read(luloc,*,err=2000,end=200,iostat=ios) itype

if (itype.eq.-998.or.itype.eq.-999.or.itype.eq.-997) then
    goto 100
endif

backspace(luloc,err=2000,iostat=ios)
read(luloc,*,err=2000,end=2000,iostat=ios) itype,n

if (itype.ne.2) then
    write(stderr,3)
 3        format(1x,'rdsym: only data type 2 implemented')
    call exitp(1)
endif

l = l + 1
SYMB(1,l) = DBLE(itype)
SYMB(2,l) = DBLE(n)
read(luloc,*,err=2000,end=2000,iostat=ios) &
((SYMB(i,j+l),i=1,2),j=1,n)
l = l + n

if (l.gt.MAXSYM) then
    write(stderr,4) numsy
 4        format(1x,'rdsym: overflow for symbol ',I2)
    call exitp(1)
endif

goto 100
200    call closef(luloc)

if (iecho.ge.4) then
    write(stderr,5) l,numsy
 5       format(1x,'rdsym: read ',I2,' entries for symbol number ', &
     I2)
endif

NSYM = l
return

1000    write(stderr,1001) ios,numsy
1001    format(1x,'rdsymb: internal read error ',I3,' while reading ', I7)
 call prerr('rdsymb')
 call exitp(1)

2000    write(stderr,2001) ios
2001    format(1x,'rdsymb: input error ',I3)
 call prerr('rdsymb')
 call exitp(1)

end subroutine
!-------------------------------------------------------------------------------
logical function RDTOO (IPRINT,ISTEP,NPRT,LMODE)

! determines if output is on file
implicit none

!-pass
integer :: ISTEP,NPRT
logical :: LMODE
integer :: IPRINT(*)
!-local
integer :: iprt
!-init

if (.not. LMODE) then
    RDTOO = .false.
    goto 300
endif

if (ISTEP.eq.0) then
    RDTOO = .true.
    goto 300
endif

RDTOO = .false.

if (NPRT.eq.0) then
    goto 300
endif

if (NPRT.lt.0) then
    if (MOD(ISTEP,IABS(NPRT)).eq.0) then
        RDTOO = .true.
    endif
    goto 300
endif

!    check if ISTEP occurs in IPRINT
IPRT = 1

100    if (IPRT.le.NPRT) then
    if (IPRINT(IPRT).eq.ISTEP) then
        goto 200
    endif
    IPRT = IPRT + 1
    goto 100
endif

goto 300

200    RDTOO = .true.
300    return
end function
!------------------------------------------------------------------------------
subroutine sclxtr (A,N,AMIN,AMAX)

use constants, only: dfmax

! Routine for determining bounds on scalar data in array A.

implicit none
!-pass
integer N
double precision A(N)
double precision amin, amax
!-local
integer :: i


!-init
AMIN =  dfmax
AMAX = -dfmax

do i=1,N
    if (A(i).gt.AMAX) then
        AMAX = A(i)
    endif
    if (A(i).lt.AMIN) then
        AMIN = A(i)
    endif
enddo

return
end
!-----------------------------------------------------------------------
subroutine setio()

USE MODELCTX,    only: part_name
use debugmodule, only: lustat
use iomodule,    only: outfil

implicit none



OUTFIL(1)=.FALSE.
OUTFIL(2)=.FALSE.

lustat = 0
return
end
!-------------------------------------------------------------------------------
subroutine OPNDSK (openlu,io)

use iomodule,           only: outfil, ascout, stderr, fedsk, luout
!
implicit none
!-pass
logical openlu
integer :: io

!-local
integer, external :: lnblk

if (.not. OUTFIL(2)) then
    return
endif

if (openlu) then
    if (ascout) then
        call openf (luout(2),FEDSK,'unknown')
    else
        call fbopen (luout(2),FEDSK,'unknown')
    endif
endif
!
!        if (ascout) then
!            write(luout(2),1,err=100,iostat=io) VERSION
!    1       format(1x,'TECTON ',A15)
!        else
!            write(luout(2),err=200,iostat=io) VERSION
!        endif
!
return
!

100     write(stderr,101) io,'"',FEDSK(1:lnblk(FEDSK)),'"'
  101   format(1x,'formatted write error ',I5,' on ',258a)
if (.not.openlu) return
call xit(1," ")

200     write(stderr,201) io,'"',FEDSK(1:lnblk(FEDSK)),'"'
  201   format(1x,'unformatted write error ',I5,' on ',258a)
if (.not.openlu) return
call xit(1," ")

end
!-------------------------------------------------------------------------------
subroutine shell(prgname,ierr)

use debugmodule,     only: iecho
use errors,          only: ok, syntax
use iomodule,        only: stderr, fcmd, fstack

implicit none
 
!-pass
integer ierr
character prgname*(*)



!-locl
#include "lenchar.i"
character(len=LENCHAR) :: arg
integer, external      :: lnblk
integer                :: iarg, narg

ierr = OK

narg  = command_argument_count()
fcmd  = ' '
 1    format(1x,80A)
call get_command_argument(0,prgname)
if (narg.gt.0) then
    do iarg=1,narg

        call get_command_argument(iarg,arg)

        if (arg.eq.'echo') then
            iecho = 1
        elseif (arg.eq.'noecho') then
            iecho = 0
        else
            if (lnblk(fcmd).le.0) then
                fcmd = arg
            else
                write(stderr,1) 'Usage: "', &
                 prgname(1:lnblk(prgname)), &
                 ' [no[echo]] [commandfile]'
                ierr = SYNTAX
                fcmd = ' '
                goto 20
            endif
        endif

    enddo
endif

20    fstack(0) = fcmd

if (lnblk(fcmd).le.0) then
    fstack(0) = 'stdin'
endif

return
end

!-------------------------------------------------------------------------------

subroutine tecout (T,lufo)

USE MESHDATAMODULE, only: meshdatactx

!
! writes "tecton" output file (for initial temperatures)
!
implicit none
!-pass
integer lufo
double precision :: T
dimension T(*)

integer :: n
!
do n=1,meshdatactx%nvglobal
    write(lufo,1) n,0,T(n)
 1        format(2I5,1PG14.6)
enddo

write(lufo,2)
 2    format('end initial temperatures')
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine tst1 (A,N,istep)

use iomodule,        only: stderr      

implicit none

integer N
double precision A(N)
integer :: istep

integer :: ii, jj

!integer N
! dimension A(N)
write(stderr,'(''A('',I3,'')'')') N
jj = MAX0(istep,1)

do ii=1,N
    if (MOD(ii,jj).eq.0) then
        write(stderr,'(I4,1PG12.4)') ii,A(ii)
    endif
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine tst2 (A,N,M,istep)

use iomodule,        only: stderr

implicit none

integer          :: n, m, istep
double precision :: A(N,M)
integer          :: ii,jj,nn

write(stderr,'(''A('',I3,'','',I3,'')'')') N,M
jj = MAX0(istep,1)

do ii=1,M
    if (MOD(ii,jj).eq.0) then
         write(stderr,'(I4,4(1PG12.4))') ii,(A(nn,ii),nn=1,N)
    endif
enddo

return
end subroutine

!-------------------------------------------------------------------------------
subroutine tst2Integers (A,N,M,istep)

use iomodule,        only: stderr

implicit none

integer          :: n, m, istep
integer          :: A(N,M)
integer          :: ii,jj,nn

write(stderr,'(''A('',I3,'','',I3,'')'')') N,M
jj = MAX0(istep,1)

do ii=1,M
    if (MOD(ii,jj).eq.0) then
         write(stderr,'(I4,4(1PG12.4))') ii,(A(nn,ii),nn=1,N)
    endif
enddo

return
end




!-------------------------------------------------------------------------------
subroutine tstinp (ITST)

USE MESHDATAMODULE
use modeldatamodule, only: modeldatactx
use materialsmodule, only: grav
use debugmodule,     only: iecho
use modeldefinition, only: mode, itmode, numsn, numnp, numel
use modeltopology,   only: ndof, nen, nsd, nstr
use iomodule,        only: stderr

! subroutine to write input data
!    ITST =  0: write FE-input file data
!    ITST <> 0: write FE-output file data

 implicit none
!-pass
integer :: itst



!-local
integer :: increm

!-init
 increm = 1

if (iecho.ne.2) then
    return
endif

if (ITST.eq.0) then
    write(stderr,'(''FE-input file contents'')')
    write(stderr,'(''nodal coordinates:'')')
    call tst2(meshdatactx%X, &
          NSD,NUMNP,increm)

    write(stderr,'(/''boundary condition codes:'')')
    ! could do with a proper interface, but seems not important enough right now.
    call tst2integers(modeldatactx%IBOND,NDOF,NUMNP,increm)
    write(stderr,'(/''boundary conditions:'')')
    call tst2(modeldatactx%BOND,NDOF,NUMNP,increm)
    write(stderr,'(/''gravity:'')')
    call tst1(grav,NDOF,1)
    write(stderr,'(/''element definitions:'')')
    call itst2(meshdatactx%IEN,NEN,NUMEL,increm)

    if (MODE.ge.3) then
        if (ITMODE.eq.0) then
            write(stderr,'(/''initial temperatures:'')')
            call tst1(modeldatactx%T,NUMNP,increm)
        endif
        write(stderr,'(/''thermal boundary condition codes:'')')
        call itst1(modeldatactx%ITBC,NUMNP,increm)
        write(stderr,'(/''thermal boundary conditions:'')')
        call tst1(modeldatactx%TBC,NUMNP,increm)
    endif
else
    if (MODE.eq.2 .or. MODE.eq.6) then
        write(stderr,'(''Mechanical FE-output files content'')')
        write(stderr,'(''displacements:'')')
        call tst2 (modeldatactx%D,NDOF,NUMNP,increm)
        write(stderr,'(/''differential displacements:'')')

        if (NUMSN.gt.0) then
            call tst2 (modeldatactx%DXE,NDOF*NEN,NUMEL,1)
        endif
        write(stderr,'(/''stresses:'')')
        call tst2 (modeldatactx%STN,NSTR,NUMEL,increm)
    endif

    if (MODE.eq.5 .or. MODE.eq.6) then
        write(stderr,'(''Thermal FE-output files contents'')')
        write(stderr,'(''temperatures:'')')
        call tst1 (modeldatactx%T,NUMNP,increm)
        write(stderr,'(/''heat fluxes:'')')
        call tst2 (modeldatactx%ThermalFlux,NSD,NUMEL,increm)
    endif
endif
return
end
!-------------------------------------------------------------------------------
 subroutine vecxtr (A,N,M,AMIN,AMAX,IDIF)

use constants, only: zero, dfmax

! Routine for determining bounds on vector data in array A.
!         AMIN(1)      <= positive A's   <=   AMAX(1)
!         AMIN(2)      <= negative A's   <=   AMAX(2)

implicit none
!-pass
integer N,M,IDIF
double precision :: A(N,M),AMAX(2),AMIN(2)

!-locl
logical pos,neg
double precision :: B(6)
integer :: i, j
!-init
AMIN(1) =  dfmax
AMAX(1) =  zero
AMAX(2) = -dfmax
AMIN(2) =  zero
pos = .false.
neg = .false.


do j=1,M
    if (IDIF.eq.0) then
        do i=1,N
            B(i) = A(i,j)
        enddo
    else
        call DEVSTN (A,B)
    endif

    do i=1,N
        if (B(i).ge.zero) then
            pos = .true.
            if (B(i).gt.AMAX(1)) AMAX(1) = B(i)
            if (B(i).lt.AMIN(1)) AMIN(1) = B(i)
        else
            neg = .true.
            if (B(i).gt.AMAX(2)) AMAX(2) = B(i)
            if (B(i).lt.AMIN(2)) AMIN(2) = B(i)
        endif
    enddo
enddo

if (.not.pos) AMIN(1) = zero
if (.not.neg) AMAX(2) = zero

return
end
!-------------------------------------------------------------------------------
subroutine WrFt (BTOT,PREQ,NEQ,lufo,ierr)

USE MODELDEFINITION, only: MODE
use modeltopology,   only: ndof
use iomodule,        only: ascout
use errors,          only: ok

! writes disk file footer


implicit none
!-pass
integer :: NEQ,lufo,ierr
double precision :: BTOT(NEQ),PREQ(NDOF)


!-init
integer :: i

ierr = OK

if (MODE.ne.2 .and. MODE.ne.6) then
    return
endif

if (ascout) then
    write(lufo,10,err=1000) (i,BTOT(i),i=1,NEQ)
10        format(I5,1x,E15.6)
else
    write(lufo,err=1000) (BTOT(i),i=1,NEQ)
endif

return

1000    call prerr('WrFt')
call exitp(1)
end subroutine
!-------------------------------------------------------------------------------
! subroutine WrHd (ierr)


!*************************************
! NOTICE:
!
! if this subroutine is ever reanimated,
! the dependence on the A and AA array will have to
! be fixed, as they no longer exist.

!*************************************


!use AAAB,            only: A, AA
!use timestepmodule,  only: NMPRT, NTPRT
!use modeldefinition, only: mode
!use iomodule,        only: ascout, lufo

! Routine to write header info

!implicit none

!-pass
! integer ierr
!
!
!-local
!integer :: MMPRT, MTPRT

!-init
! MMPRT = 0
! MTPRT = 0
! ierr = OK

! if (MODE.eq.2 .or. MODE.eq.6) MMPRT = NMPRT
! if (MODE.eq.5 .or. MODE.eq.6) MTPRT = NTPRT

! if (ascout) then
!write(lufo,10,err=100) NINTG,MMPRT,MTPRT
!10        format(3I5)
! else
!write(lufo,err=100) NINTG,MMPRT,MTPRT
! endif

!call WritHd (A(IPMSTP),imestepdat%delt,A(IPMPRT),AA(IPTPRT), &
!             NMPRT,NTPRT,NINTG,MODE,lufo,ierr)

! if (ierr.ne.OK) goto 50

!50    return

!100    call prerr('WrHd')
! ierr = EXEC
! goto 50
! end
!-------------------------------------------------------------------------------
!subroutine WritHd(MAXSTP,DELT,IMPRINT,ITPRINT,NMPRT,NTPRT,NINTG, &
!                  MODE,lufo,ierr)

!use iomodule,           only: ascout

! Routine for writing disk header

! implicit none
!-pass
! integer NMPRT,NTPRT,NINTG,MODE,lufo
!double precision :: maxstp, delt
!integer          :: imprint, itprint
!dimension MAXSTP(*),DELT(*),IMPRINT(*),ITPRINT(*)
!integer :: ierr

!
!-local
!integer :: i

! ierr = OK

! if (NINTG.le.0) return

! if (ascout) then
!write(lufo,1,err=100) (MAXSTP(i),i=1,NINTG)
! 1        format(100(I5))
!write(lufo,2,err=100) (DELT(i),i=1,NINTG)
! 2        format(100(E15.6))
!if ((MODE.eq.2 .or. MODE.eq.6).and.(NMPRT.gt.0)) &
!       write(lufo,1,err=100) (IMPRINT(i),i=1,NMPRT)
!if ((MODE.eq.5 .or. MODE.eq.6).and.(NTPRT.gt.0)) &
!       write(lufo,1,err=100) (ITPRINT(i),i=1,NTPRT)
! else
!write(lufo,err=100) (MAXSTP(i),i=1,NINTG)
!write(lufo,err=100) (DELT(i),i=1,NINTG)
!if ((MODE.eq.2 .or. MODE.eq.6).and.(NMPRT.gt.0)) &
!       write(lufo,err=100) (IMPRINT(i),i=1,NMPRT)
!if ((MODE.eq.5 .or. MODE.eq.6).and.(NTPRT.gt.0)) &
!       write(lufo,err=100) (ITPRINT(i),i=1,NTPRT)
! endif

!50    return

!100    call prerr('WritHd')
! ierr = EXEC
! goto 50
! end
!-------------------------------------------------------------------------------
