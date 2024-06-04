module subduction

!	NREF(2,NUMREF)          Element# and side# of reference elements
!                               Ordering in up-wind direction assumed.
!	NOTH(2,NUMOTH)          Element# and side# of other side elements
!                               Ordering in up-wind direction assumed.

! used to be subduct.i

implicit none

integer, parameter :: MAXREF=50
integer, parameter :: MAXOTH=80

integer            :: NUMREF,NUMOTH
integer            :: NREF(2,MAXREF),NOTH(2,MAXOTH)

end module



program main

USE MODELDEFINITION
USE MODELDATAMODULE
USE TIMESTEPMODULE,  only: time, maxitime, DELTP,ALFAP,TALFP, &
                           NINTG,MAXIT, &
                           NSTEP,ITWORK,IMWORK,IMTWRK, &
                           NMPRT,IMPRT,NTPRT,ITPRT,NPRMAT,IPRMAT, &
                           timestepdatactx, &
                           outputcontroldatactx
use debugmodule,     only: lustat, itest
use materialsmodule, only: lmat, matdatactx, grav, prpsed
use meshdatamodule,  only: meshdatactx
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: stderr
use plotControl,     only: jShell

implicit none

integer      :: ierr

integer,          allocatable :: nien(:) ! new Ien
double precision, allocatable :: dNew(:), dxNew(:), xNew(:)
integer,          allocatable :: nnsl(:) ! new slip


#ifdef EXPERIMENTAL_ReducedSlipAlloc
stop "Regridding with reduced memory not supported"
#endif


!
!-init
IWORD    = 2
ierr     = 0
NSTEP    = 0
TIME     = 0d0
ISTART   = 0
MAXITIME = 0
JSHELL   = 0
!
!    establish I/O unit numbers and open in- and output files
call setio()
!
write(stderr,1)
 1    format(1x,'regrid not compiled with -DNZ precompiler flag')

!    establish the basic geometric dimensions and element type
call INTCON ()

!    read in system data, establish the grid and define mechanical problem.
call READFE (NSD,NDOF,NSTR,NEN,ierr)
if (ierr.ne.0) then
    stop "failed to read FE input"
endif

call SYDATA (NDOF,NEN)

!    read data relating to the thermal problem
call RDTHRM (.true.,ierr)
if (ierr.ne.0) then
    stop "failed to read thermal input"
endif

call RESTIO (modeldatactx%D, &
             modeldatactx%DXE,modeldatactx%stn,modeldatactx%str,modeldatactx%BTOT, &
             modeldatactx%T,modeldatactx%FTOT,modeldatactx%VPSTR, &
             NDOF,NEN,NSTR,1)
!
write(*,2) NSTEP
 2    format(1x,'Completed reading "RESTART.DAT" (NSTEP=',I6,')')
!
allocate(nien(NEN*NUMEL))
allocate(nnsl(5*NUMEL))
allocate(dnew(NUMNP*NDOF))
allocate(dxnew(NUMEL*NEN*NDOF))
allocate(xnew(NUMNP*NSD))

!subroutine REGRID (IEN,MAT,NSLIP,X,D,DXE,NEWIEN,NEWSLP,DNEW, &
!                   DXNEW,NSELSD,XNEW,BTOT,ID,IWELM,IWTIME,IWINK,WINK,FMAXW,SKEW, &
!                   ISURF,NSD,NDOF,NEN)

call REGRID (meshdatactx%ien, &
             meshdatactx%mat, &
            modeldatactx%nslip, &
             meshdatactx%X, &
            modeldatactx%D, &
            modeldatactx%DXE, &
             nien, nnsl, dnew, dxnew, &
            modeldatactx%NSELSD, &
             XNEW, &
            modeldatactx%BTOT, &
            modeldatactx%ID, &
            modeldatactx%IWELM, &
            modeldatactx%IWTIME, &
            modeldatactx%IWINK, &
            modeldatactx%WINK, &
            modeldatactx%FMAXW, &
            modeldatactx%SKEW, &
            modeldatactx%ISURF, &
             NSD,NDOF,NEN)

!subroutine WRTECIN (IMPRINT,IMATPR,X,IBOND,BOND,LINK,IWINK,WINK, &
!                    SKEW,MAXSTP,DELT,ALFA,PRPMAT,PRPCRP,PRPDIF,PRPLAS,GRAV,PRPSED, &
!                    MAT,IEN,ISELM,ISTIME,STN0,IELNO,ISIDE,PRES,IELSTR,ISSIDE, &
!                    ISTR,STRS,IWELM,IWSIDE,IWTIME,WPRES,NFAULT,FAULT, &
!                    NSLIP,DIFORC,NSELSD,IWINX,WINX,ISURF,FMAXW,NSD,NDOF, &
!                    NSTR,NEN)

call WRTECIN (outputcontroldatactx%IMPRINT, &
              outputcontroldatactx%IMATPR, &
              xnew, &
              modeldatactx%IBOND, &
              modeldatactx%BOND, &
              modeldatactx%LINK, &
              modeldatactx%IWINK, &
              modeldatactx%WINK, &
              modeldatactx%SKEW, &
           timestepdatactx%MAXSTP, &
           timestepdatactx%DELT, &
           timestepdatactx%ALFA, &
                matdatactx%prpmat, &
                matdatactx%prpcrp, &
                matdatactx%prpdif, &
                matdatactx%PRPLAS, &
                           GRAV, &  ! in materials module, but not part of matdatactx
                           PRPsed, &  ! ,,            ,,            ,,
               meshdatactx%mat, &
                           nien, &
              modeldatactx%iselm, &
              modeldatactx%ISTIME, &
              modeldatactx%STN0, &
              modeldatactx%IELNO, &
              modeldatactx%ISIDE, &
              modeldatactx%PRES, &
              modeldatactx%IELSTR, &
              modeldatactx%ISSIDE, &
              modeldatactx%ISTR, &
              modeldatactx%STRS, &
              modeldatactx%IWELM, &
              modeldatactx%IWSIDE, &
              modeldatactx%IWTIME, &
              modeldatactx%WPRES, &
              modeldatactx%NFAULT, &
              modeldatactx%TFAULT, &
                           nnsl, &
              modeldatactx%DIFORC, &
              modeldatactx%NSELSD, &
              modeldatactx%IWINX, &
              modeldatactx%WINX, &
              modeldatactx%ISURF, &
              modeldatactx%FMAXW, &
                           NSD, NDOF, NSTR, NEN)
!
call RESTIO (dnew, dxnew, &
              modeldatactx%stn, &
              modeldatactx%str, &
              modeldatactx%btot, &
              modeldatactx%T, &
              modeldatactx%FTOT, &
              modeldatactx%VPSTR, &
               NDOF,NEN,NSTR,0)

end program
!-------------------------------------------------------------------------------
subroutine REGRID (IEN,MAT,NSLIP,X,D,DXE,NEWIEN,NEWSLP,DNEW, &
                   DXNEW,NSELSD,XNEW,BTOT,ID,IWELM,IWTIME,IWINK,WINK,FMAXW,SKEW, &
                   ISURF,NSD,NDOF,NEN)

USE ALGEBRA,         only: clear, equate
USE MODELDEFINITION, only: numslp, NSLSKEW, nSurf, numel, numnp
use constants,       only: half, one
use subduction,      only: noth, nref, numoth, numref
use iomodule,        only: stderr
!
! Fault regridding program, closely related to routine TSTRGR
! It is assumed that the the lower plate moves in the x direction.
! Reference elements (overriding plate) - ordered in up-wind direction
! (obduction direction).
! Slippery elements on other side of fault (i.e. subducting plate),
! ordered in up-wind direction (subduction direction).
!
! A node is removed upon creation of a slippery node on the down-wind side of
! the other plate. A node is created on the upwind side when a node is no
! longer slippery. The surface node number is used at the upwind side.
!
!    NREF(2,NUMREF)        Element# and side# of reference elements
!                Ordering in up-wind direction assumed.
!    NOTH(2,NUMOTH)        Element# and side# of other side elements
!                Ordering in up-wind direction assumed.
!
implicit none
!-pass

integer :: NSD,NDOF,NEN
integer ::    IEN(NEN,NUMEL)
integer :: NEWIEN(NEN,NUMEL)

integer ::  NSLIP(5,NUMEL)
integer :: NEWSLP(5,NUMEL)

double precision :: DNEW(NDOF,NUMNP)
double precision :: XNEW(NSD, NUMNP)

double precision :: D(NDOF,NUMNP)
double precision :: X(NSD, NUMNP)

double precision :: DXNEW(NDOF, NEN, NUMEL)


integer :: NSELSD(2,*),MAT(*),ID(NDOF,*), IWELM(*),IWTIME(*), &
           IWINK(NDOF,*),ISURF(*)

double precision :: DXE(NDOF,NEN,*), &
                    BTOT(*), &
                    FMAXW(NDOF,*),SKEW(2,*), &
                    WINK(NDOF,*)
!-locl
logical          :: modify,TOX
integer          :: lien(2)
double precision :: DL1(2),DL2(2),xl(2,4),dl(2,4)
integer          :: i, j, k, l, m, n
integer          :: ii, jj, kk, ll
integer          :: iref, iel
integer          :: MUMSLP, node, nnode, nsls
integer          :: IWU, IWL
!
! Application dependencies:
IWU = 1          ! dofs weight of upper plate
TOX = .FALSE.    ! motion direction of lower plate

!     derived
IWL = -IWU
!
! LWI
call EQUATE (NEWIEN,IEN,NUMEL*NEN)
call CLEAR  (NEWSLP,5*NUMEL,"newslp")
call EQUATE (NEWSLP,NSLIP,NUMSLP*5)

call EQUATE (DNEW,D,NUMNP*NDOF)
call CLEAR  (DXNEW,NUMEL*NEN*NDOF,"dxnew")
call EQUATE (XNEW,X,NUMNP*NSD)

MUMSLP = NUMSLP
nsls = 0
!
!    A non-slippery node has started to move towards the fault.
!    Find the new neighbor element of the last reference element.
n = NREF(1,NUMREF)

call SIDENP (NREF(2,NUMREF),IEN(3,n).ne.IEN(4,n),lien)

if (TOX) then
    i = lien(2)
    j = lien(1)
else
    i = lien(1)
    j = lien(2)
endif

m = 1

do while (m.le.NUMOTH)
    k = NOTH(1,m)
    call SIDENP (NOTH(2,m),IEN(3,k).ne.IEN(4,k),lien)
    if (TOX) then
        if (IEN(lien(1),k).eq.IEN(i,n) .and. &
            IEN(lien(2),k).eq.IEN(j,n)) then
            goto 1000
        endif
    else
        if (IEN(lien(1),k).eq.IEN(j,n) .and. &
            IEN(lien(2),k).eq.IEN(i,n)) then
            goto 1000
        endif
    endif
    m = m + 1
enddo

write(stderr,30)
30    format(1x,'Overriding plate trench element has no neighbor')
stop "Overriding plate trench element has no neighbor"
!
1000    m = m - 1    ! set to index of element that is not part of fault

k = NOTH(1,m)
call SIDENP (NOTH(2,m),IEN(3,k).ne.IEN(4,k),lien)

if (TOX) then
    l = lien(2)
else
    l = lien(1)
endif

node = IEN(l,k)
write(*,40) node
40    format(1x,'set node ',I12,' to become slippery')
write(*,50) NREF(2,NUMREF),n,NOTH(2,m),k
50    format(1x,'set side ',I1,' of element ',I12,' to bound side ',I1, &
   ' of element ',I12)

if (NSLSKEW.gt.0) then
    nsls = nsls + 1
    NSELSD(1,nsls) = n
    NSELSD(2,nsls) = NREF(2,NUMREF)
    if (m.gt.1) then
        nsls = nsls + 1
        NSELSD(1,nsls) = NOTH(1,m-1)
        NSELSD(2,nsls) = NOTH(2,m-1)
        nsls = nsls + 1
        NSELSD(1,nsls) = NOTH(1,m-1)
        NSELSD(2,nsls) = NOTH(2,m-1)
    endif
    nsls = nsls + 1
    NSELSD(1,nsls) = k
    NSELSD(2,nsls) = NOTH(2,m)
endif

!    define new D and DL for node IEN(j,n)
call SETDDX (IEN,X,D,DXE,n,j,k,l,DL1,DL2)

do ii=1,NSD
    DNEW(ii,IEN(j,n)) = (DL1(ii)+DL2(ii))*HALF
enddo

write(*,54) 1,IEN(j,n),(DL1(ii),ii=1,NSD)
54    format(1x,'D',I1,'(',I12,')=',2(1PG14.6))
write(*,54) 2,IEN(l,k),(DL2(ii),ii=1,NSD)
write(*,55) IEN(j,n),(D(ii,IEN(j,n)),ii=1,NDOF), &
   (DNEW(ii,IEN(j,n)),ii=1,NDOF)
55    format(1x,'D(',I12,'): ',2(1PG14.6),' -> ',2(1PG14.6))
write(*,60) node,IEN(j,n)
60    format(1x,'convert references to node ',I12,' to node ',I12)

!    convert references to NODE in NEWIEN arrays
do ii=1,NUMEL
    modify = .false.
    do jj=1,NEN
        if (IEN(jj,ii).eq.node) then
            modify = .true.
            NEWIEN(jj,ii) = IEN(j,n)
            do kk=1,NSD
                DXNEW(kk,jj,ii) = (DL2(kk)-DL1(kk))*HALF
            enddo
            write(*,65) jj,ii,(DXE(kk,jj,ii),kk=1,NDOF), &
                       (DXNEW(kk,jj,ii),kk=1,NDOF)
65            format(1x,'DXE(',I1,',',I12,'): ',2(1PG14.6),' -> ', &
           2(1PG14.6))
        endif
    enddo

    if (modify) then
        MUMSLP = MUMSLP + 1
        NEWSLP(1,MUMSLP) = ii
        NEWSLP(2,MUMSLP) = IEN(j,n)
        NEWSLP(3,MUMSLP) = IWL
        write(*,70) ii,(IEN(jj,ii),jj=1,4), &
           (NEWIEN(jj,ii),jj=1,4)
70        format(1x,' NIEN(',I12,'): ',4I12,' -> ',4I12)
        write(*,80) MUMSLP,(NEWSLP(jj,MUMSLP),jj=1,5)
80        format(1x,'New slippery node entry ',I12,': ',5I12)
    endif
enddo

!    modify differential displacements on reference plate

do ii=1,NUMSLP
    if (NSLIP(2,ii).eq.IEN(j,n).and.NSLIP(3,ii).eq.IWU) then
!        determine local elno of entry
        iel = NSLIP(1,ii)
        do kk=1,NEN
            if (IEN(kk,iel).eq.IEN(j,n)) then
                do ll=1,NSD
                    DXNEW(ll,kk,iel)=(DL1(ll)-DL2(ll))*HALF
                enddo
                write(*,65) kk,iel,(DXE(ll,kk,iel),ll=1,NDOF), &
                           (DXNEW(ll,kk,iel),ll=1,NDOF)
            endif
        enddo
    endif
enddo
!
!     Loop reference elements in down-dip direction
do iref=NUMREF-1,1,-1
    n = NREF(1,iref)
    call SIDENP (NREF(2,iref),IEN(3,n).ne.IEN(4,n),lien)

    if (TOX) then
        i = lien(2)
        j = lien(1)
    else
        i = lien(1)
        j = lien(2)
    endif

    m = m + 1
    k = NOTH(1,m)

    call SIDENP (NOTH(2,m),IEN(3,k).ne.IEN(4,k),lien)

    if (TOX) then
        l = lien(2)
    else
        l = lien(1)
    endif

    write(*,50) NREF(2,iref),n,NOTH(2,m),k

    if (NSLSKEW.gt.0) then
        nsls = nsls + 1
        NSELSD(1,nsls) = n
        NSELSD(2,nsls) = NREF(2,iref)
        nsls = nsls + 1
        NSELSD(1,nsls) = NOTH(1,m-1)
        NSELSD(2,nsls) = NOTH(2,m-1)
        nsls = nsls + 1
        NSELSD(1,nsls) = k
        NSELSD(2,nsls) = NOTH(2,m)
    endif

!        define new D and DL for node IEN(j,n)
    call SETDDX (IEN,X,D,DXE,n,j,k,l,DL1,DL2)

    do ii=1,NSD
        DNEW(ii,IEN(j,n)) = (DL1(ii)+DL2(ii))*HALF
    enddo

    write(*,54) 1,IEN(j,n),(DL1(ii),ii=1,NSD)
    write(*,54) 2,IEN(l,k),(DL2(ii),ii=1,NSD)
    write(*,55) IEN(j,n),(D(ii,IEN(j,n)),ii=1,NDOF), &
               (DNEW(ii,IEN(j,n)),ii=1,NDOF)
    write(*,60) IEN(l,k),IEN(j,n)

!        convert references to IEN(l,k) in NEWSLP and NEWIEN arrays

    do ii=1,NUMSLP
        iel = NSLIP(1,ii)

        if (NSLIP(2,ii).eq.IEN(l,k).and.NSLIP(3,ii).eq.IWL) then
!            other side of slippery fault
            NEWSLP(2,ii) = IEN(j,n)
            write(*,90) ii,(NSLIP(jj,ii),jj=1,4), &
               (NEWSLP(jj,ii),jj=1,4)
90            format(1x,'NSLIP(',I12,'): ',4I12,' -> ',4I12)
            do jj=1,NEN
                if (IEN(jj,iel).eq.IEN(l,k)) then
                     NEWIEN(jj,iel) = IEN(j,n)
                    do kk=1,NSD
                        DXNEW(kk,jj,iel)=(DL2(kk)-DL1(kk))*HALF
                    enddo
                    write(*,65) jj,iel,(DXE(kk,jj,iel), kk=1,NDOF),  &
                               (DXNEW(kk,jj,iel),kk=1,NDOF)
                endif
            enddo
            write(*,70) iel,(IEN(jj,iel),jj=1,4), &
                       (NEWIEN(jj,iel),jj=1,4)
        endif

        if (NSLIP(2,ii).eq.IEN(j,n).and.NSLIP(3,ii).eq.IWU) then
!            reference side of slippery fault
            do jj=1,NEN
                if (IEN(jj,iel).eq.IEN(j,n)) then
                    do kk=1,NSD
                        DXNEW(kk,jj,iel)=(DL1(kk)-DL2(kk))*HALF
                    enddo
                    write(*,65) jj,iel,(DXE(kk,jj,iel), kk=1,NDOF), &
                                (DXNEW(kk,jj,iel),kk=1,NDOF)
                endif
            enddo
        endif
    enddo

enddo

if (m.lt.NUMOTH) then
    nsls = nsls + 1
    NSELSD(1,nsls) = NOTH(1,m+1)
    NSELSD(2,nsls) = NOTH(2,m+1)
endif

if (TOX) then
    l = lien(1)
else
    l = lien(2)
endif

call SETDDX (IEN,X,D,DXE,n,i,k,l,DL1,DL2)

do ii=1,NSD
    DNEW(ii,IEN(i,n)) = (DL1(ii)+DL2(ii))*HALF
enddo

write(*,54) 1,IEN(j,n),(DL1(ii),ii=1,NSD)
write(*,54) 2,IEN(l,k),(DL2(ii),ii=1,NSD)
write(*,55) IEN(i,n),(D(ii,IEN(i,n)),ii=1,NDOF), &
           (DNEW(ii,IEN(i,n)),ii=1,NDOF)
write(*,60) IEN(l,k),IEN(i,n)
!    convert references to IEN(l,k) in NEWSLP and NEWIEN arrays

do ii=1,NUMSLP
    iel = NSLIP(1,ii)
    if (NSLIP(2,ii).eq.IEN(l,k).and.NSLIP(3,ii).eq.IWL) then
!        other side of slippery fault
        NEWSLP(2,ii) = IEN(i,n)
        write(*,90) ii,(NSLIP(jj,ii),jj=1,4), &
                   (NEWSLP(jj,ii),jj=1,4)
        do jj=1,NEN
            if (IEN(jj,iel).eq.IEN(l,k)) then
                NEWIEN(jj,iel) = IEN(i,n)
                do kk=1,NSD
                    DXNEW(kk,jj,iel)=(DL2(kk)-DL1(kk))*HALF
                enddo
                write(*,65) jj,iel,(DXE(kk,jj,iel),kk=1,NDOF), &
                           (DXNEW(kk,jj,iel),kk=1,NDOF)
            endif
        enddo
        write(*,70) iel,(IEN(jj,iel),jj=1,4), &
                   (NEWIEN(jj,iel),jj=1,4)
    endif

    if (NSLIP(2,ii).eq.IEN(i,n).and.NSLIP(3,ii).eq.IWU) then
!        reference side of slippery fault
        do jj=1,NEN
            if (IEN(jj,iel).eq.IEN(i,n)) then
                do kk=1,NSD
                    DXNEW(kk,jj,iel) = (DL1(kk)-DL2(kk))*HALF
                enddo
                write(*,65) jj,iel,(DXE(kk,jj,iel),kk=1,NDOF), &
                           (DXNEW(kk,jj,iel),kk=1,NDOF)
            endif
        enddo
    endif
enddo

m = m + 1
k = NOTH(1,m)

call SIDENP (NOTH(2,m),IEN(3,k).ne.IEN(4,k),lien)

if (TOX) then
    l = lien(1)
else
    l = lien(2)
endif

! todo this is is version dependent. 
! used only for pln, now.
call LCOORD (X,xl,IEN(1,k),0)
call LDISP  (dl,D,IEN(1,k),NDOF,NEN)
call ADDSNE (dl,DXE(1,1,k),NDOF,NEN)
call REZONE (xl,dl,ONE)
call EQUATE (XNEW(1:NSD,node),xl(1:NSD,l),NSD)
call CLEAR  (DNEW(1:NSD,node),NSD,"NSD")

nnode = IEN(l,k)
write(*,100) nnode
  100    format(1x,'unset node ',I3,' on other plate of being slippery')
write(*,60) nnode,node
!    convert references to IEN(l,k) in NEWSLP and NEWIEN arrays
do ii=1,NUMSLP
!        select other side of slippery fault
    if (NSLIP(3,ii).eq.IWL) then
        if (NSLIP(2,ii).eq.nnode) then
            NEWSLP(1,ii) = -999
            write(*,90) ii,(NSLIP(jj,ii),jj=1,4), &
                       (NEWSLP(jj,ii),jj=1,4)
            iel = NSLIP(1,ii)
            do jj=1,NEN
                if (IEN(jj,iel).eq.nnode) then
                    NEWIEN(jj,iel) = node
                endif
            enddo
            write(*,70) iel,(IEN(jj,iel),jj=1,4), &
                       (NEWIEN(jj,iel),jj=1,4)
        endif
    endif
enddo

write(*,110) node
  110    format(1x,'Nullify global loads on node ',I12)

do ii=1,NDOF
    jj = ID(ii,node)
    if (jj.gt.0) then
        BTOT(jj) = 0d0
    endif
enddo

write(*,120) node
  120    format(1x,'Remove entry for node ',I12,' from surface array')

ii = 1
do while (node.ne.ISURF(ii) .and. ii.le.NSURF)
    ii = ii + 1
enddo

if (ii.gt.NSURF) then
    stop 'surface node entry not found'
endif

do while (ii.lt.NSURF)
    ISURF(ii) = ISURF(ii+1)
    ii = ii + 1
enddo

NSURF = NSURF - 1
!    remove spurious entries from NEWSLP
ii = 1
do while (ii.le.MUMSLP)
    if (NEWSLP(1,ii).eq.-999) then
        call SHIFT (NEWSLP(1,ii),5,MUMSLP-ii+1)
        MUMSLP = MUMSLP - 1
    else
        ii = ii + 1
    endif
enddo
!
NUMSLP  = MUMSLP
NSLSKEW = nsls
!
!    Set Winkler pressures
!        n = NREF(1,NUMREF)
!       call SIDENP (NREF(2,NUMREF),IEN(3,n).ne.IEN(4,n),lien)
!       i = lien(1)
!       j = lien(2)
!       m = 1
!       do while (m.le.NUMOTH)
!           k = NOTH(1,m)
!           call SIDENP (NOTH(2,m),IEN(3,k).ne.IEN(4,k),lien)
!           if (IEN(lien(1),k).eq.IEN(j,n) .and. &
!            IEN(lien(2),k).eq.IEN(i,n)) goto 2000
!           m = m + 1
!       enddo
!2000    m = m - 1
!    do ii=m,NUMOTH
!        k = NOTH(1,ii)
!        do l=1,NUMWNK
!        if (IWELM(l).eq.k) IWTIME(l) = 0
!        enddo
!    enddo
!*
!    call MAXPRES (NEWIEN,XNEW,DNEW,DXNEW,IWINK,WINK,FMAXW,SKEW,NEN,
!     >     NSD,NDOF)
!
return
end subroutine
!-------------------------------------------------------------------------------
subroutine MAXPRES (IEN,X,D,DXE,IWINK,WINK,FMAXW,SKEW,NEN,NSD, &
                    NDOF)

USE MODELDEFINITION
USE CONSTANTS,       only: pi, halfpi
use subduction,      only: NUMOTH, noth, nref
!
implicit none
!-pass
integer          :: NEN,NSD, &
                    NDOF
integer          :: IEN(NEN,*),IWINK(NDOF,*)
double precision :: X(NSD,*),D(NDOF,*), &
                    WINK(NDOF,*),FMAXW(NDOF,*),SKEW(2,*),DXE(NDOF,NEN,*)

!-locl
integer          :: lien(2)
double precision :: xl(24),dl(24),xn(3)
integer          :: j, k, n, node, nodec
double precision :: angle, presn, relax, area

if (NWINK.eq.0) then
    return
endif

!    Set Winkler/Maxwell boundary condition on top of subducted plate
k = NREF(1,1)

call SIDENP (NREF(2,1),IEN(3,k).ne.IEN(4,k),lien)

nodec = IEN(lien(1),k)
n = NUMOTH
j = lien(2)

do while (n.ge.1 .and. IEN(j,k).ne.nodec)
    k = NOTH(1,n)

    call SIDENP (NOTH(2,n),IEN(3,k).ne.IEN(4,k),lien)

    j = lien(1)

    call LCOORD (X,xl,IEN(1,k),0)
    call LDISP  (dl,D,IEN(1,k),NDOF,NEN)
    call ADDSNE (dl,DXE(1,1,k),NDOF,NEN)
    call REZONE (xl,dl,1d0)
    call FCGEOM (XL,lien,area,xn)

    angle = HALFPI - ATAN2(xn(2),xn(1))
    presn = (area*5d10)/2d0
    relax = presn/2d21

    if (IEN(j,k).ne.nodec) then
        IWINK(1,IEN(j,k)) = 0
        IWINK(2,IEN(j,k)) = 1
        WINK(2,IEN(j,k)) = presn
        FMAXW(2,IEN(j,k)) = relax
        SKEW(1,IEN(j,k)) = angle
    endif

    j = lien(2)

    if (n.eq.NUMOTH) then
        if (IWINK(2,IEN(j,k)).ne.0) then
            WINK(2,IEN(j,k)) = presn
            if (NMAXW.gt.0) then
                FMAXW(2,IEN(j,k)) = relax
            endif
        endif
    else
        if (IWINK(2,IEN(j,k)).ne.0) then
            WINK(2,IEN(j,k)) = WINK(2,IEN(j,k)) + presn
            if (NMAXW.gt.0) then
                FMAXW(2,IEN(j,k)) = &
                   FMAXW(2,IEN(j,k)) + relax
            endif
        endif
        SKEW(1,IEN(j,k))=(SKEW(1,IEN(j,k))+angle)/2d0
    endif

    n = n - 1
    j = lien(1)

enddo

n = NUMOTH
j = lien(2)

do while (n.ge.1 .and. IEN(j,k).ne.nodec)
    k = NOTH(1,n)

    call SIDENP (NOTH(2,n),IEN(3,k).ne.IEN(4,k),lien)

    j = lien(2)
    node = IEN(j,k)

    write(*,10) node,SKEW(1,node),WINK(2,node),FMAXW(2,node)
10        format(1x,'node ',I3,': SKEW=',1PG14.6,' WINK=',1PG14.6, &
       ' FMAXW=',1PG14.6)

    n = n - 1
    j = lien(1)
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine SHIFT(IA,N,M)

implicit none

!-pass
integer          :: n, m
integer          :: IA(n,m)
!-local
integer          :: i, j

do j=2,M
    do i=1,N
        IA(i,j-1) = IA(i,j)
    enddo
enddo

return
end subroutine
!-------------------------------------------------------------------------------
! todo for the love of Knuth, wat happened here...
subroutine WRTECIN (IMPRINT,IMATPR,X,IBOND,BOND,LINK,IWINK,WINK, &
                    SKEW,MAXSTP,DELT,ALFA,PRPMAT,PRPCRP,PRPDIF,PRPLAS,GRAV,PRPSED, &
                    MAT,IEN,ISELM,ISTIME,STN0,IELNO,ISIDE,PRES,IELSTR,ISSIDE, &
                    ISTR,STRS,IWELM,IWSIDE,IWTIME,WPRES,NFAULT,FAULT, &
                    NSLIP,DIFORC,NSELSD,IWINX,WINX,ISURF,FMAXW,NSD,NDOF, &
                    NSTR,NEN)

USE MODELDEFINITION
USE CONSTANTS,       only: pi, rad2deg, year
use timestepmodule,  only: DELTP,ALFAP,TALFP, NINTG,MAXIT, &
                           NSTEP,ITWORK,IMWORK,IMTWRK,REFACT, &
                           NMPRT,IMPRT,NTPRT,ITPRT,NPRMAT,IPRMAT
use materialsmodule, only: lmat

!

implicit none

!-pass
integer :: NSD,NDOF, NSTR,NEN

integer :: IMPRINT(*),IMATPR(*),IBOND(NDOF,*), &
           LINK(3,*),IWINK(NDOF,*), MAXSTP(*), &
           MAT(*),IEN(NEN,*), &
           ISELM(*),ISTIME(*),IELNO(*),ISIDE(*), &
           IELSTR(*),ISSIDE(*),ISTR(2,*), &
           IWELM(*),IWSIDE(*),IWTIME(*), &
           NFAULT(3,*), NSLIP(5,*), &
           NSELSD(2,*),IWINX(NDOF,*),ISURF(*)

double precision :: X(NSD,*), &
                    BOND(NDOF,*),WINK(NDOF,*),SKEW(2,*), &
                    DELT(*),ALFA(*),PRPMAT(6,*),PRPCRP(6,*),PRPLAS(9,*), &
                    PRPDIF(9,*),GRAV(*),PRPSED(2), &
                    STN0(NSTR,*), PRES(*), &
                    STRS(NSTR,*),WPRES(*), &
                    FAULT(NDOF,*),DIFORC(NDOF,*), &
                    WINX(NDOF,*), &
                    FMAXW(NDOF,*)

!-locl
character(len=5)  :: unit(50)
integer           :: i,n,lu
logical           :: constraint
integer, external :: nextlu
!
if (NMPRT.lt.0) then
    ISTART = NMPRT
else
    ISTART = -500
endif

lu = nextlu(0)
call openf(lu,'tecin.dat_regrid','unknown')

write(lu,10) NUMNP,MODE,NINTG,NMPRT,MAXIT,NLINK,ICVIS,ISTART, &
             NPRMAT
10    format('automatically generated by "regrid"'/I12/8I12)

if (NMPRT.le.0) then
    write(lu,20) 0
20        format(51I5)
else
    write(lu,20) (IMPRINT(i),i=1,NMPRT)
endif

if (NPRMAT.gt.0) then
    write(lu,20) (IMATPR(i),i=1,NPRMAT)
endif

do n=1,NUMNP
    write(lu,30) n,0,(X(i,n),i=1,NSD)
30        format(2I12,3(1PE14.6))
enddo

write(lu,40)
40    format('end nodal coordinates')

do n=1,NUMNP
    constraint = .false.
    do i=1,NDOF
        constraint = (constraint .or. IBOND(i,n).ne.0)
    enddo
    if (constraint) then
        write(lu,20) n,0,(IBOND(i,n),i=1,NDOF)
    endif
enddo

write(lu,50)
50    format('end bc codes')

do n=1,NUMNP
    constraint = .false.
    do i=1,NDOF
        constraint = (constraint .or. BOND(i,n).ne.0d0)
    enddo
    if (constraint) then
        write(lu,30) n,0,(BOND(i,n),i=1,NDOF)
    endif
enddo

write(lu,60)
60    format('end bc magnitudes')
if (NLINK.gt.0) then
    do i=1,NLINK
        write(lu,20) LINK(1,i),LINK(2,i),LINK(3,i)
    enddo
    write(lu,70)
70        format('end linked nodes')
endif

if (NWINK.gt.0) then
    do n=1,NUMNP
        constraint = .false.
        do i=1,NDOF
            constraint = (constraint .or. IWINK(i,n).ne.0)
        enddo
        if (constraint) then
            write(lu,20) n,0,(IWINK(i,n),i=1,NDOF)
        endif
    enddo
endif

write(lu,80)
80    format('end nodal Winkler codes')

if (NWINK.gt.0) then
    do n=1,NUMNP
        constraint = .false.
        do i=1,NDOF
            constraint = (constraint .or. WINK(i,n).ne.0d0)
        enddo
        if (constraint) then
            write(lu,30) n,0,(WINK(i,n),i=1,NDOF)
        endif
    enddo
endif

write(lu,90)
90    format('end nodal Winkler magnitudes')

if (NMAXW.gt.0) then
    do n=1,NUMNP
        constraint = .false.
        do i=1,NDOF
            constraint = (constraint .or. FMAXW(i,n).ne.0d0)
        enddo
        if (constraint) then
            write(lu,30) n,0,(FMAXW(i,n),i=1,NDOF)
        endif
    enddo
endif

if (NWINK.gt.0) then
    write(lu,92)
92    format('end Maxwell bcs')
endif

if (NUMROT.gt.0) then
    do n=1,NUMNP
        constraint = .false.
        SKEW(1,n) = SKEW(1,n) * rad2deg
        SKEW(2,n) = SKEW(2,n) * rad2deg
        do i=1,2
            constraint = (constraint .or. SKEW(i,n).ne.0d0)
        enddo
        if (constraint) then
            write(lu,30) n,0,(SKEW(i,n),i=1,2)
        endif
    enddo
endif

write(lu,100)
  100    format('end Euler angles')

if (NINTG.gt.0) then
    write(lu,110) (MAXSTP(i),i=1,NINTG)
  110        format(50I5)
    do i=1,NINTG
        DELT(i) = DELT(i)/year
        unit(i) = ' year'
        if (DELT(i).lt.1.0) then
            stop 'DELTi error'
        endif
        if (DELT(i).ge.1d6) then
            DELT(i) = DELT(i)*1d-6
            unit(i) = '   Ma'
        endif
    enddo

    write(lu,110) (NINT(DELT(i)),i=1,NINTG)
    write(lu,120) (UNIT(i),i=1,NINTG)
  120        format(50A5)
    write(lu,130) (ALFA(i),i=1,NINTG)
  130        format(50F5.1)
else
    write(lu,110) 0
    write(lu,110) 0
    write(lu,120) ' '
    write(lu,130) 0.0
endif

write(lu,20) IOPT,IPOINT,LGDEF,IRESDU,IGRAV,IVLIM,INCOMP, &
             NOCOMPR,NSED,ISHELL
write(lu,20) NUMEL,NUMAT,NSURF
write(lu,20) NUMPR,NUMSTR,NUMWNK,NUMSLP,NSLSKEW,NUMFN,NPRE

do n=1,NUMAT
    write(lu,140) n,0,(PRPMAT(i,n),i=1,6)
  140        format(2I12,7(1PE14.6))
enddo

write(lu,150)
  150    format('end basic material property data')

if (NINTG.gt.0.and.ICVIS.ge.1) then
    do n=1,NUMAT
        write(lu,140) n,0,(PRPCRP(i,n),i=1,6)
    enddo
    write(lu,160)
  160        format('end powerlaw creep data')
endif

if (NINTG.gt.0.and.ICVIS.eq.3) then
    do n=1,NUMAT
        write(lu,140) n,0,(PRPDIF(i,n),i=1,7)
    enddo
    write(lu,170)
  170        format('end diffusion creep data')
endif

if (NINTG.gt.0.and.ICVIS.ge.2) then
    do n=1,NUMAT
    write(lu,180) n,0,NINT(PRPLAS(1,n)),NINT(PRPLAS(6,n)), &
       PRPLAS(2,n),PRPLAS(3,n),PRPLAS(4,n),PRPLAS(5,n), &
       PRPLAS(7,n)
  180            format(4I12,5(1PE14.6))
    enddo
    write(lu,190)
  190        format('end plasticity data')
endif

if (NINTG.gt.0.and.NSED.ne.0) then
    write(lu,200) (PRPSED(i),i=1,2)
  200    format(3(1PE14.6))
endif

if (ISHELL.eq.1) then
    write(lu,200) RADIUS
endif

write(lu,200) (GRAV(i),i=1,NDOF)

do n=1,NUMEL
    write(lu,20) n,MAT(n),(IEN(i,n),i=1,NEN)
enddo

write(lu,207)
 207    format('end element definitions')

if (NPRE.ne.0) then
    do n=1,NPRE
        write(lu,210) ISELM(n),0,ISTIME(n),(STN0(i,n),i=1,NSTR)
  210        format(2I12,I6,6(1PE14.6))
    enddo
    write(lu,215)
  215        format('end initial stress data')
endif

if (NUMPR.ne.0) then
    do n=1,IABS(NUMPR)
        write(lu,220) IELNO(n),0,ISIDE(n),PRES(n)
  220        format(3I12,1PE14.6)
    enddo
    write(lu,225)
  225        format('end pressure boundary condition data')
endif

if (NUMSTR.gt.0) then
    do n=1,NUMSTR
        write(lu,230) IELSTR(n),0,ISSIDE(n),ISTR(1,n),ISTR(2,n), &
                       (STRS(i,n),i=1,NSTR)
  230        format(3I12,2I6,6(1PE14.6))
    enddo
    write(lu,235)
  235        format('end stress boundary condition data')
endif

if (NUMWNK.ne.0) then
    do n=1,IABS(NUMWNK)
        write(lu,240) IWELM(n),0,IWSIDE(n),IWTIME(n),WPRES(n)
  240        format(3I12,I6,6(1PE14.6))
    enddo
    write(lu,250)
  250        format('end Winkler pressure data')
endif

if (NUMFN.gt.0) then
    do n=1,NUMFN
            write(lu,260) (NFAULT(i,n),i=1,3),(FAULT(i,n),i=1,NDOF)
  260        format(3I12,3(1PG10.2))
    enddo
    write(lu,270)
  270        format('end split nodes data')
endif

if (NUMSLP.gt.0) then
    do n=1,NUMSLP
        write(lu,280) (NSLIP(i,n),i=1,5),(DIFORC(i,n),i=1,NDOF)
  280        format(5I12,3(1PE10.2))
    enddo
    write(lu,290)
  290        format('end slippery nodes data')
    if (NSLSKEW.gt.0) then
        do n=1,NSLSKEW
            write(lu,20) (NSELSD(i,n),i=1,2)
        enddo
        write(lu,300)
  300        format('end fault parallel elements data')
    endif

    do n=1,NUMNP
        constraint = .false.
        do i=1,NDOF
            constraint = (constraint .or. IWINX(i,n).ne.0)
        enddo
        if (constraint) then
            write(lu,20) n,0,(IWINX(i,n),i=1,NDOF)
        endif
    enddo

    write(lu,310)
  310        format('end slippery Winkler codes')

    do n=1,NUMNP
        constraint = .false.
        do i=1,NDOF
            constraint = (constraint .or. WINX(i,n).ne.0d0)
        enddo
        if (constraint) then
            write(lu,320) n,0,(WINX(i,n),i=1,NDOF)
  320        format(2I12,3(1PE14.6))
        endif
    enddo

    write(lu,330)
  330        format('end slippery Winkler magnitudes')
endif

if (NINTG.gt.0.and.NSURF.gt.0) then
    do n=1,NSURF
        write(lu,20) ISURF(n)
    enddo
    write(lu,340)
  340        format('end surface nodes')
 endif

#ifdef NZ
write(lu,350) NUMREF,NUMOTH
  350    format(2I12)

do n=1,NUMREF
    write(lu,350) NREF(1,n),NREF(2,n)
enddo

write(lu,360)
  360    format('end overriding plate elements')

do n=1,NUMOTH
    write(lu,350) NOTH(1,n),NOTH(2,n)
enddo

write(lu,370)
  370    format('end subducting plate surface elements')
#endif

close(lu)

return
end subroutine



!-------------------------------------------------------------------------------
subroutine SETDDX (IEN,X,D,DXE,NREF,IREF,NOTH,IOTH,DL1,DL2)
!
! computes displacement vector D1 relative to local node IREF of
! element NREF. D2 is the displacement vector of local node IOTH of
! element NOTH relative to node,element IREF,NREF.
!

use constants,     only: one
use modeltopology, only: nsd, ndof, nen

implicit none

!-pass
integer          :: NREF,IREF,NOTH,IOTH
integer          :: IEN(NEN,*)
double precision :: X(NSD,*), D(NDOF,*),DXE(NDOF,NEN,*), &
                    DL1(*),DL2(*)
!-locl
integer          :: i
double precision :: dl(NSD,NEN),xoth(NSD,NEN),xref(NSD,NEN)

call LDISP  (dl,D,IEN(1,NREF),NDOF,NEN)
call ADDSNE (dl,DXE(1,1,NREF),NDOF,NEN)

do i=1,NSD
    DL1(i) = dl(i,IREF)
enddo

!    deformed coordinates of other element
call LCOORD (X,xoth,IEN(1,NOTH),0)
call LDISP  (dl,D,IEN(1,NOTH),NDOF,NEN)
call ADDSNE (dl,DXE(1,1,NOTH),NDOF,NEN)
call REZONE (xoth,dl,ONE)
call LCOORD (X,xref,IEN(1,NREF),0)

do i=1,NSD
    DL2(i) = xoth(i,IOTH) - xref(i,IREF)
enddo

return
end subroutine

!-------------------------------------------------------------------------------

subroutine testit (IWELM,IWTIME,NUMWNK)

implicit none

!-pass
integer :: NUMWNK
integer :: IWELM(*),IWTIME(*)
!-local
integer :: i

do i=1,IABS(NUMWNK)
    write(*,901) IWELM(i),IWTIME(i)
  901        format(3(I12,1x))
enddo

return
end subroutine
