subroutine wreslv(IEN,X,D,LMF,TFAULT,DXE,TENSOR,name,ISIDE,ierr)

use modeltopology,   only: NSD,NDOF,NSTR,NEN
use iomodule,        only: verbose, stderr, lumon
use monitorVars,     only: elmon, fmon, nnel, monel, monf
use errors,          only: exec

!
! Echoes components of tensor TENSOR, resolved on face ISIDE
!
implicit none
!-pass
integer :: iside
character(len=*) name

integer :: ien, lmf
double precision :: x, d, tfault, dxe, tensor
dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN,*), &
   TFAULT(NDOF,*),DXE(NDOF,NEN,*),TENSOR(NSTR,*)

integer :: ierr

!-local
integer :: n, nn, jside, ios
double precision :: sn, tn
!
if (.not.elmon) then
    return
endif
!
do nn=1,nnel
    n = monel(nn)
    jside = iside

    if (fmon) then
        jside=monf(nn)
    endif

    call resolve (IEN(1,n),X,D,LMF(1,1,n),TFAULT,n, &
                   TENSOR,jside,SN,TN)

    if (VERBOSE) then
        write(luMon,1,err=1100,iostat=ios) jside,n,name,SN,name, TN
 1           format(1x,'On face ',I1,' of element ',I5,': ', &
         'normal ',A11,' = ',1PG15.6/30x,'shear ',A11,' = ', &
         1PG15.6)
    else
        write(luMon,*,err=1100,iostat=ios) SN,TN
    endif
enddo
!
return
!
1100    write(stderr,10) ios,luMon
10   format(1x,'wreslv: write err ',I3,' on monitor logical unit ', &
    I3)
call prerr('wreslv')
ierr = EXEC
return
end subroutine
!-------------------------------------------------------------------------------
subroutine resolve (IEN,X,D,LMF,TFAULT,elemID,TENSOR,ISIDE,SN,TN)

USE ALGEBRA
USE MODELDEFINITION
use modeltopology,   only: NSD,NDOF,NSTR,NEN
use modeldatamodule, only: modeldatactx
use constants,       only: one, useDXE
!
! Computes resolved normal and shear tensor magnitudes SN and TN on side ISIDE
!
implicit none

!-pass
integer :: elemID, iside
double precision :: sn, tn

integer :: ien, lmf
double precision :: x,d,tfault, tensor
dimension IEN(NEN),X(NSD,*),D(NDOF,*),LMF(NDOF,NEN), &
   TFAULT(NDOF,*),TENSOR(NSTR)

!-locl
logical          :: QUAD
integer          :: IIEN(4)
double precision :: DL(3,8),XL(3,8),XN(3),T(3)
double precision :: area
!
call LCOORD (X,XL,IEN)
call LDISP  (DL,D,IEN,NDOF,NEN)
call AddFaultDisplacement (DL,LMF,TFAULT,NDOF,NEN,101)

call ADDSNE (DL,elemID,NDOF,NEN,useDXE)

call REZONE (XL,DL,ONE)
QUAD = (NTYPE.eq.5 .and. IEN(3).ne.IEN(4))
!       Get nodal point numbers of the side
call sidenp (ISIDE,QUAD,IIEN)
!       Calculate surface area and face normal
call FCGEOM (XL,IIEN,AREA,XN)
!
T(1) = TENSOR(1)*XN(1)+TENSOR(4)*XN(2)+TENSOR(5)*XN(3)
T(2) = TENSOR(4)*XN(1)+TENSOR(2)*XN(2)+TENSOR(6)*XN(3)
T(3) = TENSOR(5)*XN(1)+TENSOR(6)*XN(2)+TENSOR(3)*XN(3)
!
SN = DOT(T,XN,NDOF)
T(1) = T(1) - SN*XN(1)
T(2) = T(2) - SN*XN(2)
T(3) = T(3) - SN*XN(3)
TN = SQRT(T(1)*T(1)+T(2)*T(2)+T(3)*T(3))
!
return
end subroutine
