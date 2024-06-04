subroutine CloserElement(CurrentElement, NextElement, xsample, ySample, zSample)

! during a sample run, we are in element n
! n does not contain the sample coordinates
! this subroutine find the neighbour of element n
! which is closest to the sample point.
! As such, the sampling process can 'walk' toward 
! the sample coordinates.

use meshdatamodule, only: meshdatactx
use algebra,        only: twoPointsOnSameSideOfPlane

implicit none

! arguments
integer :: CurrentElement
double precision :: xsample, ySample, zSample

! return
integer :: NextElement

! local
double precision :: x1, y1, z1
double precision :: x2, y2, z2
double precision :: x3, y3, z3
double precision :: x4, y4, z4

double precision :: xCenterElem, yCenterElem, zCenterElem

double precision :: xCenterPlane123, yCenterPlane123, zCenterPlane123
double precision :: xCenterPlane124, yCenterPlane124, zCenterPlane124
double precision :: xCenterPlane134, yCenterPlane134, zCenterPlane134
double precision :: xCenterPlane234, yCenterPlane234, zCenterPlane234

double precision :: VecToSamplePointx, VecToSamplePointy, VecToSamplePointz

double precision :: Span1x, Span1y, Span1z
double precision :: Span2x, Span2y, Span2z
double precision :: Span3x, Span3y, Span3z
double precision :: Span4x, Span4y, Span4z

logical :: matchA123, matchB123, matchC123
logical :: matchA124, matchB124, matchC124
logical :: matchA134, matchB134, matchC134
logical :: matchA234, matchB234, matchC234

double precision :: r ! random number


integer :: side

integer :: v1,v2,v3,v4

side = 0
NextElement = 0

! abbreviate its vertices 
v1 = meshdatactx%IEN(1,CurrentElement)
v2 = meshdatactx%IEN(2,CurrentElement)
v3 = meshdatactx%IEN(3,CurrentElement)
v4 = meshdatactx%IEN(4,CurrentElement)

! and abbreviate coordinates

x1 = meshdatactx%X(1,v1)
y1 = meshdatactx%X(2,v1)
z1 = meshdatactx%X(3,v1)

x2 = meshdatactx%X(1,v2)
y2 = meshdatactx%X(2,v2)
z2 = meshdatactx%X(3,v2)

x3 = meshdatactx%X(1,v3)
y3 = meshdatactx%X(2,v3)
z3 = meshdatactx%X(3,v3)

x4 = meshdatactx%X(1,v4)
y4 = meshdatactx%X(2,v4)
z4 = meshdatactx%X(3,v4)


xCenterElem = (x1 + x2 + x3 + x4) / 4.0d0
yCenterElem = (y1 + y2 + y3 + y4) / 4.0d0
zCenterElem = (z1 + z2 + z3 + z4) / 4.0d0

! vectors from center of element to center of planes
xCenterPlane123 = (x1 + x2 + x3) / 3.0d0 - xCenterElem
yCenterPlane123 = (y1 + y2 + y3) / 3.0d0 - yCenterElem
zCenterPlane123 = (z1 + z2 + z3) / 3.0d0 - zCenterElem

xCenterPlane124 = (x1 + x2 + x4) / 3.0d0 - xCenterElem
yCenterPlane124 = (y1 + y2 + y4) / 3.0d0 - yCenterElem
zCenterPlane124 = (z1 + z2 + z4) / 3.0d0 - zCenterElem

xCenterPlane134 = (x1 + x3 + x4) / 3.0d0 - xCenterElem
yCenterPlane134 = (y1 + y3 + y4) / 3.0d0 - yCenterElem
zCenterPlane134 = (z1 + z3 + z4) / 3.0d0 - zCenterElem

xCenterPlane234 = (x2 + x3 + x4) / 3.0d0 - xCenterElem
yCenterPlane234 = (y2 + y3 + y4) / 3.0d0 - yCenterElem
zCenterPlane234 = (z2 + z3 + z4) / 3.0d0 - zCenterElem

! element center to sample point
VecToSamplePointx = xsample - xCenterElem
VecToSamplePointy = ysample - yCenterElem
VecToSamplePointz = zsample - zCenterElem

! from center to element corners
Span1x = x1 - xCenterElem
Span1y = y1 - yCenterElem
Span1z = z1 - zCenterElem

Span2x = x2 - xCenterElem
Span2y = y2 - yCenterElem
Span2z = z2 - zCenterElem

Span3x = x3 - xCenterElem
Span3y = y3 - yCenterElem
Span3z = z3 - zCenterElem

Span4x = x4 - xCenterElem
Span4y = y4 - yCenterElem
Span4z = z4 - zCenterElem

! Now we have all the vectors, check on which side the element
! the sample point is located

! 1 2 3

matchA123 = twoPointsOnSameSideOfPlane(Span1x,Span1y,Span1z,&
                                       Span2x,Span2y,Span2z,&
                   xCenterPlane123,yCenterPlane123,zCenterPlane123,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)

matchB123 = twoPointsOnSameSideOfPlane(Span1x,Span1y,Span1z,&
                                       Span3x,Span3y,Span3z,&
                   xCenterPlane123,yCenterPlane123,zCenterPlane123,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)

matchC123 = twoPointsOnSameSideOfPlane(Span3x,Span3y,Span3z,&
                                       Span2x,Span2y,Span2z,&
                   xCenterPlane123,yCenterPlane123,zCenterPlane123,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)


! 1 2 4

matchA124 = twoPointsOnSameSideOfPlane(Span1x,Span1y,Span1z,&
                               Span2x,Span2y,Span2z,&
                   xCenterPlane124,yCenterPlane124,zCenterPlane124,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)

matchB124 = twoPointsOnSameSideOfPlane(Span1x,Span1y,Span1z,&
                               Span4x,Span4y,Span4z,&
                   xCenterPlane124,yCenterPlane124,zCenterPlane124,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)

matchC124 = twoPointsOnSameSideOfPlane(Span4x,Span4y,Span4z,&
                               Span2x,Span2y,Span2z,&
                   xCenterPlane124,yCenterPlane124,zCenterPlane124,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)

! 1 3 4

matchA134 = twoPointsOnSameSideOfPlane(Span1x,Span1y,Span1z,&
                               Span3x,Span3y,Span3z,&
                   xCenterPlane134,yCenterPlane134,zCenterPlane134,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)

matchB134 = twoPointsOnSameSideOfPlane(Span1x,Span1y,Span1z,&
                               Span4x,Span4y,Span4z,&
                   xCenterPlane134,yCenterPlane134,zCenterPlane134,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)

matchC134 = twoPointsOnSameSideOfPlane(Span4x,Span4y,Span4z,&
                               Span3x,Span3y,Span3z,&
                   xCenterPlane134,yCenterPlane134,zCenterPlane134,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)


! 2 3 4

matchA234 = twoPointsOnSameSideOfPlane(Span3x,Span3y,Span3z,&
                               Span2x,Span2y,Span2z,&
                   xCenterPlane234,yCenterPlane234,zCenterPlane234,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)

matchB234 = twoPointsOnSameSideOfPlane(Span3x,Span3y,Span3z,&
                               Span4x,Span4y,Span4z,&
                   xCenterPlane234,yCenterPlane234,zCenterPlane234,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)

matchC234 = twoPointsOnSameSideOfPlane(Span4x,Span4y,Span4z,&
                               Span2x,Span2y,Span2z,&
                   xCenterPlane234,yCenterPlane234,zCenterPlane234,&
                  VecToSamplePointx,VecToSamplePointy,VecToSamplePointz)


! side 1: spanned by 1, 2 and 4
if      (matchA124 .and. matchB124 .and. matchC124) then
    side = 1
!    return
! side 2: spanned by 1, 2 and 3
else if (matchA123 .and. matchB123 .and. matchC123) then
    side = 2
else if (matchA134 .and. matchB134 .and. matchC134) then
    side = 3
else if (matchA234 .and. matchB234 .and. matchC234) then
    side = 4
else
    write(*,*) "Something wrong. This is not supposed to happen"
    write(*,*) "The sample point does not appear in any quadrant."
endif

if (side.eq.0) then
    write(*,*) "something went wrong with the sample approach."
    write(*,*) "No side has been selected"
    return
endif

NextElement = meshdatactx%ElementNeighbours(currentElement,side)

if (NextElement.eq.0) then
    write(*,*) "--- WARNING ---------------------------------------"
    write(*,*) "WARNING: something went wrong with the sample approach."
    write(*,*) "During the walk to the sample point, GTECTON is now in element", currentElement
    write(*,*) "The point with sample coordinates:"
    write(*,*) "x:",xSample
    write(*,*) "y:",ySample
    write(*,*) "z:",zSample
    write(*,*) "is the closest to side ", side, "of this element."
    write(*,*) "However, there is no other element on that side!"
    write(*,*) "Stepping algorithm does not where to go."
    write(*,*) "Possibly, the sample walk has run into a hole of the domain."
    write(*,*) "To circumvent this, we will continue from a random new point"
    call random_number(r)
    NextElement = int(floor(meshdatactx%neglobal * r))
    ! can be evaluated to 0, when r is very small.
    if (NextElement.eq.0) then
        NextElement=1
    endif
    write(*,*) "Start a new sample line from element", NextElement
    return
endif

return

end subroutine


!-------------------------------------------------------------------------

subroutine nsmpl (IEN,X,D,LMF,TFAULT, &
   SKEW,XR,YR,ZR,QNTY,M,QS,deformed, idata,lastSampledElement, samplePointFound)

use modeldatamodule, only: modeldatactx
USE MODELDEFINITION
use constants,       only: eps, useDXE
use meshdatamodule,  only: meshdatactx
use debugmodule,     only: NewSample, iecho
use iomodule,        only: stderr
use errors,          only: ok

!
! Routine to sample nodal point quantities on a coordinate
!
! Sample QNTY, over the coordinates XR,YR and ZR, and store the sampled
! data in QS.


!
implicit none


integer, parameter :: NEN=4,NSD=3,NDOF=3
!-pass
integer          :: M
logical          :: deformed
double precision :: XR,YR,ZR  ! the coordinates to be sampled at.
double precision :: QS(M)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif

integer          :: ien, lmf
double precision :: x,d,tfault,skew, qnty
dimension        :: X(NSD,*),IEN(NEN,*),QNTY(M,*), &
D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*)
integer          :: idata
integer          :: lastSampledElement



!-locl
double precision :: vol,w1,w2,w3,w4
double precision :: XL(NSD,NEN),DL(NDOF,NEN),QL(NDOF,NEN)
integer          :: elementID, oldElementID
integer          :: i,k,INSIDE, ierr
double precision :: dist
integer          :: numberOfTries
integer :: iDispVel
logical          :: samplePointFound

!-init
do elementID=1,M
    QS(elementID) = -999d0
enddo

! keep previous samples element, because point is quite possibly close
! to previously sampled element.

elementID = lastSampledElement

numberOfTries = 0

iDispVel = -1
if (idata.eq.1) then
    iDispVel = 0
elseif (idata.eq.2) then
    iDispVel = 1
endif

100    if (elementID.gt.NUMEL) then
    write(*,*) "Trying to sample element ", elementID
    write(*,*) "but the mesh has only ", NUMEL, "elements."
    write(*,*) "Perhaps sample point ", xr, yr, zr," outside domain?"
    return
endif

!write(*,*) "ntries", numberOfTries

if(numberOfTries.gt.1000) then
    write(*,*) "tried to sample point", xr, yr, zr
    write(*,*) "Could not find it in the domain. Giving up..."
    return
else
    numberOfTries = numberOfTries + 1
endif


!write(*,*) "nsmpl calling lcoord"

!localize coordinates from X into Xl
call LCOORD (meshdatactx%X,XL,IEN(1,elementID), elementID)

!write(*,*) "nsmpl finished lcoord"


if (deformed) then
    call LDISP  (DL,D,IEN(1,elementID),NDOF,NEN)
    call AddFaultDisplacement (DL,elementID,TFAULT,NDOF,NEN,103)

    call ADDSNE (DL,elementID,NDOF,NEN,iDispVel)

    call REZONE (XL,DL,1d0)
endif

INSIDE = -1
i = 1
w1  = 0d0
w2  = 0d0
w3  = 0d0
w4  = 0d0
vol = 1d0

!write(*,*) "sampling point", xr, yr, zr

! check whether the sample point is within eps-range of
! a vertex. If so, take the value of this point.
do while (INSIDE.eq.-1 .and. i.le.NEN)
    dist = SQRT((XR-XL(1,i))**2+(YR-XL(2,i))**2+(ZR-XL(3,i))**2)
!    write(*,*) "point" , xr, yr, zr, "close to", XL(:,i)

    if (dist.le.EPS) then
        INSIDE = 0
        if (i.eq.1) w1 = 1d0
        if (i.eq.2) w2 = 1d0
        if (i.eq.3) w3 = 1d0
        if (i.eq.4) w4 = 1d0
    endif
    i = i + 1
enddo

! if it is within an element, return the weights of the values at the vertices,
! to form the interpolated value.

!write(*,*) "checking elementID", elementID

if (INSIDE.eq.-1) then

    call PNPOLY(XL,XR,YR,ZR,vol,w1,w2,w3,w4,INSIDE,IERR,elementID)

    if (IERR.ne.OK) then
        write(stderr,5) elementID
 5        format(1x,'sample: error in element ',I8)
        return
    endif
endif

! it is neither to this point, nor in this element, check the next element.
if (INSIDE.eq.-1) then

    if (NewSample) then
        ! check which side of the element is close to the sample point
        oldElementID = elementID
        call CloserElement(oldElementID,elementID,XR,YR,ZR)
    else
        elementID = elementID + 1
    endif 

!    write(*,*) "try next element ", n, XR,YR,ZR


goto 100
endif

!write(*,*) "found next element! ", elementID
lastSampledElement = elementID
samplePointFound = .true.

w1 = w1/vol
w2 = w2/vol
w3 = w3/vol
w4 = w4/vol

if (iecho.gt.1) then
    if (INSIDE.eq.1) then
        write(stderr,10) elementID,vol,w1,w2,w3,w4
10   format(1x,'sample: coordinates in element ',I8/9x, &
    'element volume = ',1PG16.6/9x,'weights are ',4(1PG16.6))
    elseif (INSIDE.eq.0) then
        write(stderr,20) elementID,vol,w1,w2,w3,w4
20   format(1x,'sample: coordinates along element ',I8/9x, &
    'element volume = ',1PG16.6/9x,'weights are ',4(1PG16.6))
    endif
endif

do i=1,M
    k = IEN(1,elementID)
    QL(i,1) = QNTY(i,k)
    k = IEN(2,elementID)
    QL(i,2) = QNTY(i,k)
    k = IEN(3,elementID)
    QL(i,3) = QNTY(i,k)
    k = IEN(4,elementID)
    QL(i,4) = QNTY(i,k)
enddo

!write(*,*) "nsmpl extra stuff"


if (idata.eq.1) then
    ! call addfault etc
    ! now we have local coordinates, update with faulted displacements
    call AddFaultDisplacement (QL,elementID,TFAULT,NDOF,NEN,103)

    ! not tested: slippery node displacements
    ! LWITODO  addsn or addsne?
    call ADDSNE (QL,elementID,NDOF,NEN,iDispVel)

endif

! now we have local coordinates, update with faulted displacements
! call AddFaultDisplacement (QL,elementID,TFAULT,NDOF,NEN,103)
! not tested: slippery node displacements

!    call ADDSNE (QL,elementID,NDOF,NEN,iDispVel)

!write(*,*) "nsmpl extra stuff done"


do i=1,M
    QS(i) = w1*QL(i,1)+w2*QL(i,2)+w3*QL(i,3)+w4*QL(i,4)
enddo

return
end
!-----------------------------------------------------------------------
subroutine esmpl (IEN,X,D,LMF,TFAULT, &
                  SKEW,XR,YR,ZR,QNTY,M,QS,deformed,lastSampledElement, samplePointFound)

USE MODELDEFINITION
use constants,       only: eps,useDXE
use debugmodule,     only: NewSample, iecho
use meshdatamodule,  only: meshdatactx
use modeltopology,   only: nen, nsd, ndof
use iomodule,        only: stderr
use errors,          only: ok

!
! Routine to sample element quantities on a coordinate
!
implicit none

!-pass
 integer M
 logical deformed
 double precision XR,YR,ZR,QS(M)

#ifdef EXPERIMENTAL_ReducedSlipAlloc
double precision :: DXE(NDOF,NEN,nElemsWithSlipperyNodes)
#else
double precision :: DXE(NDOF,NEN,meshdatactx%neglobal)
#endif


integer :: lastSampledElement
integer :: ien,qnty,lmf
double precision :: x,d,tfault,skew
 dimension X(NSD,*),IEN(NEN,*),QNTY(*), &
 D(NDOF,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*), &
   SKEW(2,*)
logical :: samplePointFound



!-locl
double precision :: XL(NSD,NEN),DL(NDOF,NEN),QL(NDOF,NEN)
 integer n,nold,i,INSIDE,k, ierr
double precision :: dist
double precision :: vol, w1,w2,w3,w4
integer          :: numberOfTries
!-init
do i=1,M
    QS(i) = -999d0
enddo

numberOfTries = 0

n = lastSampledElement

100    if (n.gt.NUMEL) then
    return
endif

if(numberOfTries.gt.1000) then
    write(*,*) "tried to sample point", xr, yr, zr
    write(*,*) "Could not find it in the domain. Giving up..." 
    return
else
    numberOfTries =    numberOfTries + 1
endif


call LCOORD (X,XL,IEN(1,n),lastSampledElement)

if (deformed) then
    call LDISP  (DL,D,IEN(1,n),NDOF,NEN)
    call AddFaultDisplacement (DL,n,TFAULT,NDOF,NEN,104)

    call ADDSNE (DL,n,NDOF,NEN,useDXE)

    call REZONE (XL,DL,1d0)
endif

INSIDE = -1
i = 1

do while (INSIDE.eq.-1 .and. i.le.NEN)
    dist = SQRT((XR-XL(1,i))**2+(YR-XL(2,i))**2+(ZR-XL(3,i))**2)
    if (dist.le.EPS) then
        INSIDE = 0
    endif
    i = i + 1
enddo

if (INSIDE.eq.-1) then
    call PNPOLY(XL,XR,YR,ZR,vol,w1,w2,w3,w4,INSIDE,IERR, i)
    if (IERR.ne.OK) then
        write(stderr,*) "PNPOLY produced error ", ierr
        return
    endif
endif


if (INSIDE.eq.-1) then

    if (NewSample) then
        ! check which side of the element is close to the sample point
        nold = n
        call CloserElement(nold,n,XR,YR,ZR)
    else
        n = n + 1
    endif

    goto 100
endif

!write(*,*) "found it! ", n
lastSampledElement = n
samplePointFound = .true.

if (iecho.gt.1) then
    if (INSIDE.eq.1) then
        write(stderr,10) n
10           format(1x,'sample: coordinate in element ',I8)
    else
        write(stderr,20) n
20           format(1x,'sample: coordinate along element ',I8)
    endif
endif
!
do i=1,M
    k = i + (n-1)*M
    QS(i) = QNTY(k)
enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine PNPOLY (X,XR,YR,ZR,D0,D1,D2,D3,D4,INOUT,IERR, elemID)

use constants,       only: eps,zero,one
use algebra,         only: determinant4x4
use iomodule,        only: stderr
use errors,          only: exec, ok

!
implicit none

integer :: elemID
integer, parameter :: NSD=3,NEN=4
!-pass
integer          :: INOUT,IERR
double precision :: X(NSD,NEN)
double precision :: XR,YR,ZR,D0,D1,D2,D3,D4


!-locl
integer ::      i0,i1,i2,i3,i4,NINT
double precision :: A(4,4)
external     :: NINT
!-init

 IERR = OK
 INOUT = -1
!
! Let the tetrahedron have vertices
!
!         V1 = (x1, y1, z1)
!         V2 = (x2, y2, z2)
!         V3 = (x3, y3, z3)
!         V4 = (x4, y4, z4)
!
! and your test point be
!
!         P = (x, y, z)
!
! Then the point P is in the tetrahedron if following five determinants all
! have the same sign.
!
!              |x1 y1 z1 1|
!         D0 = |x2 y2 z2 1|
!              |x3 y3 z3 1|
!              |x4 y4 z4 1|
!
!              |x  y  z  1|
!         D1 = |x2 y2 z2 1|
!              |x3 y3 z3 1|
!              |x4 y4 z4 1|
!
!              |x1 y1 z1 1|
!         D2 = |x  y  z  1|
!              |x3 y3 z3 1|
!              |x4 y4 z4 1|
!
!              |x1 y1 z1 1|
!         D3 = |x2 y2 z2 1|
!              |x  y  z  1|
!              |x4 y4 z4 1|
!
!              |x1 y1 z1 1|
!         D4 = |x2 y2 z2 1|
!              |x3 y3 z3 1|
!              |x  y  z  1|
!
! Some additional notes:
!
!     * If by chance the D0=0, then your tetrahedron is degenerate (the points
!    are coplanar).
!
!     * If any other Di=0, then P lies on boundary i (boundary i being that
!    boundary formed by the three points other than Vi).
!
!     * If the sign of any Di differs from that of D0 then P is outside
!    boundary i.
!
!     * If the sign of any Di equals that of D0 then P is inside boundary i.
!
!     * If P is inside all 4 boundaries, then it is inside the tetrahedron.
!
!     * As a check, it must be that D0 = D1+D2+D3+D4.
!
!

! *** check determinant of tetrahedron 
 A(1,1) = X(1,1)
 A(1,2) = X(2,1)
 A(1,3) = X(3,1)
 A(1,4) = ONE
 A(2,1) = X(1,2)
 A(2,2) = X(2,2)
 A(2,3) = X(3,2)
 A(2,4) = ONE
 A(3,1) = X(1,3)
 A(3,2) = X(2,3)
 A(3,3) = X(3,3)
 A(3,4) = ONE
 A(4,1) = X(1,4)
 A(4,2) = X(2,4)
 A(4,3) = X(3,4)
 A(4,4) = ONE

 D0 = determinant4x4(A)

if (ABS(D0).le.EPS) then
    write(stderr,10)
10        format(1x,'sample: degenerate tetrahedron')
    write(stderr,*) "All vertices of Element ", elemID, "are in the same plane"
    write(stderr,*) "Vertex 1: ", X(1,1), X(2,1), X(3,1)
    write(stderr,*) "Vertex 2: ", X(1,2), X(2,2), X(3,2)
    write(stderr,*) "Vertex 3: ", X(1,3), X(2,3), X(3,3)
    write(stderr,*) "Vertex 4: ", X(1,4), X(2,4), X(3,4)
    ierr = EXEC
    return
endif

 A(1,1) = XR
 A(1,2) = YR
 A(1,3) = ZR
 D1 = determinant4x4(A)

 if (D1.eq.ZERO) then
INOUT = 0
 endif

 A(1,1) = X(1,1)
 A(1,2) = X(2,1)
 A(1,3) = X(3,1)
 A(2,1) = XR
 A(2,2) = YR
 A(2,3) = ZR
 D2 = determinant4x4(A)


 if (D2.eq.ZERO) then
INOUT = 0
 endif

 A(2,1) = X(1,2)
 A(2,2) = X(2,2)
 A(2,3) = X(3,2)
 A(3,1) = XR
 A(3,2) = YR
 A(3,3) = ZR
 D3 = determinant4x4(A)

 if (D3.eq.ZERO) then
INOUT = 0
 endif

 A(3,1) = X(1,3)
 A(3,2) = X(2,3)
 A(3,3) = X(3,3)
 A(4,1) = XR
 A(4,2) = YR
 A(4,3) = ZR

 D4 = determinant4x4(A)


 if (D4.eq.ZERO) then
INOUT = 0
 endif
!

 i0 = INT(SIGN(1.1d0,D0))

if (D1.ne.ZERO) then
    i1 = INT(SIGN(1.1d0,D1))
else
    i1 = i0
endif

if (D2.ne.ZERO) then
    i2 = INT(SIGN(1.1d0,D2))
else
    i2 = i0
endif

if (D3.ne.ZERO) then
    i3 = INT(SIGN(1.1d0,D3))
else
    i3 = i0
endif

if (D4.ne.ZERO) then
    i4 = INT(SIGN(1.1d0,D4))
else
    i4 = i0
endif

if (i0*i1.gt.0 .and. &
    i0*i2.gt.0 .and. &
    i0*i3.gt.0 .and. &
    i0*i4.gt.0)     then
    if (INOUT.eq.-1) then
        INOUT = 1
    endif
else
    INOUT = -1
endif
!
return
end
!-------------------------------------------------------------------------------
