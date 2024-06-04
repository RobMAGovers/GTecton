 program main

 implicit double precision (a-h,o-z)


 dimension X(2,3)

 double precision :: w1, w2, w3

 OK = 0
 EXEC = 1

 open(unit=10,file='pnpoly.in')
 rewind(10)
 read(10,*) X(1,1),X(2,1)
 read(10,*) X(1,2),X(2,2)
 read(10,*) X(1,3),X(2,3)
 read(10,*) xr,yr
 close(10)

 call PNPOLY(X,xr,yr,det,w1,w2,w3,INSIDE,ierr,1)

 if (IERR.ne.OK) stop 'pnpoly'

 w1=w1/det
 w2=w2/det
 w3=w3/det

 print*,'det,w1,w2,w3=',det,w1,w2,w3
 print*,'INSIDE=',INSIDE

 end
!-------------------------------------------------------------------------------
 subroutine PNPOLY (X,XR,YR,D0,D1,D2,D3,INOUT,IERR,elemID)

 use constants, only: eps,zero,one

!
 implicit none

integer :: elemID
 integer NSD,NEN
 parameter (NSD=2,NEN=3)
!-pass
 integer           :: INOUT,IERR
 double precision  :: X(NSD,NEN),XR,YR,D0,D1,D2,D3


!-locl
 integer           :: i0,i1,i2,i3,NINT
 double precision  :: A(3,3)
 external NINT
!-init
 IERR = OK
 INOUT = -1
!
! Let the triangle have vertices
!
!         V1 = (x1, y1)
!         V2 = (x2, y2)
!         V3 = (x3, y3)
!
! and your test point be
!
!         P = (x, y, z)
!
! Then the point P is in the triangle if following four determinants all
! have the same sign.
!
!              |x1 y1 1|
!         D0 = |x2 y2 1|
!              |x3 y3 1|
!
!              |x  y  1|
!         D1 = |x2 y2 1|
!              |x3 y3 1|
!
!              |x1 y1 1|
!         D2 = |x  y  1|
!              |x3 y3 1|
!
!              |x1 y1 1|
!         D3 = |x2 y2 1|
!              |x  y  1|
!
! Some additional notes:
!
!     * If by chance the D0=0, then your triangle is degenerate (the points
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
!     * As a check, it must be that D0 = D1+D2+D3.
!
!
 A(1,1) = X(1,1)
 A(1,2) = X(2,1)
 A(1,3) = ONE
 A(2,1) = X(1,2)
 A(2,2) = X(2,2)
 A(2,3) = ONE
 A(3,1) = X(1,3)
 A(3,2) = X(2,3)
 A(3,3) = ONE

 D0 = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
    - A(2,1)*(A(1,2)*A(3,3)-A(1,3)*A(3,2)) &
    + A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2))

 if (ABS(D0).le.EPS) then
write(stderr,10)
10        format(1x,'sample: degenerate triangle')
ierr = EXEC
goto 100
 endif

 A(1,1) = XR
 A(1,2) = YR

 D1 = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
    - A(2,1)*(A(1,2)*A(3,3)-A(1,3)*A(3,2)) &
    + A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2))

 if (D1.eq.ZERO) INOUT = 0

 A(1,1) = X(1,1)
 A(1,2) = X(2,1)
 A(2,1) = XR
 A(2,2) = YR

 D2 = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
    - A(2,1)*(A(1,2)*A(3,3)-A(1,3)*A(3,2)) &
    + A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2))

 if (D2.eq.ZERO) INOUT = 0

 A(2,1) = X(1,2)
 A(2,2) = X(2,2)
 A(3,1) = XR
 A(3,2) = YR

 D3 = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
    - A(2,1)*(A(1,2)*A(3,3)-A(1,3)*A(3,2)) &
    + A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2))

 if (D3.eq.ZERO) INOUT = 0
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

 if (i0*i1.gt.0 .and.  &
i0*i2.gt.0 .and.  &
i0*i3.gt.0)          then
if (INOUT.eq.-1) INOUT = 1
 else
INOUT = -1
 endif

!
100    return
 end
