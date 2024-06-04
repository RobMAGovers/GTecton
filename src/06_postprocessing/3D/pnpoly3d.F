!    program test

 implicit none

 use constants, only: dfmin,dlogmin,dfmax,dlogmax, eps



!-locl
 double precision d0,d1,d2,d3,d4,x(3,4)
!-init
 OK = 0
 EXEC = 1
!
 open(unit=1,file='pnpoly.in')
 rewind(1)
 read(1,*) x(1,1),x(2,1),x(3,1)
 read(1,*) x(1,2),x(2,2),x(3,2)
 read(1,*) x(1,3),x(2,3),x(3,3)
 read(1,*) x(1,4),x(2,4),x(3,4)
 read(1,*) xr,yr,zr
 close(1)
 call PNPOLY (x,xr,yr,zr,d0,d1,d2,d3,d4,inout,ierr,1)
 if (ierr.eq.OK) then
if (inout) 10,20,30
10        write(*,*) 'outside'
goto 40
20        write(*,*) 'edge or rim'
goto 40
30        write(*,*) 'inside'
goto 40
 endif
40    stop
!
 end
!-------------------------------------------------------------------------------
 subroutine PNPOLY (X,XR,YR,ZR,D0,D1,D2,D3,D4,INOUT,IERR, elemID)
!
 use constants, only: eps,zero,one

 implicit none

integer :: elemID
 integer NSD,NEN
 parameter (NSD=3,NEN=4)
!-pass
 integer INOUT,IERR
 double precision X(NSD,NEN),XR,YR,ZR,D0,D1,D2,D3,D4


!-locl
 integer i0,i1,i2,i3,i4
 double precision A(4,4),DET44
!-init
 ierr = OK
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
 D0 = DET44(A)
 if (ABS(D0).le.EPS) then
write(stderr,10)
10        format(1x,'sample: degenerate tetrahedron')
ierr = EXEC
goto 100
 endif
 A(1,1) = XR
 A(1,2) = YR
 A(1,3) = ZR
 D1 = DET44(A)
 if (D1.eq.ZERO) then
INOUT = 0
goto 100
 endif
 A(1,1) = X(1,1)
 A(1,2) = X(2,1)
 A(1,3) = X(3,1)
 A(2,1) = XR
 A(2,2) = YR
 A(2,3) = ZR
 D2 = DET44(A)
 if (D2.eq.ZERO) then
INOUT = 0
goto 100
 endif
 A(2,1) = X(1,2)
 A(2,2) = X(2,2)
 A(2,3) = X(3,2)
 A(3,1) = XR
 A(3,2) = YR
 A(3,3) = ZR
 D3 = DET44(A)
 if (D3.eq.ZERO) then
INOUT = 0
goto 100
 endif
 A(3,1) = X(1,3)
 A(3,2) = X(2,3)
 A(3,3) = X(3,3)
 A(4,1) = XR
 A(4,2) = YR
 A(4,3) = ZR
 D4 = DET44(A)
 if (D4.eq.ZERO) then
INOUT = 0
goto 100
 endif
 if (ABS(D0-D1-D2-D3-D4).gt.10*EPS) then
write(stderr,20)
20        format(1x,'sample program error:')
   write(stderr,'(5(1PE16.6,1X))') D0,D1,D2,D3,D4
write(stderr,'(1PE16.6)') D0-D1-D2-D3-D4
ierr = EXEC
goto 100
 endif
!
 i0 = INT(SIGN(1.1d0,D0))
 i1 = INT(SIGN(1.1d0,D1))
 i2 = INT(SIGN(1.1d0,D2))
 i3 = INT(SIGN(1.1d0,D3))
 i4 = INT(SIGN(1.1d0,D4))
 if (i0*i1.lt.0.or.i0*i2.lt.0.or.i0*i3.lt.0.or.i0*i4.lt.0) then
INOUT = -1
 else
INOUT =  1
 endif
!
100    return
 end
!-------------------------------------------------------------------------------
 double precision function DET44(A)
!
!-pass
 double precision A(4,4)
!-local
 double precision a11,a12,a13,a21,a22,a23,a31,a32,a33,D33,sum
!
 a11 = A(2,2)
 a12 = A(2,3)
 a13 = A(2,4)
 a21 = A(3,2)
 a22 = A(3,3)
 a23 = A(3,4)
 a31 = A(4,2)
 a32 = A(4,3)
 a33 = A(4,4)
 D33 = a11*(a22*a33-a23*a32)-a21*(a12*a33-a13*a32) &
   +a31*(a12*a23-a13*a22)
 sum = A(1,1)*D33
 a11 = A(1,2)
 a12 = A(1,3)
 a13 = A(1,4)
!    a21 = A(3,2)
!    a22 = A(3,3)
!    a23 = A(3,4)
!    a31 = A(4,2)
!    a32 = A(4,3)
!    a33 = A(4,4)
 D33 = a11*(a22*a33-a23*a32)-a21*(a12*a33-a13*a32) &
   +a31*(a12*a23-a13*a22)
 sum = sum - A(2,1)*D33
!    a11 = A(1,2)
!    a12 = A(1,3)
!    a13 = A(1,4)
 a21 = A(2,2)
 a22 = A(2,3)
 a23 = A(2,4)
!    a31 = A(4,2)
!    a32 = A(4,3)
!    a33 = A(4,4)
 D33 = a11*(a22*a33-a23*a32)-a21*(a12*a33-a13*a32) &
   +a31*(a12*a23-a13*a22)
 sum = sum + A(3,1)*D33
!    a11 = A(1,2)
!    a12 = A(1,3)
!    a13 = A(1,4)
!    a21 = A(2,2)
!    a22 = A(2,3)
!    a23 = A(2,4)
 a31 = A(3,2)
 a32 = A(3,3)
 a33 = A(3,4)
 D33 = a11*(a22*a33-a23*a32)-a21*(a12*a33-a13*a32) &
   +a31*(a12*a23-a13*a22)
 DET44 = sum - A(4,1)*D33
!
 return
 end
