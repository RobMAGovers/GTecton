 program main

 USE MODELDEFINITION

 implicit double precision (a-h,o-z)
 parameter (NEN=4,NSD=2,NROW=20,NCOL=20)
 parameter (ZERO=0d0)
!-locl
 dimension XL(NSD,NEN),SH(3,NEN),XS(NSD,NSD),ROT(3,3),g(NCOL)
 dimension IEN(4,1)
 iecho = 2
!
 NUMNP = 4
 F = 0.5d0
 angle = 0d0
 XL(1,1) = -F
 XL(2,1) = -F
 XL(1,2) =  F
 XL(2,2) = -F
 XL(1,3) =  F
 XL(2,3) =  F
 XL(1,4) = -F
 XL(2,4) =  F
!
 write(*,10) (j,(XL(i,j),i=1,NSD),j=1,NEN)
10    format(1x,'Input coordinates'/ &
   4(' XL(',I1,')=',2(F12.6)/))
 if (angle.ne.ZERO) then
write(*,20) angle
20        format(1x,'Rotate XL over ',F4.0,' degrees along Z-axis')
call rotmat(ROT,angle)
do i=1,NEN
    call rotvec (XL(1,i),ROT)
enddo
write(*,30) (j,(XL(i,j),i=1,NSD),j=1,NEN)
30        format(4(' XL(',I1,')=',2(F12.6)/))
 endif
 NUMEL = 1
 IEN(1,1) = 1
 IEN(2,1) = 2
 IEN(3,1) = 3
 IEN(4,1) = 3
 call ELMTEST (XL,IEN)

 open(unit=10,file='tst.dat')
 rewind(unit=10)
 dt = 2d0/DBLE(NROW-1)
 ds = 2d0/DBLE(NCOL-1)
!-N1
 write(10,*) '-999 N1'
 write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
 do it=NROW,1,-1
t = -1d0 + DBLE(it-1)*dt
do is=1,NCOL
s = -1d0 + DBLE(is-1)*ds
call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
    g(is) = SH(3,1)
enddo
write(10,*) (g(is),is=1,NCOL)
 enddo
 write(10,*) '-999 dN1/dx'
 write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
 do it=NROW,1,-1
t = -1d0 + DBLE(it-1)*dt
do is=1,NCOL
s = -1d0 + DBLE(is-1)*ds
call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
    g(is) = SH(1,1)
enddo
write(10,*) (g(is),is=1,NCOL)
 enddo
 write(10,*) '-999 dN1/dy'
 write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
 do it=NROW,1,-1
t = -1d0 + DBLE(it-1)*dt
do is=1,NCOL
s = -1d0 + DBLE(is-1)*ds
call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
    g(is) = SH(2,1)
enddo
write(10,*) (g(is),is=1,NCOL)
 enddo
!-N2
 write(10,*) '-999 N2'
 write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
 do it=NROW,1,-1
t = -1d0 + DBLE(it-1)*dt
do is=1,NCOL
s = -1d0 + DBLE(is-1)*ds
call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
    g(is) = SH(3,2)
enddo
write(10,*) (g(is),is=1,NCOL)
 enddo
 write(10,*) '-999 dN2/dx'
 write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
 do it=NROW,1,-1
t = -1d0 + DBLE(it-1)*dt
do is=1,NCOL
s = -1d0 + DBLE(is-1)*ds
call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
    g(is) = SH(1,2)
enddo
write(10,*) (g(is),is=1,NCOL)
 enddo
 write(10,*) '-999 dN2/dy'
 write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
 do it=NROW,1,-1
t = -1d0 + DBLE(it-1)*dt
do is=1,NCOL
s = -1d0 + DBLE(is-1)*ds
call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
    g(is) = SH(2,2)
enddo
write(10,*) (g(is),is=1,NCOL)
 enddo
!-N3
 write(10,*) '-999 N3'
 write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
 do it=NROW,1,-1
t = -1d0 + DBLE(it-1)*dt
do is=1,NCOL
s = -1d0 + DBLE(is-1)*ds
call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
    g(is) = SH(3,3)
enddo
write(10,*) (g(is),is=1,NCOL)
 enddo
 write(10,*) '-999 dN3/dx'
 write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
 do it=NROW,1,-1
t = -1d0 + DBLE(it-1)*dt
do is=1,NCOL
s = -1d0 + DBLE(is-1)*ds
call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
    g(is) = SH(1,3)
enddo
write(10,*) (g(is),is=1,NCOL)
 enddo
 write(10,*) '-999 dN3/dy'
 write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
 do it=NROW,1,-1
t = -1d0 + DBLE(it-1)*dt
do is=1,NCOL
s = -1d0 + DBLE(is-1)*ds
call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
    g(is) = SH(2,3)
enddo
write(10,*) (g(is),is=1,NCOL)
 enddo
!-N4
 if (IEN(3,1).ne.IEN(4,1)) then
write(10,*) '-999 N4'
write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
do it=NROW,1,-1
    t = -1d0 + DBLE(it-1)*dt
    do is=1,NCOL
    s = -1d0 + DBLE(is-1)*ds
    call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
        g(is) = SH(3,4)
    enddo
    write(10,*) (g(is),is=1,NCOL)
enddo
write(10,*) '-999 dN4/dx'
write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
do it=NROW,1,-1
    t = -1d0 + DBLE(it-1)*dt
    do is=1,NCOL
    s = -1d0 + DBLE(is-1)*ds
    call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
        g(is) = SH(1,4)
    enddo
    write(10,*) (g(is),is=1,NCOL)
enddo
write(10,*) '-999 dN4/dy'
write(10,*) 3,NCOL,NROW,1,-1,0.5,1,-1,0.5,1
do it=NROW,1,-1
    t = -1d0 + DBLE(it-1)*dt
    do is=1,NCOL
    s = -1d0 + DBLE(is-1)*ds
    call SHAP20 (s,t,XL,DET,SH,XS,.true.,ierr)
        g(is) = SH(2,4)
    enddo
    write(10,*) (g(is),is=1,NCOL)
enddo
 endif
 close(unit=10)
!
  end
!-------------------------------------------------------------------------------
 SUBROUTINE SHAP20 (S,T,X,DET,SH,XS,QUAD,ierr)

! Program to compute shape functions for quadrilateral or triangle

!    S,T,                = natural coordinates
!    SH(1,NEN),SH(2,NEN) = x and y derivatives of shape functions
!    SH(3,NEN)           = shape functions
!    XS(NSD,NSD)         = inverse jacobian matrix
!    DET                 = Jacobian matrix determinant
!    X(2,NEN)            = global coordinates

 parameter (NEN=4, NSD=2)
! SNGL    implicit real (a-h,o-z)
 implicit double precision (A-H,O-Z)
!-pass
 logical QUAD
 integer ierr
 dimension SH(3,NEN),X(NSD,NEN),XS(NSD,NSD)
!-locl
 dimension SA(NEN),TA(NEN)
 save SA,TA
!-init
 data SA/-0.5d0, 0.5d0,0.5d0,-0.5d0/
 data TA/-0.5d0,-0.5d0,0.5d0, 0.5d0/
 half = 0.5d0
 zero = 0d0
 ierr = 0

!    Calculate shape functions and s- and t-derivatives of shape functions
 do 10 i=1,NEN
SH(1,i) = SA(i)*(half+TA(i)*T)
SH(2,i) = TA(i)*(half+SA(i)*S)
SH(3,i) = (half+SA(i)*S)*(half+TA(i)*T)
10    continue
 if (.not.QUAD) then
!        adjust shape function of third node for triangular elements
do 20 i=1,3
SH(i,3) = SH(i,3)+SH(i,4)
SH(i,4) = zero
20        continue
 endif
!    Jacobian matrix calculation
 do 50 i=1,2
do 40 j=1,2
XS(i,j) = zero
do 30 k=1,NEN
    XS(i,j) = XS(i,j) + SH(i,k)*X(j,k)
30        continue
40        continue
50    continue
 DET = XS(1,1)*XS(2,2)-XS(1,2)*XS(2,1)
 if (DET.le.zero) then 
write (0,*) 'TempSrc/pln/tst.F says:'
go to 1000
 endif
!    Calculatate inverse Jacobian determinant (Cramers rule)
 TEMP    =  XS(1,1)/DET
 XS(1,1) =  XS(2,2)/DET
 XS(2,2) =  TEMP
 XS(1,2) = -XS(1,2)/DET
 XS(2,1) = -XS(2,1)/DET
 do 60 i=1,NEN
TEMP    = XS(1,1)*SH(1,i) + XS(1,2)*SH(2,i)
SH(2,i) = XS(2,1)*SH(1,i) + XS(2,2)*SH(2,i)
SH(1,i) = TEMP
60    continue
70    return
 1000    write(stderr,2000) DET
 2000    format(///1x,'Shape function fails! Determinant is ',1PE20.4)
 ierr = 1
 goto 70
 end
!-------------------------------------------------------------------------------
subroutine rotmat(R,rangle)

USE CONSTANTS, only: pi

implicit none

double precision PI,R,rangle
dimension R(3,3)

double precision sinA, sinB

angle = rangle *PI/1.8d2
SINA = SIN(angle)
COSA = COS(angle)

R(1,1) = COSA
R(2,2) = COSA
R(1,2) = -SINA
R(2,1) = SINA
R(3,3) = 1d0
R(1,3) = 0d0
R(2,3) = 0d0
R(3,1) = 0d0
R(3,2) = 0d0

return
end
!-------------------------------------------------------------------------------
subroutine rotvec (a,R)

implicit none

double precision :: a(2),R(3,3)
double precision :: dimension c(2)

integer :: i

do i=1,2
    c(i) = 0d0
    do j=1,2
        c(i) = c(i) + R(i,j)*a(j)
    enddo
enddo

do i=1,2
    a(i) = c(i)
enddo

return
end
!-------------------------------------------------------------------------------
 SUBROUTINE ELMTEST (X,IEN)

 USE MODELDEFINITION
 use constants, only: eps

!
! tests element definitions
!
 implicit none
 integer NEN,NSD
 parameter (NEN=4,NSD=2)
 double precision ZERO,HALF
 parameter (ZERO=0d0,HALF=5D-1)
!-pass
 integer IEN
 double precision X
 dimension X(NSD,*),IEN(NEN,*)

!-locl
 logical quad
 integer i,k
 double precision xl,area1,area2,area,a,b
 dimension xl(NSD,NEN),a(NSD),b(NSD)
!
 if (NUMEL.lt.1) then
write(stderr,10)
10        format(1x,'NUMEL < 1')
   call exitp(1)
 endif
!
 do i=1,NUMEL
k = IEN(1,i)
if (k.lt.1 .or. k.gt.NUMNP) then
write(stderr,20) i,k
20        format(1x,'element ',I6,': reference to node number ', &
       I6)
call exitp(1)
endif
XL(1,1) = X(1,k)
XL(2,1) = X(2,k)
k = IEN(2,i)
if (k.lt.1 .or. k.gt.NUMNP) then
write(stderr,20) i,k
call exitp(1)
endif
XL(1,2) = X(1,k)
XL(2,2) = X(2,k)
k = IEN(3,i)
if (k.lt.1 .or. k.gt.NUMNP) then
write(stderr,20) i,k
call exitp(1)
endif
XL(1,3) = X(1,k)
XL(2,3) = X(2,k)
k = IEN(4,i)
if (k.lt.1 .or. k.gt.NUMNP) then
write(stderr,20) i,k
call exitp(1)
endif
XL(1,4) = X(1,k)
XL(2,4) = X(2,k)
quad = (IEN(3,i).ne.IEN(4,i))
!
    a(1) = XL(1,2) - XL(1,1)
    a(2) = XL(2,2) - XL(2,1)
b(1) = XL(1,3) - XL(1,2)
b(2) = XL(2,3) - XL(2,2)
area1 = HALF*(a(1)*b(2)-a(2)*b(1))
if (quad) then
        a(1) = XL(1,4) - XL(1,3)
        a(2) = XL(2,4) - XL(2,3)
b(1) = XL(1,1) - XL(1,4)
b(2) = XL(2,1) - XL(2,4)
        area2 = HALF*(a(1)*b(2)-a(2)*b(1))
    area = area1 + area2
else
area = area1
endif
if (ABS(area).le.EPS) then
write(stderr,30) i
30        format(1x,'degenerate element ',I6)
call exitp(1)
endif
if (area.lt.ZERO) then
write(stderr,40) i
40        format(1x,'element ',I6,' is negatively oriented (node', &
       ' order is wrong)')
call exitp(1)
endif
if (ABS(area1).le.EPS) then
write(stderr,50) i
50        format(1x,'WARNING: degenerate triangle 123 in element', &
       1X,I6)
call exitp(1)
endif
if (area1.lt.ZERO) then
write(stderr,60) i
60        format(1x,'WARNING: negatively oriented triangle 123', &
       ' in element ',I6)
call exitp(1)
endif
if (quad .and. ABS(area2).le.EPS) then
write(stderr,70) i
70        format(1x,'WARNING: degenerate triangle 341 in element', &
       1X,I6)
call exitp(1)
endif
if (quad .and. area2.lt.ZERO) then
write(stderr,80) i
80        format(1x,'WARNING: negatively oriented triangle 341', &
       ' in element ',I6)
call exitp(1)
endif
 enddo
!
 if (iecho.eq.2) write(stderr,90)
90    format(1x,'Elements tested ok')
 return
 end
