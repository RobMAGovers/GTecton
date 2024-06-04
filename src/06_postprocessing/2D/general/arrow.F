integer function iarrow (mode,X,Y,len,ang,headlen,headang,arrw, &
    narrw,MAXRRW,both)

 USE CONSTANTS, only: pi, deg2rad
use debugmodule, only: iecho
use iomodule,        only: stderr

! calculates coordinates of arrow with center (X, Y) length LEN and
! angle ANG (deg).
! MODE   -1     -    no arrow head
!      0    -    open arrow heads at far ends (default)
!      1    -    open arrow heads near center
!      2    -    contour arrow with heads at far ends
!      3    -    contour arrow with heads at center
! ARANG = top hat angle [deg]
! ARLEN = length of arrow head

 implicit none

!passed
integer :: mode
 double precision len
logical :: both
double precision :: x, y
double precision :: ang, headlen, headang, arrw
integer :: narrw, maxrrw
!local
double precision, parameter :: UP=3D0, DOWN=2D0
 dimension arrw(3,*)

double precision :: sina, cosa
double precision :: sint, cost
double precision :: Xplus, Xmin
double precision :: Yplus, Ymin
double precision :: halflen

logical :: overfl
integer :: i, n0

 iarrow   = 0


 if (mode.lt.-1 .or. mode.gt.3) then
write(stderr,'(''arrow: mode error'')')
goto 100
 endif
 if (iecho.eq.5) then
write(stderr,1) mode,X,Y,len,ang,headlen,headang
 1       format(1x,'iarrow: mode=',i2,', X=',F10.4,', Y=',F10.4, &
        ', length=',F8.2,', angle=',F8.2,' deg'/1x, &
        'iarrow: unit arrow length=',F10.4,', arrow hat angle=', &
        F10.4)
 endif

 halflen = ABS(len) / 2d0

 cosa = COS(ang * deg2rad)
 sina = SIN(ang * deg2rad)
 cost = COS(headang * deg2rad)
 sint = SIN(headang * deg2rad)
 Xplus = headlen * (cosa * cost - sina * sint)
 Xmin  = headlen * (cosa * cost + sina * sint)
 Yplus = headlen * (sina * cost + sint * cosa)
 Ymin  = headlen * (sina * cost - sint * cosa)
 if (halflen .lt. headlen) then
   Xplus = Xplus  *  halflen / headlen
   Xmin  = Xmin   *  halflen / headlen
   Yplus = Yplus  *  halflen / headlen
   Ymin  = Ymin   *  halflen / headlen
 endif

 if (mode.eq.-1 .or. mode.eq.0 .or. mode.eq.2) then
!        no arrow head, open head or contour head at ends
if ( overfl(narrw,MAXRRW) ) goto 100
if (both) then
    arrw(1,narrw) = X - cosa * halflen
    arrw(2,narrw) = Y - sina * halflen
else
    arrw(1,narrw) = X
    arrw(2,narrw) = Y
endif
arrw(3,narrw) = up
    n0 = narrw
if ( overfl(narrw,MAXRRW) ) goto 100
    if (both) then
    arrw(1,narrw) = X + cosa * halflen
    arrw(2,narrw) = Y + sina * halflen
    else
    arrw(1,narrw) = X + cosa * len
    arrw(2,narrw) = Y + sina * len
    endif
arrw(3,narrw) = down
if (mode.eq.0 .or. mode.eq.2) then
!            arrow end plotting
    if ( overfl(narrw,MAXRRW) ) goto 100
    arrw(1,narrw) = arrw(1,n0+1) - Xplus
    arrw(2,narrw) = arrw(2,n0+1) - Yplus
    arrw(3,narrw) = down
    if (mode.eq.0) then
        if ( overfl(narrw,MAXRRW) ) goto 100
        arrw(1,narrw) = arrw(1,n0+1)
        arrw(2,narrw) = arrw(2,n0+1)
        arrw(3,narrw) = down
    endif
    if ( overfl(narrw,MAXRRW) ) goto 100
    arrw(1,narrw) = arrw(1,n0+1) - Xmin
    arrw(2,narrw) = arrw(2,n0+1) - Ymin
    arrw(3,narrw) = down
    if ( overfl(narrw,MAXRRW) ) goto 100
    arrw(1,narrw) = arrw(1,n0+1)
    arrw(2,narrw) = arrw(2,n0+1)
    arrw(3,narrw) = down
        if (both) then
        if ( overfl(narrw,MAXRRW) ) goto 100
        arrw(1,narrw) = arrw(1,n0)
        arrw(2,narrw) = arrw(2,n0)
        arrw(3,narrw) = down
        if ( overfl(narrw,MAXRRW) ) goto 100
        arrw(1,narrw) = arrw(1,n0) + Xplus
        arrw(2,narrw) = arrw(2,n0) + Yplus
        arrw(3,narrw) = down
        if (mode.eq.0) then
            if ( overfl(narrw,MAXRRW) ) goto 100
            arrw(1,narrw) = arrw(1,n0)
            arrw(2,narrw) = arrw(2,n0)
            arrw(3,narrw) = down
        endif
        if ( overfl(narrw,MAXRRW) ) goto 100
        arrw(1,narrw) = arrw(1,n0) + Xmin
        arrw(2,narrw) = arrw(2,n0) + Ymin
        arrw(3,narrw) = down
        if (mode.eq.2) then
            if ( overfl(narrw,MAXRRW) ) goto 100
            arrw(1,narrw) = arrw(1,n0)
            arrw(2,narrw) = arrw(2,n0)
            arrw(3,narrw) = down
        endif
        endif
    endif
 else if (mode.eq.1 .or. mode.eq.3) then
!           arrow head near center (X,Y)
if ( overfl(narrw,MAXRRW) ) goto 100
    if (both) then
    arrw(1,narrw) = X + cosa * halflen
    arrw(2,narrw) = Y + sina * halflen
else
    arrw(1,narrw) = X + cosa * len
    arrw(2,narrw) = Y + sina * len
    endif
arrw(3,narrw) = up
    n0 = narrw
if ( overfl(narrw,MAXRRW) ) goto 100
arrw(1,narrw) = X + cosa * halflen * 0.1d0
arrw(2,narrw) = Y + sina * halflen * 0.1d0
arrw(3,narrw) = down
if ( overfl(narrw,MAXRRW) ) goto 100
arrw(1,narrw) = arrw(1,n0+1) + Xplus
arrw(2,narrw) = arrw(2,n0+1) + Yplus
arrw(3,narrw) = down
if (mode.eq.1) then
    if ( overfl(narrw,MAXRRW) ) goto 100
    arrw(1,narrw) = arrw(1,n0+1)
    arrw(2,narrw) = arrw(2,n0+1)
    arrw(3,narrw) = down
endif
if ( overfl(narrw,MAXRRW) ) goto 100
arrw(1,narrw) = arrw(1,n0+1) + Xmin
arrw(2,narrw) = arrw(2,n0+1) + Ymin
arrw(3,narrw) = down
if (mode.eq.3) then
    if ( overfl(narrw,MAXRRW) ) goto 100
    arrw(1,narrw) = arrw(1,n0+1)
    arrw(2,narrw) = arrw(2,n0+1)
    arrw(3,narrw) = down
endif
    if (both) then
    if ( overfl(narrw,MAXRRW) ) goto 100
    arrw(1,narrw) = X - cosa * halflen
    arrw(2,narrw) = Y - sina * halflen
    arrw(3,narrw) = up
    if ( overfl(narrw,MAXRRW) ) goto 100
    n0 = narrw
    arrw(1,narrw) = X - cosa * halflen * 0.1d0
    arrw(2,narrw) = Y - sina * halflen * 0.1d0
    arrw(3,narrw) = down
    if ( overfl(narrw,MAXRRW) ) goto 100
    arrw(1,narrw) = arrw(1,n0) - Xplus
    arrw(2,narrw) = arrw(2,n0) - Yplus
    arrw(3,narrw) = down
    if (mode.eq.1) then
        if ( overfl(narrw,MAXRRW) ) goto 100
        arrw(1,narrw) = arrw(1,n0)
        arrw(2,narrw) = arrw(2,n0)
        arrw(3,narrw) = down
    endif
    if ( overfl(narrw,MAXRRW) ) goto 100
    arrw(1,narrw) = arrw(1,n0) - Xmin
    arrw(2,narrw) = arrw(2,n0) - Ymin
    arrw(3,narrw) = down
    if (mode.eq.3) then
        if ( overfl(narrw,MAXRRW) ) goto 100
        arrw(1,narrw) = arrw(1,n0)
        arrw(2,narrw) = arrw(2,n0)
        arrw(3,narrw) = down
    endif
    endif
 else
write(stderr,'(''arrow: mode error'')')
iarrow = 1
 endif
10    return
100    iarrow = 1
 goto 10
 end
!-------------------------------------------------------------------------------
integer function ipltrw (EVAL,EVEC,CENTR,factr,hatsize,hatangle, &
                 arrow,narrw,MAXRRW,icont,both)

USE ALGEBRA,         only: clear
USE CONSTANTS,       only: pi
use errors,          only: ok
use plotControl,     only: sysc

implicit none
!-pass
integer narrw,MAXRRW,icont

double precision :: arrow(3, maxrrw), centr
dimension        :: CENTR(*)

double precision :: EVAL(3),EVEC(3,3)
logical          :: both
double precision :: factr
double precision :: hatsize,hatangle
!-locl
double precision :: dlen
integer          :: itype
double precision :: angle
integer          :: ierr, ivec
double precision :: Xorg, Yorg
integer, external :: iarrow
double precision, external :: arctan
!integer :: ipltrw

!-init
ierr  = OK
narrw = 0

call clear(arrow,3*MAXRRW,"arrow")

Xorg = CENTR(1)
Yorg = CENTR(2)

do  ivec=1,2
    if (EVAL(ivec).ne.0d0) then

        angle = ARCTAN(EVEC(1,ivec),EVEC(2,ivec))
        if (angle.gt.-900d0) then
            angle = angle * 1.8D2 / pi
!           set arrow type according to sign of Eigenvalue

            if (EVAL(ivec).gt.0d0) then
!               tension
                if (icont.eq.0) then
                    itype = 0
                else
                    itype = 2
                endif
            else
!               compression
                if (icont.eq.0) then
                    itype = -1
                else
                    itype = -1
                endif
            endif

!           convert Eigenvalue to plot coordinate difference
            if (sysc) then
                dlen = ABS(EVAL(ivec))*factr
            else
                dlen = factr
            endif

            ierr = iarrow (itype,Xorg,Yorg,dlen,angle, &
                           hatsize,hatangle,arrow,narrw,MAXRRW,both)

            if (ierr.ne.OK) then
                goto 30
            endif

        else
            ipltrw=1
            goto 30
        endif

    endif

enddo

30    ipltrw = ierr
return
end function
!-------------------------------------------------------------------------------
logical function overfl (i,max)

implicit none
integer :: i, max

if (i.lt.max) then
    i = i + 1
    overfl = .false.
else
    overfl = .true.
endif

return
end
!-------------------------------------------------------------------------------
 subroutine princip ( stress,n,np,d,v,nrot )
!
! Computes all eigenvalues and eigenvectors of a real symmetric matrix
! A, (in this case the stress matrix)
! which is of size N by N, stored in a physical NP by NP array. D returns
! the eigenvalues of A in its first N elements. V is a matrix with the same
! logical and physical dimensions as A whose columns contain, on output,
! the normalized eigenvectors of A. NROT returns the number of rotations
! which were required. Ref: "Numerical recipies", Press, W.H., Flannery, &
! B.P , Teukolsky, S.A. and Vetterling, W.T. .
!

 USE CONSTANTS, only: zero, one, dfmin

 implicit none
!-pass

integer :: n, np, nrot
double precision :: stress(np,np),d(np),v(np,np)
!-locl
 character(len=1) answer
double precision :: b(6),z(6),a(6,6)
integer          :: i, ip, iq, j
double precision :: c, g, h, s, t, tau, sm
double precision :: theta, tresh

!

!    write(*,*) 'princip says: ', np
!    do i=1,np
!        write(*,*) 'stress: ', stress(i,1:6)
!    enddo

!     Initialize to the identity matrix and copy b into a
 do ip=1,n
do iq=1,n
    v(ip,iq) = ZERO
    a(ip,iq) = stress(ip,iq)
enddo
v(ip,ip) = ONE
 enddo

!    Initialize B and D to the diagonal of A
 do ip=1,n
b(ip) = a(ip,ip)
d(ip) = b(ip)
!        Vector Z  will accumulate terms of the form ta(p,q) as in
!        equation (11.1.14)
z(ip) = ZERO
 enddo

 nrot = 0
 do i=1,50

!        Sum off-diagonal elements
sm = ZERO
do ip=1,n-1
    do iq=ip+1,n
        sm = sm + ABS(a(ip,iq))
    enddo
enddo

!        Case that sum of off-diagonal elements zero branch to end.
if ( sm.le.dfmin ) then
    goto 100 ! computation completed. Return
else
!            write(*,*) 'residue: ', sm, dfmin
endif

if ( i.lt.4 ) then
!            On the first three sweeps (see notes)
    tresh = 0.2d0 * sm / (n*n)
else
    tresh = ZERO
endif

do ip=1,n-1
    do iq=ip+1,n
        g = 1d3*ABS( a(ip,iq) )
!                After four sweeps, skip the rotation if the off-
!                diagonal element is small.
        if ( (i.gt.4).and.(g.le.dfmin) ) then
            a(ip,iq) = ZERO
        else if ( ABS(a(ip,iq)).gt.tresh ) then
            h = d(iq) - d(ip)
            if ( g.le.dfmin ) then
                t = a(ip,iq) / h
            else
!                        Equation (11.1.10)
                theta = 0.5d0 * h / a(ip,iq)
                t = ONE/(ABS(theta)+SQRT(ONE+theta*theta))
                if ( theta.lt.ZERO ) then 
                    t = -t
                endif
            endif
            c        = ONE/SQRT(ONE+t*t)
            s        = t * c
            tau      = s/(ONE+c)
            h        = t * a(ip,iq)
            z(ip)    = z(ip) - h
            z(iq)    = z(iq) + h
            d(ip)    = d(ip) - h
            d(iq)    = d(iq) + h
            a(ip,iq) = ZERO

!                    Case of rotations 1 <= j < p
            do j=1,ip-1
                g = a(j,ip)
                h = a(j,iq)
                a(j,ip) = g - s * (h + g * tau)
                a(j,iq) = h + s * (g - h * tau)
            enddo

!                    Case of rotations p < j < q
            do j=ip+1,iq-1
                g = a(ip,j)
                h = a(j,iq)
                a(ip,j) = g - s * (h + g * tau)
                a(j,iq) = h + s * (g - h * tau)
            enddo

!                    Case of rotations q < j <= N
            do j=iq+1,n
                g = a(ip,j)
                h = a(iq,j)
                a(ip,j) = g - s * (h + g * tau)
                a(iq,j) = h + s * (g - h * tau)
            enddo
            do j=1,n
                g = v(j,ip)
                h = v(j,iq)
                v(j,ip) = g - s * (h + g * tau)
                v(j,iq) = h + s * (g - h * tau)
            enddo
            nrot = nrot + 1
        endif
    enddo
enddo

!        Update D with the sum of ta(p,q) and re-initialize Z.
do ip=1,n
    b(ip) = b(ip) + z(ip)
    d(ip) = b(ip)
    z(ip) = ZERO
enddo
 enddo

 write(*,*) 'princip: 50 iterations should never happen'

 read(*,*) answer
!
100    return
 end
