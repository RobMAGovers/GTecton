c-----------------------------------------------------------------------
	subroutine SHLCRD (XL,QUAD)
	USE ALGEBRA
	USE MODELDEFINITION, only: RADIUS,ISHELL
	USE CONSTANTS, only: pi, deg2rad, zero, half, one,eps
*
c SNGL  implicit real (a-h,o-z)
implicit double precision (A-H,O-Z)
	parameter (NEN=4,NSD=2)
c-pass
dimension XL(NSD,*)
	logical QUAD
c-comm
c-locl
	double precision north
	dimension v(3,NEN),a(3),b(3),xn(3),north(3),east(3)
*
c	Derive local coordinate directions
	xlon = XL(1,1)*DEG2RAD
	ylat = XL(2,1)*DEG2RAD
	v(1,1) = RADIUS*COS(ylat)
	v(2,1) = v(1,1)*SIN(xlon)
	v(1,1) = v(1,1)*COS(xlon)
	v(3,1) = RADIUS*SIN(ylat)
	xlon = XL(1,2)*DEG2RAD
	ylat = XL(2,2)*DEG2RAD
	v(1,2) = RADIUS*COS(ylat)
	v(2,2) = v(1,2)*SIN(xlon)
	v(1,2) = v(1,2)*COS(xlon)
	v(3,2) = RADIUS*SIN(ylat)
	xlon = XL(1,3)*DEG2RAD
	ylat = XL(2,3)*DEG2RAD
	v(1,3) = RADIUS*COS(ylat)
	v(2,3) = v(1,3)*SIN(xlon)
	v(1,3) = v(1,3)*COS(xlon)
	v(3,3) = RADIUS*SIN(ylat)
	if (QUAD) then
	    xlon = XL(1,4)*DEG2RAD
	    ylat = XL(2,4)*DEG2RAD
	    v(1,4) = RADIUS*COS(ylat)
	    v(2,4) = v(1,4)*SIN(xlon)
	    v(1,4) = v(1,4)*COS(xlon)
	    v(3,4) = RADIUS*SIN(ylat)
	endif
	a(1) = v(1,2) - v(1,1)
	a(2) = v(2,2) - v(2,1)
	a(3) = v(3,2) - v(3,1)
	b(1) = v(1,3) - v(1,1)
	b(2) = v(2,3) - v(2,1)
	b(3) = v(3,3) - v(3,1)
	call CROSSP(a,b,xn)
	if (QUAD) then
	    a(1) = v(1,4) - v(1,1)
	    a(2) = v(2,4) - v(2,1)
	    a(3) = v(3,4) - v(3,1)
	    call CROSSP(b,a,a)
	    xn(1) = (xn(1) + a(1))*HALF
	    xn(2) = (xn(2) + a(2))*HALF
	    xn(3) = (xn(3) + a(3))*HALF
	endif
	s = SQRT(DOT(xn,xn,3))
	if (s.lt.EPS) then
	    write(stderr,*) 'SHCRD: zero normal'
	    call stoper()
	endif
	call BMULT (xn,3,ONE/s)
	east(1) = -xn(2)
	east(2) =  xn(1)
	east(3) = ZERO
	s = SQRT(DOT(east,east,2))
	if (s.le.EPS) then
	    write(stderr,1)
 1	    format(1X,'SHCRD: east direction cannot be derived')
	    call stoper()
	endif
	call BMULT (east,2,ONE/s)
	call CROSSP(xn,east,north)
	s = SQRT(DOT(north,north,3))
	if (s.le.EPS) then
	    write(stderr,2)
 2	    format(1x,'SHCRD: north direction cannot be derived')
	    call stoper()
	endif
	call BMULT (north,3,ONE/s)
*
c	Project node coordinates on local east and north directions
	XL(1,1) = DOT(v(1,1),east,2)
	XL(2,1) = DOT(v(1,1),north,3)
	XL(1,2) = DOT(v(1,2),east,2)
	XL(2,2) = DOT(v(1,2),north,3)
	XL(1,3) = DOT(v(1,3),east,2)
	XL(2,3) = DOT(v(1,3),north,3)
	if (QUAD) then
	    XL(1,4) = DOT(v(1,4),east,2)
	    XL(2,4) = DOT(v(1,4),north,3)
	else
	    XL(1,4) = XL(1,3)
	    XL(2,4) = XL(2,3)
	endif
*
return
end
