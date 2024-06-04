	program main

	USE ALGEBRA, only: rotz, rotx, rotvec, clear

	implicit double precision (a-h,o-z)
	parameter (NSD=3,NEN=4,NDOF=3,NSTR=6)
	logical QUAD
	dimension XL(NSD,NEN),PL(NDOF,NEN),GRAV(NDOF)
	dimension ROT(3,3)
*
	data WT /10./
	data GRAV /0d0,0d0,-9.8d0/
F = 0.5d0*6d0**0.3333333333333333333d0
XL(1,1) = -F
XL(2,1) = -F
XL(3,1) = -F
XL(1,2) =  F
XL(2,2) = -F
XL(3,2) = -F
XL(1,3) = -F
XL(2,3) =  F
XL(3,3) = -F
XL(1,4) = -F
XL(2,4) = -F
XL(3,4) =  F
*
	write(*,10) (j,(XL(i,j),i=1,NSD),j=1,NEN)
10	format(1x,'Input data'/
  >	 4(' XL(',I1,')=',3(F12.6)/))
	write(*,25) WT
25	format(1x,'DENSITY=',F12.6)
	write(*,26) (GRAV(i),i=1,NDOF)
26	format(1x,'GRAV= ',3(F12.6,1x))

	angle = 20.
	write(*,29) angle
29	format(1x,'Rotate XL over ',F4.0,' degrees along X-axis')
	call rotx(ROT,angle)
	do i=1,NEN
	    call rotvec (XL(1,i),ROT)
	enddo
	write(*,30) (j,(XL(i,j),i=1,NSD),j=1,NEN)
30	format(4(' XL(',I1,')=',3(F12.6)/))
*
	call CLEAR (PL,NDOF*NEN)
	CALL GRAVLD (PL,XL,GRAV,DUM,WT,QUAD,ierr)
	if (ierr.ne.0) stop 'determinant error'

	write(*,50) (j,(PL(i,j),i=1,NDOF),j=1,NEN)
50	format(1x,'Output data'/
  >   4(' PL(',I1,')=',3(F12.6,1x)/))
	end
