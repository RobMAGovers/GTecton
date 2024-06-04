	program main

	USE ALGEBRA, only: rotz, rotx, rotvec, clear


	implicit double precision (a-h,o-z)
	parameter (NSD=3,NEN=4,NDOF=3,NSTR=6,NEE=NDOF*NEN)
	logical QUAD
	dimension XL(NSD,NEN),DMAT(NSTR,NSTR),SL(NEE,NEE),SH(4,NEN),
  >	 B(NSTR,NEE)
	dimension ROT(3,3)
*
	data E,POIS/1d0,0.25d0/

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

	anglex = 10d0
	anglez = 20d0
	write(*,20) anglex
20	format(1x,'Rotate XL over ',F4.0,' degrees along X-axis')
	call rotx(ROT,anglex)
	do i=1,NEN
	    call rotvec (XL(1,i),ROT)
	enddo
	write(*,29) anglez
29	format(1x,'Rotate XL over ',F4.0,' degrees along Z-axis')
	call rotz(ROT,anglez)
	do i=1,NEN
	    call rotvec (XL(1,i),ROT)
	enddo
	write(*,30) (j,(XL(i,j),i=1,NSD),j=1,NEN)
30	format(4(' XL(',I1,')=',3(F12.6)/))
*
	CALL MATERL(DMAT,E,POIS)
	write(*,40) ((DMAT(i,j),j=1,NSTR),i=1,NSTR)
40	format(1x,'Material matrix'/6(F12.6,1x)/)
*
	call SHPTET (0d0,0d0,0d0,XL,DET,SH,IERR,1)
	if (IERR.ne.0) stop 'determinant error'
c	Assemble B-matrix
	call BMATRIX (B,SH)
	write(*,45) ((B(i,j),j=1,NEE),i=1,NSTR)
45	format(1x,'B matrix'/12(F12.6,1x)/)

	call CLEAR (SL,NEE*NEE)
	CALL STIFF (DMAT,XL,SL,DUM1,QUAD,ierr)
	if (ierr.ne.0) stop 'determinant error'

	write(*,50) ((SL(i,j),j=1,NEE),i=1,NEE)
50	format(1x,'Stiffness matrix'/12(F12.6,1x)/)
	end

