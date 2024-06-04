	program main

	USE ALGEBRA, only: rotz, rotx, rotvec, clear

	implicit double precision (a-h,o-z)
	parameter (NSD=3,NEN=4)
	dimension XL(NSD,NEN),LIEN(3)
	dimension ROT(3,3),XN(3)
*
F = 0.5d0
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

	anglex = 0d0
	anglez = 45d0
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
100	write(*,40,advance='no')
40	format(/1x,'side number > ')
	read(*,*) is
	call SIDENP (is,.true.,lien)
	write(*,50) is,lien(1),lien(2),lien(3)
50	format(1x,'side ',I1,'= nodes',3(1X,I1))
	call FCGEOM (XL,lien,AREA,XN)
	write(*,60) AREA,XN(1),XN(2),XN(3)
60	FORMAT(1x,'AREA = ',1F12.6,', normal=',3(1X,F12.6))
	goto 100
	end

c--
	subroutine xit(i)
	stop
	end
