	program main

	USE ALGEBRA, only: rotz, rotx, rotvec, clear

	implicit double precision (a-h,o-z)
	parameter (NSD=3,NEN=4,NDOF=3,NSTR=6,NEE=NEN*NDOF)
	logical QUAD
	dimension XL(NSD,NEN),PL(NDOF,NEN),GRAV(NDOF)
	dimension ROT(3,3),STN(6),SH(4,NEN),B(NSTR,NEE)
*
	F = 0.5d0*6d0**0.3333333333333333333d0 ! makes volume 1
	F = 1.0d0
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
	STN(1) = 0.0
	STN(2) = 0.0
	STN(3) = 1.0
	STN(4) = 0.0
	STN(5) = 0.0
	STN(6) = 0.0
*
	write(*,10) (j,(XL(i,j),i=1,NSD),j=1,NEN)
10	format(1x,'Input data'/
  >	 4(' XL(',I1,')=',3(F12.6)/))

	anglex = 90d0
	anglez = 0d0
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

c       Assemble B-matrix
	call SHPTET (0d0,0d0,0d0,XL,DET,SH,IERR,1)
	print*,'VOLUME = ',DET
call BMATRIX (B,SH)
write(*,45) ((B(i,j),j=1,NEE),i=1,NSTR)
45   format(1x,'B matrix'/12(F12.6,1x)/)
*
	call CLEAR (PL,NDOF*NEN)
	call EFORCE (STN,PL,XL,DUM1,QUAD,IERR)
	if (ierr.ne.0) stop 'determinant error'

	write(*,50) (j,(PL(i,j),i=1,NDOF),j=1,NEN)
50	format(1x,'Output data'/
  >   4(' PL(',I1,')=',3(F12.6,1x)/))
	end


