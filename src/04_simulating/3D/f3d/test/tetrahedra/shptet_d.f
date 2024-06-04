	program main

	USE ALGEBRA, only: rotvec, rotz

	implicit none
	integer NEN,NSD
	parameter (NEN=4,NSD=3)
	double precision F,XL,DET,SH,TL,GRAD,TINT,XC,YC,ZC
	double precision angle,ROT
	integer IERR,i,j
	dimension XL(NSD,NEN),SH(4,NEN),TL(NEN),GRAD(NSD),ROT(3,3)
*
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
10	format(1x,'Input coordinates'/
  >	 4(' XL(',I1,')=',3(F12.6)/))
  	angle = 30.
	write(*,29) angle
29	format(1x,'Rotate XL over ',F4.0,' degrees along Z-axis')
!	call rotmat(ROT,angle)
	call rotz(ROT,angle)
	do i=1,NEN
	    call rotvec (XL(1,i),ROT)
	enddo
	write(*,30) (j,(XL(i,j),i=1,NSD),j=1,NEN)
30	format(4(' XL(',I1,')=',3(F12.6)/))
*
	call SHPTET (0d0,0d0,0d0,XL,DET,SH,IERR,1)
	print*,'DET=',DET
	if (ierr.ne.0) stop 'determinant error'
*
	TL(1) = 1d0
	TL(2) = 1d0
	TL(3) = 1d0
	TL(4) = 2d0
	GRAD(1)=TL(1)*SH(1,1)+TL(2)*SH(1,2)+TL(3)*SH(1,3)+TL(4)*SH(1,4)
	GRAD(2)=TL(1)*SH(2,1)+TL(2)*SH(2,2)+TL(3)*SH(2,3)+TL(4)*SH(2,4)
	GRAD(3)=TL(1)*SH(3,1)+TL(2)*SH(3,2)+TL(3)*SH(3,3)+TL(4)*SH(3,4)
	print*,'temperature gradient:'
	write(*,1) (GRAD(i),i=1,3)
 1	format(4(1X,1PG14.6))
*
	TINT = (TL(1)+TL(2)+TL(3)+TL(4))*DET*0.25D0
	print*,'Volume integral of T = ',TINT
*
	print*,'nodal shape function values:'
	do i=1,NEN
	    call SHPTET (XL(1,i),XL(2,i),XL(3,i),XL,DET,SH,IERR,2)
	    write(*,1) (SH(4,j),j=1,4)
	enddo
	XC = (XL(1,1)+XL(1,2)+XL(1,3)+XL(1,4))*0.25D0
	YC = (XL(2,1)+XL(2,2)+XL(2,3)+XL(2,4))*0.25D0
	ZC = (XL(3,1)+XL(3,2)+XL(3,3)+XL(3,4))*0.25D0
	print*,'element center and shape function values:'
	write(*,1) XC,YC,ZC
	call SHPTET (XC,YC,ZC,XL,DET,SH,IERR,2)
	write(*,1) (SH(4,j),j=1,4)
 	end
