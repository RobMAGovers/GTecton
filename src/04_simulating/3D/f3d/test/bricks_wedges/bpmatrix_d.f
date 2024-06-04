	program main

	implicit double precision (a-h,o-z)
	parameter (NSD=3,NEN=8,NDOF=3,NSTR=6)
	parameter (NEE=NSTR*NEN)
	logical QUAD
	dimension XL(NSD,NEN),SH(4,NEN),XS(NSD,NSD),B(NSTR,NEE)
	dimension ROT(3,3)

	data R,S,T /0,0,0/
	data XL /-.5,-.5,-.5,
  >	          .5,-.5,-.5,
  >	          .5, .5,-.5,
  >	         -.5, .5,-.5,
  >	         -.5,-.5, .5,
  >	          .5,-.5, .5,
  >	          .5, .5, .5,
  >	         -.5, .5, .5/
	data QUAD /.true./

	write(*,10) (j,(XL(i,j),i=1,NSD),j=1,NEN)
10	format(1x,'Input data'/
  >	 8(' XL(',I1,')=',3(F12.6)/))
	write(*,20) QUAD,R,S,T
20	format(1x,'QUAD=',L1/1x,'(R,S,T)=',3(F12.6))

	call rotmat(ROT)
	do i=1,NEN
	    call rotvec (XL(1,i),ROT)
	enddo
	write(*,30) (j,(XL(i,j),i=1,NSD),j=1,NEN)
30	format(8(' XL(',I1,')=',3(F12.6)/))

	CALL SHAP30 (R,S,T,XL,DET,SH,XS,QUAD,IERR,3)
	if (ierr.ne.0) stop 'determinant error'
	write(*,40) (i,SH(1,i),i=1,NEN)
40	format(1x,
  >	 4('dSH/dx(',I1,')=',F12.6,1x)/1x,4('dSH/dx(',I1,')=',F12.6,1x))
	write(*,50) (i,SH(2,i),i=1,NEN)
50	format(1x,
  >	 4('dSH/dy(',I1,')=',F12.6,1x)/1x,4('dSH/dy(',I1,')=',F12.6,1x))
	write(*,60) (i,SH(3,i),i=1,NEN)
60	format(1x,
  >	 4('dSH/dz(',I1,')=',F12.6,1x)/1x,4('dSH/dz(',I1,')=',F12.6,1x))

	CALL BMATRIX (B,SH)

	do k=1,NEN
	    write(*,70) k,((B(i,j),j=k*NDOF-2,k*NDOF),i=1,NSTR)
70	    format(1x,'B',I1/6(3(1x,F12.6)/))
	enddo

	end
c-------------------------------------------------------------------------------
	subroutine rotmat(R)

	USE CONSTANTS, only deg2rad

	implicit double precision (a-h,o-z)
	dimension R(3,3)
	angle = 30

	write(*,10) angle
10	format(1x,'Rotate XL over ',F4.0,' degrees along Z-axis')

	angle = angle * deg2rad
	SINA = SIN(angle)
	COSA = COS(angle)

	R(1,1) = COSA
	R(2,2) = COSA
	R(1,2) = -SINA
	R(2,1) = SINA
	R(3,3) = 1.
	R(1,3) = 0.
	R(2,3) = 0.
	R(3,1) = 0.
	R(3,2) = 0.

	return
	end
c-------------------------------------------------------------------------------
	subroutine rotvec (a,R)

	implicit double precision (a-h,o-z)
	dimension a(3),R(3,3)
	dimension c(3)

	do i=1,3
	    c(i) = 0.
	    do j=1,3
		c(i) = c(i) + R(i,j)*a(j)
	    enddo
	enddo
	do i=1,3
	    a(i) = c(i)
	enddo
	return
	end
