	program main

	USE MODELDEFINITION

	implicit double precision (a-h,o-z)
	parameter (NSD=3,NEN=8,NDOF=3,NSTR=6)
	logical QUAD
	dimension XL(NSD,NEN),PL(NDOF,NEN),STN(NSTR)
	dimension ROT(3,3)

	data XL /-.5,-.5,-.5,
  >	          .5,-.5,-.5,
  >	          .5, .5,-.5,
  >	         -.5, .5,-.5,
  >	         -.5,-.5, .5,
  >	          .5,-.5, .5,
  >	          .5, .5, .5,
  >	         -.5, .5, .5/
	data QUAD /.false./
	data STN / 0.,0.,1.,0.,0.,0./
	angle =  0
	IPOINT = 1

	if (.not.QUAD) then
	    do i=1,NSD
		XL(i,4) = XL(i,3)
		XL(i,8) = XL(i,7)
	    enddo
	endif

	write(*,"(80('-'))")
	write(*,10) (j,(XL(i,j),i=1,NSD),j=1,NEN)
10	format(1x,'Input data'/
  >	 8(' XL(',I1,')=',3(F12.6)/))
	write(*,20) QUAD,IPOINT
20	format(1x,'QUAD=',L1,'  IPOINT=',I1)
	write(*,26) (STN(i),i=1,NSTR)
26	format(1x,'STN= ',6(F12.6,1x))

	write(*,29) angle
29	format(1x,'Rotate XL over ',F4.0,' degrees along Z-axis')
	call rotmat(ROT,angle)
	do i=1,NEN
	    call rotvec (XL(1,i),ROT)
	enddo

	CALL EFORCE (STN,PL,XL,DUM1,QUAD,ierr)
	if (ierr.ne.0) stop 'determinant error'

	write(*,50) (j,(PL(i,j),i=1,NDOF),j=1,NEN)
50	format(1x,'Output data'/
  >   8(' PL(',I1,')=',3(F12.6,1x)/))
	end
c-------------------------------------------------------------------------------
	subroutine clear (A,N)
	implicit double precision (a-h,o-z)
	dimension A(*)
	do 100 i=1,N
100	    A(i) =0.
	return
	end
c-------------------------------------------------------------------------------
	subroutine rotmat(R,rangle)

	USE CONSTANTS, only: deg2rad

	implicit double precision (a-h,o-z)
	dimension R(3,3)

	angle = rangle * deg2rad
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
