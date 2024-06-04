	subroutine TSTIFF (XL,QUAD,SL,C,DUM1,IERR)
*
c computes the local conduction stiffness matrix from integration at Barlow
c points. Solely calculates upper-diagonal matrix.
*
	implicit none
	integer NEN,NSD,NDOF,NSTR,NEE
	parameter (NEN=4, NSD=3, NDOF=3, NSTR=6, NEE=NEN*NDOF)
c-pass
	integer IERR
	logical QUAD
	double precision C,XL,SL,DUM1
	dimension XL(NSD,NEN),C(NSD),SL(NEE,NEE)
c-locl
	double precision ZERO
	parameter (ZERO=0d0)
	integer j,i
	double precision SH,DET
	dimension SH(4,NEN)
*
c	Compute spatial derivatives of shape functions and determinant
	call SHPTET (ZERO,ZERO,ZERO,XL,DET,SH,IERR,1)
	if (IERR.ne.0) return
	do j=1,NEN
	    do i=1,j
		SL(i,j) = SL(i,j) + DET*(C(1)*SH(1,i)*SH(1,j)
  >		 +C(2)*SH(2,i)*SH(2,j)+C(3)*SH(3,i)*SH(3,j))
  	    enddo
	enddo
*
	return
	end
