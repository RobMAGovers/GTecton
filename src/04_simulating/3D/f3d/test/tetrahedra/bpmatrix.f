c-------------------------------------------------------------------------------
	subroutine BPMATRIX (XL,DL,QUAD,BP,IERR)
*
c Sets up derivatives matrix BP
*
	implicit none
	integer NDOF,NSD,NEN
	parameter (NDOF=3, NSD=3, NEN=4)
	double precision ZERO
	parameter (ZERO=0d0)
c-pass
	logical QUAD
	integer IERR
	double precision XL,DL,BP
	dimension XL(NSD,NEN),DL(NDOF,NEN),BP(NDOF,NSD)
c-locl
	double precision det,sh
	dimension sh(4,NEN)
*
c	Compute spatial derivatives of shape functions
	call SHPTET (ZERO,ZERO,ZERO,XL,det,sh,IERR,1)
	if (IERR.ne.0) return
*
	BP(1,1) = SH(1,1)*DL(1,1)+SH(1,2)*DL(1,2)+SH(1,3)*DL(1,3)+
  >	 SH(1,4)*DL(1,4)
	BP(2,1) = SH(1,1)*DL(2,1)+SH(1,2)*DL(2,2)+SH(1,3)*DL(2,3)+
  >	 SH(1,4)*DL(2,4)
	BP(3,1) = SH(1,1)*DL(3,1)+SH(1,2)*DL(3,2)+SH(1,3)*DL(3,3)+
  >	 SH(1,4)*DL(3,4)
	BP(1,2) = SH(2,1)*DL(1,1)+SH(2,2)*DL(1,2)+SH(2,3)*DL(1,3)+
  >	 SH(2,4)*DL(1,4)
	BP(2,2) = SH(2,1)*DL(2,1)+SH(2,2)*DL(2,2)+SH(2,3)*DL(2,3)+
  >	 SH(2,4)*DL(2,4)
	BP(3,2) = SH(2,1)*DL(3,1)+SH(2,2)*DL(3,2)+SH(2,3)*DL(3,3)+
  >	 SH(2,4)*DL(3,4)
	BP(1,3) = SH(3,1)*DL(1,1)+SH(3,2)*DL(1,2)+SH(3,3)*DL(1,3)+
  >	 SH(3,4)*DL(1,4)
	BP(2,3) = SH(3,1)*DL(2,1)+SH(3,2)*DL(2,2)+SH(3,3)*DL(2,3)+
  >	 SH(3,4)*DL(2,4)
	BP(3,3) = SH(3,1)*DL(3,1)+SH(3,2)*DL(3,2)+SH(3,3)*DL(3,3)+
  >	 SH(3,4)*DL(3,4)
*
	return
	end
