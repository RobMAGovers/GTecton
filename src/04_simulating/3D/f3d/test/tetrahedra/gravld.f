	subroutine GRAVLD (PL,XL,GRAV,DUM1,DENS,QUAD,IERR)
*
c Computes the contribution to the load from gravity body forces
*
	implicit none
	integer NEN,NSD,NDOF
	parameter (NEN=4, NSD=3, NDOF=3)
	double precision ZERO,QUART
	parameter (ZERO=0D0,QUART=0.25D0)
c-pass
	logical QUAD
	integer IERR
	double precision PL,XL,GRAV,DUM1,DENS
	dimension PL(NDOF,NEN),XL(NSD,NEN),GRAV(NDOF)
c-locl
	double precision DET,GRAV1,GRAV2,GRAV3,SH
	dimension SH(4,NEN)
*
c	Compute Jacobian determinant
	call SHPTET (ZERO,ZERO,ZERO,XL,DET,SH,IERR,1)
	if (IERR.ne.0) return
	GRAV1 = GRAV(1)*DENS*DET*QUART
	GRAV2 = GRAV(2)*DENS*DET*QUART
	GRAV3 = GRAV(3)*DENS*DET*QUART
	PL(1,1) = PL(1,1) + GRAV1
	PL(1,2) = PL(1,2) + GRAV1
	PL(1,3) = PL(1,3) + GRAV1
	PL(1,4) = PL(1,4) + GRAV1
	PL(2,1) = PL(2,1) + GRAV2
	PL(2,2) = PL(2,2) + GRAV2
	PL(2,3) = PL(2,3) + GRAV2
	PL(2,4) = PL(2,4) + GRAV2
	PL(3,1) = PL(3,1) + GRAV3
	PL(3,2) = PL(3,2) + GRAV3
	PL(3,3) = PL(3,3) + GRAV3
	PL(3,4) = PL(3,4) + GRAV3
*
	return
	end
