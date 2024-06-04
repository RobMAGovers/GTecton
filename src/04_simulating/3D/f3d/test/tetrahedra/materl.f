	subroutine MATERL (DMAT,E,POIS)
*
c Constructs the material matrix relating stress and strain.
*
	implicit none
	integer NSTR
	parameter (NSTR=6)
	double precision ONE,TWO,HALF
	parameter (ONE=1d0,TWO=2d0,HALF=.5d0)
c-pass
	double precision DMAT,E,POIS
	dimension DMAT(NSTR,NSTR)
c-locl
	double precision AM,AL
*
	call CLEAR(DMAT,36)
c	Compute Lame parameters AM(=MU) and AL(=LAMBDA)
	AM=E/(ONE+POIS)
	AL=AM*POIS/(ONE-TWO*POIS)
	AM=HALF*AM
	DMAT(1,1)=TWO*AM+AL
	DMAT(1,2)=AL
	DMAT(1,3)=AL
	DMAT(2,1)=AL
	DMAT(2,2)=DMAT(1,1)
	DMAT(2,3)=AL
	DMAT(3,1)=AL
	DMAT(3,2)=AL
	DMAT(3,3)=DMAT(2,2)
	DMAT(4,4)=AM
	DMAT(5,5)=AM
	DMAT(6,6)=AM
*
	return
	end
