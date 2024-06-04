	subroutine SHPTET (X,Y,Z,XL,DET,SH,IERR,IFLAG)
	USE CONSTANTS, only: sixth, one
*
c Computes shape functions and their derivatives for linear tetrahedra (simplex
c element).
*
	implicit none
	integer NEN,NSD
	parameter (NEN=4,NSD=3)
c-pass
	integer IERR,IFLAG
	double precision X,Y,Z,XL,SH,DET
	dimension XL(NSD,NEN),SH(4,NEN)
c-comm
c-locl
	double precision A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3,
  >	 D4,RDET
c-init
	IERR = 0
*
	A1 = XL(1,2)*(XL(2,3)*XL(3,4)-XL(3,3)*XL(2,4))-
  >	     XL(1,3)*(XL(2,2)*XL(3,4)-XL(3,2)*XL(2,4))+
  >	     XL(1,4)*(XL(2,2)*XL(3,3)-XL(3,2)*XL(2,3))
	A2 = XL(1,1)*(XL(3,3)*XL(2,4)-XL(2,3)*XL(3,4))-
  >	     XL(1,3)*(XL(3,1)*XL(2,4)-XL(2,1)*XL(3,4))+
  >	     XL(1,4)*(XL(3,1)*XL(2,3)-XL(2,1)*XL(3,3))
	A3 = XL(1,1)*(XL(2,2)*XL(3,4)-XL(3,2)*XL(2,4))-
  >	     XL(1,2)*(XL(2,1)*XL(3,4)-XL(3,1)*XL(2,4))+
  >	     XL(1,4)*(XL(2,1)*XL(3,2)-XL(3,1)*XL(2,2))
	A4 = XL(1,1)*(XL(3,2)*XL(2,3)-XL(2,2)*XL(3,3))-
  >	     XL(1,2)*(XL(3,1)*XL(2,3)-XL(2,1)*XL(3,3))+
  >	     XL(1,3)*(XL(3,1)*XL(2,2)-XL(2,1)*XL(3,2))
	B1 = XL(3,3)*XL(2,4)-XL(2,3)*XL(3,4) -
  >	     XL(3,2)*XL(2,4)+XL(2,2)*XL(3,4) +
  >	     XL(3,2)*XL(2,3)-XL(2,2)*XL(3,3)
	B2 = XL(2,3)*XL(3,4)-XL(3,3)*XL(2,4) -
  >	     XL(2,1)*XL(3,4)+XL(3,1)*XL(2,4) +
  >	     XL(2,1)*XL(3,3)-XL(3,1)*XL(2,3)
	B3 = XL(3,2)*XL(2,4)-XL(2,2)*XL(3,4) -
  >	     XL(3,1)*XL(2,4)+XL(2,1)*XL(3,4) +
  >	     XL(3,1)*XL(2,2)-XL(2,1)*XL(3,2)
	B4 = XL(2,2)*XL(3,3)-XL(3,2)*XL(2,3) -
  >	     XL(2,1)*XL(3,3)+XL(3,1)*XL(2,3) +
  >	     XL(2,1)*XL(3,2)-XL(3,1)*XL(2,2)
	C1 = XL(1,3)*XL(3,4)-XL(3,3)*XL(1,4) -
  >	     XL(1,2)*XL(3,4)+XL(3,2)*XL(1,4) +
  >	     XL(1,2)*XL(3,3)-XL(3,2)*XL(1,3)
	C2 = XL(3,3)*XL(1,4)-XL(1,3)*XL(3,4) -
  >	     XL(3,1)*XL(1,4)+XL(1,1)*XL(3,4) +
  >	     XL(3,1)*XL(1,3)-XL(1,1)*XL(3,3)
	C3 = XL(1,2)*XL(3,4)-XL(3,2)*XL(1,4) -
  >	     XL(1,1)*XL(3,4)+XL(3,1)*XL(1,4) +
  >	     XL(1,1)*XL(3,2)-XL(3,1)*XL(1,2)
	C4 = XL(3,2)*XL(1,3)-XL(1,2)*XL(3,3) -
  >	     XL(3,1)*XL(1,3)+XL(1,1)*XL(3,3) +
  >	     XL(3,1)*XL(1,2)-XL(1,1)*XL(3,2)
	D1 = XL(2,3)*XL(1,4)-XL(1,3)*XL(2,4) -
  >	     XL(2,2)*XL(1,4)+XL(1,2)*XL(2,4) +
  >	     XL(2,2)*XL(1,3)-XL(1,2)*XL(2,3)
	D2 = XL(1,3)*XL(2,4)-XL(2,3)*XL(1,4) -
  >	     XL(1,1)*XL(2,4)+XL(2,1)*XL(1,4) +
  >	     XL(1,1)*XL(2,3)-XL(2,1)*XL(1,3)
	D3 = XL(2,2)*XL(1,4)-XL(1,2)*XL(2,4) -
  >	     XL(2,1)*XL(1,4)+XL(1,1)*XL(2,4) +
  >	     XL(2,1)*XL(1,2)-XL(1,1)*XL(2,2)
	D4 = XL(1,1)*XL(2,2)-XL(2,1)*XL(1,2) -
  >	     XL(1,1)*XL(2,3)+XL(2,1)*XL(1,3) +
  >	     XL(1,2)*XL(2,3)-XL(2,2)*XL(1,3)
*
	DET = A1+A2+A3+A4
	if (DET.le.0d0) goto 1000
	RDET = ONE/DET
*
	if (IFLAG.eq.1) then
c	    Compute x-, y-, and z-derivatives of shape functions
	    SH(1,1) = B1*RDET
	    SH(1,2) = B2*RDET
	    SH(1,3) = B3*RDET
	    SH(1,4) = B4*RDET
	    SH(2,1) = C1*RDET
	    SH(2,2) = C2*RDET
	    SH(2,3) = C3*RDET
	    SH(2,4) = C4*RDET
	    SH(3,1) = D1*RDET
	    SH(3,2) = D2*RDET
	    SH(3,3) = D3*RDET
	    SH(3,4) = D4*RDET
	else
	    SH(4,1) = (A1+B1*X+C1*Y+D1*Z)*RDET
	    SH(4,2) = (A2+B2*X+C2*Y+D2*Z)*RDET
	    SH(4,3) = (A3+B3*X+C3*Y+D3*Z)*RDET
	    SH(4,4) = (A4+B4*X+C4*Y+D4*Z)*RDET
	endif
*
c	compute volume
	DET = DET*SIXTH
*
	return
*
1000	write(0,1) DET
 1	format(///1x,'Shape function fails! Determinant is ',1PE20.4)
	IERR = 1
	return
	end
