	subroutine STIFF (DMAT,XL,SL,DUM1,QUAD,IERR)
*
c Computes the super-diagonal part of the local stiffness matrix.
*
	implicit none
	integer NEN,NSD,NDOF,NSTR,NEE
	parameter (NEN=4, NSD=3, NDOF=3, NSTR=6, NEE=NEN*NDOF)
c-pass
	integer IERR
	logical QUAD
	double precision DMAT,XL,SL,DUM1
	dimension XL(NSD,NEN),DMAT(NSTR,NSTR),SL(NEE,NEE)
c-locl
	double precision ZERO
	parameter (ZERO=0d0)
	integer j,i
	double precision SH,B,DET,DB1,DB2,DB3,DB4,DB5,DB6
	dimension SH(4,NEN),B(NSTR,NEE)
*
c	Compute spatial derivatives of shape functions and determinant
	call SHPTET (ZERO,ZERO,ZERO,XL,DET,SH,IERR,1)
	if (IERR.ne.0) return
c	Assemble B-matrix
	call BMATRIX (B,SH)
	do j=1,NEE
c	    Calculate column J of DB-matrix
	    DB1=DMAT(1,1)*B(1,j)+DMAT(1,2)*B(2,j)+DMAT(1,3)*B(3,j)
  >	     +DMAT(1,4)*B(4,j)+DMAT(1,5)*B(5,j)+DMAT(1,6)*B(6,j)
	    DB2=DMAT(2,1)*B(1,j)+DMAT(2,2)*B(2,j)+DMAT(2,3)*B(3,j)
  >	     +DMAT(2,4)*B(4,j)+DMAT(2,5)*B(5,j)+DMAT(2,6)*B(6,j)
	    DB3=DMAT(3,1)*B(1,j)+DMAT(3,2)*B(2,j)+DMAT(3,3)*B(3,j)
  >	     +DMAT(3,4)*B(4,j)+DMAT(3,5)*B(5,j)+DMAT(3,6)*B(6,j)
	    DB4=DMAT(4,1)*B(1,j)+DMAT(4,2)*B(2,j)+DMAT(4,3)*B(3,j)
  >	     +DMAT(4,4)*B(4,j)+DMAT(4,5)*B(5,j)+DMAT(4,6)*B(6,j)
	    DB5=DMAT(5,1)*B(1,j)+DMAT(5,2)*B(2,j)+DMAT(5,3)*B(3,j)
  >	     +DMAT(5,4)*B(4,j)+DMAT(5,5)*B(5,j)+DMAT(5,6)*B(6,j)
	    DB6=DMAT(6,1)*B(1,j)+DMAT(6,2)*B(2,j)+DMAT(6,3)*B(3,j)
  >	     +DMAT(6,4)*B(4,j)+DMAT(6,5)*B(5,j)+DMAT(6,6)*B(6,j)
	    do i=1,j
		SL(i,j)=SL(i,j)+DET*(B(1,i)*DB1+B(2,i)*DB2+B(3,i)*DB3
  >		 +B(4,i)*DB4+B(5,i)*DB5+B(6,i)*DB6)
  	    enddo
	enddo
*
	return
	end
