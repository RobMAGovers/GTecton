	subroutine BMATRIX (B,SH)
*
C Defines the strain-displacement matrix B
*
	implicit none
	integer NEN,NSTR,NDOF,NEE
	parameter (NEN=4, NSTR=6, NDOF=3, NEE=NEN*NDOF)
	double precision ZERO
	parameter (ZERO=0d0)
c-pass
	double precision B,SH
	dimension B(NSTR,NEE),SH(4,NEN)
c-locl
	integer k,i
*
c	Column number k
	k = 0
	do i=1,NEN
	    k = k + 1
	    B(1,k) = SH(1,i)
	    B(2,k) = ZERO
	    B(3,k) = ZERO
	    B(4,k) = SH(2,i)
	    B(5,k) = SH(3,i)
	    B(6,k) = ZERO
	    k = k + 1
	    B(1,k) = ZERO
	    B(2,k) = SH(2,i)
	    B(3,k) = ZERO
	    B(4,k) = SH(1,i)
	    B(5,k) = ZERO
	    B(6,k) = SH(3,i)
	    k = k + 1
	    B(1,k) = ZERO
	    B(2,k) = ZERO
	    B(3,k) = SH(3,i)
	    B(4,k) = ZERO
	    B(5,k) = SH(1,i)
	    B(6,k) = SH(2,i)
	enddo
*
	return
	end
