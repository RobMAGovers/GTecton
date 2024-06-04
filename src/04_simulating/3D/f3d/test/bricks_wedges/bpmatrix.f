c-------------------------------------------------------------------------------
* RG 30/1/'94: not checked
*
	subroutine BPMATRIX (XL,DL,QUAD,BP,ierr)
*
c Sets up derivatives matrix BP
*
c SNGL	implicit real (a-h,o-z)
	implicit double precision (a-h,o-z)
	parameter (NDOF=3, NSD=3, NEN=8)
	parameter (ZERO=0., THIRDM=-1d0/3d0)
c-pass
	logical QUAD
	dimension XL(NSD,NEN),DL(NDOF,NEN),BP(NDOF,NSD)
c-locl
	dimension SH(4,NEN),XS(NSD,NSD)
*
c	Compute spatial derivatives of shape functions
	sn = ZERO
if (.not.QUAD) sn = THIRDM
call SHAP30(ZERO,sn,ZERO,XL,DET,SH,XS,QUAD,ierr,3)
if (ierr.ne.0) return
*
	call CLEAR(BP,NDOF*NSD)
	do 200 k=1,NEN
	    do 100 j=1,NSD
		BP(1,j) = BP(1,j)+SH(j,k)*DL(1,k)
		BP(2,j) = BP(2,j)+SH(j,k)*DL(2,k)
		BP(3,j) = BP(3,j)+SH(j,k)*DL(3,k)
100	    continue
200	continue
*
	return
	end
