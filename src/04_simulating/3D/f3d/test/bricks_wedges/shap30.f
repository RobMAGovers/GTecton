C-------------------------------------------------------------------------------
* RG 30/1/'94: not checked
*
	SUBROUTINE SHAP30 (R,S,T,X,DET,SH,XS,QUAD,IERR,IFLAG)
*
c Program to compute shape functions for quadrilateral or triangle
*
c				= 1, calculate r-, s-, t- derivatives
c				  and shape functions only
c				= 2, calculate r-, s-, t- derivatives
c				  and shape functions and Jacobian
c IFLAG				= 3, full calculation
c R,S,T,			= natural coordinates
c X(NSD,NEN)			= global coordinates
c QUAD				= TRUE if brick element
c				= FALSE if prism element
c SH(1,NEN)			= x- or r- derivatives of shape functions
c SH(2,NEN)			= y- or s- derivatives of shape functions
c SH(3,NEN)			= z- or t- derivatives of shape functions
c SH(4,NEN)			= shape functions
c XS(NSD,NSD)			= inverse jacobian matrix
c DET				= Jacobian matrix determinant
c IERR				= 0, succesful calculation
c				= 1, Jacobian determinant error
*
c SNGL	implicit real (a-h,o-z)
	implicit double precision (A-H,O-Z)
	parameter (NEN=8, NSD=3)
	parameter (ZERO=0., ONE=1., F=.125)
c-pass
	logical QUAD
	integer ierr
	dimension SH(4,NEN),X(NSD,NEN),XS(NSD,NSD)
c-locl
	dimension RA(NEN),SA(NEN),TA(NEN)
	save RA,SA,TA
c-comm
c-init
	data RA/-1., 1., 1.,-1.,-1., 1., 1.,-1./
	data SA/-1.,-1., 1., 1.,-1.,-1., 1., 1./
	data TA/-1.,-1.,-1.,-1., 1., 1., 1., 1./
	ierr = 0
*
c	Calculate shape functions and r-, s- and t-derivatives of shape functions
	do 10 i=1,NEN
	    RX = (ONE+RA(i)*R)
	    SX = (ONE+SA(i)*S)
	    TX = (ONE+TA(i)*T)
	    SH(4,i) = F*RX*SX*TX
	    SH(1,i) = F*RA(i)*SX*TX
	    SH(2,i) = F*SA(i)*RX*TX
	    SH(3,i) = F*TA(i)*RX*SX
10	continue
	if (.not.QUAD) then
c	    Collapse of brick element on prism.
	    do 20 i=1,4
		SH(i,3) = SH(i,3)+SH(i,4)
		SH(i,7) = SH(i,7)+SH(i,8)
		SH(i,4) = ZERO
		SH(i,8) = ZERO
20	    continue
	endif
c	Jacobian matrix calculation
	do 40 i=1,3
	    do 30 j=1,3
		XS(i,j) = SH(i,1)*X(j,1)+SH(i,2)*X(j,2)+SH(i,3)*X(j,3)
  >		 +SH(i,4)*X(j,4)+SH(i,5)*X(j,5)+SH(i,6)*X(j,6)
  >		 +SH(i,7)*X(j,7)+SH(i,8)*X(j,8)
30	    continue
40	continue
	DET = XS(1,1)*XS(2,2)*XS(3,3)-XS(1,1)*XS(2,3)*XS(3,2)
  >	     -XS(2,1)*XS(1,2)*XS(3,3)+XS(2,1)*XS(1,3)*XS(3,2)
  >	     +XS(3,1)*XS(1,2)*XS(2,3)-XS(3,1)*XS(1,3)*XS(2,2)
	if (DET.le.zero) go to 1000
*
	if (IFLAG.eq.1) return
*
c	Calculatate inverse Jacobian determinant and store it back in XS
	RDET=ONE/DET
	Xi11=(XS(2,2)*XS(3,3)-XS(3,2)*XS(2,3))*RDET
	Xi21=(XS(2,3)*XS(3,1)-XS(2,1)*XS(3,3))*RDET
	Xi31=(XS(2,1)*XS(3,2)-XS(2,2)*XS(3,1))*RDET
	Xi12=(XS(1,3)*XS(3,2)-XS(1,2)*XS(3,3))*RDET
	Xi22=(XS(1,1)*XS(3,3)-XS(1,3)*XS(3,1))*RDET
	Xi32=(XS(1,2)*XS(3,1)-XS(1,1)*XS(3,2))*RDET
	Xi13=(XS(1,2)*XS(2,3)-XS(1,3)*XS(2,2))*RDET
	Xi23=(XS(1,3)*XS(2,1)-XS(1,1)*XS(2,3))*RDET
	Xi33=(XS(1,1)*XS(2,2)-XS(1,2)*XS(2,1))*RDET
	XS(1,1)=Xi11
	XS(2,1)=Xi21
	XS(3,1)=Xi31
	XS(1,2)=Xi12
	XS(2,2)=Xi22
	XS(3,2)=Xi32
	XS(1,3)=Xi13
	XS(2,3)=Xi23
	XS(3,3)=Xi33
*
	if (IFLAG.eq.2) return
*
c	Calculate x-, y- and z-derivatives of shape functions
	do 50 i=1,NEN
	    dSHdR = SH(1,i)
	    dSHdS = SH(2,i)
	    dSHdT = SH(3,i)
	    SH(1,i) = XS(1,1)*dSHdR+XS(1,2)*dSHdS+XS(1,3)*dSHdT
	    SH(2,i) = XS(2,1)*dSHdR+XS(2,2)*dSHdS+XS(2,3)*dSHdT
	    SH(3,i) = XS(3,1)*dSHdR+XS(3,2)*dSHdS+XS(3,3)*dSHdT
50	continue
*
	return
*
 1000	write(stderr,2000) DET
 2000	format(///1x,'Shape function fails! Determinant is ',1PE20.4)
	ierr = 1
	return
	end
