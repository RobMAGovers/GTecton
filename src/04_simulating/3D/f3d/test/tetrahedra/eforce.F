 subroutine EFORCE (EVP,PL,XL,DUM1,QUAD,IERR)
!
! Computes stress-equivalent nodal forces PL
!
 implicit none
 integer NEN,NSD,NDOF,NSTR,NEE
 parameter (NEN=4,NSD=3,NDOF=3,NSTR=6,NEE=NEN*NDOF)
!-pass
 integer IERR
 logical QUAD
 double precision EVP,PL,XL,DUM1
 dimension XL(NSD,NEN),PL(NEE),EVP(NSTR)

!-locl
 integer i
 double precision SH,B,ZERO,DET
 parameter (ZERO=0d0)
 dimension SH(4,NEN),B(NSTR,NEE)
!
!    Compute Jacobian determinant and shape function derivatives
 call SHPTET (ZERO,ZERO,ZERO,XL,DET,SH,IERR,1)
 if (IERR.ne.0) return
 call BMATRIX (B,SH)
 do i=1,NEE
PL(i) = PL(i) + DET*(B(1,i)*EVP(1)+B(2,i)*EVP(2)+ &
       B(3,i)*EVP(3)+B(4,i)*EVP(4)+B(5,i)*EVP(5)+B(6,i)*EVP(6))
 enddo
!
 return
 end
