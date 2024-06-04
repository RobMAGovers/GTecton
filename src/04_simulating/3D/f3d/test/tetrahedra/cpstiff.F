 subroutine CPSTIFF (XL,QUAD,S,RHOCP,DUM1,IERR)
!
! computes the upper local heat capacity stiffness matrix
!
 implicit none
 integer NEN,NSD
  parameter (NEN=4,NSD=3)
 double precision ZERO,TWO
 parameter (ZERO=0d0,TWO=2D0)
!-pass
 integer IERR
 logical QUAD
 double precision XL,S,RHOCP,DUM1
 dimension XL(NSD,NEN),S(NEN,NEN)
!-locl
 integer j,i
 double precision det,sh,w
 dimension sh(4,NEN)
!
 call SHPTET (ZERO,ZERO,ZERO,XL,det,sh,IERR,1)
 if (IERR.ne.0) return
 w = det * 0.05D0 * RHOCP
 S(1,1) = S(1,1) + w*TWO
 S(1,2) = S(1,2) + w
 S(1,3) = S(1,3) + w
 S(1,4) = S(1,4) + w
 S(2,2) = S(2,2) + w*TWO
 S(2,3) = S(2,3) + w
 S(2,4) = S(2,4) + w
 S(3,3) = S(3,3) + w*TWO
 S(3,4) = S(3,4) + w
 S(4,4) = S(4,4) + w*TWO
!
 return
 end
