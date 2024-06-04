! not tested yet
 subroutine ADVSTF (S,XL,VL,PL,COND,RHOCP,TH,HEAT,QUAD,IERR)
!
! convection (advection) update of local conduction stiffness matrix for
! steady convection-conduction solution
!
 implicit none
 integer NUMOPT
 parameter (NUMOPT=1)
!
! NUMOPT = 0: no upwinding
! NUMOPT = 1: Mizukami upwind scheme
!
 integer NEN,NDOF,NSD
  parameter (NEN=4, NDOF=3, NSD=3)
 double precision ZERO
 parameter (ZERO=0d0)
!-pass
 logical QUAD
 integer IERR
 double precision S,XL,VL,PL,COND,RHOCP,TH,HEAT
 dimension XL(NSD,NEN),VL(NDOF,NEN),S(NEN,NEN),COND(NSD),PL(NEN)

!-locl
 integer j,i
 double precision sh,sm,sa,sd,det,w,vx,vy,vz,samax,cplump, &
   tau,qi
 dimension sh(4,NEN),sm(NEN,NEN),sa(NEN,NEN),sd(NEN,NEN)
 logical firstentry
save firstentry
!-init
 data firstentry /.true./
!
 if (firstentry) then
if (iecho.ne.0) then
    if (NUMOPT.eq.0) then
    write(stderr,5)
 5            format(1x,'No upwinding')
else
    write(stderr,10)
10            format(1x,'Upwind: Mizukami scheme')
endif
endif
firstentry = .false.
 endif
!
 call CLEAR(sa,16,"sa")
 if (NUMOPT.ne.0) then
!        consistent mass matrix SM
call CLEAR(sm,16,"sm")
call CPSTIFF (XL,QUAD,sm,RHOCP,TH,IERR)
if (IERR.ne.0) return
!        conduction stiffness matrix SD
call CLEAR(sd,16,"sd")
call TSTIFF (XL,QUAD,sd,COND,TH,IERR)
 endif
!
 call SHPTET (ZERO,ZERO,ZERO,XL,det,sh,IERR,1)
 if (IERR.ne.0) return
 w = det * 5.0D-2 * RHOCP
 vx = VL(1,1)+VL(1,2)+VL(1,3)+VL(1,4)
 vy = VL(2,1)+VL(2,2)+VL(2,3)+VL(2,4)
 vz = VL(3,1)+VL(3,2)+VL(3,3)+VL(3,4)
 do j=1,NEN
do i=1,NEN
sa(i,j) = w*( (vx+VL(1,i))*sh(1,j) + &
       (vy+VL(2,i))*sh(2,j) + (vz+VL(3,i))*sh(3,j) )
enddo
 enddo
 do i=1,NEN
do j=1,NEN
S(i,j)  = S(i,j) + sa(i,j)
sa(i,j) = ABS(sa(i,j))
enddo
 enddo
 if (NUMOPT.ne.0) then
call SYMMET(sm,NEN)
do i=1,NEN
samax = MAX(sa(i,1),sa(i,2),sa(i,3),sa(i,4))
!        element ii of lumped mass matrix
cplump = sm(i,1)+sm(i,2)+sm(i,3)+sm(i,4)
tau = 5.0D-1*(cplump/(samax+sd(i,i)))
qi = tau*(vx*sh(1,i)+vy*sh(2,i)+vz*sh(3,i))
PL(i) = PL(i) + HEAT*det*qi
w = det * 0.25D0 * RHOCP *qi
    do j=1,NEN
    S(i,j)=S(i,j)+w*(vx*sh(1,j)+vy*sh(2,j)+vz*sh(3,j))
enddo
enddo
 endif
!
 return
 end
