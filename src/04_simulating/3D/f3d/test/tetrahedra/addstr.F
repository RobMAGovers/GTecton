 subroutine ADDSTR (B,X,D,IEN,LM,DUM1,DUM2,IELSTR,ISSIDE,STRS, &
   ISTR,LMF,TFAULT,SKEW,DXE)

 USE MODELDEFINITION
 USE CONSTANTS, only: third,eps
use timestepmodule,  only: NSTEP,ITWORK,IMWORK,IMTWRK,REFACT

!
! Subroutine to add stresses to load vector
!
 implicit none
 integer NSD,NDOF,NSTR,NEN
 parameter (NSD=3,NDOF=3,NSTR=6,NEN=4)
!-pass
 integer IEN,LM,IELSTR,ISSIDE,ISTR,LMF
 double precision B,X,D,DUM1,DUM2,STRS,TFAULT,SKEW,DXE
 dimension DUM1(*),IEN(NEN,*),X(NSD,*),D(NDOF,*),DUM2(*),B(*), &
   LM(NDOF,NEN,*),STRS(NSTR,*),IELSTR(*),ISSIDE(*),ISTR(2,*), &
   LMF(NDOF,NEN,*),TFAULT(NDOF,*),SKEW(2,*), &
   DXE(NDOF,NEN,*)
!-locl
 integer numadd,n,k,it0,it1,lien,i,np,l,ii
 double precision xl,dl,area,xn,angle,rot,tmp,sxx,syy,szz,sxy, &
   sxz,syz,f
 dimension xl(NSD,NEN),dl(NDOF,NEN),lien(3),xn(NSD),angle(2), &
   rot(3,3),tmp(NSD),f(NDOF)
!-init
 numadd = 0
!
 do 20 n=1,NUMSTR
k = IELSTR(n)
if (k.le.0) then
write(stderr,1) n
 1        format(1x,'ADDSTR: fatal IELSTR error. Index=',I5)
call stoper()
endif
it0 = ISTR(1,n)
it1 = ISTR(2,n)
if (DBLE(NSTEP-it0)*DBLE(NSTEP-it1).gt.0d0) goto 20
numadd = numadd + 1
call LCOORD (X,xl,IEN(1,k))
if (LGDEF.ne.0) then
    call LDISP  (dl,D,IEN(1,k),NDOF,NEN)
call AddFaultDisplacement (dl,LMF(1,1,k),TFAULT,NDOF,NEN)
    call ADDSNE (dl,DXE(1,1,k),NDOF,NEN)
call REZONE (xl,dl,1d0)
endif
!        Get nodal point numbers of the side on which pressure works
call SIDENP(ISSIDE(n),.true.,lien)
!        Calculate surface area and normal vector
call FCGEOM (xl,lien,area,xn)
!        Rotate face into local direction
angle(1) = 0d0
angle(2) = 0d0
do i=1,3
np = IEN(lien(i),k)
angle(1) = angle(1) - SKEW(1,np)
angle(2) = angle(2) - SKEW(2,np)
enddo
angle(1) = angle(1)*THIRD
angle(2) = angle(2)*THIRD
if (ABS(angle(1)).gt.EPS.or.ABS(angle(2)).gt.EPS) then
call FORMRT(angle,rot,NDOF)
call EQUATE (tmp,xn,NDOF)
call VMPRD (rot,tmp,xn,NDOF,3)
endif
 sxx = STRS(1,n)
 syy = STRS(2,n)
 szz = STRS(3,n)
 sxy = STRS(4,n)
 sxz = STRS(5,n)
 syz = STRS(6,n)
    f(1)=(sxx*xn(1)+sxy*xn(2)+sxz*xn(3))*area*THIRD
    f(2)=(sxy*xn(1)+syy*xn(2)+syz*xn(3))*area*THIRD
    f(3)=(sxz*xn(1)+syz*xn(2)+szz*xn(3))*area*THIRD
do l=1,NSD
ii=LM(l,lien(1),k)
      if(ii.ne.0) B(ii)=B(ii)+f(l)
      ii=LM(l,lien(2),k)
      if(ii.ne.0) B(ii)=B(ii)+f(l)
      ii=LM(l,lien(3),k)
      if(ii.ne.0) B(ii)=B(ii)+f(l)
enddo
20    continue
!
 if (iecho.eq.1 .and. numadd.gt.0) write(stdout,2)
 2    format(1x,'Global stress loads are being added')
!
 return
 end
