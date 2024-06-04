 subroutine addneqbc (SL,IEN,S,IA,IBOND,NEN,NEE,NDOF,MEQ)
!
! adds elements of local stiffness matrix to the system of equations which
! do not actually become part of the stiffness matrix
!
 implicit double precision (a-h,o-z)
!-pass
 dimension SL(NEE,NEE),IEN(NEN),S(MEQ,*),IA(*),IBOND(NDOF,*)
!-locl
 integer nn,n,i,leq,geq,mm,m,j,meq,heq
 ieq(i,n) = i +(n-1)*NDOF

 do nn=1,NEN
n = IEN(nn)
do i=1,NDOF
if (IBOND(i,n).eq.1 .or. IBOND(i,n).eq.2 .or. &
         IBOND(i,n).eq.5) then
leq = ieq(i,nn)
geq = IA(ieq(i,n))
do mm=1,NEN
    m = IEN(mm)
    do j=1,NDOF
    meq = ieq(j,mm)
    heq = ieq(j,m)
    S(geq,heq) = S(geq,heq) + SL(leq,meq)
    enddo
enddo
enddo
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine setneqbc(IBOND,IA,NDOF,NUMNP,MEQ)
!
! identifies equation numbers that are not actually generated (bcs)
!
!-pass
 integer NDOF,NUMNP,MEQ
 dimension IBOND(NDOF,*),IA(*)
!-locl
 integer n,i,j
!
 MEQ = 0
 do n=1,NUMNP
do i=1,NDOF
if (IBOND(i,n).eq.1 .or. IBOND(i,n).eq.2 .or. &
         IBOND(i,n).eq.5) then
    MEQ = MEQ + 1
    j = i+(n-1)*NDOF
    IA(j) = MEQ
endif
enddo
 enddo
!
 return
 end
!
!    IPIA +  NDOF*NUMNP
!    IPS  +  MEQ*NDOF*NUMNP*IWORD
!    IPFBC + MEQ*IWORD
!-------------------------------------------------------------------------------
 subroutine getneqbc (IBOND,S,IA,DELD,FBC,NDOF,NUMNP,MEQ)

 implicit double precision (a-h,o-z)
!-pass
 dimension S(MEQ,*),IA(*),DELD(NDOF,*),IBOND(NDOF,*),FBC(MEQ)
!-locl
 integer n,i,m,j,geq,heq
 ieq(i,n) = i +(n-1)*NDOF

 do n=1,NUMNP
do i=1,NDOF
if (IBOND(i,n).eq.1 .or. IBOND(i,n).eq.2 .or. &
         IBOND(i,n).eq.5) then
    geq = IA(ieq(i,n))
    FBC(geq) = 0D0
    do m=1,NUMNP
    do j=1,NDOF
        heq = ieq(j,m)
        FBC(geq) = FBC(geq) + S(geq,heq)*DELD(j,m)
    enddo
    enddo
endif
enddo
 enddo
!
 return
 end
!-------------------------------------------------------------------------------
 subroutine prtneqbc (IBOND,IA,FBC,NDOF,NUMNP)

 implicit double precision (a-h,o-z)
!-pass
 dimension IA(*),IBOND(NDOF,*),FBC(*)
!-locl
 integer n,i,geq
 ieq(i,n) = i +(n-1)*NDOF

 do n=1,NUMNP
do i=1,NDOF
if (IBOND(i,n).eq.1 .or. IBOND(i,n).eq.2 .or. &
         IBOND(i,n).eq.5) then
    geq = IA(ieq(i,n))
    write(747,1) i,n,FBC(geq)
 1            format(1x,'B(',I1,',',I5,') = ',1PG15.6)
endif
enddo
 enddo
!
 return
 end
