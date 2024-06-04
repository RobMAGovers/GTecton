 subroutine CHKCCY (A,IDIAG,D,NEQ,IFLAG)
 implicit double precision (a-h,o-z)
 dimension A(*),IDIAG(*),B(*)
 save B
 dimension B(2000),RES(2000)
!
 if (iflag.eq.0) then
do i=1,NEQ
B(i) = D(i)
enddo
return
 endif
!
 sumsq = 0d0
 do i=1,NEQ
    res(i) = 0d0
    do j=1,NEQ
k = IJINDX(IDIAG,NEQ,i,j)
if (k.eq.0) then
    res(i) = res(i) - B(j)
else if (k.gt.0) then
    res(i) = res(i) + A(k)*D(j) - B(j)
else
    write(0,*) 'FATAL: negative IJINDX'
    call exitp(1)
endif
enddo
sumsq = sumsq + res(i)*res(i)
 enddo
!
 sumsq = SQRT(sumsq)
 write(0,*) 'CHKCCY: solver residu = ',sumsq
 return
 end
!------------------------------------------------------------------------------
 integer function IjIndx (IDIAG,NEQ,ii,jj)

! function returning index ij for skyline storage model

 integer ii,jj,NEQ
 dimension IDIAG(NEQ)

 if (ii.lt.1 .or. ii.gt.NEQ) then
write(0,1) NEQ
 1        format("IjIndx fatal: 1 <= i <= ",I5)
call xit(1," ")
 endif
 if (jj.lt.1 .or. jj.gt.NEQ) then
write(0,2) NEQ
 2        format("IjIndx fatal: 1 <= j <= ",I5)
call xit(1," ")
 endif

 if (ii.gt.jj) then
i = jj
j = ii
 else
i = ii
j = jj
 endif

 if (j.eq.1) then
IjIndx = 1
 else
jh = IDIAG(j)-IDIAG(j-1)
imin = j - jh + 1
if (i.lt.imin) then
IjIndx = 0
else
IjIndx = IDIAG(j) - (j-i)
endif
 endif
 return
 end
