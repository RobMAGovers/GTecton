 program skyline

! program for testing properties of the skyline storage model (akin, '82)
! for symmetric matrices

 dimension A(21),idiag(8)

 A( 1) = 11
 A( 2) = 12
 A( 3) = 22
 A( 4) = 23
 A( 5) = 33
 A( 6) = 14
 A( 7) = 24
 A( 8) = 34
 A( 9) = 44
 A(10) = 45
 A(11) = 55
 A(12) = 36
 A(13) = 46
 A(14) = 56
 A(15) = 66
 A(16) = 67
 A(17) = 77
 A(18) = 58
 A(19) = 68
 A(20) = 78
 A(21) = 88
 idiag(1) = 1
 idiag(2) = 3
 idiag(3) = 5
 idiag(4) = 9
 idiag(5) = 11
 idiag(6) = 15
 idiag(7) = 17
 idiag(8) = 21

100    write(*,1,advance='no')
 1    format("i,j > ")
 read(*,*) i,j

 indx = IjIndx (idiag,8,i,j)
 x = 0.
 if (indx.gt.0) x = A(indx)
 write(*,*) INT(x)
 goto 100
 end
!-------------------------------------------------------------------------------
 integer function IjIndx (IDIAG,NEQ,ii,jj)

! function returning index ij for skyline storage model

 integer ii,jj,NEQ
 dimension IDIAG(NEQ)
 lue=7

 if (ii.lt.1 .or. ii.gt.NEQ) then
write(lue,1) NEQ
 1        format("IjIndx fatal: 1 <= i <= ",I5)
call xit(1," ")
 endif
 if (jj.lt.1 .or. jj.gt.NEQ) then
write(lue,2) NEQ
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
