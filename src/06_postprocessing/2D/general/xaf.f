subroutine str(STN)
*
c does not resolve 3d issue
*
parameter (NSTR=3)
dimension STN(NSTR)
s1 = 5d-1*(STN(1)-STN(2))
s1 = SQRT(s1*s1+STN(3)*STN(3))
aa = 5d-1*(STN(1)+STN(2))
s2 = aa - s1
s1 = aa + s1
