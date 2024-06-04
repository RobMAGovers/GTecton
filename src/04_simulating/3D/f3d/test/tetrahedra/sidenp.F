 subroutine SIDENP (is,QUAD,lien)
!
! routine sets (local) nodal point numbers belonging to side IS.
! NPs in lien are positively oriented.
!
 implicit none
!-pass
 integer is,lien
 logical QUAD
 dimension lien(3)
!
 if (is.eq.1) then
lien(1)=1
lien(2)=2
lien(3)=4
 elseif (is.eq.2) then
lien(1)=2
lien(2)=3
lien(3)=4
 elseif (is.eq.3) then
lien(1)=3
lien(2)=1
lien(3)=4
 elseif (is.eq.4) then
lien(1)=3
lien(2)=2
lien(3)=1
 else
write(stderr,1) is
 1        format(1x,'SIDENP: a tetrahedra has 4 sides, not ',I2)
call xit(1," ")
 endif
!
 return
 end
