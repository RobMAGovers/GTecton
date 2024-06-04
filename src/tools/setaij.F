 subroutine SETAIJ (I,J,IA,JA,NEQ,MAXJA)
!
! builds IA and JA.
!
 implicit none
!-pass
 integer I,J,IA,JA,NEQ,MAXJA
 dimension IA(*),JA(MAXJA)

!-locl
 integer k,m
 logical FirstEntry
 save FirstEntry
 data FirstEntry/.true./
!
 if (FirstEntry) then
do k=1,NEQ+1
IA(k) = k
enddo
do k=1,NEQ
JA(k) = k
enddo
FirstEntry = .false.
 endif
!
 if (I.ge.J) return
!
 if (IA(I+1).gt.IA(I)) then
k = IA(I)
do while (k.le.IA(I+1)-1)
if (JA(k).eq.J) return
k = k + 1
enddo
if (J.lt.JA(IA(I))) then
if (IA(NEQ+1).ge.MAXJA) goto 100
do k=IA(NEQ+1),IA(I),-1
    if (k.gt.0) JA(k+1) = JA(k)
enddo
JA(IA(I)) = J
else if (J.gt.IA(I+1)-1) then
if (IA(NEQ+1).ge.MAXJA) goto 100
do k=IA(NEQ+1),IA(I+1),-1
    if (k.gt.0) JA(k+1) = JA(k)
enddo
JA(IA(I+1)) = J
else
m = IA(I)
do while (m.lt.IA(I+1)-1.and.(J-JA(m))*(J-JA(m+1)).gt.0)
    m = m + 1
enddo
if (IA(NEQ+1).ge.MAXJA) goto 100
do k=IA(NEQ+1),m+1,-1
    if (k.gt.0) JA(k+1) = JA(k)
enddo
JA(m+1) = J
endif
 else
if (IA(NEQ+1).ge.MAXJA) goto 100
do k=IA(NEQ+1),IA(I),-1
if (k.gt.0) JA(k+1) = JA(k)
enddo
JA(IA(I)) = J
 endif
 do k=I+1,NEQ+1
IA(k) = IA(k) + 1
 enddo
!
 return
!
100    write(stderr,1) MAXJA
 1    format(/1x,'SETNZP: need more temporary memory than ',I9)
 call xit(1," ")
 end
