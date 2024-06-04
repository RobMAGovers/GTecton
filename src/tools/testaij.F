 program main

 parameter (NJ=9,NEQ=5,NI=NEQ+1,MAX=NEQ*NEQ)
 dimension IA(NI),JA(MAX)
!-init
 do i=1,MAX
JA(i) = 0
 enddo
 call setaij (5,5,IA,JA,NEQ,MAX)
 write(*,1) (i,i=1,NJ)
 write(*,2) (IA(i),i=1,NI)
 write(*,3) (JA(i),i=1,NJ)
 call setaij (3,4,IA,JA,NEQ,MAX)
 call setaij (4,3,IA,JA,NEQ,MAX)
 write(*,1) (i,i=1,NJ)
 write(*,2) (IA(i),i=1,NI)
 write(*,3) (JA(i),i=1,NJ)
 call setaij (1,1,IA,JA,NEQ,MAX)
 write(*,1) (i,i=1,NJ)
 write(*,2) (IA(i),i=1,NI)
 write(*,3) (JA(i),i=1,NJ)
 call setaij (2,2,IA,JA,NEQ,MAX)
 write(*,1) (i,i=1,NJ)
 write(*,2) (IA(i),i=1,NI)
 write(*,3) (JA(i),i=1,NJ)
 call setaij (5,4,IA,JA,NEQ,MAX)
 call setaij (4,5,IA,JA,NEQ,MAX)
 write(*,1) (i,i=1,NJ)
 write(*,2) (IA(i),i=1,NI)
 write(*,3) (JA(i),i=1,NJ)
 call setaij (3,1,IA,JA,NEQ,MAX)
 call setaij (1,3,IA,JA,NEQ,MAX)
 write(*,1) (i,i=1,NJ)
 write(*,2) (IA(i),i=1,NI)
 write(*,3) (JA(i),i=1,NJ)
 call setaij (1,4,IA,JA,NEQ,MAX)
 call setaij (4,1,IA,JA,NEQ,MAX)
 write(*,1) (i,i=1,NJ)
 write(*,2) (IA(i),i=1,NI)
 write(*,3) (JA(i),i=1,NJ)
 call setaij (3,3,IA,JA,NEQ,MAX)
 write(*,1) (i,i=1,NJ)
 write(*,2) (IA(i),i=1,NI)
 write(*,3) (JA(i),i=1,NJ)
 call setaij (4,4,IA,JA,NEQ,MAX)
 write(*,1) (i,i=1,NJ)
 write(*,2) (IA(i),i=1,NI)
 write(*,3) (JA(i),i=1,NJ)
 call setaij (4,5,IA,JA,NEQ,MAX)
 call setaij (5,4,IA,JA,NEQ,MAX)
 write(*,1) (i,i=1,NJ)
 write(*,2) (IA(i),i=1,NI)
 write(*,3) (JA(i),i=1,NJ)

 1    format(/4x,'- ',9(1X,I2)/ &
   '---------------------------------------------')
 2    format(1X,'IA - ',6(1X,I2))
 3    format(1X,'JA - ',9(1X,I2))

 NSIZEA = 2*(IA(NEQ+1)-(NEQ+1)) + NEQ
 NSIZEA = 2*IA(NEQ+1) - NEQ - 2
 do i=1,NEQ
JA(IA(i)) = 0
 enddo
 do i=1,NEQ-1
do jj=IA(i)+1,IA(i+1)-1
j = JA(jj)
JA(IA(j)) = JA(IA(j)) + 1
enddo
 enddo
 write(*,*) 'NSIZEA = ',NSIZEA
 do i=1,NEQ
IA(i) = JA(IA(i)) + IA(i+1)-IA(i)
write(*,4) i,IA(i)
 4        format(1x,'NNZ(',I2,') = ',I2)
 enddo

 end
!-------------------------------------------------------------------------------
 subroutine xit(i)
 stop
 end
