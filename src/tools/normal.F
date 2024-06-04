 program main
!
 implicit double precision (a-h,o-z)
 parameter (NSD=2, NDOF=2, NEN=4, NSTR=3)
 parameter (EPS=2.22D-16,TWO=2d0)
!-pass
 dimension DL(NDOF,NEN),XL(NSD,NEN),STN(NSTR)
dimension NSLIP(5),DIFORC(NDOF)
!-locl
 double precision n(2),tau(2),d(2)
!
10    write(*,1,advance='no')
 1    format(1x,'Sxx, Syy, Sxy > ')
 read(*,*,err=10,end=40) STN(1),STN(2),STN(3)
20    write(*,2,advance='no')
 2    format(1x,'Face coordinate 1 (x,y) > ')
 read(*,*,err=20,end=40) XL(1,1),XL(2,1)
30    write(*,3,advance='no')
 3    format(1x,'Face coordinate 2 (x,y) > ')
 read(*,*,err=30,end=40) XL(1,2),XL(2,2)
!
 dx = XL(1,2)-XL(1,1)
 dy = XL(2,2)-XL(2,1)
 s = SQRT(dx*dx+dy*dy)
 if (s.lt.EPS) stop 'identical coordinates'
 n(1) = dy/s
 n(2) = -dx/s
 d(1) = dx/s
d(2) = dy/s
 sn = STN(1)*n(1)*n(1)+STN(2)*n(2)*n(2)+TWO*STN(3)*n(1)*n(2)
 tn = STN(1)*n(1)*d(1)+STN(3)*n(2)*d(1)+STN(3)*n(1)*d(2)+ &
   STN(2)*n(2)*d(2)
 write(*,4) sn,tn
 4    format(1x,'normal stress = ',1PG12.4,', shear stress = ', &
   1PG12.4/)
 goto 10
40    end
