	program main
*
	implicit none
	double precision E,v,eta,mhu,maxwell,year
	year = 3.1556736D+07
*
	write(*,1)
   1	format(1x,'Young''s modulus (Pa) > ',$)
	read(*,*,end=100,err=100) E
	write(*,2)
   2	format(1x,'Poisson''s ratio > ',$)
	read(*,*,end=100,err=100) v
	write(*,3)
   3	format(1x,'Viscosity (Poise) > ',$)
	read(*,*,end=100,err=100) eta
*
	mhu = E/(2d0*(1d0+v))
	maxwell = eta/mhu
	write(*,4) maxwell,maxwell/year
   4	format(1x,'Maxwell time = ',1PG14.6,' sec, ',1PG14.6,' year')
100	end
