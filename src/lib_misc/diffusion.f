    program main
*
    implicit none
    double precision k,rho,Cp,kappa,L,t
    character answer*1
*
10    write(*,1)
    1    format(1x,'Thermal conductivity [W m^-1 K^-1] > ',$)
    read(*,*,err=10,end=1000) k
20    write(*,2)
    2    format(1x,'Mass density [kg m^-3] > ',$)
    read(*,*,err=20,end=1000) rho
30    write(*,3)
    3    format(1x,'Specific heat [J kg^-1 K^-1] > ',$)
    read(*,*,err=30,end=1000) Cp
*
    kappa = k / (rho*Cp)
40    write(*,4) kappa
    4    format(1x,'Thermal diffusivity = ',1PE12.2,' [m^2 s^-1]'/
     >   1X,'Compute length (l) or time (t) scale > ',$)
    read(*,'(a1)',err=40,end=1000) answer
    if (answer.eq.'t') then
100        write(*,5)
        read(*,*,err=100,end=1000) L
    5        format(1x,'Length [m] > ',$)
        t = L**2/kappa
        write(*,6) t,t/3.1556736E+07
    6        format(1x,'t = ',1PE12.2,' sec, or ',1PG12.2,' year')
        goto 100
    else if (answer.eq.'l') then
200        write(*,7)
    7        format(1x,'Time [sec] > ',$)
        read(*,*,err=200,end=1000) t
        L = SQRT(t*kappa)
        write(*,8) L
    8        format(1x,'L = ',1PE12.2,' m')
        goto 200
    endif
*
1000    end
