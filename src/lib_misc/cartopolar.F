    subroutine cartopolar (x,y,z,R,phi,theta)

    implicit none

    double precision PI,ZERO
    parameter (PI=3.1415926535898d0,ZERO=0d0)
    double precision, parameter :: eps = 2.22d-16

!-pass
    double precision x,y,z,R,phi,theta
!-locl

    phi   = ZERO
    theta = ZERO
    R     = SQRT(x*x+y*y+z*z)
    if (R.le.eps) then
        return
    endif
    theta = ATAN2(SQRT(x*x+y*y),z)
    phi   = ATAN2(y,x)

    return
    end
