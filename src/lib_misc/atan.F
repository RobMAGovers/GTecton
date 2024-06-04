double precision function ARCTAN (x,y)

! ATAN(y/x), returns [0,2pi>

implicit none

! tools may not depend on the constants module...
! hence explicitly copied
double precision, parameter :: pi = 3.141592653589793d0
double precision, parameter :: zero = 0d0
double precision, parameter :: two = 2d0
double precision, parameter :: eps = 1e-30


!-pass
double precision :: x,y
!-local
double precision :: atn
integer          :: lue
!external
integer, external          :: iflu
double precision, external :: sgn

!    eps = FLMIN()
lue = iflu('stderr')

if (ABS(x).le.eps) then
    if (ABS(y).le.eps) then

        write(*,*) 'arctan: 0/0 error'

        write(lue,1)
    1        format(1x,'ARCTAN: atan(0/0) error')
        ARCTAN = ZERO
        return
    else
        atn = sgn(y,eps) * pi / two
        if (atn.lt.ZERO) then
            atn = two * pi + atn
        endif
    endif
else
!        ATAN yields angle e [-pi/2,pi/2]
    atn = ATAN(y/x)
    if (ABS(y).lt.eps) then
        if (sgn(x,eps).lt.ZERO) then
            atn = atn + pi
        endif
    else
        if (sgn(y,eps).lt.ZERO) then
            if (sgn(x,eps).lt.ZERO) then
                atn = atn + pi
            else
                atn = two * pi + atn
            endif
        else
            if (sgn(x,eps).lt.ZERO) atn = atn + pi
        endif
    endif
endif

ARCTAN = ABS(atn)

return
end
!-------------------------------------------------------------------------------
double precision function sgn(x,eps)


implicit none

! tools may not depend on the constants module...
double precision, parameter :: zero = 0d0

double precision  :: x, eps
integer           :: lue
integer, external :: iflu

if (ABS(x).gt.eps) then
    sgn = x/ABS(x)
else
    lue = iflu('stderr')
    write(*,*) 'arctan: division by zero.'
    write(lue,1)
1        format(1x,'sgn: division by zero error')
    sgn = ZERO
endif

return
end
