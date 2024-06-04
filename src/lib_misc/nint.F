integer function nint(x)

implicit none

double precision :: x

double precision, parameter :: HALF=0.5D0
double precision, parameter :: XNMAX = 2147483647D0

if (x.ge.0.) then
    if (x.gt.XNMAX) then
        write(0,*) 'nint overflow: ',x
        call exitp(1)
    endif
    nint = INT(x+HALF)
else
    if (x.lt.-XNMAX) then
        write(0,*) 'nint overflow: ',x
        call exitp(1)
    endif
    nint = INT(x-HALF)
endif

return
end
