double precision function temp0 (X,dummy,NSD)
! used by rdthrm.F

implicit none

!double precision, parameter :: one=1d0
!-pass
integer NSD
double precision :: X(NSD)
double precision :: dummy

#include "coord.i"
!-init
!temp0 = interpolated
temp0 = 0d0

return
end
