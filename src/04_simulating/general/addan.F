double precision function AddAn (X,D,NDOF,NSD)

USE ALGEBRA,   only: EQUATE, UPDATE
USE CONSTANTS, only: pi

implicit none

!-pass
integer          :: NDOF,NSD
double precision :: X(NSD),D(NDOF)
!-locl
double precision :: XL(3)
!-init
double precision, parameter :: h   = -105D+3
double precision, parameter :: a   =  70D+3
double precision, parameter :: dT0 =  50D0

call EQUATE (XL,X,NSD)
call UPDATE (XL,D,NDOF)

if (XL(2).gt.h) then
    AddAn = dT0 * EXP(-(XL(1)/a)**2) * &
     1.3D0*SIN(XL(2)/h*pi)*SIN((h-XL(2))*0.5D0/h*pi) &
       + dT0 * EXP(-((350D3-XL(1))/a)**2) * &
     1.3D0*SIN(XL(2)/h*pi)*SIN((h-XL(2))*0.5D0/h*pi)
else
    AddAn = 0D0
endif

return
end
