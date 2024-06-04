double precision function effs (S)
!
! Computes square root of second invariant of deviatoric part of tensor S
!
use constants, only: dfmin

implicit none
!-pass
double precision :: S
dimension S(*)
!-locl
double precision, external :: AJ2

effs = AJ2(S)

if (effs.gt.dfmin) then
    effs = SQRT(effs)
endif

return
end
!-------------------------------------------------------------------------------
double precision function etaLin (effstr,emhu,anpwr)
!
! calculates linear viscosity equivalent
!
use constants, only: dfmin, dfmax

implicit none
!-pass
double precision :: emhu,anpwr,effstr

!-locl
double precision :: xpow, eta

xpow   = anpwr - 1d0
if (xpow.gt.1E-6) then
    if (effstr.gt.dfmin) then
        eta = effstr**(-xpow) * emhu**anpwr
    else
        eta = dfmax*1d-1
    endif
else
    eta = emhu
endif

etaLin = eta

return
end
