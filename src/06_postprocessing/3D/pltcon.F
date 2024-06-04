subroutine pltcon ()

USE MODELTOPOLOGY,   only: NSD,NDOF,NSTR,NEN
use modeldefinition, only: ntype
use plotControl,     only: nplt, npldof

implicit none

call INTCON (NSD,NDOF,NSTR,NEN)
NPLT   = 3
NPLDOF = 3

return
end subroutine
