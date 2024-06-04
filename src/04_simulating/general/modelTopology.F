MODULE MODELTOPOLOGY
    
!    NSD            = Number of space dimensions
!    NDOF            = Number of degrees of freedom at each node
!    NEC = NSD*NEN           = Number of element coordinates
!    NEE = NDOF*NEN          = Number of mechanical element equations
!    NEN            = Number of nodes in an element

!    NSTR            = number of independent stress/strain tensor elements.
!                      in 2D, that is 3
!                      in 3D, that is 6
!
!    NSTRTAY        = number of independent Taylor stress tensor elements
!                   if scalar, that is 1
!                   if tensor
!                      in 2D, that is 3
!                      in 3D, that is 6
!
!    NSUR           = number of surface nodes


integer :: NSD,NDOF,NSTR,NEN, NSUR, nee !, NEC
integer :: nSidesPerElem

integer :: elementType

#ifdef EXPERIMENTALTRANSIENTDISCREEP
integer :: NSTRTAY 
#endif

END MODULE MODELTOPOLOGY
