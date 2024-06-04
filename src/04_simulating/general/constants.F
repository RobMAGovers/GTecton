MODULE CONSTANTS
! physical and mathematical constants

integer,          parameter :: MB = 2**20 ! bytes in a megabyte

double precision, parameter :: YEAR = 31556736d0
   ! (60*60*24*365.24)


double precision, parameter :: re = 6378136d0 ! radius of the earth, in meters

double precision, parameter :: dfmin = 2.23d-308
double precision, parameter :: dfmax = 1.79d+308
double precision, parameter :: dlogmin = log10(dfmin)
double precision, parameter :: dlogmax = log10(dfmax)
double precision, parameter :: eps = 2.22d-16

integer,          parameter :: maxint = 2147483647   ! 2^32 - 1
integer,          parameter :: intmax = 2147483647   ! 2^32 - 1


double precision, parameter :: pi = 3.141592653589793d0
double precision, parameter :: twopi = 2d0 * pi
double precision, parameter :: halfpi = 0.5d0 * pi 
double precision, parameter :: deg2rad = pi/1.8d2
double precision, parameter :: rad2deg = 1.8d2/pi

double precision, parameter :: THF = -2.59807621135331601d0

! G indicates position of the support points of the
! shape functions. (not yet completely clear; TODO)
double precision, parameter :: G = 0.57735026918963d0

double precision, parameter :: grav_const = 6.6732D-11

! Temperatue in Kelvin at 0 C
double precision, parameter :: TK0 = 273.15

! universal gas constant, Joule  / (Mol Kelvin)
double precision, parameter :: R = 8.3144

! ten log of e
double precision, parameter :: TenLogE = .43429448d0


double precision, parameter :: root3 = 1.73205080756887719d0

integer, parameter :: six_int = 6

double precision, parameter :: zero  = 0d0
double precision, parameter :: half  = 0.5d0
double precision, parameter :: one   = 1d0
double precision, parameter :: two   = 2d0
double precision, parameter :: three = 3d0
double precision, parameter :: four  = 4d0
double precision, parameter :: six   = 6d0
double precision, parameter :: eight = 8d0

double precision, parameter :: third = one/three
double precision, parameter :: quart = one/four
double precision, parameter :: sixth = one/six
double precision, parameter :: eighth = one/eight
double precision, parameter :: threeq = three/four

double precision, parameter :: sqrtthree = sqrt(three)

integer,          parameter :: NNeighborsMax2D = 30
integer,          parameter :: NNeighborsMax3D = 100


!*************  enumerates *************

! to clarify the arguments to ADDSNE all over the code
integer,          parameter :: useDXE               = 0
integer,          parameter :: useDELDXE            = 1

! enumerates for the plot programs:
integer,          parameter :: SampleDisplacement   = 1
integer,          parameter :: SampleVelocity       = 2
integer,          parameter :: SampleTemperature    = 3
integer,          parameter :: SampleStress         = 4
integer,          parameter :: SampleStrain         = 5
integer,          parameter :: SampleRotation       = 6
integer,          parameter :: SampleStrainRate     = 7
integer,          parameter :: SampleHeatFlux       = 8
integer,          parameter :: SampleViscosity      = 9

! the mulitple ways of sampling, based on the arguments
integer,          parameter :: SampleOnCubicGrid    = 1
integer,          parameter :: SampleOnAnyLine      = 2
integer,          parameter :: SampleFromFile       = 3

! enumerates for pln itself, to see which variety is being used
integer,          parameter :: isPln                = 1
integer,          parameter :: isOpn                = 2
integer,          parameter :: isAxi                = 3
integer,          parameter :: isF3d                = 4

! enumerates for plnplt, to see which calls the documentation module,
! because parameters are different for each version? Think for example
! of stress xy, or rz
integer,          parameter :: isPlnplt             = 1
integer,          parameter :: isOpnplt             = 2
integer,          parameter :: isAxiplt             = 3
integer,          parameter :: isPlt3d              = 4

! due keep the calls to volumetric strain and strain rate apart
integer,          parameter :: strain               = 1
integer,          parameter :: strainRate           = 2

! to keep variable types apart
integer,          parameter :: int                  = 1
integer,          parameter :: float                = 2
integer,          parameter :: single               = 1
integer,          parameter :: double               = 2

! for plt3d, to differentiate between the different boundary conditions
integer,          parameter :: dispBC               = 1
integer,          parameter :: veloBC               = 2
integer,          parameter :: forceBC              = 3
integer,          parameter :: slideInPlaneBC       = 11
integer,          parameter :: slideOnLineBC        = 12
integer,          parameter :: fixedPointBC         = 13
integer,          parameter :: EulerAnglesBC        = 21

! for plt3d, to distinguish between 3 different types quantifications
integer,          parameter :: elemQualitySide      = 1
integer,          parameter :: elemQualitySurface   = 2
integer,          parameter :: elemQualityVolume    = 3

! the Gmsh constants
! (taken from: http://gmsh.info/doc/texinfo/gmsh.html)
integer,          parameter :: elemIsPoint          = 15
integer,          parameter :: elemIsLine           = 1
integer,          parameter :: elemIsTriangle       = 2
integer,          parameter :: elemIsTetrahedron    = 4

! element type
integer,          parameter :: elementTypeP1P0      = 10 ! default
integer,          parameter :: elementTypeP2P1      = 21

! enumerate for calling restart
integer,          parameter :: readRestart          = 0
integer,          parameter :: writeRestart         = 1


! displacement or velocity
integer,          parameter :: disp                 = 1
integer,          parameter :: velo                 = 2
! which kind of displacement to return in the plot programs
integer,          parameter :: total                = 1
integer,          parameter :: normal               = 2
integer,          parameter :: slip                 = 3
integer,          parameter :: split                = 4
! which coordinate
integer,          parameter :: coordx               = 1
integer,          parameter :: coordy               = 2
integer,          parameter :: coordz               = 3
integer,          parameter :: coordr               = 4 ! for Axi
integer,          parameter :: effective            = 5 ! Pythagoras

END MODULE
