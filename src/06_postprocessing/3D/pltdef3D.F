subroutine pltdef

! Routine for setting program defaults
use iomodule
use filemodule,      only: BaseFedsk
use timestepmodule,  only: nstep
use postProcessData, only: deffac,setria,triang,deformed,opn
use monitorVars,     only: maxMon, elmon, nnel, nnp, npmon
use plotControl          ! a whole lot of it
use symbols,         only: IBcFix, IBcFor, IBcSld, IBcStn, IBcVel, IBcWnk

implicit none

#include "lenchar.i"



integer, external :: iflu, lnblk

ilu       = 0
stderr    = iflu('stderr')
stdin     = iflu('stdin')
stdout    = iflu('stdout')
luin      = 10
lucmd(ilu)= 11
luout(1)  = 12
luout(2)  = 13
lubuf     = 70
luMon     = stdout
lufo      = luout(1)
setria    = .true.
ifeout    = 0
append    = .false.
bufedit   = .false.
NSTEP     = 0
FEIN      = 'TECIN.DAT'
TEMPIN    = 'TEMPIN.DAT'
FEDSK     = 'FEDSK.DAT'
FEDSKRT   = lnblk(FEDSK)
BaseFedsk = trim(FEDSK)
readsy    = .true.
deformed  = .false.
stndev    = .true.
grnum     = .false.
npnum     = .false.
npbased   = .true.
viscnu    = .false.
ascin     = .false.
ascout    = .false.
sbhor     = .false.
visraw    = .false.
isvsrf    = -1
itemp     = 0
irelt     = 0
vislog    = .true.
ratlog    = .false.
verbose   = .true.
difstn    = .false.
BritDt    = .true.
Sscale    = 1d0
Escale    = 1d0
Fscale    = 1d0
Rscale    = 1d0
Dscale    = 1d0
Vscale    = 1d0
deffac    = 0d0
Ssb       = 1d0
Esb       = 1d0
Rsb       = 1d0
Dsb       = 1d0
Fsb       = 1d0
hatangle  = 40d0
hatsize   = 0.17d0
sysi      = 1.0d0
sysc      = .false.
ielbas    = 1
ielstp    = 1
inpbas    = 1
inpstp    = 1
iarcont   = 0
iEffD     = 1
iDisp     = 0
npmon     = .false.
elmon     = .false.
nnel      = 0
nnp       = 0
triang    = .false.
irate     = 0
viscli    = .false.
iDtfl     = 0
nstn      = 1
dispmn(1) = 0.d0
dispmn(2) = 0.d0
dispmx(1) = 0.d0
dispmx(2) = 0.d0
strsmn(1) = 0.d0
strsmn(2) = 0.d0
strsmx(1) = 0.d0
strsmx(2) = 0.d0
ratemn(1) = 0.d0
ratemn(2) = 0.d0
ratemx(1) = 0.d0
ratemx(2) = 0.d0
fluxmn(1) = 0.d0
fluxmn(2) = 0.d0
fluxmx(1) = 0.d0
fluxmx(2) = 0.d0
crdmin(1) = 0.d0
crdmin(2) = 0.d0
crdmin(3) = 0.d0
crdmax(1) = 0.d0
crdmax(2) = 0.d0
crdmax(3) = 0.d0
tempmn    = 0.d0
tempmx    = 0.d0
IBcFix    =  8
IBcSld    =  9
IBcWnk    = 13
IBcStn    = 16
IBcVel    = 15
IBcFor    = 17
Sunit     = 'Pa'
Eunit     = ' '
Funit     = 'W/m**2'
Runit     = '/s'
Dunit     = 'm'
g0        = 0.d0
icompens  = 0
Ygref     = 0
rhoa      = 3253.1d0
shrinkFactor = 0.9d0
opn       = .false.
T0calc    = .true.
maxMon    = 3000000


! gravity settings; courtesy of Taco
detpot        = .true.  ! compute gravitational potential
detgrav       = .false. ! compute gravity
detgrad       = .false. ! compute gravity gradients
coordtype     = 'sphe'  ! field point coordinate type: spherical or cartesian
detdilatation = .true.  ! determine contribution from density changes in
OpGravChange  = .true.  ! determine change in gravity (as opposed to absolute gravity)


return
end

