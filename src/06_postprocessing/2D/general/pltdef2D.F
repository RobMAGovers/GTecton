subroutine pltdef

! Routine for setting program defaults

use iomodule
use filemodule,      only: BaseFedsk
use timestepmodule,  only: nstep
use postProcessData, only: GMT,connout, &
                           deffac,setria,triang,deformed,CST, opn
use monitorVars,     only: maxMon, elmon, nnel, npmon, nnp, nsid, sidemon
use plotControl,         ! almost everything
use materialsModule, only: effsign, strinc
use symbols,         only: IBcFix, IBcFor, IBcOpn, IBcSh, IBcShr, &
                           IBcSld, IBcStn, IBcVel, IBcWnk


implicit none

#include "lenchar.i"


!-locl
integer, external :: iflu  ! from tools
integer, external :: lnblk ! from tools

ilu         = 0
stderr      = iflu('stderr')
stdin       = iflu('stdin')
stdout      = iflu('stdout')
luin        = 10
lucmd(ilu)  = 11  ! todo entry 0 of array is set, but we count from 1...
luout(1)    = 12
luout(2)    = 13
lubuf       = 70
luMon       = stdout
lufo        = 14
setria      = .true.
ifein       = 0
ifeout      = 0
append      = .false.
bufedit     = .false.
NSTEP       = 0
FEIN        = 'TECIN.DAT'
TEMPIN      = 'TEMPIN.DAT'
FEDSK       = 'FEDSK.DAT'
FEDSKRT     = lnblk(FEDSK)
BaseFedsk   = trim(FEDSK)
readsy      = .true.
deformed    = .false.
GMT         = .true.
connout     = .false.
sdev        = .true.
grnum       = .false.
npnum       = .false.
mattsym     = 1
viscnu      = .false.
ascin       = .false.
ascout      = .true. ! combines with GMT
sbhor       = .false.
visraw      = .false.
opn         = .false.
isvsrf      = -1
itemp       = 0
irelt       = 0
vislog      = .true.
ratlog      = .false.
verbose     = .true.
difstn      = .false.
BritDt      = .true.
Sscale      = 1d0
Escale      = 1d0
Fscale      = 1d0
Rscale      = 1d0
Dscale      = 1d0
Vscale      = 1d0
Velsc       = 1d0
deffac      = 1d0
Ssb         = 1D0
Esb         = 1D0
Rsb         = 1D0
Dsb         = 1D0
Vsb         = 1D0
Fsb         = 1d0
hatangle    = 40d0
hatsize     = 0.2d0
sysi        = 10d0
sysc        = .false.
ielbas      = 1
ielstp      = 1
inpbas      = 1
inpstp      = 1
iarcont     = 0
iDisp       = 0
npmon       = .false.
elmon       = .false.
sidemon     = .false.
nnel        = 0
nnp         = 0
nsid        = 0
triang      = .false.
irate       = 10
viscli      = .false.
iDtfl       = 0
nstn        = 1
dispmn(1:2) = 0d0
dispmx(1:2) = 0d0
strsmn(1:2) = 0d0
strsmx(1:2) = 0d0
ratemn(1:2) = 0d0
ratemx(1:2) = 0d0
fluxmn(1:2) = 0d0
fluxmx(1:2) = 0d0
crdmin(1:3) = 0d0
crdmax(1:3) = 0d0
tempmn      = 0d0
tempmx      = 0d0
IBcFix      =  8
IBcSld      =  9
IBcWnk      = 13
IBcStn      = 16
IBcVel      = 15
IBcFor      = 17
IBcSh       = 10
IBcOpn      = 12
IBcShr      = 45
Sunit       = 'Pa'
Eunit       = ' '
Funit       = 'W/m**2'
Runit       = '/s'
Dunit       = 'm'
Vunit       = 'm/s'
g0          = 0d0
icompens    = 0
Ygref       = 0
rhoa        = 3253.1d0
shrinkFactor = 0.8d0
CST         = .true.
opnscale    = 1.d0
opnaxnu     = 0.05d0
iopnflag    = 1
NUMZ        = 4
T0calc      = .true.
effsign     = .false.
strinc      = .false.
JSHELL      = 0
project     = .false.
gplon       = 0d0
gplat       = 0d0
gpaz        = 90d0
maxMon      = 100000

return
end subroutine
