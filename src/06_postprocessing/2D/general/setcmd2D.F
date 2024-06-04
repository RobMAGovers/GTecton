subroutine setcmd (cmd,ncmd,ierr)

USE MODELDEFINITION, only: NUMEL
USE CONVERSIONS,     only: string2int
use debugmodule,     only: iecho
use timestepmodule,  only: chktim
use filemodule,      only: do_soelim
use iomodule,        only: verbose, ifein, ifeout, ascout, stderr, append
use modeltopology,   only: nen
use postProcessData, only: GMT, cst, deffac, deformed, opn, setria, triang
use plotControl     !only: a lot
use errors,          only: exec, ok, syntax
use materialsModule, only: effsign, strinc
use symbols,         only: IBcFix, IBcFor, IBcOpn, IBcSh, IBcShr, &
                           IBcSld, IBcStn, IBcVel, IBcWnk, maxSym, &
                           mfix, mforc, mslid, mstrs, mvelo, mwink, &
                           syfix, syforc, syslid, systrs, syvelo, sywink

implicit none
!-pass
 integer ncmd,ierr
 character(len=*) cmd
 dimension cmd(*)

#include "lenchar.i"


!-locl
logical, external :: numeric
logical :: Mech,Temp
integer, external :: lnblk
double precision, external :: chreal

integer :: i, ireq, irold, iside, itmp
integer :: l4, l7, n, m

real    :: tmp, u1, u2, x

!-init

!----------------------------------------------------------
!#    se[t] parameter [value [value]]
 1    format(1x,'Usage: se[t] parameter [value [value]]')
 if (ncmd.eq.1) then
write(stderr,1)
ierr = SYNTAX
goto 1000
 endif

!----------------------------------------------------------
!#    se[t] sy[mbolsize] size
if (cmd(2)(1:2).eq.'sy') then
 5        format(1x,'Usage: se[t] sy[mbolsize] size')
    if (ncmd.lt.3) then
        write(stderr,5)
        ierr = SYNTAX
goto 1000
endif
if (.not.numeric(cmd(3))) then
write(stderr,5)
ierr = SYNTAX
goto 1000
endif
tmp = chreal(cmd(3))
if (tmp.le.0.0) then
write(stderr,10)
10        format(1x,'positive symbol size required')
ierr = SYNTAX
goto 1000
endif
sysi=tmp
goto 1000
 endif

!----------------------------------------------------------
!#    se[t] stre[ss] value1 unit1 = value2 unit2
!#    se[t] stra[in] value1 unit1 = value2 unit2
!#    se[t] ra[te] value1 unit1 = value2 unit2
!#    se[t] di[splacement] value1 unit1 = value2 unit2
!#    se[t] fl[ux] value1 unit1 = value2 unit2
!#    se[t] vel[ocity] value1 unit1 = value2 unit2
 n = 0
 if (cmd(2)(1:4).eq.'stre') then
n = 1
 else if (cmd(2)(1:4).eq.'stra') then
n = 2
 else if (cmd(2).eq.'ra' .or. cmd(2)(1:3).eq.'rat') then
n = 3
 else if (cmd(2)(1:2).eq.'di') then
n = 4
 else if (cmd(2)(1:2).eq.'fl') then
n = 5
 else if (cmd(2)(1:3).eq.'vel') then
n = 6
 endif
15   format(1x,'Usage: se[t] stre[ss] value1 unit1 = value2',1x, &
 'unit2')
16   format(1x,'Usage: se[t] stre[ss] factor')
20   format(1x,'Usage: se[t] stra[in] value1 unit1 = value2',1x, &
 'unit2')
21   format(1x,'Usage: se[t] stra[in] factor')
25   format(1x,'Usage: se[t] ra[te] value1 unit1 = value2',1x, &
 'unit2')
26   format(1x,'Usage: se[t] ra[te] factor')
30   format(1x,'Usage: se[t] di[sp] value1 unit1 = value2',1x, &
 'unit2')
31   format(1x,'Usage: se[t] di[splacement] factor')
32   format(1x,'Usage: se[t] fl[ux] value1 unit1 = value2',1x, &
 'unit2')
33   format(1x,'Usage: se[t] fl[ux] factor')
34   format(1x,'Usage: se[t] vel[ocity] value1 unit1 = value2 ', &
 'unit2')
35   format(1x,'Usage: se[t] vel[ocity] factor')
if (n.gt.0) then
   if (GMT) then
if (ncmd.lt.3) then
  if (n.eq.1) write(stderr,16)
  if (n.eq.2) write(stderr,21)
  if (n.eq.3) write(stderr,26)
  if (n.eq.4) write(stderr,31)
  if (n.eq.5) write(stderr,33)
  if (n.eq.6) write(stderr,35)
  ierr = SYNTAX
  goto 1000
endif
if (.not. numeric(cmd(3)) ) then
  if (n.eq.1) write(stderr,16)
  if (n.eq.2) write(stderr,21)
  if (n.eq.3) write(stderr,26)
  if (n.eq.4) write(stderr,31)
  if (n.eq.5) write(stderr,33)
  if (n.eq.6) write(stderr,35)
  ierr = SYNTAX
  goto 1000
endif
u1 = chreal(cmd(3))
if (u1.lt.1D-20) then
  write(stderr,36)
36          format(1x,'set: use scaling factor > 0')
  ierr = SYNTAX
  goto 1000
endif
if (n.eq.1) then
  Sscale = u1
else if (n.eq.2) then
  Escale = u1
else if (n.eq.3) then
  Rscale = u1
else if (n.eq.4) then
  Dscale = u1
else if (n.eq.5) then
  Fscale = u1
else if (n.eq.6) then
  Velsc = u1
endif
   else ! no GMT
if (ncmd.lt.7) then
if (n.eq.1) write(stderr,15)
if (n.eq.2) write(stderr,20)
if (n.eq.3) write(stderr,25)
if (n.eq.4) write(stderr,30)
if (n.eq.5) write(stderr,32)
if (n.eq.6) write(stderr,34)
ierr = SYNTAX
goto 1000
endif
if (.not. numeric(cmd(3)) ) then
if (n.eq.1) write(stderr,15)
if (n.eq.2) write(stderr,20)
if (n.eq.3) write(stderr,25)
if (n.eq.4) write(stderr,30)
if (n.eq.5) write(stderr,32)
if (n.eq.6) write(stderr,34)
ierr = SYNTAX
goto 1000
endif
u1 = chreal(cmd(3))
if (u1.lt.1D-20) then
write(stderr,36)
ierr = SYNTAX
goto 1000
endif
if (.not. numeric(cmd(6)) ) then
if (n.eq.1) write(stderr,15)
if (n.eq.2) write(stderr,20)
if (n.eq.3) write(stderr,25)
if (n.eq.4) write(stderr,30)
if (n.eq.5) write(stderr,32)
if (n.eq.6) write(stderr,34)
ierr = SYNTAX
goto 1000
endif
u2 = chreal(cmd(6))
if (u2.lt.1D-20) then
write(stderr,36)
ierr = SYNTAX
goto 1000
endif
if (n.eq.1) then
Sscale = u2/u1
Ssb    = u1
Sunit  = cmd(4)
else if (n.eq.2) then
Escale = u2/u1
Esb    = u1
Eunit  = cmd(4)
else if (n.eq.3) then
Rscale = u2/u1
Rsb    = u1
Runit  = cmd(4)
else if (n.eq.4) then
Dscale = u2/u1
Dsb    = u1
Dunit  = cmd(4)
else if (n.eq.5) then
Fscale = u2/u1
Fsb    = u1
Funit  = cmd(4)
else if (n.eq.6) then
Velsc = u2/u1
Vsb    = u1
Vunit  = cmd(4)
endif
if (iecho.eq.1) then
l4 = lnblk(cmd(4))
l7 = lnblk(cmd(7))

if (n.eq.1) then
    write(stderr,39) Sscale,cmd(7)(1:l7),'/', &
           cmd(4)(1:l4)
39            format(1x,'stress scale factor = ',1PG12.6,1X,80a)
    write(stderr,40) Ssb,cmd(4)(1:l4)
40            format(1x,'stress scale bar    = ',1PG12.6,1X,80a)
else if (n.eq.2) then
    write(stderr,45) Escale,cmd(7)(1:l7),'/', &
           cmd(4)(1:l4)
45            format(1x,'strain scale factor = ',1PG12.6,1X,80a)
    write(stderr,50) Esb,cmd(4)(1:l4)
50            format(1x,'strain scale bar    = ',1PG12.6,1X,80a)
else if (n.eq.3) then
    write(stderr,55) Rscale,cmd(7)(1:l7),'/', &
           cmd(4)(1:l4)
55            format(1x,'rate scale factor = ',1PG12.6,1X,80a)
    write(stderr,60) Rsb,cmd(4)(1:l4)
60            format(1x,'rate scale bar    = ',1PG12.6,1X,80a)
else if (n.eq.4) then
    write(stderr,65) Dscale,cmd(7)(1:l7),'/', &
           cmd(4)(1:l4)
65            format(1x,'displ. scale factor = ',1PG12.6,1X,80a)
    write(stderr,70) Dsb,cmd(4)(1:l4)
70            format(1x,'displ. scale bar    = ',1PG12.6,1X,80a)
else if (n.eq.5) then
    write(stderr,71) Fscale,cmd(7)(1:l7),'/', &
           cmd(4)(1:l4)
71            format(1x,'flux scale factor = ',1PG12.6,1X,80a)
    write(stderr,72) Fsb,cmd(4)(1:l4)
72            format(1x,'flux scale bar    = ',1PG12.6,1X,80a)
else if (n.eq.6) then
    write(stderr,73) Velsc,cmd(7)(1:l7),'/', &
           cmd(4)(1:l4)
73            format(1x,'velocity scale factor = ',1PG12.6,1X,80a)
    write(stderr,74) Vsb,cmd(4)(1:l4)
74            format(1x,'velocity scale bar    = ',1PG12.6,1X,80a)
endif
endif
   endif
   goto 1000
 endif

!----------------------------------------------------------
!#    se[t] ar[rowhead] size angle [no[co[nt]]]
 if (cmd(2)(1:2).eq.'ar') then
75       format(1x,'Usage: se[t] ar[rowhead] size angle',1x, &
     '[no[co[nt]]]')

if (ncmd.lt.4) then
    ierr = SYNTAX
    write(stderr,75)
    goto 1000
endif

if (.not.numeric(cmd(3))) then
    ierr = SYNTAX
    write(stderr,75)
    goto 1000
endif

tmp = chreal(cmd(3))

if (tmp.le.0.E0) then
    ierr = SYNTAX
    write(stderr,75)
    goto 1000
endif

hatsize = tmp

if (.not.numeric(cmd(4))) then
    ierr = SYNTAX
    write(stderr,75)
    goto 1000
endif

tmp = chreal(cmd(4))

if (tmp.le.0.E0 .or. tmp.ge.9.E1) then
    ierr = SYNTAX
    write(stderr,80)
80        format(1x,'0 < angle < 90 degrees')
    goto 1000
endif

hatangle = tmp

if (ncmd.gt.4) then
    if (cmd(5)(1:2).eq.'co') then
        iarcont = 1
    else if (cmd(5)(1:4).eq.'noco') then
        iarcont = 0
    else
        write(stderr,75)
        ierr = SYNTAX
        goto 1000
    endif
endif

goto 1000
 endif

!----------------------------------------------------------
!# se[t] [no]ec[ho] [number]
 if ( cmd(2)(1:2).eq.'ec' ) then
85        format(1x,'Usage: se[t] [no]ec[ho] [number]')

if ( ncmd.gt.2 ) then
    if (.not. numeric(cmd(3)) ) then
        ierr = SYNTAX
        write(stderr,85)
        goto 1000
    endif
    itmp = string2int(cmd(3))
!        itmp = INT( chreal(cmd(3)) )
    iecho = itmp
!        iecho=3 is magic to tecin.f
    if (iecho.eq.3) then
        iecho = 4
    endif
else
    iecho = 1
endif
goto 1000
 elseif ( cmd(2)(1:4).eq.'noec' ) then
iecho = 0
goto 1000
 endif

!----------------------------------------------------------
!#    se[t] whether to do soelim or not
 if (cmd(2)(1:4).eq.'noso') then
! noso for no soelim
if (ncmd.gt.2) then
    ierr = SYNTAX
    write(stderr,75)
    goto 1000
endif

do_soelim = .false.
goto 1000

 endif

 if (cmd(2)(1:4).eq.'doso') then
! doso for do soelim

if (ncmd.gt.2) then
    ierr = SYNTAX
    write(stderr,75)
    goto 1000
endif

do_soelim = .true.
goto 1000

 endif

!----------------------------------------------------------
!# se[t] it[ime] number
 if ( cmd(2)(1:2).eq.'it' ) then

90        format(1x,'Usage: se[t] it[ime] number|increment')

if (ncmd.lt.3) then
    ierr = SYNTAX
    write(stderr,90)
    goto 1000
endif

if (ifein.eq.0) then
    write(stderr,92)
92        format(1x,'set itime: read FE-input file first')
    ierr = EXEC
    goto 1000
endif

if (ifeout.eq.0) then
    write(stderr,93)
93        format(1x,'set itime: read FE-output file first')
    ierr = EXEC
    goto 1000
endif

if (numeric(cmd(3))) then

    IREQ = string2int(cmd(3))

    call chktim (IREQ,Mech, Temp)

    if ((.not.Mech) .and. (.not.Temp)) then
        write(stderr,95) IREQ
95            format(1x,'set itime: no output at itime = ',I4)
        ierr = EXEC
        goto 1000
    endif

else

    if (cmd(3).eq.'++') then
        call timinc (IREQ, .true.,ierr)
    else if (cmd(3).eq.'--') then
        call timinc (IREQ, .false.,ierr)
    else
        write(stderr,90)
        ierr = SYNTAX
    endif

    if (ierr.ne.OK) goto 1000

endif

call FOUTRD (IREQ,ierr)

goto 1000

 endif

!----------------------------------------------------------
!# se[t] [no]gmt
 if (cmd(2)(1:3).eq.'gmt') then
GMT = .true.
ascout = .true.
goto 1000
 elseif (cmd(2)(1:5).eq.'nogmt') then
GMT = .false.
ascout = .false.
goto 1000
 endif

!----------------------------------------------------------
!# se[t] [no]mattsym [g|c] 
 if (cmd(2)(1:7).eq.'mattsym') then
   if (ncmd.gt.2) then
      if (cmd(3)(1:1).eq.'g') then
         mattsym = 2
      else if (cmd(3)(1:1).eq.'c') then
         mattsym = 1
      else
         write(stderr,96)
96       format('Usage: set mattsym [g|c]')
         ierr = SYNTAX
       endif
   else
       mattsym = 1
   endif
   goto 1000
 elseif (cmd(2)(1:9).eq.'nomattsym') then
mattsym = 0
goto 1000
 endif

!----------------------------------------------------------
!# se[t] [no]de[formed] [amplification]
 if (cmd(2)(1:3).eq.'de ' .or. cmd(2)(1:3).eq.'def') then
97        format(1x,'Usage: se[t] [no]de[formed] [amplification]')
deformed = .true.
setria = .true.
if (ncmd.gt.2) then
if (.not.numeric(cmd(3))) then
    write(stderr,97)
    ierr = SYNTAX
    goto 1000
endif
x = chreal(cmd(3))
if (x.le.0. .or. x.gt.1e30) then
    write(stderr,98)
98               format(1x,'set deformed: amplification out of',1x, &
             'range')
    ierr = SYNTAX
    goto 1000
endif
deffac = x
endif
if (JSHELL.ne.0) write(stderr,99)
99       format(1x,'set deformed: spherical shell option on:'/1x, &
     'geographic coordinates and Carthesian displacements will', &
     ' be added')
goto 1000
 elseif (cmd(2)(1:4).eq.'node') then
deformed = .false.
goto 1000
 endif

!----------------------------------------------------------
!# se[t] [no]as[cii]
 if (cmd(2)(1:2).eq.'as') then
ascout = .true.
goto 1000
 elseif (cmd(2)(1:4).eq.'noas') then
if (GMT) then
write(stderr,101)
  101        format(1x,'set noascii: not for GMT')
else
    ascout = .false.
endif
goto 1000
 endif

!----------------------------------------------------------
!# se[t] [no]bi[nary]
 if (cmd(2)(1:2).eq.'bi') then
if (GMT) then
write(stderr,102)
  102        format(1x,'set binary: not for GMT')
else
    ascout = .false.
endif
goto 1000
 elseif (cmd(2)(1:4).eq.'nobi') then
ascout = .true.
goto 1000
 endif

!----------------------------------------------------------
!# se[t] elmi[nterval] [base] number
 if (cmd(2)(1:4).eq.'elmi') then
  100        format(1x,'Usage: se[t] elmi[nterval] [base] number')

if (ncmd.lt.3) then
    write(stderr,100)
    ierr = SYNTAX
    goto 1000
endif

if (.not.numeric(cmd(3))) then
    write(stderr,100)
    ierr = SYNTAX
    goto 1000
endif

n = string2int(cmd(3))

if (ncmd.gt.3) then
    if (.not.numeric(cmd(4))) then
        write(stderr,100)
        ierr = SYNTAX
        goto 1000
    endif
    m = string2int(cmd(4))
else
    m = n
    n = ielbas
endif

if (n.lt.1 .or. n.gt.NUMEL) then
    write(stderr,105)
  105        format(1x,'set elmi: base number out of range')
    ierr = EXEC
    goto 1000
endif

if (m.lt.1 .or. m.gt.NUMEL) then
    write(stderr,110)
  110        format(1x,'set elmi: increment out of range')
    ierr = EXEC
    goto 1000
endif

ielbas = n
ielstp = m
goto 1000
 endif

!----------------------------------------------------------
!# se[t] npi[nterval] [base] number

 if (cmd(2)(1:3).eq.'npi') then
  115        format(1x,'Usage: se[t] npi[nterval] [base] number')

if (ncmd.lt.3) then
    write(stderr,115)
    ierr = SYNTAX
    goto 1000
endif

if (.not.numeric(cmd(3))) then
    write(stderr,115)
    ierr = SYNTAX
    goto 1000
endif

n = string2int(cmd(3))

if (ncmd.gt.3) then
    write(*,*) 'ncmd > 3'
    if (.not.numeric(cmd(4))) then
        write(stderr,115)
        ierr = SYNTAX
        goto 1000
    endif
    m = string2int(cmd(4))
else
    m = n
    n = ielbas
endif


if (n.lt.1 .or. n.gt.NUMEL) then
    write(stderr,120)
  120        format(1x,'set npi: base number out of range')
    ierr = EXEC
    goto 1000
endif

if (m.lt.1 .or. m.gt.NUMEL) then
    write(stderr,125)
  125        format(1x,'set npi: increment out of range')
    ierr = EXEC
    goto 1000
endif

inpbas = n
inpstp = m
goto 1000
 endif

 if (cmd(2)(1:4).eq.'visc') then
  130        format(1x,'Usage: se[t] visc[ale] factor')

if (ncmd.lt.3) then
    write(stderr,130)
    ierr = SYNTAX
    goto 1000
endif

if (.not. numeric(cmd(3))) then
    write(stderr,130)
    ierr = SYNTAX
    goto 1000
endif

x = chreal(cmd(3))

if (x.le.0. .or. x.ge.1E30) then
    write(stderr,135)
  135        format(1x,'set viscale: scaling factor out of range')
    ierr = SYNTAX
    goto 1000
endif

Vscale = x
goto 1000
 endif

 if (cmd(2)(1:4).eq.'bcsy') then
  140        format(1x,'Usage: se[t] bcsy[mbol] type number')

if (ncmd.lt.4) then
    write(stderr,140)
    ierr = SYNTAX
    goto 1000
endif

if (.not. numeric(cmd(4)) ) then
    write(stderr,140)
    ierr = SYNTAX
    goto 1000
endif

n = string2int(cmd(4))

if (n.lt.1 .or. n.gt.40) then
    write(stderr,145)
  145        format(1x,'set bcsy: symbol number out of range')
    ierr = EXEC
    goto 1000
endif

if (cmd(3)(1:2).eq.'fi') then
IBcFix = n
if (.not.readsy) call rdsymb(IBcFix, SYFIX, MFIX,MAXSYM)
elseif (cmd(3)(1:2).eq.'sl') then
IBcSld = n
if (.not.readsy) call rdsymb(IBcSld,SYSLID,MSLID,MAXSYM)
elseif (cmd(3)(1:2).eq.'wi') then
IBcWnk = n
if (.not.readsy) call rdsymb(IBcWnk,SYWINK,MWINK,MAXSYM)
elseif (cmd(3)(1:2).eq.'ve') then
IBcVel = n
if (.not.readsy) call rdsymb(IBcVel,SYVELO,MVELO,MAXSYM)
elseif (cmd(3)(1:2).eq.'st') then
IBcStn = n
if (.not.readsy) call rdsymb(IBcVel,SYSTRS,MSTRS,MAXSYM)
elseif (cmd(3)(1:2).eq.'fo') then
IBcFor = n
if (.not.readsy) call rdsymb(IBcFor,SYFORC,MFORC,MAXSYM)
else
write(stderr,150)
  150        format(1x,'set bcsy: unknown type')
ierr = SYNTAX
endif
goto 1000
 endif

 if (cmd(2)(1:2).eq.'te') then
  155       format(1x,'Usage: se[t] te[mperature]',1x, &
     'i[ncrement]|a[bsolute]|r[elative] [itime]')
if (ncmd.lt.3) then
write(stderr,155)
ierr = SYNTAX
goto 1000
endif
if (cmd(3)(1:1).eq.'a') then
itemp = 0
else if (cmd(3)(1:1).eq.'i') then
itemp = 1
else if (cmd(3)(1:1).eq.'r') then
if (ncmd.gt.3) then
    if (.not. numeric(cmd(4)) ) then
    write(stderr,155)
    ierr = SYNTAX
    goto 1000
    endif
    irold = irelt
    irelt = string2int(cmd(4))

    call chktim (irelt, Mech,Temp)

    if ((.not.Mech) .and. (.not.Temp)) then
    write(stderr,156) irelt
  156                   format(1x,'set temperature relative:',1x, &
                 'no output at itime = ',I4)
    ierr = EXEC
    irelt = irold
    goto 1000
    endif
    if (irelt.ne.irold) write(stderr,157) irelt
  157               format(1x,'Warning: (re)read FE output for',1x, &
             'itime=',I4,' for this option to become effective')
endif
itemp = 2
else
write(stderr,155)
ierr = SYNTAX
endif
goto 1000
 endif

!----------------------------------------------------------
 if (cmd(2)(1:2).eq.'ex') then
  160       format(1x,'Usage: se[t] ex[pand] x|y SideNumber',1x, &
     'ElementNumbers')
if (ncmd.lt.5) then
write(stderr,160)
ierr = SYNTAX
goto 1000
endif
if (cmd(3).eq.'x') then
nxpnd = 0
call thenum(cmd,5,ncmd,ixpnd,nxpnd,2,MAXPND,NUMEL,ierr)
if (ierr.ne.OK) nxpnd=0
else if (cmd(3).eq.'y') then
nypnd = 0
call thenum(cmd,5,ncmd,iypnd,nypnd,2,MAXPND,NUMEL,ierr)
if (ierr.ne.OK) nypnd=0
else
write(stderr,160)
ierr = SYNTAX
goto 1000
endif
if (numeric(cmd(4))) then
iside = string2int(cmd(4))
if (iside.lt.1 .or. iside.gt.NEN) then
    write(stderr,165)
  165            format(1x,'set expand: side number out of range')
    ierr = SYNTAX
    goto 1000
endif
if (cmd(3).eq.'x'.and.nxpnd.gt.0) then
    do i=1,nxpnd
    ixpnd(1,i) = iside
    enddo
endif
if (cmd(3).eq.'y'.and.nypnd.gt.0) then
    do i=1,nypnd
      iypnd(1,i) = iside
    enddo
endif
else
write(stderr,160)
ierr = SYNTAX
goto 1000
endif
goto 1000
 endif

!----------------------------------------------------------
 if (cmd(2)(1:3).eq.'den') then
  170        format(1x,'Usage: se[t] den[sity] number')
if (ncmd.lt.3) then
write(stderr,170)
ierr = SYNTAX
goto 1000
endif
if (numeric(cmd(3))) then
RHOA = chreal(cmd(3))
else
write(stderr,170)
ierr = SYNTAX
goto 1000
endif
goto 1000
 endif

!----------------------------------------------------------
 if (cmd(2)(1:2).eq.'sh') then
  180        format(1x,'Usage: se[t] sh[rink] [factor]')
shrink = .true.
if (ncmd.eq.2) goto 1000
if (numeric(cmd(3))) then
shrinkFactor = chreal(cmd(3))
else
write(stderr,180)
ierr = SYNTAX
endif
goto 1000
 else if (cmd(2)(1:4).eq.'nosh') then
shrink = .false.
goto 1000
 endif

!----------------------------------------------------------
 if (cmd(2).eq.'raw') then
icompens = 0
goto 1000
 else if (cmd(2)(1:3).eq.'fre') then
  190        format(1x,'Usage: se[t] fre|bou [Yref]')
icompens = 1
if (ncmd.eq.2) goto 1000
if (numeric(cmd(3))) then
Ygref = chreal(cmd(3))
else
write(stderr,190)
ierr = SYNTAX
endif
goto 1000
 elseif (cmd(2)(1:3).eq.'bou') then
icompens = 2
if (ncmd.eq.2) goto 1000
if (numeric(cmd(3))) then
Ygref = chreal(cmd(3))
else
write(stderr,190)
ierr = SYNTAX
endif
goto 1000
 endif

!----------------------------------------------------------
 if (cmd(2)(1:2).eq.'sc') then
sysc = .true.
goto 1000
 elseif (cmd(2)(1:4).eq.'nosc') then
sysc = .false.
goto 1000
 endif

!----------------------------------------------------------
 if (cmd(2)(1:2).eq.'ap') then
append = .true.
goto 1000
 elseif (cmd(2)(1:4).eq.'noap') then
append = .false.
goto 1000
 endif

!----------------------------------------------------------
 if (cmd(2).eq.'cst') then
CST = .true.
goto 1000
 elseif (cmd(2).eq.'nocst') then
CST = .false.
!           generate element ij-arrays
    call elmsys (ierr)
setria = .true.
goto 1000
 endif

!----------------------------------------------------------
 if (cmd(2)(1:2).eq.'tr') then
triang = .true.
goto 1000
 elseif (cmd(2)(1:4).eq.'notr') then
triang = .false.
goto 1000
 endif

!----------------------------------------------------------
 if (cmd(2)(1:3).eq.'ve ' .or. cmd(2)(1:3).eq.'ver') then
verbose = .true.
goto 1000
 elseif (cmd(2)(1:4).eq.'nove') then
verbose = .false.
goto 1000
 endif

!----------------------------------------------------------
if (cmd(2).eq.'opn') then
    opn = .true.
if (ncmd.eq.2) goto 1000
  210        format(1x,'Usage: se[t] opn [flag [NUMZ [scale [axnu]]]]')
if (cmd(3).eq.'x') then
iopnflag = 0
else if (cmd(3).eq.'d') then
iopnflag = 1
else
write(stderr,211)
  211        format(1x,'set opn: flag should be one of "x","d"')
ierr = SYNTAX
goto 1000
endif
if (ncmd.eq.3) goto 1000
if (.not.numeric(cmd(4))) then
write(stderr,210)
ierr = SYNTAX
goto 1000
endif
NUMZ = nint(chreal(cmd(4)))
if (NUMZ.lt.2 .or. NUMZ.gt.100) then
write(stderr,212)
  212        format(1x,'set opn: NUMNZ should be in range [2,100]')
ierr = SYNTAX
goto 1000
endif
if (ncmd.eq.4) goto 1000
if (.not.numeric(cmd(5))) then
write(stderr,210)
ierr = SYNTAX
goto 1000
endif
opnscale = chreal(cmd(5))
if (ncmd.eq.4) goto 1000
if (.not.numeric(cmd(5))) then
write(stderr,210)
ierr = SYNTAX
goto 1000
endif
opnaxnu = chreal(cmd(5))
goto 1000
elseif (cmd(2).eq.'noopn') then
    opn = .false.
    goto 1000
endif
!
!----------------------------------------------------------
if (cmd(2)(1:3).eq.'T0c') then
    T0calc = .true.
    goto 1000
else if (cmd(2)(1:5).eq.'noT0c') then
    T0calc = .false.
    goto 1000
endif
!
!----------------------------------------------------------
 if (cmd(2).eq.'sign') then
effsign = .true.
goto 1000
 else if (cmd(2).eq.'nosign') then
effsign = .false.
goto 1000
 endif
!
!----------------------------------------------------------
 if (cmd(2)(1:4).eq.'inco') then
strinc = .true.
goto 1000
 else if (cmd(2)(1:4).eq.'comp') then
strinc = .false.
goto 1000
 endif
!
!----------------------------------------------------------
 if (cmd(2)(1:2).eq.'pr') then
project = .true.
goto 1000
 else if (cmd(2)(1:4).eq.'nopr') then
project = .false.
goto 1000
 endif


!----------------------------------------------------------

 ierr = SYNTAX
 write(stderr,995) '"',cmd(2)(1:lnblk(cmd(2))),'"'
  995    format(1x,'set: unknown option ',80A)
 write(stderr,996)
  996    format(1x,'type "? set" for a list of options')

1000    return
 end
