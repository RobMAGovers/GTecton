subroutine setcmd (cmd,ncmd,ierr)

use meshdatamodule,  only: meshdatactx
use conversions,     only: string2int
use vtkmodule,       only: plotthis, outputTags, nOutputTags
use timestepmodule,  only: chktim, faultvelo
use debugmodule,     only: iecho,itest
use errors,          only: exec, syntax, ok
use filemodule,      only: do_soelim
use modeldefinition, only: numfnglobal, NUMSLPglobal, NUMSLP, numel
use modeltopology,   only: nen
use modelctx,        only: getsize
use iomodule,        only: ifein, ifeout, verbose, ascout, stderr
use postProcessData, only: deffac,setria,triang,deformed,CST,opn
use plotControl,     only: coordtype, detdilatation, detgrad, &
                           detgrav, detpot, icompens, ielbas, &
                           ielstp, inpbas, inpstp, iopnflag, &
                           irelt, itemp, numz, opgravchange, &
                           opnaxnu, opnscale, readsy, sysc, &
                           T0calc, vscale, Ygref, &
                           Sunit,Eunit,Funit,Runit,Dunit, &
                           shrink, shrinkFactor, &
                           ixpnd, iypnd, maxpnd, nxpnd, nypnd
use symbols !        only: everything

implicit none
!    implicit none
!-pass
integer ncmd,ierr
character(len=*) cmd
dimension cmd(*)

#include "lenchar.i"


!-locl
logical           :: numeric,Mech,Temp, append
integer, external :: lnblk
integer :: i, ireq, irold, iside, itmp
integer :: l4, l7
integer :: m, n
double precision, external :: chreal
real           :: rhoA, tmp
real           :: u1, u2, x

integer        :: newTag

! integer, parameter ::  MAXPND=80
!-init

!----------------------------------------------------------
!#    se[t] parameter [value [value]]
 1    format(1x,'Usage: se[t] parameter [value [value]]')
 if (ncmd.eq.1) then
! set what?
write(stderr,1)
ierr = SYNTAX
goto 1000
 endif


!----------------------------------------------------------
!#    se[t] coordinate type for gravity calculations
! spherical (in case model is spherical with surface at Re)
! cartesian (in case model has surface at z=0)
 if (cmd(2)(1:10).eq.'coordinate') then
 300     format(1x,'Usage: se[t] coordinate type([cart]esian or [sphe]rical)')
  if (ncmd.lt.3) then
     write(stderr,300)
     ierr = SYNTAX
     goto 1000
  endif
  if (cmd(3)(1:4).eq.'sphe') then
   coordtype = 'sphe'
  elseif (cmd(3)(1:4).eq.'cart') then
   coordtype = 'cart'
  else
   write(stderr,300)
   ierr = SYNTAX
   goto 1000
  endif
  goto 1000
 endif

!----------------------------------------------------------
!#    se[t] gravity field functional (potential, gravity, gravity gradient) to
! be calculated
 if (cmd(2)(1:4).eq.'grav') then
 301     format(1x,'Usage: se[t] [grav]functional type([pote]ntial or', &
        ' [grav]ity or [grad]ient)')
   if (ncmd.lt.4) then
     write(stderr,301)
     ierr = SYNTAX
     goto 1000
   endif
   if (.not.numeric(cmd(4))) then
     ierr = SYNTAX
     write(stderr,301)
     goto 1000
   endif
   if (cmd(3)(1:4).eq.'pote') then
       if (cmd(4).eq.'1') then    
         detpot = .true.
       elseif (cmd(4).eq.'0') then
         detpot = .false.
       else
         ierr = SYNTAX
         write(stderr,301)
         goto 1000
       endif
   elseif (cmd(3)(1:4).eq.'grav') then
       if (cmd(4).eq.'1') then             
         detgrav = .true.
       elseif (cmd(4).eq.'0') then
         detgrav = .false.
       else
         ierr = SYNTAX
         write(stderr,301)
         goto 1000      
       endif
   elseif (cmd(3)(1:4).eq.'grad') then
       if (cmd(4).eq.'1') then
         detgrad = .true.
       elseif (cmd(4).eq.'0') then
         detgrad = .false.
       else                     
         ierr = SYNTAX                           
         write(stderr,301)                                
         goto 1000                                                 
      endif   
   elseif (cmd(3)(1:4).eq.'dila') then
      if (cmd(4).eq.'1') then
         detdilatation = .true.
      elseif (cmd(4).eq.'0') then
         detdilatation = .false.
      else
         ierr = SYNTAX
          write(stderr,301)
         goto 1000
      endif
    elseif (cmd(3)(1:4).eq.'chan') then
      if (cmd(4).eq.'1') then
        OpGravChange = .true.
      elseif (cmd(4).eq.'0') then
        OpGravChange = .false.
      else
        ierr = SYNTAX
        write(stderr,301)
        goto 1000
       endif
   else
      ierr = SYNTAX                           
      write(stderr,301)                                
      goto 1000
   endif
!   write(stderr,*) 'detpot,detgrav,detgrad'
!   write(stderr,*) detpot,detgrav,detgrad
   goto 1000
  endif





!----------------------------------------------------------
!#    se[t] sy[mbolsize] size
!  if (cmd(2)(1:2).eq.'sy') then
!  5        format(1x,'Usage: se[t] sy[mbolsize] size')
! if (ncmd.lt.3) then
!     write(stderr,5)
!     ierr = SYNTAX
!     goto 1000
! endif
! if (.not.numeric(cmd(3))) then
!     write(stderr,5)
!     ierr = SYNTAX
!     goto 1000
! endif
! tmp = chreal(cmd(3))
! if (tmp.le.0.0) then
!     write(stderr,10)
! 10        format(1x,'positive symbol size required')
!     ierr = SYNTAX
!     goto 1000
! endif
! sysi=tmp
! goto 1000
!  endif


!----------------------------------------------------------
!#    se[t] stre[ss] value1 unit1 = value2 unit2
!#    se[t] stra[in] value1 unit1 = value2 unit2
!#    se[t] ra[te] value1 unit1 = value2 unit2
!#    se[t] di[splacement] value1 unit1 = value2 unit2
!#    se[t] fl[ux] value1 unit1 = value2 unit2
! 
!  n = 0
! 
!  if (cmd(2)(1:4).eq.'stre') then
! n = 1
!  else if (cmd(2)(1:4).eq.'stra') then
! n = 2
!  else if (cmd(2).eq.'ra' .or. cmd(2)(1:3).eq.'rat') then
! n = 3
!  else if (cmd(2)(1:2).eq.'di') then
! n = 4
!  else if (cmd(2)(1:2).eq.'fl') then
! n = 5
!  endif
! 
!  if (n.gt.0) then
! 15       format(1x,'Usage: se[t] stre[ss] value1 unit1 = value2',1x, &
!      'unit2')
! 20       format(1x,'Usage: se[t] stra[ss] value1 unit1 = value2',1x, &
!      'unit2')
! 25       format(1x,'Usage: se[t] ra[te] value1 unit1 = value2',1x, &
!      'unit2')
! 30       format(1x,'Usage: se[t] di[sp] value1 unit1 = value2',1x, &
!      'unit2')
! 32       format(1x,'Usage: se[t] fl[ux] value1 unit1 = value2',1x, &
!      'unit2')
! if (ncmd.lt.7) then
!     if (n.eq.1) write(stderr,15)
!     if (n.eq.2) write(stderr,20)
!     if (n.eq.3) write(stderr,25)
!     if (n.eq.4) write(stderr,30)
!     if (n.eq.5) write(stderr,32)
!     ierr = SYNTAX
!     goto 1000
! endif
! 
! if (.not. numeric(cmd(3)) ) then
!     if (n.eq.1) write(stderr,15)
!     if (n.eq.2) write(stderr,20)
!     if (n.eq.3) write(stderr,25)
!     if (n.eq.4) write(stderr,30)
!     if (n.eq.5) write(stderr,32)
!     ierr = SYNTAX
!     goto 1000
! endif
! 
! u1 = chreal(cmd(3))
! 
! if (u1.lt.1D-20) then
!     if (n.eq.1) write(stderr,15)
!     if (n.eq.2) write(stderr,20)
!     if (n.eq.3) write(stderr,25)
!     if (n.eq.4) write(stderr,30)
!     if (n.eq.5) write(stderr,32)
!     ierr = SYNTAX
!     goto 1000
! endif
! 
! if (.not. numeric(cmd(6)) ) then
!     if (n.eq.1) write(stderr,15)
!     if (n.eq.2) write(stderr,20)
!     if (n.eq.3) write(stderr,25)
!     if (n.eq.4) write(stderr,30)
!     if (n.eq.5) write(stderr,32)
!     ierr = SYNTAX
!     goto 1000
! endif
! 
! u2 = chreal(cmd(6))
! if (u2.lt.1D-20) then
!     if (n.eq.1) write(stderr,15)
!     if (n.eq.2) write(stderr,20)
!     if (n.eq.3) write(stderr,25)
!     if (n.eq.4) write(stderr,30)
!     if (n.eq.5) write(stderr,32)
!     ierr = SYNTAX
!     goto 1000
! endif
! 
! if (n.eq.1) then
!     Sscale = u2/u1
!     Ssb    = u1
!     Sunit  = cmd(4)
! else if (n.eq.2) then
!     Escale = u2/u1
!     Esb    = u1
!     Eunit  = cmd(4)
! else if (n.eq.3) then
!     Rscale = u2/u1
!     Rsb    = u1
!     Runit  = cmd(4)
! else if (n.eq.4) then
!     Dscale = u2/u1
!     Dsb    = u1
!     Dunit  = cmd(4)
! else if (n.eq.5) then
!     Fscale = u2/u1
!     Fsb    = u1
!     Funit  = cmd(4)
! endif
! 
! if (iecho.eq.1) then
!     l4 = lnblk(cmd(4))
!     l7 = lnblk(cmd(7))
! 
!     if (n.eq.1) then
!         write(stderr,35) Sscale,cmd(7)(1:l7),'/', &
!            cmd(4)(1:l4)
! 35            format(1x,'stress scale factor = ',1PG12.6,1X,80a)
!         write(stderr,40) Ssb,cmd(4)(1:l4)
! 40            format(1x,'stress scale bar    = ',1PG12.6,1X,80a)
!     else if (n.eq.2) then
!         write(stderr,45) Escale,cmd(7)(1:l7),'/', &
!            cmd(4)(1:l4)
! 45            format(1x,'strain scale factor = ',1PG12.6,1X,80a)
!         write(stderr,50) Esb,cmd(4)(1:l4)
! 50            format(1x,'strain scale bar    = ',1PG12.6,1X,80a)
!     else if (n.eq.3) then
!         write(stderr,55) Rscale,cmd(7)(1:l7),'/', &
!            cmd(4)(1:l4)
! 55            format(1x,'rate scale factor = ',1PG12.6,1X,80a)
!         write(stderr,60) Rsb,cmd(4)(1:l4)
! 60            format(1x,'rate scale bar    = ',1PG12.6,1X,80a)
!     else if (n.eq.4) then
!         write(stderr,65) Dscale,cmd(7)(1:l7),'/', &
!            cmd(4)(1:l4)
! 65            format(1x,'displ. scale factor = ',1PG12.6,1X,80a)
!         write(stderr,70) Dsb,cmd(4)(1:l4)
! 70            format(1x,'displ. scale bar    = ',1PG12.6,1X,80a)
!     else if (n.eq.5) then
!         write(stderr,72) Fscale,cmd(7)(1:l7),'/', &
!            cmd(4)(1:l4)
! 72            format(1x,'flux scale factor = ',1PG12.6,1X,80a)
!         write(stderr,74) Fsb,cmd(4)(1:l4)
! 74            format(1x,'flux scale bar    = ',1PG12.6,1X,80a)
! 
!     endif
! endif
! goto 1000
!  endif

!----------------------------------------------------------
!#    se[t] number of piece files in the VTK output.
! TODO; this is not yet fully implemented.
 if (cmd(2)(1:6).eq.'nprocs') then

if (ncmd.gt.2) then
    ierr = SYNTAX
    write(stderr,*) 'Usage: set nprocs number'
    goto 1000
endif

!        nprocs = cmd(2)

write(*,*) "Understood, will write", getsize(),"piece files."

 endif


!----------------------------------------------------------
!#    se[t] whether to do soelim or not
 if (cmd(2)(1:4).eq.'noso') then
    ! noso for no soelim
    if (ncmd.gt.2) then
        ierr = SYNTAX
        write(stderr,*) 'Usage: set nosoelim'
        goto 1000
    endif

    do_soelim = .false.
    goto 1000
 endif
 if (cmd(2)(1:4).eq.'doso') then
    ! doso for do soelim

    if (ncmd.gt.2) then
        ierr = SYNTAX  
        write(stderr,*) 'Usage: set dosoelim'
        goto 1000
    endif
    do_soelim = .true.
    goto 1000
 endif

!----------------------------------------------------------
!#    se[t] ar[rowhead] size angle [no[co[nt]]]
!  if (cmd(2)(1:2).eq.'ar') then
! 75       format(1x,'Usage: se[t] ar[rowhead] size angle',1x, &
!      '[no[co[nt]]]')
! 
! if (ncmd.lt.4) then
!     ierr = SYNTAX
!     write(stderr,75)
!     goto 1000
! endif
! 
! if (.not.numeric(cmd(3))) then
!     ierr = SYNTAX
!     write(stderr,75)
!     goto 1000
! endif
! 
! tmp = chreal(cmd(3))
! if (tmp.le.0.E0) then
!     ierr = SYNTAX
!     write(stderr,75)
!    goto 1000
! endif
! 
! hatsize = tmp
! if (.not.numeric(cmd(4))) then
!     ierr = SYNTAX
!     write(stderr,75)
!     goto 1000
! endif
! 
! tmp = chreal(cmd(4))
! if (tmp.le.0.E0 .or. tmp.ge.9.E1) then
!     ierr = SYNTAX
!     write(stderr,80)
! 80        format(1x,'0 < angle < 90 degrees')
!     goto 1000
! endif
! 
! hatangle = tmp
! if (ncmd.gt.4) then
!     if (cmd(5)(1:2).eq.'co') then
!         iarcont = 1
!     else if (cmd(5)(1:4).eq.'noco') then
!         iarcont = 0
!     else
!         write(stderr,75)
!         ierr = SYNTAX
!         goto 1000
!     endif
! endif
! goto 1000
!  endif

!----------------------------------------------------------
!# se[t] [no]ec[ho] [number]
if ( cmd(2)(1:2).eq.'ec' ) then

!    write(*,*) "1setting echo: ", cmd(3)
!    write(*,*) "2setting echo: ", trim(cmd(3))
!    write(*,*) "3setting echo: ", chreal(cmd(3))
!    write(*,*) "4setting echo: ", int(chreal(cmd(3)))


85        format(1x,'Usage: se[t] [no]ec[ho] [number]')
    if ( ncmd.gt.2 ) then
        if (.not. numeric(cmd(3)) ) then
            ierr = SYNTAX
            write(stderr,85)
            goto 1000
        endif
!        itmp = INT( chreal(cmd(3)) )
        read(cmd(3),*) itmp
        iecho = itmp
!        iecho=3 is magic to tecin.f
        if (iecho.eq.3) then
            iecho = 4
        endif
    else
        iecho = 1
    endif
    write(*,*) "set echo to", iecho
    goto 1000
elseif ( cmd(2)(1:4).eq.'noec' ) then
    iecho = 0
    goto 1000
endif

!----------------------------------------------------------
!# se[t] tag


if ( cmd(2)(1:3).eq.'tag' ) then
    if (ncmd.lt.3) then
        ierr = SYNTAX
        write(stderr,*) "Usage: se[t] tag number"
        goto 1000
    endif
    
    nOutputTags = nOutputTags + 1

    newTag = INT(chreal(cmd(3)))

    if (newTag.eq.0) then
        ! reset tags
        nOutputTags = 0
        outputTags = 0    
    else
        ! append a new tag to the list
        if (numeric(cmd(3))) then
            outputTags(nOutputTags) = INT(chreal(cmd(3)))
        else
            write(*,*) "Third argument of 'set tag xxx' must be numeric"
        endif
    endif
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

    IREQ = INT(chreal(cmd(3)))

!            write(*,*) 'found IREQ: ', IREQ

    call chktim (IREQ,Mech, Temp)

!            write(*,*) 'with mech: ', mech

    if ((.not.Mech) .and. (.not.Temp)) then
        write(stderr,95) IREQ
95            format(1x,'set itime: no output at itime = ',I6)
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

    if (x.lt.0. .or. x.gt.1e30) then
        write(stderr,98)
98               format(1x,'set deformed: amplification out of',1x, &
             'range')
        ierr = SYNTAX
        goto 1000
    endif

    deffac = x
else
    deffac = 1.d0
endif

goto 1000
 elseif (cmd(2)(1:4).eq.'node') then
deformed = .false.
deffac = 0.d0
goto 1000
 endif

!----------------------------------------------------------
!# se[t] [no]as[cii]
 if (cmd(2)(1:2).eq.'as') then
ascout = .true.
goto 1000
 elseif (cmd(2)(1:4).eq.'noas') then
ascout = .false.
goto 1000
 endif

!----------------------------------------------------------
!# se[t] [no]bi[nary]
 if (cmd(2)(1:2).eq.'bi') then
ascout = .false.
goto 1000
 elseif (cmd(2)(1:4).eq.'nobi') then
ascout = .true.
goto 1000
 endif

!---------------------------------------------------------------
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
    m = INT(chreal(cmd(4)))
else
    m = n
    n = ielbas
endif

if (n.lt.1 .or. n.gt.meshdatactx%neglobal) then
    write(stderr,105)
  105        format(1x,'set elmi: base number out of range')
    ierr = EXEC
    goto 1000
endif

if (m.lt.1 .or. m.gt.meshdatactx%neglobal) then
    write(stderr,110)
  110        format(1x,'set elmi: increment out of range')
    ierr = EXEC
    goto 1000
endif

ielbas = n
ielstp = m
goto 1000
 endif

!--------------------------------------------------------------------
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
    if (.not.numeric(cmd(4))) then
        write(stderr,115)
        ierr = SYNTAX
        goto 1000
    endif
    m = INT(chreal(cmd(4)))
else
    m = n
    n = ielbas
endif


if (n.lt.1 .or. n.gt.meshdatactx%nvglobal) then
    write(stderr,120)
  120        format(1x,'set npi: base number out of range')

    ierr = EXEC
    goto 1000
endif

if (m.lt.1 .or. m.gt.meshdatactx%nvglobal) then
    write(stderr,125)
  125        format(1x,'set npi: increment out of range')
    ierr = EXEC
    goto 1000
endif

inpbas = n
inpstp = m
goto 1000
 endif

!-------------------------------------------------------------
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

!-------------------------------------------------------------
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
n = INT(chreal(cmd(4)))
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

!-------------------------------------------------------------
 if (cmd(2)(1:3).eq.'tem') then
  155       format(1x,'Usage: se[t] tem[perature]',1x, &
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
    irelt = INT(chreal(cmd(4)))

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

!-------------------------------------------------------------
if (cmd(2).eq.'test') then
    itest = 1
    goto 1000
endif
if (cmd(2).eq.'notest') then
    itest = 0
    goto 1000
endif
!-------------------------------------------------------------
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
        if (ierr.ne.OK) then
            nxpnd=0
        endif
    else if (cmd(3).eq.'y') then
        nypnd = 0
        call thenum(cmd,5,ncmd,iypnd,nypnd,2,MAXPND,NUMEL,ierr)
        if (ierr.ne.OK) then
            nypnd=0
        endif
    else
        write(stderr,160)
        ierr = SYNTAX
        goto 1000
    endif

    if (numeric(cmd(4))) then
        iside = INT(chreal(cmd(4)))
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
shrinkfactor = chreal(cmd(3))
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
 if (cmd(2)(1:2).eq.'tr') then
triang = .true.
goto 1000
 elseif (cmd(2)(1:4).eq.'notr') then
triang = .false.
goto 1000
 endif

!----------------------------------------------------------
 if (cmd(2)(1:2).eq.'ve') then
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
  210       format(1x,'Usage: se[t] opn [flag [NUMZ [scale [axnu]]]]')
    if (cmd(3).eq.'x') then
        iopnflag = 0
    else if (cmd(3).eq.'d') then
        iopnflag = 1
    else
        write(stderr,211)
  211           format(1x,'set opn: flag should be one of "x","d"')
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
  212           format(1x,'set opn: NUMNZ should be in range [2,100]')
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
if (cmd(2)(1:3).eq.'vtu') then

if (ncmd.eq.2) then
    write(*,*) 'set what?'
   	call writeAllOptions()


    if (plotthis%precision.eq.1) then
        write(*,*) 'precision;       currently set to single precision'
    else if (plotthis%precision.eq.2) then
        write(*,*) 'precision;       currently set to double precision'
    endif
    write(*,*)
else if (ncmd.eq.3) then
    write(*,*) 'supply whether to set it ''on'' or ''off''.'
    write(*,*) ' ''1'' and ''0'' also accepted for ''on'' and ''off''.'
else if (ncmd.eq.4) then
    !------- prestress ----------
    if (cmd(3).eq.'prestress') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%prestress = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%prestress = .false.
        else
         write(*,*) 'set prestress with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- materials ----------
    else if (cmd(3).eq.'materials') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%materials = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%materials = .false.
        else
         write(*,*) 'set materials with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- velocity BC ----------
    else if (cmd(3).eq.'velocityBC') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
            plotthis%velocityBC = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
            plotthis%velocityBC = .false.
        else
            write(*,*) 'set velocityBC with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- Euler Angles BC ----------
    else if (cmd(3).eq.'Euler' .or. cmd(3).eq.'euler') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
            plotthis%EulerAnglesBC = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
            plotthis%EulerAnglesBC = .false.
        else
            write(*,*) 'set euler with ''0'',''1'',''on'' or ''off'' '
        endif

    !------- velocity ----------
    else if (cmd(3).eq.'velocity') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%velocity = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%velocity = .false.
        else
            write(*,*) 'set velocity with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- slippery markers ----------
    else if (cmd(3).eq.'slipperymarkers') then
        write(*,*) "Warning, slipperymarkers not yet implemented", NUMSLP, NUMSLPglobal
        if (NUMSLPglobal.eq.0) then
            write(*,*) "Warning: Trying to change slippery markers setting, but the model has no slippery nodes"
            write(*,*) "Leaving slipperymarkers .false."
        else
            if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
                   plotthis%slipperymarkers = .true.
            else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
                plotthis%slipperymarkers = .false.
            else
                write(*,*) 'set slipperymarkers with ''0'',''1'',''on'' or ''off'' '
            endif
        endif
    !------- slippery vectors ----------
    else if (cmd(3).eq.'slipperyvectors') then
        if (NUMSLPglobal.eq.0) then
            write(*,*) "Warning: Trying to change slippery vectors setting, but the model has no slippery nodes"
            write(*,*) "Leaving slipperyvectors .false."
        else
            if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
                 plotthis%slipperyvectors = .true.
            else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
                 plotthis%slipperyvectors = .false.
            else
             write(*,*) 'set slipperyvectors with ''0'',''1'',''on'' or ''off'' '
            endif
        endif
    !------- slippery displacement ----------
    else if (cmd(3).eq.'slipperydisplacement' .or. &
             cmd(3).eq.'slipperyDisplacement') then
        if (NUMSLPglobal.eq.0) then
            write(*,*) "Warning: Trying to change slippery displacement setting, but the model has no slippery nodes."
            write(*,*) "Leaving slipperydisplacement .false."
        else
            if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
                 plotthis%slipperydisplacement = .true.
            else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
                 plotthis%slipperydisplacement = .false.
            else
             write(*,*) 'set slipperydisplacement with ''0'',''1'',''on'' or ''off'' '
            endif
        endif
    !------- slippery velocity ----------
    else if (cmd(3).eq.'slipperyvelocity' .or. &
             cmd(3).eq.'slipperyVelocity') then
        if (NUMSLPglobal.eq.0) then
            write(*,*) "Warning: Trying to change slippery velocity setting, but the model has no slippery nodes"
            write(*,*) "Leaving slipperyvelocity .false."
        else
            if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
                 plotthis%slipperyVelocity = .true.
            else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
                 plotthis%slipperyVelocity = .false.
            else
             write(*,*) 'set slipperyvelocity with ''0'',''1'',''on'' or ''off'' '
            endif
        endif
    !------- slippery weights ----------
    else if (cmd(3).eq.'slipperyweights' .or. &
             cmd(3).eq.'slipperyWeights') then
        if (NUMSLPglobal.eq.0) then
            write(*,*) "Warning: Trying to change slippery weights setting, but the model has no slippery nodes"
            write(*,*) "Leaving slipperyweights .false."
        else
            if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
                 plotthis%slipperyweights = .true.
            else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
                 plotthis%slipperyweights = .false.
            else
                 write(*,*) 'set slipperyweights with ''0'',''1'',''on'' or ''off'' '
            endif
        endif
    !------- faulted displacement ----------
    else if (cmd(3).eq.'faulteddisplacement' .or. &
             cmd(3).eq.'faultedDisplacement') then

        if (NUMFNglobal.eq.0) then
            write(*,*) "Warning: Trying to change faulted displacement setting, but the model has no faulted nodes"
            write(*,*) "Leaving faultedDisplacement .false."
        else
            if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
                 plotthis%faultedDisplacement = .true.
            else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
                 plotthis%faultedDisplacement = .false.
            else
                 write(*,*) 'set faultedDisplacement with ''0'',''1'',''on'' or ''off'' '
            endif
        endif
    !------- faulted velocity ----------
    else if (cmd(3).eq.'faultedvelocity' .or. &
             cmd(3).eq.'faultedVelocity') then

        if (NUMFNglobal.eq.0) then
            write(*,*) "Warning: Trying to change faulted velocity setting, but the model has no faults."
            write(*,*) "Leaving faultedVelocity .false."
        else
            if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
                 plotthis%faultedVelocity = .true.
            else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
                 plotthis%faultedVelocity = .false.
            else
                 write(*,*) 'set faultedvelocity with ''0'',''1'',''on'' or ''off'' '
            endif
        endif
    !------- displacement ----------
    else if (cmd(3).eq.'displacement') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%displacement = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%displacement = .false.
        else
         write(*,*) 'set displacement with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- displacementBC ----------
    else if (cmd(3).eq.'displacementBC') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%displacementBC = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%displacementBC = .false.
        else
         write(*,*) 'set displacementBC with ''0'',''1'',''on'' or ''off'' '
        endif

    !------- slideInPlaneBC ----------
    else if (cmd(3).eq.'slideInPlaneBC') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%slideInPlaneBC = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%slideInPlaneBC = .false.
        else
         write(*,*) 'set slideInPlaneBC with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- slideOnLineBC ----------
    else if (cmd(3).eq.'slideOnLineBC') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%slideOnLineBC = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%slideOnLineBC = .false.
        else
         write(*,*) 'set slideOnLineBC with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- fixedPointBC ----------
    else if (cmd(3).eq.'fixedPointBC') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%fixedPointBC = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%fixedPointBC = .false.
        else
         write(*,*) 'set fixedPointBC with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- stress ----------
    else if (cmd(3).eq.'stress') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%stress = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%stress = .false.
        else
         write(*,*) 'set stress with ''0'',''1'',''on'' or ''off'' '
        endif
#ifdef EXPERIMENTALTRANSIENTDISCREEP
    !------- Taylor stress ----------
    else if (cmd(3).eq.'taylorstress') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%taylorstress = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%taylorstress = .false.
        else
        write(*,*) 'set taylor stress with ''0'',''1'',''on'' or ''off'' '
    endif
#endif      
    !------- viscosity ----------
    else if (cmd(3).eq.'viscosity') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%viscosity = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%viscosity = .false.
        else
         write(*,*) 'set viscosity with ''0'',''1'',''on'' or ''off'' '
    endif
    !------- strain ----------
    else if (cmd(3).eq.'strain') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%strain = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%strain = .false.
        else
         write(*,*) 'set strain with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- strain rate ----------
    else if (cmd(3).eq.'strainrate') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%strainrate = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%strainrate = .false.
        else 
         write(*,*) 'set strainrate with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- element volume ----------
    else if (cmd(3).eq.'volume') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%volume = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%volume = .false.
        else
         write(*,*) 'set volume with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- elempartition ----------
    else if (cmd(3).eq.'elempartition') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%elempartition = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%elempartition = .false.
        else
         write(*,*) 'set elempartition with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- vertexpartition ----------
    else if (cmd(3).eq.'vertexpartition') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%vertexpartition = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%vertexpartition = .false.
        else
         write(*,*) 'set vertexpartition with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- vertexpartition ----------
    else if (cmd(3).eq.'vertexmarker') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%vertexmarker = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%vertexmarker = .false.
        else
         write(*,*) 'set vertexmarker with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- ascii/binary writing ----------
    else if (cmd(3).eq.'writeASCII') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%writeASCII = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%writeASCII = .false.
        else
         write(*,*) 'set writeASCII with ''0'',''1'',''on'' or ''off'' '
        endif

    !------- Three element qualities. 1: sides  ----------
    else if (cmd(3).eq.'qualitySide') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%elemQualitySide = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%elemQualitySide = .false.
        else
         write(*,*) 'set qualitySide with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- Three element qualities. 2: surfaces  ----------
    else if (cmd(3).eq.'qualitySurface') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%elemQualitySurface = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%elemQualitySurface = .false.
        else
         write(*,*) 'set qualitySurface with ''0'',''1'',''on'' or ''off'' '
        endif
    !------- Three element qualities. 3: volumes  ----------
    else if (cmd(3).eq.'qualityVolume') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%elemQualityVolume = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%elemQualityVolume = .false.
        else
         write(*,*) 'set qualityVolume with ''0'',''1'',''on'' or ''off'' '
        endif

    !------- precision ----------
    else if (cmd(3).eq.'precision') then
        if (plotthis%writeASCII .eqv. .true.) then
            write(*,*) 'warning: precision only relevant for binary writing'
            write(*,*) 'see: set vtu writeASCII [on/off]'
        endif

        if ((cmd(4).eq.'double').or.(cmd(4).eq.'2')) then
             plotthis%precision = 2
        else if ((cmd(4).eq.'single').or.(cmd(4).eq.'1')) then
             plotthis%precision = 1
        else
         write(*,*) 'set precision with ''1'',''2'',''single'' or ''double'' '
        endif
    !------- error  ----------
    else if (cmd(3).eq.'error') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
             plotthis%error = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
             plotthis%error = .false.
        else
         write(*,*) 'set error plotting with ''0'',''1'',''on'' or ''off'' '
        endif

    else if (cmd(3).eq.'all') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
            plotthis%prestress = .true.
            plotthis%materials = .true.
            plotthis%slipperymarkers = .true.
            plotthis%slipperyvectors = .true.
            plotthis%slipperyweights = .true.
            plotthis%faultedDisplacement = .true.
            plotthis%velocity = .true.
            plotthis%displacement = .true.
            plotthis%stress = .true.
            plotthis%strain = .true.
            plotthis%strainrate = .true.
            plotthis%volume = .true.
            plotthis%elempartition = .true.
            plotthis%vertexpartition = .true.
            plotthis%vertexmarker = .true.
            plotthis%elemQualitySide  = .true.
            plotthis%elemQualitySurface  = .true.
            plotthis%elemQualityVolume = .true.
            plotthis%velocityBC      = .true.
            plotthis%slideInPlaneBC  =.true.
            plotthis%slideOnLineBC   =.true.
            plotthis%fixedPointBC    =.true.
            plotthis%displacementBC  =.true.
            plotthis%forceBC         =.true.
            plotthis%EulerAnglesBC   =.true.
            plotthis%error           =.true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
            plotthis%prestress = .false.
            plotthis%materials = .false.
            plotthis%slipperymarkers = .false.
            plotthis%slipperyvectors = .false.
            plotthis%slipperyweights = .false.
            plotthis%faultedDisplacement = .false.
            plotthis%velocity = .false.
            plotthis%displacement = .false.
            plotthis%stress = .false.
            plotthis%strain = .false.
            plotthis%strainrate = .false.
            plotthis%volume = .false.
            plotthis%elempartition = .false.
            plotthis%vertexpartition = .false.
            plotthis%vertexmarker = .false.
            plotthis%elemQualitySide  = .false.
            plotthis%elemQualitySurface  = .false.
            plotthis%elemQualityVolume = .false.
            plotthis%velocityBC      = .false.
            plotthis%slideInPlaneBC  =.false.
            plotthis%slideOnLineBC   =.false.
            plotthis%fixedPointBC    =.false.
            plotthis%displacementBC  =.false.
            plotthis%forceBC         =.false.
            plotthis%EulerAnglesBC   =.false.
            plotthis%error           =.false.
        endif
    else if (cmd(3).eq.'allBCs') then
        if ((cmd(4).eq.'on').or.(cmd(4).eq.'1')) then
            plotthis%velocityBC      = .true.
            plotthis%slideInPlaneBC  = .true.
            plotthis%slideOnLineBC   = .true.
            plotthis%fixedPointBC    = .true.
            plotthis%displacementBC  = .true.
            plotthis%forceBC         = .true.
            plotthis%EulerAnglesBC   = .true.
        else if ((cmd(4).eq.'off').or.(cmd(4).eq.'0')) then
            plotthis%velocityBC      = .false.
            plotthis%slideInPlaneBC  = .false.
            plotthis%slideOnLineBC   = .false.
            plotthis%fixedPointBC    = .false.
            plotthis%displacementBC  = .false.
            plotthis%forceBC         = .false.
            plotthis%EulerAnglesBC   = .false.
        endif
    else
        !------- we have a woopsie ----------
        write(*,*) 'unknown VTU option: ', trim(cmd(3))
        call writeAllOptions()
    endif
else
    write(*,*) 'Too many or too few arguments. Use four.'
    write(*,*) 'example: set vtu [property] [on/off]'
endif
goto 1000
 endif


!----------------------------------------------------------

 ierr = SYNTAX
 write(stderr,995) '"',cmd(2)(1:lnblk(cmd(2))),'"'
  995    format(1x,'set: unknown option ',80A)
 write(stderr,996)
  996    format(1x,'type "? set" for a list of options')

1000    return
 end subroutine

!----------------------------------------------------------


subroutine writeAllOptions()

use vtkmodule,       only: plotthis

implicit none

    write(*,*) 'choose from: '
    write(*,*)
    write(*,*) 'all;             turn everything on or off'
    write(*,*) 'allBCs;          turn all boundary conditions on or off'
    write(*,*)
    write(*,*) 'Boundary conditions on nodal points:'
    write(*,*)
    write(*,*) 'velocityBC;	 currently set to', plotthis%velocityBC
    write(*,*) 'slideInPlaneBC;  currently set to', plotthis%slideInPlaneBC
    write(*,*) 'slideOnLineBC;   currently set to', plotthis%slideOnLineBC
    write(*,*) 'fixedPointBC;    currently set to', plotthis%fixedPointBC
    write(*,*) 'displacementBC;  currently set to', plotthis%displacementBC
    write(*,*) 'forceBC;         currently set to', plotthis%forceBC
    write(*,*) 'Euler angles;    currently set to', plotthis%EulerAnglesBC
    write(*,*)
    write(*,*) 'Boundary conditions on elements:'
    write(*,*)
    write(*,*) 'prestress;	 currently set to', plotthis%prestress
    write(*,*) 'materials;	 currently set to', plotthis%materials
    write(*,*) 'slipperymarkers; currently set to', plotthis%slipperyMarkers
    write(*,*) 'slipperyvectors; currently set to', plotthis%slipperyVectors
    write(*,*) 'slipperyweights; currently set to', plotthis%slipperyWeights
    write(*,*)
    write(*,*) 'Output data:'
    write(*,*)
    write(*,*) 'displacement;    currently set to', plotthis%displacement
    write(*,*) 'faultedDisplacement;  currently set to', plotthis%faultedDisplacement
    write(*,*) 'velocity;        currently set to', plotthis%velocity
    write(*,*) 'stress;          currently set to', plotthis%stress
#ifdef EXPERIMENTALTRANSIENTDISCREEP
    write(*,*) 'taylor stress;   currently set to', plotthis%taylorstress
#endif
    write(*,*) 'strain;          currently set to', plotthis%strain
    write(*,*) 'strainrate;  	 currently set to', plotthis%strainrate
    write(*,*) 'viscosity;       currently set to', plotthis%viscosity
    write(*,*) 'element volume;  currently set to', plotthis%volume
    write(*,*) 'elempartition;   currently set to', plotthis%elempartition
    write(*,*) 'vertexpartition; currently set to', plotthis%vertexpartition
    write(*,*) 'writeASCII;	     currently set to', plotthis%writeASCII
    write(*,*)
    write(*,*) 'quality of elements, based on property:'
    write(*,*)
    write(*,*) 'qualitySide;     currently set to', plotthis%elemQualitySide
    write(*,*) 'qualitySurface;  currently set to', plotthis%elemQualitySurface
    write(*,*) 'qualityVolume;   currently set to', plotthis%elemQualityVolume
    write(*,*)
    write(*,*) 'quality of solution, based on difference of [strain + rotation] . normalvector'
    write(*,*) 'across element boundaries'
    write(*,*)
    write(*,*) 'error;           currently set to', plotthis%error

    stop "Leaving Plt3d"

end subroutine

