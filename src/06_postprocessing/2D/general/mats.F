subroutine matplt (MAT,PRPMAT,PRPCRP,PRPLAS,PRPTEM, &
                   LUMON,IFLAG,VERBOSE)

USE MODELDEFINITION
use materialsmodule, only: MaterialIDOfElement
use monitorVars,     only: elmon, monel, nnel

! subroutine to write material numbers to output file

implicit none
!-pass
integer LUMON,IFLAG
logical VERBOSE
integer   :: mat
double precision :: prpmat, prpcrp, prptem, prplas
dimension MAT(*),PRPMAT(6,*), &
          PRPCRP(6,*),PRPTEM(6,*),PRPLAS(9,*)

!-locl
#include "lenchar.i"
character(len=LENCHAR) string
real D2REAL
integer           :: m, n, nn, nflow, nptype
integer, external :: lnblk
!
if (.not.elmon) then
    return
endif
!
do nn=1,nnel

    n = nn
    n = monel(nn)
    m = MaterialIDOfElement(n)

    if (iflag.eq.0) then

        if (VERBOSE) then
            write(LUMON,10) n,m
10            format(1x,'element ',I5,': material property ',I5)
        else
            write(LUMON,*) n,m
        endif

    elseif (iflag.eq. 1) then

        if (VERBOSE) then
            write(LUMON,20) n,PRPMAT(1,m)
20               format(1x,'element ',I5,': Young''s modulus ', &
             1PG12.2)
        else
            write(LUMON,*) n,PRPMAT(1,m)
        endif

    elseif (iflag.eq.2) then

        if (VERBOSE) then
            write(LUMON,30) n,PRPMAT(2,m)
30            format(1x,'element ',I5,': Poison ratio ',F5.2)
        else
            write(LUMON,*) n,PRPMAT(2,m)
        endif

    elseif (iflag.eq. 3) then

        if (VERBOSE) then
            write(LUMON,40) n,PRPMAT(3,m)
40            format(1x,'element ',I5,': viscosity ',1PG14.6)
            write(LUMON,*) n,PRPMAT(3,m)
        endif

    elseif (iflag.eq. 4) then

        if (VERBOSE) then
            write(LUMON,50) n,PRPMAT(4,m)
50            format(1x,'element ',I5,': powerlaw index ',F4.2)
        else
            write(LUMON,*) n,PRPMAT(4,m)
        endif

    elseif (iflag.eq. 5) then

        if (VERBOSE) then
            write(LUMON,60) n,PRPMAT(5,m)
60            format(1x,'element ',I5,': density ',1PG12.2)
        else
            write(LUMON,*) n,PRPMAT(5,m)
        endif

    elseif (iflag.eq. 6) then

        if (VERBOSE) then
            write(LUMON,70) n,PRPMAT(6,m)
70            format(1x,'element ',I5,': thickness ',1PG12.2)
        else
            write(LUMON,*) n,PRPMAT(6,m)
        endif

    elseif (iflag.eq.6+1) then

        if (VERBOSE) then
            write(LUMON,80) n,PRPCRP(1,m)
80               format(1x,'element ',I5, &
             ': powerlaw activation energy ',1PG14.6)
        else
            write(LUMON,*) n,PRPCRP(1,m)
        endif

    elseif (iflag.eq.6+2) then

        if (VERBOSE) then
            write(LUMON,90) n,PRPCRP(2,m)
90               format(1x,'element ',I5, &
             ': powerlaw pre-exponential factor ',1PG14.6)
        else
            write(LUMON,*) n,PRPCRP(2,m)
        endif

    elseif (iflag.eq.6+3) then

        if (VERBOSE) then
            write(LUMON,100) n,PRPCRP(3,m)
  100            format(1x,'element ',I5,': powerlaw power ',1PG14.6)
        else
            write(LUMON,*) n,PRPCRP(3,m)
        endif

    elseif (iflag.eq.6+6) then

        if (VERBOSE) then
            write(LUMON,110) n,PRPCRP(6,m)
  110               format(1x,'element ',I5, &
             ': powerlaw activation volume ',1PG14.6)
        else
            write(LUMON,*) n,PRPCRP(6,m)
        endif

    elseif (iflag.eq.12+1) then

        nptype=INT(PRPLAS(1,m))

        if (nptype.eq.0) then
            string = 'none'
        elseif (nptype.eq.1) then
            string = 'Tresca'
        elseif (nptype.eq.2) then
            string = 'Von Mises'
        elseif (nptype.eq.3) then
            string = 'Mohr-Coulomb'
        elseif (nptype.eq.4) then
            string = 'Drucker-Prager'
        else
            string = 'unknown'
        endif

        if (VERBOSE) then
            write(LUMON,120) n,string(1:lnblk(string))
  120            format(1x,'element ',I5,': plasticity type ',12A)
        else
            write(LUMON,*) n,string(1:lnblk(string))
        endif

    elseif (iflag.eq.12+2) then

        if (VERBOSE) then
            write(LUMON,130) n,PRPLAS(2,m)
  130            format(1x,'element ',I5,': friction angle ',1PG12.6)
        else
            write(LUMON,*) n,PRPLAS(2,m)
        endif

    elseif (iflag.eq.12+3) then

        if (VERBOSE) then
            write(LUMON,140) n,PRPLAS(3,m)
  140               format(1x,'element ',I5,': (equivalent) uniaxial ', &
             'yield stress ',1PG12.6)
        else
            write(LUMON,*) n,PRPLAS(3,m)
        endif

    elseif (iflag.eq.12+4) then

        if (VERBOSE) then
            write(LUMON,150) n,PRPLAS(4,m)
  150            format(1x,'element ',I5,': fluidity ',1PG12.6)
        else
            write(LUMON,*) n,PRPLAS(4,m)
        endif

    elseif (iflag.eq.12+5) then

        if (VERBOSE) then
            write(LUMON,160) n,PRPLAS(5,m)
  160            format(1x,'element ',I5,': yield power ',1PG12.6)
        else
            write(LUMON,*) n,PRPLAS(5,m)
        endif

    elseif (iflag.eq.12+6) then

        nflow=INT(PRPLAS(6,m))

        if (nflow.eq.0) then
            string = 'Exponential'
        elseif (nflow.eq.1) then
            string = 'Powerlaw'
        else
            string = 'unknown'
        endif

        if (VERBOSE) then
            write(LUMON,170) n,string(1:lnblk(string))
  170            format(1x,'element ',I5,': plastic flow law ',12A)
        else
            write(LUMON,*) n,string(1:lnblk(string))
        endif

    elseif (iflag.eq.12+7) then

        if (VERBOSE) then
            write(LUMON,180) n,PRPLAS(7,m)
  180               format(1x,'element ',I5,': strain hardening ', &
             1PG12.6)
        else
            write(LUMON,*) n,PRPLAS(7,m)
        endif

    elseif (iflag.eq.21+1) then

        if (VERBOSE) then
            write(LUMON,190) n,PRPTEM(1,m)
  190               format(1x,'element ',I5,': heat production ', &
             1PG12.6)
        else
            write(LUMON,*) n,PRPTEM(1,m)
        endif

    elseif (iflag.eq.21+2) then

        if (VERBOSE) then
            write(LUMON,200) n,PRPTEM(2,m)
  200               format(1x,'element ',I5,': specific heat Cp ', &
                 1PG12.6)
        else
            write(LUMON,*) n,PRPTEM(2,m)
        endif

    elseif (iflag.eq.21+3) then

        if (VERBOSE) then
            write(LUMON,210) n,PRPTEM(3,m)
  210               format(1x,'element ',I5, &
             ': volumetric expansion coefficient ',1PG12.6)
        else
            write(LUMON,*) n,PRPTEM(3,m)
        endif

    elseif (iflag.eq.21+4) then

        if (VERBOSE) then
            write(LUMON,220) n,PRPTEM(4,m)
  220               format(1x,'element ',I5,': conductivity Kx ', &
             1PG12.6)
        else
            write(LUMON,*) n,PRPTEM(4,m)
        endif

    elseif (iflag.eq.21+5) then

        if (VERBOSE) then
            write(LUMON,230) n,PRPTEM(5,m)
  230               format(1x,'element ',I5,': conductivity Ky ', &
             1PG12.6)
        else
            write(LUMON,*) n,PRPTEM(5,m)
        endif

    endif
enddo

return
end subroutine
!-------------------------------------------------------------------------------
subroutine nummat (IEN,X,D,MAT,PRPMAT,STN,LMF,TFAULT,DXE, &
                   SKEW,NEN,NSD,NDOF,NSTR,DEFORM,deffac,lu,ierr, &
                   luMon,verbose)

USE ALGEBRA
USE MODELDEFINITION
use materialsmodule, only: MaterialIDOfElement
use constants,       only: zero, half
use timestepmodule,  only: nintg, maxit
use materialsmodule, only: lmat
use iomodule,        only: stderr
use monitorVars,     only: elmon, monel, nnel
use plotControl,     only: Vscale, viscli, vislog
use errors,          only: ok
!
implicit none

!-pass
logical DEFORM,verbose
double precision :: deffac
integer NEN,NSD,NDOF,NSTR,lu,ierr, lumon
integer          :: ien, mat, lmf
double precision :: x, d, prpmat, stn, tfault, dxe, skew
dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),MAT(*),PRPMAT(6,*), &
          STN(NSTR,*),LMF(NDOF,NEN,*),TFAULT(NDOF,*),DXE(NDOF,NEN,*), &
          SKEW(2,*)



!-locl
double precision :: def, anpwr, eff, emhu, value
integer :: m, n, nn, i
double precision :: CENTR(3)
double precision, external :: effstn, etaLin

!-init
ierr = OK
def  = ZERO
if (DEFORM) then
    def = deffac
endif

if (elmon) then
    do nn=1,nnel

        n = monel(nn)
        m = MaterialIDOfElement(n)

        emhu  = PRPMAT(3,m)
        anpwr = PRPMAT(4,m)

        if (viscli) then
            eff  = effstn (STN(1,n),HALF,.false.)
            value = etaLin (eff,emhu,anpwr)
        else
            value = emhu
        endif

        if (Vscale.gt.ZERO) then
            value = value*Vscale
        endif

        if (vislog) then
            value = TENLOG(ABS(value))
        endif

        if (verbose) then
            if (viscli) then
                if (vislog) then
                    write(luMon,1) n,value
 1                       format(1x,'element ',I5,': linearized LOG ', &
                     1PG15.6)
                else
                    write(luMon,2) n,value
 2                       format(1x,'element ',I5,': linearized ', &
                     1PG15.6)
                endif
            else
                if (vislog) then
                    write(luMon,3) n,value
 3                format(1x,'element ',I5,': LOG ',1PG15.6)
                else
                    write(luMon,4) n,value
 4                format(1x,'element ',I5,': ',1PG15.6)
                endif
            endif
        else
            write(luMon,*) value
        endif
    enddo
else

    if (viscli) then
        if (vislog) then
            write(lu,5)
 5            format(1x,'* Linearized Logarithmic viscosities *')
        else
            write(lu,6)
 6            format(1x,'*** Linearized viscosities ***')
        endif
    else
        if (vislog) then
            write(lu,7)
 7            format(1x,'*** Logarithmic function values ***')
        else
            write(lu,8)
 8            format(1x,'*** straight function values ***')
        endif
    endif

    do nn=1,NUMEL

        n = nn
        m = MaterialIDOfElement(n)

        emhu  = PRPMAT(3,m)
        anpwr = PRPMAT(4,m)

!             calculate element center coordinates
        call ELCTR (CENTR,IEN(1,n),X,D,n,TFAULT, &
                    SKEW, def,ierr)
        if (ierr.ne.0) then
            write(stderr,13) n
13                format(1x,'nummat: error in element ',I5)
            call stoper()
        endif

        if (viscli) then
            eff   = effstn (STN(1,n),HALF,.false.)
            value = etaLin (eff,emhu,anpwr)
        else
            value = EMHU
        endif

        if (Vscale.gt.ZERO) then
            value = value*Vscale
        endif

        if (vislog) then
            value = TENLOG(ABS(value))
            write(lu,9,advance='no') (CENTR(i),i=1,NSD)
 9            format(1x,'pt -c',1PG15.6,1x,1PG15.6,1x)
            write(lu,10) value
10            format('@10 sup ',0PF4.1,'@')
        else
            write(lu,9) (CENTR(i),i=1,NSD)

            if (ABS(value).lt.1E2) then
                write(lu,11) value
11            format(0PF4.1)
            else
                write(lu,12) value
12            format(1PE7.1)
            endif
        endif
    enddo
endif

return
end subroutine
