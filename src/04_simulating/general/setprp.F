subroutine setprp (IMPRINT,NSD,NDOF,NEN,NSTR)

USE MESHDATAMODULE
USE MATERIALSMODULE
USE MODELDEFINITION
USE MODELDATAMODULE
USE MODELCTX
use algebra,        only: idig
use debugmodule,    only: xit

!
! Routine for calculating material properties
!
implicit none
!-pass
integer NSD,NDOF,NEN,NSTR
integer :: imprint
dimension IMPRINT(*)

!       Determine linearized viscosity minimum
!-local
integer :: iElem
double precision :: strainRate(NSTR)

if (IVLIM.eq.2) then
    if (LMAT) then
        call vismin (matdatactx%PRPMAT,NUMEL)
    else
        call vismin (matdatactx%PRPMAT,NUMAT)
    endif
endif
!
if (MODE.eq.6 .and. IDIG(ICVIS,4,1).ne.0) then
!           Calculate non-linear viscosities from temperatures
! can be power-law or fully fledged non-linear viscosity
    call setvisc (meshdatactx%IEN,meshdatactx%X,modeldatactx%D, &
       modeldatactx%T,modeldatactx%TX, &
     modeldatactx%STN,modeldatactx%STR,modeldatactx%LMF, &
     modeldatactx%TFAULT,modeldatactx%DXE, &
     modeldatactx%LMTX,modeldatactx%SKEW,NSD,NSTR,NDOF, &
     NEN)
endif
!
if (IDIG(ICVIS,4,3).eq.1) then
!           Compute viscous strain weakening
    call VWEAK (modeldatactx%STR,NSTR)
else if (IDIG(ICVIS,4,3).eq.2) then
!           Compute viscous strain rate weakening

#ifdef EXPERIMENTAL_NORATESTORAGE
    do iElem=1,meshdatactx%neglobal
        call StrainrateOfElem (meshdatactx%X,&
                              modeldatactx%D,&
                              modeldatactx%DELD, &
                               strainRate,&
                               meshdatactx%IEN,&
                              modeldatactx%LMF, &
                              modeldatactx%TFAULT,&
                              modeldatactx%DFAULT,&
                              modeldatactx%LMX, &
                              modeldatactx%DX, &
                              modeldatactx%SKEW,&
                              modeldatactx%DXE, &
                               iElem)

        call VWEAKStrainRate (strainRate,NSTR)
    enddo
#else
!    call VWEAKStrainRate (modeldatactx%RATE,NSTR)
    call VWEAK (modeldatactx%RATE,NSTR)
#endif

endif
!
if (IVLIM.ne.0) then
    call VisLim (meshdatactx%mat,matdatactx%PRPMAT, &
  matdatactx%PRPCRP, &
    modeldatactx%STN,NSTR)
endif
!
if (IDIG(ICVIS,4,4).eq.1) then
!   Compute yield weakening from visco-plastic strain
    call CWEAK (meshdatactx%mat,matdatactx%PRPMAT, &
                matdatactx%PRPLAS, &
                matdatactx%PRPDIF,modeldatactx%VPSTR,NSTR)
else if (IDIG(ICVIS,4,4).eq.2) then
!   Compute yield weakening from visco-plastic strain rate
    call CWEAK (meshdatactx%mat,matdatactx%PRPMAT, &
                matdatactx%PRPLAS, &
                matdatactx%PRPDIF,modeldatactx%BETB,NSTR)
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine setvisc (IEN,X,D,T,TX,STN, &
                    STR,LMF,TFAULT,DXE,LMTX,SKEW,NSD,NSTR, &
                    NDOF,NEN)

USE MATERIALSMODULE, only: LMAT, &
                           matdatactx, &
                           MaterialIDOfElement
#ifdef EXPERIMENTALTRANSIENTDISCREEP                          
use modeldatamodule, only: modeldatactx
#endif
USE MODELCTX
USE MODELDEFINITION
use constants,       only: dfmin,dlogmin,dfmax,dlogmax, half, useDXE
use iomodule
use algebra,         only: idig

! calculates viscosities as function of temperature and pressure
! for each element

implicit none
!-pass-
integer NSD,NEN,NDOF,NSTR

integer :: ien, lmf, lmtx
double precision :: x, d, t, stn, str, tfault, tx, skew, dxe

dimension IEN(NEN,*),X(NSD,*),D(NDOF,*),T(*), &
 STN(NSTR,*),STR(NSTR,*), &
 LMF(NDOF,NEN,*),TFAULT(NDOF,*),TX(*),LMTX(NEN,*), &
 SKEW(2,*),DXE(NDOF,NEN,*)

double precision :: xl, dl

integer :: nn,n,m, ierr
double precision :: pois, Tc, Estrs, P

double precision, external :: EFFSTN, PRESSURE, Tcentr

!-locl-
dimension xl(24),dl(24)
!-init
ierr = 0
! loop over elements

do nn=1,NUMEL
    n = nn
    m = MaterialIDOfElement(n)
!
!           localize coordinate data
    call LCOORD (X,xl,IEN(1,n),nn)
    if (LGDEF.ne.0) then
        call LDISP  (dl,D,IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (dl,n,TFAULT,NDOF,NEN,12)
        call ADDSNE (dl,n,NDOF,NEN,useDXE)
        call REZONE (xl,dl,1d0)
    endif

!           determine temperature [C] in center of element
    Tc = Tcentr (T,TX,XL,IEN(1,n),LMTX(1,n),ierr)
    if (ierr.ne.0) then
        write(stderr,1) n
 1           format(1x,'setvisc: error in element ',I5)
        call stoper()
    endif
!
!           Effective stress from deviatoric incompressible tensor
    Estrs = EFFSTN (STN(1,n),HALF,.false.)
!
!           Compute pressure from incompressible stress tensor
    P = PRESSURE(STN(1,n),HALF)


#ifdef EXPERIMENTALTRANSIENTDISCREEP
if (IDIG(ICVIS,4,1).eq.2) then
! Calculate viscosity with composite transient/steady-state dislocation
! creep and steady-state diffusion creep
    
! compute viscosity
!    call SetCompositeVisc (matdatactx%PRPDISLC(1:10,m),matdatactx%PRPDIFFC(1:5,m),\
!    matdatactx%PRPMAT(1:6,m),modeldatactx%STNTAY(1,m),modeldatactx%STNTHRES(1,m),\
!    Tc,P,Estrs,modeldatactx%InvEffViscDisl(m),modeldatactx%InvEffViscDiff(m))
        call SetCompositeVisc (matdatactx%PRPDISLC(1:10,m),matdatactx%PRPDIFFC(1:5,m),\
        matdatactx%PRPMAT(1:6,m),modeldatactx%STNTAY(1,m),modeldatactx%STNTHRES(m),\
    Tc,P,Estrs,modeldatactx%InvEffViscDisl(m),modeldatactx%InvEffViscDiff(m))
    
else
   !           Calculate effective viscosity with power-law rheology
    call SetMhu (matdatactx%PRPCRP(1,m),matdatactx%PRPMAT(1,m),Tc,P,Estrs)
endif
#else
!           Calculate effective viscosity with power-law rheology
    call SetMhu (matdatactx%PRPCRP(1,m),matdatactx%PRPMAT(1,m),Tc,P,Estrs)
#endif 

enddo
!
return
end
!-------------------------------------------------------------------------------
subroutine CWEAK (MAT,PRPMAT,PRPLAS,PRPDIF,TENSOR,NSTR)
!
! Adjusts uniaxial yield stress for strain/strain rate weakening
! TENSOR is either the strain tensor or the strain rate tensor
!
USE MATERIALSMODULE, only : LMAT, &
                            matdatactx, &
                            MaterialIDOfElement
USE MODELCTX
USE MODELDEFINITION
use constants, only: dfmin,dlogmin,dfmax,dlogmax

implicit none
!-pass-
integer NSTR, mat
double precision :: prpmat, prplas, prpdif, tensor
dimension MAT(*),PRPMAT(6,*),PRPLAS(9,*),PRPDIF(9,*),TENSOR(NSTR,*)
!-locl
logical :: firstentry
save    :: firstentry
integer :: n,m
double precision :: pois, eff
double precision, external :: EFFEE
integer :: ierr
!-init
data firstentry /.true./
ierr = 0
 
do n=1,NUMEL
    m = MaterialIDOfElement(n)
    pois   = matdatactx%PRPMAT(2,m)
    eff = EFFEE (TENSOR(1,n),pois,.false.)
    if (firstentry) then
!       backup (input) uniaxial yield stress
        PRPDIF(8,m) = matdatactx%PRPLAS(3,m)
    else
!       retrieve (input) uniaxial yield stress
        PRPLAS(3,m) = matdatactx%PRPDIF(8,m)
    endif
    call YLDWEAK (matdatactx%PRPDIF(1,m),matdatactx%PRPLAS(1,m),eff)
enddo
firstentry = .false.
 
return
end
!-------------------------------------------------------------------------------

#ifdef EXPERIMENTAL_NORATESTORAGE

subroutine VWEAKStrainRate (TENSOR,NSTR, elemID)

USE MATERIALSMODULE, only : LMAT, &
                            matdatactx, &
                            MaterialIDOfElement
USE MODELDEFINITION, only : mode, icvis, numel
use algebra,         only : idig
use constants, only: dfmin,dlogmin,dfmax,dlogmax


! Adjusts viscosities for strain/strain rate weakening
! TENSOR is either the strain tensor or the strain rate tensor
! depending on the caller of this sub.

implicit none
!-pass-
integer NSTR

double precision TENSOR
dimension TENSOR(NSTR)
integer :: elemID
!-locl
logical :: firstentry
save    :: firstentry
integer :: ierr
integer :: nn,n,m
double precision :: pois, eff
double precision, external :: effee
!-init
data firstentry /.true./

ierr = 0

m = MaterialIDOfElement(elemID)

pois = matdatactx%PRPMAT(2,m)
eff  = EFFEE (TENSOR(1),pois,.false.)

if (MODE.ne.6 .or. IDIG(ICVIS,4,1).eq.0) then
    if (firstentry) then
        matdatactx%PRPDIF(9,m) = matdatactx%PRPMAT(3,m)
    else
        matdatactx%PRPMAT(3,m) = matdatactx%PRPDIF(9,m)
    endif
endif

call VISWEAK (matdatactx%PRPDIF(1,m),matdatactx%PRPMAT(1,m),eff)

firstentry = .false.
!
return
end  


#endif

!-------------------------------------------------------------------------------

subroutine VWEAK (TENSOR,NSTR)

USE MATERIALSMODULE, only : LMAT, &
                            matdatactx, &
                            MaterialIDOfElement
USE MODELDEFINITION, only : mode, icvis, numel
use algebra,         only : idig
use constants, only: dfmin,dlogmin,dfmax,dlogmax


! Adjusts viscosities for strain/strain rate weakening
! TENSOR is either the strain tensor or the strain rate tensor
! depending on the caller of this sub.

implicit none
!-pass-
integer NSTR

double precision TENSOR
dimension TENSOR(NSTR,*)
!-locl
logical :: firstentry
save    :: firstentry
integer :: ierr
integer :: nn,n,m
double precision :: pois, eff
double precision, external :: effee
!-init
data firstentry /.true./
ierr = 0

do nn=1,NUMEL
    n = nn

    m = MaterialIDOfElement(n)

    pois = matdatactx%PRPMAT(2,m)
    eff  = EFFEE (TENSOR(1,n),pois,.false.)

    if (MODE.ne.6 .or. IDIG(ICVIS,4,1).eq.0) then
        if (firstentry) then
            matdatactx%PRPDIF(9,m) = matdatactx%PRPMAT(3,m)
        else
            matdatactx%PRPMAT(3,m) = matdatactx%PRPDIF(9,m)
        endif
    endif

    call VISWEAK (matdatactx%PRPDIF(1,m),matdatactx%PRPMAT(1,m),eff)
enddo
firstentry = .false.
!
return
end

!-------------------------------------------------------------------------------
subroutine SetMhu (PRPCRP,PRPMAT,Tc,P,Estrs)

use constants, only: TK0

! calculates effective viscosity EMHU by powerlaw creep
!
!       Qpl     = Powerlaw activation energy in kJ/mole
!       Apl     = Powerlaw Pre-exponential constant [Pa^-n.s^-1]
!       Vpl     = Powerlaw activation volume [m^3/mole]
!       ANPL    = Powerlaw power for stress
!       Tc      = Temperature in Celsius
!       P       = Total pressure
!       Estrs   = Effective stress

implicit none
!-pass
dimension PRPCRP(6),PRPMAT(6)
double precision :: prpmat, prpcrp
double precision :: Tc, P, Estrs
!-local
double precision :: Qpl, Vpl, Apl, ANPL, TK
double precision :: EMHU, ANPWR
double precision, external :: visco

!-init
Qpl  = PRPCRP(1)*1d3
Vpl  = PRPCRP(6)
Apl  = PRPCRP(2)
ANPL = PRPCRP(3)
TK   = Tc + TK0

#ifdef LINPL
EMHU = visco (Qpl,Vpl,Apl,ANPL,Estrs,TK,P,.true.)
ANPWR = 1d0
#else
EMHU = visco (Qpl,Vpl,Apl,ANPL,Estrs,TK,P,.false.)
ANPWR = ANPL
#endif

PRPMAT(3) = EMHU
PRPMAT(4) = ANPWR
!
return
end
!-------------------------------------------------------------------------------


#ifdef EXPERIMENTALTRANSIENTDISCREEP
subroutine SetCompositeVisc (PRPDISLC,PRPDIFFC,PRPMAT,StnTay,StnThres,Tc,P,Estrs,\
InvEffViscDisl,InvEffViscDiff)
use constants, only: TK0, one
use modeltopology, only: NSTRTAY
! calculates effective viscosity for both dislocation creep as well as diffusion creep

implicit none
!-pass
dimension PRPDISLC(10),PRPDIFFC(5),PRPMAT(6)
dimension StnTay(NSTRTAY)
double precision :: prpmat, prpdislc, prpdiffc
double precision :: Tc, P, Estrs ! temperature in C, pressure, effective stress
double precision :: InvEffViscDisl,InvEffViscDiff ! effective viscosities for dislocation and diffusion creep, respectively
double precision :: StnTay,StnThres
!-local
double precision :: QDislC, ADislC, SigmaP, TK
double precision :: M, R0GBD, R0Pipe, VDislC,d
double precision :: QDiffC, ADiffC, Pexp, VDiffC
double precision :: EffVisc
double precision, external :: InvViscoDisl,InvViscoDiff

!-init
! 
! TRANSIENT/STEADY-STATE DISLOCATION CREEP MATERIAL PROPERTIES
!
    QDislc=PRPDISLC(1)       ! Activation energy.
    ADislc=PRPDISLC(2)       ! Pre-exponent constant.
    SigmaP=PRPDISLC(3)       ! Peierls stress
    M = PRPDISLC(6)          ! M coefficient Taylor stress.
    R0GBD = PRPDISLC(7)      ! Rate coefficient dislocation change due to grain boundary diffusion.
    R0Pipe = PRPDISLC(8)     ! Rate coefficient dislocation change due to pipe diffusion
    VDislc = PRPDISLC(9)     ! Activation volume

! STEADY-STATE DIFFUSION CREEP MATERIAL PROPERTIES

    QDiffc = PRPDIFFC(1)     ! Diffusion creep activation energy.
    ADiffc = PRPDIFFC(2)     ! Diffusion creep scaling constant.
    pExp = PRPDIFFC(3)       ! Grain size exponent.
    VDiffC = PRPDIFFC(4)     ! Diffusion creep activation volume.
    d = PRPDIFFC(5)          ! grain size


! convert temperature to kelvin
TK   = Tc + TK0

! compute inverse viscosity for dislocation creep
InvEffViscDisl = InvViscoDisl (QDislc,ADislc,VDislc,SigmaP,StnTay(1),StnThres,\
                 TK,P,Estrs)

! compute inverse viscosity for diffusion creep
InvEffViscDiff = InvViscoDiff (QDiffc,ADiffc,VDiffc,d,pExp,TK,P)

! compute composite viscosity
EffVisc = one/(InvEffViscDisl+InvEffViscDiff)


! write back to PRPMAT
PRPMAT(3) = EffVisc
PRPMAT(4) = one ! set power to one

!
return
end


#endif
!-------------------------------------------------------------------------------
#ifdef EXPERIMENTALTRANSIENTDISCREEP
double precision function InvViscoDisl (Q,A,V,SigmaP,StnTay,StnThres,\
                 TK,P,Estrs)
! transient+steady-state dislocation creep, which is dependent on internal stress (i.e. Taylor stress) related to variable
! dislocation density
! based on Breithaupt et al. 2023 https://doi.org/10.1073/pnas.2203448120

! computes inverse viscosity, i.e. 1/viscosity
use constants, only: R,sqrtthree,zero
implicit none
!-pass
double precision :: Q,V,A,SigmaP,StnTay,StnThres,TK,P,Estrs
! Q : activation energy 
! V : activation volume
! A : constant
! SigmaP : peierls stress
! StnTay : taylor stress
! StnThres : threshold stress
! TK : temperature in Kelvin
! P : pressure
! Estrs : effective stress


! inverse viscosity
if (sqrtthree*Estrs.gt.StnThres) then
! pressure has not yet been included
   invviscodisl = sqrtthree*A*exp(-Q/(R*TK))*StnTay**2*sinh(Q/(R*TK)*(sqrtthree*Estrs-StnTay-StnThres)/SigmaP)/Estrs
else
! if effective stress (stress invariant) is less than threshold stress, return zero inverse viscosity
   invviscodisl = 1e-40
endif

return 

end
#endif
!-------------------------------------------------------------------------------
#ifdef EXPERIMENTALTRANSIENTDISCREEP
double precision function InvViscoDiff (Q,A,V,d,pExp,TK,P)
! steady-state diffusion creep. 

! computes inverse viscosity, i.e. 1/viscosity
use constants, only: R,three
implicit none
!-pass
double precision :: Q,V,A,d,pExp,TK,P,Estrs
! Q : activation energy 
! V : activation volume
! A : constant
! d : grain size
! pExp : grain size exponent
! TK : temperature in Kelvin
! P : pressure

! pressure has not yet been included

! inverse viscosity
InvViscoDiff = three*A*d**(-pExp)*exp(-Q/(R*TK))

return 

end
#endif
!-------------------------------------------------------------------------------
#ifdef EXPERIMENTALTRANSIENTDISCREEP
subroutine SetTaylorStress(STNTAY,STNTAYRATE,\
    RATEDISLEFF,STNTHRES,PRPDISLC,PRPDIFFC,T)

! calculates increment of Taylor stress, the internal stress resulting from dislocation density
! based on Breithaupt et al. 2023 https://doi.org/10.1073/pnas.2203448120

use constants, only: TK0,R,sqrtthree,half,useDXE
USE MODELDEFINITION, only: LGDEF
use timestepmodule,  only: DELTP
use modeltopology, only: NSTRTAY ,NDOF,NEN
use modeldatamodule, only: modeldatactx
USE MESHDATAMODULE,  only: meshdatactx
use materialsmodule, only: MaterialIDOfElement
use iomodule
! 
implicit none
!-pass
dimension PRPDISLC(10,*),STNTAY(NSTRTAY,*),STNTAYRATE(NSTRTAY,*),STNTHRES(*)
dimension RATEDISLEFF(*),T(*),PRPDIFFC(5,*)
double precision :: PRPDISLC,STNTAY,STNTAYRATE,RATEDISLEFF,STNTHRES,PRPDIFFC
double precision :: T
! local
double precision :: EffStress,signEffStress,ExpminQRT
double precision :: Q,MTay,R0GBD,R0Pipe,E,d
double precision :: Tc,TK,P,Estrs
double precision :: xl(24),dl(24)
integer :: ierr,n,nn,m

double precision, external :: EFFSTN, PRESSURE, Tcentr  

if (NSTRTAY.ne.1) then
       write(stderr,1) 
 2           format(1x,'SetTaylorStress: error dimensions of Taylor stress, only scalar stress is implemented')
        call stoper()
endif
! loop over elements
do nn=1,meshdatactx%nelocal
    n = nn
    m = MaterialIDOfElement(n)

    !localize coordinate data
    call LCOORD (meshdatactx%X,xl,meshdatactx%IEN(1,n),nn)
    if (LGDEF.ne.0) then
        call LDISP  (dl,modeldatactx%D,meshdatactx%IEN(1,n),NDOF,NEN)
        call AddFaultDisplacement (dl,n,modeldatactx%TFAULT,NDOF,NEN,12)
        call ADDSNE (dl,n,NDOF,NEN,useDXE)
        call REZONE (xl,dl,1d0)
    endif

    !   determine temperature [C] in center of element
    Tc = Tcentr (T,modeldatactx%TX,XL,meshdatactx%IEN(1,n),modeldatactx%LMTX(1,n),ierr)
    if (ierr.ne.0) then
        write(stderr,1) n
 1           format(1x,'SetTaylorStress: error in element ',I5)
        call stoper()
    endif
!
!   Effective stress from deviatoric incompressible tensor
    Estrs = EFFSTN (modeldatactx%STN(1,n),HALF,.false.)
!
!   Compute pressure from incompressible stress tensor
    P = PRESSURE(modeldatactx%STN(1,n),HALF)
    

    ! convert temperature to Kelvin
    TK = Tc + TK0

    ! get parameters from PRPDISLC
    Q=PRPDISLC(1,m)! activation energy in [J/mol]
    MTay=PRPDISLC(6,m)! constant for Taylor stress change [Pa]
    R0GBD=PRPDISLC(7,m)! coefficient for dislocation change due to grain boundary diffusion [Pa**4]
    R0Pipe=PRPDISLC(8,m)! coefficient for dislocation change due to pipe diffusion [Pa**5]
    E=PRPDISLC(9,m) ! not yet used, activation volume in [m3/mol]

    !d=PRPDIFFC(5,m) ! grain size [in meters]

    ExpminQRT=exp(-Q/(R*TK))
    ! effective stress
    EffStress=sqrtthree*Estrs-STNTAY(1,m)-STNTHRES(m)
    ! sign of effective stress
    signEffStress=abs(EffStress)/EffStress

    ! rate of Taylor stress
    STNTAYRATE(1,m)=MTay*(signEffStress*(STNTAY(1,m)+STNTHRES(m))/STNTAY(1,m)*2/sqrtthree*RATEDISLEFF(m) \
    -R0GBD*ExpminQRT*STNTAY(1,m)**3*STNTHRES(m)\
    -R0Pipe*ExpminQRT*STNTAY(1,m)**5)

    ! update of Taylor stress, by integrating the rate of stress change
    STNTAY(1,m)=STNTAY(1,m)+STNTAYRATE(1,m)*DELTP

  !  write(*,*) 'Q, MTAY,R0GBD,R0Pipe,StayRate,STay'
  !  write(*,*) Q,MTay,R0GBD,R0Pipe,STNTAYRATE(1,m),STNTAY(1,m)
enddo

return
end


#endif
!-------------------------------------------------------------------------------
#ifdef EXPERIMENTALTRANSIENTDISCREEP
subroutine InitTaylorStess(STNTAY)
! this subroutine sets an initial Taylor stress
! The Taylor stress is related to dislocation density as sigma = alpha * mu * b * sqrt(rho)
! with alpha a coefficient, mu shear modulus, b Burgers length and rho dislocation density
! without initial dislocation density, there is no dislocation motion possible. 
! Taylor stress (and dislocation density) thus needs an intial value
! If it is too low, it takes a long time before Taylor stress reaches its steady state value
! If it is too high (higher than the element effective stress) it will drive deformation (instead of the environment stress due to loading)
! Therefore we set it initially at a user-supplied fraction of 



USE MATERIALSMODULE, only : matdatactx,MaterialIDOfElement
use meshdatamodule, only: meshdatactx
USE MODELDATAMODULE, only: modeldatactx
use modeltopology, only: NSTRTAY
use constants, only: Half,sqrtthree

implicit none

! - pass
double precision :: STNTAY
dimension STNTAY(NSTRTAY,*)
! - local
integer :: nn,n,m
double precision :: Estrs,FractionElastic
double precision, external :: EFFSTN

do nn=1,meshdatactx%nelocal
    n = nn
    m = MaterialIDOfElement(n)

    ! fraction to elastic stress
    FractionElastic = matdatactx%PRPDISLC(10,m)
    ! Effective stress from deviatoric incompressible tensor
    Estrs = EFFSTN (modeldatactx%STN(1,n),HALF,.false.) ! half and false are not used for 3D
    ! use fraction of effective stress
    ! square root of three comes from the factor sqrt(3) in the constitutive equation when stress is written as effective deviatoric stress
    STNTAY(1,m) = sqrtthree*Estrs*FractionElastic
enddo

end

#endif
!-------------------------------------------------------------------------------
#ifdef EXPERIMENTALTRANSIENTDISCREEP

subroutine SetThresholdStress(STNTHRES)
! set the grain size dependent threshold stress
! threshold stress functions as a minimum stress that has to be overcome before
! dislocation creep can be active

use meshdatamodule, only: meshdatactx
USE MATERIALSMODULE, only : matdatactx, MaterialIDOfElement
use constants, only: one, half
                           
implicit none
!-pass
double precision :: STNTHRES
dimension STNTHRES(*)

! local
integer :: nn,n,m
double precision :: d,Beta,Mu,E,Pois,b


do nn=1,meshdatactx%nelocal
    n = nn
    m = MaterialIDOfElement(n)
    
    d = matdatactx%PRPDIFFC(5,m) ! grain size
    b = matdatactx%PRPDISLC(4,m) ! burgers length
    Beta = matdatactx%PRPDISLC(5,m) ! Beta Threshold Stress Coefficient
    E = matdatactx%PRPMAT(1,m) ! Young's modulus
    Pois = matdatactx%PRPMAT(2,n) ! Poisson ratio
    !    Compute Lame parameter Mu = shear modulus
    Mu=E/(One+Pois)*half
    ! threshold stress
    STNTHRES(m)=Beta*Mu*b/d
    
enddo
    
return
end
#endif

!-------------------------------------------------------------------------------



subroutine YLDWEAK (PRPDIF,PRPLAS,EFF)

use constants, only: zero, one, two
!
! calculates yield stress weakening
!
implicit none
!-pass
double precision :: EFF
double precision :: PRPDIF, PRPLAS
dimension PRPDIF(9),PRPLAS(9)
!-local
double precision :: weak, sigma, arg, fac

!-init
weak = PRPDIF(3)
if (weak.gt.ZERO .and. weak.ne.ONE) then
    sigma = PRPDIF(4)
    arg = EFF/sigma
    arg = arg*arg/TWO
    fac = weak+(ONE-weak)*DEXP(-arg)
    PRPLAS(3) = MAX(PRPLAS(3)*fac,5d6)
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine VISWEAK (PRPDIF,PRPMAT,EFF)

! calculates viscous weakening
!
use constants, only: zero, one, two

implicit none
!-pass
double precision :: EFF
double precision :: PRPDIF, PRPMAT
dimension PRPDIF(9),PRPMAT(6)
!-local
double precision :: weak, sigma, arg, fac

!-init
weak = PRPDIF(1)
if (weak.gt.ZERO) then
    sigma = PRPDIF(2)
    arg = EFF/sigma
    arg = arg*arg/TWO
    fac = weak+(ONE-weak)*DEXP(-arg)
    PRPMAT(3) = PRPMAT(3)*fac
endif
!
return
end
!-------------------------------------------------------------------------------
subroutine vismin (PRPMAT,NUMAT)

use materialsmodule, only: matdatactx
use timestepmodule,  only: DELTP,ALFAP,TALFP
use constants,       only: dfmin,dlogmin,dfmax,dlogmax, &
                           one, two, three, six
use viscoLimits,     only: viscoMin, viscoMax

!
! Routine for setting the minimum to (non-linear) viscosity, based upon
! Oscillation, Error and Stability criteria at a given time step size.
!
implicit none

double precision, parameter :: TOLER=0.01d0
!-pass
integer NUMAT
double precision :: PRPMAT
dimension PRPMAT(6,*)
!-local
double precision :: c, acc, osc, stb, eigmin
double precision :: E, Pois
integer :: n

!
!    Approximate minimum eigenvalue for local truncation error
c = (three * ALFAP * (ALFAP-one)+one)/six * DELTP*DELTP
acc = (ABS(TOLER)/c)**(one/three)
!
!    Minimum for Oscillation limit
c = (one-ALFAP)
if (c.gt.dfmin) then
    osc = one/(c*DELTP)
else
    osc = dfmax
endif
!
!    Minimum for stability
c = one - two*ALFAP
if (c.gt.dfmin) then
    stb = two/(c*DELTP)
else
    stb = dfmax
endif
!
eigmin = acc
if (eigmin.gt.osc) then
    eigmin = osc
endif
if (eigmin.gt.stb) then
    eigmin = stb
endif
!
!    find minimum Young modulus and Poison ratio
E = dfmax
POIS = -dfmax
do n=1,NUMAT
    if (matdatactx%PRPMAT(1,n).lt.E) then
        E = matdatactx%PRPMAT(1,n)
    endif
    if (matdatactx%PRPMAT(2,n).gt.POIS) then
        POIS = matdatactx%PRPMAT(2,n)
    endif
enddo
!
!    linearized viscosity minimum
! todo these values are set but never used...?
viscoMin = E/(one+POIS)/(one-two*POIS)/eigmin
viscoMax = viscoMin * 1d12
!
return
end
!-----------------------------------------------------------------------
subroutine WrtCrp (PRPMAT,NUMAT)

use materialsmodule, only: matdatactx
use debugmodule,     only: xit
use iomodule

! Writes calculated viscosities and powers to file

implicit none
!-pass
integer NUMAT
double precision :: PRPMAT
dimension PRPMAT(6,*)
!-locl
 character(len=1) quote
integer, external :: lnblk
integer :: i, irec, ios
double precision :: emhu, anpwr

!-init
irec=51
 quote=CHAR(34)

if (OUTFIL(1)) then
    do i=1,NUMAT
        EMHU  = matdatactx%PRPMAT(3,i)
        ANPWR = matdatactx%PRPMAT(4,i)
        if (irec.ge.51) then
            irec = 0
            write(luout(1),2,err=100,iostat=ios)
 2               format(///1X,'V I S C O S I T Y'//5X, &
             'SET',6X,'VISCOSITY',4X,'POWER'/1X,'NUMBER'/)
        endif
        irec = irec + 1
        write (luout(1),3,err=100,iostat=ios) i,EMHU,ANPWR
 3           format(I12,5X,1PG10.2,5X,0PF4.2)
    enddo
endif

if (OUTFIL(2)) then
!             write (*,*) 'Going to write powers and viscosities'
    if (ascout) then
        do i=1,NUMAT
            EMHU  = matdatactx%PRPMAT(3,i)
            ANPWR = matdatactx%PRPMAT(4,i)
            write(luout(2),4,err=110,iostat=ios) i,EMHU,ANPWR
 4               format(I12,1X,1PG15.6,1X,0PF4.2)
        enddo
    else
        write(luout(2),err=110,iostat=ios) &
         (PRPMAT(3,i),PRPMAT(4,i),i=1,NUMAT)
    endif
endif

return

100     write(stderr,101) ios,quote,FEOUT(1:lnblk(FEOUT)),quote
101     format(1x,'WrtCrp: fatal write error ',I3,' on ',80a)
call prerr('WrtCrp')
call xit(1," ")

110     write(stderr,111) ios,quote,FEDSK(1:lnblk(FEDSK)),quote
111     format(1x,'WrtCrp: fatal write error ',I3,' on ',80a)
call prerr('WrtCrp')
call xit(1," ")

 end
!-------------------------------------------------------------------------------
subroutine WrtPls (NUMAT)

! Writes calculated uniaxial yield stress and fluidity to file

! Note that NUMAT can also be the number of elements, when
! LMAT is true, and every element has its own material entry.
! The caller arranges this.

use modeldefinition, only: icvis
use materialsmodule, only: matdatactx
use algebra,         only: idig
use debugmodule,     only: xit
use iomodule

implicit none
!-pass
integer          :: NUMAT
!-locl
character(len=1)  :: quote
integer, external :: lnblk
integer           :: i, irec, ios
double precision  :: uniyld, fluidy
logical           :: hasPlasticity

hasPlasticity = (IDIG(ICVIS,4,2).ne.0)
if (.not. hasPlasticity) then
    return
endif

!-init
irec=51
quote=CHAR(34)

if (OUTFIL(1)) then
!    do 10 i=1,NUMAT
    do i=1,NUMAT
!        write(*,*) "enjoying numat i", i, matdatactx%PRPLAS
!        if (INT(PRPLAS(1,i)).eq.0) goto 10


        if (INT(matdatactx%PRPLAS(1,i)).ne.0) then

            uniyld = matdatactx%PRPLAS(3,i)
            FLUIDY = matdatactx%PRPLAS(4,i)
            if (irec.ge.51) then
                irec = 0
                write(luout(1),2,err=100,iostat=ios)
 2               format(///1X,'P L A S T I C I T Y'//5X, &
             'SET',7X,' UNIYLD ',6X,'FLUIDITY'/2X,'NO.'/)
            endif
            irec = irec + 1
            write (luout(1),3,err=100,iostat=ios) i,uniyld,FLUIDY
 3           format(1X,I12,5X,1PG10.2,4X,1PG10.2)
        endif ! added for positive selection


    enddo
!10          continue
endif

if (OUTFIL(2)) then
!    write (*,*) 'Going to write uniaxial yield and fluidity to FEDSK'
    if (ascout) then
        do i=1,NUMAT
            uniyld  = matdatactx%PRPLAS(3,i)
            FLUIDY  = matdatactx%PRPLAS(4,i)
            write(luout(2),4,err=110,iostat=ios) i,uniyld,FLUIDY
 4               format(I12,1X,1PG15.6,1X,1PG15.6)
        enddo
    else
        write(luout(2),err=110,iostat=ios) &
         (matdatactx%PRPLAS(3,i),matdatactx%PRPLAS(4,i),i=1,NUMAT)
    endif
endif

return

100     write(stderr,101) ios,quote,FEOUT(1:lnblk(FEOUT)),quote
101     format(1x,'WrtPls: fatal write error ',I3,' on ',80a)
call prerr('WrtPls')
call xit(1,"Error while writing FEOUT.DAT")

110     write(stderr,111) ios,quote,FEDSK(1:lnblk(FEDSK)),quote
111     format(1x,'WrtPls: fatal write error ',I3,' on ',80a)
call prerr('WrtPls')
call xit(1,"Error while writing fedsk file")

 end
!-------------------------------------------------------------------------------
subroutine VisLim (MAT,PRPMAT,PRPCRP,STN,NSTR)

USE MATERIALSMODULE, only: LMAT, &
                           matdatactx, &
                           MaterialIDOfElement
USE MODELDEFINITION
use meshdatamodule,  only: meshdatactx
use constants,       only: dfmin,dlogmin,dfmax,dlogmax, half, one
use viscoLimits,     only: viscoMin, viscoMax
use iomodule
!
implicit none
!-pass
integer :: nstr, mat
double precision :: prpmat, prpcrp, stn
dimension MAT(*),PRPMAT(6,*),PRPCRP(6,*),STN(NSTR,*)
!-local
integer :: nn,n,m
double precision :: vmn, vmx, Estrs, anpwr, amhu
double precision, external :: effstn
double precision :: emhu

if (viscoMax.lt.0d0) then
do nn=1,NUMEL
    n = nn
    m = MaterialIDOfElement(n)
    if (IVLIM.eq.1) then
        VMN = PRPCRP(4,m)
        VMX = PRPCRP(5,m)
    else if (IVLIM.eq.2) then
        VMN = viscoMin
        VMX = viscoMax
        if (VMN.lt.PRPCRP(4,m)) VMN = PRPCRP(4,m)
        if (VMX.gt.PRPCRP(5,m)) VMX = PRPCRP(5,m)
    endif
    anpwr = PRPMAT(4,m)
    if (anpwr.ne.ONE) then
!       effective stress from deviatoric incompressible tensor
        Estrs = EFFSTN(STN(1,n),HALF,.false.)
        if (Estrs.gt.dfmin) then
            VMN = (VMN*Estrs**(anpwr-ONE))**(ONE/anpwr)
            VMX = (VMX*Estrs**(anpwr-ONE))**(ONE/anpwr)
        else
            VMN = VMN**(ONE/anpwr)
            VMX = VMX**(ONE/anpwr)
        endif
        emhu = matdatactx%PRPMAT(3,m)
        emhu = SQRT(emhu*emhu+VMN*VMN)
        matdatactx%PRPMAT(3,m) = MIN(emhu,VMX)
    else
        emhu = matdatactx%PRPMAT(3,m)
        emhu = SQRT(emhu*emhu+VMN*VMN)
        matdatactx%PRPMAT(3,m) = MIN(emhu,VMX)
    endif
enddo
else
   write(stderr,1) n
 1 format(1x,'no viscosity limits set while IVLIM ne 0')
   call stoper()
   write(*,*) ''
endif
!
return
end

!-------------------------------------------------------------------------------

double precision function Finc(X,Y)

 use constants, only: dfmin,dlogmin,dfmax,dlogmax
use iomodule

 implicit none

double precision :: X,Y
double precision :: a



 if (X.lt.dfmin) then
write(stderr,1) X
 1        format(1x,'Finc: EMHU < dfmin error: EMHU = ',1PG15.6)
call stoper()
 endif

 a = Y/X

 if (a.lt.dfmin) then
write(stderr,2) X,Y
 2       format(1x,'Finc: Y/X < dfmin error: X,Y = ',1PG15.6,1X, &
     G15.6)
call stoper()
 endif

 a = LOG10(a) / 2d0

 if (a.gt.dlogmax) then
write(stderr,3) X,Y
 3       format(1x,'Finc: LOG10(Y/X)/2 > dlogmax error'/ &
     1x,'Finc: X,Y = ',1PG15.6,1X,G15.6)
call stoper()
 endif

 a = X * 1D1**a

 if (a.gt.dfmax) then
write(stderr,4) X,Y
 4       format(1x,'Finc: X*10**(LOG10(Y/X)/2) > dlogmax error'/ &
     1x,'Finc: X,Y = ',1PG15.6,1X,G15.6)
call stoper()
 endif

 Finc = a

 return

 end
!-----------------------------------------------------------------------
function visco (Q,V,A,N,S,TK,P,LINEAR)

use constants, only: dfmin,dlogmin,dfmax,dlogmax, &
                     R, &
                     one, two, &
                     TenLogE

implicit none
!-pass
logical LINEAR
double precision :: Q,V,A,N,S,TK,P

!-local
double precision :: f1, f2, f3, f4, pow, visco

if (LINEAR) then

    if (LOG10(TWO*A).le.dlogmin) then
        visco = dfmax
        return
    else if (LOG10(TWO*A).ge.dlogmax) then
        visco = dfmin
        return
    endif

    f1 = ONE/(TWO*A)
    f2 = (Q+ABS(P)*V)/(R*TK)

    if (f2.le.dlogmin/TENLOGE) then
        visco = dfmin
        return
    else if (f2.ge.dlogmax/TENLOGE) then
        visco = dfmax
        return
    endif

    f2 = EXP(f2)
    pow = ONE-N

    if (pow.ge.-dfmin) then
        f3 = 1d0
    else
        if (S.gt.dfmin) then
            if (pow*LOG10(S).le.dlogmin) then
!                high stress
                visco = dfmin
                return
            endif
            if (pow*LOG10(S).ge.dlogmax) then
!                low stress
                visco = dfmax
                return
            endif
        else
            visco = dfmax
            return
        endif
        f3 = S**pow
    endif
!
    if (LOG10(f1)+LOG10(f2).le.dlogmin) then
        visco = dfmin
        return
    else if (LOG10(f1)+LOG10(f2).ge.dlogmax) then
        visco = dfmax
        return
    endif

    f4 = f1*f2

    if (LOG10(f4)+LOG10(f3).le.dlogmin) then
        visco = 0d0
        return
    else if (LOG10(f4)+LOG10(f3).ge.dlogmax) then
        visco = dfmax
        return
    endif

    visco = f4*f3

else
    f1 = (TWO*A)**(-ONE/N)
    f2 = (Q+ABS(P)*V)/(N*R*TK)
    if (f2.le.dlogmin/TENLOGE) then
        visco = dfmin
        return
    else if (f2.ge.dlogmax/TENLOGE) then
        visco = dfmax
        return
    endif
    f2 = EXP(f2)
!
    if (LOG10(f1)+LOG10(f2).le.dlogmin) then
        visco = dfmin
        return
    else if (LOG10(f1)+LOG10(f2).ge.dlogmax) then
        visco = dfmax
        return
    endif
    visco = f1*f2
endif
!
return
end
