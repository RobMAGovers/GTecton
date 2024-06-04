module postProcessData

implicit none

double precision  :: deffac
logical           :: setria,triang,deformed,CST,opn

logical           :: GMT, connout

end module

!-----------------------------------------------------------------------
module monitorVars

implicit none

integer :: MAXMON  ! set in PltSrc/pltdef2D or Plt3D/pltdef3D
logical :: elmon,npmon,fmon, sidemon
integer :: nnp, nnel, nncol, nnrow, nsid
integer, allocatable :: monel(:), monnp(:), monf(:), monsid(:,:)


contains

subroutine alloc_elem_monitor_data()
  use modeldefinition, only: NUMEL
  implicit none
  if (.not. allocated(monel)) then
      allocate(monel(NUMEL))
  endif
  return
end subroutine alloc_elem_monitor_data

subroutine dealloc_elem_monitor_data()
  implicit none
  if (allocated(monel)) then
      deallocate(monel)
  endif
  return
end subroutine dealloc_elem_monitor_data

subroutine alloc_node_monitor_data()
  use modeldefinition, only: NUMNP
  implicit none
  if (.not. allocated(monnp)) then
      allocate(monnp(NUMNP))
  endif
  return
end subroutine alloc_node_monitor_data

subroutine dealloc_node_monitor_data()
  implicit none
  if (allocated(monnp)) then
      deallocate(monnp)
  endif
  return
end subroutine dealloc_node_monitor_data

! RG-2/8/2023 subroutine resizeMonitorData()
! RG-2/8/2023 
! RG-2/8/2023 implicit none
! RG-2/8/2023 
! RG-2/8/2023 if (allocated(monel)) then
    ! RG-2/8/2023 deallocate(monel)
! RG-2/8/2023 endif
! RG-2/8/2023 if (allocated(monnp)) then
    ! RG-2/8/2023 deallocate(monnp)
! RG-2/8/2023 endif
! RG-2/8/2023 if (allocated(monf)) then
    ! RG-2/8/2023 deallocate(monf)
! RG-2/8/2023 endif
! RG-2/8/2023 if (allocated(monsid)) then
    ! RG-2/8/2023 deallocate(monsid)
! RG-2/8/2023 endif
! RG-2/8/2023 
! RG-2/8/2023 allocate(monsid(2,maxMon))
! RG-2/8/2023 allocate(monel(maxMon))
! RG-2/8/2023 allocate(monnp(maxMon))
! RG-2/8/2023 allocate(monf(maxMon))
! RG-2/8/2023 
! RG-2/8/2023 end subroutine

end module

!-----------------------------------------------------------------------

module label

implicit none

!todo: only labelS is used... throw away the rest?

character(len=4)  :: LABELC(3),LABELD(6),LABELS(6),LABELR(6)
character(len=10) :: MATLBL(6),CRPLBL(5),TEMLBL(6)

end module

!-----------------------------------------------------------------------

module plotControl

! used to be pltcom.i

implicit none

! nodal points, used to be /NPVAR/
logical :: npnum
integer :: inpbas, inpstp
logical :: npbased

! elements data, used to be /ELMVAR/
logical :: grnum
integer :: ielbas, ielstp

! stress data
logical          :: difstn
double precision :: Sscale,Ssb

! strain data
logical          :: sdev
double precision :: Escale,Esb

! flux data
double precision :: Fscale,Fsb

! strain rate
logical          :: ratlog
double precision :: Rscale,Rsb
integer          :: irate

! displacement
double precision :: Dscale, Dsb
integer          :: iDisp

type displacementProperties
integer          :: dispOrVelo
integer          :: dispType
integer          :: coord
end type

type(displacementProperties) :: printThisDisp


! velocity
double precision :: Velsc,Vsb

! viscosity
logical          :: viscnu,viscli,vislog,visraw
double precision :: Vscale

! thermals... what?  todo
logical          :: BritDt
double precision :: idtfl

! force
integer          :: nstn

! temperature
logical          :: T0calc
integer          :: itemp, irelt

! opn
integer          :: iopnflag, NUMZ
double precision :: opnscale, opnaxnu

! gravity
! for all
integer          :: icompens
double precision :: Ygref
! for 2D
double precision :: g0,rhoa
! for 3d
! detpot  = determine (gravity) potential
! detgrav = determine (gravity) acceleration
logical          :: detpot,detgrav,detgrad,detdilatation, OpGravChange
character(len=4) :: coordtype
integer          :: iEffD



! surface var
integer          :: isvsrf

! plot dimensions
integer          :: NPLT,NPLDOF

! generic 
logical          :: readsy, sbhor, sysc
double precision :: hatangle,hatsize,sysi,iarcont
integer          :: mattsym

! projection
logical          :: project
double precision :: gplon,gplat,gpaz

! strings containing the units of various physical quantities
character(len=10000) :: Sunit,Eunit,Funit,Runit,Dunit,Vunit


logical          :: shrink
integer          :: jshell
double precision :: shrinkFactor

integer :: MEQ
logical :: showfbc

double precision ::  dispmn(2),dispmx(2),strsmn(2),strsmx(2), &
   ratemn(2),ratemx(2),tempmn,tempmx,fluxmn(2),fluxmx(2), &
   crdmin(3),crdmax(3)

integer, parameter :: MAXPND=80
integer            :: ixpnd(2,MAXPND), iypnd(2,MAXPND), nxpnd,nypnd

logical            :: stndev

end module

module symbols

integer, parameter :: MAXSYM=100

double precision :: SYFIX(2,MAXSYM),SYSLID(2,MAXSYM),SYWINK(2,MAXSYM),&
SYSTRS(2,MAXSYM),SYVELO(2,MAXSYM),SYFORC(2,MAXSYM),&
  SYSHEAR(2,MAXSYM),SYOPEN(2,MAXSYM),SYSHR(2,MAXSYM)

integer :: IBcFix,IBcSld,IBcWnk,IBcStn,IBcVel,IBcFor,IBcSh,&
  IBcOpn,IBcShr,MFIX,MSLID,MWINK,MSTRS,MVELO,MFORC,MSHEAR,MOPEN,&
  MSHR

end module
