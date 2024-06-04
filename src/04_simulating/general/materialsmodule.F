   MODULE MATERIALSMODULE
#define PETSC_AVOID_DECLARATIONS
#include "petsc/finclude/petsc.h"
#undef PETSC_AVOID_DECLARATIONS
!     module containing all mechanical and thermal material properties
!     replaces common blocks EPOINT (only mat pointers), PTPROP, and VISDAT (only LMAT)
!     also contains sediment related properties and spherical geometry
!     properites

   logical :: LMAT
   PetscReal, ALLOCATABLE :: PRPSED(:)
!   PetscReal              :: RADIUS    ! radius of planet, in case ishell = 1
   PetscReal, ALLOCATABLE :: GRAV(:)
!      PetscInt :: NUMAT


logical :: effsign,strinc


   TYPE matdata
!#define PETSC_AVOID_DECLARATIONS
!#include "petsc/finclude/petsc.h"
!#undef PETSC_AVOID_DECLARATIONS
  SEQUENCE
  logical :: hascrp,hasdif,hasplas
  PetscReal, ALLOCATABLE :: PRPMAT(:,:)
  PetscReal, ALLOCATABLE :: PRPCRP(:,:)
  PetscReal, ALLOCATABLE :: PRPDIF(:,:)
  PetscReal, ALLOCATABLE :: PRPLAS(:,:)
  PetscReal, ALLOCATABLE :: PRPTEM(:,:)
#ifdef EXPERIMENTALTRANSIENTDISCREEP
  ! material parameters for transient dislocation creep
  PetscReal, ALLOCATABLE :: PRPDISLC(:,:)
  ! material parameters for steady state diffusion creep
  PetscReal, ALLOCATABLE :: PRPDIFFC(:,:)
#endif

   END TYPE matdata

type materialFlags
!logical 


end type materialFlags



   type (matdata) :: matdatactx

   contains


!********************************************************************
function MaterialIDOfElement(elementID)

use meshdatamodule, only: meshdatactx

implicit none

integer :: MaterialIDOfElement, elementID

if (elementID.le.0) then
    write(*,*) "MaterialIDOfElement should be >0 but is", elementID
    stop "You may have found a bug. Please contact model support. Leaving GTecton"
else if (elementID.gt.size(meshdatactx%mat,1)) then
    write(*,*) "MaterialIDOfElement requests material ID of elem", elementID
    write(*,*) "but the mesh has only ", size(meshdatactx%mat,1), "elements"
    stop "You may have found a bug. Please contact model support. Leaving GTecton"
endif


if (LMAT) then
    MaterialIDOfElement = elementID
else
    MaterialIDOfElement = meshdatactx%mat(elementID)
endif

end function
!********************************************************************






!*********************************************************************
   subroutine allocatematdatamech(numats,ndof,docrp,dodif,doplas)
   use debugmodule, only: AllocateError
!*********************************************************************
   PetscInt numats, ndof, ERROR
   logical docrp,dodif,doplas
!*********************************************************************
if (allocated(matdatactx%PRPMAT)) then
    deallocate(matdatactx%PRPMAT)
endif
   ALLOCATE(matdatactx%PRPMAT(6,numats), STAT=ERROR)

   if (docrp) then
if (allocated(matdatactx%PRPCRP)) then
    deallocate(matdatactx%PRPCRP)
endif
  ALLOCATE(matdatactx%PRPCRP(6,numats), STAT=ERROR)
  call AllocateError("PRPCRP", ERROR)
  matdatactx%hascrp = docrp
   endif

   if (dodif) then
if (allocated(matdatactx%PRPDIF)) then
    deallocate(matdatactx%PRPDIF)
endif
  ALLOCATE(matdatactx%PRPDIF(9,numats), STAT=ERROR)
  call AllocateError("PRPDIF", ERROR)
  matdatactx%hasdif = dodif
   endif

if (doplas) then
    if (allocated(matdatactx%PRPLAS)) then
        deallocate(matdatactx%PRPLAS)
    endif
    ALLOCATE(matdatactx%PRPLAS(9,numats), STAT=ERROR)
    call AllocateError("PRPLAS", ERROR)
    matdatactx%hasplas = doplas
endif

if (allocated(PRPSED)) then
    deallocate(PRPSED)
endif
   ALLOCATE(PRPSED(2), STAT=ERROR)
   call AllocateError("PRPSED", ERROR)

if (allocated(GRAV)) then
    deallocate(GRAV)
endif
   ALLOCATE(GRAV(ndof), STAT=ERROR)
   call AllocateError("Grav", ERROR)

#ifdef EXPERIMENTALTRANSIENTDISCREEP
  if (allocated(matdatactx%PRPDISLC)) then
    deallocate(matdatactx%PRPDISLC)
  endif
  ALLOCATE(matdatactx%PRPDISLC(10,numats), STAT=ERROR)
  if (allocated(matdatactx%PRPDIFFC)) then
    deallocate(matdatactx%PRPDIFFC)
  endif
  ALLOCATE(matdatactx%PRPDIFFC(5,numats), STAT=ERROR)

#endif

   end subroutine allocatematdatamech
!*********************************************************************

!*********************************************************************
   subroutine allocatematdatatemp(numats)
   PetscInt numats, ERROR

if (allocated(matdatactx%PRPTEM)) then
    deallocate(matdatactx%PRPTEM)
endif
   ALLOCATE(matdatactx%PRPTEM(6,numats), STAT=ERROR)

   end subroutine allocatematdatatemp
!*********************************************************************

!*********************************************************************
   subroutine deallocatematdatamech()
   PetscInt ERROR

if (allocated(matdatactx%PRPMAT)) then
   DEALLOCATE(matdatactx%PRPMAT, STAT=ERROR)
endif

if (matdatactx%hascrp) then
    if (allocated(matdatactx%PRPCRP)) then
        DEALLOCATE(matdatactx%PRPCRP, STAT=ERROR)
    endif
endif

if (matdatactx%hasdif) then
    if (allocated(matdatactx%PRPDIF)) then
        DEALLOCATE(matdatactx%PRPDIF, STAT=ERROR)
    endif
endif

if (matdatactx%hasplas) then
    if (allocated(matdatactx%PRPLAS)) then
        DEALLOCATE(matdatactx%PRPLAS, STAT=ERROR)
    endif
endif


if (allocated(PRPSED)) then
   DEALLOCATE(PRPSED, STAT=ERROR)
endif

if (allocated(GRAV)) then
   DEALLOCATE(GRAV, STAT=ERROR)
endif


#ifdef EXPERIMENTALTRANSIENTDISCREEP
if (allocated(matdatactx%PRPDISLC)) then
   DEALLOCATE(matdatactx%PRPDISLC, STAT=ERROR)
endif

if (allocated(matdatactx%PRPDIFFC)) then
   DEALLOCATE(matdatactx%PRPDIFFC, STAT=ERROR)
endif
#endif

   end subroutine deallocatematdatamech
!*********************************************************************

!*********************************************************************
   subroutine deallocatematdatatemp()
   PetscInt ERROR

if (allocated(matdatactx%PRPTEM)) then
   DEALLOCATE(matdatactx%PRPTEM, STAT=ERROR)
endif

   end subroutine deallocatematdatatemp
!*********************************************************************

!*********************************************************************
   subroutine materialgen (PROP,NDIM,NPROP,NUMAT,IFLAG)

!        call materialgen (matdatactx%PRPLAS,9,7,NUMAT,1)

!*********************************************************************
! Routine to read and generate material data
!
! This routine can read both mechanical element data (typically in TECIN.DAT)
! as well as temperature data (typically in TEMPIN.DAT)
!
!       NPROP   = number of material properties
!       PROP    = material properties array
!                  can be one of:
!                  matdatactx%PRPMAT
!                  matdatactx%PRPCRP
!                  matdatactx%PRPLAS
!                  matdatactx%PRPDIF
!                  matdatactx%PRPDISLC
!                  matdatactx%PRPDIFFC
!                  matdatactx%PRPTEM
!                 depending on the caller
!       NUMAT   = number of material properties
!       IFLAG   = 0,1: first two items float or integer
!       ng      = element generation increment
!*********************************************************************
use modelctx,    only: getrank
use debugmodule, only: iecho, xit, debug
use iomodule,    only: luin, stderr

implicit none

integer MAXPRP
#ifdef EXPERIMENTALTRANSIENTDISCREEP
parameter (MAXPRP=10)
#else
parameter (MAXPRP=7)
#endif
!-pass
integer NDIM,NPROP,NUMAT,IFLAG
double precision,  dimension(NDIM,*), intent(inout) :: PROP

!-locl
integer if1, if2, nn, mmat, nmat, ngen, mgen, numgen, ios,i
character(len=256) record
character(len=23) frm
double precision f(MAXPRP)
integer lnblk
external lnblk
!*********************************************************************
!-init
mgen = 0
mmat = 0
numgen = 0



if (NPROP.lt.1 .or. NPROP.gt.MAXPRP) then
    write(stderr,1)
 1       format(1x,'MATGEN: NPROP out of range')
    call xit(1," ")
endif

if (IFLAG.eq.0) then
    write(frm,2,err=700) MAXPRP
 2       format('(I12,',I2,'(G14.0))')
else
    write(frm,3,err=700) MAXPRP-2
 3       format('(I12,2I5,',I2,'(G14.0))')
endif

if (iecho.eq.2) then
     write(stderr,*) ' *** read material data ***'
     write(stderr,*) ' nread = ',NUMAT, ' using format ',frm
endif

if (NUMAT.le.0) then
    write(stderr,4)
 4       format(1x,'MATGEN: no material properties ?')
    call xit(1," ")
endif

!##### read the material lines ##########

100     read(luin,"(a)",err=700,end=710,iostat=ios) record
nn = lnblk(record)

if (debug) then
    write(*,*) 'material gen read record: ', record
endif

if (record(1:3).eq.'end' .or. nn.le.0) then
    if (debug) then
        write(*,*) "End detected. Finished reading block of material data."
    endif
    goto 600
endif

if (record(12:12).eq."") then
    write(*,*) "Input reader says: Something wrong with material number", numgen+1
    write(*,*) "position 12, where the material number should be, is empty."
endif

if (IFLAG.eq.0) then
!            write(*,*) 'iflag 0 (mechanical?)', NPROP, frm
    read (record,fmt=frm,err=700,end=710,iostat=ios) nmat, (f(i),i=1,NPROP)
    ! those properties are:
!            E(n)         2. Elastic Young's modulus
!            POIS(n)      4. Elastic Poisson's ratio
!            EMHU(n)      5. Effective viscosity
!            ANPWR(n)     6. Power
!            WT(n)        7. Mass density (if MODE>=3 interpreted
!                            as zero temperature mass density)
!            TH(n)



    if (debug) then
        write(*,*) "interpreting with flag 0"
        write(*,*) "interpreting with format: ", frm
        write(*,*) "interpreted: ", nmat, (f(i),i=1,NPROP)
    endif
else
!            write(*,*) 'iflag not 0 (thermal?)'
    read (record,fmt=frm,err=700,end=710,iostat=ios) nmat, if1,if2,f(2),f(3),f(4),f(5),f(7)
    f(1) = DBLE(if1)
    f(6) = DBLE(if2)
    if (debug) then
        write(*,*) "interpreting with flag 1"
        write(*,*) 'interpreted: ', nmat, if1,if2,f(2),f(3),f(4),f(5),f(7)
    endif
endif

if (nmat.eq.0) then
    goto 600
endif

 8               format(1x,'matgen: no of material inputs exceeds ', I12)


numgen = numgen + 1

if (numgen.gt.NUMAT) then
    write(stderr,8) NUMAT
    call xit(1," ")
endif

do i=1,NPROP
    PROP(i,nmat) = f(i)
enddo


mmat = nmat
mgen = ngen

goto 100
!####################



600     if (numgen.ne.NUMAT) then
    write(stderr,9) numgen, NUMAT
 9       format(1x,'materialgen: no of generated materials ',I12, ' not equal to NUMAT in TECIN.DAT ', I12)
    call xit(1," ")
endif

return

700     write(stderr,10) ios
10   format(1x,'materialgen: fatal input error ',i8)
call prerr('materialgen')
call xit(1," ")

710     write(stderr,20) ios
20   format(1x,'materialgen: reached premature end ',i8)
call prerr('materialgen')
call xit(1," ")


end subroutine materialgen
!*********************************************************************

END MODULE MATERIALSMODULE

module viscoLimits

implicit none

double precision :: viscoMin, viscoMax

end module
