#ifdef SPARSE
module spetscmodule
!#else
!module petscmodule
!#endif


use fileModule,  only: myTmpDir
!use debugmodule, only: iecho, &
!                       iAmAtTimestep
use petscksp

! this module contains all the PETSc interfacing.
! All of the subroutine names begin with 'PETSc_', for clarity
! The surname equals the name of the PETSc call.

! The idea of this module is to put all the PETSc stuff
! in a single module, so that PETSc version upgrades can
! be handled more efficiently and transparently.
! Also, it is clearly possible to distinguish between
! PETSc calls and internal calls.

! More over, it facilitates logging of the traffic
! between PETSc and GTecton.

! On a side note, it contains the variables that
! used to be in the included files:
! petsc_init.i  petsc_m.i  petsc_t.i

#include "petsc/finclude/petscmat.h" 
! for the Mat type

#include "petsc/finclude/petscvec.h" 
! for the Vec type

#include "petsc/finclude/petscsys.h" 

#include "petsc/finclude/petscksp.h" 
! for the KSP type

implicit none

logical, parameter :: logpetsc     = .false.
integer, parameter :: petsclogfile = 74

! Mechanical variables
Mat :: M_A
Vec :: M_B,M_X
Vec :: M_Bseq,M_Xseq
KSP :: M_SLES

! Thermal variables
Mat :: T_A
Vec :: T_B,T_X
Vec :: T_Bseq,T_Xseq
KSP :: T_SLES


! positions of non zeros.
! not sure of still necessary with newer PETSc
!integer, allocatable :: mechanicalNonZeros(:)
!integer, allocatable ::    thermalNonZeros(:)






! from PETSc_init
logical :: PetSc_Init
logical :: M_Init
logical :: T_Init


! logging
public :: PETSc_OpenLogFile 
public :: PETSc_CloseLogFile
! general
public :: PETSc_PetscInitialize
public :: PETSc_PetscFinalize
! vectors
public :: PETSc_VecCreate
public :: PETSc_VecCreateMPI
public :: PETSc_VecCreateSEQ
public :: PETSc_VecDuplicate
public :: PETSc_VecSetSizes
public :: PETSc_VecSetFromOptions
public :: PETSc_VecAssemblyBegin
public :: PETSc_VecAssemblyEnd
public :: PETSc_VecGetArrayF90
public :: PETSc_VecRestoreArrayF90
public :: PETSc_VecDestroy
! matrices
public :: PETSc_MatCreate
public :: PETSc_MatDuplicate
public :: PETSc_MatSetSizes
public :: PETSc_MatMPIAIJSetPreallocation
public :: PETSc_MatSeqAIJSetPreallocation
public :: PETSc_MatSetType
public :: PETSc_MatSetOption
public :: PETSc_MatSetFromOptions
public :: PETSc_MatSetValuesLocal
public :: PETSc_MatGetOwnershipRange
public :: PETSc_MatAssemblyBegin
public :: PETSc_MatAssemblyEnd
public :: PETSc_MatDestroy
! mappings (IS = Index Set)
public :: PETSc_ISLocalToGlobalMappingCreate
public :: PETSc_ISGetIndicesF90
public :: PETSc_ISRestoreIndicesF90
public :: PETSc_ISCreateBlock
public :: PETSc_ISCreateStride
public :: PETSc_MatSetLocalToGlobalMapping
public :: PETSc_VecSetLocalToGlobalMapping
public :: PETSc_ISLocalToGlobalMappingDestroy
public :: PETSc_ISDestroy
! Orderings
public :: PETSc_AOCreateBasic
public :: PETSc_AODestroy
public :: PETSc_AOApplicationToPetsc
! Scattering
public :: PETSc_VecScatterCreate
public :: PETSc_VecScatterCreateToAll
public :: PETSc_VecScatterBegin
public :: PETSc_VecScatterEnd
public :: PETSc_VecScatterDestroy
! solving
public :: PETSc_KSPCreate
public :: PETSc_KSPSetTolerances
public :: PETSc_KSPSetFromOptions
public :: PETSc_KSPSetOperators
public :: PETSc_KSPSolve
public :: PETSc_KSPGetConvergedReason
public :: PETSc_KSPGetIterationNumber
public :: PETSc_KSPDestroy
!Preconditioning
public :: PETSc_KSPGetPC
public :: PETSc_PCSetType
public :: PETSc_PCFactorSetMatSolverPackage
!Miscellaneous
public :: PETSC_error ! to handle non-zero returns from PETSc calls

#ifdef SPARSE
public :: UserFinalize
#endif

contains 

!-------------------------------------------------------------------------

#ifdef SPARSE
subroutine UserFinalize ()

use petscksp

implicit none


!-locl
PetscErrorCode ierr

if (.not.PetSc_Init) then
    return
endif

call PetscFinalize(ierr)

call MPI_FINALIZE(ierr)

return
end subroutine

#endif


!-------------------------------------------------------------------------
! logging subs

subroutine PETSc_openlogfile()
 use filemodule,      only: workpath

 implicit none

!    write(*,*) 'opening petsc_interface.log in', trim(workpath)

 if (logpetsc) then
open(unit=petsclogfile, file=trim(myTmpDir)//'/'//'petsc_interface.log')
 endif

 write(petsclogfile,*) 'opened file; Start logging PETSc commands'
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_closelogfile()

 implicit none

 write(petsclogfile,*) 'Finished logging PETSc commands; closing file'

 if (logpetsc) then
close(unit=petsclogfile)
 endif

end subroutine

!-------------------------------------------------------------------------

!general
subroutine PETSc_PetscInitialize
implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling my PetscInitialize with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_PetscFinalize
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling PetscFinalize with'
endif
end subroutine 

!-------------------------------------------------------------------------


! vectors
subroutine PETSc_VecCreate(vectorHandle)
    implicit none

!#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscvec.h"

    Vec     :: vectorHandle
!    integer :: ierr

!    write(petsclogfile,*) 'calling VecCreate with'

!    call VecCreate(MPI_COMM_WORLD, vectorHandle, ierr)

!    write(petsclogfile,*) 'VecCreate returned error ', ierr

end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecCreateMPI
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecCreateSEQ
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecDuplicate
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecSetSizes
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecSetFromOptions
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecAssemblyBegin
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecAssemblyEnd
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecGetArrayF90
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecRestoreArrayF90
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecDestroy
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

! matrices
subroutine PETSc_MatCreate
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatDuplicate
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatSetSizes
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatMPIAIJSetPreallocation
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling MatMPIAIJSetPreallocation with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatSeqAIJSetPreallocation
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling MatSeqAIJSetPreallocation with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatSetType
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatSetOption
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatSetFromOptions
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatSetValuesLocal(matrix,k,m,v,callerID,globalElemID,SLcol,SLrow, ierr)

use constants,      only: dfmin, dfmax
!use debugmodule,    only: iecho
use meshdatamodule, only: meshdatactx
use modelctx,       only: getrank
!#ifdef PLOTNONZEROES
!use debugmodule,    only: debug, matrixEntryNonZero
!#else
!use debugmodule,    only: debug
!#endif
use modeldefinition,only: NEQGlobal
!#ifdef SPARSE
!use AOmodule,       only: equations
!#endif
use petscksp
use iomodule

! This function is an interface between PETSC and and GTecton.
! This inserts values into the stiffness matrix
implicit none


! for PETSC_VIEWER_ASCII_INFO
#include "petsc/finclude/petscviewer.h"


integer          :: k,m, ierr  !  k and m extracted from LM by caller.
integer          :: SLcol,SLrow, globalElemID, callerID ! used to identify the call.
                                                  ! Every combination of these values
                                                  ! should add only a single
                                                  ! value to the stiffness matrix.
                                                  ! During debugging, it is possible to sort
                                                  ! on these values, for examples, to see omissions.
                                                  ! SL is local stiffness matrix entry.
integer          :: entryID
double precision :: v
Mat              :: matrix ! should either correspond to M_A of T_A (mech or temp)


!#ifdef SPARSE
if (getrank().eq.0) then
!if (iAmAtTimestep.eq.11) then
!if (callerID.eq.1 .and. globalElemID.eq.1) then
!if ((k.eq.8 .and. m.eq.7) .or. (k.eq.7 .and. m.eq.8))  then
    entryID = callerID * 1000000 + globalElemID * 10000 + SLcol * 100 + SLrow
    ! positioning start counting at 1, but arguments to MatSetValues starts counting at 0
!    write(*,*) "rank",getrank(),"entryID ", entryID, "is value ", v, " to loc pos: ", k, m
endif
!endif
!#endif



! shooting values into the PETSc stiffness matrix.
! Notice the 'local' in the call 'MatSetValuesLocal' below.
! The k and m are not the global positions, and they are not 
! the global equation numbers, but the local numbers.

! The relationship between the two is set up in
! subroutine createlocaltoglobalmapping.
! For reference, it is possible to find the global 
! equation number by finding equations(m) (or k).



if (v.ne.v) then
    write(*,*) "GTECTON tried to insert a NaN in position", k, m, &
               "of the stiffness matrix. Leaving GTecton..."
    stop "Leaving GTecton..."
else if (v.gt.dfmax) then
    write(*,*) "GTECTON tried to insert a +infinity in position", k, m, &
               "of the stiffness matrix. Leaving GTecton..."
    stop "Leaving GTecton..."
else if (v.lt.-dfmax) then
    write(*,*) "GTECTON tried to insert a -infinity in position", k, m, &
         "of the stiffness matrix. Leaving GTecton..."
    stop "Leaving GTecton..."
endif


call MatSetValuesLocal(matrix,1,k,1,m,v,ADD_VALUES,ierr)


#ifdef PLOTNONZEROES
#ifdef SPARSE
matrixEntryNonZero(k+1,m+1) = .true.
#endif
#endif


if (ierr.ne.0) then
    write(stderr,*) 'PETSc_MatSetValuesLocal, called by', callerID
    write(stderr,*) 'MatSetValuesLocal error ', ierr
    write(stderr,*) 'writing value ',v,' to pos ', k,m, ' in mat ', M_A
    write(stderr,*) 'please contact model support'
    !write(stderr,*) 'PETSC_VIEWER_ASCII_INFO ', PETSC_VIEWER_ASCII_INFO
    
!#    call matview(M_A,PETSC_VIEWER_ASCII_INFO,ierr)
    !call stoper()
endif

end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatGetOwnershipRange
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling MatGetOwnershipRange with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatAssemblyBegin
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatAssemblyEnd
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatDestroy
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

! mappings (IS = Index Set)
subroutine PETSc_ISLocalToGlobalMappingCreate
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_ISGetIndicesF90
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_ISRestoreIndicesF90
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_ISCreateBlock
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_ISCreateStride
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_MatSetLocalToGlobalMapping
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecSetLocalToGlobalMapping
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_ISLocalToGlobalMappingDestroy
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_ISDestroy
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

! Application orderings.
subroutine PETSc_AOCreateBasic
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling AOCreateBasic with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_AOApplicationToPetsc
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling AOApplicationToPetsc with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_AODestroy
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling AODestroy with'
endif
end subroutine

!-------------------------------------------------------------------------

! Scattering
subroutine PETSc_VecScatterCreate
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecScatterCreateToAll
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecScatterBegin
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecScatterEnd
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_VecScatterDestroy
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

! solving
subroutine PETSc_KSPCreate
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_KSPSetTolerances
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling KSPSetTolerances with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_KSPSetFromOptions
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling KSPSetFromOptions with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_KSPSetOperators
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling KSPSetOperators with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_KSPSolve
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_KSPGetConvergedReason
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling KSPGetConvergedReason with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_KSPGetIterationNumber
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling KSPGetIterationNumber with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_KSPDestroy
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling with'
endif
end subroutine

!-------------------------------------------------------------------------

! preconditioning

subroutine PETSc_KSPGetPC
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling KSPGetPC with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_PCSetType
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling PCSetType with'
endif
end subroutine

!-------------------------------------------------------------------------

subroutine PETSc_PCFactorSetMatSolverPackage
 implicit none
if (logpetsc) then
 write(petsclogfile,*) 'calling PCFactorSetMatSolverPackage with'
endif
end subroutine

!---------------------------------------------------------------
subroutine PETSC_error(callingSub, PETScCall, ierror)

    use modelctx,    only: getrank

    implicit none

    integer          :: ierror
    character(len=*) :: callingSub
    character(len=*) :: PETScCall

    if (ierror.eq.0) then
!        write(*,*)        "caller ", trim(callingSub), &
!                        " on rank ", getrank(), &
!                         " called ", trim(PETScCall), " and did great!"
    else
        write(*,*)        "caller ", trim(callingSub), &
                         " called ", trim(PETScCall),  &
        " which failed with error ", ierror
        stop "Leaving GTECTON..."
    endif


end subroutine

!#ifdef SPARSE
end module spetscmodule
!#else
!end module petscmodule
#endif
