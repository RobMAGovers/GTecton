   MODULE MODELCTX
!     this module contains global information related to the model, such
!     as version number, title, rank, number of cpus

#define PETSC_AVOID_DECLARATIONS
#include "petsc/finclude/petsc.h"
#undef PETSC_AVOID_DECLARATIONS

  character(len=1024) TITLE 
  character(len=1024) TTITLE
  integer FILE_outputf
  logical parout
  character(len=1024) part_name

  logical hasthermal
  logical hasdiffthermal

  PetscMPIInt  :: myrank, mysize
  PetscErrorCode :: pierr

  PRIVATE :: myrank, mysize
  PUBLIC :: getrank, getsize
  contains

  subroutine initializemodelctx()

    use petscksp, only: PETSC_COMM_WORLD


    implicit none

#ifdef SPARSE
      call MPI_Comm_rank (PETSC_COMM_WORLD,myrank,pierr)
      call MPI_Comm_size (PETSC_COMM_WORLD,mysize,pierr)
#endif
      parout         = .false.
      hasthermal     = .false.
      hasdiffthermal = .false.

  end subroutine initializemodelctx

  PetscMPIInt function getrank()
      implicit none
      getrank = myrank
  end function getrank

  PetscMPIInt function getsize()
      implicit none
      getsize = mysize
  end function getsize

  subroutine setrank(setrnk)
      implicit none
      integer setrnk
      myrank = setrnk
  end subroutine setrank

  subroutine setsize(setsze)
      implicit none
      integer setsze
      mysize = setsze
  end subroutine setsize


   END MODULE MODELCTX
