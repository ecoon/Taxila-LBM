!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_distribution_function.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            14:06:07 MDT
!!!     last modified:   28 March 2011
!!!       at:            16:03:02 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module LBM_Distribution_Function_module
  use petsc
  use LBM_Info_module
  use LBM_Discretization_Type_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: distribution_type
     MPI_Comm comm
     PetscInt s
     PetscInt b
     type(info_type),pointer:: info
     type(discretization_type),pointer:: disc
     Vec fi
     Vec fi_g
     PetscScalar,pointer:: fi_a(:)
     PetscScalar,pointer:: fi_eq(:)
     DM,pointer :: da
  end type distribution_type

  public:: DistributionCreate, &
       DistributionDestroy, &
       DistributionSetInfo, &
       DistributionSetDiscretization, &
       DistributionSetSizes, &
       DistributionSetDA, &
       DistributionLocalToLocal

contains
  function DistributionCreate(comm) result(distribution)
    MPI_Comm comm
    type(distribution_type),pointer:: distribution
    allocate(distribution)
    distribution%comm = comm
    distribution%s = -1
    distribution%b = -1
    nullify(distribution%info)
    nullify(distribution%disc)
    distribution%fi = 0
    distribution%fi_g = 0
    nullify(distribution%fi_a)
    nullify(distribution%fi_eq)
    nullify(distribution%da)
  end function DistributionCreate

  subroutine DistributionDestroy(distribution, ierr)
    type(distribution_type) distribution
    PetscErrorCode ierr
    
    if (distribution%fi /= 0) call VecDestroy(distribution%fi, ierr)
    if (distribution%fi_g /= 0) call VecDestroy(distribution%fi_g, ierr)
    if (associated(distribution%fi_eq)) deallocate(distribution%fi_eq)
  end subroutine DistributionDestroy

  subroutine DistributionSetInfo(distribution, info)
    type(distribution_type) distribution
    type(info_type),pointer:: info
    distribution%info => info
  end subroutine DistributionSetInfo

  subroutine DistributionSetDiscretization(distribution, disc)
    type(distribution_type) distribution
    type(discretization_type),pointer:: disc
    distribution%disc => disc
    distribution%b = disc%b
  end subroutine DistributionSetDiscretization

  subroutine DistributionSetSizes(distribution, s)
    type(distribution_type) distribution
    PetscInt s
    distribution%s = s
    allocate(distribution%fi_eq(distribution%s*(distribution%b+1)*distribution%info%gxyzl))
  end subroutine DistributionSetSizes
    
  subroutine DistributionSetDA(distribution, da)
    type(distribution_type) distribution
    DM,target:: da
    distribution%da => da
  end subroutine DistributionSetDA

  subroutine DistributionLocalToLocal(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DMDAVecRestoreArrayF90(distribution%da, distribution%fi, distribution%fi_a, ierr)
    call DMDALocalToLocalBegin(distribution%da, distribution%fi, INSERT_VALUES, &
         distribution%fi, ierr)
    call DMDALocalToLocalEnd(distribution%da, distribution%fi, INSERT_VALUES, &
         distribution%fi, ierr)
    call DMDAVecGetArrayF90(distribution%da, distribution%fi, distribution%fi_a, ierr)
  end subroutine DistributionLocalToLocal
end module LBM_Distribution_Function_module
    
