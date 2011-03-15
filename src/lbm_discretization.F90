!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_discretization.F90
!!!     version:         
!!!     created:         14 March 2011
!!!       on:            16:33:56 MDT
!!!     last modified:   15 March 2011
!!!       at:            17:40:32 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

! generic class for discretizations
module LBM_Discretization_module
  use petsc
  use LBM_Discretization_Type_module
  use LBM_Discretization_D3Q19_module
  use LBM_Discretization_D2Q9_module
  implicit none
  
  private
#include "lbm_definitions.h"

  public:: DiscretizationCreate, &
       DiscretizationDestroy

contains
  function DiscretizationCreate(comm) result(disc)
    MPI_Comm comm
    type(discretization_type),pointer:: disc

    allocate(disc)
    disc%comm = comm
    disc%name = NULL_DISCRETIZATION
    disc%ndims = -1
    disc%b = -1
    nullify(disc%ci)
    nullify(disc%weights)
  end function DiscretizationCreate

  subroutine DiscretizationDestroy(disc, ierr)
    type(discretization_type) disc
    PetscErrorCode ierr

    if (associated(disc%ci)) deallocate(disc%ci)
    if (associated(disc%weights)) deallocate(disc%weights)
  end subroutine DiscretizationDestroy

  subroutine DiscretizationSetup(disc, name)
    type(discretization_type) disc
    PetscInt name
    PetscErrorCode ierr

    select case(name)
    case(D3Q19_DISCRETIZATION)
       call DiscretizationSetup_D3Q19(disc)
    case(D2Q9_DISCRETIZATION)
       call DiscretizationSetup_D2Q9(disc)
    case DEFAULT
       if (disc%ndims < 0 ) SETERRQ(1, 1, 'Invalid Discretization', ierr)
    end select
  end subroutine DiscretizationSetup

  subroutine DiscretizationSetupConstants(disc, constants)
    use LBM_Constants_module
    type(discretization_type) disc
    type(constants_type) constants
    PetscErrorCode ierr
    
    select case(disc%name)
    case(D3Q19_DISCRETIZATION)
       call DiscretizationSetupConstants_D3Q19(disc, constants)
    case(D2Q9_DISCRETIZATION)
       call DiscretizationSetupConstants_D2Q9(disc, constants)
    case DEFAULT
       SETERRQ(1,1,'invalid discretization in LBM',ierr)
    end select
  end subroutine DiscretizationSetupConstants
end module LBM_Discretization_module
