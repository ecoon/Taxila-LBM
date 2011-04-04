!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_discretization.F90
!!!     version:         
!!!     created:         14 March 2011
!!!       on:            16:33:56 MDT
!!!     last modified:   04 April 2011
!!!       at:            12:11:45 MDT
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
       DiscretizationDestroy, &
       DiscretizationSetType, &
       DiscretizationSetUp, &
       DiscretizationSetUpPhase

contains
  function DiscretizationCreate(comm) result(disc)
    MPI_Comm comm
    type(discretization_type),pointer:: disc

    allocate(disc)
    disc%comm = comm
    disc%name = NULL_DISCRETIZATION
    disc%ndims = -1
    disc%b = -1
    nullify(disc%weights)
  end function DiscretizationCreate

  subroutine DiscretizationDestroy(disc, ierr)
    type(discretization_type) disc
    PetscErrorCode ierr

    if (associated(disc%ci)) deallocate(disc%ci)
    if (associated(disc%weights)) deallocate(disc%weights)
    if (associated(disc%mt_mrt)) deallocate(disc%mt_mrt)
    if (associated(disc%mmt_mrt)) deallocate(disc%mmt_mrt)
  end subroutine DiscretizationDestroy
  
  subroutine DiscretizationSetType(disc, name)
    type(discretization_type) disc
    PetscInt name
    disc%name = name
  end subroutine DiscretizationSetType

  subroutine DiscretizationSetUp(disc)
    type(discretization_type) disc

    PetscErrorCode ierr
    select case(disc%name)
    case(D3Q19_DISCRETIZATION)
       call DiscretizationSetUp_D3Q19(disc)
    case(D2Q9_DISCRETIZATION)
       call DiscretizationSetUp_D2Q9(disc)
    case DEFAULT
       if (disc%ndims < 0 ) SETERRQ(1, 1, 'Invalid Discretization', ierr)
    end select
  end subroutine DiscretizationSetUp

  subroutine DiscretizationSetUpPhase(disc, phase)
    use LBM_Phase_module
    type(discretization_type) disc
    type(phase_type) phase
    PetscErrorCode ierr
    
    select case(disc%name)
    case(D3Q19_DISCRETIZATION)
       call DiscretizationSetUpPhase_D3Q19(disc, phase)
    case(D2Q9_DISCRETIZATION)
       call DiscretizationSetUpPhase_D2Q9(disc, phase)
    case DEFAULT
       SETERRQ(1,1,'invalid discretization in LBM',ierr)
    end select
  end subroutine DiscretizationSetUpPhase
end module LBM_Discretization_module
