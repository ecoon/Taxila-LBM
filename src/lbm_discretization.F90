!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_discretization.F90
!!!     version:         
!!!     created:         14 March 2011
!!!       on:            16:33:56 MDT
!!!     last modified:   13 April 2011
!!!       at:            11:17:43 MDT
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
       DiscretizationSetSizes, &
       DiscretizationSetUp, &
       DiscretizationSetUpRelax, &
       DiscretizationEquilf

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
    if (associated(disc%opposites)) deallocate(disc%opposites)
    if (associated(disc%mt_mrt)) deallocate(disc%mt_mrt)
    if (associated(disc%mmt_mrt)) deallocate(disc%mmt_mrt)
    if (associated(disc%ffw)) deallocate(disc%ffw)
  end subroutine DiscretizationDestroy
  
  subroutine DiscretizationSetSizes(disc, stencil_size)
    type(discretization_type) disc
    PetscInt stencil_size
    disc%stencil_size = stencil_size
  end subroutine DiscretizationSetSizes

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

  subroutine DiscretizationSetUpRelax(disc, relax)
    use LBM_Relaxation_module
    type(discretization_type) disc
    type(relaxation_type) relax
    PetscErrorCode ierr
    
    select case(disc%name)
    case(D3Q19_DISCRETIZATION)
       call DiscretizationSetUpRelax_D3Q19(disc, relax)
    case(D2Q9_DISCRETIZATION)
       call DiscretizationSetUpRelax_D2Q9(disc, relax)
    case DEFAULT
       SETERRQ(1,1,'invalid discretization in LBM',ierr)
    end select
  end subroutine DiscretizationSetUpRelax

  subroutine DiscretizationEquilf(disc, rho, u, walls, feq, relax, dist)
    use LBM_Distribution_Function_type_module
    use LBM_Relaxation_module
    type(discretization_type) disc
    type(distribution_type) dist
    type(relaxation_type) relax

    PetscScalar,dimension(0:dist%b,dist%info%gxyzl):: feq
    PetscScalar,dimension(dist%info%gxyzl):: rho
    PetscScalar,dimension(dist%info%ndims,dist%info%gxyzl):: u
    PetscScalar,dimension(dist%info%gxyzl):: walls

    PetscErrorCode ierr

    select case(disc%name)
    case(D3Q19_DISCRETIZATION)
       call DiscretizationEquilf_D3Q19(disc, rho, u, walls, feq, relax, dist)
    case(D2Q9_DISCRETIZATION)
       call DiscretizationEquilf_D2Q9(disc, rho, u, walls, feq, relax, dist)
    case DEFAULT 
       SETERRQ(1,1,'invalid discretization in LBM',ierr)
    end select
  end subroutine DiscretizationEquilf
end module LBM_Discretization_module
