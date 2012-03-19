!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_discretization.F90
!!!     version:
!!!     created:         14 March 2011
!!!       on:            16:33:56 MDT
!!!     last modified:   19 October 2011
!!!       at:            10:07:17 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

! generic class for discretizations
module LBM_Discretization_module
  use petsc
  use LBM_Error_module
  use LBM_Discretization_Type_module
  use LBM_Discretization_D3Q19_module
  use LBM_Discretization_D2Q9_module
  implicit none

  private
#include "lbm_definitions.h"

  public:: DiscretizationCreate, &
       DiscretizationDestroy, &
       DiscretizationSetType, &
       DiscretizationSetDerivOrder, &
       DiscretizationSetUp, &
       DiscretizationSetUpRelax, &
       DiscretizationEquilf, &
       DiscSetLocalDirections, &
       DiscApplyBCDirichletToBoundary, &
       DiscApplyBCFluxToBoundary, &
       DiscApplyBCVelocityToBoundary

contains
  function DiscretizationCreate(comm) result(disc)
    MPI_Comm comm
    type(discretization_type),pointer:: disc

    allocate(disc)
    disc%comm = comm
    disc%name = NULL_DISCRETIZATION
    disc%ndims = -1
    disc%b = -1
    disc%c_0 = 0
    nullify(disc%ci)
    nullify(disc%weights)
    nullify(disc%opposites)
    nullify(disc%reflect_x)
    nullify(disc%reflect_y)
    nullify(disc%reflect_z)
    nullify(disc%mt_mrt)
    nullify(disc%mmt_mrt)
    nullify(disc%ffw)
  end function DiscretizationCreate

  subroutine DiscretizationDestroy(disc, ierr)
    type(discretization_type) disc
    PetscErrorCode ierr

    if (associated(disc%ci)) deallocate(disc%ci)
    if (associated(disc%weights)) deallocate(disc%weights)
    if (associated(disc%opposites)) deallocate(disc%opposites)
    if (associated(disc%reflect_x)) deallocate(disc%reflect_x)
    if (associated(disc%reflect_y)) deallocate(disc%reflect_y)
    if (associated(disc%reflect_z)) deallocate(disc%reflect_z)
    if (associated(disc%mt_mrt)) deallocate(disc%mt_mrt)
    if (associated(disc%mmt_mrt)) deallocate(disc%mmt_mrt)
    if (associated(disc%ffw)) deallocate(disc%ffw)
  end subroutine DiscretizationDestroy

  subroutine DiscretizationSetDerivOrder(disc, isotropy_order)
    type(discretization_type) disc
    PetscInt isotropy_order
    PetscErrorCode ierr

    if ((isotropy_order.eq.0).or.(isotropy_order.eq.4).or. &
         (isotropy_order.eq.8).or.(isotropy_order.eq.10)) then
       disc%isotropy_order = isotropy_order
    else
       call LBMError(PETSC_COMM_SELF, 1, 'Invalid derivative order', ierr)
    end if
  end subroutine DiscretizationSetDerivOrder

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
       if (disc%ndims < 0 ) call LBMError(PETSC_COMM_SELF, 1, 'Invalid Discretization', ierr)
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
       call LBMError(PETSC_COMM_SELF,1,'invalid discretization in LBM',ierr)
    end select
  end subroutine DiscretizationSetUpRelax

  subroutine DiscretizationEquilf(disc, rho, u, walls, feq, m, relax, dist)
    use LBM_Distribution_Function_type_module
    use LBM_Relaxation_module
    type(discretization_type) disc
    type(distribution_type) dist
    type(relaxation_type) relax

    PetscScalar,dimension(dist%s,0:dist%b,dist%info%gxyzl):: feq
    PetscScalar,dimension(dist%s,dist%info%rgxyzl):: rho
    PetscScalar,dimension(dist%s,dist%info%ndims,dist%info%gxyzl):: u
    PetscScalar,dimension(dist%info%rgxyzl):: walls
    PetscInt m

    PetscErrorCode ierr

    select case(disc%name)
    case(D3Q19_DISCRETIZATION)
       call DiscretizationEquilf_D3Q19(disc, rho, u, walls, feq, m, relax, dist)
    case(D2Q9_DISCRETIZATION)
       call DiscretizationEquilf_D2Q9(disc, rho, u, walls, feq, m, relax, dist)
    case DEFAULT
       call LBMError(PETSC_COMM_SELF,1,'invalid discretization in LBM',ierr)
    end select
  end subroutine DiscretizationEquilf

  subroutine DiscSetLocalDirections(disc, boundary, directions, cardinals)
    type(discretization_type) disc
    PetscInt,intent(in):: boundary
    PetscInt,intent(out),dimension(0:disc%b) :: directions
    PetscInt,intent(out),dimension(1:disc%ndims) :: cardinals
    PetscErrorCode ierr

    select case(disc%name)
    case(D3Q19_DISCRETIZATION)
       call DiscSetLocalDirections_D3Q19(disc, boundary, directions, cardinals)
    case(D2Q9_DISCRETIZATION)
       call DiscSetLocalDirections_D2Q9(disc, boundary, directions, cardinals)
    case DEFAULT
       call LBMError(PETSC_COMM_SELF,1,'invalid discretization in LBM',ierr)
    end select
  end subroutine DiscSetLocalDirections

  subroutine DiscApplyBCDirichletToBoundary(disc, fi, forces, vals, directions, cardinals, dist)
    use LBM_Distribution_Function_type_module
    type(discretization_type) disc
    type(distribution_type) dist
    PetscInt,intent(in),dimension(0:dist%b):: directions
    PetscInt,intent(in),dimension(1:dist%info%ndims):: cardinals
    PetscScalar,intent(inout),dimension(1:dist%s, 0:dist%b):: fi
    PetscScalar,intent(in),dimension(dist%s, dist%info%ndims):: forces
    PetscScalar,intent(in),dimension(dist%s,dist%info%ndims):: vals
    PetscErrorCode ierr

    select case(disc%name)
    case(D3Q19_DISCRETIZATION)
       call DiscApplyBCDirichletToBoundary_D3Q19(disc, fi, forces, vals, directions, cardinals, dist)
    case(D2Q9_DISCRETIZATION)
       call DiscApplyBCDirichletToBoundary_D2Q9(disc, fi, forces, vals, directions, cardinals, dist)
    case DEFAULT
       call LBMError(PETSC_COMM_SELF,1,'invalid discretization in LBM',ierr)
    end select
  end subroutine DiscApplyBCDirichletToBoundary

  subroutine DiscApplyBCFluxToBoundary(disc, fi, forces, vals, directions, cardinals, dist)
    use LBM_Distribution_Function_type_module
    type(discretization_type) disc
    type(distribution_type) dist
    PetscInt,intent(in),dimension(0:disc%b):: directions
    PetscInt,intent(in),dimension(0:disc%ndims):: cardinals
    PetscScalar,intent(inout),dimension(1:dist%s, 0:disc%b):: fi
    PetscScalar,intent(in),dimension(dist%s, dist%info%ndims):: forces
    PetscScalar,intent(in),dimension(dist%s,dist%info%ndims):: vals
    PetscErrorCode ierr

    select case(disc%name)
    case(D3Q19_DISCRETIZATION)
       call DiscApplyBCFluxToBoundary_D3Q19(disc, fi, forces, vals, directions, cardinals, &
            dist)
    case(D2Q9_DISCRETIZATION)
       call DiscApplyBCFluxToBoundary_D2Q9(disc, fi, forces, vals, directions, cardinals, &
            dist)
    case DEFAULT
       call LBMError(PETSC_COMM_SELF,1,'invalid discretization in LBM',ierr)
    end select
  end subroutine DiscApplyBCFluxToBoundary

  subroutine DiscApplyBCVelocityToBoundary(disc, fi, forces, vals, directions, cardinals, &
       dist)
    use LBM_Distribution_Function_type_module
    type(discretization_type) disc
    type(distribution_type) dist
    PetscInt,intent(in),dimension(0:disc%b):: directions
    PetscInt,intent(in),dimension(0:disc%ndims):: cardinals
    PetscScalar,intent(inout),dimension(1:dist%s, 0:disc%b):: fi
    PetscScalar,intent(in),dimension(dist%s, dist%info%ndims):: forces
    PetscScalar,intent(in),dimension(dist%s, dist%info%ndims):: vals
    PetscErrorCode ierr

    select case(disc%name)
    case(D3Q19_DISCRETIZATION)
       call DiscApplyBCVelocityToBoundary_D3Q19(disc, fi, forces, vals, directions, cardinals, &
            dist)
    case(D2Q9_DISCRETIZATION)
       call DiscApplyBCVelocityToBoundary_D2Q9(disc, fi, forces, vals, directions, cardinals, &
            dist)
    case DEFAULT
       call LBMError(PETSC_COMM_SELF,1,'invalid discretization in LBM',ierr)
    end select
  end subroutine DiscApplyBCVelocityToBoundary
end module LBM_Discretization_module
