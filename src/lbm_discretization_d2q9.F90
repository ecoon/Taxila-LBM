!!!==================================================================
!!! Fortran-file
!!!    author:        Ethan T. Coon
!!!    filename:      c_constants.f
!!!    version:
!!!    created:       05 November 2010
!!!      on:          12:16:11 MDT
!!!    last modified:  05 November 2010
!!!      at:          12:16:11 MDT
!!!    URL:           http://www.ldeo.columbia.edu/~ecoon/
!!!    email:         ecoon _at_ ldeo.columbia.edu
!!!
!!!==================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

module LBM_Discretization_Directions_D2Q9_module
  implicit none
  
  ! discretization-specific enum directions
  PetscInt, public, parameter :: ORIGIN = 0
  PetscInt, public, parameter :: EAST = 1
  PetscInt, public, parameter :: NORTH = 2
  PetscInt, public, parameter :: WEST = 3
  PetscInt, public, parameter :: SOUTH = 4
  PetscInt, public, parameter :: NORTHEAST = 5
  PetscInt, public, parameter :: NORTHWEST = 6
  PetscInt, public, parameter :: SOUTHWEST = 7
  PetscInt, public, parameter :: SOUTHEAST = 8

  PetscInt, public, parameter :: discretization_dims = 2
  PetscInt, public, parameter :: discretization_directions = 8
end module LBM_Discretization_Directions_D2Q9_module
  
module LBM_Discretization_D2Q9_module
  use LBM_Discretization_Type_module
  use LBM_Discretization_Directions_D2Q9_module
  implicit none

  private
#include "lbm_definitions.h"

  public:: DiscretizationSetup_D2Q9

contains
  subroutine DiscretizationSetup_D2Q9(disc)
    type(discretization_type) disc
    disc%name = D2Q9_DISCRETIZATION
    disc%ndims = discretization_dims
    disc%b = discretization_directions
    allocate(disc%ci(0:disc%b,1:disc%ndims))
    allocate(disc%weights(0:disc%b))

    disc%ci(:,X_DIRECTION) = (/ 0, 1, 0,-1, 0, 1,-1,-1, 1/)
    disc%ci(:,Y_DIRECTION) = (/ 0, 0, 1, 0,-1, 1, 1,-1,-1/)

    disc%weights = (/ 4.d0/9.d0, &
       1.d0/9.d0,  1.d0/9.d0,  1.d0/9.d0,  1.d0/9.d0, &
       1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0 /)
  end subroutine DiscretizationSetup_D2Q9

  subroutine DiscretizationSetupConstants_D2Q9(disc, constants)
    use LBM_Constants_module
    type(discretization_type) disc
    type(constants_type) constants
    
    constants%alpha_0 = (1.d0 + constants%d_k*5.d0)/6.d0
    constants%alpha_1 = -2.d0/3.d0
  end subroutine DiscretizationSetupConstants_D2Q9
end module LBM_Discretization_D2Q9_module

