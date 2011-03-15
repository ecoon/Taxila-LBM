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
  
module LBM_Discretization_Directions_D3Q19_module
  implicit none

  PetscInt, parameter, public :: ORIGIN = 0
  PetscInt, parameter, public :: EAST = 1
  PetscInt, parameter, public :: NORTH = 2
  PetscInt, parameter, public :: WEST = 3
  PetscInt, parameter, public :: SOUTH = 4
  PetscInt, parameter, public :: UP = 5
  PetscInt, parameter, public :: DOWN = 6
  PetscInt, parameter, public :: NORTHEAST = 7
  PetscInt, parameter, public :: NORTHWEST = 8
  PetscInt, parameter, public :: SOUTHWEST = 9
  PetscInt, parameter, public :: SOUTHEAST = 10
  PetscInt, parameter, public :: EASTUP = 11
  PetscInt, parameter, public :: WESTUP = 12
  PetscInt, parameter, public :: WESTDOWN = 13
  PetscInt, parameter, public :: EASTDOWN = 14
  PetscInt, parameter, public :: NORTHUP = 15
  PetscInt, parameter, public :: SOUTHUP = 16
  PetscInt, parameter, public :: SOUTHDOWN = 17
  PetscInt, parameter, public :: NORTHDOWN = 18

  PetscInt, public, parameter :: discretization_dims = 3
  PetscInt, public, parameter :: discretization_directions = 18
end module LBM_Discretization_Directions_D3Q19_module

module LBM_Discretization_D3Q19_module
  use LBM_Discretization_Type_module
  use LBM_Discretization_Directions_D3Q19_module
  implicit none

  private
#include "lbm_definitions.h"

  public:: DiscretizationSetup_D3Q19

contains
  subroutine DiscretizationSetup_D3Q19(disc)
    type(discretization_type) disc
    disc%name = D3Q19_DISCRETIZATION
    disc%ndims = 3
    disc%b = 18
    allocate(disc%ci(0:disc%b,1:disc%ndims))
    allocate(disc%weights(0:disc%b))

    disc%ci(:,X_DIRECTION) = (/ 0, 1, 0,-1, 0, 0, 0, 1,-1,-1, &
         1, 1,-1,-1, 1, 0, 0, 0, 0/)
    disc%ci(:,Y_DIRECTION) = (/ 0, 0, 1, 0,-1, 0, 0, 1, 1,-1, &
        -1, 0, 0, 0, 0, 1,-1,-1, 1/)
    disc%ci(:,Z_DIRECTION) = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, &
         0, 1, 1,-1,-1, 1, 1,-1,-1/)

    disc%weights = (/ 1.d0/3.d0, &
       1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, &
       1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, &
       1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0/)
  end subroutine DiscretizationSetup_D3Q19

  subroutine DiscretizationSetupConstants_D3Q19(disc, constants)
    use LBM_Constants_module
    type(discretization_type) disc
    type(constants_type) constants
    
    constants%alpha_0 = constants%d_k
    constants%alpha_1 = -1.d0/2.d0
  end subroutine DiscretizationSetupConstants_D3Q19
end module LBM_Discretization_D3Q19_module

