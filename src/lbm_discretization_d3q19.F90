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

  public:: DiscretizationSetup_D3Q19, &
       DiscretizationSetupConstants_D3Q19
contains
  subroutine DiscretizationSetup_D3Q19(disc)
    type(discretization_type) disc
    disc%name = D3Q19_DISCRETIZATION
    disc%ndims = 3
    disc%b = 18
    allocate(disc%ci(0:disc%b,1:disc%ndims))
    allocate(disc%weights(0:disc%b))
    allocate(disc%m(0:disc%b,0:disc%b))                 ! transformation matrix
    allocate(disc%mt(0:disc%b,0:disc%b))                ! transpose of M
    allocate(disc%mmt(0:disc%b))                        ! diagonal M dot MT matrix 
    !allocate(disc%S(0:disc%b))                         ! diagonal relaxation matrix

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

    disc%m(:, 0) = (/  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /)
    disc%m(:, 1) = (/-30,-11,-11,-11,-11,-11,-11,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8 /)
    disc%m(:, 2) = (/ 12, -4, -4, -4, -4, -4, -4,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /)
    disc%m(:, 3) = (/  0,  1,  0, -1,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,  0,  0,  0,  0 /)
    disc%m(:, 4) = (/  0, -4,  0,  4,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,  0,  0,  0,  0 /)
    disc%m(:, 5) = (/  0,  0,  1,  0, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1, -1,  1 /)
    disc%m(:, 6) = (/  0,  0, -4,  0,  4,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1, -1,  1 /)
    disc%m(:, 7) = (/  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1 /)
    disc%m(:, 8) = (/  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1 /)
    disc%m(:, 9) = (/  0,  2, -1,  2, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1, -2, -2, -2, -2 /)
    disc%m(:,10) = (/  0, -4,  2, -4,  2,  2,  2,  1,  1,  1,  1,  1,  1,  1,  1, -2, -2, -2, -2 /)
    disc%m(:,11) = (/  0,  0,  1,  0,  1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0 /)
    disc%m(:,12) = (/  0,  0, -2,  0, -2,  2,  2,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0 /)
    disc%m(:,13) = (/  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0 /)
    disc%m(:,14) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1 /)
    disc%m(:,15) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0 /)
    disc%m(:,16) = (/  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1, -1,  1,  1, -1,  0,  0,  0,  0 /)
    disc%m(:,17) = (/  0,  0,  0,  0,  0,  0,  0, -1, -1,  1,  1,  0,  0,  0,  0,  1, -1, -1,  1 /)
    disc%m(:,18) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1, -1, -1, -1,  1,  1 /)

    disc%mt(:, 0) = (/ 1,-30, 12,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /) 
    disc%mt(:, 1) = (/ 1,-11, -4,  1, -4,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,  0,  0,  0 /)
    disc%mt(:, 2) = (/ 1,-11, -4,  0,  0,  1, -4,  0,  0, -1,  2,  1, -2,  0,  0,  0,  0,  0,  0 /)
    disc%mt(:, 3) = (/ 1,-11, -4, -1,  4,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,  0,  0,  0 /)
    disc%mt(:, 4) = (/ 1,-11, -4,  0,  0, -1,  4,  0,  0, -1,  2,  1, -2,  0,  0,  0,  0,  0,  0 /)
    disc%mt(:, 5) = (/ 1,-11, -4,  0,  0,  0,  0,  1, -4, -1,  2, -1,  2,  0,  0,  0,  0,  0,  0 /)
    disc%mt(:, 6) = (/ 1,-11, -4,  0,  0,  0,  0, -1,  4, -1,  2, -1,  2,  0,  0,  0,  0,  0,  0 /)
    disc%mt(:, 7) = (/ 1,  8,  1,  1,  1,  1,  1,  0,  0,  1,  1,  1,  1,  1,  0,  0,  1, -1,  0 /)
    disc%mt(:, 8) = (/ 1,  8,  1, -1, -1,  1,  1,  0,  0,  1,  1,  1,  1, -1,  0,  0, -1, -1,  0 /)
    disc%mt(:, 9) = (/ 1,  8,  1, -1, -1, -1, -1,  0,  0,  1,  1,  1,  1,  1,  0,  0, -1,  1,  0 /) 
    disc%mt(:,11) = (/ 1,  8,  1,  1,  1, -1, -1,  0,  0,  1,  1,  1,  1, -1,  0,  0,  1,  1,  0 /)
    disc%mt(:,12) = (/ 1,  8,  1,  1,  1,  0,  0,  1,  1,  1,  1, -1, -1,  0,  0,  1, -1,  0,  1 /)
    disc%mt(:,13) = (/ 1,  8,  1, -1, -1,  0,  0,  1,  1,  1,  1, -1, -1,  0,  0, -1,  1,  0,  1 /)
    disc%mt(:,14) = (/ 1,  8,  1, -1, -1,  0,  0, -1, -1,  1,  1, -1, -1,  0,  0,  1,  1,  0, -1 /)
    disc%mt(:,15) = (/ 1,  8,  1,  1,  1,  0,  0, -1, -1,  1,  1, -1, -1,  0,  0, -1, -1,  0, -1 /)
    disc%mt(:,16) = (/ 1,  8,  1,  0,  0,  1,  1,  1,  1, -2, -2,  0,  0,  0,  1,  0,  0,  1, -1 /)
    disc%mt(:,17) = (/ 1,  8,  1,  0,  0, -1, -1,  1,  1, -2, -2,  0,  0,  0, -1,  0,  0, -1, -1 /)
    disc%mt(:,18) = (/ 1,  8,  1,  0,  0, -1, -1, -1, -1, -2, -2,  0,  0,  0,  1,  0,  0, -1,  1 /)
    disc%mt(:,19) = (/ 1,  8,  1,  0,  0,  1,  1, -1, -1, -2, -2,  0,  0,  0, -1,  0,  0,  1,  1 /)

    disc%mmt = (/ 19, 2394, 252, 10, 40, 10, 40, 10, 40, 36, 72, 12, 24, 4, 4, 4, 8, 8, 8 /) 

    !disc%S = (/ /)
    
  end subroutine DiscretizationSetup_D3Q19

  subroutine DiscretizationSetupConstants_D3Q19(disc, constants)
    use LBM_Constants_module
    type(discretization_type) disc
    type(constants_type) constants
    
    constants%alpha_0 = constants%d_k
    constants%alpha_1 = -1.d0/2.d0
  end subroutine DiscretizationSetupConstants_D3Q19
end module LBM_Discretization_D3Q19_module

