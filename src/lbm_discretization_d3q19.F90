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
  
module LBM_Discretization_D3Q19_module
  implicit none
  
  PetscInt, parameter :: discretization_dims = 3
  PetscInt, parameter :: discretization_directions = 18

  ! directions
  PetscInt, parameter :: ORIGIN = 0
  PetscInt, parameter :: EAST = 1
  PetscInt, parameter :: NORTH = 2
  PetscInt, parameter :: WEST = 3
  PetscInt, parameter :: SOUTH = 4
  PetscInt, parameter :: UP = 5
  PetscInt, parameter :: DOWN = 6
  PetscInt, parameter :: NORTHEAST = 7
  PetscInt, parameter :: NORTHWEST = 8
  PetscInt, parameter :: SOUTHWEST = 9
  PetscInt, parameter :: SOUTHEAST = 10
  PetscInt, parameter :: EASTUP = 11
  PetscInt, parameter :: WESTUP = 12
  PetscInt, parameter :: WESTDOWN = 13
  PetscInt, parameter :: EASTDOWN = 14
  PetscInt, parameter :: NORTHUP = 15
  PetscInt, parameter :: SOUTHUP = 16
  PetscInt, parameter :: SOUTHDOWN = 17
  PetscInt, parameter :: NORTHDOWN = 18

  PetscScalar,parameter,dimension(0:discretization_directions):: &
       cix = (/0.d0, & 
               1.d0, &
               0.d0, &
              -1.d0, &
               0.d0, &
               0.d0, &
               0.d0, &
               1.d0, &
              -1.d0, &
              -1.d0, &
               1.d0, &
               1.d0, &
              -1.d0, &
              -1.d0, &
               1.d0, &
               0.d0, &
               0.d0, &
               0.d0, &
               0.d0/)

  PetscScalar,parameter,dimension(0:discretization_directions):: &
       ciy = (/0.d0, &
               0.d0, &
               1.d0, &
               0.d0, &
              -1.d0, &
               0.d0, &
               0.d0, &
               1.d0, &
               1.d0, &
              -1.d0, &
              -1.d0, &
               0.d0, &
               0.d0, &
               0.d0, &
               0.d0, &
               1.d0, &
              -1.d0, &
              -1.d0, &
               1.d0/)

  PetscScalar,parameter,dimension(0:discretization_directions):: &
       ciz = (/0.d0, &
               0.d0, &
               0.d0, &
               0.d0, &
               0.d0, &
               1.d0, &
              -1.d0, &
               0.d0, &
               0.d0, &
               0.d0, &
               0.d0, &
               1.d0, &
               1.d0, &
              -1.d0, &
              -1.d0, &
               1.d0, &
               1.d0, &
              -1.d0, &
              -1.d0/)

  PetscScalar,parameter,dimension(0:discretization_directions):: &
       weights = (/ 1.d0/3.d0, &
       1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, &
       1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, &
       1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0/)

  PetscScalar,dimension(0:discretization_directions, 1:discretization_dims):: ci

contains
  subroutine LBMDiscretizationSetup()
    ci(:,1) = cix
    ci(:,2) = ciy
    ci(:,3) = ciz
  end subroutine LBMDiscretizationSetup
end module LBM_Discretization_D3Q19_module

