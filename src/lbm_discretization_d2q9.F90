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
  
module LBM_Discretization_D2Q9_module
  implicit none
  
  PetscInt, parameter :: discretization_dims = 2
  PetscInt, parameter :: discretization_directions = 8

  ! directions
  PetscInt, parameter :: ORIGIN = 0
  PetscInt, parameter :: EAST = 1
  PetscInt, parameter :: NORTH = 2
  PetscInt, parameter :: WEST = 3
  PetscInt, parameter :: SOUTH = 4
  PetscInt, parameter :: NORTHEAST = 5
  PetscInt, parameter :: NORTHWEST = 6
  PetscInt, parameter :: SOUTHWEST = 7
  PetscInt, parameter :: SOUTHEAST = 8

  PetscScalar,parameter,dimension(0:discretization_directions):: &
       cix = (/ 0.d0, & 
                1.d0, &
                0.d0, &
               -1.d0, &
                0.d0, &
                1.d0, &
               -1.d0, &
               -1.d0, &
                1.d0/)
  PetscScalar,parameter,dimension(0:discretization_directions):: &
       ciy = (/ 0.d0, &
                0.d0, &
                1.d0, &
                0.d0, &
               -1.d0, &
                1.d0, &
                1.d0, &
               -1.d0, &
               -1.d0/)

  PetscScalar,parameter,dimension(0:discretization_directions):: &
       weights = (/ 4.d0/9.d0, &
       1.d0/9.d0, 1.d0/9.d0,1.d0/9.d0, 1.d0/9.d0, &
       1.d0/36.d0, 1.d0/36.d0,1.d0/36.d0, 1.d0/36.d0 /)
  
  PetscScalar,dimension(0:discretization_directions, 1:discretization_dims):: ci

contains
  subroutine LBMDiscretizationSetup()
    ci(:,1) = cix
    ci(:,2) = ciy
  end subroutine LBMDiscretizationSetup
 
end module LBM_Discretization_D2Q9_module

