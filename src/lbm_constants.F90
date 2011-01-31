!!!==================================================================
!!! Fortran-file
!!!    author:        Ethan T. Coon
!!!    filename:      constants.f
!!!    version:
!!!    created:       10 November 2010
!!!      on:          13:04:46 MST
!!!    last modified:  10 November 2010
!!!      at:          13:04:46 MST
!!!    URL:           http://www.ldeo.columbia.edu/~ecoon/
!!!    email:         ecoon _at_ ldeo.columbia.edu
!!!
!!!==================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"
  
module LBM_Constants_module
  implicit none

  private 

  type, public:: constants_type
     PetscScalar g,g11,g22          ! constants for mobility forces
     PetscScalar,pointer,dimension(:):: tau  ! relaxation times
     PetscScalar,pointer,dimension(:):: gvt  ! gravity (i.e. body forces)
     PetscScalar,pointer,dimension(:):: gw   ! fluid-solid interaction forces
     PetscScalar,pointer,dimension(:):: mm   ! mass per number (rho = mm*n)
     PetscScalar,pointer,dimension(:):: alf   ! mass per number (rho = mm*n)
     PetscScalar,pointer,dimension(:):: rho1, rho2         ! left and right fluid densities?
  end type constants_type

  public :: ConstantsCreate, &
       ConstantsSetFromOptions, &
       ConstantsDestroy

contains
  function ConstantsCreate() result(constants)
    type(constants_type),pointer:: constants

    allocate(constants)
    constants%g = 0
    constants%g11 = 0
    constants%g22 = 0

    nullify(constants%tau)
    nullify(constants%gvt)
    nullify(constants%gw)
    nullify(constants%mm)
    nullify(constants%alf)
    nullify(constants%rho1)
    nullify(constants%rho2)

  end function ConstantsCreate
    
  subroutine ConstantsSetSizes(constants, s)
    type(constants_type) constants
    PetscInt s

    allocate(tau(1:s))
    allocate(gvt(1:s))
    allocate(gw(1:s))
    allocate(mm(1:s))
    allocate(alf(1:s))
    allocate(rho1(1:s))
    allocate(rho2(1:s))
  end subroutine ConstantsSetSizes

  subroutine ConstantsSetFromOptions(constants, options, ierr)
    use LBM_Options_module
    type(constants_type) constants
    type(options_type) options
    PetscErrorCode ierr

    PetscInt nmax
    
    call ConstantsSetSizes(constants, options%s)

    call PetscOptionsGetReal(options%my_prefix, '-g', constants%g, flag, ierr)
    call PetscOptionsGetReal(options%my_prefix, '-g11', constants%g, flag, ierr)
    call PetscOptionsGetReal(options%my_prefix, '-g22', constants%g, flag, ierr)

    nmax = options%s
    call PetscOptionsGetRealArray(options%my_prefix, '-tau', constants%tau, nmax, flag, ierr)
    nmax = options%s
    call PetscOptionsGetRealArray(options%my_prefix, '-gvt', constants%gvt, nmax, flag, ierr)
    nmax = options%s
    call PetscOptionsGetRealArray(options%my_prefix, '-gw', constants%gw, nmax, flag, ierr)
    nmax = options%s
    call PetscOptionsGetRealArray(options%my_prefix, '-mm', constants%mm, nmax, flag, ierr)
    nmax = options%s
    call PetscOptionsGetRealArray(options%my_prefix, '-rho1', constants%rho1, nmax, flag, ierr)
    nmax = options%s
    call PetscOptionsGetRealArray(options%my_prefix, '-rho2', constants%rho2, nmax, flag, ierr)

    constants%alf = 1.-0.555555555/mm

  end subroutine ConstantsSetFromOptions
  
  subroutine ConstantsDestroy(constants, ierr)
    type(constants_type) constants
    PetscErrorCode ierr
    
    if (associated(constants%tau)) deallocate(constants%tau)
    if (associated(constants%gvt)) deallocate(constants%gvt)
    if (associated(constants%gw)) deallocate(constants%gw)
    if (associated(constants%mm)) deallocate(constants%mm)
    if (associated(constants%rho1)) deallocate(constants%rho1)
    if (associated(constants%rho2)) deallocate(constants%rho2)
    
    deallocate(constants)
  end subroutine ConstantsDestroy
end module LBM_Constants_module
  
