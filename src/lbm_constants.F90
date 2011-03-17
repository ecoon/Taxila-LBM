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
  
module LBM_Constants_module
  implicit none

  private 
#include "lbm_definitions.h"

  type, public:: constants_type
     PetscInt s
     PetscScalar g,g11,g22          ! constants for mobility forces
     PetscScalar,pointer,dimension(:):: tau  ! relaxation times
     PetscScalar,pointer,dimension(:,:):: gvt  ! gravity (i.e. body forces)
     PetscScalar,pointer,dimension(:):: gw   ! fluid-solid interaction forces
     PetscScalar,pointer,dimension(:):: mm   ! mass per number (rho = mm*n)
     PetscScalar c_s2   ! sound speed squared
     PetscScalar,pointer,dimension(:):: d_k   ! 
     PetscScalar,pointer,dimension(:):: alpha_0   ! 
     PetscScalar alpha_1
     PetscBool mrt
  end type constants_type

  public :: ConstantsCreate, &
       ConstantsSetFromOptions, &
       ConstantsView, &
       ConstantsDestroy

contains
  function ConstantsCreate() result(constants)
    type(constants_type),pointer:: constants

    allocate(constants)
    constants%g = 0
    constants%g11 = 0
    constants%g22 = 0
    constants%c_s2 = 0
    constants%mrt = .FALSE.

    nullify(constants%tau)
    nullify(constants%gvt)
    nullify(constants%gw)
    nullify(constants%mm)
    nullify(constants%d_k)
    nullify(constants%alpha_0)
  end function ConstantsCreate
    
  subroutine ConstantsSetFromOptions(constants, options, ierr)
    use LBM_Options_module
    
    type(constants_type) constants
    type(options_type) options
    PetscErrorCode ierr

    PetscInt nmax
    PetscBool flag

    constants%s = options%nphases
    allocate(constants%tau(1:options%nphases))
    allocate(constants%gvt(1:options%nphases,1:options%ndims))
    allocate(constants%gw(1:options%nphases))
    allocate(constants%mm(1:options%nphases))
    allocate(constants%d_k(1:options%nphases))
    allocate(constants%alpha_0(1:options%nphases))

    ! defaults
    constants%tau = 1.d0
    constants%gvt = 0.d0
    constants%gw = 0.d0
    constants%mm = 1.d0
    constants%c_s2 = 1.d0/3.d0

    if (options%nphases > 1) then
       call PetscOptionsGetReal(options%my_prefix, '-g', constants%g, flag, ierr)
       call PetscOptionsGetReal(options%my_prefix, '-g11', constants%g11, flag, ierr)
       call PetscOptionsGetReal(options%my_prefix, '-g22', constants%g22, flag, ierr)
    end if

    nmax = constants%s
    call PetscOptionsGetRealArray(options%my_prefix, '-tau', constants%tau, nmax, &
         flag, ierr)

    nmax = constants%s
    call PetscOptionsGetRealArray(options%my_prefix, '-gvt_x', &
         constants%gvt(:,X_DIRECTION), nmax, flag, ierr)
    nmax = constants%s
    call PetscOptionsGetRealArray(options%my_prefix, '-gvt_y', &
         constants%gvt(:,Y_DIRECTION), nmax, flag, ierr)

    if (options%ndims > 2) then
       nmax = constants%s
       call PetscOptionsGetRealArray(options%my_prefix, '-gvt_z', &
            constants%gvt(:,Z_DIRECTION), nmax, flag, ierr)
       if (.not.flag) then
          nmax = constants%s
          call PetscOptionsGetRealArray(options%my_prefix, '-gvt', &
               constants%gvt(:,Z_DIRECTION), nmax, flag, ierr)
       end if
    end if

    nmax = constants%s
    call PetscOptionsGetRealArray(options%my_prefix, '-gw', constants%gw, nmax, &
         flag, ierr)
    nmax = constants%s
    call PetscOptionsGetRealArray(options%my_prefix, '-mm', constants%mm, nmax, &
         flag, ierr)
    constants%d_k = 1.d0 - 2.d0/(3.d0*constants%mm)

    call PetscOptionsGetBool(options%my_prefix, '-flow_mrt', &
         constants%mrt, flag, ierr)
  end subroutine ConstantsSetFromOptions
  
  subroutine ConstantsView(constants)
    type(constants_type) constants
    
    print*, ' Physics:'
    print*, '  gw =', constants%gw
    print*, '  tau =', constants%tau
    print*, '  gvt =', constants%gvt
    print*, '  mm =', constants%mm
    if (constants%s > 1) then
       print*, '  g =', constants%g
       print*, '  g11 =', constants%g11
       print*, '  g22 =', constants%g22
    end if
  end subroutine ConstantsView

  subroutine ConstantsDestroy(constants, ierr)
    type(constants_type) constants
    PetscErrorCode ierr
    
    if (associated(constants%tau)) deallocate(constants%tau)
    if (associated(constants%gvt)) deallocate(constants%gvt)
    if (associated(constants%gw)) deallocate(constants%gw)
    if (associated(constants%mm)) deallocate(constants%mm)
    if (associated(constants%d_k)) deallocate(constants%d_k)
    if (associated(constants%alpha_0)) deallocate(constants%alpha_0)
  end subroutine ConstantsDestroy
end module LBM_Constants_module
  
