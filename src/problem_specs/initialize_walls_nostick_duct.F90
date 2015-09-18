!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_walls.F90
!!!     version:
!!!     created:         14 January 2011
!!!       on:            17:26:22 MST
!!!     last modified:   01 October 2011
!!!       at:            15:53:25 PDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petscdmdef.h"

  subroutine initialize_walls(walls, info, options)
    use petsc
    use LBM_Options_module
    use LBM_Info_module
    implicit none

#include "lbm_definitions.h"
!   input variables
    type(
info_type) info
    type(options_type) options
    PetscScalar,dimension(info%rgxyzl) :: walls

    select case(info%ndims)
    case(3)
       call initialize_walls_d3(walls, info, options)
    case(2)
       call initialize_walls_d2(walls, info, options)
    end select
  end subroutine initialize_walls

  subroutine initialize_walls_d3(walls, info, options)
    use petsc
    use LBM_Info_module
    implicit none

#include "lbm_definitions.h"
!   input variables
    type(info_type) info
    type(options_type) options
    PetscScalar,dimension(info%rgxs:info%rgxe, &
         info%rgys:info%rgye, &
         info%rgzs:info%rgze):: walls

    PetscBool flag
    PetscInt normal
    PetscErrorCode ierr

    walls=0
    normal = Y_DIRECTION
    call OptionsGetInt(options, "-duct_normal_direction", "normal direction", &
         normal, flag, ierr)

    select case(normal)
    case(X_DIRECTION)
      if (info%xs.eq.1) walls(1,:,:) = WALL_NORMAL_X
      if (info%xe.eq.info%NX) walls(info%NX,:,:) = WALL_NORMAL_X
    case(Y_DIRECTION)
      if (info%ys.eq.1) walls(:,1,:) = WALL_NORMAL_Y
      if (info%ye.eq.info%NY) walls(:,info%NY,:) = WALL_NORMAL_Y
    case(Z_DIRECTION)
      if (info%zs.eq.1) walls(:,:,1) = WALL_NORMAL_Z
      if (info%ze.eq.info%NY) walls(:,:,info%NZ) = WALL_NORMAL_Z
    end select
    return
  end subroutine initialize_walls_d3

  subroutine initialize_walls_d2(walls, filename, info)
    use petsc
    use LBM_Info_module
    implicit none

#include "lbm_definitions.h"
!   input variables
    type(info_type) info
    type(options_type) options
    PetscScalar,dimension(info%rgxs:info%rgxe, &
         info%rgys:info%rgye):: walls

    walls=0
    PetscBool flag
    PetscInt normal
    PetscErrorCode ierr

    walls=0

    normal = Y_DIRECTION
    call OptionsGetInt(options, "-duct_normal_direction", "normal direction", &
         normal, flag, ierr)

    select case(normal)
    case(X_DIRECTION)
      if (info%xs.eq.1) walls(1,:) = WALL_NORMAL_X
      if (info%xe.eq.info%NX) walls(info%NX,:) = WALL_NORMAL_X
    case(Y_DIRECTION)
      if (info%ys.eq.1) walls(:,1) = WALL_NORMAL_Y
      if (info%ye.eq.info%NY) walls(:,info%NY) = WALL_NORMAL_Y
    end select
    return
  end subroutine initialize_walls_d3
