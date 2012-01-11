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
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

  subroutine initialize_walls(walls, filename, info)
    use petsc
    use LBM_Info_module
    implicit none

#include "lbm_definitions.h"
!   input variables
    type(info_type) info
    PetscScalar,dimension(info%rgxyzl) :: walls
    character(len=MAXSTRINGLENGTH) filename

    select case(info%ndims)
    case(3)
       call initialize_walls_d3(walls, filename, info)
    case(2)
       call initialize_walls_d2(walls, filename, info)
    end select
  end subroutine initialize_walls

  subroutine initialize_walls_d3(walls, filename, info)
    use petsc
    use LBM_Info_module
    implicit none

#include "lbm_definitions.h"
!   input variables
    type(info_type) info
    PetscScalar,dimension(info%rgxs:info%rgxe, &
         info%rgys:info%rgye, &
         info%rgzs:info%rgze):: walls
    character(len=MAXSTRINGLENGTH) filename

    walls=0

    if (info%ys.eq.1) walls(:,1,:) = WALL_NORMAL_Y
    if (info%ye.eq.info%NY) walls(:,info%NY,:) = WALL_NORMAL_Y

    if (info%zs.eq.1) walls(:,:,1) = WALL_NORMAL_Z
    if (info%ze.eq.info%NZ) walls(:,:,info%NZ) = WALL_NORMAL_Z
    return
  end subroutine initialize_walls_d3

  subroutine initialize_walls_d2(walls, filename, info)
    use petsc
    use LBM_Info_module
    implicit none

#include "lbm_definitions.h"
!   input variables
    type(info_type) info
    PetscScalar,dimension(info%rgxs:info%rgxe, &
         info%rgys:info%rgye):: walls
    character(len=MAXSTRINGLENGTH) filename

    walls=0

    if (info%ys.eq.1) walls(:,1) = WALL_NORMAL_Y
    if (info%ye.eq.info%NY) walls(:,info%NY) = WALL_NORMAL_Y
    return
  end subroutine initialize_walls_d2