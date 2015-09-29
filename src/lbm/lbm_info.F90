!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        grid_info.F90
!!!     version:         
!!!     created:         06 December 2010
!!!       on:            15:19:22 MST
!!!     last modified:   18 October 2011
!!!       at:            13:50:14 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ ldeo.columbia.edu
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "petsc/finclude/petscsysdef.h"

module LBM_Info_module
  use petsc
  implicit none

  private
#include "lbm_definitions.h"

  type, public:: info_type
     MPI_Comm comm
     PetscInt ndims

     PetscInt NX,NY,NZ
     PetscBool,pointer,dimension(:) :: periodic
     PetscScalar,pointer,dimension(:,:) :: corners
     PetscInt stencil_size_rho
     PetscInt stencil_size
     PetscInt stencil_type

     ! dependent parameters
     PetscInt xs,xe,xl,gxs,gxe,gxl,rgxs,rgxe,rgxl
     PetscInt ys,ye,yl,gys,gye,gyl,rgys,rgye,rgyl
     PetscInt zs,ze,zl,gzs,gze,gzl,rgzs,rgze,rgzl
     PetscInt xyzl,gxyzl,rgxyzl
     PetscInt nproc_x,nproc_y,nproc_z
     PetscMPIInt rank, nprocs
     PetscReal,pointer,dimension(:) :: gridsize
     PetscInt,pointer,dimension(:):: ownership_x,ownership_y,ownership_z
  end type info_type

  public :: InfoCreate, &
       InfoDestroy, &
       InfoSetFromOptions, &
       InfoView

contains
  function InfoCreate(comm) result(info)
    MPI_Comm comm
    type(info_type),pointer:: info
    allocate(info)
    info%comm = comm

    ! intialize
    info%NX = -1
    info%NY = -1
    info%NZ = -1
    info%stencil_size_rho = -1
    info%stencil_size = 1
    info%stencil_type = DMDA_STENCIL_BOX
    info%ndims = -1

    info%xs = -1
    info%xe = -1
    info%xl = -1
    info%ys = -1
    info%ye = -1
    info%yl = -1
    info%zs = -1
    info%ze = -1
    info%zl = -1

    info%xyzl = -1
    info%gxyzl = -1

    info%gxs = -1
    info%gxe = -1
    info%gxl = -1
    info%gys = -1
    info%gye = -1
    info%gyl = -1
    info%gzs = -1
    info%gze = -1
    info%gzl = -1

    info%rgxs = -1
    info%rgxe = -1
    info%rgxl = -1
    info%rgys = -1
    info%rgye = -1
    info%rgyl = -1
    info%rgzs = -1
    info%rgze = -1
    info%rgzl = -1

    info%rank = -1
    info%nprocs = -1
    info%nproc_x = PETSC_DECIDE
    info%nproc_y = PETSC_DECIDE
    info%nproc_z = PETSC_DECIDE

    nullify(info%ownership_x)
    nullify(info%ownership_y)
    nullify(info%ownership_z)
    nullify(info%periodic)
    nullify(info%gridsize)
    nullify(info%corners)
  end function InfoCreate

  subroutine InfoSetFromOptions(info, options, ierr)
    use LBM_Options_module

    type(info_type) info
    type(options_type) options
    PetscBool flag
    PetscErrorCode ierr

    ! allocate sizes
    info%ndims = options%ndims
    allocate(info%gridsize(info%ndims))
    info%gridsize = -1.d0
    allocate(info%corners(info%ndims,2))
    info%corners = 0.d0
    allocate(info%periodic(info%ndims))
    info%periodic = PETSC_FALSE

    ! get options
    call OptionsGroupHeader(options, " Grid Info Options", ierr)
    ! -- grid size
    call OptionsGetInt(options, "-NX", "grid size in X", info%NX, flag, ierr)
    call OptionsGetInt(options, "-NY", "grid size in Y", info%NY, flag, ierr)
    if (info%ndims > 2) then
      call OptionsGetInt(options, "-NZ", "grid size in Z", info%NZ, flag, ierr)
    else
      info%NZ = 1
    end if

    ! -- stencil info
    call OptionsGetInt(options, "-stencil_size", "number of grid points in the stencil", &
         info%stencil_size, flag, ierr)
    call OptionsGetInt(options, "-stencil_size_rho", &
         "number of grid points in the stencil for rho -- defaults to match derivative order", &
         info%stencil_size_rho, flag, ierr)
    call OptionsGetInt(options, "-stencil_type", &
         "stencil type: 0=STAR (i.e. D2Q5), 1=BOX (i.e. D2Q9)", &
         info%stencil_type, flag, ierr)

    ! -- grid perioidicty
    call OptionsGetBool(options, "-bc_periodic_x", "periodic in x-direction", &
         info%periodic(X_DIRECTION), flag, ierr)
    call OptionsGetBool(options, "-bc_periodic_y", "periodic in y-direction", &
         info%periodic(Y_DIRECTION), flag, ierr)
    if (info%ndims > 2) then
      call OptionsGetBool(options, "-bc_periodic_z", "periodic in z-direction", &
           info%periodic(Z_DIRECTION), flag, ierr)
    end if

    ! -- grid corners
    if (info%periodic(X_DIRECTION)) then
      info%corners(X_DIRECTION,2) = info%NX
    else
      info%corners(X_DIRECTION,2) = info%NX-1
    end if
    call OptionsGetReal(options, "-x_start", "lower x coordinate", &
         info%corners(X_DIRECTION,1), flag, ierr)
    call OptionsGetReal(options, "-x_end", "upper x coordinate", &
         info%corners(X_DIRECTION,2), flag, ierr)

    if (info%periodic(Y_DIRECTION)) then
      info%corners(Y_DIRECTION,2) = info%NY
    else
      info%corners(Y_DIRECTION,2) = info%NY-1
    end if
    call OptionsGetReal(options, "-y_start", "lower y coordinate", &
         info%corners(Y_DIRECTION,1), flag, ierr)
    call OptionsGetReal(options, "-y_end", "upper y coordinate", &
         info%corners(Y_DIRECTION,2), flag, ierr)

    if (info%ndims > 2) then
      if (info%periodic(Z_DIRECTION)) then
        info%corners(Z_DIRECTION,2) = info%NZ
      else
        info%corners(Z_DIRECTION,2) = info%NZ-1
      end if
      call OptionsGetReal(options, "-z_start", "lower z coordinate", &
           info%corners(Z_DIRECTION,1), flag, ierr)
      call OptionsGetReal(options, "-z_end", "upper z coordinate", &
           info%corners(Z_DIRECTION,2), flag, ierr)
    end if

    ! calculate grid spacing
    if (info%periodic(X_DIRECTION)) then
       info%gridsize(X_DIRECTION) = (info%corners(X_DIRECTION,2) - &
            info%corners(X_DIRECTION, 1))/info%NX
    else
       info%gridsize(X_DIRECTION) = (info%corners(X_DIRECTION,2) - &
            info%corners(X_DIRECTION, 1))/(info%NX-1)
    end if
    if (info%periodic(Y_DIRECTION)) then
       info%gridsize(Y_DIRECTION) = (info%corners(Y_DIRECTION,2) - &
            info%corners(Y_DIRECTION, 1))/(info%NY)
    else
       info%gridsize(Y_DIRECTION) = (info%corners(Y_DIRECTION,2) - &
            info%corners(Y_DIRECTION, 1))/(info%NY-1)
    end if
    if (info%ndims > 2) then
       if (info%periodic(Z_DIRECTION)) then
          info%gridsize(Z_DIRECTION) = (info%corners(Z_DIRECTION,2) - &
               info%corners(Z_DIRECTION, 1))/(info%NZ)
       else
          info%gridsize(Z_DIRECTION) = (info%corners(Z_DIRECTION,2) - &
               info%corners(Z_DIRECTION, 1))/(info%NZ-1)
       end if
    end if

    ! nullify z-things for a 2D DA
    if (info%ndims.eq.2) then
       info%zs = PETSC_NULL_INTEGER
       info%ze = PETSC_NULL_INTEGER
       info%zl = 1
       info%NZ = 1
       info%gzs = PETSC_NULL_INTEGER
       info%gze = PETSC_NULL_INTEGER
       info%gzl = 1
       info%nproc_z = PETSC_DECIDE
    end if

    call MPI_Comm_rank(info%comm, info%rank, ierr)
    call MPI_Comm_size(info%comm, info%nprocs, ierr)
  end subroutine InfoSetFromOptions

  subroutine InfoView(info)
    type(info_type) info
    print*, ' Domain size =', info%NX,info%NY,info%NZ
  end subroutine InfoView

  subroutine InfoDestroy(info, ierr)
    type(info_type) info
    PetscErrorCode ierr

    if (associated(info%gridsize)) deallocate(info%gridsize)
    if (associated(info%ownership_x)) deallocate(info%ownership_x)
    if (associated(info%ownership_y)) deallocate(info%ownership_y)
    if (associated(info%ownership_z)) deallocate(info%ownership_z)
    if (associated(info%corners)) deallocate(info%corners)
    if (associated(info%periodic)) deallocate(info%periodic)
  end subroutine InfoDestroy
end module LBM_Info_module
