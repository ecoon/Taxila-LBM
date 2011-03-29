!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_grid.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            09:24:24 MDT
!!!     last modified:   29 March 2011
!!!       at:            16:16:00 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscdmdef.h"

module LBM_Grid_module
  use petsc
  use LBM_Info_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public:: grid_type
     MPI_Comm comm
     DM,pointer,dimension(:) :: da
     type(info_type), pointer:: info
     PetscInt,pointer,dimension(:) :: da_sizes
     PetscInt nda
     character(len=MAXWORDLENGTH) name       
  end type grid_type

  public:: GridCreate, &
       GridDestroy, &
       GridSetName, &
       GridSetFromOptions

contains
  function GridCreate(comm) result(grid)
    MPI_Comm comm
    character(len=MAXWORDLENGTH) name       
    type(grid_type),pointer :: grid
    allocate(grid)
    grid%comm = comm
    nullify(grid%da)
    nullify(grid%da_sizes)
    grid%info => InfoCreate(comm)
    grid%nda = 0
    grid%name = ''
  end function GridCreate

  subroutine GridDestroy(grid, ierr)
    type(grid_type) grid
    PetscErrorCode ierr
    integer lcv

    if (associated(grid%da)) then
       do lcv=1,grid%nda
          call DMDestroy(grid%da(lcv),ierr)
       end do
    end if
    if (associated(grid%da_sizes)) deallocate(grid%da_sizes)
    if (associated(grid%info)) call InfoDestroy(grid%info, ierr)
  end subroutine GridDestroy

  subroutine GridSetName(grid, name)
    type(grid_type) grid
    character(len=MAXWORDLENGTH) name       
    grid%name = name
  end subroutine GridSetName

  subroutine GridSetFromOptions(grid, options, ierr)
    use LBM_Options_module

    type(grid_type) grid
    type(options_type) options
    PetscErrorCode ierr

    call InfoSetFromOptions(grid%info, options, ierr)
    
    if (options%transport_disc /= NULL_DISCRETIZATION) then
       grid%nda = 6
    else
       grid%nda = 4
    end if

    allocate(grid%da(1:grid%nda))
    allocate(grid%da_sizes(1:grid%nda))
  end subroutine GridSetFromOptions

  subroutine GridSetUp(grid)
    type(grid_type) grid

    PetscErrorCode ierr
    PetscInt,dimension(3):: btype
    PetscInt xs,gxs
    PetscInt ys,gys
    PetscInt zs,gzs
    PetscInt,allocatable,dimension(:):: ownership_x, ownership_y, ownership_z
    integer lcv

    btype(:) = DMDA_BOUNDARY_GHOSTED
    if (grid%info%periodic(X_DIRECTION)) btype(X_DIRECTION) = DMDA_BOUNDARY_PERIODIC
    if (grid%info%periodic(Y_DIRECTION)) btype(Y_DIRECTION) = DMDA_BOUNDARY_PERIODIC
    if (grid%info%ndims > 2) then
       if (grid%info%periodic(Z_DIRECTION)) btype(Z_DIRECTION) = DMDA_BOUNDARY_PERIODIC
    else
       btype(Z_DIRECTION) = PETSC_NULL_INTEGER
    end if

    call DMDACreate3d(grid%comm, btype(X_DIRECTION), btype(Y_DIRECTION), &
         btype(Z_DIRECTION), grid%info%stencil_type, grid%info%NX, grid%info%NY, &
         grid%info%NZ, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, &
         grid%da_sizes(ONEDOF), grid%info%stencil_size, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, grid%da(ONEDOF), ierr)

    call DMDAGetCorners(grid%da(ONEDOF), xs, ys, zs, grid%info%xl, grid%info%yl, grid%info%zl, ierr)
    call DMDAGetGhostCorners(grid%da(ONEDOF), gxs, gys, gzs, grid%info%gxl, grid%info%gyl, &
         grid%info%gzl, ierr)
    call DMDAGetInfo(grid%da(ONEDOF), PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, grid%info%nproc_x, &
         grid%info%nproc_y, grid%info%nproc_z, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr)
    allocate(ownership_x(1:grid%info%nproc_x))
    allocate(ownership_y(1:grid%info%nproc_y))
    allocate(ownership_z(1:grid%info%nproc_z))
    ownership_x = 0
    ownership_y = 0
    ownership_z = 0
    call DMDAGetOwnershipRanges(grid%da(ONEDOF),ownership_x,ownership_y,ownership_z,ierr)
    
    ! set grid%info including corners
    grid%info%xs = xs+1
    grid%info%gxs = gxs+1
    grid%info%xe = grid%info%xs+grid%info%xl-1
    grid%info%gxe = grid%info%gxs+grid%info%gxl-1
    
    grid%info%ys = ys+1
    grid%info%gys = gys+1
    grid%info%ye = grid%info%ys+grid%info%yl-1
    grid%info%gye = grid%info%gys+grid%info%gyl-1
    
    if (grid%info%ndims > 2) then
       grid%info%zs = zs+1
       grid%info%gzs = gzs+1
       grid%info%ze = grid%info%zs+grid%info%zl-1
       grid%info%gze = grid%info%gzs+grid%info%gzl-1
    end if
    
    grid%info%xyzl = grid%info%xl*grid%info%yl*grid%info%zl
    grid%info%gxyzl = grid%info%gxl*grid%info%gyl*grid%info%gzl
    
    do lcv=2,grid%nda
       call DMDACreate3d(grid%comm, btype(X_DIRECTION), btype(Y_DIRECTION), &
            btype(Z_DIRECTION), grid%info%stencil_type, &
            grid%info%NX, grid%info%NY, grid%info%NZ, &
            grid%info%nproc_x, grid%info%nproc_y, grid%info%nproc_z, &
            grid%da_sizes(lcv), grid%info%stencil_size, &
            ownership_x, ownership_y, ownership_z, &
            grid%da(lcv), ierr)
    end do
    deallocate(ownership_x)
    deallocate(ownership_y)
    deallocate(ownership_z)

    ! set the names
    call PetscObjectSetName(grid%da(ONEDOF), trim(grid%name)//'DA_one', ierr)
    call PetscObjectSetName(grid%da(NPHASEDOF), trim(grid%name)//'DA_NPHASE', ierr)
    call PetscObjectSetName(grid%da(NPHASEXBDOF), trim(grid%name)//'DA_NPHASEXB', ierr)
    call PetscObjectSetName(grid%da(NFLOWDOF), trim(grid%name)//'DA_NFLOW', ierr)
    if (grid%nda > 4) then
       call PetscObjectSetName(grid%da(NCOMPONENTDOF), &
            trim(grid%name)//'DA_NCOMPONENT', ierr)
       call PetscObjectSetName(grid%da(NCOMPONENTXBDOF), &
            trim(grid%name)//'DA_NCOMPONENTXB', ierr)
    end if
    
    ! set the coordinates, but only on the one-dof (no need for extra memory)
    if (grid%info%ndims == 2) then
       call DMDASetUniformCoordinates(grid%da(ONEDOF), &
            grid%info%corners(X_DIRECTION,1), grid%info%corners(X_DIRECTION,2), &
            grid%info%corners(Y_DIRECTION,1), grid%info%corners(Y_DIRECTION,2), &
            grid%info%corners(Z_DIRECTION,1), grid%info%corners(Z_DIRECTION,2), ierr)
    else
       call DMDASetUniformCoordinates(grid%da(ONEDOF), &
            grid%info%corners(X_DIRECTION,1), grid%info%corners(X_DIRECTION,2), &
            grid%info%corners(Y_DIRECTION,1), grid%info%corners(Y_DIRECTION,2), &
            PETSC_NULL_SCALAR, PETSC_NULL_SCALAR, ierr)
    end if
  end subroutine GridSetUp

  subroutine GridViewCoordinates(grid, io)
    use LBM_IO_module
    type(grid_type) grid
    type(io_type) io
    Vec coords
    PetscErrorCode ierr

    call DMDAGetCoordinates(grid%da(ONEDOF), coords, ierr)
    call IOView(io, coords, 'coords')
  end subroutine GridViewCoordinates
end module LBM_Grid_module
