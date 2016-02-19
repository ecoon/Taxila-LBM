!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_grid.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            09:24:24 MDT
!!!     last modified:   18 October 2011
!!!       at:            13:47:21 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscdmdef.h"

module LBM_Grid_module
  use petsc
  use LBM_Error_module
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
     PetscScalar length_scale
     character(len=MAXWORDLENGTH) name       
  end type grid_type

  public:: GridCreate, &
       GridDestroy, &
       GridSetName, &
       GridSetFromOptions, &
       GridSetPhysicalScales, &
       GridSetUp, &
       GridViewCoordinates

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
    PetscInt lcv

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

  subroutine GridSetPhysicalScales(grid, ierr)
    type(grid_type) grid
    PetscErrorCode ierr
    PetscScalar,parameter:: eps=1.e-8

    grid%length_scale = grid%info%gridsize(X_DIRECTION)
    if (abs(grid%info%gridsize(Y_DIRECTION)-grid%length_scale) > eps) then
       call LBMError(PETSC_COMM_WORLD, 1, &
            'Lattice sizes must be equal, check grid parameters.', ierr)
    end if
    if (grid%info%ndims > 2) then
       if (abs(grid%info%gridsize(Z_DIRECTION)-grid%length_scale) > eps) then
          call LBMError(PETSC_COMM_WORLD,1,'Lattice sizes must be equal, check grid parameters.', ierr)
       end if
    end if
  end subroutine GridSetPhysicalScales

  subroutine GridSetFromOptions(grid, options, ierr)
    use LBM_Options_module

    type(grid_type) grid
    type(options_type) options
    PetscErrorCode ierr

    call InfoSetFromOptions(grid%info, options, ierr)

    ! update a few things
    if (grid%info%stencil_size_rho < 0) then
       select case(options%isotropy_order)
       case (0)
          grid%info%stencil_size_rho = 0
       case (4)
          grid%info%stencil_size_rho = 1
       case (8)
          grid%info%stencil_size_rho = 2
       case (10)
          grid%info%stencil_size_rho = 3
       case default
          grid%info%stencil_size_rho = 1
       end select
    end if

    if (options%transport_disc /= NULL_DISCRETIZATION) then
       grid%nda = 6
    else
       grid%nda = 4
    end if

    allocate(grid%da(1:grid%nda))
    allocate(grid%da_sizes(1:grid%nda))
    grid%da_sizes(:) = 0
    grid%da(:) = 0
  end subroutine GridSetFromOptions

  subroutine GridSetUp(grid)
    type(grid_type) grid

    PetscErrorCode ierr
    PetscInt,dimension(3):: btype
    PetscInt xs,gxs,rgxs
    PetscInt ys,gys,rgys
    PetscInt zs,gzs,rgzs

    PetscInt lcv
    PetscReal zero
    PetscReal one

    zero = 0
    one = 1

    btype(:) = DM_BOUNDARY_GHOSTED
    if (grid%info%periodic(X_DIRECTION)) btype(X_DIRECTION) = DM_BOUNDARY_PERIODIC
    if (grid%info%periodic(Y_DIRECTION)) btype(Y_DIRECTION) = DM_BOUNDARY_PERIODIC
    if (grid%info%ndims > 2) then
       if (grid%info%periodic(Z_DIRECTION)) btype(Z_DIRECTION) = DM_BOUNDARY_PERIODIC
    else
       btype(Z_DIRECTION) = PETSC_NULL_INTEGER
    end if

    ! create the 1-dof DA, and set up ghost/local info
    if (associated(grid%info%ownership_x)) then
       ! most of this is pre-calculated or provided in options
       call DMDACreate3d(grid%comm, btype(X_DIRECTION), btype(Y_DIRECTION), &
            btype(Z_DIRECTION), grid%info%stencil_type, grid%info%NX, grid%info%NY, &
            grid%info%NZ, grid%info%nproc_x, grid%info%nproc_y, grid%info%nproc_z, &
            grid%da_sizes(ONEDOF), grid%info%stencil_size_rho, grid%info%ownership_x, &
            grid%info%ownership_y, grid%info%ownership_z,grid%da(ONEDOF), ierr)
    else
      ! petsc does a bad job of guessing a partitioning, because it thinks 
      ! we're in 3D, not in 2D, so give it a hint
       if (grid%info%ndims.eq.2) grid%info%nproc_z = 1

       call DMDACreate3d(grid%comm, btype(X_DIRECTION), btype(Y_DIRECTION), &
            btype(Z_DIRECTION), grid%info%stencil_type, grid%info%NX, grid%info%NY, &
            grid%info%NZ, grid%info%nproc_x, grid%info%nproc_y, grid%info%nproc_z, &
            grid%da_sizes(ONEDOF), grid%info%stencil_size_rho, PETSC_NULL_INTEGER, &
            PETSC_NULL_INTEGER, PETSC_NULL_INTEGER,grid%da(ONEDOF), ierr)
       
       call DMDAGetInfo(grid%da(ONEDOF), PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
            PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, grid%info%nproc_x, &
            grid%info%nproc_y, grid%info%nproc_z, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
            PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
            PETSC_NULL_INTEGER, ierr)
       allocate(grid%info%ownership_x(1:grid%info%nproc_x))
       allocate(grid%info%ownership_y(1:grid%info%nproc_y))
       allocate(grid%info%ownership_z(1:grid%info%nproc_z))
       grid%info%ownership_x = 0
       grid%info%ownership_y = 0
       grid%info%ownership_z = 0
       call DMDAGetOwnershipRanges(grid%da(ONEDOF),grid%info%ownership_x, &
            grid%info%ownership_y,grid%info%ownership_z,ierr)
    end if

    ! create the rho DA seperately, since the stencil size (and therefore
    ! the ghost start/end/length) are potentially different
    call DMDACreate3d(grid%comm, btype(X_DIRECTION), btype(Y_DIRECTION), &
         btype(Z_DIRECTION), grid%info%stencil_type, &
         grid%info%NX, grid%info%NY, grid%info%NZ, &
         grid%info%nproc_x, grid%info%nproc_y, grid%info%nproc_z, &
         grid%da_sizes(NCOMPONENTDOF), grid%info%stencil_size_rho, &
         grid%info%ownership_x, grid%info%ownership_y, grid%info%ownership_z, &
         grid%da(NCOMPONENTDOF), ierr)

    ! create all other DAs
    do lcv=3,grid%nda
       call DMDACreate3d(grid%comm, btype(X_DIRECTION), btype(Y_DIRECTION), &
            btype(Z_DIRECTION), grid%info%stencil_type, &
            grid%info%NX, grid%info%NY, grid%info%NZ, &
            grid%info%nproc_x, grid%info%nproc_y, grid%info%nproc_z, &
            grid%da_sizes(lcv), grid%info%stencil_size, &
            grid%info%ownership_x, grid%info%ownership_y, grid%info%ownership_z, &
            grid%da(lcv), ierr)
    end do

    ! set grid%info including corners
    call DMDAGetCorners(grid%da(NCOMPONENTXBDOF), xs, ys, zs, grid%info%xl, grid%info%yl, &
         grid%info%zl, ierr)
    call DMDAGetGhostCorners(grid%da(NCOMPONENTXBDOF), gxs, gys, gzs, grid%info%gxl, grid%info%gyl,&
         grid%info%gzl, ierr)

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

    call DMDAGetGhostCorners(grid%da(NCOMPONENTDOF), rgxs, rgys, rgzs, &
         grid%info%rgxl, grid%info%rgyl, grid%info%rgzl, ierr)
    grid%info%rgxs = rgxs+1
    grid%info%rgxe = grid%info%rgxs+grid%info%rgxl-1

    grid%info%rgys = rgys+1
    grid%info%rgye = grid%info%rgys+grid%info%rgyl-1

    if (grid%info%ndims > 2) then
      grid%info%rgzs = rgzs+1
      grid%info%rgze = grid%info%rgzs+grid%info%rgzl-1
    end if

    grid%info%rgxyzl = grid%info%rgxl*grid%info%rgyl*grid%info%rgzl

    ! set the names
    call PetscObjectSetName(grid%da(ONEDOF), trim(grid%name)//'DA_one', ierr)
    call PetscObjectSetName(grid%da(NCOMPONENTDOF), trim(grid%name)//'DA_NCOMPONENT', ierr)
    call PetscObjectSetName(grid%da(NCOMPONENTXBDOF), trim(grid%name)//'DA_NCOMPONENTXB', ierr)
    call PetscObjectSetName(grid%da(NFLOWDOF), trim(grid%name)//'DA_NFLOW', ierr)
    if (grid%nda > 4) then
       call PetscObjectSetName(grid%da(NSPECIEDOF), &
            trim(grid%name)//'DA_NSPECIE', ierr)
       call PetscObjectSetName(grid%da(NSPECIEXBDOF), &
            trim(grid%name)//'DA_NSPECIEXB', ierr)
    end if

    ! set the coordinates, but only on the one-dof (no need for extra memory)
    if (grid%info%ndims > 2) then
       call DMDASetUniformCoordinates(grid%da(ONEDOF), &
            grid%info%corners(X_DIRECTION,1), grid%info%corners(X_DIRECTION,2), &
            grid%info%corners(Y_DIRECTION,1), grid%info%corners(Y_DIRECTION,2), &
            grid%info%corners(Z_DIRECTION,1), grid%info%corners(Z_DIRECTION,2), ierr)
       call DMDASetUniformCoordinates(grid%da(NFLOWDOF), &
            grid%info%corners(X_DIRECTION,1), grid%info%corners(X_DIRECTION,2), &
            grid%info%corners(Y_DIRECTION,1), grid%info%corners(Y_DIRECTION,2), &
            grid%info%corners(Z_DIRECTION,1), grid%info%corners(Z_DIRECTION,2), ierr)
    else
       call DMDASetUniformCoordinates(grid%da(ONEDOF), &
            grid%info%corners(X_DIRECTION,1), grid%info%corners(X_DIRECTION,2), &
            grid%info%corners(Y_DIRECTION,1), grid%info%corners(Y_DIRECTION,2), &
            zero, one, ierr)
       call DMDASetUniformCoordinates(grid%da(NFLOWDOF), &
            grid%info%corners(X_DIRECTION,1), grid%info%corners(X_DIRECTION,2), &
            grid%info%corners(Y_DIRECTION,1), grid%info%corners(Y_DIRECTION,2), &
            zero, one, ierr)
    end if
    CHKERRQ(ierr)
  end subroutine GridSetUp

  subroutine GridViewCoordinates(grid, io)
    use LBM_IO_module
    type(grid_type) grid
    type(io_type) io
    Vec coords
    PetscErrorCode ierr

    call DMGetCoordinates(grid%da(ONEDOF), coords, ierr)
    call IOView(io, coords, 'coords')
  end subroutine GridViewCoordinates
end module LBM_Grid_module
