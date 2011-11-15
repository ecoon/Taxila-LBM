!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_walls.F90
!!!     version:         
!!!     created:         21 June 2011
!!!       on:            10:07:45 MDT
!!!     last modified:   14 November 2011
!!!       at:            18:22:56 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"

module LBM_Walls_module
  use petsc
  use LBM_Mineral_module
  use LBM_Grid_module
  use LBM_Info_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public:: walls_type
    MPI_Comm comm
    character(len=MAXWORDLENGTH):: name       

    PetscInt ndims
    type(grid_type),pointer:: grid

    PetscInt nminerals
    type(mineral_type),pointer,dimension(:):: minerals
    
    Vec walls
    Vec walls_g
    PetscScalar,pointer:: walls_a(:)
    character(len=MAXSTRINGLENGTH):: filename
    PetscInt walls_type
  end type walls_type

  public :: WallsCreate, &
       WallsDestroy, &
       WallsSetGrid, &
       WallsSetFromOptions, &
       WallsSetUp, &
       WallsCommunicate, &
       WallsOutputDiagnostics
  
contains
  function WallsCreate(comm) result(walls)
    MPI_Comm comm
    type(walls_type),pointer :: walls
  
    allocate(walls)
    walls%comm = comm
    walls%ndims = -1
    walls%nminerals = -1
    nullify(walls%minerals)
    nullify(walls%grid)
    walls%walls = 0
    walls%walls_g = 0
    nullify(walls%walls_a)
    walls%name = ''
    walls%filename = 'geometry.dat'
    walls%walls_type = WALLS_TYPE_PETSC
  end function WallsCreate

  subroutine WallsDestroy(walls, ierr)
    type(walls_type) walls
    PetscErrorCode ierr
    PetscInt lcv

    if (associated(walls%minerals)) then
      do lcv=1,walls%nminerals
        call MineralDestroy(walls%minerals(lcv), ierr)
      end do
    end if
    if (walls%walls /= 0) call VecDestroy(walls%walls, ierr)
    if (walls%walls_g /= 0) call VecDestroy(walls%walls_g, ierr)
  end subroutine WallsDestroy

  subroutine WallsSetGrid(walls, grid)
    type(walls_type) walls
    type(grid_type),pointer:: grid
    walls%grid => grid
  end subroutine WallsSetGrid
    
  subroutine WallsSetName(walls, name) 
    type(walls_type) walls 
    character(len=MAXWORDLENGTH):: name       
    walls%name = name
  end subroutine WallsSetName

  subroutine WallsSetFromOptions(walls, options, ierr)
    use LBM_Options_module
    type(walls_type) walls
    type(options_type) options
    PetscErrorCode ierr

    PetscScalar,parameter:: eps=1.e-15 ! slightly larger than machine epsilon
    PetscInt lcv, lcv2
    PetscBool help
    PetscBool flag

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    walls%ndims = options%ndims
    walls%nminerals = options%nminerals
    walls%minerals => MineralCreate(walls%comm, walls%nminerals)
    
    do lcv=1,walls%nminerals
      call MineralSetID(walls%minerals(lcv), lcv)
      call MineralSetFromOptions(walls%minerals(lcv), options, ierr)
      do lcv2=1,options%ncomponents
        if (ABS(walls%minerals(lcv)%gw(lcv2)) > eps) then
          options%flow_fluidsolid_forces = PETSC_TRUE
        end if
      end do
    end do

    if (help) then
      call PetscPrintf(options%comm, "  -walls_file=<geometry.dat>: filename for porescale walls \n", ierr)
    end if
    call PetscOptionsGetString(options%my_prefix, '-walls_file', walls%filename, &
         flag, ierr)

    if (help) call PetscPrintf(options%comm, "  -walls_type <1>: enum: (1) use "// &
         "PETSc .dat file, (2) from initialize_walls subroutine\n", ierr)
    call PetscOptionsGetInt(options%my_prefix, '-walls_type', walls%walls_type, flag, &
         ierr)
  end subroutine WallsSetFromOptions

  subroutine WallsSetUp(walls)
    type(walls_type) walls
    PetscInt lcv
    PetscErrorCode ierr
    PetscViewer viewer
    external initialize_walls

    call DMCreateGlobalVector(walls%grid%da(ONEDOF), walls%walls_g, ierr)
    call DMCreateLocalVector(walls%grid%da(ONEDOF), walls%walls, ierr)
    call VecSet(walls%walls_g, 0.d0, ierr)
    call VecSet(walls%walls, 0.d0, ierr)
    call PetscObjectSetName(walls%walls_g, trim(walls%name)//'walls', ierr)

    if (walls%walls_type == WALLS_TYPE_PETSC) then
      call PetscViewerBinaryOpen(walls%comm, walls%filename, FILE_MODE_READ, viewer, ierr)
      call VecLoad(walls%walls_g, viewer, ierr)
      call PetscViewerDestroy(viewer, ierr)
      call DMGlobalToLocalBegin(walls%grid%da(ONEDOF), walls%walls_g, INSERT_VALUES, &
           walls%walls, ierr)
      call DMGlobalToLocalEnd(walls%grid%da(ONEDOF), walls%walls_g, INSERT_VALUES, &
           walls%walls, ierr)
      call DMDAVecGetArrayF90(walls%grid%da(ONEDOF), walls%walls, walls%walls_a, ierr)
      call WallsSetGhostNodes(walls, walls%walls_a)
    else
      call DMDAVecGetArrayF90(walls%grid%da(ONEDOF), walls%walls, walls%walls_a, ierr)
      call initialize_walls(walls%walls_a, walls%filename, walls%grid%info)
      call DMDAVecRestoreArrayF90(walls%grid%da(ONEDOF), walls%walls, &
           walls%walls_a, ierr)
      call DMLocalToGlobalBegin(walls%grid%da(ONEDOF), walls%walls, INSERT_VALUES, &
           walls%walls_g, ierr)
      call DMLocalToGlobalEnd(walls%grid%da(ONEDOF), walls%walls, INSERT_VALUES, &
           walls%walls_g, ierr)
      call DMDAVecGetArrayF90(walls%grid%da(ONEDOF), walls%walls, walls%walls_a, ierr)
      call WallsSetGhostNodes(walls, walls%walls_a)
    end if
  end subroutine WallsSetUp

  subroutine WallsSetGhostNodes(walls, data)
    type(walls_type) walls
    PetscScalar,dimension(walls%grid%info%rgxyzl):: data
    select case(walls%ndims)
    case (2)
      call WallsSetGhostNodesD2(walls, data, walls%grid%info)
    case (3)
      call WallsSetGhostNodesD3(walls, data, walls%grid%info)
    end select
  end subroutine WallsSetGhostNodes
  
  subroutine WallsSetGhostNodesD2(walls, data, info)
    type(walls_type) walls
    type(info_type) info
    PetscScalar,dimension(info%rgxs:info%rgxe,info%rgys:info%rgye):: data
    if ((info%xs.eq.1).and.(.not.info%periodic(X_DIRECTION))) then
      data(info%rgxs:info%xs-1,:) = 999.
    end if
    if ((info%xe.eq.info%NX).and.(.not.info%periodic(X_DIRECTION))) then
      data(info%xe+1:info%rgxe,:) = 999.
    end if
    if ((info%ys.eq.1).and.(.not.info%periodic(Y_DIRECTION))) then
      data(:,info%rgys:info%ys-1) = 999.
    end if
    if ((info%ye.eq.info%NY).and.(.not.info%periodic(Y_DIRECTION))) then
      data(:,info%ye+1:info%rgye) = 999.
    end if
  end subroutine WallsSetGhostNodesD2
  
  subroutine WallsSetGhostNodesD3(walls, data, info)
    type(walls_type) walls
    type(info_type) info
    PetscScalar,dimension(info%rgxs:info%rgxe,info%rgys:info%rgye, &
         info%rgzs:info%rgze):: data
    if ((info%xs.eq.1).and.(.not.info%periodic(X_DIRECTION))) then
      data(info%rgxs:info%xs-1,:,:) = 999.
    end if
    if ((info%xe.eq.info%NX).and.(.not.info%periodic(X_DIRECTION))) then
      data(info%xe+1:info%rgxe,:,:) = 999.
    end if
    if ((info%ys.eq.1).and.(.not.info%periodic(Y_DIRECTION))) then
      data(:,info%rgys:info%ys-1,:) = 999.
    end if
    if ((info%ye.eq.info%NY).and.(.not.info%periodic(Y_DIRECTION))) then
      data(:,info%ye+1:info%rgye,:) = 999.
    end if
    if ((info%zs.eq.1).and.(.not.info%periodic(Z_DIRECTION))) then
      data(:,:,info%rgzs:info%zs-1) = 999.
    end if
    if ((info%ze.eq.info%NZ).and.(.not.info%periodic(Z_DIRECTION))) then
      data(:,:,info%ze+1:info%rgze) = 999.
    end if
  end subroutine WallsSetGhostNodesD3

  subroutine WallsCommunicate(walls)
    type(walls_type) walls

    PetscErrorCode ierr

    call DMDAVecRestoreArrayF90(walls%grid%da(ONEDOF), walls%walls, walls%walls_a, ierr)
    call DMDALocalToLocalBegin(walls%grid%da(ONEDOF), walls%walls, INSERT_VALUES, &
         walls%walls, ierr)
    call DMDALocalToLocalEnd(walls%grid%da(ONEDOF), walls%walls, INSERT_VALUES, &
           walls%walls, ierr)
    call DMDAVecGetArrayF90(walls%grid%da(ONEDOF), walls%walls, walls%walls_a, ierr)
  end subroutine WallsCommunicate

  subroutine WallsOutputDiagnostics(walls, io)
    use LBM_IO_module
    type(walls_type) walls
    type(io_type) io
    PetscErrorCode ierr

    call DMDAVecRestoreArrayF90(walls%grid%da(ONEDOF), walls%walls, &
         walls%walls_a, ierr)
    call DMLocalToGlobalBegin(walls%grid%da(ONEDOF), walls%walls, INSERT_VALUES, &
         walls%walls_g, ierr)
    call DMLocalToGlobalEnd(walls%grid%da(ONEDOF), walls%walls, INSERT_VALUES, &
         walls%walls_g, ierr)
    call IOView(io, walls%walls_g, 'walls')
    call DMDAVecGetArrayF90(walls%grid%da(ONEDOF), walls%walls, walls%walls_a, ierr)
  end subroutine WallsOutputDiagnostics
end module LBM_Walls_module
