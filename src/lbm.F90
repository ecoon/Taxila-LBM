!!!==================================================================
!!! Fortran-file
!!!    author:        Ethan T. Coon
!!!    filename:      lbregion.f
!!!    version:
!!!    created:       17 November 2010
!!!      on:          12:29:33 MST
!!!    last modified:  17 November 2010
!!!      at:          12:29:33 MST
!!!    URL:           http://www.ldeo.columbia.edu/~ecoon/
!!!    email:         ecoon _at_ ldeo.columbia.edu
!!!
!!!==================================================================
  ! module for region on which LBM will be used
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"
  
  module LBM_module
    use petsc
    use Timing_module
    use LBM_Info_module
    use LBM_Constants_module
    use LBM_BC_module
    use LBM_Options_module
    implicit none

    private

#include "lbm_definitions.h"
    type, public:: lbm_type
       MPI_Comm comm
       type(info_type),pointer:: info
       type(bc_type),pointer:: bc
       type(options_type),pointer:: options
       type(constants_type),pointer:: constants

       DM,pointer:: da_one ! pressure, rhot, etc
       DM,pointer:: da_s   ! rho -- #dofs = s = # of components
       DM,pointer:: da_sb  ! fi -- #dofs = s*b = components * directions
       DM,pointer:: da_flow ! ut -- #dofs = 3
       PetscInt dm_index_to_ndof(4)

       Vec fi
       Vec rho
       Vec ut
       Vec prs
       Vec rhot
       Vec walls

       Vec fi_g
       Vec rho_g
       Vec ut_g
       Vec uyt_g
       Vec uzt_g
       Vec prs_g
       Vec rhot_g
       Vec walls_g

       PetscScalar,pointer:: ut_a(:)
       PetscScalar,pointer:: prs_a(:)
       PetscScalar,pointer:: walls_a(:)
       PetscScalar,pointer:: rho_a(:)
       PetscScalar,pointer:: fi_a(:)
       PetscScalar,pointer:: rhot_a(:)

       PetscScalar,pointer,dimension(:,:,:):: vel
       PetscScalar,pointer,dimension(:,:,:):: vel_eq
       PetscScalar,pointer,dimension(:,:,:):: forces
       character(len=MAXSTRINGLENGTH) name       
    end type lbm_type

    interface LBMCreate
       module procedure LBMCreate_Comm
       module procedure LBMCreate_CommName
    end interface
    
    public :: LBMCreate, &
         LBMSetFromOptions, &
         LBMDestroy, &
         LBMSetDomain, &
         LBMRun, &
         LBMLocalToGlobal, &
         LBMInitializeWalls, &
         LBMInitializeWallsPetsc, &
         LBMInitializeState, &
         LBMGetDMByIndex, &
         LBMGetCorners, &
         LBMView, &
         LBMOutput, &
         LBMInput

  contains
    function LBMCreate_Comm(comm) result(lbm)
      type(lbm_type),pointer:: lbm
      MPI_Comm comm

      allocate(lbm)
      allocate(lbm%da_one)
      allocate(lbm%da_s)
      allocate(lbm%da_sb)
      allocate(lbm%da_flow)

      lbm%comm = 0
      lbm%da_one = 0
      lbm%da_s = 0
      lbm%da_sb = 0
      lbm%da_flow = 0
      lbm%info => InfoCreate()
      lbm%bc => BCCreate()
      lbm%options => OptionsCreate(comm)
      lbm%constants => ConstantsCreate()
      lbm%comm = comm

      lbm%fi = 0
      lbm%rho = 0
      lbm%ut = 0
      lbm%prs = 0
      lbm%rhot = 0
      lbm%walls = 0

      lbm%fi_g = 0
      lbm%rho_g = 0
      lbm%ut_g = 0
      lbm%prs_g = 0
      lbm%rhot_g = 0
      lbm%walls_g = 0

      nullify(lbm%ut_a)
      nullify(lbm%prs_a)
      nullify(lbm%walls_a)
      nullify(lbm%rho_a)
      nullify(lbm%fi_a)
      nullify(lbm%rhot_a)

      nullify(lbm%vel)
      nullify(lbm%vel_eq)
      nullify(lbm%forces)

      lbm%name = ''
    end function LBMCreate_Comm

    function LBMCreate_CommName(comm, name) result(lbm)
      type(lbm_type),pointer:: lbm 
      MPI_Comm comm
      character(len=*):: name       

      lbm => LBMCreate_Comm(comm)
      lbm%name = name
    end function LBMCreate_CommName

    ! --- set up LB method
    subroutine LBMSetFromOptions(lbm, options)
      type(lbm_type) lbm
      type(options_type) options

      ! local
      PetscErrorCode ierr
      PetscInt xs,gxs
      PetscInt ys,gys
      PetscInt zs,gzs
      PetscInt,allocatable,dimension(:):: ownership_x, ownership_y, ownership_z
      PetscScalar,dimension(3,2):: corners
      PetscScalar,dimension(3):: tmpcorners
      PetscInt nmax
      PetscBool flag, flag2
      integer charlen
      Vec coords
      character(len=MAXSTRINGLENGTH) coordfile
      PetscViewer viewer
      PetscScalar zero
      PetscInt periodicity
        
      call mpi_comm_rank(lbm%comm,lbm%info%id,ierr)
      call mpi_comm_size(lbm%comm,lbm%info%nproc, ierr)

      call InfoSetFromOptions(lbm%info, options, ierr)

      lbm%dm_index_to_ndof(ONEDOF) = 1
      lbm%dm_index_to_ndof(NPHASEDOF) = lbm%info%s
      lbm%dm_index_to_ndof(NPHASEXBDOF) = lbm%info%s*(lbm%info%flow_disc%b+1)
      lbm%dm_index_to_ndof(NFLOWDOF) = lbm%info%ndims

      ! create DAs
      periodicity = DMDA_NONPERIODIC
      if (lbm%info%periodic(X_DIRECTION)) then
         periodicity = IOR(periodicity,DMDA_XPERIODIC)
      else
         periodicity = IOR(periodicity,DMDA_XGHOSTED)
      end if

      if (lbm%info%periodic(Y_DIRECTION)) then
         periodicity = IOR(periodicity,DMDA_YPERIODIC)
      else
         periodicity = IOR(periodicity,DMDA_YGHOSTED)
      end if

      if (lbm%info%ndims > 2) then
         if (lbm%info%periodic(Z_DIRECTION)) then
            periodicity = IOR(periodicity,DMDA_ZPERIODIC)
         else
            periodicity = IOR(periodicity,DMDA_ZGHOSTED)
         end if
      end if

      call DMDACreate3d(lbm%comm, periodicity, DMDA_STENCIL_BOX, &
           lbm%info%NX, lbm%info%NY, lbm%info%NZ, &
           PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, lbm%dm_index_to_ndof(ONEDOF), 1, &
           PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, lbm%da_one, ierr)
      CHKERRQ(ierr)
      call PetscObjectSetName(lbm%da_one, trim(lbm%name)//'DA_one_dof', ierr)
!      call DMView(lbm%da_one, PETSC_VIEWER_STDOUT_WORLD, ierr)

      call DMDAGetCorners(lbm%da_one, xs, ys, zs, lbm%info%xl, lbm%info%yl, lbm%info%zl, ierr)
      call DMDAGetGhostCorners(lbm%da_one, gxs, gys, gzs, lbm%info%gxl, lbm%info%gyl, &
           lbm%info%gzl, ierr)
      call DMDAGetInfo(lbm%da_one, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
           PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, lbm%info%nproc_x, &
           lbm%info%nproc_y, lbm%info%nproc_z, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
           PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
      allocate(ownership_x(1:lbm%info%nproc_x))
      allocate(ownership_y(1:lbm%info%nproc_y))
      allocate(ownership_z(1:lbm%info%nproc_z))
      ownership_x = 0
      ownership_y = 0
      ownership_z = 0
      call DMDAGetOwnershipRanges(lbm%da_one, ownership_x, ownership_y, ownership_z, ierr)

      ! set lbm%info including corners
      lbm%info%xs = xs+1
      lbm%info%gxs = gxs+1
      lbm%info%xe = lbm%info%xs+lbm%info%xl-1
      lbm%info%gxe = lbm%info%gxs+lbm%info%gxl-1

      lbm%info%ys = ys+1
      lbm%info%gys = gys+1
      lbm%info%ye = lbm%info%ys+lbm%info%yl-1
      lbm%info%gye = lbm%info%gys+lbm%info%gyl-1

      if (lbm%info%ndims > 2) then
         lbm%info%zs = zs+1
         lbm%info%gzs = gzs+1
         lbm%info%ze = lbm%info%zs+lbm%info%zl-1
         lbm%info%gze = lbm%info%gzs+lbm%info%gzl-1
      end if

      lbm%info%xyzl = lbm%info%xl*lbm%info%yl*lbm%info%zl
      lbm%info%gxyzl = lbm%info%gxl*lbm%info%gyl*lbm%info%gzl

      call DMDACreate3d(lbm%comm, periodicity, DMDA_STENCIL_BOX, &
           lbm%info%NX, lbm%info%NY, lbm%info%NZ, &
           lbm%info%nproc_x, lbm%info%nproc_y, lbm%info%nproc_z, &
           lbm%dm_index_to_ndof(NPHASEDOF), 1, &
           ownership_x, ownership_y, ownership_z, &
           lbm%da_s, ierr)
      call PetscObjectSetName(lbm%da_s, trim(lbm%name)//'DA_NPHASE', ierr)

      call DMDACreate3d(lbm%comm, periodicity, DMDA_STENCIL_BOX, &
           lbm%info%NX, lbm%info%NY, lbm%info%NZ, &
           lbm%info%nproc_x, lbm%info%nproc_y, lbm%info%nproc_z, &
           lbm%dm_index_to_ndof(NPHASEXBDOF), 1, &
           ownership_x, ownership_y, ownership_z, &
           lbm%da_sb, ierr)
      call PetscObjectSetName(lbm%da_sb, trim(lbm%name)//'DA_fi', ierr)

      call DMDACreate3d(lbm%comm, periodicity, DMDA_STENCIL_BOX, &
           lbm%info%NX, lbm%info%NY, lbm%info%NZ, &
           lbm%info%nproc_x, lbm%info%nproc_y, lbm%info%nproc_z, &
           lbm%dm_index_to_ndof(NFLOWDOF), 1, &
           ownership_x, ownership_y, ownership_z, &
           lbm%da_flow, ierr)
      call PetscObjectSetName(lbm%da_flow, trim(lbm%name)//'DA_nflow', ierr)

      deallocate(ownership_x)
      deallocate(ownership_y)
      deallocate(ownership_z)

      ! set up constants, bcs, periodicity
      call ConstantsSetFromOptions(lbm%constants, lbm%info, options, ierr)
      call BCSetFromOptions(lbm%bc, lbm%info, options, ierr)

      ! coordinates
      tmpcorners = 0.d0
      corners = 0.d0

      nmax = lbm%info%ndims
      call PetscOptionsGetRealArray(options%my_prefix, '-corner0', tmpcorners, nmax, flag, ierr)
      if (flag) corners(:,1) = tmpcorners

      tmpcorners = 0.d0
      nmax = lbm%info%ndims
      call PetscOptionsGetRealArray(options%my_prefix, '-corner1', tmpcorners, nmax, flag2, ierr)
      if (flag2) corners(:,2) = tmpcorners

      if (flag .and. flag2) then
         call LBMSetDomain(lbm, corners)
         call DMDAGetCoordinates(lbm%da_one, coords, ierr)
         charlen = LEN_TRIM(lbm%options%output_prefix)
         coordfile = options%output_prefix(1:charlen)//'coords.dat'
         call PetscViewerBinaryOpen(PETSC_COMM_WORLD, coordfile, FILE_MODE_WRITE, &
              viewer, ierr)
         call VecView(coords, viewer, ierr)
         call PetscViewerDestroy(viewer,ierr)
      end if



      ! allocate, associate workspace
      allocate(lbm%vel(1:lbm%info%s, 1:lbm%info%ndims, 1:lbm%info%gxyzl))
      allocate(lbm%vel_eq(1:lbm%info%s, 1:lbm%info%ndims, 1:lbm%info%gxyzl))
      allocate(lbm%forces(1:lbm%info%s, 1:lbm%info%ndims, 1:lbm%info%gxyzl))

      lbm%vel = 0
      lbm%vel_eq = 0
      lbm%forces = 0

      ! get vectors
      call DMCreateLocalVector(lbm%da_one, lbm%prs, ierr)
      call VecDuplicate(lbm%prs, lbm%rhot, ierr)
      call VecDuplicate(lbm%prs, lbm%walls, ierr)

      call DMCreateLocalVector(lbm%da_s, lbm%rho, ierr)
      call DMCreateLocalVector(lbm%da_sb, lbm%fi, ierr)
      call DMCreateLocalVector(lbm%da_flow, lbm%ut, ierr)

      call DMCreateGlobalVector(lbm%da_one, lbm%prs_g, ierr)
      call VecDuplicate(lbm%prs_g, lbm%walls_g, ierr)
      call VecDuplicate(lbm%prs_g, lbm%rhot_g, ierr)
      call DMCreateGlobalVector(lbm%da_s, lbm%rho_g, ierr)
      call DMCreateGlobalVector(lbm%da_sb, lbm%fi_g, ierr)
      call DMCreateGlobalVector(lbm%da_flow, lbm%ut_g, ierr)

      zero = 0.d0
      call VecSet(lbm%prs_g, zero, ierr)
      call VecSet(lbm%walls_g, zero, ierr)
      call VecSet(lbm%rhot_g, zero, ierr)
      call VecSet(lbm%rho_g, zero, ierr)
      call VecSet(lbm%fi_g, zero, ierr)
      call VecSet(lbm%ut_g, zero, ierr)
      call VecSet(lbm%prs, zero, ierr)
      call VecSet(lbm%walls, zero, ierr)
      call VecSet(lbm%rhot, zero, ierr)
      call VecSet(lbm%rho, zero, ierr)
      call VecSet(lbm%fi, zero, ierr)
      call VecSet(lbm%ut, zero, ierr)

      call PetscObjectSetName(lbm%prs_g, trim(lbm%name)//'prs', ierr)
      call PetscObjectSetName(lbm%walls_g, trim(lbm%name)//'walls', ierr)
      call PetscObjectSetName(lbm%rhot_g, trim(lbm%name)//'rhot', ierr)
      call PetscObjectSetName(lbm%rho_g, trim(lbm%name)//'rho', ierr)
      call PetscObjectSetName(lbm%fi_g, trim(lbm%name)//'fi', ierr)
      call PetscObjectSetName(lbm%ut_g, trim(lbm%name)//'ut', ierr)

      CHKERRQ(ierr)
      CHKMEMQ
      return
    end subroutine LBMSetFromOptions

    subroutine LBMView(lbm)
      type(lbm_type) lbm
      
      call OptionsView(lbm%options)
      call InfoView(lbm%info)
      call ConstantsView(lbm%constants)
    end subroutine LBMView

    ! --- destroy things
    subroutine LBMDestroy(lbm, ierr)
      type(lbm_type) lbm
      PetscErrorCode ierr

      if (lbm%ut /= 0) call VecDestroy(lbm%ut,ierr)
      if (lbm%prs /= 0) call VecDestroy(lbm%prs,ierr)
      if (lbm%rhot /= 0) call VecDestroy(lbm%rhot,ierr)
      if (lbm%rho /= 0) call VecDestroy(lbm%rho,ierr)
      if (lbm%fi /= 0) call VecDestroy(lbm%fi,ierr)
      if (lbm%walls /= 0) call VecDestroy(lbm%walls,ierr)
      if (lbm%ut_g /= 0) call VecDestroy(lbm%ut_g,ierr)
      if (lbm%prs_g /= 0) call VecDestroy(lbm%prs_g,ierr)
      if (lbm%rhot_g /= 0) call VecDestroy(lbm%rhot_g,ierr)
      if (lbm%rho_g /= 0) call VecDestroy(lbm%rho_g,ierr)
      if (lbm%walls_g /= 0) call VecDestroy(lbm%walls_g,ierr)

      if (associated(lbm%vel)) deallocate(lbm%vel)
      if (associated(lbm%vel_eq)) deallocate(lbm%vel_eq)
      if (associated(lbm%forces)) deallocate(lbm%forces)

      if (lbm%da_one /= 0) call DMDestroy(lbm%da_one, ierr)
      if (lbm%da_s /= 0) call DMDestroy(lbm%da_s, ierr)
      if (lbm%da_sb /= 0) call DMDestroy(lbm%da_sb, ierr)
      if (lbm%da_flow /= 0) call DMDestroy(lbm%da_flow, ierr)
      deallocate(lbm%da_one)
      deallocate(lbm%da_s)
      deallocate(lbm%da_sb)
      deallocate(lbm%da_flow)

      call BCDestroy(lbm%bc, ierr)
      call InfoDestroy(lbm%info, ierr)
      call ConstantsDestroy(lbm%constants, ierr)
      return
    end subroutine LBMDestroy

    ! --- do things
    subroutine LBMSetDomain(lbm, corners)
      type(lbm_type) lbm
      PetscScalar,dimension(lbm%info%ndims,2):: corners
      PetscErrorCode ierr
      PetscScalar deltacoord
      PetscScalar newcoord_x
      PetscScalar newcoord_y
      PetscScalar newcoord_z

      ! uniform coordinates DA
      if (.not.lbm%info%periodic(X_DIRECTION)) then
         lbm%info%gridsize(X_DIRECTION) = (corners(X_DIRECTION,2) - &
              corners(X_DIRECTION,1))/dble(lbm%info%NX-1)
      else
         lbm%info%gridsize(X_DIRECTION) = (corners(X_DIRECTION,2) - &
              corners(X_DIRECTION,1))/dble(lbm%info%NX)
      endif

      if (.not.lbm%info%periodic(Y_DIRECTION)) then
         lbm%info%gridsize(Y_DIRECTION) = (corners(Y_DIRECTION,2) - &
              corners(Y_DIRECTION,1))/dble(lbm%info%NY-1)
      else
         lbm%info%gridsize(Y_DIRECTION) = (corners(Y_DIRECTION,2) - &
              corners(Y_DIRECTION,1))/dble(lbm%info%NY)
      endif

      if (lbm%info%ndims > 2) then
         if (.not.lbm%info%periodic(Z_DIRECTION)) then
            lbm%info%gridsize(Z_DIRECTION) = (corners(Z_DIRECTION,2) - &
                 corners(Z_DIRECTION,1))/dble(lbm%info%NZ-1)
         else
            lbm%info%gridsize(Z_DIRECTION) = (corners(Z_DIRECTION,2) - &
                 corners(Z_DIRECTION,1))/dble(lbm%info%NZ)
         endif

         call DMDASetUniformCoordinates(lbm%da_one, corners(X_DIRECTION,1), &
              corners(X_DIRECTION, 2), corners(Y_DIRECTION,1), corners(Y_DIRECTION, 2), &
              corners(Z_DIRECTION,1), corners(Z_DIRECTION, 2), ierr)
         call DMDASetUniformCoordinates(lbm%da_s, corners(X_DIRECTION,1), &
              corners(X_DIRECTION, 2), corners(Y_DIRECTION,1), corners(Y_DIRECTION, 2), &
              corners(Z_DIRECTION,1), corners(Z_DIRECTION, 2), ierr)
         call DMDASetUniformCoordinates(lbm%da_sb, corners(X_DIRECTION,1), &
              corners(X_DIRECTION, 2), corners(Y_DIRECTION,1), corners(Y_DIRECTION, 2), &
              corners(Z_DIRECTION,1), corners(Z_DIRECTION, 2), ierr)
         call DMDASetUniformCoordinates(lbm%da_flow, corners(X_DIRECTION,1), &
              corners(X_DIRECTION, 2), corners(Y_DIRECTION,1), corners(Y_DIRECTION, 2), &
              corners(Z_DIRECTION,1), corners(Z_DIRECTION, 2), ierr)
      else 
         call DMDASetUniformCoordinates(lbm%da_one, corners(X_DIRECTION,1), &
              corners(X_DIRECTION, 2), corners(Y_DIRECTION,1), corners(Y_DIRECTION, 2), &
              PETSC_NULL_SCALAR, PETSC_NULL_SCALAR, ierr)

         call DMDASetUniformCoordinates(lbm%da_s, corners(X_DIRECTION,1), &
              corners(X_DIRECTION, 2), corners(Y_DIRECTION,1), corners(Y_DIRECTION, 2), &
              PETSC_NULL_SCALAR, PETSC_NULL_SCALAR, ierr)

         call DMDASetUniformCoordinates(lbm%da_sb, corners(X_DIRECTION,1), &
              corners(X_DIRECTION, 2), corners(Y_DIRECTION,1), corners(Y_DIRECTION, 2), &
              PETSC_NULL_SCALAR, PETSC_NULL_SCALAR, ierr)

         call DMDASetUniformCoordinates(lbm%da_flow, corners(X_DIRECTION,1), &
              corners(X_DIRECTION, 2), corners(Y_DIRECTION,1), corners(Y_DIRECTION, 2), &
              PETSC_NULL_SCALAR, PETSC_NULL_SCALAR, ierr)
      end if

      CHKERRQ(ierr)
      CHKMEMQ

      lbm%info%corners = corners
      return
    end subroutine LBMSetDomain

    subroutine LBMRun(lbm, istep, kstep, kwrite)
      use LBM_Updates_module
      use LBM_Streaming_module
      use LBM_Bounceback_module
      use LBM_Forcing_module
      use LBM_Collision_module

      ! input
      type(lbm_type) lbm
      integer istep
      integer kstep
      integer kwrite

      ! local
      PetscErrorCode ierr
      logical,dimension(0:10):: bcs_done        ! flag for whether boundary condition
      integer lcv_sides, lcv_step
      type(timing_type),pointer:: timer1, timer2, timer3
      character(len=MAXWORDLENGTH) timerunits
      character(len=MAXWORDLENGTH) timername

      ! communicate to initialize
      call DMDALocalToLocalBegin(lbm%da_one, lbm%walls, INSERT_VALUES, lbm%walls, ierr)
      call DMDALocalToLocalEnd(lbm%da_one, lbm%walls, INSERT_VALUES, lbm%walls, ierr)

      call DMDALocalToLocalBegin(lbm%da_sb, lbm%fi, INSERT_VALUES, lbm%fi, ierr)
      call DMDALocalToLocalEnd(lbm%da_sb, lbm%fi, INSERT_VALUES, lbm%fi, ierr)

      if (istep.eq.0) then
         ! get arrays
         call DMDAVecGetArrayF90(lbm%da_one, lbm%rhot, lbm%rhot_a, ierr)
         call DMDAVecGetArrayF90(lbm%da_one, lbm%prs, lbm%prs_a, ierr)
         call DMDAVecGetArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
         call DMDAVecGetArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
         call DMDAVecGetArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
         call DMDAVecGetArrayF90(lbm%da_flow, lbm%ut, lbm%ut_a, ierr)
         
         ! update values for zero time i/o
         call LBMUpdateMoments(lbm%fi_a,lbm%rho_a, lbm%vel, lbm%walls_a, lbm%info)
         call LBMUpdateDiagnostics(lbm%rho_a, lbm%vel, lbm%walls_a, lbm%ut_a, &
              lbm%rhot_a, lbm%prs_a, lbm%forces, lbm%info, lbm%constants)
         
         ! --  --  restore arrays
         call DMDAVecRestoreArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
         call DMDAVecRestoreArrayF90(lbm%da_one, lbm%prs, lbm%prs_a, ierr)
         call DMDAVecRestoreArrayF90(lbm%da_one, lbm%rhot, lbm%rhot_a, ierr)
         call DMDAVecRestoreArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
         call DMDAVecRestoreArrayF90(lbm%da_flow, lbm%ut, lbm%ut_a, ierr)
         call DMDAVecRestoreArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)

         call LBMLocalToGlobal(lbm, ierr)
         ! output at zero time
         call LBMOutput(lbm, istep, kwrite)
      else
         call LBMLocalToGlobal(lbm, ierr)
      endif

      ! get arrays
      call DMDAVecGetArrayF90(lbm%da_one, lbm%rhot, lbm%rhot_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_one, lbm%prs, lbm%prs_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_flow, lbm%ut, lbm%ut_a, ierr)

      call BCGetArrays(lbm%bc, ierr)


      timername = trim(lbm%name)//'Simulation'
      timer1 => TimingCreate(lbm%comm, timername)
      do lcv_step = istep+1,kstep
!         timername = 'Streaming'
!         timer2 => TimingCreate(lbm%comm, timername)
         call LBMStreaming(lbm%fi_a, lbm%info)
!         call TimingEnd(timer2)
!         timer2%name = 'Streaming + Bounceback'
         call LBMBounceback(lbm%fi_a, lbm%walls_a, lbm%info)
!         call TimingEnd(timer2)
!         timername = 'BCs + moments'
!         timer2 => TimingCreate(lbm%comm, timername)

         bcs_done=.FALSE.
         bcs_done(0) = .TRUE.   ! periodic done by default
         do lcv_sides = 1,6
            if (.not.bcs_done(lbm%bc%flags(lcv_sides))) then
               select case (lbm%bc%flags(lcv_sides))
               case (BC_PSEUDOPERIODIC)         ! pseudo-periodic
                  call BCPseudoperiodic(lbm%bc, lbm%fi_a, lbm%walls_a, lbm%info)

               case (BC_FLUX)         ! flux
                  call BCFlux(lbm%bc, lbm%fi_a, lbm%walls_a, lbm%info)

               case (BC_PRESSURE)         ! pressure
                  call BCPressure(lbm%bc, lbm%fi_a, lbm%walls_a, lbm%info)
               end select
               bcs_done(lbm%bc%flags(lcv_sides)) = .TRUE. ! only do each bc type once
            endif
         enddo

         call LBMUpdateMoments(lbm%fi_a,lbm%rho_a, lbm%vel, lbm%walls_a, lbm%info)

!         call TimingEnd(timer2)
!         timername = 'communication'
!         timer2 => TimingCreate(lbm%comm, timername)

         ! update rho ghosts values
         call DMDAVecRestoreArrayF90(lbm%da_s,lbm%rho, lbm%rho_a, ierr)
         call DMDALocalToLocalBegin(lbm%da_s, lbm%rho, INSERT_VALUES, lbm%rho, ierr)
         call DMDALocalToLocalEnd(lbm%da_s, lbm%rho, INSERT_VALUES, lbm%rho, ierr)
         call DMDAVecGetArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
!         call TimingEnd(timer2)
!         timer2%name = 'communication + forces'
!         timername = 'forces'
!         timer3 => TimingCreate(lbm%comm, timername)


         !calculate forces
         lbm%forces=0.d0

         if (lbm%info%s > 1) then
            if ((lbm%constants%g /= 0.0d0).or.(lbm%constants%g11 /= 0.0d0).or. &
                 (lbm%constants%g22 /= 0.0d0)) then
               call LBMAddFluidFluidForces(lbm%rho_a, lbm%forces, &
                    lbm%walls_a, lbm%info, lbm%constants)
            end if

            if ((lbm%constants%gw(1) /= 0.0d0).or.(lbm%constants%gw(2) /= 0.0d0)) then
               call LBMAddFluidSolidForces(lbm%rho_a, lbm%forces, &
                    lbm%walls_a, lbm%info, lbm%constants)
            end if
         endif
         call LBMAddBodyForces(lbm%rho_a, lbm%forces, lbm%walls_a, lbm%info, lbm%constants)
         call LBMZeroBoundaryForces(lbm%bc%flags, lbm%forces, lbm%info%ndims, &
              lbm%info, lbm%constants)

!         call TimingEnd(timer2)
!         call TimingEnd(timer3)
!         timername = 'equilibrium'
!         timer2 => TimingCreate(lbm%comm, timername)

         ! calculate u_equilibrium
         call LBMUpdateUEquilibrium(lbm%fi_a, lbm%rho_a, lbm%vel, lbm%walls_a, &
              lbm%vel_eq, lbm%rhot_a, lbm%forces, lbm%info, lbm%constants)
!         call TimingEnd(timer2)
!         timername = 'collision'
!         timer2 => TimingCreate(lbm%comm, timername)

         ! collision
         call LBMCollision(lbm%fi_a, lbm%rho_a, lbm%vel_eq, lbm%walls_a, &
              lbm%info, lbm%constants)
!         call TimingEnd(timer2)

         ! communicate, update fi
         call DMDAVecRestoreArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
         call DMDALocalToLocalBegin(lbm%da_sb, lbm%fi, INSERT_VALUES, lbm%fi, ierr)
         call DMDALocalToLocalEnd(lbm%da_sb, lbm%fi, INSERT_VALUES, lbm%fi, ierr)

         ! check for output?
         if(mod(lcv_step,kwrite).eq.0) then
            ! --  --  update diagnostics
            call LBMUpdateDiagnostics(lbm%rho_a, lbm%vel, lbm%walls_a, lbm%ut_a, &
                 lbm%rhot_a, lbm%prs_a, lbm%forces, lbm%info, lbm%constants)

            ! --  --  restore arrays
            call DMDAVecRestoreArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%da_one, lbm%prs, lbm%prs_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%da_one, lbm%rhot, lbm%rhot_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%da_flow, lbm%ut, lbm%ut_a, ierr)

            ! --  --  output
            call LBMLocalToGlobal(lbm, ierr)
            call LBMOutput(lbm, lcv_step, kwrite)

            ! --  --  reopen arrays
            call DMDAVecGetArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
            call DMDAVecGetArrayF90(lbm%da_one, lbm%prs, lbm%prs_a, ierr)
            call DMDAVecGetArrayF90(lbm%da_one, lbm%rhot, lbm%rhot_a, ierr)
            call DMDAVecGetArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
            call DMDAVecGetArrayF90(lbm%da_flow, lbm%ut, lbm%ut_a, ierr)
         endif

         call DMDAVecGetArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
      end do

      timerunits = 'timestep'
      call TimingEndPerUnit(timer1, (kstep-istep+1), timerunits)
      call TimingDestroy(timer1)

      ! restore arrays in prep for communication
      call BCRestoreArrays(lbm%bc, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_one, lbm%prs, lbm%prs_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_one, lbm%rhot, lbm%rhot_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_flow, lbm%ut, lbm%ut_a, ierr)

      ! communicate local to global
      call LBMLocalToGlobal(lbm, ierr)
      return
    end subroutine LBMRun

    subroutine LBMLocalToGlobal(lbm, ierr)
      type(lbm_type) lbm
      PetscErrorCode ierr

      call DMLocalToGlobalBegin(lbm%da_flow, lbm%ut, INSERT_VALUES, lbm%ut_g, ierr)
      call DMLocalToGlobalEnd(lbm%da_flow, lbm%ut, INSERT_VALUES, lbm%ut_g, ierr)

      call DMLocalToGlobalBegin(lbm%da_one, lbm%prs, INSERT_VALUES, lbm%prs_g, ierr)
      call DMLocalToGlobalEnd(lbm%da_one, lbm%prs, INSERT_VALUES, lbm%prs_g, ierr)

      call DMLocalToGlobalBegin(lbm%da_one, lbm%walls, INSERT_VALUES, lbm%walls_g, ierr)
      call DMLocalToGlobalEnd(lbm%da_one, lbm%walls, INSERT_VALUES, lbm%walls_g, ierr)

      call DMLocalToGlobalBegin(lbm%da_one, lbm%rhot, INSERT_VALUES, lbm%rhot_g, ierr)
      call DMLocalToGlobalEnd(lbm%da_one, lbm%rhot, INSERT_VALUES, lbm%rhot_g, ierr)

      call DMLocalToGlobalBegin(lbm%da_s, lbm%rho, INSERT_VALUES, lbm%rho_g, ierr)
      call DMLocalToGlobalEnd(lbm%da_s, lbm%rho, INSERT_VALUES, lbm%rho_g, ierr)

      call DMLocalToGlobalBegin(lbm%da_sb, lbm%fi, INSERT_VALUES, lbm%fi_g, ierr)
      call DMLocalToGlobalEnd(lbm%da_sb, lbm%fi, INSERT_VALUES, lbm%fi_g, ierr)
      return
    end subroutine LBMLocalToGlobal

    subroutine LBMInitializeWalls(lbm, init_subroutine)
      type(lbm_type) lbm
      !      interface
      !         subroutine init_subroutine(walls, info)
      !           use LBM_Info_module
      !           type(info_type) info
      !           PetscScalar walls(:) ! problem is that this does not work, as we want
      !                                  to specify an explicit shape within the subroutine
      !         end subroutine init_subroutine
      !      end interface
      external :: init_subroutine
      PetscErrorCode ierr
      PetscInt vsize

      call DMDAVecGetArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
      call init_subroutine(lbm%walls_a, lbm%options%walls_file, lbm%info)
      call DMDAVecRestoreArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
      return
    end subroutine LBMInitializeWalls

    subroutine LBMInitializeWallsPetsc(lbm, filename)
      type(lbm_type) lbm
      character(len=MAXSTRINGLENGTH) filename
      
      PetscViewer viewer
      PetscErrorCode ierr
      call PetscViewerBinaryOpen(lbm%comm, filename, FILE_MODE_READ, viewer, ierr)
      call VecLoad(lbm%walls_g, viewer, ierr)
      call PetscViewerDestroy(viewer, ierr)
      call DMGlobalToLocalBegin(lbm%da_one, lbm%walls_g, INSERT_VALUES, lbm%walls, ierr)
      call DMGlobalToLocalEnd(lbm%da_one, lbm%walls_g, INSERT_VALUES, lbm%walls, ierr)
      return
    end subroutine LBMInitializeWallsPetsc

    subroutine LBMInitializeState(lbm, init_subroutine)
      type(lbm_type) lbm
      !      interface
      !         subroutine init_subroutine(fi, rho, vel, walls, info)
      !           use LBM_Info_module
      !           type(info_type) info
      !           PetscScalar fi(:)
      !           PetscScalar rho(:)
      !           PetscScalar vel(:,:,:,:,:)
      !           PetscScalar walls(:)
      !         end subroutine init_subroutine
      !      end interface
      external :: init_subroutine
      PetscErrorCode ierr
      call DMDAVecGetArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
      call init_subroutine(lbm%fi_a, lbm%rho_a, lbm%vel, lbm%walls_a, lbm%info, lbm%constants)
      call DMDAVecRestoreArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
      return
    end subroutine LBMInitializeState

    function LBMGetDMByIndex( lbm, dm_index ) result(dm)
      type(lbm_type) lbm
      PetscInt dm_index
      DM,pointer:: dm
      PetscErrorCode ierr
      nullify(dm)
      select case(dm_index)
      case (ONEDOF)
         dm => lbm%da_one
      case (NPHASEDOF)
         dm => lbm%da_s
      case (NPHASEXBDOF)
         dm => lbm%da_sb
      case (NFLOWDOF)
         dm => lbm%da_flow
      end select
    end function LBMGetDMByIndex

    ! get the corner coordinates of the domain
    subroutine LBMGetCorners(lbm, corners)
      type(lbm_type) lbm
      PetscScalar,dimension(3,2):: corners
      
      corners = lbm%info%corners
    end subroutine LBMGetCorners

  subroutine LBMOutput(lbm, istep, kwrite)
    ! input variables
    type(lbm_type) lbm
    PetscInt istep, kwrite

    ! local variables
    PetscViewer viewer
    PetscErrorCode ierr
    character(len=MAXIODIGITS):: outnum
    character(len=MAXSTRINGLENGTH):: stringformat
    character(len=MAXSTRINGLENGTH):: flnm1,flnm2,flnm3,flnm4,flnm5,flnm6,flnm7,flnm8
    integer charlen

    charlen = LEN_TRIM(lbm%options%output_prefix)

    if (lbm%info%id.eq.0) then
       write(*,*) 'outputing (istep, kwrite)', istep, kwrite
    endif

    charlen = LEN_TRIM(lbm%options%output_prefix)

    write(stringformat, '("(I0.",I1,")")'), MAXIODIGITS
    write(outnum, stringformat) istep/kwrite
    flnm1=lbm%options%output_prefix(1:charlen)//'fi'//outnum//'.dat'
    flnm2=lbm%options%output_prefix(1:charlen)//'rho'//outnum//'.dat'
    flnm3=lbm%options%output_prefix(1:charlen)//'u'//outnum//'.dat'
    flnm6=lbm%options%output_prefix(1:charlen)//'walls'//outnum//'.dat'
    flnm7=lbm%options%output_prefix(1:charlen)//'rhot'//outnum//'.dat'
    flnm8=lbm%options%output_prefix(1:charlen)//'prs'//outnum//'.dat'

    call PetscViewerCreate(lbm%comm, viewer, ierr)
    call PetscViewerSetType(viewer, PETSCVIEWERBINARY, ierr)
    call PetscViewerBinarySetMPIIO(viewer, ierr)
    call PetscViewerFileSetMode(viewer, FILE_MODE_WRITE, ierr)
    call PetscViewerFileSetName(viewer, flnm1, ierr)
    call VecView(lbm%fi_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)

    call PetscViewerCreate(lbm%comm, viewer, ierr)
    call PetscViewerSetType(viewer, PETSCVIEWERBINARY, ierr)
    call PetscViewerBinarySetMPIIO(viewer, ierr)
    call PetscViewerFileSetMode(viewer, FILE_MODE_WRITE, ierr)
    call PetscViewerFileSetName(viewer, flnm2, ierr)
    call VecView(lbm%rho_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)

    call PetscViewerCreate(lbm%comm, viewer, ierr)
    call PetscViewerSetType(viewer, PETSCVIEWERBINARY, ierr)
    call PetscViewerBinarySetMPIIO(viewer, ierr)
    call PetscViewerFileSetMode(viewer, FILE_MODE_WRITE, ierr)
    call PetscViewerFileSetName(viewer, flnm3, ierr)
    call VecView(lbm%ut_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)

    call PetscViewerCreate(lbm%comm, viewer, ierr)
    call PetscViewerSetType(viewer, PETSCVIEWERBINARY, ierr)
    call PetscViewerBinarySetMPIIO(viewer, ierr)
    call PetscViewerFileSetMode(viewer, FILE_MODE_WRITE, ierr)
    call PetscViewerFileSetName(viewer, flnm6, ierr)
    call VecView(lbm%walls_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)

    call PetscViewerCreate(lbm%comm, viewer, ierr)
    call PetscViewerSetType(viewer, PETSCVIEWERBINARY, ierr)
    call PetscViewerBinarySetMPIIO(viewer, ierr)
    call PetscViewerFileSetMode(viewer, FILE_MODE_WRITE, ierr)
    call PetscViewerFileSetName(viewer, flnm7, ierr)
    call VecView(lbm%rhot_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)

    call PetscViewerCreate(lbm%comm, viewer, ierr)
    call PetscViewerSetType(viewer, PETSCVIEWERBINARY, ierr)
    call PetscViewerBinarySetMPIIO(viewer, ierr)
    call PetscViewerFileSetMode(viewer, FILE_MODE_WRITE, ierr)
    call PetscViewerFileSetName(viewer, flnm8, ierr)
    call VecView(lbm%prs_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)

    return
  end subroutine LBMOutput

  subroutine LBMInput(lbm, istep, kwrite)

    ! input variables
    type(lbm_type) lbm
    integer istep, kwrite

    ! local variables
    PetscViewer viewer
    PetscErrorCode ierr
    character(len=MAXIODIGITS):: outnum
    character(len=MAXSTRINGLENGTH):: stringformat
    character(len=MAXSTRINGLENGTH):: flnm1,flnm2,flnm3,flnm4,flnm5,flnm6,flnm7,flnm8
    integer charlen

    charlen = LEN_TRIM(lbm%options%output_prefix)

    if (lbm%info%id.eq.0) then
       write(*,*) 'outputing (istep, kwrite)', istep, kwrite
    endif

    charlen = LEN_TRIM(lbm%options%output_prefix)

    write(stringformat, '("(I0.",I1)'), MAXIODIGITS
    write(outnum, stringformat) istep/kwrite
    flnm1=lbm%options%output_prefix(1:charlen)//'fi'//outnum//'.dat'
    flnm2=lbm%options%output_prefix(1:charlen)//'rho'//outnum//'.dat'
    flnm3=lbm%options%output_prefix(1:charlen)//'u'//outnum//'.dat'
    flnm6=lbm%options%output_prefix(1:charlen)//'walls'//outnum//'.dat'
    flnm7=lbm%options%output_prefix(1:charlen)//'rhot'//outnum//'.dat'
    flnm8=lbm%options%output_prefix(1:charlen)//'prs'//outnum//'.dat'

    call PetscViewerBinaryOpen(lbm%comm, flnm1, FILE_MODE_READ, viewer, ierr)
    call VecLoad(lbm%fi_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)

    call PetscViewerBinaryOpen(lbm%comm, flnm2, FILE_MODE_READ, viewer, ierr)
    call VecLoad(lbm%rho_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)

    call PetscViewerBinaryOpen(lbm%comm, flnm3, FILE_MODE_READ, viewer, ierr)
    call VecLoad(lbm%ut_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)

    call PetscViewerBinaryOpen(lbm%comm, flnm6, FILE_MODE_READ, viewer, ierr)
    call VecLoad(lbm%walls_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)

    call PetscViewerBinaryOpen(lbm%comm, flnm7, FILE_MODE_READ, viewer, ierr)
    call VecLoad(lbm%rhot_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)

    call PetscViewerBinaryOpen(lbm%comm, flnm8, FILE_MODE_READ, viewer, ierr)
    call VecLoad(lbm%prs_g, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)
  end subroutine LBMInput
end module LBM_module

