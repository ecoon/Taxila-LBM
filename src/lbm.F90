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
    use Info_module
    use BC_module
    use Timing_module
    use Options_module
    implicit none

    private

#include "lbm_definitions.h"
    type, public:: lbm_type
       MPI_Comm comm
       type(info_type),pointer:: info
       type(bc_type),pointer:: bc
       type(options_type),pointer:: options

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

       PetscScalar,pointer,dimension(:,:,:,:):: ux,uy,uz
       PetscScalar,pointer,dimension(:,:,:,:):: uxe,uye,uze
       PetscScalar,pointer,dimension(:,:,:,:):: Fx,Fy,Fz

       PetscScalar,dimension(3,2):: corners
    end type lbm_type

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
         LBMGetCorners

  contains
    function LBMCreate(comm) result(lbm)
      implicit none
      type(lbm_type),pointer:: lbm
      MPI_Comm comm
      PetscErrorCode ierr

      allocate(lbm)
      allocate(lbm%da_one)
      allocate(lbm%da_s)
      allocate(lbm%da_sb)
      allocate(lbm%da_flow)

      lbm%da_one = 0
      lbm%da_s = 0
      lbm%da_sb = 0
      lbm%da_flow = 0
      lbm%info => InfoCreate()
      lbm%bc => BCCreate()
      lbm%options => OptionsCreate()
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

      nullify(lbm%ux)
      nullify(lbm%uy)
      nullify(lbm%uz)

      nullify(lbm%uxe)
      nullify(lbm%uye)
      nullify(lbm%uze)

      nullify(lbm%Fx)
      nullify(lbm%Fy)
      nullify(lbm%Fz)
    end function LBMCreate

    ! --- set up LB method
    subroutine LBMSetFromOptions(lbm, options)
      implicit none

      ! input
      type(lbm_type) lbm
      type(options_type) options

      ! local
      integer xs,ys,zs,gxs,gys,gzs
      PetscErrorCode ierr

      call mpi_comm_rank(lbm%comm,lbm%info%id,ierr)
      call mpi_comm_size(lbm%comm,lbm%info%nproc, ierr)

      lbm%info%dim = 3
      lbm%info%s = options%s
      lbm%info%b = options%b
      lbm%dm_index_to_ndof(ONEDOF) = 1
      lbm%dm_index_to_ndof(NPHASEDOF) = options%s
      lbm%dm_index_to_ndof(NPHASEXBDOF) = options%s*(options%b+1)
      lbm%dm_index_to_ndof(NFLOWDOF) = lbm%info%dim

      lbm%info%NX = options%NX
      lbm%info%NY = options%NY
      lbm%info%NZ = options%NZ

      lbm%info%options_prefix = options%my_prefix
      lbm%bc%flags(:) = options%bc_flags(:)

      ! create DAs
      call DMDACreate3d(lbm%comm, DMDA_XYZPERIODIC, DMDA_STENCIL_BOX, &
           options%NX, options%NY, options%NZ, &
           PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, lbm%dm_index_to_ndof(ONEDOF), 1, &
           PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, lbm%da_one, ierr)
      CHKERRQ(ierr)

      call DMDACreate3d(lbm%comm, DMDA_XYZPERIODIC, DMDA_STENCIL_BOX, &
           options%NX, options%NY, options%NZ, &
           PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, lbm%dm_index_to_ndof(NPHASEDOF), 1, &
           PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, lbm%da_s, ierr)

      call DMDACreate3d(lbm%comm, DMDA_XYZPERIODIC, DMDA_STENCIL_BOX, &
           options%NX, options%NY, options%NZ, &
           PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, lbm%dm_index_to_ndof(NPHASEXBDOF), 1,&
           PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, lbm%da_sb, ierr)

      call DMDACreate3d(lbm%comm, DMDA_XYZPERIODIC, DMDA_STENCIL_BOX, &
           options%NX, options%NY, options%NZ, &
           PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, lbm%dm_index_to_ndof(NFLOWDOF), 1, &
           PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, lbm%da_flow, ierr)

      call DMDAGetCorners(lbm%da_one, xs, ys, zs, lbm%info%xl, lbm%info%yl, lbm%info%zl, ierr)
      call DMDAGetGhostCorners(lbm%da_one, gxs, gys, gzs, lbm%info%gxl, lbm%info%gyl, lbm%info%gzl, ierr)

      ! set lbm%info including corners
      lbm%info%xs = xs+1
      lbm%info%ys = ys+1
      lbm%info%zs = zs+1
      lbm%info%gxs = gxs+1
      lbm%info%gys = gys+1
      lbm%info%gzs = gzs+1

      lbm%info%xe = lbm%info%xs+lbm%info%xl-1
      lbm%info%ye = lbm%info%ys+lbm%info%yl-1
      lbm%info%ze = lbm%info%zs+lbm%info%zl-1
      lbm%info%gxe = lbm%info%gxs+lbm%info%gxl-1
      lbm%info%gye = lbm%info%gys+lbm%info%gyl-1
      lbm%info%gze = lbm%info%gzs+lbm%info%gzl-1

      ! allocate, associate workspace
      allocate(lbm%ux(1:lbm%info%s,lbm%info%gxs:lbm%info%gxe, &
           lbm%info%gys:lbm%info%gye,lbm%info%gzs:lbm%info%gze))
      allocate(lbm%uy(1:lbm%info%s,lbm%info%gxs:lbm%info%gxe, &
           lbm%info%gys:lbm%info%gye,lbm%info%gzs:lbm%info%gze))
      allocate(lbm%uz(1:lbm%info%s,lbm%info%gxs:lbm%info%gxe, &
           lbm%info%gys:lbm%info%gye,lbm%info%gzs:lbm%info%gze))
      lbm%ux = 0
      lbm%uy = 0
      lbm%uz = 0

      allocate(lbm%uxe(1:lbm%info%s,lbm%info%gxs:lbm%info%gxe, &
           lbm%info%gys:lbm%info%gye,lbm%info%gzs:lbm%info%gze))
      allocate(lbm%uye(1:lbm%info%s,lbm%info%gxs:lbm%info%gxe, &
           lbm%info%gys:lbm%info%gye,lbm%info%gzs:lbm%info%gze))
      allocate(lbm%uze(1:lbm%info%s,lbm%info%gxs:lbm%info%gxe, &
           lbm%info%gys:lbm%info%gye,lbm%info%gzs:lbm%info%gze))
      lbm%uxe = 0
      lbm%uye = 0
      lbm%uze = 0

      allocate(lbm%Fx(1:lbm%info%s,lbm%info%gxs:lbm%info%gxe, &
           lbm%info%gys:lbm%info%gye,lbm%info%gzs:lbm%info%gze))
      allocate(lbm%Fy(1:lbm%info%s,lbm%info%gxs:lbm%info%gxe, &
           lbm%info%gys:lbm%info%gye,lbm%info%gzs:lbm%info%gze))
      allocate(lbm%Fz(1:lbm%info%s,lbm%info%gxs:lbm%info%gxe, &
           lbm%info%gys:lbm%info%gye,lbm%info%gzs:lbm%info%gze))
      lbm%Fx = 0
      lbm%Fy = 0
      lbm%Fz = 0

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

      call PetscObjectSetName(lbm%prs_g, 'prs', ierr)
      call PetscObjectSetName(lbm%walls_g, 'walls', ierr)
      call PetscObjectSetName(lbm%rhot_g, 'rhot', ierr)
      call PetscObjectSetName(lbm%rho_g, 'rho', ierr)
      call PetscObjectSetName(lbm%fi_g, 'fi', ierr)
      call PetscObjectSetName(lbm%ut_g, 'ut', ierr)

      ! set up BC vectors
      call BCSetSizes(lbm%bc, lbm%comm, lbm%info)

      CHKERRQ(ierr)
      CHKMEMQ
      return
    end subroutine LBMSetFromOptions

    ! --- destroy things
    subroutine LBMDestroy(lbm, ierr)
      implicit none
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

      if (associated(lbm%ux)) deallocate(lbm%ux)
      if (associated(lbm%uy)) deallocate(lbm%uy)
      if (associated(lbm%uz)) deallocate(lbm%uz)
      if (associated(lbm%uxe)) deallocate(lbm%uxe)
      if (associated(lbm%uye)) deallocate(lbm%uye)
      if (associated(lbm%uze)) deallocate(lbm%uze)
      if (associated(lbm%Fx)) deallocate(lbm%Fx)
      if (associated(lbm%Fy)) deallocate(lbm%Fy)
      if (associated(lbm%Fz)) deallocate(lbm%Fz)

      if (lbm%da_one /= 0) call DMDestroy(lbm%da_one, ierr)
      if (lbm%da_s /= 0) call DMDestroy(lbm%da_s, ierr)
      if (lbm%da_sb /= 0) call DMDestroy(lbm%da_sb, ierr)
      if (lbm%da_flow /= 0) call DMDestroy(lbm%da_flow, ierr)

      call BCDestroy(lbm%bc, ierr)
      call InfoDestroy(lbm%info, ierr)
      return
    end subroutine LBMDestroy

    ! --- do things
    subroutine LBMSetDomain(lbm, corners)
      implicit none
      type(lbm_type) lbm
      PetscScalar,dimension(3,2):: corners
      PetscErrorCode ierr
      PetscScalar deltacoord
      PetscScalar newcoord_x
      PetscScalar newcoord_y
      PetscScalar newcoord_z

      lbm%corners = corners
      ! note, because DAs are all periodic, it really screws up the corners in a 
      ! uniform coordinates DA.  We must adjust these in the nonperiodic cases.
      if (lbm%bc%flags(1) /= 0) then
         deltacoord = (corners(1,2) - corners(1,1))/(lbm%info%NX-1)
         newcoord_x = corners(1,2) + deltacoord
      else
         newcoord_x = corners(1,2)
      endif

      if (lbm%bc%flags(3) /= 0) then
         deltacoord = (corners(2,2) - corners(2,1))/(lbm%info%NY-1)
         newcoord_y = corners(2,2) + deltacoord
      else
         newcoord_y = corners(2,2)
      endif

      if (lbm%bc%flags(5) /= 0) then
         deltacoord = (corners(3,2) - corners(3,1))/(lbm%info%NZ-1)
         newcoord_z = corners(3,2) + deltacoord
      else
         newcoord_z = corners(3,2)
      endif

      call DMDASetUniformCoordinates(lbm%da_one, corners(1,1), newcoord_x, corners(2,1), &
           newcoord_y, corners(3,1), newcoord_z, ierr)
      call DMDASetUniformCoordinates(lbm%da_s, corners(1,1), newcoord_x, corners(2,1), &
           newcoord_y, corners(3,1), newcoord_z, ierr)
      call DMDASetUniformCoordinates(lbm%da_sb, corners(1,1), newcoord_x, corners(2,1), &
           newcoord_y, corners(3,1), newcoord_z, ierr)
      call DMDASetUniformCoordinates(lbm%da_flow, corners(1,1), newcoord_x, corners(2,1), &
           newcoord_y, corners(3,1), newcoord_z, ierr)
      CHKERRQ(ierr)
      CHKMEMQ
      return
    end subroutine LBMSetDomain

    subroutine LBMRun(lbm, istep, kstep, kwrite)
      implicit none

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

      ! output at zero time
      call LBMLocalToGlobal(lbm)
      call LBMOutput(lbm, istep-1, kwrite)

      ! get arrays
      call DMDAVecGetArrayF90(lbm%da_one, lbm%rhot, lbm%rhot_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_one, lbm%prs, lbm%prs_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_flow, lbm%ut, lbm%ut_a, ierr)

      call BCGetArrays(lbm%bc, ierr)


      timername = 'Simulation'
      timer1 => TimingCreate(lbm%comm, timername)
      do lcv_step = istep,kstep
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
               case (1)         ! pseudo-periodic
                  call BCPseudoperiodic(lbm%fi_a, lbm%walls_a, lbm%bc%flags, lbm%info%dim, &
                       lbm%bc%xm_a, lbm%bc%xp_a, lbm%bc%ym_a, lbm%bc%yp_a, lbm%bc%zm_a, &
                       lbm%bc%zp_a, lbm%info)
               case (2)         ! flux
                  call BCFlux(lbm%fi_a, lbm%walls_a, lbm%bc%flags, lbm%info%dim, lbm%bc%xm_a, &
                       lbm%bc%xp_a, lbm%bc%ym_a, lbm%bc%yp_a, lbm%bc%zm_a, lbm%bc%zp_a, lbm%info)
               case (3)         ! pressure
                  call BCPressure(lbm%fi_a, lbm%walls_a, lbm%bc%flags, lbm%info%dim, lbm%bc%xm_a, &
                       lbm%bc%xp_a, lbm%bc%ym_a, lbm%bc%yp_a, lbm%bc%zm_a, lbm%bc%zp_a, lbm%info)
               end select
               bcs_done(lbm%bc%flags(lcv_sides)) = .TRUE. ! only do each bc type once
            endif
         enddo

         call LBMUpdateMoments(lbm%fi_a,lbm%rho_a, lbm%ux,lbm%uy,lbm%uz,lbm%walls_a, lbm%info)

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
         lbm%Fx=0.
         lbm%Fy=0.
         lbm%Fz=0.

         if (lbm%info%s.eq.2) then
            call LBMAddFluidFluidForces(lbm%rho_a, lbm%Fx, lbm%Fy, lbm%Fz, lbm%walls_a, lbm%info)
         endif
         call LBMAddBodyForces(lbm%rho_a, lbm%Fx, lbm%Fy, lbm%Fz, lbm%walls_a, lbm%info)
         call LBMAddFluidSolidForces(lbm%rho_a, lbm%Fx, lbm%Fy, lbm%Fz, lbm%walls_a, lbm%info)
         call LBMZeroBoundaryForces(lbm%bc%flags, lbm%Fx, lbm%Fy, lbm%Fz, lbm%info%dim, lbm%info)

!         call TimingEnd(timer2)
!         call TimingEnd(timer3)
!         timername = 'equilibrium'
!         timer2 => TimingCreate(lbm%comm, timername)

         ! calculate u_equilibrium
         call LBMUpdateUEquilibrium(lbm%fi_a, lbm%rho_a, lbm%ux, lbm%uy, lbm%uz, lbm%walls_a, &
              lbm%uxe, lbm%uye, lbm%uze, lbm%rhot_a, lbm%Fx, lbm%Fy, lbm%Fz, lbm%info)
!         call TimingEnd(timer2)
!         timername = 'collision'
!         timer2 => TimingCreate(lbm%comm, timername)

         ! collision
         call LBMCollision(lbm%fi_a, lbm%rho_a, lbm%uxe, lbm%uye, lbm%uze, lbm%walls_a, lbm%info)
!         call TimingEnd(timer2)

         ! communicate, update fi
         call DMDAVecRestoreArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
         call DMDALocalToLocalBegin(lbm%da_sb, lbm%fi, INSERT_VALUES, lbm%fi, ierr)
         call DMDALocalToLocalEnd(lbm%da_sb, lbm%fi, INSERT_VALUES, lbm%fi, ierr)

         ! check for output?
         if(mod(lcv_step,kwrite).eq.0) then
            ! --  --  update diagnostics
            call LBMUpdateDiagnostics(lbm%rho_a, lbm%ux, lbm%uy, lbm%uz, lbm%walls_a, lbm%ut_a, &
                 lbm%rhot_a, lbm%prs_a, lbm%Fx, lbm%Fy, lbm%Fz, lbm%info)

            ! --  --  restore arrays
            call DMDAVecRestoreArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%da_one, lbm%prs, lbm%prs_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%da_one, lbm%rhot, lbm%rhot_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%da_flow, lbm%ut, lbm%ut_a, ierr)

            ! --  --  output
            call LBMLocalToGlobal(lbm)
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

      ! restore arrays in prep for communication
      call BCRestoreArrays(lbm%bc, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_one, lbm%prs, lbm%prs_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_one, lbm%rhot, lbm%rhot_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_flow, lbm%ut, lbm%ut_a, ierr)

      ! communicate local to global
      call LBMLocalToGlobal(lbm)
      return
    end subroutine LBMRun

    subroutine LBMLocalToGlobal(lbm)
      implicit none
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
      implicit none
      type(lbm_type) lbm
      !      interface
      !         subroutine init_subroutine(walls, info)
      !           use Info_module
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
      implicit none
      type(lbm_type) lbm
      character(len=MAXSTRINGLENGTH) filename
      
      PetscViewer viewer
      PetscErrorCode ierr
      call PetscViewerBinaryOpen(lbm%comm, filename, FILE_MODE_READ, viewer, ierr)
      call VecLoad(lbm%walls_g, viewer, ierr)
      call PetscViewerDestroy(viewer, ierr)
      call DMGlobalToLocalBegin(lbm%da_one, lbm%walls_g, lbm%walls, ierr)
      call DMGlobalToLocalEnd(lbm%da_one, lbm%walls_g, lbm%walls, ierr)
      return
    end subroutine LBMInitializeWallsPetsc

    subroutine LBMInitializeState(lbm, init_subroutine)
      implicit none
      type(lbm_type) lbm
      !      interface
      !         subroutine init_subroutine(fi, rho, ux, uy, uz, walls, info)
      !           use Info_module
      !           type(info_type) info
      !           PetscScalar fi(:)
      !           PetscScalar rho(:)
      !           PetscScalar ux(:,:,:,:)
      !           PetscScalar uy(:,:,:,:)
      !           PetscScalar uz(:,:,:,:)
      !           PetscScalar walls(:)
      !         end subroutine init_subroutine
      !      end interface
      external :: init_subroutine
      PetscErrorCode ierr
      call DMDAVecGetArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
      call DMDAVecGetArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
      call init_subroutine(lbm%fi_a, lbm%rho_a, lbm%ux, lbm%uy, lbm%uz, lbm%walls_a, lbm%info)
      call DMDAVecRestoreArrayF90(lbm%da_one, lbm%walls, lbm%walls_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_sb, lbm%fi, lbm%fi_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%da_s, lbm%rho, lbm%rho_a, ierr)
      return
    end subroutine LBMInitializeState

    function LBMGetDMByIndex( lbm, dm_index ) result(dm)
      implicit none
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
      
      corners = lbm%corners
    end subroutine LBMGetCorners

    subroutine LBMPrintAFew(rho, fi, walls, uxe, uye, uze, info)
      implicit none

      type(info_type) info
      PetscScalar,dimension(1:info%s,info%gxs:info%gxe,info%gys:info%gye,info%gzs:info%gze)::rho
      PetscScalar,dimension(1:info%s,0:info%b,info%gxs:info%gxe,info%gys:info%gye, info%gzs:info%gze)::fi
      PetscScalar,dimension(info%gxs:info%gxe,info%gys:info%gye,info%gzs:info%gze)::walls
      PetscScalar,dimension(1:info%s,info%gxs:info%gxe,info%gys:info%gye,info%gzs:info%gze)::uxe,uye,uze

      ! write(*,*) 'walls check:'
      ! write(*,*) 'walls(1,0,13):', walls(2,1,14)
      ! write(*,*) 'walls(1,1,13):', walls(2,2,14)
      ! write(*,*) 'walls(1,NY,13):', walls(2,16,14)
      ! write(*,*) 'walls(1,NY-1,13):', walls(2,15,14)

      ! write(*,*) 'walls(1,0,13):', walls(1,1,14)
      ! write(*,*) 'walls(1,1,13):', walls(1,2,14)
      ! write(*,*) 'walls(1,NY,13):', walls(1,16,14)
      ! write(*,*) 'walls(1,NY-1,13):', walls(1,15,14)

      ! write(*,*) 'walls(1,0,13):', walls(3,1,14)
      ! write(*,*) 'walls(1,1,13):', walls(3,2,14)
      ! write(*,*) 'walls(1,NY,13):', walls(3,16,14)
      ! write(*,*) 'walls(1,NY-1,13):', walls(3,15,14)

      if (info%id.eq.0) then
         write(*,*) '---------------------------------------'
         write(*,*) 'output: rho(0,0,8,10):', rho(1,1,2,11)
         write(*,*) 'output: uxe(0,1,8,10):', uxe(1,1,2,11)
         write(*,*) 'output: uye(0,1,8,10):', uye(1,1,2,11)
         write(*,*) 'output: uze(0,1,8,10):', uze(1,1,2,11)
         write(*,*) 'output: fi(0,0,0,8,10):', fi(1,0,1,2,11)
         write(*,*) 'output: fi(0,1,0,8,10):', fi(1,2,1,2,11)
         write(*,*) 'output: fi(0,5,0,8,10):', fi(1,4,1,2,11)
      endif
      if (info%id.eq.info%nproc-1) then
         write(*,*) '---------------------------------------'
         print *, 'walls:', walls(1,25,99)
         print *, 'walls:', walls(2,25,99)
         print *, 'walls:', walls(3,25,99)
         write(*,*) 'output: walls(0,8,72):', walls(1,26,99)
         write(*,*) 'output: walls(0,8,72):', walls(1,25,99)
         write(*,*) 'output: rho(0,0,8,72):', rho(1,1,25,99)
         write(*,*) 'output: rho(0,0,8,72):', rho(2,1,25,99)
         write(*,*) 'output: rho(0,0,8,72):', rho(1,1,23,99)
         write(*,*) 'output: rho(0,0,8,72):', rho(1,1,25,99)
         write(*,*) 'output: uxe(0,1,8,72):', uxe(1,1,25,99)
         write(*,*) 'output: uye(0,1,8,72):', uye(1,1,25,99)
         write(*,*) 'output: uze(0,1,8,72):', uze(1,1,25,99)
         write(*,*) 'output: fi(0,0,0,8,72):', fi(1,0,1,25,99)
         write(*,*) 'output: fi(0,1,0,8,72):', fi(1,2,1,25,99)
         write(*,*) 'output: fi(0,5,0,8,72):', fi(1,4,1,25,99)
         write(*,*) '---------------------------------------'
         write(*,*) 'output: rho(0,0,8,90):', rho(1,1,2,91)
         write(*,*) 'output: rho(1,0,8,90):', rho(2,1,2,91)
         write(*,*) 'output: uxe(0,1,8,90):', uxe(1,1,2,91)
         write(*,*) 'output: uye(0,1,8,90):', uye(1,1,2,91)
         write(*,*) 'output: uze(0,1,8,90):', uze(1,1,2,91)
         write(*,*) 'output: fi(0,0,0,8,90):', fi(1,0,1,2,91)
         write(*,*) 'output: fi(0,1,0,8,90):', fi(1,2,1,2,91)
         write(*,*) 'output: fi(0,5,0,8,90):', fi(1,4,1,2,91)
         write(*,*) '========================================'
      endif
      return
    end subroutine LBMPrintAFew
  end module LBM_module

