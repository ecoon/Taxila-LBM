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
    use LBM_Discretization_module
    use LBM_Grid_module
    use LBM_BC_module
    use LBM_Flow_module
!    use LBM_Transport_module
    use LBM_Options_module
    implicit none

    private

#include "lbm_definitions.h"
    type, public:: lbm_type
       MPI_Comm comm
       type(options_type),pointer:: options
       type(grid_type),pointer:: grid
       type(bc_type),pointer:: bc
       type(flow_type),pointer:: flow
!       type(transport_type),pointer:: transport

       Vec walls
       Vec walls_g
       PetscScalar,pointer:: walls_a(:)

       character(len=MAXSTRINGLENGTH) name       
    end type lbm_type
    
    public :: LBMCreate, &
         LBMDestroy, &
         LBMSetName, &
         LBMSetFromOptions, &
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
    function LBMCreate(comm) result(lbm)
      type(lbm_type),pointer:: lbm
      MPI_Comm comm

      allocate(lbm)
      lbm%comm = comm

      lbm%options => OptionsCreate(comm)
      lbm%grid => GridCreate(comm)
      lbm%bc => BCCreate(comm)
      nullify(lbm%flow)
!      nullify(lbm%transport)

      lbm%walls = 0
      lbm%walls_g = 0
      nullify(lbm%walls_a)

      lbm%name = ''
    end function LBMCreate

    ! --- destroy things
    subroutine LBMDestroy(lbm, ierr)
      type(lbm_type) lbm
      PetscErrorCode ierr

      if (lbm%walls /= 0) call VecDestroy(lbm%walls,ierr)
      if (lbm%walls_g /= 0) call VecDestroy(lbm%walls_g,ierr)

      call BCDestroy(lbm%bc, ierr)
      if (associated(lbm%flow)) call FlowDestroy(lbm%flow,ierr)
!      if (associated(lbm%transport)) call TransportDestroy(lbm%transport,ierr)
      return
    end subroutine LBMDestroy

    subroutine LBMSetName(lbm, name) 
      type(lbm_type) lbm 
      character(len=MAXWORDLENGTH):: name       
      lbm%name = name
    end subroutine LBMSetName

    subroutine LBMSetFromOptions(lbm, options, ierr)
      type(lbm_type) lbm
      type(options_type) options
      PetscErrorCode ierr
      
      call GridSetFromOptions(lbm%grid, options, ierr)
      call GridSetName(lbm%grid, lbm%name)
      
      lbm%flow => FlowCreate(lbm%comm)
      call FlowSetFromOptions(lbm%flow, options)
!      if (options%transport_disc /= NULL_DISCRETIZATION) then
!         lbm%transport => TransportCreate(lbm%comm)
!         call TransportSetFromOptions(lbm%transport, options)
!      end if

      call BCSetFromOptions(lbm%bc,options)
    end subroutine LBMSetFromOptions

    subroutine LBMSetUp(lbm)
      PetscScalar zero

      call GridSetUp(lbm%grid)
      call FlowSetGrid(lbm%flow, lbm%grid)
      call FlowSetUp(lbm%flow)
      ! if (associated(lbm%transport)) then
      !    call TransportSetGrid(lbm%transport, lbm%grid)
      !    call TransportSetUp(lbm%transport)
      ! end if
      call BCSetGrid(lbm%bc, lbm%grid)
      call BCSetUp(lbm%bc)

      ! get vectors
      call DMCreateLocalVector(lbm%grid%da(ONEDOF), lbm%walls, ierr)
      call DMCreateGlobalVector(lbm%grid%da(ONEDOF), lbm%walls_g, ierr)

      zero = 0.d0
      call VecSet(lbm%walls_g, zero, ierr)
      call VecSet(lbm%walls, zero, ierr)

      call PetscObjectSetName(lbm%walls_g, trim(lbm%name)//'walls', ierr)

      CHKERRQ(ierr)
      return
    end subroutine LBMSetUp

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
      call DMDALocalToLocalBegin(lbm%grid%da(ONEDOF), lbm%walls, INSERT_VALUES, &
           lbm%walls, ierr)
      call DMDALocalToLocalEnd(lbm%grid%da(ONEDOF), lbm%walls, INSERT_VALUES, &
           lbm%walls, ierr)

      call DMDALocalToLocalBegin(lbm%grid%da(NPHASEXBDOF), lbm%flow%distribution%fi, &
           INSERT_VALUES, lbm%flow%distribution%fi, ierr)
      call DMDALocalToLocalEnd(lbm%grid%da(NPHASEXBDOF), lbm%flow%distribution%fi, &
           INSERT_VALUES, lbm%flow%distribution%fi, ierr)

      if (associated(lbm%transport)) then
         call DMDALocalToLocalBegin(lbm%grid%da(NCOMPONENTXBDOF), &
              lbm%transport%distribution%fi, INSERT_VALUES, &
              lbm%transport%distribution%fi, ierr)
         call DMDALocalToLocalEnd(lbm%grid%da(NCOMPONENTXBDOF), &
              lbm%transport%distribution%fi, INSERT_VALUES, &
              lbm%transport%distribution%fi, ierr)
      end if

      if (istep.eq.0) then
         ! get arrays
         call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
         call FlowGetArrays(lbm%flow)

         ! update values for zero time i/o
         call FlowUpdateMoments(lbm%fi_a, lbm%rho_a, lbm%vel, lbm%walls_a, lbm%info)

!!!!!!!!!!!!! here!!!!!!!!

         call LBMUpdateUEquilibrium(lbm%fi_a, lbm%rho_a, lbm%vel, lbm%walls_a, &
              lbm%vel_eq, lbm%rhot_a, lbm%forces, lbm%info, lbm%constants)
         call LBMUpdateDiagnostics(lbm%rho_a, lbm%vel, lbm%walls_a, lbm%velt_a, &
              lbm%rhot_a, lbm%prs_a, lbm%forces, lbm%info, lbm%constants)
         
         ! --  --  restore arrays
         call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
         call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%prs, lbm%prs_a, ierr)
         call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%rhot, lbm%rhot_a, ierr)
         call DMDAVecRestoreArrayF90(lbm%grid%da(NPHASEDOF), lbm%rho, lbm%rho_a, ierr)
         call DMDAVecRestoreArrayF90(lbm%grid%da(NFLOWDOF), lbm%velt, lbm%velt_a, ierr)
         call DMDAVecRestoreArrayF90(lbm%grid%da(NPHASEXBDOF), lbm%fi, lbm%fi_a, ierr)

         call LBMLocalToGlobal(lbm, ierr)
         ! output at zero time
         call LBMOutput(lbm, istep, kwrite)
      else
         call LBMLocalToGlobal(lbm, ierr)
      endif

      ! get arrays
      call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%rhot, lbm%rhot_a, ierr)
      call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%prs, lbm%prs_a, ierr)
      call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
      call DMDAVecGetArrayF90(lbm%grid%da(NPHASEXBDOF), lbm%fi, lbm%fi_a, ierr)
      call DMDAVecGetArrayF90(lbm%grid%da(NPHASEDOF), lbm%rho, lbm%rho_a, ierr)
      call DMDAVecGetArrayF90(lbm%grid%da(NFLOWDOF), lbm%velt, lbm%velt_a, ierr)

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
         call DMDAVecRestoreArrayF90(lbm%grid%da(NPHASEDOF),lbm%rho, lbm%rho_a, ierr)
         call DMDALocalToLocalBegin(lbm%grid%da(NPHASEDOF), lbm%rho, INSERT_VALUES, lbm%rho, ierr)
         call DMDALocalToLocalEnd(lbm%grid%da(NPHASEDOF), lbm%rho, INSERT_VALUES, lbm%rho, ierr)
         call DMDAVecGetArrayF90(lbm%grid%da(NPHASEDOF), lbm%rho, lbm%rho_a, ierr)
!         call TimingEnd(timer2)
!         timer2%name = 'communication + forces'
!         timername = 'forces'
!         timer3 => TimingCreate(lbm%comm, timername)


         !calculate forces
         lbm%forces=0.d0

         if (lbm%info%nphases > 1) then
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
         call DMDAVecRestoreArrayF90(lbm%grid%da(NPHASEXBDOF), lbm%fi, lbm%fi_a, ierr)
         call DMDALocalToLocalBegin(lbm%grid%da(NPHASEXBDOF), lbm%fi, INSERT_VALUES, lbm%fi, ierr)
         call DMDALocalToLocalEnd(lbm%grid%da(NPHASEXBDOF), lbm%fi, INSERT_VALUES, lbm%fi, ierr)

         ! check for oveltput?
         if(mod(lcv_step,kwrite).eq.0) then
            ! --  --  update diagnostics
            call LBMUpdateDiagnostics(lbm%rho_a, lbm%vel, lbm%walls_a, lbm%velt_a, &
                 lbm%rhot_a, lbm%prs_a, lbm%forces, lbm%info, lbm%constants)

            ! --  --  restore arrays
            call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%prs, lbm%prs_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%rhot, lbm%rhot_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%grid%da(NPHASEDOF), lbm%rho, lbm%rho_a, ierr)
            call DMDAVecRestoreArrayF90(lbm%grid%da(NFLOWDOF), lbm%velt, lbm%velt_a, ierr)

            ! --  --  output
            call LBMLocalToGlobal(lbm, ierr)
            call LBMOutput(lbm, lcv_step, kwrite)

            ! --  --  reopen arrays
            call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
            call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%prs, lbm%prs_a, ierr)
            call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%rhot, lbm%rhot_a, ierr)
            call DMDAVecGetArrayF90(lbm%grid%da(NPHASEDOF), lbm%rho, lbm%rho_a, ierr)
            call DMDAVecGetArrayF90(lbm%grid%da(NFLOWDOF), lbm%velt, lbm%velt_a, ierr)
         endif

         call DMDAVecGetArrayF90(lbm%grid%da(NPHASEXBDOF), lbm%fi, lbm%fi_a, ierr)
      end do

      timerunits = 'timestep'
      call TimingEndPerUnit(timer1, (kstep-istep+1), timerunits)
      call TimingDestroy(timer1)

      ! restore arrays in prep for communication
      call BCRestoreArrays(lbm%bc, ierr)
      call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%prs, lbm%prs_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%rhot, lbm%rhot_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%grid%da(NPHASEXBDOF), lbm%fi, lbm%fi_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%grid%da(NPHASEDOF), lbm%rho, lbm%rho_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%grid%da(NFLOWDOF), lbm%velt, lbm%velt_a, ierr)

      ! communicate local to global
      call LBMLocalToGlobal(lbm, ierr)
      return
    end subroutine LBMRun

    subroutine LBMLocalToGlobal(lbm, ierr)
      type(lbm_type) lbm
      PetscErrorCode ierr

      call DMLocalToGlobalBegin(lbm%grid%da(NFLOWDOF), lbm%velt, INSERT_VALUES, &
           lbm%velt_g, ierr)
      call DMLocalToGlobalEnd(lbm%grid%da(NFLOWDOF), lbm%velt, INSERT_VALUES, &
           lbm%velt_g, ierr)

      call DMLocalToGlobalBegin(lbm%grid%da(ONEDOF), lbm%prs, INSERT_VALUES, &
           lbm%prs_g, ierr)
      call DMLocalToGlobalEnd(lbm%grid%da(ONEDOF), lbm%prs, INSERT_VALUES, &
           lbm%prs_g, ierr)

      call DMLocalToGlobalBegin(lbm%grid%da(ONEDOF), lbm%walls, INSERT_VALUES, &
           lbm%walls_g, ierr)
      call DMLocalToGlobalEnd(lbm%grid%da(ONEDOF), lbm%walls, INSERT_VALUES, &
           lbm%walls_g, ierr)

      call DMLocalToGlobalBegin(lbm%grid%da(ONEDOF), lbm%rhot, INSERT_VALUES, &
           lbm%rhot_g, ierr)
      call DMLocalToGlobalEnd(lbm%grid%da(ONEDOF), lbm%rhot, INSERT_VALUES, &
           lbm%rhot_g, ierr)

      call DMLocalToGlobalBegin(lbm%grid%da(NPHASEDOF), lbm%rho, INSERT_VALUES, &
           lbm%rho_g, ierr)
      call DMLocalToGlobalEnd(lbm%grid%da(NPHASEDOF), lbm%rho, INSERT_VALUES, &
           lbm%rho_g, ierr)

      call DMLocalToGlobalBegin(lbm%grid%da(NPHASEXBDOF), lbm%fi, INSERT_VALUES, &
           lbm%fi_g, ierr)
      call DMLocalToGlobalEnd(lbm%grid%da(NPHASEXBDOF), lbm%fi, INSERT_VALUES, &
           lbm%fi_g, ierr)
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

      call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
      call init_subroutine(lbm%walls_a, lbm%options%walls_file, lbm%info)
      call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
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
      call DMGlobalToLocalBegin(lbm%grid%da(ONEDOF), lbm%walls_g, INSERT_VALUES, lbm%walls, ierr)
      call DMGlobalToLocalEnd(lbm%grid%da(ONEDOF), lbm%walls_g, INSERT_VALUES, lbm%walls, ierr)
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
      call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
      call DMDAVecGetArrayF90(lbm%grid%da(NPHASEXBDOF), lbm%fi, lbm%fi_a, ierr)
      call DMDAVecGetArrayF90(lbm%grid%da(NPHASEDOF), lbm%rho, lbm%rho_a, ierr)
      call init_subroutine(lbm%fi_a, lbm%rho_a, lbm%vel, lbm%walls_a, lbm%info, lbm%constants)
      call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%grid%da(NPHASEXBDOF), lbm%fi, lbm%fi_a, ierr)
      call DMDAVecRestoreArrayF90(lbm%grid%da(NPHASEDOF), lbm%rho, lbm%rho_a, ierr)
      return
    end subroutine LBMInitializeState

    ! function LBMGetDMByIndex( lbm, dm_index ) result(dm)
    !   type(lbm_type) lbm
    !   PetscInt dm_index
    !   DM,pointer:: dm
    !   PetscErrorCode ierr
    !   nullify(dm)
    !   select case(dm_index)
    !   case (ONEDOF)
    !      dm => lbm%grid%da(ONEDOF)
    !   case (NPHASEDOF)
    !      dm => lbm%grid%da(NPHASEDOF)
    !   case (NPHASEXBDOF)
    !      dm => lbm%grid%da(NPHASEXBDOF)
    !   case (NFLOWDOF)
    !      dm => lbm%grid%da(NFLOWDOF)
    !   end select
    ! end function LBMGetDMByIndex

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
    character(len=MAXSTRINGLENGTH):: flnm0,flnm1,flnm2,flnm3,flnm6,flnm7,flnm8
    integer charlen
    Vec coords

    charlen = LEN_TRIM(lbm%options%output_prefix)

    if (lbm%info%id.eq.0) then
       write(*,*) 'outputing (istep, kwrite)', istep, kwrite
    endif

    charlen = LEN_TRIM(lbm%options%output_prefix)

    write(stringformat, '("(I0.",I1,")")'), MAXIODIGITS
    write(outnum, stringformat) istep/kwrite

    if (istep == 0) then
       flnm0=lbm%options%output_prefix(1:charlen)//'coords'//outnum//'.dat'
       call DMDAGetCoordinates(lbm%grid%da(ONEDOF), coords, ierr)
       call PetscViewerCreate(lbm%comm, viewer, ierr)
       call PetscViewerSetType(viewer, PETSCVIEWERBINARY, ierr)
       call PetscViewerBinarySetMPIIO(viewer, ierr)
       call PetscViewerFileSetMode(viewer, FILE_MODE_WRITE, ierr)
       call PetscViewerFileSetName(viewer, flnm0, ierr)
       call VecView(coords, viewer, ierr)
       call PetscViewerDestroy(viewer,ierr)
    end if

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
    call VecView(lbm%velt_g, viewer, ierr)
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
    call VecLoad(lbm%velt_g, viewer, ierr)
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

