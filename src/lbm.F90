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
    use LBM_Distribution_Function_module
    use LBM_Walls_module
    use LBM_Flow_module
    use LBM_Transport_module
    use LBM_IO_module
    use LBM_Options_module
    use LBM_BC_module
    implicit none

    private

#include "lbm_definitions.h"
    type, public:: lbm_type
       MPI_Comm comm
       type(options_type),pointer:: options
       type(grid_type),pointer:: grid
       type(walls_type), pointer :: walls
       type(flow_type),pointer:: flow
       type(transport_type),pointer:: transport
       type(io_type),pointer :: io

       character(len=MAXWORDLENGTH) name       
    end type lbm_type

    interface LBMInitializeState
       module procedure LBMInitializeState_Flow
       module procedure LBMInitializeState_FlowTransport
    end interface
    
    interface LBMInit
       module procedure LBMInit1
       module procedure LBMInit2
    end interface

    interface LBMRun
       module procedure LBMRun1
       module procedure LBMRun2
    end interface

    public :: LBMCreate, &
         LBMDestroy, &
         LBMSetName, &
         LBMSetFromOptions, &
         LBMSetUp, &
         LBMInit, &
         LBMRun, &
         LBMOutput, &
         LBMInitializeState, &
         LBMInitializeStateRestarted, &
         LBMLoadSteadyStateFlow, &
         LBMGetCorners

  contains
    function LBMCreate(comm) result(lbm)
      type(lbm_type),pointer:: lbm
      MPI_Comm comm

      allocate(lbm)
      lbm%comm = comm

      lbm%options => OptionsCreate(comm)
      lbm%grid => GridCreate(comm)
      lbm%io => IOCreate(comm)
      lbm%walls => WallsCreate(lbm%comm)
      lbm%flow => FlowCreate(lbm%comm)
      nullify(lbm%transport)

      lbm%name = ''
    end function LBMCreate

    ! --- destroy things
    subroutine LBMDestroy(lbm, ierr)
      type(lbm_type) lbm
      PetscErrorCode ierr

      if (associated(lbm%flow)) call FlowDestroy(lbm%flow,ierr)
      if (associated(lbm%transport)) call TransportDestroy(lbm%transport,ierr)
      call GridDestroy(lbm%grid, ierr)
      call IODestroy(lbm%io, ierr)
      if (associated(lbm%options%waypoints)) deallocate(lbm%options%waypoints)
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

      call IOSetFromOptions(lbm%io, options, ierr);CHKERRQ(ierr)
      call GridSetFromOptions(lbm%grid, options, ierr);CHKERRQ(ierr)
      call GridSetName(lbm%grid, lbm%name)
      call GridSetPhysicalScales(lbm%grid, ierr);CHKERRQ(ierr)

      call WallsSetGrid(lbm%walls, lbm%grid)
      call WallsSetFromOptions(lbm%walls, options, ierr);CHKERRQ(ierr)
      
      call FlowSetGrid(lbm%flow, lbm%grid)
      call FlowSetFromOptions(lbm%flow, options, ierr);CHKERRQ(ierr)
      call FlowSetPhysicalScales(lbm%flow,ierr);CHKERRQ(ierr)
      if (options%transport_disc /= NULL_DISCRETIZATION) then
         lbm%transport => TransportCreate(lbm%comm)
         call TransportSetGrid(lbm%transport, lbm%grid)
         call TransportSetFromOptions(lbm%transport, options, ierr);CHKERRQ(ierr)
      end if
    end subroutine LBMSetFromOptions

    subroutine LBMSetUp(lbm)
      type(lbm_type) lbm
      PetscErrorCode ierr
      PetscScalar zero
      zero = 0.d0

      ! set up the DA sizes
      print*, 'onedof:', ONEDOF
      print*, 'nx:', lbm%grid%info%NX
      print*, 'SHAPE:', SHAPE(lbm%grid%da_sizes)
      lbm%grid%da_sizes(ONEDOF) = 1
      lbm%grid%da_sizes(NPHASEDOF) = lbm%flow%nphases
      lbm%grid%da_sizes(NPHASEXBDOF) = lbm%flow%nphases*(lbm%flow%disc%b+1)
      lbm%grid%da_sizes(NFLOWDOF) = lbm%flow%ndims
      if (associated(lbm%transport)) then
         lbm%grid%da_sizes(NSPECIEDOF) = lbm%transport%nspecies
         lbm%grid%da_sizes(NSPECIEXBDOF) = lbm%transport%nspecies*(lbm%transport%disc%b+1)
      end if

      call GridSetUp(lbm%grid)
      call WallsSetUp(lbm%walls)
      call FlowSetUp(lbm%flow)
      if (associated(lbm%transport)) then
         call TransportSetUp(lbm%transport)
      end if

      CHKERRQ(ierr)
      return
    end subroutine LBMSetUp
    
    subroutine LBMInit1(lbm, istep)
      type(lbm_type) lbm
      PetscInt istep
      call LBMInit2(lbm, istep, PETSC_FALSE)
    end subroutine LBMInit1

    subroutine LBMInit2(lbm, istep, supress_output)
      type(lbm_type) lbm
      PetscInt istep
      PetscErrorCode ierr
      PetscBool supress_output

      call WallsCommunicate(lbm%walls)

      if (istep.eq.0) then
         call DMDALocalToLocalBegin(lbm%grid%da(NPHASEXBDOF), lbm%flow%distribution%fi, &
              INSERT_VALUES, lbm%flow%distribution%fi, ierr)
         call DMDALocalToLocalEnd(lbm%grid%da(NPHASEXBDOF), lbm%flow%distribution%fi, &
              INSERT_VALUES, lbm%flow%distribution%fi, ierr)
         if (associated(lbm%transport)) then
            call DMDALocalToLocalBegin(lbm%grid%da(NSPECIEXBDOF), &
                 lbm%transport%distribution%fi, INSERT_VALUES, &
                 lbm%transport%distribution%fi, ierr)
            call DMDALocalToLocalEnd(lbm%grid%da(NSPECIEXBDOF), &
                 lbm%transport%distribution%fi, INSERT_VALUES, &
                 lbm%transport%distribution%fi, ierr)
         end if

         ! get arrays
         call FlowGetArrays(lbm%flow, ierr)
         if (associated(lbm%transport)) then
            call TransportGetArrays(lbm%transport, ierr)
         end if

         ! update values and view
         if ((.not.supress_output).and.(lbm%grid%info%rank.eq.0)) then
            write(*,*) 'outputing step', istep, 'to file', lbm%io%counter
         endif

         if ((.not.lbm%options%steadystate).or. &
              (lbm%options%steadystate_rampup_steps > 0)) then
            call FlowUpdateMoments(lbm%flow, lbm%walls%walls_a)
            call FlowUpdateDiagnostics(lbm%flow, lbm%walls%walls_a)
         end if
         if (.not.supress_output) then
            call FlowOutputDiagnostics(lbm%flow, lbm%walls%walls_a, lbm%io)
         end if

         if (associated(lbm%transport)) then
            call TransportUpdateMoments(lbm%transport, lbm%walls%walls_a)
            call TransportUpdateDiagnostics(lbm%transport, lbm%walls%walls_a)
            if (.not.supress_output) then
               call TransportOutputDiagnostics(lbm%transport, lbm%walls%walls_a, lbm%io)
            end if
         end if

         if (.not.supress_output) then
           call WallsOutputDiagnostics(lbm%walls, lbm%io)
           ! view grid coordinates
           call GridViewCoordinates(lbm%grid, lbm%io)
         end if
      else
         ! just get arrays
         call FlowGetArrays(lbm%flow, ierr)
      endif
      call IOIncrementCounter(lbm%io)
    end subroutine LBMInit2

    subroutine LBMRun1(lbm, istep, kstep, kwrite)
      type(lbm_type) lbm
      PetscInt istep
      PetscInt kstep
      PetscInt kwrite
      call LBMRun2(lbm, istep, kstep, kwrite, PETSC_FALSE)
    end subroutine LBMRun1

    subroutine LBMRun2(lbm, istep, kstep, kwrite, supress_output)
      ! input
      type(lbm_type) lbm
      PetscInt istep
      PetscInt kstep
      PetscInt kwrite
      PetscBool supress_output

      ! local
      PetscErrorCode ierr
      PetscInt lcv_step
      type(timing_type),pointer:: timer1, timer2, timer3
      character(len=MAXWORDLENGTH) timerunits
      character(len=MAXWORDLENGTH) timername

      timername = trim(lbm%name)//'Simulation'
      timer1 => TimingCreate(lbm%comm, timername)
      do lcv_step = istep+1,kstep
         ! streaming
         if ((.not.lbm%options%steadystate).or. &
              (lcv_step < lbm%options%steadystate_rampup_steps)) then
            call FlowStream(lbm%flow)
         end if

         if (associated(lbm%transport)) then
            call TransportStream(lbm%transport)
         end if

         ! internal fluid-solid boundary conditions
         if ((.not.lbm%options%steadystate).or. &
              (lcv_step < lbm%options%steadystate_rampup_steps)) then
            call FlowBounceback(lbm%flow, lbm%walls%walls_a)
         end if

         if (associated(lbm%transport)) then
            call TransportReactWithWalls(lbm%transport, lbm%walls%walls_a)
         end if

         ! external boundary conditions
         if ((.not.lbm%options%steadystate).or. &
              (lcv_step < lbm%options%steadystate_rampup_steps)) then
            call FlowApplyBCs(lbm%flow, lbm%walls%walls_a)
         end if

         if (associated(lbm%transport)) then
            call TransportApplyBCs(lbm%transport, lbm%walls%walls_a)
         end if

         ! update moments for rho, psi
         if ((.not.lbm%options%steadystate).or. &
              (lcv_step < lbm%options%steadystate_rampup_steps)) then
            call FlowUpdateMoments(lbm%flow, lbm%walls%walls_a)
         end if

         if (associated(lbm%transport)) then
            call TransportUpdateMoments(lbm%transport, lbm%walls%walls_a)
         end if

         ! add in momentum forcing terms
         if ((.not.lbm%options%steadystate).or. &
              (lcv_step < lbm%options%steadystate_rampup_steps)) then
            call DistributionCommunicateDensity(lbm%flow%distribution)
            call FlowCalcForces(lbm%flow, lbm%walls%walls_a)
            call BCZeroForces(lbm%flow%bc, lbm%flow%forces, lbm%flow%distribution)
         ! reaction?

         ! collision
            call FlowCollision(lbm%flow, lbm%walls%walls_a)
         end if

         if (associated(lbm%transport)) then
            if ((.not.lbm%options%steadystate).or. &
                 (lcv_step < lbm%options%steadystate_rampup_steps)) then
               call FlowUpdateDiagnostics(lbm%flow, lbm%walls%walls_a)
            end if
            call TransportCollision(lbm%transport, lbm%walls%walls_a, lbm%flow)
         end if

         ! communicate, update fi
         if ((.not.lbm%options%steadystate).or. &
              (lcv_step < lbm%options%steadystate_rampup_steps)) then
            call DistributionCommunicateFi(lbm%flow%distribution)
         end if

         if (associated(lbm%transport)) then
            call DistributionCommunicateFi(lbm%transport%distribution)
         end if

         
         ! check for output
         if (lbm%options%current_waypoint > 0) then
            if (lbm%options%waypoints(lbm%options%current_waypoint) .eq. lcv_step) then
               call LBMOutput(lbm, lcv_step)
               lbm%options%current_waypoint = lbm%options%current_waypoint + 1
            else if((kwrite > 0) .and. (mod(lcv_step,kwrite).eq.0)) then
               call LBMOutput(lbm, lcv_step)
            endif
         else if((kwrite > 0) .and. (mod(lcv_step,kwrite).eq.0)) then
            call LBMOutput(lbm, lcv_step)
         endif
      end do

      timerunits = 'timestep'
      call TimingEndPerUnit(timer1, (kstep-istep+1), timerunits, supress_output)
      call TimingDestroy(timer1)
      return
    end subroutine LBMRun2
    
    subroutine LBMOutput(lbm, istep)
      type(lbm_type) lbm
      PetscInt istep

      if (lbm%grid%info%rank.eq.0) then
         write(*,*) 'outputing step', istep, 'to file', lbm%io%counter
      endif
      if ((.not.lbm%options%steadystate).or. &
           (istep < lbm%options%steadystate_rampup_steps)) then
         call FlowUpdateDiagnostics(lbm%flow, lbm%walls%walls_a)
         call FlowOutputDiagnostics(lbm%flow, lbm%walls%walls_a, lbm%io)
      end if
      
      if (associated(lbm%transport)) then
         call TransportUpdateDiagnostics(lbm%transport, lbm%walls%walls_a)
         call TransportOutputDiagnostics(lbm%transport, lbm%walls%walls_a, lbm%io)
      end if
      call IOIncrementCounter(lbm%io)
    end subroutine LBMOutput

    subroutine LBMInitializeState_Flow(lbm, init_subroutine)
      type(lbm_type) lbm
      external :: init_subroutine
      PetscErrorCode ierr

      call FlowGetArrays(lbm%flow, ierr)
      call init_subroutine(lbm%flow%distribution%fi_a, lbm%flow%distribution%rho_a, &
           lbm%flow%distribution%flux, lbm%walls%walls_a, lbm%flow%distribution, &
           lbm%flow%phases, lbm%options)
      call FlowRestoreArrays(lbm%flow, ierr)
      return
    end subroutine LBMInitializeState_Flow

    subroutine LBMInitializeState_FlowTransport(lbm, flow_subroutine, trans_subroutine)
      type(lbm_type) lbm
      external :: flow_subroutine, trans_subroutine
      PetscErrorCode ierr

      ! flow init
      call FlowGetArrays(lbm%flow, ierr)
      call flow_subroutine(lbm%flow%distribution%fi_a, lbm%flow%distribution%rho_a, &
           lbm%flow%distribution%flux, lbm%walls%walls_a, lbm%flow%distribution, &
           lbm%flow%phases, lbm%options)
      call FlowRestoreArrays(lbm%flow, ierr)

      ! transport init
      call TransportGetArrays(lbm%transport, ierr)
      call trans_subroutine(lbm%transport%distribution%fi_a, &
           lbm%transport%distribution%rho_a, lbm%transport%distribution%flux, &
           lbm%walls%walls_a, lbm%transport%distribution, lbm%transport%species, lbm%options)
      call TransportRestoreArrays(lbm%transport, ierr)
    end subroutine LBMInitializeState_FlowTransport

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
      
      corners = lbm%grid%info%corners
    end subroutine LBMGetCorners

    subroutine LBMInitializeStateRestarted(lbm, istep, kwrite)
      use petsc
      type(lbm_type) lbm
      PetscInt istep, kwrite
      PetscErrorCode ierr

      if (lbm%options%restart_counter > 0 ) then
         lbm%io%counter = lbm%options%restart_counter
      else if (kwrite > 0) then
         lbm%io%counter = istep/kwrite
      else
         SETERRQ(PETSC_COMM_WORLD, 1, &
              'must set -restart_counter to define which file to restart from', ierr)
      end if

      if (lbm%grid%info%rank.eq.0) then
         write(*,*) 'reading step', istep, 'from file', lbm%io%counter
      endif
      call IOLoad(lbm%io, lbm%flow%distribution%fi_g, 'fi')
      
      call DMGlobalToLocalBegin(lbm%grid%da(NPHASEXBDOF), lbm%flow%distribution%fi_g, &
           INSERT_VALUES, lbm%flow%distribution%fi, ierr)
      call DMGlobalToLocalEnd(lbm%grid%da(NPHASEXBDOF), lbm%flow%distribution%fi_g, &
           INSERT_VALUES, lbm%flow%distribution%fi, ierr)
      return
    end subroutine LBMInitializeStateRestarted
    
    subroutine LBMLoadSteadyStateFlow(lbm, filename)
      type(lbm_type) lbm 
      character(len=MAXSTRINGLENGTH) filename
      PetscViewer viewer
      PetscErrorCode ierr
      call PetscViewerBinaryOpen(lbm%comm, trim(filename), FILE_MODE_READ, viewer, ierr)
      CHKERRQ(ierr)
      call VecLoad(lbm%flow%velt_g, viewer, ierr)
      CHKERRQ(ierr)
      call PetscViewerDestroy(viewer, ierr)
    end subroutine LBMLoadSteadyStateFlow
  end module LBM_module

