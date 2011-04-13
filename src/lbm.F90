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
    use LBM_Distribution_Function_module
    use LBM_Flow_module
    use LBM_Transport_module
    use LBM_IO_module
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
       type(transport_type),pointer:: transport
       type(io_type),pointer :: io

       Vec walls
       Vec walls_g
       PetscScalar,pointer:: walls_a(:)

       character(len=MAXWORDLENGTH) name       
    end type lbm_type
    
    public :: LBMCreate, &
         LBMDestroy, &
         LBMSetName, &
         LBMSetFromOptions, &
         LBMSetUp, &
         LBMRun, &
         LBMInitializeWalls, &
         LBMInitializeWallsPetsc, &
         LBMInitializeState, &
         LBMInitializeStateRestarted, &
         LBMGetCorners

  contains
    function LBMCreate(comm) result(lbm)
      type(lbm_type),pointer:: lbm
      MPI_Comm comm

      allocate(lbm)
      lbm%comm = comm

      lbm%options => OptionsCreate(comm)
      lbm%grid => GridCreate(comm)
      lbm%bc => BCCreate(comm)
      lbm%io => IOCreate(comm)
      lbm%flow => FlowCreate(lbm%comm)
      nullify(lbm%transport)

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
      if (associated(lbm%transport)) call TransportDestroy(lbm%transport,ierr)
      call GridDestroy(lbm%grid, ierr)
      call IODestroy(lbm%io, ierr)
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

      call IOSetFromOptions(lbm%io, options, ierr)
      call GridSetFromOptions(lbm%grid, options, ierr)
      call GridSetName(lbm%grid, lbm%name)
      call FlowSetFromOptions(lbm%flow, options, ierr)
      if (options%transport_disc /= NULL_DISCRETIZATION) then
         lbm%transport => TransportCreate(lbm%comm)
         call TransportSetFromOptions(lbm%transport, options)
      end if

      call BCSetFromOptions(lbm%bc, options, ierr)
    end subroutine LBMSetFromOptions

    subroutine LBMSetUp(lbm)
      type(lbm_type) lbm
      PetscErrorCode ierr
      PetscScalar zero
      zero = 0.d0

      ! set up the DA sizes
      lbm%grid%da_sizes(ONEDOF) = 1
      lbm%grid%da_sizes(NPHASEDOF) = lbm%flow%nphases
      lbm%grid%da_sizes(NPHASEXBDOF) = lbm%flow%nphases*(lbm%flow%disc%b+1)
      lbm%grid%da_sizes(NFLOWDOF) = lbm%flow%ndims
      if (associated(lbm%transport)) then
         lbm%grid%da_sizes(NCOMPONENTDOF) = lbm%transport%ncomponents
         lbm%grid%da_sizes(NCOMPONENTXBDOF) = lbm%transport%ncomponents*(lbm%transport%disc%b+1)
      end if

      call GridSetUp(lbm%grid)
      call FlowSetGrid(lbm%flow, lbm%grid)
      call FlowSetUp(lbm%flow)
      if (associated(lbm%transport)) then
         call TransportSetGrid(lbm%transport, lbm%grid)
         call TransportSetUp(lbm%transport)
      end if
      call BCSetGrid(lbm%bc, lbm%grid)
      call BCSetUp(lbm%bc)

      ! get vectors
      call DMCreateLocalVector(lbm%grid%da(ONEDOF), lbm%walls, ierr)
      call DMCreateGlobalVector(lbm%grid%da(ONEDOF), lbm%walls_g, ierr)

      call VecSet(lbm%walls_g, zero, ierr)
      call VecSet(lbm%walls, zero, ierr)

      call PetscObjectSetName(lbm%walls_g, trim(lbm%name)//'walls', ierr)

      CHKERRQ(ierr)
      return
    end subroutine LBMSetUp

    subroutine LBMRun(lbm, istep, kstep, kwrite)
      ! input
      type(lbm_type) lbm
      integer istep
      integer kstep
      integer kwrite

      ! local
      PetscErrorCode ierr
      integer lcv_step
      type(timing_type),pointer:: timer1, timer2, timer3
      character(len=MAXWORDLENGTH) timerunits
      character(len=MAXWORDLENGTH) timername

      call DMDALocalToLocalBegin(lbm%grid%da(ONEDOF), lbm%walls, INSERT_VALUES, &
           lbm%walls, ierr)
      call DMDALocalToLocalEnd(lbm%grid%da(ONEDOF), lbm%walls, INSERT_VALUES, &
           lbm%walls, ierr)

      if (istep.eq.0) then
         call DMDALocalToLocalBegin(lbm%grid%da(NPHASEXBDOF), lbm%flow%distribution%fi, &
              INSERT_VALUES, lbm%flow%distribution%fi, ierr)
         call DMDALocalToLocalEnd(lbm%grid%da(NPHASEXBDOF), lbm%flow%distribution%fi, &
              INSERT_VALUES, lbm%flow%distribution%fi, ierr)

         ! if (associated(lbm%transport)) then
         !    call DMDALocalToLocalBegin(lbm%grid%da(NCOMPONENTXBDOF), &
         !         lbm%transport%distribution%fi, INSERT_VALUES, &
         !         lbm%transport%distribution%fi, ierr)
         !    call DMDALocalToLocalEnd(lbm%grid%da(NCOMPONENTXBDOF), &
         !         lbm%transport%distribution%fi, INSERT_VALUES, &
         !         lbm%transport%distribution%fi, ierr)
         ! end if

         ! get arrays
         call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
         call FlowGetArrays(lbm%flow)
         call BCGetArrays(lbm%bc, ierr)

         ! update values and view
         if (lbm%grid%info%rank.eq.0) then
            write(*,*) 'outputing step', istep, 'to file', lbm%io%counter
         endif
         call FlowUpdateMoments(lbm%flow, lbm%walls_a)
         call FlowOutputDiagnostics(lbm%flow, lbm%walls_a, lbm%io)

         ! view walls
         call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
         call DMLocalToGlobalBegin(lbm%grid%da(ONEDOF), lbm%walls, INSERT_VALUES, lbm%walls_g, ierr)
         call DMLocalToGlobalEnd(lbm%grid%da(ONEDOF), lbm%walls, INSERT_VALUES, lbm%walls_g, ierr)
         call IOView(lbm%io, lbm%walls_g, 'walls')
         call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)

         ! view grid coordinates
         call GridViewCoordinates(lbm%grid, lbm%io)
      else
         ! just get arrays
         call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
         call FlowGetArrays(lbm%flow)
         call BCGetArrays(lbm%bc, ierr)
      endif
      call IOIncrementCounter(lbm%io)

      timername = trim(lbm%name)//'Simulation'
      timer1 => TimingCreate(lbm%comm, timername)
      do lcv_step = istep+1,kstep
         call FlowStream(lbm%flow)
         call FlowBounceback(lbm%flow, lbm%walls_a)

         call BCApplyFlow(lbm%bc, lbm%walls_a, lbm%flow%distribution)
         call FlowUpdateMoments(lbm%flow, lbm%walls_a)

         ! update rho ghosts values
         call DistributionCommunicateDensity(lbm%flow%distribution)

         ! calculate forces
         call FlowCalcForces(lbm%flow, lbm%walls_a)
         call BCZeroForces(lbm%bc, lbm%flow%forces, lbm%flow%distribution)

         ! collision
         call FlowCollision(lbm%flow, lbm%walls_a)

         ! communicate, update fi
         call DistributionCommunicateFi(lbm%flow%distribution)

         ! check for output
         if(mod(lcv_step,kwrite).eq.0) then
            if (lbm%grid%info%rank.eq.0) then
               write(*,*) 'outputing step', lcv_step, 'to file', lbm%io%counter
            endif
            call FlowOutputDiagnostics(lbm%flow, lbm%walls_a, lbm%io)
            call IOIncrementCounter(lbm%io)
         endif
      end do

      timerunits = 'timestep'
      call TimingEndPerUnit(timer1, (kstep-istep+1), timerunits)
      call TimingDestroy(timer1)

      return
    end subroutine LBMRun

    subroutine LBMInitializeWalls(lbm, init_subroutine)
      type(lbm_type) lbm
      external :: init_subroutine
      PetscErrorCode ierr
      PetscInt vsize

      call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
      call init_subroutine(lbm%walls_a, lbm%options%walls_file, lbm%grid%info)
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
      external :: init_subroutine
      PetscErrorCode ierr

      call DMDAVecGetArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
      call FlowGetArrays(lbm%flow)
      call init_subroutine(lbm%flow%distribution%fi_a, lbm%flow%distribution%rho_a, &
           lbm%flow%distribution%flux, lbm%walls_a, lbm%flow%distribution, &
           lbm%flow%phases, lbm%options)
      call DMDAVecRestoreArrayF90(lbm%grid%da(ONEDOF), lbm%walls, lbm%walls_a, ierr)
      call FlowRestoreArrays(lbm%flow)
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
      
      corners = lbm%grid%info%corners
    end subroutine LBMGetCorners

  subroutine LBMInitializeStateRestarted(lbm, istep, kwrite)
    use petsc
    type(lbm_type) lbm
    PetscInt istep, kwrite
    PetscErrorCode ierr

    lbm%io%counter = istep/kwrite
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


end module LBM_module

