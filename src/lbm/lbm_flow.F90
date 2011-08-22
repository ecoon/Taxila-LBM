!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_flow.F90
!!!     version:         
!!!     created:         17 March 2011
!!!       on:            17:58:06 MDT
!!!     last modified:   22 August 2011
!!!       at:            15:32:20 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"

module LBM_Flow_module
  use petsc
  use LBM_Relaxation_module
  use LBM_Phase_module
  use LBM_Info_module
  use LBM_Discretization_Type_module
  use LBM_Discretization_module
  use LBM_Grid_module
  use LBM_Distribution_Function_type_module
  use LBM_Distribution_Function_module
  use LBM_BC_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public:: flow_type
     MPI_Comm comm
     PetscInt ndims
     PetscInt nphases
     type(discretization_type),pointer :: disc
     type(phase_type),pointer,dimension(:):: phases
     type(grid_type),pointer:: grid
     type(distribution_type),pointer:: distribution
     type(bc_type),pointer:: bc
     
     PetscBool io_fi
     PetscBool io_rho
     Vec rhot_g
     PetscScalar,pointer:: rhot_a(:)
     PetscBool io_rhot

     Vec prs_g
     PetscScalar,pointer:: prs_a(:)
     PetscBool io_prs

     Vec velt_g
     PetscScalar,pointer:: velt_a(:)
     PetscBool io_velt

     PetscScalar,pointer,dimension(:,:,:):: vel_eq
     PetscScalar,pointer,dimension(:,:,:):: fi_eq
     PetscScalar,pointer,dimension(:,:,:):: forces
     PetscScalar,pointer,dimension(:,:):: psi_of_rho
     character(len=MAXWORDLENGTH) name       

     PetscBool :: body_forces
     PetscBool :: fluidfluid_forces
     PetscBool :: use_nonideal_eos
     PetscInt :: use_nonideal_eos_eqn
     PetscBool :: fluidsolid_forces

     PetscScalar,pointer,dimension(:) :: gvt
     PetscScalar velocity_scale
     PetscScalar time_scale
     PetscScalar mass_scale
  end type flow_type

  public :: FlowCreate, &
       FlowDestroy, &
       FlowSetName, &
       FlowSetFromOptions, &
       FlowSetPhysicalScales, &
       FlowSetGrid, &
       FlowSetUp, &
       FlowGetArrays, &
       FlowRestoreArrays, &
       FlowUpdateMoments, &
       FlowUpdateDiagnostics, &
       FlowCalcForces, &
       FlowStream, &
       FlowBounceback, &
       FlowCollision, &
       FlowApplyBCs, &
       FlowOutputDiagnostics, &
       print_a_few

contains
  function FlowCreate(comm) result(flow)
    MPI_Comm comm
    type(flow_type),pointer:: flow

    allocate(flow)
    flow%comm = comm
    flow%nphases = -1
    flow%ndims = -1
    nullify(flow%phases)
    nullify(flow%grid)
    flow%disc => DiscretizationCreate(flow%comm)
    flow%distribution => DistributionCreate(flow%comm)
    flow%bc => BCCreate(flow%comm)

    flow%velt_g = 0
    flow%prs_g = 0
    flow%rhot_g = 0
    
    nullify(flow%velt_a)
    nullify(flow%prs_a)
    nullify(flow%rhot_a)
    
    nullify(flow%vel_eq)
    nullify(flow%forces)
    nullify(flow%fi_eq)
    nullify(flow%psi_of_rho)

    flow%io_prs = PETSC_TRUE
    flow%io_rhot = PETSC_FALSE
    flow%io_rho = PETSC_TRUE
    flow%io_fi = PETSC_TRUE
    flow%io_velt = PETSC_TRUE

    flow%body_forces = PETSC_FALSE
    flow%fluidfluid_forces = PETSC_FALSE
    flow%use_nonideal_eos = PETSC_FALSE
    flow%use_nonideal_eos_eqn = EOS_DENSITY
    flow%fluidsolid_forces = PETSC_FALSE
    nullify(flow%gvt)
  end function FlowCreate

  subroutine FlowDestroy(flow, ierr)
    type(flow_type) flow
    PetscErrorCode ierr
    integer lcv

    if (associated(flow%disc)) call DiscretizationDestroy(flow%disc, ierr)
    if (associated(flow%phases)) then
       do lcv=1,flow%nphases
          call PhaseDestroy(flow%phases(lcv),ierr)
       end do
    end if
    if (associated(flow%distribution)) call DistributionDestroy(flow%distribution, ierr)
    if (associated(flow%bc)) call BCDestroy(flow%bc, ierr)

    if (flow%velt_g /= 0) call VecDestroy(flow%velt_g,ierr)
    if (flow%prs_g /= 0) call VecDestroy(flow%prs_g,ierr)
    if (flow%rhot_g /= 0) call VecDestroy(flow%rhot_g,ierr)
    if (associated(flow%vel_eq)) deallocate(flow%vel_eq)
    if (associated(flow%forces)) deallocate(flow%forces)
    if (associated(flow%fi_eq)) deallocate(flow%fi_eq)
    if (associated(flow%psi_of_rho)) deallocate(flow%psi_of_rho)
    if (associated(flow%gvt)) deallocate(flow%gvt)
  end subroutine FlowDestroy

  subroutine FlowSetName(flow, name) 
    type(flow_type) flow 
    character(len=MAXWORDLENGTH):: name       
    flow%name = name
    call DistributionSetName(flow%distribution, name)
  end subroutine FlowSetName
    
  subroutine FlowSetFromOptions(flow, options, ierr)
    use String_module
    use LBM_Options_module
    type(flow_type) flow
    type(options_type) options
    PetscErrorCode ierr
    integer lcv, lcv2
    PetscBool help
    PetscInt nmax
    PetscBool bcvalue
    PetscScalar gravity(options%ndims)
    PetscScalar,parameter:: eps=1.e-15 ! slightly larger than machine epsilon
    PetscBool flag
    character(len=MAXWORDLENGTH) eos_name
    character(len=MAXWORDLENGTH) test_eos_name

    flag = PETSC_FALSE
    flow%nphases = options%nphases
    flow%ndims = options%ndims
    flow%phases => PhaseCreate(flow%comm, flow%nphases)

    call DiscretizationSetType(flow%disc, options%flow_disc)
    call DiscretizationSetSizes(flow%disc, flow%grid%info%stencil_size)
    call DiscretizationSetUp(flow%disc)
    
    do lcv=1,flow%nphases
       call PhaseSetSizes(flow%phases(lcv), flow%nphases, flow%disc%b)
       call PhaseSetID(flow%phases(lcv), lcv)
       call PhaseSetFromOptions(flow%phases(lcv), options, ierr)
    end do

    ! set up control for forcing terms
    gravity(:) = 0.d0
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)
    if (help) call PetscPrintf(options%comm, "-gvt=<0,0,0>: gravity\n", ierr)
    nmax = flow%ndims
    call PetscOptionsGetRealArray(options%my_prefix, '-gvt', gravity, &
         nmax, flow%body_forces, ierr)
    if (flow%body_forces) then
      allocate(flow%gvt(flow%ndims))
      flow%gvt = gravity
    end if

    do lcv=1,flow%nphases
      do lcv2=1,flow%nphases
        if (ABS(flow%phases(lcv)%gf(lcv2)) > eps) then
          flow%fluidfluid_forces = PETSC_TRUE
        end if
      end do
    end do

    if (help) call PetscPrintf(options%comm, "-flow_use_nonideal_eos <eos-type>", &
         ierr)
    call PetscOptionsGetString(options%my_prefix, '-flow_use_nonideal_eos', &
         eos_name, flag, ierr)
    if (flag) then
      flow%use_nonideal_eos = PETSC_TRUE

      test_eos_name = 'sc'
      if (StringCompareIgnoreCase(test_eos_name, eos_name, 2)) then
         flow%use_nonideal_eos_eqn = EOS_SC
      end if

      test_eos_name = 'pr'
      if (StringCompareIgnoreCase(test_eos_name, eos_name, 2)) then
         flow%use_nonideal_eos_eqn = EOS_PR
      end if
    end if

    flow%fluidsolid_forces = options%flow_fluidsolid_forces
    
    ! set up the vectors for holding boundary data
    ! dimension 
    call BCSetSizes(flow%bc, options%ndims*options%nphases)
    call BCSetFromOptions(flow%bc, options, ierr)

    ! parse all flow boundary conditions
    if (help) call PetscPrintf(options%comm, "-bc_pressure_{xyz}{mp}: use pressure bcs\n",&
         ierr)
    if (help) call PetscPrintf(options%comm, "-bc_velocity_{xyz}{mp}: use total"//&
         " velocity bcs\n", ierr)
    if (help) call PetscPrintf(options%comm, "-bc_velocity_{xyz}{mp}_poiseuille:"// &
         " use total velocity bcs with a poiseuille profile\n", ierr)
    if (help) call PetscPrintf(options%comm, "-bc_flux_{xyz}{mp}: use VOLUMETRIC flux"//&
         " bcs (NOT mass flux!)\n", ierr)

    ! xm boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_xm', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_XM) = BC_DIRICHLET

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_xm', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_XM) = BC_FLUX

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_xm', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_XM) = BC_VELOCITY

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_xm_poiseuille', bcvalue, &
         flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_XM) = BC_VELOCITY

    ! xp boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_xp', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_XP) = BC_DIRICHLET

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_xp', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_XP) = BC_FLUX

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_xp', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_XP) = BC_VELOCITY

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_xp_poiseuille', bcvalue, &
         flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_XP) = BC_VELOCITY

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_x', bcvalue, flag, ierr)
    if (bcvalue) then
       flow%bc%flags(BOUNDARY_XM) = BC_PERIODIC
       flow%bc%flags(BOUNDARY_XP) = BC_PERIODIC
    end if

    ! ym boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_ym', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_YM) = BC_DIRICHLET

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_ym', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_YM) = BC_FLUX

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_ym', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_YM) = BC_VELOCITY

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_ym_poiseuille', bcvalue, &
         flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_YM) = BC_VELOCITY

    ! yp boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_yp', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_YP) = BC_DIRICHLET
    
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_yp', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_YP) = BC_FLUX

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_yp', bcvalue, flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_YP) = BC_VELOCITY

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_yp_poiseuille', bcvalue, &
         flag, ierr)
    if (bcvalue) flow%bc%flags(BOUNDARY_YP) = BC_VELOCITY

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_y', bcvalue, flag, ierr)
    if (bcvalue) then
       flow%bc%flags(BOUNDARY_YM) = BC_PERIODIC
       flow%bc%flags(BOUNDARY_YP) = BC_PERIODIC
    end if

    if (options%ndims > 2) then
       ! zm boundary
       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_zm', bcvalue, flag, ierr)
       if (bcvalue) flow%bc%flags(BOUNDARY_ZM) = BC_DIRICHLET

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_flux_zm', bcvalue, flag, ierr)
       if (bcvalue) flow%bc%flags(BOUNDARY_ZM) = BC_FLUX

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_zm', bcvalue, flag, ierr)
       if (bcvalue) flow%bc%flags(BOUNDARY_ZM) = BC_VELOCITY

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_zm_poiseuille', bcvalue, &
            flag, ierr)
       if (bcvalue) flow%bc%flags(BOUNDARY_ZM) = BC_VELOCITY

       ! zp boundary
       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_zp', bcvalue, flag, ierr)
       if (bcvalue) flow%bc%flags(BOUNDARY_ZP) = BC_DIRICHLET

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_flux_zp', bcvalue, flag, ierr)
       if (bcvalue) flow%bc%flags(BOUNDARY_ZP) = BC_FLUX

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_zp', bcvalue, flag, ierr)
       if (bcvalue) flow%bc%flags(BOUNDARY_ZP) = BC_VELOCITY

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_velocity_zp_poiseuille', bcvalue, &
            flag, ierr)
       if (bcvalue) flow%bc%flags(BOUNDARY_ZP) = BC_VELOCITY

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_z', bcvalue, flag, ierr)
       if (bcvalue) then
          flow%bc%flags(BOUNDARY_ZM) = BC_PERIODIC
          flow%bc%flags(BOUNDARY_ZP) = BC_PERIODIC
       end if
    end if
  end subroutine FlowSetFromOptions

  subroutine FlowSetGrid(flow, grid)
    use LBM_Grid_module
    type(flow_type) flow
    type (grid_type),pointer:: grid
    flow%grid => grid
    call BCSetGrid(flow%bc, flow%grid)
  end subroutine FlowSetGrid

  subroutine FlowSetPhysicalScales(flow, ierr)
    type(flow_type) flow
    PetscScalar consistent_tau, consistent_mm
    PetscErrorCode ierr
    PetscScalar,parameter:: eps=1.e-8
    PetscInt lcv

    ! deal with consistency.  Assume phase 1's tau, mm are correct
    ! viscosity/time scale
    if (flow%phases(1)%viscosity < -990.) then
       ! tau is correct, viscosity is not given
       flow%phases(1)%time_scale = 1.d0
       flow%phases(1)%viscosity = flow%grid%length_scale * &
         flow%grid%length_scale * (flow%phases(1)%relax%tau - 0.5d0)/3.d0
    else        
       flow%phases(1)%time_scale = flow%grid%length_scale * &
            flow%grid%length_scale * (flow%phases(1)%relax%tau - 0.5d0)/ &
            (flow%phases(1)%viscosity*3.d0)
    end if
    flow%time_scale = flow%phases(1)%time_scale
    flow%velocity_scale = flow%grid%length_scale/flow%time_scale

    ! density
    if (flow%phases(1)%density < -990.) then
       flow%phases(1)%density = 1.d0
    end if
    flow%mass_scale = flow%phases(1)%mm*flow%phases(1)%density*(flow%grid%length_scale**3)

    if (flow%nphases > 1) then
      do lcv=2,flow%nphases
        ! Assert correct viscosity ratios
        if (flow%phases(lcv)%viscosity < -990.d0) then
           flow%phases(lcv)%viscosity = flow%grid%length_scale * &
                flow%grid%length_scale * &
                (flow%phases(1)%relax%tau - 0.5d0)/(3.d0*flow%time_scale)
        else 
          consistent_tau = (3.d0*flow%time_scale*flow%phases(lcv)%viscosity) &
                /(flow%grid%length_scale* &
                flow%grid%length_scale) + 0.5d0
          if (abs(flow%phases(lcv)%relax%tau - 1.d0) < eps) then
             ! tau not specified, so set it
             flow%phases(lcv)%relax%tau = consistent_tau
             if (flow%grid%info%rank .eq. 0) then
                print*, '  Setting phase', TRIM(flow%phases(lcv)%name), 'tau = ', &
                     consistent_tau
             end if
          else if (abs(flow%phases(lcv)%relax%tau - consistent_tau) > eps) then
             SETERRQ(PETSC_COMM_WORLD, 1, 'Viscosities and relaxation times specified '// &
                  'for phases are not consistent.', ierr)
          else
             ! all ok
          end if
        end if

        ! Assert correct density ratios
        if (flow%phases(lcv)%density < -990.d0) then
          flow%phases(lcv)%density = flow%mass_scale/flow%phases(1)%mm / &
               (flow%grid%length_scale**3)
        else
          consistent_mm = flow%mass_scale / flow%phases(1)%density / &
               (flow%grid%length_scale**3)
          if (abs(flow%phases(lcv)%mm - 1.d0) < eps) then
             ! mm not specified, so set it
             flow%phases(lcv)%mm = consistent_mm
             if (flow%grid%info%rank .eq. 0) then
                print*, '  Setting molecular mass', TRIM(flow%phases(lcv)%name), 'mm = ', &
                     consistent_mm
             end if
          else if (abs(flow%phases(lcv)%mm - consistent_mm) > eps) then
             SETERRQ(PETSC_COMM_WORLD, 1, 'Densities and molecular masses specified '// &
                  'for phases are not consistent.', ierr)
          else
             ! all ok
          end if
        end if
      end do
    end if

    
    if (flow%grid%info%rank .eq. 0) then
       print*, 'Scale report from LBM:'
       print*, '  length scale [m]:', flow%grid%length_scale
       print*, '  time scale [s]:', flow%time_scale
       print*, '  velocity scale [m]:', flow%velocity_scale
       print*, '  mass scale [kg]:', flow%mass_scale
    end if
  end subroutine FlowSetPhysicalScales

  subroutine FlowSetUp(flow) 
    type(flow_type) flow
    integer lcv
    PetscErrorCode ierr
    PetscScalar zero
    zero = 0.d0

    do lcv=1,flow%nphases
       call DiscretizationSetUpRelax(flow%disc, flow%phases(lcv)%relax)
    end do
    call DistributionSetInfo(flow%distribution, flow%grid%info)
    call DistributionSetDiscretization(flow%distribution, flow%disc)
    call DistributionSetSizes(flow%distribution, flow%nphases)
    call DistributionSetDAs(flow%distribution, flow%grid%da(NPHASEXBDOF), &
         flow%grid%da(NPHASEDOF))
    call DistributionSetUp(flow%distribution)
    call BCSetUp(flow%bc)

    ! allocate, initialize workspace
    allocate(flow%vel_eq(1:flow%nphases, 1:flow%ndims, &
         1:flow%grid%info%gxyzl))
    allocate(flow%forces(1:flow%nphases, 1:flow%ndims, &
         1:flow%grid%info%gxyzl))
    allocate(flow%fi_eq(1:flow%nphases, 0:flow%disc%b, &
         1:flow%grid%info%gxyzl))
    allocate(flow%psi_of_rho(flow%nphases,flow%grid%info%rgxyzl))
    flow%vel_eq = 0.d0
    flow%forces = 0.d0
    flow%fi_eq = 0.d0

    call DMCreateGlobalVector(flow%grid%da(ONEDOF), flow%prs_g, ierr)
    call DMCreateGlobalVector(flow%grid%da(ONEDOF), flow%rhot_g, ierr)
    call DMCreateGlobalVector(flow%grid%da(NFLOWDOF), flow%velt_g, ierr)

    call VecSet(flow%prs_g, zero, ierr)
    call VecSet(flow%rhot_g, zero, ierr)
    call VecSet(flow%velt_g, zero, ierr)

    call PetscObjectSetName(flow%prs_g, trim(flow%name)//'prs', ierr)
    call PetscObjectSetName(flow%rhot_g, trim(flow%name)//'rhot', ierr)
    call PetscObjectSetName(flow%velt_g, trim(flow%name)//'velt', ierr)
   end subroutine FlowSetUp

  subroutine FlowGetArrays(flow, ierr)
    type(flow_type) flow
    PetscErrorCode ierr
    call DistributionGetArrays(flow%distribution, ierr)
    call BCGetArrays(flow%bc, ierr)
  end subroutine FlowGetArrays

  subroutine FlowRestoreArrays(flow, ierr)
    type(flow_type) flow
    PetscErrorCode ierr
    call DistributionRestoreArrays(flow%distribution, ierr)
    call BCRestoreArrays(flow%bc, ierr)
  end subroutine FlowRestoreArrays

  subroutine FlowUpdateMoments(flow, walls)
    type(flow_type) flow
    PetscScalar,dimension(1:flow%grid%info%gxyzl):: walls
    call DistributionCalcDensity(flow%distribution, walls)
    call DistributionCalcFlux(flow%distribution, walls)
  end subroutine FlowUpdateMoments

  subroutine FlowOutputDiagnostics(flow, walls, io)
    use LBM_IO_module
    type(flow_type) flow
    type(io_type) io
    PetscScalar,dimension(1:flow%grid%info%gxyzl):: walls
    PetscErrorCode ierr

    if (flow%io_fi) then
       call DMDAVecRestoreArrayF90(flow%grid%da(NPHASEXBDOF),flow%distribution%fi, &
            flow%distribution%fi_a, ierr)
       call DMLocalToGlobalBegin(flow%grid%da(NPHASEXBDOF),flow%distribution%fi, &
            INSERT_VALUES, flow%distribution%fi_g, ierr)
       call DMLocalToGlobalEnd(flow%grid%da(NPHASEXBDOF),flow%distribution%fi, &
            INSERT_VALUES, flow%distribution%fi_g, ierr)
       call DMDAVecGetArrayF90(flow%grid%da(NPHASEXBDOF),flow%distribution%fi, &
            flow%distribution%fi_a, ierr)
       call IOView(io, flow%distribution%fi_g, 'fi')
    end if
    if (flow%io_rho) then
       call DMDAVecRestoreArrayF90(flow%grid%da(NPHASEDOF),flow%distribution%rho, &
            flow%distribution%rho_a, ierr)
       call DMLocalToGlobalBegin(flow%grid%da(NPHASEDOF),flow%distribution%rho, &
            INSERT_VALUES, flow%distribution%rho_g, ierr)
       call DMLocalToGlobalEnd(flow%grid%da(NPHASEDOF),flow%distribution%rho, &
            INSERT_VALUES, flow%distribution%rho_g, ierr)
       call VecScale(flow%distribution%rho_g, &
            flow%mass_scale/(flow%grid%length_scale**3), ierr)
       call DMDAVecGetArrayF90(flow%grid%da(NPHASEDOF),flow%distribution%rho, &
            flow%distribution%rho_a, ierr)
       call IOView(io, flow%distribution%rho_g, 'rho')
       call VecScale(flow%distribution%rho_g, &
            (flow%grid%length_scale**3)/flow%mass_scale, ierr)
    end if

    call VecScale(flow%velt_g, flow%velocity_scale, ierr)
    if (flow%io_velt) call IOView(io, flow%velt_g, 'u')
    call VecScale(flow%velt_g, 1.d0/flow%velocity_scale, ierr)

    if (flow%io_rhot) call IOView(io, flow%rhot_g, 'rhot')
    if (flow%io_prs) call IOView(io, flow%prs_g, 'prs')
  end subroutine FlowOutputDiagnostics

  subroutine FlowUpdateDiagnostics(flow, walls)
    type(flow_type) flow
    PetscScalar,dimension(1:flow%grid%info%gxyzl):: walls
    PetscErrorCode ierr

    call DMDAVecGetArrayF90(flow%grid%da(NFLOWDOF),flow%velt_g,flow%velt_a,ierr)
    call DMDAVecGetArrayF90(flow%grid%da(ONEDOF), flow%rhot_g, flow%rhot_a,ierr)
    call DMDAVecGetArrayF90(flow%grid%da(ONEDOF), flow%prs_g, flow%prs_a, ierr)

    if (flow%use_nonideal_eos) then
      call FlowCalcPsiOfRho(flow, flow%distribution%rho_a)
    end if

    select case(flow%ndims)
    case(2)
      if (flow%use_nonideal_eos) then
        call FlowUpdateDiagnosticsD2(flow, flow%distribution%rho_a, flow%psi_of_rho, &
             flow%distribution%flux, flow%forces, walls, flow%rhot_a, &
             flow%prs_a, flow%velt_a)
      else
        call FlowUpdateDiagnosticsD2(flow, flow%distribution%rho_a, flow%distribution%rho_a, &
             flow%distribution%flux, flow%forces, walls, flow%rhot_a, &
             flow%prs_a, flow%velt_a)
      end if
    case(3)
      if (flow%use_nonideal_eos) then
        call FlowUpdateDiagnosticsD3(flow, flow%distribution%rho_a, flow%psi_of_rho, &
             flow%distribution%flux, flow%forces, walls, flow%rhot_a, &
             flow%prs_a, flow%velt_a)
      else
        call FlowUpdateDiagnosticsD3(flow, flow%distribution%rho_a, flow%distribution%rho_a, &
             flow%distribution%flux, flow%forces, walls, flow%rhot_a, &
             flow%prs_a, flow%velt_a)
      end if
    end select
    call DMDAVecRestoreArrayF90(flow%grid%da(NFLOWDOF),flow%velt_g,flow%velt_a,ierr)
    call DMDAVecRestoreArrayF90(flow%grid%da(ONEDOF), flow%rhot_g, flow%rhot_a, ierr)
    call DMDAVecRestoreArrayF90(flow%grid%da(ONEDOF), flow%prs_g, flow%prs_a, ierr)

  end subroutine FlowUpdateDiagnostics

  subroutine FlowUpdateDiagnosticsD3(flow, rho, psi, vel, forces, walls, rhot, prs, velt)
    type(flow_type) flow
    PetscScalar,dimension(flow%grid%info%gxs:flow%grid%info%gxe, &
         flow%grid%info%gys:flow%grid%info%gye, &
         flow%grid%info%gzs:flow%grid%info%gze):: walls
    PetscScalar,dimension(flow%nphases, &
         flow%grid%info%rgxs:flow%grid%info%rgxe, &
         flow%grid%info%rgys:flow%grid%info%rgye, &
         flow%grid%info%rgzs:flow%grid%info%rgze):: rho, psi
    PetscScalar,dimension(flow%nphases, flow%ndims, &
         flow%grid%info%gxs:flow%grid%info%gxe, &
         flow%grid%info%gys:flow%grid%info%gye, &
         flow%grid%info%gzs:flow%grid%info%gze):: vel,forces
    PetscScalar,dimension(flow%grid%info%xs:flow%grid%info%xe, &
         flow%grid%info%ys:flow%grid%info%ye, &
         flow%grid%info%zs:flow%grid%info%ze):: rhot,prs
    PetscScalar,dimension(flow%ndims, &
         flow%grid%info%xs:flow%grid%info%xe, &
         flow%grid%info%ys:flow%grid%info%ye, &
         flow%grid%info%zs:flow%grid%info%ze):: velt

    PetscInt i,j,k,m,d
    PetscScalar mm(1:flow%nphases)
    
    do m=1,flow%nphases
       mm(m) = flow%phases(m)%mm
    end do

    do k=flow%grid%info%zs,flow%grid%info%ze
    do j=flow%grid%info%ys,flow%grid%info%ye
    do i=flow%grid%info%xs,flow%grid%info%xe
    if (walls(i,j,k).eq.0) then
       rhot(i,j,k) = sum(rho(:,i,j,k)*mm,1)
       do d=1,flow%ndims
          velt(d,i,j,k) = (sum(vel(:,d,i,j,k)*mm,1) + 0.5*sum(forces(:,d,i,j,k),1))/rhot(i,j,k)
       end do
       prs(i,j,k) = rhot(i,j,k)/3.d0
       if (flow%nphases > 1) then
          do m=1,flow%nphases
            ! In Qinjun's original code the pressure was calculated as:
            ! prs(i,j) = prs(i,j) + 3*rho(m,i,j)*sum(rho(:,i,j)*flow%phases(m)%gf,1)
            ! However, with the HOD for fluid-fluid the new gf is 1/6 of the old gf, so
            ! I have accounted for it in the new pressure calculation.  This equation
            ! has been verified by Laplace's law for viscosity ratios of 1 and 5. 
             prs(i,j,k) = prs(i,j,k) + flow%disc%c_0/2.d0*psi(m,i,j,k) &
                  *sum(flow%phases(m)%gf*psi(:,i,j,k),1)
          end do
       end if
    end if
    end do
    end do
    end do
  end subroutine FlowUpdateDiagnosticsD3

  subroutine FlowUpdateDiagnosticsD2(flow, rho, psi, vel, forces, walls, rhot, prs, velt)
    type(flow_type) flow
    PetscScalar,dimension(flow%grid%info%gxs:flow%grid%info%gxe, &
         flow%grid%info%gys:flow%grid%info%gye):: walls
    PetscScalar,dimension(flow%nphases, &
         flow%grid%info%rgxs:flow%grid%info%rgxe, &
         flow%grid%info%rgys:flow%grid%info%rgye):: rho, psi
    PetscScalar,dimension(flow%nphases, flow%ndims, &
         flow%grid%info%gxs:flow%grid%info%gxe, &
         flow%grid%info%gys:flow%grid%info%gye):: vel,forces
    PetscScalar,dimension(flow%grid%info%xs:flow%grid%info%xe, &
         flow%grid%info%ys:flow%grid%info%ye):: rhot,prs
    PetscScalar,dimension(flow%ndims, &
         flow%grid%info%xs:flow%grid%info%xe, &
         flow%grid%info%ys:flow%grid%info%ye):: velt

    PetscInt i,j,m,d
    PetscScalar mm(1:flow%nphases)
    
    do m=1,flow%nphases
       mm(m) = flow%phases(m)%mm
    end do
    do j=flow%grid%info%ys,flow%grid%info%ye
    do i=flow%grid%info%xs,flow%grid%info%xe
    if (walls(i,j).eq.0) then
       rhot(i,j) = sum(rho(:,i,j)*mm,1)
       do d=1,flow%ndims
          velt(d,i,j) = (sum(vel(:,d,i,j)*mm,1) + 0.5*sum(forces(:,d,i,j),1))/rhot(i,j)
       end do
       prs(i,j) = rhot(i,j)/3.d0
       if (flow%nphases > 1) then
          do m=1,flow%nphases
            ! In Qinjun's original code the pressure was calculated as:
            ! prs(i,j) = prs(i,j) + 1.5*rho(m,i,j)*sum(rho(:,i,j)*flow%phases(m)%gf,1)
            ! However, with the hod for fluid-fluid the new gf is 1/3 of the old gf, so
            ! I have accounted for it in the new pressure calculation.  This equation
            ! has been verified by Laplace's law for viscosity ratios of 1 and 5.        
            prs(i,j) = prs(i,j) + flow%disc%c_0/2.d0*psi(m,i,j) &
                  *sum(flow%phases(m)%gf*psi(:,i,j),1)
          end do
       end if
    end if
    end do
    end do
  end subroutine FlowUpdateDiagnosticsD2

  subroutine FlowCalcForces(flow, walls)
    use LBM_Walls_module
    use LBM_Forcing_module
    use LBM_Logging_module
    type(flow_type) flow
    type(walls_type) walls
    PetscErrorCode ierr
    
    flow%forces = 0.d0
    if (flow%fluidfluid_forces) then
      call PetscLogEventBegin(logger%event_forcing_fluidfluid,ierr)
      if (flow%use_nonideal_eos) then
        call FlowCalcPsiOfRho(flow, flow%distribution%rho_a)
        call LBMAddFluidFluidForces(flow%distribution, flow%phases, &
             flow%psi_of_rho, walls%walls_a, flow%forces)
      else
        call LBMAddFluidFluidForces(flow%distribution, flow%phases, &
             flow%distribution%rho_a, walls%walls_a, flow%forces)
      end if
      call PetscLogEventEnd(logger%event_forcing_fluidfluid,ierr)
    end if

    if (flow%fluidsolid_forces) then
      call PetscLogEventBegin(logger%event_forcing_fluidsolid,ierr)
      call LBMAddFluidSolidForces(flow%distribution, flow%phases, walls, &
           flow%distribution%rho_a, flow%forces)
      call PetscLogEventEnd(logger%event_forcing_fluidsolid,ierr)
    end if

    if (flow%body_forces) then
      call PetscLogEventBegin(logger%event_forcing_body,ierr)
      call LBMAddBodyForces(flow%distribution, flow%phases, flow%gvt, &
           flow%distribution%rho_a, walls%walls_a, flow%forces)
      call PetscLogEventEnd(logger%event_forcing_body,ierr)
    end if
  end subroutine FlowCalcForces

  subroutine FlowStream(flow)
    type(flow_type) flow
    call DistributionStream(flow%distribution)
  end subroutine FlowStream

  subroutine FlowBounceback(flow, walls)
    type(flow_type) flow
    PetscScalar,dimension(flow%grid%info%gxyzl):: walls
    call DistributionBounceback(flow%distribution, walls)
  end subroutine FlowBounceback

  subroutine FlowCollision(flow, walls)
    type(flow_type) flow
    PetscScalar,dimension(flow%grid%info%gxyzl):: walls
    select case(flow%ndims)
    case(2)
       call FlowCollisionD2(flow, flow%distribution%fi_a, flow%distribution%rho_a, &
            flow%distribution%flux, flow%forces, walls, flow%fi_eq, flow%distribution)
    case(3)
       call FlowCollisionD3(flow, flow%distribution%fi_a, flow%distribution%rho_a, &
            flow%distribution%flux, flow%forces, walls, flow%fi_eq, flow%distribution)
    end select
  end subroutine FlowCollision

  subroutine FlowCollisionD3(flow, fi, rho, u, forces, walls, fi_eq, dist)
    use LBM_Logging_module
    type(flow_type) flow
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(flow%nphases, 0:flow%disc%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze)::fi,fi_eq
    PetscScalar,dimension(1:dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: u,forces
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls

    PetscScalar,dimension(dist%info%ndims,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze ):: up
    PetscScalar,dimension(dist%s,dist%info%ndims,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: ue

    PetscInt m,i,j,k,d
    PetscScalar,parameter:: eps=1.e-15 ! slightly larger than machine epsilon
    PetscScalar,dimension(1:dist%s) :: mmot, mm
    PetscErrorCode ierr

    call PetscLogEventBegin(logger%event_collision_precalc,ierr)
    do m=1,dist%s
       mmot(m) = flow%phases(m)%mm/flow%phases(m)%relax%tau
       mm(m) = flow%phases(m)%mm
    end do

    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j,k).eq.0) then
       do d=1,dist%info%ndims
          up(d,i,j,k) = sum(u(:,d,i,j,k)*mmot,1)/sum(rho(:,i,j,k)*mmot,1)
       end do
       do m=1,dist%s
          ue(m,:,i,j,k) = up(:,i,j,k) + forces(m,:,i,j,k)/(rho(m,i,j,k)*mmot(m)+eps)
       end do
    end if
    end do
    end do
    end do
    call PetscLogEventEnd(logger%event_collision_precalc,ierr)

    do m=1,dist%s
      call PetscLogEventBegin(logger%event_collision_feq,ierr)
      call DiscretizationEquilf(flow%disc, rho, ue, &
           walls, fi_eq, m, flow%phases(m)%relax, dist)
      call PetscLogEventEnd(logger%event_collision_feq,ierr)

      call PetscLogEventBegin(logger%event_collision_relax,ierr)
      call RelaxationCollide(flow%phases(m)%relax, fi, fi_eq, &
           walls, m, dist)
      call PetscLogEventEnd(logger%event_collision_relax,ierr)
    end do
  end subroutine FlowCollisionD3

  subroutine FlowCollisionD2(flow, fi, rho, u, forces, walls, fi_eq, dist)
    use LBM_Logging_module
    type(flow_type) flow
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(flow%nphases, 0:flow%disc%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye)::fi,fi_eq
    PetscScalar,dimension(1:dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u,forces
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls

    PetscScalar,dimension(dist%info%ndims,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye ):: up
    PetscScalar,dimension(dist%s,dist%info%ndims,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: ue

    PetscInt m,i,j,d
    PetscScalar,parameter:: eps=1.e-15 ! slightly larger than machine epsilon
    PetscScalar,dimension(1:dist%s) :: mmot, mm
    PetscErrorCode ierr

    call PetscLogEventBegin(logger%event_collision_precalc,ierr)
    do m=1,dist%s
       mmot(m) = flow%phases(m)%mm/flow%phases(m)%relax%tau
       mm(m) = flow%phases(m)%mm
    end do

    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j).eq.0) then
       do d=1,dist%info%ndims
          up(d,i,j) = sum(u(:,d,i,j)*mmot,1)/sum(rho(:,i,j)*mmot,1)
       end do
       do m=1,dist%s
          ue(m,:,i,j) = up(:,i,j) + forces(m,:,i,j)/(rho(m,i,j)*mmot(m)+eps)
       end do
    end if
    end do
    end do
    call PetscLogEventEnd(logger%event_collision_precalc,ierr)

    do m=1,dist%s
      call PetscLogEventBegin(logger%event_collision_feq,ierr)
      call DiscretizationEquilf(flow%disc, rho, ue, walls, &
           fi_eq, m, flow%phases(m)%relax, dist)
      call PetscLogEventEnd(logger%event_collision_feq,ierr)

      call PetscLogEventBegin(logger%event_collision_relax,ierr)
      call RelaxationCollide(flow%phases(m)%relax, fi, fi_eq, &
           walls, m, dist)
      call PetscLogEventEnd(logger%event_collision_relax,ierr)
    end do
  end subroutine FlowCollisionD2

  subroutine FlowApplyBCs(flow, walls)
    type(flow_type) flow
    PetscScalar,dimension(flow%grid%info%gxyzl):: walls
    call BCApply(flow%bc, walls, flow%distribution)
  end subroutine FlowApplyBCs

  subroutine FlowCalcPsiOfRho(flow, rho)
    type(flow_type) flow
    PetscScalar,dimension(flow%distribution%s,flow%grid%info%rgxyzl):: rho
    PetscErrorCode ierr
    PetscInt m,i

    select case(flow%use_nonideal_eos_eqn)
    case(EOS_SC)
      do i=1,flow%grid%info%gxyzl
        do m=1,flow%distribution%s
          flow%psi_of_rho(m,i) = 1 - EXP(-rho(m,i))
        end do
      end do

    case(EOS_PR)
      do i=1,flow%grid%info%gxyzl
        do m=1,flow%distribution%s

        end do
      end do
      
    case DEFAULT
      SETERRQ(1, 1, 'invalid equation of state', ierr)
    end select


  end subroutine FlowCalcPsiOfRho

  subroutine print_a_few(fi, u, forces, dist,istep)
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s,0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: u,forces

    PetscInt j,k,istep

    print*, 'step:', istep
    print*, '---------------------'
    j=63
    k=19
    print*, 'outer:'
    print*, 'fi(1,:):', fi(1,:,1,j,k)
    print*, 'fi(2,:):', fi(2,:,1,j,k)
    ! print*, 'fi(1,0):', fi(1,0,1,j,k)
    ! print*, 'fi(1,5):', fi(1,5,1,j,k)
    ! print*, 'fi(2,0):', fi(2,0,1,j,k)
    ! print*, 'fi(2,5):', fi(2,5,1,j,k)
    print*, 'u(1,z):', u(1,3,1,j,k)
    print*, 'u(2,z):', u(2,3,1,j,k)
    print*, 'forces(1,z):', forces(1,3,1,j,k)
    print*, 'forces(2,z):', forces(2,3,1,j,k)
    print*, '---------------------'

    j=63
    k=20
    print*, 'inner:'
    print*, 'fi(1,:):', fi(1,:,1,j,k)
    print*, 'fi(2,:):', fi(2,:,1,j,k)
    ! print*, 'fi(1,0):', fi(1,0,1,j,k)
    ! print*, 'fi(1,5):', fi(1,5,1,j,k)
    ! print*, 'fi(2,0):', fi(2,0,1,j,k)
    ! print*, 'fi(2,5):', fi(2,5,1,j,k)
    print*, 'u(1,z):', u(1,3,1,j,k)
    print*, 'u(2,z):', u(2,3,1,j,k)
    print*, 'forces(1,z):', forces(1,3,1,j,k)
    print*, 'forces(2,z):', forces(2,3,1,j,k)
    print*, '========================'
  end subroutine print_a_few
end module LBM_Flow_module
