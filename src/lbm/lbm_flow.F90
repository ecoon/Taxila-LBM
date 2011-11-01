!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_flow.F90
!!!     version:         
!!!     created:         17 March 2011
!!!       on:            17:58:06 MDT
!!!     last modified:   01 November 2011
!!!       at:            15:40:10 MDT
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
  use LBM_EOS_module
  use LBM_Component_module
  use LBM_Info_module
  use LBM_Discretization_Type_module
  use LBM_Discretization_module
  use LBM_Grid_module
  use LBM_Distribution_Function_type_module
  use LBM_Distribution_Function_module
  use LBM_Walls_module
  use LBM_BC_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public:: flow_type
     MPI_Comm comm
     PetscInt ndims
     PetscInt ncomponents
     type(discretization_type),pointer :: disc
     type(component_type),pointer,dimension(:):: components
     type(grid_type),pointer:: grid
     type(distribution_type),pointer:: distribution
     type(bc_type),pointer:: bc
     
     PetscBool io_fi
     PetscBool io_last_fi
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
       FlowStream, &
       FlowBounceback, &
       FlowCollision, &
       FlowApplyBCs, &
       FlowOutputDiagnostics, &
       FlowFiInit, &
       print_a_few

contains
  function FlowCreate(comm) result(flow)
    MPI_Comm comm
    type(flow_type),pointer:: flow

    allocate(flow)
    flow%comm = comm
    flow%ncomponents = -1
    flow%ndims = -1
    nullify(flow%components)
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
    flow%io_fi = PETSC_FALSE
    flow%io_last_fi = PETSC_FALSE
    flow%io_velt = PETSC_TRUE

    flow%body_forces = PETSC_FALSE
    flow%fluidfluid_forces = PETSC_FALSE
    flow%use_nonideal_eos = PETSC_FALSE
    flow%fluidsolid_forces = PETSC_FALSE
    nullify(flow%gvt)
  end function FlowCreate

  subroutine FlowDestroy(flow, ierr)
    type(flow_type) flow
    PetscErrorCode ierr
    PetscInt lcv

    if (associated(flow%disc)) call DiscretizationDestroy(flow%disc, ierr)
    if (associated(flow%components)) then
       do lcv=1,flow%ncomponents
          call ComponentDestroy(flow%components(lcv),ierr)
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
    use LBM_Options_module
    type(flow_type) flow
    type(options_type) options
    PetscErrorCode ierr
    PetscInt lcv, lcv2
    PetscBool help
    PetscInt nmax
    PetscBool bcvalue
    PetscScalar gravity(options%ndims)
    PetscScalar,parameter:: eps=1.e-15 ! slightly larger than machine epsilon
    PetscBool flag

    flag = PETSC_FALSE
    flow%ncomponents = options%ncomponents
    flow%ndims = options%ndims
    flow%components => ComponentCreate(flow%comm, flow%ncomponents)

    call DiscretizationSetType(flow%disc, options%flow_disc)
    call DiscretizationSetDerivOrder(flow%disc, options%deriv_order)
    call DiscretizationSetUp(flow%disc)
    
    flow%use_nonideal_eos = options%flow_use_nonideal_eos
    do lcv=1,flow%ncomponents
       call ComponentSetSizes(flow%components(lcv), flow%ncomponents, flow%disc%b)
       if (flow%use_nonideal_eos) then
         flow%components(lcv)%eos => EOSCreate(flow%comm)
       end if
       call ComponentSetID(flow%components(lcv), lcv)
       call ComponentSetFromOptions(flow%components(lcv), options, ierr)
    end do

    ! set up control for forcing terms
    gravity(:) = 0.
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)
    if (help) call PetscPrintf(options%comm, "-gvt=<0,0,0>: gravity\n", ierr)
    nmax = flow%ndims
    call PetscOptionsGetRealArray(options%my_prefix, '-gvt', gravity, &
         nmax, flow%body_forces, ierr)
    if (flow%body_forces) then
      allocate(flow%gvt(flow%ndims))
      flow%gvt = gravity
    end if

    do lcv=1,flow%ncomponents
      do lcv2=1,flow%ncomponents
        if (ABS(flow%components(lcv)%gf(lcv2)) > eps) then
          flow%fluidfluid_forces = PETSC_TRUE
        end if
      end do
    end do

    flow%fluidsolid_forces = options%flow_fluidsolid_forces
    
    ! set up the vectors for holding boundary data
    ! dimension 
    call BCSetSizes(flow%bc, options%ndims*options%ncomponents)
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

    ! check io options
    if (help) call PetscPrintf(options%comm, &
         "  -output_flow_fi <FALSE>: output distribution functions for flow\n", ierr)
    call PetscOptionsGetBool(options%my_prefix, '-output_flow_fi', flow%io_fi, flag, ierr)
    if (help) call PetscPrintf(options%comm, &
         "  -output_flow_last_fi <FALSE>: output distribution functions at the last timestep\n", ierr)
    call PetscOptionsGetBool(options%my_prefix, '-output_flow_last_fi', flow%io_last_fi, flag, ierr)
    if (help) call PetscPrintf(options%comm, "  -output_flow_rho <TRUE>: output densities\n", ierr)
    call PetscOptionsGetBool(options%my_prefix, '-output_flow_rho', flow%io_rho, flag, ierr)
    if (help) call PetscPrintf(options%comm, "  -output_flow_velt <TRUE>: output total velocity\n", ierr)
    call PetscOptionsGetBool(options%my_prefix, '-output_flow_velt', flow%io_velt, flag, ierr)
    if (help) call PetscPrintf(options%comm, "  -output_flow_rhot <FALSE>: output total density\n", ierr)
    call PetscOptionsGetBool(options%my_prefix, '-output_flow_rhot', flow%io_rhot, flag, ierr)
    if (help) call PetscPrintf(options%comm, "  -output_flow_prs <TRUE>: output pressure\n", ierr)
    call PetscOptionsGetBool(options%my_prefix, '-output_flow_prs', flow%io_prs, flag, ierr)
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

    ! deal with consistency.  Assume component 1's tau, mm are correct
    ! viscosity/time scale
    if (flow%components(1)%viscosity < -990.) then
       ! tau is correct, viscosity is not given
       flow%components(1)%time_scale = 1.
       flow%components(1)%viscosity = flow%grid%length_scale * &
         flow%grid%length_scale * (flow%components(1)%relax%tau - 0.5)/3.
    else        
       flow%components(1)%time_scale = flow%grid%length_scale * &
            flow%grid%length_scale * (flow%components(1)%relax%tau - 0.5)/ &
            (flow%components(1)%viscosity*3.)
    end if
    flow%time_scale = flow%components(1)%time_scale
    flow%velocity_scale = flow%grid%length_scale/flow%time_scale

    ! density
    if (flow%components(1)%density < -990.) then
       flow%components(1)%density = 1.
    end if
    flow%mass_scale = flow%components(1)%mm*flow%components(1)%density*(flow%grid%length_scale**3)

    if (flow%ncomponents > 1) then
      do lcv=2,flow%ncomponents
        ! Assert correct viscosity ratios
        if (flow%components(lcv)%viscosity < -990.) then
           flow%components(lcv)%viscosity = flow%grid%length_scale * &
                flow%grid%length_scale * &
                (flow%components(1)%relax%tau - 0.5)/(3.*flow%time_scale)
        else 
          consistent_tau = (3.*flow%time_scale*flow%components(lcv)%viscosity) &
                /(flow%grid%length_scale* &
                flow%grid%length_scale) + 0.5
          if (abs(flow%components(lcv)%relax%tau - 1.) < eps) then
             ! tau not specified, so set it
             flow%components(lcv)%relax%tau = consistent_tau
             if (flow%grid%info%rank .eq. 0) then
                print*, '  Setting component', TRIM(flow%components(lcv)%name), 'tau = ', &
                     consistent_tau
             end if
          else if (abs(flow%components(lcv)%relax%tau - consistent_tau) > eps) then
             SETERRQ(PETSC_COMM_WORLD, 1, 'Viscosities and relaxation times specified '// &
                  'for components are not consistent.', ierr)
          else
             ! all ok
          end if
        end if

        ! Assert correct density ratios
        if (flow%components(lcv)%density < -990.) then
          flow%components(lcv)%density = flow%mass_scale/flow%components(1)%mm / &
               (flow%grid%length_scale**3)
        else
          consistent_mm = flow%mass_scale / flow%components(1)%density / &
               (flow%grid%length_scale**3)
          if (abs(flow%components(lcv)%mm - 1.) < eps) then
             ! mm not specified, so set it
             flow%components(lcv)%mm = consistent_mm
             if (flow%grid%info%rank .eq. 0) then
                print*, '  Setting molecular mass', TRIM(flow%components(lcv)%name), 'mm = ', &
                     consistent_mm
             end if
          else if (abs(flow%components(lcv)%mm - consistent_mm) > eps) then
             SETERRQ(PETSC_COMM_WORLD, 1, 'Densities and molecular masses specified '// &
                  'for components are not consistent.', ierr)
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
    PetscInt lcv
    PetscErrorCode ierr
    PetscScalar zero
    zero = 0.

    do lcv=1,flow%ncomponents
       call DiscretizationSetUpRelax(flow%disc, flow%components(lcv)%relax)
    end do
    call DistributionSetInfo(flow%distribution, flow%grid%info)
    call DistributionSetDiscretization(flow%distribution, flow%disc)
    call DistributionSetSizes(flow%distribution, flow%ncomponents)
    call DistributionSetDAs(flow%distribution, flow%grid%da(NCOMPONENTXBDOF), &
         flow%grid%da(NCOMPONENTDOF))
    call DistributionSetUp(flow%distribution)
    call BCSetUp(flow%bc)

    ! allocate, initialize workspace
    allocate(flow%vel_eq(1:flow%ncomponents, 1:flow%ndims, &
         1:flow%grid%info%gxyzl))
    allocate(flow%forces(1:flow%ncomponents, 1:flow%ndims, &
         1:flow%grid%info%gxyzl))
    allocate(flow%fi_eq(1:flow%ncomponents, 0:flow%disc%b, &
         1:flow%grid%info%gxyzl))
    allocate(flow%psi_of_rho(flow%ncomponents,flow%grid%info%rgxyzl))
    flow%vel_eq = zero
    flow%forces = zero
    flow%fi_eq = zero

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
    type(walls_type) walls

    call DistributionCalcDensity(flow%distribution, walls%walls_a)
    call DistributionCommunicateDensityBegin(flow%distribution)
    call DistributionCalcFlux(flow%distribution, walls%walls_a)
    call DistributionCommunicateDensityEnd(flow%distribution)
    call FlowCalcForces(flow, walls)
    call BCZeroForces(flow%bc, flow%forces, flow%distribution)
    call FlowUpdateU(flow, walls%walls_a)
  end subroutine FlowUpdateMoments

  subroutine FlowUpdateU(flow, walls)
    type(flow_type) flow
    PetscScalar,dimension(1:flow%grid%info%rgxyzl):: walls

    select case(flow%ndims)
    case(2)
       call FlowUpdateUD2(flow, flow%distribution%rho_a, &
            flow%distribution%flux, flow%forces, walls, flow%distribution)
    case(3)
       call FlowUpdateUD3(flow, flow%distribution%rho_a, &
            flow%distribution%flux, flow%forces, walls, flow%distribution)
    end select
  end subroutine FlowUpdateU

  subroutine FlowUpdateUD3(flow, rho, u, forces, walls, dist)
    type(flow_type) flow
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(flow%ncomponents, &
         dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, &
         dist%info%rgzs:dist%info%rgze):: rho
    PetscScalar,dimension(flow%ncomponents, flow%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: u,forces

    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, &
         dist%info%rgzs:dist%info%rgze):: walls

    PetscInt i,j,k,m

    do k=flow%grid%info%zs,flow%grid%info%ze
    do j=flow%grid%info%ys,flow%grid%info%ye
    do i=flow%grid%info%xs,flow%grid%info%xe
      if (walls(i,j,k).eq.0) then
        do m=1,flow%ncomponents
          u(m,:,i,j,k) = (u(m,:,i,j,k) + .5*forces(m,:,i,j,k))/rho(m,i,j,k)
        end do
      end if
    end do
    end do
    end do
  end subroutine FlowUpdateUD3

  subroutine FlowUpdateUD2(flow, rho, u, forces, walls, dist)
    type(flow_type) flow
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(1:dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u,forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls

    PetscInt i,j,m

    do j=flow%grid%info%ys,flow%grid%info%ye
    do i=flow%grid%info%xs,flow%grid%info%xe
      if (walls(i,j).eq.0) then
        do m=1,flow%ncomponents
          u(m,:,i,j) = (u(m,:,i,j) + .5*forces(m,:,i,j))/rho(m,i,j)
        end do
      end if
    end do
    end do
  end subroutine FlowUpdateUD2

  subroutine FlowOutputDiagnostics(flow, io)
    use LBM_IO_module
    type(flow_type) flow
    type(io_type) io
    PetscErrorCode ierr

    if (flow%io_fi) then
       call DMDAVecRestoreArrayF90(flow%grid%da(NCOMPONENTXBDOF),flow%distribution%fi, &
            flow%distribution%fi_a, ierr)
       call DMLocalToGlobalBegin(flow%grid%da(NCOMPONENTXBDOF),flow%distribution%fi, &
            INSERT_VALUES, flow%distribution%fi_g, ierr)
       call DMLocalToGlobalEnd(flow%grid%da(NCOMPONENTXBDOF),flow%distribution%fi, &
            INSERT_VALUES, flow%distribution%fi_g, ierr)
       call DMDAVecGetArrayF90(flow%grid%da(NCOMPONENTXBDOF),flow%distribution%fi, &
            flow%distribution%fi_a, ierr)
       call IOView(io, flow%distribution%fi_g, 'fi')
    end if
    if (flow%io_rho) then
       call DMDAVecRestoreArrayF90(flow%grid%da(NCOMPONENTDOF),flow%distribution%rho, &
            flow%distribution%rho_a, ierr)
       call DMLocalToGlobalBegin(flow%grid%da(NCOMPONENTDOF),flow%distribution%rho, &
            INSERT_VALUES, flow%distribution%rho_g, ierr)
       call DMLocalToGlobalEnd(flow%grid%da(NCOMPONENTDOF),flow%distribution%rho, &
            INSERT_VALUES, flow%distribution%rho_g, ierr)
       call VecScale(flow%distribution%rho_g, &
            flow%mass_scale/(flow%grid%length_scale**3), ierr)
       call DMDAVecGetArrayF90(flow%grid%da(NCOMPONENTDOF),flow%distribution%rho, &
            flow%distribution%rho_a, ierr)
       call IOView(io, flow%distribution%rho_g, 'rho')
       call VecScale(flow%distribution%rho_g, &
            (flow%grid%length_scale**3)/flow%mass_scale, ierr)
    end if

    call VecScale(flow%velt_g, flow%velocity_scale, ierr)
    if (flow%io_velt) call IOView(io, flow%velt_g, 'u')
    call VecScale(flow%velt_g, 1./flow%velocity_scale, ierr)

    if (flow%io_rhot) call IOView(io, flow%rhot_g, 'rhot')
    if (flow%io_prs) call IOView(io, flow%prs_g, 'prs')
  end subroutine FlowOutputDiagnostics

  subroutine FlowUpdateDiagnostics(flow, walls)
    type(flow_type) flow
    type(walls_type) walls
    PetscErrorCode ierr
    
    PetscInt m

    call DMDAVecGetArrayF90(flow%grid%da(NFLOWDOF),flow%velt_g,flow%velt_a,ierr)
    call DMDAVecGetArrayF90(flow%grid%da(ONEDOF), flow%rhot_g, flow%rhot_a,ierr)
    call DMDAVecGetArrayF90(flow%grid%da(ONEDOF), flow%prs_g, flow%prs_a, ierr)

    if (flow%use_nonideal_eos) then
      do m=1,flow%distribution%s
        call EOSApply(flow%components(m)%eos, flow%distribution%rho_a, &
             flow%psi_of_rho, flow%components(m)%gf(m), m, flow%distribution)
      end do
    end if

    select case(flow%ndims)
    case(2)
      if (flow%use_nonideal_eos) then
        call FlowUpdateDiagnosticsD2(flow, flow%distribution%rho_a, flow%psi_of_rho, &
             flow%distribution%flux, flow%forces, walls%walls_a, flow%rhot_a, &
             flow%prs_a, flow%velt_a)
      else
        call FlowUpdateDiagnosticsD2(flow, flow%distribution%rho_a, flow%distribution%rho_a, &
             flow%distribution%flux, flow%forces, walls%walls_a, flow%rhot_a, &
             flow%prs_a, flow%velt_a)
      end if
    case(3)
      if (flow%use_nonideal_eos) then
        call FlowUpdateDiagnosticsD3(flow, flow%distribution%rho_a, flow%psi_of_rho, &
             flow%distribution%flux, flow%forces, walls%walls_a, flow%rhot_a, &
             flow%prs_a, flow%velt_a)
      else
        call FlowUpdateDiagnosticsD3(flow, flow%distribution%rho_a, flow%distribution%rho_a, &
             flow%distribution%flux, flow%forces, walls%walls_a, flow%rhot_a, &
             flow%prs_a, flow%velt_a)
      end if
    end select
    call DMDAVecRestoreArrayF90(flow%grid%da(NFLOWDOF),flow%velt_g,flow%velt_a,ierr)
    call DMDAVecRestoreArrayF90(flow%grid%da(ONEDOF), flow%rhot_g, flow%rhot_a, ierr)
    call DMDAVecRestoreArrayF90(flow%grid%da(ONEDOF), flow%prs_g, flow%prs_a, ierr)

  end subroutine FlowUpdateDiagnostics

  subroutine FlowUpdateDiagnosticsD3(flow, rho, psi, u, forces, walls, rhot, prs, velt)
    type(flow_type) flow
    PetscScalar,dimension(flow%grid%info%rgxs:flow%grid%info%rgxe, &
         flow%grid%info%rgys:flow%grid%info%rgye, &
         flow%grid%info%rgzs:flow%grid%info%rgze):: walls
    PetscScalar,dimension(flow%ncomponents, &
         flow%grid%info%rgxs:flow%grid%info%rgxe, &
         flow%grid%info%rgys:flow%grid%info%rgye, &
         flow%grid%info%rgzs:flow%grid%info%rgze):: rho, psi
    PetscScalar,dimension(flow%ncomponents, flow%ndims, &
         flow%grid%info%gxs:flow%grid%info%gxe, &
         flow%grid%info%gys:flow%grid%info%gye, &
         flow%grid%info%gzs:flow%grid%info%gze):: u,forces
    PetscScalar,dimension(flow%grid%info%xs:flow%grid%info%xe, &
         flow%grid%info%ys:flow%grid%info%ye, &
         flow%grid%info%zs:flow%grid%info%ze):: rhot,prs
    PetscScalar,dimension(flow%ndims, &
         flow%grid%info%xs:flow%grid%info%xe, &
         flow%grid%info%ys:flow%grid%info%ye, &
         flow%grid%info%zs:flow%grid%info%ze):: velt

    PetscInt i,j,k,m,d
    PetscScalar mm(1:flow%ncomponents)
    
    do m=1,flow%ncomponents
       mm(m) = flow%components(m)%mm
    end do

    do k=flow%grid%info%zs,flow%grid%info%ze
    do j=flow%grid%info%ys,flow%grid%info%ye
    do i=flow%grid%info%xs,flow%grid%info%xe
    if (walls(i,j,k).eq.0) then
       rhot(i,j,k) = sum(rho(:,i,j,k)*mm,1)
       do d=1,flow%ndims
          velt(d,i,j,k) = sum(u(:,d,i,j,k)*mm(:)*rho(:,i,j,k))/rhot(i,j,k)
       end do
       prs(i,j,k) = rhot(i,j,k)/3.
       if (flow%use_nonideal_eos .or. (flow%ncomponents > 1)) then
          do m=1,flow%ncomponents
             prs(i,j,k) = prs(i,j,k) + flow%disc%c_0/2.*psi(m,i,j,k) &
                  *sum(flow%components(m)%gf*psi(:,i,j,k),1)
          end do
       end if
    end if
    end do
    end do
    end do
  end subroutine FlowUpdateDiagnosticsD3

  subroutine FlowUpdateDiagnosticsD2(flow, rho, psi, u, forces, walls, rhot, prs, velt)
    type(flow_type) flow
    PetscScalar,dimension(flow%grid%info%rgxs:flow%grid%info%rgxe, &
         flow%grid%info%rgys:flow%grid%info%rgye):: walls
    PetscScalar,dimension(flow%ncomponents, &
         flow%grid%info%rgxs:flow%grid%info%rgxe, &
         flow%grid%info%rgys:flow%grid%info%rgye):: rho, psi
    PetscScalar,dimension(flow%ncomponents, flow%ndims, &
         flow%grid%info%gxs:flow%grid%info%gxe, &
         flow%grid%info%gys:flow%grid%info%gye):: u,forces
    PetscScalar,dimension(flow%grid%info%xs:flow%grid%info%xe, &
         flow%grid%info%ys:flow%grid%info%ye):: rhot,prs
    PetscScalar,dimension(flow%ndims, &
         flow%grid%info%xs:flow%grid%info%xe, &
         flow%grid%info%ys:flow%grid%info%ye):: velt

    PetscInt i,j,m,d
    PetscScalar mm(1:flow%ncomponents)
    
    do m=1,flow%ncomponents
       mm(m) = flow%components(m)%mm
    end do

    do j=flow%grid%info%ys,flow%grid%info%ye
    do i=flow%grid%info%xs,flow%grid%info%xe
    if (walls(i,j).eq.0) then
       rhot(i,j) = sum(rho(:,i,j)*mm,1)
       do d=1,flow%ndims
          velt(d,i,j) = sum(u(:,d,i,j)*mm(:)*rho(:,i,j))/rhot(i,j)
       end do
       prs(i,j) = rhot(i,j)/3.
       if (flow%use_nonideal_eos .or. (flow%ncomponents > 1)) then
          do m=1,flow%ncomponents
             prs(i,j) = prs(i,j) + flow%disc%c_0/2.*psi(m,i,j) &
                  *sum(flow%components(m)%gf*psi(:,i,j),1)
          end do
       end if
    end if
    end do
    end do
  end subroutine FlowUpdateDiagnosticsD2

  subroutine FlowCalcForces(flow, walls)
    use LBM_Forcing_module
    use LBM_Logging_module
    type(flow_type) flow
    type(walls_type) walls

    PetscErrorCode ierr
    PetscInt m

    flow%forces = 0.
    if (flow%fluidfluid_forces) then
      call PetscLogEventBegin(logger%event_forcing_fluidfluid,ierr)
      if (flow%use_nonideal_eos) then
        do m=1,flow%distribution%s
          call EOSApply(flow%components(m)%eos, flow%distribution%rho_a, &
               flow%psi_of_rho, flow%components(m)%gf(m), m, flow%distribution)
        end do
        call LBMAddFluidFluidForces(flow%distribution, flow%components, &
             flow%psi_of_rho, walls%walls_a, flow%forces)
      else
        call LBMAddFluidFluidForces(flow%distribution, flow%components, &
             flow%distribution%rho_a, walls%walls_a, flow%forces)
      end if
      call PetscLogEventEnd(logger%event_forcing_fluidfluid,ierr)
    end if

    if (flow%fluidsolid_forces) then
      call PetscLogEventBegin(logger%event_forcing_fluidsolid,ierr)
      call LBMAddFluidSolidForces(flow%distribution, flow%components, walls, &
           flow%distribution%rho_a, flow%forces)
      call PetscLogEventEnd(logger%event_forcing_fluidsolid,ierr)
    end if

    if (flow%body_forces) then
      call PetscLogEventBegin(logger%event_forcing_body,ierr)
      call LBMAddBodyForces(flow%distribution, flow%components, flow%gvt, &
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
    type(walls_type) walls
    call DistributionBounceback(flow%distribution, walls%walls_a)
  end subroutine FlowBounceback

  subroutine FlowUpdateFeq(flow, walls)
    type(flow_type) flow
    PetscScalar,dimension(flow%grid%info%rgxyzl):: walls

    PetscInt m
    PetscErrorCode ierr

    do m=1,flow%distribution%s
      call DiscretizationEquilf(flow%disc, flow%distribution%rho_a, &
           flow%distribution%flux, walls, flow%fi_eq, m, flow%components(m)%relax, &
           flow%distribution)
    end do
  end subroutine FlowUpdateFeq

  subroutine FlowFiBarEqPrefactor(flow, rho, u, forces, prefactor, dist)
    type(flow_type) flow
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(1:dist%s):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims):: u,forces
    PetscScalar,dimension(1:dist%s,0:flow%disc%b):: prefactor

    PetscInt m,n

    do n=0,flow%disc%b
      do m=1,dist%s
        prefactor(m,n) = sum(forces(m,:)*(flow%disc%ci(n,:)-u(m,:))) / &
             (rho(m)*flow%components(m)%relax%c_s2)
      end do
    end do
  end subroutine FlowFiBarEqPrefactor


  subroutine FlowFiBarInit(flow, walls)
    type(flow_type) flow
    PetscScalar,dimension(flow%grid%info%rgxyzl):: walls

    select case(flow%ndims)
    case(2)
       call FlowFeqBarD2(flow, flow%distribution%fi_a, flow%distribution%rho_a, &
            flow%distribution%flux, flow%forces, walls, flow%fi_eq, flow%distribution)
    case(3)
       call FlowFeqBarD3(flow, flow%distribution%fi_a, flow%distribution%rho_a, &
            flow%distribution%flux, flow%forces, walls, flow%fi_eq, flow%distribution)
     end select
   end subroutine FlowFiBarInit

  subroutine FlowFeqBarD3(flow, fi_eq_bar, rho, u, forces, walls, fi_eq, dist)
    type(flow_type) flow
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(flow%ncomponents, 0:flow%disc%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze)::fi_eq_bar, fi_eq
    PetscScalar,dimension(1:dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: u,forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: walls

    PetscScalar,dimension(1:dist%s,0:flow%disc%b):: prefactor
    PetscInt i,j,k,m,n

    ! \bar{f^eq} = (1 + dt/2 * (F*(e_i - u))/(rho*RT))*f^eq
    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j,k).eq.0) then
      call FlowFiBarEqPrefactor(flow, rho(:,i,j,k), u(:,:,i,j,k), forces(:,:,i,j,k), &
           prefactor, dist)
      fi_eq_bar(:,:,i,j,k) = (1. - 0.5*prefactor)*fi_eq(:,:,i,j,k)
    end if
    end do
    end do
    end do
  end subroutine FlowFeqBarD3

  subroutine FlowFeqBarD2(flow, fi_eq_bar, rho, u, forces, walls, fi_eq, dist)
    type(flow_type) flow
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(flow%ncomponents, 0:flow%disc%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye)::fi_eq_bar, fi_eq
    PetscScalar,dimension(1:dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u,forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls

    PetscScalar,dimension(1:dist%s,0:flow%disc%b):: prefactor
    PetscInt i,j,m,n

    ! \bar{f^eq} = (1 + dt/2 * (F*(e_i - u))/(rho*RT))*f^eq
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j).eq.0) then
      call FlowFiBarEqPrefactor(flow, rho(:,i,j), u(:,:,i,j), forces(:,:,i,j), &
           prefactor, dist)
      fi_eq_bar(:,:,i,j) = (1. - 0.5*prefactor)*fi_eq(:,:,i,j)
    end if
    end do
    end do
  end subroutine FlowFeqBarD2

  subroutine FlowFiInit(flow, walls)
    type(flow_type) flow
    type(walls_type) walls
    
    ! Note that the initialization process is calculated without walls.
    ! This is because we must fake that a bounceback has occured in a
    ! previous step.  Wall nodes, which should hold the to-be-bounced 
    ! back fi, will instead have zeros.  
    PetscScalar,dimension(1:flow%grid%info%rgxyzl):: tmp_no_walls
    tmp_no_walls = 0.

    call DistributionCommunicateDensity(flow%distribution)
    call FlowCalcForces(flow, walls) ! not sure if this is ok or not!  
                                     ! Does this need no-walls as well?
    call FlowUpdateFeq(flow, tmp_no_walls)
    call FlowFiBarInit(flow, tmp_no_walls)
  end subroutine FlowFiInit

  subroutine FlowCollision(flow, walls)
    use LBM_Logging_module
    type(flow_type) flow
    type(walls_type) walls
    PetscErrorCode ierr

    call PetscLogEventBegin(logger%event_collision_feq,ierr)
    call FlowUpdateFeq(flow, walls%walls_a)
    call PetscLogEventEnd(logger%event_collision_feq,ierr)

    call PetscLogEventBegin(logger%event_collision_relax,ierr)
    select case(flow%ndims)
    case(2)
       call FlowCollisionD2(flow, flow%distribution%fi_a, flow%distribution%rho_a, &
            flow%distribution%flux, flow%forces, walls%walls_a, flow%fi_eq, &
            flow%distribution)
    case(3)
       call FlowCollisionD3(flow, flow%distribution%fi_a, flow%distribution%rho_a, &
            flow%distribution%flux, flow%forces, walls%walls_a, flow%fi_eq, &
            flow%distribution)
    end select
    call PetscLogEventEnd(logger%event_collision_relax,ierr)
  end subroutine FlowCollision

  subroutine FlowCollisionD3(flow, fi, rho, u, forces, walls, fi_eq, dist)
    type(flow_type) flow
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(flow%ncomponents, 0:flow%disc%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze)::fi,fi_eq
    PetscScalar,dimension(1:dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: u,forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: walls

    PetscInt m,i,j,k,d
    PetscScalar,dimension(1:dist%s,0:flow%disc%b):: prefactor
    PetscScalar,dimension(1:dist%s,0:flow%disc%b):: fi_eq_bar
    PetscErrorCode ierr

    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j,k).eq.0) then
      call FlowFiBarEqPrefactor(flow, rho(:,i,j,k), u(:,:,i,j,k), forces(:,:,i,j,k), &
           prefactor, dist)
      fi_eq_bar = (1. - .5*prefactor(:,:))*fi_eq(:,:,i,j,k)

      do m=1,dist%s
        call RelaxationCollide(flow%components(m)%relax, fi(:,:,i,j,k), fi_eq_bar, m, dist)
      end do
      fi(:,:,i,j,k) = fi(:,:,i,j,k) + prefactor(:,:)*fi_eq(:,:,i,j,k)
    end if
    end do
    end do
    end do
  end subroutine FlowCollisionD3

  subroutine FlowCollisionD2(flow, fi, rho, u, forces, walls, fi_eq, dist)
    type(flow_type) flow
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(flow%ncomponents, 0:flow%disc%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye)::fi,fi_eq
    PetscScalar,dimension(1:dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u,forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls

    PetscScalar,dimension(1:dist%s,0:flow%disc%b):: prefactor
    PetscScalar,dimension(1:dist%s,0:flow%disc%b):: fi_eq_bar

    PetscInt m,i,j,d
    PetscScalar,parameter:: eps=1.e-15 ! slightly larger than machine epsilon
    PetscScalar,dimension(1:dist%s) :: mm
    PetscErrorCode ierr

    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j).eq.0) then
      call FlowFiBarEqPrefactor(flow, rho(:,i,j), u(:,:,i,j), forces(:,:,i,j), &
           prefactor, dist)
      fi_eq_bar = (1. - .5*prefactor(:,:))*fi_eq(:,:,i,j)

      do m=1,dist%s
        call RelaxationCollide(flow%components(m)%relax, fi(:,:,i,j), fi_eq_bar, m, dist)
      end do
      fi(:,:,i,j) = fi(:,:,i,j) + prefactor(:,:)*fi_eq(:,:,i,j)
    end if
    end do
    end do
  end subroutine FlowCollisionD2

  subroutine FlowApplyBCs(flow, walls)
    type(flow_type) flow
    type(walls_type) walls

    call BCApply(flow%bc, walls%walls_a, flow%distribution)
  end subroutine FlowApplyBCs

  subroutine print_a_few(fi, rho, u, forces, walls, dist, istep)
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s,0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u,forces
    PetscScalar,dimension(1:dist%s, dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls

    PetscInt i,j,istep

    print*, 'step:', istep
    print*, '---------------------'
    i=37
    j=37
    print*, 'walls:', walls(i,j)
    print*, 'outer:'
    print*, '  fi(1,:):', fi(1,:,i,j)
    print*, '  rho(1,x):', rho(1,i,j)
    print*, '  u(1,x):', u(1,:,i,j)
    print*, '  forces(1,x):', forces(1,:,i,j)
    print*, 'inner:'
    print*, '  fi(2,:):', fi(2,:,i,j)
    print*, '  rho(2,x):', rho(2,i,j)
    print*, '  u(1,x):', u(2,:,i,j)
    print*, '  forces(1,x):', forces(2,:,i,j)
    print*, '---------------------'

    ! i=1
    ! j=100
    ! print*, 'inner:'
    ! print*, 'walls:', walls(i,j)
    ! print*, 'outer:'
    ! print*, 'fi(1,:):', fi(1,:,i,j)
    ! print*, 'rho(1,x):', rho(:,i,j)
    ! print*, 'u(1,x):', u(1,:,i,j)
    ! print*, 'forces(1,x):', forces(1,:,i,j)
    print*, '========================'
  end subroutine print_a_few
end module LBM_Flow_module
