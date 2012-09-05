!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_flow.F90
!!!     version:         
!!!     created:         17 March 2011
!!!       on:            17:58:06 MDT
!!!     last modified:   17 November 2011
!!!       at:            12:03:10 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"

module LBM_Flow_module
  use petsc
  use LBM_Error_module
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
     PetscInt, dimension(6):: bc_flags
     PetscBool, dimension(15):: bc_done
     PetscScalar,pointer,dimension(:,:) :: bc_data

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
     PetscScalar null_pressure ! the value assigned to wall cells for pressure
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
       FlowUpdateFlux, &
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
    flow%name = ''
    flow%ncomponents = -1
    flow%ndims = -1
    nullify(flow%components)
    nullify(flow%grid)
    flow%disc => DiscretizationCreate(flow%comm)
    flow%distribution => DistributionCreate(flow%comm)
    flow%bc => BCCreate(flow%comm)

    flow%bc_flags = BC_NULL
    flow%bc_done = PETSC_FALSE
    nullify(flow%bc_data)

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

    flow%velocity_scale = 1.
    flow%time_scale = 1.
    flow%mass_scale = 1.
    flow%null_pressure = 0.
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
    if (associated(flow%bc_data)) deallocate(flow%bc_data)
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
    character(len=3) boundaryname

    flag = PETSC_FALSE
    flow%ncomponents = options%ncomponents
    flow%ndims = options%ndims
    flow%components => ComponentCreate(flow%comm, flow%ncomponents)

    call DiscretizationSetType(flow%disc, options%flow_disc)
    call DiscretizationSetDerivOrder(flow%disc, options%isotropy_order)
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
    if (help) call PetscPrintf(flow%comm, "-gvt=<0,0,0>: gravity\n", ierr)
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
    allocate(flow%bc_data(flow%bc%nbcs,2*flow%ndims))
    flow%bc_data(:,:) = 0.d0

    ! parse flow boundary conditions
    boundaryname = "xm"
    call FlowParseBC(flow, options, BOUNDARY_XM, boundaryname)
    boundaryname = "xp"
    call FlowParseBC(flow, options, BOUNDARY_XP, boundaryname)
    boundaryname = "ym"
    call FlowParseBC(flow, options, BOUNDARY_YM, boundaryname)
    boundaryname = "yp"
    call FlowParseBC(flow, options, BOUNDARY_YP, boundaryname)
    if (flow%ndims > 2) then
      boundaryname = "zm"
      call FlowParseBC(flow, options, BOUNDARY_ZM, boundaryname)
      boundaryname = "zp"
      call FlowParseBC(flow, options, BOUNDARY_ZP, boundaryname)
    endif

    ! check io options
    call OptionsGetBool(options, "-output_flow_fi", &
         "ouput distribution function for flow", flow%io_fi, flag, ierr)
    call OptionsGetBool(options, "-output_flow_last_fi", &
         "ouput distribution function for flow at the last timestep", &
         flow%io_last_fi, flag, ierr)
    call OptionsGetBool(options, "-output_rho", "output density", &
         flow%io_rho, flag, ierr)
    call OptionsGetBool(options, "-output_velocity", "output total velocity", &
         flow%io_velt, flag, ierr)
    call OptionsGetBool(options, "-output_rhot", "output total density", &
         flow%io_rhot, flag, ierr)
    call OptionsGetBool(options, "-output_pressure", "output pressure", &
         flow%io_prs, flag, ierr)

    ! steady state checking
    call DistributionSetTrackOld(flow%distribution, options%steadystate_field_rho, &
         options%steadystate_field_fi)
  end subroutine FlowSetFromOptions

  subroutine FlowSetGrid(flow, grid)
    use LBM_Grid_module
    type(flow_type) flow
    type (grid_type),pointer:: grid
    flow%grid => grid
    call BCSetGrid(flow%bc, flow%grid)
  end subroutine FlowSetGrid

  ! This currently only makes sense for SRT -- should be updated at some point.
  subroutine FlowSetPhysicalScales(flow, ierr)
    type(flow_type) flow
    PetscScalar consistent_tau, consistent_mm
    PetscErrorCode ierr
    PetscScalar,parameter:: eps=1.e-8
    PetscInt lcv

    ! deal with consistency.  Assume component 1's tau, mm are correct
    ! viscosity/time scale
    if (flow%components(1)%viscosity < -990.d0) then
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

          else if (abs(flow%components(lcv)%relax%tau - consistent_tau) > eps) then
             call LBMError(flow%comm, 1, &
                  'Viscosities and relaxation times specified '// &
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

          else if (abs(flow%components(lcv)%mm - consistent_mm) > eps) then
             call LBMError(flow%comm, 1, &
                  'Densities and molecular masses specified '// &
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
       print*, '  velocity scale [m/s]:', flow%velocity_scale
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

    ! set up boundary conditions
    call BCSetUp(flow%bc)
    select case(flow%ndims)
    case(2)
      call FlowSetUpBCsD2(flow, flow%bc%xm_a, flow%bc%xp_a, &
           flow%bc%ym_a, flow%bc%yp_a, flow%distribution)
    case(3)
      call FlowSetUpBCsD3(flow, flow%bc%xm_a, flow%bc%xp_a, flow%bc%ym_a, &
           flow%bc%yp_a, flow%bc%zm_a, flow%bc%zp_a, flow%distribution)
    end select

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

  subroutine FlowCalcRhoForces(flow, walls)
    type(flow_type) flow
    type(walls_type) walls

    call DistributionCalcDensity(flow%distribution, walls%walls_a)
    call BCApplyDirichletToRho(flow%bc, walls%walls_a, flow%distribution)
    call DistributionCommunicateDensity(flow%distribution)
    call FlowCalcForces(flow, walls)
  end subroutine FlowCalcRhoForces

  subroutine FlowUpdateFlux(flow, walls)
    type(flow_type) flow
    type(walls_type) walls

    call DistributionCalcFlux(flow%distribution, walls%walls_a)
    call FlowUpdateUE(flow, walls%walls_a)
  end subroutine FlowUpdateFlux

  subroutine FlowUpdateMoments(flow, walls)
    type(flow_type) flow
    type(walls_type) walls

    call DistributionCalcDensity(flow%distribution, walls%walls_a)
    call DistributionCommunicateDensityBegin(flow%distribution)
    call DistributionCalcFlux(flow%distribution, walls%walls_a)
    call DistributionCommunicateDensityEnd(flow%distribution)
    call FlowCalcForces(flow, walls)
    call FlowUpdateUE(flow, walls%walls_a)
  end subroutine FlowUpdateMoments

  subroutine FlowUpdateUE(flow, walls)
    type(flow_type) flow
    PetscScalar,dimension(1:flow%grid%info%rgxyzl):: walls

    select case(flow%ndims)
    case(2)
       call FlowUpdateUED2(flow, flow%distribution%rho_a, &
            flow%distribution%flux, flow%forces, walls, flow%distribution)
    case(3)
       call FlowUpdateUED3(flow, flow%distribution%rho_a, &
            flow%distribution%flux, flow%forces, walls, flow%distribution)
    end select
  end subroutine FlowUpdateUE

  subroutine FlowUpdateUED3(flow, rho, ue, forces, walls, dist)
    type(flow_type) flow
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(flow%ncomponents, &
         dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, &
         dist%info%rgzs:dist%info%rgze):: rho
    PetscScalar,dimension(flow%ncomponents, flow%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: ue,forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, &
         dist%info%rgzs:dist%info%rgze):: walls

    PetscScalar,dimension(flow%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: up
    PetscInt i,j,k,m,d
    PetscScalar,dimension(1:dist%s) :: mmot

    do m=1,dist%s
       mmot(m) = flow%components(m)%mm*flow%components(m)%relax%s_c
    end do

    do k=flow%grid%info%zs,flow%grid%info%ze
    do j=flow%grid%info%ys,flow%grid%info%ye
    do i=flow%grid%info%xs,flow%grid%info%xe
      if (walls(i,j,k).eq.0) then
        ue(:,:,i,j,k) = ue(:,:,i,j,k) + .5*forces(:,:,i,j,k)

        do d=1,dist%info%ndims
          up(d,i,j,k) = sum(ue(:,d,i,j,k)*mmot,1)/sum(rho(:,i,j,k)*mmot,1)
        end do

        do m=1,flow%ncomponents
          ue(m,:,i,j,k) = up(:,i,j,k)
        end do
      end if
    end do
    end do
    end do
  end subroutine FlowUpdateUED3

  subroutine FlowUpdateUED2(flow, rho, ue, forces, walls, dist)
    type(flow_type) flow
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(1:dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: ue,forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls

    PetscScalar,dimension(flow%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: up
    PetscInt i,j,m,d
    PetscScalar,dimension(1:dist%s) :: mmot

    do m=1,dist%s
       mmot(m) = flow%components(m)%mm*flow%components(m)%relax%s_c
    end do

    do j=flow%grid%info%ys,flow%grid%info%ye
    do i=flow%grid%info%xs,flow%grid%info%xe
      if (walls(i,j).eq.0) then
        ue(:,:,i,j) = ue(:,:,i,j) + .5*forces(:,:,i,j)

        do d=1,dist%info%ndims
          up(d,i,j) = sum(ue(:,d,i,j)*mmot,1)/sum(rho(:,i,j)*mmot,1)
        end do

        do m=1,flow%ncomponents
          ue(m,:,i,j) = up(:,i,j)
        end do
      end if
    end do
    end do
  end subroutine FlowUpdateUED2

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

    PetscScalar,dimension(flow%ncomponents, flow%ndims, &
         flow%grid%info%gxyzl):: u
    PetscInt m

    call DistributionCalcFlux(flow%distribution, walls%walls_a, u)

    call DMDAVecGetArrayF90(flow%grid%da(NFLOWDOF),flow%velt_g,flow%velt_a,ierr)
    call DMDAVecGetArrayF90(flow%grid%da(ONEDOF), flow%rhot_g, flow%rhot_a,ierr)
    call DMDAVecGetArrayF90(flow%grid%da(ONEDOF), flow%prs_g, flow%prs_a, ierr)

    if (flow%use_nonideal_eos) then
      do m=1,flow%ncomponents
        call EOSApply(flow%components(m)%eos, flow%distribution%rho_a, &
             flow%psi_of_rho, flow%components(m)%gf(m), m, flow%distribution)
      end do
    end if


    select case(flow%ndims)
    case(2)
      if (flow%use_nonideal_eos) then
        call FlowUpdateDiagnosticsD2(flow, flow%distribution%rho_a, flow%psi_of_rho, &
             u, flow%forces, walls%walls_a, flow%rhot_a, &
             flow%prs_a, flow%velt_a)
      else
        call FlowUpdateDiagnosticsD2(flow, flow%distribution%rho_a, flow%distribution%rho_a, &
             u, flow%forces, walls%walls_a, flow%rhot_a, &
             flow%prs_a, flow%velt_a)
      end if
    case(3)
      if (flow%use_nonideal_eos) then
        call FlowUpdateDiagnosticsD3(flow, flow%distribution%rho_a, flow%psi_of_rho, &
             u, flow%forces, walls%walls_a, flow%rhot_a, &
             flow%prs_a, flow%velt_a)
      else
        call FlowUpdateDiagnosticsD3(flow, flow%distribution%rho_a, flow%distribution%rho_a, &
             u, flow%forces, walls%walls_a, flow%rhot_a, &
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
       prs(i,j,k) = rhot(i,j,k)/3.
       if (flow%use_nonideal_eos .or. (flow%ncomponents > 1)) then
          do m=1,flow%ncomponents
             prs(i,j,k) = prs(i,j,k) + flow%disc%c_0/2.*psi(m,i,j,k) &
                  *sum(flow%components(m)%gf*psi(:,i,j,k),1)
          end do
       end if

       do m=1,flow%ncomponents
         u(m,:,i,j,k) = u(m,:,i,j,k) + .5*forces(m,:,i,j,k)
       end do

       do d=1,flow%ndims
         velt(d,i,j,k) = sum(u(:,d,i,j,k)*mm(:))/rhot(i,j,k)
       end do
    else
      prs(i,j,k) = flow%null_pressure
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
       prs(i,j) = rhot(i,j)/3.
       if (flow%use_nonideal_eos .or. (flow%ncomponents > 1)) then
          do m=1,flow%ncomponents
             prs(i,j) = prs(i,j) + flow%disc%c_0/2.*psi(m,i,j) &
                  *sum(flow%components(m)%gf*psi(:,i,j),1)
          end do
       end if

       do m=1,flow%ncomponents
         u(m,:,i,j) = u(m,:,i,j) + .5*forces(m,:,i,j)
       end do

       do d=1,flow%ndims
          velt(d,i,j) = sum(u(:,d,i,j)*mm(:))/rhot(i,j)
       end do
    else
      prs(i,j) = flow%null_pressure
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
        do m=1,flow%ncomponents
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
    call BCPreStream(flow%bc,flow%distribution)
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

    do m=1,flow%ncomponents
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

    call DistributionCommunicateDensity(flow%distribution)
    call FlowCalcForces(flow, walls)
    call FlowUpdateFeq(flow, walls%walls_a)
    call FlowFiBarInit(flow, walls%walls_a)
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

  subroutine FlowParseBC(flow, options, boundary, bname)
    use LBM_Options_module

    type(flow_type) flow
    type(options_type) options
    PetscInt boundary
    character(len=3):: bname

    PetscBool done, bcvalue, flag
    PetscInt bccount, m
    PetscErrorCode ierr

    done = PETSC_FALSE

    ! periodic boundary condition
    bcvalue = PETSC_FALSE
    call OptionsGetBool(options, "-bc_periodic_"//bname(1:1), "periodic in x", &
         bcvalue, flag, ierr)
    if (bcvalue) then
      done = PETSC_TRUE
      flow%bc%flags(BOUNDARY) = BC_PERIODIC
      flow%bc_flags(BOUNDARY) = BC_PERIODIC
    end if

    ! reflecting boundary condition... much like a "no stick" BC
    bcvalue = PETSC_FALSE
    call OptionsGetBool(options, "-bc_reflecting_"//trim(bname), "reflecting bc", &
         bcvalue, flag, ierr)
    if (bcvalue) then
      done = PETSC_TRUE
      flow%bc%flags(BOUNDARY) = BC_REFLECTING
      flow%bc_flags(BOUNDARY) = BC_REFLECTING
    end if

    ! density boundary condition
    bcvalue = PETSC_FALSE
    call OptionsGetBool(options, "-bc_density_"//trim(bname),"provide density", &
         bcvalue, flag, ierr)
    if (bcvalue) then
      if (done) call LBMError(flow%comm, 1, &
           "Multiple BCs provided for boundary "//trim(bname), ierr)
      done = PETSC_TRUE
      flow%bc%flags(boundary) = BC_DIRICHLET
      flow%bc_flags(boundary) = BC_DENSITY
      flow%bc_data(:,boundary) = -1.d0

      if (flow%ncomponents.eq.1) then
        call OptionsGetReal(options, "-bc_density_"//trim(bname)//"_value", &
             "density", flow%bc_data(1,boundary), flag, ierr)
        if (.not.flag) call LBMWarn(flow%comm, "Density BC on boundary " &
             //trim(bname)//" not set in input file!", ierr)
      else
        do m=1,flow%ncomponents
          call OptionsGetReal(options, "-bc_density_"//trim(bname)//"_"// &
               trim(flow%components(m)%name), "density", flow%bc_data(m,boundary), &
               flag, ierr)
          if (.not.flag) then
            call LBMWarn(flow%comm, "Density BC on boundary "// &
                 trim(bname)//" not set in input file!", ierr)
          endif
        end do
      end if
    endif

    ! single phase secondary boundary conditions
    if (flow%ncomponents.eq.1) then
      bcvalue = PETSC_FALSE
      call OptionsGetBool(options,"-bc_pressure_"//trim(bname), &
           "provide pressure", bcvalue, flag, ierr)
      if (bcvalue) then
        if (done) call LBMError(flow%comm, 1, &
             "Multiple BCs provided for boundary "//trim(bname), ierr)
        done = PETSC_TRUE
        flow%bc%flags(boundary) = BC_DIRICHLET
        flow%bc_flags(boundary) = BC_PRESSURE
        flow%bc_data(:,boundary) = -1.d0

        call OptionsGetReal(options, "-bc_pressure_"//trim(bname)//"_value", &
             "pressure", flow%bc_data(1,boundary), flag, ierr)
        if (.not.flag) then
          call LBMWarn(flow%comm, "Pressure BC on boundary "//trim(bname)//&
               " not set in input file!", ierr)
        end if
      end if

      ! flux inlet provides total flux and rho1/rho_total
      bcvalue = PETSC_FALSE
      call OptionsGetBool(options,"-bc_flux_"//trim(bname), &
           "provide total flux", bcvalue, flag, ierr)
      if (bcvalue) then
        if (done) call LBMError(flow%comm, 1, &
             "Multiple BCs provided for boundary "//trim(bname), ierr)
        done = PETSC_TRUE
        flow%bc%flags(boundary) = BC_NEUMANN
        flow%bc_flags(boundary) = BC_FLUX
        flow%bc_data(:,boundary) = 0.d0

        call OptionsGetReal(options, "-bc_flux_"//trim(bname)//"_value", &
             "flux in the inward-normal direction", flow%bc_data(1,boundary), flag, ierr)

        if (.not.flag) then
          call LBMWarn(flow%comm, "Flux BC on boundary "//trim(bname)// &
               " not set in input file!", ierr)
        end if
      end if
    endif

    ! pressure on the inlet provides pressure and rho1/rho_total
    bcvalue = PETSC_FALSE
    call OptionsGetBool(options,"-bc_pressure_inlet_"//trim(bname), &
         "provide pressure (and rho_1/rho_t if two-component)", &
         bcvalue, flag, ierr)
    if (bcvalue) then
      if (done) call LBMError(flow%comm, 1, &
           "Multiple BCs provided for boundary "//trim(bname), ierr)
      done = PETSC_TRUE
      flow%bc%flags(boundary) = BC_DIRICHLET
      flow%bc_flags(boundary) = BC_PRESSURE_INLET
      flow%bc_data(:,boundary) = -1.d0

      call OptionsGetReal(options, "-bc_pressure_"//trim(bname)//"_value", &
           "pressure", flow%bc_data(1,boundary), flag, ierr)
      if (.not.flag) then
        call LBMWarn(flow%comm, "Pressure BC on boundary "//trim(bname)//&
           " not set in input file!", ierr)
      end if
      if (flow%ncomponents > 1) then
        call OptionsGetReal(options, "-bc_rho1_fraction_"//trim(bname)// &
             "_value", "mass fraction of rho1, rho1/sum(rho)", &
             flow%bc_data(2,boundary), flag, ierr)
        if (.not.flag) then
          call LBMWarn(flow%comm, &
               "Mass fraction for pressure BC on boundary " &
               //trim(bname)//" not set in input file!", ierr)
        end if
      endif
    end if

    ! pressure on the outlet provides pressure only
    bcvalue = PETSC_FALSE
    call OptionsGetBool(options,"-bc_pressure_outlet_"//trim(bname), &
         "provide pressure only", &
         bcvalue, flag, ierr)
    if (bcvalue) then
      if (done) call LBMError(flow%comm, 1, &
           "Multiple BCs provided for boundary "//trim(bname), ierr)
      done = PETSC_TRUE
      flow%bc%flags(boundary) = BC_DIRICHLET
      flow%bc_flags(boundary) = BC_PRESSURE_OUTLET
      flow%bc_data(:,boundary) = -1.d0

      call OptionsGetReal(options, "-bc_pressure_"//trim(bname)//"_value", &
           "pressure", flow%bc_data(1,boundary), flag, ierr)
      if (.not.flag) then
        call LBMWarn(flow%comm, "Pressure BC on boundary "//trim(bname)// &
           " not set in input file!", ierr)
      end if
    end if

    ! momentum boundary condition
    bcvalue = PETSC_FALSE
    call OptionsGetBool(options, "-bc_momentum_"//trim(bname), &
         "provide all components of rho*u", bcvalue, flag, ierr)
    if (bcvalue) then
      if (done) call LBMError(flow%comm, 1, &
           "Multiple BCs provided for boundary "//trim(bname), ierr)
      done = PETSC_TRUE
      flow%bc%flags(boundary) = BC_NEUMANN
      flow%bc_flags(boundary) = BC_MOMENTUM
      flow%bc_data(:,boundary) = -1.d0

      if (flow%ncomponents.eq.1) then
        bccount = flow%ndims
        call OptionsGetRealArray(options,"-bc_momentum_"//trim(bname)//"_value",&
             "list of momentums in x,y,(z) directions", &
             flow%bc_data(:,boundary), bccount, flag, ierr)
        if (.not.flag) then
          call LBMWarn(flow%comm, "Density BC on boundary "//trim(bname)// &
               " not set in input file!", ierr)
        else if (bccount /= flow%ndims) then
          call LBMError(flow%comm,1,"Momentum BC on boundary "//trim(bname)// &
               " was provided the wrong number of components.", ierr)
        end if
      else
        do m=1,flow%ncomponents
          bccount = flow%ndims
          call OptionsGetRealArray(options,"-bc_momentum_"//trim(bname)// &
               "_"//trim(flow%components(m)%name), &
               "list of momentums in x,y,(z) directions", &
               flow%bc_data((m-1)*flow%ndims+1:m*flow%ndims,boundary), &
               bccount, flag, ierr)
          if (.not.flag) then
            call LBMWarn(flow%comm, "Momentum BC on boundary " &
                 //trim(bname)//" not set in input file!", ierr)
          else if (bccount /= flow%ndims) then
            call LBMError(flow%comm, 1, "Momentum BC on boundary "// &
                 trim(bname)//" for component "//trim(flow%components(m)%name)&
                 //" was provided the wrong number of components.", ierr)
          end if
        end do
      end if
    endif

    ! flux inlet provides total flux and rho1/rho_total
    bcvalue = PETSC_FALSE
    call OptionsGetBool(options,"-bc_flux_inlet_"//trim(bname), &
         "provide total flux (and rho_1/rho_t if two-component)", &
         bcvalue, flag, ierr)
    if (bcvalue) then
      if (done) call LBMError(flow%comm, 1, &
           "Multiple BCs provided for boundary "//trim(bname), ierr)
      done = PETSC_TRUE
      flow%bc%flags(boundary) = BC_NEUMANN
      flow%bc_flags(boundary) = BC_FLUX_INLET
      flow%bc_data(:,boundary) = -1.d0

      call OptionsGetReal(options, "-bc_flux_inlet_"//trim(bname)//"_value", &
           "flux", flow%bc_data(1,boundary), flag, ierr)
      if (.not.flag) then
        call LBMWarn(flow%comm, "Flux inlet BC on boundary "//trim(bname)// &
           " not set in input file!", ierr)
      end if

      if (flow%ncomponents > 1) then
        call OptionsGetReal(options, "-bc_rho1_fraction_"//trim(bname)// &
             "_value", "mass fraction of rho1, rho1/sum(rho)", &
             flow%bc_data(2,boundary), flag, ierr)
        if (.not.flag) then
          call LBMWarn(flow%comm, "Mass fraction for flux BC on boundary " &
             //trim(bname)//" not set in input file!", ierr)
        end if
      endif
    end if

    ! flux outlet provides total flux
    bcvalue = PETSC_FALSE
    call OptionsGetBool(options,"-bc_flux_outlet_"//trim(bname), &
         "provide total flux", bcvalue, flag, ierr)
    if (bcvalue) then
      if (done) call LBMError(flow%comm, 1, &
           "Multiple BCs provided for boundary "//trim(bname), ierr)
      done = PETSC_TRUE
      flow%bc%flags(boundary) = BC_NEUMANN
      flow%bc_flags(boundary) = BC_FLUX_OUTLET
      flow%bc_data(:,boundary) = -1.d0

      call OptionsGetReal(options, "-bc_flux_outlet_"//trim(bname)//"_value", &
           "flux", flow%bc_data(1,boundary), flag, ierr)
      if (.not.flag) then
        call LBMWarn(flow%comm, "Flux outlet BC on boundary "//trim(bname)// &
             " not set in input file!", ierr)
      end if
    end if

    ! velocity provides total velocities
    bcvalue = PETSC_FALSE
    call OptionsGetBool(options,"-bc_velocity_"//trim(bname), &
         "provide velocity", bcvalue, flag, ierr)
    if (bcvalue) then
      if (done) call LBMError(flow%comm, 1, &
           "Multiple BCs provided for boundary "//trim(bname), ierr)
      done = PETSC_TRUE
      flow%bc%flags(boundary) = BC_VELOCITY
      flow%bc_flags(boundary) = BC_VELOCITY
      flow%bc_data(:,boundary) = -1.d0

      bccount = flow%ndims
      call OptionsGetRealArray(options, "-bc_velocity_"//trim(bname)//"_values",&
           "velocity in x,y,(z) directions", flow%bc_data(:,boundary), bccount,&
           flag, ierr)
      if (.not.flag) then
        call LBMWarn(flow%comm, "Velocity BC on boundary "//trim(bname)// &
             " not set in input file!", ierr)
      else if (bccount /= flow%ndims) then
        call LBMError(flow%comm, 1, "Velocity BC on boundary "//trim(bname)// &
             " was provided the wrong number of components.", ierr)
      end if
    endif
  end subroutine FlowParseBC

  subroutine FlowSetUpBCsD3(flow, xm_vals, xp_vals, ym_vals, yp_vals, &
       zm_vals, zp_vals, dist)
    type(flow_type) flow
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,flow%ndims,dist%info%ys:dist%info%ye, &
         dist%info%zs:dist%info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(dist%s,flow%ndims,dist%info%xs:dist%info%xe, &
         dist%info%zs:dist%info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(dist%s,flow%ndims,dist%info%xs:dist%info%xe, &
         dist%info%ys:dist%info%ye):: zm_vals, zp_vals

    PetscScalar density(dist%s)
    PetscInt i,j,k,m
    PetscErrorCode ierr

    ! xm boundary
    if (dist%info%xs.eq.1) then
      select case(flow%bc_flags(BOUNDARY_XM))
      case(BC_DENSITY)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          xm_vals(:,1,j,k) = flow%bc_data(1:dist%s,BOUNDARY_XM)
        end do
        end do
      case(BC_PRESSURE_INLET, BC_PRESSURE)
        density(:) = 0.d0
        call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_XM), &
             flow%bc_data(2,BOUNDARY_XM), density)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          xm_vals(:,1,j,k) = density
        end do
        end do
      case(BC_PRESSURE_OUTLET)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          xm_vals(1,1,j,k) = flow%bc_data(1,BOUNDARY_XM)
        end do
        end do
      case(BC_MOMENTUM)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          do m=1,dist%s
            xm_vals(m,:,j,k) = flow%bc_data((m-1)*flow%ndims+1:m*flow%ndims, &
                 BOUNDARY_XM)
          end do
        end do
        end do
      case(BC_FLUX_INLET,BC_FLUX)
        if (dist%s.eq.1) then
          do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
            xm_vals(1,X_DIRECTION,j,k) = flow%bc_data(1,BOUNDARY_XM)
          end do
          end do
        else if (dist%s.eq.2) then
          do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
            xm_vals(1,X_DIRECTION,j,k) = flow%bc_data(1,BOUNDARY_XM) &
                 * flow%bc_data(2,BOUNDARY_XM)
            xm_vals(2,X_DIRECTION,j,k) = flow%bc_data(1,BOUNDARY_XM) &
                 * (1.d0 - flow%bc_data(2,BOUNDARY_XM))
          end do
          end do
        else
          call LBMError(flow%comm, 1, &
               "BC flux inlet only implemented for 1 or 2 components.", ierr)
        end if
      case(BC_FLUX_OUTLET)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          xm_vals(1,X_DIRECTION,j,k) = flow%bc_data(1,BOUNDARY_XM)
        end do
        end do
      case(BC_VELOCITY)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          xm_vals(1,:,j,k) = flow%bc_data(1:flow%ndims,BOUNDARY_XM)
        end do
        end do
      end select
    end if

    ! xp boundary
    if (dist%info%xe.eq.dist%info%NX) then
      select case(flow%bc_flags(BOUNDARY_XP))
      case(BC_DENSITY)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          xp_vals(:,1,j,k) = flow%bc_data(1:dist%s,BOUNDARY_XP)
        end do
        end do
      case(BC_PRESSURE_INLET,BC_PRESSURE)
        density(:) = 0.d0
        call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_XP), &
             flow%bc_data(2,BOUNDARY_XP), density)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          xp_vals(:,1,j,k) = density
        end do
        end do
      case(BC_PRESSURE_OUTLET)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          xp_vals(1,1,j,k) = flow%bc_data(1,BOUNDARY_XP)
        end do
        end do
      case(BC_MOMENTUM)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          do m=1,dist%s
            xp_vals(m,:,j,k) = flow%bc_data((m-1)*flow%ndims+1:m*flow%ndims, &
                 BOUNDARY_XP)
          end do
        end do
        end do
      case(BC_FLUX_INLET,BC_FLUX)
        if (dist%s.eq.1) then
          do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
            xp_vals(1,X_DIRECTION,j,k) = -flow%bc_data(1,BOUNDARY_XP)
          end do
          end do
        else if (dist%s.eq.2) then
          do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
            xp_vals(1,X_DIRECTION,j,k) = flow%bc_data(1,BOUNDARY_XP) &
                 * flow%bc_data(2,BOUNDARY_XP)
            xp_vals(2,X_DIRECTION,j,k) = flow%bc_data(1,BOUNDARY_XP) &
                 * (1.d0 - flow%bc_data(2,BOUNDARY_XP))
          end do
          end do
        else
          call LBMError(flow%comm, 1, &
               "BC flux inlet only implemented for 1 or 2 components.", ierr)
        end if
      case(BC_FLUX_OUTLET)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          xp_vals(1,X_DIRECTION,j,k) = -flow%bc_data(1,BOUNDARY_XP)
        end do
        end do
      case(BC_VELOCITY)
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          xp_vals(1,:,j,k) = flow%bc_data(1:flow%ndims,BOUNDARY_XP)
        end do
        end do
      end select
    end if

    ! ym boundary
    if (dist%info%ys.eq.1) then
      select case(flow%bc_flags(BOUNDARY_YM))
      case(BC_DENSITY)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          ym_vals(:,1,i,k) = flow%bc_data(1:dist%s,BOUNDARY_YM)
        end do
        end do
      case(BC_PRESSURE_INLET,BC_PRESSURE)
        density(:) = 0.d0
        call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_YM), &
             flow%bc_data(2,BOUNDARY_YM), density)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          ym_vals(:,1,i,k) = density
        end do
        end do
      case(BC_PRESSURE_OUTLET)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          ym_vals(1,1,i,k) = flow%bc_data(1,BOUNDARY_YM)
        end do
        end do
      case(BC_MOMENTUM)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          do m=1,dist%s
            ym_vals(m,:,i,k) = flow%bc_data((m-1)*flow%ndims+1:m*flow%ndims, &
                 BOUNDARY_YM)
          end do
        end do
        end do
      case(BC_FLUX_INLET,BC_FLUX)
        if (dist%s.eq.1) then
          do k=dist%info%zs,dist%info%ze
          do i=dist%info%xs,dist%info%xe
            ym_vals(1,Y_DIRECTION,i,k) = flow%bc_data(1,BOUNDARY_YM)
          end do
          end do
        else if (dist%s.eq.2) then
          do k=dist%info%zs,dist%info%ze
          do i=dist%info%xs,dist%info%xe
            ym_vals(1,Y_DIRECTION,i,k) = flow%bc_data(1,BOUNDARY_YM) &
                 * flow%bc_data(2,BOUNDARY_YM)
            ym_vals(2,Y_DIRECTION,i,k) = flow%bc_data(1,BOUNDARY_YM) &
                 * (1.d0 - flow%bc_data(2,BOUNDARY_YM))
          end do
          end do
        else
          call LBMError(flow%comm, 1, &
               "BC flux inlet only implemented for 1 or 2 components.", ierr)
        end if
      case(BC_FLUX_OUTLET)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          ym_vals(1,Y_DIRECTION,i,k) = flow%bc_data(1,BOUNDARY_YM)
        end do
        end do
      case(BC_VELOCITY)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          ym_vals(1,:,i,k) = flow%bc_data(1:flow%ndims,BOUNDARY_YM)
        end do
        end do
      end select
    end if

    ! yp boundary
    if (dist%info%ye.eq.dist%info%NY) then
      select case(flow%bc_flags(BOUNDARY_YP))
      case(BC_DENSITY)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          yp_vals(:,1,i,k) = flow%bc_data(1:dist%s,BOUNDARY_YP)
        end do
        end do
      case(BC_PRESSURE_INLET,BC_PRESSURE)
        density(:) = 0.d0
        call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_YP), &
             flow%bc_data(2,BOUNDARY_YP), density)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          yp_vals(:,1,i,k) = density
        end do
        end do
      case(BC_PRESSURE_OUTLET)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          yp_vals(1,1,i,k) = flow%bc_data(1,BOUNDARY_YP)
        end do
        end do
      case(BC_MOMENTUM)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          do m=1,dist%s
            yp_vals(m,:,i,k) = flow%bc_data((m-1)*flow%ndims+1:m*flow%ndims, &
                 BOUNDARY_YP)
          end do
        end do
        end do
      case(BC_FLUX_INLET,BC_FLUX)
        if (dist%s.eq.1) then
          do k=dist%info%zs,dist%info%ze
          do i=dist%info%xs,dist%info%xe
            yp_vals(1,Y_DIRECTION,i,k) = -flow%bc_data(1,BOUNDARY_YP)
          end do
          end do
        else if (dist%s.eq.2) then
          do k=dist%info%zs,dist%info%ze
          do i=dist%info%xs,dist%info%xe
            yp_vals(1,Y_DIRECTION,i,k) = flow%bc_data(1,BOUNDARY_YP) &
                 * flow%bc_data(2,BOUNDARY_YP)
            yp_vals(2,Y_DIRECTION,i,k) = flow%bc_data(1,BOUNDARY_YP) &
                 * (1.d0 - flow%bc_data(2,BOUNDARY_YP))
          end do
          end do
        else
          call LBMError(flow%comm, 1, &
               "BC flux inlet only implemented for 1 or 2 components.", ierr)
        end if
      case(BC_FLUX_OUTLET)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          yp_vals(1,Y_DIRECTION,i,k) = -flow%bc_data(1,BOUNDARY_YP)
        end do
        end do
      case(BC_VELOCITY)
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          yp_vals(1,:,i,k) = flow%bc_data(1:flow%ndims,BOUNDARY_YP)
        end do
        end do
      end select
    end if

    ! zm boundary
    if (dist%info%zs.eq.1) then
      select case(flow%bc_flags(BOUNDARY_ZM))
      case(BC_DENSITY)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          zm_vals(:,1,i,j) = flow%bc_data(1:dist%s,BOUNDARY_ZM)
        end do
        end do
      case(BC_PRESSURE_INLET,BC_PRESSURE)
        density(:) = 0.d0
        call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_ZM), &
             flow%bc_data(2,BOUNDARY_ZM), density)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          zm_vals(:,1,i,j) = density
        end do
        end do
      case(BC_PRESSURE_OUTLET)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          zm_vals(1,1,i,j) = flow%bc_data(1,BOUNDARY_ZM)
        end do
        end do
      case(BC_MOMENTUM)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          do m=1,dist%s
            zm_vals(m,:,i,j) = flow%bc_data((m-1)*flow%ndims+1:m*flow%ndims, &
                 BOUNDARY_ZM)
          end do
        end do
        end do
      case(BC_FLUX_INLET,BC_FLUX)
        if (dist%s.eq.1) then
          do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
            zm_vals(1,Z_DIRECTION,i,j) = flow%bc_data(1,BOUNDARY_ZM)
          end do
          end do
        else if (dist%s.eq.2) then
          do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
            zm_vals(1,Z_DIRECTION,i,j) = flow%bc_data(1,BOUNDARY_ZM) &
                 * flow%bc_data(2,BOUNDARY_ZM)
            zm_vals(2,Z_DIRECTION,i,j) = flow%bc_data(1,BOUNDARY_ZM) &
                 * (1.d0 - flow%bc_data(2,BOUNDARY_ZM))
          end do
          end do
        else
          call LBMError(flow%comm, 1, &
               "BC flux inlet only implemented for 1 or 2 components.", ierr)
        end if
      case(BC_FLUX_OUTLET)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          zm_vals(1,Z_DIRECTION,i,j) = flow%bc_data(1,BOUNDARY_ZM)
        end do
        end do
      case(BC_VELOCITY)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          zm_vals(1,:,i,j) = flow%bc_data(1:flow%ndims,BOUNDARY_ZM)
        end do
        end do
      end select
    end if

    ! zp boundary
    if (dist%info%ze.eq.dist%info%NZ) then
      select case(flow%bc_flags(BOUNDARY_ZP))
      case(BC_DENSITY)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          zp_vals(:,1,i,j) = flow%bc_data(1:dist%s,BOUNDARY_ZP)
        end do
        end do
      case(BC_PRESSURE_INLET,BC_PRESSURE)
        density(:) = 0.d0
        call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_ZP), &
             flow%bc_data(2,BOUNDARY_ZP), density)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          zp_vals(:,1,i,j) = density
        end do
        end do
      case(BC_PRESSURE_OUTLET)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          zp_vals(1,1,i,j) = flow%bc_data(1,BOUNDARY_ZP)
        end do
        end do
      case(BC_MOMENTUM)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          do m=1,dist%s
            zp_vals(m,:,i,j) = flow%bc_data((m-1)*flow%ndims+1:m*flow%ndims, &
                 BOUNDARY_ZP)
          end do
        end do
        end do
      case(BC_FLUX_INLET,BC_FLUX)
        if (dist%s.eq.1) then
          do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
            zp_vals(1,Z_DIRECTION,i,j) = -flow%bc_data(1,BOUNDARY_ZP)
          end do
          end do
        else if (dist%s.eq.2) then
          do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
            zp_vals(1,Z_DIRECTION,i,j) = flow%bc_data(1,BOUNDARY_ZP) &
                 * flow%bc_data(2,BOUNDARY_ZP)
            zp_vals(2,Z_DIRECTION,i,j) = flow%bc_data(1,BOUNDARY_ZP) &
                 * (1.d0 - flow%bc_data(2,BOUNDARY_ZP))
          end do
          end do
        else
          call LBMError(flow%comm, 1, &
               "BC flux inlet only implemented for 1 or 2 components.", ierr)
        end if
      case(BC_FLUX_OUTLET)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          zp_vals(1,Z_DIRECTION,i,j) = -flow%bc_data(1,BOUNDARY_ZP)
        end do
        end do
      case(BC_VELOCITY)
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          zp_vals(1,:,i,j) = flow%bc_data(1:flow%ndims,BOUNDARY_ZP)
        end do
        end do
      end select
    end if
  end subroutine FlowSetUpBCsD3

  subroutine FlowSetUpBCsD2(flow, xm_vals, xp_vals, ym_vals, yp_vals, dist)
    type(flow_type) flow
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,flow%ndims, &
         dist%info%ys:dist%info%ye):: xm_vals, xp_vals
    PetscScalar,dimension(dist%s,flow%ndims, &
         dist%info%xs:dist%info%xe):: ym_vals, yp_vals

    PetscScalar density(dist%s)
    PetscInt i,j,m
    PetscErrorCode ierr

    ! xm boundary
    if (dist%info%xs.eq.1) then
      select case(flow%bc_flags(BOUNDARY_XM))
      case(BC_DENSITY)
        do j=dist%info%ys,dist%info%ye
          xm_vals(:,1,j) = flow%bc_data(1:dist%s,BOUNDARY_XM)
        end do
      case(BC_PRESSURE_INLET,BC_PRESSURE)
        density(:) = 0.d0
        call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_XM), &
             flow%bc_data(2,BOUNDARY_XM), density)
        do j=dist%info%ys,dist%info%ye
          xm_vals(:,1,j) = density
        end do
      case(BC_PRESSURE_OUTLET)
        do j=dist%info%ys,dist%info%ye
          xm_vals(1,1,j) = flow%bc_data(1,BOUNDARY_XM)
        end do
      case(BC_MOMENTUM)
        do j=dist%info%ys,dist%info%ye
          do m=1,dist%s
            xm_vals(m,:,j) = flow%bc_data((m-1)*flow%ndims+1:m*flow%ndims, &
                 BOUNDARY_XM)
          end do
        end do
      case(BC_FLUX_INLET,BC_FLUX)
        if (dist%s.eq.1) then
          do j=dist%info%ys,dist%info%ye
            xm_vals(1,X_DIRECTION,j) = flow%bc_data(1,BOUNDARY_XM)
          end do
        else if (dist%s.eq.2) then
          do j=dist%info%ys,dist%info%ye
            xm_vals(1,X_DIRECTION,j) = flow%bc_data(1,BOUNDARY_XM) &
                 * flow%bc_data(2,BOUNDARY_XM)
            xm_vals(2,X_DIRECTION,j) = flow%bc_data(1,BOUNDARY_XM) &
                 * (1.d0 - flow%bc_data(2,BOUNDARY_XM))
          end do
        else
          call LBMError(flow%comm, 1, &
               "BC flux inlet only implemented for 1 or 2 components.", ierr)
        end if
      case(BC_FLUX_OUTLET)
        do j=dist%info%ys,dist%info%ye
          xm_vals(1,X_DIRECTION,j) = flow%bc_data(1,BOUNDARY_XM)
        end do
      case(BC_VELOCITY)
        do j=dist%info%ys,dist%info%ye
          xm_vals(1,:,j) = flow%bc_data(1:flow%ndims,BOUNDARY_XM)
        end do
      end select
    end if

    ! xp boundary
    if (dist%info%xe.eq.dist%info%NX) then
      select case(flow%bc_flags(BOUNDARY_XP))
      case(BC_DENSITY)
        do j=dist%info%ys,dist%info%ye
          xp_vals(:,1,j) = flow%bc_data(1:dist%s,BOUNDARY_XP)
        end do
      case(BC_PRESSURE_INLET,BC_PRESSURE)
        density(:) = 0.d0
        call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_XP), &
             flow%bc_data(2,BOUNDARY_XP), density)
        do j=dist%info%ys,dist%info%ye
          xp_vals(:,1,j) = density
        end do
      case(BC_PRESSURE_OUTLET)
        do j=dist%info%ys,dist%info%ye
          xp_vals(1,1,j) = flow%bc_data(1,BOUNDARY_XP)
        end do
      case(BC_MOMENTUM)
        do j=dist%info%ys,dist%info%ye
          do m=1,dist%s
            xp_vals(m,:,j) = flow%bc_data((m-1)*flow%ndims+1:m*flow%ndims, &
                 BOUNDARY_XP)
          end do
        end do
      case(BC_FLUX_INLET,BC_FLUX)
        if (dist%s.eq.1) then
          do j=dist%info%ys,dist%info%ye
            ! flux given in the inward normal direction, but uses
            ! bc_neumann which takes things in x,y,z coordinates
            xp_vals(1,X_DIRECTION,j) = -flow%bc_data(1,BOUNDARY_XP)
          end do
        else if (dist%s.eq.2) then
          do j=dist%info%ys,dist%info%ye
            xp_vals(1,X_DIRECTION,j) = flow%bc_data(1,BOUNDARY_XP) &
                 * flow%bc_data(2,BOUNDARY_XP)
            xp_vals(2,X_DIRECTION,j) = flow%bc_data(1,BOUNDARY_XP) &
                 * (1.d0 - flow%bc_data(2,BOUNDARY_XP))
          end do
        else
          call LBMError(flow%comm, 1, &
               "BC flux inlet only implemented for 1 or 2 components.", ierr)
        end if
      case(BC_FLUX_OUTLET)
        do j=dist%info%ys,dist%info%ye
          xp_vals(1,X_DIRECTION,j) = -flow%bc_data(1,BOUNDARY_XP)
        end do
      case(BC_VELOCITY)
        do j=dist%info%ys,dist%info%ye
          xp_vals(1,:,j) = flow%bc_data(1:flow%ndims,BOUNDARY_XP)
        end do
      end select
    end if

    ! ym boundary
    if (dist%info%ys.eq.1) then
      select case(flow%bc_flags(BOUNDARY_YM))
      case(BC_DENSITY)
        do i=dist%info%xs,dist%info%xe
          ym_vals(:,1,i) = flow%bc_data(1:dist%s,BOUNDARY_YM)
        end do
      case(BC_PRESSURE_INLET,BC_PRESSURE)
        density(:) = 0.d0
        call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_YM), &
             flow%bc_data(2,BOUNDARY_YM), density)
        do i=dist%info%xs,dist%info%xe
          ym_vals(:,1,i) = density
        end do
      case(BC_PRESSURE_OUTLET)
        do i=dist%info%xs,dist%info%xe
          ym_vals(1,1,i) = flow%bc_data(1,BOUNDARY_YM)
        end do
      case(BC_MOMENTUM)
        do i=dist%info%xs,dist%info%xe
          do m=1,dist%s
            ym_vals(m,:,i) = flow%bc_data((m-1)*flow%ndims+1:m*flow%ndims, &
                 BOUNDARY_YM)
          end do
        end do
      case(BC_FLUX_INLET,BC_FLUX)
        if (dist%s.eq.1) then
          do i=dist%info%xs,dist%info%xe
            ym_vals(1,Y_DIRECTION,i) = flow%bc_data(1,BOUNDARY_YM)
          end do
        else if (dist%s.eq.2) then
          do i=dist%info%xs,dist%info%xe
            ym_vals(1,Y_DIRECTION,i) = flow%bc_data(1,BOUNDARY_YM) &
                 * flow%bc_data(2,BOUNDARY_YM)
            ym_vals(2,Y_DIRECTION,i) = flow%bc_data(1,BOUNDARY_YM) &
                 * (1.d0 - flow%bc_data(2,BOUNDARY_YM))
          end do
        else
          call LBMError(flow%comm, 1, &
               "BC flux inlet only implemented for 1 or 2 components.", ierr)
        end if
      case(BC_FLUX_OUTLET)
        do i=dist%info%xs,dist%info%xe
          ym_vals(1,Y_DIRECTION,i) = flow%bc_data(1,BOUNDARY_YM)
        end do
      case(BC_VELOCITY)
        do i=dist%info%xs,dist%info%xe
          ym_vals(1,:,i) = flow%bc_data(1:flow%ndims,BOUNDARY_YM)
        end do
      end select
    end if

    ! yp boundary
    if (dist%info%ye.eq.dist%info%NY) then
      select case(flow%bc_flags(BOUNDARY_YP))
      case(BC_DENSITY)
        do i=dist%info%xs,dist%info%xe
          yp_vals(:,1,i) = flow%bc_data(1:dist%s,BOUNDARY_YP)
        end do
      case(BC_PRESSURE_INLET,BC_PRESSURE)
        density(:) = 0.d0
        call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_YP), &
             flow%bc_data(2,BOUNDARY_YP), density)
        do i=dist%info%xs,dist%info%xe
          yp_vals(:,1,i) = density
        end do
      case(BC_PRESSURE_OUTLET)
        do i=dist%info%xs,dist%info%xe
          yp_vals(1,1,i) = flow%bc_data(1,BOUNDARY_YP)
        end do
      case(BC_MOMENTUM)
        do i=dist%info%xs,dist%info%xe
          do m=1,dist%s
            yp_vals(m,:,i) = flow%bc_data((m-1)*flow%ndims+1:m*flow%ndims, &
                 BOUNDARY_YP)
          end do
        end do
      case(BC_FLUX_INLET,BC_FLUX)
        if (dist%s.eq.1) then
          do i=dist%info%xs,dist%info%xe
            yp_vals(1,Y_DIRECTION,i) = -flow%bc_data(1,BOUNDARY_YP)
          end do
        else if (dist%s.eq.2) then
          do i=dist%info%xs,dist%info%xe
            yp_vals(1,Y_DIRECTION,i) = flow%bc_data(1,BOUNDARY_YP) &
                 * flow%bc_data(2,BOUNDARY_YP)
            yp_vals(2,Y_DIRECTION,i) = flow%bc_data(1,BOUNDARY_YP) &
                 * (1.d0 - flow%bc_data(2,BOUNDARY_YP))
          end do
        else
          call LBMError(flow%comm, 1, &
               "BC flux inlet only implemented for 1 or 2 components.", ierr)
        end if
      case(BC_FLUX_OUTLET)
        do i=dist%info%xs,dist%info%xe
          yp_vals(1,Y_DIRECTION,i) = -flow%bc_data(1,BOUNDARY_YP)
        end do
      case(BC_VELOCITY)
        do i=dist%info%xs,dist%info%xe
          yp_vals(1,:,i) = flow%bc_data(1:flow%ndims,BOUNDARY_YP)
        end do
      end select
    end if
  end subroutine FlowSetUpBCsD2

  subroutine FlowApplyBCs(flow, walls)
    type(flow_type) flow
    type(walls_type) walls

    PetscInt lcv_side

    ! some BCs are special cases and need preprocessing
    do lcv_side=1,flow%ndims*2
      if ((flow%bc_flags(lcv_side).eq.BC_PRESSURE_OUTLET).and. &
           .not.flow%bc_done(BC_PRESSURE_OUTLET)) then
        if (flow%ncomponents /= 1) then
          call FlowUpdateBCPressureOutlet(flow, walls)
        end if
        flow%bc_done(BC_PRESSURE_OUTLET) = PETSC_TRUE
      else if ((flow%bc_flags(lcv_side).eq.BC_FLUX_OUTLET).and. &
           .not.flow%bc_done(BC_FLUX_OUTLET)) then
        call FlowUpdateBCFluxOutlet(flow, walls)
        flow%bc_done(BC_FLUX_OUTLET) = PETSC_TRUE
      end if
    end do
    flow%bc_done(BC_FLUX_OUTLET) = PETSC_FALSE
    flow%bc_done(BC_PRESSURE_OUTLET) = PETSC_FALSE

    ! Calculate rho, using the streamed values for internal, the rho
    ! boundary condition for Dirichlet BCs, and f^* as the previous
    ! timestep value.  Then calculate forces from that rho.
    call FlowCalcRhoForces(flow, walls)

    ! apply boundary conditions, using the calculated forces
    call BCApply(flow%bc, flow%forces, walls%walls_a, flow%distribution)

    ! update rho on the boundary using the new values
    call BCUpdateRho(flow%bc, walls%walls_a, flow%distribution)
  end subroutine FlowApplyBCs

  subroutine FlowUpdateBCPressureOutlet(flow, walls)
    type(flow_type) flow
    type(walls_type) walls

    PetscErrorCode ierr
    PetscScalar,parameter:: eps=1.e-10 ! slightly larger than machine epsilon

    if ((DABS(flow%components(1)%gf(1)) > eps) .or. &
        (DABS(flow%components(1)%gf(1)) > eps) .or. &
        flow%use_nonideal_eos) then
      call LBMError(PETSC_COMM_SELF, 1, &
           "Pressure outlet not implemented for non-ideal EOS or g11 or g22 /= 0", ierr)
    end if

    select case(flow%ndims)
    case (2)
      call FlowUpdateBCPressureOutletD2(flow, walls%walls_a, &
           flow%distribution%rho_a, flow%bc%xm_a, flow%bc%xp_a, &
           flow%bc%ym_a, flow%bc%yp_a, flow%distribution)
    case (3)
      call FlowUpdateBCPressureOutletD3(flow, walls%walls_a, &
           flow%distribution%rho_a, flow%bc%xm_a, flow%bc%xp_a, &
           flow%bc%ym_a, flow%bc%yp_a, flow%bc%zm_a, flow%bc%zp_a, &
           flow%distribution)
    end select
  end subroutine FlowUpdateBCPressureOutlet

  subroutine FlowUpdateBCPressureOutletD2(flow, walls, rho, xm_vals, xp_vals, &
           ym_vals, yp_vals, dist)
    type(flow_type) flow
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls
    PetscScalar,dimension(flow%bc%nbcs,dist%info%ys:dist%info%ye):: xm_vals, xp_vals
    PetscScalar,dimension(flow%bc%nbcs,dist%info%xs:dist%info%xe):: ym_vals, yp_vals

    PetscInt i,j
    PetscScalar rho1frac

    ! Enforce the condition that rho1/rho_tot is fixed

    ! xm boundary
    if ((flow%bc_flags(BOUNDARY_XM).eq.BC_PRESSURE_OUTLET).and.(dist%info%xs.eq.1)) then
      i = 1
      do j=dist%info%ys,dist%info%ye
        if (walls(i,j).eq.0) then
          if (walls(i+1,j).eq.0) then
            rho1frac = rho(1,i+1,j)/sum(rho(:,i+1,j))
          else
            rho1frac = rho(1,i,j)/sum(rho(:,i,j))
          endif
          call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_XM), &
               rho1frac, xm_vals(:,j))
        endif
      end do
    endif

    ! xp boundary
    if ((flow%bc_flags(BOUNDARY_XP).eq.BC_PRESSURE_OUTLET).and. &
         (dist%info%xe.eq.dist%info%NX)) then
      i = dist%info%NX
      do j=dist%info%ys,dist%info%ye
        if (walls(i,j).eq.0) then
          if (walls(i-1,j).eq.0) then
            rho1frac = rho(1,i-1,j)/sum(rho(:,i-1,j))
          else
            rho1frac = rho(1,i,j)/sum(rho(:,i,j))
          endif
          call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_XP), &
               rho1frac, xp_vals(:,j))
        end if
      end do
    endif

    ! ym boundary
    if ((flow%bc_flags(BOUNDARY_YM).eq.BC_PRESSURE_OUTLET).and.(dist%info%ys.eq.1)) then
      j = 1
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j).eq.0) then
          if (walls(i,j+1).eq.0) then
            rho1frac = rho(1,i,j+1)/sum(rho(:,i,j+1))
          else
            rho1frac = rho(1,i,j)/sum(rho(:,i,j))
          endif
          call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_YM), &
               rho1frac, ym_vals(:,i))
        end if
      end do
    endif

    ! yp boundary
    if ((flow%bc_flags(BOUNDARY_YP).eq.BC_PRESSURE_OUTLET).and. &
         (dist%info%ye.eq.dist%info%NY)) then
      j = dist%info%NY
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j).eq.0) then
          if (walls(i,j-1).eq.0) then
            rho1frac = rho(1,i,j-1)/sum(rho(:,i,j-1))
          else
            rho1frac = rho(1,i,j)/sum(rho(:,i,j))
          endif
          call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_YP), &
               rho1frac, yp_vals(:,i))
        endif
      end do
    endif
  end subroutine FlowUpdateBCPressureOutletD2

  subroutine FlowUpdateBCPressureOutletD3(flow, walls, rho, xm_vals, xp_vals, &
           ym_vals, yp_vals, zm_vals, zp_vals, dist)
    type(flow_type) flow
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye,dist%info%rgzs:dist%info%rgze):: rho
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye,dist%info%rgzs:dist%info%rgze):: walls
    PetscScalar,dimension(flow%bc%nbcs,dist%info%ys:dist%info%ye, &
         dist%info%zs:dist%info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(flow%bc%nbcs,dist%info%xs:dist%info%xe, &
         dist%info%zs:dist%info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(flow%bc%nbcs,dist%info%xs:dist%info%xe, &
         dist%info%ys:dist%info%ye):: zm_vals, zp_vals

    PetscInt i,j,k
    PetscScalar rho1frac

    ! Enforce the condition that rho1/rho_tot is fixed

    ! xm boundary
    if ((flow%bc_flags(BOUNDARY_XM).eq.BC_PRESSURE_OUTLET).and.(dist%info%xs.eq.1)) then
      i = 1
      do k=dist%info%zs,dist%info%ze
      do j=dist%info%ys,dist%info%ye
        if (walls(i,j,k).eq.0) then
          if (walls(i+1,j,k).eq.0) then
            rho1frac = rho(1,i+1,j,k)/sum(rho(:,i+1,j,k))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_XM), &
               rho1frac, xm_vals(:,j,k))
        endif
      end do
      end do
    endif

    ! xp boundary
    if ((flow%bc_flags(BOUNDARY_XP).eq.BC_PRESSURE_OUTLET).and. &
         (dist%info%xe.eq.dist%info%NX)) then
      i = dist%info%NX
      do k=dist%info%zs,dist%info%ze
      do j=dist%info%ys,dist%info%ye
        if (walls(i,j,k).eq.0) then
          if (walls(i-1,j,k).eq.0) then
            rho1frac = rho(1,i-1,j,k)/sum(rho(:,i-1,j,k))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_XP), &
               rho1frac, xp_vals(:,j,k))
        end if
      end do
      end do
    endif

    ! ym boundary
    if ((flow%bc_flags(BOUNDARY_YM).eq.BC_PRESSURE_OUTLET).and.(dist%info%ys.eq.1)) then
      j = 1
      do k=dist%info%zs,dist%info%ze
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j,k).eq.0) then
          if (walls(i,j+1,k).eq.0) then
            rho1frac = rho(1,i,j+1,k)/sum(rho(:,i,j+1,k))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_YM), &
               rho1frac, ym_vals(:,i,k))
        end if
      end do
      end do
    endif

    ! yp boundary
    if ((flow%bc_flags(BOUNDARY_YP).eq.BC_PRESSURE_OUTLET).and. &
         (dist%info%ye.eq.dist%info%NY)) then
      j = dist%info%NY
      do k=dist%info%zs,dist%info%ze
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j,k).eq.0) then
          if (walls(i,j-1,k).eq.0) then
            rho1frac = rho(1,i,j-1,k)/sum(rho(:,i,j-1,k))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_YP), &
               rho1frac, yp_vals(:,i,k))
        endif
      end do
      end do
    endif

    ! zm boundary
    if ((flow%bc_flags(BOUNDARY_ZM).eq.BC_PRESSURE_OUTLET).and.(dist%info%zs.eq.1)) then
      k = 1
      do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j,k).eq.0) then
          if (walls(i,j,k+1).eq.0) then
            rho1frac = rho(1,i,j,k+1)/sum(rho(:,i,j,k+1))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_ZM), &
               rho1frac, zm_vals(:,i,j))
        end if
      end do
      end do
    endif

    ! zp boundary
    if ((flow%bc_flags(BOUNDARY_ZP).eq.BC_PRESSURE_OUTLET).and. &
         (dist%info%ze.eq.dist%info%NZ)) then
      k = dist%info%NZ
      do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j,k).eq.0) then
          if (walls(i,j,k-1).eq.0) then
            rho1frac = rho(1,i,j,k-1)/sum(rho(:,i,j,k-1))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          call FlowUpdateDensityFromPressure(flow, flow%bc_data(1,BOUNDARY_ZP), &
               rho1frac, zp_vals(:,i,j))
        endif
      end do
      end do
    endif
  end subroutine FlowUpdateBCPressureOutletD3

  subroutine FlowUpdateDensityFromPressure(flow, pressure, rho1frac, rho)
    type(flow_type) flow
    PetscScalar pressure, rho1frac
    PetscScalar,dimension(flow%ncomponents) :: rho
    PetscScalar alpha
    PetscErrorCode ierr
    PetscScalar :: eps = 1.d-10

    if (flow%ncomponents.eq.1) then
      rho(1) = pressure*3.d0
    else if (flow%ncomponents.eq.2) then
      if (rho1frac < eps) then
        rho(1) = 0.d0
        rho(2) = pressure*3.d0
      else if (rho1frac > 1-eps) then
        rho(1) = pressure*3.d0
        rho(2) = 0.d0
      else
        alpha = 1.d0/(1.d0/rho1frac - 1.d0)

        rho(1) = (-(1.d0 + alpha)/3.d0 + DSQRT((1.d0 + alpha)/3.d0*(1.d0 + alpha)/3.d0&
             + 4.*flow%disc%c_0 *flow%components(2)%gf(1)*alpha*pressure)) / &
             (2*flow%disc%c_0*flow%components(2)%gf(1))
        rho(2) = rho(1)/alpha
      end if
    else
      call LBMError(flow%comm, 1, "Invalid number of components for flow", ierr)
    end if
  end subroutine FlowUpdateDensityFromPressure

  subroutine FlowUpdateBCFluxOutlet(flow, walls)
    type(flow_type) flow
    type(walls_type) walls

    select case(flow%ndims)
    case (2)
      call FlowUpdateBCFluxOutletD2(flow, walls%walls_a, &
           flow%distribution%rho_a, flow%bc%xm_a, flow%bc%xp_a, &
           flow%bc%ym_a, flow%bc%yp_a, flow%distribution)
    case (3)
      call FlowUpdateBCFluxOutletD3(flow, walls%walls_a, &
           flow%distribution%rho_a, flow%bc%xm_a, flow%bc%xp_a, &
           flow%bc%ym_a, flow%bc%yp_a, flow%bc%zm_a, flow%bc%zp_a, &
           flow%distribution)
    end select
  end subroutine FlowUpdateBCFluxOutlet

  subroutine FlowUpdateBCFluxOutletD2(flow, walls, rho, xm_vals, xp_vals, &
           ym_vals, yp_vals, dist)
    type(flow_type) flow
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls
    PetscScalar,dimension(dist%s,flow%ndims, &
         dist%info%ys:dist%info%ye):: xm_vals, xp_vals
    PetscScalar,dimension(dist%s,flow%ndims, &
         dist%info%xs:dist%info%xe):: ym_vals, yp_vals

    PetscInt i,j
    PetscScalar rho1frac

    ! Enforce the condition that rho1/rho_tot is fixed

    ! xm boundary
    if ((flow%bc_flags(BOUNDARY_XM).eq.BC_PRESSURE_OUTLET).and.(dist%info%xs.eq.1)) then
      i = 1
      do j=dist%info%ys,dist%info%ye
        if (walls(i,j).eq.0) then
          if (walls(i+1,j).eq.0) then
            rho1frac = rho(1,i+1,j)/sum(rho(:,i+1,j))
          else
            rho1frac = rho(1,i,j)/sum(rho(:,i,j))
          endif
          xm_vals(1,X_DIRECTION,j) = rho1frac*flow%bc_data(1,BOUNDARY_XM)
          xm_vals(2,X_DIRECTION,j) = (1.d0-rho1frac)*flow%bc_data(1,BOUNDARY_XM)
        endif
      end do
    endif

    ! xp boundary
    if ((flow%bc_flags(BOUNDARY_XP).eq.BC_PRESSURE_OUTLET).and. &
         (dist%info%xe.eq.dist%info%NX)) then
      i = dist%info%NX
      do j=dist%info%ys,dist%info%ye
        if (walls(i,j).eq.0) then
          if (walls(i-1,j).eq.0) then
            rho1frac = rho(1,i-1,j)/sum(rho(:,i-1,j))
          else
            rho1frac = rho(1,i,j)/sum(rho(:,i,j))
          endif
          xp_vals(1,X_DIRECTION,j) = rho1frac*flow%bc_data(1,BOUNDARY_XP)
          xp_vals(2,X_DIRECTION,j) = (1.d0-rho1frac)*flow%bc_data(1,BOUNDARY_XP)
        end if
      end do
    endif

    ! ym boundary
    if ((flow%bc_flags(BOUNDARY_YM).eq.BC_PRESSURE_OUTLET).and.(dist%info%ys.eq.1)) then
      j = 1
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j).eq.0) then
          if (walls(i,j+1).eq.0) then
            rho1frac = rho(1,i,j+1)/sum(rho(:,i,j+1))
          else
            rho1frac = rho(1,i,j)/sum(rho(:,i,j))
          endif
          ym_vals(1,Y_DIRECTION,i) = rho1frac*flow%bc_data(1,BOUNDARY_YM)
          ym_vals(2,Y_DIRECTION,i) = (1.d0-rho1frac)*flow%bc_data(1,BOUNDARY_YM)
        end if
      end do
    endif

    ! yp boundary
    if ((flow%bc_flags(BOUNDARY_YP).eq.BC_PRESSURE_OUTLET).and. &
         (dist%info%ye.eq.dist%info%NY)) then
      j = dist%info%NY
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j).eq.0) then
          if (walls(i,j-1).eq.0) then
            rho1frac = rho(1,i,j-1)/sum(rho(:,i,j-1))
          else
            rho1frac = rho(1,i,j)/sum(rho(:,i,j))
          endif
          yp_vals(1,Y_DIRECTION,i) = rho1frac*flow%bc_data(1,BOUNDARY_YP)
          yp_vals(2,Y_DIRECTION,i) = (1.d0-rho1frac)*flow%bc_data(1,BOUNDARY_YP)
        endif
      end do
    endif
  end subroutine FlowUpdateBCFluxOutletD2

  subroutine FlowUpdateBCFluxOutletD3(flow, walls, rho, xm_vals, xp_vals, &
           ym_vals, yp_vals, zm_vals, zp_vals, dist)
    type(flow_type) flow
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye,dist%info%rgzs:dist%info%rgze):: rho
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye,dist%info%rgzs:dist%info%rgze):: walls
    PetscScalar,dimension(dist%s,flow%ndims,dist%info%ys:dist%info%ye, &
         dist%info%zs:dist%info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(dist%s,flow%ndims,dist%info%xs:dist%info%xe, &
         dist%info%zs:dist%info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(dist%s,flow%ndims,dist%info%xs:dist%info%xe, &
         dist%info%ys:dist%info%ye):: zm_vals, zp_vals

    PetscInt i,j,k
    PetscScalar rho1frac

    ! Enforce the condition that rho1/rho_tot is fixed

    ! xm boundary
    if ((flow%bc_flags(BOUNDARY_XM).eq.BC_PRESSURE_OUTLET).and.(dist%info%xs.eq.1)) then
      i = 1
      do k=dist%info%zs,dist%info%ze
      do j=dist%info%ys,dist%info%ye
        if (walls(i,j,k).eq.0) then
          if (walls(i+1,j,k).eq.0) then
            rho1frac = rho(1,i+1,j,k)/sum(rho(:,i+1,j,k))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          xm_vals(1,X_DIRECTION,j,k) = rho1frac*flow%bc_data(1,BOUNDARY_XM)
          xm_vals(2,X_DIRECTION,j,k) = (1.d0-rho1frac)*flow%bc_data(1,BOUNDARY_XM)
        endif
      end do
      end do
    endif

    ! xp boundary
    if ((flow%bc_flags(BOUNDARY_XP).eq.BC_PRESSURE_OUTLET).and. &
         (dist%info%xe.eq.dist%info%NX)) then
      i = dist%info%NX
      do k=dist%info%zs,dist%info%ze
      do j=dist%info%ys,dist%info%ye
        if (walls(i,j,k).eq.0) then
          if (walls(i-1,j,k).eq.0) then
            rho1frac = rho(1,i-1,j,k)/sum(rho(:,i-1,j,k))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          xp_vals(1,X_DIRECTION,j,k) = rho1frac*flow%bc_data(1,BOUNDARY_XP)
          xp_vals(2,X_DIRECTION,j,k) = (1.d0-rho1frac)*flow%bc_data(1,BOUNDARY_XP)
        end if
      end do
      end do
    endif

    ! ym boundary
    if ((flow%bc_flags(BOUNDARY_YM).eq.BC_PRESSURE_OUTLET).and.(dist%info%ys.eq.1)) then
      j = 1
      do k=dist%info%zs,dist%info%ze
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j,k).eq.0) then
          if (walls(i,j+1,k).eq.0) then
            rho1frac = rho(1,i,j+1,k)/sum(rho(:,i,j+1,k))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          ym_vals(1,Y_DIRECTION,i,k) = rho1frac*flow%bc_data(1,BOUNDARY_YM)
          ym_vals(2,Y_DIRECTION,i,k) = (1.d0-rho1frac)*flow%bc_data(1,BOUNDARY_YM)
        end if
      end do
      end do
    endif

    ! yp boundary
    if ((flow%bc_flags(BOUNDARY_YP).eq.BC_PRESSURE_OUTLET).and. &
         (dist%info%ye.eq.dist%info%NY)) then
      j = dist%info%NY
      do k=dist%info%zs,dist%info%ze
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j,k).eq.0) then
          if (walls(i,j-1,k).eq.0) then
            rho1frac = rho(1,i,j-1,k)/sum(rho(:,i,j-1,k))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          yp_vals(1,Y_DIRECTION,i,k) = rho1frac*flow%bc_data(1,BOUNDARY_YP)
          yp_vals(2,Y_DIRECTION,i,k) = (1.d0-rho1frac)*flow%bc_data(1,BOUNDARY_YP)
        endif
      end do
      end do
    endif

    ! zm boundary
    if ((flow%bc_flags(BOUNDARY_ZM).eq.BC_PRESSURE_OUTLET).and.(dist%info%zs.eq.1)) then
      k = 1
      do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j,k).eq.0) then
          if (walls(i,j,k+1).eq.0) then
            rho1frac = rho(1,i,j,k+1)/sum(rho(:,i,j,k+1))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          zm_vals(1,Z_DIRECTION,i,j) = rho1frac*flow%bc_data(1,BOUNDARY_ZM)
          zm_vals(2,Z_DIRECTION,i,j) = (1.d0-rho1frac)*flow%bc_data(1,BOUNDARY_ZM)
        end if
      end do
      end do
    endif

    ! zp boundary
    if ((flow%bc_flags(BOUNDARY_ZP).eq.BC_PRESSURE_OUTLET).and. &
         (dist%info%ze.eq.dist%info%NZ)) then
      k = dist%info%NZ
      do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j,k).eq.0) then
          if (walls(i,j,k-1).eq.0) then
            rho1frac = rho(1,i,j,k-1)/sum(rho(:,i,j,k-1))
          else
            rho1frac = rho(1,i,j,k)/sum(rho(:,i,j,k))
          end if
          zp_vals(1,Z_DIRECTION,i,j) = rho1frac*flow%bc_data(1,BOUNDARY_ZP)
          zp_vals(2,Z_DIRECTION,i,j) = (1.d0-rho1frac)*flow%bc_data(1,BOUNDARY_ZP)
        endif
      end do
      end do
    endif
  end subroutine FlowUpdateBCFluxOutletD3

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
