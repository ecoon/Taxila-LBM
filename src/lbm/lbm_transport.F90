!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_transport.F90
!!!     version:         
!!!     created:         04 April 2011
!!!       on:            14:35:39 MDT
!!!     last modified:   31 October 2011
!!!       at:            16:04:36 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"

module LBM_Transport_module
  use petsc
  use LBM_Specie_module
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

  type, public:: transport_type
     MPI_Comm comm
     PetscInt ndims
     PetscInt nspecies
     PetscInt nsec_species
     PetscInt nminerals
     type(discretization_type),pointer :: disc
     type(specie_type),pointer,dimension(:):: species
     type(grid_type),pointer:: grid
     type(distribution_type),pointer:: distribution
     type(bc_type),pointer:: bc
     
     PetscBool io_fi
     PetscBool io_rho

     Vec solidmass_g, solidmass
     PetscScalar,pointer:: solidmass_a(:)
     PetscBool io_solidmass
     PetscBool reactive_matrix
     
     PetscScalar,pointer,dimension(:,:,:):: fi_eq
     character(len=MAXWORDLENGTH) name       
  end type transport_type

  public :: TransportCreate, &
       TransportDestroy, &
       TransportSetName, &
       TransportSetFromOptions, &
       TransportSetGrid, &
       TransportSetUp, &
       TransportGetArrays, &
       TransportRestoreArrays, &
       TransportUpdateMoments, &
       TransportStream, &
       TransportReactWithWalls, &
       TransportFiInit, &
       TransportCollision, &
       TransportApplyBCs, &
       TransportUpdateDiagnostics, &
       TransportOutputDiagnostics
contains

  function TransportCreate(comm) result(transport)
    MPI_Comm comm
    type(transport_type),pointer:: transport

    allocate(transport)
    transport%comm = comm
    transport%nspecies = -1
    transport%ndims = -1
    nullify(transport%species)
    nullify(transport%grid)
    transport%disc => DiscretizationCreate(transport%comm)
    transport%distribution => DistributionCreate(transport%comm)
    transport%bc => BCCreate(transport%comm)

    transport%solidmass = 0
    transport%solidmass_g = 0
    nullify(transport%solidmass_a)
    transport%reactive_matrix = PETSC_FALSE

    nullify(transport%fi_eq)

    transport%io_fi = PETSC_TRUE
    transport%io_rho = PETSC_TRUE
    transport%io_solidmass = PETSC_TRUE
  end function TransportCreate

  subroutine TransportDestroy(transport, ierr)
    type(transport_type),pointer:: transport
    PetscErrorCode ierr
    PetscInt lcv

    if (associated(transport%disc)) call DiscretizationDestroy(transport%disc, ierr)
    if (associated(transport%species)) then
       do lcv=1,transport%nspecies
          call SpecieDestroy(transport%species(lcv),ierr)
       end do
    end if
    if (associated(transport%distribution)) then
       call DistributionDestroy(transport%distribution, ierr)
    end if
    if (associated(transport%bc)) call BCDestroy(transport%bc, ierr)

    if (transport%solidmass /= 0) call VecDestroy(transport%solidmass, ierr)
    if (transport%solidmass_g /= 0) call VecDestroy(transport%solidmass_g, ierr)
    if (associated(transport%fi_eq)) deallocate(transport%fi_eq)
  end subroutine TransportDestroy

  subroutine TransportSetName(transport, name) 
    type(transport_type) transport 
    character(len=MAXWORDLENGTH):: name       
    transport%name = name
    call DistributionSetName(transport%distribution, name)
  end subroutine TransportSetName
    
  subroutine TransportSetFromOptions(transport, options, ierr)
    use LBM_Options_module
    type(transport_type) transport
    type(options_type) options
    PetscErrorCode ierr
    PetscInt lcv
    PetscBool help
    PetscBool flag
    PetscInt nmax
    PetscBool bcvalue

    transport%nspecies = options%nspecies
    transport%ndims = options%ndims
    transport%species => SpecieCreate(transport%comm, transport%nspecies)
    transport%reactive_matrix = options%transport_reactive_matrix

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    call DiscretizationSetType(transport%disc, options%transport_disc)
    call DiscretizationSetDerivOrder(transport%disc, options%deriv_order)
    call DiscretizationSetUp(transport%disc)
    
    do lcv=1,transport%nspecies
       call SpecieSetSizes(transport%species(lcv), transport%nspecies, transport%disc%b)
       call SpecieSetID(transport%species(lcv), lcv)
       call SpecieSetFromOptions(transport%species(lcv), options, ierr)
    end do

    ! set up the vectors for holding boundary data
    ! dimension 
    call BCSetSizes(transport%bc, options%ndims*options%nspecies)
    call BCSetFromOptions(transport%bc, options, ierr)

    ! parse all flow boundary conditions
    if (help) call PetscPrintf(options%comm, "-bc_conc_{xyz}{mp}: use concentration bcs\n",&
         ierr)
    if (help) call PetscPrintf(options%comm, "-bc_conc_flux_{xyz}{mp}: use concentration"//&
         " flux bcs\n", ierr)

    ! xm boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix,'-bc_conc_xm', bcvalue, flag, ierr)
    if (bcvalue) transport%bc%flags(BOUNDARY_XM) = BC_DIRICHLET

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix,'-bc_conc_flux_xm', bcvalue, flag, ierr)
    if (bcvalue) transport%bc%flags(BOUNDARY_XM) = BC_NEUMANN

    ! xp boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix,'-bc_conc_xp', bcvalue, flag, ierr)
    if (bcvalue) transport%bc%flags(BOUNDARY_XP) = BC_DIRICHLET

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix,'-bc_conc_flux_xp', bcvalue, flag, ierr)
    if (bcvalue) transport%bc%flags(BOUNDARY_XP) = BC_NEUMANN

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_x', bcvalue, flag, ierr)
    if (bcvalue) then
       transport%bc%flags(BOUNDARY_XM) = BC_PERIODIC
       transport%bc%flags(BOUNDARY_XP) = BC_PERIODIC
    end if

    ! ym boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix,'-bc_conc_ym', bcvalue, flag, ierr)
    if (bcvalue) transport%bc%flags(BOUNDARY_YM) = BC_DIRICHLET

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix,'-bc_conc_flux_ym', bcvalue, flag, ierr)
    if (bcvalue) transport%bc%flags(BOUNDARY_YM) = BC_NEUMANN

    ! yp boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix,'-bc_conc_yp', bcvalue, flag, ierr)
    if (bcvalue) transport%bc%flags(BOUNDARY_YP) = BC_DIRICHLET

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix,'-bc_conc_flux_yp', bcvalue, flag, ierr)
    if (bcvalue) transport%bc%flags(BOUNDARY_YP) = BC_NEUMANN

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_y', bcvalue, flag, ierr)
    if (bcvalue) then
       transport%bc%flags(BOUNDARY_YM) = BC_PERIODIC
       transport%bc%flags(BOUNDARY_YP) = BC_PERIODIC
    end if

    if (options%ndims > 2) then
       ! zm boundary
       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix,'-bc_conc_zm',bcvalue,flag,ierr)
       if (bcvalue) transport%bc%flags(BOUNDARY_ZM) = BC_DIRICHLET

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix,'-bc_conc_flux_zm', bcvalue, flag, ierr)
       if (bcvalue) transport%bc%flags(BOUNDARY_ZM) = BC_NEUMANN

       ! zp boundary
       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix,'-bc_conc_zp',bcvalue,flag,ierr)
       if (bcvalue) transport%bc%flags(BOUNDARY_ZP) = BC_DIRICHLET

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix,'-bc_conc_flux_zp', bcvalue, flag, ierr)
       if (bcvalue) transport%bc%flags(BOUNDARY_ZP) = BC_NEUMANN

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_z', bcvalue, flag, ierr)
       if (bcvalue) then
          transport%bc%flags(BOUNDARY_ZM) = BC_PERIODIC
          transport%bc%flags(BOUNDARY_ZP) = BC_PERIODIC
       end if
    end if
  end subroutine TransportSetFromOptions

  subroutine TransportSetGrid(transport, grid)
    type(transport_type) transport
    type (grid_type),pointer:: grid
    transport%grid => grid
    call BCSetGrid(transport%bc, transport%grid)
  end subroutine TransportSetGrid

  subroutine TransportSetUp(transport) 
    type(transport_type) transport
    PetscInt lcv
    PetscErrorCode ierr
    PetscScalar zero
    zero = 0.d0

    do lcv=1,transport%nspecies
       call DiscretizationSetUpRelax(transport%disc, transport%species(lcv)%relax)
    end do
    call DistributionSetInfo(transport%distribution, transport%grid%info)
    call DistributionSetDiscretization(transport%distribution, transport%disc)
    call DistributionSetSizes(transport%distribution, transport%nspecies)
    call DistributionSetDAs(transport%distribution, transport%grid%da(NSPECIEXBDOF), &
         transport%grid%da(NSPECIEDOF))
    call DistributionSetUp(transport%distribution)
    call BCSetUp(transport%bc)

    ! allocate, initialize workspace
    allocate(transport%fi_eq(1:transport%nspecies, 0:transport%disc%b, &
         1:transport%grid%info%gxyzl))
    transport%fi_eq = 0.d0

    if (transport%reactive_matrix) then
       call DMCreateGlobalVector(transport%grid%da(NSPECIEXBDOF), &
            transport%solidmass_g, ierr)
       call DMCreateLocalVector(transport%grid%da(NSPECIEXBDOF), transport%solidmass, ierr)
       call VecSet(transport%solidmass_g, zero, ierr)
       call VecSet(transport%solidmass, zero, ierr)
       call PetscObjectSetName(transport%solidmass_g, &
            trim(transport%name)//'solidmass', ierr)
    end if
  end subroutine TransportSetUp

  subroutine TransportGetArrays(transport, ierr)
    type(transport_type) transport
    PetscErrorCode ierr

    call DistributionGetArrays(transport%distribution, ierr)
    if (transport%solidmass /= 0) then
       call DMDAVecGetArrayF90(transport%grid%da(NSPECIEXBDOF), transport%solidmass, &
            transport%solidmass_a, ierr)
    end if
    call BCGetArrays(transport%bc, ierr)
  end subroutine TransportGetArrays

  subroutine TransportRestoreArrays(transport, ierr)
    type(transport_type) transport
    PetscErrorCode ierr

    call DistributionRestoreArrays(transport%distribution, ierr)
    if (transport%solidmass /= 0) then
       call DMDAVecRestoreArrayF90(transport%grid%da(NSPECIEXBDOF), transport%solidmass, &
            transport%solidmass_a, ierr)
    end if
    call BCRestoreArrays(transport%bc, ierr)
  end subroutine TransportRestoreArrays

  subroutine TransportUpdateMoments(transport, walls)
    type(transport_type) transport
    type(walls_type) walls

    call DistributionCalcDensity(transport%distribution, walls%walls_a)
    call DistributionCalcFlux(transport%distribution, walls%walls_a)
  end subroutine TransportUpdateMoments

  subroutine TransportFiInit(transport, walls)
    type(transport_type) transport
    type(walls_type) walls
    
    call TransportUpdateFeq(transport, walls%walls_a)
    call TransportFiInit_(transport, transport%distribution%fi_a)
  end subroutine TransportFiInit

  subroutine TransportFiInit_(transport, fi)
    type(transport_type) transport
    PetscScalar,dimension(transport%nspecies, &
         0:transport%disc%b,transport%distribution%info%gxyzl):: fi

    fi = transport%fi_eq
  end subroutine TransportFiInit_

  subroutine TransportUpdateDiagnostics(transport, walls)
    type(transport_type) transport
    type(walls_type) walls
    ! nothing to do
  end subroutine TransportUpdateDiagnostics

  subroutine TransportOutputDiagnostics(transport, io)
    use LBM_IO_module
    type(transport_type) transport
    type(io_type) io
    PetscErrorCode ierr

    if (transport%io_fi) then
       call DMDAVecRestoreArrayF90(transport%grid%da(NSPECIEXBDOF), &
            transport%distribution%fi, transport%distribution%fi_a, ierr)
       call DMLocalToGlobalBegin(transport%grid%da(NSPECIEXBDOF), &
            transport%distribution%fi, INSERT_VALUES, &
            transport%distribution%fi_g, ierr)
       call DMLocalToGlobalEnd(transport%grid%da(NSPECIEXBDOF), &
            transport%distribution%fi, INSERT_VALUES, &
            transport%distribution%fi_g, ierr)
       call DMDAVecGetArrayF90(transport%grid%da(NSPECIEXBDOF), &
            transport%distribution%fi, transport%distribution%fi_a, ierr)
       call IOView(io, transport%distribution%fi_g, 'gi')
    end if
    if (transport%io_rho) then
       call DMDAVecRestoreArrayF90(transport%grid%da(NSPECIEDOF), &
            transport%distribution%rho, transport%distribution%rho_a, ierr)
       call DMLocalToGlobalBegin(transport%grid%da(NSPECIEDOF), &
            transport%distribution%rho, INSERT_VALUES, &
            transport%distribution%rho_g, ierr)
       call DMLocalToGlobalEnd(transport%grid%da(NSPECIEDOF), &
            transport%distribution%rho, INSERT_VALUES, &
            transport%distribution%rho_g, ierr)
       call DMDAVecGetArrayF90(transport%grid%da(NSPECIEDOF), &
            transport%distribution%rho, transport%distribution%rho_a, ierr)
       call IOView(io, transport%distribution%rho_g, 'psi')
    end if
  end subroutine TransportOutputDiagnostics

  subroutine TransportStream(transport)
    type(transport_type) transport
    call DistributionStream(transport%distribution)
  end subroutine TransportStream

  subroutine TransportUpdateFeq(transport, walls)
    use LBM_Logging_module
    type(transport_type) transport
    PetscScalar,dimension(transport%grid%info%rgxyzl):: walls

    PetscInt m
    PetscErrorCode ierr

    call PetscLogEventBegin(logger%event_collision_feq,ierr)
    do m=1,transport%distribution%s
      call DiscretizationEquilf(transport%disc, transport%distribution%rho_a, &
           transport%distribution%flux, walls, transport%fi_eq, m, &
           transport%species(m)%relax, transport%distribution)
    end do
    call PetscLogEventEnd(logger%event_collision_feq,ierr)
  end subroutine TransportUpdateFeq


  subroutine TransportCollision(transport, walls, flow)
    use LBM_Flow_module
    type(transport_type) transport
    type(flow_type) flow
    type(walls_type) walls
    PetscErrorCode ierr

    call TransportUpdateFeq(transport, walls%walls_a)
    call DMDAVecGetArrayF90(flow%grid%da(NFLOWDOF),flow%velt_g,flow%velt_a,ierr)
    select case(transport%ndims)
    case(2)
       call TransportCollisionD2(transport, transport%distribution%fi_a, &
            transport%distribution%rho_a, flow%velt_a, walls%walls_a, transport%fi_eq, &
            transport%distribution)
    case(3)
       call TransportCollisionD3(transport, transport%distribution%fi_a, &
            transport%distribution%rho_a, flow%velt_a, walls%walls_a, transport%fi_eq, &
            transport%distribution)
    end select
    call DMDAVecRestoreArrayF90(flow%grid%da(NFLOWDOF),flow%velt_g,flow%velt_a,ierr)
  end subroutine TransportCollision

  subroutine TransportCollisionD3(transport, fi, rho, u, walls, fi_eq, dist)
    use LBM_Relaxation_module
    type(transport_type) transport
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(transport%nspecies, &
         0:transport%disc%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze)::fi,fi_eq
    PetscScalar,dimension(1:dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: rho
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: walls
    PetscScalar,dimension(1:dist%info%ndims, dist%info%xs:dist%info%xe, &
         dist%info%ys:dist%info%ye, dist%info%zs:dist%info%ze):: u

    PetscInt m,i,j,k

    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j,k).eq.0) then
      do m=1,dist%s
        call RelaxationCollide(transport%species(m)%relax, fi(:,:,i,j,k), &
            fi_eq(:,:,i,j,k), m, dist)
      end do
    end if
    end do
    end do
    end do
  end subroutine TransportCollisionD3

  subroutine TransportCollisionD2(transport, fi, rho, u, walls, fi_eq, dist)
    use LBM_Relaxation_module
    type(transport_type) transport
    type(distribution_type) dist ! just for convenience
    PetscScalar,dimension(transport%nspecies, 0:transport%disc%b, &
         dist%info%gxs:dist%info%gxe, dist%info%gys:dist%info%gye)::fi,fi_eq
    PetscScalar,dimension(1:dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe,&
         dist%info%rgys:dist%info%rgye):: walls
    PetscScalar,dimension(1:dist%info%ndims, &
         dist%info%xs:dist%info%xe, dist%info%ys:dist%info%ye):: u

    PetscInt m,i,j

    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j) > 0) then
      do m=1,dist%s
        call RelaxationCollide(transport%species(m)%relax, fi(:,:,i,j), fi_eq(:,:,i,j), &
             m, dist)
      end do
    end if
    end do
    end do
  end subroutine TransportCollisionD2

  ! bounceback if nonreactive walls, otherwise call reaction algorithm
  subroutine TransportReactWithWalls(transport, walls)
    type(transport_type) transport
    type(walls_type) walls
    
    if (transport%reactive_matrix) then
       select case(transport%ndims)
       case(2)
          call TransportReactWithWallsD2(transport, transport%distribution%fi_a, &
               transport%distribution%rho_a, walls%walls_a, transport%distribution)
       case(3)
          call TransportReactWithWallsD3(transport, transport%distribution%fi_a, &
               transport%distribution%rho_a, walls%walls_a, transport%distribution)
       end select
    else 
       call DistributionBounceback(transport%distribution, walls%walls_a)
    end if
  end subroutine TransportReactWithWalls

  subroutine TransportReactWithWallsD3(transport, gi, psi, walls, dist)
    type(transport_type) transport
    type(distribution_type) dist ! for convenience
    PetscScalar,dimension(transport%nspecies, &
         0:transport%disc%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze)::gi
    PetscScalar,dimension(1:dist%s,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: psi
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: walls

    PetscScalar gi_opp(transport%nspecies)

    PetscInt i,j,k,m,n
    do k=transport%distribution%info%zs,transport%distribution%info%ze
    do j=transport%distribution%info%ys,transport%distribution%info%ye
    do i=transport%distribution%info%xs,transport%distribution%info%xe
    if (walls(i,j,k) > 0) then
       do n=1,dist%b
          if (walls(i+transport%disc%ci(n,X_DIRECTION), j+transport%disc%ci(n,Y_DIRECTION),&
               k+transport%disc%ci(n,Z_DIRECTION)).eq.0) then
             ! must specify gi(b), as the node in that direction is fluid
             gi_opp = gi(:,transport%disc%opposites(n),i,j,k)
          end if
       end do
    end if
    end do
    end do
    end do
  end subroutine TransportReactWithWallsD3

  subroutine TransportReactWithWallsD2(transport, gi, psi, walls, dist)
    type(transport_type) transport
    type(distribution_type) dist ! for convenience
    PetscScalar,dimension(transport%nspecies, &
         0:transport%disc%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye)::gi
    PetscScalar,dimension(1:dist%s,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: psi
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls

    PetscScalar gi_opp(transport%nspecies)

    PetscInt i,j,m,n
    do j=transport%distribution%info%ys,transport%distribution%info%ye
    do i=transport%distribution%info%xs,transport%distribution%info%xe
    if (walls(i,j) > 0) then
       do n=1,dist%b
          if (walls(i+transport%disc%ci(n,X_DIRECTION), &
               j+transport%disc%ci(n,Y_DIRECTION)).eq.0) then
             ! must specify gi(b), as the node in that direction is fluid
             gi_opp = gi(:,transport%disc%opposites(n),i,j)
          end if
       end do
    end if
    end do
    end do
  end subroutine TransportReactWithWallsD2

  subroutine TransportApplyBCs(transport, walls)
    type(transport_type) transport
    type(walls_type) walls
    PetscScalar,dimension(transport%distribution%s, &
         transport%distribution%disc%ndims, &
         transport%distribution%info%gxyzl) :: zeros
    call BCApply(transport%bc, zeros, walls%walls_a, transport%distribution)
  end subroutine TransportApplyBCs
end module LBM_Transport_module

