!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        bc.F90
!!!     version:         
!!!     created:         06 December 2010
!!!       on:            09:03:18 MST
!!!     last modified:   05 August 2011
!!!       at:            11:05:11 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ ldeo.columbia.edu
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module LBM_BC_module
  use LBM_Info_module
  use LBM_Grid_module
  use LBM_Distribution_Function_type_module
  use LBM_Discretization_module
  use LBM_Options_module
  use petsc
  implicit none
  private
#include "lbm_definitions.h"

  type, public:: bc_type
     MPI_Comm comm
     type(grid_type),pointer:: grid
     PetscInt nbcs ! number of boundary conditions per face
     PetscInt, dimension(6):: flags ! enum for boundary conditions
     PetscInt xml,xpl,yml,ypl,zml,zpl ! sizes
     PetscScalar,pointer,dimension(:):: xm_a, xp_a
     PetscScalar,pointer,dimension(:):: ym_a, yp_a
     PetscScalar,pointer,dimension(:):: zm_a, zp_a
     Vec xm
     Vec xp
     Vec ym
     Vec yp
     Vec zm
     Vec zp
  end type bc_type

  public :: BCCreate, &
       BCDestroy, &
       BCSetGrid, &
       BCSetSizes, &
       BCSetFromOptions, &
       BCSetUp, &
       BCSetValues, &
       BCZeroForces, &
       BCGetArrays, &
       BCRestoreArrays, &
       BCApply, &
       BCApplyPseudoperiodic, &
       BCApplyDirichlet, &
       BCApplyFlux, &
       BCApplyVelocity, &
       BCApplyZeroGradient
contains

  ! constructor
  function BCCreate(comm) result(bc)
    MPI_Comm comm
    type(bc_type),pointer :: bc
    allocate(bc)
    bc%comm = comm
    nullify(bc%grid)
    nullify(bc%xm_a)
    nullify(bc%xp_a)
    nullify(bc%ym_a)
    nullify(bc%yp_a)
    nullify(bc%zm_a)
    nullify(bc%zp_a)
    bc%flags = 0
    bc%nbcs = 0
    
    bc%xm = 0
    bc%xp = 0
    bc%ym = 0
    bc%yp = 0
    bc%zm = 0
    bc%zp = 0

    bc%xml = 0
    bc%xpl = 0
    bc%yml = 0
    bc%ypl = 0
    bc%zml = 0
    bc%zpl = 0
  end function BCCreate
  
  ! destructor
  subroutine BCDestroy(bc, ierr)
    type(bc_type) bc
    PetscErrorCode ierr
    if (bc%xm /= 0) call VecDestroy(bc%xm, ierr)
    if (bc%xp /= 0) call VecDestroy(bc%xp, ierr)
    if (bc%ym /= 0) call VecDestroy(bc%ym, ierr)
    if (bc%yp /= 0) call VecDestroy(bc%yp, ierr)
    if (bc%zm /= 0) call VecDestroy(bc%zm, ierr)
    if (bc%zp /= 0) call VecDestroy(bc%zp, ierr)
  end subroutine BCDestroy
  
  subroutine BCSetGrid(bc, grid)
    type(bc_type) bc
    type(grid_type),pointer:: grid
    bc%grid => grid
  end subroutine BCSetGrid

  subroutine BCSetSizes(bc, nbcs)
    type(bc_type) bc
    PetscInt nbcs
    bc%nbcs = nbcs
  end subroutine BCSetSizes

  ! set up the vectors for holding boundary data
  subroutine BCSetFromOptions(bc, options, ierr)
    type(bc_type) bc
    type(options_type) options
    PetscErrorCode ierr

    ! locals
    PetscBool flag
    PetscInt nmax

    ! flags and constant values from options
    nmax = 6
    call PetscOptionsGetIntArray(options%my_prefix, '-bc_flags', bc%flags, &
         nmax, flag, ierr)
  end subroutine BCSetFromOptions
 
  subroutine BCSetUp(bc)
    type(bc_type) bc
    PetscInt locn
    PetscErrorCode ierr

    type(info_type),pointer:: info 
    info => bc%grid%info

    ! now make the boundary vecs/arrays
    ! x boundaries
    if (info%xs.eq.1) then
       bc%xml = info%yl*info%zl 
    else
       bc%xml = 0
    endif
    locn = bc%xml*bc%nbcs
    call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%xm, ierr)
    call VecSetBlockSize(bc%xm, bc%nbcs, ierr)
    call PetscObjectSetName(bc%xm, 'xm_bc', ierr)
    
    if (info%xe.eq.info%NX) then
       bc%xpl = info%yl*info%zl 
    else
       bc%xpl = 0
    endif
    locn = bc%xpl*bc%nbcs
    call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%xp, ierr)
    call VecSetBlockSize(bc%xp, bc%nbcs, ierr)
    call PetscObjectSetName(bc%xp, 'xp_bc', ierr)
    
    ! y boundaries
    if (info%ys.eq.1) then
       bc%yml = info%xl*info%zl
    else
       bc%yml = 0
    endif
    locn = bc%yml*bc%nbcs
    call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%ym, ierr)
    call VecSetBlockSize(bc%ym, bc%nbcs, ierr)
    call PetscObjectSetName(bc%ym, 'ym_bc', ierr)
    
    if (info%ye.eq.info%NY) then
       bc%ypl = info%xl*info%zl
    else
       bc%ypl = 0
    endif
    locn = bc%ypl*bc%nbcs
    call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%yp, ierr)
    call VecSetBlockSize(bc%yp, bc%nbcs, ierr)
    call PetscObjectSetName(bc%yp, 'yp_bc', ierr)
    
    if ((info%ndims > 2).and.(info%zs.eq.1)) then
       bc%zml = info%xl*info%yl
    else
       bc%zml = 0
    endif
    locn = bc%zml*bc%nbcs
    call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%zm, ierr)
    call VecSetBlockSize(bc%zm, bc%nbcs, ierr)
    call PetscObjectSetName(bc%zm, 'zm_bc', ierr)
       
    if ((info%ndims > 2).and.(info%ze.eq.info%NZ)) then
       bc%zpl = info%xl*info%yl
    else
       bc%zpl = 0
    endif
    locn = bc%zpl*bc%nbcs
    call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%zp, ierr)
    call VecSetBlockSize(bc%zp, bc%nbcs, ierr)
    call PetscObjectSetName(bc%zp, 'zp_bc', ierr)
 
    call BCGetArrays(bc, ierr)
  end subroutine BCSetUp

  ! call initialize
  subroutine BCSetValues(bc, dist, options, bc_subroutine)
    type(bc_type) bc
    type(distribution_type) dist
    type(options_type) options
    external bc_subroutine
    
    PetscErrorCode ierr
    call BCGetArrays(bc, ierr)
    call bc_subroutine(bc%flags, bc%xm_a, bc%xp_a, &
         bc%ym_a, bc%yp_a, bc%zm_a, bc%zp_a, bc%nbcs, dist, options)
    CHKERRQ(ierr)
    CHKMEMQ
    call BCRestoreArrays(bc, ierr)
  end subroutine BCSetValues
    
  subroutine BCGetArrays(bc, ierr)
    type(bc_type) bc
    PetscErrorCode ierr
    call VecGetArrayF90(bc%xm, bc%xm_a, ierr)
    call VecGetArrayF90(bc%xp, bc%xp_a, ierr)
    call VecGetArrayF90(bc%ym, bc%ym_a, ierr)
    call VecGetArrayF90(bc%yp, bc%yp_a, ierr)
    if (bc%zm /= 0) call VecGetArrayF90(bc%zm, bc%zm_a, ierr)
    if (bc%zp /= 0) call VecGetArrayF90(bc%zp, bc%zp_a, ierr)
  end subroutine BCGetArrays

  subroutine BCRestoreArrays(bc, ierr)
    type(bc_type) bc
    PetscErrorCode ierr
    call VecRestoreArrayF90(bc%xm, bc%xm_a, ierr)
    call VecRestoreArrayF90(bc%xp, bc%xp_a, ierr)
    call VecRestoreArrayF90(bc%ym, bc%ym_a, ierr)
    call VecRestoreArrayF90(bc%yp, bc%yp_a, ierr)
    if (bc%zm /= 0) call VecRestoreArrayF90(bc%zm, bc%zm_a, ierr)
    if (bc%zp /= 0) call VecRestoreArrayF90(bc%zp, bc%zp_a, ierr)
  end subroutine BCRestoreArrays
  
  subroutine BCZeroForces(bc, forces, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 1:dist%info%ndims, 1:dist%info%gxyzl):: forces
    select case(dist%info%ndims)
    case(2)
       call BCZeroForcesD2(bc, forces, dist)
    case(3)
       call BCZeroForcesD3(bc, forces, dist)
    end select
  end subroutine BCZeroForces

  subroutine BCZeroForcesD3(bc, forces, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: forces

    ! -- x
    if ((dist%info%xs.eq.1).and.((bc%flags(BOUNDARY_XM).eq.BC_FLUX).or.&
         (bc%flags(BOUNDARY_XM).eq.BC_VELOCITY).or.&
         (bc%flags(BOUNDARY_XM).eq.BC_DIRICHLET))) then
       forces(:,:,1,:,:) = 0
    endif

    if ((dist%info%xe.eq.dist%info%NX).and.((bc%flags(BOUNDARY_XP).eq.BC_FLUX).or. &
         (bc%flags(BOUNDARY_XP).eq.BC_VELOCITY).or.&
         (bc%flags(BOUNDARY_XP).eq.BC_DIRICHLET))) then
       forces(:,:,dist%info%NX,:,:) = 0
    endif

    ! -- y
    if ((dist%info%ys.eq.1).and.((bc%flags(BOUNDARY_YM).eq.BC_FLUX).or.&
         (bc%flags(BOUNDARY_YM).eq.BC_VELOCITY).or.&
         (bc%flags(BOUNDARY_YM).eq.BC_DIRICHLET))) then
       forces(:,:,:,1,:) = 0
    endif

    if ((dist%info%ye.eq.dist%info%NY).and.((bc%flags(BOUNDARY_YP).eq.BC_FLUX).or. &
         (bc%flags(BOUNDARY_YP).eq.BC_VELOCITY).or.&
         (bc%flags(BOUNDARY_YP).eq.BC_DIRICHLET))) then
       forces(:,:,:,dist%info%NY,:) = 0
    endif

    ! -- z
    if ((dist%info%zs.eq.1).and.((bc%flags(BOUNDARY_ZM).eq.BC_FLUX).or.&
         (bc%flags(BOUNDARY_ZM).eq.BC_VELOCITY).or.&
         (bc%flags(BOUNDARY_ZM).eq.BC_DIRICHLET))) then
       forces(:,:,:,:,1) = 0
    endif

    if ((dist%info%ze.eq.dist%info%NZ).and.((bc%flags(BOUNDARY_ZP).eq.BC_FLUX).or. &
         (bc%flags(BOUNDARY_ZP).eq.BC_VELOCITY).or.&
         (bc%flags(BOUNDARY_ZP).eq.BC_DIRICHLET))) then
       forces(:,:,:,:,dist%info%NZ) = 0
    endif

  end subroutine BCZeroForcesD3

  subroutine BCZeroForcesD2(bc, forces, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: forces

    ! -- x
    if ((dist%info%xs.eq.1).and.((bc%flags(BOUNDARY_XM).eq.BC_FLUX).or.&
         (bc%flags(BOUNDARY_XM).eq.BC_VELOCITY).or.&
         (bc%flags(BOUNDARY_XM).eq.BC_DIRICHLET))) then
       forces(:,:,1,:) = 0
    endif

    if ((dist%info%xe.eq.dist%info%NX).and.((bc%flags(BOUNDARY_XP).eq.BC_FLUX).or. &
         (bc%flags(BOUNDARY_XP).eq.BC_VELOCITY).or.&
         (bc%flags(BOUNDARY_XP).eq.BC_DIRICHLET))) then
       forces(:,:,dist%info%NX,:) = 0
    endif

    ! -- y
    if ((dist%info%ys.eq.1).and.((bc%flags(BOUNDARY_YM).eq.BC_FLUX).or.&
         (bc%flags(BOUNDARY_YM).eq.BC_VELOCITY).or.&
         (bc%flags(BOUNDARY_YM).eq.BC_DIRICHLET))) then
       forces(:,:,:,1) = 0
    endif

    if ((dist%info%ye.eq.dist%info%NY).and.((bc%flags(BOUNDARY_YP).eq.BC_FLUX).or. &
         (bc%flags(BOUNDARY_YP).eq.BC_VELOCITY).or.&
         (bc%flags(BOUNDARY_YP).eq.BC_DIRICHLET))) then
       forces(:,:,:,dist%info%NY) = 0
    endif
  end subroutine BCZeroForcesD2

  subroutine BCApply(bc, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%info%gxyzl):: walls

    logical,dimension(0:10):: bcs_done    
    PetscInt lcv_sides
    bcs_done=.FALSE.
    bcs_done(BC_PERIODIC) = .TRUE.   ! periodic done by default

    do lcv_sides = 1,6
       if (.not.bcs_done(bc%flags(lcv_sides))) then
          select case (bc%flags(lcv_sides))
          case (BC_PSEUDOPERIODIC)         ! pseudo-periodic
             call BCApplyPseudoperiodic(bc, walls, dist)
          case (BC_FLUX)         ! flux
             call BCApplyFlux(bc, walls, dist)
          case (BC_VELOCITY)         ! flux
             call BCApplyVelocity(bc, walls, dist)
          case (BC_DIRICHLET)         ! dirichlet conc/pressure
             call BCApplyDirichlet(bc, walls, dist)
          case (BC_ZERO_GRADIENT)         ! dirichlet conc/pressure
             call BCApplyZeroGradient(bc, walls, dist)
          end select
          bcs_done(bc%flags(lcv_sides)) = .TRUE. ! only do each bc type once
       endif
    enddo
  end subroutine BCApply

  subroutine BCApplyPseudoperiodic(bc, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%info%gxyzl):: walls
    PetscErrorCode ierr

    ! first, check to make sure number of phases = 2
    if (dist%s /= 2) then
       SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic only makes sense for two-phase', ierr)
       return
    end if

    select case(dist%info%ndims)
    case(3)
       call BCApplyPseudoperiodicD3(bc, dist%fi_a, walls, dist)
    case(D2Q9_DISCRETIZATION)
       call BCApplyPseudoperiodicD2(bc, dist%fi_a, walls, dist)
    case DEFAULT
       SETERRQ(PETSC_COMM_SELF, 1, 'invalid discretization in LBM', ierr)
    end select
  end subroutine BCApplyPseudoperiodic

  subroutine BCApplyDirichlet(bc, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%info%gxyzl):: walls

    select case(dist%info%ndims)
    case(2)
       call BCApplyDirichletD2(bc, dist%fi_a, walls, bc%xm_a, bc%xp_a, &
            bc%ym_a, bc%yp_a, dist)
    case(3)
       call BCApplyDirichletD3(bc, dist%fi_a, walls, bc%xm_a, bc%xp_a, &
            bc%ym_a, bc%yp_a, bc%zm_a, bc%zp_a, dist)
    end select
  end subroutine BCApplyDirichlet

  subroutine BCApplyFlux(bc, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%info%gxyzl):: walls

    select case(dist%info%ndims)
    case(2)
       call BCApplyFluxD2(bc, dist%fi_a, walls, bc%xm_a, bc%xp_a, &
            bc%ym_a, bc%yp_a, dist)
    case(3)
       call BCApplyFluxD3(bc, dist%fi_a, walls, bc%xm_a, bc%xp_a, &
            bc%ym_a, bc%yp_a, bc%zm_a, bc%zp_a, dist)
    end select
  end subroutine BCApplyFlux

  subroutine BCApplyVelocity(bc, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%info%gxyzl):: walls

    select case(dist%info%ndims)
    case(2)
       call BCApplyVelocityD2(bc, dist%fi_a, walls, bc%xm_a, bc%xp_a, &
            bc%ym_a, bc%yp_a, dist)
    case(3)
       call BCApplyVelocityD3(bc, dist%fi_a, walls, bc%xm_a, bc%xp_a, &
            bc%ym_a, bc%yp_a, bc%zm_a, bc%zp_a, dist)
    end select
  end subroutine BCApplyVelocity

  subroutine BCApplyZeroGradient(bc, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%info%gxyzl):: walls
    PetscErrorCode ierr

    SETERRQ(PETSC_COMM_SELF, 1, 'zero gradient BC not implemented yet', ierr)
    return
    select case(dist%info%ndims)
    case(2)
       call BCApplyZeroGradientD2(bc, dist%fi_a, walls, dist)
    case(3)
       call BCApplyZeroGradientD3(bc, dist%fi_a, walls, dist)
    end select
  end subroutine BCApplyZeroGradient

  subroutine BCApplyDirichletD3(bc, fi, walls, xm_vals, xp_vals, &
       ym_vals, yp_vals, zm_vals, zp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls
    PetscScalar,dimension(bc%nbcs,dist%info%ys:dist%info%ye, &
         dist%info%zs:dist%info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe, &
         dist%info%zs:dist%info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe, &
         dist%info%ys:dist%info%ye):: zm_vals, zp_vals

    PetscInt i,j,k
    PetscInt directions(0:dist%b)
    PetscInt cardinals(1:dist%info%ndims)

    directions(:) = 0
    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_DIRICHLET).and.(dist%info%xs.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XM, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
             i = 1
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCDirichletToBoundary(dist%disc, fi(:,:,i,j,k), &
                     xm_vals(:,j,k), directions, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_DIRICHLET).and.(dist%info%xe.eq.dist%info%NX)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XP, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
             i = dist%info%NX
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCDirichletToBoundary(dist%disc, fi(:,:,i,j,k), &
                     xp_vals(:,j,k), directions, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_DIRICHLET).and.(dist%info%ys.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YM, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do i=dist%info%xs,dist%info%xe
             j = 1
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCDirichletToBoundary(dist%disc, fi(:,:,i,j,k), &
                     ym_vals(:,i,k), directions, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_DIRICHLET).and.(dist%info%ye.eq.dist%info%NY)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YP, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do i=dist%info%xs,dist%info%xe
             j = dist%info%NY
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCDirichletToBoundary(dist%disc, fi(:,:,i,j,k), &
                     yp_vals(:,i,k), directions, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! ZM BOUNDARY
    if ((bc%flags(BOUNDARY_ZM).eq.BC_DIRICHLET).and.(dist%info%zs.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_ZM, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
             k = 1
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCDirichletToBoundary(dist%disc, fi(:,:,i,j,k), &
                     zm_vals(:,i,j), directions, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! ZP BOUNDARY
    if ((bc%flags(BOUNDARY_ZP).eq.BC_DIRICHLET).and.(dist%info%ze.eq.dist%info%NZ)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_ZP, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
             k = dist%info%NZ
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCDirichletToBoundary(dist%disc, fi(:,:,i,j,k), &
                     zp_vals(:,i,j), directions, dist)
             end if
          end do
       end do
    endif
    return
  end subroutine BCApplyDirichletD3

  subroutine BCApplyDirichletD2(bc, fi, walls, xm_vals, xp_vals, ym_vals, yp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls
    PetscScalar,dimension(bc%nbcs,dist%info%ys:dist%info%ye):: xm_vals, xp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe):: ym_vals, yp_vals

    PetscInt i,j
    PetscInt directions(0:dist%b)
    PetscInt cardinals(1:dist%info%ndims)
    directions(:) = 0

    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_DIRICHLET).and.(dist%info%xs.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XM, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          i = 1
          if (walls(i,j).eq.0) then
             call DiscApplyBCDirichletToBoundary(dist%disc, fi(:,:,i,j), &
                  xm_vals(:,j), directions, dist)
          end if
       end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_DIRICHLET).and.(dist%info%xe.eq.dist%info%NX)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XP, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          i = dist%info%NX
          if (walls(i,j).eq.0) then
             call DiscApplyBCDirichletToBoundary(dist%disc, fi(:,:,i,j), &
                  xp_vals(:,j), directions, dist)
          end if
       end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_DIRICHLET).and.(dist%info%ys.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YM, directions, cardinals)
       do i=dist%info%xs,dist%info%xe
          j = 1
          if (walls(i,j).eq.0) then
             call DiscApplyBCDirichletToBoundary(dist%disc, fi(:,:,i,j), &
                  ym_vals(:,i), directions, dist)
          end if
       end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_DIRICHLET).and.(dist%info%ye.eq.dist%info%NY)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YP, directions, cardinals)
       do i=dist%info%xs,dist%info%xe
          j = dist%info%NY
          if (walls(i,j).eq.0) then
             call DiscApplyBCDirichletToBoundary(dist%disc, fi(:,:,i,j), &
                  yp_vals(:,i), directions, dist)
          end if
       end do
    endif
  end subroutine BCApplyDirichletD2

  subroutine BCApplyFluxD3(bc, fi, walls, xm_vals, xp_vals, &
       ym_vals, yp_vals, zm_vals, zp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls
    PetscScalar,dimension(bc%nbcs,dist%info%ys:dist%info%ye, &
         dist%info%zs:dist%info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe, &
         dist%info%zs:dist%info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe, &
         dist%info%ys:dist%info%ye):: zm_vals, zp_vals

    PetscInt i,j,k
    PetscInt directions(0:dist%b)
    PetscInt cardinals(1:dist%info%ndims)

    directions(:) = 0
    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_FLUX).and.(dist%info%xs.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XM, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
             i = 1
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCFluxToBoundary(dist%disc, fi(:,:,i,j,k), &
                     xm_vals(:,j,k), directions, cardinals, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_FLUX).and.(dist%info%xe.eq.dist%info%NX)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XP, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
             i = dist%info%NX
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCFluxToBoundary(dist%disc, fi(:,:,i,j,k), &
                     xp_vals(:,j,k), directions, cardinals, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_FLUX).and.(dist%info%ys.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YM, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do i=dist%info%xs,dist%info%xe
             j = 1
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCFluxToBoundary(dist%disc, fi(:,:,i,j,k), &
                     ym_vals(:,i,k), directions, cardinals, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_FLUX).and.(dist%info%ye.eq.dist%info%NY)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YP, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do i=dist%info%xs,dist%info%xe
             j = dist%info%NY
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCFluxToBoundary(dist%disc, fi(:,:,i,j,k), &
                     yp_vals(:,i,k), directions, cardinals, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! ZM BOUNDARY
    if ((bc%flags(BOUNDARY_ZM).eq.BC_FLUX).and.(dist%info%zs.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_ZM, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
             k = 1
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCFluxToBoundary(dist%disc, fi(:,:,i,j,k), &
                     zm_vals(:,i,j), directions, cardinals, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! ZP BOUNDARY
    if ((bc%flags(BOUNDARY_ZP).eq.BC_FLUX).and.(dist%info%ze.eq.dist%info%NZ)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_ZP, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
             k = dist%info%NZ
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCFluxToBoundary(dist%disc, fi(:,:,i,j,k), &
                     zp_vals(:,i,j), directions, cardinals, dist)
             end if
          end do
       end do
    endif
    return
  end subroutine BCApplyFluxD3

  subroutine BCApplyFluxD2(bc, fi, walls, xm_vals, xp_vals, ym_vals, yp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls
    PetscScalar,dimension(bc%nbcs,dist%info%ys:dist%info%ye):: xm_vals, xp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe):: ym_vals, yp_vals

    PetscInt i,j
    PetscInt directions(0:dist%b)
    PetscInt cardinals(1:dist%info%ndims)
    directions(:) = 0

    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_FLUX).and.(dist%info%xs.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XM, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          i = 1
          if (walls(i,j).eq.0) then
             call DiscApplyBCFluxToBoundary(dist%disc, fi(:,:,i,j), &
                  xm_vals(:,j), directions, cardinals, dist)
          end if
       end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_FLUX).and.(dist%info%xe.eq.dist%info%NX)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XP, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          i = dist%info%NX
          if (walls(i,j).eq.0) then
             call DiscApplyBCFluxToBoundary(dist%disc, fi(:,:,i,j), &
                  xp_vals(:,j), directions, cardinals, dist)
          end if
       end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_FLUX).and.(dist%info%ys.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YM, directions, cardinals)
       do i=dist%info%xs,dist%info%xe
          j = 1
          if (walls(i,j).eq.0) then
             call DiscApplyBCFluxToBoundary(dist%disc, fi(:,:,i,j), &
                  ym_vals(:,i), directions, cardinals, dist)
          end if
       end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_FLUX).and.(dist%info%ye.eq.dist%info%NY)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YP, directions, cardinals)
       do i=dist%info%xs,dist%info%xe
          j = dist%info%NY
          if (walls(i,j).eq.0) then
             call DiscApplyBCFluxToBoundary(dist%disc, fi(:,:,i,j), &
                  yp_vals(:,i), directions, cardinals, dist)
          end if
       end do
    endif
  end subroutine BCApplyFluxD2

  subroutine BCApplyVelocityD3(bc, fi, walls, xm_vals, xp_vals, &
       ym_vals, yp_vals, zm_vals, zp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls
    PetscScalar,dimension(bc%nbcs,dist%info%ys:dist%info%ye, &
         dist%info%zs:dist%info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe, &
         dist%info%zs:dist%info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe, &
         dist%info%ys:dist%info%ye):: zm_vals, zp_vals

    PetscInt i,j,k
    PetscInt directions(0:dist%b)
    PetscInt cardinals(1:dist%info%ndims)

    directions(:) = 0
    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_VELOCITY).and.(dist%info%xs.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XM, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
             i = 1
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCVelocityToBoundary(dist%disc, fi(:,:,i,j,k), &
                     xm_vals(:,j,k), directions, cardinals, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_VELOCITY).and.(dist%info%xe.eq.dist%info%NX)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XP, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
             i = dist%info%NX
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCVelocityToBoundary(dist%disc, fi(:,:,i,j,k), &
                     xp_vals(:,j,k), directions, cardinals, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_VELOCITY).and.(dist%info%ys.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YM, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do i=dist%info%xs,dist%info%xe
             j = 1
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCVelocityToBoundary(dist%disc, fi(:,:,i,j,k), &
                     ym_vals(:,i,k), directions, cardinals, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_VELOCITY).and.(dist%info%ye.eq.dist%info%NY)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YP, directions, cardinals)
       do k=dist%info%zs,dist%info%ze
          do i=dist%info%xs,dist%info%xe
             j = dist%info%NY
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCVelocityToBoundary(dist%disc, fi(:,:,i,j,k), &
                     yp_vals(:,i,k), directions, cardinals, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! ZM BOUNDARY
    if ((bc%flags(BOUNDARY_ZM).eq.BC_VELOCITY).and.(dist%info%zs.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_ZM, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
             k = 1
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCVelocityToBoundary(dist%disc, fi(:,:,i,j,k), &
                     zm_vals(:,i,j), directions, cardinals, dist)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! ZP BOUNDARY
    if ((bc%flags(BOUNDARY_ZP).eq.BC_VELOCITY).and.(dist%info%ze.eq.dist%info%NZ)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_ZP, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
             k = dist%info%NZ
             if (walls(i,j,k).eq.0) then
                call DiscApplyBCVelocityToBoundary(dist%disc, fi(:,:,i,j,k), &
                     zp_vals(:,i,j), directions, cardinals, dist)
             end if
          end do
       end do
    endif
    return
  end subroutine BCApplyVelocityD3

  subroutine BCApplyVelocityD2(bc, fi, walls, xm_vals, xp_vals, ym_vals, yp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls
    PetscScalar,dimension(bc%nbcs,dist%info%ys:dist%info%ye):: xm_vals, xp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe):: ym_vals, yp_vals

    PetscInt i,j
    PetscInt directions(0:dist%b)
    PetscInt cardinals(1:dist%info%ndims)
    directions(:) = 0

    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_VELOCITY).and.(dist%info%xs.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XM, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          i = 1
          if (walls(i,j).eq.0) then
             call DiscApplyBCVelocityToBoundary(dist%disc, fi(:,:,i,j), &
                  xm_vals(:,j), directions, cardinals, dist)
          end if
       end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_VELOCITY).and.(dist%info%xe.eq.dist%info%NX)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_XP, directions, cardinals)
       do j=dist%info%ys,dist%info%ye
          i = dist%info%NX
          if (walls(i,j).eq.0) then
             call DiscApplyBCVelocityToBoundary(dist%disc, fi(:,:,i,j), &
                  xp_vals(:,j), directions, cardinals, dist)
          end if
       end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_VELOCITY).and.(dist%info%ys.eq.1)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YM, directions, cardinals)
       do i=dist%info%xs,dist%info%xe
          j = 1
          if (walls(i,j).eq.0) then
             call DiscApplyBCVelocityToBoundary(dist%disc, fi(:,:,i,j), &
                  ym_vals(:,i), directions, cardinals, dist)
          end if
       end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_VELOCITY).and.(dist%info%ye.eq.dist%info%NY)) then
       call DiscSetLocalDirections(dist%disc, BOUNDARY_YP, directions, cardinals)
       do i=dist%info%xs,dist%info%xe
          j = dist%info%NY
          if (walls(i,j).eq.0) then
             call DiscApplyBCVelocityToBoundary(dist%disc, fi(:,:,i,j), &
                  yp_vals(:,i), directions, cardinals, dist)
          end if
       end do
    endif
  end subroutine BCApplyVelocityD2

  subroutine BCApplyPseudoperiodicD3(bc, fi, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls
    PetscErrorCode ierr

    PetscInt i,j,k,n
    PetscScalar tmp

    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_PSEUDOPERIODIC).and.(dist%info%xs.eq.1)) then
       if (.not.dist%info%periodic(X_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_x', ierr)
          return
       end if
       do k=dist%info%zs,dist%info%ze
       do j=dist%info%ys,dist%info%ye
          i = 1
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,X_DIRECTION) > 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if

    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_PSEUDOPERIODIC).and.(dist%info%xe.eq.dist%info%NX)) &
         then
       if (.not.dist%info%periodic(X_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_x', ierr)
          return
       end if
       do k=dist%info%zs,dist%info%ze
       do j=dist%info%ys,dist%info%ye
          i = dist%info%NX
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,X_DIRECTION) < 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if

    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_PSEUDOPERIODIC).and.(dist%info%ys.eq.1)) then
       if (.not.dist%info%periodic(Y_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_y', ierr)
          return
       end if
       do k=dist%info%zs,dist%info%ze
       do i=dist%info%xs,dist%info%xe
          j = 1
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Y_DIRECTION) > 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if

    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_PSEUDOPERIODIC).and.(dist%info%ye.eq.dist%info%NY)) &
         then
       if (.not.dist%info%periodic(Y_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_y', ierr)
          return
       end if
       do k=dist%info%zs,dist%info%ze
       do i=dist%info%xs,dist%info%xe
          j = dist%info%NY
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Y_DIRECTION) < 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if

    ! ZM BOUNDARY
    if ((bc%flags(BOUNDARY_ZM).eq.BC_PSEUDOPERIODIC).and.(dist%info%zs.eq.1)) then
       if (.not.dist%info%periodic(Z_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_z', ierr)
          return
       end if
       do j=dist%info%ys,dist%info%ye
       do i=dist%info%xs,dist%info%xe
          k = 1
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Z_DIRECTION) > 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if

    ! ZM BOUNDARY
    if ((bc%flags(BOUNDARY_ZP).eq.BC_PSEUDOPERIODIC).and.(dist%info%ze.eq.dist%info%NZ)) &
         then
       if (.not.dist%info%periodic(Z_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_z', ierr)
          return
       end if
       do j=dist%info%ys,dist%info%ye
       do i=dist%info%xs,dist%info%xe
          k = dist%info%NZ
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Z_DIRECTION) < 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if
  end subroutine BCApplyPseudoperiodicD3

  subroutine BCApplyPseudoperiodicD2(bc, fi, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls
    PetscErrorCode ierr

    PetscInt i,j,n
    PetscScalar tmp

    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_PSEUDOPERIODIC).and.(dist%info%xs.eq.1)) then
       if (.not.dist%info%periodic(X_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_x', ierr)
          return
       end if
       do j=dist%info%ys,dist%info%ye
          i = 1
          if (walls(i,j).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,X_DIRECTION) > 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j)
                   fi(1,n,i,j) = fi(2,n,i,j)
                   fi(2,n,i,j) = tmp
                end if
             end do
          end if
       end do
    end if

    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_PSEUDOPERIODIC).and.(dist%info%xe.eq.dist%info%NX)) &
         then
       if (.not.dist%info%periodic(X_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_x', ierr)
          return
       end if
       do j=dist%info%ys,dist%info%ye
          i = dist%info%NX
          if (walls(i,j).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,X_DIRECTION) < 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j)
                   fi(1,n,i,j) = fi(2,n,i,j)
                   fi(2,n,i,j) = tmp
                end if
             end do
          end if
       end do
    end if

    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_PSEUDOPERIODIC).and.(dist%info%ys.eq.1)) then
       if (.not.dist%info%periodic(Y_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_y', ierr)
          return
       end if
       do i=dist%info%xs,dist%info%xe
          j = 1
          if (walls(i,j).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Y_DIRECTION) > 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j)
                   fi(1,n,i,j) = fi(2,n,i,j)
                   fi(2,n,i,j) = tmp
                end if
             end do
          end if
       end do
    end if

    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_PSEUDOPERIODIC).and.(dist%info%ye.eq.dist%info%NY)) &
         then
       if (.not.dist%info%periodic(Y_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_y', ierr)
          return
       end if
       do i=dist%info%xs,dist%info%xe
          j = dist%info%NY
          if (walls(i,j).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Y_DIRECTION) < 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j)
                   fi(1,n,i,j) = fi(2,n,i,j)
                   fi(2,n,i,j) = tmp
                end if
             end do
          end if
       end do
    end if
  end subroutine BCApplyPseudoperiodicD2

  subroutine BCApplyZeroGradientD3(bc, fi, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls

    PetscInt i,j,k,n
    PetscScalar tmp
    PetscErrorCode ierr

    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_PSEUDOPERIODIC).and.(dist%info%xs.eq.1)) then
       if (.not.dist%info%periodic(X_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_x', ierr)
          return
       end if
       do k=dist%info%zs,dist%info%ze
       do j=dist%info%ys,dist%info%ye
          i = 1
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,X_DIRECTION) > 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if

    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_PSEUDOPERIODIC).and.(dist%info%xe.eq.dist%info%NX)) &
         then
       if (.not.dist%info%periodic(X_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_x', ierr)
          return
       end if
       do k=dist%info%zs,dist%info%ze
       do j=dist%info%ys,dist%info%ye
          i = dist%info%NX
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,X_DIRECTION) < 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if

    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_PSEUDOPERIODIC).and.(dist%info%ys.eq.1)) then
       if (.not.dist%info%periodic(Y_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_y', ierr)
          return
       end if
       do k=dist%info%zs,dist%info%ze
       do i=dist%info%xs,dist%info%xe
          j = 1
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Y_DIRECTION) > 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if

    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_PSEUDOPERIODIC).and.(dist%info%ye.eq.dist%info%NY)) &
         then
       if (.not.dist%info%periodic(Y_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_y', ierr)
          return
       end if
       do k=dist%info%zs,dist%info%ze
       do i=dist%info%xs,dist%info%xe
          j = dist%info%NY
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Y_DIRECTION) < 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if

    ! ZM BOUNDARY
    if ((bc%flags(BOUNDARY_ZM).eq.BC_PSEUDOPERIODIC).and.(dist%info%zs.eq.1)) then
       if (.not.dist%info%periodic(Z_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_z', ierr)
          return
       end if
       do j=dist%info%ys,dist%info%ye
       do i=dist%info%xs,dist%info%xe
          k = 1
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Z_DIRECTION) > 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if

    ! ZM BOUNDARY
    if ((bc%flags(BOUNDARY_ZP).eq.BC_PSEUDOPERIODIC).and.(dist%info%ze.eq.dist%info%NZ)) &
         then
       if (.not.dist%info%periodic(Z_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_z', ierr)
          return
       end if
       do j=dist%info%ys,dist%info%ye
       do i=dist%info%xs,dist%info%xe
          k = dist%info%NZ
          if (walls(i,j,k).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Z_DIRECTION) < 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j,k)
                   fi(1,n,i,j,k) = fi(2,n,i,j,k)
                   fi(2,n,i,j,k) = tmp
                end if
             end do
          end if
       end do
       end do
    end if
  end subroutine BCApplyZeroGradientD3

  subroutine BCApplyZeroGradientD2(bc, fi, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls

    PetscInt i,j,n
    PetscScalar tmp
    PetscErrorCode ierr

    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_PSEUDOPERIODIC).and.(dist%info%xs.eq.1)) then
       if (.not.dist%info%periodic(X_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_x', ierr)
          return
       end if
       do j=dist%info%ys,dist%info%ye
          i = 1
          if (walls(i,j).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,X_DIRECTION) > 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j)
                   fi(1,n,i,j) = fi(2,n,i,j)
                   fi(2,n,i,j) = tmp
                end if
             end do
          end if
       end do
    end if

    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_PSEUDOPERIODIC).and.(dist%info%xe.eq.dist%info%NX)) &
         then
       if (.not.dist%info%periodic(X_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_x', ierr)
          return
       end if
       do j=dist%info%ys,dist%info%ye
          i = dist%info%NX
          if (walls(i,j).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,X_DIRECTION) < 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j)
                   fi(1,n,i,j) = fi(2,n,i,j)
                   fi(2,n,i,j) = tmp
                end if
             end do
          end if
       end do
    end if

    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_PSEUDOPERIODIC).and.(dist%info%ys.eq.1)) then
       if (.not.dist%info%periodic(Y_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_y', ierr)
          return
       end if
       do i=dist%info%xs,dist%info%xe
          j = 1
          if (walls(i,j).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Y_DIRECTION) > 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j)
                   fi(1,n,i,j) = fi(2,n,i,j)
                   fi(2,n,i,j) = tmp
                end if
             end do
          end if
       end do
    end if

    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_PSEUDOPERIODIC).and.(dist%info%ye.eq.dist%info%NY)) &
         then
       if (.not.dist%info%periodic(Y_DIRECTION)) then
          SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic also must get -bc_periodic_y', ierr)
          return
       end if
       do i=dist%info%xs,dist%info%xe
          j = dist%info%NY
          if (walls(i,j).eq.0) then
             do n=1,dist%b
                if (dist%disc%ci(n,Y_DIRECTION) < 0) then
                   ! has a component in the inward-normal direction
                   tmp = fi(1,n,i,j)
                   fi(1,n,i,j) = fi(2,n,i,j)
                   fi(2,n,i,j) = tmp
                end if
             end do
          end if
       end do
    end if
  end subroutine BCApplyZeroGradientD2
end module LBM_BC_module
