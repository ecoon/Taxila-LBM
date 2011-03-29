!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        bc.F90
!!!     version:         
!!!     created:         06 December 2010
!!!       on:            09:03:18 MST
!!!     last modified:   28 March 2011
!!!       at:            13:26:01 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ ldeo.columbia.edu
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module LBM_BC_module
  use LBM_Options_module
  use LBM_Info_module
  use LBM_Flow_module
  use LBM_Grid_module
  use petsc
  implicit none
  private
#include "lbm_definitions.h"

  type, public:: bc_type
     MPI_Comm comm
     type(grid_type),pointer:: grid
     PetscInt, dimension(6):: flags ! enum for boundary conditions
     PetscInt dim
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
       BCSetFromOptions, &
       BCSetUp, &
       BCSetValues, &
       BCGetArrays, &
       BCRestoreArrays, &
       BCFlux, &
       BCPressure, &
       BCPseudoperiodic
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
    bc%dim = 0
    
    bc%xm = 0
    bc%xp = 0
    bc%ym = 0
    bc%yp = 0
    bc%zm = 0
    bc%zp = 0
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

  ! set up the vectors for holding boundary data
  subroutine BCSetFromOptions(bc, options, ierr)
    type(bc_type) bc
    type(options_type) options
    PetscErrorCode ierr

    ! locals
    PetscBool flag
    PetscInt nmax
    PetscBool bcvalue

    ! lots of internals for flags
    PetscScalar xp3_ave, xm3_ave, xp3_max, xm3_max
    PetscScalar,dimension(1:options%nphases):: xp3_ave_p, xm3_ave_p
    PetscScalar yp3_ave, ym3_ave, yp3_max, ym3_max
    PetscScalar,dimension(1:options%nphases):: yp3_ave_p, ym3_ave_p
    PetscScalar zp3_ave, zm3_ave, zp3_max, zm3_max
    PetscScalar,dimension(1:options%nphases):: zp3_ave_p, zm3_ave_p

    ! dimension 
    bc%dim = MAX(options%ndims, options%nphases)

    ! flags and constant values from options
    nmax = 6
    call PetscOptionsGetIntArray(options%my_prefix, '-bc_flags', bc%flags, &
         nmax, flag, ierr)

    ! parse all potential boundary conditions
    ! xm boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_xm', bcvalue, flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_XM) = BC_PRESSURE

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_xm', bcvalue, flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_XM) = BC_FLUX

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_xm_poiseuille', bcvalue, &
         flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_XM) = BC_FLUX

    ! xp boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_xp', bcvalue, flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_XP) = BC_PRESSURE

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_xp', bcvalue, flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_XP) = BC_FLUX

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_xp_poiseuille', bcvalue, &
         flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_XP) = BC_FLUX

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_x', bcvalue, flag, ierr)
    if (bcvalue) then
       bc%flags(BOUNDARY_XM) = BC_PERIODIC
       bc%flags(BOUNDARY_XP) = BC_PERIODIC
    end if

    ! ym boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_ym', bcvalue, flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_YM) = BC_PRESSURE

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_ym', bcvalue, flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_YM) = BC_FLUX

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_ym_poiseuille', bcvalue, &
         flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_YM) = BC_FLUX

    ! yp boundary
    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_yp', bcvalue, flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_YP) = BC_PRESSURE
    

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_yp', bcvalue, flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_YP) = BC_FLUX

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_flux_yp_poiseuille', bcvalue, &
         flag, ierr)
    if (bcvalue) bc%flags(BOUNDARY_YP) = BC_FLUX

    bcvalue = PETSC_FALSE
    call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_y', bcvalue, flag, ierr)
    if (bcvalue) then
       bc%flags(BOUNDARY_YM) = BC_PERIODIC
       bc%flags(BOUNDARY_YP) = BC_PERIODIC
    end if

    if (options%ndims > 2) then
       ! zm boundary
       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_zm', bcvalue, flag, ierr)
       if (bcvalue) bc%flags(BOUNDARY_ZM) = BC_PRESSURE

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_flux_zm', bcvalue, flag, ierr)
       if (bcvalue) bc%flags(BOUNDARY_ZM) = BC_FLUX

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_flux_zm_poiseuille', bcvalue, &
            flag, ierr)
       if (bcvalue) bc%flags(BOUNDARY_ZM) = BC_FLUX

       ! zp boundary
       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_pressure_zp', bcvalue, flag, ierr)
       if (bcvalue) bc%flags(BOUNDARY_ZP) = BC_PRESSURE

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_flux_zp', bcvalue, flag, ierr)
       if (bcvalue) bc%flags(BOUNDARY_ZP) = BC_FLUX

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_flux_zp_poiseuille', bcvalue, &
            flag, ierr)
       if (bcvalue) bc%flags(BOUNDARY_ZP) = BC_FLUX

       bcvalue = PETSC_FALSE
       call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_z', bcvalue, flag, ierr)
       if (bcvalue) then
          bc%flags(BOUNDARY_ZM) = BC_PERIODIC
          bc%flags(BOUNDARY_ZP) = BC_PERIODIC
       end if
    end if

    call BCRestoreArrays(bc, ierr)
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
       locn = info%yl*info%zl*bc%dim
    else
       locn = 0
    endif
    call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%xm, ierr)
    call VecSetBlockSize(bc%xm, bc%dim, ierr)
    call PetscObjectSetName(bc%xm, 'xm_bc', ierr)
    
    if (info%xe.eq.info%NX) then
       locn = info%yl*info%zl*bc%dim
    else
       locn = 0
    endif
    call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%xp, ierr)
    call VecSetBlockSize(bc%xp, bc%dim, ierr)
    call PetscObjectSetName(bc%xp, 'xp_bc', ierr)
    
    ! y boundaries
    if (info%ys.eq.1) then
       locn = info%xl*info%zl*bc%dim
    else
       locn = 0
    endif
    call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%ym, ierr)
    call VecSetBlockSize(bc%ym, bc%dim, ierr)
    call PetscObjectSetName(bc%ym, 'ym_bc', ierr)
    
    if (info%ye.eq.info%NY) then
       locn = info%xl*info%zl*bc%dim
    else
       locn = 0
    endif
    call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%yp, ierr)
    call VecSetBlockSize(bc%yp, bc%dim, ierr)
    call PetscObjectSetName(bc%yp, 'yp_bc', ierr)
    
    if (info%ndims > 2) then
       ! z boundaries
       if (info%zs.eq.1) then
          locn = info%xl*info%yl*bc%dim
       else
          locn = 0
       endif
       call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%zm, ierr)
       call VecSetBlockSize(bc%zm, bc%dim, ierr)
       call PetscObjectSetName(bc%zm, 'zm_bc', ierr)
       
       if (info%ze.eq.info%NZ) then
          locn = info%xl*info%yl*bc%dim
       else
          locn = 0
       endif
       call VecCreateMPI(bc%comm, locn, PETSC_DETERMINE, bc%zp, ierr)
       call VecSetBlockSize(bc%zp, bc%dim, ierr)
       call PetscObjectSetName(bc%zp, 'zp_bc', ierr)
    end if

    call BCGetArrays(bc, ierr)
  end subroutine BCSetUp


  ! call initialize
  subroutine BCSetValues(bc, info, bc_subroutine)
    type(bc_type) bc
    type(info_type) info
    
    !  interface 
    !     subroutine bc_subroutine(xm, xp, ym, yp, zm, zp, dim, info)
    !       use LBM_Info_module
    !       PetscScalar xm(:)
    !       PetscScalar xp(:)
    !       PetscScalar ym(:)
    !       PetscScalar yp(:)
    !       PetscScalar zm(:)
    !       PetscScalar zp(:)
    !       PetscInt dim
    !       type(info_type) info
    !     end subroutine bc_subroutine
    !  end interface
    external bc_subroutine
    
    PetscErrorCode ierr

    call BCGetArrays(bc, ierr)

    call bc_subroutine(bc%flags, bc%xm_a, bc%xp_a, &
         bc%ym_a, bc%yp_a, bc%zm_a, bc%zp_a, bc%dim, info)
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
  
  subroutine BCPressure(bc, fi, walls, flow)
    type(bc_type) bc
    type(flow_type) flow
    
    PetscScalar,dimension(1:flow%nphases, 0:flow%disc%b, 1:flow%grid%info%gxyzl):: fi
    PetscScalar,dimension(1:flow%grid%info%gxyzl):: walls
    PetscErrorCode ierr

    select case(flow%disc%name)
    case(D3Q19_DISCRETIZATION)
       call BCPressureD3Q19(bc,fi,walls,bc%xm_a,bc%xp_a,bc%ym_a,bc%yp_a, &
            bc%zm_a,bc%zp_a,flow%grid%info)
    case(D2Q9_DISCRETIZATION)
       call BCPressureD2Q9(bc,fi,walls,bc%xm_a,bc%xp_a,bc%ym_a,bc%yp_a,flow%grid%info)
    case DEFAULT
       SETERRQ(1, 1, 'invalid discretization in LBM', ierr)
    end select
  end subroutine BCPressure

  subroutine BCFlux(bc, fi, walls, flow)
    type(bc_type) bc
    type(flow_type) flow
    
    PetscScalar,dimension(1:flow%nphases, 0:flow%disc%b, 1:flow%grid%info%gxyzl):: fi
    PetscScalar,dimension(1:flow%grid%info%gxyzl):: walls
    PetscErrorCode ierr

    select case(flow%disc%name)
    case(D3Q19_DISCRETIZATION)
       call BCFluxD3Q19(bc, fi, walls, bc%xm_a, bc%xp_a, bc%ym_a, bc%yp_a, &
            bc%zm_a, bc%zp_a, flow%grid%info)
    case(D2Q9_DISCRETIZATION)
       call BCFluxD2Q9(bc, fi, walls, bc%xm_a, bc%xp_a, bc%ym_a, bc%yp_a, flow%grid%info)
    case DEFAULT
       SETERRQ(1, 1, 'invalid discretization in LBM', ierr)
    end select
  end subroutine BCFlux

  subroutine BCPseudoperiodic(bc, fi, walls, flow)
    type(bc_type) bc
    type(flow_type) flow
    
    PetscScalar,dimension(1:flow%nphases, 0:flow%disc%b, 1:flow%grid%info%gxyzl):: fi
    PetscScalar,dimension(1:flow%grid%info%gxyzl):: walls
    PetscErrorCode ierr

    select case(flow%disc%name)
    case(D3Q19_DISCRETIZATION)
       call BCPseudoperiodicD3Q19(bc, fi, walls, bc%xm_a, bc%xp_a, bc%ym_a, bc%yp_a, &
            bc%zm_a, bc%zp_a, flow%grid%info)
    case(D2Q9_DISCRETIZATION)
       SETERRQ(PETSC_COMM_SELF, 1, 'pseudoperiodic bcs in 2d are not implemented', ierr)
    case DEFAULT
       SETERRQ(PETSC_COMM_SELF, 1, 'invalid discretization in LBM', ierr)
    end select
  end subroutine BCPseudoperiodic

  subroutine BCPressureD3Q19(bc, fi, walls, xm_vals, xp_vals, &
       ym_vals, yp_vals, zm_vals, zp_vals, info)

    type(bc_type) bc
    type(info_type) info

    PetscScalar,dimension(1:info%nphases, 0:info%flow_b, info%gxs:info%gxe, &
         info%gys:info%gye, info%gzs:info%gze):: fi
    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: walls
    PetscScalar,dimension(bc%dim,info%ys:info%ye,info%zs:info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(bc%dim,info%xs:info%xe,info%zs:info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(bc%dim,info%xs:info%xe,info%ys:info%ye):: zm_vals, zp_vals

    PetscInt i,j,k
    PetscInt directions(0:info%flow_b)
    PetscInt cardinals(1:info%ndims)

    directions(:) = 0
    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_PRESSURE).and.(info%xs.eq.1)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_XM, directions, cardinals)
       do k=info%zs,info%ze
          do j=info%ys,info%ye
             i = 1
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyD3Q19(bc, fi(:,:,i,j,k), xm_vals(:,j,k), &
                     directions, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_PRESSURE).and.(info%xe.eq.info%NX)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_XP, directions, cardinals)
       do k=info%zs,info%ze
          do j=info%ys,info%ye
             i = info%NX
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyD3Q19(bc, fi(:,:,i,j,k), xp_vals(:,j,k), &
                     directions, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_PRESSURE).and.(info%ys.eq.1)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_YM, directions, cardinals)
       do k=info%zs,info%ze
          do i=info%xs,info%xe
             j = 1
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyD3Q19(bc, fi(:,:,i,j,k), ym_vals(:,i,k), &
                     directions, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_PRESSURE).and.(info%ye.eq.info%NY)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_YP, directions, cardinals)
       do k=info%zs,info%ze
          do i=info%xs,info%xe
             j = info%NY
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyD3Q19(bc, fi(:,:,i,j,k), yp_vals(:,i,k), &
                     directions, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! ZM BOUNDARY
    if ((bc%flags(BOUNDARY_ZM).eq.BC_PRESSURE).and.(info%zs.eq.1)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_ZM, directions, cardinals)
       do j=info%ys,info%ye
          do i=info%xs,info%xe
             k = 1
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyD3Q19(bc, fi(:,:,i,j,k), zm_vals(:,i,j), &
                     directions, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! ZP BOUNDARY
    if ((bc%flags(BOUNDARY_ZP).eq.BC_PRESSURE).and.(info%ze.eq.info%NZ)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_ZP, directions, cardinals)
       do j=info%ys,info%ye
          do i=info%xs,info%xe
             k = info%NZ
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyD3Q19(bc, fi(:,:,i,j,k), zp_vals(:,i,j), &
                     directions, info)
             end if
          end do
       end do
    endif

    return
  end subroutine BCPressureD3Q19

  subroutine BCPressureApplyD3Q19(bc, fi, pvals, directions, info)
    use LBM_Discretization_Directions_D3Q19_module

    type(bc_type) bc
    type(info_type) info
    PetscInt,intent(in),dimension(0:info%flow_b):: directions
    PetscScalar,intent(inout),dimension(1:info%nphases, 0:info%flow_b):: fi
    PetscScalar,intent(in),dimension(bc%dim):: pvals

!    PetscScalar utmp, vtmp
    PetscScalar wtmp
    PetscScalar,dimension(0:info%flow_b)::ftmp
    integer m
    
!    utmp = 0
!    vtmp = 0
    wtmp = 0

    do m=1,info%nphases
       ftmp = 0.0
       wtmp = fi(m,directions(ORIGIN)) + fi(m,directions(EAST)) &
            + fi(m,directions(NORTH)) + fi(m,directions(WEST)) &
            + fi(m,directions(SOUTH)) + fi(m,directions(NORTHEAST)) &
            + fi(m,directions(NORTHWEST)) + fi(m,directions(SOUTHWEST)) &
            + fi(m,directions(SOUTHEAST)) + 2.*(fi(m,directions(DOWN)) &
            + fi(m,directions(WESTDOWN)) + fi(m,directions(EASTDOWN)) &
            + fi(m,directions(SOUTHDOWN)) + fi(m,directions(NORTHDOWN)))
       wtmp = 1.0-wtmp/pvals(m)
       
       
       ! Choice should not affect the momentum significantly
       ftmp(directions(UP)) = fi(m,directions(DOWN))
       ftmp(directions(EASTUP)) = fi(m,directions(WESTDOWN))
       ftmp(directions(WESTUP)) = fi(m,directions(EASTDOWN))
       ftmp(directions(NORTHUP)) = fi(m,directions(SOUTHDOWN))
       ftmp(directions(SOUTHUP)) = fi(m,directions(NORTHDOWN))
       
       fi(m,directions(UP)) = 2./3.*ftmp(directions(UP)) &
            + 1./3.*pvals(m)*wtmp &
            - 1./3.*(ftmp(directions(EASTUP))+ ftmp(directions(WESTUP)) &
            + ftmp(directions(NORTHUP)) + ftmp(directions(SOUTHUP))) &
            + 1./3.*(fi(m,directions(DOWN)) &
            + fi(m,directions(WESTDOWN)) + fi(m,directions(EASTDOWN)) &
            + fi(m,directions(SOUTHDOWN)) + fi(m,directions(NORTHDOWN)))
       
       fi(m,directions(EASTUP)) = 1./3.*ftmp(directions(EASTUP)) &
!            + 1./2.*pvals(m)*utmp &
            + 1./6.*pvals(m)*wtmp &
            - 1./2.*(fi(m,directions(EAST)) &
            - fi(m,directions(WEST)) + fi(m,directions(NORTHEAST)) &
            - fi(m,directions(NORTHWEST)) - fi(m,directions(SOUTHWEST)) &
            + fi(m,directions(SOUTHEAST))) &
            - 1./6.*(ftmp(directions(UP)) - fi(m,directions(DOWN)) &
            + ftmp(directions(NORTHUP)) + ftmp(directions(SOUTHUP)) &
            - fi(m,directions(SOUTHDOWN)) - fi(m,directions(NORTHDOWN))) &
            + 1./3.*(ftmp(directions(WESTUP)) - fi(m,directions(EASTDOWN))) &
            + 2./3.*fi(m,directions(WESTDOWN))
       
       fi(m,directions(WESTUP)) = 1./3.*ftmp(directions(WESTUP)) &
!            - 1./2.*pvals(m)*utmp &
            + 1./6.*pvals(m)*wtmp &
            + 1./2.*(fi(m,directions(EAST)) &
            - fi(m,directions(WEST)) + fi(m,directions(NORTHEAST)) &
            - fi(m,directions(NORTHWEST)) - fi(m,directions(SOUTHWEST)) &
            + fi(m,directions(SOUTHEAST))) &
            - 1./6.*(ftmp(directions(UP)) - fi(m,directions(DOWN)) &
            + ftmp(directions(NORTHUP)) + ftmp(directions(SOUTHUP)) &
            - fi(m,directions(SOUTHDOWN)) - fi(m,directions(NORTHDOWN))) &
            + 1./3.*(ftmp(directions(EASTUP)) - fi(m,directions(WESTDOWN))) &
            + 2./3.*fi(m,directions(EASTDOWN))
       
       fi(m,directions(NORTHUP)) = 1./3.*ftmp(directions(NORTHUP)) &
!            + 1./2.*pvals(m)*vtmp &
            + 1./6.*pvals(m)*wtmp &
            - 1./2.*(fi(m,directions(NORTH)) &
            - fi(m,directions(SOUTH)) + fi(m,directions(NORTHEAST)) &
            + fi(m,directions(NORTHWEST)) - fi(m,directions(SOUTHWEST)) &
            - fi(m,directions(SOUTHEAST))) &
            - 1./6.*(ftmp(directions(UP)) - fi(m,directions(DOWN)) &
            + ftmp(directions(EASTUP)) + ftmp(directions(WESTUP)) &
            - fi(m,directions(WESTDOWN)) - fi(m,directions(EASTDOWN))) &
            + 1./3.*(ftmp(directions(SOUTHUP)) - fi(m,directions(NORTHDOWN))) &
            + 2./3.*fi(m,directions(SOUTHDOWN))
       
       fi(m,directions(SOUTHUP)) = 1./3.*ftmp(directions(SOUTHUP)) &
!            - 1./2.*pvals(m)*vtmp &
            + 1./6.*pvals(m)*wtmp &
            + 1./2.*(fi(m,directions(NORTH)) &
            - fi(m,directions(SOUTH)) + fi(m,directions(NORTHEAST)) &
            + fi(m,directions(NORTHWEST)) - fi(m,directions(SOUTHWEST)) &
            - fi(m,directions(SOUTHEAST))) &
            - 1./6.*(ftmp(directions(UP)) - fi(m,directions(DOWN)) &
            + ftmp(directions(EASTUP)) + ftmp(directions(WESTUP)) &
            - fi(m,directions(WESTDOWN)) - fi(m,directions(EASTDOWN))) &
            + 1./3.*(ftmp(directions(NORTHUP)) - fi(m,directions(SOUTHDOWN))) &
            + 2./3.*fi(m,directions(NORTHDOWN))
    enddo
    return
  end subroutine BCPressureApplyD3Q19

  subroutine BCFluxD3Q19(bc, fi, walls, xm_vals, xp_vals, &
       ym_vals, yp_vals, zm_vals, zp_vals, info)

    type(bc_type) bc
    type(info_type) info

    PetscScalar,dimension(1:info%nphases, 0:info%flow_b, info%gxs:info%gxe, &
         info%gys:info%gye, info%gzs:info%gze):: fi
    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: walls
    PetscScalar,dimension(bc%dim,info%ys:info%ye,info%zs:info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(bc%dim,info%xs:info%xe,info%zs:info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(bc%dim,info%xs:info%xe,info%ys:info%ye):: zm_vals, zp_vals

    PetscInt i,j,k
    PetscInt directions(0:info%flow_b)
    PetscInt cardinals(1:info%ndims)

    directions(:) = 0
    cardinals(:) = 0
    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_FLUX).and.(info%xs.eq.1)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_XM, directions, cardinals)
       do k=info%zs,info%ze
          do j=info%ys,info%ye
             i = 1
             if (walls(i,j,k).eq.0) then
                call BCFluxApplyD3Q19(bc, fi(:,:,i,j,k), xm_vals(:,j,k), &
                     directions, cardinals, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    cardinals(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_FLUX).and.(info%xe.eq.info%NX)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_XP, directions, cardinals)
       do k=info%zs,info%ze
          do j=info%ys,info%ye
             i = info%NX
             if (walls(i,j,k).eq.0) then
                call BCFluxApplyD3Q19(bc, fi(:,:,i,j,k), xp_vals(:,j,k), &
                     directions, cardinals, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    cardinals(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_FLUX).and.(info%ys.eq.1)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_YM, directions, cardinals)
       do k=info%zs,info%ze
          do i=info%xs,info%xe
             j = 1
             if (walls(i,j,k).eq.0) then
                call BCFluxApplyD3Q19(bc, fi(:,:,i,j,k), ym_vals(:,i,k), &
                     directions, cardinals, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    cardinals(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_FLUX).and.(info%ye.eq.info%NY)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_YP, directions, cardinals)
       do k=info%zs,info%ze
          do i=info%xs,info%xe
             j = info%NY
             if (walls(i,j,k).eq.0) then
                call BCFluxApplyD3Q19(bc, fi(:,:,i,j,k), yp_vals(:,i,k), &
                     directions, cardinals, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    cardinals(:) = 0
    ! ZM BOUNDARY
    if ((bc%flags(BOUNDARY_ZM).eq.BC_FLUX).and.(info%zs.eq.1)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_ZM, directions, cardinals)
       do j=info%ys,info%ye
          do i=info%xs,info%xe
             k = 1
             if (walls(i,j,k).eq.0) then
                call BCFluxApplyD3Q19(bc, fi(:,:,i,j,k), zm_vals(:,i,j), &
                     directions, cardinals, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    cardinals(:) = 0
    ! ZP BOUNDARY
    if ((bc%flags(BOUNDARY_ZP).eq.BC_FLUX).and.(info%ze.eq.info%NZ)) then
       call BCSetLocalDirectionsD3Q19(BOUNDARY_ZP, directions, cardinals)
       do j=info%ys,info%ye
          do i=info%xs,info%xe
             k = info%NZ
             if (walls(i,j,k).eq.0) then
                call BCFluxApplyD3Q19(bc, fi(:,:,i,j,k), zp_vals(:,i,j), &
                     directions, cardinals, info)
             end if
          end do
       end do
    endif

    return
  end subroutine BCFluxD3Q19
  
  subroutine BCFluxApplyD3Q19(bc, fi, fvals, directions, cardinals, info)
    use LBM_Discretization_Directions_D3Q19_module

    type(bc_type) bc
    type(info_type) info

    PetscScalar,intent(inout),dimension(1:info%nphases, 0:info%flow_b):: fi
    PetscScalar,intent(in),dimension(bc%dim):: fvals
    PetscInt,intent(in),dimension(0:info%flow_b):: directions
    PetscInt,intent(in),dimension(1:info%ndims):: cardinals

    PetscScalar rhotmp
    PetscScalar,dimension(0:info%flow_b)::ftmp
    PetscInt m

    ftmp = 0.0
    rhotmp = 0.0

    do m=1,info%nphases
       rhotmp = fi(m,directions(ORIGIN)) + fi(m,directions(EAST)) &
            + fi(m,directions(NORTH)) + fi(m,directions(WEST)) &
            + fi(m,directions(SOUTH)) + fi(m,directions(NORTHEAST)) &
            + fi(m,directions(NORTHWEST)) + fi(m,directions(SOUTHWEST)) &
            + fi(m,directions(SOUTHEAST)) + 2.*(fi(m,directions(DOWN)) &
            + fi(m,directions(WESTDOWN)) + fi(m,directions(EASTDOWN)) &
            + fi(m,directions(SOUTHDOWN)) + fi(m,directions(NORTHDOWN))) 
       rhotmp = rhotmp/(1. - fvals(cardinals(CARDINAL_NORMAL)))
       
       ! Choice should not affect the momentum significantly
       ftmp(directions(UP)) = fi(m,directions(DOWN))
       ftmp(directions(EASTUP)) = fi(m,directions(WESTDOWN))
       ftmp(directions(WESTUP)) = fi(m,directions(EASTDOWN))
       ftmp(directions(NORTHUP)) = fi(m,directions(SOUTHDOWN))
       ftmp(directions(SOUTHUP)) = fi(m,directions(NORTHDOWN))

       fi(m,directions(UP)) = 2./3.*ftmp(directions(UP)) &
            + 1./3.*rhotmp*fvals(cardinals(CARDINAL_NORMAL)) &
            - 1./3.*(ftmp(directions(EASTUP))+ ftmp(directions(WESTUP)) &
            + ftmp(directions(NORTHUP)) + ftmp(directions(SOUTHUP))) &
            + 1./3.*(fi(m,directions(DOWN)) &
            + fi(m,directions(WESTDOWN)) + fi(m,directions(EASTDOWN)) &
            + fi(m,directions(SOUTHDOWN)) + fi(m,directions(NORTHDOWN)))
       
       fi(m,directions(EASTUP)) = 1./3.*ftmp(directions(EASTUP)) &
            + 1./2.*rhotmp*fvals(cardinals(CARDINAL_CROSS)) &
            + 1./6.*rhotmp*fvals(cardinals(CARDINAL_NORMAL)) &
            - 1./2.*(fi(m,directions(EAST)) &
            - fi(m,directions(WEST)) + fi(m,directions(NORTHEAST)) &
            - fi(m,directions(NORTHWEST)) - fi(m,directions(SOUTHWEST)) &
            + fi(m,directions(SOUTHEAST))) &
            - 1./6.*(ftmp(directions(UP)) - fi(m,directions(DOWN)) &
            + ftmp(directions(NORTHUP)) + ftmp(directions(SOUTHUP)) &
            - fi(m,directions(SOUTHDOWN)) - fi(m,directions(NORTHDOWN))) &
            + 1./3.*(ftmp(directions(WESTUP)) - fi(m,directions(EASTDOWN))) &
            + 2./3.*fi(m,directions(WESTDOWN)) 

       fi(m,directions(WESTUP)) = 1./3.*ftmp(directions(WESTUP)) &
            - 1./2.*rhotmp*fvals(cardinals(CARDINAL_CROSS)) &
            + 1./6.*rhotmp*fvals(cardinals(CARDINAL_NORMAL)) &
            + 1./2.*(fi(m,directions(EAST)) &
            - fi(m,directions(WEST)) + fi(m,directions(NORTHEAST)) &
            - fi(m,directions(NORTHWEST)) - fi(m,directions(SOUTHWEST)) &
            + fi(m,directions(SOUTHEAST))) &
            - 1./6.*(ftmp(directions(UP)) - fi(m,directions(DOWN)) &
            + ftmp(directions(NORTHUP)) + ftmp(directions(SOUTHUP)) &
            - fi(m,directions(SOUTHDOWN)) - fi(m,directions(NORTHDOWN))) &
            + 1./3.*(ftmp(directions(EASTUP)) - fi(m,directions(WESTDOWN))) &
            + 2./3.*fi(m,directions(EASTDOWN))
       
       fi(m,directions(NORTHUP)) = 1./3.*ftmp(directions(NORTHUP)) &
            + 1./2.*rhotmp*fvals(cardinals(CARDINAL_RESULTANT)) &
            + 1./6.*rhotmp*fvals(cardinals(CARDINAL_NORMAL)) &
            - 1./2.*(fi(m,directions(NORTH)) &
            - fi(m,directions(SOUTH)) + fi(m,directions(NORTHEAST)) &
            + fi(m,directions(NORTHWEST)) - fi(m,directions(SOUTHWEST)) &
            - fi(m,directions(SOUTHEAST))) &
            - 1./6.*(ftmp(directions(UP)) - fi(m,directions(DOWN)) &
            + ftmp(directions(EASTUP)) + ftmp(directions(WESTUP)) &
            - fi(m,directions(WESTDOWN)) - fi(m,directions(EASTDOWN))) &
            + 1./3.*(ftmp(directions(SOUTHUP)) - fi(m,directions(NORTHDOWN))) &
            + 2./3.*fi(m,directions(SOUTHDOWN))
       
       fi(m,directions(SOUTHUP)) = 1./3.*ftmp(directions(SOUTHUP)) &
            - 1./2.*rhotmp*fvals(cardinals(CARDINAL_RESULTANT)) &
            + 1./6.*rhotmp*fvals(cardinals(CARDINAL_NORMAL)) &
            + 1./2.*(fi(m,directions(NORTH)) &
            - fi(m,directions(SOUTH)) + fi(m,directions(NORTHEAST)) &
            + fi(m,directions(NORTHWEST)) - fi(m,directions(SOUTHWEST)) &
            - fi(m,directions(SOUTHEAST))) &
            - 1./6.*(ftmp(directions(UP)) - fi(m,directions(DOWN)) &
            + ftmp(directions(EASTUP)) + ftmp(directions(WESTUP)) &
            - fi(m,directions(WESTDOWN)) - fi(m,directions(EASTDOWN))) &
            + 1./3.*(ftmp(directions(NORTHUP)) - fi(m,directions(SOUTHDOWN))) &
            + 2./3.*fi(m,directions(NORTHDOWN))
    enddo


    return
  end subroutine BCFluxApplyD3Q19

  subroutine BCSetLocalDirectionsD3Q19(boundary, directions, cardinals)
    use LBM_Discretization_Directions_D3Q19_module
    PetscInt,intent(in):: boundary
    PetscInt,intent(out),dimension(0:discretization_directions) :: directions
    PetscInt,intent(out),dimension(1:discretization_dims) :: cardinals

    PetscInt i
    
    select case(boundary)
    case (BOUNDARY_ZM)
       ! identity mapping
       do i=0,discretization_directions
          directions(i) = i
       end do
       cardinals(CARDINAL_NORMAL) = Z_DIRECTION
       cardinals(CARDINAL_CROSS) = X_DIRECTION
       cardinals(CARDINAL_RESULTANT) = Y_DIRECTION

    case (BOUNDARY_ZP)
       ! inverted mapping around the origin
       directions(ORIGIN) = ORIGIN
       directions(EAST) = WEST
       directions(WEST) = EAST
       directions(NORTH) = SOUTH
       directions(SOUTH) = NORTH
       directions(UP) = DOWN
       directions(DOWN) = UP
       directions(NORTHEAST) = SOUTHWEST
       directions(SOUTHWEST) = NORTHEAST
       directions(NORTHWEST) = SOUTHEAST
       directions(SOUTHEAST) = NORTHWEST
       directions(NORTHUP) = SOUTHDOWN
       directions(SOUTHDOWN) = NORTHUP
       directions(NORTHDOWN) = SOUTHUP
       directions(SOUTHUP) = NORTHDOWN
       directions(EASTUP) = WESTDOWN
       directions(WESTDOWN) = EASTUP
       directions(EASTDOWN) = WESTUP
       directions(WESTUP) = EASTDOWN

       cardinals(CARDINAL_NORMAL) = Z_DIRECTION
       cardinals(CARDINAL_CROSS) = X_DIRECTION
       cardinals(CARDINAL_RESULTANT) = Y_DIRECTION

    case (BOUNDARY_XM)
       ! map z -> x -> y
       directions(ORIGIN) = ORIGIN
       directions(EAST) = NORTH
       directions(WEST) = SOUTH
       directions(NORTH) = UP
       directions(SOUTH) = DOWN
       directions(UP) = EAST
       directions(DOWN) = WEST
       directions(NORTHEAST) = NORTHUP
       directions(SOUTHWEST) = SOUTHDOWN
       directions(NORTHWEST) = SOUTHUP
       directions(SOUTHEAST) = NORTHDOWN
       directions(NORTHUP) = EASTUP
       directions(SOUTHDOWN) = WESTDOWN
       directions(NORTHDOWN) = WESTUP
       directions(SOUTHUP) = EASTDOWN
       directions(EASTUP) = NORTHEAST
       directions(WESTDOWN) = SOUTHWEST
       directions(EASTDOWN) = NORTHWEST
       directions(WESTUP) = SOUTHEAST

       cardinals(CARDINAL_NORMAL) = X_DIRECTION
       cardinals(CARDINAL_CROSS) = Y_DIRECTION
       cardinals(CARDINAL_RESULTANT) = Z_DIRECTION

    case (BOUNDARY_XP)
       ! map z -> x -> y and invert
       directions(ORIGIN) = ORIGIN
       directions(EAST) = SOUTH
       directions(WEST) = NORTH
       directions(NORTH) = DOWN
       directions(SOUTH) = UP
       directions(UP) = WEST
       directions(DOWN) = EAST
       directions(NORTHEAST) = SOUTHDOWN
       directions(SOUTHWEST) = NORTHUP
       directions(NORTHWEST) = NORTHDOWN
       directions(SOUTHEAST) = SOUTHUP
       directions(NORTHUP) = WESTDOWN
       directions(SOUTHDOWN) = EASTUP
       directions(NORTHDOWN) = EASTDOWN
       directions(SOUTHUP) = WESTUP
       directions(EASTUP) = SOUTHWEST
       directions(WESTDOWN) = NORTHEAST
       directions(EASTDOWN) = SOUTHEAST
       directions(WESTUP) = NORTHWEST

       cardinals(CARDINAL_NORMAL) = X_DIRECTION
       cardinals(CARDINAL_CROSS) = Y_DIRECTION
       cardinals(CARDINAL_RESULTANT) = Z_DIRECTION

    case (BOUNDARY_YM)       
       ! cycle x <-- y <-- z
       directions(ORIGIN) = ORIGIN
       directions(EAST) = UP
       directions(WEST) = DOWN
       directions(NORTH) = EAST
       directions(SOUTH) = WEST
       directions(UP) = NORTH
       directions(DOWN) = SOUTH
       directions(NORTHEAST) = EASTUP
       directions(SOUTHWEST) = WESTDOWN
       directions(NORTHWEST) = EASTDOWN
       directions(SOUTHEAST) = WESTUP
       directions(NORTHUP) = NORTHEAST
       directions(SOUTHDOWN) = SOUTHWEST
       directions(NORTHDOWN) = SOUTHEAST
       directions(SOUTHUP) = NORTHWEST
       directions(EASTUP) = NORTHUP
       directions(WESTDOWN) = SOUTHDOWN
       directions(EASTDOWN) = SOUTHUP
       directions(WESTUP) = NORTHDOWN

       cardinals(CARDINAL_NORMAL) = Y_DIRECTION
       cardinals(CARDINAL_CROSS) = Z_DIRECTION
       cardinals(CARDINAL_RESULTANT) = X_DIRECTION

    case (BOUNDARY_YP)       
       ! cycle x <-- y <-- z, then invert
       directions(ORIGIN) = ORIGIN
       directions(EAST) = DOWN
       directions(WEST) = UP
       directions(NORTH) = WEST
       directions(SOUTH) = EAST
       directions(UP) = SOUTH
       directions(DOWN) = NORTH
       directions(NORTHEAST) = WESTDOWN
       directions(SOUTHWEST) = EASTUP
       directions(NORTHWEST) = WESTUP
       directions(SOUTHEAST) = EASTDOWN
       directions(NORTHUP) = SOUTHWEST
       directions(SOUTHDOWN) = NORTHEAST
       directions(NORTHDOWN) = NORTHWEST
       directions(SOUTHUP) = SOUTHEAST
       directions(EASTUP) = SOUTHDOWN
       directions(WESTDOWN) = NORTHUP
       directions(EASTDOWN) = NORTHDOWN
       directions(WESTUP) = SOUTHUP

       cardinals(CARDINAL_NORMAL) = Y_DIRECTION
       cardinals(CARDINAL_CROSS) = Z_DIRECTION
       cardinals(CARDINAL_RESULTANT) = X_DIRECTION

    end select
  end subroutine BCSetLocalDirectionsD3Q19

  subroutine BCPseudoperiodicD3Q19(bc, fi, walls, xm_vals, xp_vals, &
       ym_vals, yp_vals, zm_vals, zp_vals, info)

    type(bc_type) bc
    type(info_type) info

    PetscScalar,dimension(1:info%nphases, 0:info%flow_b, info%gxs:info%gxe, &
         info%gys:info%gye, info%gzs:info%gze):: fi
    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: walls

    PetscScalar,dimension(bc%dim,info%ys:info%ye,info%zs:info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(bc%dim,info%xs:info%xe,info%zs:info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(bc%dim,info%xs:info%xe,info%ys:info%ye):: zm_vals, zp_vals

    integer i,j,m
    PetscScalar ftmp

    ! note that pseudo-periodic is a requirement on both plus and minus
    ! boundaries.  Only one is checked, and it is assumed that the other
    ! one is correctly flagged

    ! In the y-direction
    if ((bc%flags(BOUNDARY_YP).eq.BC_PSEUDOPERIODIC).and.(info%ye.eq.info%NY)) then
       do i=info%xs,info%xe
          do j=info%zs,info%ze
             !   at the exit
             ftmp = fi(1,4,i,info%NY,j)
             fi(1,4,i,info%NY,j)=fi(2,4,i,info%NY,j)
             fi(2,4,i,info%NY,j)= ftmp

             ftmp = fi(1,9,i,info%NY,j)
             fi(1,9,i,info%NY,j)=fi(2,9,i,info%NY,j)
             fi(2,9,i,info%NY,j)= ftmp

             ftmp = fi(1,10,i,info%NY,j)
             fi(1,10,i,info%NY,j)=fi(2,10,i,info%NY,j)
             fi(2,10,i,info%NY,j)= ftmp

             ftmp = fi(1,16,i,info%NY,j)
             fi(1,16,i,info%NY,j)=fi(2,16,i,info%NY,j)
             fi(2,16,i,info%NY,j)= ftmp

             ftmp = fi(1,17,i,info%NY,j)
             fi(1,17,i,info%NY,j)=fi(2,17,i,info%NY,j)
             fi(2,17,i,info%NY,j)= ftmp
          enddo
       enddo
    endif

    if ((bc%flags(BOUNDARY_YM).eq.BC_PSEUDOPERIODIC).and.(info%ys.eq.1)) then
       do i=info%xs,info%xe
          do j=info%zs,info%ze
             !    at the entrance
             ftmp = fi(1,2,i,1,j)
             fi(1,2,i,1,j)=fi(2,2,i,1,j)
             fi(2,2,i,1,j)= ftmp

             ftmp = fi(1,7,i,1,j)
             fi(1,7,i,1,j)=fi(2,7,i,1,j)
             fi(2,7,i,1,j)= ftmp

             ftmp = fi(1,8,i,1,j)
             fi(1,8,i,1,j)=fi(2,8,i,1,j)
             fi(2,8,i,1,j)= ftmp

             ftmp = fi(1,15,i,1,j)
             fi(1,15,i,1,j)=fi(2,15,i,1,j)
             fi(2,15,i,1,j)= ftmp

             ftmp = fi(1,18,i,1,j)
             fi(1,18,i,1,j)=fi(2,18,i,1,j)
             fi(2,18,i,1,j)= ftmp

          end do
       end do
    endif

    ! In the z-direction (parallelized)
    if ((bc%flags(BOUNDARY_ZP).eq.BC_PSEUDOPERIODIC).and.(info%ze.eq.info%NZ)) then
       do i=info%xs,info%xe
          do j=info%ys,info%ye
             !  at the exit
             ftmp = fi(1,6,i,j,info%NZ)
             fi(1,6,i,j,info%NZ)=fi(2,6,i,j,info%NZ)
             fi(2,6,i,j,info%NZ)= ftmp

             ftmp = fi(1,13,i,j,info%NZ)
             fi(1,13,i,j,info%NZ)=fi(2,13,i,j,info%NZ)
             fi(2,13,i,j,info%NZ)= ftmp

             ftmp = fi(1,14,i,j,info%NZ)
             fi(1,14,i,j,info%NZ)=fi(2,14,i,j,info%NZ)
             fi(2,14,i,j,info%NZ)= ftmp

             ftmp = fi(1,17,i,j,info%NZ)
             fi(1,17,i,j,info%NZ)=fi(2,17,i,j,info%NZ)
             fi(2,17,i,j,info%NZ)= ftmp

             ftmp = fi(1,18,i,j,info%NZ)
             fi(1,18,i,j,info%NZ)=fi(2,18,i,j,info%NZ)
             fi(2,18,i,j,info%NZ)= ftmp
          end do
       end do
    endif

    if ((bc%flags(BOUNDARY_ZM).eq.BC_PSEUDOPERIODIC).and.(info%zs.eq.1)) then
       do i=info%xs,info%xe
          do j=info%ys,info%ye

             !  at the entrance
             ftmp = fi(1,5,i,j,1)
             fi(1,5,i,j,1)=fi(2,5,i,j,1)
             fi(2,5,i,j,1)= ftmp

             ftmp = fi(1,11,i,j,1)
             fi(1,11,i,j,1)=fi(2,11,i,j,1)
             fi(2,11,i,j,1)= ftmp

             ftmp = fi(1,12,i,j,1)
             fi(1,12,i,j,1)=fi(2,12,i,j,1)
             fi(2,12,i,j,1)= ftmp

             ftmp = fi(1,15,i,j,1)
             fi(1,15,i,j,1)=fi(2,15,i,j,1)
             fi(2,15,i,j,1)= ftmp

             ftmp = fi(1,16,i,j,1)
             fi(1,16,i,j,1)=fi(2,16,i,j,1)
             fi(2,16,i,j,1)= ftmp

          end do
       end do
    end if
    return
  end subroutine BCPseudoperiodicD3Q19

  subroutine BCPressureD2Q9(bc, fi, walls, xm_vals, xp_vals, &
       ym_vals, yp_vals, info)

    type(bc_type) bc
    type(info_type) info

    PetscScalar,dimension(1:info%nphases, 0:info%flow_b, info%gxs:info%gxe, info%gys:info%gye):: fi
    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye):: walls
    PetscScalar,dimension(bc%dim,info%ys:info%ye):: xm_vals, xp_vals
    PetscScalar,dimension(bc%dim,info%xs:info%xe):: ym_vals, yp_vals

    PetscInt i,j
    PetscInt directions(0:info%flow_b)
    PetscInt cardinals(1:info%ndims)

    directions(:) = 0
    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_PRESSURE).and.(info%xs.eq.1)) then
       call BCSetLocalDirectionsD2Q9(BOUNDARY_XM, directions, cardinals)
       do j=info%ys,info%ye
          i = 1
          if (walls(i,j).eq.0) then
             call BCPressureApplyD2Q9(bc, fi(:,:,i,j), xm_vals(:,j), directions, info)
          end if
       end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_PRESSURE).and.(info%xe.eq.info%NX)) then
       call BCSetLocalDirectionsD2Q9(BOUNDARY_XP, directions, cardinals)
       do j=info%ys,info%ye
          i = info%NX
          if (walls(i,j).eq.0) then
             call BCPressureApplyD2Q9(bc, fi(:,:,i,j), xp_vals(:,j), directions, info)
          end if
       end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_PRESSURE).and.(info%ys.eq.1)) then
       call BCSetLocalDirectionsD2Q9(BOUNDARY_YM, directions, cardinals)
       do i=info%xs,info%xe
          j = 1
          if (walls(i,j).eq.0) then
             call BCPressureApplyD2Q9(bc, fi(:,:,i,j), ym_vals(:,i), &
                  directions, info)
          end if
       end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_PRESSURE).and.(info%ye.eq.info%NY)) then
       call BCSetLocalDirectionsD2Q9(BOUNDARY_YP, directions, cardinals)
       do i=info%xs,info%xe
          j = info%NY
          if (walls(i,j).eq.0) then
             call BCPressureApplyD2Q9(bc, fi(:,:,i,j), yp_vals(:,i), &
                  directions, info)
          end if
       end do
    endif
    return
  end subroutine BCPressureD2Q9

  subroutine BCFluxD2Q9(bc, fi, walls, xm_vals, xp_vals, &
       ym_vals, yp_vals, info)

    type(bc_type) bc
    type(info_type) info

    PetscScalar,dimension(1:info%nphases, 0:info%flow_b, info%gxs:info%gxe, info%gys:info%gye):: fi
    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye):: walls
    PetscScalar,dimension(bc%dim,info%ys:info%ye):: xm_vals, xp_vals
    PetscScalar,dimension(bc%dim,info%xs:info%xe):: ym_vals, yp_vals

    PetscInt i,j
    PetscInt directions(0:info%flow_b)
    PetscInt cardinals(1:info%ndims)

    directions(:) = 0
    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_PRESSURE).and.(info%xs.eq.1)) then
       call BCSetLocalDirectionsD2Q9(BOUNDARY_XM, directions, cardinals)
       do j=info%ys,info%ye
          i = 1
          if (walls(i,j).eq.0) then
             call BCFluxApplyD2Q9(bc, fi(:,:,i,j), xm_vals(:,j), &
                  directions, cardinals, info)
          end if
       end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_PRESSURE).and.(info%xe.eq.info%NX)) then
       call BCSetLocalDirectionsD2Q9(BOUNDARY_XP, directions, cardinals)
       do j=info%ys,info%ye
          i = info%NX
          if (walls(i,j).eq.0) then
             call BCFluxApplyD2Q9(bc, fi(:,:,i,j), xp_vals(:,j), &
                  directions, cardinals, info)
          end if
       end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_PRESSURE).and.(info%ys.eq.1)) then
       call BCSetLocalDirectionsD2Q9(BOUNDARY_YM, directions, cardinals)
       do i=info%xs,info%xe
          j = 1
          if (walls(i,j).eq.0) then
             call BCFluxApplyD2Q9(bc, fi(:,:,i,j), ym_vals(:,i), &
                  directions, cardinals, info)
          end if
       end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_PRESSURE).and.(info%ye.eq.info%NY)) then
       call BCSetLocalDirectionsD2Q9(BOUNDARY_YP, directions, cardinals)
       do i=info%xs,info%xe
          j = info%NY
          if (walls(i,j).eq.0) then
             call BCFluxApplyD2Q9(bc, fi(:,:,i,j), yp_vals(:,i), &
                  directions, cardinals, info)
          end if
       end do
    endif
    return
  end subroutine BCFluxD2Q9

  subroutine BCSetLocalDirectionsD2Q9(boundary, directions, cardinals)
    use LBM_Discretization_Directions_D2Q9_module
    PetscInt,intent(in):: boundary
    PetscInt,intent(out),dimension(0:discretization_directions) :: directions
    PetscInt,intent(out),dimension(1:discretization_dims) :: cardinals

    PetscInt i
    
    select case(boundary)
    case (BOUNDARY_XP)
       ! identity mapping
       do i=0,discretization_directions
          directions(i) = i
       end do
       cardinals(CARDINAL_NORMAL) = X_DIRECTION
       cardinals(CARDINAL_CROSS) = Y_DIRECTION

    case (BOUNDARY_XM)
       ! inverted mapping around the origin
       directions(ORIGIN) = ORIGIN
       directions(EAST) = WEST
       directions(WEST) = EAST
       directions(NORTH) = SOUTH
       directions(SOUTH) = NORTH
       directions(NORTHEAST) = SOUTHWEST
       directions(SOUTHWEST) = NORTHEAST
       directions(NORTHWEST) = SOUTHEAST
       directions(SOUTHEAST) = NORTHWEST

       cardinals(CARDINAL_NORMAL) = X_DIRECTION
       cardinals(CARDINAL_CROSS) = Y_DIRECTION

    case (BOUNDARY_YP)
       ! map x -> y
       directions(ORIGIN) = ORIGIN
       directions(EAST) = NORTH
       directions(WEST) = SOUTH
       directions(NORTH) = EAST
       directions(SOUTH) = WEST
       directions(NORTHEAST) = NORTHEAST
       directions(SOUTHWEST) = SOUTHWEST
       directions(NORTHWEST) = SOUTHEAST
       directions(SOUTHEAST) = NORTHWEST

       cardinals(CARDINAL_NORMAL) = Y_DIRECTION
       cardinals(CARDINAL_CROSS) = X_DIRECTION

    case (BOUNDARY_YM)
       ! map x -> y
       directions(ORIGIN) = ORIGIN
       directions(EAST) = SOUTH
       directions(WEST) = NORTH
       directions(NORTH) = WEST
       directions(SOUTH) = EAST
       directions(NORTHEAST) = SOUTHWEST
       directions(SOUTHWEST) = NORTHEAST
       directions(NORTHWEST) = NORTHWEST
       directions(SOUTHEAST) = SOUTHEAST

       cardinals(CARDINAL_NORMAL) = Y_DIRECTION
       cardinals(CARDINAL_CROSS) = X_DIRECTION

    end select
  end subroutine BCSetLocalDirectionsD2Q9

!!!!! MY 2D Blux BC Additions !!!!!

  subroutine BCFluxApplyD2Q9(bc, fi, fvals, directions, cardinals, info)
    use LBM_Discretization_Directions_D2Q9_module

    type(bc_type) bc
    type(info_type) info

    PetscScalar,intent(inout),dimension(1:info%nphases, 0:info%flow_b):: fi
    PetscScalar,intent(in),dimension(bc%dim):: fvals
    PetscInt,intent(in),dimension(0:info%flow_b):: directions
    PetscInt,intent(in),dimension(1:info%ndims):: cardinals

    PetscScalar rhotmp
    PetscScalar,dimension(0:info%flow_b)::ftmp
    PetscInt m

    ftmp = 0.0
    rhotmp = 0.0

    !!!!! Ethan, this written for the NORTH boundary.

    do m=1,info%nphases
       rhotmp = fi(m,directions(ORIGIN)) + fi(m,directions(EAST)) + fi(m,directions(WEST)) &
            + 2.*(fi(m,directions(NORTH)) + fi(m,directions(NORTHEAST)) + fi(m,directions(NORTHWEST))) 
       rhotmp = rhotmp/(1. + fvals(cardinals(CARDINAL_NORMAL)))
       
       ! Choice should not affect the momentum significantly
       ftmp(directions(SOUTH)) = fi(m,directions(NORTH))
       ftmp(directions(SOUTHEAST)) = fi(m,directions(NORTHWEST))
       ftmp(directions(SOUTHWEST)) = fi(m,directions(NORTHEAST))

       fi(m,directions(SOUTH)) = 1./3.*ftmp(directions(SOUTH)) &
            - 2./3.*rhotmp*fvals(cardinals(CARDINAL_NORMAL)) &
            - 2./3.*(ftmp(directions(SOUTHEAST)) + ftmp(directions(SOUTHWEST))) &
            + 2./3.*(fi(m,directions(NORTH)) + fi(m,directions(NORTHEAST)) + fi(m,directions(NORTHWEST)))

       fi(m,directions(SOUTHWEST)) = 1./3.*ftmp(directions(SOUTHWEST)) &
            + 1./2.*rhotmp*fvals(cardinals(CARDINAL_CROSS)) &
            - 1./6.*rhotmp*fvals(cardinals(CARDINAL_NORMAL)) &
            + 1./2.*(fi(m,directions(EAST)) - fi(m,directions(WEST))) &
            + 1./6.*(fi(m,directions(NORTH)) - ftmp(directions(SOUTH))) &
            - 1./3.*(fi(m,directions(NORTHWEST)) - ftmp(directions(SOUTHEAST))) &
            + 2./3.*fi(m,directions(NORTHEAST))

       fi(m,directions(SOUTHEAST)) = 1./3.*ftmp(directions(SOUTHEAST)) &
            - 1./2.*rhotmp*fvals(cardinals(CARDINAL_CROSS)) &
            - 1./6.*rhotmp*fvals(cardinals(CARDINAL_NORMAL)) &
            - 1./2.*(fi(m,directions(EAST)) - fi(m,directions(WEST))) &
            + 1./6.*(fi(m,directions(NORTH)) - ftmp(directions(SOUTH))) &
            - 1./3.*(fi(m,directions(NORTHEAST)) - ftmp(directions(SOUTHWEST))) &
            + 2./3.*fi(m,directions(NORTHWEST)) 

    enddo


    return
  end subroutine BCFluxApplyD2Q9

  
  subroutine BCPressureApplyD2Q9(bc, fi, pvals, directions, info)
    use LBM_Discretization_Directions_D2Q9_module

    type(bc_type) bc
    type(info_type) info
    PetscInt,intent(in),dimension(0:info%flow_b):: directions
    PetscScalar,intent(inout),dimension(1:info%nphases, 0:info%flow_b):: fi
    PetscScalar,intent(in),dimension(bc%dim):: pvals

!    PetscScalar utmp
    PetscScalar vtmp
    PetscScalar,dimension(0:info%flow_b)::ftmp
    integer m
    
!    utmp = 0
    vtmp = 0

    !!!!! Ethan, this written for the NORTH boundary.

    do m=1,info%nphases
       ftmp = 0.0
       vtmp = fi(m,directions(ORIGIN)) + fi(m,directions(EAST)) + fi(m,directions(WEST)) &
            + 2.*(fi(m,directions(NORTH)) + fi(m,directions(NORTHEAST)) + fi(m,directions(NORTHWEST)))
       vtmp = vtmp/pvals(m)-1.0       
       
       ! Choice should not affect the momentum significantly
       ftmp(directions(SOUTH)) = fi(m,directions(NORTH))
       ftmp(directions(SOUTHWEST)) = fi(m,directions(NORTHEAST))
       ftmp(directions(SOUTHEAST)) = fi(m,directions(NORTHWEST))
              
       fi(m,directions(SOUTH)) = 1./3.*ftmp(directions(SOUTH)) &
            - 2./3.*pvals(m)*vtmp &
            - 2./3.*(ftmp(directions(SOUTHWEST)) + ftmp(directions(SOUTHEAST))) &
            + 2./3.*(fi(m,directions(NORTH)) + fi(m,directions(NORTHEAST)) + fi(m,directions(NORTHWEST)))
       
       fi(m,directions(SOUTHWEST)) = 1./3.*ftmp(directions(SOUTHWEST)) &
!            - 1./2.*pvals(m)*utmp &
            - 1./6.*pvals(m)*vtmp &
            + 1./2.*(fi(m,directions(EAST)) - fi(m,directions(WEST))) &
            + 1./6.*(fi(m,directions(NORTH)) - ftmp(directions(SOUTH))) &
            - 1./3.*(fi(m,directions(NORTHWEST)) - ftmp(directions(SOUTHEAST))) &
            + 2./3.*fi(m,directions(NORTHEAST))
       
       fi(m,directions(SOUTHEAST)) = 1./3.*ftmp(directions(SOUTHEAST)) &
!            + 1./2.*pvals(m)*utmp &
            - 1./6.*pvals(m)*vtmp &
            - 1./2.*(fi(m,directions(EAST)) - fi(m,directions(WEST))) &
            + 1./6.*(fi(m,directions(NORTH)) - ftmp(directions(SOUTH))) &
            - 1./3.*(fi(m,directions(NORTHEAST)) - ftmp(directions(SOUTHWEST))) &
            + 2./3.*fi(m,directions(NORTHWEST))
       
    enddo
    
    return
  end subroutine BCPressureApplyD2Q9

end module LBM_BC_module
