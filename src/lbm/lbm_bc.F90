!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        bc.F90
!!!     version:         
!!!     created:         06 December 2010
!!!       on:            09:03:18 MST
!!!     last modified:   14 September 2011
!!!       at:            12:32:20 PDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ ldeo.columbia.edu
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module LBM_BC_module
  use LBM_Error_module
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
       BCGetArrays, &
       BCRestoreArrays, &
       BCUpdateRho, &
       BCPreStream, &
       BCApplyDirichletToRho, &
       BCApply
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
    bc%flags = BC_NULL
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

  subroutine BCSetFromOptions(bc, options, ierr)
    type(bc_type) bc
    type(options_type) options
    PetscErrorCode ierr

    ! nothing to do, BCs get set from their containing physics module
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

    ! z boundaries
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

  subroutine BCApplyDirichletToRho(bc, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%info%rgxyzl):: walls

    PetscInt lcv_sides
    do lcv_sides = 1,6
      if (bc%flags(lcv_sides).eq.BC_DIRICHLET) then
        select case(dist%info%ndims)
        case (2)
          call BCApplyDirichletToRho_D2(bc, walls, dist%rho_a, dist, lcv_sides, &
               bc%xm_a, bc%xp_a, bc%ym_a, bc%yp_a)
        case (3)
          call BCApplyDirichletToRho_D3(bc, walls, dist%rho_a, dist, lcv_sides, &
               bc%xm_a, bc%xp_a, bc%ym_a, bc%yp_a, bc%zm_a, bc%zp_a)
        end select
      end if
    end do
  end subroutine BCApplyDirichletToRho

  subroutine BCApplyDirichletToRho_D2(bc, walls, rho, dist, boundary, &
       xm_vals, xp_vals, ym_vals, yp_vals)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls
    PetscScalar,dimension(dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscInt boundary
    PetscScalar,dimension(bc%nbcs,dist%info%ys:dist%info%ye):: xm_vals, xp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe):: ym_vals, yp_vals

    PetscInt m

    ! note this requires m is fastest varying in nbcs
    select case(boundary)
    case(BOUNDARY_XM)
      if (dist%info%xs.eq.1) then
        do m=1,dist%s
          rho(m,dist%info%xs,dist%info%ys:dist%info%ye) = &
               xm_vals(m,dist%info%ys:dist%info%ye)
        end do
      end if
    case(BOUNDARY_XP)
      if (dist%info%xe.eq.dist%info%NX) then
        do m=1,dist%s
          rho(m,dist%info%xe,dist%info%ys:dist%info%ye) = &
               xp_vals(m,dist%info%ys:dist%info%ye)
        end do
      end if
    case(BOUNDARY_YM)
      if (dist%info%ys.eq.1) then
        do m=1,dist%s
          rho(m,dist%info%xs:dist%info%xe,dist%info%ys) = &
               ym_vals(m,dist%info%xs:dist%info%xe)
        end do
      end if
    case(BOUNDARY_YP)
      if (dist%info%ye.eq.dist%info%NY) then
        do m=1,dist%s
          rho(m,dist%info%xs:dist%info%xe,dist%info%ye) = &
               yp_vals(m,dist%info%xs:dist%info%xe)
        end do
      end if
    end select
  end subroutine BCApplyDirichletToRho_D2

  subroutine BCApplyDirichletToRho_D3(bc, walls, rho, dist, boundary, &
       xm_vals, xp_vals, ym_vals, yp_vals, zm_vals, zp_vals)
    type(bc_type) bc
    type(distribution_type) dist

    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: walls
    PetscScalar,dimension(dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: rho
    PetscInt boundary
    PetscScalar,dimension(bc%nbcs,dist%info%ys:dist%info%ye, &
         dist%info%zs:dist%info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe, &
         dist%info%zs:dist%info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%ys:dist%info%ye, &
         dist%info%zs:dist%info%ze):: zm_vals, zp_vals

    PetscInt m

    select case(boundary)
    case(BOUNDARY_XM)
      if (dist%info%xs.eq.1) then
        do m=1,dist%s
          rho(m,dist%info%xs,dist%info%ys:dist%info%ye,dist%info%zs:dist%info%ze) = &
               xm_vals(m,dist%info%ys:dist%info%ye,dist%info%zs:dist%info%ze)
        end do
      end if
    case(BOUNDARY_XP)
      if (dist%info%xe.eq.dist%info%NX) then
        do m=1,dist%s
          rho(m,dist%info%xe,dist%info%ys:dist%info%ye,dist%info%zs:dist%info%ze) = &
               xp_vals(m,dist%info%ys:dist%info%ye,dist%info%zs:dist%info%ze)
        end do
      end if
    case(BOUNDARY_YM)
      if (dist%info%ys.eq.1) then
        do m=1,dist%s
          rho(m,dist%info%xs:dist%info%xe,dist%info%ys,dist%info%zs:dist%info%ze) = &
               ym_vals(m,dist%info%xs:dist%info%xe,dist%info%zs:dist%info%ze)
        end do
      end if
    case(BOUNDARY_YP)
      if (dist%info%ye.eq.dist%info%NY) then
        do m=1,dist%s
          rho(m,dist%info%xs:dist%info%xe,dist%info%ye,dist%info%zs:dist%info%ze) = &
               yp_vals(m,dist%info%xs:dist%info%xe,dist%info%zs:dist%info%ze)
        end do
      end if
    case(BOUNDARY_ZM)
      if (dist%info%zs.eq.1) then
        do m=1,dist%s
          rho(m,dist%info%xs:dist%info%xe,dist%info%ys:dist%info%ye,dist%info%zs) = &
               ym_vals(m,dist%info%xs:dist%info%xe,dist%info%ys:dist%info%ye)
        end do
      end if
    case(BOUNDARY_ZP)
      if (dist%info%ze.eq.dist%info%NZ) then
        do m=1,dist%s
          rho(m,dist%info%xs:dist%info%xe,dist%info%ys:dist%info%ye,dist%info%ze) = &
               yp_vals(m,dist%info%xs:dist%info%xe,dist%info%ys:dist%info%ye)
        end do
      end if
    end select
  end subroutine BCApplyDirichletToRho_D3

  subroutine BCUpdateRho(bc, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%info%rgxyzl):: walls

    PetscInt lcv_sides

    do lcv_sides = 1,6
      if (bc%flags(lcv_sides).eq.BC_NEUMANN .or. &
           bc%flags(lcv_sides).eq.BC_VELOCITY .or. &
           bc%flags(lcv_sides).eq.BC_DIRICHLET) then

        select case(dist%info%ndims)
        case (2)
          call BCUpdateRho_D2(bc, walls, dist%fi_a, dist%rho_a, dist, lcv_sides)
        case (3)
          call BCUpdateRho_D3(bc, walls, dist%fi_a, dist%rho_a, dist, lcv_sides)
        end select
      end if
    enddo
  end subroutine BCUpdateRho

  subroutine BCUpdateRho_D2(bc, walls, fi, rho, dist, boundary)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe,dist%info%rgys:dist%info%rgye):: walls
    PetscScalar,dimension(dist%s,0:dist%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscInt boundary

    PetscInt i,j,m

    select case(boundary)
    case(BOUNDARY_XM)
      if (dist%info%xs.eq.1) then
        i = 1
        do j=dist%info%ys,dist%info%ye
          if (walls(i,j).eq.0.d0) then
            do m=1,dist%s
              rho(m,i,j) = sum(fi(m,:,i,j))
            end do
          end if
        end do
      end if
    case(BOUNDARY_XP)
      if (dist%info%xe.eq.dist%info%NX) then
        i = dist%info%NX
        do j=dist%info%ys,dist%info%ye
          if (walls(i,j).eq.0.d0) then
            do m=1,dist%s
              rho(m,i,j) = sum(fi(m,:,i,j))
            end do
          end if
        end do
      end if
    case(BOUNDARY_YM)
      if (dist%info%ys.eq.1) then
        j = 1
        do i=dist%info%xs,dist%info%xe
          if (walls(i,j).eq.0.d0) then
            do m=1,dist%s
              rho(m,i,j) = sum(fi(m,:,i,j))
            end do
          end if
        end do
      end if
    case(BOUNDARY_YP)
      if (dist%info%ye.eq.dist%info%NY) then
        j = dist%info%NY
        do i=dist%info%xs,dist%info%xe
          if (walls(i,j).eq.0.d0) then
            do m=1,dist%s
              rho(m,i,j) = sum(fi(m,:,i,j))
            end do
          end if
        end do
      end if
    end select
  end subroutine BCUpdateRho_D2

  subroutine BCUpdateRho_D3(bc, walls, fi, rho, dist, boundary)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye,dist%info%rgzs:dist%info%rgze):: walls
    PetscScalar,dimension(dist%s,0:dist%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye,dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye,dist%info%rgzs:dist%info%rgze):: rho
    PetscInt boundary

    PetscInt i,j,k,m

    select case(boundary)
    case(BOUNDARY_XM)
      if (dist%info%xs.eq.1) then
        i = 1
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          if (walls(i,j,k).eq.0.d0) then
            do m=1,dist%s
              rho(m,i,j,k) = sum(fi(m,:,i,j,k))
            end do
          end if
        end do
        end do
      end if
    case(BOUNDARY_XP)
      if (dist%info%xe.eq.dist%info%NX) then
        i = dist%info%NX
        do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          if (walls(i,j,k).eq.0.d0) then
            do m=1,dist%s
              rho(m,i,j,k) = sum(fi(m,:,i,j,k))
            end do
          end if
        end do
        end do
      end if
    case(BOUNDARY_YM)
      if (dist%info%ys.eq.1) then
        j = 1
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          if (walls(i,j,k).eq.0.d0) then
            do m=1,dist%s
              rho(m,i,j,k) = sum(fi(m,:,i,j,k))
            end do
          end if
        end do
        end do
      end if
    case(BOUNDARY_YP)
      if (dist%info%ye.eq.dist%info%NY) then
        j = dist%info%NY
        do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          if (walls(i,j,k).eq.0.d0) then
            do m=1,dist%s
              rho(m,i,j,k) = sum(fi(m,:,i,j,k))
            end do
          end if
        end do
        end do
      end if
    case(BOUNDARY_ZM)
      if (dist%info%zs.eq.1) then
        k = 1
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          if (walls(i,j,k).eq.0.d0) then
            do m=1,dist%s
              rho(m,i,j,k) = sum(fi(m,:,i,j,k))
            end do
          end if
        end do
        end do
      end if
    case(BOUNDARY_ZP)
      if (dist%info%ze.eq.dist%info%NZ) then
        k = dist%info%NZ
        do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          if (walls(i,j,k).eq.0.d0) then
            do m=1,dist%s
              rho(m,i,j,k) = sum(fi(m,:,i,j,k))
            end do
          end if
        end do
        end do
      end if
    end select
  end subroutine BCUpdateRho_D3

  subroutine BCPreStream(bc, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscInt lcv_sides

    do lcv_sides = 1,6
      if (bc%flags(lcv_sides).eq.BC_NEUMANN .or. &
           bc%flags(lcv_sides).eq.BC_VELOCITY .or. &
           bc%flags(lcv_sides).eq.BC_DIRICHLET) then

        select case(dist%info%ndims)
        case (2)
          call BCPreStream_D2(bc, dist%fi_a, dist, lcv_sides)
        case (3)
          call BCPreStream_D3(bc, dist%fi_a, dist, lcv_sides)
        end select
      end if
    enddo
  end subroutine BCPreStream

  subroutine BCPreStream_D2(bc, fi, dist, boundary)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,0:dist%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscInt boundary

    PetscInt n

    select case(boundary)
    case(BOUNDARY_XM)
      if (dist%info%xs.eq.1) then
        do n=1,dist%b
          if (dist%disc%ci(n,X_DIRECTION) > 0) then
            fi(:,n,dist%info%xs-dist%disc%ci(n,X_DIRECTION), &
                 dist%info%ys-dist%disc%ci(n,Y_DIRECTION):&
                 dist%info%ye-dist%disc%ci(n,Y_DIRECTION)) = &
                 fi(:,n,dist%info%xs,dist%info%ys:dist%info%ye)
          end if
        end do
      end if
    case(BOUNDARY_XP)
      if (dist%info%xe.eq.dist%info%NX) then
        do n=1,dist%b
          if (dist%disc%ci(n,X_DIRECTION) < 0) then
            fi(:,n,dist%info%xe-dist%disc%ci(n,X_DIRECTION), &
                 dist%info%ys-dist%disc%ci(n,Y_DIRECTION):&
                 dist%info%ye-dist%disc%ci(n,Y_DIRECTION)) = &
                 fi(:,n,dist%info%xe,dist%info%ys:dist%info%ye)
          end if
        end do
      end if
    case(BOUNDARY_YM)
      if (dist%info%ys.eq.1) then
        do n=1,dist%b
          if (dist%disc%ci(n,Y_DIRECTION) > 0) then
            fi(:,n,dist%info%xs-dist%disc%ci(n,X_DIRECTION):&
                 dist%info%xe-dist%disc%ci(n,X_DIRECTION), &
                 dist%info%ys-dist%disc%ci(n,Y_DIRECTION)) = &
                 fi(:,n,dist%info%xs:dist%info%xe,dist%info%ys)
          end if
        end do
      end if
    case(BOUNDARY_YP)
      if (dist%info%ye.eq.dist%info%NY) then
        do n=1,dist%b
          if (dist%disc%ci(n,Y_DIRECTION) < 0) then
            fi(:,n,dist%info%xs-dist%disc%ci(n,X_DIRECTION):&
                 dist%info%xe-dist%disc%ci(n,X_DIRECTION), &
                 dist%info%ye-dist%disc%ci(n,Y_DIRECTION)) = &
                 fi(:,n,dist%info%xs:dist%info%xe,dist%info%ye)
          end if
        end do
      end if
    end select
  end subroutine BCPreStream_D2

  subroutine BCPreStream_D3(bc, fi, dist, boundary)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,0:dist%b,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye,dist%info%gzs:dist%info%gze):: fi
    PetscInt boundary

    PetscInt n

    select case(boundary)
    case(BOUNDARY_XM)
      if (dist%info%xs.eq.1) then
        do n=1,dist%b
          if (dist%disc%ci(n,X_DIRECTION) > 0) then
            fi(:,n,dist%info%xs-dist%disc%ci(n,X_DIRECTION), &
                 dist%info%ys-dist%disc%ci(n,Y_DIRECTION):&
                 dist%info%ye-dist%disc%ci(n,Y_DIRECTION),&
                 dist%info%zs-dist%disc%ci(n,Z_DIRECTION):&
                 dist%info%ze-dist%disc%ci(n,Z_DIRECTION)) = &
                 fi(:,n,dist%info%xs,dist%info%ys:dist%info%ye,dist%info%zs:dist%info%ze)
          end if
        end do
      end if
    case(BOUNDARY_XP)
      if (dist%info%xe.eq.dist%info%NX) then
        do n=1,dist%b
          if (dist%disc%ci(n,X_DIRECTION) < 0) then
            fi(:,n,dist%info%xe-dist%disc%ci(n,X_DIRECTION), &
                 dist%info%ys-dist%disc%ci(n,Y_DIRECTION):&
                 dist%info%ye-dist%disc%ci(n,Y_DIRECTION),&
                 dist%info%zs-dist%disc%ci(n,Z_DIRECTION):&
                 dist%info%ze-dist%disc%ci(n,Z_DIRECTION)) = &
                 fi(:,n,dist%info%xe,dist%info%ys:dist%info%ye,dist%info%zs:dist%info%ze)
          end if
        end do
      end if
    case(BOUNDARY_YM)
      if (dist%info%ys.eq.1) then
        do n=1,dist%b
          if (dist%disc%ci(n,Y_DIRECTION) > 0) then
            fi(:,n,dist%info%xs-dist%disc%ci(n,X_DIRECTION):&
                 dist%info%xe-dist%disc%ci(n,X_DIRECTION), &
                 dist%info%ys-dist%disc%ci(n,Y_DIRECTION), &
                 dist%info%zs-dist%disc%ci(n,Z_DIRECTION):&
                 dist%info%ze-dist%disc%ci(n,Z_DIRECTION)) = &
                 fi(:,n,dist%info%xs:dist%info%xe,dist%info%ys,dist%info%zs:dist%info%ze)
          end if
        end do
      end if
    case(BOUNDARY_YP)
      if (dist%info%ye.eq.dist%info%NY) then
        do n=1,dist%b
          if (dist%disc%ci(n,Y_DIRECTION) < 0) then
            fi(:,n,dist%info%xs-dist%disc%ci(n,X_DIRECTION):&
                 dist%info%xe-dist%disc%ci(n,X_DIRECTION), &
                 dist%info%ye-dist%disc%ci(n,Y_DIRECTION), &
                 dist%info%zs-dist%disc%ci(n,Z_DIRECTION):&
                 dist%info%ze-dist%disc%ci(n,Z_DIRECTION)) = &
                 fi(:,n,dist%info%xs:dist%info%xe,dist%info%ye,dist%info%zs:dist%info%ze)
          end if
        end do
      end if
    case(BOUNDARY_ZM)
      if (dist%info%zs.eq.1) then
        do n=1,dist%b
          if (dist%disc%ci(n,Z_DIRECTION) > 0) then
            fi(:,n,dist%info%xs-dist%disc%ci(n,X_DIRECTION):&
                 dist%info%xe-dist%disc%ci(n,X_DIRECTION), &
                 dist%info%ys-dist%disc%ci(n,Y_DIRECTION):&
                 dist%info%ye-dist%disc%ci(n,Y_DIRECTION), &
                 dist%info%zs-dist%disc%ci(n,Z_DIRECTION)) = &
                 fi(:,n,dist%info%xs:dist%info%xe,dist%info%ys:dist%info%ye,dist%info%zs)
          end if
        end do
      end if
    case(BOUNDARY_ZP)
      if (dist%info%ze.eq.dist%info%NZ) then
        do n=1,dist%b
          if (dist%disc%ci(n,Z_DIRECTION) < 0) then
            fi(:,n,dist%info%xs-dist%disc%ci(n,X_DIRECTION):&
                 dist%info%xe-dist%disc%ci(n,X_DIRECTION), &
                 dist%info%ys-dist%disc%ci(n,Y_DIRECTION):&
                 dist%info%ye-dist%disc%ci(n,Y_DIRECTION), &
                 dist%info%ze-dist%disc%ci(n,Z_DIRECTION)) = &
                 fi(:,n,dist%info%xs:dist%info%xe,dist%info%ys:dist%info%ye,dist%info%ze)
          end if
        end do
      end if
    end select
  end subroutine BCPreStream_D3

  subroutine BCApply(bc, forces, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%info%rgxyzl):: walls
    PetscScalar,dimension(dist%s,dist%info%ndims,dist%info%gxyzl):: forces

    logical,dimension(0:15):: bcs_done
    PetscInt lcv_sides
    bcs_done=.FALSE.
    bcs_done(BC_PERIODIC) = .TRUE.   ! periodic done by default

    do lcv_sides = 1,6
      if (.not.bcs_done(bc%flags(lcv_sides))) then
        select case (bc%flags(lcv_sides))
        case (BC_REFLECTING)
          call BCApplyReflecting(bc, walls, dist)
        case (BC_DIRICHLET)
          call BCApplyDirichlet(bc, forces, walls, dist)
        case (BC_NEUMANN)
          call BCApplyNeumann(bc, forces, walls, dist)
        case (BC_VELOCITY)
          call BCApplyVelocity(bc, forces, walls, dist)
        end select
        bcs_done(bc%flags(lcv_sides)) = .TRUE. ! only do each bc type once
      endif
    enddo
  end subroutine BCApply

  subroutine BCApplyReflecting(bc, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%info%rgxyzl):: walls
    PetscErrorCode ierr

    select case(dist%info%ndims)
    case(3)
      call BCApplyReflectingD3(bc, dist%fi_a, walls, dist)
    case(D2Q9_DISCRETIZATION)
      call BCApplyReflectingD2(bc, dist%fi_a, walls, dist)
    case DEFAULT
      call LBMError(PETSC_COMM_SELF, 1, 'invalid discretization in LBM', ierr)
    end select
  end subroutine BCApplyReflecting

  subroutine BCApplyReflectingD3(bc, fi, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: walls
    PetscErrorCode ierr

    PetscInt i,j,k,n,p
    PetscScalar tmp

    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_REFLECTING).and.(dist%info%xs.eq.1)) then
      do k=dist%info%zs,dist%info%ze
      do j=dist%info%ys,dist%info%ye
        i = 1
        if (walls(i,j,k).eq.0) then
          do n=1,dist%b
            if (dist%disc%ci(n,X_DIRECTION) > 0) then
              ! has a component in the inward-normal direction
              do p=1,dist%b
                if ((dist%disc%ci(n,Y_DIRECTION).eq.dist%disc%ci(p,Y_DIRECTION)).and. &
                    (dist%disc%ci(n,Z_DIRECTION).eq.dist%disc%ci(p,Z_DIRECTION)).and. &
                    (dist%disc%ci(n,X_DIRECTION).eq.-dist%disc%ci(p,Z_DIRECTION))) then
                  fi(:,n,i,j,k) = fi(:,p,i,j,k)
                end if
              end do
            end if
          end do
        end if
      end do
      end do
    end if

    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_REFLECTING).and. &
         (dist%info%xe.eq.dist%info%NX)) then
      do k=dist%info%zs,dist%info%ze
      do j=dist%info%ys,dist%info%ye
        i = dist%info%NX
        if (walls(i,j,k).eq.0) then
          do n=1,dist%b
            if (dist%disc%ci(n,X_DIRECTION) < 0) then
              ! has a component in the inward-normal direction
              do p=1,dist%b
                if ((dist%disc%ci(n,Y_DIRECTION).eq.dist%disc%ci(p,Y_DIRECTION)).and. &
                    (dist%disc%ci(n,Z_DIRECTION).eq.dist%disc%ci(p,Z_DIRECTION)).and. &
                    (dist%disc%ci(n,X_DIRECTION).eq.-dist%disc%ci(p,X_DIRECTION))) then
                  fi(:,n,i,j,k) = fi(:,p,i,j,k)
                end if
              end do
            end if
          end do
        end if
      end do
      end do
    end if

    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_REFLECTING).and.(dist%info%ys.eq.1)) then
      do k=dist%info%zs,dist%info%ze
      do i=dist%info%xs,dist%info%xe
        j = 1
        if (walls(i,j,k).eq.0) then
          do n=1,dist%b
            if (dist%disc%ci(n,Y_DIRECTION) > 0) then
              ! has a component in the inward-normal direction
              do p=1,dist%b
                if ((dist%disc%ci(n,X_DIRECTION).eq.dist%disc%ci(p,X_DIRECTION)).and. &
                    (dist%disc%ci(n,Z_DIRECTION).eq.dist%disc%ci(p,Z_DIRECTION)).and. &
                    (dist%disc%ci(n,Y_DIRECTION).eq.-dist%disc%ci(p,Y_DIRECTION))) then
                  fi(:,n,i,j,k) = fi(:,p,i,j,k)
                end if
              end do
            end if
          end do
        end if
      end do
      end do
    end if

    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_REFLECTING).and. &
         (dist%info%ye.eq.dist%info%NY)) then
      do k=dist%info%zs,dist%info%ze
      do i=dist%info%xs,dist%info%xe
        j = dist%info%NY
        if (walls(i,j,k).eq.0) then
          do n=1,dist%b
            if (dist%disc%ci(n,Y_DIRECTION) < 0) then
              ! has a component in the inward-normal direction
              do p=1,dist%b
                if ((dist%disc%ci(n,X_DIRECTION).eq.dist%disc%ci(p,X_DIRECTION)).and. &
                    (dist%disc%ci(n,Z_DIRECTION).eq.dist%disc%ci(p,Z_DIRECTION)).and. &
                    (dist%disc%ci(n,Y_DIRECTION).eq.-dist%disc%ci(p,Y_DIRECTION))) then
                  fi(:,n,i,j,k) = fi(:,p,i,j,k)
                end if
              end do
            end if
          end do
        end if
      end do
      end do
    end if

    ! ZM BOUNDARY
    if ((bc%flags(BOUNDARY_ZM).eq.BC_REFLECTING).and.(dist%info%zs.eq.1)) then
      do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe
        k = 1
        if (walls(i,j,k).eq.0) then
          do n=1,dist%b
            if (dist%disc%ci(n,Z_DIRECTION) > 0) then
              ! has a component in the inward-normal direction
              do p=1,dist%b
                if ((dist%disc%ci(n,X_DIRECTION).eq.dist%disc%ci(p,X_DIRECTION)).and. &
                    (dist%disc%ci(n,Y_DIRECTION).eq.dist%disc%ci(p,Y_DIRECTION)).and. &
                    (dist%disc%ci(n,Z_DIRECTION).eq.-dist%disc%ci(p,Z_DIRECTION))) then
                  fi(:,n,i,j,k) = fi(:,p,i,j,k)
                end if
              end do
            end if
          end do
        end if
      end do
      end do
    end if

    ! ZP BOUNDARY
    if ((bc%flags(BOUNDARY_ZP).eq.BC_REFLECTING).and. &
         (dist%info%ze.eq.dist%info%NZ)) then
      do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe
        k = dist%info%NZ
        if (walls(i,j,k).eq.0) then
          do n=1,dist%b
            if (dist%disc%ci(n,Z_DIRECTION) < 0) then
              ! has a component in the inward-normal direction
              do p=1,dist%b
                if ((dist%disc%ci(n,X_DIRECTION).eq.dist%disc%ci(p,X_DIRECTION)).and. &
                    (dist%disc%ci(n,Y_DIRECTION).eq.dist%disc%ci(p,Y_DIRECTION)).and. &
                    (dist%disc%ci(n,Z_DIRECTION).eq.-dist%disc%ci(p,Z_DIRECTION))) then
                  fi(:,n,i,j,k) = fi(:,p,i,j,k)
                end if
              end do
            end if
          end do
        end if
      end do
      end do
    end if
  end subroutine BCApplyReflectingD3

  subroutine BCApplyReflectingD2(bc, fi, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls
    PetscErrorCode ierr

    PetscInt i,j,n,p
    PetscScalar tmp


    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_REFLECTING).and.(dist%info%xs.eq.1)) then
      do j=dist%info%ys,dist%info%ye
        i = 1
        if (walls(i,j).eq.0) then
          do n=1,dist%b
            if (dist%disc%ci(n,X_DIRECTION) > 0) then
              ! has a component in the inward-normal direction
              do p=1,dist%b
                if ((dist%disc%ci(n,Y_DIRECTION).eq.dist%disc%ci(p,Y_DIRECTION)).and. &
                    (dist%disc%ci(n,X_DIRECTION).eq.-dist%disc%ci(p,Z_DIRECTION))) then
                  fi(:,n,i,j) = fi(:,p,i,j)
                end if
              end do
            end if
          end do
        end if
      end do
    end if

    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_REFLECTING).and. &
         (dist%info%xe.eq.dist%info%NX)) then
      do j=dist%info%ys,dist%info%ye
        i = dist%info%NX
        if (walls(i,j).eq.0) then
          do n=1,dist%b
            if (dist%disc%ci(n,X_DIRECTION) < 0) then
              ! has a component in the inward-normal direction
              do p=1,dist%b
                if ((dist%disc%ci(n,Y_DIRECTION).eq.dist%disc%ci(p,Y_DIRECTION)).and. &
                    (dist%disc%ci(n,X_DIRECTION).eq.-dist%disc%ci(p,X_DIRECTION))) then
                  fi(:,n,i,j) = fi(:,p,i,j)
                end if
              end do
            end if
          end do
        end if
      end do
    end if

    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_REFLECTING).and.(dist%info%ys.eq.1)) then
      do i=dist%info%xs,dist%info%xe
        j = 1
        if (walls(i,j).eq.0) then
          do n=1,dist%b
            if (dist%disc%ci(n,Y_DIRECTION) > 0) then
              ! has a component in the inward-normal direction
              do p=1,dist%b
                if ((dist%disc%ci(n,X_DIRECTION).eq.dist%disc%ci(p,X_DIRECTION)).and. &
                    (dist%disc%ci(n,Y_DIRECTION).eq.-dist%disc%ci(p,Y_DIRECTION))) then
                  fi(:,n,i,j) = fi(:,p,i,j)
                end if
              end do
            end if
          end do
        end if
      end do
    end if

    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_REFLECTING).and. &
         (dist%info%ye.eq.dist%info%NY)) then
      do i=dist%info%xs,dist%info%xe
        j = dist%info%NY
        if (walls(i,j).eq.0) then
          do n=1,dist%b
            if (dist%disc%ci(n,Y_DIRECTION) < 0) then
              ! has a component in the inward-normal direction
              do p=1,dist%b
                if ((dist%disc%ci(n,X_DIRECTION).eq.dist%disc%ci(p,X_DIRECTION)).and. &
                    (dist%disc%ci(n,Y_DIRECTION).eq.-dist%disc%ci(p,Y_DIRECTION))) then
                  fi(:,n,i,j) = fi(:,p,i,j)
                end if
              end do
            end if
          end do
        end if
      end do
    end if
  end subroutine BCApplyReflectingD2

  subroutine BCApplyDirichlet(bc, forces, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%info%rgxyzl):: walls
    PetscScalar,dimension(dist%s,dist%info%ndims,dist%info%gxyzl):: forces

    select case(dist%info%ndims)
    case(2)
      call BCApplyDirichletD2(bc, dist%fi_a, forces, walls, bc%xm_a, bc%xp_a, &
           bc%ym_a, bc%yp_a, dist)
    case(3)
      call BCApplyDirichletD3(bc, dist%fi_a, forces, walls, bc%xm_a, bc%xp_a, &
           bc%ym_a, bc%yp_a, bc%zm_a, bc%zp_a, dist)
    end select
  end subroutine BCApplyDirichlet

  subroutine BCApplyDirichletD3(bc, fi, forces, walls, xm_vals, xp_vals, &
       ym_vals, yp_vals, zm_vals, zp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%s,dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: walls
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
            call BCApplyDirichletNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
                 xm_vals(:,j,k), directions, cardinals, dist)
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
            call BCApplyDirichletNode(bc, fi(:,:,i,j,k), &
                 forces(:,:,i,j,k), xp_vals(:,j,k), directions, cardinals, dist)
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
            call BCApplyDirichletNode(bc, fi(:,:,i,j,k), &
                 forces(:,:,i,j,k), ym_vals(:,i,k), directions, cardinals, dist)
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
            call BCApplyDirichletNode(bc, fi(:,:,i,j,k), &
                 forces(:,:,i,j,k), yp_vals(:,i,k), directions, cardinals, dist)
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
            call BCApplyDirichletNode(bc, fi(:,:,i,j,k), &
                 forces(:,:,i,j,k), zm_vals(:,i,j), directions, cardinals, dist)
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
            call BCApplyDirichletNode(bc, fi(:,:,i,j,k), &
                 forces(:,:,i,j,k), zp_vals(:,i,j), directions, cardinals, dist)
          end if
        end do
      end do
    endif
    return
  end subroutine BCApplyDirichletD3

  subroutine BCApplyDirichletD2(bc, fi, forces, walls, xm_vals, xp_vals, ym_vals, yp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%s,dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls
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
          call BCApplyDirichletNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
               xm_vals(:,j), directions, cardinals, dist)
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
          call BCApplyDirichletNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
               xp_vals(:,j), directions, cardinals, dist)
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
          call BCApplyDirichletNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
               ym_vals(:,i), directions, cardinals, dist)
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
          call BCApplyDirichletNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
               yp_vals(:,i), directions, cardinals, dist)
        end if
      end do
    endif
  end subroutine BCApplyDirichletD2

  subroutine BCApplyDirichletNode(bc, fi, forces, pvals, directions, cardinals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,intent(inout),dimension(1:dist%s, 0:dist%disc%b):: fi
    PetscScalar,intent(in),dimension(1:dist%s, dist%info%ndims):: forces
    PetscScalar,intent(in),dimension(dist%s,dist%info%ndims):: pvals
    PetscInt,intent(in),dimension(0:dist%disc%b):: directions
    PetscInt,intent(in),dimension(1:dist%info%ndims):: cardinals

    PetscScalar, dimension(dist%info%ndims) :: Q, momentum, weightsum
    PetscInt m,n,p
    PetscInt local_normal

    ! NOTE ON ASSUMPTIONS:
    ! -- f*_i(x,t) = f_i(x,t-dt) (see Chang, Liu, and Lin, Comp & Math with Apps 2009)
    ! -- weights are isotropic (for instance, that weights(12) - weights(11) = 0)
    ! -- rho*u = F_x/2, i.e. that tangential momentum is due to forces only

    ! calculation is based upon the positive z-direction boundary, and
    ! all others are rotated to match via directions/cardinals.
    local_normal = directions(dist%disc%local_normal)

    do m=1,dist%s
      ! note fi has fi in the interior/defined directions and f* in
      ! the incoming/unknown directions
      weightsum = 0.d0
      do n=1,dist%disc%b
        if (dist%disc%ci(local_normal,cardinals(CARDINAL_NORMAL)) &
             * dist%disc%ci(n,cardinals(CARDINAL_NORMAL)).eq.1) then
          ! incoming f_i, as c_i is positive in the inward-normal direction
          weightsum(cardinals(CARDINAL_NORMAL)) = weightsum(cardinals(CARDINAL_NORMAL)) &
               + dist%disc%weights(n)
        end if
      end do

      Q(cardinals(CARDINAL_NORMAL)) = dist%disc%ci(local_normal,cardinals(CARDINAL_NORMAL))&
           * (pvals(m,1) - sum(fi(m,:))) / weightsum(cardinals(CARDINAL_NORMAL))

      momentum = 0.d0
      do p=2,dist%info%ndims
        do n=1,dist%disc%b
          momentum(cardinals(p)) = momentum(cardinals(p)) + & ! disc%ci(directions(dist%disc%LOCAL_NORMAL),cardinals(CARDINAL_NORMAL)) *
               fi(m,n)*dist%disc%ci(n, cardinals(p))
          if ((dist%disc%ci(local_normal,cardinals(CARDINAL_NORMAL)) &
               * dist%disc%ci(n,cardinals(CARDINAL_NORMAL)).eq.1).and. &
               (dist%disc%ci(n,cardinals(p)).ne.0)) then
            ! pointed inward and has a component in the p-direction
            weightsum(cardinals(p)) = weightsum(cardinals(p)) + dist%disc%weights(n)
          end if
        end do
        Q(cardinals(p)) = - momentum(cardinals(p)) / weightsum(cardinals(p))
      end do

      do n=1,dist%disc%b
        if (dist%disc%ci(local_normal,cardinals(CARDINAL_NORMAL)) &
             * dist%disc%ci(n,cardinals(CARDINAL_NORMAL)).eq.1) then
          fi(m,n) = fi(m,n) + dist%disc%weights(n)*sum(dist%disc%ci(n,:)*Q(:))
        end if
      end do
    end do
  end subroutine BCApplyDirichletNode

  subroutine BCApplyNeumann(bc, forces, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,dist%info%ndims,dist%info%gxyzl):: forces
    PetscScalar,dimension(dist%info%rgxyzl):: walls

    select case(dist%info%ndims)
    case(2)
      call BCApplyNeumannD2(bc, dist%fi_a, forces, walls, bc%xm_a, bc%xp_a, &
           bc%ym_a, bc%yp_a, dist)
    case(3)
      call BCApplyNeumannD3(bc, dist%fi_a, forces, walls, bc%xm_a, bc%xp_a, &
           bc%ym_a, bc%yp_a, bc%zm_a, bc%zp_a, dist)
    end select
  end subroutine BCApplyNeumann

  subroutine BCApplyNeumannD3(bc, fi, forces, walls, xm_vals, xp_vals, &
       ym_vals, yp_vals, zm_vals, zp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%s,dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: walls
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
    if ((bc%flags(BOUNDARY_XM).eq.BC_NEUMANN).and.(dist%info%xs.eq.1)) then
      call DiscSetLocalDirections(dist%disc, BOUNDARY_XM, directions, cardinals)
      do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          i = 1
          if (walls(i,j,k).eq.0) then
            call BCApplyNeumannNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
                 xm_vals(:,j,k), directions, cardinals, dist)
          end if
        end do
      end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_NEUMANN).and.(dist%info%xe.eq.dist%info%NX)) then
      call DiscSetLocalDirections(dist%disc, BOUNDARY_XP, directions, cardinals)
      do k=dist%info%zs,dist%info%ze
        do j=dist%info%ys,dist%info%ye
          i = dist%info%NX
          if (walls(i,j,k).eq.0) then
            call BCApplyNeumannNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
                 xp_vals(:,j,k), directions, cardinals, dist)
          end if
        end do
      end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_NEUMANN).and.(dist%info%ys.eq.1)) then
      call DiscSetLocalDirections(dist%disc, BOUNDARY_YM, directions, cardinals)
      do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          j = 1
          if (walls(i,j,k).eq.0) then
            call BCApplyNeumannNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
                 ym_vals(:,i,k), directions, cardinals, dist)
          end if
        end do
      end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_NEUMANN).and.(dist%info%ye.eq.dist%info%NY)) then
      call DiscSetLocalDirections(dist%disc, BOUNDARY_YP, directions, cardinals)
      do k=dist%info%zs,dist%info%ze
        do i=dist%info%xs,dist%info%xe
          j = dist%info%NY
          if (walls(i,j,k).eq.0) then
            call BCApplyNeumannNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
                 yp_vals(:,i,k), directions, cardinals, dist)
          end if
        end do
      end do
    endif

    directions(:) = 0
    ! ZM BOUNDARY
    if ((bc%flags(BOUNDARY_ZM).eq.BC_NEUMANN).and.(dist%info%zs.eq.1)) then
      call DiscSetLocalDirections(dist%disc, BOUNDARY_ZM, directions, cardinals)
      do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          k = 1
          if (walls(i,j,k).eq.0) then
            call BCApplyNeumannNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
                 zm_vals(:,i,j), directions, cardinals, dist)
          end if
        end do
      end do
    endif

    directions(:) = 0
    ! ZP BOUNDARY
    if ((bc%flags(BOUNDARY_ZP).eq.BC_NEUMANN).and.(dist%info%ze.eq.dist%info%NZ)) then
      call DiscSetLocalDirections(dist%disc, BOUNDARY_ZP, directions, cardinals)
      do j=dist%info%ys,dist%info%ye
        do i=dist%info%xs,dist%info%xe
          k = dist%info%NZ
          if (walls(i,j,k).eq.0) then
            call BCApplyNeumannNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
                 zp_vals(:,i,j), directions, cardinals, dist)
          end if
        end do
      end do
    endif
    return
  end subroutine BCApplyNeumannD3

  subroutine BCApplyNeumannD2(bc, fi, forces, walls, xm_vals, xp_vals, ym_vals, yp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%s,dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls
    PetscScalar,dimension(bc%nbcs,dist%info%ys:dist%info%ye):: xm_vals, xp_vals
    PetscScalar,dimension(bc%nbcs,dist%info%xs:dist%info%xe):: ym_vals, yp_vals

    PetscInt i,j
    PetscInt directions(0:dist%b)
    PetscInt cardinals(1:dist%info%ndims)
    directions(:) = 0

    ! XM BOUNDARY
    if ((bc%flags(BOUNDARY_XM).eq.BC_NEUMANN).and.(dist%info%xs.eq.1)) then
      call DiscSetLocalDirections(dist%disc, BOUNDARY_XM, directions, cardinals)
      do j=dist%info%ys,dist%info%ye
        i = 1
        if (walls(i,j).eq.0) then
          call BCApplyNeumannNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
               xm_vals(:,j), directions, cardinals, dist)
        end if
      end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc%flags(BOUNDARY_XP).eq.BC_NEUMANN).and.(dist%info%xe.eq.dist%info%NX)) then
      call DiscSetLocalDirections(dist%disc, BOUNDARY_XP, directions, cardinals)
      do j=dist%info%ys,dist%info%ye
        i = dist%info%NX
        if (walls(i,j).eq.0) then
          call BCApplyNeumannNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
               xp_vals(:,j), directions, cardinals, dist)
        end if
      end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc%flags(BOUNDARY_YM).eq.BC_NEUMANN).and.(dist%info%ys.eq.1)) then
      call DiscSetLocalDirections(dist%disc, BOUNDARY_YM, directions, cardinals)
      do i=dist%info%xs,dist%info%xe
        j = 1
        if (walls(i,j).eq.0) then
          call BCApplyNeumannNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
               ym_vals(:,i), directions, cardinals, dist)
        end if
      end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc%flags(BOUNDARY_YP).eq.BC_NEUMANN).and.(dist%info%ye.eq.dist%info%NY)) then
      call DiscSetLocalDirections(dist%disc, BOUNDARY_YP, directions, cardinals)
      do i=dist%info%xs,dist%info%xe
        j = dist%info%NY
        if (walls(i,j).eq.0) then
          call BCApplyNeumannNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
               yp_vals(:,i), directions, cardinals, dist)
        end if
      end do
    endif
  end subroutine BCApplyNeumannD2

  subroutine BCApplyNeumannNode(bc, fi, forces, mvals, directions, cardinals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,intent(inout),dimension(1:dist%s, 0:dist%disc%b):: fi
    PetscScalar,intent(in),dimension(1:dist%s, dist%info%ndims):: forces
    PetscScalar,intent(in),dimension(dist%s,dist%info%ndims):: mvals
    PetscInt,intent(in),dimension(0:dist%disc%b):: directions
    PetscInt,intent(in),dimension(1:dist%info%ndims):: cardinals

    PetscScalar, dimension(dist%info%ndims) :: Q, momentum, weightsum
    PetscInt m,n,p
    PetscInt local_normal

    ! NOTE ON ASSUMPTIONS:
    ! -- f*_i(x,t) = f_i(x,t-dt) (see Chang, Liu, and Lin, Comp & Math with Apps 2009)
    ! -- weights are isotropic (for instance, that weights(12) - weights(11) = 0)

    ! calculation is based upon the positive z-direction boundary, and
    ! all others are rotated to match via directions/cardinals.
    local_normal = directions(dist%disc%local_normal)

    do m=1,dist%s
      ! print*, 'momentum in:', mvals(m,:)

      ! note fi has fi in the interior/defined directions and f* in
      ! the incoming/unknown directions
      weightsum = 0.d0
      momentum = 0.d0
      do p=1,dist%info%ndims
        do n=1,dist%disc%b
          momentum(p) = momentum(p) &
               + fi(m,n)*dist%disc%ci(n, p)
          if ((dist%disc%ci(local_normal,cardinals(CARDINAL_NORMAL)) &
               * dist%disc%ci(n,cardinals(CARDINAL_NORMAL)).eq.1).and. &
               (dist%disc%ci(n,p).ne.0)) then
            ! points inward and has a component in the p direction
            weightsum(p) = weightsum(p) + dist%disc%weights(n)
          end if
        end do
        Q(p) = (mvals(m,p) - forces(m,p)/2.d0 - momentum(p)) / weightsum(p)
      end do

      do n=1,dist%disc%b
        if (dist%disc%ci(local_normal,cardinals(CARDINAL_NORMAL)) &
             * dist%disc%ci(n,cardinals(CARDINAL_NORMAL)).eq.1) then
          fi(m,n) = fi(m,n) + dist%disc%weights(n)*sum(dist%disc%ci(n,:)*Q(:))
        end if
      end do

      ! ! recalc and affirm it worked
      ! momentum(:) = 0.d0
      ! do p=1,dist%info%ndims
      !   do n=1,dist%disc%b
      !     momentum(p) = momentum(p) &
      !          + fi(m,n)*dist%disc%ci(n, p)
      !   end do
      ! end do
      ! momentum = momentum + forces(m,:)/2.d0
      ! print*, 'momentum out:', momentum
    end do
  end subroutine BCApplyNeumannNode

  subroutine BCApplyVelocity(bc, forces, walls, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,dist%info%ndims,dist%info%gxyzl):: forces
    PetscScalar,dimension(dist%info%rgxyzl):: walls

    select case(dist%info%ndims)
    case(2)
      call BCApplyVelocityD2(bc, dist%fi_a, forces, walls, bc%xm_a, bc%xp_a, &
           bc%ym_a, bc%yp_a, dist)
    case(3)
      call BCApplyVelocityD3(bc, dist%fi_a, forces, walls, bc%xm_a, bc%xp_a, &
           bc%ym_a, bc%yp_a, bc%zm_a, bc%zp_a, dist)
    end select
  end subroutine BCApplyVelocity

  subroutine BCApplyVelocityD3(bc, fi, forces, walls, xm_vals, xp_vals, &
       ym_vals, yp_vals, zm_vals, zp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%s,dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: walls
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
            call BCApplyVelocityNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
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
            call BCApplyVelocityNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
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
            call BCApplyVelocityNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
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
            call BCApplyVelocityNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
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
            call BCApplyVelocityNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
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
            call BCApplyVelocityNode(bc, fi(:,:,i,j,k), forces(:,:,i,j,k), &
                 zp_vals(:,i,j), directions, cardinals, dist)
          end if
        end do
      end do
    endif
    return
  end subroutine BCApplyVelocityD3

  subroutine BCApplyVelocityD2(bc, fi, forces, walls, xm_vals, xp_vals, ym_vals, yp_vals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,dimension(1:dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%s,dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: forces
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls
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
          call BCApplyVelocityNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
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
          call BCApplyVelocityNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
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
          call BCApplyVelocityNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
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
          call BCApplyVelocityNode(bc, fi(:,:,i,j), forces(:,:,i,j), &
               yp_vals(:,i), directions, cardinals, dist)
        end if
      end do
    endif
  end subroutine BCApplyVelocityD2

  subroutine BCApplyVelocityNode(bc, fi, forces, uvals, directions, cardinals, dist)
    type(bc_type) bc
    type(distribution_type) dist
    PetscScalar,intent(inout),dimension(1:dist%s, 0:dist%disc%b):: fi
    PetscScalar,intent(in),dimension(1:dist%s, dist%info%ndims):: forces
    PetscScalar,intent(in),dimension(dist%s,dist%info%ndims):: uvals
    PetscInt,intent(in),dimension(0:dist%disc%b):: directions
    PetscInt,intent(in),dimension(1:dist%info%ndims):: cardinals

    PetscScalar, dimension(dist%info%ndims) :: Q, momentum, weightsum
    PetscScalar rho
    PetscInt m,n,p
    PetscInt local_normal

    ! NOTE ON ASSUMPTIONS:
    ! -- f*_i(x,t) = f_i(x,t-dt) (see Chang, Liu, and Lin, Comp & Math with Apps 2009)
    ! -- weights are isotropic (for instance, that weights(12) - weights(11) = 0)

    ! calculation is based upon the positive z-direction boundary, and
    ! all others are rotated to match via directions/cardinals.
    local_normal = directions(dist%disc%local_normal)



    do m=1,dist%s
      ! note fi has fi in the interior/defined directions and f* in
      ! the incoming/unknown directions
      weightsum = 0.d0
      momentum = 0.d0
      do n=1,dist%disc%b
        momentum(cardinals(CARDINAL_NORMAL)) = momentum(cardinals(CARDINAL_NORMAL)) &
             + fi(m,n)*dist%disc%ci(n, cardinals(CARDINAL_NORMAL))
        if (dist%disc%ci(local_normal,cardinals(CARDINAL_NORMAL)) &
             * dist%disc%ci(n,cardinals(CARDINAL_NORMAL)).eq.1) then
          ! incoming f_i, as c_i is positive in the inward-normal direction
          weightsum(cardinals(CARDINAL_NORMAL)) = weightsum(cardinals(CARDINAL_NORMAL)) &
               + dist%disc%weights(n)
        end if
      end do

      Q(cardinals(CARDINAL_NORMAL)) = &
           (sum(fi(m,:))*uvals(1,cardinals(CARDINAL_NORMAL)) &
           - momentum(cardinals(CARDINAL_NORMAL)) - forces(m,cardinals(CARDINAL_NORMAL))/2.d0) &
           / (1.d0 - dist%disc%ci(local_normal, cardinals(CARDINAL_NORMAL))*uvals(1,cardinals(CARDINAL_NORMAL))) &
           / weightsum(cardinals(CARDINAL_NORMAL))

      rho = sum(fi(m,:)) &
           + weightsum(cardinals(CARDINAL_NORMAL))*Q(cardinals(CARDINAL_NORMAL))

      do p=2,dist%info%ndims
        do n=1,dist%disc%b
          momentum(cardinals(p)) = momentum(cardinals(p)) &
               + fi(m,n)*dist%disc%ci(n, cardinals(p))
          if ((dist%disc%ci(local_normal,cardinals(CARDINAL_NORMAL)) &
               * dist%disc%ci(n,cardinals(CARDINAL_NORMAL)).eq.1).and. &
               (dist%disc%ci(n,cardinals(p)).ne.0)) then
            ! pointed inward and has a component in the p-direction
            weightsum(cardinals(p)) = weightsum(cardinals(p)) + dist%disc%weights(n)
          end if
        end do
        Q(cardinals(p)) = (rho*uvals(1,cardinals(p)) - forces(m,cardinals(p))/2.d0 &
             - momentum(cardinals(p))) / weightsum(cardinals(p))
      end do

      do n=1,dist%disc%b
        if (dist%disc%ci(local_normal,cardinals(CARDINAL_NORMAL)) &
             * dist%disc%ci(n,cardinals(CARDINAL_NORMAL)).eq.1) then
          fi(m,n) = fi(m,n) + dist%disc%weights(n)*sum(dist%disc%ci(n,:)*Q(:))
        end if
      end do
    end do
  end subroutine BCApplyVelocityNode
end module LBM_BC_module
