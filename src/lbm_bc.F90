!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        bc.F90
!!!     version:         
!!!     created:         06 December 2010
!!!       on:            09:03:18 MST
!!!     last modified:   02 February 2011
!!!       at:            16:40:23 MST
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
  use petsc
  implicit none
  private
#include "lbm_definitions.h"

  type, public:: bc_type
     integer, dimension(6):: flags ! enum for boundary conditions
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
       BCSetSizes, &
       BCDestroy, &
       BCSetValues, &
       BCGetArrays, &
       BCRestoreArrays, &
       BCFlux, &
       BCPressure, &
       BCPseudoperiodic
contains

! constructor
function BCCreate() result(bc)
  type(bc_type),pointer :: bc
  allocate(bc)
  nullify(bc%xm_a)
  nullify(bc%xp_a)
  nullify(bc%ym_a)
  nullify(bc%yp_a)
  nullify(bc%zm_a)
  nullify(bc%zp_a)
  bc%flags = 0

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

! set up the vectors for holding boundary data
subroutine BCSetSizes(bc, comm, info)
  type(bc_type) bc
  MPI_Comm comm
  type(info_type) info
  PetscInt locn
  PetscErrorCode ierr

  ! x boundaries
  if (info%xs.eq.1) then
     locn = info%yl*info%zl*info%dim
  else
     locn = 0
  endif
  call VecCreateMPI(comm, locn, info%NY*info%NZ*info%dim, bc%xm, ierr)

  if (info%xe.eq.info%NX) then
     locn = info%yl*info%zl*info%dim
  else
     locn = 0
  endif
  call VecCreateMPI(comm, locn, info%NY*info%NZ*info%dim, bc%xp, ierr)

  ! y boundaries
  if (info%ys.eq.1) then
     locn = info%xl*info%zl*info%dim
  else
     locn = 0
  endif
  call VecCreateMPI(comm, locn, info%NX*info%NZ*info%dim, bc%ym, ierr)

  if (info%ye.eq.info%NY) then
     locn = info%xl*info%zl*info%dim
  else
     locn = 0
  endif
  call VecCreateMPI(comm, locn, info%NX*info%NZ*info%dim, bc%yp, ierr)

  ! z boundaries
  if (info%zs.eq.1) then
     locn = info%xl*info%yl*info%dim
  else
     locn = 0
  endif
  call VecCreateMPI(comm, locn, info%NX*info%NY*info%dim, bc%zm, ierr)

  if (info%ze.eq.info%NZ) then
     locn = info%xl*info%yl*info%dim
  else
     locn = 0
  endif
  call VecCreateMPI(comm, locn, info%NX*info%NY*info%dim, bc%zp, ierr)
end subroutine BCSetSizes

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
  call bc_subroutine(bc%xm_a, bc%xp_a, bc%ym_a, bc%yp_a, bc%zm_a, bc%zp_a, info%dim, info)
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
  call VecGetArrayF90(bc%zm, bc%zm_a, ierr)
  call VecGetArrayF90(bc%zp, bc%zp_a, ierr)
end subroutine BCGetArrays

subroutine BCRestoreArrays(bc, ierr)
  type(bc_type) bc
  PetscErrorCode ierr
  call VecRestoreArrayF90(bc%xm, bc%xm_a, ierr)
  call VecRestoreArrayF90(bc%xp, bc%xp_a, ierr)
  call VecRestoreArrayF90(bc%ym, bc%ym_a, ierr)
  call VecRestoreArrayF90(bc%yp, bc%yp_a, ierr)
  call VecRestoreArrayF90(bc%zm, bc%zm_a, ierr)
  call VecRestoreArrayF90(bc%zp, bc%zp_a, ierr)
end subroutine BCRestoreArrays
  
  subroutine BCFlux(fi, walls, bc_flags, bc_dim, &
       xm_vals, xp_vals, ym_vals, yp_vals, zm_vals, zp_vals, info)
    type(info_type) info
    integer bc_dim


    PetscInt,dimension(1:6)::bc_flags
    PetscScalar,dimension(1:info%s, 0:info%b, info%gxs:info%gxe, &
         info%gys:info%gye, info%gzs:info%gze):: fi
    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: walls
    PetscScalar,dimension(bc_dim,info%ys:info%ye,info%zs:info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(bc_dim,info%xs:info%xe,info%zs:info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(bc_dim,info%xs:info%xe,info%ys:info%ye):: zm_vals, zp_vals

    PetscScalar rhotmp
    PetscScalar,dimension(0:info%b)::ftmp
    integer i,j,k,m

    ftmp = 0.0

    ! In the z-direction only for now
    ! --- zp side
    if (bc_flags(6).eq.2) then
       if(info%ze.eq.info%NZ) then
          do m=1,info%s
             do i=info%xs,info%xe
                do j=info%ys,info%ye
                   k = info%NZ
                   if (walls(i,j,k).eq.0) then
                      rhotmp = fi(m,0,i,j,k) + fi(m,1,i,j,k) &
                           + fi(m,2,i,j,k) + fi(m,3,i,j,k) &
                           + fi(m,4,i,j,k) + fi(m,7,i,j,k) &
                           + fi(m,8,i,j,k) + fi(m,9,i,j,k) &
                           + fi(m,10,i,j,k) + 2.*(fi(m,5,i,j,k) &
                           + fi(m,11,i,j,k) + fi(m,12,i,j,k) &
                           + fi(m,15,i,j,k) + fi(m,16,i,j,k)) 
                      rhotmp = rhotmp/(1. + zp_vals(3,i,j))

                      ! Choice should not affect the momentum significantly
                      ftmp(6) = fi(m,5,i,j,k)
                      ftmp(13) = fi(m,11,i,j,k)
                      ftmp(14) = fi(m,12,i,j,k)
                      ftmp(17) = fi(m,15,i,j,k)
                      ftmp(18) = fi(m,16,i,j,k)
                      
                      ! This choice does not work for pressure bc
                      !ftmp(6) = fi(m,6,i,j,k-1)
                      !ftmp(13) = fi(m,13,i-1,j,k-1)
                      !ftmp(14) = fi(m,14,i+1,j,k-1)
                      !ftmp(17) = fi(m,17,i,j-1,k-1)
                      !ftmp(18) = fi(m,18,i,j+1,k-1)
                      
                      fi(m,6,i,j,k) = 2./3.*ftmp(6) &
                           - 1./3.*rhotmp*zp_vals(3,i,j) &
                           + 1./3.*(fi(m,5,i,j,k) &
                           + fi(m,11,i,j,k) + fi(m,12,i,j,k) &
                           + fi(m,15,i,j,k) + fi(m,16,i,j,k)) &
                           - 1./3.*(ftmp(13) + ftmp(14) &
                           + ftmp(17) + ftmp(18))

                      fi(m,13,i,j,k) = 1./3.*ftmp(13) &
                           - 1./2.*rhotmp*zp_vals(1,i,j) &
                           - 1./6.*rhotmp*zp_vals(3,i,j) &
                           + 1./2.*(fi(m,1,i,j,k) & 
                           - fi(m,3,i,j,k) + fi(m,7,i,j,k) &
                           - fi(m,8,i,j,k) - fi(m,9,i,j,k) &
                           + fi(m,10,i,j,k)) &
                           + 1./6.*( fi(m,5,i,j,k) - ftmp(6) &
                           + fi(m,15,i,j,k) + fi(m,16,i,j,k) &
                           - ftmp(17) - ftmp(18) ) &
                           - 1./3.*(fi(m,12,i,j,k) - ftmp(14)) &
                           + 2./3.*fi(m,11,i,j,k) 
                      
                      fi(m,14,i,j,k) = 1./3.*ftmp(14) &
                           + 1./2.*rhotmp*zp_vals(1,i,j) &
                           - 1./6.*rhotmp*zp_vals(3,i,j) &
                           - 1./2.*(fi(m,1,i,j,k) &
                           - fi(m,3,i,j,k) + fi(m,7,i,j,k) &
                           - fi(m,8,i,j,k) - fi(m,9,i,j,k) &
                           + fi(m,10,i,j,k)) &
                           + 1./6.*(fi(m,5,i,j,k) - ftmp(6) &
                           + fi(m,15,i,j,k) + fi(m,16,i,j,k) &
                           - ftmp(17) - ftmp(18)) &
                           - 1./3.*(fi(m,11,i,j,k) - ftmp(13)) &
                           + 2./3.*fi(m,12,i,j,k)

                      fi(m,17,i,j,k) = 1./3.*ftmp(17) &
                           - 1./2.*rhotmp*zp_vals(2,i,j) &
                           - 1./6.*rhotmp*zp_vals(3,i,j) &
                           + 1./2.*(fi(m,2,i,j,k) &
                           - fi(m,4,i,j,k) + fi(m,7,i,j,k) &
                           + fi(m,8,i,j,k) - fi(m,9,i,j,k) &
                           - fi(m,10,i,j,k)) &
                           + 1./6.*(fi(m,5,i,j,k) - ftmp(6) &
                           + fi(m,11,i,j,k) + fi(m,12,i,j,k) &
                           - ftmp(13) - ftmp(14)) &
                           - 1./3.*(fi(m,16,i,j,k) - ftmp(18)) &
                           + 2./3.*fi(m,15,i,j,k) 

                      fi(m,18,i,j,k) = 1./3.*ftmp(18) &
                           + 1./2.*rhotmp*zp_vals(2,i,j) &
                           - 1./6.*rhotmp*zp_vals(3,i,j) &
                           - 1./2.*(fi(m,2,i,j,k) &
                           - fi(m,4,i,j,k) + fi(m,7,i,j,k) &
                           + fi(m,8,i,j,k) - fi(m,9,i,j,k) &
                           - fi(m,10,i,j,k)) &
                           + 1./6.*(fi(m,5,i,j,k) - ftmp(6) &
                           + fi(m,11,i,j,k) + fi(m,12,i,j,k) &
                           - ftmp(13) - ftmp(14)) &
                           - 1./3.*(fi(m,15,i,j,k) - ftmp(17)) &
                           + 2./3.*fi(m,16,i,j,k)
                   endif
                enddo
             enddo
          enddo
       endif
    endif

    ! --- z minus side
    if (bc_flags(5).eq.2) then
       if(info%zs.eq.1) then
          do m=1,info%s
             do i=info%xs,info%xe
                do j=info%ys,info%ye
                   k = 1
                   if (walls(i,j,k).eq.0) then
                      rhotmp = fi(m,0,i,j,k) + fi(m,1,i,j,k) &
                           + fi(m,2,i,j,k) + fi(m,3,i,j,k) &
                           + fi(m,4,i,j,k) + fi(m,7,i,j,k) &
                           + fi(m,8,i,j,k) + fi(m,9,i,j,k) &
                           + fi(m,10,i,j,k) + 2.*(fi(m,6,i,j,k) &
                           + fi(m,13,i,j,k) + fi(m,14,i,j,k) &
                           + fi(m,17,i,j,k) + fi(m,18,i,j,k)) 
                      rhotmp = rhotmp/(1. - zm_vals(3,i,j))
                      
                      ! Choice should not affect the momentum significantly
                      ftmp(5) = fi(m,6,i,j,k);
                      ftmp(11) = fi(m,13,i,j,k)
                      ftmp(12) = fi(m,14,i,j,k)
                      ftmp(15) = fi(m,17,i,j,k)
                      ftmp(16) = fi(m,18,i,j,k)

                      ! This choice does not work for pressure bc
                      !ftmp(5) = fi(m,5,i,j,2)
                      !ftmp(11) = fi(m,11,i+1,j,2)
                      !ftmp(12) = fi(m,12,i-1,j,2)
                      !ftmp(15) = fi(m,15,i,j+1,2)
                      !ftmp(16) = fi(m,16,i,j-1,2)

                      fi(m,5,i,j,k) = 2./3.*ftmp(5) &
                           + 1./3.*rhotmp*zm_vals(3,i,j) &
                           - 1./3.*(ftmp(11)+ ftmp(12) &
                           + ftmp(15) + ftmp(16)) &
                           + 1./3.*(fi(m,6,i,j,k) &
                           + fi(m,13,i,j,k) + fi(m,14,i,j,k) &
                           + fi(m,17,i,j,k) + fi(m,18,i,j,k))
                      
                      fi(m,11,i,j,k) = 1./3.*ftmp(11) &
                           + 1./2.*rhotmp*zm_vals(1,i,j) &
                           + 1./6.*rhotmp*zm_vals(3,i,j) &
                           - 1./2.*(fi(m,1,i,j,k) &
                           - fi(m,3,i,j,k) + fi(m,7,i,j,k) &
                           - fi(m,8,i,j,k) - fi(m,9,i,j,k) &
                           + fi(m,10,i,j,k)) &
                           - 1./6.*(ftmp(5) - fi(m,6,i,j,k) &
                           + ftmp(15) + ftmp(16) &
                           - fi(m,17,i,j,k) - fi(m,18,i,j,k)) &
                           + 1./3.*(ftmp(12) - fi(m,14,i,j,k)) &
                           + 2./3.*fi(m,13,i,j,k) 

                      fi(m,12,i,j,k) = 1./3.*ftmp(12) &
                           - 1./2.*rhotmp*zm_vals(1,i,j) &
                           + 1./6.*rhotmp*zm_vals(3,i,j) &
                           + 1./2.*(fi(m,1,i,j,k) &
                           - fi(m,3,i,j,k) + fi(m,7,i,j,k) &
                           - fi(m,8,i,j,k) - fi(m,9,i,j,k) &
                           + fi(m,10,i,j,k)) &
                           - 1./6.*(ftmp(5) - fi(m,6,i,j,k) &
                           + ftmp(15) + ftmp(16) &
                           - fi(m,17,i,j,k) - fi(m,18,i,j,k)) &
                           + 1./3.*(ftmp(11) - fi(m,13,i,j,k)) &
                           + 2./3.*fi(m,14,i,j,k)

                      fi(m,15,i,j,k) = 1./3.*ftmp(15) &
                           + 1./2.*rhotmp*zm_vals(2,i,j) &
                           + 1./6.*rhotmp*zm_vals(3,i,j) &
                           - 1./2.*(fi(m,2,i,j,k) &
                           - fi(m,4,i,j,k) + fi(m,7,i,j,k) &
                           + fi(m,8,i,j,k) - fi(m,9,i,j,k) &
                           - fi(m,10,i,j,k)) &
                           - 1./6.*(ftmp(5) - fi(m,6,i,j,k) &
                           + ftmp(11) + ftmp(12) &
                           - fi(m,13,i,j,k) - fi(m,14,i,j,k)) &
                           + 1./3.*(ftmp(16) - fi(m,18,i,j,k)) &
                           + 2./3.*fi(m,17,i,j,k)
                      
                      fi(m,16,i,j,k) = 1./3.*ftmp(16) &
                           - 1./2.*rhotmp*zm_vals(2,i,j) &
                           + 1./6.*rhotmp*zm_vals(3,i,j) &
                           + 1./2.*(fi(m,2,i,j,k) &
                           - fi(m,4,i,j,k) + fi(m,7,i,j,k) &
                           + fi(m,8,i,j,k) - fi(m,9,i,j,k) &
                           - fi(m,10,i,j,k)) &
                           - 1./6.*(ftmp(5) - fi(m,6,i,j,k) &
                           + ftmp(11) + ftmp(12) &
                           - fi(m,13,i,j,k) - fi(m,14,i,j,k)) &
                           + 1./3.*(ftmp(15) - fi(m,17,i,j,k)) &
                           + 2./3.*fi(m,18,i,j,k)
                   endif
                enddo
             enddo
          enddo
       endif
    endif


    ! --- other sides not implemented
    if (bc_flags(4).eq.2) then
       write(*,*) 'FLUX BCS IN THE NON-Z-DIRECTIONS ARE NOT IMPLEMENTED -- you get periodic now!'
    endif
    if (bc_flags(3).eq.2) then
       write(*,*) 'FLUX BCS IN THE NON-Z-DIRECTIONS ARE NOT IMPLEMENTED -- you get periodic now!'
    endif
    if (bc_flags(2).eq.2) then
       write(*,*) 'FLUX BCS IN THE NON-Z-DIRECTIONS ARE NOT IMPLEMENTED -- you get periodic now!'
    endif
    if (bc_flags(1).eq.2) then
       write(*,*) 'FLUX BCS IN THE NON-Z-DIRECTIONS ARE NOT IMPLEMENTED -- you get periodic now!'
    endif

    return
  end subroutine BCFlux

  subroutine BCPressure(fi, walls, bc_flags, bc_dim, &
       xm_vals, xp_vals, ym_vals, yp_vals, zm_vals, zp_vals, info)

    type(info_type) info
    PetscInt bc_dim

    PetscInt,dimension(1:6)::bc_flags
    PetscScalar,dimension(1:info%s, 0:info%b, info%gxs:info%gxe, &
         info%gys:info%gye, info%gzs:info%gze):: fi
    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: walls
    PetscScalar,dimension(bc_dim,info%ys:info%ye,info%zs:info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(bc_dim,info%xs:info%xe,info%zs:info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(bc_dim,info%xs:info%xe,info%ys:info%ye):: zm_vals, zp_vals

    PetscInt i,j,k
    PetscInt directions(0:info%b)

    directions(:) = 0
    ! XM BOUNDARY
    if ((bc_flags(BOUNDARY_XM).eq.BC_PRESSURE).and.(info%xs.eq.1)) then
       call BCSetLocalDirections(BOUNDARY_XM, directions)
       do k=info%zs,info%ze
          do j=info%ys,info%ye
             i = 1
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyBoundary(fi(:,:,i,j,k), bc_dim, xm_vals(:,j,k), directions, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! XP BOUNDARY
    if ((bc_flags(BOUNDARY_XP).eq.BC_PRESSURE).and.(info%xe.eq.info%NX)) then
       call BCSetLocalDirections(BOUNDARY_XP, directions)
       do k=info%zs,info%ze
          do j=info%ys,info%ye
             i = info%NX
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyBoundary(fi(:,:,i,j,k), bc_dim, xp_vals(:,j,k), directions, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! YM BOUNDARY
    if ((bc_flags(BOUNDARY_YM).eq.BC_PRESSURE).and.(info%ys.eq.1)) then
       call BCSetLocalDirections(BOUNDARY_YM, directions)
       do k=info%zs,info%ze
          do i=info%xs,info%xe
             j = 1
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyBoundary(fi(:,:,i,j,k), bc_dim, ym_vals(:,i,k), directions, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! YP BOUNDARY
    if ((bc_flags(BOUNDARY_YP).eq.BC_PRESSURE).and.(info%ye.eq.info%NY)) then
       call BCSetLocalDirections(BOUNDARY_YP, directions)
       do k=info%zs,info%ze
          do i=info%xs,info%xe
             j = info%NY
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyBoundary(fi(:,:,i,j,k), bc_dim, yp_vals(:,i,k), directions, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! ZM BOUNDARY
    if ((bc_flags(BOUNDARY_ZM).eq.BC_PRESSURE).and.(info%zs.eq.1)) then
       call BCSetLocalDirections(BOUNDARY_ZM, directions)
       do j=info%ys,info%ye
          do i=info%xs,info%xe
             k = 1
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyBoundary(fi(:,:,i,j,k), bc_dim, zm_vals(:,i,j), directions, info)
             end if
          end do
       end do
    endif

    directions(:) = 0
    ! ZP BOUNDARY
    if ((bc_flags(BOUNDARY_ZP).eq.BC_PRESSURE).and.(info%ze.eq.info%NZ)) then
       call BCSetLocalDirections(BOUNDARY_ZP, directions)
       do j=info%ys,info%ye
          do i=info%xs,info%xe
             k = info%NZ
             if (walls(i,j,k).eq.0) then
                call BCPressureApplyBoundary(fi(:,:,i,j,k), bc_dim, zp_vals(:,i,j), directions, info)
             end if
          end do
       end do
    endif

    return
  end subroutine BCPressure

  subroutine BCPressureApplyBoundary(fi, bc_dim, pvals, directions, info)
    use LBM_Discretization_D3Q19_module
          ! assume x- and y-velocities equal zero

    type(info_type) info
    PetscInt bc_dim
    PetscInt,intent(in),dimension(0:info%b):: directions
    PetscScalar,intent(inout),dimension(1:info%s, 0:info%b):: fi
    PetscScalar,intent(in),dimension(bc_dim):: pvals

!    PetscScalar utmp, vtmp
    PetscScalar wtmp
    PetscScalar,dimension(0:info%b)::ftmp
    integer m
    
!    utmp = 0
!    vtmp = 0
    wtmp = 0

    do m=1,info%s
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
  end subroutine BCPressureApplyBoundary

  subroutine BCSetLocalDirections(boundary, directions)
    use LBM_Discretization_D3Q19_module
    PetscInt,intent(in):: boundary
    PetscInt,intent(out),dimension(0:discretization_directions) :: directions

    PetscInt i
    
    select case(boundary)
    case (BOUNDARY_ZM)
       ! identity mapping
       do i=0,discretization_directions
          directions(i) = i
       end do
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
    end select
  end subroutine BCSetLocalDirections

  subroutine BCPseudoperiodic(fi, walls, bc_flags, bc_dim, &
       xm_vals, xp_vals, ym_vals, yp_vals, zm_vals, zp_vals, info)

    type(info_type) info
    integer bc_dim

    PetscScalar,dimension(1:info%s, 0:info%b, info%gxs:info%gxe, &
         info%gys:info%gye, info%gzs:info%gze):: fi
    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: walls

    PetscInt,dimension(1:6)::bc_flags

    PetscScalar,dimension(bc_dim,info%ys:info%ye,info%zs:info%ze):: xm_vals, xp_vals
    PetscScalar,dimension(bc_dim,info%xs:info%xe,info%zs:info%ze):: ym_vals, yp_vals
    PetscScalar,dimension(bc_dim,info%xs:info%xe,info%ys:info%ye):: zm_vals, zp_vals

    integer i,j,m
    PetscScalar ftmp

    ! note that pseudo-periodic is a requirement on both plus and minus
    ! boundaries.  Only one is checked, and it is assumed that the other
    ! one is correctly flagged

    ! In the y-direction
    if ((info%ye.eq.info%NY).and.(bc_flags(4).eq.1)) then
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


    if ((info%ys.eq.1).and.(bc_flags(3).eq.1)) then
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
    if ((info%ze.eq.info%NZ).and.(bc_flags(6).eq.1)) then
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

    if ((info%zs.eq.1).and.(bc_flags(5).eq.1)) then
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
  end subroutine BCPseudoperiodic
end module LBM_BC_module
