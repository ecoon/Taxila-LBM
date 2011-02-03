!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        bc.F90
!!!     version:         
!!!     created:         06 December 2010
!!!       on:            09:03:18 MST
!!!     last modified:   03 February 2011
!!!       at:            10:36:38 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ ldeo.columbia.edu
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"

module LBM_BC_module
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
  implicit none
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
  implicit none
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
  use LBM_Info_module
  implicit none
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
  use LBM_Info_module
  implicit none
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
  implicit none
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
  implicit none
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
    use LBM_Info_module
    implicit none
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
    use LBM_Info_module
    implicit none
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

    PetscScalar utmp, vtmp, wtmp
    PetscScalar,dimension(0:info%b)::ftmp
    integer i,j,k,m

    ftmp = 0.0

    ! In the z-direction
    ! TOP BOUNDARY
    if (bc_flags(6).eq.3) then
       if(info%ze.eq.info%NZ) then
          ! assume x- and y-velocities equal zero
          utmp = 0.0
          vtmp = 0.0
          do m=1,info%s
             do i=info%xs,info%xe
                do j=info%ys,info%ye
                   k = info%NZ
                   if (walls(i,j,k).eq.0) then
                      wtmp = fi(m,0,i,j,k) + fi(m,1,i,j,k) &
                           + fi(m,2,i,j,k) + fi(m,3,i,j,k) &
                           + fi(m,4,i,j,k) + fi(m,7,i,j,k) &
                           + fi(m,8,i,j,k) + fi(m,9,i,j,k) &
                           + fi(m,10,i,j,k) + 2.*(fi(m,5,i,j,k) &
                           + fi(m,11,i,j,k) + fi(m,12,i,j,k) &
                           + fi(m,15,i,j,k) + fi(m,16,i,j,k))
                      wtmp = wtmp/zp_vals(m,i,j)-1.0

                      ! Choice should not affect the momentum significantly

                      ftmp(6) = fi(m,5,i,j,k)
                      ftmp(13) = fi(m,11,i,j,k)
                      ftmp(14) = fi(m,12,i,j,k)
                      ftmp(17) = fi(m,15,i,j,k)
                      ftmp(18) = fi(m,16,i,j,k)

                      ! This choice does not work for Pressure bc
                      !ftmp(6) = fi(m,6,i,j,k-1)
                      !ftmp(13) = fi(m,13,i-1,j,k-1)
                      !ftmp(14) = fi(m,14,i+1,j,k-1)
                      !ftmp(17) = fi(m,17,i,j-1,k-1)
                      !ftmp(18) = fi(m,18,i,j+1,k-1)

                      fi(m,6,i,j,k) = 2./3.*ftmp(6) &
                           - 1./3.*zp_vals(m,i,j)*wtmp &
                           + 1./3.*(fi(m,5,i,j,k) &
                           + fi(m,11,i,j,k) + fi(m,12,i,j,k) &
                           + fi(m,15,i,j,k) + fi(m,16,i,j,k)) &
                           - 1./3.*(ftmp(13) + ftmp(14) &
                           + ftmp(17) + ftmp(18))
                           
                      fi(m,13,i,j,k) = 1./3.*ftmp(13) &
                           - 1./2.*zp_vals(m,i,j)*utmp &
                           - 1./6.*zp_vals(m,i,j)*wtmp &
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
                           + 1./2.*zp_vals(m,i,j)*utmp &
                           - 1./6.*zp_vals(m,i,j)*wtmp &
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
                           - 1./2.*zp_vals(m,i,j)*vtmp &
                           - 1./6.*zp_vals(m,i,j)*wtmp &
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
                           + 1./2.*zp_vals(m,i,j)*vtmp &
                           - 1./6.*zp_vals(m,i,j)*wtmp &
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


    ! BOTTOM BOUNDARY
    if (bc_flags(5).eq.3) then
       if(info%zs.eq.1) then
          ! assume x- and y-velocities equal zero
          utmp = 0.0
          vtmp = 0.0
          do m=1,info%s
             do i=info%xs,info%xe
                do j=info%ys,info%ye
                   k = 1
                   if (walls(i,j,k).eq.0) then
                      wtmp = fi(m,0,i,j,k) + fi(m,1,i,j,k) &
                           + fi(m,2,i,j,k) + fi(m,3,i,j,k) &
                           + fi(m,4,i,j,k) + fi(m,7,i,j,k) &
                           + fi(m,8,i,j,k) + fi(m,9,i,j,k) &
                           + fi(m,10,i,j,k) + 2.*(fi(m,6,i,j,k) &
                           + fi(m,13,i,j,k) + fi(m,14,i,j,k) &
                           + fi(m,17,i,j,k) + fi(m,18,i,j,k))
                      wtmp = 1.0-wtmp/zm_vals(m,i,j)


                      ! Choice should not affect the momentum significantly
                      ftmp(5) = fi(m,6,i,j,k)
                      ftmp(11) = fi(m,13,i,j,k)
                      ftmp(12) = fi(m,14,i,j,k)
                      ftmp(15) = fi(m,17,i,j,k)
                      ftmp(16) = fi(m,18,i,j,k)

                      ! This choice does not work for Pressure bc
                      !ftmp(5) = fi(m,5,i,j,2)
                      !ftmp(11) = fi(m,11,i+1,j,2)
                      !ftmp(12) = fi(m,12,i-1,j,2)
                      !ftmp(15) = fi(m,15,i,j+1,2)
                      !ftmp(16) = fi(m,16,i,j-1,2)

                      fi(m,5,i,j,k) = 2./3.*ftmp(5) &
                           + 1./3.*zm_vals(m,i,j)*wtmp &
                           - 1./3.*(ftmp(11)+ ftmp(12) &
                           + ftmp(15) + ftmp(16)) &
                           + 1./3.*(fi(m,6,i,j,k) &
                           + fi(m,13,i,j,k) + fi(m,14,i,j,k) &
                           + fi(m,17,i,j,k) + fi(m,18,i,j,k))
                           
                      fi(m,11,i,j,k) = 1./3.*ftmp(11) &
                           + 1./2.*zm_vals(m,i,j)*utmp &
                           + 1./6.*zm_vals(m,i,j)*wtmp &
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
                           - 1./2.*zm_vals(m,i,j)*utmp &
                           + 1./6.*zm_vals(m,i,j)*wtmp &
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
                           + 1./2.*zm_vals(m,i,j)*vtmp &
                           + 1./6.*zm_vals(m,i,j)*wtmp &
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
                           - 1./2.*zm_vals(m,i,j)*vtmp &
                           + 1./6.*zm_vals(m,i,j)*wtmp &
                           + 1./2.*(fi(m,2,i,j,k) &
                           - fi(m,4,i,j,k) + fi(m,7,i,j,k) &
                           + fi(m,8,i,j,k) - fi(m,9,i,j,k) &
                           - fi(m,10,i,j,k)) &
                           - 1./6.*(ftmp(5) - fi(m,6,i,j,k) &
                           + ftmp(11) + ftmp(12) &
                           - fi(m,13,i,j,k) - fi(m,14,i,j,k)) &
                           + 1./3.*(ftmp(15) - fi(m,17,i,j,k)) &
                           + 2./3.*fi(m,18,i,j,k)
                   end if
                enddo
             enddo
          enddo
       endif
    endif

    return
  end subroutine BCPressure

!#include "bc_pressure_dev.F90"

  subroutine BCPseudoperiodic(fi, walls, bc_flags, bc_dim, &
       xm_vals, xp_vals, ym_vals, yp_vals, zm_vals, zp_vals, info)
    use LBM_Info_module
    implicit none

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
