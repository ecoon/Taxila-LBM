!!!==================================================================
!!! Fortran-file
!!!    author:        Ethan T. Coon
!!!    filename:      c_constants.f
!!!    version:
!!!    created:       05 November 2010
!!!      on:          12:16:11 MDT
!!!    last modified:  05 November 2010
!!!      at:          12:16:11 MDT
!!!    URL:           http://www.ldeo.columbia.edu/~ecoon/
!!!    email:         ecoon _at_ ldeo.columbia.edu
!!!
!!!==================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
  
module LBM_Discretization_Directions_D3Q19_module
  implicit none

  PetscInt, parameter, public :: ORIGIN = 0
  PetscInt, parameter, public :: EAST = 1
  PetscInt, parameter, public :: NORTH = 2
  PetscInt, parameter, public :: WEST = 3
  PetscInt, parameter, public :: SOUTH = 4
  PetscInt, parameter, public :: UP = 5
  PetscInt, parameter, public :: DOWN = 6
  PetscInt, parameter, public :: NORTHEAST = 7
  PetscInt, parameter, public :: NORTHWEST = 8
  PetscInt, parameter, public :: SOUTHWEST = 9
  PetscInt, parameter, public :: SOUTHEAST = 10
  PetscInt, parameter, public :: EASTUP = 11
  PetscInt, parameter, public :: WESTUP = 12
  PetscInt, parameter, public :: WESTDOWN = 13
  PetscInt, parameter, public :: EASTDOWN = 14
  PetscInt, parameter, public :: NORTHUP = 15
  PetscInt, parameter, public :: SOUTHUP = 16
  PetscInt, parameter, public :: SOUTHDOWN = 17
  PetscInt, parameter, public :: NORTHDOWN = 18

  PetscInt, public, parameter :: discretization_dims = 3
  PetscInt, public, parameter :: discretization_directions = 18
end module LBM_Discretization_Directions_D3Q19_module

module LBM_Discretization_D3Q19_module
  use LBM_Discretization_Type_module
  use LBM_Discretization_Directions_D3Q19_module
  implicit none

  private
#include "lbm_definitions.h"

  public:: DiscretizationSetup_D3Q19, &
       DiscretizationSetupRelax_D3Q19, &
       DiscretizationEquilf_D3Q19

contains
  subroutine DiscretizationSetup_D3Q19(disc)
    type(discretization_type) disc
    disc%name = D3Q19_DISCRETIZATION
    disc%ndims = 3
    disc%b = 18
    allocate(disc%ci(0:disc%b,1:disc%ndims))
    allocate(disc%weights(0:disc%b))
    allocate(disc%opposites(0:disc%b))
    allocate(disc%mt_mrt(0:disc%b,0:disc%b))                ! transpose of M
    allocate(disc%mmt_mrt(0:disc%b))                        ! diagonal M dot MT matrix 

    disc%opposites(ORIGIN) = ORIGIN
    disc%opposites(EAST) = WEST
    disc%opposites(WEST) = EAST
    disc%opposites(NORTH) = SOUTH
    disc%opposites(SOUTH) = NORTH
    disc%opposites(UP) = DOWN
    disc%opposites(DOWN) = UP
    disc%opposites(NORTHEAST) = SOUTHWEST
    disc%opposites(SOUTHWEST) = NORTHEAST
    disc%opposites(SOUTHEAST) = NORTHWEST
    disc%opposites(NORTHWEST) = SOUTHEAST
    disc%opposites(NORTHUP) = SOUTHDOWN
    disc%opposites(SOUTHDOWN) = NORTHUP
    disc%opposites(SOUTHUP) = NORTHDOWN
    disc%opposites(NORTHDOWN) = SOUTHUP
    disc%opposites(EASTUP) = WESTDOWN
    disc%opposites(WESTDOWN) = EASTUP
    disc%opposites(WESTUP) = EASTDOWN
    disc%opposites(EASTDOWN) = WESTUP

    disc%ci(:,X_DIRECTION) = (/ 0, 1, 0,-1, 0, 0, 0, 1,-1,-1, &
         1, 1,-1,-1, 1, 0, 0, 0, 0/)
    disc%ci(:,Y_DIRECTION) = (/ 0, 0, 1, 0,-1, 0, 0, 1, 1,-1, &
        -1, 0, 0, 0, 0, 1,-1,-1, 1/)
    disc%ci(:,Z_DIRECTION) = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, &
         0, 1, 1,-1,-1, 1, 1,-1,-1/)

    disc%weights = (/ 1.d0/3.d0, &
       1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, 1.d0/18.d0, &
       1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, &
       1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0/)
    
    disc%mt_mrt(:, 0) = (/  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /)
    disc%mt_mrt(:, 1) = (/-30,-11,-11,-11,-11,-11,-11,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8 /)
    disc%mt_mrt(:, 2) = (/ 12, -4, -4, -4, -4, -4, -4,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /)
    disc%mt_mrt(:, 3) = disc%ci(:,X_DIRECTION)
    disc%mt_mrt(:, 4) = (/  0, -4,  0,  4,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,  0,  0,  0,  0 /)
    disc%mt_mrt(:, 5) = disc%ci(:,Y_DIRECTION)
    disc%mt_mrt(:, 6) = (/  0,  0, -4,  0,  4,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1, -1,  1 /)
    disc%mt_mrt(:, 7) = disc%ci(:,Z_DIRECTION)
    disc%mt_mrt(:, 8) = (/  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1 /)
    disc%mt_mrt(:, 9) = (/  0,  2, -1,  2, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1, -2, -2, -2, -2 /)
    disc%mt_mrt(:,10) = (/  0, -4,  2, -4,  2,  2,  2,  1,  1,  1,  1,  1,  1,  1,  1, -2, -2, -2, -2 /)
    disc%mt_mrt(:,11) = (/  0,  0,  1,  0,  1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0 /)
    disc%mt_mrt(:,12) = (/  0,  0, -2,  0, -2,  2,  2,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0 /)
    disc%mt_mrt(:,13) = (/  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0 /)
    disc%mt_mrt(:,14) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1 /)
    disc%mt_mrt(:,15) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0 /)
    disc%mt_mrt(:,16) = (/  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1, -1,  1,  1, -1,  0,  0,  0,  0 /)
    disc%mt_mrt(:,17) = (/  0,  0,  0,  0,  0,  0,  0, -1, -1,  1,  1,  0,  0,  0,  0,  1, -1, -1,  1 /)
    disc%mt_mrt(:,18) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1, -1, -1, -1,  1,  1 /)

    disc%mmt_mrt = (/ 19, 2394, 252, 10, 40, 10, 40, 10, 40, 36, 72, 12, 24, 4, 4, 4, 8, 8, 8 /) 
  end subroutine DiscretizationSetup_D3Q19

  subroutine DiscretizationSetupRelax_D3Q19(disc, relax)
    use LBM_Relaxation_module
    type(discretization_type) disc
    type(relaxation_type) relax
    PetscScalar oneontau
    
    oneontau = 1.d0/relax%tau

    if (relax%mode .eq. RELAXATION_MODE_MRT) then
       !! Curently following the code Qinjun gave me.
       relax%tau_mrt(0) = oneontau
       relax%tau_mrt(1) = oneontau
       relax%tau_mrt(2) = oneontau
       relax%tau_mrt(3) = oneontau  ! 0.d0
       relax%tau_mrt(4) = 8.d0*(2.d0-oneontau)/(8.d0 - oneontau)
       relax%tau_mrt(5) = oneontau  ! 0.d0
       relax%tau_mrt(6) = 8.d0*(2.d0-oneontau)/(8.d0 - oneontau)
       relax%tau_mrt(7) = oneontau  !0.d0
       relax%tau_mrt(8) = 8.d0*(2.d0-oneontau)/(8.d0 - oneontau)
       relax%tau_mrt(9) = oneontau
       relax%tau_mrt(10) = oneontau
       relax%tau_mrt(11) = oneontau
       relax%tau_mrt(12) = oneontau
       relax%tau_mrt(13) = oneontau
       relax%tau_mrt(14) = oneontau
       relax%tau_mrt(15) = oneontau
       relax%tau_mrt(16) = 8.d0*(2.d0-oneontau)/(8.d0 - oneontau)
       relax%tau_mrt(17) = 8.d0*(2.d0-oneontau)/(8.d0 - oneontau)
       relax%tau_mrt(18) = 8.d0*(2.d0-oneontau)/(8.d0 - oneontau)
    end if
  end subroutine DiscretizationSetupRelax_D3Q19

  subroutine DiscretizationEquilf_D3Q19(disc, rho, u, walls, feq, relax, dist)
    use LBM_Distribution_Function_type_module
    use LBM_Relaxation_module
    type(discretization_type) disc
    type(distribution_type) dist
    type(relaxation_type) relax

    PetscScalar,dimension(0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: feq
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: rho
    PetscScalar,dimension(1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: u
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls

    PetscInt i,j,k,d,n,m

    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: usqr
    PetscScalar udote

    usqr = sum(u*u,1)

    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j,k).eq.0) then
       feq(0,i,j,k) = rho(i,j,k)*(relax%d_k - usqr(i,j,k)/2.d0)
       do n=1,dist%b
          udote = sum(disc%ci(n,:)*u(:,i,j,k), 1)
          feq(n,i,j,k)= disc%weights(n)*rho(i,j,k)* &
               (1.5d0*(1.d0-relax%d_k) + udote/relax%c_s2 &
               + udote*udote/(2.d0*relax%c_s2*relax%c_s2) &
               - usqr(i,j,k)/(2.d0*relax%c_s2))
       end do
    end if
    end do
    end do
    end do
  end subroutine DiscretizationEquilf_D3Q19
end module LBM_Discretization_D3Q19_module

