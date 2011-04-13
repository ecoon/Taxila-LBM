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

module LBM_Discretization_Directions_D2Q9_module
  implicit none
  
  ! discretization-specific enum directions
  PetscInt, public, parameter :: ORIGIN = 0
  PetscInt, public, parameter :: EAST = 1
  PetscInt, public, parameter :: NORTH = 2
  PetscInt, public, parameter :: WEST = 3
  PetscInt, public, parameter :: SOUTH = 4
  PetscInt, public, parameter :: NORTHEAST = 5
  PetscInt, public, parameter :: NORTHWEST = 6
  PetscInt, public, parameter :: SOUTHWEST = 7
  PetscInt, public, parameter :: SOUTHEAST = 8

  PetscInt, public, parameter :: discretization_dims = 2
  PetscInt, public, parameter :: discretization_directions = 8
end module LBM_Discretization_Directions_D2Q9_module
  
module LBM_Discretization_D2Q9_module
  use LBM_Discretization_Type_module
  use LBM_Discretization_Directions_D2Q9_module
  implicit none

  private
#include "lbm_definitions.h"

  public:: DiscretizationSetUp_D2Q9, &
       DiscretizationSetUpRelax_D2Q9, &
       DiscretizationEquilf_D2Q9

contains
  subroutine DiscretizationSetUp_D2Q9(disc)
    type(discretization_type) disc
    disc%name = D2Q9_DISCRETIZATION
    disc%ndims = discretization_dims
    disc%b = discretization_directions
    allocate(disc%ci(0:disc%b,1:disc%ndims))
    allocate(disc%weights(0:disc%b))
    allocate(disc%opposites(0:disc%b))
    allocate(disc%mt_mrt(0:disc%b,0:disc%b))         ! transpose of M
    allocate(disc%mmt_mrt(0:disc%b))                 ! diagonal M dot MT matrix
    allocate(disc%ffw(0:2,0:2))  

    disc%opposites(ORIGIN) = ORIGIN
    disc%opposites(EAST) = WEST
    disc%opposites(WEST) = EAST
    disc%opposites(NORTH) = SOUTH
    disc%opposites(SOUTH) = NORTH
    disc%opposites(NORTHEAST) = SOUTHWEST
    disc%opposites(SOUTHWEST) = NORTHEAST
    disc%opposites(SOUTHEAST) = NORTHWEST
    disc%opposites(NORTHWEST) = SOUTHEAST

    disc%ci(:,X_DIRECTION) = (/ 0, 1, 0,-1, 0, 1,-1,-1, 1/)
    disc%ci(:,Y_DIRECTION) = (/ 0, 0, 1, 0,-1, 1, 1,-1,-1/)

    disc%weights = (/ 4.d0/9.d0, &
       1.d0/9.d0,  1.d0/9.d0,  1.d0/9.d0,  1.d0/9.d0, &
       1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0, 1.d0/36.d0 /)

    disc%mt_mrt(:, 0) = (/  1,  1,  1,  1,  1,  1,  1,  1,  1 /)
    disc%mt_mrt(:, 1) = (/ -4, -1, -1, -1, -1,  2,  2,  2,  2 /)
    disc%mt_mrt(:, 2) = (/  4, -2, -2, -2, -2,  1,  1,  1,  1 /)
    disc%mt_mrt(:, 3) = disc%ci(:,X_DIRECTION)
    disc%mt_mrt(:, 4) = (/  0, -2,  0,  2,  0,  1, -1, -1,  1 /)
    disc%mt_mrt(:, 5) = disc%ci(:,Y_DIRECTION)
    disc%mt_mrt(:, 6) = (/  0,  0, -2,  0,  2,  1,  1, -1, -1 /)
    disc%mt_mrt(:, 7) = (/  0,  1, -1,  1, -1,  0,  0,  0,  0 /)
    disc%mt_mrt(:, 8) = (/  0,  0,  0,  0,  0,  1, -1,  1, -1 /)

    disc%mmt_mrt = (/ 9, 36, 36, 6, 12, 6, 12, 4, 4 /)

!!$    disc%ffw = 0.0
!!$    disc%ffw( 1, 0) = 4./21.
!!$    disc%ffw( 0, 1) = 4./21.
!!$    disc%ffw(-1, 0) = 4./21.
!!$    disc%ffw( 0,-1) = 4./21.
!!$    disc%ffw( 1, 1) = 4./45.
!!$    disc%ffw(-1, 1) = 4./45.
!!$    disc%ffw(-1,-1) = 4./45.
!!$    disc%ffw( 1,-1) = 4./45.
!!$    disc%ffw( 2, 0) = 1./60.
!!$    disc%ffw( 0, 2) = 1./60.
!!$    disc%ffw(-2, 0) = 1./60.
!!$    disc%ffw( 0,-2) = 1./60.
!!$    disc%ffw( 2, 1) = 2./315.
!!$    disc%ffw( 1, 2) = 2./315.
!!$    disc%ffw(-1, 2) = 2./315.
!!$    disc%ffw(-2, 1) = 2./315.
!!$    disc%ffw(-2,-1) = 2./315.
!!$    disc%ffw(-1,-2) = 2./315.
!!$    disc%ffw( 1,-2) = 2./315.
!!$    disc%ffw( 2,-1) = 2./315.
!!$    disc%ffw( 2, 2) = 1./5040.
!!$    disc%ffw(-2, 2) = 1./5040.
!!$    disc%ffw(-2,-2) = 1./5040.
!!$    disc%ffw( 2,-2) = 1./5040.

  end subroutine DiscretizationSetUp_D2Q9


  subroutine DiscretizationSetUpRelax_D2Q9(disc, relax)
    use LBM_Relaxation_module
    type(discretization_type) disc
    type(relaxation_type) relax
    PetscScalar oneontau
     
    oneontau = 1.d0/relax%tau
    if (relax%mode .eq. RELAXATION_MODE_MRT) then
       !! Curently following the code Qinjun gave me.
       relax%tau_mrt(0) = 1.d0
       relax%tau_mrt(1) = 0.1d0
       relax%tau_mrt(2) = 1.5d0
       relax%tau_mrt(3) = 1.d0
       relax%tau_mrt(4) = 0.7d0
       relax%tau_mrt(5) = 1.d0
       relax%tau_mrt(6) = 0.7d0
       relax%tau_mrt(7) = oneontau
       relax%tau_mrt(8) = oneontau
    end if

  end subroutine DiscretizationSetUpRelax_D2Q9

  subroutine DiscretizationEquilf_D2Q9(disc, rho, u, walls, feq, relax, dist)
    use LBM_Distribution_Function_type_module
    use LBM_Relaxation_module
    type(discretization_type) disc
    type(distribution_type) dist
    type(relaxation_type) relax

    PetscScalar,dimension(0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: feq
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: rho
    PetscScalar,dimension(1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls

    PetscInt i,j,d,n
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: usqr
    PetscScalar udote

    usqr = sum(u*u,1)

    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j).eq.0) then
       feq(0,i,j) = rho(i,j)*((1.d0 + relax%d_k*5.d0)/6.d0 - 2.d0*usqr(i,j)/3.d0)
       do n=1,dist%b
          udote = sum(disc%ci(n,:)*u(:,i,j), 1)
          feq(n,i,j)= disc%weights(n)*rho(i,j)* &
               (1.5d0*(1.d0-relax%d_k) + udote/relax%c_s2 &
               + udote*udote/(2.d0*relax%c_s2*relax%c_s2) &
               - usqr(i,j)/(2.d0*relax%c_s2))
       end do
    end if
    end do
    end do
  end subroutine DiscretizationEquilf_D2Q9
end module LBM_Discretization_D2Q9_module

