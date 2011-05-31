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
  use LBM_Distribution_Function_type_module
  implicit none

  private
#include "lbm_definitions.h"

  public:: DiscretizationSetup_D3Q19, &
       DiscretizationSetupRelax_D3Q19, &
       DiscretizationEquilf_D3Q19, &
       DiscApplyBCDirichletToBoundary_D3Q19, &
       DiscApplyBCVelocityToBoundary_D3Q19, &
       DiscApplyBCFluxToBoundary_D3Q19, &
       DiscSetLocalDirections_D3Q19

contains
  subroutine DiscretizationSetup_D3Q19(disc)
    type(discretization_type) disc
    disc%name = D3Q19_DISCRETIZATION
    disc%ndims = 3
    disc%b = 18
    allocate(disc%ci(0:disc%b,1:disc%ndims))
    allocate(disc%weights(0:disc%b))
    allocate(disc%opposites(0:disc%b))
    allocate(disc%mt_mrt(0:disc%b,0:disc%b))         ! transpose of M
    allocate(disc%mmt_mrt(0:disc%b))                 ! diagonal M dot MT matrix 
    allocate(disc%ffw(1:4*disc%stencil_size))        ! slightly larger than needed in all cases

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

    disc%ffw = 0.0
    if(disc%stencil_size.eq.1) then
      disc%ffw(1) = 1./6.
      disc%ffw(2) = 1./12.
    end if

    if(disc%stencil_size.eq.2) then 
      disc%ffw(1) = 4./45.
      disc%ffw(2) = 1./21.
      disc%ffw(3) = 2./105.
      disc%ffw(4) = 5./504.
      disc%ffw(5) = 1./315.
      disc%ffw(6) = 1./630.
      disc%ffw(8) = 1./5040.
    end if

    ! Not implemented yet
    if(disc%stencil_size.eq.3) then 
      disc%ffw( 1) = 352./5355.
      disc%ffw( 2) = 38./1071.
      disc%ffw( 3) = 271./14280.
      disc%ffw( 4) = 139./14280.
      disc%ffw( 5) = 53./10710.
      disc%ffw( 6) = 5./2142.
      disc%ffw( 7) = 1./4284.   ! w_{221}(9) in Sbragaglia, Phys Rev E (2007)
      disc%ffw( 8) = 41./85680.
      disc%ffw( 9) = 1./5355.   ! w_{300}(9) in Sbragaglia, Phys Rev E (2007)
      disc%ffw(10) = 1./10710.
      disc%ffw(11) = 1./42840.    
    end if

  end subroutine DiscretizationSetup_D3Q19

  subroutine DiscretizationSetupRelax_D3Q19(disc, relax)
    use LBM_Relaxation_module
    type(discretization_type) disc
    type(relaxation_type) relax
    PetscScalar oneontau
    
    oneontau = 1.d0/relax%tau

    if (relax%mode .eq. RELAXATION_MODE_MRT) then
       !! Not tested yet
       relax%tau_mrt(0) = oneontau
       relax%tau_mrt(1) = 0.6 !oneontau  !s_e
       relax%tau_mrt(2) = 1.4 !oneontau  !s_ep
       relax%tau_mrt(3) = oneontau 
       relax%tau_mrt(4) = 1.99 !8.d0*(2.d0-oneontau)/(8.d0-oneontau) !s_q
       relax%tau_mrt(5) = oneontau 
       relax%tau_mrt(6) = 1.99 !8.d0*(2.d0-oneontau)/(8.d0-oneontau) !s_q
       relax%tau_mrt(7) = oneontau 
       relax%tau_mrt(8) = 1.99  !8.d0*(2.d0-oneontau)/(8.d0-oneontau) !s_q
       relax%tau_mrt(9) = oneontau  !s_v
       relax%tau_mrt(10) = 0.1 !oneontau !s_pi
       relax%tau_mrt(11) = oneontau !s_v
       relax%tau_mrt(12) = 0.1 !oneontau !s_pi
       relax%tau_mrt(13) = oneontau !s_v
       relax%tau_mrt(14) = oneontau !s_v
       relax%tau_mrt(15) = oneontau !s_v
       relax%tau_mrt(16) = 0.3 !8.d0*(2.d0-oneontau)/(8.d0-oneontau) !s_m
       relax%tau_mrt(17) = 0.3 !8.d0*(2.d0-oneontau)/(8.d0-oneontau) !s_m
       relax%tau_mrt(18) = 0.3 !8.d0*(2.d0-oneontau)/(8.d0-oneontau) !s_m
    end if
  end subroutine DiscretizationSetupRelax_D3Q19

  subroutine DiscretizationEquilf_D3Q19(disc, rho, u, walls, feq, relax, dist)
    use LBM_Relaxation_module
    type(discretization_type) disc
    type(distribution_type) dist
    type(relaxation_type) relax

    PetscScalar,dimension(0:disc%b, dist%info%gxs:dist%info%gxe, &
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
       do n=1,disc%b
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

  subroutine DiscApplyBCDirichletToBoundary_D3Q19(disc, fi, pvals, directions, dist)
    type(discretization_type) disc
    type(distribution_type) dist
    PetscInt,intent(in),dimension(0:disc%b):: directions
    PetscScalar,intent(inout),dimension(1:dist%s, 0:disc%b):: fi
    PetscScalar,intent(in),dimension(dist%s,dist%info%ndims):: pvals

    PetscScalar wtmp
    PetscScalar,dimension(0:disc%b)::ftmp
    integer m
    
    wtmp = 0

    do m=1,dist%s
       ftmp = 0.0
       wtmp = fi(m,directions(ORIGIN)) + fi(m,directions(EAST)) &
            + fi(m,directions(NORTH)) + fi(m,directions(WEST)) &
            + fi(m,directions(SOUTH)) + fi(m,directions(NORTHEAST)) &
            + fi(m,directions(NORTHWEST)) + fi(m,directions(SOUTHWEST)) &
            + fi(m,directions(SOUTHEAST)) + 2.*(fi(m,directions(DOWN)) &
            + fi(m,directions(WESTDOWN)) + fi(m,directions(EASTDOWN)) &
            + fi(m,directions(SOUTHDOWN)) + fi(m,directions(NORTHDOWN)))
       wtmp = wtmp - pvals(m,1)
       
       ! Choice should not affect the momentum significantly
       ftmp(directions(UP)) = fi(m,directions(DOWN))
       ftmp(directions(EASTUP)) = fi(m,directions(WESTDOWN))
       ftmp(directions(WESTUP)) = fi(m,directions(EASTDOWN))
       ftmp(directions(NORTHUP)) = fi(m,directions(SOUTHDOWN))
       ftmp(directions(SOUTHUP)) = fi(m,directions(NORTHDOWN))
       
       fi(m,directions(UP)) = 2./3.*ftmp(directions(UP)) &
            + 1./3.*wtmp &
            - 1./3.*(ftmp(directions(EASTUP))+ ftmp(directions(WESTUP)) &
            + ftmp(directions(NORTHUP)) + ftmp(directions(SOUTHUP))) &
            + 1./3.*(fi(m,directions(DOWN)) &
            + fi(m,directions(WESTDOWN)) + fi(m,directions(EASTDOWN)) &
            + fi(m,directions(SOUTHDOWN)) + fi(m,directions(NORTHDOWN)))
       
       fi(m,directions(EASTUP)) = 1./3.*ftmp(directions(EASTUP)) &
            + 1./6.*wtmp &
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
            + 1./6.*wtmp &
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
            + 1./6.*wtmp &
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
            + 1./6.*wtmp &
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
  end subroutine DiscApplyBCDirichletToBoundary_D3Q19

  subroutine DiscApplyBCVelocityToBoundary_D3Q19(disc, fi, fvals, directions, cardinals, &
       dist)
    type(discretization_type) disc
    type(distribution_type) dist
    PetscScalar,intent(inout),dimension(1:dist%s, 0:disc%b):: fi
    PetscScalar,intent(in),dimension(dist%s,dist%info%ndims):: fvals
    PetscInt,intent(in),dimension(0:disc%b):: directions
    PetscInt,intent(in),dimension(1:dist%info%ndims):: cardinals

    PetscScalar rhotmp
    PetscScalar,dimension(0:disc%b)::ftmp
    PetscInt m

    ftmp = 0.0
    rhotmp = 0.0

    do m=1,dist%s
       rhotmp = fi(m,directions(ORIGIN)) + fi(m,directions(EAST)) &
            + fi(m,directions(NORTH)) + fi(m,directions(WEST)) &
            + fi(m,directions(SOUTH)) + fi(m,directions(NORTHEAST)) &
            + fi(m,directions(NORTHWEST)) + fi(m,directions(SOUTHWEST)) &
            + fi(m,directions(SOUTHEAST)) + 2.*(fi(m,directions(DOWN)) &
            + fi(m,directions(WESTDOWN)) + fi(m,directions(EASTDOWN)) &
            + fi(m,directions(SOUTHDOWN)) + fi(m,directions(NORTHDOWN))) 
       rhotmp = rhotmp/(1. - fvals(1,cardinals(CARDINAL_NORMAL)))
       
       ! Choice should not affect the momentum significantly
       ftmp(directions(UP)) = fi(m,directions(DOWN))
       ftmp(directions(EASTUP)) = fi(m,directions(WESTDOWN))
       ftmp(directions(WESTUP)) = fi(m,directions(EASTDOWN))
       ftmp(directions(NORTHUP)) = fi(m,directions(SOUTHDOWN))
       ftmp(directions(SOUTHUP)) = fi(m,directions(NORTHDOWN))

       fi(m,directions(UP)) = 2./3.*ftmp(directions(UP)) &
            + 1./3.*rhotmp*fvals(1,cardinals(CARDINAL_NORMAL)) &
            - 1./3.*(ftmp(directions(EASTUP))+ ftmp(directions(WESTUP)) &
            + ftmp(directions(NORTHUP)) + ftmp(directions(SOUTHUP))) &
            + 1./3.*(fi(m,directions(DOWN)) &
            + fi(m,directions(WESTDOWN)) + fi(m,directions(EASTDOWN)) &
            + fi(m,directions(SOUTHDOWN)) + fi(m,directions(NORTHDOWN)))
       
       fi(m,directions(EASTUP)) = 1./3.*ftmp(directions(EASTUP)) &
            + 1./2.*rhotmp*fvals(1,cardinals(CARDINAL_CROSS)) &
            + 1./6.*rhotmp*fvals(1,cardinals(CARDINAL_NORMAL)) &
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
            - 1./2.*rhotmp*fvals(1,cardinals(CARDINAL_CROSS)) &
            + 1./6.*rhotmp*fvals(1,cardinals(CARDINAL_NORMAL)) &
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
            + 1./2.*rhotmp*fvals(1,cardinals(CARDINAL_RESULTANT)) &
            + 1./6.*rhotmp*fvals(1,cardinals(CARDINAL_NORMAL)) &
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
            - 1./2.*rhotmp*fvals(1,cardinals(CARDINAL_RESULTANT)) &
            + 1./6.*rhotmp*fvals(1,cardinals(CARDINAL_NORMAL)) &
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
  end subroutine DiscApplyBCVelocityToBoundary_D3Q19

  subroutine DiscApplyBCFluxToBoundary_D3Q19(disc, fi, fvals, directions, cardinals, &
       dist)
    type(discretization_type) disc
    type(distribution_type) dist
    PetscScalar,intent(inout),dimension(1:dist%s, 0:disc%b):: fi
    PetscScalar,intent(in),dimension(dist%s,dist%info%ndims):: fvals
    PetscInt,intent(in),dimension(0:disc%b):: directions
    PetscInt,intent(in),dimension(1:dist%info%ndims):: cardinals

    PetscScalar rhotmp
    PetscScalar,dimension(0:disc%b)::ftmp
    PetscInt m

    ftmp = 0.0
    rhotmp = 0.0

    do m=1,dist%s
       ! Choice should not affect the momentum significantly
       ftmp(directions(UP)) = fi(m,directions(DOWN))
       ftmp(directions(EASTUP)) = fi(m,directions(WESTDOWN))
       ftmp(directions(WESTUP)) = fi(m,directions(EASTDOWN))
       ftmp(directions(NORTHUP)) = fi(m,directions(SOUTHDOWN))
       ftmp(directions(SOUTHUP)) = fi(m,directions(NORTHDOWN))

       fi(m,directions(UP)) = 2./3.*ftmp(directions(UP)) &
            + 1./3.*fvals(m,cardinals(CARDINAL_NORMAL)) &
            - 1./3.*(ftmp(directions(EASTUP))+ ftmp(directions(WESTUP)) &
            + ftmp(directions(NORTHUP)) + ftmp(directions(SOUTHUP))) &
            + 1./3.*(fi(m,directions(DOWN)) &
            + fi(m,directions(WESTDOWN)) + fi(m,directions(EASTDOWN)) &
            + fi(m,directions(SOUTHDOWN)) + fi(m,directions(NORTHDOWN)))
       
       fi(m,directions(EASTUP)) = 1./3.*ftmp(directions(EASTUP)) &
            + 1./2.*fvals(m,cardinals(CARDINAL_CROSS)) &
            + 1./6.*fvals(m,cardinals(CARDINAL_NORMAL)) &
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
            - 1./2.*fvals(m,cardinals(CARDINAL_CROSS)) &
            + 1./6.*fvals(m,cardinals(CARDINAL_NORMAL)) &
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
            + 1./2.*fvals(m,cardinals(CARDINAL_RESULTANT)) &
            + 1./6.*fvals(m,cardinals(CARDINAL_NORMAL)) &
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
            - 1./2.*fvals(m,cardinals(CARDINAL_RESULTANT)) &
            + 1./6.*fvals(m,cardinals(CARDINAL_NORMAL)) &
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
  end subroutine DiscApplyBCFluxToBoundary_D3Q19

  subroutine DiscSetLocalDirections_D3Q19(disc, boundary, directions, cardinals)
    type(discretization_type) disc
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
  end subroutine DiscSetLocalDirections_D3Q19
end module LBM_Discretization_D3Q19_module

