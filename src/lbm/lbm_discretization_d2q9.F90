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
  use LBM_Distribution_Function_type_module
  implicit none

  private
#include "lbm_definitions.h"

  public:: DiscretizationSetUp_D2Q9, &
       DiscretizationSetUpRelax_D2Q9, &
       DiscretizationEquilf_D2Q9, &
       DiscApplyBCDirichletToBoundary_D2Q9, &
       DiscApplyBCFluxToBoundary_D2Q9, &
       DiscApplyBCVelocityToBoundary_D2Q9, &
       DiscSetLocalDirections_D2Q9

contains
  subroutine DiscretizationSetUp_D2Q9(disc)
    type(discretization_type) disc
    disc%name = D2Q9_DISCRETIZATION
    disc%ndims = discretization_dims
    disc%b = discretization_directions
    allocate(disc%ci(0:disc%b,1:disc%ndims))
    allocate(disc%weights(0:disc%b))
    allocate(disc%opposites(0:disc%b))              ! direction reflected in all directions
    allocate(disc%reflect_x(0:disc%b))              ! direction reflected in x-direction
    allocate(disc%reflect_y(0:disc%b))              ! direction reflected in y-direction
    allocate(disc%mt_mrt(0:disc%b,0:disc%b))        ! transpose of M
    allocate(disc%mmt_mrt(0:disc%b))                ! diagonal M dot MT matrix
    allocate(disc%ffw(1:4*disc%isotropy_order))        ! slightly larger than needed in all cases

    disc%local_normal = SOUTH

    disc%opposites(ORIGIN) = ORIGIN
    disc%opposites(EAST) = WEST
    disc%opposites(WEST) = EAST
    disc%opposites(NORTH) = SOUTH
    disc%opposites(SOUTH) = NORTH
    disc%opposites(NORTHEAST) = SOUTHWEST
    disc%opposites(SOUTHWEST) = NORTHEAST
    disc%opposites(SOUTHEAST) = NORTHWEST
    disc%opposites(NORTHWEST) = SOUTHEAST

    disc%reflect_x(ORIGIN) = ORIGIN
    disc%reflect_x(EAST) = WEST
    disc%reflect_x(WEST) = EAST
    disc%reflect_x(NORTH) = NORTH
    disc%reflect_x(SOUTH) = SOUTH
    disc%reflect_x(NORTHEAST) = NORTHWEST
    disc%reflect_x(SOUTHWEST) = SOUTHEAST
    disc%reflect_x(SOUTHEAST) = SOUTHWEST
    disc%reflect_x(NORTHWEST) = NORTHEAST

    disc%reflect_y(ORIGIN) = ORIGIN
    disc%reflect_y(EAST) = EAST
    disc%reflect_y(WEST) = WEST
    disc%reflect_y(NORTH) = SOUTH
    disc%reflect_y(SOUTH) = NORTH
    disc%reflect_y(NORTHEAST) = SOUTHEAST
    disc%reflect_y(SOUTHWEST) = NORTHWEST
    disc%reflect_y(SOUTHEAST) = NORTHEAST
    disc%reflect_y(NORTHWEST) = SOUTHWEST

    disc%ci(:,X_DIRECTION) = (/ 0, 1, 0,-1, 0, 1,-1,-1, 1/)
    disc%ci(:,Y_DIRECTION) = (/ 0, 0, 1, 0,-1, 1, 1,-1,-1/)

    disc%c_0 = 6.d0

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

    disc%ffw = 0.0
    if(disc%isotropy_order.eq.4) then
      disc%ffw(1) = 1.d0/3.
      disc%ffw(2) = 1.d0/12.
    end if
    
    if(disc%isotropy_order.eq.8) then 
      disc%ffw(1) = 4.d0/21.
      disc%ffw(2) = 4.d0/45.
      disc%ffw(4) = 1.d0/60.
      disc%ffw(5) = 2.d0/315.
      disc%ffw(8) = 1.d0/5040.
    end if

    if(disc%isotropy_order.eq.10) then 
      disc%ffw( 1) = 262.d0/1785.
      disc%ffw( 2) = 93.d0/1190.
      disc%ffw( 4) = 7.d0/340.
      disc%ffw( 5) = 6.d0/595.
      disc%ffw( 8) = 9.d0/9520.
      disc%ffw( 9) = 2.d0/5355.
      disc%ffw(10) = 1.d0/7140.
    end if

  end subroutine DiscretizationSetUp_D2Q9


  subroutine DiscretizationSetUpRelax_D2Q9(disc, relax)
    use LBM_Relaxation_module
    type(discretization_type) disc
    type(relaxation_type) relax
    PetscScalar oneontau
     
    !oneontau = 1.d0/relax%tau
    if (relax%mode .eq. RELAXATION_MODE_MRT) then
       relax%tau_mrt(0) = relax%s_c 
       relax%tau_mrt(1) = relax%s_e 
       relax%tau_mrt(2) = relax%s_e2 
       relax%tau_mrt(3) = relax%s_c 
       relax%tau_mrt(4) = relax%s_q 
       relax%tau_mrt(5) = relax%s_c
       relax%tau_mrt(6) = relax%s_q 
       relax%tau_mrt(7) = relax%s_nu
       relax%tau_mrt(8) = relax%s_nu
    end if

  end subroutine DiscretizationSetUpRelax_D2Q9

  subroutine DiscretizationEquilf_D2Q9(disc, rho, u, walls, feq, m, relax, dist)
    use LBM_Relaxation_module
    type(discretization_type) disc
    type(distribution_type) dist
    type(relaxation_type) relax

    PetscScalar,dimension(dist%s,0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: feq
    PetscScalar,dimension(dist%s,dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls
    PetscInt m

    PetscInt i,j,d,n
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: usqr
    PetscScalar udote

    usqr = sum(u(m,:,:,:)*u(m,:,:,:),1)

    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j).eq.0) then
       feq(m,0,i,j) = rho(m,i,j)*((1. + relax%d_k*5.)/6. - 2.*usqr(i,j)/3.)
       do n=1,dist%b
          udote = sum(disc%ci(n,:)*u(m,:,i,j), 1)
          feq(m,n,i,j)= disc%weights(n)*rho(m,i,j)* &
               (1.5*(1.-relax%d_k) + udote/relax%c_s2 &
               + udote*udote/(2.*relax%c_s2*relax%c_s2) &
               - usqr(i,j)/(2.*relax%c_s2))
       end do
    end if
    end do
    end do
  end subroutine DiscretizationEquilf_D2Q9

  subroutine DiscApplyBCDirichletToBoundary_D2Q9(disc, fi, forces, pvals, &
       directions, cardinals, dist)
    type(discretization_type) disc
    type(distribution_type) dist
    PetscScalar,intent(inout),dimension(1:dist%s, 0:dist%b):: fi
    PetscScalar,intent(in),dimension(1:dist%s, dist%info%ndims):: forces
    PetscInt,intent(in),dimension(0:dist%b):: directions
    PetscInt,intent(in),dimension(1:dist%info%ndims):: cardinals
    PetscScalar,intent(in),dimension(dist%s,dist%info%ndims):: pvals

    PetscScalar rhovtmp
    PetscScalar,dimension(0:dist%b)::ftmp
    PetscInt m
    rhovtmp = 0

    ! written for the NORTH boundary.
    do m=1,dist%s
       ftmp = 0.0
       rhovtmp = fi(m,directions(ORIGIN)) + fi(m,directions(EAST)) &
            + fi(m,directions(WEST)) + 2.*(fi(m,directions(NORTH)) &
            + fi(m,directions(NORTHEAST)) + fi(m,directions(NORTHWEST))) &
            + forces(m,cardinals(CARDINAL_NORMAL))/2.d0
       rhovtmp = rhovtmp - pvals(m,1)

       ! Choice should not affect the momentum significantly
       ftmp(directions(SOUTH)) = fi(m,directions(NORTH))
       ftmp(directions(SOUTHWEST)) = fi(m,directions(NORTHEAST))
       ftmp(directions(SOUTHEAST)) = fi(m,directions(NORTHWEST))

       fi(m,directions(SOUTH)) = 1./3.*ftmp(directions(SOUTH)) &
            - 2./3.*(rhovtmp - forces(m,cardinals(CARDINAL_NORMAL))/2.d0) &
            - 2./3.*(ftmp(directions(SOUTHWEST)) + ftmp(directions(SOUTHEAST))) &
            + 2./3.*(fi(m,directions(NORTH)) + fi(m,directions(NORTHEAST)) + fi(m,directions(NORTHWEST)))

       fi(m,directions(SOUTHWEST)) = 1./3.*ftmp(directions(SOUTHWEST)) &
            - 1./6.*(rhovtmp - forces(m,cardinals(CARDINAL_NORMAL))/2.d0) &
            - 1./2.*(-forces(m,cardinals(CARDINAL_CROSS))/2.0) &
            + 1./2.*(fi(m,directions(EAST)) - fi(m,directions(WEST))) &
            + 1./6.*(fi(m,directions(NORTH)) - ftmp(directions(SOUTH))) &
            - 1./3.*(fi(m,directions(NORTHWEST)) - ftmp(directions(SOUTHEAST))) &
            + 2./3.*fi(m,directions(NORTHEAST))

       fi(m,directions(SOUTHEAST)) = 1./3.*ftmp(directions(SOUTHEAST)) &
            - 1./6.*(rhovtmp - forces(m,cardinals(CARDINAL_NORMAL))/2.d0) &
            + 1./2.*(-forces(m,cardinals(CARDINAL_CROSS))/2.0) &
            - 1./2.*(fi(m,directions(EAST)) - fi(m,directions(WEST))) &
            + 1./6.*(fi(m,directions(NORTH)) - ftmp(directions(SOUTH))) &
            - 1./3.*(fi(m,directions(NORTHEAST)) - ftmp(directions(SOUTHWEST))) &
            + 2./3.*fi(m,directions(NORTHWEST))
    enddo
  end subroutine DiscApplyBCDirichletToBoundary_D2Q9

  subroutine DiscApplyBCVelocityToBoundary_D2Q9(disc, fi, forces, fvals, &
       directions, cardinals, dist)
    type(discretization_type) disc
    type(distribution_type) dist
    PetscScalar,intent(inout),dimension(1:dist%s, 0:dist%b):: fi
    PetscScalar,intent(in),dimension(1:dist%s, dist%info%ndims):: forces
    PetscInt,intent(in),dimension(0:dist%b):: directions
    PetscScalar,intent(in),dimension(dist%s,dist%info%ndims):: fvals
    PetscInt,intent(in),dimension(1:dist%info%ndims):: cardinals

    PetscScalar rhotmp
    PetscScalar,dimension(0:dist%b)::ftmp
    PetscInt m
    ftmp = 0.0
    rhotmp = 0.0

    ! written for the NORTH boundary.
    do m=1,dist%s
       rhotmp = fi(m,directions(ORIGIN)) + fi(m,directions(EAST)) &
            + fi(m,directions(WEST)) + 2.*(fi(m,directions(NORTH)) &
            + fi(m,directions(NORTHEAST)) + fi(m,directions(NORTHWEST))) &
            + forces(m,cardinals(CARDINAL_NORMAL))/2.0
       rhotmp = rhotmp/(1. - fvals(1,cardinals(CARDINAL_NORMAL)))

       ! Choice should not affect the momentum significantly
       ftmp(directions(SOUTH)) = fi(m,directions(NORTH))
       ftmp(directions(SOUTHEAST)) = fi(m,directions(NORTHWEST))
       ftmp(directions(SOUTHWEST)) = fi(m,directions(NORTHEAST))

       fi(m,directions(SOUTH)) = 1./3.*ftmp(directions(SOUTH)) &
            + 2./3.*(rhotmp*fvals(1,cardinals(CARDINAL_NORMAL)) - forces(m,cardinals(CARDINAL_NORMAL))/2.) &
            - 2./3.*(ftmp(directions(SOUTHEAST)) + ftmp(directions(SOUTHWEST))) &
            + 2./3.*(fi(m,directions(NORTH)) + fi(m,directions(NORTHEAST)) + fi(m,directions(NORTHWEST)))

       fi(m,directions(SOUTHWEST)) = 1./3.*ftmp(directions(SOUTHWEST)) &
            - 1./2.*(rhotmp*fvals(1,cardinals(CARDINAL_CROSS)) - forces(m,cardinals(CARDINAL_CROSS))/2.) &
            + 1./6.*(rhotmp*fvals(1,cardinals(CARDINAL_NORMAL)) - forces(m,cardinals(CARDINAL_NORMAL))/2.) &
            + 1./2.*(fi(m,directions(EAST)) - fi(m,directions(WEST))) &
            + 1./6.*(fi(m,directions(NORTH)) - ftmp(directions(SOUTH))) &
            - 1./3.*(fi(m,directions(NORTHWEST)) - ftmp(directions(SOUTHEAST))) &
            + 2./3.*fi(m,directions(NORTHEAST))

       fi(m,directions(SOUTHEAST)) = 1./3.*ftmp(directions(SOUTHEAST)) &
            + 1./2.*(rhotmp*fvals(1,cardinals(CARDINAL_CROSS)) - forces(m,cardinals(CARDINAL_CROSS))/2.) &
            + 1./6.*(rhotmp*fvals(1,cardinals(CARDINAL_NORMAL)) - forces(m,cardinals(CARDINAL_NORMAL))/2.) &
            - 1./2.*(fi(m,directions(EAST)) - fi(m,directions(WEST))) &
            + 1./6.*(fi(m,directions(NORTH)) - ftmp(directions(SOUTH))) &
            - 1./3.*(fi(m,directions(NORTHEAST)) - ftmp(directions(SOUTHWEST))) &
            + 2./3.*fi(m,directions(NORTHWEST))

    enddo
    return
  end subroutine DiscApplyBCVelocityToBoundary_D2Q9

  subroutine DiscApplyBCFluxToBoundary_D2Q9(disc, fi, forces, fvals, &
       directions, cardinals, dist)
    type(discretization_type) disc
    type(distribution_type) dist
    PetscScalar,intent(inout),dimension(1:dist%s, 0:dist%b):: fi
    PetscScalar,intent(in),dimension(1:dist%s, dist%info%ndims):: forces
    PetscInt,intent(in),dimension(0:dist%b):: directions
    PetscScalar,intent(in),dimension(dist%s,dist%info%ndims):: fvals
    PetscInt,intent(in),dimension(1:dist%info%ndims):: cardinals

    PetscScalar rhotmp
    PetscScalar,dimension(0:dist%b)::ftmp
    PetscInt m
    PetscScalar :: twothirds,onethird,onesixth,onehalf

    ftmp = 0.0
    rhotmp = 0.0
    twothirds = 2.d0/3.d0
    onethird = 1.d0/3.d0
    onesixth = 1.d0/6.d0
    onehalf = 0.5d0

    ! written for the NORTH boundary.
    do m=1,dist%s
       ! Choice should not affect the momentum significantly
       ftmp(directions(SOUTH)) = fi(m,directions(NORTH))
       ftmp(directions(SOUTHEAST)) = fi(m,directions(NORTHWEST))
       ftmp(directions(SOUTHWEST)) = fi(m,directions(NORTHEAST))

       ! Chang et al '09, eqn 11
       fi(m,directions(SOUTH)) = ftmp(directions(SOUTH)) &
            + twothirds*(fvals(m,cardinals(CARDINAL_NORMAL)) &
                          - forces(m,cardinals(CARDINAL_NORMAL))/2.0) &
            + twothirds*(  fi(m,directions(NORTH)) - ftmp(directions(SOUTH)) &
                         + fi(m,directions(NORTHEAST)) - ftmp(directions(SOUTHWEST)) &
                         + fi(m,directions(NORTHWEST)) - ftmp(directions(SOUTHEAST)))

       ! Chang et al '09, eqn 13
       fi(m,directions(SOUTHWEST)) = ftmp(directions(SOUTHWEST)) &
            - onehalf*(fvals(m,cardinals(CARDINAL_CROSS)) - forces(m,cardinals(CARDINAL_CROSS))/2.0) &
            + onesixth*(fvals(m,cardinals(CARDINAL_NORMAL)) - forces(m,cardinals(CARDINAL_NORMAL))/2.0) &
            + onehalf*(fi(m,directions(EAST)) - fi(m,directions(WEST))) &
            + onesixth*(fi(m,directions(NORTH)) - ftmp(directions(SOUTH))) &
            - onethird*(fi(m,directions(NORTHWEST)) - ftmp(directions(SOUTHEAST))) &
            + twothirds*(fi(m,directions(NORTHEAST)) - ftmp(directions(SOUTHWEST)))

       ! Chang et al '09, eqn 12
       fi(m,directions(SOUTHEAST)) = ftmp(directions(SOUTHEAST)) &
            + onehalf*(fvals(m,cardinals(CARDINAL_CROSS)) - forces(m,cardinals(CARDINAL_CROSS))/2.0) &
            + onesixth*(fvals(m,cardinals(CARDINAL_NORMAL)) - forces(m,cardinals(CARDINAL_NORMAL))/2.0) &
            - onehalf*(fi(m,directions(EAST)) - fi(m,directions(WEST))) &
            + onesixth*(fi(m,directions(NORTH)) - ftmp(directions(SOUTH))) &
            + twothirds*(fi(m,directions(NORTHWEST)) - ftmp(directions(SOUTHEAST))) &
            - onethird*(fi(m,directions(NORTHEAST)) - ftmp(directions(SOUTHWEST)))

    enddo
    return
  end subroutine DiscApplyBCFluxToBoundary_D2Q9

  subroutine DiscSetLocalDirections_D2Q9(disc, boundary, directions, cardinals)
    type(discretization_type) disc
    PetscInt,intent(in):: boundary
    PetscInt,intent(out),dimension(0:discretization_directions) :: directions
    PetscInt,intent(out),dimension(1:discretization_dims) :: cardinals

    PetscInt i

    select case(boundary)
    case (BOUNDARY_YP)
       ! identity mapping
       do i=0,discretization_directions
          directions(i) = i
       end do
       cardinals(CARDINAL_NORMAL) = Y_DIRECTION
       cardinals(CARDINAL_CROSS) = X_DIRECTION

    case (BOUNDARY_YM)
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

       cardinals(CARDINAL_NORMAL) = Y_DIRECTION
       cardinals(CARDINAL_CROSS) = X_DIRECTION

    case (BOUNDARY_XP)
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

       cardinals(CARDINAL_NORMAL) = X_DIRECTION
       cardinals(CARDINAL_CROSS) = Y_DIRECTION

    case (BOUNDARY_XM)
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

       cardinals(CARDINAL_NORMAL) = X_DIRECTION
       cardinals(CARDINAL_CROSS) = Y_DIRECTION
    end select
  end subroutine DiscSetLocalDirections_D2Q9
end module LBM_Discretization_D2Q9_module

