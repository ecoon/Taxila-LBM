!!!==================================================================
!!! Fortran-file
!!!    author:        Ethan T. Coon
!!!    filename:      get_forces.f
!!!    version:
!!!    created:       09 November 2010
!!!      on:          16:27:53 MST
!!!    last modified:  09 November 2010
!!!      at:          16:27:53 MST
!!!    URL:           http://www.ldeo.columbia.edu/~ecoon/
!!!    email:         ecoon _at_ ldeo.columbia.edu
!!!
!!!==================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

module LBM_Forcing_module
  use LBM_Distribution_Function_type_module
  use LBM_Distribution_Function_module
  use LBM_Phase_module
  use petsc
  implicit none

  private 
#include "lbm_definitions.h"

  public:: LBMAddFluidFluidForces, &
       LBMAddFluidSolidForces, &
       LBMAddBodyForces

contains
  ! --- Fluid-fluid interaction forces, from
  ! ---  (Kang 2002 Eq. 6)
  subroutine LBMAddFluidFluidForces(dist, phases, rho, walls, forces)
    !     NONLOCAL IN RHO
    type(distribution_type) dist
    type(phase_type) phases(dist%s)
    PetscScalar,dimension(1:dist%s, 1:dist%info%gxyzl):: rho
    PetscScalar,dimension(1:dist%s, 1:dist%info%ndims, 1:dist%info%gxyzl):: forces
    PetscScalar,dimension(1:dist%info%gxyzl):: walls

    select case(dist%info%ndims)
    case(2)
       call LBMAddFluidFluidForcesD2(dist, phases, rho, walls, forces)
    case(3)
       call LBMAddFluidFluidForcesD3(dist, phases, rho, walls, forces)
    end select
  end subroutine LBMAddFluidFluidForces

  subroutine LBMAddFluidFluidForcesD3(dist, phases, rho, walls, forces)
    type(distribution_type) dist
    type(phase_type) phases(dist%s)
    PetscScalar,dimension(dist%s, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: rho
    PetscScalar,dimension(1:dist%s, 1:dist%info%ndims,  dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: forces
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls

    ! local
    PetscInt i,j,k,m,n,d
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: tmp
    PetscScalar,dimension(dist%s, dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: gradrho
    PetscScalar,dimension(dist%info%ndims):: weightsum
    PetscErrorCode ierr

    gradrho = 0.d0
    weightsum = 0.d0

    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if (walls(i,j,k).eq.0) then
          weightsum = 0.d0

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!! 4th order isotropy !!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! N/S/E/W sqared length 1 
          if (walls(i+1,j,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(1)*(rho(:,i+1,j,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(1)
          end if
          if (walls(i-1,j,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(1)*(rho(:,i-1,j,k)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(1)
          end if
          if (walls(i,j+1,k).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(1)*(rho(:,i,j+1,k)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(1)
          end if
          if (walls(i,j-1,k).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(1)*(rho(:,i,j-1,k)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(1)
          end if
          if (walls(i,j,k+1).eq.0) then
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(1)*(rho(:,i,j,k+1)-rho(:,i,j,k)) 
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(1)
          end if
          if (walls(i,j,k-1).eq.0) then
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(1)*(rho(:,i,j,k-1)-rho(:,i,j,k)) 
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(1)
          end if
   
          ! diagonals squared length 2 (x and y)
          if (walls(i+1,j+1,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i+1,j+1,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i+1,j+1,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i-1,j+1,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i-1,j+1,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i-1,j+1,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i+1,j-1,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i+1,j-1,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i+1,j-1,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i-1,j-1,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i-1,j-1,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i-1,j-1,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
          end if 

          ! diagonals squared length 2 (x and z)
          if (walls(i+1,j,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i+1,j,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i+1,j,k+1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i-1,j,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i-1,j,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i-1,j,k+1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i+1,j,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i+1,j,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i+1,j,k-1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i-1,j,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i-1,j,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i-1,j,k-1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(2)
          end if 

          ! diagonals squared length 2 (y and z)
          if (walls(i,j+1,k+1).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i,j+1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i,j+1,k+1)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i,j-1,k+1).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i,j-1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i,j-1,k+1)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i,j+1,k-1).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(2)*(rho(:,i,j+1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i,j+1,k-1)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i,j-1,k-1).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i,j-1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(2)*(rho(:,i,j-1,k-1)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(2)
          end if 

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!! 8th order isotropy !!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(dist%info%stencil_size > 1) then

          ! diagonals squared length 3
          if (walls(i+1,j+1,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i+1,j+1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i+1,j+1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i+1,j+1,k+1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(3) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(3)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(3)
          end if
          if (walls(i-1,j+1,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i-1,j+1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i-1,j+1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i-1,j+1,k+1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(3) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(3)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(3)
          end if 
          if (walls(i-1,j-1,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i-1,j-1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i-1,j-1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i-1,j-1,k+1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(3) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(3)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(3)
          end if
          if (walls(i+1,j-1,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i+1,j-1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i+1,j-1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i+1,j-1,k+1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(3) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(3)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(3)
          end if
          if (walls(i+1,j+1,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i+1,j+1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i+1,j+1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i+1,j+1,k-1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(3) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(3)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(3)
          end if
          if (walls(i-1,j+1,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i-1,j+1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i-1,j+1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i-1,j+1,k-1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(3) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(3)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(3)
          end if
          if (walls(i-1,j-1,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i-1,j-1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i-1,j-1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i-1,j-1,k-1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(3) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(3)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(3)
          end if
          if (walls(i+1,j-1,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(3)*(rho(:,i+1,j-1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i+1,j-1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(3)*(rho(:,i+1,j-1,k-1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(3) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(3)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(3)
          end if

          ! diagonals squared length 4
          if (walls(i+2,j,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(4)*(rho(:,i+2,j,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(4)*4. 
           end if
           if (walls(i-2,j,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(4)*(rho(:,i-2,j,k)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(4)*4.
          end if
          if (walls(i,j+2,k).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(4)*(rho(:,i,j+2,k)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(4)*4.
          end if
          if (walls(i,j-2,k).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(4)*(rho(:,i,j-2,k)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(4)*4.
          end if
          if (walls(i,j,k+2).eq.0) then
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(4)*(rho(:,i,j,k+2)-rho(:,i,j,k)) 
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(4)*4.
          end if
          if (walls(i,j,k-2).eq.0) then
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(4)*(rho(:,i,j,k-2)-rho(:,i,j,k)) 
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(4)*4.
          end if

          ! diagonals squared length 5 (x: 2 and -2)
          if (walls(i+2,j+1,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i+2,j+1,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i+2,j+1,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
          end if
          if (walls(i-2,j+1,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i-2,j+1,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i-2,j+1,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
          end if
          if (walls(i-2,j-1,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i-2,j-1,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i-2,j-1,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
          end if
          if (walls(i+2,j-1,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i+2,j-1,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i+2,j-1,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
          end if
          if (walls(i+2,j,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i+2,j,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i+2,j,k+1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)
          end if
          if (walls(i-2,j,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i-2,j,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i-2,j,k+1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)
          end if
          if (walls(i-2,j,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i-2,j,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i-2,j,k-1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)
          end if
          if (walls(i+2,j,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i+2,j,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i+2,j,k-1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)
          end if

          ! diagonals squared length 5 (y: 2 and -2)
          if (walls(i+1,j+2,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i+1,j+2,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i+1,j+2,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
          end if
          if (walls(i-1,j+2,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i-1,j+2,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i-1,j+2,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
          end if
          if (walls(i-1,j-2,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i-1,j-2,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i-1,j-2,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
          end if
          if (walls(i+1,j-2,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i+1,j-2,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i+1,j-2,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
          end if
          if (walls(i,j+2,k+1).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i,j+2,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i,j+2,k+1)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)
          end if
          if (walls(i,j-2,k+1).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i,j-2,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i,j-2,k+1)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)
          end if
          if (walls(i,j+2,k-1).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i,j+2,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i,j+2,k-1)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)
          end if
          if (walls(i,j-2,k-1).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i,j-2,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i,j-2,k-1)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)
          end if

          ! diagonals squared length 5 (z: 2 and -2)
          if (walls(i+1,j,k+2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i+1,j,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i+1,j,k+2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)*4.
          end if
          if (walls(i,j+1,k+2).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i,j+1,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i,j+1,k+2)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)*4.
          end if
          if (walls(i-1,j,k+2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i-1,j,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i-1,j,k+2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)*4.
          end if
          if (walls(i,j-1,k+2).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i,j-1,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(5)*(rho(:,i,j-1,k+2)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)*4.
          end if

          if (walls(i+1,j,k-2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i+1,j,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i+1,j,k-2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)*4.
          end if
          if (walls(i,j+1,k-2).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(5)*(rho(:,i,j+1,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i,j+1,k-2)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)*4.
          end if
          if (walls(i-1,j,k-2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i-1,j,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i-1,j,k-2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)*4.
          end if
          if (walls(i,j-1,k-2).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(5)*(rho(:,i,j-1,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(5)*(rho(:,i,j-1,k-2)-rho(:,i,j,k)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(5)*4.
          end if

          ! diagonals squared length 6 (x: 2 and -2)
          if (walls(i+2,j+1,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i+2,j+1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+2,j+1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+2,j+1,k+1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i-2,j+1,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i-2,j+1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i-2,j+1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i-2,j+1,k+1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i-2,j-1,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i-2,j-1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-2,j-1,k+1)-rho(:,i,j,k))
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i-2,j-1,k+1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i+2,j-1,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i+2,j-1,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i+2,j-1,k+1)-rho(:,i,j,k))
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+2,j-1,k+1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i+2,j+1,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i+2,j+1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+2,j+1,k-1)-rho(:,i,j,k))
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i+2,j+1,k-1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i-2,j+1,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i-2,j+1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i-2,j+1,k-1)-rho(:,i,j,k))
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-2,j+1,k-1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i-2,j-1,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i-2,j-1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-2,j-1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-2,j-1,k-1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i+2,j-1,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i+2,j-1,k-1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i+2,j-1,k-1)-rho(:,i,j,k))
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i+2,j-1,k-1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if

          ! diagonals squared length 6 (y: 2 and -2)
          if (walls(i+1,j+2,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j+2,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i+1,j+2,k+1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j+2,k+1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i-1,j+2,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j+2,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i-1,j+2,k+1)-rho(:,i,j,k))
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i-1,j+2,k+1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i-1,j-2,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j-2,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i-1,j-2,k+1)-rho(:,i,j,k))
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i-1,j-2,k+1)-rho(:,i,j,k))  
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i+1,j-2,k+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j-2,k+1)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i+1,j-2,k+1)-rho(:,i,j,k))
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j-2,k+1)-rho(:,i,j,k))  
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i+1,j+2,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j+2,k-1)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i+1,j+2,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i+1,j+2,k-1)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i-1,j-2,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j-2,k-1)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i-1,j-2,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j-2,k-1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i-1,j+2,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j+2,k-1)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i-1,j+2,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j+2,k-1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if
          if (walls(i+1,j-2,k-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j-2,k-1)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i+1,j-2,k-1)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i+1,j-2,k-1)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)
          end if

          ! diagonals squared length 6 (z: 2 and -2)
          if (walls(i+1,j+1,k+2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j+1,k+2)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j+1,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i+1,j+1,k+2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)*4.
          end if
          if (walls(i-1,j+1,k+2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,x_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j+1,k+2)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i-1,j+1,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i-1,j+1,k+2)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)*4.
          end if
          if (walls(i-1,j-1,k+2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j-1,k+2)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j-1,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i-1,j-1,k+2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)*4.
          end if
          if (walls(i+1,j-1,k+2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j-1,k+2)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i+1,j-1,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(6)*(rho(:,i+1,j-1,k+2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)*4.
          end if
          if (walls(i+1,j+1,k-2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j+1,k-2)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j+1,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i+1,j+1,k-2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)*4.
          end if
          if (walls(i-1,j+1,k-2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j+1,k-2)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i-1,j+1,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i-1,j+1,k-2)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)*4.
          end if
          if (walls(i-1,j-1,k-2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j-1,k-2)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i-1,j-1,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i-1,j-1,k-2)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)*4.
          end if
          if (walls(i+1,j-1,k-2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + dist%disc%ffw(6)*(rho(:,i+1,j-1,k-2)-rho(:,i,j,k))
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - dist%disc%ffw(6)*(rho(:,i+1,j-1,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(6)*(rho(:,i+1,j-1,k-2)-rho(:,i,j,k))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(6) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(6)
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(6)*4.
          end if

          ! diagonals squared length 8 (x and y)
          if (walls(i+2,j+2,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i+2,j+2,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i+2,j+2,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
          end if
          if (walls(i-2,j+2,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i-2,j+2,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i-2,j+2,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
          end if
          if (walls(i-2,j-2,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i-2,j-2,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i-2,j-2,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
          end if
          if (walls(i+2,j-2,k).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i+2,j-2,k)-rho(:,i,j,k)) 
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i+2,j-2,k)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
          end if

          ! diagonals squared length 8 (x and z)
          if (walls(i+2,j,k+2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i+2,j,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i+2,j,k+2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(8)*4.
          end if
          if (walls(i-2,j,k+2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i-2,j,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i-2,j,k+2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(8)*4.
          end if
          if (walls(i-2,j,k-2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i-2,j,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i-2,j,k-2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(8)*4.
          end if
          if (walls(i+2,j,k-2).eq.0) then
             gradrho(:,X_DIRECTION,i,j,k) = gradrho(:,X_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i+2,j,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i+2,j,k-2)-rho(:,i,j,k)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(8)*4.
          end if

          ! diagonals squared length 8 (y and z)
          if (walls(i,j+2,k+2).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i,j+2,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i,j+2,k+2)-rho(:,i,j,k))
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(8)*4.
          end if
          if (walls(i,j-2,k+2).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i,j-2,k+2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i,j-2,k+2)-rho(:,i,j,k))
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(8)*4.
          end if
          if (walls(i,j-2,k-2).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i,j-2,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i,j-2,k-2)-rho(:,i,j,k))
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(8)*4.
          end if
          if (walls(i,j+2,k-2).eq.0) then
             gradrho(:,Y_DIRECTION,i,j,k) = gradrho(:,Y_DIRECTION,i,j,k) &
                  + 2.*dist%disc%ffw(8)*(rho(:,i,j+2,k-2)-rho(:,i,j,k)) 
             gradrho(:,Z_DIRECTION,i,j,k) = gradrho(:,Z_DIRECTION,i,j,k) &
                  - 2.*dist%disc%ffw(8)*(rho(:,i,j+2,k-2)-rho(:,i,j,k))
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
             weightsum(Z_DIRECTION) = weightsum(Z_DIRECTION) + dist%disc%ffw(8)*4.
          end if

          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!! Calc. FF force using gradrho !!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do m=1,dist%s
          do d=1,dist%info%ndims
             forces(m,d,i,j,k) = -rho(m,i,j,k)*sum(phases(m)%gf &
                  *(gradrho(:,d,i,j,k)/weightsum(d)),1)
          end do
          end do
       end if
    
    end do
    end do
    end do
    
  end subroutine LBMAddFluidFluidForcesD3

  subroutine LBMAddFluidFluidForcesD2(dist, phases, rho, walls, forces)
    !     NONLOCAL IN RHO

    ! input
    type(distribution_type) dist
    type(phase_type) phases(dist%s)
    PetscScalar,dimension(dist%s, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: rho
    PetscScalar,dimension(1:dist%s, 1:dist%info%ndims,  dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: forces
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls

    ! local
    PetscInt i,j,m,n,d
    PetscScalar,dimension(dist%s, dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: gradrho
    PetscScalar,dimension(dist%info%ndims):: weightsum
    PetscErrorCode ierr

    gradrho = 0.d0
    weightsum = 0.d0

    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if (walls(i,j).eq.0) then
          weightsum = 0.d0

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!! 4th order isotropy !!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! N/S/E/W sqared length 1 
          if (walls(i+1,j).eq.0) then
             gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                  + dist%disc%ffw(1)*(rho(:,i+1,j)-rho(:,i,j)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(1)
          end if
          if (walls(i-1,j).eq.0) then
             gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                  - dist%disc%ffw(1)*(rho(:,i-1,j)-rho(:,i,j))
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(1)
          end if
          if (walls(i,j+1).eq.0) then
             gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                  + dist%disc%ffw(1)*(rho(:,i,j+1)-rho(:,i,j)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(1)
          end if
          if (walls(i,j-1).eq.0) then
             gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                  - dist%disc%ffw(1)*(rho(:,i,j-1)-rho(:,i,j)) 
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(1)
          end if

          ! diagonals squared length 2
          if (walls(i+1,j+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                  + dist%disc%ffw(2)*(rho(:,i+1,j+1)-rho(:,i,j)) 
             gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                  + dist%disc%ffw(2)*(rho(:,i+1,j+1)-rho(:,i,j)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i-1,j+1).eq.0) then
             gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                  - dist%disc%ffw(2)*(rho(:,i-1,j+1)-rho(:,i,j)) 
             gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                  + dist%disc%ffw(2)*(rho(:,i-1,j+1)-rho(:,i,j)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i+1,j-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                  + dist%disc%ffw(2)*(rho(:,i+1,j-1)-rho(:,i,j)) 
             gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                  - dist%disc%ffw(2)*(rho(:,i+1,j-1)-rho(:,i,j)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
          end if
          if (walls(i-1,j-1).eq.0) then
             gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                  - dist%disc%ffw(2)*(rho(:,i-1,j-1)-rho(:,i,j)) 
             gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                  - dist%disc%ffw(2)*(rho(:,i-1,j-1)-rho(:,i,j)) 
             weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(2)
             weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(2)
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!! 8th order isotropy !!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(dist%info%stencil_size > 1) then
             ! N/S/E/W squared length 4 
             if (walls(i+2,j).eq.0.and.walls(i+1,j).eq.0) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + 2.*dist%disc%ffw(4)*(rho(:,i+2,j)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(4)*4.
             end if
             if (walls(i-2,j).eq.0.and.walls(i-1,j).eq.0) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - 2.*dist%disc%ffw(4)*(rho(:,i-2,j)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(4)*4.
             end if
             if (walls(i,j+2).eq.0.and.walls(i,j+1).eq.0) then
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + 2.*dist%disc%ffw(4)*(rho(:,i,j+2)-rho(:,i,j)) 
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(4)*4.
             end if
             if (walls(i,j-2).eq.0.and.walls(i,j-1).eq.0) then
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - 2.*dist%disc%ffw(4)*(rho(:,i,j-2)-rho(:,i,j)) 
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(4)*4.
             end if

             ! short diagonals, length squared 5
             if (walls(i+2,j+1).eq.0.and.(walls(i+1,j+1).eq.0.or.walls(i+1,j).eq.0)) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + 2.*dist%disc%ffw(5)*(rho(:,i+2,j+1)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + dist%disc%ffw(5)*(rho(:,i+2,j+1)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
             end if
             if (walls(i+2,j-1).eq.0.and.(walls(i+1,j-1).eq.0.or.walls(i+1,j).eq.0)) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + 2.*dist%disc%ffw(5)*(rho(:,i+2,j-1)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - dist%disc%ffw(5)*(rho(:,i+2,j-1)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
             end if
             if (walls(i-2,j+1).eq.0.and.(walls(i-1,j+1).eq.0.or.walls(i-1,j).eq.0)) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - 2.*dist%disc%ffw(5)*(rho(:,i-2,j+1)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + dist%disc%ffw(5)*(rho(:,i-2,j+1)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
             end if
             if (walls(i-2,j-1).eq.0.and.(walls(i-1,j-1).eq.0.or.walls(i-1,j).eq.0)) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - 2.*dist%disc%ffw(5)*(rho(:,i-2,j-1)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - dist%disc%ffw(5)*(rho(:,i-2,j-1)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)*4.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)
             end if
             if (walls(i+1,j+2).eq.0.and.(walls(i+1,j+1).eq.0.or.walls(i,j+1).eq.0)) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + dist%disc%ffw(5)*(rho(:,i+1,j+2)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + 2.*dist%disc%ffw(5)*(rho(:,i+1,j+2)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
             end if
             if (walls(i+1,j-2).eq.0.and.(walls(i+1,j-1).eq.0.or.walls(i,j-1).eq.0)) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + dist%disc%ffw(5)*(rho(:,i+1,j-2)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - 2.*dist%disc%ffw(5)*(rho(:,i+1,j-2)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
             end if
             if (walls(i-1,j+2).eq.0.and.(walls(i-1,j+1).eq.0.or.walls(i,j+1).eq.0)) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - dist%disc%ffw(5)*(rho(:,i-1,j+2)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + 2.*dist%disc%ffw(5)*(rho(:,i-1,j+2)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
             end if
             if (walls(i-1,j-2).eq.0.and.(walls(i-1,j-1).eq.0.or.walls(i,j-1).eq.0)) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - dist%disc%ffw(5)*(rho(:,i-1,j-2)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - 2.*dist%disc%ffw(5)*(rho(:,i-1,j-2)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(5)
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(5)*4.
             end if

             ! long diagonals, length squared = 8
             if (walls(i+2,j+2).eq.0.and.walls(i+1,j+1).eq.0) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + 2.*dist%disc%ffw(8)*(rho(:,i+2,j+2)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + 2.*dist%disc%ffw(8)*(rho(:,i+2,j+2)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
             end if
             if (walls(i-2,j+2).eq.0.and.walls(i-1,j+1).eq.0) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - 2.*dist%disc%ffw(8)*(rho(:,i-2,j+2)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + 2.*dist%disc%ffw(8)*(rho(:,i-2,j+2)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
             end if
             if (walls(i+2,j-2).eq.0.and.walls(i+1,j-1).eq.0) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + 2.*dist%disc%ffw(8)*(rho(:,i+2,j-2)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - 2.*dist%disc%ffw(8)*(rho(:,i+2,j-2)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
             end if
             if (walls(i-2,j-2).eq.0.and.walls(i-1,j-1).eq.0) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - 2.*dist%disc%ffw(8)*(rho(:,i-2,j-2)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - 2.*dist%disc%ffw(8)*(rho(:,i-2,j-2)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(8)*4.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(8)*4.
             end if
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!! 10th order isotropy !!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(dist%info%stencil_size > 2) then
             ! N/S/E/W squared length 9 
             if (walls(i+3,j).eq.0.and.walls(i+2,j).eq.0.and.walls(i+1,j).eq.0) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + 3.*dist%disc%ffw(9)*(rho(:,i+3,j)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(9)*9.
             end if
             if (walls(i-3,j).eq.0.and.walls(i-2,j).eq.0.and.walls(i-1,j).eq.0) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - 3.*dist%disc%ffw(9)*(rho(:,i-3,j)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(9)*9.
             end if
             if (walls(i,j+3).eq.0.and.walls(i,j+2).eq.0.and.walls(i,j+1).eq.0) then
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + 3.*dist%disc%ffw(9)*(rho(:,i,j+3)-rho(:,i,j)) 
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(9)*9.
             end if
             if (walls(i,j-3).eq.0.and.walls(i,j-2).eq.0.and.walls(i,j-1).eq.0) then
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - 3.*dist%disc%ffw(9)*(rho(:,i,j-3)-rho(:,i,j)) 
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(9)*9.
             end if

             ! short diagonals, length squared 10
             if (walls(i+3,j+1).eq.0.and.((walls(i+1,j).eq.0.and.walls(i+2,j).eq.0) &
                                      .or.(walls(i+1,j).eq.0.and.walls(i+2,j).eq.0) &
                                      .or.(walls(i+1,j+1).eq.0.and.walls(i+2,j+1).eq.0))) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + 3.*dist%disc%ffw(10)*(rho(:,i+3,j+1)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + dist%disc%ffw(10)*(rho(:,i+3,j+1)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(10)*9.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(10)
             end if
             if (walls(i+3,j-1).eq.0.and.((walls(i+1,j).eq.0.and.walls(i+2,j).eq.0) &
                                      .or.(walls(i+1,j).eq.0.and.walls(i+2,j).eq.0) &
                                      .or.(walls(i+1,j-1).eq.0.and.walls(i+2,j-1).eq.0))) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + 3.*dist%disc%ffw(10)*(rho(:,i+3,j-1)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - dist%disc%ffw(10)*(rho(:,i+3,j-1)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(10)*9.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(10)
             end if
             if (walls(i-3,j+1).eq.0.and.((walls(i-1,j).eq.0.and.walls(i-2,j).eq.0) &
                                      .or.(walls(i-1,j).eq.0.and.walls(i-2,j).eq.0) &
                                      .or.(walls(i-1,j+1).eq.0.and.walls(i-2,j+1).eq.0))) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - 3.*dist%disc%ffw(10)*(rho(:,i-3,j+1)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + dist%disc%ffw(10)*(rho(:,i-3,j+1)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(10)*9.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(10)
             end if
             if (walls(i-3,j-1).eq.0.and.((walls(i-1,j).eq.0.and.walls(i-2,j).eq.0) &
                                      .or.(walls(i-1,j).eq.0.and.walls(i-2,j).eq.0) &
                                      .or.(walls(i-1,j-1).eq.0.and.walls(i-2,j-1).eq.0))) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - 3.*dist%disc%ffw(10)*(rho(:,i-3,j-1)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - dist%disc%ffw(10)*(rho(:,i-3,j-1)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(10)*9.
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(10)
             end if
             if (walls(i+1,j+3).eq.0.and.((walls(i,j+1).eq.0.and.walls(i,j+2).eq.0) &
                                      .or.(walls(i,j+1).eq.0.and.walls(i,j+2).eq.0) &
                                      .or.(walls(i+1,j+1).eq.0.and.walls(i+1,j+2).eq.0))) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + dist%disc%ffw(10)*(rho(:,i+1,j+3)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + 3.*dist%disc%ffw(10)*(rho(:,i+1,j+3)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(10)
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(10)*9.
             end if
             if (walls(i+1,j-3).eq.0.and.((walls(i,j-1).eq.0.and.walls(i,j-2).eq.0) &
                                      .or.(walls(i,j-1).eq.0.and.walls(i,j-2).eq.0) &
                                      .or.(walls(i+1,j-1).eq.0.and.walls(i+1,j-2).eq.0))) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     + dist%disc%ffw(10)*(rho(:,i+1,j-3)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - 3.*dist%disc%ffw(10)*(rho(:,i+1,j-3)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(10)
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(10)*9.
             end if
             if (walls(i-1,j+3).eq.0.and.((walls(i,j+1).eq.0.and.walls(i,j+2).eq.0) &
                                      .or.(walls(i,j+1).eq.0.and.walls(i,j+2).eq.0) &
                                      .or.(walls(i-1,j+1).eq.0.and.walls(i-1,j+2).eq.0))) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - dist%disc%ffw(10)*(rho(:,i-1,j+3)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     + 3.*dist%disc%ffw(10)*(rho(:,i-1,j+3)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(10)
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(10)*9.
             end if
             if (walls(i-1,j-3).eq.0.and.((walls(i,j-1).eq.0.and.walls(i,j-2).eq.0) &
                                      .or.(walls(i,j-1).eq.0.and.walls(i,j-2).eq.0) &
                                      .or.(walls(i-1,j-1).eq.0.and.walls(i-1,j-2).eq.0))) then
                gradrho(:,X_DIRECTION,i,j) = gradrho(:,X_DIRECTION,i,j) &
                     - dist%disc%ffw(10)*(rho(:,i-1,j-3)-rho(:,i,j)) 
                gradrho(:,Y_DIRECTION,i,j) = gradrho(:,Y_DIRECTION,i,j) &
                     - 3.*dist%disc%ffw(10)*(rho(:,i-1,j-3)-rho(:,i,j)) 
                weightsum(X_DIRECTION) = weightsum(X_DIRECTION) + dist%disc%ffw(10)
                weightsum(Y_DIRECTION) = weightsum(Y_DIRECTION) + dist%disc%ffw(10)*9.
             end if
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!! Calc. FF force using gradrho !!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do m=1,dist%s
          do d=1,dist%info%ndims
             forces(m,d,i,j) = -rho(m,i,j)*sum(phases(m)%gf &
                  *(gradrho(:,d,i,j)/weightsum(d)),1)
          end do
          end do
       end if
    end do
    end do

  end subroutine LBMAddFluidFluidForcesD2

  ! --- Fluid-solid interaction forces, from
  ! ---  (Kang 2002 Eq. 8)
  ! --- NOTE -- I'm ignoring the bug that the forces will
  !             wrap even in the non-periodic case.  This seems
  !             save at this point, since it seems unlikely that
  !             there will be walls on one side but not the other
  ! New code added by MLP.  This accounts for fluid-solid forces on the 
  ! diagonals, which is not accounted for in the other version.
 
  subroutine LBMAddFluidSolidForces(dist, phases, rho, walls, forces)
    !     NONLOCAL IN WALLS
    type(distribution_type) dist
    type(phase_type) phases(dist%s)
    PetscScalar,dimension(1:dist%s, 1:dist%info%gxyzl):: rho
    PetscScalar,dimension(1:dist%s, 1:dist%info%ndims, 1:dist%info%gxyzl):: forces
    PetscScalar,dimension(1:dist%info%gxyzl):: walls

    select case(dist%info%ndims)
    case(2)
       call LBMAddFluidSolidForcesD2(dist, phases, rho, walls, forces)
    case(3)
       call LBMAddFluidSolidForcesD3(dist, phases, rho, walls, forces)
    end select
  end subroutine LBMAddFluidSolidForces
    
  subroutine LBMAddFluidSolidForcesD3(dist, phases, rho, walls, forces)
    type(distribution_type) dist
    type(phase_type) phases(dist%s)
    PetscScalar,dimension(dist%s, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: rho
    PetscScalar,dimension(1:dist%s, 1:dist%info%ndims,  dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: forces
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls

    ! local
    PetscInt i,j,k,m,n,d
    PetscScalar,dimension(0:dist%b,  dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: tmp
    PetscScalar gw(1:dist%s)
    PetscErrorCode ierr

    call DistributionGatherValueToDirection(dist, walls, tmp)

    do m=1,dist%s
       gw(m) = phases(m)%gw
    end do
    
    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j,k).eq.0) then
       do d=1,dist%info%ndims
          do n=1,2*dist%info%ndims
             forces(:,d,i,j,k) = forces(:,d,i,j,k) &
                  - rho(:,i,j,k)*tmp(n,i,j,k)*dist%disc%ci(n,d)*gw
          end do
          do n=2*dist%info%ndims+1,dist%b
             forces(:,d,i,j,k) = forces(:,d,i,j,k) &
                  - 0.5*rho(:,i,j,k)*tmp(n,i,j,k)*dist%disc%ci(n,d)*gw
          end do
       end do
    end if
    end do
    end do
    end do
  end subroutine LBMAddFluidSolidForcesD3

  subroutine LBMAddFluidSolidForcesD2(dist, phases, rho, walls, forces)
    type(distribution_type) dist
    type(phase_type) phases(dist%s)
    PetscScalar,dimension(dist%s, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: rho
    PetscScalar,dimension(1:dist%s, 1:dist%info%ndims,  dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: forces
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls

    ! local
    PetscInt i,j,m,n,d
    PetscScalar,dimension(0:dist%b,  dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: tmp
    PetscScalar gw(1:dist%s)
    PetscErrorCode ierr

    call DistributionGatherValueToDirection(dist, walls, tmp)

    do m=1,dist%s
       gw(m) = phases(m)%gw
    end do
    
    ! The weights in here were 1 for the neighbors and 0.25 for the next 
    ! nearet neighbors without the HOD in FF.  Now I am setting the weights
    ! to 1/3 and 1/12 to be consisitent with the weight in the 4th order 
    ! isotropy HOD.

    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j).eq.0) then
       do d=1,dist%info%ndims
          do n=1,2*dist%info%ndims
             forces(:,d,i,j) = forces(:,d,i,j) &
                  - 1./3.*rho(:,i,j)*tmp(n,i,j)*dist%disc%ci(n,d)*gw
          end do
          do n=2*dist%info%ndims+1,dist%b
             forces(:,d,i,j) = forces(:,d,i,j) &
                  - 1./12.*rho(:,i,j)*tmp(n,i,j)*dist%disc%ci(n,d)*gw
          end do
       end do
    end if
    end do
    end do
  end subroutine LBMAddFluidSolidForcesD2

  ! --- body forces on fluid
  subroutine LBMAddBodyForces(dist, phases, gvt, rho, walls, forces)
    type(distribution_type) dist
    type(phase_type) phases(1:dist%s)
    PetscScalar,dimension(dist%s, dist%info%gxyzl):: rho
    PetscScalar,dimension(dist%s, dist%info%ndims, 1:dist%info%gxyzl):: forces
    PetscScalar,dimension(dist%info%gxyzl):: walls
    PetscScalar,dimension(dist%info%ndims) :: gvt

    select case(dist%info%ndims)
    case(2)
       call LBMAddBodyForcesD2(dist, phases, gvt, rho, walls, forces)
    case(3)
       call LBMAddBodyForcesD3(dist, phases, gvt, rho, walls, forces)
    end select
  end subroutine LBMAddBodyForces
    
  subroutine LBMAddBodyForcesD3(dist, phases, gvt, rho, walls, forces)
    type(distribution_type) dist
    type(phase_type) phases(dist%s)
    PetscScalar,dimension(dist%s, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: rho
    PetscScalar,dimension(1:dist%s, 1:dist%info%ndims,  dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: forces
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls
    PetscScalar,dimension(dist%info%ndims) :: gvt

    ! local
    PetscInt i,j,k,m
    PetscErrorCode ierr

    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if (walls(i,j,k).eq.0) then
          do m=1,dist%s
             forces(m,:,i,j,k) = forces(m,:,i,j,k) &
                  + gvt*phases(m)%mm*rho(m,i,j,k)
          end do
       end if
    end do
    end do
    end do
    return
  end subroutine LBMAddBodyForcesD3

  subroutine LBMAddBodyForcesD2(dist, phases, gvt, rho, walls, forces)
    type(distribution_type) dist
    type(phase_type) phases(dist%s)
    PetscScalar,dimension(dist%s, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: rho
    PetscScalar,dimension(1:dist%s, 1:dist%info%ndims,  dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: forces
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls
    PetscScalar,dimension(dist%info%ndims) :: gvt

    ! local
    PetscInt i,j,m
    PetscErrorCode ierr

    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if (walls(i,j).eq.0) then
          do m=1,dist%s
             forces(m,:,i,j) = forces(m,:,i,j) &
                  + gvt*phases(m)%mm*rho(m,i,j)
          end do
       end if
    end do
    end do
    return
  end subroutine LBMAddBodyForcesD2
end module LBM_Forcing_module
