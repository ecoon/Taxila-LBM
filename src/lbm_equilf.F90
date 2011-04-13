#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module LBM_Equilibrium_module
  use LBM_Phase_module
  use LBM_Distribution_Function_module
  implicit none
  private
#include "lbm_definitions.h"
  public :: LBMEquilfFlow

contains
  subroutine LBMEquilfFlow(feq, rho, u, forces, walls, phase, dist)
    type(distribution_type) dist
    type(phase_type) phase(1:dist%s)

    PetscScalar,dimension(dist%s,0:dist%b,dist%info%gxyzl):: feq
    PetscScalar,dimension(dist%s,dist%info%gxyzl):: rho
    PetscScalar,dimension(dist%s,dist%info%ndims,dist%info%gxyzl):: u,forces
    PetscScalar,dimension(dist%info%gxyzl):: walls

    select case(dist%info%ndims)
    case(3)
       call LBMEquilfFlowD3(feq, rho, u, forces, walls, phase, dist)
    case(2)
       call LBMEquilfFlowD2(feq, rho, u, forces, walls, phase, dist)
    end select
  end subroutine LBMEquilfFlow

  subroutine LBMEquilfFlowD3(feq, rho, u, forces, walls, phase, dist)
    type(distribution_type) dist
    type(phase_type) phase(1:dist%s)

    PetscScalar,dimension(1:dist%s,0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: feq
    PetscScalar,dimension(1:dist%s,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: u,forces
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls

    PetscInt i,j,k,d,n,m
    PetscScalar,parameter:: eps=1.e-15 ! slightly larger than machine epsilon
    PetscScalar,dimension(1:dist%s) :: mmot, mm, alpha_0, alpha_1, usqr
    PetscScalar :: up(1:dist%info%ndims)
    PetscScalar :: ue(1:dist%s,1:dist%info%ndims)
    PetscScalar tmp
    do m=1,dist%s
       mmot(m) = phase(m)%mm/phase(m)%tau
       mm(m) = phase(m)%mm
       alpha_0(m) = phase(m)%alpha_0
       alpha_1(m) = phase(m)%alpha_1
    end do

    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j,k).eq.0) then
       do d=1,dist%info%ndims
          up(d) = sum(u(:,d,i,j,k)*mmot,1)/sum(rho(:,i,j,k)*mmot,1)
       end do
       do m=1,dist%s
          ue(m,:) = up + forces(m,:,i,j,k)/(rho(m,i,j,k)*mmot(m)+eps)
       end do
       
       usqr = sum(ue*ue, 2)
       feq(:,0,i,j,k) = rho(:,i,j,k)*(alpha_0 + alpha_1*usqr)
       do n=1,dist%b
       do m=1,dist%s
          tmp = sum(dist%disc%ci(n,:)*ue(m,:), 1)
          feq(m,n,i,j,k)= dist%disc%weights(n)*rho(m,i,j,k)* &
               (1.5d0*(1.d0-phase(m)%d_k) + tmp/phase(m)%c_s2 &
               + tmp*tmp/(2.d0*phase(m)%c_s2*phase(m)%c_s2) &
               - usqr(m)/(2.d0*phase(m)%c_s2))
       end do
       end do
    end if
    end do
    end do
    end do
  end subroutine LBMEquilfFlowD3

  subroutine LBMEquilfFlowD2(feq, rho, u, forces, walls, phase, dist)
    type(distribution_type) dist
    type(phase_type) phase(1:dist%s)

    PetscScalar,dimension(1:dist%s,0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: feq
    PetscScalar,dimension(1:dist%s,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u,forces
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls

    PetscInt i,j,d,n,m
    PetscScalar,parameter:: eps=1.e-15 ! slightly larger than machine epsilon
    PetscScalar,dimension(1:dist%s) :: mmot, mm, alpha_0, alpha_1, usqr
    PetscScalar :: up(1:dist%info%ndims)
    PetscScalar :: ue(1:dist%s,1:dist%info%ndims)
    PetscScalar tmp
    do m=1,dist%s
       mmot(m) = phase(m)%mm/phase(m)%tau
       mm(m) = phase(m)%mm
       alpha_0(m) = phase(m)%alpha_0
       alpha_1(m) = phase(m)%alpha_1
    end do
    
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
    if (walls(i,j).eq.0) then
       do d=1,dist%info%ndims
          up(d) = sum(u(:,d,i,j)*mmot,1)/sum(rho(:,i,j)*mmot,1)
       end do
       do m=1,dist%s
          ue(m,:) = up + forces(m,:,i,j)/(rho(m,i,j)*mmot(m)+eps)
       end do
       usqr = sum(ue*ue, 2)
       feq(:,0,i,j) = rho(:,i,j)*(alpha_0 + alpha_1*usqr)

       do n=1,dist%b
       do m=1,dist%s
          tmp = sum(dist%disc%ci(n,:)*ue(m,:), 1)
          feq(m,n,i,j)= dist%disc%weights(n)*rho(m,i,j)* &
               (1.5d0*(1.d0-phase(m)%d_k) + tmp/phase(m)%c_s2 &
               + tmp*tmp/(2.d0*phase(m)%c_s2*phase(m)%c_s2) &
               - usqr(m)/(2.d0*phase(m)%c_s2))
       end do
       end do
    end if
    end do
    end do
  end subroutine LBMEquilfFlowD2
end module LBM_Equilibrium_module
