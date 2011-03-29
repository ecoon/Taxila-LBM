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
  subroutine LBMEquilfFlow(feq, rho, ue, phase, dist)
    type(distribution_type) dist
    type(phase_type) phase(1:dist%s)

    PetscScalar,dimension(1:dist%s,0:dist%b,1:dist%info%gxyzl):: feq
    PetscScalar,dimension(1:dist%s,1:dist%info%gxyzl):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims,1:dist%info%gxyzl):: ue

    select case(dist%info%ndims)
    case(3)
       call LBMEquilfFlowD3(feq, rho, ue, phase, dist)
    case(2)
       call LBMEquilfFlowD2(feq, rho, ue, phase, dist)
    end select
  end subroutine LBMEquilfFlow

  subroutine LBMEquilfFlowD3(feq, rho, ue, phase, dist)
    type(distribution_type) dist
    type(phase_type) phase(1:dist%s)

    PetscScalar,dimension(1:dist%s,0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: feq
    PetscScalar,dimension(1:dist%s,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: ue

    PetscScalar usqr(1:dist%s)
    PetscScalar tmp
    PetscInt i,j,k,n,m
    
    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       usqr = SUM(ue(:,:,i,j,k)*ue(:,:,i,j,k), 2)
       do m=1,dist%s
          feq(m,0,i,j,k) = rho(m,i,j,k)*(phase(m)%alpha_0 + phase(m)%alpha_1*usqr(m))
       end do

       do n=1,dist%b
          do m=1,dist%s
             tmp = SUM(dist%disc%ci(n,:)*ue(m,:,i,j,k), 1)
             feq(m,n,i,j,k)= dist%disc%weights(n)*rho(m,i,j,k)* &
                  (1.5d0*(1.d0-phase(m)%d_k) + tmp/phase(m)%c_s2 &
                  + tmp*tmp/(2.d0*phase(m)%c_s2*phase(m)%c_s2) &
                  - usqr(m)/(2.d0*phase(m)%c_s2))
          end do
       end do
    end do
    end do
    end do
  end subroutine LBMEquilfFlowD3

  subroutine LBMEquilfFlowD2(feq, rho, ue, phase, dist)
    type(distribution_type) dist
    type(phase_type) phase(1:dist%s)

    PetscScalar,dimension(1:dist%s,0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: feq
    PetscScalar,dimension(1:dist%s,dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: rho
    PetscScalar,dimension(1:dist%s,1:dist%info%ndims, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: ue

    PetscScalar usqr(1:dist%s)
    PetscScalar tmp
    PetscInt i,j,n,m
    
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       usqr = SUM(ue(:,:,i,j)*ue(:,:,i,j), 2)
       do m=1,dist%s
          feq(m,0,i,j) = rho(m,i,j)*(phase(m)%alpha_0 + phase(m)%alpha_1*usqr(m))
       end do

       do n=1,dist%b
          do m=1,dist%s
             tmp = SUM(dist%disc%ci(n,:)*ue(m,:,i,j), 1)
             feq(m,n,i,j)= dist%disc%weights(n)*rho(m,i,j)* &
                  (1.5d0*(1.d0-phase(m)%d_k) + tmp/phase(m)%c_s2 &
                  + tmp*tmp/(2.d0*phase(m)%c_s2*phase(m)%c_s2) &
                  - usqr(m)/(2.d0*phase(m)%c_s2))
          end do
       end do
    end do
    end do
  end subroutine LBMEquilfFlowD2
end module LBM_Equilibrium_module
