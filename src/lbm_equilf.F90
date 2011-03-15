#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module LBM_Equilibrium_module
  use LBM_Info_module
  use LBM_Constants_module
  implicit none
  private
#include "lbm_definitions.h"
  public :: LBMEquilf

contains
  subroutine LBMEquilf(feq, rho, ue, info, constants)
    use petsc
    type(info_type) info
    type(constants_type) constants
    PetscScalar,dimension(1:info%s,0:info%flow_disc%b)::feq
    PetscScalar,dimension(1:info%s):: rho
    PetscScalar,dimension(1:info%s,1:info%ndims):: ue
    PetscErrorCode ierr
    
    PetscScalar usqr(1:info%s),tmp(1:info%s)
    PetscInt i,m
    
    usqr = SUM(ue*ue, 2)
    feq(:,0) = rho*(constants%alpha_0 + constants%alpha_1*usqr)

    do i=1,info%flow_disc%b
       do m=1,info%s
          tmp(m) = SUM(info%flow_disc%ci(i,:)*ue(m,:), 1)
       end do
       feq(:,i)= info%flow_disc%weights(i)*rho* &
            (1.5d0*(1.d0-constants%d_k) + tmp/constants%c_s2 &
            + tmp*tmp/(2.d0*constants%c_s2*constants%c_s2) &
            - usqr/(2.d0*constants%c_s2))
    end do
  end subroutine LBMEquilf
end module LBM_Equilibrium_module
