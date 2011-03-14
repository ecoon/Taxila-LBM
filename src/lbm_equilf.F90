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
    PetscScalar,dimension(0:info%flow_disc%b)::feq
    PetscScalar rho
    PetscScalar,dimension(1:info%ndims):: ue
    PetscErrorCode ierr
    
    PetscScalar usqr,tmp
    PetscInt i
    
    usqr = SUM(ue*ue, 1)
    do i=0,info%flow_disc%b
       tmp = SUM(info%flow_disc%ci(i,:)*ue, 1)
       feq(i)= info%flow_disc%weights(i)*rho* &
            (1 + tmp/constants%c_s2 + tmp*tmp/(2.d0*constants%c_s2*constants%c_s2) &
            - usqr/(2.d0*constants%c_s2))
    end do
  end subroutine LBMEquilf
end module LBM_Equilibrium_module
