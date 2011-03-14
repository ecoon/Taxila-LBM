#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module LBM_Equilibrium_module
  use LBM_Info_module
  implicit none
  private
#include "lbm_definitions.h"
  public :: LBMEquilf

contains
  subroutine LBMEquilf(feq, rho, ue, alf, info)
    use petsc
    type(info_type) info
    PetscScalar,dimension(0:info%b)::feq
    PetscScalar rho, alf
    PetscScalar,dimension(1:info%dim):: ue
    PetscErrorCode ierr
    
    select case(info%discretization)
    case(D3Q19_DISCRETIZATION)
       call LBMEquilfD3Q19(feq, rho, ue, alf, info)
    case(D2Q9_DISCRETIZATION)
       call LBMEquilfD2Q9(feq, rho, ue, alf, info)
    case DEFAULT
       SETERRQ(1, 1, 'invalid discretization in LBM', ierr)
    end select
  end subroutine LBMEquilf
    
  subroutine LBMEquilfD3Q19(feq, rho, ue, constants, info)
    use LBM_Discretization_D3Q19_module
    use LBM_Constants_module
    
    type(info_type) info
    type(constants_type) constants
    PetscScalar,dimension(0:info%b)::feq
    PetscScalar rho
    PetscScalar,dimension(1:info%dim):: ue
    
    PetscScalar usqr,tmp
    PetscInt i
    
    usqr = SUM(ue*ue, 1)
    
    do i=0,info%b
       tmp = SUM(ci(i,:)*ue, 1)
       feq(i)= weights(i)*rho* &
            (1 + tmp/constants%c_s2 + tmp*tmp/(2.d0*constants%c_s2*constants&c_s2) - usqr/(2.d0*constants%c_s2))
    end do
    
    return
  end subroutine LBMEquilfD3Q19

  subroutine LBMEquilfD2Q9(feq, rho, ue, constants, info)
    use LBM_Discretization_D2Q9_module
    
    type(info_type) info
    type(constants_type) constants
    PetscScalar,dimension(0:info%b)::feq
    PetscScalar rho
    PetscScalar,dimension(1:info%dim):: ue
    
    PetscScalar usqr,tmp
    PetscInt i
    
    usqr = SUM(ue*ue, 1)
    
    do i=0,info%b
       tmp = SUM(ci(i,:)*ue, 1)
       feq(i)= weights(i)*rho* &
            (1 + tmp/constants%c_s2 + tmp*tmp/(2.d0*constants%c_s2*constants&c_s2) - usqr/(2.d0*constants%c_s2))
    end do
    
    return
  end subroutine LBMEquilfD2Q9
end module LBM_Equilibrium_module
