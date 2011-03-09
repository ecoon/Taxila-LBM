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
    
  subroutine LBMEquilfD3Q19(feq, rho, ue, alf, info)
    use LBM_Discretization_D3Q19_module
    
    type(info_type) info
    PetscScalar,dimension(0:info%b)::feq
    PetscScalar rho, alf
    PetscScalar,dimension(1:info%dim):: ue
    
    PetscScalar usqr,tmp
    PetscInt i
    
    usqr = SUM(ue*ue, 1)
    feq(0)=(1./3.-0.5*usqr)*rho
    
    do i=1,info%dim*2
       tmp = SUM(ci(i,:)*ue, 1)
       feq(i)=(1./18.+tmp/6.+tmp*tmp/4.-usqr/12.)*rho
    end do
      
    do i=info%dim*2+1,info%b
       tmp = SUM(ci(i,:)*ue, 1)
       feq(i)=(1./36.+tmp/12.+tmp*tmp/8.-usqr/24.)*rho
    end do
    
    return
  end subroutine LBMEquilfD3Q19

  subroutine LBMEquilfD2Q9(feq, rho, ue, alf, info)
    use LBM_Discretization_D2Q9_module
    
    type(info_type) info
    PetscScalar,dimension(0:info%b)::feq
    PetscScalar rho, alf
    PetscScalar,dimension(1:info%dim):: ue
    
    PetscScalar usqr,tmp
    PetscInt i
    
    usqr = SUM(ue*ue, 1)
    feq(0)= (alf-2./3.*usqr)*rho

    do i=1,info%dim*2
       tmp = SUM(ci(i,:)*ue, 1)
       feq(i)=((1.-alf)/5.+tmp/3.+tmp*tmp/2.-usqr/6.)*rho
    end do
    
    do i=info%dim*2+1,info%b
       tmp = SUM(ci(i,:)*ue, 1)
       feq(i)=((1.-alf)/20.+tmp/12.+tmp*tmp/8.-usqr/24.)*rho
    end do
    return
  end subroutine LBMEquilfD2Q9
end module LBM_Equilibrium_module
