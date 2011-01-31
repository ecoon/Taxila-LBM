#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module LBM_Equilibrium_module
  implicit none
  private
  public :: LBMEquilf

contains

  
  subroutine LBMEquilf(feq,rho,uxe,uye,uze,alf,b)
    use c_constants
    
    integer b
    PetscScalar,dimension(0:b)::feq
    PetscScalar rho, uxe, uye, uze, alf
    
    PetscScalar usqr,tmp
    integer i
    
    usqr = uxe*uxe + uye*uye + uze*uze
    
    feq(0)=(1./3.-0.5*usqr)*rho
    
    do i=1,6
       tmp = (cix(i)*uxe + ciy(i)*uye +ciz(i)*uze )
       feq(i)=(1./18.+tmp/6.+tmp*tmp/4.-usqr/12.)*rho
    end do
      
    do i=7,18
       tmp = (cix(i)*uxe + ciy(i)*uye +ciz(i)*uze)
       feq(i)=(1./36.+tmp/12.+tmp*tmp/8.-usqr/24.)*rho
    end do
    
    return
  end subroutine LBMEquilf
end module LBM_Equilibrium_module
