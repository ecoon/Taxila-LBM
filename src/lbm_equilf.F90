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
    PetscScalar,dimension(1:info%s):: jx
    PetscScalar,dimension(1:info%s):: jy
    PetscScalar,dimension(1:info%s):: jz
    PetscScalar,dimension(1:info%s):: jsqr

    if(info%MRT) then
     if(info%ndims.eq.3) then
       

       
!! It appears that ue is the equilibrium velocity and I need u = sum(fi,ci)
!!$       jx(:) = ue(:,1)
!!$       jy(:) = ue(:,2)
!!$       jz(:) = ue(:,3)
       jx(:) = rho(:)*ue(:,1)
       jy(:) = rho(:)*ue(:,2)
       jz(:) = rho(:)*ue(:,3)
       jsqr(:) = jx(:)**2 + jy(:)**2 + jz(:)**2

       
       ! Assumes rho_0 = 1.0 from d'Humieres et al. (2002) and
       ! Pan et al. (2006).  I am following Pan et al. (2006) and
       ! I have chosen to write out meq explicitly because there 
       ! are free parameters in the derivation by d'Humieres. 
       feq(:, 0) = rho(:)
       feq(:, 1) = -11.d0*rho(:)+19.d0*jsqr(:)
       !feq(:, 2) = 0.d0*rho(:)-475.d0/63.d0*jsqr(:)
       feq(:, 2) = 3.d0*rho(:)-11.d0/2.d0*jsqr(:)
       feq(:, 3) = jx(:)
       feq(:, 4) = -2.d0/3.d0*jx(:)
       feq(:, 5) = jy(:)
       feq(:, 6) = -2.d0/3.d0*jy(:)
       feq(:, 7) = jz(:)
       feq(:, 8) = -2.d0/3.d0*jz(:)
       feq(:, 9) = 2.d0*jx(:)**2-(jy(:)**2+jz(:)**2)
       !feq(:,10) = 0.d0
       feq(:,10) = -1.d0/2.d0*feq(:,9)
       feq(:,11) = jy(:)**2-jz(:)**2
       !feq(:,12) = 0.d0 
       feq(:,12) = -1.d0/2.d0*feq(:,11)
       feq(:,13) = jx(:)*jy(:)
       feq(:,14) = jy(:)*jz(:)
       feq(:,15) = jx(:)*jz(:)
       feq(:,16) = 0.d0
       feq(:,17) = 0.d0
       feq(:,18) = 0.d0  

     endif
     ! Need to work on the 2D MRT

    else    

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
    endif

  end subroutine LBMEquilf

end module LBM_Equilibrium_module
