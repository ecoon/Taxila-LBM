!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            18:21:06 MST
!!!     last modified:   31 January 2011
!!!       at:            14:58:46 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

  subroutine initialize_state(fi, rho, ux, uy, uz, walls, info, constants)
    use petsc
    use LBM_Info_module
    use LBM_Constants_module
    use LBM_Equilibrium_module
    implicit none

    ! input variables
    type(info_type) info
    type(constants_type) constants
    PetscScalar,dimension(1:info%s,0:info%b, &
         info%gxs:info%gxe, &
         info%gys:info%gye, &
         info%gzs:info%gze):: fi
    PetscScalar,dimension(1:info%s, &
         info%gxs:info%gxe, &
         info%gys:info%gye, &
         info%gzs:info%gze):: rho
    PetscScalar,dimension(1:info%s, &
         info%gxs:info%gxe, &
         info%gys:info%gye, &
         info%gzs:info%gze):: ux,uy,uz
    PetscScalar,dimension(info%gxs:info%gxe, &
         info%gys:info%gye, &
         info%gzs:info%gze):: walls

    ! local variables
    PetscErrorCode ierr
    logical,dimension(info%gxs:info%gxe, &
         info%gys:info%gye, &
         info%gzs:info%gze):: bound

    PetscInt i,j,k,m ! local values
    
    ! initialize state
    fi=0.0
    ux=0.0
    uy=0.0
    uz=0.0

    ! flushing experiement 
    bound=.false.
    do i=info%xs,info%xe
       do j=info%ys,info%ye
          do k=info%zs,info%ze
             if(k.le.10) bound(i,j,k)=.true.
          enddo
       enddo
    enddo
    
    ! set density
    where(bound)
       rho(1,:,:,:)=constants%rho1(2)
       rho(2,:,:,:)=constants%rho2(1)
    else where
       rho(1,:,:,:)=constants%rho1(1)
       rho(2,:,:,:)=constants%rho2(2)
    end where
    
    where(walls.eq.1)
       rho(1,:,:,:) = 0
       rho(2,:,:,:) = 0
    end where
    
    ! set state at equilibrium       
    do i=info%xs,info%xe
       do j=info%ys,info%ye
          do k=info%zs,info%ze
             do m=1,info%s
                call LBMEquilf(fi(m,:,i,j,k),rho(m,i,j,k),ux(m,i,j,k), &
                     uy(m,i,j,k), uz(m,i,j,k), constants%alf(m), info%b)
             enddo
          enddo
       enddo
    enddo
    
    return
  end subroutine initialize_state
