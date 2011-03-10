!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            18:21:06 MST
!!!     last modified:   10 March 2011
!!!       at:            10:27:34 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

  subroutine initialize_state(fi, rho, u, walls, info, constants)
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
         info%gys:info%gye):: fi
    PetscScalar,dimension(1:info%s, &
         info%gxs:info%gxe, &
         info%gys:info%gye):: rho
    PetscScalar,dimension(1:info%s, 1:info%dim, &
         info%gxs:info%gxe, &
         info%gys:info%gye):: u
    PetscScalar,dimension(info%gxs:info%gxe, &
         info%gys:info%gye):: walls

    ! local variables
    PetscErrorCode ierr
    PetscBool flag
    logical,dimension(info%gxs:info%gxe, &
         info%gys:info%gye):: bound
    PetscScalar,dimension(1:info%s):: rho1, rho2         ! left and right fluid densities?
    PetscInt nmax
    PetscBool flushx, flushy

    PetscInt i,j,m ! local values

    ! input data
    nmax = constants%s
    call PetscOptionsGetRealArray(info%options_prefix, '-rho1', rho1, nmax, flag, ierr)
    nmax = constants%s
    call PetscOptionsGetRealArray(info%options_prefix, '-rho2', rho2, nmax, flag, ierr)
    
    flushy = .FALSE.
    flushx = .TRUE.
    call PetscOptionsGetBool(info%options_prefix, '-flush_direction_x', flushx, flag, ierr)
    call PetscOptionsGetBool(info%options_prefix, '-flush_direction_y', flushy, flag, ierr)
    
    ! initialize state
    fi=0.d0
    u=0.d0

    ! flushing experiement 
    if (flushx) then
       bound=.false.
       do i=info%xs,info%xe
          do j=info%ys,info%ye
             if(i.le.10) bound(i,j)=.true.
          enddo
       enddo
    else if (flushy) then
       bound=.false.
       do i=info%xs,info%xe
          do j=info%ys,info%ye
             if(j.le.10) bound(i,j)=.true.
          enddo
       enddo
    end if

    ! set density
    where(bound)
       rho(1,:,:)=rho1(2)
       rho(2,:,:)=rho2(1)
    else where
       rho(1,:,:)=rho1(1)
       rho(2,:,:)=rho2(2)
    end where
    
    where(walls.eq.1)
       rho(1,:,:) = 0
       rho(2,:,:) = 0
    end where
    
    ! set state at equilibrium       
    do i=info%xs,info%xe
       do j=info%ys,info%ye
          do m=1,info%s
             call LBMEquilf(fi(m,:,i,j),rho(m,i,j),u(m,:,i,j), &
                  constants%alf(m), info)
          enddo
       enddo
    enddo
    
    return
  end subroutine initialize_state
