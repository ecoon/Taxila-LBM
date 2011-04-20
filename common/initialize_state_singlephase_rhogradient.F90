!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            18:21:06 MST
!!!     last modified:   20 April 2011
!!!       at:            13:45:19 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

! initializes the state to a constant gradient in rho in the
! x-direction, given by the left and right rho values given as BCs
  subroutine initialize_state(fi, rho, u, walls, dist, phases, options)
    use petsc
    use LBM_Distribution_Function_type_module
    use LBM_Phase_module
    use LBM_Options_module
    use LBM_Discretization_module
    implicit none

    ! input variables
    type(distribution_type) dist
    type(phase_type) phases(dist%s)
    type(options_type) options
    PetscScalar,dimension(dist%s,0:dist%b, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%s, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: rho
    PetscScalar,dimension(dist%s, 1:dist%info%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: u
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: walls

    ! local variables
    PetscErrorCode ierr
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: nowalls

    PetscScalar :: xp3_ave_p, xm3_ave_p, t
    PetscBool flag
    PetscInt i

    ! initialize state
    fi = 0.d0
    u = 0.d0

    xp3_ave_p = 0.d0
    xm3_ave_p = 0.d0

    call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_xp_phase1', xp3_ave_p, &
         flag, ierr)
    if (.not.flag) SETERRQ(1, 1, 'invalid boundary pressure for xp_phase1', ierr)
    call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_xm_phase1', xm3_ave_p, &
         flag, ierr)
    if (.not.flag) SETERRQ(1, 1, 'invalid boundary value for xm_phase1', ierr)

    do i=dist%info%xs,dist%info%xe
       t = dble(i-1)/dble(dist%info%NX-1)
       rho(1,i,:,:) = (1-t)*xm3_ave_p + t*xp3_ave_p
    end do    
    
    ! set state at equilibrium       
    nowalls = 0.d0
    call DiscretizationEquilf(dist%disc, rho, u, nowalls, fi, phases(1)%relax, dist)    
    return
  end subroutine initialize_state
