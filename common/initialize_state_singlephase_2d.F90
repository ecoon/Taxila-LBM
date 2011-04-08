!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            18:21:06 MST
!!!     last modified:   08 April 2011
!!!       at:            11:17:40 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

  subroutine initialize_state(fi, rho, u, walls, dist, phases, options)
    use petsc
    use LBM_Distribution_Function_module
    use LBM_Phase_module
    use LBM_Options_module
    use LBM_Equilibrium_module
    implicit none

    ! input variables
    type(distribution_type) dist
    type(phase_type) phases(dist%s)
    type(options_type) options
    PetscScalar,dimension(dist%s,0:dist%b, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%s, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: rho
    PetscScalar,dimension(dist%s, 1:dist%info%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls

    ! local variables
    PetscErrorCode ierr
    PetscBool flag
    logical,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: bound
    PetscScalar,dimension(dist%s):: rho1, rho2         ! left and right fluid densities?
    PetscInt nmax
    PetscBool flushx, flushy
    PetscScalar,dimension(dist%s, 1:dist%info%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: forces

    PetscInt i,j,m ! local values
    PetscBool help

    ! input data
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    
    ! initialize state
    fi=0.d0
    u=0.d0

    ! set density
    rho(1,:,:)=1.d0
    
    where(walls.eq.1)
       rho(1,:,:) = 0
    end where
    
    ! set state at equilibrium       
    forces = 0.d0
    call LBMEquilfFlow(fi, rho, u, forces, walls, phases, dist)
    
    return
  end subroutine initialize_state
