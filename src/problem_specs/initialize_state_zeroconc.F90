!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state_zeroconc.F90
!!!     version:         
!!!     created:         20 April 2011
!!!       on:            16:55:58 MDT
!!!     last modified:   20 April 2011
!!!       at:            17:10:47 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

  subroutine initialize_state_transport(fi, rho, u, walls, dist, phases, options)
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

    ! initialize state
    fi = 0.0
    u = 0.0
    rho = 0.0
    return
  end subroutine initialize_state_transport
