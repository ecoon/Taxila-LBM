!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_bc_constant.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            17:30:22 MST
!!!     last modified:   05 August 2011
!!!       at:            11:35:32 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "petsc/finclude/petscsysdef.h"

! does nothing -- homogenous BCs are set on the command line
  subroutine initialize_bcs(bc_flags, xm_bcvals, xp_bcvals, ym_bcvals, &
       yp_bcvals, zm_bcvals, zp_bcvals, bc_dim, dist, options)
    use petsc
    use LBM_Error_module
    use LBM_Distribution_Function_type_module
    use LBM_Options_module
    implicit none

    PetscInt, dimension(6):: bc_flags ! enum for boundary conditions
    type(distribution_type) dist
    type(options_type) options
    PetscScalar :: xm_bcvals(*), xp_bcvals(*)
    PetscScalar :: ym_bcvals(*), yp_bcvals(*)
    PetscScalar :: zm_bcvals(*), zp_bcvals(*)
    PetscInt bc_dim

    return
  end subroutine initialize_bcs
