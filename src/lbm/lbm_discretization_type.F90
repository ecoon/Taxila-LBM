!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_discretization_type.F90
!!!     version:         
!!!     created:         15 March 2011
!!!       on:            17:04:20 MDT
!!!     last modified:   14 September 2011
!!!       at:            12:06:06 PDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

module LBM_Discretization_Type_module
  implicit none

  private

  type, public:: discretization_type
     MPI_Comm comm
     PetscInt name
     PetscInt ndims
     PetscInt b
     PetscInt deriv_order

     PetscScalar :: c_0 ! grid factor for fluid-fluid terms
     PetscInt,pointer:: ci(:,:)
     PetscScalar,pointer:: weights(:)
     PetscInt,pointer:: opposites(:)
     PetscInt,pointer:: reflect_x(:)
     PetscInt,pointer:: reflect_y(:)
     PetscInt,pointer:: reflect_z(:)
     PetscScalar,pointer:: mt_mrt(:,:)
     PetscScalar,pointer:: mmt_mrt(:)
     PetscScalar,pointer:: ffw(:)
  end type discretization_type

end module LBM_Discretization_Type_module
