!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_discretization_type.F90
!!!     version:         
!!!     created:         15 March 2011
!!!       on:            17:04:20 MDT
!!!     last modified:   15 March 2011
!!!       at:            17:04:47 MDT
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

     PetscInt,pointer:: ci(:,:)
     PetscScalar,pointer:: weights(:)
     PetscScalar,pointer:: m_mrt(:,:)
     PetscScalar,pointer:: mt_mrt(:,:)
     PetscScalar,pointer:: mmt_mrt(:)
     !PetscScalar,pointer:: tau_mrt(:,:)

  end type discretization_type

end module LBM_Discretization_Type_module
