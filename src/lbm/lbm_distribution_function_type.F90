!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_distribution_function_type.F90
!!!     version:         
!!!     created:         12 April 2011
!!!       on:            11:15:29 MDT
!!!     last modified:   12 April 2011
!!!       at:            11:20:32 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petscdmdef.h"

module LBM_Distribution_Function_type_module
  use petsc
  use LBM_Info_module
  use LBM_Discretization_Type_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: distribution_type
     MPI_Comm comm
     PetscInt s
     PetscInt b
     type(info_type),pointer:: info
     type(discretization_type),pointer:: disc
     DM,pointer :: da_fi, da_rho

     Vec fi
     Vec fi_g
     PetscScalar,pointer:: fi_a(:)
     Vec fi_g_old

     Vec rho
     Vec rho_g
     PetscScalar,pointer:: rho_a(:)
     Vec rho_g_old

     PetscScalar,pointer,dimension(:,:,:):: flux
     PetscBool flux_required
     character(len=MAXWORDLENGTH) name
     PetscBool track_old_rho
     PetscBool track_old_fi
  end type distribution_type
end module LBM_Distribution_Function_type_module
