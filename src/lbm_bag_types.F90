!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_bag_types.F90
!!!     version:         
!!!     created:         17 March 2011
!!!       on:            16:59:38 MDT
!!!     last modified:   17 March 2011
!!!       at:            17:34:46 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

! collection of all the stupid bag types that need to get defined ahead of time
! due to fortran stupidity
module LBM_Phase_Bag_Data_type_module
  implicit none

  ! physical parameters for either a phase phase or phase component
  type, public :: phase_bag_data_type
     PetscScalar mm ! molecular mass
     PetscScalar tau ! relaxation time
     PetscScalar,pointer,dimension(:) :: tau_mrt ! components of S vector for mrt
     PetscScalar gw ! solid affinity? for phase-wall interaction forces
     PetscScalar,pointer,dimension(:) :: gf ! phase-phase force coefs
     PetscBool mrt ! do MRT?
  end type phase_bag_data_type
end module LBM_Phase_Bag_Data_type_module

module LBM_Info_Bag_Data_type_module
  implicit none
  
  type, public :: info_bag_data_type
     PetscInt NX,NY,NZ
     PetscBool,pointer,dimension(:) :: periodic
     PetscScalar,pointer,dimension(:,:) :: corners
  end type info_bag_data_type
end module LBM_Info_Bag_Data_type_module
