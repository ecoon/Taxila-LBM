!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_bag_types.F90
!!!     version:         
!!!     created:         17 March 2011
!!!       on:            16:59:38 MDT
!!!     last modified:   22 August 2011
!!!       at:            15:55:46 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

! collection of all the stupid bag types that need to get defined ahead of time
! due to fortran stupidity
module LBM_Component_Bag_Data_type_module
  implicit none
  private
#include "lbm_definitions.h"
  ! physical parameters for a component
  type, public :: component_bag_data_type
     PetscScalar mm ! molecular mass
     PetscScalar viscosity
     PetscScalar density
     PetscScalar,dimension(NMAX_COMPONENTS) :: gf ! component-component force coefs
  end type component_bag_data_type
end module LBM_Component_Bag_Data_type_module

module LBM_Mineral_Bag_Data_type_module
  implicit none
  private
#include "lbm_definitions.h"
  ! physical parameters for a mineral
  type, public :: mineral_bag_data_type
     PetscScalar,dimension(NMAX_COMPONENTS) :: gw ! mineral-mineral force coefs
  end type mineral_bag_data_type
end module LBM_Mineral_Bag_Data_type_module

module LBM_Specie_Bag_Data_type_module
  implicit none

  ! physical parameters for a chemical specie
  type, public :: specie_bag_data_type
     PetscInt component 
  end type specie_bag_data_type
end module LBM_Specie_Bag_Data_type_module

module LBM_Info_Bag_Data_type_module
  implicit none
  
  type, public :: info_bag_data_type
     PetscInt NX,NY,NZ
     PetscInt stencil_size
     PetscInt stencil_size_rho
     PetscInt stencil_type
     PetscBool,dimension(3) :: periodic
     PetscScalar,dimension(3,2) :: corners
  end type info_bag_data_type
end module LBM_Info_Bag_Data_type_module

module LBM_Relaxation_Bag_Data_type_module
  implicit none
  private
#include "lbm_definitions.h"
  
  type, public :: relaxation_bag_data_type
     PetscScalar tau ! relaxation time
     PetscScalar s1  ! MRT relaxation time
     PetscScalar s2  ! MRT relaxation time
     PetscScalar s3  ! MRT relaxation time
     PetscScalar s4  ! MRT relaxation time (only for 3D)
     PetscScalar s5  ! MRT relaxation time (only for 3D)
  end type relaxation_bag_data_type
end module LBM_Relaxation_Bag_Data_type_module
