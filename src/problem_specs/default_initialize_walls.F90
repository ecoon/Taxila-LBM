!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_walls.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            17:26:22 MST
!!!     last modified:   08 August 2011
!!!       at:            15:41:04 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petscdmdef.h"

  subroutine initialize_walls(walls, info, options)
    use petsc
    use LBM_Options_module
    use LBM_Info_module
    implicit none

#include "lbm_definitions.h"
!   input variables
    type(info_type) info
    type(options_type) options
    PetscScalar,dimension(info%rgxyzl):: walls

    walls=0
    return
  end subroutine initialize_walls
