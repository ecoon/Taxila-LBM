!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_walls_from_petsc.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            17:15:55 MST
!!!     last modified:   14 January 2011
!!!       at:            17:17:19 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscviewerdef.h"

  subroutine initialize_walls(walls, filename, info)
    use petsc
    use Info_module
    implicit none

#include "lbm_definitions.h"
! input variables
    type(info_type) info
    PetscScalar,dimension(info%gxs:info%gxe, &
         info%gys:info%gye, &
         info%gzs:info%gze):: walls
    character(len=MAXSTRINGLENGTH) filename
    
    ! do nothing
    return
  end subroutine initialize_walls
