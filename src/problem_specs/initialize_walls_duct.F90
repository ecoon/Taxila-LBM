!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_walls.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            17:26:22 MST
!!!     last modified:   14 September 2011
!!!       at:            12:42:22 PDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

  subroutine initialize_walls(walls, filename, info)
    use petsc
    use LBM_Info_module
    implicit none

#include "lbm_definitions.h"
!   input variables
    type(info_type) info
    PetscScalar,dimension(info%rgxs:info%rgxe, &
         info%rgys:info%rgye, &
         info%rgzs:info%rgze):: walls
    character(len=MAXSTRINGLENGTH) filename

    walls=0
    
    if (info%ys.eq.1) walls(:,1,:)=1
    if (info%ye.eq.info%NY) walls(:,info%NY,:)=1
    
    if (info%zs.eq.1) walls(:,:,1)=1
    if (info%ze.eq.info%NZ) walls(:,:,info%NZ)=1
    return
  end subroutine initialize_walls
