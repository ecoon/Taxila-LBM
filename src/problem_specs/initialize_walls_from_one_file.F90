!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_walls_from_one_file.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            17:15:47 MST
!!!     last modified:   14 January 2011
!!!       at:            17:16:18 MST
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
    PetscScalar,dimension(info%gxs:info%gxe, &
         info%gys:info%gye, &
         info%gzs:info%gze):: walls
    character(len=MAXSTRINGLENGTH) filename

    logical,dimension(1:info%NX, 1:info%NY, 1:info%NZ):: tmp
    
    integer i,j,k
    
    walls=0
    
    open(11,file=filename,status='unknown',form='unformatted')
    read(11) tmp
    close(11)
    
    do i=info%xs,info%xe
       do j=info%ys,info%ye
          do k=info%zs,info%ze
             if (tmp(i,j,k)) walls(i,j,k) = 1
          enddo
       enddo
    enddo
    return
  end subroutine initialize_walls
  
