!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_walls_from_files.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            17:15:36 MST
!!!     last modified:   14 January 2011
!!!       at:            17:16:10 MST
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
! input variables
    type(info_type) info
    PetscScalar,dimension(info%gxs:info%gxe, &
         info%gys:info%gye, &
         info%gzs:info%gze):: walls 
    character(len=MAXSTRINGLENGTH):: filename
    
    integer charlen
    integer i1d, i2d, i3d
    character(len=MAXSTRINGLENGTH) flnm1
    logical,dimension(info%xs:info%xe, info%ys:info%ye, info%zs:info%ze):: tmp
    
    walls=0
    tmp = .FALSE.
    
    i1d=info%id/100
    i2d=mod(info%id,100)/10
    i3d=mod(info%id,10)   
    
    charlen = LEN_TRIM(filename)
    flnm1=filename(1:charlen)//'.'//char(i1d+48)//char(i2d+48)//char(i3d+48)

    open(11,file=flnm1,status='unknown',form='unformatted')
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
  
