!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_relaxation.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            15:15:25 MDT
!!!     last modified:   28 March 2011
!!!       at:            16:55:08 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscbagdef.h"

module LBM_Relaxation_module
  use petsc
  use LBM_Relaxation_Bag_Data_type_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: relaxation_type
     MPI_Comm comm
     PetscInt s
     PetscInt b
     PetscInt mode
     PetscScalar,pointer :: tau ! relaxation time
     PetscScalar,pointer,dimension(:) :: tau_mrt ! components of S vector for mrt

     ! bag 
     character(len=MAXWORDLENGTH):: name
     type(relaxation_bag_data_type),pointer:: data
     PetscBag bag
  end type relaxation_type
     
  interface PetscBagGetData
     subroutine PetscBagGetData(bag, data, ierr)
       use LBM_Relaxation_Bag_Data_type_module
       PetscBag bag
       type(relaxation_bag_data_type),pointer :: data
       PetscErrorCode ierr
     end subroutine PetscBagGetData
  end interface

  public :: RelaxationCreate, &
       RelaxationDestroy, &
       RelaxationSetSizes, &
       RelaxationSetName, &
       RelaxationSetFromOptions

contains
  function RelaxationCreate(comm) result(relax)
    MPI_Comm comm
    type(relaxation_type),pointer:: relax
    allocate(relax)
    relax%comm = comm
    relax%s = -1
    relax%b = -1
    nullify(relax%data)
    relax%bag = 0
  end function RelaxationCreate

  subroutine RelaxationDestroy(relax, ierr)
    type(relaxation_type) relax
    PetscErrorCode ierr
    if (relax%bag /= 0) call PetscBagDestroy(relax%bag, ierr)
  end subroutine RelaxationDestroy

  subroutine RelaxationSetSizes(relax, s, b)
    type(relaxation_type) relax
    PetscInt:: s, b
    relax%s = s
    relax%b = b
  end subroutine RelaxationSetSizes

  subroutine RelaxationSetMode(relax, mode)
    type(relaxation_type) relax
    PetscInt mode
    relax%mode = mode
  end subroutine RelaxationSetMode

  subroutine RelaxationSetName(relax, name)
    type(relaxation_type) relax
    character(len=MAXWORDLENGTH):: name
    
    relax%name = name
  end subroutine RelaxationSetName

  subroutine RelaxationSetFromOptions(relax, options, ierr)
    use LBM_Options_module
    type(relaxation_type) relax
    type(options_type) options
    PetscErrorCode ierr

    PetscInt sizeofscalar, sizeofdata
    PetscInt lcv
    character(len=MAXWORDLENGTH):: paramname

    ! create the bag
    call PetscDataTypeGetSize(PETSC_SCALAR, sizeofscalar, ierr)
    sizeofdata = (relax%b+2)*sizeofscalar 
    call PetscBagCreate(relax%comm, sizeofdata, relax%bag, ierr)
    call PetscBagSetName(relax%bag, TRIM(options%my_prefix)//relax%name, ierr)


    ! register data
    call PetscBagGetData(relax%bag, relax%data, ierr)

    call PetscBagRegisterScalar(relax%bag, relax%data%tau, 1.d0, &
         'tau', 'relaxation time', ierr)
    relax%tau => relax%data%tau

    if (relax%mode == RELAXATION_MODE_MRT) then
       allocate(relax%data%tau_mrt(0:relax%b))
       do lcv=0,relax%b
          write(paramname, '(I0.2)') lcv
          call PetscBagRegisterScalar(relax%bag, relax%data%tau_mrt(lcv), 0.d0, &
               'tau_'//paramname, 'MRT relaxation moment coefficient', ierr)
       end do
    else
       relax%data%tau_mrt = 0.d0
    end if
    relax%tau_mrt => relax%data%tau_mrt
  end subroutine RelaxationSetFromOptions
end module LBM_Relaxation_module

    
