!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_relaxation.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            15:15:25 MDT
!!!     last modified:   30 March 2011
!!!       at:            17:02:48 MDT
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
     PetscInt id
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
       RelaxationSetID, &
       RelaxationSetMode, &
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
    relax%name = ''
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

  subroutine RelaxationSetID(relax, id)
    type(relaxation_type) relax
    PetscInt id
    relax%id = id
  end subroutine RelaxationSetID

  subroutine RelaxationSetFromOptions(relax, options, ierr)
    use LBM_Options_module
    type(relaxation_type) relax
    type(options_type) options
    PetscErrorCode ierr

    PetscSizeT sizeofscalar, sizeofdata
    PetscInt lcv
    character(len=MAXWORDLENGTH):: paramname, paramname2

    ! create the bag
    call PetscDataTypeGetSize(PETSC_SCALAR, sizeofscalar, ierr)
    sizeofdata = (relax%b+2)*sizeofscalar 
    call PetscBagCreate(relax%comm, sizeofdata, relax%bag, ierr)
    call PetscBagGetData(relax%bag, relax%data, ierr)

    write(paramname, '(I1)') relax%id
    call PetscBagRegisterScalar(relax%bag, relax%data%tau, 1.d0, &
         'tau'//paramname, 'relaxation time', ierr)
    relax%tau => relax%data%tau

    if (relax%mode == RELAXATION_MODE_MRT) then
       do lcv=0,relax%b
          write(paramname2, '(I0.2)') lcv
          call PetscBagRegisterScalar(relax%bag, relax%data%tau_mrt(lcv), 0.d0, &
               'mrt'//trim(paramname)//'_s'//paramname2, 'MRT relaxation moment coefficient', ierr)
       end do
    end if
    call PetscBagSetName(relax%bag, TRIM(options%my_prefix)//relax%name, "", ierr)
    relax%tau_mrt => relax%data%tau_mrt
  end subroutine RelaxationSetFromOptions
end module LBM_Relaxation_module

    
