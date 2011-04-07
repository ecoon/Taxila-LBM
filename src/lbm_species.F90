!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_specie.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            15:34:44 MDT
!!!     last modified:   05 April 2011
!!!       at:            08:46:25 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscbagdef.h"

module LBM_Component_module
  use petsc
  use LBM_Component_Bag_Data_type_module
  use LBM_Relaxation_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: specie_type
     MPI_Comm comm
     ! sizes and identifiers (set pre-bag)
     PetscInt id

     ! bagged parameters

     ! dependent parameters
     type(relaxation_type),pointer:: relax

     ! bag 
     character(len=MAXWORDLENGTH):: name
     type(specie_bag_data_type),pointer:: data
     PetscBag bag
  end type specie_type

  interface PetscBagGetData
     subroutine PetscBagGetData(bag, data, ierr)
       use LBM_Component_Bag_Data_type_module
       PetscBag bag
       type(specie_bag_data_type),pointer :: data
       PetscErrorCode ierr
     end subroutine PetscBagGetData
  end interface

  interface ComponentCreate
     procedure ComponentCreateOne, ComponentCreateN
  end interface

  public :: ComponentCreate, &
       ComponentDestroy, &
       ComponentSetSizes, &
       ComponentSetName, &
       ComponentSetFromOptions

contains
  function ComponentCreateOne(comm) result(specie)
    MPI_Comm comm
    type(specie_type),pointer :: specie
    allocate(specie)
    specie%comm = comm
    call ComponentInitialize(specie)
    specie%relax => RelaxationCreate(comm)
  end function ComponentCreateOne

  function ComponentCreateN(comm, n) result(specie)
    MPI_Comm comm
    PetscInt n
    type(specie_type),pointer,dimension(:):: specie
    PetscInt lcv
    allocate(specie(n))

    do lcv=1,n
       specie(lcv)%comm = comm
       call ComponentInitialize(specie(lcv))
       specie(lcv)%relax => RelaxationCreate(comm)
    end do
  end function ComponentCreateN

  subroutine ComponentInitialize(specie)
    type(specie_type) specie
    specie%s = -1
    specie%b = -1
    specie%id = 0
    specie%name = ''
    nullify(specie%data)
    specie%bag = 0
  end subroutine ComponentInitialize
  
  subroutine ComponentSetSizes(specie, s, b)
    type(specie_type) :: specie
    PetscInt s,b
    specie%s = s
    call RelaxationSetSizes(specie%relax, s, b)
  end subroutine ComponentSetSizes

  subroutine ComponentSetName(specie, name)
    type(specie_type) specie
    character(len=MAXWORDLENGTH):: name
    
    specie%name = name
    call RelaxationSetName(specie%relax, name)
  end subroutine ComponentSetName

  subroutine ComponentSetID(specie, id)
    type(specie_type) specie
    PetscInt id
    specie%id = id
  end subroutine ComponentSetID

  subroutine ComponentSetFromOptions(specie, options, ierr)
    use LBM_Options_module
    type(specie_type) specie
    type(options_type) options
    PetscErrorCode ierr

    PetscInt sizeofint, sizeofscalar, sizeofbool, sizeofdata
    PetscInt lcv
    character(len=MAXWORDLENGTH):: paramname

    ! set up the data

    ! create the bag
    call PetscDataTypeGetSize(PETSC_SCALAR, sizeofscalar, ierr)
    sizeofdata = 0 !???
    call PetscBagCreate(specie%comm, sizeofdata, specie%bag, ierr)
!    call PetscBagSetName(specie%bag, TRIM(options%my_prefix)//specie%name, ierr)

    ! register data

    call RelaxationSetFromOptions(specie%relax, options, ierr)
  end subroutine ComponentSetFromOptions

  subroutine ComponentDestroy(specie, ierr)
    type(specie_type) specie
    PetscErrorCode ierr
    if (specie%bag /= 0) call PetscBagDestroy(specie%bag, ierr)
  end subroutine ComponentDestroy
end module LBM_Component_module
