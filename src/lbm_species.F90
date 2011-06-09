!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_specie.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            15:34:44 MDT
!!!     last modified:   05 May 2011
!!!       at:            10:29:53 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscbagdef.h"

module LBM_Specie_module
  use petsc
  use LBM_Specie_Bag_Data_type_module
  use LBM_Relaxation_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: specie_type
     MPI_Comm comm
     ! sizes and identifiers (set pre-bag)
     PetscInt s
     PetscInt id

     ! bagged parameters
     PetscInt,pointer:: phase

     ! dependent parameters
     type(relaxation_type),pointer:: relax

     ! bag 
     character(len=MAXWORDLENGTH):: name
     type(specie_bag_data_type),pointer:: data
     PetscBag bag
  end type specie_type

  interface PetscBagGetData
     subroutine PetscBagGetData(bag, data, ierr)
       use LBM_Specie_Bag_Data_type_module
       PetscBag bag
       type(specie_bag_data_type),pointer :: data
       PetscErrorCode ierr
     end subroutine PetscBagGetData
  end interface

  interface SpecieCreate
     module procedure SpecieCreateOne
     module procedure SpecieCreateN
  end interface

  public :: SpecieCreate, &
       SpecieDestroy, &
       SpecieSetSizes, &
       SpecieSetID, &
       SpecieSetName, &
       SpecieSetFromOptions

contains
  function SpecieCreateOne(comm) result(specie)
    MPI_Comm comm
    type(specie_type),pointer :: specie
    character(len=MAXWORDLENGTH):: name
    allocate(specie)
    specie%comm = comm
    call SpecieInitialize(specie)
    specie%relax => RelaxationCreate(comm)
    name = 'specie1'
    call SpecieSetName(specie, name)
  end function SpecieCreateOne
 
  function SpecieCreateN(comm, n) result(species)
    MPI_Comm comm
    PetscInt n
    type(specie_type),pointer,dimension(:):: species
    type(specie_type),pointer:: aspecie
    PetscInt lcv
    character(len=MAXWORDLENGTH):: name
    allocate(species(1:n))
    name = ''

    do lcv=1,n
       aspecie => species(lcv)
       aspecie%comm = comm
       call SpecieInitialize(aspecie)
       aspecie%relax => RelaxationCreate(comm)
       call SpecieSetID(aspecie, lcv)
       name = 'specie'//char(lcv+48)
       call SpecieSetName(aspecie, name)
    end do
  end function SpecieCreateN

  subroutine SpecieInitialize(specie)
    type(specie_type) specie
    specie%s = -1
    specie%id = 0
    specie%name = ''
    nullify(specie%data)
    nullify(specie%relax)
    specie%bag = 0
  end subroutine SpecieInitialize
  
  subroutine SpecieSetSizes(specie, s, b)
    type(specie_type) :: specie
    PetscInt s,b
    specie%s = s
    call RelaxationSetSizes(specie%relax, s, b)
  end subroutine SpecieSetSizes

  subroutine SpecieSetName(specie, name)
    type(specie_type) specie
    character(len=MAXWORDLENGTH):: name
    
    specie%name = name
    call RelaxationSetName(specie%relax, name)
  end subroutine SpecieSetName

  subroutine SpecieSetID(specie, id)
    type(specie_type) specie
    PetscInt id
    specie%id = id
    call RelaxationSetID(specie%relax, id)
  end subroutine SpecieSetID

  subroutine SpecieSetFromOptions(specie, options, ierr)
    use LBM_Options_module
    type(specie_type) specie
    type(options_type) options
    PetscErrorCode ierr

    PetscSizeT sizeofint, sizeofscalar, sizeofbool, sizeofdata
    PetscInt lcv
    character(len=MAXWORDLENGTH):: paramname
    PetscBool help, flag
    write(paramname, '(I1)') specie%id
    
    ! set the species name from options
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)
    if (help) call PetscPrintf(options%comm, "-specie"//trim(paramname)//"_name=<specie"// &
         trim(paramname)//">: name the specie -- for use with parameter options\n", ierr)
    call PetscOptionsGetString(options%my_prefix, "-specie"//trim(paramname)//"_name", &
         specie%name, flag, ierr)
    call RelaxationSetName(specie%relax, specie%name)

    ! create the bag
    call PetscDataTypeGetSize(PETSC_INT, sizeofint, ierr)
    sizeofdata = sizeofint
    call PetscBagCreate(specie%comm, sizeofdata, specie%bag, ierr)
    call PetscBagSetName(specie%bag, TRIM(options%my_prefix)//trim(specie%name), "", ierr)
    call PetscBagGetData(specie%bag, specie%data, ierr)

    ! register data
    call PetscBagRegisterInt(specie%bag, specie%data%phase, 1, &
         trim(options%my_prefix)//'phase_'//trim(specie%name), &
         'Phase of which specie is a component', ierr)
    specie%phase => specie%data%phase

    call RelaxationSetMode(specie%relax, options%transport_relaxation_mode)
    call RelaxationSetFromOptions(specie%relax, options, ierr)
  end subroutine SpecieSetFromOptions

  subroutine SpecieDestroy(specie, ierr)
    type(specie_type) specie
    PetscErrorCode ierr
    if (specie%bag /= 0) call PetscBagDestroy(specie%bag, ierr)
  end subroutine SpecieDestroy
end module LBM_Specie_module
