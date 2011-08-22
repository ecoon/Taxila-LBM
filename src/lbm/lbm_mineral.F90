!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_mineral.F90
!!!     version:         
!!!     created:         21 June 2011
!!!       on:            10:44:10 MDT
!!!     last modified:   21 June 2011
!!!       at:            12:43:09 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscbagdef.h"

module LBM_Mineral_module
  use petsc
  use LBM_Mineral_Bag_Data_type_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: mineral_type
    MPI_Comm comm
    PetscInt id
    character(len=MAXWORDLENGTH):: name

    ! bagged parameters
    PetscScalar,pointer,dimension(:) :: gw ! component-mineral force coefs

    type(mineral_bag_data_type),pointer:: data
    PetscBag bag
  end type mineral_type

  interface PetscBagGetData
     subroutine PetscBagGetData(bag, data, ierr)
       use LBM_Mineral_Bag_Data_type_module
       PetscBag bag
       type(mineral_bag_data_type),pointer :: data
       PetscErrorCode ierr
     end subroutine PetscBagGetData
  end interface
  
  interface MineralCreate
    module procedure MineralCreateOne
    module procedure MineralCreateN
  end interface
  
  public :: MineralCreate, &
       MineralDestroy, &
       MineralSetFromOptions, &
       MineralSetName, &
       MineralSetID
  
contains
  
  function MineralCreateOne(comm) result(mineral)
    MPI_Comm comm
    type(mineral_type),pointer :: mineral
    character(len=MAXWORDLENGTH):: name
    allocate(mineral)
    mineral%comm = comm
    mineral%id = 0
    nullify(mineral%gw)
    nullify(mineral%data)
    mineral%bag = 0
    name = 'mineral1'
    call MineralSetName(mineral, name)
  end function MineralCreateOne

  function MineralCreateN(comm, n) result(minerals)
    MPI_Comm comm
    PetscInt n
    type(mineral_type),pointer,dimension(:):: minerals
    type(mineral_type),pointer:: amineral
    character(len=MAXWORDLENGTH):: name
    PetscInt lcv
    allocate(minerals(1:n))
    name = ''

    do lcv=1,n
       amineral => minerals(lcv)
       amineral%comm = comm
       nullify(amineral%gw)
       nullify(amineral%data)
       amineral%bag = 0
       call MineralSetID(amineral, lcv)
       name = 'mineral'//char(lcv+48)
       call MineralSetName(amineral, name)
    end do
  end function MineralCreateN

  subroutine MineralDestroy(mineral, ierr)
    type(mineral_type) mineral
    PetscErrorCode ierr
  end subroutine MineralDestroy

  subroutine MineralSetName(mineral, name)
    type(mineral_type) mineral
    character(len=MAXWORDLENGTH):: name
    mineral%name = name
  end subroutine MineralSetName

  subroutine MineralSetID(mineral, id)
    type(mineral_type) mineral
    PetscInt id
    mineral%id = id
  end subroutine MineralSetID

  subroutine MineralSetFromOptions(mineral, options, ierr)
    use LBM_Options_module
    type(mineral_type) mineral
    type(options_type) options
    PetscErrorCode ierr

    PetscSizeT sizeofscalar, sizeofdata
    PetscInt lcv
    character(len=MAXWORDLENGTH):: paramname,componentname
    PetscBool help, flag
    write(paramname, '(I1)') mineral%id
    
    ! set the mineral name from options
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)
    if (help) call PetscPrintf(options%comm, "-mineral"//trim(paramname)//"_name=<mineral"// &
         trim(paramname)//">: name the mineral -- for use with parameter options\n", ierr)
    call PetscOptionsGetString(options%my_prefix, "-mineral"//trim(paramname)//"_name", &
         mineral%name, flag, ierr)

    ! create the bag
    call PetscDataTypeGetSize(PETSC_SCALAR, sizeofscalar, ierr)
    sizeofdata = (NMAX_PHASES)*sizeofscalar
    call PetscBagCreate(mineral%comm, sizeofdata, mineral%bag, ierr)
    call PetscBagSetName(mineral%bag, TRIM(options%my_prefix)//mineral%name, "", ierr)
    call PetscBagGetData(mineral%bag, mineral%data, ierr)

    ! register data
    do lcv=1,options%ncomponents
      write(paramname, '(I1)') lcv
      call PetscOptionsGetString(options%my_prefix, "-component"//trim(paramname)//"_name", &
           componentname, flag, ierr)
      paramname = 'gw_'//trim(mineral%name)//'_'//trim(componentname)
      call PetscBagRegisterScalar(mineral%bag, mineral%data%gw(lcv), 0.d0, &
           trim(options%my_prefix)//trim(paramname), &
           'mineral-component interaction potential coefficient', ierr)
    end do
    mineral%gw => mineral%data%gw(1:options%ncomponents)
  end subroutine MineralSetFromOptions
end module LBM_Mineral_module
