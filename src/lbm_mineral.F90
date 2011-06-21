!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_mineral.F90
!!!     version:         
!!!     created:         21 June 2011
!!!       on:            10:44:10 MDT
!!!     last modified:   21 June 2011
!!!       at:            10:58:45 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscbagdef.h"

module LBM_Mineral_module
  use petsc
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: mineral_type
    MPI_Comm comm
    PetscInt id
    character(len=MAXWORDLENGTH):: name
  end type mineral_type
  
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

    PetscInt lcv
    character(len=MAXWORDLENGTH):: paramname
    PetscBool help, flag
    write(paramname, '(I1)') mineral%id
    
    ! set the mineral name from options
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)
    if (help) call PetscPrintf(options%comm, "-mineral"//trim(paramname)//"_name=<mineral"// &
         trim(paramname)//">: name the mineral -- for use with parameter options\n", ierr)
    call PetscOptionsGetString(options%my_prefix, "-mineral"//trim(paramname)//"_name", &
         mineral%name, flag, ierr)
  end subroutine MineralSetFromOptions
end module LBM_Mineral_module
