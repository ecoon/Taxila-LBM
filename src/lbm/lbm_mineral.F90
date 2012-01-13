!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_mineral.F90
!!!     version:         
!!!     created:         21 June 2011
!!!       on:            10:44:10 MDT
!!!     last modified:   18 October 2011
!!!       at:            13:38:22 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

module LBM_Mineral_module
  use petsc
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: mineral_type
    MPI_Comm comm
    PetscInt id
    character(len=MAXWORDLENGTH):: name

    PetscScalar,pointer,dimension(:) :: gw ! component-mineral force coefs
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
    call MineralInitialize(mineral)

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

    do lcv=1,n
       amineral => minerals(lcv)
       amineral%comm = comm
       call MineralInitialize(amineral)
       call MineralSetID(amineral, lcv)
       name = 'mineral'//char(lcv+48)
       call MineralSetName(amineral, name)
    end do
  end function MineralCreateN

  subroutine MineralInitialize(mineral)
    type(mineral_type) mineral
    mineral%name = ''
    mineral%id = 0
    nullify(mineral%gw)
  end subroutine MineralInitialize

  subroutine MineralDestroy(mineral, ierr)
    type(mineral_type) mineral
    PetscErrorCode ierr
    if (associated(mineral%gw)) deallocate(mineral%gw)
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
    PetscBool flag
    character(len=MAXWORDLENGTH):: idstring, componentname
    write(idstring, '(I1)') mineral%id

    ! set the mineral name from options
    call OptionsGetString(options, "-"//trim(mineral%name)//"_name", &
         "name the mineral", mineral%name, flag, ierr)
    call OptionsGroupHeader(options, " "//trim(mineral%name)//" Options", ierr)

    ! register data
    allocate(mineral%gw(options%ncomponents))
    mineral%gw(:) = 0.d0
    do lcv=1,options%ncomponents
      write(idstring, '(I1)') lcv
      componentname = "component"//trim(idstring)
      call PetscOptionsGetString(options%my_prefix, "-component"//trim(idstring)//"_name", &
           componentname, flag, ierr)
      idstring = '-gw_'//trim(mineral%name)//'_'//trim(componentname)
      call OptionsGetReal(options, trim(idstring), &
           "mineral-component interaction potential coefficient", mineral%gw(lcv), &
           flag, ierr)
    end do

    call OptionsGroupFooter(options, " "//trim(mineral%name)//" Options", ierr)
  end subroutine MineralSetFromOptions
end module LBM_Mineral_module
