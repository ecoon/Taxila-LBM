!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_specie.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            15:34:44 MDT
!!!     last modified:   18 October 2011
!!!       at:            13:50:51 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

module LBM_Specie_module
  use petsc
  use LBM_Relaxation_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: specie_type
     MPI_Comm comm
     character(len=MAXWORDLENGTH):: name
     PetscInt s
     PetscInt id

     PetscInt component
     type(relaxation_type),pointer:: relax
  end type specie_type

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
    nullify(specie%relax)
    specie%component = -1
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
    PetscBool flag
    PetscErrorCode ierr

    character(len=MAXWORDLENGTH):: idstring
    write(idstring, '(I1)') specie%id

    ! set the species name from options
    call OptionsGetString(options, "-"//trim(specie%name)//"_name", &
         "name the chemical specie", specie%name, flag, ierr)
    call OptionsGroupHeader(options, " "//trim(specie%name)//" Options", ierr)
    call RelaxationSetName(specie%relax, specie%name)

    call OptionsGetInt(options, "-component_"//trim(specie%name), &
         "Component id in which specie exists", specie%component, flag, ierr)

    call RelaxationSetMode(specie%relax, options%transport_relaxation_mode)
    call RelaxationSetFromOptions(specie%relax, options, ierr)
    call OptionsGroupFooter(options, " "//trim(specie%name)//" Options", ierr)
  end subroutine SpecieSetFromOptions

  subroutine SpecieDestroy(specie, ierr)
    type(specie_type) specie
    PetscErrorCode ierr
    if (associated(specie%relax)) call RelaxationDestroy(specie%relax, ierr)
  end subroutine SpecieDestroy
end module LBM_Specie_module
