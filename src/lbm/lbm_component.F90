!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_component.F90
!!!     version:         
!!!     created:         17 March 2011
!!!       on:            13:43:00 MDT
!!!     last modified:   18 October 2011
!!!       at:            13:39:35 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

module LBM_Component_module
  use petsc
  use LBM_Relaxation_module
  use LBM_EOS_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: component_type
     MPI_Comm comm
     character(len=MAXWORDLENGTH):: name
     ! sizes and identifiers (set pre-bag)
     PetscInt s
     PetscInt id
     PetscScalar time_scale

     PetscScalar mm ! molecular mass
     PetscScalar viscosity ! kinematic viscoscity, in m^2/s
     PetscScalar density ! reference density, in kg/m^3
     PetscScalar,pointer,dimension(:) :: gf ! component-component force coefs

     type(relaxation_type),pointer:: relax
     type(eos_type),pointer :: eos
  end type component_type

  interface ComponentCreate
     module procedure ComponentCreateOne
     module procedure ComponentCreateN
  end interface

  public :: ComponentCreate, &
       ComponentDestroy, &
       ComponentSetSizes, &
       ComponentSetName, &
       ComponentSetID, &
       ComponentSetFromOptions

contains
  function ComponentCreateOne(comm) result(component)
    MPI_Comm comm
    type(component_type),pointer :: component
    character(len=MAXWORDLENGTH):: name
    allocate(component)
    component%comm = comm
    call ComponentInitialize(component)
    component%relax => RelaxationCreate(comm)
    name = 'component1'
    call ComponentSetName(component, name)
  end function ComponentCreateOne

  function ComponentCreateN(comm, n) result(components)
    MPI_Comm comm
    PetscInt n
    type(component_type),pointer,dimension(:):: components
    type(component_type),pointer:: acomponent
    character(len=MAXWORDLENGTH):: name
    PetscInt lcv
    allocate(components(1:n))

    do lcv=1,n
       acomponent => components(lcv)
       acomponent%comm = comm
       call ComponentInitialize(acomponent)
       acomponent%relax => RelaxationCreate(comm)
       call ComponentSetID(acomponent, lcv)
       name = 'component'//char(lcv+48)
       call ComponentSetName(acomponent, name)
    end do
  end function ComponentCreateN

  subroutine ComponentInitialize(component)
    type(component_type) component
    component%name = ''
    component%s = -1
    component%id = 0
    component%time_scale = 0

    component%mm = 1.d0
    component%density = -999.d0
    component%viscosity = -999.d0
    nullify(component%gf)
    nullify(component%relax)
    nullify(component%eos)

  end subroutine ComponentInitialize

  subroutine ComponentSetSizes(component, s, b)
    type(component_type) :: component
    PetscInt s,b
    component%s = s
    call RelaxationSetSizes(component%relax, s, b)
  end subroutine ComponentSetSizes

  subroutine ComponentSetName(component, name)
    type(component_type) component
    character(len=MAXWORDLENGTH):: name
    component%name = name
    call RelaxationSetName(component%relax, name)
    if (associated(component%eos)) then
      call EOSSetName(component%eos, name)
    end if
  end subroutine ComponentSetName

  subroutine ComponentSetID(component, id)
    type(component_type) component
    PetscInt id
    component%id = id
    call RelaxationSetID(component%relax, id)
    if (associated(component%eos)) then
      call EOSSetID(component%eos, id)
    end if
  end subroutine ComponentSetID

  subroutine ComponentSetFromOptions(component, options, ierr)
    use LBM_Options_module
    type(component_type) component
    type(options_type) options
    PetscErrorCode ierr

    PetscInt lcv
    PetscBool flag
    character(len=MAXWORDLENGTH):: idstring
    write(idstring, '(I1)') component%id

    ! set the component name from options
    call OptionsGetString(options, "-"//trim(component%name)//"_name", &
         "name the component", component%name, flag, ierr)
    call OptionsGroupHeader(options, " "//trim(component%name)//" Options", ierr)

    call RelaxationSetName(component%relax, component%name)
    call RelaxationSetMode(component%relax, options%flow_relaxation_mode)
    call RelaxationSetFromOptions(component%relax, options, ierr)

    if (associated(component%eos)) then
      call EOSSetName(component%eos, component%name)
      call EOSSetFromOptions(component%eos, options, ierr)
    end if

    call OptionsGetReal(options, "-mm_"//trim(component%name), &
         "molecular mass", component%mm, flag, ierr)
    component%relax%d_k = 1. - 2./(3.*component%mm) ! d_k = 1/3 for mm=1

    call OptionsGetReal(options, "-viscosity_"//trim(component%name), &
         "kinematic viscosity, defaults to nondimensional value", &
         component%viscosity, flag, ierr)
    call OptionsGetReal(options, "-density_"//trim(component%name), &
         "density scale, defaults to nondimensional value", component%density, flag, ierr)

    allocate(component%gf(component%s))
    component%gf(:) = 0.d0
    do lcv=1,component%s
       write(idstring, '(I1, I1)') component%id, lcv
       call OptionsGetReal(options, "-g_"//trim(idstring), &
            "component-component interaction potential coefficient", component%gf(lcv), &
            flag, ierr)
    end do

    call OptionsGroupFooter(options, " "//trim(component%name)//" Options", ierr)
  end subroutine ComponentSetFromOptions

  subroutine ComponentDestroy(component, ierr)
    type(component_type) component
    PetscErrorCode ierr
    if (associated(component%gf)) deallocate(component%gf)
    if (associated(component%relax)) call RelaxationDestroy(component%relax, ierr)
    if (associated(component%eos)) call EOSDestroy(component%eos, ierr)
  end subroutine ComponentDestroy
end module LBM_Component_module
