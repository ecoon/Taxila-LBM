!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_component.F90
!!!     version:         
!!!     created:         17 March 2011
!!!       on:            13:43:00 MDT
!!!     last modified:   23 August 2011
!!!       at:            16:39:11 MDT
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
  use LBM_EOS_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: component_type
     MPI_Comm comm
     ! sizes and identifiers (set pre-bag)
     PetscInt s
     PetscInt id
     PetscScalar time_scale
     
     ! bagged parameters
     PetscScalar,pointer :: mm ! molecular mass
     PetscScalar,pointer :: viscosity ! kinematic viscoscity, in m^2/s
     PetscScalar,pointer :: density ! reference density, in kg/m^3
     PetscScalar,pointer,dimension(:) :: gf ! component-component force coefs

     ! dependent parameters, for equilf and collision
     type(relaxation_type),pointer:: relax
     type(eos_type),pointer :: eos

     ! bag 
     character(len=MAXWORDLENGTH):: name
     type(component_bag_data_type),pointer:: data
     PetscBag bag
  end type component_type

  interface PetscBagGetData
     subroutine PetscBagGetData(bag, data, ierr)
       use LBM_Component_Bag_Data_type_module
       PetscBag bag
       type(component_bag_data_type),pointer :: data
       PetscErrorCode ierr
     end subroutine PetscBagGetData
  end interface

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
    name = ''

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
    component%s = -1
    component%id = 0
    component%time_scale = 0

    nullify(component%mm)
    nullify(component%gf)
    nullify(component%relax)
    nullify(component%eos)

    component%name = ''
    nullify(component%data)
    component%bag = 0
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

    character(len=MAXWORDLENGTH):: name
    PetscSizeT sizeofint, sizeofscalar, sizeofbool, sizeofdata
    PetscInt lcv
    character(len=MAXWORDLENGTH):: paramname
    PetscBool help, flag
    write(paramname, '(I1)') component%id
    
    ! set the component name from options
    name = ''
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)
    if (help) call PetscPrintf(options%comm, "-component"//trim(paramname)// &
         "_name=<component"//trim(paramname)// &
         ">: name the component -- for use with parameter options\n", ierr)
    call PetscOptionsGetString(options%my_prefix, "-component"//trim(paramname)//"_name", &
         name, flag, ierr)
    if (flag) then
      call ComponentSetName(component, name)
    end if

    call RelaxationSetMode(component%relax, options%flow_relaxation_mode)
    call RelaxationSetFromOptions(component%relax, options, ierr)
    if (associated(component%eos)) then
      call EOSSetFromOptions(component%eos, options, ierr)
    end if

    ! create the bag
    call PetscDataTypeGetSize(PETSC_SCALAR, sizeofscalar, ierr)
    sizeofdata = (4+NMAX_COMPONENTS)*sizeofscalar
    call PetscBagCreate(component%comm, sizeofdata, component%bag, ierr)
    call PetscBagSetName(component%bag, TRIM(options%my_prefix)//component%name, "", ierr)
    call PetscBagGetData(component%bag, component%data, ierr)

    call PetscBagRegisterScalar(component%bag, component%data%mm, 1.d0, &
         trim(options%my_prefix)//'mm_'//trim(component%name), 'molecular mass', ierr)
    component%mm => component%data%mm

    call PetscBagRegisterScalar(component%bag, component%data%viscosity, -999.d0, &
         trim(options%my_prefix)//'viscosity_'//trim(component%name), 'kinematic viscosity [m^2/s], defaults to ND value', &
         ierr)
    component%viscosity => component%data%viscosity

    call PetscBagRegisterScalar(component%bag, component%data%density, -999.d0, &
         trim(options%my_prefix)//'density_'//trim(component%name), &
         'density [kg/m^3], defaults to ND value', ierr)
    component%density => component%data%density

    component%relax%d_k = 1.d0 - 2.d0/(3.d0*component%mm) ! d_k = 1/3 for mm=1

    do lcv=1,component%s
       write(paramname, '(I1, I1)') component%id, lcv
       call PetscBagRegisterScalar(component%bag, component%data%gf(lcv), 0.d0, &
            trim(options%my_prefix)//'g_'//trim(paramname), 'component-component interaction potential coefficient', ierr)
    end do
    component%gf => component%data%gf(1:options%ncomponents)

  end subroutine ComponentSetFromOptions

  subroutine ComponentDestroy(component, ierr)
    type(component_type) component
    PetscErrorCode ierr
    if (associated(component%relax)) call RelaxationDestroy(component%relax, ierr)
    if (associated(component%eos)) call EOSDestroy(component%eos, ierr)
    if (component%bag /= 0) call PetscBagDestroy(component%bag, ierr)
  end subroutine ComponentDestroy
end module LBM_Component_module
