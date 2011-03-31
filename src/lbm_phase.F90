!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_phase.F90
!!!     version:         
!!!     created:         17 March 2011
!!!       on:            13:43:00 MDT
!!!     last modified:   31 March 2011
!!!       at:            09:47:37 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscbagdef.h"

module LBM_Phase_module
  use petsc
  use LBM_Phase_Bag_Data_type_module
  use LBM_Relaxation_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: phase_type
     MPI_Comm comm
     ! sizes and identifiers (set pre-bag)
     PetscInt s
     PetscInt id

     ! bagged parameters
     PetscScalar,pointer :: mm ! molecular mass
     PetscScalar,pointer :: gw ! solid affinity? for phase-wall interaction forces
     PetscScalar,pointer,dimension(:) :: gf ! phase-phase force coefs
     PetscScalar,pointer :: tau ! relaxation time

     ! dependent parameters
     PetscScalar alpha_0, alpha_1 
     PetscScalar d_k
     PetscScalar c_s2
     type(relaxation_type),pointer:: relax

     ! bag 
     character(len=MAXWORDLENGTH):: name
     type(phase_bag_data_type),pointer:: data
     PetscBag bag
  end type phase_type

  interface PetscBagGetData
     subroutine PetscBagGetData(bag, data, ierr)
       use LBM_Phase_Bag_Data_type_module
       PetscBag bag
       type(phase_bag_data_type),pointer :: data
       PetscErrorCode ierr
     end subroutine PetscBagGetData
  end interface

  interface PhaseCreate
     procedure PhaseCreateOne, PhaseCreateN
  end interface

  public :: PhaseCreate, &
       PhaseDestroy, &
       PhaseSetSizes, &
       PhaseSetName, &
       PhaseSetID, &
       PhaseSetFromOptions

contains
  function PhaseCreateOne(comm) result(phase)
    MPI_Comm comm
    type(phase_type),pointer :: phase
    allocate(phase)
    phase%comm = comm
    call PhaseInitialize(phase)
    phase%relax => RelaxationCreate(comm)
  end function PhaseCreateOne

  function PhaseCreateN(comm, n) result(phases)
    MPI_Comm comm
    PetscInt n
    type(phase_type),pointer,dimension(:):: phases
    type(phase_type),pointer:: aphase
    PetscInt lcv
    allocate(phases(1:n))

    do lcv=1,n
       aphase => phases(lcv)
       aphase%comm = comm
       call PhaseInitialize(aphase)
       aphase%relax => RelaxationCreate(comm)
       call PhaseSetID(aphase, lcv)
    end do
  end function PhaseCreateN

  subroutine PhaseInitialize(phase)
    type(phase_type) phase
    phase%s = -1
    phase%id = 0

    nullify(phase%mm)
    nullify(phase%gw)
    nullify(phase%gf)
    nullify(phase%tau)
    nullify(phase%relax)

    phase%alpha_0 = 0.
    phase%alpha_1 = 0.
    phase%d_k = 0.
    phase%c_s2 = 1.d0/3.d0

    phase%name = ''
    nullify(phase%data)
    phase%bag = 0
  end subroutine PhaseInitialize
  
  subroutine PhaseSetSizes(phase, s, b)
    type(phase_type) :: phase
    PetscInt s,b
    phase%s = s
    call RelaxationSetSizes(phase%relax, s, b)
  end subroutine PhaseSetSizes

  subroutine PhaseSetName(phase, name)
    type(phase_type) phase
    character(len=MAXWORDLENGTH):: name
    
    phase%name = name
    call RelaxationSetName(phase%relax, name)
  end subroutine PhaseSetName

  subroutine PhaseSetID(phase, id)
    type(phase_type) phase
    PetscInt id
    phase%id = id
    call RelaxationSetID(phase%relax, id)
  end subroutine PhaseSetID

  subroutine PhaseSetFromOptions(phase, options, ierr)
    use LBM_Options_module
    type(phase_type) phase
    type(options_type) options
    PetscErrorCode ierr

    PetscSizeT sizeofint, sizeofscalar, sizeofbool, sizeofdata
    PetscInt lcv
    character(len=MAXWORDLENGTH):: paramname
    write(paramname, '(I1)') phase%id

    ! create the bag
    call PetscDataTypeGetSize(PETSC_SCALAR, sizeofscalar, ierr)
    sizeofdata = (3+phase%s)*sizeofscalar
    call PetscBagCreate(phase%comm, sizeofdata, phase%bag, ierr)
    call PetscBagSetName(phase%bag, TRIM(options%my_prefix)//phase%name, "", ierr)
    call PetscBagGetData(phase%bag, phase%data, ierr)

    ! register data
    call PetscBagRegisterScalar(phase%bag, phase%data%gw, 0.d0, &
         'gw'//paramname, 'Phase-solid interaction potential coefficient', ierr)
    phase%gw => phase%data%gw
    call PetscBagRegisterScalar(phase%bag, phase%data%mm, 1.d0, &
         'mm'//paramname, 'molecular mass', ierr)
    phase%mm => phase%data%mm
    phase%d_k = 1.d0 - 2.d0/(3.d0*phase%mm)

    do lcv=1,phase%s
       write(paramname, '(I1, I1)') lcv, phase%id
       call PetscBagRegisterScalar(phase%bag, phase%data%gf(lcv), 0.d0, &
            'g_'//paramname, 'phase-phase interaction potential coefficient', ierr)
    end do
    phase%gf => phase%data%gf(1:options%nphases)

    call RelaxationSetMode(phase%relax, options%flow_relaxation_mode)
    call RelaxationSetFromOptions(phase%relax, options, ierr)
    phase%tau => phase%relax%tau
  end subroutine PhaseSetFromOptions

  subroutine PhaseDestroy(phase, ierr)
    type(phase_type) phase
    PetscErrorCode ierr
    if (phase%bag /= 0) call PetscBagDestroy(phase%bag, ierr)
  end subroutine PhaseDestroy
end module LBM_Phase_module
