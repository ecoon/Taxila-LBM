!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_phase.F90
!!!     version:         
!!!     created:         17 March 2011
!!!       on:            13:43:00 MDT
!!!     last modified:   17 March 2011
!!!       at:            17:38:21 MDT
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
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: phase_type
     MPI_Comm comm
     ! sizes and identifiers (set pre-bag)
     PetscInt s
     PetscInt b
     PetscInt id

     ! bagged parameters
     PetscScalar,pointer :: mm ! molecular mass
     PetscScalar,pointer :: tau ! relaxation time
     PetscScalar,pointer,dimension(:) :: tau_mrt ! components of S vector for mrt
     PetscScalar,pointer :: gw ! solid affinity? for phase-wall interaction forces
     PetscScalar,pointer,dimension(:) :: gf ! phase-phase force coefs
     PetscBool, pointer:: mrt ! do MRT?

     ! dependent parameters
     PetscScalar alpha_0, alpha_1 
     PetscScalar d_k
     PetscScalar c_s2

     ! bag 
     character(len=MAXWORDLENGTH):: name
     type(phase_bag_data_type),pointer:: data
     PetscBag bag
  end type phase_bag_type

  interface PetscBagGetData
     subroutine PetscBagGetData(bag, data, ierr)
       PetscBag bag
       type(phase_bag_data_type),pointer :: data
       PetscErrorCode ierr
     end subroutine PetscBagGetData
  end interface

  public :: PhaseCreate, &
       PhaseDestroy, &
       PhaseSetUp!, &
!       PhaseSetFromOptions

contains
  function PhaseCreate(comm) result(phase)
    MPI_Comm comm
    type(phase_type),pointer :: phase
    allocate(phase)
    phase%comm = comm
    phase%s = -1
    phase%b = -1
    phase%id = 0
    phase%alpha_0 = 0.
    phase%alpha_1 = 0.
    phase%d_k = 0.
    phase%c_s2 = 1.d0/3.d0
    phase%name = ''
    nullify(phase%data)
    phase%bag = 0
  end function PhaseCreate

  subroutine PhaseSetSizes(phase, s, b)
    type(phase_type),pointer :: phase
    PetscInt s,b
    phase%s = s
    phase%b = b
  end subroutine PhaseSetSizes
  
  subroutine PhaseSetName(phase, name, id)
    type(phase_type) phase
    character(len=MAXWORDLENGTH):: name
    PetscInt id
    
    phase%name = name
    phase%id = id
  end subroutine PhaseSetName

  subroutine LBMPhaseSetFromOptions(phase, options, ierr)
    use LBM_Options_module
    type(phase_type) phase
    type(options_type) options
    PetscErrorCode ierr

    PetscInt sizeofint, sizeofscalar, sizeofbool, sizeofdata
    PetscInt lcv
    character(len=MAXWORDLENGTH):: paramname

    ! set up the data
    allocate(phase%data%tau_mrt(0%phase%b))
    allocate(phase%data%gf(1%phase%s))

    ! create the bag
    call PetscDataTypeGetSize(PETSC_SCALAR, sizeofscalar, ierr)
    call PetscDataTypeGetSize(PETSC_BOOL, sizeofbool, ierr)
    sizeofdata = (3+b+1+s)*sizeofscalar + sizeofbool
    call PetscBagCreate(phase%comm, sizeofdata, phase%bag, ierr)
    call PetscBagSetName(phase%bag, TRIM(options%my_prefix)//phase%name, ierr)

    ! register data
    call PetscBagGetData(phase%bag, phase%data, ierr)
    call PetscBagRegisterBool(phase%bag, phase%data%mrt, PETSC_FALSE, 'mrt', &
         'Use MRT on this phase?', ierr)
    phase%mrt => phase%data%mrt
    
    if (phase%mrt) then
       do lcv=0,b
          write(paramname, '(I0.2)') lcv
          call PetscBagRegisterScalar(phase%bag, phase%data%tau_mrt(lcv), 0.d0, &
               'tau_'//paramname, 'MRT relaxation moment coefficient', ierr)
       end do
    else
       phase%data%tau_mrt = 0.d0
    end if
    phase%tau_mrt => phase%data%tau_mrt
    
    do lcv=1,s
       write(paramname, '(I1, I1)') lcv, phase
       call PetscBagRegisterScalar(phase%bag, phase%data%gf(lcv), 0.d0, &
            'g_'//paramname, 'phase-phase interaction potential coefficient', ierr)
    end do
    phase%gf => phase%data%gf

    call PetscBagRegisterScalar(phase%bag, phase%data%gw, 0.d0, &
         'gw', 'Phase-solid interaction potential coefficient', ierr)
    phase%gw => phase%data%gw
    call PetscBagRegisterScalar(phase%bag, phase%data%tau, 1.d0, &
         'tau', 'relaxation time', ierr)
    phase%tau => phase%data%tau
    call PetscBagRegisterScalar(phase%bag, phase%data%mm, 1.d0, &
         'mm', 'molecular mass', ierr)
    phase%mm => phase%data%mm
  end subroutine LBMPhaseSetFromOptions

  subroutine LBMPhaseDestroy(phase, ierr)
    if (phase%bag /= 0) call PetscBagDestroy(phase%bag, ierr)
  end subroutine LBMPhaseDestroy
end module LBM_Phase_module
