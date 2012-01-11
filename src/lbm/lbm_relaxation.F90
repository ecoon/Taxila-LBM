!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_relaxation.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            15:15:25 MDT
!!!     last modified:   31 October 2011
!!!       at:            15:23:36 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscbagdef.h"

module LBM_Relaxation_module
  use petsc
  use LBM_Relaxation_Bag_Data_type_module
  use LBM_Distribution_Function_type_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: relaxation_type
     MPI_Comm comm
     PetscInt id
     PetscInt s
     PetscInt b
     PetscInt mode
     
     PetscScalar,pointer :: tau ! relaxation time
     PetscScalar,pointer :: s1  ! MRT relaxation time
     PetscScalar,pointer :: s2  ! MRT relaxation time
     PetscScalar,pointer :: s3  ! MRT relaxation time
     PetscScalar,pointer :: s4  ! MRT relaxation time (only for 3D)
     PetscScalar,pointer :: s5  ! MRT relaxation time (only for 3D)
     PetscScalar,pointer,dimension(:) :: tau_mrt ! species of S vector for mrt

     ! dependents, set by discretization
     PetscScalar d_k
     PetscScalar c_s2

     ! bag 
     character(len=MAXWORDLENGTH):: name
     type(relaxation_bag_data_type),pointer:: data
     PetscBag bag
  end type relaxation_type
     
  interface PetscBagGetData
     subroutine PetscBagGetData(bag, data, ierr)
       use LBM_Relaxation_Bag_Data_type_module
       PetscBag bag
       type(relaxation_bag_data_type),pointer :: data
       PetscErrorCode ierr
     end subroutine PetscBagGetData
  end interface

  public :: RelaxationCreate, &
       RelaxationDestroy, &
       RelaxationSetSizes, &
       RelaxationSetName, &
       RelaxationSetID, &
       RelaxationSetMode, &
       RelaxationSetFromOptions, &
       RelaxationCollide

contains
  function RelaxationCreate(comm) result(relax)
    MPI_Comm comm
    type(relaxation_type),pointer:: relax
    allocate(relax)
    relax%comm = comm
    relax%id = 0
    relax%s = -1
    relax%b = -1
    relax%mode = RELAXATION_MODE_SRT

    nullify(relax%tau)
    nullify(relax%s1)
    nullify(relax%s2)
    nullify(relax%s3)
    nullify(relax%s4)
    nullify(relax%s5)
    nullify(relax%tau_mrt)

    relax%d_k = 0. 
    relax%c_s2 = 1.d0/3.

    nullify(relax%data)
    relax%bag = 0
    relax%name = ''
  end function RelaxationCreate

  subroutine RelaxationDestroy(relax, ierr)
    type(relaxation_type) relax
    PetscErrorCode ierr
    if (relax%bag /= 0) call PetscBagDestroy(relax%bag, ierr)
    if (associated(relax%tau_mrt)) deallocate(relax%tau_mrt)
  end subroutine RelaxationDestroy

  subroutine RelaxationSetSizes(relax, s, b)
    type(relaxation_type) relax
    PetscInt:: s, b
    relax%s = s
    relax%b = b
    allocate(relax%tau_mrt(0:b))
  end subroutine RelaxationSetSizes

  subroutine RelaxationSetMode(relax, mode)
    type(relaxation_type) relax
    PetscInt mode
    relax%mode = mode
  end subroutine RelaxationSetMode

  subroutine RelaxationSetName(relax, name)
    type(relaxation_type) relax
    character(len=MAXWORDLENGTH):: name
    
    relax%name = name
  end subroutine RelaxationSetName

  subroutine RelaxationSetID(relax, id)
    type(relaxation_type) relax
    PetscInt id
    relax%id = id
  end subroutine RelaxationSetID

  subroutine RelaxationSetFromOptions(relax, options, ierr)
    use LBM_Options_module
    type(relaxation_type) relax
    type(options_type) options
    PetscErrorCode ierr

    PetscSizeT sizeofscalar, sizeofdata
    PetscInt lcv
    character(len=MAXWORDLENGTH):: paramname, paramname2

    ! create the bag
    call PetscDataTypeGetSize(PETSC_SCALAR, sizeofscalar, ierr)
    sizeofdata = 6*sizeofscalar 
    sizeofdata = sizeofdata + sizeofscalar

    call PetscBagCreate(relax%comm, sizeofdata, relax%bag, ierr)
    call PetscBagGetData(relax%bag, relax%data, ierr)

    call PetscBagRegisterScalar(relax%bag, relax%data%tau, ONE_S, &
         trim(options%my_prefix)//'tau_'//trim(relax%name), 'relaxation time', ierr)
    relax%tau => relax%data%tau

    call PetscBagRegisterScalar(relax%bag, relax%data%s1, ONE_S, &
         trim(options%my_prefix)//'s1_'//trim(relax%name), 'MRT relaxation time', ierr)
    relax%s1 => relax%data%s1

    call PetscBagRegisterScalar(relax%bag, relax%data%s2, ONE_S, &
         trim(options%my_prefix)//'s2_'//trim(relax%name), 'MRT relaxation time', ierr)
    relax%s2 => relax%data%s2

    call PetscBagRegisterScalar(relax%bag, relax%data%s3, ONE_S, &
         trim(options%my_prefix)//'s3_'//trim(relax%name), 'MRT relaxation time', ierr)
    relax%s3 => relax%data%s3
   
    call PetscBagRegisterScalar(relax%bag, relax%data%s4, ONE_S, &
         trim(options%my_prefix)//'s4_'//trim(relax%name), 'MRT relaxation time', ierr)
    relax%s4 => relax%data%s4

    call PetscBagRegisterScalar(relax%bag, relax%data%s5, ONE_S, &
         trim(options%my_prefix)//'s5_'//trim(relax%name), 'MRT relaxation time', ierr)
    relax%s5 => relax%data%s5

    call PetscBagSetName(relax%bag, TRIM(options%my_prefix)//relax%name, "", ierr)
  end subroutine RelaxationSetFromOptions

  subroutine RelaxationCollide(relax, fi, fi_eq, m, dist)
    type(relaxation_type) relax
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, 0:dist%b):: fi
    PetscScalar,dimension(dist%s, 0:dist%b):: fi_eq
    PetscInt m
    PetscErrorCode ierr

    select case(relax%mode)
    case(RELAXATION_MODE_MRT)
      call RelaxationCollideMRT(relax, fi, fi_eq, m, dist)
    case(RELAXATION_MODE_SRT)
      call RelaxationCollideSRT(relax, fi, fi_eq, m, dist)
    case DEFAULT
      SETERRQ(1, 1, 'invalid relaxation mode in LBM', ierr)
    end select
  end subroutine RelaxationCollide

  subroutine RelaxationCollideSRT(relax, fi, fi_eq, m, dist)
    ! input variables
    type(relaxation_type) relax
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, 0:dist%b):: fi, fi_eq
    PetscInt m

    fi(m,:) = fi(m,:) - (fi(m,:)-fi_eq(m,:))/relax%tau
    return
  end subroutine RelaxationCollideSRT

  subroutine RelaxationCollideMRT(relax, fi, fi_eq, m, dist)
    ! input variables
    type(relaxation_type) relax
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, 0:dist%b):: fi, fi_eq
    PetscInt m

    ! local variables
    PetscInt n                  ! loop variables
    PetscScalar :: mdiff
    PetscScalar :: d_fi(0:dist%b)
   
    d_fi = fi(m,:)-fi_eq(m,:)
    do n=0,dist%b
      mdiff = sum(dist%disc%mt_mrt(:,n)*d_fi,1)
      fi(m,:) = fi(m,:) - &
           relax%tau_mrt(n)*mdiff/dist%disc%mmt_mrt(n)*dist%disc%mt_mrt(:,n)
    enddo
  end subroutine RelaxationCollideMRT
end module LBM_Relaxation_module

    
