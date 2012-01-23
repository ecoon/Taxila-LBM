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

module LBM_Relaxation_module
  use petsc
  use LBM_Distribution_Function_type_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: relaxation_type
     MPI_Comm comm
     character(len=MAXWORDLENGTH):: name
     PetscInt id
     PetscInt s
     PetscInt b
     PetscInt mode

     PetscScalar tau ! relaxation time
     PetscScalar s1  ! MRT relaxation time
     PetscScalar s2  ! MRT relaxation time
     PetscScalar s3  ! MRT relaxation time
     PetscScalar s4  ! MRT relaxation time (only for 3D)
     PetscScalar s5  ! MRT relaxation time (only for 3D)
     PetscScalar,pointer,dimension(:) :: tau_mrt ! species of S vector for mrt

     ! dependents, set by discretization
     PetscScalar d_k
     PetscScalar c_s2
  end type relaxation_type

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
    relax%name = ''
    relax%id = 0
    relax%s = -1
    relax%b = -1
    relax%mode = RELAXATION_MODE_SRT

    relax%tau = 1.d0
    relax%s1 = 1.d0
    relax%s2 = 1.d0
    relax%s3 = 1.d0
    relax%s4 = 1.d0
    relax%s5 = 1.d0
    nullify(relax%tau_mrt)

    relax%d_k = 0.d0
    relax%c_s2 = 1.d0/3.
  end function RelaxationCreate

  subroutine RelaxationDestroy(relax, ierr)
    type(relaxation_type) relax
    PetscErrorCode ierr
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
    PetscBool flag
    PetscErrorCode ierr

    PetscInt lcv

    call OptionsGroupHeader(options, "  "//trim(relax%name)//" Relaxation Options", ierr)

    ! single and multiple relaxation times
    call OptionsGetReal(options, "-tau_"//trim(relax%name), "relaxation time", &
         relax%tau, flag, ierr)

    relax%s1 = 1.d0/relax%tau
    call OptionsGetReal(options, "-s1_"//trim(relax%name), "MRT relaxation time", &
         relax%s1, flag, ierr)
    relax%s2 = 1.d0/relax%tau
    call OptionsGetReal(options, "-s2_"//trim(relax%name), "MRT relaxation time", &
         relax%s2, flag, ierr)
    relax%s3 = 1.d0/relax%tau
    call OptionsGetReal(options, "-s3_"//trim(relax%name), "MRT relaxation time", &
         relax%s3, flag, ierr)
    relax%s4 = 1.d0/relax%tau
    call OptionsGetReal(options, "-s4_"//trim(relax%name), "MRT relaxation time", &
         relax%s4, flag, ierr)
    relax%s5 = 1.d0/relax%tau
    call OptionsGetReal(options, "-s5_"//trim(relax%name), "MRT relaxation time", &
         relax%s5, flag, ierr)
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

    
