!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_eos.F90
!!!     version:         
!!!     created:         22 August 2011
!!!       on:            15:53:59 MDT
!!!     last modified:   18 October 2011
!!!       at:            13:43:19 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================


#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscbagdef.h"

module LBM_EOS_module
  use petsc
  use LBM_Error_module
  use LBM_Distribution_Function_type_module
  implicit none

private
#include "lbm_definitions.h"

  type, public :: eos_type
     MPI_Comm comm
     PetscInt id
     PetscInt eos_type

     PetscScalar :: omega
     PetscScalar :: T, T_c
     PetscScalar :: p_c
     PetscScalar :: R
     PetscScalar :: a,b
     PetscScalar :: psi0, rho0
     character(len=MAXWORDLENGTH):: name
   end type eos_type

   public :: EOSCreate, &
        EOSDestroy, &
        EOSSetName, &
        EOSSetID, &
        EOSSetFromOptions, &
        EOSApply

contains
  function EOSCreate(comm) result(eos)
    MPI_Comm comm
    type(eos_type),pointer :: eos
    allocate(eos)
    eos%comm = comm
    eos%id = 0
    eos%eos_type = EOS_NULL

    eos%omega = 0.344
    eos%T = 0.
    eos%T_c = 0.
    eos%p_c = 0.
    eos%R = 1.
    eos%a = 0.
    eos%b = 0.
    eos%psi0 = 0.
    eos%rho0 = 0.

    eos%name = ''
  end function EOSCreate

  subroutine EOSDestroy(eos, ierr)
    type(eos_type) :: eos
    PetscErrorCode ierr
  end subroutine EOSDestroy

  subroutine EOSSetName(eos, name)
    type(eos_type) eos
    character(len=MAXWORDLENGTH):: name
    
    eos%name = name
  end subroutine EOSSetName

  subroutine EOSSetID(eos, id)
    type(eos_type) eos
    PetscInt id
    
    eos%id = id
  end subroutine EOSSetID

  subroutine EOSSetFromOptions(eos, options, ierr)
    use String_module
    use LBM_Options_module
    type(eos_type) eos
    type(options_type) options
    PetscErrorCode ierr

    character(len=MAXWORDLENGTH) eos_name
    character(len=MAXWORDLENGTH) test_eos_name
    PetscBool help, flag, done

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    if (help) call PetscPrintf(options%comm, '  -eos_'//trim(eos%name)//&
         '={rho,sc,pr,thermo}: name of equation of state\n', ierr)
    call PetscOptionsGetString(options%my_prefix, '-eos_'//trim(eos%name), eos_name, &
         flag, ierr)

    done = PETSC_FALSE
    test_eos_name = 'rho'
    if (StringCompareIgnoreCase(test_eos_name, eos_name, 3)) then
      eos%eos_type = EOS_DENSITY
      call EOSSetFromOptions_Rho(eos, options, ierr)
      done = PETSC_TRUE
    end if

    if (.not.done) then
      test_eos_name = 'sc'
      if (StringCompareIgnoreCase(test_eos_name, eos_name, 2)) then
        eos%eos_type = EOS_SC
        call EOSSetFromOptions_SC(eos, options, ierr)
        done = PETSC_TRUE
      end if
    end if

    if (.not.done) then
      test_eos_name = 'thermo'
      if (StringCompareIgnoreCase(test_eos_name, eos_name, 6)) then
        eos%eos_type = EOS_THERMO
        call EOSSetFromOptions_Thermo(eos, options, ierr)
        done = PETSC_TRUE
      end if
    end if

    if (.not.done) then
      test_eos_name = 'pr'
      if (StringCompareIgnoreCase(test_eos_name, eos_name, 2)) then
        eos%eos_type = EOS_PR
        call EOSSetFromOptions_PR(eos, options, ierr)
        done = PETSC_TRUE
      end if
    end if

    if (.not.done) then
      call LBMError(PETSC_COMM_WORLD, 1, 'Invalid EOS type', ierr)
    end  if
  end subroutine EOSSetFromOptions

  subroutine EOSApply(eos, rho, psi, g_mm, m, dist)
    type(eos_type) eos
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,dist%info%rgxyzl) :: rho, psi
    PetscScalar g_mm
    PetscInt m

    PetscErrorCode ierr

    select case(eos%eos_type)
    case (EOS_DENSITY)
      call EOSApply_Rho(eos, rho, psi, g_mm, m, dist)
    case (EOS_SC)
      call EOSApply_SC(eos, rho, psi, g_mm, m, dist)
    case (EOS_Thermo)
      call EOSApply_Thermo(eos, rho, psi, g_mm, m, dist)
    case (EOS_PR)
      call EOSApply_PR(eos, rho, psi, g_mm, m, dist)
    case DEFAULT
      call LBMError(PETSC_COMM_WORLD, 1, 'Invalid EOS type', ierr)
    end select
  end subroutine EOSApply

  ! private implementations

  ! EOS_DENSITY
  ! psi = rho
  subroutine EOSSetFromOptions_Rho(eos, options, ierr)
    use LBM_Options_module
    type(eos_type) eos
    type(options_type) options
    PetscErrorCode ierr
  end subroutine EOSSetFromOptions_Rho

  subroutine EOSApply_Rho(eos, rho, psi, g_mm, m, dist)
    type(eos_type) eos
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,dist%info%rgxyzl) :: rho, psi
    PetscScalar g_mm
    PetscInt m

    PetscInt i
    do i=1,dist%info%rgxyzl
      psi(m,i) = rho(m,i)
    end do
  end subroutine EOSApply_Rho

  ! EOS_SC
  ! EOS from Shan & Chen '93 & '94
  ! psi = 1 - exp(-rho)
  subroutine EOSSetFromOptions_SC(eos, options, ierr)
    use LBM_Options_module
    type(eos_type) eos
    type(options_type) options
    PetscErrorCode ierr

    PetscBool help, flag
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    ! set defaults
    eos%rho0 = 1.

    ! get options
    if (help) call PetscPrintf(options%comm, '  -eos_sc_rho0_'//trim(eos%name)//&
         '=1:'//' coefficient\n', ierr)
    call PetscOptionsGetReal(options%my_prefix, '-eos_sc_rho0_'//trim(eos%name), &
         eos%rho0, flag, ierr)
  end subroutine EOSSetFromOptions_SC

  subroutine EOSApply_SC(eos, rho, psi, g_mm, m, dist)
    type(eos_type) eos
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,dist%info%rgxyzl) :: rho, psi
    PetscScalar g_mm
    PetscInt m

    PetscInt i
    do i=1,dist%info%rgxyzl
      psi(m,i) = eos%rho0*(1. - EXP(-rho(m,i)/eos%rho0))
    end do
  end subroutine EOSApply_SC

  ! EOS_THERMO
  ! EOS from Shan & Chen '94
  ! psi = psi0*exp(-rho0/rho)
  subroutine EOSSetFromOptions_Thermo(eos, options, ierr)
    use LBM_Options_module
    type(eos_type) eos
    type(options_type) options
    PetscErrorCode ierr
    
    PetscBool help, flag
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    ! set defaults
    eos%rho0 = 1.
    eos%psi0 = 1.

    ! get options
    if (help) call PetscPrintf(options%comm, '  -eos_thermo_psi0_'//trim(eos%name)//&
         '=1:'//' coefficient\n', ierr)
    call PetscOptionsGetReal(options%my_prefix, '-eos_thermo_psi0_'//trim(eos%name), &
         eos%psi0, flag, ierr)
    if (help) call PetscPrintf(options%comm, '  -eos_thermo_rho0_'//trim(eos%name)//&
         '=1:'//' coefficient\n', ierr)
    call PetscOptionsGetReal(options%my_prefix, '-eos_thermo_rho0_'//trim(eos%name), &
         eos%rho0, flag, ierr)
  end subroutine EOSSetFromOptions_Thermo

  subroutine EOSApply_Thermo(eos, rho, psi, g_mm, m, dist)
    type(eos_type) eos
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,dist%info%rgxyzl) :: rho, psi
    PetscScalar g_mm
    PetscInt m

    PetscInt i
    do i=1,dist%info%rgxyzl
      psi(m,i) = eos%psi0*EXP(-eos%rho0/rho(m,i))
    end do
  end subroutine EOSApply_Thermo

  ! EOS_PR
  ! EOS from Peng & Robinson
  ! psi = ...
  subroutine EOSSetFromOptions_PR(eos, options, ierr)
    use LBM_Options_module
    type(eos_type) eos
    type(options_type) options
    PetscErrorCode ierr
    PetscBool help, flag
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    ! set defaults
    eos%a = 2.d0/49.
    eos%b = 2.d0/21.

    if (help) call PetscPrintf(options%comm, '  -eos_pr_a_'//trim(eos%name)//&
         '=2/49:'//' coefficient a\n', ierr)
    call PetscOptionsGetReal(options%my_prefix, '-eos_pr_a_'//trim(eos%name), &
         eos%a, flag, ierr)
    if (help) call PetscPrintf(options%comm, '  -eos_pr_b_'//trim(eos%name)//&
         '=2/49:'//' coefficient b\n', ierr)
    call PetscOptionsGetReal(options%my_prefix, '-eos_pr_b_'//trim(eos%name), &
         eos%b, flag, ierr)

    ! from these, calculate T_c and p_c
    eos%T_c = eos%a/eos%b*0.0778/0.45724/eos%R
    eos%p_c = 0.0778*eos%R*eos%T_c/eos%b

    ! get T, or default to T/T_c = 0.9
    eos%T = 0.9*eos%T_c
    if (help) call PetscPrintf(options%comm, '  -eos_pr_T_'//trim(eos%name)//&
         '=1.:'//' temperature (LU)\n', ierr)
    call PetscOptionsGetReal(options%my_prefix, '-eos_pr_T_'//trim(eos%name), &
         eos%T, flag, ierr)

    ! get omega, or default to water
    eos%omega = 0.344
    if (help) call PetscPrintf(options%comm, '  -eos_pr_omega_'//trim(eos%name)//&
         '=0.344:'//' coefficient\n', ierr)
    call PetscOptionsGetReal(options%my_prefix, '-eos_pr_omega_'//trim(eos%name), &
         eos%omega, flag, ierr)
  end subroutine EOSSetFromOptions_PR

  subroutine EOSApply_PR(eos, rho, psi, g_mm, m, dist)
    type(eos_type) eos
    type(distribution_type) dist
    PetscScalar,dimension(dist%s,dist%info%rgxyzl) :: rho, psi
    PetscScalar g_mm
    PetscInt m

    PetscInt i
    PetscScalar alpha   

    alpha = (1. + (0.37464 + 1.54226*eos%omega - 0.26992*eos%omega**2)*(1. - & 
            SQRT(eos%T/eos%T_c)))**2
    do i=1,dist%info%rgxyzl
      psi(m,i) = SQRT( 2.*( rho(m,i)*eos%R*eos%T/(1. - eos%b*rho(m,i)) - &
                 eos%a*alpha*rho(m,i)**2/(1. + 2.*eos%b*rho(m,i) - eos%b**2*rho(m,i)**2) &
                 - rho(m,i)/3. )/(dist%disc%c_0*g_mm ) )
    end do
  end subroutine EOSApply_PR
end module LBM_EOS_module
