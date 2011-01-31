!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        options.F90
!!!     version:         
!!!     created:         09 December 2010
!!!       on:            14:16:32 MST
!!!     last modified:   31 January 2011
!!!       at:            09:53:36 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"
  module LBM_Options_module
    use petsc
    use String_module
    implicit none

    private
#include "lbm_definitions.h"

    type,public:: options_type
       character(len=MAXWORDLENGTH):: my_prefix
       PetscInt NX, NY, NZ
       PetscInt s
       PetscInt discretization
       PetscScalar g,g11,g22          ! constants for mobility forces
       PetscScalar,pointer,dimension(:)::tau  ! relaxation times
       PetscScalar,pointer,dimension(:)::gvt  ! gravity (i.e. body forces)
       PetscScalar,pointer,dimension(:)::gw   ! fluid-solid interaction forces
       PetscScalar,pointer,dimension(:)::mm   ! mass per number (rho = mm*n)
       PetscScalar,pointer,dimension(:):: rho1, rho2         ! left and right fluid densiti
       PetscInt ntimes, npasses
       PetscInt kprint, kwrite
       PetscBool new_simulation
       PetscInt istep
       PetscInt,dimension(1:6):: bc_flags

       character(len=MAXSTRINGLENGTH):: output_prefix
       character(len=MAXSTRINGLENGTH):: walls_file
       character(len=MAXSTRINGLENGTH):: old_options_file
       PetscInt walls_type
       PetscBool:: use_old_options_style
    end type options_type

    public :: OptionsCreate, &
         OptionsSetSizes, &
         OptionsInitialize, &
         OptionsReadOldFile, &
         OptionsPrint

  contains
    function OptionsCreate() result(options)
      implicit none
      type(options_type),pointer:: options

      allocate(options)
      options%NX = -1
      options%NY = -1
      options%NZ = -1
      options%discretization = NULL_DISCRETIZATION
      options%g = 0
      options%g11 = 0
      options%g22 = 0
      options%ntimes = 0
      options%npasses = 0
      options%kprint = 0
      options%kwrite = 0
      options%new_simulation = .TRUE.
      options%istep = 0
      options%bc_flags = 0

      nullify(options%tau)
      nullify(options%gvt)
      nullify(options%gw)
      nullify(options%mm)
      nullify(options%rho1)
      nullify(options%rho2)

      options%output_prefix = 'test_solution/'
      options%walls_file = 'geometry.dat'
      options%walls_type = WALLS_TYPE_PETSC
      options%old_options_file = ''
      options%use_old_options_style = .FALSE.

    end function OptionsCreate

    subroutine OptionsSetSizes(options, s)
      implicit none
      type(options_type) options
      PetscInt s

      if (.not.associated(options%tau)) allocate(options%tau(s))
      options%tau = 0

      if (.not.associated(options%gvt)) allocate(options%gvt(s))
      options%gvt = 0

      if (.not.associated(options%gw)) allocate(options%gw(s))
      options%gw = 0

      if (.not.associated(options%mm)) allocate(options%mm(s))
      options%mm = 0

      if (.not.associated(options%rho1)) allocate(options%rho1(s))
      options%rho1 = 0
      if (.not.associated(options%rho2)) allocate(options%rho2(s))
      options%rho2 = 0
    end subroutine OptionsSetSizes

    subroutine OptionsInitialize(options, prefix, ierr)
      implicit none
      type(options_type) options
      character(len=MAXWORDLENGTH):: prefix
      character(len=MAXWORDLENGTH):: discretization
      PetscBool flag
      PetscErrorCode ierr
      PetscInt nmax
      integer charlen
      

      charlen = LEN_TRIM(prefix)
      options%my_prefix = prefix(1:charlen)
     
      ! grab sizes
      call PetscOptionsGetInt(options%my_prefix,'-s',options%s,flag,ierr)
      if (.not.flag) options%s = 2

      ! allocate space
      call OptionsSetSizes(options, options%s)

      ! check for old style
      call PetscOptionsGetBool(options%my_prefix,'-use_old_options_style', &
           options%use_old_options_style,flag,ierr)
      if (flag) then
         call PetscOptionsGetString(options%my_prefix,'-old_options_file', &
              options%old_options_file,flag,ierr)
         if (.not.flag) options%old_options_file = 'input_data'
         call OptionsReadOldFile(options, options%old_options_file)
      end if

      ! update sizes
      call PetscOptionsGetInt(options%my_prefix,'-nx',options%NX,flag,ierr)
      call PetscOptionsGetInt(options%my_prefix,'-ny',options%NY,flag,ierr)
      call PetscOptionsGetInt(options%my_prefix,'-nz',options%NZ,flag,ierr)

      call PetscOptionsGetString(options%my_prefix,'-discretization',discretization,flag,ierr)
      if (.not.flag) then
         options%discretization = D3Q19
      else
         if (StringCompareIgnoreCase(discretization, 'd3q19', 6)) then
            options%discretization = D3Q19
         else if (StringCompareIgnoreCase(discretization, 'd2q9', 5)) then
            options%discretization = D2Q9
         else
            write(*,*) 'INVALID DISCRETIZATION:', discretization
            ierr = 1
            CHKERRQ(ierr)
         end if
      end if
      
      ! update defaults from input file
      call PetscOptionsGetReal(options%my_prefix, '-g', options%g, flag, ierr)
      call PetscOptionsGetReal(options%my_prefix, '-g11', options%g, flag, ierr)
      call PetscOptionsGetReal(options%my_prefix, '-g22', options%g, flag, ierr)
      call PetscOptionsGetInt(options%my_prefix, '-ntimes', options%ntimes, flag, ierr)
      call PetscOptionsGetInt(options%my_prefix, '-npasses', options%npasses, flag, ierr)
      call PetscOptionsGetInt(options%my_prefix, '-kprint', options%kprint, flag, ierr)
      call PetscOptionsGetInt(options%my_prefix, '-kwrite', options%kwrite, flag, ierr)
      call PetscOptionsGetBool(options%my_prefix, '-new_simulation', options%new_simulation, flag, ierr)
      call PetscOptionsGetInt(options%my_prefix, '-istep', options%istep, flag, ierr)

      nmax = 6
      call PetscOptionsGetIntArray(options%my_prefix, '-bc_flags', options%bc_flags, nmax, flag, ierr)

      nmax = options%s
      call PetscOptionsGetRealArray(options%my_prefix, '-tau', options%tau, nmax, flag, ierr)
      nmax = options%s
      call PetscOptionsGetRealArray(options%my_prefix, '-gvt', options%gvt, nmax, flag, ierr)
      nmax = options%s
      call PetscOptionsGetRealArray(options%my_prefix, '-gw', options%gw, nmax, flag, ierr)
      nmax = options%s
      call PetscOptionsGetRealArray(options%my_prefix, '-mm', options%mm, nmax, flag, ierr)
      nmax = options%s
      call PetscOptionsGetRealArray(options%my_prefix, '-rho1', options%rho1, nmax, flag, ierr)
      nmax = options%s
      call PetscOptionsGetRealArray(options%my_prefix, '-rho2', options%rho2, nmax, flag, ierr)

      call PetscOptionsGetString(options%my_prefix, '-output_file_prefix', options%output_prefix, flag, ierr)
      call PetscOptionsGetString(options%my_prefix, '-walls_file', options%walls_file, flag, ierr)
      call PetscOptionsGetInt(options%my_prefix, '-walls_type', options%walls_type, flag, ierr)
      return
    end subroutine OptionsInitialize

    subroutine OptionsReadOldFile(options, filename)
      implicit none
      type(options_type) options
      character(len=MAXSTRINGLENGTH):: filename

      open(90,file=filename,status='unknown')
      read(90,*) options%NX,options%NY,options%NZ
      read(90,*) options%g, options%g11, options%g22, options%gw(1), options%gw(2)
      read(90,*) options%tau(1),options%tau(2)
      read(90,*) options%gvt(1),options%gvt(2)
      read(90,*) options%mm(1),options%mm(2)
      read(90,*) options%rho1(1),options%rho1(2)
      read(90,*) options%rho2(1),options%rho2(2)
      read(90,*) options%ntimes
      read(90,*) options%npasses
      read(90,*) options%new_simulation
      read(90,*) options%istep          ! istep=0 if new=.true.
      read(90,*) options%kprint,options%kwrite
      read(90,*) options%bc_flags(1),options%bc_flags(2),options%bc_flags(3), &
           options%bc_flags(4),options%bc_flags(5),options%bc_flags(6)
      close(90)
    end subroutine OptionsReadOldFile

    subroutine OptionsPrint(options)
      implicit none
      type(options_type) options

      print*, 'Options Used:'
      print*, ' Domain size =', options%NX,options%NY,options%NZ
      print*, ' Timestepping:'
      print*, '  ntimes =',options%ntimes
      print*, '  npasses =',options%npasses
      print*, '  i/o interval =', options%kprint,options%kwrite
      print*, ' Physics:'
      print*, '  g =', options%g
      print*, '  g11,g22 =', options%g11, options%g22
      print*, '  gw =', options%gw
      print*, '  tau =', options%tau
      print*, '  gvt =', options%gvt
      print*, '  mm =', options%mm
      print*, '  rho1 =', options%rho1
      print*, '  rho2 =', options%rho2
      print*, ' BC Enums:'
      print*, options%bc_flags
    end subroutine OptionsPrint

  end module LBM_Options_module
