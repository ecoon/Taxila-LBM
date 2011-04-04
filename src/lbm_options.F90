!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        options.F90
!!!     version:         
!!!     created:         09 December 2010
!!!       on:            14:16:32 MST
!!!     last modified:   04 April 2011
!!!       at:            12:18:29 MDT
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
    implicit none

    private
#include "lbm_definitions.h"

    type,public:: options_type
       MPI_Comm comm
       character(len=MAXWORDLENGTH):: my_prefix
       PetscInt ntimes, npasses
       PetscInt kprint, kwrite
       PetscBool restart
       PetscInt istep

       character(len=MAXSTRINGLENGTH):: output_prefix
       character(len=MAXSTRINGLENGTH):: walls_file
       PetscBool mpiio
       PetscInt walls_type

       PetscInt flow_disc
       PetscInt transport_disc
       PetscInt ndims
       PetscInt nphases
       PetscInt ncomponents
       PetscInt flow_relaxation_mode
       PetscInt transport_relaxation_mode
    end type options_type

    public :: OptionsCreate, &
         OptionsSetUp, &
         OptionsSetPrefix, &
         OptionsView

  contains
    function OptionsCreate(comm) result(options)
      type(options_type),pointer:: options

      MPI_Comm comm

      allocate(options)
      options%comm = comm
      options%ntimes = 1
      options%npasses = 1
      options%kprint = 0
      options%kwrite = 1
      options%restart = PETSC_FALSE
      options%istep = 0
      options%mpiio = PETSC_FALSE

      options%output_prefix = 'test_solution/'
      options%walls_file = 'geometry.dat'
      options%walls_type = WALLS_TYPE_PETSC
      
      options%flow_disc = NULL_DISCRETIZATION
      options%transport_disc = NULL_DISCRETIZATION
      options%ndims = 0
      options%nphases = 1
      options%ncomponents = 1
      options%flow_relaxation_mode = RELAXATION_MODE_SRT
      options%transport_relaxation_mode = RELAXATION_MODE_SRT
    end function OptionsCreate

    subroutine OptionsSetPrefix(options, prefix)
      type(options_type) options
      character(len=MAXWORDLENGTH):: prefix
      integer charlen

      charlen = LEN_TRIM(prefix)
      options%my_prefix = prefix(1:charlen)
    end subroutine OptionsSetPrefix

    subroutine OptionsSetUp(options)
      use String_module

      type(options_type) options
      PetscBool help
      PetscBool flag
      PetscErrorCode ierr
      PetscInt nmax
      
      character(len=MAXWORDLENGTH):: name
      character(len=MAXWORDLENGTH):: test_discretization
      PetscInt tmpdims

      call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)
    
      ! update defaults from input file
      if (help) call PetscPrintf(options%comm, "  -ntimes=<1>: total ??? to run (defunct)\n", ierr)
      call PetscOptionsGetInt(options%my_prefix, '-ntimes', options%ntimes, flag, ierr)

      if (help) call PetscPrintf(options%comm, "  -npasses=<1>: total timesteps to run\n", ierr)
      call PetscOptionsGetInt(options%my_prefix, '-npasses', options%npasses, flag, ierr)

      if (help) call PetscPrintf(options%comm, "  -kwrite=<1>: output interval in timesteps\n", ierr)
      call PetscOptionsGetInt(options%my_prefix, '-kwrite', options%kwrite, flag, ierr)

      if (help) call PetscPrintf(options%comm, "  -restart: start from an old simulation\n", ierr)
      call PetscOptionsGetBool(options%my_prefix, '-restart', options%restart, flag, ierr)

      if (help) call PetscPrintf(options%comm, "  -istep=<0>: initial timestep\n", ierr)
      call PetscOptionsGetInt(options%my_prefix, '-istep', options%istep, flag, ierr)

      if (help) call PetscPrintf(options%comm, "  -output_file_prefix=<test_solution/>: "// &
           "prefix for solution data output files\n", ierr)
      call PetscOptionsGetString(options%my_prefix, '-output_file_prefix', options%output_prefix, flag, ierr)

      if (help) then
         call PetscPrintf(options%comm, "  -walls_file=<geometry.dat>: filename for porescale walls \n", ierr)
      end if
      call PetscOptionsGetString(options%my_prefix, '-walls_file', options%walls_file, flag, ierr)

      if (help) call PetscPrintf(options%comm, "  -mpiio: use mpiio for i/o\n", ierr)
      call PetscOptionsGetBool(options%my_prefix, '-mpiio', options%mpiio, flag, ierr)

      if (help) call PetscPrintf(options%comm, &
           "  -walls_type <1>: number of phases\n", ierr)
      call PetscOptionsGetInt(options%my_prefix, '-walls_type', options%walls_type, flag, ierr)

      if (help) call PetscPrintf(options%comm, &
           "  -flow_relaxation_mode <0>: flow relaxation as SRT=0, MRT=1\n", ierr)
      call PetscOptionsGetInt(options%my_prefix, '-flow_relaxation_mode', &
           options%flow_relaxation_mode, flag, ierr)

      ! set the flow discretization
      if (help) call PetscPrintf(options%comm, &
           "  -nphases <1>: number of phases\n", ierr)
      call PetscOptionsGetInt(options%my_prefix,'-nphases', options%nphases,flag,ierr)

      if (help) call PetscPrintf(options%comm, &
           "  -flow_discretization 'd3q19': discretization type\n", ierr)
      call PetscOptionsGetString(options%my_prefix, '-flow_discretization', &
           name, flag, ierr)
      if (.not.flag) then
         call PetscOptionsGetString(options%my_prefix, '-discretization', &
              name, flag, ierr)
      end if
      if (.not.flag) name = 'D3Q19'

      test_discretization = 'd3q19'
      if (StringCompareIgnoreCase(name, test_discretization, 6)) then
         options%flow_disc = D3Q19_DISCRETIZATION
         options%ndims = 3
      end if

      test_discretization = 'd2q9'
      if (StringCompareIgnoreCase(name, test_discretization, 5)) then
         options%flow_disc = D2Q9_DISCRETIZATION
         options%ndims = 2
      end if
      
      if (options%flow_disc == NULL_DISCRETIZATION) then
         SETERRQ(1, 1, 'Invalid Discretization', ierr)
      end if
         
      ! set the tran discretization
      options%ncomponents = 1
      if (help) call PetscPrintf(options%comm, &
           "  -ncomponents <1>: number of major components\n", ierr)
      call PetscOptionsGetInt(options%my_prefix,'-ncomponents', options%ncomponents, &
           flag,ierr)

      if (help) call PetscPrintf(options%comm, &
           "  -transport_relaxation_mode <0>: transport relaxation as SRT=0, MRT=1\n", ierr)
      call PetscOptionsGetInt(options%my_prefix, '-transport_relaxation_mode', &
           options%transport_relaxation_mode, flag, ierr)

      if (help) call PetscPrintf(options%comm, &
           "  -transport_discretization 'd3q19': discretization type\n", ierr)
      call PetscOptionsGetString(options%my_prefix, '-transport_discretization', &
           name, flag, ierr)

      if (flag) then
         test_discretization = 'd3q19'
         if (StringCompareIgnoreCase(name, test_discretization, 6)) then
            options%transport_disc = D3Q19_DISCRETIZATION
            tmpdims = 3
         end if

         test_discretization = 'd2q9'
         if (StringCompareIgnoreCase(name, test_discretization, 5)) then
            options%transport_disc = D2Q9_DISCRETIZATION
            tmpdims = 2
         end if

         if (options%transport_disc == NULL_DISCRETIZATION) then
            SETERRQ(1, 1, 'Invalid Discretization', ierr)
         else if (tmpdims /= options%ndims) then
            SETERRQ(1,1,"Discretization dimensions don't match", ierr)
         end if
      end if
      return
    end subroutine OptionsSetUp

    subroutine OptionsView(options)
      type(options_type) options

      print*, 'Options Used:'
      print*, ' Timestepping:'
      print*, '  ntimes =',options%ntimes
      print*, '  npasses =',options%npasses
      print*, '  i/o interval =', options%kprint,options%kwrite
    end subroutine OptionsView

  end module LBM_Options_module
