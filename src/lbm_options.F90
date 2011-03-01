!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        options.F90
!!!     version:         
!!!     created:         09 December 2010
!!!       on:            14:16:32 MST
!!!     last modified:   15 February 2011
!!!       at:            11:45:13 MST
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
       PetscBool new_simulation
       PetscInt istep

       character(len=MAXSTRINGLENGTH):: output_prefix
       character(len=MAXSTRINGLENGTH):: walls_file
       PetscInt walls_type
    end type options_type

    public :: OptionsCreate, &
         OptionsInitialize, &
         OptionsView

  contains
    function OptionsCreate(comm) result(options)
      implicit none
      type(options_type),pointer:: options

      MPI_Comm comm

      allocate(options)
      options%comm = comm
      options%ntimes = 0
      options%npasses = 0
      options%kprint = 0
      options%kwrite = 0
      options%new_simulation = .TRUE.
      options%istep = 0

      options%output_prefix = 'test_solution/'
      options%walls_file = 'geometry.dat'
      options%walls_type = WALLS_TYPE_PETSC

    end function OptionsCreate

    subroutine OptionsInitialize(options, prefix, ierr)
      implicit none
      type(options_type) options
      character(len=MAXWORDLENGTH):: prefix
      PetscBool flag
      PetscErrorCode ierr
      PetscInt nmax
      integer charlen
      

      charlen = LEN_TRIM(prefix)
      options%my_prefix = prefix(1:charlen)
     
      ! update defaults from input file
      call PetscOptionsGetInt(options%my_prefix, '-ntimes', options%ntimes, flag, ierr)
      call PetscOptionsGetInt(options%my_prefix, '-npasses', options%npasses, flag, ierr)
      call PetscOptionsGetInt(options%my_prefix, '-kprint', options%kprint, flag, ierr)
      call PetscOptionsGetInt(options%my_prefix, '-kwrite', options%kwrite, flag, ierr)
      call PetscOptionsGetBool(options%my_prefix, '-new_simulation', options%new_simulation, flag, ierr)
      call PetscOptionsGetInt(options%my_prefix, '-istep', options%istep, flag, ierr)

      call PetscOptionsGetString(options%my_prefix, '-output_file_prefix', options%output_prefix, flag, ierr)
      call PetscOptionsGetString(options%my_prefix, '-walls_file', options%walls_file, flag, ierr)
      call PetscOptionsGetInt(options%my_prefix, '-walls_type', options%walls_type, flag, ierr)
      return
    end subroutine OptionsInitialize

    subroutine OptionsView(options)
      implicit none
      type(options_type) options

      print*, 'Options Used:'
      print*, ' Timestepping:'
      print*, '  ntimes =',options%ntimes
      print*, '  npasses =',options%npasses
      print*, '  i/o interval =', options%kprint,options%kwrite
    end subroutine OptionsView

  end module LBM_Options_module
