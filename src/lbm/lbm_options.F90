!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        options.F90
!!!     version:
!!!     created:         09 December 2010
!!!       on:            14:16:32 MST
!!!     last modified:   15 November 2011
!!!       at:            16:31:11 MST
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
       PetscInt restart_counter
       PetscInt current_waypoint
       PetscInt, pointer:: waypoints(:)

       PetscBool ic_from_file
       character(len=MAXSTRINGLENGTH):: ic_file
       character(len=MAXSTRINGLENGTH):: output_prefix
       PetscBool mpiio
       PetscBool supress_ic_output
       PetscBool print_help

       PetscInt flow_disc
       PetscInt transport_disc

       PetscInt ndims
       PetscInt ncomponents
       PetscInt nspecies
       PetscInt nminerals

       PetscInt flow_relaxation_mode
       PetscBool flow_fluidsolid_forces
       PetscBool flow_use_nonideal_eos
       PetscInt deriv_order

       PetscBool steadystate
       PetscInt steadystate_rampup_steps
       character(len=MAXSTRINGLENGTH):: steadystate_flow_file
       PetscBool steadystate_hasfile

       PetscInt transport_relaxation_mode
       PetscBool transport_reactive_matrix
    end type options_type

    public :: OptionsCreate, &
         OptionsSetUp, &
         OptionsSetPrefix, &
         OptionsView, &
         OptionsGroupHeader, &
         OptionsGroupFooter, &
         OptionsGetBool, &
         OptionsGetInt, &
         OptionsGetIntArray, &
         OptionsGetReal, &
         OptionsGetRealArray, &
         OptionsGetString

  contains
    function OptionsCreate(comm) result(options)
      type(options_type),pointer:: options

      MPI_Comm comm

      allocate(options)
      options%comm = comm
      options%ntimes = 1
      options%npasses = 1
      options%kprint = 0
      options%kwrite = -1
      options%restart = PETSC_FALSE
      options%restart_counter = -1
      options%istep = 0
      options%mpiio = PETSC_FALSE
      options%supress_ic_output = PETSC_FALSE
      options%print_help = PETSC_FALSE
      options%ic_from_file = PETSC_FALSE
      options%ic_file = ''

      options%current_waypoint = 0
      nullify(options%waypoints)

      options%output_prefix = 'test_solution/'
      
      options%flow_disc = NULL_DISCRETIZATION
      options%transport_disc = NULL_DISCRETIZATION

      options%ndims = 0
      options%nminerals = 1
      options%ncomponents = 1
      options%nspecies = 0

      options%flow_relaxation_mode = RELAXATION_MODE_SRT
      options%flow_fluidsolid_forces = PETSC_FALSE
      options%flow_use_nonideal_eos = PETSC_FALSE
      options%deriv_order = 4

      options%steadystate = PETSC_FALSE
      options%steadystate_rampup_steps = 0
      options%steadystate_flow_file = ''
      options%steadystate_hasfile = PETSC_FALSE

      options%transport_relaxation_mode = RELAXATION_MODE_SRT
      options%transport_reactive_matrix = PETSC_FALSE
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
      PetscBool flag
      PetscErrorCode ierr
      PetscInt nmax

      character(len=MAXWORDLENGTH):: name
      character(len=MAXWORDLENGTH):: test_discretization
      PetscInt tmpdims
      PetscInt wpnum, wpdummy, lcv
      character(len=3):: wpstring
      PetscBool help

      call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", options%print_help, ierr)
      call OptionsGroupHeader(options, "Simulation Options", ierr)

      ! options for timestepping and steady state solutions
      call OptionsGroupHeader(options, " Timestepping Control Options", ierr)
      call OptionsGetInt(options, "-ntimes", "total ??? to run (defunct)", options%ntimes, &
           flag, ierr)
      call OptionsGetInt(options, "-npasses", "total timesteps to run", options%npasses, &
           flag, ierr)
      call OptionsGetBool(options, "-steadystate", "turn off flow, assuming steady state", &
           options%steadystate, flag, ierr)
      if (options%steadystate) then
        call OptionsGetString(options, "-steadystate_flow_file", &
             "set flow via pre-computed steadystate", options%steadystate_flow_file, &
             flag, ierr)
        if (.not.flag) then
          call OptionsGetInt(options, "-steadystate_rampup_steps", &
               "allow flow to ramp up if flow not set via file", &
               options%steadystate_rampup_steps, flag, ierr)
        end if
      end if

      ! options for initial condition
      call OptionsGroupHeader(options, " IC Options", ierr)
      call OptionsGetBool(options, "-restart", "restart from an old simulation", &
           options%restart, flag, ierr)
      call OptionsGetInt(options, "-restart_counter", "output file number to start from", &
           options%restart_counter, flag, ierr)
      call OptionsGetString(options, "-ic_file", &
           "start from a given IC for the distribution function", &
           options%ic_file, options%ic_from_file, ierr)
      call OptionsGetInt(options, "-istep", "intial timestep", options%istep, flag, ierr)

      ! options for i/o
      call OptionsGroupHeader(options, " I/O Options", ierr)
      call OptionsGetInt(options, "-kwrite", "output interval in timesteps", &
           options%kwrite, flag, ierr)
      call OptionsGetString(options, "-output_file_prefix", &
           "prefix for path of output files", options%output_prefix, flag, ierr)
      call OptionsGetBool(options, "-mpiio", "use mpiio for parallel i/o", &
           options%mpiio, flag, ierr)
      call OptionsGetBool(options, "-supress_ic_output", "do not output IC", &
           options%supress_ic_output, flag, ierr)

      ! waypoints and checkpointing, i/o stuff
      wpnum = 0
      wpstring = ''
      do
        wpdummy = -1
         wpnum = wpnum + 1
         if (wpnum > 99) then
            write(wpstring, '(I3)') wpnum
         else if (wpnum > 9) then
            write(wpstring, '(I2)') wpnum
         else
            write(wpstring, '(I1)') wpnum
         end if
         call OptionsGetInt(options,'-waypoint'//TRIM(wpstring), "checkpoint istep", &
              wpdummy, flag, ierr)
         if (.not.flag) exit
      end do

      if (wpnum > 1) then
         allocate(options%waypoints(wpnum))
         options%waypoints = -1
         do lcv=1,wpnum-1
            if (lcv > 99) then
               write(wpstring, '(I3)') lcv
            else if (lcv > 9) then
               write(wpstring, '(I2)') lcv
            else
               write(wpstring, '(I1)') lcv
            end if
            call PetscOptionsGetInt(options%my_prefix,'-waypoint'//TRIM(wpstring), &
                 options%waypoints(lcv), flag, ierr)
         end do
         options%current_waypoint = 1
      end if

      ! --- options that shouldn't be here and really need to be moved  ---
      ! --- into their respective modules                               ---
      ! flow stuff
      call OptionsGroupHeader(options, " Flow Options", ierr)
      call OptionsGetInt(options, "-flow_relaxation_mode", &
           "flow relaxation as SRT=0, MRT=1", options%flow_relaxation_mode, flag, ierr)
      call OptionsGetBool(options, "-flow_use_nonideal_eos", &
           "use a nonideal eos equation, set by phase", options%flow_use_nonideal_eos, &
           flag, ierr)

      call OptionsGetInt(options, "-nminerals", "number of minerals", options%nminerals, &
           flag,ierr)
      call OptionsGetInt(options, "-ncomponents", "number of components", &
           options%ncomponents,flag,ierr)
      call OptionsGetInt(options, "-derivative_order", &
           "order of fluid-fluid term derivatives", options%deriv_order,flag,ierr)

      ! flow discretization
      name = "D3Q19"
      call OptionsGetString(options, "-flow_discretization", "flow discretization type", &
           name, flag, ierr)
      if (.not.flag) then
        call OptionsGetString(options, "-discretization", "flow discretization type", &
             name, flag, ierr)
      end if

      print*, 'ndims:', options%ndims
      print*, 'disc:', name
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

      ! transport stuff
      call OptionsGroupHeader(options, " Transport Options", ierr)
      name = "D3Q19"
      call OptionsGetString(options, "-transport_discretization", &
           "transport discretization type", name, flag, ierr)

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
            SETERRQ(1,1,"Flow and transport discretization dimensions don't match", ierr)
        end if

        options%nspecies = 1
        call OptionsGetInt(options, "-nspecies", "number of species", options%nspecies, &
             flag, ierr)
        call OptionsGetInt(options, "-transport_relaxation_mode", &
             "transport relaxation as SRT=0, MRT=1", options%transport_relaxation_mode, &
             flag, ierr)
        call OptionsGetBool(options, "-reactive_matrix", &
             "allow dissolution/precipitation", options%transport_reactive_matrix, &
             flag, ierr)
      end if

      call OptionsGroupFooter(options, "Simulation Options", ierr)
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

    subroutine OptionsGroupHeader(options, group, ierr)
      type(options_type) options
      character(len=*):: group
      PetscErrorCode ierr

      if (options%print_help) then
         call PetscPrintf(options%comm, trim(group), ierr)
         call PetscPrintf(options%comm, " -------------------------\n", ierr)
      end if
    end subroutine OptionsGroupHeader

    subroutine OptionsGroupFooter(options, group, ierr)
      type(options_type) options
      character(len=*):: group
      PetscErrorCode ierr

      if (options%print_help) then
        call PetscPrintf(options%comm, "------------------------- ", ierr)
        call PetscPrintf(options%comm, trim(group)//" \n", ierr)
      end if
    end subroutine OptionsGroupFooter

    subroutine OptionsGetBool(options, name, help, val, flag, ierr)
      type(options_type) options
      character(len=*):: name
      character(len=*):: help
      PetscBool val
      PetscBool flag
      PetscErrorCode ierr

      character(len=MAXWORDLENGTH):: default = "FALSE"

      if (options%print_help) then
        if (val) default = "TRUE"
        call PetscPrintf(options%comm, "  -"//trim(options%my_prefix)// &
             trim(name(2:))//" <"//trim(default)//">: "//trim(help)//"\n", ierr)
      end if
      call PetscOptionsGetBool(options%my_prefix, name, val, flag, ierr)
    end subroutine OptionsGetBool

    subroutine OptionsGetInt(options, name, help, val, flag, ierr)
      type(options_type) options
      character(len=*):: name
      character(len=*):: help
      PetscInt val
      PetscBool flag
      PetscErrorCode ierr

      character(len=MAXWORDLENGTH):: default = ""

      if (options%print_help) then
        write(default,"(I0)") val
        call PetscPrintf(options%comm, "  -"//trim(options%my_prefix)// &
             trim(name(2:))//" <"//trim(default)//">: "//trim(help)//"\n", ierr)
      end if
      call PetscOptionsGetInt(options%my_prefix, name, val, flag, ierr)
    end subroutine OptionsGetInt

    subroutine OptionsGetIntArray(options, name, help, val, nvals, flag, ierr)
      type(options_type) options
      character(len=*):: name
      character(len=*):: help
      PetscInt nvals
      PetscInt val(nvals)
      PetscBool flag
      PetscErrorCode ierr

      PetscInt lcv

      character(len=MAXWORDLENGTH):: default = ""
      character(len=MAXWORDLENGTH):: tmpdefault = ""

      if (options%print_help) then
        do lcv=1,nvals
          write(tmpdefault,"(I0),") val(lcv)
          default = trim(default)//tmpdefault
        end do

        call PetscPrintf(options%comm, "  -"//trim(options%my_prefix)// &
             trim(name(2:))//" <"//default(1:len_trim(default)-1)//">: "//trim(help)//"\n", ierr)
      end if
      call PetscOptionsGetIntArray(options%my_prefix, name, val, nvals, flag, ierr)
    end subroutine OptionsGetIntArray

    subroutine OptionsGetReal(options, name, help, val, flag, ierr)
      type(options_type) options
      character(len=*):: name
      character(len=*):: help
      PetscReal val
      PetscBool flag
      PetscErrorCode ierr

      character(len=MAXWORDLENGTH):: default = ""

      if (options%print_help) then
        write(default,"(D1.1)") val
        call PetscPrintf(options%comm, "  -"//trim(options%my_prefix)// &
             trim(name(2:))//" <"//trim(default)//">: "//trim(help)//"\n", ierr)
      end if
      call PetscOptionsGetReal(options%my_prefix, name, val, flag, ierr)
    end subroutine OptionsGetReal

    subroutine OptionsGetRealArray(options, name, help, val, nvals, flag, ierr)
      type(options_type) options
      character(len=*):: name
      character(len=*):: help
      PetscInt nvals
      PetscReal val(nvals)
      PetscBool flag
      PetscErrorCode ierr

      PetscInt lcv

      character(len=MAXWORDLENGTH):: default = ""
      character(len=MAXWORDLENGTH):: tmpdefault = ""

      if (options%print_help) then
        do lcv=1,nvals
          write(tmpdefault,"(D1.1),") val(lcv)
          default = trim(default)//tmpdefault
        end do

        call PetscPrintf(options%comm, "  -"//trim(options%my_prefix)// &
             trim(name(2:))//" <"//default(1:len_trim(default)-1)//">: "//trim(help)//"\n", ierr)
      end if
      call PetscOptionsGetRealArray(options%my_prefix, name, val, nvals, flag, ierr)
    end subroutine OptionsGetRealArray

    subroutine OptionsGetString(options, name, help, val, flag, ierr)
      type(options_type) options
      character(len=*):: name
      character(len=*):: help
      character(len=*):: val

      PetscBool flag
      PetscErrorCode ierr

      if (options%print_help) then
        call PetscPrintf(options%comm, "  -"//trim(options%my_prefix)// &
             trim(name(2:))//" <"//trim(val)//">: "//trim(help)//"\n", ierr)
      end if
      call PetscOptionsGetString(options%my_prefix, name, val, flag, ierr)
    end subroutine OptionsGetString

  end module LBM_Options_module
