!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        main.F90
!!!     version:
!!!     created:         08 December 2010
!!!       on:            11:48:19 MST
!!!     last modified:   02 May 2011
!!!       at:            11:07:20 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

program MAIN
  use petsc
  use LBM_Options_module
  use LBM_BC_module
  use LBM_module
  implicit none
#include "lbm_definitions.h"
  
  PetscInt istep
  PetscInt ntimes, npasses
  PetscInt kwrite, kprint
  PetscErrorCode ierr
  character(len=MAXSTRINGLENGTH) infile
  character(len=MAXWORDLENGTH) prefix
  
  external initialize_bcs
  external initialize_bcs_transport
  external initialize_state
  external initialize_state_transport
  external initialize_walls
  type(lbm_type),pointer:: lbm
  type(options_type),pointer:: options
  
  
  ! --- setup environment
  call getarg(1, infile)
  call PetscInitialize(infile, ierr)
  lbm => LBMCreate(PETSC_COMM_WORLD)
  options => lbm%options
  
  ! initialize options and constants
  prefix = ''
  call OptionsSetPrefix(options, prefix)
  call OptionsSetUp(options)
  call LBMSetFromOptions(lbm, options, ierr)
  call LBMSetUp(lbm)

  ! --- initialize state
  ! walls
  if (lbm%grid%info%rank.eq.0) then
     write(*,*) 'initialization of walls'
  end if
  call LBMInitializeWalls(lbm)
  
  ! set boundary conditions
  call BCSetValues(lbm%flow%bc, lbm%flow%distribution, options, initialize_bcs)
  if (associated(lbm%transport)) then
     call BCSetValues(lbm%transport%bc, lbm%transport%distribution, &
          options, initialize_bcs_transport)
  end if

  ! set initial conditions
  if (options%restart) then
     call LBMInitializeStateRestarted(lbm, options%istep, options%kwrite)
     istep = options%istep
  else
     if (associated(lbm%transport)) then
        call LBMInitializeState(lbm, initialize_state, initialize_state_transport)
     else
        call LBMInitializeState(lbm, initialize_state)
     end if
     
     if (options%steadystate_hasfile) then
        call LBMLoadSteadyStateFlow(lbm, options%steadystate_flow_file)
     end if
     istep=0
  endif
  
  ! start lbm
  if (lbm%grid%info%rank.eq.0) then
     write(*,*) 'calling lbm from inital step', istep, 'to final step', &
          options%ntimes*options%npasses
  end if
  
  call LBMInit(lbm, istep)
  call LBMRun(lbm, istep, options%ntimes*options%npasses, options%kwrite)
  call LBMDestroy(lbm, ierr)
  call PetscFinalize(ierr)
  stop
end program main
  !----------------------------------------------------------

