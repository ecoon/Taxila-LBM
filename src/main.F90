!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        main.F90
!!!     version:
!!!     created:         08 December 2010
!!!       on:            11:48:19 MST
!!!     last modified:   29 March 2011
!!!       at:            17:46:41 MDT
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
  external initialize_state
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
  
  ! --- initialize state
  ! walls
  if (lbm%grid%info%rank.eq.0) then
     write(*,*) 'initialization of walls'
  end if
  if (lbm%options%walls_type.eq.WALLS_TYPE_PETSC) then
     call LBMInitializeWallsPetsc(lbm, options%walls_file)
  else
     call LBMInitializeWalls(lbm, initialize_walls)
  end if
  
  ! bcs
  call BCSetValues(lbm%bc, lbm%flow%distribution, options, initialize_bcs)
  
  ! fi/state
  if (options%new_simulation) then
     call LBMInitializeState(lbm, initialize_state)
     istep=0
  else
     call LBMInitializeStateRestarted(lbm, options%istep, options%kwrite)
     istep = options%istep
  endif
  
  ! start lbm
  if (lbm%grid%info%rank.eq.0) then
     write(*,*) 'calling lbm from inital step', istep, 'to final step', &
          options%ntimes*options%npasses
  end if
  
  call LBMRun(lbm, istep, options%ntimes*options%npasses, options%kwrite)
  call LBMDestroy(lbm, ierr)
  call PetscFinalize(ierr)
  stop
end program main
  !----------------------------------------------------------

