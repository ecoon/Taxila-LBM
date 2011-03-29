!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        main.F90
!!!     version:
!!!     created:         08 December 2010
!!!       on:            11:48:19 MST
!!!     last modified:   28 March 2011
!!!       at:            13:25:27 MDT
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
    type(lbm_type),pointer:: user
    type(options_type),pointer:: options


    ! --- setup environment
    call getarg(1, infile)
    call PetscInitialize(infile, ierr)
    user => LBMCreate(PETSC_COMM_WORLD)
    options => user%options

    ! initialize options and constants
    prefix = ''
    call OptionsSetPrefix(options, prefix)
    call OptionsSetUp(options)
    call LBMSetFromOptions(user, options, ierr)

!    if (user%info%id.eq.0) call LBMView(user)

    ! --- initialize state
    ! walls
    if(user%info%id.eq.0) print*,'initialization of walls'
    if (user%options%walls_type.eq.WALLS_TYPE_PETSC) then
       call LBMInitializeWallsPetsc(user, options%walls_file)
    else
       call LBMInitializeWalls(user, initialize_walls)
    end if

    ! bcs
    call BCSetValues(user%bc, user%info, initialize_bcs)

    ! fi/state
    if (options%new_simulation) then
       call LBMInitializeState(user, initialize_state)
       istep=0
    else
       call LBMInitializeStateRestarted(user, options%istep, options%kwrite)
       istep = options%istep
    endif

    ! start lbm
    if(user%info%id.eq.0) print*,'calling lbm from inital step', istep, 'to final step', options%ntimes*options%npasses

    call LBMRun(user, istep, options%ntimes*options%npasses, options%kwrite)
    call LBMDestroy(user, ierr)
    call PetscFinalize(ierr)
    stop
  end program main
  !----------------------------------------------------------

