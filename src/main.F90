!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        main.F90
!!!     version:
!!!     created:         08 December 2010
!!!       on:            11:48:19 MST
!!!     last modified:   16 December 2010
!!!       at:            16:21:02 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

  program main
    use BC_module
    use LBM_module
    use Options_module
    use constants
    use petsc
    implicit none

    PetscInt istep
    PetscInt ntimes, npasses
    PetscInt kwrite, kprint
    PetscErrorCode ierr
    character(60) infile
    character*20 prefix

    external initialize_bcs
    external initialize_state
    external initialize_walls
    type(lbm_type),pointer:: user
    type(options_type),pointer:: options

#include "lbm_definitions.h"

    ! --- setup environment
    call getarg(1, infile)
    call PetscInitialize(infile, ierr)
    user => LBMCreate(PETSC_COMM_WORLD)
    options => user%options

    ! initialize options and constants
    prefix = ''
    call OptionsInitialize(options, prefix, ierr)
    call constants_initialize(options%s)
    call constants_set_from_options(user%options)

    ! --- initialize memory
    call LBMSetSizes(user, options%NX, options%NY, options%NZ, options%s, options%b)
    user%bc%flags(:) = options%bc_flags(:)

    if (user%info%id.eq.0) call OptionsPrint(options)
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
       istep=1
    else
       call initialize_state_restarted(user%fi, user%rho, user%walls, &
            options%istep, options%kwrite, user)
       istep = options%istep
    endif

    ! start lbm
    if(user%info%id.eq.0) print*,'calling lbm'

    call LBMRun(user, istep, options%ntimes*options%npasses, options%kwrite) ! for the moment, this is crap
    call LBMDestroy(user, ierr)
    call PetscFinalize(ierr)
    stop
  end program main
  !----------------------------------------------------------

