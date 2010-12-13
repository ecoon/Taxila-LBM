!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        main.F90
!!!     version:
!!!     created:         08 December 2010
!!!       on:            11:48:19 MST
!!!     last modified:   09 December 2010
!!!       at:            15:16:54 MST
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

    integer,parameter:: s=2
    integer,parameter:: b=18

    external initialize_bcs
    external initialize_state
    external initialize_walls
    type(lbm_type),pointer:: user
    type(options_type),pointer:: options

    

    ! --- setup environment
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call constants_initialize(s)
    options => OptionsCreate(s)

    ! --- read parameters
    infile = 'input_data'
    call OptionsReadFile(options, infile)
    call constants_set_from_options(options)

    ! --- initialize memory
    user => LBMCreate(PETSC_COMM_WORLD)
    call LBMSetSizes(user, options%NX, options%NY, options%NZ, s, b)
    user%bc%flags(:) = options%bc_flags(:)

    if (user%info%id.eq.0) call OptionsPrint(options)
    ! --- initialize state
    ! walls
    if(user%info%id.eq.0) print*,'initialization of walls'
    call LBMInitializeWalls(user, initialize_walls)

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

