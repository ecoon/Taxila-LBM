!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        main.F90
!!!     version:
!!!     created:         08 December 2010
!!!       on:            11:48:19 MST
!!!     last modified:   14 December 2010
!!!       at:            12:16:41 MST
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

    ! --- setup environment
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call constants_initialize(s)
    user => LBMCreate(PETSC_COMM_WORLD)
    call OptionsSetSizes(user%options, s)

    ! --- read parameters
    infile = 'input_data'
    call OptionsReadFile(user%options, infile)
    call constants_set_from_options(user%options)

    ! --- initialize memory
    call LBMSetSizes(user, user%options%NX, user%options%NY, user%options%NZ, s, b)
    user%bc%flags(:) = user%options%bc_flags(:)

    if (user%info%id.eq.0) call OptionsPrint(user%options)
    ! --- initialize state
    ! walls
    if(user%info%id.eq.0) print*,'initialization of walls'
    call LBMInitializeWalls(user, initialize_walls)

    ! bcs
    call BCSetValues(user%bc, user%info, initialize_bcs)

    ! fi/state
    if (user%options%new_simulation) then
       call LBMInitializeState(user, initialize_state)
       istep=1
    else
       call initialize_state_restarted(user%fi, user%rho, user%walls, &
            user%options%istep, user%options%kwrite, user)
       istep = user%options%istep
    endif

    ! start lbm
    if(user%info%id.eq.0) print*,'calling lbm'

    call LBMRun(user, istep, user%options%ntimes*user%options%npasses, user%options%kwrite) ! for the moment, this is crap
    call LBMDestroy(user, ierr)
    call PetscFinalize(ierr)
    stop
  end program main
  !----------------------------------------------------------

