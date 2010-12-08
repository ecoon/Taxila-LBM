!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        main.F90
!!!     version:
!!!     created:         08 December 2010
!!!       on:            11:48:19 MST
!!!     last modified:   08 December 2010
!!!       at:            12:37:59 MST
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
    use constants
    use petsc
    implicit none

    logical new
    integer istep             ! restart disabled at the moment
    PetscErrorCode ierr

    integer NX, NY, NZ             ! global domain size

    integer,parameter:: s=2
    integer,parameter:: b=18

    integer,dimension(1:6):: bc_flags
    external initialize_bcs
    external initialize_state
    external initialize_walls
    type(lbm_type),pointer:: user


    ! --- setup environment
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call constants_initialize(s)

    ! --- read parameters from node id=0

    open(90,file='input_data',status='unknown')
    read(90,*) NX,NY,NZ
    read(90,*) g,g11,g22,gw(1),gw(2)
    read(90,*) tau(1),tau(2)
    read(90,*) gvt(1),gvt(2)
    read(90,*) mm(1),mm(2)
    read(90,*) rhol,rhor
    read(90,*) ntimes
    read(90,*) npasses
    read(90,*) new
    read(90,*) istep          ! istep=0 if new=.true.
    read(90,*) kprint,kwrite
    read(90,*) bc_flags(1),bc_flags(2),bc_flags(3),bc_flags(4),bc_flags(5),bc_flags(6)

    close(90)

    ! --- initialize memory
    user => LBMCreate(PETSC_COMM_WORLD)
    call LBMSetSizes(user, NX, NY, NZ, s, b)

    ! --- output logging info
    if(user%info%id.eq.0) then
       write(*,*) 'ntimes = ',ntimes,' npasses= ',npasses
       write(*,*) 'NX,NY,NZ'
       write(*,*)  NX,NY,NZ
       write(*,*) 'g,g11,g22,g1w,g2w'
       write(*,*)  g,g11,g22,gw(1),gw(2)
       write(*,*) 'tau1,tau2'
       write(*,*)  tau(1),tau(2)
       write(*,*) 'gvt1,gvt2'
       write(*,*)  gvt(1),gvt(2)
       write(*,*) 'mm(1),mm(2)'
       write(*,*)  mm(1),mm(2)
       write(*,*) 'rhol,rhor'
       write(*,*)  rhol,rhor
       write(*,*) 'new,kprint,kwrite'
       write(*,*)  new,kprint,kwrite
       write(*,*) 'bc flags: 0=periodic, 1=pseudo-periodic, 2=flux, 3=pressure'
       write(*,*) bc_flags
    endif
    user%bc%flags(:) = bc_flags(:)

    ! --- initialize state
    ! walls
    if(user%info%id.eq.0) print*,'initialization of walls'
    call LBMInitializeWalls(user, initialize_walls)

    ! bcs
    call BCSetValues(user%bc, user%info, initialize_bcs)

    ! fi/state
    if (new) then
       call LBMInitializeState(user, initialize_state)
       istep=1
    else
       call initialize_state_restarted(user%fi, user%rho, user%walls, &
            istep, kwrite, user)
    endif

    ! start lbm
    if(user%info%id.eq.0) print*,'calling lbm'

    call LBMRun(user, istep, ntimes*npasses, kwrite) ! for the moment, this is crap
    call LBMDestroy(user, ierr)
    call constants_finalize()
    call PetscFinalize(ierr)
    stop
  end program main
  !----------------------------------------------------------

