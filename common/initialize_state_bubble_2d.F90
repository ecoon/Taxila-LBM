!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state.F90
!!!     version:
!!!     created:         14 January 2011
!!!       on:            17:27:04 MST
!!!     last modified:   25 April 2011
!!!       at:            15:35:59 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"
  
  subroutine initialize_state(fi, rho, u, walls, dist, phases, options)
    use petsc
    use LBM_Distribution_Function_type_module
    use LBM_Phase_module
    use LBM_Options_module
    use LBM_Discretization_module
    implicit none

    ! input variables
    type(distribution_type) dist
    type(phase_type) phases(dist%s)
    type(options_type) options

    PetscScalar,dimension(dist%s,0:dist%b, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%s, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: rho
    PetscScalar,dimension(dist%s,1:dist%info%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls

    ! local variables
    PetscErrorCode ierr
    PetscBool flag
    logical,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: bound
    PetscScalar,dimension(dist%s):: rho1, rho2   ! initial region densities
    PetscInt nmax
    character(len=30) bubblestring

    PetscInt i,j,m ! local values
    PetscInt lx, rx, ly, ry

    ! input data
    PetscBool help

    ! input data
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    rho1 = 0.d0
    if (help) call PetscPrintf(options%comm, "-rho_inner=<0,0>: ???", ierr)
    nmax = dist%s
    call PetscOptionsGetRealArray(options%my_prefix, '-rho_inner', rho1, nmax, flag, ierr)

    rho2 = 0.d0
    if (help) call PetscPrintf(options%comm, "-rho_outer=<0,0>: ???", ierr)
    nmax = dist%s
    call PetscOptionsGetRealArray(options%my_prefix, '-rho_outer', rho2, nmax, flag, ierr)

    ! initialize state
    fi=0.0
    u=0.0

    !---duct with nonwetting bubble in middle ----------
    !--- rho2: nonwetting
    !--- rho1: wetting

    ! for odd length side of square
    lx = (dist%info%NX+1)/2-26
    rx = (dist%info%NX+1)/2+26
    ly = (dist%info%NY+1)/2-26
    ry = (dist%info%NY+1)/2+26

    write(bubblestring,'(a,I4,I4,I4,I4)') 'bubble format:', ly, ry
    call PetscPrintf(options%comm, bubblestring, ierr)

    ! set bound to indicate fluid
    bound=.false.
    do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe    
        if((i.ge.lx.and.i.le.rx).and.(j.ge.ly.and.j.le.ry)) then
          bound(i,j)=.true.
        endif
      enddo
    enddo

    ! set density
    where(bound)
       rho(1,:,:)=rho1(1)
       rho(2,:,:)=rho1(2)
    else where
       rho(1,:,:)=rho2(1)
       rho(2,:,:)=rho2(2)
    end where

    where(walls.eq.1)
       rho(1,:,:) = 0
       rho(2,:,:) = 0
    end where

    ! set state at equilibrium       
    do m=1,dist%s
      call DiscretizationEquilf(dist%disc, rho(m,:,:), u(m,:,:,:), &
           walls, fi(m,:,:,:), phases(m)%relax, dist)    
    end do
    return
  end subroutine initialize_state
