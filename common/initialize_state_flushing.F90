!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            18:21:06 MST
!!!     last modified:   05 May 2011
!!!       at:            10:22:10 MDT
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
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%s, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: rho
    PetscScalar,dimension(dist%s, 1:dist%info%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: u
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: walls

    ! local variables
    PetscErrorCode ierr
    PetscBool flag
    logical,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: bound
    PetscScalar,dimension(dist%s):: rho1, rho2         ! left and right fluid densities?
    PetscInt nmax
    PetscBool flushx, flushy, flushz
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: nowalls

    PetscInt i,j,k,m ! local values
    PetscBool help

    ! input data
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    rho1 = 0.d0
    if (help) call PetscPrintf(options%comm, "-rho_invading=<0,0>: phase density of the invading fluid", ierr)
    nmax = dist%s
    call PetscOptionsGetRealArray(options%my_prefix, '-rho_invading', rho1, nmax, flag, ierr)

    rho2 = 0.d0
    if (help) call PetscPrintf(options%comm, "-rho_defending=<0,0>: phase density of the defending fluid", ierr)
    nmax = dist%s
    call PetscOptionsGetRealArray(options%my_prefix, '-rho_defending', rho2, nmax, flag, ierr)
    
    flushz = .TRUE.
    flushy = .FALSE.
    flushx = .FALSE.
    if (help) call PetscPrintf(options%comm, "-flush_direction_{xyz}: set the direction"// &
         " of flushing", ierr)
    call PetscOptionsGetBool(options%my_prefix, '-flush_direction_x', flushx, flag, ierr)
    call PetscOptionsGetBool(options%my_prefix, '-flush_direction_y', flushy, flag, ierr)
    call PetscOptionsGetBool(options%my_prefix, '-flush_direction_z', flushz, flag, ierr)
    
    ! initialize state
    fi=0.0
    u=0.0

    ! flushing experiement 
    if (flushx) then
       bound=.false.
       do i=dist%info%xs,dist%info%xe
          do j=dist%info%ys,dist%info%ye
             do k=dist%info%zs,dist%info%ze
                if(i.le.10) bound(i,j,k)=.true.
             enddo
          enddo
       enddo
    else if (flushy) then
       bound=.false.
       do i=dist%info%xs,dist%info%xe
          do j=dist%info%ys,dist%info%ye
             do k=dist%info%zs,dist%info%ze
                if(j.le.10) bound(i,j,k)=.true.
             enddo
          enddo
       enddo
    else if (flushz) then
       bound=.false.
       do i=dist%info%xs,dist%info%xe
          do j=dist%info%ys,dist%info%ye
             do k=dist%info%zs,dist%info%ze
                if(k.le.10) bound(i,j,k)=.true.
             enddo
          enddo
       enddo
    end if

    ! set density
    where(bound)
       rho(1,:,:,:)=rho1(1)
       rho(2,:,:,:)=rho1(2)
    else where
       rho(1,:,:,:)=rho2(1)
       rho(2,:,:,:)=rho2(2)
    end where
    
    ! set state at equilibrium       
    nowalls = 0.d0
    do m=1,dist%s
       call DiscretizationEquilf(dist%disc, rho(m,:,:,:), u(m,:,:,:,:), &
            nowalls, fi(m,:,:,:,:), phases(m)%relax, dist)    
    end do
    return
  end subroutine initialize_state
