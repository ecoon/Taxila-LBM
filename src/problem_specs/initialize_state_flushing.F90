!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            18:21:06 MST
!!!     last modified:   14 September 2011
!!!       at:            12:39:30 PDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

  subroutine initialize_state(fi, rho, u, walls, dist, components, options)
    use petsc
    use LBM_Distribution_Function_type_module
    use LBM_Component_module
    use LBM_Options_module
    use LBM_Discretization_module
    implicit none

    ! input variables
    type(distribution_type) dist
    type(component_type) components(dist%s)
    type(options_type) options
    PetscScalar,dimension(dist%s,0:dist%b,dist%info%gxyzl) :: fi
    PetscScalar,dimension(dist%s,dist%info%rgxyzl) :: rho
    PetscScalar,dimension(dist%s, 1:dist%info%ndims, dist%info%gxyzl):: u
    PetscScalar,dimension(dist%info%rgxyzl):: walls

    select case(dist%info%ndims)
    case (2) 
      call initialize_state_d2(fi, rho, u, walls, dist, components, options)
    case (3) 
      call initialize_state_d3(fi, rho, u, walls, dist, components, options)
    end select
  end subroutine initialize_state

  subroutine initialize_state_d3(fi, rho, u, walls, dist, components, options)
    use petsc
    use LBM_Distribution_Function_type_module
    use LBM_Component_module
    use LBM_Options_module
    use LBM_Discretization_module
    implicit none

    ! input variables
    type(distribution_type) dist
    type(component_type) components(dist%s)
    type(options_type) options
    PetscScalar,dimension(dist%s,0:dist%b, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%s, &
         dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, &
         dist%info%rgzs:dist%info%rgze):: rho
    PetscScalar,dimension(dist%s, 1:dist%info%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: u
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, &
         dist%info%rgzs:dist%info%rgze):: walls

    ! local variables
    PetscErrorCode ierr
    PetscBool flag
    PetscScalar,dimension(dist%s):: rho1, rho2         ! left and right fluid densities?
    PetscInt nmax
    PetscBool flushx, flushy, flushz
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, &
         dist%info%rgzs:dist%info%rgze):: nowalls

    PetscInt i,j,k,m ! local values
    PetscBool help

    ! input data
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    rho1 = 0.d0
    if (help) call PetscPrintf(options%comm, &
         "-rho_invading=<0,0>: component density of the invading fluid", ierr)
    nmax = dist%s
    call PetscOptionsGetRealArray(options%my_prefix, '-rho_invading', &
         rho1, nmax, flag, ierr)

    rho2 = 0.d0
    if (help) call PetscPrintf(options%comm, &
         "-rho_defending=<0,0>: component density of the defending fluid", ierr)
    nmax = dist%s
    call PetscOptionsGetRealArray(options%my_prefix, '-rho_defending', &
         rho2, nmax, flag, ierr)
    
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
      do k=dist%info%zs,dist%info%ze
      do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j,k).eq.0) then
          if (i.le.10) then
            rho(:,i,j,k)=rho1(:)
          else 
            rho(:,i,j,k)=rho2(:)
          end if
        end if
      enddo
      enddo
      enddo
    else if (flushy) then
      do k=dist%info%zs,dist%info%ze
      do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j,k).eq.0) then
          if (j.le.10) then
            rho(:,i,j,k)=rho1(:)
          else 
            rho(:,i,j,k)=rho2(:)
          end if
        end if
      enddo
      enddo
      enddo
    else if (flushz) then
      do k=dist%info%zs,dist%info%ze
      do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j,k).eq.0) then
          if (k.le.10) then
            rho(:,i,j,k)=rho1(:)
          else 
            rho(:,i,j,k)=rho2(:)
          end if
        end if
      enddo
      enddo
      enddo
    end if

    ! set state at equilibrium       
    nowalls = 0.d0
    do m=1,dist%s
       call DiscretizationEquilf(dist%disc, rho, u, &
            nowalls, fi, m, components(m)%relax, dist)    
    end do
    return
  end subroutine initialize_state_d3


  subroutine initialize_state_d2(fi, rho, u, walls, dist, components, options)
    use petsc
    use LBM_Distribution_Function_type_module
    use LBM_Component_module
    use LBM_Options_module
    use LBM_Discretization_module
    implicit none

    ! input variables
    type(distribution_type) dist
    type(component_type) components(dist%s)
    type(options_type) options
    PetscScalar,dimension(dist%s,0:dist%b, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%s, &
         dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(dist%s, 1:dist%info%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls

    ! local variables
    PetscErrorCode ierr
    PetscBool flag
    logical,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: bound
    PetscScalar,dimension(dist%s):: rho1, rho2         ! left and right fluid densities?
    PetscInt nmax
    PetscBool flushx, flushy
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: nowalls

    PetscInt i,j,m ! local values
    PetscBool help

    ! input data
    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    rho1 = 0.d0
    if (help) call PetscPrintf(options%comm, "-rho_invading=<0,0>: ???", ierr)
    nmax = dist%s
    call PetscOptionsGetRealArray(options%my_prefix, '-rho_invading', rho1, nmax, flag, ierr)

    rho2 = 0.d0
    if (help) call PetscPrintf(options%comm, "-rho_defending=<0,0>: ???", ierr)
    nmax = dist%s
    call PetscOptionsGetRealArray(options%my_prefix, '-rho_defending', rho2, nmax, flag, ierr)
    
    flushy = .TRUE.
    flushx = .FALSE.
    if (help) call PetscPrintf(options%comm, "-flush_direction_{xy}: set the direction"// &
         " of flushing", ierr)
    call PetscOptionsGetBool(options%my_prefix, '-flush_direction_x', flushx, flag, ierr)
    call PetscOptionsGetBool(options%my_prefix, '-flush_direction_y', flushy, flag, ierr)
    
    ! initialize state
    fi=0.d0
    u=0.d0

    ! flushing experiement 
    if (flushx) then
      do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j).eq.0) then
          if (i.le.10) then
            rho(:,i,j)=rho1(:)
          else 
            rho(:,i,j)=rho2(:)
          end if
        end if
      enddo
      enddo
    else if (flushy) then
      do j=dist%info%ys,dist%info%ye
      do i=dist%info%xs,dist%info%xe
        if (walls(i,j).eq.0) then
          if (j.le.10) then
            rho(:,i,j)=rho1(:)
          else 
            rho(:,i,j)=rho2(:)
          end if
        end if
      enddo
      enddo
    end if
    
    ! set state at equilibrium       
    nowalls = 0.d0
    do m=1,dist%s
       call DiscretizationEquilf(dist%disc, rho, u, &
            nowalls, fi, m, components(m)%relax, dist)    
    end do
    return
  end subroutine initialize_state_d2
