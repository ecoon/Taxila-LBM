!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            18:21:06 MST
!!!     last modified:   21 October 2011
!!!       at:            18:16:35 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

! initializes the state to a constant gradient in rho in the
! x-direction, given by the left and right rho values given as BCs

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
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, &
         dist%info%rgzs:dist%info%rgze):: nowalls

    PetscScalar :: xp3_ave_p, xm3_ave_p, t
    PetscBool flag
    PetscInt i,one
    PetscInt direction

    ! initialize state
    fi = 0.
    u = 0.

    xp3_ave_p = 0.
    xm3_ave_p = 0.

    direction = 1
    call PetscOptionsGetInt(options%my_prefix, '-gradient_direction', direction, &
         flag, ierr)

    select case(direction)
    case(1)
      call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_xp_component1', xp3_ave_p, &
           flag, ierr)
      if (.not.flag) SETERRQ(1, 1, 'invalid boundary pressure for xp_component1', ierr)
      call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_xm_component1', xm3_ave_p, &
           flag, ierr)
      if (.not.flag) SETERRQ(1, 1, 'invalid boundary value for xm_component1', ierr)
      
      do i=dist%info%xs,dist%info%xe
        t = dble(i-1)/dble(dist%info%NX-1)
        rho(1,i,:,:) = (1-t)*xm3_ave_p + t*xp3_ave_p
      end do
    case(2)
      call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_yp_component1', xp3_ave_p, &
           flag, ierr)
      if (.not.flag) SETERRQ(1, 1, 'invalid boundary pressure for yp_component1', ierr)
      call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_ym_component1', xm3_ave_p, &
           flag, ierr)
      if (.not.flag) SETERRQ(1, 1, 'invalid boundary value for ym_component1', ierr)
      
      do i=dist%info%ys,dist%info%ye
        t = dble(i-1)/dble(dist%info%NY-1)
        rho(1,:,i,:) = (1-t)*xm3_ave_p + t*xp3_ave_p
      end do
    case(3)
      call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_zp_component1', xp3_ave_p, &
           flag, ierr)
      if (.not.flag) SETERRQ(1, 1, 'invalid boundary pressure for zp_component1', ierr)
      call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_zm_component1', xm3_ave_p, &
           flag, ierr)
      if (.not.flag) SETERRQ(1, 1, 'invalid boundary value for zm_component1', ierr)
      
      do i=dist%info%zs,dist%info%ze
        t = dble(i-1)/dble(dist%info%NZ-1)
        rho(1,:,:,i) = (1-t)*xm3_ave_p + t*xp3_ave_p
      end do
    end select
    
    ! set state at equilibrium       
    nowalls = 0.
    one = 1
    call DiscretizationEquilf(dist%disc, rho, u, nowalls, fi, one, components(1)%relax, dist)    
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
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: nowalls

    PetscScalar :: xp3_ave_p, xm3_ave_p, t
    PetscBool flag
    PetscInt i,one
    PetscInt direction

    ! initialize state
    fi = 0.
    u = 0.

    xp3_ave_p = 0.
    xm3_ave_p = 0.

    direction = 1
    call PetscOptionsGetInt(options%my_prefix, '-gradient_direction', direction, &
         flag, ierr)

    select case(direction)
    case(1)
      call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_xp_component1', xp3_ave_p, &
           flag, ierr)
      if (.not.flag) SETERRQ(1, 1, 'invalid boundary pressure for xp_component1', ierr)
      call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_xm_component1', xm3_ave_p, &
           flag, ierr)
      if (.not.flag) SETERRQ(1, 1, 'invalid boundary value for xm_component1', ierr)
      
      do i=dist%info%xs,dist%info%xe
        t = dble(i-1)/dble(dist%info%NX-1)
        rho(1,i,:) = (1-t)*xm3_ave_p + t*xp3_ave_p
      end do
    case(2)
      call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_yp_component1', xp3_ave_p, &
           flag, ierr)
      if (.not.flag) SETERRQ(1, 1, 'invalid boundary pressure for yp_component1', ierr)
      call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_ym_component1', xm3_ave_p, &
           flag, ierr)
      if (.not.flag) SETERRQ(1, 1, 'invalid boundary value for ym_component1', ierr)
      
      do i=dist%info%ys,dist%info%ye
        t = dble(i-1)/dble(dist%info%NY-1)
        rho(1,:,i) = (1-t)*xm3_ave_p + t*xp3_ave_p
      end do
    end select

    ! set state at equilibrium       
    nowalls = 0.
    one = 1
    call DiscretizationEquilf(dist%disc, rho, u, nowalls, fi, one, components(1)%relax, dist)    
    return
  end subroutine initialize_state_d2
