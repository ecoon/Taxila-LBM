!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            18:21:06 MST
!!!     last modified:   31 October 2011
!!!       at:            15:53:22 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petscdmdef.h"

  subroutine initialize_state(rho, u, walls, dist, components, options)
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
    PetscScalar,dimension(dist%s,dist%info%rgxyzl) :: rho
    PetscScalar,dimension(dist%s, 1:dist%info%ndims, dist%info%gxyzl):: u
    PetscScalar,dimension(dist%info%rgxyzl):: walls

    select case(dist%info%ndims)
    case (2) 
      call initialize_state_d2(rho, u, walls, dist, components, options)
    case (3) 
      call initialize_state_d3(rho, u, walls, dist, components, options)
    end select
  end subroutine initialize_state

  subroutine initialize_state_d3(rho, u, walls, dist, components, options)
    use petsc
    use LBM_Distribution_Function_type_module
    use LBM_Component_module
    use LBM_Options_module
    use LBM_Discretization_module
    implicit none

#include "lbm_definitions.h"

    ! input variables
    type(distribution_type) dist
    type(component_type) components(dist%s)
    type(options_type) options
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
    PetscReal val_left, val_right
    PetscReal,dimension(3) :: vels
    PetscBool flag, bcpresent
    PetscInt flow_direction, count
    PetscInt bc_left, bc_right
    PetscInt i,j,k
    PetscReal t

    ! assumes bc is one of pressure, flux, or velocity
    flow_direction = X_DIRECTION
    call OptionsGetInt(options, "-flow_direction", &
         "initialization for a single phase flow problem", &
         flow_direction, flag, ierr)

    val_left = 0.d0
    val_right = 0.d0
    bc_left = BC_NULL
    bc_right = BC_NULL

    select case(flow_direction)
    case(X_DIRECTION)
      flag = PETSC_FALSE
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_pressure_xm_value", &
             val_left, flag, ierr)
        if (flag) bc_left = BC_PRESSURE
      endif
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_flux_xm_value", &
             val_left, flag, ierr)
        if (flag) bc_left = BC_FLUX
      endif
      if (.not.flag) then
        count = dist%info%ndims
        vels = 0.d0
        call PetscOptionsGetRealArray(options%my_prefix, "-bc_velocity_xm_values", &
             vels, count, flag, ierr)
        if (flag) then
          bc_left = BC_VELOCITY
          val_left = vels(X_DIRECTION)
        endif
      endif

      flag = PETSC_FALSE
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_pressure_xp_value", &
             val_right, flag, ierr)
        if (flag) bc_right = BC_PRESSURE
      endif
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_flux_xp_value", &
             val_right, flag, ierr)
        if (flag) bc_right = BC_FLUX
      endif
      if (.not.flag) then
        count = dist%info%ndims
        vels = 0.d0
        call PetscOptionsGetRealArray(options%my_prefix, "-bc_velocity_xp_values", &
             vels, count, flag, ierr)
        if (flag) then
          bc_right = BC_VELOCITY
          val_right = vels(X_DIRECTION)
        endif
      endif

    case(Y_DIRECTION)
      flag = PETSC_FALSE
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_pressure_ym_value", &
             val_left, flag, ierr)
        if (flag) bc_left = BC_PRESSURE
      endif
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_flux_ym_value", &
             val_left, flag, ierr)
        if (flag) bc_left = BC_FLUX
      endif
      if (.not.flag) then
        count = dist%info%ndims
        vels = 0.d0
        call PetscOptionsGetRealArray(options%my_prefix, "-bc_velocity_ym_values", &
             vels, count, flag, ierr)
        if (flag) then
          bc_left = BC_VELOCITY
          val_left = vels(Y_DIRECTION)
        endif
      endif

      flag = PETSC_FALSE
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_pressure_yp_value", &
             val_right, flag, ierr)
        if (flag) bc_right = BC_PRESSURE
      endif
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_flux_yp_value", &
             val_right, flag, ierr)
        if (flag) bc_right = BC_FLUX
      endif
      if (.not.flag) then
        count = dist%info%ndims
        vels = 0.d0
        call PetscOptionsGetRealArray(options%my_prefix, "-bc_velocity_yp_values", &
             vels, count, flag, ierr)
        if (flag) then
          bc_right = BC_VELOCITY
          val_right = vels(Y_DIRECTION)
        endif
      endif

    case(Z_DIRECTION)
      flag = PETSC_FALSE
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_pressure_zm_value", &
             val_left, flag, ierr)
        if (flag) bc_left = BC_PRESSURE
      endif
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_flux_zm_value", &
             val_left, flag, ierr)
        if (flag) bc_left = BC_FLUX
      endif
      if (.not.flag) then
        count = dist%info%ndims
        vels = 0.d0
        call PetscOptionsGetRealArray(options%my_prefix, "-bc_velocity_zm_values", &
             vels, count, flag, ierr)
        if (flag) then
          bc_left = BC_VELOCITY
          val_left = vels(Z_DIRECTION)
        endif
      endif

      flag = PETSC_FALSE
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_pressure_zp_value", &
             val_right, flag, ierr)
        if (flag) bc_right = BC_PRESSURE
      endif
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_flux_zp_value", &
             val_right, flag, ierr)
        if (flag) bc_right = BC_FLUX
      endif
      if (.not.flag) then
        count = dist%info%ndims
        vels = 0.d0
        call PetscOptionsGetRealArray(options%my_prefix, "-bc_velocity_zp_values", &
             vels, count, flag, ierr)
        if (flag) then
          bc_right = BC_VELOCITY
          val_right = vels(Z_DIRECTION)
        endif
      endif
    end select

    u = 0.
    rho = 1.

    ! initialize state
    if (bc_left.eq.BC_PRESSURE) then
      if (bc_right.eq.BC_PRESSURE) then
        ! u = 0.d0  ! (already set)
        select case(flow_direction)
        case(X_DIRECTION)
          do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
            t = dble(i-1)/dble(dist%info%NX-1)
            rho(1,i,j,k) = (val_left*(1-t) + val_right*t)*3.d0
          end do
          end do
          end do
        case(Y_DIRECTION)
          do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
            t = dble(j-1)/dble(dist%info%NY-1)
            rho(1,i,j,k) = (val_left*(1-t) + val_right*t)*3.d0
          end do
          end do
          end do
        case(Z_DIRECTION)
          do k=dist%info%zs,dist%info%ze
          do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
            t = dble(k-1)/dble(dist%info%NZ-1)
            rho(1,i,j,k) = (val_left*(1-t) + val_right*t)*3.d0
          end do
          end do
          end do
        end select
      else if (bc_right.eq.BC_FLUX) then
        rho(1,:,:,:) = val_left*3.d0
        u(1,flow_direction,:,:,:) = -val_right/(val_left*3.d0)
      else if (bc_right.eq.BC_VELOCITY) then
        rho(1,:,:,:) = val_left*3.d0
        u(1,flow_direction,:,:,:) = val_right
      endif

    else if (bc_left.eq.BC_FLUX) then
      if (bc_right.eq.BC_PRESSURE) then
        rho(1,:,:,:) = val_right*3.d0
        u(1,flow_direction,:,:,:) = val_left/(val_right*3.d0)
      else if (bc_right.eq.BC_FLUX) then
        rho(1,:,:,:) = 1.d0
        u(1,flow_direction,:,:,:) = (val_left - val_right)/2.d0
      else if (bc_right.eq.BC_VELOCITY) then
        rho(1,:,:,:) = val_left/val_right
        u(1,flow_direction,:,:,:) = val_right
      endif

    else if (bc_left.eq.BC_VELOCITY) then
      if (bc_right.eq.BC_PRESSURE) then
        rho(1,:,:,:) = val_right*3.d0
        u(1,flow_direction,:,:,:) = val_left
      else if (bc_right.eq.BC_FLUX) then
        rho(1,:,:,:) = -val_right/val_left
        u(1,flow_direction,:,:,:) = val_left
      else if (bc_right.eq.BC_VELOCITY) then
        rho(1,:,:,:) = 1.d0
        u(1,flow_direction,:,:,:) = (val_left+val_right)/2.d0
      endif
    end if
    return
  end subroutine initialize_state_d3

  subroutine initialize_state_d2(rho, u, walls, dist, components, options)
    use petsc
    use LBM_Distribution_Function_type_module
    use LBM_Component_module
    use LBM_Options_module
    use LBM_Discretization_module
    implicit none

#include "lbm_definitions.h"

    ! input variables
    type(distribution_type) dist
    type(component_type) components(dist%s)
    type(options_type) options
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
    PetscReal val_left, val_right
    PetscReal,dimension(3) :: vels
    PetscBool flag, bcpresent
    PetscInt flow_direction, count
    PetscInt bc_left, bc_right
    PetscInt i,j,k
    PetscReal t

    ! assumes bc is one of pressure, flux, or velocity
    flow_direction = X_DIRECTION
    call OptionsGetInt(options, "-flow_direction", &
         "initialization for a single phase flow problem", &
         flow_direction, flag, ierr)

    val_left = 0.d0
    val_right = 0.d0
    bc_left = BC_NULL
    bc_right = BC_NULL

    select case(flow_direction)
    case(X_DIRECTION)
      flag = PETSC_FALSE
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_pressure_xm_value", &
             val_left, flag, ierr)
        if (flag) bc_left = BC_PRESSURE
      endif
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_flux_xm_value", &
             val_left, flag, ierr)
        if (flag) bc_left = BC_FLUX
      endif
      if (.not.flag) then
        count = dist%info%ndims
        vels = 0.d0
        call PetscOptionsGetRealArray(options%my_prefix, "-bc_velocity_xm_values", &
             vels, count, flag, ierr)
        if (flag) then
          bc_left = BC_VELOCITY
          val_left = vels(X_DIRECTION)
        endif
      endif

      flag = PETSC_FALSE
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_pressure_xp_value", &
             val_right, flag, ierr)
        if (flag) bc_right = BC_PRESSURE
      endif
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_flux_xp_value", &
             val_right, flag, ierr)
        if (flag) bc_right = BC_FLUX
      endif
      if (.not.flag) then
        count = dist%info%ndims
        vels = 0.d0
        call PetscOptionsGetRealArray(options%my_prefix, "-bc_velocity_xp_values", &
             vels, count, flag, ierr)
        if (flag) then
          bc_right = BC_VELOCITY
          val_right = vels(X_DIRECTION)
        endif
      endif

    case(Y_DIRECTION)
      flag = PETSC_FALSE
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_pressure_ym_value", &
             val_left, flag, ierr)
        if (flag) bc_left = BC_PRESSURE
      endif
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_flux_ym_value", &
             val_left, flag, ierr)
        if (flag) bc_left = BC_FLUX
      endif
      if (.not.flag) then
        count = dist%info%ndims
        vels = 0.d0
        call PetscOptionsGetRealArray(options%my_prefix, "-bc_velocity_ym_values", &
             vels, count, flag, ierr)
        if (flag) then
          bc_left = BC_VELOCITY
          val_left = vels(Y_DIRECTION)
        endif
      endif

      flag = PETSC_FALSE
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_pressure_yp_value", &
             val_right, flag, ierr)
        if (flag) bc_right = BC_PRESSURE
      endif
      if (.not.flag) then
        call PetscOptionsGetReal(options%my_prefix, "-bc_flux_yp_value", &
             val_right, flag, ierr)
        if (flag) bc_right = BC_FLUX
      endif
      if (.not.flag) then
        count = dist%info%ndims
        vels = 0.d0
        call PetscOptionsGetRealArray(options%my_prefix, "-bc_velocity_yp_values", &
             vels, count, flag, ierr)
        if (flag) then
          bc_right = BC_VELOCITY
          val_right = vels(Y_DIRECTION)
        endif
      endif
    end select

    u = 0.
    rho = 1.

    ! initialize state
    if (bc_left.eq.BC_PRESSURE) then
      if (bc_right.eq.BC_PRESSURE) then
        ! u = 0.d0  ! (already set)
        select case(flow_direction)
        case(X_DIRECTION)
          do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
            t = dble(i-1)/dble(dist%info%NX-1)
            rho(1,i,j) = (val_left*(1-t) + val_right*t)*3.d0
          end do
          end do
        case(Y_DIRECTION)
          do j=dist%info%ys,dist%info%ye
          do i=dist%info%xs,dist%info%xe
            t = dble(j-1)/dble(dist%info%NY-1)
            rho(1,i,j) = (val_left*(1-t) + val_right*t)*3.d0
          end do
          end do
        end select
      else if (bc_right.eq.BC_FLUX) then
        rho(1,:,:) = val_left*3.d0
        u(1,flow_direction,:,:) = -val_right/(val_left*3.d0)
      else if (bc_right.eq.BC_VELOCITY) then
        rho(1,:,:) = val_left*3.d0
        u(1,flow_direction,:,:) = val_right
      endif

    else if (bc_left.eq.BC_FLUX) then
      if (bc_right.eq.BC_PRESSURE) then
        rho(1,:,:) = val_right*3.d0
        u(1,flow_direction,:,:) = val_left/(val_right*3.d0)
      else if (bc_right.eq.BC_FLUX) then
        rho(1,:,:) = 1.d0
        u(1,flow_direction,:,:) = (val_left - val_right)/2.d0
      else if (bc_right.eq.BC_VELOCITY) then
        rho(1,:,:) = val_left/val_right
        u(1,flow_direction,:,:) = val_right
      endif

    else if (bc_left.eq.BC_VELOCITY) then
      if (bc_right.eq.BC_PRESSURE) then
        rho(1,:,:) = val_right*3.d0
        u(1,flow_direction,:,:) = val_left
      else if (bc_right.eq.BC_FLUX) then
        rho(1,:,:) = -val_right/val_left
        u(1,flow_direction,:,:) = val_left
      else if (bc_right.eq.BC_VELOCITY) then
        rho(1,:,:) = 1.d0
        u(1,flow_direction,:,:) = (val_left+val_right)/2.d0
      endif
    end if
    return
  end subroutine initialize_state_d2
