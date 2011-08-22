!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_bc_zfluxflux.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            17:30:22 MST
!!!     last modified:   05 May 2011
!!!       at:            10:14:37 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

! initializes the BCs from options to be a general constant flux/pressure bc
  subroutine initialize_bcs(bc_flags, xm_bcvals, xp_bcvals, ym_bcvals, & 
       yp_bcvals, zm_bcvals, zp_bcvals, bc_dim, dist, options)
    use petsc
    use LBM_Distribution_Function_type_module
    use LBM_Options_module
    implicit none

#include "lbm_definitions.h"    
    PetscInt, dimension(6):: bc_flags ! enum for boundary conditions

    type(distribution_type) dist
    type(options_type) options
    PetscScalar,dimension(dist%s, dist%info%ndims, dist%info%ys:dist%info%ye):: xm_bcvals, xp_bcvals
    PetscScalar,dimension(dist%s, dist%info%ndims, dist%info%xs:dist%info%xe):: ym_bcvals, yp_bcvals
    PetscScalar,pointer:: zm_bcvals, zp_bcvals
    PetscInt bc_dim

    ! local
    PetscScalar xp3_ave, xm3_ave, xp3_max, xm3_max
    PetscScalar,dimension(dist%s):: xp3_ave_p, xm3_ave_p
    PetscScalar,dimension(dist%s,dist%info%ndims):: xp3_ave_f, xm3_ave_f
    PetscBool pressure_xm, pressure_xp, flux_xm, flux_xp
    PetscBool velocity_xm, velocity_xp, poise_xm, poise_xp

    PetscScalar yp3_ave, ym3_ave, yp3_max, ym3_max
    PetscScalar,dimension(dist%s):: yp3_ave_p, ym3_ave_p
    PetscScalar,dimension(dist%s,dist%info%ndims):: yp3_ave_f, ym3_ave_f
    PetscBool pressure_ym, pressure_yp, flux_ym, flux_yp
    PetscBool velocity_ym, velocity_yp, poise_ym, poise_yp

    PetscInt i, j, m, d
    PetscBool flag
    PetscErrorCode ierr
    PetscBool help

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    ! get options
    pressure_xm = .FALSE.
    pressure_xp = .FALSE.
    pressure_ym = .FALSE.
    pressure_yp = .FALSE.

    velocity_xm = .FALSE.
    velocity_xp = .FALSE.
    velocity_ym = .FALSE.
    velocity_yp = .FALSE.

    flux_xm = .FALSE.
    flux_xp = .FALSE.
    flux_ym = .FALSE.
    flux_yp = .FALSE.

    poise_xm = .FALSE.
    poise_xp = .FALSE.
    poise_ym = .FALSE.
    poise_yp = .FALSE.

    xp3_ave = 3.8e-3
    xm3_ave = 3.8e-3
    xp3_max = 0.d0
    xm3_max = 0.d0
    xp3_ave_p = 0.d0
    xm3_ave_p = 0.d0
    xp3_ave_f = 0.d0
    xm3_ave_f = 0.d0

    yp3_ave = 3.8e-3
    ym3_ave = 3.8e-3
    yp3_max = 0.d0
    ym3_max = 0.d0
    yp3_ave_p = 0.d0
    ym3_ave_p = 0.d0
    yp3_ave_f = 0.d0
    ym3_ave_f = 0.d0

    if (help) call PetscPrintf(options%comm, "-bc_pressure_{xy}{mp}_component*: "// &
         "density of component * for a Dirichlet BC\n", ierr)
    if (help) call PetscPrintf(options%comm, "-bc_velocity_{xy}{mp}_avg: "// &
         "mean velocity for constant velocity, or max velocity for Poiseuille flow\n", ierr)
    if (help) call PetscPrintf(options%comm, "-bc_flux_{xy}{mp}_component*: "// &
         "normal volumetric flux of component *\n", ierr)

    ! get average values on minus edge
    if (bc_flags(BOUNDARY_XM).eq.BC_VELOCITY) then
       ! check if poiseuille velocity
       velocity_xm = .TRUE.
       call PetscOptionsGetBool(options%my_prefix,'-bc_velocity_xm_poiseuille', &
            poise_xm, flag, ierr)
       call PetscOptionsGetReal(options%my_prefix,'-bc_velocity_xm_avg', xm3_ave, &
            flag, ierr)
       if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)

    else if (bc_flags(BOUNDARY_XM).eq.BC_FLUX) then
       ! check if poiseuille velocity
       flux_xm = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_flux_xm_component'//char(m+48),&
               xm3_ave_f(m,X_DIRECTION), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) then
             call PetscOptionsGetReal(options%my_prefix,&
                  '-bc_flux_xm', xm3_ave_f(1,X_DIRECTION), flag, ierr)
          end if
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do

    else if (bc_flags(BOUNDARY_XM).eq.BC_DIRICHLET) then
       pressure_xm = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_xm_component'//char(m+48),&
               xm3_ave_p(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_pressure_xm', xm3_ave_p(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif

    ! get average values on minus edge
    if (bc_flags(BOUNDARY_XP).eq.BC_VELOCITY) then
       ! check if poiseuille velocity
       velocity_xp = .TRUE.
       call PetscOptionsGetBool(options%my_prefix,'-bc_velocity_xp_poiseuille', &
            poise_xp, flag, ierr)
       call PetscOptionsGetReal(options%my_prefix,'-bc_velocity_xp_avg', xp3_ave, &
            flag, ierr)
       if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
        
    else if (bc_flags(BOUNDARY_XP).eq.BC_FLUX) then
       ! check if poiseuille velocity
       flux_xp = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_flux_xp_component'//char(m+48),&
               xp3_ave_f(m,X_DIRECTION), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) then
             call PetscOptionsGetReal(options%my_prefix,&
                  '-bc_flux_xp', xp3_ave_f(1,X_DIRECTION), flag, ierr)
          end if
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do

    else if (bc_flags(BOUNDARY_XP).eq.BC_DIRICHLET) then
       pressure_xp = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_xp_component'//char(m+48),&
               xp3_ave_p(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_pressure_xp', xp3_ave_p(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif

    ! get average values on minus edge
    if (bc_flags(BOUNDARY_YM).eq.BC_VELOCITY) then
       ! check if poiseuille velocity
       velocity_ym = .TRUE.
       call PetscOptionsGetBool(options%my_prefix,'-bc_velocity_ym_poiseuille', &
            poise_ym, flag, ierr)
       call PetscOptionsGetReal(options%my_prefix,'-bc_velocity_ym_avg', ym3_ave, &
            flag, ierr)
       if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
        
    else if (bc_flags(BOUNDARY_YM).eq.BC_FLUX) then
       ! check if poiseuille velocity
       flux_ym = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_flux_ym_component'//char(m+48),&
               ym3_ave_f(m,Y_DIRECTION), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) then
             call PetscOptionsGetReal(options%my_prefix,&
                  '-bc_flux_ym', ym3_ave_f(1,Y_DIRECTION), flag, ierr)
          end if
       end do

    else if (bc_flags(BOUNDARY_YM).eq.BC_DIRICHLET) then
       pressure_ym = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_ym_component'//char(m+48),&
               ym3_ave_p(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_pressure_ym', ym3_ave_p(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif

    ! get average values on minus edge
    if (bc_flags(BOUNDARY_YP).eq.BC_VELOCITY) then
       velocity_yp = .TRUE.
       ! check if poiseuille velocity
       call PetscOptionsGetBool(options%my_prefix,'-bc_velocity_yp_poiseuille', &
            poise_yp, flag, ierr)
       call PetscOptionsGetReal(options%my_prefix,'-bc_velocity_yp_avg', yp3_ave, &
            flag, ierr)
       if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
        
    else if (bc_flags(BOUNDARY_YP).eq.BC_FLUX) then
       ! check if poiseuille velocity
       flux_yp = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_flux_yp_component'//char(m+48),&
               yp3_ave_f(m,Y_DIRECTION), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) then
             call PetscOptionsGetReal(options%my_prefix,&
                  '-bc_flux_yp', yp3_ave_f(1,Y_DIRECTION), flag, ierr)
          end if
       end do

    else if (bc_flags(BOUNDARY_YP).eq.BC_DIRICHLET) then
       pressure_yp = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_pressure_yp_component'//char(m+48),&
               yp3_ave_p(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_pressure_yp', yp3_ave_p(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif


! --- xm boundary
    if (dist%info%xs.eq.1) then
       xm_bcvals = 0.0        ! zero out, setting x and y velocityes on z boundary =0
       if (pressure_xm) then
          do m=1,dist%s
             xm_bcvals(m,1,:) = xm3_ave_p(m)
          end do
       else if (flux_xm) then
          do m=1,dist%s
             xm_bcvals(m,X_DIRECTION,:) = xm3_ave_f(m,X_DIRECTION)
          end do
       else if (poise_xm) then
          ! unclear what is intended here... for now just guessing, and doing poise in y
          xm3_max = 3./2.*xm3_ave
          do j=dist%info%ys,dist%info%ye
             xm_bcvals(1,X_DIRECTION,j) = 4.*xm3_max/(dist%info%NY-2.)**2*(j-1.5)*(dist%info%NY-0.5-j)
          enddo
       else if (velocity_xm) then
          xm_bcvals(1,X_DIRECTION,:) = xm3_ave
       endif
    end if

! --- xm boundary
    if (dist%info%xe.eq.dist%info%NX) then
       xp_bcvals = 0.0        ! xero out, setting x and y velocityes on x boundary =0
       if (pressure_xp) then
          do m=1,dist%s
             xp_bcvals(m,1,:) = xp3_ave_p(m)
          end do
       else if (flux_xp) then
          do m=1,dist%s
             xp_bcvals(m,X_DIRECTION,:) = xp3_ave_f(m,X_DIRECTION)
          end do
       else if (poise_xp) then
          ! unclear what is intended here... for now just guessing, and doing poise in y
          xp3_max = 3./2.*xp3_ave
          do j=dist%info%ys,dist%info%ye
             xp_bcvals(1,X_DIRECTION,j) = 4.*xp3_max/(dist%info%NY-2.)**2*(j-1.5)*(dist%info%NY-0.5-j)
          enddo
       else if (velocity_xp) then
          xp_bcvals(1,X_DIRECTION,:) = xp3_ave
       endif
    end if

! --- ym boundary
    if (dist%info%ys.eq.1) then
       ym_bcvals = 0.0        ! yero out, setting x and y velocityes on y boundary =0
       if (pressure_ym) then
          do m=1,dist%s
             ym_bcvals(m,1,:) = ym3_ave_p(m)
          end do
       else if (flux_ym) then
          do m=1,dist%s
             ym_bcvals(m,Y_DIRECTION,:) = ym3_ave_f(m,Y_DIRECTION)
          end do
       else if (poise_ym) then
          ym3_max = 3./2.*ym3_ave
          ! unclear what is intended here... for now just guessing, and doing poise in x
          do i=dist%info%xs,dist%info%xe
             ym_bcvals(1,Y_DIRECTION,i) = 4.*ym3_max/(dist%info%NX-2.)**2*(i-1.5)*(dist%info%NX-0.5-i)
          enddo
       else if (velocity_ym) then
          ym_bcvals(1,Y_DIRECTION,:) = ym3_ave
       endif
    end if

! --- ym boundary
    if (dist%info%ye.eq.dist%info%NY) then
       yp_bcvals = 0.0        ! yero out, setting x and y velocityes on y boundary =0
       if (pressure_yp) then
          do m=1,dist%s
             yp_bcvals(m,1,:) = yp3_ave_p(m)
          end do
       else if (flux_yp) then
          do m=1,dist%s
             yp_bcvals(m,Y_DIRECTION,:) = yp3_ave_f(m,Y_DIRECTION)
          end do
       else if (poise_yp) then
          yp3_max = 3./2.*yp3_ave
          ! unclear what is intended here... for now just guessing, and doing poise in x
          do i=dist%info%xs,dist%info%xe
             yp_bcvals(1,Y_DIRECTION,i) = 4.*yp3_max/(dist%info%NX-2.)**2*(i-1.5)*(dist%info%NX-0.5-i)
          enddo
       else if (velocity_yp) then
          yp_bcvals(1,Y_DIRECTION,:) = yp3_ave
       endif
    end if
    return
  end subroutine initialize_bcs
  
