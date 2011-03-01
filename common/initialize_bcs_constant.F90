!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_bc_zfluxflux.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            17:30:22 MST
!!!     last modified:   25 February 2011
!!!       at:            09:17:00 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

! initializes the BCs from options to be a general constant flux/pressure bc
  subroutine initialize_bcs(bc_flags, xm_bcvals, xp_bcvals, ym_bcvals, & 
       yp_bcvals, zm_bcvals, zp_bcvals, bc_dim, info)
    use petsc
    use LBM_Info_module
    implicit none

#include "lbm_definitions.h"    
    PetscInt, dimension(6):: bc_flags ! enum for boundary conditions

    type(info_type) info
    PetscScalar,dimension(bc_dim, info%ys:info%ye, info%zs:info%ze):: xm_bcvals, xp_bcvals
    PetscScalar,dimension(bc_dim, info%xs:info%xe, info%zs:info%ze):: ym_bcvals, yp_bcvals
    PetscScalar,dimension(bc_dim, info%xs:info%xe, info%ys:info%ye):: zm_bcvals, zp_bcvals
    PetscInt bc_dim

    ! local
    PetscScalar xp3_ave, xm3_ave, xp3_max, xm3_max
    PetscScalar,dimension(1:info%s):: xp3_ave_p, xm3_ave_p
    PetscBool pressure_xm, pressure_xp, flux_xm, flux_xp, poise_xm, poise_xp

    PetscScalar yp3_ave, ym3_ave, yp3_max, ym3_max
    PetscScalar,dimension(1:info%s):: yp3_ave_p, ym3_ave_p
    PetscBool pressure_ym, pressure_yp, flux_ym, flux_yp, poise_ym, poise_yp

    PetscScalar zp3_ave, zm3_ave, zp3_max, zm3_max
    PetscScalar,dimension(1:info%s):: zp3_ave_p, zm3_ave_p
    PetscBool pressure_zm, pressure_zp, flux_zm, flux_zp, poise_zm, poise_zp

    PetscInt i, j, k, m
    PetscBool flag
    PetscErrorCode ierr

    ! get options
    pressure_xm = .FALSE.
    pressure_xp = .FALSE.
    pressure_ym = .FALSE.
    pressure_yp = .FALSE.
    pressure_zm = .FALSE.
    pressure_zp = .FALSE.

    flux_xm = .FALSE.
    flux_xp = .FALSE.
    flux_ym = .FALSE.
    flux_yp = .FALSE.
    flux_zm = .FALSE.
    flux_zp = .FALSE.

    poise_xm = .FALSE.
    poise_xp = .FALSE.
    poise_ym = .FALSE.
    poise_yp = .FALSE.
    poise_zm = .FALSE.
    poise_zp = .FALSE.

    xp3_ave = 3.8e-3
    xm3_ave = 3.8e-3
    xp3_max = 0.d0
    xm3_max = 0.d0
    xp3_ave_p = 0.d0
    xm3_ave_p = 0.d0

    yp3_ave = 3.8e-3
    ym3_ave = 3.8e-3
    yp3_max = 0.d0
    ym3_max = 0.d0
    yp3_ave_p = 0.d0
    ym3_ave_p = 0.d0

    zp3_ave = 3.8e-3
    zm3_ave = 3.8e-3
    zp3_max = 0.d0
    zm3_max = 0.d0
    zp3_ave_p = 0.d0
    zm3_ave_p = 0.d0

    ! get average values on minus edge
    if (bc_flags(BOUNDARY_XM).eq.BC_FLUX) then
       ! check if poiseuille flux
       flux_xm = .TRUE.
       call PetscOptionsGetBool(info%options_prefix,'-bc_flux_xm_poiseuille', &
            poise_xm, flag, ierr)
       call PetscOptionsGetReal(info%options_prefix,'-bc_flux_xm_avg', xm3_ave, &
            flag, ierr)
       if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
        
    else if (bc_flags(BOUNDARY_XM).eq.BC_PRESSURE) then
       pressure_xm = .TRUE.
       do m=1,info%s
          call PetscOptionsGetReal(info%options_prefix,'-bc_pressure_xm_phase'//char(m+48),&
               xm3_ave_p(m), flag, ierr)
          if ((info%s.eq.1).and. .not.flag) call PetscOptionsGetReal(info%options_prefix,&
               '-bc_pressure_xm', xm3_ave_p(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif

    ! get average values on minus edge
    if (bc_flags(BOUNDARY_XP).eq.BC_FLUX) then
       ! check if poiseuille flux
       flux_xp = .TRUE.
       call PetscOptionsGetBool(info%options_prefix,'-bc_flux_xp_poiseuille', &
            poise_xp, flag, ierr)
       call PetscOptionsGetReal(info%options_prefix,'-bc_flux_xp_avg', xp3_ave, &
            flag, ierr)
       if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
        
    else if (bc_flags(BOUNDARY_XP).eq.BC_PRESSURE) then
       pressure_xp = .TRUE.
       do m=1,info%s
          call PetscOptionsGetReal(info%options_prefix,'-bc_pressure_xp_phase'//char(m+48),&
               xp3_ave_p(m), flag, ierr)
          if ((info%s.eq.1).and. .not.flag) call PetscOptionsGetReal(info%options_prefix,&
               '-bc_pressure_xp', xp3_ave_p(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif

    ! get average values on minus edge
    if (bc_flags(BOUNDARY_YM).eq.BC_FLUX) then
       ! check if poiseuille flux
       flux_ym = .TRUE.
       call PetscOptionsGetBool(info%options_prefix,'-bc_flux_ym_poiseuille', &
            poise_ym, flag, ierr)
       call PetscOptionsGetReal(info%options_prefix,'-bc_flux_ym_avg', ym3_ave, &
            flag, ierr)
       if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
        
    else if (bc_flags(BOUNDARY_YM).eq.BC_PRESSURE) then
       pressure_ym = .TRUE.
       do m=1,info%s
          call PetscOptionsGetReal(info%options_prefix,'-bc_pressure_ym_phase'//char(m+48),&
               ym3_ave_p(m), flag, ierr)
          if ((info%s.eq.1).and. .not.flag) call PetscOptionsGetReal(info%options_prefix,&
               '-bc_pressure_ym', ym3_ave_p(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif

    ! get average values on minus edge
    if (bc_flags(BOUNDARY_YP).eq.BC_FLUX) then
       flux_yp = .TRUE.
       ! check if poiseuille flux
       call PetscOptionsGetBool(info%options_prefix,'-bc_flux_yp_poiseuille', &
            poise_yp, flag, ierr)
       call PetscOptionsGetReal(info%options_prefix,'-bc_flux_yp_avg', yp3_ave, &
            flag, ierr)
       if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
        
    else if (bc_flags(BOUNDARY_YP).eq.BC_PRESSURE) then
       pressure_yp = .TRUE.
       do m=1,info%s
          call PetscOptionsGetReal(info%options_prefix,'-bc_pressure_yp_phase'//char(m+48),&
               yp3_ave_p(m), flag, ierr)
          if ((info%s.eq.1).and. .not.flag) call PetscOptionsGetReal(info%options_prefix,&
               '-bc_pressure_yp', yp3_ave_p(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif

    ! get average values on minus edge
    if (bc_flags(BOUNDARY_ZM).eq.BC_FLUX) then
       flux_zm = .TRUE.
       ! check if poiseuille flux
       call PetscOptionsGetBool(info%options_prefix,'-bc_flux_zm_poiseuille', &
            poise_zm, flag, ierr)
       call PetscOptionsGetReal(info%options_prefix,'-bc_flux_zm_avg', zm3_ave, &
            flag, ierr)
       if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
        
    else if (bc_flags(BOUNDARY_ZM).eq.BC_PRESSURE) then
       pressure_zm = .TRUE.
       do m=1,info%s
          call PetscOptionsGetReal(info%options_prefix,'-bc_pressure_zm_phase'//char(m+48),&
               zm3_ave_p(m), flag, ierr)
          if ((info%s.eq.1).and. .not.flag) call PetscOptionsGetReal(info%options_prefix,&
               '-bc_pressure_zm', zm3_ave_p(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif

    ! get average values on minus edge
    if (bc_flags(BOUNDARY_ZP).eq.BC_FLUX) then
       flux_zp = .TRUE.
       ! check if poiseuille flux
       call PetscOptionsGetBool(info%options_prefix,'-bc_flux_zp_poiseuille', &
            poise_zp, flag, ierr)
       call PetscOptionsGetReal(info%options_prefix,'-bc_flux_zp_avg', zp3_ave, &
            flag, ierr)
       if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
        
    else if (bc_flags(BOUNDARY_ZP).eq.BC_PRESSURE) then
       pressure_zp = .TRUE.
       do m=1,info%s
          call PetscOptionsGetReal(info%options_prefix,'-bc_pressure_zp_phase'//char(m+48),&
               zp3_ave_p(m), flag, ierr)
          if ((info%s.eq.1).and. .not.flag) call PetscOptionsGetReal(info%options_prefix,&
               '-bc_pressure_zp', zp3_ave_p(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif


! --- xm boundary
    if (info%xs.eq.1) then
       xm_bcvals = 0.0        ! zero out, setting x and y fluxes on z boundary =0
       if (pressure_xm) then
          do m=1,info%s
             xm_bcvals(m,:,:) = xm3_ave_p(m)
          end do
       else if (poise_xm) then
          ! unclear what is intended here... for now just guessing, and doing poise in y
          xm3_max = 3./2.*xm3_ave
          do k=info%zs,info%ze
             do j=info%ys,info%ye
                xm_bcvals(1,j,k) = 4.*xm3_max/(info%NY-2.)**2*(j-1.5)*(info%NY-0.5-j)
             enddo
          enddo
       else if (flux_xm) then
          xm_bcvals(1,:,:) = xm3_ave
       endif
    end if

! --- xm boundary
    if (info%xe.eq.info%NX) then
       xp_bcvals = 0.0        ! xero out, setting x and y fluxes on x boundary =0
       if (pressure_xp) then
          do m=1,info%s
             xp_bcvals(m,:,:) = xp3_ave_p(m)
          end do
       else if (poise_xp) then
          ! unclear what is intended here... for now just guessing, and doing poise in y
          xp3_max = 3./2.*xp3_ave
          do k=info%zs,info%ze
             do j=info%ys,info%ye
                xp_bcvals(1,j,k) = 4.*xp3_max/(info%NY-2.)**2*(j-1.5)*(info%NY-0.5-j)
             enddo
          enddo
       else if (flux_xp) then
          xp_bcvals(1,:,:) = xp3_ave
       endif
    end if

! --- ym boundary
    if (info%ys.eq.1) then
       ym_bcvals = 0.0        ! yero out, setting x and y fluxes on y boundary =0
       if (pressure_ym) then
          do m=1,info%s
             ym_bcvals(m,:,:) = ym3_ave_p(m)
          end do
       else if (poise_ym) then
          ym3_max = 3./2.*ym3_ave
          ! unclear what is intended here... for now just guessing, and doing poise in x
          do i=info%xs,info%xe
             do k=info%zs,info%ze
                ym_bcvals(2,i,k) = 4.*ym3_max/(info%NX-2.)**2*(i-1.5)*(info%NX-0.5-i)
             enddo
          enddo
       else if (flux_ym) then
          ym_bcvals(2,:,:) = ym3_ave
       endif
    end if

! --- ym boundary
    if (info%ye.eq.info%NY) then
       yp_bcvals = 0.0        ! yero out, setting x and y fluxes on y boundary =0
       if (pressure_yp) then
          do m=1,info%s
             yp_bcvals(m,:,:) = yp3_ave_p(m)
          end do
       else if (poise_yp) then
          yp3_max = 3./2.*yp3_ave
          ! unclear what is intended here... for now just guessing, and doing poise in x
          do i=info%xs,info%xe
             do k=info%zs,info%ze
                yp_bcvals(2,i,k) = 4.*yp3_max/(info%NX-2.)**2*(i-1.5)*(info%NX-0.5-i)
             enddo
          enddo
       else if (flux_yp) then
          yp_bcvals(2,:,:) = yp3_ave
       endif
    end if

! --- zm boundary
    if (info%zs.eq.1) then
       zm_bcvals = 0.0        ! zero out, setting x and y fluxes on z boundary =0
       if (pressure_zm) then
          do m=1,info%s
             zm_bcvals(m,:,:) = zm3_ave_p(m)
          end do
       else if (poise_zm) then
          zm3_max = 3./2.*zm3_ave
          do i=info%xs,info%xe
             do j=info%ys,info%ye
                zm_bcvals(3,i,j) = 4.*zm3_max/(info%NY-2.)**2*(j-1.5)*(info%NY-0.5-j)
             enddo
          enddo
       else if (flux_zm) then
          zm_bcvals(3,:,:) = zm3_ave
       endif
    end if

! --- zm boundary
    if (info%ze.eq.info%NZ) then
       zp_bcvals = 0.0        ! zero out, setting x and y fluxes on z boundary =0
       if (pressure_zp) then
          do m=1,info%s
             zp_bcvals(m,:,:) = zp3_ave_p(m)
          end do
       else if (poise_zp) then
          zp3_max = 3./2.*zp3_ave
          do i=info%xs,info%xe
             do j=info%ys,info%ye
                zp_bcvals(3,i,j) = 4.*zp3_max/(info%NY-2.)**2*(j-1.5)*(info%NY-0.5-j)
             enddo
          enddo
       else if (flux_zp) then
          zp_bcvals(3,:,:) = zp3_ave
       endif
    end if
    return
  end subroutine initialize_bcs
  
