!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_bc_zfluxflux.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            17:30:22 MST
!!!     last modified:   05 May 2011
!!!       at:            10:17:09 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

! initializes the BCs from options to be a general constant flux/pressure bc
  subroutine initialize_bcs_transport(bc_flags, xm_bcvals, xp_bcvals, ym_bcvals, & 
       yp_bcvals, zm_bcvals, zp_bcvals, bc_dim, dist, options)
    use petsc
    use LBM_Distribution_Function_type_module
    use LBM_Options_module
    implicit none

#include "lbm_definitions.h"    
    PetscInt, dimension(6):: bc_flags ! enum for boundary conditions

    type(distribution_type) dist
    type(options_type) options
    PetscScalar,dimension(dist%s,dist%info%ndims, dist%info%ys:dist%info%ye):: xm_bcvals, xp_bcvals
    PetscScalar,dimension(dist%s,dist%info%ndims, dist%info%xs:dist%info%xe):: ym_bcvals, yp_bcvals
    PetscScalar,pointer:: zm_bcvals, zp_bcvals
    PetscInt bc_dim

    ! local
    PetscScalar,dimension(dist%s):: xp3_ave, xm3_ave
    PetscScalar,dimension(dist%s):: yp3_ave, ym3_ave

    PetscBool conc_xm, conc_xp
    PetscBool conc_ym, conc_yp
    PetscBool flux_xm, flux_xp
    PetscBool flux_ym, flux_yp

    PetscInt i, j, m
    PetscBool flag
    PetscErrorCode ierr
    PetscBool help

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    ! get options
    conc_xm = .FALSE.
    conc_xp = .FALSE.
    conc_ym = .FALSE.
    conc_yp = .FALSE.

    flux_xm = .FALSE.
    flux_xp = .FALSE.
    flux_ym = .FALSE.
    flux_yp = .FALSE.

    xp3_ave = 0.d0
    xm3_ave = 0.d0
    yp3_ave = 0.d0
    ym3_ave = 0.d0

    if (help) call PetscPrintf(options%comm, "-bc_conc_{xy}{mp}_specie*: "// &
         "concentration of specie * for a Dirichlet BC\n", ierr)
    if (help) call PetscPrintf(options%comm, "-bc_conc_flux_{xy}{mp}_specie*: "// &
         "normal volumetric flux of specie *\n", ierr)

    ! get average values on xm edge
    if (bc_flags(BOUNDARY_XM).eq.BC_DIRICHLET) then
       conc_xm = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_conc_xm_specie'// &
               char(m+48), xm3_ave(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_conc_xm', xm3_ave(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    else if (bc_flags(BOUNDARY_XM).eq.BC_FLUX) then
       flux_xm = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_conc_flux_xm_specie'// &
               char(m+48), xm3_ave(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_conc_flux_xm', xm3_ave(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif

    ! get average values on xp edge
    if (bc_flags(BOUNDARY_XP).eq.BC_DIRICHLET) then
       conc_xp = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_conc_xp_specie'// &
               char(m+48), xp3_ave(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_conc_xp', xp3_ave(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    else if (bc_flags(BOUNDARY_XP).eq.BC_FLUX) then
       flux_xp = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_conc_flux_xp_specie'// &
               char(m+48), xp3_ave(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_conc_flux_xp', xp3_ave(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif

    ! get average values on ym edge
    if (bc_flags(BOUNDARY_YM).eq.BC_DIRICHLET) then
       conc_ym = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_conc_ym_specie'// &
               char(m+48), ym3_ave(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_conc_ym', ym3_ave(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    else if (bc_flags(BOUNDARY_YM).eq.BC_FLUX) then
       flux_ym = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_conc_flux_ym_specie'// &
               char(m+48), ym3_ave(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_conc_flux_ym', ym3_ave(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif

    ! get average values on yp edge
    if (bc_flags(BOUNDARY_YP).eq.BC_DIRICHLET) then
       conc_yp = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_conc_yp_specie'// &
               char(m+48), yp3_ave(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_conc_yp', yp3_ave(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    else if (bc_flags(BOUNDARY_YP).eq.BC_FLUX) then
       flux_yp = .TRUE.
       do m=1,dist%s
          call PetscOptionsGetReal(options%my_prefix,'-bc_conc_flux_yp_specie'// &
               char(m+48), yp3_ave(m), flag, ierr)
          if ((dist%s.eq.1).and. .not.flag) call PetscOptionsGetReal(options%my_prefix,&
               '-bc_conc_flux_yp', yp3_ave(1), flag, ierr)
          if (.not.flag) SETERRQ(1, 1, 'invalid boundary value', ierr)
       end do
    endif


! --- xm boundary
    if (dist%info%xs.eq.1) then
       xm_bcvals = 0.
       if (conc_xm) then
          do m=1,dist%s
             xm_bcvals(m,1,:) = xm3_ave(m)
          end do
       else if (flux_xm) then
          do m=1,dist%s
             xm_bcvals(m,X_DIRECTION,:) = xm3_ave(m)
          end do
       end if
    end if

! --- xp boundary
    if (dist%info%xe.eq.dist%info%NX) then
       xp_bcvals = 0.
       if (conc_xp) then
          do m=1,dist%s
             xp_bcvals(m,1,:) = xp3_ave(m)
          end do
       else if (flux_xp) then
          do m=1,dist%s
             xp_bcvals(m,X_DIRECTION,:) = xp3_ave(m)
          end do
       end if
    end if

! --- ym boundary
    if (dist%info%ys.eq.1) then
       ym_bcvals = 0.
       if (conc_ym) then
          do m=1,dist%s
             ym_bcvals(m,1,:) = ym3_ave(m)
          end do
       else if (flux_ym) then
          do m=1,dist%s
             ym_bcvals(m,Y_DIRECTION,:) = ym3_ave(m)
          end do
       end if
    end if

! --- yp boundary
    if (dist%info%ye.eq.dist%info%NY) then
       yp_bcvals = 0.
       if (conc_yp) then
          do m=1,dist%s
             yp_bcvals(m,1,:) = yp3_ave(m)
          end do
       else if (flux_yp) then
          do m=1,dist%s
             yp_bcvals(m,Y_DIRECTION,:) = yp3_ave(m)
          end do
       end if
    end if

    return
  end subroutine initialize_bcs_transport
  
