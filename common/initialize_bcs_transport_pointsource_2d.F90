!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_bcs_transport_constant_pointsource_2d.F90
!!!     version:         
!!!     created:         20 April 2011
!!!       on:            16:58:06 MDT
!!!     last modified:   21 April 2011
!!!       at:            17:38:27 MDT
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
    PetscScalar pointsource_conc
    PetscBool flag, help
    PetscErrorCode ierr
    PetscInt m

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    if (help) call PetscPrintf(options%comm, "-bc_conc_pointsource':"// &
         "set the psi flux of primary specie in the inward normal direction", ierr)
    pointsource_conc = 0.d0

    call PetscOptionsGetReal(options%my_prefix,'-bc_conc_pointsource', pointsource_conc,&
         flag, ierr)
    if (.not.flag) then
       SETERRQ(1,1,'invalid boundary value for point source concentration', ierr)
    end if

    if (dist%info%xs.eq.1) then
       xm_bcvals = 0.d0
       if ((dist%info%NY+1)/2 >= dist%info%ys .and. (dist%info%NY+1)/2 <= dist%info%ye) &
          xm_bcvals(1,X_DIRECTION,(dist%info%NY+1)/2) = pointsource_conc
    end if

    if (dist%info%xe.eq.dist%info%NX) then
       xp_bcvals = 0.d0
    end if

    if (dist%info%ys.eq.1) then
       ym_bcvals = 0.d0
    end if

    if (dist%info%ye.eq.dist%info%NY) then
       yp_bcvals = 0.d0
    end if

  end subroutine initialize_bcs_transport
    
