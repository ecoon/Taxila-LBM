!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_bcs_transport_constant_pointsource_2d.F90
!!!     version:         
!!!     created:         20 April 2011
!!!       on:            16:58:06 MDT
!!!     last modified:   20 April 2011
!!!       at:            17:12:09 MDT
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
    PetscScalar,dimension(bc_dim, dist%info%ys:dist%info%ye):: xm_bcvals, xp_bcvals
    PetscScalar,dimension(bc_dim, dist%info%xs:dist%info%xe):: ym_bcvals, yp_bcvals
    PetscScalar,pointer:: zm_bcvals, zp_bcvals
    PetscInt bc_dim

    ! local
    PetscScalar pointsource_conc
    PetscBool flag, help
    PetscErrorCode ierr
    PetscInt m

    xm_bcvals = 0.d0
    xp_bcvals = 0.d0
    ym_bcvals = 0.d0
    yp_bcvals = 0.d0

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    if (help) call PetscPrintf(options%comm, "-bc_conc_pointsource{N}':"// &
         "set the psi of primary specie N at the point source", ierr)
    pointsource_conc = 0.d0
    do m=1,dist%s
       call PetscOptionsGetReal(options%my_prefix,'-bc_conc_pointsource'// &
            char(m+48), pointsource_conc, flag, ierr)
       if (.not.flag) then
          SETERRQ(1,1,'invalid boundary value for point source concentration', ierr)
       end if
       xm_bcvals(m,dist%info%NY/2) = pointsource_conc
    end do
  end subroutine initialize_bcs_transport
    
