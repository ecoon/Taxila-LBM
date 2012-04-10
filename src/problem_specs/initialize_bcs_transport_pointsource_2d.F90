!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_bcs_transport_constant_pointsource_2d.F90
!!!     version:         
!!!     created:         20 April 2011
!!!       on:            16:58:06 MDT
!!!     last modified:   05 May 2011
!!!       at:            10:18:47 MDT
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
    use LBM_Error_module
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
    PetscInt pointsource_node
    PetscBool flag, help
    PetscErrorCode ierr
    PetscInt m

    call PetscOptionsHasName(PETSC_NULL_CHARACTER, "-help", help, ierr)

    if (help) call PetscPrintf(options%comm, "-bc_conc_pointsource <val>: "// &
         "set the psi flux of primary specie in the inward normal direction", ierr)
    if (help) call PetscPrintf(options%comm, "-bc_conc_pointsource_location <val>: "// &
         "set the location of the point source in y-grid points", ierr)
    pointsource_conc = 0.

    call PetscOptionsGetReal(options%my_prefix,'-bc_conc_pointsource', pointsource_conc,&
         flag, ierr)
    if (.not.flag) then
       call LBMError(PETSC_COMM_SELF,1,'invalid boundary value for point source concentration', ierr)
    end if
    call PetscOptionsGetInt(options%my_prefix,'-bc_conc_pointsource_location', &
         pointsource_node, flag, ierr)
    if (.not.flag) then
       call LBMError(PETSC_COMM_SELF,1,'invalid node for point source concentration', ierr)
    end if

    if (dist%info%xs.eq.1) then
       xm_bcvals = 0.
       if (pointsource_node >= dist%info%ys .and. pointsource_node <= dist%info%ye) &
          xm_bcvals(1,X_DIRECTION,pointsource_node) = pointsource_conc
    end if

    if (dist%info%xe.eq.dist%info%NX) then
       xp_bcvals = 0.
    end if

    if (dist%info%ys.eq.1) then
       ym_bcvals = 0.
    end if

    if (dist%info%ye.eq.dist%info%NY) then
       yp_bcvals = 0.
    end if

  end subroutine initialize_bcs_transport
    
