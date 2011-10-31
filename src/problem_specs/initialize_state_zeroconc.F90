!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state_zeroconc.F90
!!!     version:         
!!!     created:         20 April 2011
!!!       on:            16:55:58 MDT
!!!     last modified:   31 October 2011
!!!       at:            15:54:22 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

  subroutine initialize_state_transport(rho, u, walls, &
       dist, components, options)
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
      call initialize_state_transport_d2(rho, u, walls, &
           dist, components, options)
    case (3) 
      call initialize_state_transport_d3(rho, u, walls, &
           dist, components, options)
    end select
  end subroutine initialize_state

  subroutine initialize_state_transport_d3(rho, u, walls, dist, components, options)
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

    ! initialize state
    u = 0.0
    rho = 0.0
    return
  end subroutine initialize_state_transport_d3


  subroutine initialize_state_transport_d2(rho, u, walls, dist, components, options)
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
    PetscScalar,dimension(dist%s, &
         dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: rho
    PetscScalar,dimension(dist%s, 1:dist%info%ndims, &
         dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: u
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls

    ! initialize state
    u = 0.0
    rho = 0.0
    return
  end subroutine initialize_state_transport_d2
