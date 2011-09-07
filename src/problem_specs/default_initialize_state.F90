!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        initialize_state.F90
!!!     version:         
!!!     created:         14 January 2011
!!!       on:            18:21:06 MST
!!!     last modified:   17 August 2011
!!!       at:            17:16:42 MDT
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
    PetscScalar,dimension(dist%info%gxyzl) :: walls

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
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: walls

    ! local variables
    PetscErrorCode ierr
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, &
         dist%info%gzs:dist%info%gze):: nowalls

    ! initialize state
    fi = 0.0
    u = 0.0
    rho = 1.d0
    
    ! set state at equilibrium       
    nowalls = 0.d0
    call DiscretizationEquilf(dist%disc, rho, u, nowalls, fi, 1, components(1)%relax, dist)    
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
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls

    ! local variables
    PetscErrorCode ierr
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: nowalls

    ! initialize state
    fi = 0.d0
    u = 0.d0
    rho = 1.d0
    
    ! set state at equilibrium       
    nowalls = 0.d0
    call DiscretizationEquilf(dist%disc, rho, u, nowalls, fi, 1, components(1)%relax, dist)    
    return
  end subroutine initialize_state_d2
