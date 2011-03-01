!!!==================================================================
!!! Fortran-file
!!!    author:        Ethan T. Coon
!!!    filename:      get_forces.f
!!!    version:
!!!    created:       09 November 2010
!!!      on:          16:27:53 MST
!!!    last modified:  09 November 2010
!!!      at:          16:27:53 MST
!!!    URL:           http://www.ldeo.columbia.edu/~ecoon/
!!!    email:         ecoon _at_ ldeo.columbia.edu
!!!
!!!==================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

module LBM_Forcing_module
  use LBM_Info_module
  use LBM_Constants_module
  use petsc
  implicit none

  private 
#include "lbm_definitions.h"

  public:: LBMAddFluidFluidForces, &
       LBMAddBodyForces, &
       LBMAddFluidSolidForces, &
       LBMZeroBoundaryForces

contains
  ! --- Fluid-fluid interaction forces, from
  ! ---  (Kang 2002 Eq. 6)
  subroutine LBMAddFluidFluidForces(rho, forces, walls, info, constants)
    use LBM_Discretization_D3Q19_module
    !     NONLOCAL IN RHO

    ! input
    type(info_type) info
    type(constants_type) constants
    PetscScalar,dimension(1:info%s, 1:info%gxyzl):: rho
    PetscScalar,dimension(1:info%s, 1:info%dim, 1:info%gxyzl):: forces
    PetscScalar,dimension(1:info%gxyzl):: walls

    ! local
    PetscInt i,j,k,m,n,d
    PetscScalar,dimension(1:info%s, 0:info%b, 1:info%gxyzl)::tmp
    PetscErrorCode ierr

    do m=1,info%s
       call InfoGatherValueToDirection(info, rho(m,:), tmp(m,:,:))
    end do

    do i=1,info%gxyzl
       if (walls(i).eq.0) then
          do d=1,info%dim
             do n=1,2*info%dim
                forces(1,d,i) = forces(1,d,i) &
                     - rho(1,i)*tmp(2,n,i)*info%ci_int(n,d)*constants%g
                forces(1,d,i) = forces(1,d,i) &
                     - rho(1,i)*tmp(1,n,i)*info%ci_int(n,d)*constants%g11
                forces(2,d,i) = forces(2,d,i) &
                     - rho(2,i)*tmp(1,n,i)*info%ci_int(n,d)*constants%g
                forces(2,d,i) = forces(2,d,i) &
                     - rho(2,i)*tmp(2,n,i)*info%ci_int(n,d)*constants%g22
             enddo
             do n=2*info%dim+1,info%b
                forces(1,d,i) = forces(1,d,i) &
                     - rho(1,i)*tmp(2,n,i)*info%ci_int(n,d)*constants%g*0.5
                forces(1,d,i) = forces(1,d,i) &
                     - rho(1,i)*tmp(1,n,i)*info%ci_int(n,d)*constants%g11*0.5
                forces(2,d,i) = forces(2,d,i) &
                     - rho(2,i)*tmp(1,n,i)*info%ci_int(n,d)*constants%g*0.5
                forces(2,d,i) = forces(2,d,i) &
                     - rho(2,i)*tmp(2,n,i)*info%ci_int(n,d)*constants%g22*0.5
             enddo
          end do
       end if
    end do
  end subroutine LBMAddFluidFluidForces


  ! --- body forces on fluid
  subroutine LBMAddBodyForces(rho,forces,walls,info,constants)
    type(info_type) info
    type(constants_type) constants
    PetscScalar,dimension(1:info%s, 1:info%gxyzl):: rho
    PetscScalar,dimension(1:info%s, 1:info%dim, 1:info%gxyzl):: forces
    PetscScalar,dimension(1:info%gxyzl):: walls

    PetscInt i,m,n,d

    do m=1,info%s
       do i=1,info%gxyzl
          if (walls(i).eq.0) then
             do d=1,info%dim
                forces(m,d,i) = forces(m,d,i) &
                     + constants%gvt(m,d)*constants%mm(m)*rho(m,i)
             end do
          end if
       end do
    enddo
    return
  end subroutine LBMAddBodyForces

  ! --- Fluid-solid interaction forces, from
  ! ---  (Kang 2002 Eq. 8)
  ! --- NOTE -- I'm ignoring the bug that the forces will
  !             wrap even in the non-periodic case.  This seems
  !             save at this point, since it seems unlikely that
  !             there will be walls on one side but not the other
  subroutine LBMAddFluidSolidForces(rho, forces, walls, info, constants)
    use LBM_Discretization_D3Q19_module
    !     NONLOCAL IN WALLS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! New code added by MLP.  This accounts for fluid-solid forces on the 
    ! diagonals, which is not accounted for in the other version.
 
    ! input
    type(info_type) info
    type(constants_type) constants
    PetscScalar,dimension(1:info%s, 1:info%gxyzl):: rho
    PetscScalar,dimension(1:info%s, 1:info%dim, 1:info%gxyzl):: forces
    PetscScalar,dimension(1:info%gxyzl):: walls

    ! local
    PetscInt i,m,n,d
    PetscScalar,dimension(0:info%b,1:info%gxyzl)::tmp
    PetscErrorCode ierr

    call InfoGatherValueToDirection(info, walls, tmp)

    do m=1,info%s
      do i=1,info%gxyzl
        if (walls(i).eq.0) then
          do d=1,info%dim
            do n=1,2*info%dim
               forces(m,d,i) = forces(m,d,i) &
                    - rho(m,i)*tmp(n,i)*info%ci_int(n,d)*constants%gw(m)
            enddo
            do n=2*info%dim+1,info%b
               forces(m,d,i) = forces(m,d,i) &
                    - rho(m,i)*tmp(n,i)*info%ci_int(n,d)*constants%gw(m)*0.5
            enddo
          end do
        end if
      end do
    end do
  end subroutine LBMAddFluidSolidForces


  ! --- zero out boundaries, as they affect the fi
  subroutine LBMZeroBoundaryForces(bc_flags, forces, bc_dim, info,constants)
    type(info_type) info
    type(constants_type) constants
    PetscInt bc_dim
    PetscInt, dimension(6)::bc_flags ! enum for boundary conditions
    PetscScalar,dimension(1:info%s, 1:info%dim, 1:info%gxyzl):: forces

    if (info%dim.eq.2) then
       call LBMZeroBoundaryForcesD2(bc_flags, forces, bc_dim, info,constants)
    else if (info%dim.eq.3) then
       call LBMZeroBoundaryForcesD3(bc_flags, forces, bc_dim, info,constants)
    end if
  end subroutine LBMZeroBoundaryForces

  subroutine LBMZeroBoundaryForcesD2(bc_flags, forces, bc_dim, info,constants)
    type(info_type) info
    type(constants_type) constants
    PetscInt bc_dim
    PetscInt, dimension(6)::bc_flags ! enum for boundary conditions

    PetscScalar,dimension(1:info%s, 1:info%dim, info%gxs:info%gxe, &
         info%gys:info%gye):: forces
    
    ! -- x
    if ((info%xs.eq.1).and.((bc_flags(BOUNDARY_XM).eq.BC_FLUX).or.&
         (bc_flags(BOUNDARY_XM).eq.BC_PRESSURE))) then
       forces(:,:,1,:) = 0
    endif

    if ((info%xe.eq.info%NX).and.((bc_flags(BOUNDARY_XP).eq.BC_FLUX).or. &
         (bc_flags(BOUNDARY_XP).eq.BC_PRESSURE))) then
       forces(:,:,info%NX,:) = 0
    endif

    ! -- y
    if ((info%ys.eq.1).and.((bc_flags(BOUNDARY_YM).eq.BC_FLUX).or.&
         (bc_flags(BOUNDARY_YM).eq.BC_PRESSURE))) then
       forces(:,:,:,1) = 0
    endif

    if ((info%ye.eq.info%NY).and.((bc_flags(BOUNDARY_YP).eq.BC_FLUX).or. &
         (bc_flags(BOUNDARY_YP).eq.BC_PRESSURE))) then
       forces(:,:,:,info%NY) = 0
    endif

  end subroutine LBMZeroBoundaryForcesD2

  subroutine LBMZeroBoundaryForcesD3(bc_flags, forces, bc_dim, info,constants)
    type(info_type) info
    type(constants_type) constants
    PetscInt bc_dim
    PetscInt, dimension(6)::bc_flags ! enum for boundary conditions

    PetscScalar,dimension(1:info%s, 1:info%dim, info%gxs:info%gxe, &
         info%gys:info%gye, info%gzs:info%gze):: forces
    
    ! -- x
    if ((info%xs.eq.1).and.((bc_flags(BOUNDARY_XM).eq.BC_FLUX).or.&
         (bc_flags(BOUNDARY_XM).eq.BC_PRESSURE))) then
       forces(:,:,1,:,:) = 0
    endif

    if ((info%xe.eq.info%NX).and.((bc_flags(BOUNDARY_XP).eq.BC_FLUX).or. &
         (bc_flags(BOUNDARY_XP).eq.BC_PRESSURE))) then
       forces(:,:,info%NX,:,:) = 0
    endif

    ! -- y
    if ((info%ys.eq.1).and.((bc_flags(BOUNDARY_YM).eq.BC_FLUX).or.&
         (bc_flags(BOUNDARY_YM).eq.BC_PRESSURE))) then
       forces(:,:,:,1,:) = 0
    endif

    if ((info%ye.eq.info%NY).and.((bc_flags(BOUNDARY_YP).eq.BC_FLUX).or. &
         (bc_flags(BOUNDARY_YP).eq.BC_PRESSURE))) then
       forces(:,:,:,info%NY,:) = 0
    endif

    ! -- z
    if ((info%zs.eq.1).and.((bc_flags(BOUNDARY_ZM).eq.BC_FLUX).or.&
         (bc_flags(BOUNDARY_ZM).eq.BC_PRESSURE))) then
       forces(:,:,:,:,1) = 0
    endif

    if ((info%ze.eq.info%NZ).and.((bc_flags(BOUNDARY_ZP).eq.BC_FLUX).or. &
         (bc_flags(BOUNDARY_ZP).eq.BC_PRESSURE))) then
       forces(:,:,:,:,info%NZ) = 0
    endif
  end subroutine LBMZeroBoundaryForcesD3
end module LBM_Forcing_module
