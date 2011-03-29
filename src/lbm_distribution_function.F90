!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_distribution_function.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            14:06:07 MDT
!!!     last modified:   29 March 2011
!!!       at:            16:19:49 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"

module LBM_Distribution_Function_module
  use petsc
  use LBM_Info_module
  use LBM_Discretization_Type_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: distribution_type
     MPI_Comm comm
     PetscInt s
     PetscInt b
     type(info_type),pointer:: info
     type(discretization_type),pointer:: disc
     DM,pointer :: da_fi, da_rho

     Vec fi
     Vec fi_g
     PetscScalar,pointer:: fi_a(:)

     Vec rho
     Vec rho_g
     PetscScalar,pointer:: rho_a(:)

     PetscScalar,pointer,dimension(:,:,:):: flux
     PetscBool flux_required
     character(len=MAXWORDLENGTH) name       
  end type distribution_type

  public:: DistributionCreate, &
       DistributionDestroy, &
       DistributionSetName, &
       DistributionSetInfo, &
       DistributionSetDiscretization, &
       DistributionSetSizes, &
       DistributionSetDAs, &
       DistributionCalcDensity, &
       DistributionCalcFlux

contains
  function DistributionCreate(comm) result(distribution)
    MPI_Comm comm
    type(distribution_type),pointer:: distribution
    allocate(distribution)
    distribution%comm = comm
    distribution%s = -1
    distribution%b = -1
    nullify(distribution%info)
    nullify(distribution%disc)
    nullify(distribution%da_fi)
    nullify(distribution%da_rho)

    distribution%fi = 0
    distribution%rho = 0
    distribution%fi_g = 0
    distribution%rho_g = 0
    nullify(distribution%fi_a)
    nullify(distribution%rho_a)
    nullify(distribution%flux)
    distribution%flux_required = PETSC_TRUE
  end function DistributionCreate

  subroutine DistributionDestroy(distribution, ierr)
    type(distribution_type) distribution
    PetscErrorCode ierr
    if (distribution%fi /= 0) call VecDestroy(distribution%fi,ierr)
    if (distribution%rho /= 0) call VecDestroy(distribution%rho,ierr)
    if (distribution%fi_g /= 0) call VecDestroy(distribution%fi_g,ierr)
    if (distribution%rho_g /= 0) call VecDestroy(distribution%rho_g,ierr)
    if (associated(distribution%flux)) deallocate(distribution%flux)
  end subroutine DistributionDestroy

  subroutine DistributionSetName(distribution, name) 
    type(distribution_type) distribution 
    character(len=MAXWORDLENGTH):: name       
    distribution%name = name
  end subroutine DistributionSetName

  subroutine DistributionSetInfo(distribution, info)
    type(distribution_type) distribution
    type(info_type),pointer:: info
    distribution%info => info
  end subroutine DistributionSetInfo

  subroutine DistributionSetDiscretization(distribution, disc)
    type(distribution_type) distribution
    type(discretization_type),pointer:: disc
    distribution%disc => disc
    distribution%b = disc%b
  end subroutine DistributionSetDiscretization

  subroutine DistributionSetSizes(distribution, s)
    type(distribution_type) distribution
    PetscInt s
    distribution%s = s
  end subroutine DistributionSetSizes
    
  subroutine DistributionSetDAs(distribution, da_fi, da_rho)
    type(distribution_type) distribution
    DM,target:: da_fi, da_rho
    distribution%da_fi => da_fi
    distribution%da_rho => da_rho
  end subroutine DistributionSetDAs

  subroutine DistributionSetUp(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    PetscScalar zero
    zero = 0.d0

    if (distribution%flux_required) then
       allocate(distribution%flux(1:distribution%s, 1:distribution%info%ndims, &
            1:distribution%info%gxyzl))
       distribution%flux = 0.d0
    end if

    call DMCreateLocalVector(distribution%da_fi, distribution%fi, ierr)
    call VecSet(distribution%fi, zero, ierr)

    call DMCreateGlobalVector(distribution%da_fi, distribution%fi_g, ierr)
    call VecSet(distribution%fi_g, zero, ierr)
    call PetscObjectSetName(distribution%rho_g, trim(distribution%name)//'fi', ierr)

    call DMCreateLocalVector(distribution%da_rho, distribution%rho, ierr)
    call VecSet(distribution%rho, zero, ierr)

    call DMCreateGlobalVector(distribution%da_rho, distribution%rho_g, ierr)
    call VecSet(distribution%rho_g, zero, ierr)
    call PetscObjectSetName(distribution%rho_g, trim(distribution%name)//'rho', ierr)
  end subroutine DistributionSetUp

  subroutine DistributionGetArrays(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DMDAVecGetArrayF90(distribution%da_rho, distribution%rho, distribution%rho_a, ierr)
    call DMDAVecGetArrayF90(distribution%da_fi, distribution%fi, distribution%fi_a, ierr)
  end subroutine DistributionGetArrays

  subroutine DistributionRestoreArrays(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DMDAVecRestoreArrayF90(distribution%da_rho,distribution%rho,distribution%rho_a,ierr)
    call DMDAVecRestoreArrayF90(distribution%da_fi,distribution%fi,distribution%fi_a,ierr)
  end subroutine DistributionRestoreArrays

  subroutine DistributionCommunicateAll(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DistributionRestoreArrays(distribution)
    call DMDALocalToLocalBegin(distribution%da_fi, distribution%fi, INSERT_VALUES, &
         distribution%fi, ierr)
    call DMDALocalToLocalEnd(distribution%da_fi, distribution%fi, INSERT_VALUES, &
         distribution%fi, ierr)
    call DMDALocalToLocalBegin(distribution%da_rho, distribution%rho, INSERT_VALUES, &
         distribution%rho, ierr)
    call DMDALocalToLocalEnd(distribution%da_rho, distribution%rho, INSERT_VALUES, &
         distribution%rho, ierr)
    call DistributionGetArrays(distribution)
  end subroutine DistributionCommunicateAll

  subroutine DistributionCommunicateFi(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DMDAVecRestoreArrayF90(distribution%da_fi, distribution%fi, &
         distribution%fi_a, ierr)
    call DMDALocalToLocalBegin(distribution%da_fi, distribution%fi, INSERT_VALUES, &
         distribution%fi, ierr)
    call DMDALocalToLocalEnd(distribution%da_fi, distribution%fi, INSERT_VALUES, &
         distribution%fi, ierr)
    call DMDAVecGetArrayF90(distribution%da_fi, distribution%fi, distribution%fi_a, ierr)
  end subroutine DistributionCommunicateFi

  subroutine DistributionCommunicateRho(distribution)
    type(distribution_type) distribution
    PetscErrorCode ierr
    call DMDAVecRestoreArrayF90(distribution%da_rho, distribution%rho, &
         distribution%rho_a, ierr)
    call DMDALocalToLocalBegin(distribution%da_rho, distribution%rho, INSERT_VALUES, &
         distribution%rho, ierr)
    call DMDALocalToLocalEnd(distribution%da_rho, distribution%rho, INSERT_VALUES, &
         distribution%rho, ierr)
    call DMDAVecGetArrayF90(distribution%da_rho, distribution%rho, distribution%rho_a, ierr)
  end subroutine DistributionCommunicateRho

  subroutine DistributionCalcDensity(distribution, walls)
    type(distribution_type) distribution
    PetscScalar,dimension(1:distribution%info%gxyzl):: walls

    select case(distribution%info%ndims)
    case(2)
       call DistributionCalcDensityD2(distribution, distribution%fi_a, walls, &
            distribution%rho_a)
    case(3)
       call DistributionCalcDensityD3(distribution, distribution%fi_a, walls, & 
            distribution%rho_a)
    end select
  end  subroutine DistributionCalcDensity

  subroutine DistributionCalcDensityD2(distribution, fi, walls, rho)
    type(distribution_type) distribution
    PetscScalar,dimension(1:distribution%s, 0:distribution%b, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye):: fi
    PetscScalar,dimension(distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye):: walls
    PetscScalar,dimension(distribution%s, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye):: rho
    
    PetscInt i,j,m
    do j=distribution%info%ys,distribution%info%ye
    do i=distribution%info%xs,distribution%info%xe
       if (walls(i,j).eq.0) then
          rho(:,i,j) = sum(fi(:,:,i,j),2)
       end if
    end do
    end do
  end  subroutine DistributionCalcDensityD2

  subroutine DistributionCalcDensityD3(distribution, fi, walls, rho)
    type(distribution_type) distribution
    PetscScalar,dimension(1:distribution%s, 0:distribution%b, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye, &
         distribution%info%gzs:distribution%info%gze):: fi
    PetscScalar,dimension(distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye, &
         distribution%info%gzs:distribution%info%gze):: walls
    PetscScalar,dimension(distribution%s, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye, &
         distribution%info%gzs:distribution%info%gze):: rho
    
    PetscInt i,j,k,m
    do k=distribution%info%zs,distribution%info%ze
    do j=distribution%info%ys,distribution%info%ye
    do i=distribution%info%xs,distribution%info%xe
       if (walls(i,j,k).eq.0) then
          rho(:,i,j,k) = sum(fi(:,:,i,j,k),2)
       else
          rho(:,i,j,k) = 0.d0
       end if
    end do
    end do
    end do
  end  subroutine DistributionCalcDensityD3

  subroutine DistributionCalcFlux(distribution, walls)
    type(distribution_type) distribution
    PetscScalar,dimension(distribution%info%gxyzl):: walls

    select case(distribution%info%ndims)
    case(2)
       call DistributionCalcFluxD2(distribution,distribution%fi_a,walls,distribution%flux)
    case(3)
       call DistributionCalcFluxD3(distribution,distribution%fi_a,walls,distribution%flux)
    end select
  end  subroutine DistributionCalcFlux

  subroutine DistributionCalcFluxD3(distribution, fi, walls, u)
    type(distribution_type) distribution
    PetscScalar,dimension(distribution%s, 0:distribution%b, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye, &
         distribution%info%gzs:distribution%info%gze):: fi
    PetscScalar,dimension(distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye, &
         distribution%info%gzs:distribution%info%gze):: walls
    PetscScalar,dimension(distribution%s, distribution%info%ndims, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye, &
         distribution%info%gzs:distribution%info%gze):: u
    
    PetscInt i,j,k,d,m
    do k=distribution%info%zs,distribution%info%ze
    do j=distribution%info%ys,distribution%info%ye
    do i=distribution%info%xs,distribution%info%xe
       if (walls(i,j,k).eq.0) then
          do m=1,distribution%s
          do d=1,distribution%info%ndims
             u(m,d,i,j,k) = sum(dble(distribution%disc%ci(:,d))*fi(m,:,i,j,k),1)
          end do
          end do
       else
          u(:,:,i,j,k) = 0.d0
       end if
    end do
    end do
    end do
  end subroutine DistributionCalcFluxD3

  subroutine DistributionCalcFluxD2(distribution, fi, walls, u)
    type(distribution_type) distribution
    PetscScalar,dimension(distribution%s, 0:distribution%b, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye):: fi
    PetscScalar,dimension(distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye):: walls
    PetscScalar,dimension(distribution%s, distribution%info%ndims, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye):: u
    
    PetscInt i,j,d,m
    do j=distribution%info%ys,distribution%info%ye
    do i=distribution%info%xs,distribution%info%xe
       if (walls(i,j).eq.0) then
          do m=1,distribution%s
          do d=1,distribution%info%ndims
             u(m,d,i,j) = sum(dble(distribution%disc%ci(:,d))*fi(m,:,i,j),1)
          end do
          end do
       else
          u(:,:,i,j) = 0.d0
       end if
    end do
    end do
  end subroutine DistributionCalcFluxD2

  subroutine DistributionGatherValueToDirection(distribution, val, out)
    type(distribution_type) distribution
    PetscScalar,intent(in),dimension(1:distribution%info%gxyzl):: val
    PetscScalar,intent(out),dimension(0:distribution%b, 1:distribution%info%gxyzl):: out
    PetscErrorCode ierr
    
    if (distribution%info%ndims.eq.2) then
       call DistributionGatherValueToDirectionD2(distribution, val, out)
    else if (distribution%info%ndims.eq.3) then
       call DistributionGatherValueToDirectionD3(distribution, val, out)
    else 
       SETERRQ(1, 1, 'invalid ndims in LBM', ierr)
    end if
  end subroutine DistributionGatherValueToDirection
  
  subroutine DistributionGatherValueToDirectionD2(distribution, val, out)
    type(distribution_type) distribution
    PetscScalar,intent(in),dimension(distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye):: val
    PetscScalar,intent(out),dimension(0:distribution%b, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye):: out

    PetscInt n
    do n=0,distribution%b
       out(n,distribution%info%xs:distribution%info%xe, &
            distribution%info%ys:distribution%info%ye) = &
              val(distribution%info%xs + distribution%disc%ci(n,X_DIRECTION): &
                  distribution%info%xe + distribution%disc%ci(n,X_DIRECTION), &
                  distribution%info%ys + distribution%disc%ci(n,Y_DIRECTION): &
                  distribution%info%ye + distribution%disc%ci(n,Y_DIRECTION))
    end do
  end subroutine DistributionGatherValueToDirectionD2

  subroutine DistributionGatherValueToDirectionD3(distribution, val, out)
    type(distribution_type) distribution
    PetscScalar,intent(in),dimension(distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye, distribution%info%gzs:distribution%info%gze):: val
    PetscScalar,intent(out),dimension(0:distribution%b, &
         distribution%info%gxs:distribution%info%gxe, &
         distribution%info%gys:distribution%info%gye, &
         distribution%info%gzs:distribution%info%gze):: out
    
    PetscInt n
    do n=0,distribution%b
       out(n,distribution%info%xs:distribution%info%xe, &
            distribution%info%ys:distribution%info%ye, &
            distribution%info%zs:distribution%info%ze) = &
              val(distribution%info%xs + distribution%disc%ci(n,X_DIRECTION): &
                  distribution%info%xe + distribution%disc%ci(n,X_DIRECTION), &
                  distribution%info%ys + distribution%disc%ci(n,Y_DIRECTION): &
                  distribution%info%ye + distribution%disc%ci(n,Y_DIRECTION), &
                  distribution%info%zs + distribution%disc%ci(n,Z_DIRECTION): &
                  distribution%info%ze + distribution%disc%ci(n,Z_DIRECTION))
    end do
  end subroutine DistributionGatherValueToDirectionD3
end module LBM_Distribution_Function_module
    
