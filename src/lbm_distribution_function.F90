!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_distribution_function.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            14:06:07 MDT
!!!     last modified:   29 March 2011
!!!       at:            08:12:49 PDT
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
     DM,pointer :: da
  end type distribution_type

  public:: DistributionCreate, &
       DistributionDestroy, &
       DistributionSetInfo, &
       DistributionSetDiscretization, &
       DistributionSetSizes, &
       DistributionSetDA, &
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
    nullify(distribution%da)
  end function DistributionCreate

  subroutine DistributionDestroy(distribution, ierr)
    type(distribution_type) distribution
    PetscErrorCode ierr
  end subroutine DistributionDestroy

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
    
  subroutine DistributionSetDA(distribution, da)
    type(distribution_type) distribution
    DM,target:: da
    distribution%da => da
  end subroutine DistributionSetDA

  subroutine DistributionCalcDensity(distribution, fi, walls, rho)
    type(distribution_type) distribution
    PetscScalar,dimension(1:distribution%s, 0:distribution%b, &
         1:distribution%info%gxyzl):: fi
    PetscScalar,dimension(1:distribution%s, 1:distribution%info%gxyzl):: rho
    PetscScalar,dimension(1:distribution%info%gxyzl):: walls

    select case(distribution%info%ndims)
    case(2)
       call DistributionCalcDensityD2(distribution, fi, walls, rho)
    case(3)
       call DistributionCalcDensityD3(distribution, fi, walls, rho)
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
       else:
          rho(:,i,j,k) = 0.d0
       end if
    end do
    end do
    end do
  end  subroutine DistributionCalcDensityD3

  subroutine DistributionCalcFlux(distribution, fi, walls, u)
    type(distribution_type) distribution
    PetscScalar,dimension(distribution%s, 0:distribution%b, &
         1:distribution%info%gxyzl):: fi
    PetscScalar,dimension(distribution%s, distribution%info%ndims, &
         1:distribution%info%gxyzl):: u
    PetscScalar,dimension(distribution%info%gxyzl):: walls

    select case(distribution%info%ndims)
    case(2)
       call DistributionCalcDensityD2(distribution, fi, walls, u)
    case(3)
       call DistributionCalcDensityD3(distribution, fi, walls, u)
    end select
  end  subroutine DistributionCalcDensity

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

end module LBM_Distribution_Function_module
    
