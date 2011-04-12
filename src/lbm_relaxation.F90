!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_relaxation.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            15:15:25 MDT
!!!     last modified:   12 April 2011
!!!       at:            12:16:28 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscbagdef.h"

module LBM_Relaxation_module
  use petsc
  use LBM_Relaxation_Bag_Data_type_module
  use LBM_Distribution_Function_type_module
  implicit none

  private
#include "lbm_definitions.h"

  type, public :: relaxation_type
     MPI_Comm comm
     PetscInt id
     PetscInt s
     PetscInt b
     PetscInt mode
     PetscScalar,pointer :: tau ! relaxation time
     PetscScalar,pointer,dimension(:) :: tau_mrt ! species of S vector for mrt

     ! dependents, set by discretization
     PetscScalar d_k
     PetscScalar c_s2

     ! bag 
     character(len=MAXWORDLENGTH):: name
     type(relaxation_bag_data_type),pointer:: data
     PetscBag bag
  end type relaxation_type
     
  interface PetscBagGetData
     subroutine PetscBagGetData(bag, data, ierr)
       use LBM_Relaxation_Bag_Data_type_module
       PetscBag bag
       type(relaxation_bag_data_type),pointer :: data
       PetscErrorCode ierr
     end subroutine PetscBagGetData
  end interface

  public :: RelaxationCreate, &
       RelaxationDestroy, &
       RelaxationSetSizes, &
       RelaxationSetName, &
       RelaxationSetID, &
       RelaxationSetMode, &
       RelaxationSetFromOptions, &
       RelaxationCollide

contains
  function RelaxationCreate(comm) result(relax)
    MPI_Comm comm
    type(relaxation_type),pointer:: relax
    allocate(relax)
    relax%comm = comm
    relax%s = -1
    relax%b = -1
    relax%d_k = 0.
    relax%c_s2 = 1.d0/3.d0

    nullify(relax%data)
    relax%bag = 0
    relax%name = ''
  end function RelaxationCreate

  subroutine RelaxationDestroy(relax, ierr)
    type(relaxation_type) relax
    PetscErrorCode ierr
    if (relax%bag /= 0) call PetscBagDestroy(relax%bag, ierr)
    if (associated(relax%tau_mrt)) deallocate(relax%tau_mrt)
  end subroutine RelaxationDestroy

  subroutine RelaxationSetSizes(relax, s, b)
    type(relaxation_type) relax
    PetscInt:: s, b
    relax%s = s
    relax%b = b
    allocate(relax%tau_mrt(0:b))
  end subroutine RelaxationSetSizes

  subroutine RelaxationSetMode(relax, mode)
    type(relaxation_type) relax
    PetscInt mode
    relax%mode = mode
  end subroutine RelaxationSetMode

  subroutine RelaxationSetName(relax, name)
    type(relaxation_type) relax
    character(len=MAXWORDLENGTH):: name
    
    relax%name = name
  end subroutine RelaxationSetName

  subroutine RelaxationSetID(relax, id)
    type(relaxation_type) relax
    PetscInt id
    relax%id = id
  end subroutine RelaxationSetID

  subroutine RelaxationSetFromOptions(relax, options, ierr)
    use LBM_Options_module
    type(relaxation_type) relax
    type(options_type) options
    PetscErrorCode ierr

    PetscSizeT sizeofscalar, sizeofdata
    PetscInt lcv
    character(len=MAXWORDLENGTH):: paramname, paramname2

    ! create the bag
    call PetscDataTypeGetSize(PETSC_SCALAR, sizeofscalar, ierr)
    sizeofdata = sizeofscalar 
    call PetscBagCreate(relax%comm, sizeofdata, relax%bag, ierr)
    call PetscBagGetData(relax%bag, relax%data, ierr)

    write(paramname, '(I1)') relax%id
    call PetscBagRegisterScalar(relax%bag, relax%data%tau, 1.d0, &
         trim(options%my_prefix)//'tau'//paramname, 'relaxation time', ierr)
    relax%tau => relax%data%tau

    call PetscBagSetName(relax%bag, TRIM(options%my_prefix)//relax%name, "", ierr)
  end subroutine RelaxationSetFromOptions

  subroutine RelaxationCollide(relax, fi, fi_eq, walls, dist)
    type(distribution_type) dist
    PetscScalar,dimension(0:dist%b, 1:dist%info%gxyzl):: fi
    PetscScalar,dimension(0:dist%b, 1:dist%info%gxyzl):: fi_eq
    PetscScalar,dimension(1:dist%info%gxyzl):: walls
    type(relaxation_type) relax
    PetscErrorCode ierr

    select case(dist%info%ndims)
    case(3)
       select case(relax%mode)
       case(RELAXATION_MODE_MRT)
          call RelaxationCollideMRTD3(relax, fi, fi_eq, walls, dist)
       case(RELAXATION_MODE_SRT)
          call RelaxationCollideSRTD3(relax, fi, fi_eq, walls, dist)
       case DEFAULT
          SETERRQ(1, 1, 'invalid relaxation mode in LBM', ierr)
       end select
    case(2)
       select case(relax%mode)
       case(RELAXATION_MODE_MRT)
          call RelaxationCollideMRTD2(relax, fi, fi_eq, walls, dist)
       case(RELAXATION_MODE_SRT)
          call RelaxationCollideSRTD2(relax, fi, fi_eq, walls, dist)
       case DEFAULT
          SETERRQ(1, 1, 'invalid relaxation mode in LBM', ierr)
       end select
    case DEFAULT
       SETERRQ(1, 1, 'invalid discretization in LBM', ierr)
    end select
  end subroutine RelaxationCollide

  subroutine RelaxationCollideSRTD3(relax, fi, fi_eq, walls, dist)
    ! input variables
    type(distribution_type) dist
    PetscScalar,dimension(0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi_eq
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls
    type(relaxation_type) relax

    ! local variables
    integer i,j,k,m                  ! loop variables

    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if(walls(i,j,k).eq.0) then
          fi(:,i,j,k) = fi(:,i,j,k) - (fi(:,i,j,k)-fi_eq(:,i,j,k))/relax%tau
       endif
    enddo
    enddo
    enddo
    return
  end subroutine RelaxationCollideSRTD3

  subroutine RelaxationCollideSRTD2(relax, fi, fi_eq, walls, dist)
    ! input variables
    type(distribution_type) dist
    PetscScalar,dimension(0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi_eq
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls
    type(relaxation_type) relax

    ! local variables
    integer i,j,k,m                  ! loop variables

    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if(walls(i,j).eq.0) then
          fi(:,i,j) = fi(:,i,j) - (fi(:,i,j)-fi_eq(:,i,j))/relax%tau
       endif
    enddo
    enddo
    return
  end subroutine RelaxationCollideSRTD2

  subroutine RelaxationCollideMRTD3(relax, fi, fi_eq, walls, dist)
    ! input variables
    type(distribution_type) dist
    PetscScalar,dimension(0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi_eq
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: walls
    type(relaxation_type) relax

    ! local variables
    integer i,j,k,m,n                  ! loop variables
    PetscScalar,dimension(0:dist%b):: mdiff
    PetscScalar,dimension(0:dist%b):: rhs
    
    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if(walls(i,j,k).eq.0) then
          do n=0,dist%b
             mdiff(n) = dot_product(dist%disc%mt_mrt(:,n), fi(:,i,j,k)-fi_eq(:,i,j,k))
          end do
          rhs(:)=0.d0
          do n=0,dist%b
             rhs(:) = rhs(:) - (relax%tau_mrt(n)*(mdiff(n)))/dist%disc%mmt_mrt(n)*dist%disc%mt_mrt(:,n)
          enddo
          fi(:,i,j,k)=fi(:,i,j,k) + rhs(:)
       endif
    enddo
    enddo
    enddo
  end subroutine RelaxationCollideMRTD3

  subroutine RelaxationCollideMRTD2(relax, fi, fi_eq, walls, dist)
    ! input variables
    type(distribution_type) dist
    PetscScalar,dimension(0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi_eq
    PetscScalar,dimension(dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: walls
    type(relaxation_type) relax

    ! local variables
    integer i,j,m,n                  ! loop variables
    PetscScalar,dimension(0:dist%b):: mdiff
    PetscScalar,dimension(0:dist%b):: rhs
    
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if(walls(i,j).eq.0) then
          do n=0,dist%b
             mdiff(n) = dot_product(dist%disc%mt_mrt(:,n), fi(:,i,j)-fi_eq(:,i,j))
          end do
          rhs(:)=0.d0
          do n=0,dist%b
             rhs(:) = rhs(:) - (relax%tau_mrt(n)*(mdiff(n)))/dist%disc%mmt_mrt(n)*dist%disc%mt_mrt(:,n)
          enddo
          fi(:,i,j)=fi(:,i,j) + rhs(:)
       endif
    enddo
    enddo
  end subroutine RelaxationCollideMRTD2
end module LBM_Relaxation_module

    
