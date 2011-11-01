!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_relaxation.F90
!!!     version:         
!!!     created:         28 March 2011
!!!       on:            15:15:25 MDT
!!!     last modified:   18 October 2011
!!!       at:            13:48:12 MDT
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
     PetscScalar,pointer :: s1  ! MRT relaxation time
     PetscScalar,pointer :: s2  ! MRT relaxation time
     PetscScalar,pointer :: s3  ! MRT relaxation time
     PetscScalar,pointer :: s4  ! MRT relaxation time (only for 3D)
     PetscScalar,pointer :: s5  ! MRT relaxation time (only for 3D)
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
    relax%id = 0
    relax%s = -1
    relax%b = -1
    relax%mode = RELAXATION_MODE_SRT

    nullify(relax%tau)
    nullify(relax%s1)
    nullify(relax%s2)
    nullify(relax%s3)
    nullify(relax%s4)
    nullify(relax%s5)
    nullify(relax%tau_mrt)

    relax%d_k = 0. 
    relax%c_s2 = 1.d0/3.

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
    sizeofdata = 6*sizeofscalar 
    call PetscBagCreate(relax%comm, sizeofdata, relax%bag, ierr)
    call PetscBagGetData(relax%bag, relax%data, ierr)

    call PetscBagRegisterScalar(relax%bag, relax%data%tau, ONE_S, &
         trim(options%my_prefix)//'tau_'//trim(relax%name), 'relaxation time', ierr)
    relax%tau => relax%data%tau

    call PetscBagRegisterScalar(relax%bag, relax%data%s1, ONE_S, &
         trim(options%my_prefix)//'s1_'//trim(relax%name), 'MRT relaxation time', ierr)
    relax%s1 => relax%data%s1

    call PetscBagRegisterScalar(relax%bag, relax%data%s2, ONE_S, &
         trim(options%my_prefix)//'s2_'//trim(relax%name), 'MRT relaxation time', ierr)
    relax%s2 => relax%data%s2

    call PetscBagRegisterScalar(relax%bag, relax%data%s3, ONE_S, &
         trim(options%my_prefix)//'s3_'//trim(relax%name), 'MRT relaxation time', ierr)
    relax%s3 => relax%data%s3
   
    call PetscBagRegisterScalar(relax%bag, relax%data%s4, ONE_S, &
         trim(options%my_prefix)//'s4_'//trim(relax%name), 'MRT relaxation time', ierr)
    relax%s4 => relax%data%s4

    call PetscBagRegisterScalar(relax%bag, relax%data%s5, ONE_S, &
         trim(options%my_prefix)//'s5_'//trim(relax%name), 'MRT relaxation time', ierr)
    relax%s5 => relax%data%s5

    call PetscBagSetName(relax%bag, TRIM(options%my_prefix)//relax%name, "", ierr)
  end subroutine RelaxationSetFromOptions

  subroutine RelaxationCollide(relax, fi, fi_eq, walls, m, dist)
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, 0:dist%b, 1:dist%info%gxyzl):: fi
    PetscScalar,dimension(dist%s, 0:dist%b, 1:dist%info%gxyzl):: fi_eq
    PetscScalar,dimension(1:dist%info%rgxyzl):: walls
    PetscInt m
    type(relaxation_type) relax
    PetscErrorCode ierr

    select case(dist%info%ndims)
    case(3)
       select case(relax%mode)
       case(RELAXATION_MODE_MRT)
          call RelaxationCollideMRTD3(relax, fi, fi_eq, walls, m, dist)
       case(RELAXATION_MODE_SRT)
          call RelaxationCollideSRTD3(relax, fi, fi_eq, walls, m, dist)
       case DEFAULT
          SETERRQ(1, 1, 'invalid relaxation mode in LBM', ierr)
       end select
    case(2)
       select case(relax%mode)
       case(RELAXATION_MODE_MRT)
          call RelaxationCollideMRTD2(relax, fi, fi_eq, walls, m, dist)
       case(RELAXATION_MODE_SRT)
          call RelaxationCollideSRTD2(relax, fi, fi_eq, walls, m, dist)
       case DEFAULT
          SETERRQ(1, 1, 'invalid relaxation mode in LBM', ierr)
       end select
    case DEFAULT
       SETERRQ(1, 1, 'invalid discretization in LBM', ierr)
    end select
  end subroutine RelaxationCollide

  subroutine RelaxationCollideSRTD3(relax, fi, fi_eq, walls, m, dist)
    ! input variables
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi_eq
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: walls
    type(relaxation_type) relax
    PetscInt m

    ! local variables
    PetscInt i,j,k                  ! loop variables

    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if(walls(i,j,k).eq.0) then
          fi(m,:,i,j,k) = fi(m,:,i,j,k) - (fi(m,:,i,j,k)-fi_eq(m,:,i,j,k))/relax%tau
       endif
    enddo
    enddo
    enddo
    return
  end subroutine RelaxationCollideSRTD3

  subroutine RelaxationCollideSRTD2(relax, fi, fi_eq, walls, m, dist)
    ! input variables
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi_eq
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls
    type(relaxation_type) relax
    PetscInt m

    ! local variables
    PetscInt i,j,k                  ! loop variables

    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if(walls(i,j).eq.0) then
          fi(m,:,i,j) = fi(m,:,i,j) - (fi(m,:,i,j)-fi_eq(m,:,i,j))/relax%tau
       endif
    enddo
    enddo
    return
  end subroutine RelaxationCollideSRTD2

  subroutine RelaxationCollideMRTD3(relax, fi, fi_eq, walls, m, dist)
    ! input variables
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye, dist%info%gzs:dist%info%gze):: fi_eq
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye, dist%info%rgzs:dist%info%rgze):: walls
    type(relaxation_type) relax
    PetscInt m

    ! local variables
    PetscInt i,j,k,n                  ! loop variables
    PetscScalar :: mdiff
    PetscScalar :: d_fi(0:dist%b)
    
    do k=dist%info%zs,dist%info%ze
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if(walls(i,j,k).eq.0) then
          d_fi = fi(m,:,i,j,k)-fi_eq(m,:,i,j,k)
          do n=0,dist%b
             mdiff = sum(dist%disc%mt_mrt(:,n)*d_fi,1)
             fi(m,:,i,j,k) = fi(m,:,i,j,k) - &
                (relax%tau_mrt(n)*(mdiff))/dist%disc%mmt_mrt(n)*dist%disc%mt_mrt(:,n)
          enddo
       endif
    enddo
    enddo
    enddo
  end subroutine RelaxationCollideMRTD3

  subroutine RelaxationCollideMRTD2(relax, fi, fi_eq, walls, m, dist)
    ! input variables
    type(distribution_type) dist
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi
    PetscScalar,dimension(dist%s, 0:dist%b, dist%info%gxs:dist%info%gxe, &
         dist%info%gys:dist%info%gye):: fi_eq
    PetscScalar,dimension(dist%info%rgxs:dist%info%rgxe, &
         dist%info%rgys:dist%info%rgye):: walls
    type(relaxation_type) relax
    PetscInt m

    ! local variables
    PetscInt i,j,n                  ! loop variables
    PetscScalar :: mdiff
    PetscScalar :: d_fi(0:dist%b)
   
    do j=dist%info%ys,dist%info%ye
    do i=dist%info%xs,dist%info%xe
       if(walls(i,j).eq.0) then
          d_fi = fi(m,:,i,j)-fi_eq(m,:,i,j)
          do n=0,dist%b
             mdiff = sum(dist%disc%mt_mrt(:,n)*d_fi,1)
             fi(m,:,i,j) = fi(m,:,i,j) - &
                  relax%tau_mrt(n)*mdiff/dist%disc%mmt_mrt(n)*dist%disc%mt_mrt(:,n)
          enddo
       endif
    enddo
    enddo
  end subroutine RelaxationCollideMRTD2
end module LBM_Relaxation_module

    
