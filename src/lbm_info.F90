!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        grid_info.F90
!!!     version:         
!!!     created:         06 December 2010
!!!       on:            15:19:22 MST
!!!     last modified:   01 February 2011
!!!       at:            18:06:37 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ ldeo.columbia.edu
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

  module LBM_Info_module
    use petsc
    implicit none

    private
#include "lbm_definitions.h"

    type, public:: info_type
       PetscInt xs,xe,xl,gxs,gxe,gxl
       PetscInt ys,ye,yl,gys,gye,gyl
       PetscInt zs,ze,zl,gzs,gze,gzl
       PetscInt NX,NY,NZ
       PetscInt nproc_x, nproc_y, nproc_z
       PetscInt id, nproc
       PetscInt s
       PetscInt b
       PetscInt dim
       PetscInt discretization
       character(len=MAXWORDLENGTH):: options_prefix
    end type info_type
    
    public :: InfoCreate, &
         InfoSetFromOptions, &
         InfoView, &
         InfoDestroy

  contains
    function InfoCreate() result(info)
      type(info_type),pointer:: info
      allocate(info)

      info%xs = -1
      info%xe = -1
      info%xl = -1
      info%ys = -1
      info%ye = -1
      info%yl = -1
      info%zs = -1
      info%ze = -1
      info%zl = -1

      info%gxs = -1
      info%gxe = -1
      info%gxl = -1
      info%gys = -1
      info%gye = -1
      info%gyl = -1
      info%gzs = -1
      info%gze = -1
      info%gzl = -1

      info%NX = -1
      info%NY = -1
      info%NZ = -1

      info%s = -1
      info%b = -1
      info%dim = -1
      info%discretization = NULL_DISCRETIZATION

      info%id = -1
      info%nproc = -1

      info%options_prefix = ''

    end function InfoCreate

    subroutine InfoSetFromOptions(info, options, ierr)
      use LBM_Options_module
      use String_module

      type(info_type) info
      type(options_type) options
      PetscErrorCode ierr

      PetscBool flag
      character(len=MAXWORDLENGTH):: discretization
      character(len=MAXWORDLENGTH):: test_discretization
      
      info%options_prefix = options%my_prefix

      ! grab sizes
      info%s = 1
      call PetscOptionsGetInt(options%my_prefix,'-s', info%s,flag,ierr)

      ! update sizes
      call PetscOptionsGetInt(options%my_prefix,'-nx',info%NX,flag,ierr)
      call PetscOptionsGetInt(options%my_prefix,'-ny',info%NY,flag,ierr)
      call PetscOptionsGetInt(options%my_prefix,'-nz',info%NZ,flag,ierr)

      call PetscOptionsGetString(options%my_prefix,'-discretization',discretization,flag,ierr)
      if (.not.flag) discretization = 'D3Q19'

      test_discretization = 'd3q19'
      if (StringCompareIgnoreCase(discretization, test_discretization, 6)) then
         call InfoSetDiscretizationD3Q19(info)
      endif
      
      test_discretization = 'd2q9'
      if (StringCompareIgnoreCase(discretization, test_discretization, 5)) then
         call InfoSetDiscretizationD2Q9(info)
      end if

      ! make sure we have an answer
      if (info%discretization .eq. NULL_DISCRETIZATION) then
         write(*,*) 'INVALID DISCRETIZATION:', discretization
         ierr = 1
         CHKERRQ(ierr)
      end if

    end subroutine InfoSetFromOptions

    subroutine InfoSetDiscretizationD3Q19(info)
      use LBM_Discretization_D3Q19_module
      type(info_type) info
      
      info%discretization = D3Q19_DISCRETIZATION
      info%dim = discretization_dims
      info%b = discretization_directions
    end subroutine InfoSetDiscretizationD3Q19

    subroutine InfoSetDiscretizationD2Q9(info)
      use LBM_Discretization_D2Q9_module
      use LBM_Discretization_D2Q9_module
      type(info_type) info
      
      info%discretization = D2Q9_DISCRETIZATION
      info%dim = 2
      info%b = 8
    end subroutine InfoSetDiscretizationD2Q9

    subroutine InfoView(info)
      type(info_type) info

      print*, ' Domain size =', info%NX,info%NY,info%NZ
    end subroutine InfoView

    subroutine InfoDestroy(info, ierr)
      type(info_type) info
      PetscErrorCode ierr
    end subroutine InfoDestroy

  end module LBM_Info_module
