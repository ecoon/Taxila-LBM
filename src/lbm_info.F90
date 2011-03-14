!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        grid_info.F90
!!!     version:         
!!!     created:         06 December 2010
!!!       on:            15:19:22 MST
!!!     last modified:   14 March 2011
!!!       at:            17:55:42 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ ldeo.columbia.edu
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

  module LBM_Info_module
    use petsc
    use LBM_Discretization_module
    implicit none

    private
#include "lbm_definitions.h"

    type, public:: info_type
       PetscInt xs,xe,xl,gxs,gxe,gxl
       PetscInt ys,ye,yl,gys,gye,gyl
       PetscInt zs,ze,zl,gzs,gze,gzl
       PetscInt xyzl, gxyzl
       PetscInt NX,NY,NZ
       PetscInt nproc_x, nproc_y, nproc_z
       PetscInt id, nproc
       PetscInt s
       PetscInt ndims
       PetscBool,pointer:: periodic(:)
       PetscReal,pointer:: gridsize(:)
       PetscReal,pointer:: corners(:,:)
       character(len=MAXWORDLENGTH):: options_prefix
       type(discretization_type),pointer :: flow_disc
       type(discretization_type),pointer :: tran_disc

    end type info_type
    
    public :: InfoCreate, &
         InfoDestroy, &
         InfoSetFromOptions, &
         InfoView, &
         InfoGatherValueToDirection

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

      info%xyzl = -1
      info%gxyzl = -1

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
      info%ndims = -1
      nullify(info%flow_disc)
      nullify(info%tran_disc)

      info%id = -1
      info%nproc = -1

      info%nproc_x = -1
      info%nproc_y = -1
      info%nproc_z = -1

      info%options_prefix = ''

      nullify(info%periodic)
      nullify(info%gridsize)
      nullify(info%corners)
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
      PetscInt tmp_dim
      
      info%options_prefix = options%my_prefix

      ! grab sizes
      info%s = 1
      call PetscOptionsGetInt(options%my_prefix,'-s', info%s,flag,ierr)

      ! update sizes
      call PetscOptionsGetInt(options%my_prefix,'-nx',info%NX,flag,ierr)
      call PetscOptionsGetInt(options%my_prefix,'-ny',info%NY,flag,ierr)
      call PetscOptionsGetInt(options%my_prefix,'-nz',info%NZ,flag,ierr)

      ! set the flow discretization
      call PetscOptionsGetString(options%my_prefix, '-discretization', &
           discretization, flag, ierr)
      if (.not.flag) then
         call PetscOptionsGetString(options%my_prefix, '-flow_discretization', &
              discretization, flag, ierr)
      end if
      if (.not.flag) discretization = 'D3Q19'

      test_discretization = 'd3q19'
      if (StringCompareIgnoreCase(discretization, test_discretization, 6)) then
         call InfoSetDiscretizationD3Q19(info, info%flow_disc)
      endif
      
      test_discretization = 'd2q9'
      if (StringCompareIgnoreCase(discretization, test_discretization, 5)) then
         call InfoSetDiscretizationD2Q9(info, info%flow_disc)
      end if
      
      ! make sure we have an answer
      if (.not.associated(info%flow_disc)) SETERRQ(1, 1, 'Invalid Discretization', ierr)
      tmp_dim = info%ndims
      
      ! set the tran discretization
      call PetscOptionsGetString(options%my_prefix, '-transport_discretization', &
           discretization, flag, ierr)

      if (flag) then
         test_discretization = 'd3q19'
         if (StringCompareIgnoreCase(discretization, test_discretization, 6)) then
            call InfoSetDiscretizationD3Q19(info, info%tran_disc)
         endif

         test_discretization = 'd2q9'
         if (StringCompareIgnoreCase(discretization, test_discretization, 5)) then
            call InfoSetDiscretizationD2Q9(info, info%tran_disc)
         end if
         
         ! make sure we have an answer
         if (.not.associated(info%tran_disc)) SETERRQ(1, 1, 'Invalid Discretization', ierr)
         if (tmp_dim /= info%ndims) SETERRQ(1,1,"Discretization dimensions don't match", ierr)
      end if

      ! set up the grid
      allocate(info%corners(info%ndims,2))
      allocate(info%gridsize(info%ndims))
      allocate(info%periodic(info%ndims))
      info%periodic = PETSC_FALSE
      
      ! nullify z-things for a 2D DA
      if (info%ndims.eq.2) then
         info%zs = PETSC_NULL_INTEGER
         info%ze = PETSC_NULL_INTEGER
         info%zl = 1
         info%NZ = 1
         info%gzs = PETSC_NULL_INTEGER
         info%gze = PETSC_NULL_INTEGER
         info%gzl = 1
         info%nproc_z = PETSC_NULL_INTEGER
      end if

      call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_x', &
           info%periodic(X_DIRECTION), flag,ierr)
      call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_y', &
           info%periodic(Y_DIRECTION), flag,ierr)
      if (info%ndims > 2) then
         call PetscOptionsGetBool(options%my_prefix, '-bc_periodic_z', &
              info%periodic(Z_DIRECTION), flag,ierr)
      end if
      
    end subroutine InfoSetFromOptions

    subroutine InfoSetDiscretizationD3Q19(info, disc)
      use LBM_Discretization_D3Q19_module
      type(info_type) info
      type(discretization_type),pointer:: disc

      allocate(disc)
      call DiscretizationSetup_D3Q19(disc)
      info%ndims = disc%ndims
    end subroutine InfoSetDiscretizationD3Q19

    subroutine InfoSetDiscretizationD2Q9(info, disc)
      use LBM_Discretization_D2Q9_module
      type(info_type) info
      type(discretization_type),pointer:: disc

      allocate(disc)
      call DiscretizationSetup_D2Q9(disc)
      info%ndims = disc%ndims
    end subroutine InfoSetDiscretizationD2Q9

    subroutine InfoView(info)
      type(info_type) info

      print*, ' Domain size =', info%NX,info%NY,info%NZ
    end subroutine InfoView

    subroutine InfoDestroy(info, ierr)
      type(info_type) info
      PetscErrorCode ierr

      if (associated(info%flow_disc)) call DiscretizationDestroy(info%flow_disc, ierr)
      if (associated(info%tran_disc)) call DiscretizationDestroy(info%tran_disc, ierr)
      if (associated(info%corners)) deallocate(info%corners)
      if (associated(info%gridsize)) deallocate(info%gridsize)
      if (associated(info%periodic)) deallocate(info%periodic)
    end subroutine InfoDestroy

    subroutine InfoGatherValueToDirection(info, val, out)
      type(info_type) info
      PetscScalar,intent(in),dimension(1:info%gxyzl):: val
      PetscScalar,intent(out),dimension(0:info%flow_disc%b, 1:info%gxyzl):: out
      PetscErrorCode ierr

      if (info%ndims.eq.2) then
         call InfoGatherValueToDirection2D(info, val, out)
      else if (info%ndims.eq.3) then
         call InfoGatherValueToDirection3D(info, val, out)
      else 
         SETERRQ(1, 1, 'invalid ndims in LBM', ierr)
      end if
    end subroutine InfoGatherValueToDirection


    subroutine InfoGatherValueToDirection2D(info, val, out)
      type(info_type) info
      PetscScalar,intent(in),dimension(info%gxs:info%gxe, &
           info%gys:info%gye):: val
      PetscScalar,intent(out),dimension(0:info%flow_disc%b, &
           info%gxs:info%gxe, info%gys:info%gye):: out

      PetscInt n
      do n=0,info%flow_disc%b
         out(n,info%xs:info%xe,info%ys:info%ye) = val( &
              info%xs+info%flow_disc%ci(n,X_DIRECTION): &
              info%xe+info%flow_disc%ci(n,X_DIRECTION), &
              info%ys+info%flow_disc%ci(n,Y_DIRECTION): &
              info%ye+info%flow_disc%ci(n,Y_DIRECTION))
      end do
    end subroutine InfoGatherValueToDirection2D

    subroutine InfoGatherValueToDirection3D(info, val, out)
      type(info_type) info
      PetscScalar,intent(in),dimension(info%gxs:info%gxe, &
           info%gys:info%gye, info%gzs:info%gze):: val
      PetscScalar,intent(out),dimension(0:info%flow_disc%b, &
           info%gxs:info%gxe, info%gys:info%gye, info%gzs:info%gze):: out

      PetscInt n
      do n=0,info%flow_disc%b
         out(n,info%xs:info%xe,info%ys:info%ye,info%zs:info%ze) = val( &
              info%xs+info%flow_disc%ci(n,X_DIRECTION): &
              info%xe+info%flow_disc%ci(n,X_DIRECTION), &
              info%ys+info%flow_disc%ci(n,Y_DIRECTION): &
              info%ye+info%flow_disc%ci(n,Y_DIRECTION), &
              info%zs+info%flow_disc%ci(n,Z_DIRECTION): &
              info%ze+info%flow_disc%ci(n,Z_DIRECTION))
      end do
    end subroutine InfoGatherValueToDirection3D
  end module LBM_Info_module
