!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        grid_info.F90
!!!     version:         
!!!     created:         06 December 2010
!!!       on:            15:19:22 MST
!!!     last modified:   17 March 2011
!!!       at:            17:55:22 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ ldeo.columbia.edu
!!!  
!!! ====================================================================

#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscbagdef.h"

module LBM_Info_module
  use petsc
  use LBM_Info_Bag_Data_type_module
  implicit none

  private
#include "lbm_definitions.h"
  
  type, public:: info_type
     MPI_Comm comm
     PetscInt ndims

     ! bagged params
     PetscInt,pointer:: NX,NY,NZ
     PetscBool,pointer,dimension(:) :: periodic
     PetscScalar,pointer,dimension(:,:) :: corners

     ! dependent parameters
     PetscInt xs,xe,xl,gxs,gxe,gxl
     PetscInt ys,ye,yl,gys,gye,gyl
     PetscInt zs,ze,zl,gzs,gze,gzl
     PetscInt xyzl, gxyzl
     PetscInt nproc_x, nproc_y, nproc_z, nprocs
     PetscInt rank
     PetscReal,pointer,dimension(:) :: gridsize

     ! bag
     character(len=MAXWORDLENGTH):: name
     type(phase_bag_data_type),pointer:: data
     PetscBag bag
  end type info_type

  interface PetscBagGetData
     subroutine PetscBagGetData(bag, data, ierr)
       PetscBag bag
       type(phase_bag_data_type),pointer :: data
       PetscErrorCode ierr
     end subroutine PetscBagGetData
  end interface

  public :: InfoCreate, &
       InfoDestroy, &
       InfoSetFromOptions, &
       InfoView, &
       InfoGatherValueToDirection

contains
  function InfoCreate(comm) result(info)
    MPI_Comm comm
    type(info_type),pointer:: info
    allocate(info)
    info%comm = comm
    info%name = 'grid'

    ! intialize
    nullify(info%NX)
    nullify(info%NY)
    nullify(info%NZ)
    info%ndims = -1

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

    info%id = -1
    info%nproc = -1
    info%nproc_x = -1
    info%nproc_y = -1
    info%nproc_z = -1

    nullify(info%periodic)
    nullify(info%gridsize)
    nullify(info%corners)
    nullify(info%data)
    info%bag = 0
  end function InfoCreate

  subroutine InfoSetFromOptions(info, options, ierr)
    use LBM_Options_module

    type(info_type) info
    type(options_type) options
    PetscErrorCode ierr
      
    PetscInt sizeofint, sizeofscalar, sizeofbool, sizeofdata

    ! grab, allocate sizes
    info%ndims = options%ndims
    allocate(info%data%periodic(1:info%ndims))
    allocate(info%data%corners(1:info%ndims,2))
    allocate(info%gridsize(1:info%ndims))

    ! create the bag
    call PetscDataTypeGetSize(PETSC_SCALAR, sizeofscalar, ierr)
    call PetscDataTypeGetSize(PETSC_BOOL, sizeofbool, ierr)
    call PetscDataTypeGetSize(PETSC_INT, sizeofint, ierr)
    sizeofdata = info%ndims*2*sizeofscalar + info%ndims*sizeofbool + 3*sizeofint
    call PetscBagCreate(info%comm, sizeofdata, info%bag, ierr)
    call PetscBagSetName(info%bag, TRIM(options%my_prefix)//info%name, ierr)

    ! register data
    ! -- grid size
    call PetscBagGetData(info%bag, info%data, ierr)
    call PetscBagRegisterInt(info%bag, info%data%NX, 0, 'NX', 'grid size in X', ierr)
    call PetscBagRegisterInt(info%bag, info%data%NY, 0, 'NY', 'grid size in Y', ierr)
    if (info%ndims > 2) then
       call PetscBagRegisterInt(info%bag, info%data%NZ, 0, 'NZ', 'grid size in Z', ierr)
    else 
       info%data%NZ = 1
    end if
    info%NX => info%data%NX
    info%NY => info%data%NY
    info%NZ => info%data%NZ

    ! -- grid perioidicty
    call PetscBagRegisterBool(info%bag, info%data%periodic(X_DIRECTION), PETSC_FALSE, &
         'bc_x_periodic', 'x-direction periodic?', ierr)
    call PetscBagRegisterBool(info%bag, info%data%periodic(Y_DIRECTION), PETSC_FALSE, &
         'bc_y_periodic', 'y-direction periodic?', ierr)
    if (info%ndims > 2) then
       call PetscBagRegisterBool(info%bag, info%data%periodic(Z_DIRECTION), PETSC_FALSE, &
            'bc_z_periodic', 'z-direction periodic?', ierr)
    end if
    info%periodic => info%data%periodic

    ! -- grid corners
    call PetscBagRegisterScalar(info%bag, info%data%corners(X_DIRECTION,1), 0.d0, &
         'x_start', 'lower x coordinate', ierr)
    call PetscBagRegisterScalar(info%bag, info%data%corners(X_DIRECTION,2), 1.d0, &
         'x_end', 'upper x coordinate', ierr)
    call PetscBagRegisterScalar(info%bag, info%data%corners(Y_DIRECTION,1), 0.d0, &
         'y_start', 'lower y coordinate', ierr)
    call PetscBagRegisterScalar(info%bag, info%data%corners(Y_DIRECTION,2), 1.d0, &
         'y_end', 'upper y coordinate', ierr)
    if (info%ndims > 2) then
       call PetscBagRegisterScalar(info%bag, info%data%corners(Z_DIRECTION,1), 0.d0, &
            'z_start', 'lower z coordinate', ierr)
       call PetscBagRegisterScalar(info%bag, info%data%corners(Z_DIRECTION,2), 1.d0, &
            'z_end', 'upper z coordinate', ierr)
    end if
    info%corners => info%data%corners
      
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

    ! calculate grid spacing
    if (info%periodic(X_DIRECTION)) then
       info%gridsize(X_DIRECTION) = (info%corners(X_DIRECTION,2) - &
            info%corners(X_DIRECTION, 1))/dble(info%NX)
    else 
       info%gridsize(X_DIRECTION) = (info%corners(X_DIRECTION,2) - &
            info%corners(X_DIRECTION, 1))/dble(info%NX-1)
    end if
    if (info%periodic(Y_DIRECTION)) then
       info%gridsize(Y_DIRECTION) = (info%corners(Y_DIRECTION,2) - &
            info%corners(Y_DIRECTION, 1))/dble(info%NY)
    else 
       info%gridsize(Y_DIRECTION) = (info%corners(Y_DIRECTION,2) - &
            info%corners(Y_DIRECTION, 1))/dble(info%NY-1)
    end if
    if (info%ndims > 2) then
       if (info%periodic(Z_DIRECTION)) then
          info%gridsize(Z_DIRECTION) = (info%corners(Z_DIRECTION,2) - &
               info%corners(Z_DIRECTION, 1))/dble(info%NZ)
       else 
          info%gridsize(Z_DIRECTION) = (info%corners(Z_DIRECTION,2) - &
               info%corners(Z_DIRECTION, 1))/dble(info%NZ-1)
       end if
    else
       info%gridsize(Z_DIRECTION)= 1.d0
    end if
  end subroutine InfoSetFromOptions

  subroutine InfoView(info)
    type(info_type) info
    
    print*, ' Domain size =', info%NX,info%NY,info%NZ
  end subroutine InfoView
  
  subroutine InfoDestroy(info, ierr)
    type(info_type) info
    PetscErrorCode ierr
    
    if (associated(info%gridsize)) deallocate(info%gridsize)
    if (info%bag /= 0) call PetscBagDestroy(info%bag, ierr)
  end subroutine InfoDestroy
  
  subroutine InfoGatherValueToDirection(info, val, out, disc)
    use LBM_Discretization_type_module
    type(info_type) info
    type(disc_type) disc
    PetscScalar,intent(in),dimension(1:info%gxyzl):: val
    PetscScalar,intent(out),dimension(0:disc%b, 1:info%gxyzl):: out
    PetscErrorCode ierr
    
    if (info%ndims.eq.2) then
       call InfoGatherValueToDirection2D(info, val, out, disc)
    else if (info%ndims.eq.3) then
       call InfoGatherValueToDirection3D(info, val, out, disc)
    else 
       SETERRQ(1, 1, 'invalid ndims in LBM', ierr)
    end if
  end subroutine InfoGatherValueToDirection

  
  subroutine InfoGatherValueToDirection2D(info, val, out, disc)
    type(info_type) info
    type(disc_type) disc
    PetscScalar,intent(in),dimension(info%gxs:info%gxe, &
         info%gys:info%gye):: val
    PetscScalar,intent(out),dimension(0:disc%b, &
         info%gxs:info%gxe, info%gys:info%gye):: out

    PetscInt n
    do n=0,disc%b
       out(n,info%xs:info%xe,info%ys:info%ye) = val( &
            info%xs+info%flow_disc%ci(n,X_DIRECTION): &
            info%xe+info%flow_disc%ci(n,X_DIRECTION), &
            info%ys+info%flow_disc%ci(n,Y_DIRECTION): &
            info%ye+info%flow_disc%ci(n,Y_DIRECTION))
    end do
  end subroutine InfoGatherValueToDirection2D

  subroutine InfoGatherValueToDirection3D(info, val, out, disc)
    type(info_type) info
    type(disc_type) disc
    PetscScalar,intent(in),dimension(info%gxs:info%gxe, &
         info%gys:info%gye, info%gzs:info%gze):: val
    PetscScalar,intent(out),dimension(0:disc%b, &
         info%gxs:info%gxe, info%gys:info%gye, info%gzs:info%gze):: out
    
    PetscInt n
    do n=0,disc%b
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
