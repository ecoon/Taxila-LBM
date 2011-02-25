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
  implicit none

  private 

  public:: LBMAddFluidFluidForces, &
       LBMAddBodyForces, &
       LBMAddFluidSolidForces, &
       LBMZeroBoundaryForces

contains
  ! --- Fluid-fluid interaction forces, from
  ! ---  (Kang 2002 Eq. 6)
  subroutine LBMAddFluidFluidForces(rho,Fx,Fy,Fz,walls,info,constants)
    !     NONLOCAL IN RHO
    use LBM_Discretization_D3Q19_module

    ! input
    type(info_type) info
    type(constants_type) constants
    PetscScalar,dimension(1:info%s,info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: rho
    PetscScalar,dimension(1:info%s, info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: Fx,Fy,Fz
    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: walls

    ! local
    integer i,j,k,m,n
    PetscScalar,dimension(1:info%s,0:info%b,info%gxs:info%gxe, &
         info%gys:info%gye,info%gzs:info%gze)::tmp
!!$    PetscScalar,dimension(info%gxs:info%gxe,info%gys:info%gye, &
!!$         info%gzs:info%gze):: fluidfluid

    do m=1,info%s
       tmp(m, 1,:,:,:)= (cshift(rho(m,:,:,:),shift= 1,dim=1))
       tmp(m, 2,:,:,:)= (cshift(rho(m,:,:,:),shift= 1,dim=2))
       tmp(m, 3,:,:,:)= (cshift(rho(m,:,:,:),shift=-1,dim=1))
       tmp(m, 4,:,:,:)= (cshift(rho(m,:,:,:),shift=-1,dim=2))
       tmp(m, 5,:,:,:)= (cshift(rho(m,:,:,:),shift= 1,dim=3))
       tmp(m, 6,:,:,:)= (cshift(rho(m,:,:,:),shift=-1,dim=3))

       tmp(m,7 ,:,:,:)= (cshift(tmp(m,1,:,:,:),shift= 1,dim=2))
       tmp(m,8 ,:,:,:)= (cshift(tmp(m,3,:,:,:),shift= 1,dim=2))
       tmp(m,9 ,:,:,:)= (cshift(tmp(m,3,:,:,:),shift=-1,dim=2))
       tmp(m,10,:,:,:)= (cshift(tmp(m,1,:,:,:),shift=-1,dim=2))

       tmp(m,11,:,:,:)= (cshift(tmp(m,5,:,:,:),shift= 1,dim=1))
       tmp(m,12,:,:,:)= (cshift(tmp(m,5,:,:,:),shift=-1,dim=1))
       tmp(m,13,:,:,:)= (cshift(tmp(m,6,:,:,:),shift=-1,dim=1))
       tmp(m,14,:,:,:)= (cshift(tmp(m,6,:,:,:),shift= 1,dim=1))

       tmp(m,15,:,:,:)= (cshift(tmp(m,5,:,:,:),shift= 1,dim=2))
       tmp(m,16,:,:,:)= (cshift(tmp(m,5,:,:,:),shift=-1,dim=2))
       tmp(m,17,:,:,:)= (cshift(tmp(m,6,:,:,:),shift=-1,dim=2))
       tmp(m,18,:,:,:)= (cshift(tmp(m,6,:,:,:),shift= 1,dim=2))
    enddo

!!$  MLP: I don't think we need the fluidfluid and I think we can vectorize the loops
    !fluidfluid = (cshift(walls, shift= 1, dim=1) + cshift(walls, shift=-1, dim=1) + walls)
    do i=info%xs,info%xe 
       do j=info%ys,info%ye 
          do k=info%zs,info%ze
             if (walls(i,j,k).eq.0) then
                do n=1,6
                   Fx(1,i,j,k)=Fx(1,i,j,k)-rho(1,i,j,k)*tmp(2,n,i,j,k)*cix(n)*constants%g
                   Fx(1,i,j,k)=Fx(1,i,j,k)-rho(1,i,j,k)*tmp(1,n,i,j,k)*cix(n)*constants%g11
                   Fx(2,i,j,k)=Fx(2,i,j,k)-rho(2,i,j,k)*tmp(1,n,i,j,k)*cix(n)*constants%g
                   Fx(2,i,j,k)=Fx(2,i,j,k)-rho(2,i,j,k)*tmp(2,n,i,j,k)*cix(n)*constants%g22
                   Fy(1,i,j,k)=Fy(1,i,j,k)-rho(1,i,j,k)*tmp(2,n,i,j,k)*ciy(n)*constants%g
                   Fy(1,i,j,k)=Fy(1,i,j,k)-rho(1,i,j,k)*tmp(1,n,i,j,k)*ciy(n)*constants%g11
                   Fy(2,i,j,k)=Fy(2,i,j,k)-rho(2,i,j,k)*tmp(1,n,i,j,k)*ciy(n)*constants%g
                   Fy(2,i,j,k)=Fy(2,i,j,k)-rho(2,i,j,k)*tmp(2,n,i,j,k)*ciy(n)*constants%g22
                   Fz(1,i,j,k)=Fz(1,i,j,k)-rho(1,i,j,k)*tmp(2,n,i,j,k)*ciz(n)*constants%g
                   Fz(1,i,j,k)=Fz(1,i,j,k)-rho(1,i,j,k)*tmp(1,n,i,j,k)*ciz(n)*constants%g11
                   Fz(2,i,j,k)=Fz(2,i,j,k)-rho(2,i,j,k)*tmp(1,n,i,j,k)*ciz(n)*constants%g
                   Fz(2,i,j,k)=Fz(2,i,j,k)-rho(2,i,j,k)*tmp(2,n,i,j,k)*ciz(n)*constants%g22
                enddo
                do n=7,18
                   Fx(1,i,j,k)=Fx(1,i,j,k)-rho(1,i,j,k)*tmp(2,n,i,j,k)*cix(n)*constants%g*0.5
                   Fx(1,i,j,k)=Fx(1,i,j,k)-rho(1,i,j,k)*tmp(1,n,i,j,k)*cix(n)*constants%g11*0.5
                   Fx(2,i,j,k)=Fx(2,i,j,k)-rho(2,i,j,k)*tmp(1,n,i,j,k)*cix(n)*constants%g*0.5
                   Fx(2,i,j,k)=Fx(2,i,j,k)-rho(2,i,j,k)*tmp(2,n,i,j,k)*cix(n)*constants%g22*0.5
                   Fy(1,i,j,k)=Fy(1,i,j,k)-rho(1,i,j,k)*tmp(2,n,i,j,k)*ciy(n)*constants%g*0.5
                   Fy(1,i,j,k)=Fy(1,i,j,k)-rho(1,i,j,k)*tmp(1,n,i,j,k)*ciy(n)*constants%g11*0.5
                   Fy(2,i,j,k)=Fy(2,i,j,k)-rho(2,i,j,k)*tmp(1,n,i,j,k)*ciy(n)*constants%g*0.5
                   Fy(2,i,j,k)=Fy(2,i,j,k)-rho(2,i,j,k)*tmp(2,n,i,j,k)*ciy(n)*constants%g22*0.5
                   Fz(1,i,j,k)=Fz(1,i,j,k)-rho(1,i,j,k)*tmp(2,n,i,j,k)*ciz(n)*constants%g*0.5
                   Fz(1,i,j,k)=Fz(1,i,j,k)-rho(1,i,j,k)*tmp(1,n,i,j,k)*ciz(n)*constants%g11*0.5
                   Fz(2,i,j,k)=Fz(2,i,j,k)-rho(2,i,j,k)*tmp(1,n,i,j,k)*ciz(n)*constants%g*0.5
                   Fz(2,i,j,k)=Fz(2,i,j,k)-rho(2,i,j,k)*tmp(2,n,i,j,k)*ciz(n)*constants%g22*0.5
                enddo
             end if
          end do
       end do
    end do

!!$    !fluidfluid = cshift(walls, shift= 1, dim=2) + cshift(walls, shift=-1, dim=2) + walls
!!$    do i=info%xs,info%xe 
!!$       do j=info%ys,info%ye 
!!$          do k=info%zs,info%ze
!!$             if (walls(i,j,k).eq.0) then
!!$                ! where i'm fluid and fluid on both sides of me
!!$                do n=1,6
!!$                   Fy(1,i,j,k)=Fy(1,i,j,k)+tmp(2,n,i,j,k)*ciy(n)*constants%g
!!$                   Fy(1,i,j,k)=Fy(1,i,j,k)+tmp(1,n,i,j,k)*ciy(n)*constants%g11
!!$                   Fy(2,i,j,k)=Fy(2,i,j,k)+tmp(1,n,i,j,k)*ciy(n)*constants%g
!!$                   Fy(2,i,j,k)=Fy(2,i,j,k)+tmp(2,n,i,j,k)*ciy(n)*constants%g22
!!$                enddo
!!$                do n=7,18
!!$                   Fy(1,i,j,k)=Fy(1,i,j,k)+tmp(2,n,i,j,k)*ciy(n)*constants%g*0.5
!!$                   Fy(1,i,j,k)=Fy(1,i,j,k)+tmp(1,n,i,j,k)*ciy(n)*constants%g11*0.5
!!$                   Fy(2,i,j,k)=Fy(2,i,j,k)+tmp(1,n,i,j,k)*ciy(n)*constants%g*0.5
!!$                   Fy(2,i,j,k)=Fy(2,i,j,k)+tmp(2,n,i,j,k)*ciy(n)*constants%g22*0.5
!!$                enddo
!!$             end if
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    !fluidfluid = cshift(walls, shift= 1, dim=3) + cshift(walls, shift=-1, dim=3) + walls
!!$    do i=info%xs,info%xe 
!!$       do j=info%ys,info%ye 
!!$          do k=info%zs,info%ze
!!$             if (walls(i,j,k).eq.0) then
!!$                ! where i'm fluid and fluid on both sides of me
!!$                do n=1,6
!!$                   Fz(1,i,j,k)=Fz(1,i,j,k)+tmp(2,n,i,j,k)*ciz(n)*constants%g
!!$                   Fz(1,i,j,k)=Fz(1,i,j,k)+tmp(1,n,i,j,k)*ciz(n)*constants%g11
!!$                   Fz(2,i,j,k)=Fz(2,i,j,k)+tmp(1,n,i,j,k)*ciz(n)*constants%g
!!$                   Fz(2,i,j,k)=Fz(2,i,j,k)+tmp(2,n,i,j,k)*ciz(n)*constants%g22
!!$                enddo
!!$                do n=7,18
!!$                   Fz(1,i,j,k)=Fz(1,i,j,k)+tmp(2,n,i,j,k)*ciz(n)*constants%g*0.5
!!$                   Fz(1,i,j,k)=Fz(1,i,j,k)+tmp(1,n,i,j,k)*ciz(n)*constants%g11*0.5
!!$                   Fz(2,i,j,k)=Fz(2,i,j,k)+tmp(1,n,i,j,k)*ciz(n)*constants%g*0.5
!!$                   Fz(2,i,j,k)=Fz(2,i,j,k)+tmp(2,n,i,j,k)*ciz(n)*constants%g22*0.5
!!$                enddo
!!$             end if
!!$          end do
!!$       end do
!!$    end do

!!$    ! code that replaces the code above
!!$    do i=info%xs,info%xe 
!!$      do j=info%ys,info%ye 
!!$        do k=info%zs,info%ze
!!$          if (walls(i,j,k).eq.0) then
!!$            do n=1,6
!!$              Fx(1,i,j,k)=Fx(1,i,j,k)+tmp(2,n,i,j,k)*cix(i)*constants%g
!!$              Fx(1,i,j,k)=Fx(1,i,j,k)+tmp(1,n,i,j,k)*cix(i)*constants%g11
!!$              Fy(1,i,j,k)=Fy(1,i,j,k)+tmp(2,n,i,j,k)*ciy(i)*constants%g
!!$              Fy(1,i,j,k)=Fy(1,i,j,k)+tmp(1,n,i,j,k)*ciy(i)*constants%g11
!!$              Fz(1,i,j,k)=Fz(1,i,j,k)+tmp(2,n,i,j,k)*ciz(i)*constants%g
!!$              Fz(1,i,j,k)=Fz(1,i,j,k)+tmp(1,n,i,j,k)*ciz(i)*constants%g11
!!$              Fx(2,i,j,k)=Fx(2,i,j,k)+tmp(1,n,i,j,k)*cix(i)*constants%g
!!$              Fx(2,i,j,k)=Fx(2,i,j,k)+tmp(2,n,i,j,k)*cix(i)*constants%g22
!!$              Fy(2,i,j,k)=Fy(2,i,j,k)+tmp(1,n,i,j,k)*ciy(i)*constants%g
!!$              Fy(2,i,j,k)=Fy(2,i,j,k)+tmp(2,n,i,j,k)*ciy(i)*constants%g22
!!$              Fz(2,i,j,k)=Fz(2,i,j,k)+tmp(1,n,i,j,k)*ciz(i)*constants%g
!!$              Fz(2,i,j,k)=Fz(2,i,j,k)+tmp(2,n,i,j,k)*ciz(i)*constants%g22 
!!$            enddo
!!$            do n=7,18
!!$              Fx(1,i,j,k)=Fx(1,i,j,k)+tmp(2,n,i,j,k)*cix(i)*constants%g*0.5
!!$              Fx(1,i,j,k)=Fx(1,i,j,k)+tmp(1,n,i,j,k)*cix(i)*constants%g11*0.5
!!$              Fy(1,i,j,k)=Fy(1,i,j,k)+tmp(2,n,i,j,k)*ciy(i)*constants%g*0.5
!!$              Fy(1,i,j,k)=Fy(1,i,j,k)+tmp(1,n,i,j,k)*ciy(i)*constants%g11*0.5
!!$              Fz(1,i,j,k)=Fz(1,i,j,k)+tmp(2,n,i,j,k)*ciz(i)*constants%g*0.5
!!$              Fz(1,i,j,k)=Fz(1,i,j,k)+tmp(1,n,i,j,k)*ciz(i)*constants%g11*0.5
!!$              Fx(2,i,j,k)=Fx(2,i,j,k)+tmp(1,n,i,j,k)*cix(i)*constants%g*0.5
!!$              Fx(2,i,j,k)=Fx(2,i,j,k)+tmp(2,n,i,j,k)*cix(i)*constants%g22*0.5
!!$              Fy(2,i,j,k)=Fy(2,i,j,k)+tmp(1,n,i,j,k)*ciy(i)*constants%g*0.5
!!$              Fy(2,i,j,k)=Fy(2,i,j,k)+tmp(2,n,i,j,k)*ciy(i)*constants%g22*0.5
!!$              Fz(2,i,j,k)=Fz(2,i,j,k)+tmp(1,n,i,j,k)*ciz(i)*constants%g*0.5
!!$              Fz(2,i,j,k)=Fz(2,i,j,k)+tmp(2,n,i,j,k)*ciz(i)*constants%g22*0.5
!!$            enddo
!!$          end if
!!$        end do
!!$      end do
!!$    end do


    ! this code is way more efficient than the existing code,
    ! but we'll go by the "make it work, then make it correct,
    ! then make it fast" approach.  Note that similar code would
    ! have to be added for Fy,Fz

    !       do i=info%xs,info%xe
    !         do j=info%ys,info%ye
    !           do k=info%zs,info%ze
    !             if (walls(i,j,k).eq.0) then
    !                 Fx(1,i,j,k)=Fx(1,i,j,k)
    !      &              + rho(1,i+1,j,k)*cix(1)*g11
    !      &              + rho(1,i-1,j,k)*cix(3)*g11
    !      &              + rho(2,i+1,j,k)*cix(1)*g
    !      &              + rho(2,i-1,j,k)*cix(3)*g
    !      &              + rho(1,i+1,j+1,k)*cix(7)*g11*.5
    !      &              + rho(1,i-1,j+1,k)*cix(8)*g11*.5
    !      &              + rho(1,i-1,j-1,k)*cix(9)*g11*.5
    !      &              + rho(1,i+1,j-1,k)*cix(10)*g11*.5
    !      &              + rho(2,i+1,j+1,k)*cix(7)*g*.5
    !      &              + rho(2,i-1,j+1,k)*cix(8)*g*.5
    !      &              + rho(2,i-1,j-1,k)*cix(9)*g*.5
    !      &              + rho(2,i+1,j-1,k)*cix(10)*g*.5
    !      &              + rho(1,i+1,j,k+1)*cix(11)*g11*.5
    !      &              + rho(1,i-1,j,k+1)*cix(12)*g11*.5
    !      &              + rho(1,i-1,j,k-1)*cix(13)*g11*.5
    !      &              + rho(1,i+1,j,k-1)*cix(14)*g11*.5
    !      &              + rho(2,i+1,j,k+1)*cix(11)*g11*.5
    !      &              + rho(2,i-1,j,k+1)*cix(12)*g11*.5
    !      &              + rho(2,i-1,j,k-1)*cix(13)*g11*.5
    !      &              + rho(2,i+1,j,k-1)*cix(14)*g11*.5

    !                 Fx(2,i,j,k)=Fx(2,i,j,k)
    !      &              + rho(2,i+1,j,k)*cix(1)*g22
    !      &              + rho(2,i-1,j,k)*cix(3)*g22
    !      &              + rho(1,i+1,j,k)*cix(1)*g
    !      &              + rho(1,i-1,j,k)*cix(3)*g
    !      &              + rho(2,i+1,j+1,k)*cix(7)*g22*.5
    !      &              + rho(2,i-1,j+1,k)*cix(8)*g22*.5
    !      &              + rho(2,i-1,j-1,k)*cix(9)*g22*.5
    !      &              + rho(2,i+1,j-1,k)*cix(10)*g22*.5
    !      &              + rho(1,i+1,j+1,k)*cix(7)*g*.5
    !      &              + rho(1,i-1,j+1,k)*cix(8)*g*.5
    !      &              + rho(1,i-1,j-1,k)*cix(9)*g*.5
    !      &              + rho(1,i+1,j-1,k)*cix(10)*g*.5
    !      &              + rho(2,i+1,j,k+1)*cix(11)*g22*.5
    !      &              + rho(2,i-1,j,k+1)*cix(12)*g22*.5
    !      &              + rho(2,i-1,j,k-1)*cix(13)*g22*.5
    !      &              + rho(2,i+1,j,k-1)*cix(14)*g22*.5
    !      &              + rho(1,i+1,j,k+1)*cix(11)*g22*.5
    !      &              + rho(1,i-1,j,k+1)*cix(12)*g22*.5
    !      &              + rho(1,i-1,j,k-1)*cix(13)*g22*.5
    !      &              + rho(1,i+1,j,k-1)*cix(14)*g22*.5
    return
  end subroutine LBMAddFluidFluidForces


  ! --- body forces on fluid
  subroutine LBMAddBodyForces(rho,Fx,Fy,Fz,walls,info,constants)
    type(info_type) info
    type(constants_type) constants
    PetscScalar,dimension(1:info%s,info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: rho
    PetscScalar,dimension(1:info%s, info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: Fx,Fy,Fz
    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: walls

    integer i,j,k,m,n

    do m=1,info%s
       where( (walls.eq.0) )
          Fz(m,:,:,:) = Fz(m,:,:,:) + constants%gvt(m)*constants%mm(m)*rho(m,:,:,:)
       end where
    enddo
    return
  end subroutine LBMAddBodyForces

  ! --- Fluid-solid interaction forces, from
  ! ---  (Kang 2002 Eq. 8)
  ! --- NOTE -- I'm ignoring the bug that the forces will
  !             wrap even in the non-periodic case.  This seems
  !             save at this point, since it seems unlikely that
  !             there will be walls on one side but not the other
  subroutine LBMAddFluidSolidForces(rho,Fx,Fy,Fz,walls,info,constants)
    !     NONLOCAL IN WALLS

!!$    ! input
!!$    type(info_type) info
!!$    type(constants_type) constants
!!$    PetscScalar,dimension(1:info%s,info%gxs:info%gxe, info%gys:info%gye, &
!!$         info%gzs:info%gze):: rho
!!$    PetscScalar,dimension(1:info%s, info%gxs:info%gxe, info%gys:info%gye, &
!!$         info%gzs:info%gze):: Fx,Fy,Fz
!!$    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye, &
!!$         info%gzs:info%gze):: walls
!!$
!!$    ! local
!!$    integer i,j,k,m,n
!!$
!!$    do m=1,info%s
!!$       where( (walls.eq.0) )
!!$          Fx(m,:,:,:)=Fx(m,:,:,:) - rho(m,:,:,:)*constants%gw(m) &
!!$               *(cshift(walls, shift= 1, dim=1) - cshift(walls, shift=-1, dim=1))
!!$
!!$          Fy(m,:,:,:)=Fy(m,:,:,:) - rho(m,:,:,:)*constants%gw(m) &
!!$               *(cshift(walls, shift= 1, dim=2) - cshift(walls, shift=-1, dim=2))
!!$
!!$          Fz(m,:,:,:)=Fz(m,:,:,:) - rho(m,:,:,:)*constants%gw(m) &
!!$               *(cshift(walls, shift= 1, dim=3) - cshift(walls, shift=-1, dim=3))
!!$       end where
!!$    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! New code added by MLP.  This accounts for fluid-solid forces on the diagonals, which is
    ! not accounted for in the other version.
 
    use LBM_Discretization_D3Q19_module

    ! input
    type(info_type) info
    type(constants_type) constants
    PetscScalar,dimension(1:info%s,info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: rho
    PetscScalar,dimension(1:info%s, info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: Fx,Fy,Fz
    PetscScalar,dimension(info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: walls


    ! local
    integer i,j,k,m,n
    PetscScalar,dimension(0:info%b,info%gxs:info%gxe, &
         info%gys:info%gye,info%gzs:info%gze)::tmp
        
    tmp(1,:,:,:)= (cshift(walls(:,:,:),shift= 1,dim=1))
    tmp(2,:,:,:)= (cshift(walls(:,:,:),shift= 1,dim=2))
    tmp(3,:,:,:)= (cshift(walls(:,:,:),shift=-1,dim=1))
    tmp(4,:,:,:)= (cshift(walls(:,:,:),shift=-1,dim=2))
    tmp(5,:,:,:)= (cshift(walls(:,:,:),shift= 1,dim=3))
    tmp(6,:,:,:)= (cshift(walls(:,:,:),shift=-1,dim=3))

    tmp(7 ,:,:,:)= (cshift(tmp(1,:,:,:),shift= 1,dim=2))
    tmp(8 ,:,:,:)= (cshift(tmp(3,:,:,:),shift= 1,dim=2))
    tmp(9 ,:,:,:)= (cshift(tmp(3,:,:,:),shift=-1,dim=2))
    tmp(10,:,:,:)= (cshift(tmp(1,:,:,:),shift=-1,dim=2))

    tmp(11,:,:,:)= (cshift(tmp(5,:,:,:),shift= 1,dim=1))
    tmp(12,:,:,:)= (cshift(tmp(5,:,:,:),shift=-1,dim=1))
    tmp(13,:,:,:)= (cshift(tmp(6,:,:,:),shift=-1,dim=1))
    tmp(14,:,:,:)= (cshift(tmp(6,:,:,:),shift= 1,dim=1))

    tmp(15,:,:,:)= (cshift(tmp(5,:,:,:),shift= 1,dim=2))
    tmp(16,:,:,:)= (cshift(tmp(5,:,:,:),shift=-1,dim=2))
    tmp(17,:,:,:)= (cshift(tmp(6,:,:,:),shift=-1,dim=2))
    tmp(18,:,:,:)= (cshift(tmp(6,:,:,:),shift= 1,dim=2))


    do m=1,info%s
      do i=info%xs,info%xe 
        do j=info%ys,info%ye 
          do k=info%zs,info%ze
            if (walls(i,j,k).eq.0) then
              do n=1,6
                Fx(m,i,j,k)=Fx(m,i,j,k)-rho(m,i,j,k)*tmp(n,i,j,k)*cix(n)*constants%gw(m)
                Fy(m,i,j,k)=Fy(m,i,j,k)-rho(m,i,j,k)*tmp(n,i,j,k)*ciy(n)*constants%gw(m)
                Fz(m,i,j,k)=Fz(m,i,j,k)-rho(m,i,j,k)*tmp(n,i,j,k)*ciz(n)*constants%gw(m)
              enddo
              do n=7,18
                Fx(m,i,j,k)=Fx(m,i,j,k)-rho(m,i,j,k)*tmp(n,i,j,k)*cix(n)*constants%gw(m)*0.5
                Fy(m,i,j,k)=Fy(m,i,j,k)-rho(m,i,j,k)*tmp(n,i,j,k)*ciy(n)*constants%gw(m)*0.5
	      	Fz(m,i,j,k)=Fz(m,i,j,k)-rho(m,i,j,k)*tmp(n,i,j,k)*ciz(n)*constants%gw(m)*0.5
              enddo
            end if
            !write(*,*) Fz(1,i,j,k),Fz(2,i,j,k)
          end do
        end do
      end do
    end do

     
!!$    write(*,*) i-1,j-1,k-1
!!$    write(*,*) Fx(1,i-1,j-1,k-1),Fx(2,i-1,j-1,k-1)
!!$    write(*,*) Fy(1,i-1,j-1,k-1),Fy(2,i-1,j-1,k-1)
!!$    write(*,*) Fz(1,i-1,j-1,k-1),Fz(2,i-1,j-1,k-1)
!!$    write(*,*) tmp(1,i-1,j-1,k-1)
    return
  end subroutine LBMAddFluidSolidForces


  ! --- zero out boundaries, as they affect the fi
  subroutine LBMZeroBoundaryForces(bc_flags, Fx, Fy, Fz, bc_dim, info,constants)
    type(info_type) info
    type(constants_type) constants
    integer bc_dim
    integer, dimension(6)::bc_flags ! enum for boundary conditions
    PetscScalar,dimension(1:info%s, info%gxs:info%gxe, info%gys:info%gye, &
         info%gzs:info%gze):: Fx,Fy,Fz

    ! -- x
    if ((info%xs.eq.1).and.((bc_flags(1).eq.2).or.(bc_flags(1).eq.3))) then
       Fx(:,1,:,:) = 0
       Fy(:,1,:,:) = 0
       Fz(:,1,:,:) = 0
    endif

    if ((info%xe.eq.info%NX).and.((bc_flags(2).eq.2).or.(bc_flags(2).eq.3))) then
       Fx(:,info%NX,:,:) = 0
       Fy(:,info%NX,:,:) = 0
       Fz(:,info%NX,:,:) = 0
    endif

    ! -- y
    if ((info%ys.eq.1).and.((bc_flags(3).eq.2).or.(bc_flags(3).eq.3))) then
       Fx(:,:,1,:) = 0
       Fy(:,:,1,:) = 0
       Fz(:,:,1,:) = 0
    endif

    if ((info%ye.eq.info%NY).and.((bc_flags(4).eq.2).or.(bc_flags(4).eq.3))) then
       Fx(:,:,info%NY,:) = 0
       Fy(:,:,info%NY,:) = 0
       Fz(:,:,info%NY,:) = 0
    endif

    ! -- z
    if ((info%zs.eq.1).and.((bc_flags(5).eq.2).or.(bc_flags(5).eq.3))) then
       Fx(:,:,:,1) = 0
       Fy(:,:,:,1) = 0
       Fz(:,:,:,1) = 0
    endif

    if ((info%ze.eq.info%NZ).and.((bc_flags(6).eq.2).or.(bc_flags(6).eq.3))) then
       Fx(:,:,:,info%NZ) = 0
       Fy(:,:,:,info%NZ) = 0
       Fz(:,:,:,info%NZ) = 0
    endif

    return
  end subroutine LBMZeroBoundaryForces
end module LBM_Forcing_module
