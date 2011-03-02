!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        CreateGeom3DParallel.f90
!!!     version:
!!!     created:         16 December 2010
!!!       on:            14:14:04 MST
!!!     last modified:   01 March 2011
!!!       at:            13:26:28 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscdmdef.h"
#include "finclude/petscviewerdef.h"

program geometry
  use petsc
  implicit none

#include "lbm_definitions.h"

  PetscInt NX, NY, NZ
  PetscScalar,pointer,dimension(:,:,:):: walls
  PetscInt, allocatable, dimension(:,:,:)::pores
  character(len=MAXSTRINGLENGTH) flnm1, flnm2
  character(len=MAXSTRINGLENGTH) optionsfile
  PetscBool flag
  PetscInt i,j,k

  PetscScalar poro

  DM da
  Vec walls_g
  PetscErrorCode ierr
  PetscViewer viewer

  call getarg(1, optionsfile)
  call PetscInitialize(optionsfile, ierr)

  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-nx',NX,flag,ierr)
  if (.not.flag) NX=3
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-ny',NY,flag,ierr)
  if (.not.flag) NY=65
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-nz',NZ,flag,ierr)
  if (.not.flag)   NZ=200 !default values, may be altered


  call DMDACreate3d(PETSC_COMM_SELF, DMDA_XYZPERIODIC, DMDA_STENCIL_BOX, NX, NY, NZ, &
       PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 0, &
       PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, da, ierr)
  call DMSetFromOptions(da, ierr)
  call DMCreateGlobalVector(da, walls_g, ierr)
  call DMDAVecGetArrayF90(da, walls_g, walls, ierr)

  call DMDAGetInfo(da, PETSC_NULL_INTEGER, NX, NY, NZ, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
       PETSC_NULL_INTEGER, ierr)

  walls = 0.0
  allocate(pores(NX,NY,NZ))

  call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-input_filename',flnm1,flag,ierr)
  if (.not.flag) flnm1 = 'StaggCyl.txt'
  call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-output_filename',flnm2,flag,ierr)
  if (.not.flag) flnm2 = 'mygeom.000'

99405 FORMAT(1000I2)
  open(9,file=flnm1)

  do i=1,NX
     do j=1,NY
        do k=1,NZ
           read(9,99405) pores(i,j,k)
        enddo
     enddo
  enddo
  close(9)

  where (pores(:,:,:).eq.1) 
     walls(:,:,:) = 1.0
  end where

  poro=sum(pores)/float(NX*NY*NZ)
  write(*,*) 'poro=', 1-poro

  call DMDAVecRestoreArrayF90(da, walls_g, walls, ierr)

  call PetscViewerBinaryOpen(PETSC_COMM_SELF, flnm2, FILE_MODE_WRITE, viewer, ierr)
  call VecView(walls_g, viewer, ierr)
  call PetscViewerDestroy(viewer, ierr)

  call VecDestroy(walls_g, ierr)
  call DMDestroy(da, ierr)
  call PetscFinalize(ierr)

end program geometry
