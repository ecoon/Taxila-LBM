!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_error.F90
!!!     version:
!!!     created:         23 January 2012
!!!       on:            14:16:32 MST
!!!     created:         23 January 2012
!!!       at:            16:31:11 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petscdmdef.h"

module LBM_Error_module
  use petsc
  implicit none

  private
#include "lbm_definitions.h"

  public LBMError, &
       LBMWarn

contains
  subroutine LBMError(comm, code, message, ierr)
    MPI_Comm comm
    character(len=*):: message
    PetscInt code
    PetscErrorCode ierr

    PetscMPIInt rank, size
    character(len=MAXSTRINGLENGTH):: fullmessage

    call mpi_comm_rank(comm, rank, ierr)
    call mpi_comm_size(comm, size, ierr)
    write(fullmessage, "(a,I0,a,I0,a)"), "ERROR [", rank, "/", size, "]:"//message
    print*, fullmessage
    SETERRQ(comm, code, message, ierr)
    return
  end subroutine LBMError

  subroutine LBMWarn(comm, message, ierr)
    MPI_Comm comm
    character(len=*):: message
    PetscErrorCode ierr

    PetscMPIInt rank, size
    character(len=MAXSTRINGLENGTH):: fullmessage

    call mpi_comm_rank(comm, rank, ierr)
    if (rank.eq.0) then
      write(fullmessage, "(a,I0,a,I0,a)"), "WARNING:"//message
      print*, fullmessage
    end if
    return
  end subroutine LBMWarn
end module LBM_Error_module
