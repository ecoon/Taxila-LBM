!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_io.F90
!!!     version:         
!!!     created:         29 March 2011
!!!       on:            15:38:36 MDT
!!!     last modified:   30 March 2011
!!!       at:            12:40:01 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscviewerdef.h"

module LBM_IO_module
  use petsc
  implicit none

  private
#include "lbm_definitions.h"

  type, public:: io_type
     MPI_Comm comm
     character(len=MAXSTRINGLENGTH) prefix
     PetscBool mpiio
     PetscInt kwrite
     PetscInt counter
     PetscInt prefix_len
  end type io_type

  public:: IOCreate, &
       IODestroy, &
       IOSetFromOptions, &
       IOView, &
       IOLoad, &
       IOIncrementCounter

contains
  function IOCreate(comm) result(io)
    type(io_type),pointer:: io
    MPI_Comm comm
    allocate(io)
    io%comm = comm
    io%counter = 0
    io%prefix_len = -1
    io%prefix = ''
  end function IOCreate

  subroutine IODestroy(io, ierr)
    type(io_type) io
    PetscErrorCode ierr
  end subroutine IODestroy

  subroutine IOSetFromOptions(io, options, ierr)
    use LBM_Options_module
    type(io_type) io
    type(options_type) options
    PetscErrorCode ierr

    io%kwrite = options%kwrite
    io%prefix = options%output_prefix
    io%prefix_len = LEN_TRIM(io%prefix)
    io%mpiio = options%mpiio
  end subroutine IOSetFromOptions

  subroutine IOView(io, vec, name)
    type(io_type) io
    Vec vec
    character(len=*) name
    character(len=MAXIODIGITS):: outnum
    character(len=MAXSTRINGLENGTH):: stringformat
    character(len=MAXSTRINGLENGTH):: filename
    PetscErrorCode ierr
    PetscViewer viewer

    write(stringformat, '("(I0.",I1,")")'), MAXIODIGITS
    write(outnum, stringformat) io%counter
    filename = io%prefix(1:io%prefix_len)//trim(name)//outnum//'.dat'

    call PetscViewerCreate(io%comm, viewer, ierr)
    call PetscViewerSetType(viewer, PETSCVIEWERBINARY, ierr)
    if (io%mpiio) call PetscViewerBinarySetMPIIO(viewer, ierr)
    call PetscViewerFileSetMode(viewer, FILE_MODE_WRITE, ierr)
    call PetscViewerFileSetName(viewer, filename, ierr)
    call VecView(vec, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)
    CHKERRQ(ierr)
  end subroutine IOView

  subroutine IOLoad(io, vec, name)
    type(io_type) io
    Vec vec
    character(len=*) name
    character(len=MAXIODIGITS):: outnum
    character(len=MAXSTRINGLENGTH):: stringformat
    character(len=MAXSTRINGLENGTH):: filename
    PetscErrorCode ierr
    PetscViewer viewer

    write(stringformat, '("(I0.",I1,")")'), MAXIODIGITS
    write(outnum, stringformat) io%counter
    filename = io%prefix(1:io%prefix_len)//trim(name)//outnum//'.dat'

    call PetscViewerCreate(io%comm, viewer, ierr)
    call PetscViewerSetType(viewer, PETSCVIEWERBINARY, ierr)
    if (io%mpiio) call PetscViewerBinarySetMPIIO(viewer, ierr)
    call PetscViewerFileSetMode(viewer, FILE_MODE_READ, ierr)
    call PetscViewerFileSetName(viewer, filename, ierr)
    call VecLoad(vec, viewer, ierr)
    call PetscViewerDestroy(viewer,ierr)
    CHKERRQ(ierr)
  end subroutine IOLoad

  subroutine IOIncrementCounter(io)
    type(io_type) io
    io%counter = io%counter + 1
  end subroutine IOIncrementCounter
end module LBM_IO_module
