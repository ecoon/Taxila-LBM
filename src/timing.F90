!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        timing.F90
!!!     version:         
!!!     created:         08 December 2010
!!!       on:            14:35:16 MST
!!!     last modified:   08 December 2010
!!!       at:            15:11:08 MST
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
  
  module Timing_module
    use petsc
    implicit none

    type timing_type
       MPI_Comm comm
       integer time_start
       integer time_end
       integer rate
       integer max
       real time
       character name(60)
    end type timing_type

  contains

    function TimingCreate(comm, name) result(timing)
      implicit none
      MPI_Comm comm
      type(timing_type),pointer:: timing
      character name(60)

      allocate(timing)
      timing%comm = comm
      timing%name = name
      timing%time_end = 0
      timing%time = 0

      call system_clock (timing%time_start, timing%rate, timing%max)
    end function TimingCreate

    subroutine TimingEnd(timing)
      implicit none
      type(timing_type) timing
      integer id, nprocs
      integer ierr
      real mean_time

      call system_clock (timing%time_end, timing%rate, timing%max)
      timing%time = real(timing%time_end - timing%time_start)/real(timing%rate)
      call MPI_Comm_Rank(timing%comm, id, ierr)
      
      if (id.eq.0) then
         call MPI_Comm_Size(timing%comm, nprocs, ierr)
         call MPI_Reduce(timing%time, mean_time, 1, MPI_REAL, MPI_SUM, id, timing%comm, ierr)
         write(*,*) 'Timing:', timing%name
         write(*,*) '  took (s):', mean_time/real(nprocs)
      end if
    end subroutine TimingEnd
  end module Timing_module
