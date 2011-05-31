!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        timing.F90
!!!     version:         
!!!     created:         08 December 2010
!!!       on:            14:35:16 MST
!!!     last modified:   25 May 2011
!!!       at:            15:31:15 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"
  
module Timing_module
  use petsc
  implicit none
  
  private
#include "lbm_definitions.h"
  
  type,public:: timing_type
     MPI_Comm comm
     integer time_start
     integer time_end
     integer rate
     integer max
     real time
     character(len=MAXWORDLENGTH) name
  end type timing_type
  
  interface TimingEndPerUnit
     module procedure TimingEndPerUnit1
     module procedure TimingEndPerUnit2
  end interface
     
  public :: TimingCreate, &
          TimingDestroy, &
         TimingEnd, &
         TimingEndPerUnit

  contains

    function TimingCreate(comm, name) result(timing)
      MPI_Comm comm
      type(timing_type),pointer:: timing
      character(len=MAXWORDLENGTH) name

      allocate(timing)
      timing%comm = comm
      timing%name = name
      timing%time_end = 0
      timing%time = 0

      call system_clock (timing%time_start, timing%rate, timing%max)
    end function TimingCreate

    subroutine TimingDestroy(timing)
      type(timing_type),pointer:: timing
      deallocate(timing)
    end subroutine TimingDestroy

    subroutine TimingEnd(timing)
      type(timing_type) timing
      integer id, nprocs
      integer ierr
      real mean_time

      call system_clock (timing%time_end, timing%rate, timing%max)
      timing%time = real(timing%time_end - timing%time_start)/real(timing%rate)
      call MPI_Comm_Rank(timing%comm, id, ierr)
      call MPI_Reduce(timing%time, mean_time, 1, MPI_REAL, MPI_SUM, 0, &
           timing%comm, ierr)
      if (id.eq.0) then
         call MPI_Comm_Size(timing%comm, nprocs, ierr)
         write(*,*) 'Timing:', timing%name
         write(*,*) '  took (s):', mean_time/real(nprocs)
      end if
    end subroutine TimingEnd

    subroutine TimingEndPerUnit1(timing, numunits, unitname)
      type(timing_type) timing
      integer numunits
      character(len=MAXWORDLENGTH) unitname
      call TimingEndPerUnit2(timing, numunits, unitname, .FALSE.)
    end subroutine TimingEndPerUnit1

    subroutine TimingEndPerUnit2(timing, numunits, unitname, supress_output)
      ! same as TimingEnd, this just divides by a factor of the
      ! number of units, for instance the number of timesteps
      type(timing_type) timing
      integer numunits
      character(len=MAXWORDLENGTH) unitname
      logical supress_output

      integer id, nprocs
      integer ierr
      real mean_time
      integer charlen

      call system_clock (timing%time_end, timing%rate, timing%max)
      timing%time = real(timing%time_end - timing%time_start)/real(timing%rate)
      call MPI_Comm_Rank(timing%comm, id, ierr)
      
      call MPI_Reduce(timing%time, mean_time, 1, MPI_REAL, MPI_SUM, 0, timing%comm, ierr)
      if (id.eq.0) then
         call MPI_Comm_Size(timing%comm, nprocs, ierr)
         charlen = LEN_TRIM(unitname)
         if (.not.supress_output) then
            write(*,*) 'Timing:', timing%name
            write(*,*) '  took (per unit', unitname(1:charlen), ') (s):', mean_time/real(nprocs)/real(numunits)
         end if
      end if
    end subroutine TimingEndPerUnit2

  end module Timing_module
