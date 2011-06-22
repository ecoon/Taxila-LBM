!!! ====================================================================
!!!  Fortran-90-file
!!!     author:          Ethan T. Coon
!!!     filename:        lbm_logging.F90
!!!     version:         
!!!     created:         22 June 2011
!!!       on:            09:37:52 MDT
!!!     last modified:   22 June 2011
!!!       at:            10:52:39 MDT
!!!     URL:             http://www.ldeo.columbia.edu/~ecoon/
!!!     email:           ecoon _at_ lanl.gov
!!!  
!!! ====================================================================
#define PETSC_USE_FORTRAN_MODULES 1
#include "finclude/petscsysdef.h"

module LBM_Logging_module
  implicit none

  private
#include "lbm_definitions.h"

  PetscInt, parameter, public :: INIT_STAGE = 1
  PetscInt, parameter, public :: STREAM_STAGE = 2
  PetscInt, parameter, public :: BC_STAGE = 3
  PetscInt, parameter, public :: FORCING_STAGE = 4
  PetscInt, parameter, public :: COLLISION_STAGE = 5
  PetscInt, parameter, public :: OUTPUT_STAGE = 6
  PetscInt, parameter, public :: COMMUNICATION_STAGE = 7
  PetscInt, parameter, public :: DESTROY_STAGE = 8
  
  type, public :: log_type
    PetscLogStage :: stage(10)
    PetscClassId :: class_lbm

    PetscLogEvent :: event_init_options
    PetscLogEvent :: event_init_setfromoptions
    PetscLogEvent :: event_init_gridsetup
    PetscLogEvent :: event_init_wallssetup
    PetscLogEvent :: event_init_flowsetup
    PetscLogEvent :: event_init_transetup
    PetscLogEvent :: event_init_bcsetup
    PetscLogEvent :: event_init_icsetup

    PetscLogEvent :: event_moments
    PetscLogEvent :: event_diagnostics

    PetscLogEvent :: event_stream_flow
    PetscLogEvent :: event_stream_tran

    PetscLogEvent :: event_bc_bounceback
    PetscLogEvent :: event_bc_tranwallreact
    PetscLogEvent :: event_bc_flow
    PetscLogEvent :: event_bc_tran
    
    PetscLogEvent :: event_forcing_fluidfluid
    PetscLogEvent :: event_forcing_fluidsolid
    PetscLogEvent :: event_forcing_body

    PetscLogEvent :: event_collision_flow
    PetscLogEvent :: event_collision_tran

    PetscLogEvent :: event_communicate_fi
    PetscLogEvent :: event_communicate_rho

    PetscLogEvent :: event_output
    PetscLogEvent :: event_output_grid
    PetscLogEvent :: event_output_walls
    PetscLogEvent :: event_cleanup

    
  end type log_type

  type(log_type), pointer, public :: logger

  public :: LoggerCreate, &
       LoggerDestroy

contains
  subroutine LoggerCreate()
    PetscErrorCode ierr
    allocate(logger)
    call PetscLogStageRegister('Init Stage', logger%stage(INIT_STAGE), ierr)
    call PetscLogStageRegister('Streaming Stage', logger%stage(STREAM_STAGE), ierr)
    call PetscLogStageRegister('BC Stage', logger%stage(BC_STAGE), ierr)
    call PetscLogStageRegister('Forcing Stage', logger%stage(FORCING_STAGE), ierr)
    call PetscLogStageRegister('Collision Stage', logger%stage(COLLISION_STAGE), ierr)
    call PetscLogStageRegister('Comm Stage', &
         logger%stage(COMMUNICATION_STAGE), ierr)
    call PetscLogStageRegister('Output Stage', logger%stage(OUTPUT_STAGE), ierr)
    call PetscLogStageRegister('Destroy Stage', logger%stage(DESTROY_STAGE), ierr)

    call PetscClassIdRegister('LBM', logger%class_lbm, ierr)

    call PetscLogEventRegister('Init options', logger%class_lbm, &
         logger%event_init_options, ierr)
    call PetscLogEventRegister('Init SetFromOptions()', logger%class_lbm, &
         logger%event_init_setfromoptions, ierr)
    call PetscLogEventRegister('Init Grid', logger%class_lbm, &
         logger%event_init_gridsetup, ierr)
    call PetscLogEventRegister('Init Walls', logger%class_lbm, &
         logger%event_init_wallssetup, ierr)
    call PetscLogEventRegister('Init Flow', logger%class_lbm, &
         logger%event_init_flowsetup, ierr)
    call PetscLogEventRegister('Init Transport', logger%class_lbm, &
         logger%event_init_transetup, ierr)
    call PetscLogEventRegister('Init BCs', logger%class_lbm, &
         logger%event_init_bcsetup, ierr)
    call PetscLogEventRegister('Init ICs', logger%class_lbm, &
         logger%event_init_icsetup, ierr)

    call PetscLogEventRegister('Update Moments', logger%class_lbm, &
         logger%event_moments, ierr)
    call PetscLogEventRegister('Update Soln', logger%class_lbm, &
         logger%event_diagnostics, ierr)

    call PetscLogEventRegister('Stream Flow', logger%class_lbm, &
         logger%event_stream_flow, ierr)
    call PetscLogEventRegister('Stream Transport', logger%class_lbm, &
         logger%event_stream_tran, ierr)

    call PetscLogEventRegister('BC Flow Bounceback', logger%class_lbm, &
         logger%event_bc_bounceback, ierr)
    call PetscLogEventRegister('BC Transport React with Walls', logger%class_lbm, &
         logger%event_bc_tranwallreact, ierr)
    call PetscLogEventRegister('BC Exterior Flow', logger%class_lbm, &
         logger%event_bc_flow, ierr)
    call PetscLogEventRegister('BC Exterior Transport', logger%class_lbm, &
         logger%event_bc_tran, ierr)

    call PetscLogEventRegister('Fluid-Fluid Forcing', logger%class_lbm, &
         logger%event_forcing_fluidfluid, ierr)
    call PetscLogEventRegister('Fluid-Solid Forcing', logger%class_lbm, &
         logger%event_forcing_fluidsolid, ierr)
    call PetscLogEventRegister('Body Forces', logger%class_lbm, &
         logger%event_forcing_body, ierr)

    call PetscLogEventRegister('Collision Flow', logger%class_lbm, &
         logger%event_collision_flow, ierr)
    call PetscLogEventRegister('Collision Transport', logger%class_lbm, &
         logger%event_collision_tran, ierr)

    call PetscLogEventRegister('Communicate f_i', logger%class_lbm, &
         logger%event_communicate_fi, ierr)
    call PetscLogEventRegister('Communicate density', logger%class_lbm, &
         logger%event_communicate_rho, ierr)

    call PetscLogEventRegister('Output Soln', logger%class_lbm, &
         logger%event_output, ierr)
    call PetscLogEventRegister('Output Grid', logger%class_lbm, &
         logger%event_output_grid, ierr)
    call PetscLogEventRegister('Output Walls', logger%class_lbm, &
         logger%event_output_walls, ierr)

    call PetscLogEventRegister('Destruction', logger%class_lbm, &
         logger%event_cleanup, ierr)
  end subroutine LoggerCreate
  
  subroutine LoggerDestroy()
    deallocate(logger)
    nullify(logger)
  end subroutine LoggerDestroy
end module LBM_Logging_module
