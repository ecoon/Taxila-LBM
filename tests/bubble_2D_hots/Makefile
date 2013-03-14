ifndef LBM_DIR
  LBM_DIR:=$(shell pwd)/../..
endif

ifndef COMP
  HOST:= $(shell hostname | cut -d"-" -f1)
  ifeq ($(HOST),jaguar)
    COMP:=ftn
  else
    COMP:= $(shell basename $(shell mpif90 -show | cut -d" " -f 1))
  endif
endif

LBM_SRC_DIR= $(LBM_DIR)/src/lbm
LBM_PROBLEM_DIR= $(LBM_DIR)/src/problem_specs
LBM_INCLUDE_DIR= $(LBM_DIR)/include
LBM_LIB_DIR= $(LBM_DIR)/lib
LBM_MOD_DIR= $(LBM_LIB_DIR)/modules
LBM_LIB_NAME= lbm.a

FFLAGS = -I$(LBM_INCLUDE_DIR) 

# module location
FFLAGS += -I$(LBM_MOD_DIR)

# optimization flags
ifndef LBM_OPT
  ifeq ($(COMP),gfortran)
    FFLAGS += -fbounds-check 
  else ifeq ($(COMP),pgf90)
    FFLAGS += -Mbounds
  else ifeq ($(COMP),pgf95)
    FFLAGS += -Mbounds
  endif
else
  FFLAGS += -O3
endif

MAIN = $(LBM_SRC_DIR)/main.o 


ifeq ($(strip $(wildcard initialize_state.F90)),)
  FORTRANOBJS = $(LBM_PROBLEM_DIR)/default_initialize_state.o
else
  FORTRANOBJS = initialize_state.o
endif
ifeq ($(strip $(wildcard initialize_walls.F90)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_walls.o
else
  FORTRANOBJS += initialize_walls.o
endif
ifeq ($(strip $(wildcard initialize_bcs.F90)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_bcs.o
else
  FORTRANOBJS += initialize_bcs.o
endif
ifeq ($(strip $(wildcard initialize_bcs_transport.F90)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_bcs_transport.o
else
  FORTRANOBJS += initialize_bcs_transport.o
endif
ifeq ($(strip $(wildcard initialize_state_transport.F90)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_state_transport.o
else
  FORTRANOBJS += initialize_state_transport.o
endif

fortran: BEGIN source $(FORTRANOBJS) chkopts
	@echo ""
	@echo "Linking executable."
	@echo "-------------------"
	-${FLINKER} -o runSimulation $(FORTRANOBJS) $(MAIN) $(LBM_LIB_DIR)/$(LBM_LIB_NAME) ${PETSC_DM_LIB} 

BEGIN:
	@echo "Building lbm simulation in $(LBM_DIR) with $(COMP)" 
	@echo ""

allclean: 
	@cd $(LBM_SRC_DIR); make allclean
	@cd $(LBM_PROBLEM_DIR); make allclean
	rm -f ./*.o
	rm -f runSimulation

source: 
	@cd $(LBM_SRC_DIR); make fortran

test: fortran
	@echo ""
	@echo "Running simulation"
	@echo "------------------"
	mpiexec -n 8 ./runSimulation input_data 
	@echo ""
	@echo "Comparing simulation solution to reference solution"
	@echo "---------------------------------------------------"
	python $(LBM_DIR)/src/testing/check_solution.py input_data

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

