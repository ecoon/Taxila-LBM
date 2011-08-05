ifndef LBM_DIR
  LBM_DIR:=$(shell pwd)/../..
endif

COMP:= $(shell basename $(shell mpif90 -show | cut -d" " -f 1))

ifndef MAKE
  MAKE=make
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
  ifeq ($(COMP), gfortran)
    FFLAGS += -fbounds-check 
  else ifeq ($(COMP), pgf90)
    FFLAGS += -Mbounds
  else ifeq ($(COMP), pgf95)
    FFLAGS += -Mbounds
  endif
else
  FFLAGS += -O3
endif

MAIN = $(LBM_SRC_DIR)/main.o 


ifeq ($(strip $(wildcard initialize_state.o)),)
  FORTRANOBJS = $(LBM_PROBLEM_DIR)/default_initialize_state.o
else
  FORTRANOBJS = initialize_state.o
endif
ifeq ($(strip $(wildcard initialize_walls.o)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_walls.o
else
  FORTRANOBJS += initialize_walls.o
endif
ifeq ($(strip $(wildcard initialize_bcs.o)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_bcs.o
else
  FORTRANOBJS += initialize_bcs.o
endif
ifeq ($(strip $(wildcard initialize_bcs_transport.o)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_bcs_transport.o
else
  FORTRANOBJS += initialize_bcs_transport.o
endif
ifeq ($(strip $(wildcard initialize_state_transport.o)),)
  FORTRANOBJS += $(LBM_PROBLEM_DIR)/default_initialize_state_transport.o
else
  FORTRANOBJS += initialize_state_transport.o
endif

fortran: BEGIN $(FORTRANOBJS) chkopts
	echo $(LBM_LIB_NAME)
	-${FLINKER} -o runSimulation $(FORTRANOBJS) $(MAIN) $(LBM_LIB_DIR)/$(LBM_LIB_NAME) ${PETSC_DM_LIB} 

BEGIN:
	@echo "building lbm simulation in $(LBM_DIR) with $(COMP)" 
	@echo "and objects: ${FORTRANOBJS}"


allclean: 
	cd $(LBM_SRC_DIR); $(MAKE) allclean
	cd $(LBM_PROBLEM_DIR); $(MAKE) allclean
	rm -f ./*.o
	rm -f runSimulation runLBMSimulation

source: 
	cd $(LBM_SRC_DIR); $(MAKE) fortran

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

