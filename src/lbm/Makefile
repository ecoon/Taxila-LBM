ifndef LBM_DIR
  LBM_DIR:=$(shell pwd)/..
endif

COMP:= $(shell basename $(shell mpif90 -show | cut -d" " -f 1))

ifndef MAKE
  MAKE=make
endif

LBM_SRC_DIR= $(LBM_DIR)/src
LBM_INCLUDE_DIR= $(LBM_DIR)/include
LBM_LIB_DIR= $(LBM_DIR)/lib
LBM_MOD_DIR= $(LBM_LIB_DIR)/modules
LBM_LIB_NAME= lbm.a

FFLAGS += -I$(LBM_INCLUDE_DIR) 

# module location
ifeq ($(COMP), gfortran)
  FFLAGS += -J$(LBM_MOD_DIR) -I$(LBM_MOD_DIR)
else ifeq ($(COMP), pgf90)
  FFLAGS += -module $(LBM_MOD_DIR) -I$(LBM_MOD_DIR)
else ifeq ($(COMP), pgf95)
  FFLAGS += -module $(LBM_MOD_DIR) -I$(LBM_MOD_DIR)
endif

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

ifndef AR
  AR = ar
endif

ifndef ARFLAGS
  ARFLAGS = csrv
endif

MAIN = main.o 

MODULES = timing.o \
	string.o \
	lbm_logging.o \
	lbm_bag_types.o \
	lbm_options.o \
	lbm_io.o \
	lbm_info.o \
	lbm_grid.o \
	lbm_discretization_type.o \
	lbm_distribution_function_type.o \
	lbm_relaxation.o \
	lbm_phase.o \
	lbm_species.o \
	lbm_mineral.o \
	lbm_discretization_d3q19.o \
	lbm_discretization_d2q9.o \
	lbm_discretization.o \
	lbm_distribution_function.o \
	lbm_walls.o \
	lbm_forcing.o \
        lbm_bc.o \
	lbm_flow.o \
	lbm_transport.o \
	lbm.o 

fortran: library BEGINEX $(MAIN)

BEGINEX:
	@echo "building lbm $(MAIN)"
	@echo " with compiler family $(COMP)"

library: BEGINLIB $(foreach obj, $(MODULES), $(LBM_LIB_DIR)/$(LBM_LIB_NAME)($(obj)))

BEGINLIB:
	@echo "building lbm objects and library $(LBM_LIB_NAME) ... "
	@echo " with compiler family $(COMP)"

allclean: 
	rm -f ./*.mod
	rm -f ./*.o
	rm -f $(LBM_LIB_DIR)/$(LBM_LIB_NAME)
	rm -f $(LBM_MOD_DIR)/*.mod

tags:
	rm -f TAGS
	etags *.F90

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
