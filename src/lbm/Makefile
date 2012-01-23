ifndef LBM_DIR
  LBM_DIR:=$(shell pwd)/../..
endif

HOST:= $(shell hostname | cut -d"-" -f1)
ifeq ($(HOST),jaguar)
  COMP:=ftn
else
  COMP:= $(shell basename $(shell mpif90 -show | cut -d" " -f 1))
endif

LBM_SRC_DIR= $(LBM_DIR)/src/lbm
LBM_UTIL_DIR= $(LBM_DIR)/src/utils
LBM_PROBLEM_DIR= $(LBM_DIR)/src/problem_specs
LBM_INCLUDE_DIR= $(LBM_DIR)/include
LBM_LIB_DIR= $(LBM_DIR)/lib
LBM_MOD_DIR= $(LBM_LIB_DIR)/modules
LBM_LIB_NAME= lbm.a

FFLAGS += -I$(LBM_INCLUDE_DIR) 

# module location
ifeq ($(COMP),gfortran)
  FFLAGS += -J$(LBM_MOD_DIR) -I$(LBM_MOD_DIR)
else ifeq ($(COMP),pgf90)
  FFLAGS += -module $(LBM_MOD_DIR) -I$(LBM_MOD_DIR)
else ifeq ($(COMP),pgf95)
  FFLAGS += -module $(LBM_MOD_DIR) -I$(LBM_MOD_DIR)
endif

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

ifndef AR
  AR = ar
endif

ifndef ARFLAGS
  ARFLAGS = csrv
endif

MAIN = main.o 

UTILS = timing.o \
	string.o

MODULES = lbm_logging.o \
	lbm_error.o \
	lbm_options.o \
	lbm_io.o \
	lbm_info.o \
	lbm_grid.o \
	lbm_discretization_type.o \
	lbm_distribution_function_type.o \
	lbm_eos.o \
	lbm_relaxation.o \
	lbm_component.o \
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

fortran: BEGINLIB library $(MAIN)

library: BEGINLIB utils $(foreach obj, $(MODULES), $(LBM_LIB_DIR)/$(LBM_LIB_NAME)($(obj)))
BEGINLIB:
	@echo "building lbm objects and library $(LBM_LIB_NAME) ... "
	@echo " with compiler family $(COMP)"
	@echo ""
allclean: 
	@cd $(LBM_UTIL_DIR); make allclean
	@cd $(LBM_PROBLEM_DIR); make allclean
	rm -f ./*.mod
	rm -f ./*.o
	rm -f $(LBM_LIB_DIR)/$(LBM_LIB_NAME)
	rm -f $(LBM_MOD_DIR)/*.mod
utils:
	@cd $(LBM_UTIL_DIR); make fortran

tags:
	rm -f TAGS
	etags *.F90

# these dependencies should likely be made explicit, but I'm not exactly sure what the
# best way to do this without forcing the remake of .o files that already exist within 
# the library (and have therefore been removed).
# lbm_options.o: utils
# lbm_io.o: lbm_options.o
# lbm_info.o: lbm_bag_types.o lbm_options.o
# lbm_grid.o: lbm_info.o lbm_options.o
# lbm_distribution_function_type.o: lbm_bag_types.o lbm_info.o lbm_discretization_type.o 
# lbm_relaxation.o: lbm_bag_types.o lbm_options.o lbm_distribution_function_type.o
# lbm_component.o: lbm_bag_types.o lbm_options.o lbm_relaxation.o
# lbm_species.o: lbm_bag_types.o lbm_options.o lbm_relaxation.o
# lbm_mineral.o: lbm_bag_types.o lbm_options.o
# lbm_discretization_d3q19.o: lbm_discretization_type.o lbm_distribution_function_type.o \
# 	lbm_relaxation.o
# lbm_discretization_d2q9.o: lbm_discretization_type.o lbm_distribution_function_type.o \
# 	lbm_relaxation.o
# lbm_discretization.o: lbm_discretization_type.o lbm_distribution_function_type.o \
# 	lbm_relaxation.o lbm_discretization_d3q19.o lbm_discretization_d2q9.o
# lbm_distribution_function.o: lbm_info.o lbm_distribution_function_type.o \
# 	lbm_discretization_type.o
# lbm_walls.o: lbm_options.o lbm_info.o lbm_grid.o lbm_mineral.o 
# lbm_forcing.o: lbm_distribution_function_type.o lbm_component.o lbm_distribution_function.o \
# 	lbm_walls.o 
# lbm_bc.o: lbm_options.o lbm_info.o lbm_grid.o lbm_distribution_function_type.o \
# 	lbm_discretization.o
# lbm_flow.o: lbm_logging.o lbm_io.o lbm_options.o lbm_info.o lbm_grid.o \
# 	lbm_distribution_function_type.o lbm_discretization_type.o lbm_relaxation.o \
# 	lbm_component.o lbm_discretization.o lbm_distribution_function.o lbm_walls.o \
# 	lbm_forcing.o lbm_bc.o 
# lbm_transport.o: lbm_logging.o lbm_io.o lbm_options.o lbm_info.o lbm_grid.o \
# 	lbm_distribution_function_type.o lbm_discretization_type.o lbm_relaxation.o \
# 	lbm_component.o lbm_discretization.o lbm_distribution_function.o lbm_walls.o \
# 	lbm_forcing.o lbm_bc.o 
# lbm.o: utils lbm_logging.o lbm_options.o lbm_io.o lbm_info.o lbm_grid.o \
# 	lbm_distribution_function_type.o lbm_discretization_type.o lbm_relaxation.o \
# 	lbm_component.o lbm_discretization.o lbm_distribution_function.o lbm_walls.o \
# 	lbm_forcing.o lbm_bc.o lbm_flow.o lbm_transport.o
# main.o: lbm_logging.o lbm_options.o lbm_bc.o lbm.o

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
