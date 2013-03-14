ifndef LBM_DIR
  LBM_DIR:=$(shell pwd)/..
endif

ifndef COMP
  HOST:= $(shell hostname | cut -d"-" -f1)
  ifeq ($(HOST),jaguar)
    COMP:=ftn
  else
    COMP:= $(shell basename $(shell mpif90 -show | cut -d" " -f 1))
  endif
endif

LBM_INCLUDE_DIR= $(LBM_DIR)/include
FFLAGS = -I$(LBM_INCLUDE_DIR) 

geometry: CreateGeom3DPetsc.o CreateGeom2DPetsc.o
	-${FLINKER} -o createGeom3DPetsc CreateGeom3DPetsc.o  ${PETSC_DM_LIB}
	-${FLINKER} -o createGeom2DPetsc CreateGeom2DPetsc.o  ${PETSC_DM_LIB}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
