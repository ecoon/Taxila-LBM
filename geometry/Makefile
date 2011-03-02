COMP = ${PETSC_DIR}/${PETSC_ARCH}/bin/mpif90
FFLAGS = -I../src

petsc: CreateGeom3DPetsc.o
	-${FLINKER} -o createGeom3DPetsc CreateGeom3DPetsc.o  ${PETSC_DM_LIB}


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
