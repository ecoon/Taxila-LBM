Geom: CreateGeom3DParallel.o
	mpif90 -o Geom CreateGeom3DParallel.o
CreateGeom3DParallel.o: CreateGeom3DParallel.f90
	mpif90 -c CreateGeom3DParallel.f90
clean:
	rm CreateGeom3DParallel.o Geom *~
