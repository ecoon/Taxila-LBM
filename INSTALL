========================================================================
    author:              Ethan T. Coon
    filename:            INSTALL
    version:
    created:             12 July 2011
      on:                12:16:19 MDT
    last modified:   08 August 2011
      at:            15:34:23 MDT
    URL:                 http://www.ldeo.columbia.edu/~ecoon/
    email:               ecoon _at_ lanl.gov

========================================================================

INSTALL instructions for Taxila LBM

1. Install some flavor of mpi and PETSc.  See
http://www.mcs.anl.gov/petsc/petsc-as/developers/index.html

The development version of PETSc is required.  Shared libraries are
recommended for use with petsc4py, the preferred mode of loading
output for plotting Taxila LBM output in matplotlib/python.  Note that
PETSc can install MPICH2 for you automatically using the
--download-mpich option.

Recommended (quick start) configure method for PETSc:

./configure --with-shared-libararies=1 --download-mpich=1 --with-debugging=1


2. Set the environmental variable LBM_DIR to directory containing the
Taxila installation (likely the directory containing this file!)

3. Build the base objects by running ``make'' in the src/lbm directory.

4. Run a test by by running ``make test'' in one of the test
directories in $LBM_DIR/tests
