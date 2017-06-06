Installation instructions for Taxila LBM
==================================================

Note that Taxila uses Make, and has no method for out-of-source
builds.  If you use multiple configurations (most likely a Debug
version and a Release version, as Release builds are much faster than
Debug builds but result in little help if something goes wrong),
simply download two copies of the source.

Step 1: Install PETSc and MPI
--------------------------------------------------

Install some flavor of mpi and PETSc.  See [PETSc's
website](http://www.mcs.anl.gov/petsc/) for downloads and installation
information.

Taxila attempts to stay current with the current release version of
PETSc.  Sometimes this doesn't happen too quickly, and may use an
older release.  The current version of PETSc that Taxila builds
against is:

**PETSc version 3.6**

Shared libraries are recommended for use with petsc4py, the preferred
mode of loading output for plotting Taxila LBM output in
matplotlib/python.  Note that PETSc can install MPICH2 for you
automatically using the --download-mpich option if you are building
for a small workstation without it.

Recommended (quick start) configure method for PETSc:

./configure --with-fortran-interfaces --with-shared-libraries --with-debugging --download-mpich --download-petsc4py --download-hdf5


Ensure that PETSc is installed and working, then set an environment
variable PETSC_DIR to point to your PETSc installation.


Step 2: Install Taxila LBM
--------------------------------------------------

Set the environmental variable TAXILA_DIR to the directory containing
the Taxila installation (likely the directory containing this file!)
Build the base library.

    export TAXILA_DIR=/path/to/my/taxila-lbm
    cd ${TAXILA_DIR}
    make


Step 3: Run a test
--------------------------------------------------

Run a test by by running ``make test'' in one of the test
directories in ${TAXILA_DIR}/tests
