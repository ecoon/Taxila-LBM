Taxila LBM
========================================================================

NOTE: This code is no longer supported, and has reached the end of its software lifecycle.  Any are welcome to use any or all of this code according to the license below, but it is provided as is.


About
-------

Taxila LBM is a parallel implementation of the Lattice Boltzmann Method for simulation of flow in porous and geometrically complex media.

Taxila LBM uses the LGPL license as seen in LICENSE

The implementation:

* solves both single and multi-phase systems, with or without passive transport of multiple species.
* is capable of solving D2Q9, D3Q19, and other mesh dependencies, on 2D or 3D grids. It is easily extended to other models of connectivity.
* uses the Shan and Chen Lattice Boltzmann Method (ref:).
* handles multi-phase systems with different phase viscosities and/or molecular masses.
* includes the ability to use (or not use) higher order derivatives (ref:) or multiple relaxation times (ref:) to improve stability at large viscosity ratios.
* enables multiple mineral/wall materials, allowing for different wettabilities and contact angles on each mineral.
* supports arbitrary, heterogeneous boundary conditions. 

* is massively parallel, showing excellent strong and weak scaling to X cores (update with Jaguar numbers!)
* is modular, "object-oriented" Fortran-90, enabling extension for new features.
* leverages  PETSc, the Portable, Extensible Toolkit for Scientific Computation for data structures, communication, and (parallel) I/O.
* may be called from  PFloTran, a massively parallel solver for reactive transport in porous media. This allows for micro-scale solution of reactive transport, using Taxila LBM for flow. 

* currently only allows structured, regular meshes.
* currently is not implemented on GPUs (coming soon?)
* currently does not handle curved geometries or other non-bounceback style interior boundary conditions (but we're accepting contributions!) 

Authors
-------
Ethan Coon
ecoon _at_ ornl.gov

Mark Porter

Qinjun Kang
qkang _at_ lanl.gov


Quick start links
------------------

[Getting Taxila LBM](http://github.com/ecoon/taxila-lbm)

[Installing Taxila LBM](INSTALL.md)

Testing installation:

    $>  cd tests/{test_name}/
    $>  make test

Using Taxila LBM: 
  To use, see [the doc page](docs/README) on how to make the manual and quickstart guide.
