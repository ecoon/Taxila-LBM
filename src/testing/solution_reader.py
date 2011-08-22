### ====================================================================
###  Python-file
###     author:          Ethan T. Coon
###     filename:        solution_reader.py
###     version:         
###     created:         25 January 2011
###       on:            10:53:33 MST
###     last modified:   30 March 2011
###       at:            12:30:14 MDT
###     URL:             http://www.ldeo.columbia.edu/~ecoon/
###     email:           ecoon _at_ lanl.gov
###  
### ====================================================================

import os
import numpy as np
import pyvtk

def generate_commandline(filename):
    infile = open(filename, 'r')
    args = [line[:-1] for line in infile if not line.strip().startswith('#')]
    infile.close()

    if any('use_old_options_style' in line for line in args):
        raise NotImplementedError("Old style options can't be used here... generate a PETSc-style options file and try again")
    else:
        import petsc4py
        petsc4py.init(' '.join(args))
        from petsc4py import PETSc
    return

class SolutionReader(object):
    def __init__( self, prefix ):
        from petsc4py import PETSc
        self._prefix = prefix
        self._opts = PETSc.Options(prefix)
        self._discretization = self._opts.getString('discretization', default='d3q19')
        if  self._discretization == 'd2q9':
            self._size = (self._opts.getInt('NX'), self._opts.getInt('NY'))
            self._dim = 2
            print 'got d2q9'
        elif self._discretization == 'd3q19':
            self._size = (self._opts.getInt('NX'), self._opts.getInt('NY'),
                          self._opts.getInt('NZ'))
            self._dim = 3
            print 'got d3q19'
        self._s = self._opts.getInt('ncomponents',default=2)
        lsize = list(self._size)
        lsize.reverse()
        self._size_r = tuple(lsize)
        self._file_prefix = self._opts.getString('output_file_prefix', default='test_solution/')
        self._vecs = dict()
        self._scalefactor = self._opts.getReal('velocity_scalefactor', default=1.0)

    def loadVec( self, name, ndofs=1 ):
        from petsc4py import PETSc
        length = np.array(self._size).prod()
        print 'loading', self._file_prefix+name
        vec = PETSc.Vec().createSeq(length*ndofs)
        vec.setBlockSize(ndofs)
        viewer = PETSc.Viewer().createBinary(self._file_prefix+name, PETSc.Viewer.Mode.R)
        vec.load(viewer)
        if self._dim == 3:
            npvec = vec[...].reshape(self._size_r+(ndofs,)).transpose((2,1,0,3))[:]
        elif self._dim == 2:
            npvec = vec[...].reshape(self._size_r+(ndofs,)).transpose((1,0,2))[:]
        viewer.destroy()
        vec.destroy()
        del viewer
        del vec
        return npvec

    def loadVecToVTK( self, name, ndofs=1 ):
        from petsc4py import PETSc
        length = np.array(self._size).prod()
        print 'loading', self._file_prefix+name
        vec = PETSc.Vec().createSeq(length*ndofs)
        viewer = PETSc.Viewer().createBinary(self._file_prefix+name, PETSc.Viewer.Mode.R)
        vec.load(viewer)

        if self._size[2] == 1:
            length = 2*length
            npvec = vec[...].reshape((self._size_r+(ndofs,)))[:]
            npvec = np.repeat(npvec,2,axis=0)
            npvec = npvec.reshape((length,ndofs))
        else:
            npvec = vec[...].reshape((length,ndofs))[:]

        # scale
        if name.startswith('u'):
            npvec = npvec*self._scalefactor
            print np.where(np.abs(npvec)/np.abs(npvec).mean() > 1e2)[0]
            npvec = np.where(np.abs(npvec)/np.abs(npvec).mean() > 1e2, 0., npvec)

        if ndofs == 1:
            data = pyvtk.Scalars(npvec, self._prefix.strip('_')+' '+name[:-7], 'default')
        else:
            data = pyvtk.Vectors([tuple(npvec[i,:]) for i in range(length)], self._prefix.strip('_')+' '+name[:-7])
        viewer.destroy()
        vec.destroy()
        del viewer
        del vec
        return data

    def solnToVTK( self ):
        done = False
        lcv = 0
        
        coords = self.loadVec('coords.dat', 3)
        dims = list(coords.shape[:-1])
        try:
            dx = coords[1,1,1]-coords[0,0,0]
        except IndexError:
            try:
                dx = coords[0,1,1] - coords[0,0,0]
            except IndexError:
                try:
                    dx = coords[1,0,1] - coords[0,0,0]
                except IndexError:
                    dx = coords[1,1,0] - coords[0,0,0]
            
        dx = np.where(dx==0., 0.1, dx)
        if dims[2] == 1:
            dims[2] = 2
        dims = tuple(dims)
        print dims
        dx = tuple(dx)
        vtkgrid = pyvtk.StructuredPoints(dims, coords[0,0,0], dx)

        while not done:
            try:
                if not os.path.exists(self._file_prefix+'prs%03d.dat'%lcv):
                    raise IOError('Nonexistent file')
            except IOError:
                done = True
                print 'Read %d timesteps'%lcv
            else:
                prs_data = self.loadVecToVTK('prs%03d.dat'%lcv, 1)
                vel_data = self.loadVecToVTK('u%03d.dat'%lcv, 3)
                wall_data = self.loadVecToVTK('walls%03d.dat'%lcv, 1)

                pointdata = pyvtk.PointData(prs_data, vel_data, wall_data)
                data = pyvtk.VtkData(vtkgrid, self._prefix.strip('_')+' step %d'%lcv, pointdata)
                data.tofile(self._file_prefix+'soln_%03d.vtk'%lcv)
                lcv += 1
        return

