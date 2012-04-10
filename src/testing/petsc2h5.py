### ====================================================================
###  Python-file
###     author:          Ethan T. Coon
###     filename:        petsc2h5.py
###     version:         
###     created:         15 November 2011
###       on:            10:07:25 MST
###     last modified:   15 November 2011
###       at:            10:26:02 MST
###     URL:             http://www.ldeo.columbia.edu/~ecoon/
###     email:           ecoon _at_ lanl.gov
###  
### ====================================================================

import sys,os
sys.path.append(os.path.join(os.environ['LBM_DIR'], 'src', 'testing'))
import solution_reader3
import h5py
import numpy as np

def _create_coords(out, shape):
    coords = out.create_group('Coordinates')
    coords.create_dataset(name='X [m]', shape=(shape[0]+1,), dtype=np.float, data=np.arange(shape[0]+1))
    coords.create_dataset(name='Y [m]', shape=(shape[1]+1,), dtype=np.float, data=np.arange(shape[1]+1))
    coords.create_dataset(name='Z [m]', shape=(shape[2]+1,), dtype=np.float, data=np.arange(shape[2]+1))


def petsc2h5(directory, infile='input_data', outfile='soln.h5'):
    sr = solution_reader3.SolutionReader(infile)
    out = h5py.File(outfile, 'w')

    # make the coordinates, which are exteriors of cells if pflotran's h5 format
    if len(sr._size) < 3:
        _create_coords(out, sr._size+(1,))
    else:
        _create_coords(out, sr._size)

    for i in range(int(sr._options['npasses'])/int(sr._options['kwrite'])+1):
        if len(sr._size) < 3:
            prs = np.array(sr.loadVec('prs%03d.dat'%i)[:,:,0], dtype=np.float)
            rho_raw = np.array(sr.loadVec('rho%03d.dat'%i), dtype=np.float)
            rho = [rho[:,:,i] for i in range(rho.shape[2])]
            u = np.array(sr.loadVec('u%03d.dat'%i), dtype=np.float)
            u_x = u[:,:,0]
            u_y = u[:,:,1]
            u_z = np.zeros(u_x.shape, u_x.dtype)
            shape = sr._size+(1,)
        else:
            prs = np.array(sr.loadVec('prs%03d.dat'%i)[:,:,:,0], dtype=np.float)
            rho_raw = np.array(sr.loadVec('rho%03d.dat'%i), dtype=np.float)
            rho = [rho[:,:,:,i] for i in range(rho.shape[3])]
            u = np.array(sr.loadVec('u%03d.dat'%i),dtype=np.float)
            u_x = u[:,:,:,0]
            u_y = u[:,:,:,1]
            u_z = u[:,:,:,2]
            shape = sr._size

        groupname = 'Time:  %d.0000E+00 s'%i
        group = out.create_group(groupname)
        group.create_dataset(name='Pressure', shape=shape, data=prs)
        group.create_dataset(name='Liquid X-Velocity', shape=shape, data=u_x)
        group.create_dataset(name='Liquid Y-Velocity', shape=shape, data=u_y)
        group.create_dataset(name='Liquid Z-Velocity', shape=shape, data=u_z)
        for i,rho in enumerate(rho):
            group.create_dataset(name='Component %d Density'%i, shape=shape,
                                 data=rho)

    out.close()

if __name__ == '__main__':
    directory = sys.argv[-1]
    petsc2h5(directory)
