### ====================================================================
###  Python-file
###     author:          Ethan T. Coon
###     filename:        petsc2tec.py
###     version:         
###     created:         14 July 2011
###       on:            14:17:38 MDT
###     last modified:   14 July 2011
###       at:            15:23:37 MDT
###     URL:             http://www.ldeo.columbia.edu/~ecoon/
###     email:           ecoon _at_ lanl.gov
###  
### ====================================================================

import numpy as np
import solution_reader2
import types

def write_tecplot(soln, outfile):
    fh = open(outfile, 'w')

    try:
        rhos = ','.join(['RHO%d'%m for m in range(soln._s)])
        rhods = ' '.join(['%g' for m in range(soln._s)])
    
        if soln._dim == 2:
            fh.write('VARIABLES=X,Y,U,V,%s,PRS\n'%rhos)
            fh.write('ZONE T="ONLY",I=%d,J=%d\n'%soln._size)
            i = 1
            j = 1

            data = '\n'.join([('%d %d %g %g '+rhods+' %g')%((i+1,j+1, \
                                                          soln._vecs['u'][i,j,0],soln._vecs['u'][i,j,1],) \
                                  +tuple(soln._vecs['rho'][i,j,m] for m in range(soln._s)) \
                                  +(soln._vecs['prs'][i,j,0],)) for j in range(soln._size[1]) for i in range(soln._size[0])])
            fh.write(data+'\n')

        elif soln._dim == 3:
            fh.write('VARIABLES=X,Y,Z,U,V,W,%s,PRS\n'%rhos)
            fh.write('ZONE T="ONLY",I=%d,J=%d,K=%d\n'%soln._size)
            
            data = '\n'.join([('%d %d %d %g %g %g '+rhods+' %g')%((i+1,j+1,k+1, \
                                                                soln._vecs['u'][i,j,k,0],soln._vecs['u'][i,j,k,1],soln._vecs['u'][i,j,k,2],) \
                                  +tuple(soln._vecs['rho'][i,j,k,m] for m in range(soln._s)) \
                                  +(soln._vecs['prs'][i,j,k,0],)) for k in range(soln._size[2]) for j in range(soln._size[1]) for i in range(soln._size[0])])
            fh.write(data+'\n')
    finally:
        fh.close()
    return
        

if __name__ == '__main__':
    import sys
    import optparse

    opts = optparse.OptionParser()
    opts.add_option('-n', '--stepnum',
                    dest='stepnum', default='1',
                    help='step number of files to read')
    opts.add_option('-o', '--outfile',
                    dest='outfile', default='lbm_output.tec',
                    help='file name for output')
    opts.add_option('-p', '--prefix',
                    dest='prefix', default='',
                    help='prefix for input file')
    opts.set_usage('petsc2tec.py [OPTIONS] input_filename')
    options, others = opts.parse_args(sys.argv[1:])
    if len(others) != 1:
        opts.print_usage()
        sys.exit(1)
    else:
        infile = others[0]

    sr = solution_reader2.SolutionReader(infile, options.prefix)
    sr.loadSoln(int(options.stepnum))
    write_tecplot(sr, options.outfile)
