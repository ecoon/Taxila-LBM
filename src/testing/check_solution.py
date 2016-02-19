#!/usr/bin/env python2.6
### ====================================================================
###  Python-file
###     author:          Ethan T. Coon
###     filename:        check_solution.py
###     version:         
###     created:         28 January 2011
###       on:            10:57:45 MST
###     last modified:   08 August 2011
###       at:            15:29:50 MDT
###     URL:             http://www.ldeo.columbia.edu/~ecoon/
###     email:           ecoon _at_ lanl.gov
###  
### ====================================================================

import sys, os
import solution_reader2 as solution_reader
import optparse
import numpy



if __name__ == '__main__':
    infile = sys.argv.pop()
    
    # get options
    opts = optparse.OptionParser()
    opts.add_option('-n','--nsteps',
                    dest='nsteps', default='1',
                    help='number of timesteps to compare')
    opts.add_option('-p','--prefix',
                    dest='prefix', default='',
                    help='optional PETSc file prefix')
    opts.add_option('-e', '--eps',
                    dest='eps', default='1.e-5')
    opts.add_option('--rotate',
                    action='store_true', dest='rotate', default=False)
    opts.add_option('--double-rotate',
                    action='store_true', dest='doublerotate', default=False)

    options, others = opts.parse_args(sys.argv[1:])

    directory = infile.split('/')[:-1]
    if directory:
        if len(directory) > 0:
            directory = os.path.join(*directory)

    # set up the readers
    test = solution_reader.SolutionReader(infile, options.prefix)
    truth = solution_reader.SolutionReader(infile, options.prefix)
    truth._file_prefix = truth._file_prefix.replace('test_solution', 'reference_solution')
    if options.rotate:
        truth._size = (test._size[1], test._size[2], test._size[0])
        truth._size_r = (test._size[0], test._size[2], test._size[1])
        truth._file_prefix = truth._file_prefix.replace('_rot', '')
        print truth._size
        print truth._size_r

    if options.doublerotate:
        truth._size = (test._size[2], test._size[0], test._size[1])
        truth._size_r = (test._size[1], test._size[0], test._size[2])
        truth._file_prefix = truth._file_prefix.replace('_rotrot', '')

    # load and compare
    ecode = 0
    eps = numpy.double(options.eps)
    print 'Testing:', test._file_prefix

    if test._dim == 2:
        tests = [('fi',9*test._s)]
    elif test._dim == 3:
        tests = [('fi',19*test._s)]
    else:
        raise RuntimeError('invalid dims')
        
    for i in range(1,int(options.nsteps)+1):
        for vecname, ndofs in  tests:
            testdata = test.loadVec(vecname+'%03d.dat'%i)
            truedata = truth.loadVec(vecname+'%03d.dat'%i)
            if options.rotate:
                truedata = truedata.transpose((2,0,1,3))
                if vecname == 'u':
                    # vectorized, must cycle dimensions
                    truedata = numpy.roll(truedata, 1, axis=3)
            if options.doublerotate:
                truedata = truedata.transpose((1,2,0,3))
                if vecname == 'u':
                    # vectorized, must cycle dimensions
                    truedata = numpy.roll(truedata, 2, axis=3)

            print truedata.shape
            print testdata.shape
            if numpy.linalg.norm((truedata - testdata).ravel(),ord=numpy.inf) > eps:
                ecode = 1
                print ' ',vecname+'%03d.dat'%i,'failed, L_inf norm:', numpy.linalg.norm((truedata - testdata).ravel(),ord=numpy.inf)

                from matplotlib import pyplot as plt
                plt.figure()
                ndofs = testdata.shape[-1]
                for m in range(ndofs):
                    plt.subplot(ndofs,1,m)
                    if test._dim == 2:
                        plt.imshow((testdata-truedata)[:,:,m],origin='lower')
                    else:
                        plt.imshow((testdata-truedata)[:,:,testdata.shape[2]/2,m],origin='lower')
                    plt.colorbar()
                plt.title(vecname+'%03d'%i)
                plt.show()

            else:
                print ' ',vecname+'%03d.dat'%i,'passed'
                
    print 'ecode:', ecode
    sys.exit(ecode)

    
