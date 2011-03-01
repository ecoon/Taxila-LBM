#!/usr/bin/env python2.6
### ====================================================================
###  Python-file
###     author:          Ethan T. Coon
###     filename:        check_solution.py
###     version:         
###     created:         28 January 2011
###       on:            10:57:45 MST
###     last modified:   03 February 2011
###       at:            11:50:01 MST
###     URL:             http://www.ldeo.columbia.edu/~ecoon/
###     email:           ecoon _at_ lanl.gov
###  
### ====================================================================

import sys, os
import solution_reader
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
                    dest='eps', default='1.e-7')
    opts.add_option('--rotate',
                    action='store_true', dest='rotate', default=False)
    opts.add_option('--double-rotate',
                    action='store_true', dest='doublerotate', default=False)

    options, others = opts.parse_args(sys.argv[1:])

    solution_reader.generate_commandline(infile)
    directory = infile.split('/')[:-1]
    if directory:
        if len(directory) > 0:
            directory = os.path.join(*directory)
            os.chdir(directory)

    # set up the readers
    test = solution_reader.SolutionReader(options.prefix)
    truth = solution_reader.SolutionReader(options.prefix)
    truth._file_prefix = truth._file_prefix.replace('test_solution', 'solution')
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

    for i in range(int(options.nsteps)+1):
        for vecname, ndofs in  [('prs',1), ('u',3), ('rho',test._s)]:
            testdata = test.loadVec(vecname+'%03d.dat'%i, ndofs)
            truedata = truth.loadVec(vecname+'%03d.dat'%i, ndofs)
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
            if numpy.linalg.norm((truedata - testdata).ravel()) > eps:
                ecode = 1
                print ' ',vecname+'%03d.dat'%i,'failed, norm:', numpy.linalg.norm((truedata - testdata).ravel())

                from matplotlib import pyplot as plt
                plt.figure()
                ndofs = testdata.shape[-1]
                for m in range(ndofs):
                    plt.subplot(ndofs,1,m)
                    plt.imshow((testdata-truedata)[0,:,:,m],origin='lower')
                    plt.colorbar()
                plt.title(vecname+'%03d'%i)
                plt.show()

            else:
                print ' ',vecname+'%03d.dat'%i,'passed'
                
    print 'ecode:', ecode
    sys.exit(ecode)

    
