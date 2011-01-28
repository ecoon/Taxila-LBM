#!/usr/bin/env python2.6
### ====================================================================
###  Python-file
###     author:          Ethan T. Coon
###     filename:        check_solution.py
###     version:         
###     created:         28 January 2011
###       on:            10:57:45 MST
###     last modified:   28 January 2011
###       at:            12:44:58 MST
###     URL:             http://www.ldeo.columbia.edu/~ecoon/
###     email:           ecoon _at_ lanl.gov
###  
### ====================================================================

import sys, os
sys.path.append('../../pysrc')
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

    options, others = opts.parse_args(sys.argv[1:])

    solution_reader.generate_commandline(infile)
    directory = infile.split('/')[:-1]
    if directory:
        if len(directory > 0):
            directory = os.path.join(*directory)
            os.chdir(directory)

    # set up the readers
    testsoln = solution_reader.SolutionReader(options.prefix)
    truth = solution_reader.SolutionReader(options.prefix)
    truth._file_prefix = truth._file_prefix.replace('test_solution', 'solution')

    # load and compare
    ecode = 0
    eps = numpy.double(options.eps)
    print 'Testing:', testsoln._file_prefix

    for i in range(int(options.nsteps)+1):
        for vecname, ndofs in  [('prs',1), ('u',3), ('rho',testsoln._s)]:
            testdata = testsoln.loadVec(vecname+'%03d.dat'%i, ndofs)
            truedata = truth.loadVec(vecname+'%03d.dat'%i, ndofs)
            if numpy.linalg.norm((truedata - testdata).ravel()) > eps:
                ecode = 1
                print ' ',vecname+'%03d.dat'%i,'failed, norm:', numpy.linalg.norm((truedata - testdata).ravel())
            else:
                print ' ',vecname+'%03d.dat'%i,'passed'
                
    print 'ecode:', ecode
    sys.exit(ecode)

    
