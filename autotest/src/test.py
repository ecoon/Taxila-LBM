### ====================================================================
###  Python-file
###     author:          Ethan T. Coon
###     filename:        test.py
###     version:         
###     created:         28 January 2011
###       on:            09:52:47 MST
###     last modified:   17 March 2011
###       at:            10:47:19 MDT
###     URL:             http://www.ldeo.columbia.edu/~ecoon/
###     email:           ecoon _at_ lanl.gov
###  
### ====================================================================


import unittest
import subprocess
import os
import sys
import solution_reader

_startdir = os.getcwd()
_np = '1'


class LBMTest(unittest.TestCase):
    def __init__( self, directory, postfix, doc ):
        self._directory = os.path.join(_startdir, '..', directory)
        self._executable = os.path.join(self._directory, 'runSimulation')
        self._postfix = postfix
        self._doc = doc + ' in ../%s, input_data_%s'%(directory, postfix)
        super(LBMTest, self).__init__('test')

    def setUp( self ):
        os.chdir(self._directory)
        self._stdout = file(os.path.join(self._directory, 'test_solution_'+self._postfix, 'output.log'), 'w')
        self._stderr = file(os.path.join(self._directory, 'test_solution_'+self._postfix, 'output.err'), 'w')

    def test( self ):
        # make the executable
        self._stdout.write('make\n\n')
        self._stdout.flush()
        make = subprocess.Popen(['make', 'allclean'], stdout=self._stdout, stderr=self._stderr)
        make.wait()
        self.failUnlessEqual(make.returncode, 0, msg='clean failed')
        make = subprocess.Popen(['make'], stdout=self._stdout, stderr=self._stderr)
        make.wait()
        self.failUnlessEqual(make.returncode, 0, msg='compilation failed')

        # run the executable
        self._stdout.write('\n\n-------------------------------------------\n')
        self._stdout.write('mpiexec -n %s ./runSimulation input_data_%s\n\n'%(_np, self._postfix))
        self._stdout.flush()
        run = subprocess.Popen(['mpiexec', '-n', _np, './runSimulation', 'input_data_'+self._postfix], stdout=self._stdout, stderr=self._stderr)
        run.wait()
        self.failUnlessEqual(run.returncode, 0, msg='execution failed')

        # test the result
        self._stdout.write('\n\n-------------------------------------------\n')
        self._stdout.write('python ../src/check_solution.py input_data_%s\n\n'%(self._postfix))
        self._stdout.flush()
        test = subprocess.Popen(['python', '../src/check_solution.py', 'input_data_'+self._postfix], stdout=self._stdout, stderr=self._stderr)
        test.wait()
        self.failUnlessEqual(test.returncode, 0, msg='regression failed')

    def tearDown( self ):
        self._stdout.close()        
        self._stderr.close()
        os.chdir(_startdir)
        
    def shortDescription( self ):
        return self._doc


class LBMRestartTest(LBMTest):
    def test( self ):
        # make the executable
        self._stdout.write('make\n\n')
        self._stdout.flush()
        make = subprocess.Popen(['make'], stdout=self._stdout, stderr=self._stderr)
        make.wait()
        self.failUnlessEqual(make.returncode, 0, msg='compilation failed')
        
        # run the executable
        self._stdout.write('\n\n-------------------------------------------\n')
        self._stdout.write('mpiexec -n %s ./runSimulation input_data_%s\n\n'%(_np, self._postfix))
        self._stdout.flush()
        run = subprocess.Popen(['mpiexec', '-n', _np, './runSimulation', 'input_data_'+self._postfix], stdout=self._stdout, stderr=self._stderr)
        run.wait()
        self.failUnlessEqual(run.returncode, 0, msg='execution failed')

        # run the restart
        self._stdout.write('\n\n-------------------------------------------\n')
        self._stdout.write('mpiexec -n %s ./runSimulation input_data_%s_restart\n\n'%(_np, self._postfix))
        self._stdout.flush()
        run = subprocess.Popen(['mpiexec', '-n', _np, './runSimulation', 'input_data_'+self._postfix+'_restart'], stdout=self._stdout, stderr=self._stderr)
        run.wait()
        self.failUnlessEqual(run.returncode, 0, msg='restart execution failed')

        # test the result
        self._stdout.write('\n\n-------------------------------------------\n')
        self._stdout.write('python ../src/check_solution.py --nsteps=2 input_data_%s\n\n'%(self._postfix))
        self._stdout.flush()
        test = subprocess.Popen(['python', '../src/check_solution.py', '--nsteps=2', 'input_data_'+self._postfix], stdout=self._stdout, stderr=self._stderr)
        test.wait()
        self.failUnlessEqual(test.returncode, 0, msg='regression failed')
        

class LBMRotatedTest(LBMTest):
    def test( self ):
        # make the executable
        self._stdout.write('make\n\n')
        self._stdout.flush()
        make = subprocess.Popen(['make', 'allclean'], stdout=self._stdout, stderr=self._stderr)
        make.wait()
        self.failUnlessEqual(make.returncode, 0, msg='clean failed')
        make = subprocess.Popen(['make'], stdout=self._stdout, stderr=self._stderr)
        make.wait()
        self.failUnlessEqual(make.returncode, 0, msg='compilation failed')

        # run the executable
        self._stdout.write('\n\n-------------------------------------------\n')
        self._stdout.write('mpiexec -n %s ./runSimulation input_data_%s\n\n'%(_np, self._postfix))
        self._stdout.flush()
        run = subprocess.Popen(['mpiexec', '-n', _np, './runSimulation', 'input_data_'+self._postfix], stdout=self._stdout, stderr=self._stderr)
        run.wait()
        self.failUnlessEqual(run.returncode, 0, msg='execution failed')

        # test the result
        self._stdout.write('\n\n-------------------------------------------\n')
        self._stdout.write('python ../src/check_solution.py --rotate input_data_%s\n\n'%(self._postfix))
        self._stdout.flush()
        test = subprocess.Popen(['python', '../src/check_solution.py', '--rotate', 'input_data_'+self._postfix], stdout=self._stdout, stderr=self._stderr)
        test.wait()
        self.failUnlessEqual(test.returncode, 0, msg='regression failed')

class LBMDoubleRotatedTest(LBMTest):
    def test( self ):
        # make the executable
        self._stdout.write('make\n\n')
        self._stdout.flush()
        make = subprocess.Popen(['make', 'allclean'], stdout=self._stdout, stderr=self._stderr)
        make.wait()
        self.failUnlessEqual(make.returncode, 0, msg='clean failed')
        make = subprocess.Popen(['make'], stdout=self._stdout, stderr=self._stderr)
        make.wait()
        self.failUnlessEqual(make.returncode, 0, msg='compilation failed')

        # run the executable
        self._stdout.write('\n\n-------------------------------------------\n')
        self._stdout.write('mpiexec -n %s ./runSimulation input_data_%s\n\n'%(_np, self._postfix))
        self._stdout.flush()
        run = subprocess.Popen(['mpiexec', '-n', _np, './runSimulation', 'input_data_'+self._postfix], stdout=self._stdout, stderr=self._stderr)
        run.wait()
        self.failUnlessEqual(run.returncode, 0, msg='execution failed')

        # test the result
        self._stdout.write('\n\n-------------------------------------------\n')
        self._stdout.write('python ../src/check_solution.py --double-rotate input_data_%s\n\n'%(self._postfix))
        self._stdout.flush()
        test = subprocess.Popen(['python', '../src/check_solution.py', '--double-rotate', 'input_data_'+self._postfix], stdout=self._stdout, stderr=self._stderr)
        test.wait()
        self.failUnlessEqual(test.returncode, 0, msg='regression failed')


def suite():
    thesuite = unittest.TestSuite()

    thesuite.addTest(LBMTest('bubble', '3D', '3D bubble test'))
    thesuite.addTest(LBMTest('bubble', '3D_singledensity', '3D bubble test with mm=1'))
    thesuite.addTest(LBMTest('flushing', '3D_fluxpressure', '3D flux and pressure BC test'))
    thesuite.addTest(LBMTest('flushing', '3D_pressure', '3D pressure BC test'))
    thesuite.addTest(LBMTest('flushing', '3D_poiseuille', '3D, nonconstant bc test'))
    thesuite.addTest(LBMRestartTest('flushing', '3D_flux', '3D flux, restart test'))
    thesuite.addTest(LBMRotatedTest('flushing', '3D_pressure_rot', '3D pressure BC test, rotated to x'))
    thesuite.addTest(LBMDoubleRotatedTest('flushing', '3D_pressure_rotrot', '3D pressure BC test, rotated to y'))
    thesuite.addTest(LBMRotatedTest('flushing', '3D_flux_rot', '3D flux BC test, rotated to x'))
    thesuite.addTest(LBMDoubleRotatedTest('flushing', '3D_flux_rotrot', '3D flux BC test, rotated to y'))
    return thesuite


if __name__ == '__main__':
    import optparse
    opts = optparse.OptionParser()
    opts.add_option('-n', 
                    dest='np', default='1',
                    help='number of processors on which to run tests')
    options, others = opts.parse_args(sys.argv[1:])
    globals()['_np'] = options.np
    testsuite = suite()
    unittest.TextTestRunner().run(testsuite)
    
