import glob
import os

Import("*")

# Note that this script is executed from within the build/ folder
curdir = os.getcwd()

# Get our top-level directory
toplevel_dir = os.path.basename(os.path.dirname(curdir))

# Look for .cpp files within the src folder
os.chdir('..') # This is so .o files are built in `toplevel_dir'/build/
files = []
for dirpath, dirnames, filenames in os.walk('src'):
  for filename in filenames:
    if filename[-4:] == '.cpp':
      files.append(os.path.join(dirpath, filename))
os.chdir(curdir)

# Look for files containing a test suite
# A list of test suites to run will be found in a test/<name>TestPack.txt
# file, one per line.
# Alternatively, a single test suite may have been specified on the command
# line
testfiles = []
if single_test_suite:
  if single_test_suite_dir == toplevel_dir:
    testfiles = [single_test_suite]
    # Remove any old test output file to force a re-run
    try:
      os.remove(single_test_suite[:-4] + '.log')
    except OSError:
      print "Hmm"
else:
  for testpack in build.TestPacks():
    try:
      packfile = file('../test/'+testpack+'TestPack.txt', 'r')
      for testfile in map(lambda s: s.strip(), packfile.readlines()):
        # Ignore blank lines and repeated tests.
        if testfile and not testfile in testfiles:
          testfiles.append(testfile)
      packfile.close()
    except IOError:
      pass


# Look for source files that tests depend on in test/.
_testsource = os.listdir('../test')
testsource = []
for file in _testsource:
  if file[-4:] == '.cpp':
    testsource.append('test/' + file)
del _testsource

#print files, testfiles, testsource


petsc_libs = ['petscts', 'petscsnes', 'petscksp', 'petscdm', 
              'petscmat', 'petscvec', 'petsc']
chaste_libs = ['global', 'io', 'ode', 'pde', 'coupled', 'maths', 'mesh']

all_libs = petsc_libs + chaste_libs + ['test'+toplevel_dir]

opt = Environment(ENV = {'PATH' : os.environ['PATH']})
opt.Append(CCFLAGS = petsc_incs+extra_flags)
opt.Append(LINKFLAGS = link_flags)
opt.Append(BOPT = 'g_c++')
opt.Replace(CXX = mpicxx)
opt.Replace(AR = ar)
opt.Replace(CPPPATH = cpppath)

test = Builder(action = 'cxxtest/cxxtestgen.pl --error-printer -o $TARGET $SOURCES')
runtests = Builder(action = 'python/TestRunner.py $SOURCE $TARGET ' +
                   build_type + ' ' + build.GetTestReportDir())
#runtests = Builder(action = './$SOURCE | tee $TARGET')
#runparalleltests = Builder(action = mpirun + ' -np 2 ./$SOURCE | tee $TARGET')

opt['BUILDERS']['Test'] = test
opt['BUILDERS']['RunTests'] = runtests
#opt['BUILDERS']['RunParallelTests'] = runparalleltests

opt['ENV']['LD_LIBRARY_PATH'] = petsc_base+'lib/libg_c++/linux-gnu/'
opt.Library(toplevel_dir, files)
opt.Install('../../lib', 'lib'+toplevel_dir+'.a')
opt.Library('test'+toplevel_dir, testsource)

for testfile in testfiles:
  prefix = testfile[:-4]
  opt.Test(prefix+'Runner.cpp', 'test/' + testfile) 
  opt.Program(testfile[:-4]+'Runner', [prefix+'Runner.cpp'],
              LIBS = all_libs,
              LIBPATH = ['../../lib', '.', petsc_libpath])
  opt.RunTests(prefix+'.log', prefix+'Runner')


# Parallel tests
#if not ARGUMENTS.get('no_parallel', 0):
#  opt.Test('parallelrunner.cpp', paralleltestfiles)
#  opt.Program('paralleltestrunner', ['parallelrunner.cpp'],
#              LIBS=all_libs,
#              LIBPATH=['../../datawriters/build', '.', petsc_libpath])
#  opt.RunParallelTests('parbuild.log', 'paralleltestrunner')

