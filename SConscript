import glob
import os

Import("*")

# Note that this script is executed from within the build/<something>/ folder
curdir = os.getcwd()

# Get our top-level directory
toplevel_dir = os.path.basename(os.path.dirname(os.path.dirname(curdir)))

# Look for .cpp files within the src folder
os.chdir('../..') # This is so .o files are built in `toplevel_dir'/build/
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
# line.
testfiles = set()
if single_test_suite:
  if single_test_suite_dir == toplevel_dir:
    testfiles.add(single_test_suite)
    # Remove any old test output file to force a re-run
    try:
      os.remove(single_test_suite[:-4] + '.log')
    except OSError:
      pass
else:
  packfiles = []
  if all_tests:
    for packfile in glob.glob('../../test/*TestPack.txt'):
      try:
        packfiles.append(file(packfile, 'r'))
      except IOError:
        pass
  else:
    for testpack in build.TestPacks():
      try:
        packfile = '../../test/'+testpack+'TestPack.txt'
        packfiles.append(file(packfile, 'r'))
      except IOError:
        pass
  for packfile in packfiles:
    try:
      for testfile in map(lambda s: s.strip(), packfile.readlines()):
        # Ignore blank lines and repeated tests.
        if testfile and not testfile in testfiles:
          testfiles.add(testfile)
      packfile.close()
    except IOError:
      pass


# Look for source files that tests depend on in test/.
testsource = []
for file in os.listdir('../../test'):
  if file[-4:] == '.cpp':
    testsource.append('test/' + file)

#print files, testfiles, testsource


petsc_libs = ['petscts', 'petscsnes', 'petscksp', 'petscdm', 
              'petscmat', 'petscvec', 'petsc']
chaste_libs = ['io', 'ode', 'pde', 'coupled', 'linalg', 'mesh', 'models', 'global']

all_libs = chaste_libs + petsc_libs + blas_libs + other_libs + ['test'+toplevel_dir]

opt = Environment(ENV = {'PATH': os.environ['PATH'],
                         'USER': os.environ['USER'],
                         'CHASTE_TEST_OUTPUT': os.environ.get('CHASTE_TEST_OUTPUT',
                                                              '/tmp/'+os.environ['USER']+'/testoutput/')})
opt.Append(CCFLAGS = petsc_incs+extra_flags)
opt.Append(LINKFLAGS = link_flags)
opt.Append(BOPT = 'g_c++')
opt.Replace(CXX = mpicxx)
opt.Replace(AR = ar)
opt.Replace(CPPPATH = cpppath)

test = Builder(action = 'cxxtest/cxxtestgen.py --error-printer -o $TARGET $SOURCES')
runtests = Builder(action = 'python/TestRunner.py $SOURCE $TARGET ' +
                   build_type + ' ' + build.GetTestReportDir() + 
                   ' ' + run_time_flags)

opt['BUILDERS']['Test'] = test
opt['BUILDERS']['RunTests'] = runtests

opt['ENV']['LD_LIBRARY_PATH'] = petsc_base+'lib/libg_c++/linux-gnu/'
opt.Library(toplevel_dir, files)
opt.Install('../../../lib', 'lib'+toplevel_dir+'.a')
opt.Library('test'+toplevel_dir, testsource)

for testfile in testfiles:
  prefix = testfile[:-4]
  opt.Test(prefix+'Runner.cpp', 'test/' + testfile) 
  opt.Program(testfile[:-4]+'Runner', [prefix+'Runner.cpp'],
              LIBS = all_libs,
              LIBPATH = ['../../../lib', '.', petsc_libpath, blas_libpath] + other_libpaths)
  if not compile_only:
    opt.RunTests(prefix+'.log', prefix+'Runner')

