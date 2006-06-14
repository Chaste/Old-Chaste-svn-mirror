# Controlling scons build script for Chaste

import sys, os, glob
sys.path.append('python')
import BuildTypes

build_type = ARGUMENTS.get('build', 'default')
build = BuildTypes.GetBuildType(build_type)
build.SetRevision(ARGUMENTS.get('revision', ''))
Export('build', 'build_type')

# Specify test_summary=0 to scons to *NOT* generate a summary html page
test_summary = ARGUMENTS.get('test_summary', 1)

# Specify system_name=finarfin to scons to change default paths
system_name = ARGUMENTS.get('system_name', '')

# Specifying extra run-time flags
run_time_flags = ARGUMENTS.get('run_time_flags', '')
Export('run_time_flags')

# Specify compile_only=1 to not run any tests.
compile_only = ARGUMENTS.get('compile_only', 0)
Export('compile_only')

# To run a single test suite only, give its path (relative to the Chaste
# root) as the test_suite=<path> argument.
# This will force the test suite to be run even if the source is unchanged.
single_test_suite = ARGUMENTS.get('test_suite', '')
if single_test_suite:
  single_test_suite = single_test_suite.split(os.path.sep)
  single_test_suite_dir = single_test_suite[0]
  single_test_suite = single_test_suite[-1]
  #print single_test_suite, single_test_suite_dir
else:
  single_test_suite_dir = ''
Export('single_test_suite', 'single_test_suite_dir')


## PETSc library paths
if system_name == 'finarfin':
  # Finarfin (Debian sarge)
  petsc_base = '/home/jonc/work/dtc/courses/IB/petsc-2.2.1/'
  petsc_inc = '-I'+petsc_base+'include '
  petsc_bmake = '-I'+petsc_base+'bmake/linux-gnu '
  petsc_mpi = '-I'+petsc_base+'include/mpiuni '
  petsc_incs = petsc_inc+petsc_bmake+petsc_mpi
  
  petsc_libpath = petsc_base+'lib/libg_c++/linux-gnu/'
  blas_libs = ['f2clapack', 'f2cblas']
elif system_name == 'zuse':
  petsc_base = '/home/zuse/system/software/petsc-2.2.1/'
  petsc_inc = '-I'+petsc_base+'include '
  petsc_bmake = '-I'+petsc_base+'bmake/linux-mpich-gnu-mkl '
  #petsc_mpi = '-I'+petsc_base+'include/mpiuni '
  petsc_mpi = ' '
  boost = '-I/home/zuse/system/joe '
  petsc_incs = petsc_inc+petsc_bmake+petsc_mpi+boost
  
  petsc_libpath = petsc_base+'lib/libg_c++/linux-mpich-gnu-mkl/'
  blas_libpath = '/opt/intel/mkl/8.0/lib/em64t/'
  blas_libs = ['']

else:
  # Default
  petsc_base = '../../../petsc-2.3.1-p13/'
  petsc_inc = '-I'+petsc_base+'include '
  petsc_bmake = '-I'+petsc_base+'bmake/linux-gnu '
  # petsc_mpi = '-I'+petsc_base+'include/mpiuni '
  petsc_mpi = ''
  petsc_incs = petsc_inc+petsc_bmake+petsc_mpi
  blas_libpath = '#'+petsc_base+'externalpackages/f2cblaslapack/linux-gnu'
  blas_libs = ['f2clapack', 'f2cblas']
  
  petsc_libpath = '#'+petsc_base+'lib/linux-gnu/'

Export("petsc_base", "petsc_inc", "petsc_bmake", "petsc_mpi", "petsc_incs", "petsc_libpath", "blas_libpath", "blas_libs")


## C++ build tools & MPI runner
if system_name == 'finarfin':
  mpirun = '/usr/bin/mpirun'
  if build.CompilerType() == 'intel':
    # Use intel compiler
    mpicxx = '/usr/bin/mpicxx -CC=icpc'
    cxx    = '/opt/intel_cc_80/bin/icpc'
    ar     = '/opt/intel_cc_80/bin/xiar'
  else:
    # Use gcc
    mpicxx = '/usr/bin/mpicxx'
    cxx    = '/usr/bin/g++'
    ar     = '/usr/bin/ar'
elif system_name == 'zuse':
   mpicxx = '/home/zuse/system/software/mpich-gcc/bin/mpicxx'
   mpirun = '/home/zuse/system/software/mpich-gcc/bin/mpirun'
   cxx = '/usr/bin/g++'
   ar = '/usr/bin/ar'
   
else:
  mpicxx = 'mpicxx'
  mpirun = 'mpirun'
  cxx = '/usr/bin/g++'
  ar = '/usr/bin/ar'

Export("mpicxx", "mpirun", "cxx", "ar")


## Any extra CCFLAGS and LINKFLAGS
extra_flags = build.CcFlags()
link_flags  = build.LinkFlags()

Export("extra_flags", "link_flags")

# Search path for #includes
import glob
cpppath = ['#/', '#/cxxtest']
src_folders = glob.glob('*/src')
for src_folder in src_folders:
  cpppath.append('#/'+src_folder)
  for dirpath, dirnames, filenames in os.walk(src_folder):
    for dirname in dirnames[:]:
      if dirname == '.svn':
        dirnames.remove(dirname)
      else:
        cpppath.append('#/'+os.path.join(dirpath, dirname))
Export("cpppath")

# Check for orphaned test files
os.system('python/TestRunner.py python/CheckForOrphanedTests.py ' +
          'testoutput/OrphanedTests.log ' + build_type + ' ' +
          build.GetTestReportDir() + ' --no-stdout')
# Check for duplicate file names in multiple directories
os.system('python/TestRunner.py python/CheckForDuplicateFileNames.py ' +
          'testoutput/DuplicateFileNames.log ' + build_type + ' ' +
          build.GetTestReportDir() + ' --no-stdout')

build_dir = build.build_dir
for toplevel_dir in ['linalg', 'mesh', 'global', 'io', 'ode', 'pde', 'coupled']:
    bld_dir = toplevel_dir + '/build/' + build_dir
    if not os.path.exists(bld_dir):
        os.mkdir(bld_dir)
    SConscript(toplevel_dir + '/SConscript', build_dir=bld_dir, duplicate=0)


# Remove the contents of testoutput/ on a clean build
test_output_files = glob.glob('testoutput/*')
Clean('.', test_output_files)



# Test summary generation
if test_summary and not compile_only:
  import socket, time
  # Touch a file, which we use as source for the summary target, so the summary
  # is done on every build.
  fp = file('buildtime.txt', 'w')
  print >>fp, time.asctime()
  fp.close()
  # Get the directory to put results & summary in
  machine = socket.getfqdn()
  output_dir = os.path.join(build.GetTestReportDir(), machine+'.'+build_type)
  # Remove old results. Note that this command gets run before anything is built.
  #for oldfile in os.listdir(output_dir):
  #  os.remove(os.path.join(output_dir, oldfile))
  # Add a summary generator to the list of things for scons to do
  if build_type == 'Coverage':
    # Remove old .gcda files before running more tests
    # First, find appropriate build directories
    build_dirs = glob.glob('*/build/' + build.build_dir)
    # Now find & remove .gcda files within there.
    # Also remove .log files so tests are re-run
    for build_dir in build_dirs:
      for dirpath, dirnames, filenames in os.walk(build_dir):
        for filename in filenames:
          if filename[-5:] == '.gcda' or filename[-4:] == '.log':
            os.remove(os.path.join(dirpath, filename))
    # For a Coverage build, run gcov & summarise instead
    summary = Builder(action = 'python python/DisplayCoverage.py ' + output_dir+' '+build_type)
  else:
    summary = Builder(action = 'python python/DisplayTests.py '+output_dir+' '+build_type)
  opt = Environment(ENV = {'PATH': os.environ['PATH'],
                           'PYTHONPATH': os.environ.get('PYTHONPATH', ''),
                           'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH', '')})
  opt['BUILDERS']['TestSummary'] = summary
  opt.TestSummary(os.path.join(output_dir, 'index.html'),
                  'buildtime.txt')
