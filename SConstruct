# Controlling scons build script for Chaste

# This script is executed within the root Chaste source directory

import sys
import os
import glob
import socket 

sys.path.append('python')
import BuildTypes

# The type of build to perform (see python/BuildTypes.py for options)
build_type = ARGUMENTS.get('build', 'default')
build = BuildTypes.GetBuildType(build_type)
build.SetRevision(ARGUMENTS.get('revision', ''))
Export('build', 'build_type')

# Specify test_summary=0 to scons to *NOT* generate a summary html page
test_summary = ARGUMENTS.get('test_summary', 1)

# Allow the system_name to be derived automatically
machine_fqdn = socket.getfqdn()
if machine_fqdn in ["userpc30.comlab.ox.ac.uk", "userpc33.comlab.ox.ac.uk"]:
    system_name = 'joe'
elif machine_fqdn in ["userpc44.comlab.ox.ac.uk", "userpc60.comlab.ox.ac.uk",
                      "userpc58.comlab.ox.ac.uk", "userpc59.comlab.ox.ac.uk"]:
    system_name = 'new_chaste'
elif machine_fqdn == "zuse.osc.ox.ac.uk":
    system_name = 'zuse'
elif machine_fqdn.endswith(".comlab.ox.ac.uk"):
    system_name = 'chaste'
elif machine_fqdn.startswith('finarfin'):
    system_name = 'finarfin'
elif machine_fqdn.endswith(".maths.nottingham.ac.uk"):
    system_name = 'Nottingham'
elif machine_fqdn.endswith(".maths.ox.ac.uk"):
    system_name = 'maths'
else:
    system_name = ''

# Specify system_name=<whatever> to scons to change default paths
system_name = ARGUMENTS.get('system_name', system_name)

# Specifying extra run-time flags
run_time_flags = ARGUMENTS.get('run_time_flags', '')
Export('run_time_flags')

# Specify all_tests=1 to select all tests for running (useful with
# compile_only=1)
all_tests = ARGUMENTS.get('all_tests', 0)
Export('all_tests')

# Specify compile_only=1 to not run any tests
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

# To run tests of only a single component, specify it with the
# test_component=<component> argument.
test_component = ARGUMENTS.get('test_component', '')
Export('test_component')


# Chaste components (top level dirs).
# We hard code dependencies between them, and use this to work out the
# order to link them in.  Each one is linked against just its dependencies,
# in the order given here.
comp_deps = {'models': ['ode', 'mesh', 'linalg', 'io', 'global'],
			 'dealii': ['models', 'coupled', 'pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
			 'coupled': ['pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
			 'pde': ['mesh', 'linalg', 'io', 'global'],
			 'mesh': ['linalg', 'global'],
			 'linalg': ['global'],
			 'ode': ['linalg', 'io', 'global'],
			 'io': ['global'],
			 'global': []}
components = ['models', 'coupled', 'pde', 'ode',
               'mesh', 'linalg', 'io', 'global']
if build.using_dealii:
    components = ['dealii'] + components
Export('components', 'comp_deps')


# Set extra paths to search for libraries and include files.
# Paths to PETSc, and any other external libraries, should be set here.
# The three variables exported are:
#  other_libs: names of libraries to link against.  Do not include PETSc
#              libraries.  Do not include the 'lib' prefix or any suffixes
#              in the library name (e.g. f2cblas not libf2cblas.a)
#  other_libpaths: paths to search for libraries to link against.  These
#                  should all be absolute paths, for safety.  Do include
#                  the path to PETSc libraries (unless they're in a standard
#                  system location).
#  other_includepaths: paths to search for header files.  Do include the
#                      path to PETSc headers (unless it's standard).
if system_name == 'finarfin':
  # Finarfin (Debian etch)
  petsc_base = '/home/jonc/work/dphil/petsc-2.3.1-p19/'
  if build.is_optimised:
      petsc_libpath = os.path.abspath(petsc_base+'lib/linux-gnu-opt/')
      petsc_bmake = petsc_base+'bmake/linux-gnu-opt'
  else:
      petsc_libpath = os.path.abspath(petsc_base+'lib/linux-gnu/')
      petsc_bmake = petsc_base+'bmake/linux-gnu'
  petsc_inc = petsc_base+'include'

  other_libs = ['lapack', 'blas', 'boost_serialization', 'xerces-c']
  other_libpaths = [petsc_libpath, '/usr/lib/atlas/sse2']
  other_includepaths = [petsc_inc, petsc_bmake]
elif system_name == 'maths':
  # Oxford uni maths inst
  petsc_base = '/scratch/chaste/petsc-2.3.2-p4/'
  petsc_inc = petsc_base+'include'
  petsc_bmake = petsc_base+'bmake/linux-gnu'
  petsc_mpi = petsc_base+'include/mpiuni'
  petsc_incs = [petsc_inc, petsc_bmake]
  petsc_libpath = petsc_base+'lib/linux-gnu/'

  other_libs = ['lapack', 'blas-3']
  other_libpaths = [petsc_libpath]
  other_includepaths = petsc_incs
elif system_name == 'joe' or system_name == 'intel_joe':
  # Joe Pitt-Francis userpc30 (Suse 9.3) and userpc33 (Ubuntu 6.06 Dapper Drake)
  petsc_base = '/home/jmpf/petsc-2.3.1-p16/'
  petsc_inc = petsc_base+'include'
  petsc_bmake = petsc_base+'bmake/linux-gnu'
  boost = '/home/jmpf'
  other_includepaths = [petsc_inc, petsc_bmake, boost]
  petsc_libpath = petsc_base+'lib/linux-gnu/'
  blas_libpath=  petsc_base+'externalpackages/f2cblaslapack/linux-gnu/'
  intel_libpath = '/opt/intel/cc/9.1.039/lib'
  other_libs = ['f2clapack', 'f2cblas','boost_serialization', 'xerces-c']
  other_libpaths = [petsc_libpath, blas_libpath, intel_libpath]
  xsd_inc = '/home/jmpf/xsd-2.3.1-i686-linux-gnu/libxsd'
  other_includepaths.extend([petsc_inc, petsc_bmake, xsd_inc])
elif system_name == 'zuse':
  petsc_base = '/home/zuse/system/software/petsc-2.2.1/'
  petsc_inc = petsc_base+'include'
  petsc_bmake = petsc_base+'bmake/linux-mpich-gnu-mkl'
  #petsc_mpi = petsc_base+'include/mpiuni'
  petsc_mpi = ''
  boost = '/home/zuse/system/joe'
  other_includepaths = [petsc_inc, petsc_bmake, boost]
  
  petsc_libpath = petsc_base+'lib/libg_c++/linux-mpich-gnu-mkl/'
  blas_libpath = '/opt/intel/mkl/8.0/lib/em64t/'

  other_libpaths = [petsc_libpath, blas_libpath]
  other_libs = []
elif system_name == 'zuse_opt':
  petsc_base = '/home/zuse/system/software/petsc-2.2.1/'
  petsc_inc = petsc_base+'include'
  petsc_bmake = petsc_base+'bmake/linux-mpich-gnu-mkl'
  #petsc_mpi = petsc_base+'include/mpiuni'
  petsc_mpi = ' '
  boost = '/home/zuse/system/joe'
  other_includepaths = [petsc_inc, petsc_bmake, boost]
  
  petsc_libpath = petsc_base+'lib/libg_c++/linux-mpich-gnu-mkl/'
  blas_libpath = '/opt/intel/mkl/8.0/lib/em64t/'
  other_libpaths = [petsc_libpath, blas_libpath,
                      '/home/zuse/system/software/opt/opt/lib/',
                      '/home/zuse/system/software/opt/opt-deps/gsoap/lib/',
                      '/home/zuse/system/software/opt/opt-deps/papi/lib/',
                      '/home/zuse/system/software/opt/opt-deps/libunwind/lib',
                      '/home/zuse/system/software/opt/opt-deps/papi/lib64']
  other_libs = ['opt', 'gsoap', 'stdc++', 'dl', 'papi','unwind-x86_64', 'unwind', 'perfctr']
elif system_name == 'chaste':
  # Chaste machines in comlab
  petsc_base = '../../../petsc-2.3.1-p13/'
  petsc_inc = petsc_base+'include'
  petsc_bmake = petsc_base+'bmake/linux-gnu'
  # petsc_mpi = petsc_base+'include/mpiuni'
  petsc_mpi = ''
  other_includepaths = [petsc_inc, petsc_bmake]
  blas_libpath = os.path.abspath(petsc_base+'externalpackages/f2cblaslapack/linux-gnu')
  petsc_libpath = os.path.abspath(petsc_base+'lib/linux-gnu/')
  other_libpaths = [petsc_libpath, blas_libpath]
  other_libs = ['f2clapack', 'f2cblas', 'boost_serialization']
elif (system_name == 'new_chaste' or system_name == 'intel_chaste'):
  # New Chaste machines in comlab
  other_includepaths = []
  if build.using_dealii:
    petsc_base = '../../../petsc-2.2.1/'
    dealii_base = '../../../deal.II/'
    Export('dealii_base')
    petsc_bmake = petsc_base+'bmake/linux-gnu'
    if build.is_optimised:
        petsc_libpath = os.path.abspath(petsc_base+'lib/libO_c++/linux-gnu/')
    else:
        petsc_libpath = os.path.abspath(petsc_base+'lib/libg_c++/linux-gnu/')
    dealii_libpath = os.path.abspath(dealii_base+'lib/')
    metis_libpath = os.path.abspath('../../../metis-4.0/')
    other_libpaths = [petsc_libpath, dealii_libpath, metis_libpath]
    metis_includepath = metis_libpath + '/Lib' # Bizarre, I know!
    dealii_includepaths = ['base/include', 'lac/include', 'deal.II/include']
    other_includepaths.extend(map(lambda s: dealii_base + s, dealii_includepaths))
    other_libs = build.GetDealiiLibraries(dealii_base) + ['blas', 'lapack', 'boost_serialization', 'xerces-c']
  else:
    petsc_base = '../../../petsc-2.3.2-p4/'
    if build.is_optimised:
        petsc_libpath = os.path.abspath(petsc_base+'lib/linux-gnu-opt/')
        petsc_bmake = petsc_base+'bmake/linux-gnu-opt'
    else:
        petsc_libpath = os.path.abspath(petsc_base+'lib/linux-gnu/')
        petsc_bmake = petsc_base+'bmake/linux-gnu'
    blas_libpath = os.path.abspath(petsc_base+'externalpackages/f2cblaslapack/linux-gnu')
    intel_libpath = '/opt/intel/cc/9.1.039/lib'
 
    other_libpaths = [petsc_libpath, blas_libpath, intel_libpath]
    other_libs = ['f2clapack', 'f2cblas', 'boost_serialization', 'xerces-c']
  petsc_inc = petsc_base+'include'
  # TODO: Make sure Chaste paths come first in the -I list.
  xsd_inc = '../../../xsd-2.3.1-i686-linux-gnu/libxsd'
  other_includepaths.extend([petsc_inc, petsc_bmake, xsd_inc])
elif system_name == 'Nottingham':
  # Gary, Alex and Helen's machines in Nottingham
  petsc_base = '/opt/petsc-2.2.1-with-mpi/'
  petsc_inc = petsc_base+'include'
  petsc_bmake = petsc_base+'bmake/linux-gnu'
  boost_path = '/opt/boost/include/boost-1_33_1'
  other_includepaths = [petsc_inc, petsc_bmake,  boost_path]
  other_libs = ['lapack', 'blas', 'boost_serialization-gcc','xerces-c']
  other_libpaths = [petsc_base+'lib/libO_c++/linux-gnu/',
                    '/opt/boost/lib/']
else:
  # Default for cancer course in the DTC
  petsc_base = '/usr/local/petsc-2.3.1-p15/'
  petsc_inc = petsc_base+'include'
  petsc_bmake = petsc_base+'bmake/linux-gnu'
  other_includepaths = [petsc_inc, petsc_bmake]
  other_libs = ['lapack', 'blas']
  other_libpaths = [petsc_base+'lib/linux-gnu/']

Export("other_includepaths", "other_libpaths", "other_libs")


## C++ build tools & MPI runner
if system_name == 'finarfin':
  mpirun = 'mpirun'
  if build.CompilerType() == 'intel':
    # Use intel compiler
    mpicxx = '/usr/bin/mpicxx -CC=icpc'
    cxx    = '/opt/intel_cc_80/bin/icpc'
    ar     = '/opt/intel_cc_80/bin/xiar'
  else:
    # Use gcc
    mpicxx = 'mpicxx'
    cxx    = 'g++'
    ar     = 'ar'
elif system_name == 'maths':
  mpicxx = 'mpicxx'
  mpirun = 'mpirun'
  cxx = 'g++'
  ar = 'ar'
elif system_name == 'zuse':
   mpicxx = '/home/zuse/system/software/mpich-gcc/bin/mpicxx'
   mpirun = '/home/zuse/system/software/mpich-gcc/bin/mpirun'
   cxx = '/usr/bin/g++'
   ar = '/usr/bin/ar'
elif system_name == 'zuse_opt':
   mpicxx = '/home/zuse/system/software/mpich-gcc/bin/mpicxx'
   mpirun = '/home/zuse/system/software/mpich-gcc/bin/mpirun'
   cxx = '/usr/bin/g++'
   ar = '/usr/bin/ar'
elif system_name == 'joe':
  mpicxx = '/home/jmpf/mpi/bin/mpicxx'
  mpirun = '/home/jmpf/mpi/bin/mpirun'
  cxx = '/usr/bin/g++'
  ar = '/usr/bin/ar'
elif system_name == 'intel_joe':
  mpicxx = 'mpicxx -CC=icpc'
  mpirun = 'mpirun'
  cxx = '/opt/intel/cc/9.1.039/bin/icpc'
  ar = ' /opt/intel/cc/9.1.039/bin/xiar'
elif system_name == 'chaste':
  mpicxx = 'mpicxx'
  mpirun = 'mpirun'
  cxx = '/usr/bin/g++'
  ar = '/usr/bin/ar'
elif system_name == 'new_chaste':
    mpirun = 'mpirun'
    if build.CompilerType() == 'intel':
        mpicxx = 'mpicxx -CC=icpc'
        cxx = '/opt/intel/cc/9.1.039/bin/icpc'
        ar = ' /opt/intel/cc/9.1.039/bin/xiar'
    else:    
        mpicxx = 'mpicxx'
        cxx = '/usr/bin/g++'
        ar = '/usr/bin/ar'
elif system_name == 'Nottingham':
  mpicxx = '/opt/mpi/bin/mpicxx'
  mpirun = '/opt/mpi/bin/mpirun'
  cxx = '/usr/bin/g++'
  ar = '/usr/bin/ar'
else:
  # DTC cancer course defaults
  mpicxx = '/usr/local/mpi/bin/mpicxx'
  mpirun = '/usr/local/mpi/bin/mpirun'
  cxx = '/usr/bin/g++'
  ar = '/usr/bin/ar'

build.tools['mpicxx'] = mpicxx
build.tools['mpirun'] = mpirun
build.tools['cxx'] = cxx
build.tools['ar'] = ar

# Find full path to valgrind, as parallel memory testing needs it to be
# given explicitly.
# We search on os.environ['PATH'] for now.  When #258 is done use the environment.
vg_path = WhereIs(build.tools['valgrind'])
if vg_path:
    build.tools['valgrind'] = vg_path
del vg_path


## Any extra CCFLAGS and LINKFLAGS
extra_flags = build.CcFlags()
link_flags  = build.LinkFlags()

# Hack to get around Debian sarge strangeness
if system_name in ['maths']:
    extra_flags = extra_flags + " -DCWD_HACK "

if system_name == 'Nottingham':
    extra_flags = "-isystem " + boost_path + " " + extra_flags

Export("extra_flags", "link_flags")

# Search path for Chaste #includes
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
          build.GetTestReportDir() + 'OrphanedTests.log ' + build_type + ' --no-stdout')
# Check for duplicate file names in multiple directories
os.system('python/TestRunner.py python/CheckForDuplicateFileNames.py ' +
          build.GetTestReportDir() + 'DuplicateFileNames.log ' + build_type + ' --no-stdout')

build_dir = build.build_dir
test_depends = [File(build.GetTestReportDir() + 'OrphanedTests.log'),
                File(build.GetTestReportDir() + 'DuplicateFileNames.log')]
for toplevel_dir in components:
    bld_dir = os.path.join(toplevel_dir, 'build', build_dir)
    if not os.path.exists(bld_dir):
        os.mkdir(bld_dir)
    test_depends.append(SConscript('SConscript', src_dir=toplevel_dir, build_dir=bld_dir,
                                   duplicate=0))


# Remove the contents of build.GetTestReportDir() on a clean build
test_output_files = glob.glob(build.GetTestReportDir() + '*')
Clean('.', test_output_files)
# Also remove the entire build.build_dir for each component, so we
# don't have stale tests, etc. still present
for toplevel_dir in components:
    Clean('.', os.path.join(toplevel_dir, 'build', build_dir))
# Also make sure we remove any libraries still hanging around, just in case
for lib in glob.glob('lib/*'):
    Clean('.', lib)
for lib in glob.glob('linklib/*'):
    Clean('.', lib)



# Test summary generation
if test_summary and not compile_only:
  # Get the directory to put results & summary in
  output_dir = build.output_dir
  # Remove old results. Note that this command gets run before anything is built.
  #for oldfile in os.listdir(output_dir):
  #  os.remove(os.path.join(output_dir, oldfile))
  # Add a summary generator to the list of things for scons to do
  if isinstance(build, BuildTypes.Coverage):
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
  def output_dir_lister(target, source, env, output_dir=output_dir):
    """Create a file containing a directory listing."""
    files = filter(lambda f: f[0] != '.', os.listdir(output_dir))
    files.sort()
    fp = file(str(target[0]), 'w')
    for f in files:
      fp.write(f + '\n')
    fp.close()
    return None # Successful build
  lister = Action(output_dir_lister,
                  lambda ts, ss, env: "Generating file list %s" % ts[0])
  opt = Environment(ENV = {'PATH': os.environ['PATH'],
                           'PYTHONPATH': os.environ.get('PYTHONPATH', ''),
                           'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH', '')})
  opt['BUILDERS']['TestSummary'] = summary
  # Make sure we use the file list contents as its signature
  opt.TargetSignatures('content')
  file_list = File(os.path.join(output_dir, '.filelist'))
  opt.AlwaysBuild(file_list)
  opt.Command(file_list, test_depends, lister)
  opt.TestSummary(os.path.join(output_dir, 'index.html'), [file_list, test_depends])
