"""Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.
"""

# Controlling SCons build script for Chaste.

# This script is executed within the root Chaste source directory.
# We need at least Python 2.3.
EnsurePythonVersion(2,3)

# We're also no longer compatible with SCons 0.96
EnsureSConsVersion(0,97)

Help("""
  Type: 'scons -c .' to remove all the compiled files (clean build),
        'scons' to do a default build,
        'scons test_suite=<Path from chaste folder>' to run a single test,
        'scons <component>' to build and test a single component.
  
  For other options, such as profiling, optimised builds and 
  memory testing please refer to:
  
  https://chaste.comlab.ox.ac.uk/cgi-bin/trac.cgi/wiki/ChasteGuides/BuildGuide

""")

import sys
import os
import glob
import socket
import time

import SCons

sys.path[0:0] = ['python', 'python/hostconfig']
import BuildTypes
import SConsTools
Export('SConsTools')

import hostconfig

# If building a loadable module at run-time
dyn_libs_only = int(ARGUMENTS.get('dyn_libs_only', 0))
Export('dyn_libs_only')
if dyn_libs_only:
    #print sys.argv
    # Set some other options
    ARGUMENTS['test_summary'] = 0
    ARGUMENTS['do_inf_tests'] = 0
    # Note what folder is being built
    dyn_folder = os.path.join(Dir('#').abspath, COMMAND_LINE_TARGETS[0])
    Export('dyn_folder')

# Turn on some build-script debugging?
debug = int(ARGUMENTS.get('debug', 0))
Export('debug')

# The type of build to perform (see python/BuildTypes.py for options)
build_type = ARGUMENTS.get('build', ARGUMENTS.get('b','default'))
build = BuildTypes.GetBuildType(build_type)
build.SetRevision(ARGUMENTS.get('revision', ''))
build.debug = debug
Export('build')

# Whether to use static or shared libraries
static_libs = int(ARGUMENTS.get('static', 0))
if build.is_profile:
    static_libs = 1
Export('static_libs')

# Whether to build Chaste libraries, or link tests against object files directly
use_chaste_libs = int(ARGUMENTS.get('chaste_libs',  ARGUMENTS.get('cl', 0)))
Export('use_chaste_libs')

# Specify test_summary=0 to scons to *NOT* generate a summary html page
test_summary = int(ARGUMENTS.get('test_summary', 1))

# Used by the automated build system
run_infrastructure_tests = int(ARGUMENTS.get('do_inf_tests', 1))
check_failing_tests = int(ARGUMENTS.get('check_failing_tests', 0))

# Specifying extra run-time flags
run_time_flags = ARGUMENTS.get('run_time_flags', '')

# Specify all_tests=1 to select all tests for running (useful with
# compile_only=1)
all_tests = int(ARGUMENTS.get('all_tests', ARGUMENTS.get('at', 0)))
Export('all_tests')

# Specify compile_only=1 to not run any tests
compile_only = int(ARGUMENTS.get('compile_only', ARGUMENTS.get('co', 0)))
Export('compile_only')

# To run a single test suite only, give its path (relative to the Chaste
# root) as the test_suite=<path> argument.
# This will force the test suite to be run even if the source is unchanged.
single_test_suite = ARGUMENTS.get('test_suite', ARGUMENTS.get('ts', ''))
if single_test_suite:
    single_test_suite = single_test_suite.split(os.path.sep)
    if (len(single_test_suite)<2):
        raise ValueError('Path to test suite is too short')
    for i in [-2, -3, -4]:
        if single_test_suite[i] == 'test':
            single_test_suite_dir = single_test_suite[i-1]
            single_test_suite = os.path.sep.join(single_test_suite[i+1:])
            break
    else:
        raise ValueError('Test suite is not in a test folder')
    #print single_test_suite, single_test_suite_dir
else:
    single_test_suite_dir = ''
Export('single_test_suite', 'single_test_suite_dir')

# Force re-running of all (selected) tests even if the source is unchanged.
force_test_runs = int(ARGUMENTS.get('force_test_runs', 0))
Export('force_test_runs')

# Don't update the provenance information (Version.cpp file).
update_provenance = int(ARGUMENTS.get('update_provenance', ARGUMENTS.get('up', 1)))

# Whether to kill the tests if they run too long (limit in seconds, 0 means don't kill)
test_time_limit = int(ARGUMENTS.get('test_time_limit', 0))

# Whether to build executables, or just tests
build_exes = int(ARGUMENTS.get('exe', 0))
Export('build_exes')
if build_exes:
    assert use_chaste_libs, "Cannot build executables unless building Chaste libraries"


# Experimental support for installing Chaste as a normal collection of
# libraries and headers.
install_prefix = ARGUMENTS.get('install_prefix', '/usr/local')
Export('install_prefix')
if 'install' in BUILD_TARGETS:
    assert use_chaste_libs, "Cannot install unless building Chaste libraries"

# Check for an easy mistake, where the user forgets the 'test_suite='.
for target in BUILD_TARGETS:
    if target.endswith('.hpp'):
        raise ValueError('Unexpected target ' + target +
                         '; did you forget "test_suite="?')

# To run tests of only a single component, specify it with the
# test_component=<component> argument (deprecated).
test_component = ARGUMENTS.get('test_component', '')
Export('test_component')


# If building static libraries, get rid of any old shared libraries,
# in order to stop the automatic dependency algorithm getting confused.
if use_chaste_libs and static_libs:
    for lib in glob.glob('lib/lib*.so'):
        Execute(Delete(lib))


# Use a single file to store signatures.
# Forwards-compatible with SCons 0.97, and nicer for svn ignore.
if not dyn_libs_only:
    if sys.platform == 'cygwin':
        SConsignFile('.sconsign-cygwin')
    else:
        SConsignFile('.sconsign')
else:
    # Use a .sconsign file in the folder we're building to avoid conflicts.
    assert(len(COMMAND_LINE_TARGETS) == 1)
    SConsignFile(os.path.join(COMMAND_LINE_TARGETS[0], '.sconsign'))

# Chaste components (top level dirs).
# We hard code dependencies between them, and use this to work out the
# order to link them in.  Each one is linked against just its dependencies,
# in the order given here.
comp_deps = {'cell_based': ['pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'crypt': ['cell_based', 'pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'notforrelease': ['heart', 'pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'notforrelease_cell_based': ['crypt', 'cell_based', 'pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'heart': ['pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'pde': ['mesh', 'linalg', 'io', 'global'],
             'mesh': ['linalg', 'global'],
             'linalg': ['global'],
             'ode': ['linalg', 'io', 'global'],
             'io': ['global'],
             'global': [],
             'core': ['pde', 'ode', 'mesh', 'linalg', 'io', 'global']}
SConsTools.comp_deps = comp_deps
components = ['python', 'global', 'io', 'linalg', 'mesh', 'ode', 'pde',
              'heart', 'cell_based', 'crypt', 'notforrelease', 'notforrelease_cell_based']
# Ignore non-existent components
# e.g. notforrelease wont appear in a release version
for comp in components[:]:
    if not os.path.isdir(comp):
        components.remove(comp)
Export('components', 'comp_deps')

Alias('core', Split('global io linalg mesh ode pde'))

# Set extra paths to search for libraries and include files.
# Paths to PETSc, and any other external libraries, should be set here.
# The three variables exported are:
#   other_libs: names of libraries to link against.
#   other_libpaths: paths to search for libraries to link against.
#   other_includepaths: paths to search for header files.
# This is now done by the hostconfig subsystem.
hostconfig.Configure(build)
other_libs = hostconfig.libraries
other_libpaths = hostconfig.libpaths
other_includepaths = hostconfig.incpaths
if isinstance(build, BuildTypes.CovTool):
    build.UseCovTool(other_includepaths, other_libs)

Export("other_libpaths", "other_libs")


# Set up the environment to use for building.
other_libpaths.append(os.path.abspath('lib'))
if os.environ.get('XTPE_COMPILE_TARGET', ''):
    env = Environment(ENV = os.environ)
else:
    env = Environment(
        #tools = ['g++', 'gnulink', 'gas', 'ar', 'g77'],
        ENV={'PATH': '.:' + os.environ['PATH'],
             'PYTHONPATH': os.environ.get('PYTHONPATH', ''),
             'USER': os.environ['USER'],
             'INTEL_LICENSE_FILE': '28518@lic1.osc.ox.ac.uk:' +
                                   os.environ.get('INTEL_LICENSE_FILE', '.'),
             'CHASTE_TEST_OUTPUT':
             os.environ.get('CHASTE_TEST_OUTPUT',
                            '/tmp/'+os.environ['USER']+'/testoutput/'),
             'LD_LIBRARY_PATH': ':'.join(other_libpaths),
             'HOME': os.environ['HOME']
            })
env.Append(BOPT = 'g_c++') # Needed for some versions of PETSc?
env.Replace(CXX = build.tools['mpicxx'])
env.Replace(AR = build.tools['ar'])
env.Replace(CXXFILESUFFIX = '.cpp')
env['INSTALL_PREFIX'] = install_prefix

if int(ARGUMENTS.get('br', ARGUMENTS.get('brief', 0))):
    env.Replace(CXXCOMSTR = '$CXX -o $TARGET -c <flags etc. omitted> $SOURCES')
    env.Replace(SHCXXCOMSTR = '$SHCXX -o $TARGET -c <flags etc. omitted> $SOURCES')

# Any extra CCFLAGS and LINKFLAGS
extra_flags = build.CcFlags() + ' ' + hostconfig.CcFlags()
link_flags  = build.LinkFlags() + ' ' + hostconfig.LdFlags()
include_flag = ' ' + build.IncludeFlag() + ' '
env.Append(CCFLAGS = include_flag + include_flag.join(other_includepaths)
           + ' ' + extra_flags)
env.Append(LINKFLAGS = link_flags)
env.Append(CPPDEFINES = hostconfig.CppDefines() + ['TRILIBRARY', 'TETLIBRARY', 'ANSI_DECLARATORS'])

# Search path for Chaste #includes
cpppath = [Dir('.'), Dir('cxxtest')]
for component in components:
    src_folder = os.path.join(component, 'src')
    src_dirs = SConsTools.FindSourceFiles(env, src_folder, dirsOnly=True, includeRoot=True)
    bld_dir = os.path.join(component, 'build', build.build_dir)
    clen = len(component)
    for d in src_dirs:
        cpppath.extend([Dir(d), Dir(bld_dir + d[clen:])])
cpppath = map(lambda p: Dir(p), cpppath)
env.Replace(CPPPATH = cpppath)

# Some state needed by our build system
env['build'] = build
env['buildsig'] = build.GetSignature()
env['CHASTE_COMPONENTS'] = components + ['projects']
env['CHASTE_OBJECTS'] = {}
env['UPDATE_CHASTE_PROVENANCE'] = update_provenance

if not single_test_suite:
    # Default is to build all components, but not user projects
    Default(components)


# Create Builders for generating test .cpp files, and running test executables
test = Builder(action='cxxtest/cxxtestgen.py --error-printer -o $TARGET $SOURCES')
import TestRunner
def TestDescription(target, source, env):
    return "Running '%s'" % (source[0])
test_action = Action(TestRunner.get_build_function(build, run_time_flags, test_time_limit),
                     TestDescription, varlist=['buildsig'])
runtest = Builder(action=test_action)
env['BUILDERS']['Test'] = test
env['BUILDERS']['RunTest'] = runtest
env['BUILDERS']['BuildTest'] = Builder(action=SConsTools.BuildTest)

# Faster builds of shared libraries
import fasterSharedLibrary
fasterSharedLibrary.Copy = Copy # Bit of a hack this!
env['BUILDERS']['OriginalSharedLibrary'] = env['BUILDERS']['SharedLibrary']
env['BUILDERS']['SharedLibrary'] = fasterSharedLibrary.fasterSharedLibrary

# Builder for generating C++ code from XML Schema files
SConsTools.CreateXsdBuilder(build, env, dyn_libs_only)

# Builder for generating C++ code from CellML files
SConsTools.CreatePyCmlBuilder(build, env)

# Find full path to valgrind, as parallel memory testing needs it to be
# given explicitly.
vg_path = env.WhereIs(build.tools['valgrind'])
if vg_path:
    build.tools['valgrind'] = vg_path
del vg_path

# Record key build info for the provenance system
SConsTools.RecordBuildInfo(env, build_type, static_libs, use_chaste_libs)

# Allow hostconfig scripts to modify the build environment
if hasattr(hostconfig.conf, 'ModifyEnv') and callable(hostconfig.conf.ModifyEnv):
    hostconfig.conf.ModifyEnv(env)

# We need different linker flags when compiling dynamically loadable modules
dynenv = env.Clone()
env.Append(LINKFLAGS=' '+build.rdynamic_link_flag)
# Try to avoid likely conflicts
for path in [p for p in dynenv['CPPPATH'] if 'cellml' in str(p)]:
    dynenv['CPPPATH'].remove(path)

# Export the build environment to SConscript files
Export('env', 'dynenv')

# Test log files to summarise
test_log_files = []

if run_infrastructure_tests and not GetOption('clean'):
    # Check for orphaned test files
    out = File(build.GetTestReportDir() + 'OrphanedTests.log')
    os.system('python/TestRunner.py python/CheckForOrphanedTests.py ' +
              str(out) + ' ' + build_type + ' --no-stdout')
    test_log_files.append(out)
    # Check for duplicate file names in multiple directories
    out = File(build.GetTestReportDir() + 'DuplicateFileNames.log')
    os.system('python/TestRunner.py python/CheckForDuplicateFileNames.py ' +
              str(out) + ' ' + build_type + ' --no-stdout')
    test_log_files.append(out)
    # Check for duplicate file names in multiple directories
    out = File(build.GetTestReportDir() + 'Copyrights.log')
    os.system('python/TestRunner.py python/CheckForCopyrights.py ' +
              str(out) + ' ' + build_type + ' --no-stdout')
    test_log_files.append(out)
    ## Do not check for stale semaphores - it's only important on MPICH with Gnu Linux
    #out = File(build.GetTestReportDir() + 'Semaphores.log')
    #os.system('python/TestRunner.py python/CheckSemaphores.py ' +
    #          str(out) + ' ' + build_type + ' --no-stdout')
    #test_log_files.append(out)
    # Check for stray schemas
    out = File(build.GetTestReportDir() + 'Schemas.log')
    os.system('python/TestRunner.py python/CheckSchemas.py ' +
              str(out) + ' ' + build_type + ' --no-stdout')
    test_log_files.append(out)
if check_failing_tests:
    out = File(build.GetTestReportDir() + 'FailingTests.log')
    os.system('python/TestRunner.py python/CheckForFailingTests.py ' +
              str(out) + ' ' + build_type + ' --no-stdout')
    test_log_files.append(out)

build_dir = build.build_dir
if not isinstance(build, BuildTypes.DoxygenCoverage):
    # Build each component.
    for toplevel_dir in components:
        bld_dir = os.path.join(toplevel_dir, 'build', build_dir)
        if not os.path.exists(bld_dir):
            os.mkdir(bld_dir)
        script = os.path.join(toplevel_dir, 'SConscript')
        (test_logs, lib) = SConscript(script, src_dir=toplevel_dir, build_dir=bld_dir,
                                      duplicate=0)
        if not lib:
            # This component hasn't created a library file, so don't try to link
            # against it.  Happens if the component consists only of headers.
            for v in comp_deps.itervalues():
                try:
                    v.remove(toplevel_dir)
                except ValueError:
                    pass
        test_log_files.append(test_logs)
    
    # Any user projects?
    for project in glob.glob('projects/[_a-zA-Z]*'):
        if not os.path.isdir(project):
            if debug:
                print "Found non-dir", project, "in projects folder"
            continue
        if not os.path.exists(os.path.join(project, 'SConscript')):
            print >>sys.stderr, "Unexpected folder", project, "in projects folder."
            continue
        if ('install' in BUILD_TARGETS and not
            (project in BUILD_TARGETS or project+os.sep in BUILD_TARGETS)):
            # Only install projects if explicitly requested
            continue
        bld_dir = os.path.join(project, 'build', build_dir)
        if not os.path.exists(bld_dir):
            os.makedirs(bld_dir)
        script = os.path.join(project, 'SConscript')
        test_log_files.append(SConscript(script, src_dir=project, build_dir=bld_dir,
                                         duplicate=0))
    
    # Make sure test executables get built if compile_only=1
    env.Default(env.Alias('test_exes'))
    # Work around an error on some SCons versions, by ensuring the alias has a builder
    import SCons.Environment
    env.Alias('test_exes')[0].builder_set(SCons.Environment.AliasBuilder)


# Remove the contents of build.GetTestReportDir() on a clean build
test_output_files = glob.glob(build.GetTestReportDir() + '*')
Clean('.', test_output_files)
# Also remove the entire build.build_dir for each component, so we
# don't have stale tests, etc. still present
for toplevel_dir in components:
    Clean('.', os.path.join(toplevel_dir, 'build', build_dir))
# Also make sure we remove any libraries still hanging around, just in case
Clean('.', glob.glob('lib/*'))
Clean('.', glob.glob('linklib/*'))


def RequestedProjects():
    """Return a list of projects explicitly mentioned on the command line."""
    projects = []
    for targ in COMMAND_LINE_TARGETS:
        if str(targ).startswith('projects'):
            projects.append(str(targ))
    return projects

# Test summary generation
if test_summary and not compile_only:
    # Copy the build env, since we change TargetSigs
    senv = env.Clone()
    # Get the directory to put results & summary in
    output_dir = build.output_dir
    Execute(Mkdir(output_dir))
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
        summary_action = ('python python/DisplayCoverage.py ' + output_dir + ' ' + build_type
                          + ' ' + ' '.join(RequestedProjects()))
    elif isinstance(build, BuildTypes.DoxygenCoverage):
        # Run Doxygen and parse the output
        doxy_conf = ['cat Doxyfile', 'echo "PROJECT_NUMBER=Build::r%s"' % build._revision]
        # Include projects?
        project_inputs = RequestedProjects()
        if project_inputs:
            doxy_conf.append('echo "INPUT += %s"' % (' '.join(map(lambda p: os.path.join(p, 'src'), project_inputs))))
        build.ExtendDoxygenConfig(doxy_conf)
        cmd = ('( ' + ' ; '.join(doxy_conf) + ' ) '
               + '| doxygen - 2>doxygen-error.log 1>doxygen-output.log')
        summary_action = cmd + '; python python/ParseDoxygen.py doxygen-output.log doxygen-error.log ' + output_dir
    else:
        summary_action = 'python python/DisplayTests.py '+output_dir+' '+build_type
  
    summary_index = os.path.join(output_dir, 'index.html')
    senv.Command(summary_index, Flatten(test_log_files), summary_action)
    # Avoid circular dependencies
    senv.Ignore(summary_index, summary_index)
    senv.Ignore(Dir(output_dir), summary_index)
    # Make sure the summary is always required by any build targets requested
    # explicitly
    for targ in COMMAND_LINE_TARGETS:
        senv.Depends(targ, summary_index)
    senv.Default(summary_index)
    # Also allow other code to add dependencies, by making the summary depend on an Alias
    senv.Depends(summary_index, senv.Alias('test_summary_dependencies'))
    # Try to ensure it runs even if SCons thinks it's up-to-date, just to
    # re-assure the user
    senv.AlwaysBuild(summary_index)


if build_exes:
    SConsTools.BuildExes(build, env, 'apps',
                         components=['heart']+comp_deps['heart'],
                         otherVars=globals())
