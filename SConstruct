"""Copyright (C) University of Oxford, 2005-2009

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
  Type: 'scons -c' to remove all the compiled files (clean build),
        'scons' to do a default build,
        'scons test_suite=<Path from chaste folder>' to run a single test,
        'scons <component>' to build and test a single component.
  
  For other options, such as profiling, optimised builds and 
  memory testing please refer to:
  
  https://chaste.comlab.ox.ac.uk/cgi-bin/trac.cgi/wiki/BuildGuide

""")

import sys
import os
import glob
import socket
import time

sys.path.append('python')
import BuildTypes
import SConsTools
Export('SConsTools')

sys.path.append('python/hostconfig')
import hostconfig

# The type of build to perform (see python/BuildTypes.py for options)
build_type = ARGUMENTS.get('build', 'default')
build = BuildTypes.GetBuildType(build_type)
build.SetRevision(ARGUMENTS.get('revision', ''))
Export('build')

# Whether to use static or shared libraries
static_libs = int(ARGUMENTS.get('static', 0))
if build.is_profile:
    static_libs = 1
Export('static_libs')

# Whether to build Chaste libraries, or link tests against object files directly
use_chaste_libs = int(ARGUMENTS.get('chaste_libs', 0))
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
all_tests = int(ARGUMENTS.get('all_tests', 0))
Export('all_tests')

# Specify compile_only=1 to not run any tests
compile_only = int(ARGUMENTS.get('compile_only', 0))
Export('compile_only')

# To run a single test suite only, give its path (relative to the Chaste
# root) as the test_suite=<path> argument.
# This will force the test suite to be run even if the source is unchanged.
single_test_suite = ARGUMENTS.get('test_suite', '')
if single_test_suite:
    single_test_suite = single_test_suite.split(os.path.sep)
    if (len(single_test_suite)<2):
        raise ValueError('Path to test suite is too short')
    for i in [-2, -3]:
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
force_test_runs = bool(ARGUMENTS.get('force_test_runs', 0))
Export('force_test_runs')

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
if sys.platform == 'cygwin':
    SConsignFile('.sconsign-cygwin')
else:
    SConsignFile('.sconsign')


# Chaste components (top level dirs).
# We hard code dependencies between them, and use this to work out the
# order to link them in.  Each one is linked against just its dependencies,
# in the order given here.
comp_deps = {'cancer': ['pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'notforrelease': ['heart', 'pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'notforrelease_cancer': ['cancer', 'pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'heart': ['pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'pde': ['mesh', 'linalg', 'io', 'global'],
             'mesh': ['linalg', 'global'],
             'linalg': ['global'],
             'ode': ['linalg', 'io', 'global'],
             'io': ['global'],
             'global': [],
             'core': ['pde', 'ode', 'mesh', 'linalg', 'io', 'global']}
SConsTools.comp_deps = comp_deps
components = ['global', 'io', 'linalg', 'mesh', 'ode', 'pde',
              'heart', 'cancer', 'notforrelease', 'notforrelease_cancer']
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
#  other_libs: names of libraries to link against.  Do not include PETSc
#              libraries.  Do not include the 'lib' prefix or any suffixes
#              in the library name (e.g. f2cblas not libf2cblas.a)
#  other_libpaths: paths to search for libraries to link against.  These
#                  should all be absolute paths, for safety.  Do include
#                  the path to PETSc libraries (unless they're in a standard
#                  system location).
#  other_includepaths: paths to search for header files.  Do include the
#                      path to PETSc headers (unless it's standard).
# This is now done by the hostconfig subsystem.
hostconfig.configure(build)
other_libs = hostconfig.libraries
other_libpaths = hostconfig.libpaths
other_includepaths = hostconfig.incpaths
if isinstance(build, BuildTypes.CovTool):
    build.UseCovTool(other_includepaths, other_libs)

Export("other_libpaths", "other_libs")


# Any extra CCFLAGS and LINKFLAGS
extra_flags = build.CcFlags() + ' ' + hostconfig.ccflags() \
              + ' -DTRILIBRARY -DANSI_DECLARATORS '
link_flags  = build.LinkFlags() + ' ' + hostconfig.ldflags()

# Search path for Chaste #includes
cpppath = ['.', 'cxxtest']
src_folders = glob.glob('*/src')
for src_folder in src_folders:
    cpppath.extend(SConsTools.FindSourceFiles(src_folder, dirsOnly=True, includeRoot=True))
cpppath = map(lambda p: '#/'+p, cpppath)


# Set up the environment to use for building.
other_libpaths.append(os.path.abspath('lib'))
env = Environment(
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
env.Append(CCFLAGS = '-isystem ' + ' -isystem '.join(other_includepaths)
           + ' ' + extra_flags)
env.Append(LINKFLAGS = link_flags)
env.Append(BOPT = 'g_c++')
env.Replace(CXX = build.tools['mpicxx'])
env.Replace(AR = build.tools['ar'])
env.Replace(CPPPATH = cpppath)
env['buildsig'] = build.GetSignature()
env['CHASTE_COMPONENTS'] = components + ['projects']
env['CHASTE_OBJECTS'] = {}


if not single_test_suite:
    # Default is to build all components, but not user projects
    Default(components)


# Create Builders for generating test .cpp files, and running test executables
test = Builder(action='cxxtest/cxxtestgen.py --error-printer -o $TARGET $SOURCES')
import TestRunner
def test_description(target, source, env):
    return "running '%s'" % (source[0])
test_action = Action(TestRunner.get_build_function(build, run_time_flags),
                     test_description, varlist=['buildsig'])
runtest = Builder(action=test_action)
env['BUILDERS']['Test'] = test
env['BUILDERS']['RunTest'] = runtest
env['BUILDERS']['BuildTest'] = Builder(action=SConsTools.BuildTest)

# Faster builds of shared libraries
import fasterSharedLibrary
fasterSharedLibrary.Copy = Copy # Bit of a hack this!
env['BUILDERS']['OriginalSharedLibrary'] = env['BUILDERS']['SharedLibrary']
env['BUILDERS']['SharedLibrary'] = fasterSharedLibrary.fasterSharedLibrary

# 'Builder' for running xsd to generate parser code from an XML schema.
# Getting this to integrate with the build is non-trivial, so we cheat.
def run_xsd(schema_file):
    output_dir = os.path.dirname(schema_file)
    command = ' '.join([build.tools['xsd'], 'cxx-tree',
                        '--generate-serialization',
                        '--output-dir', output_dir,
                        '--hxx-suffix', '.hpp', '--cxx-suffix', '.cpp',
                        '--prologue-file', 'heart/src/io/XsdPrologue.txt',
                        '--epilogue-file', 'heart/src/io/XsdEpilogue.txt',
                        '--namespace-regex', "'X.* $Xchaste::parametersX'",
                        '--namespace-regex', "'X.* https://chaste.comlab.ox.ac.uk/nss/parameters/(.+)Xchaste::parameters::v$1X'",
                        schema_file])
    print "Running xsd on", schema_file
    os.system(command)
# Check if we need to run XSD (and that 'xsd' is really xsd...)
command = build.tools['xsd'] + ' version 2>&1'
xsd_version_string = os.popen(command).readline().strip()
if xsd_version_string.startswith('XML Schema Definition Compiler'):
    xsd_version = 2
elif xsd_version_string.startswith('CodeSynthesis XSD XML Schema to C++ compiler'):
    xsd_version = 3
else:
    print "Unexpected XSD program found:"
    print xsd_version_string
    sys.exit(1)
# If it's the old version, always run XSD; otherwise only run if generated code is out of date
for schema_file in glob.glob('heart/src/io/ChasteParameters*.xsd'):
    cpp_file = schema_file[:-3] + 'cpp'
    if (xsd_version == 2 or
        not os.path.exists(cpp_file) or
        os.stat(schema_file).st_mtime > os.stat(cpp_file).st_mtime):
        run_xsd(schema_file)

# Find full path to valgrind, as parallel memory testing needs it to be
# given explicitly.
vg_path = env.WhereIs(build.tools['valgrind'])
if vg_path:
    build.tools['valgrind'] = vg_path
del vg_path

# Export the build environment to SConscript files
Export('env')

# Test log files to summarise
test_log_files = []

if run_infrastructure_tests:
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
    for project in glob.glob('projects/[_a-zA-z]*'):
	if not os.path.exists(os.path.join(project, 'SConscript')):
            print >>sys.stderr, "Unexpected folder", project, "in projects folder."
            continue
        bld_dir = os.path.join(project, 'build', build_dir)
        if not os.path.exists(bld_dir):
            os.mkdir(bld_dir)
        script = os.path.join(project, 'SConscript')
        test_log_files.append(SConscript(script, src_dir=project, build_dir=bld_dir,
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
    # Copy the build env, since we change TargetSigs
    senv = env.Copy()
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
        summary_action = 'python python/DisplayCoverage.py ' + output_dir+' '+build_type
    elif isinstance(build, BuildTypes.DoxygenCoverage):
        # Run Doxygen and parse the output
        cmd = '( cat Doxyfile ; echo "PROJECT_NUMBER=Build::r%s" ) ' % build._revision \
            + '| doxygen - 2>doxygen-error.log 1>doxygen-output.log'
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
    # Try to ensure it runs even if SCons thinks it's up-to-date, just to
    # re-assure the user
    senv.AlwaysBuild(summary_index)


if ARGUMENTS.get('exe', 0):
    assert use_chaste_libs
    env = env.Copy()
    # Build information to supply to the executable
    svn_rev = os.popen("svnversion").read().strip()
    uname = ' '.join(os.uname()).replace(' ', '-')
    
    env.Append(CCFLAGS=' -DUNAME=\'"'+uname+'"\' -DBUILD_TYPE=\'"'+build_type+'"\' ')

    if static_libs:
        libpath = '#lib'
        env.Append(LINKFLAGS=' -static ')
        if sys.platform != 'cygwin':
            env.Append(LINKFLAGS='-pthread ')
    else:
        libpath = '#linklib'
    env.Replace(LIBPATH=[libpath] + other_libpaths)

    exes = []
    for main_cpp in glob.glob('apps/src/*.cpp'):
        exes.append(env.Program(main_cpp,
                                LIBS=['heart'] + comp_deps['heart'] + other_libs))

    if not compile_only:
        import re
        texttest_result_line = re.compile(r'<H2>.*: 1 tests: 1 (FAILED|succeeded) </H2>')
        def texttest_parse_function(target, source, env):
            """Parse results from texttest to figure out if acceptance tests passed.
            target is a dummy file, since we don't know what we'll output until we're done.
            source is the texttest results file.
            """
            fp = open(str(source[0]))
            fails, succs = 0, 0
            for line in fp:
                m = texttest_result_line.match(line)
                if m:
                    result = m.group(1)
                    if result == 'FAILED':
                        fails += 1
                    else:
                        succs += 1
            fp.close()
            if fails == 0 and succs == 0:
                status = 'unknown'
            elif fails == 0:
                status = 'OK'
            else:
                status = '%d_%d' % (fails, fails+succs)
            to_file_name = build.output_dir + '/AcceptanceTests.' + status + '.0'
            # Remove any old copies of results from this test
            oldfiles = glob.glob(os.path.join(build.output_dir, 'AcceptanceTests.*'))
            for oldfile in oldfiles:
                os.remove(oldfile)
            # Copy results and update summary dependencies
            env.Execute(Copy(to_file_name, str(source[0])))
            if test_summary:
                env.Depends(summary_index, to_file_name)
            return None
        parse_builder = Builder(action=texttest_parse_function)
        env['BUILDERS']['ParseTexttest'] = parse_builder
        # Run acceptance tests
        print "Running acceptance tests", map(str, exes)
        checkout_dir = Dir('#').abspath
        texttest = build.tools['texttest'] + ' -d ' + checkout_dir + '/apps/texttest/chaste'
        texttest_output_dir = env['ENV']['CHASTE_TEST_OUTPUT']+'/texttest_reports/chaste'
        time_eight_hours_ago = time.time() - 8*60*60
        canonical_test_date = time.strftime("%d%b%Y", time.localtime(time_eight_hours_ago))
        todays_file = os.path.join(texttest_output_dir, 'test__' + canonical_test_date + '.html')
        # The next 2 lines make sure the acceptance tests will get run, and the right results stored
        env.Execute(Delete(todays_file))
        env.Execute(Delete(os.path.join(env['ENV']['CHASTE_TEST_OUTPUT'], 'texttest_output')))
        env.Command(todays_file, exes,
                    [texttest + ' -b -c ' + checkout_dir,
                     texttest + ' -c ' + checkout_dir + ' -s batch.GenerateHistoricalReport default'])
        dummy = build.output_dir + '/dummy.texttest'
        env.ParseTexttest(dummy, todays_file)
        env.Depends('apps', [todays_file, dummy])
        if test_summary:
            env.Depends(summary_index, dummy)
        
