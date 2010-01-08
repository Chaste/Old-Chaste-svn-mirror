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

# SCons build script for core Chaste components.

import glob
import os

Import("*")

# Note that this script is executed from within the <component>/build/<something>/ folder
curdir = os.getcwd()

# Get our top-level directory
toplevel_dir = os.path.basename(os.path.dirname(os.path.dirname(curdir)))

#print curdir, toplevel_dir

# Look for .cpp files within the src folder
os.chdir('../..') # This is so .o files are built in `toplevel_dir'/build/<something>/
files, _ = SConsTools.FindSourceFiles(env, 'src')

# Look for source files that tests depend on under test/.
# We also need to add any subfolders to the CPPPATH, so they are searched
# for #includes.
testsource, test_cpppath = SConsTools.FindSourceFiles(env, 'test',
                                                      ignoreDirs=['data'],
                                                      includeRoot=True)

os.chdir(curdir)

# Look for files containing a test suite
# A list of test suites to run will be found in a test/<name>TestPack.txt
# file, one per line.
# Alternatively, a single test suite may have been specified on the command
# line.
testfiles = SConsTools.FindTestsToRun(
    build, BUILD_TARGETS,
    single_test_suite, single_test_suite_dir, all_tests,
    component=toplevel_dir)

#print test_cpppath, testsource
#print files, testfiles, testsource

# Add test folders to CPPPATH only for this component
if test_cpppath:
    newenv = env.Copy()
    newenv.Prepend(CPPPATH=test_cpppath)
    # Make sure both envs reference the same dict *object*,
    # so updates in one env are reflected in all.
    newenv['CHASTE_OBJECTS'] = env['CHASTE_OBJECTS']
    env = newenv

# Build and install the library for this component
if use_chaste_libs:
    if static_libs:
        lib = env.Library(toplevel_dir, files)
        lib = env.Install('#lib', lib)
        libpath = '#lib'
    else:
        if files:
            lib = env.SharedLibrary(toplevel_dir, files)
        else:
            lib = None
        libpath = '#linklib'
    # Build the test library for this component
    env.Library('test'+toplevel_dir, testsource)
else:
    # Don't build libraries - tests will link against object files directly
    lib = None
    for source_file in files + testsource:
        obj = env.StaticObject(source_file)
        key = os.path.join(toplevel_dir, str(source_file))
        #print toplevel_dir, "source", key
        env['CHASTE_OBJECTS'][key] = obj[0]
    

# Determine libraries to link against.
# Note that order does matter!
if lib:
	chaste_libs = [toplevel_dir] + comp_deps[toplevel_dir]
else:
	chaste_libs = comp_deps[toplevel_dir]
all_libs = ['test'+toplevel_dir] + chaste_libs + other_libs



# Make test output depend on shared libraries, so if implementation changes
# then tests are re-run.  Choose which line according to taste.
#lib_deps = map(lambda lib: '#lib/lib%s.so' % lib, chaste_libs) # all libs
lib_deps = lib # only this lib
#linklib_deps = map(lambda lib: '#linklib/lib%s.so' % lib, chaste_libs)

# Collect a list of test log files to use as dependencies for the test
# summary generation
test_log_files = []

# Build and run tests of this component
if not use_chaste_libs:
    env['TestBuilder'] = \
        lambda target, source: env.Program(target, source,
                    LIBS=other_libs,
                    LIBPATH=other_libpaths)

for testfile in testfiles:
    prefix = os.path.splitext(testfile)[0]
    #print toplevel_dir, 'test', prefix
    test_hpp = os.path.join('test', testfile)
    runner_cpp = env.Test(prefix+'Runner.cpp', test_hpp)
    runner_exe = File(prefix+'Runner').abspath
    if use_chaste_libs:
        runner_exe = env.Program(runner_exe, runner_cpp,
                                 LIBS = all_libs,
                                 LIBPATH = [libpath, '.'] + other_libpaths)
    else:
        runner_obj = env.StaticObject(runner_cpp)
        runner_dummy = runner_exe+'.dummy'
        runner_exe = File(SConsTools.ExeName(env, runner_exe))
        env.BuildTest(runner_dummy, runner_obj, RUNNER_EXE=runner_exe)
        env.AlwaysBuild(runner_dummy)
        env.Depends(runner_exe, runner_dummy)
    # Make sure we build the test unless the user says otherwise
    Default(runner_exe)
    if not compile_only:
        log_file = env.File(prefix+'.log')
        if use_chaste_libs:
            env.Depends(log_file, lib_deps)
        else:
            env.Depends(log_file, runner_dummy)
        test_log_files.append(log_file)
        env.RunTest(log_file, runner_exe)
        if force_test_runs:
            env.AlwaysBuild(log_file)

return_value = (test_log_files, lib)
Return("return_value")
