
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


"""Useful functions for use by the build system."""

import os
import sys
import glob
import time
import subprocess

from SCons.Script import Command, Dir, Value
import SCons.Action
import SCons.Tool
import SCons.Script

# Compatability with Python 2.3
try:
    set = set
except NameError:
    import sets
    set = sets.Set

# Possible extensions for source files in Chaste
chaste_source_exts = ['.cpp', '.xsd', '.cellml']

def FindSourceFiles(env, rootDir, ignoreDirs=[], dirsOnly=False, includeRoot=False,
                    source_exts=None):
    """Look for source files under rootDir.
    
    Returns 2 lists: the first of source (.cpp, .xsd) files, and the second
    of the directories in which they may be found.
    
    Optionally:
     * specify ignoreDirs to not search within particular folder names
     * set dirsOnly to True to only find source directories.  In this case
       only a single list is returned
     * set includeRoot to True to include the rootDir in the returned folder list
    """
    source_files = []
    source_dirs = []
    if source_exts is None:
        source_exts = chaste_source_exts
    ignoreDirs.append('.svn')
    if includeRoot:
        source_dirs.append(rootDir)
    for dirpath, dirnames, filenames in os.walk(rootDir):
        for dirname in dirnames[:]:
            if dirname in ignoreDirs:
                dirnames.remove(dirname)
            else:
                source_dirs.append(os.path.join(dirpath, dirname))
        if not dirsOnly:
            for filename in filenames:
                if os.path.splitext(filename)[1] in source_exts:
                    filepath = os.path.join(dirpath, filename)
                    source_files.append(filepath)
    if dirsOnly:
        return source_dirs
    else:
        component = os.path.basename(os.path.dirname(os.path.abspath(rootDir)))
        if component == 'global' and rootDir == 'src':
            # Special-case the version info files.
            file_name = os.path.join('src', 'Version.cpp')
            file_node = env.File(file_name)
            if (not env['UPDATE_CHASTE_PROVENANCE'] and
                os.path.exists(file_node.abspath)):
                # Don't update provenance info - just use the existing file
                version_value = Value(open(file_node.abspath).read())
            else:
                version_value = Value(GetVersionCpp(file_name + '.in', env))
            file_node = env.Command(file_name, [version_value], GenerateCppFromValue)[0]
            source_files.append(file_node)
            # This one just contains the path to Chaste
            source_files.append(env.Command(os.path.join('src', 'ChasteBuildRoot.cpp'),
                                            [Value(GetChasteBuildRootCpp(env))],
                                            GenerateCppFromValue)[0])
        return source_files, source_dirs

def BuildTest(target, source, env):
    """A builder for test executables.

    Takes as a single source the object file compiled for the test,
    and uses its implicit dependencies to work out which other object
    files must be linked against.  For each header file, if it is a
    Chaste header and has a corresponding source .cpp file, then we
    should link with the associated object file.  This analysis is
    recursive - we analyse each object file in the same way.

    It requires a few attributes of the environment:
     * env['CHASTE_COMPONENTS']
       A list of the Chaste components.
     * env['CHASTE_OBJECTS']
       A dictionary mapping cpp file paths (relative to the Chaste root)
       to object file nodes.
     * env['RUNNER_EXE']
       The test runner executable (absolute) path.
     * env['TestBuilder']
       A callable for use in building the test runner executable.
       Should take keyword parameters target (filled by RUNNER_EXE)
       and source (will be given the required object files).
    """
    header_files = set()
    objects = []
    #print sorted(env['CHASTE_OBJECTS'].keys())

    def process(o):
        """Process an object file as described in BuildTest.__doc__"""
        o.scan() # Needed to ensure scons dependencies are set up
        for d in o.implicit:
            hdr = str(d)
            if hdr not in header_files:
                header_files.add(hdr)
                # Is this a Chaste header?
                parts = hdr.split(os.path.sep)
                component = parts[0]
                if component in env['CHASTE_COMPONENTS']:
                    # Does it have a source file?
                    base, ext = os.path.splitext(hdr)
                    if base in ['global/src/Version', 'global/src/ChasteBuildRoot']:
                        # Special cases
                        has_source = True
                        source_filename = base + '.cpp'
                    else:
                        if parts[1] == 'build':
                            # It's a generated source file
                            parts = [component] + parts[3:]
                            base = os.path.join(*parts)[:-len(ext)]
                        for ext in chaste_source_exts:
                            source_filename = base + ext
                            has_source = source_filename in env['CHASTE_OBJECTS']
                            if has_source:
                                break
                    if has_source:
                        # Find the object file(s) and analyse it/them
                        objs = env['CHASTE_OBJECTS'][source_filename]
                        objects.extend(objs)
                        for obj in objs:
                            process(obj)

    for o in source:
        process(o)
    # Build the test itself
    runner = env['RUNNER_EXE']
    #print "Building", runner, "from", pns(source+objects)
    env['TestBuilder'](target=runner, source=source+objects)
    return None

def pns(nodes):
  """Pretty-print nodes for debugging"""
  return map(str, nodes)


def FindTestsToRun(build, BUILD_TARGETS,
                   singleTestSuite, singleTestSuiteDir, allTests,
                   component=None, project=None):
    """Find header files defining tests to run.

    One of component or project must be specified; this says which Chaste
    component or user project to hunt for tests in.

    If singleTestSuite is set, and singleTestSuiteDir is equal to component
    (or project), then just run the requested test.
    If instead allTests is True, then find all tests listed in test packs
    in this component/project.
    Otherwise, if this component (or project) is being built (determined by
    checking BUILD_TARGETS) then search for all tests listed in the test packs
    specified by build.TestPacks().

    Returns an iterable of header file leaf paths (relative to the test folder).
    """
    testfiles = set()
    # Check arguments
    assert component or project
    if component:
        assert project is None
    else:
        component = project
    # Check for a single test
    if singleTestSuite:
        #print singleTestSuite, singleTestSuiteDir
        if singleTestSuiteDir == component:
            testfiles.add(singleTestSuite)
            # Remove any old test output file to force a re-run
            try:
                base = os.path.splitext(singleTestSuite)[0]
                os.remove(base + '.log')
            except OSError:
                pass
    else:
        # Are we building this component/project?
        test_this_comp = False
        root_dir = Dir('#').abspath
        this_comp_targets = ['.', root_dir]
        if not project and component in comp_deps['core']:
            this_comp_targets.append('core')
        if project:
            this_comp_targets.extend(
                [os.path.join('projects', project),
                 os.path.join(root_dir, 'projects', project)])
        else:
            this_comp_targets.extend([component,
                                      os.path.join(root_dir, component)])
        #print map(str, BUILD_TARGETS)
        #print component, project, this_comp_targets
        for targ in BUILD_TARGETS:
            if str(targ).endswith(os.sep):
                # Allow users to specify (e.g.) "global/" as a target
                # (handy for use with tab completion).
                targ = str(targ)[:-len(os.sep)]
            if str(targ) in this_comp_targets:
                test_this_comp = True
                break
        if test_this_comp:
            # Find appropriate test pack files
            packfiles = []
            if allTests:
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
            # Find tests in those test pack files
            for packfile in packfiles:
                try:
                    for testfile in map(lambda s: s.strip(), packfile.readlines()):
                        # Ignore blank lines and repeated tests.
                        if testfile and not testfile in testfiles:
                            testfiles.add(testfile)
                    packfile.close()
                except IOError:
                    pass
    return testfiles


def ExeName(env, exePath):
    """Figure out the real name of an executable.
    
    Given the Linux-style path, this works out what an executable is actually
    called, so that we can run on cygwin.
    """
    pre = env.subst('$PROGPREFIX')
    suf = env.subst('$PROGSUFFIX')
    dirpath = os.path.dirname(exePath)
    name = os.path.basename(exePath)
    return os.path.join(dirpath, pre+name+suf)


def GetVersionCpp(templateFilePath, env):
    """Return the contents of the Version.cpp source file."""
    chaste_root = Dir('#').abspath
    version_file = os.path.join(chaste_root, 'ReleaseVersion.txt')
    if os.path.exists(version_file):
        # Extract just the revision number from the file.
        full_version = open(version_file).read().strip()
        chaste_revision = int(full_version[1+full_version.rfind('.'):])
        wc_modified = False
    else:
        version_pipe = os.popen("svnversion")
        chaste_revision = version_pipe.read().strip()
        if version_pipe.close():
            chaste_revision = 'UINT_MAX'
        else:
            # Extract upper end of range, and store modified flag
            wc_modified = chaste_revision[-1] == 'M'
            if wc_modified:
                chaste_revision = chaste_revision[:-1]
            chaste_revision = int(chaste_revision[1+chaste_revision.rfind(':'):])
    time_format = "%a, %d %b %Y %H:%M:%S +0000"
    build_time = time.strftime(time_format, time.gmtime())
    subst = {'example': '%(example)s',
             'chaste_root': chaste_root,
             'revision': chaste_revision,
             'wc_modified': str(wc_modified).lower(),
             'time_format': time_format,
             'time_size': len(build_time)+1,
             'build_time': build_time,
             'uname': ' '.join(os.uname()),
             'build_type': env['build'].build_type,
             'build_dir': env['build'].build_dir,
             'build_info': env['CHASTE_BUILD_INFO']}
    return open(templateFilePath).read() % subst

def GetChasteBuildRootCpp(env):
    """Return the contents of the ChasteBuildRoot.cpp source file."""
    subst = {'chaste_root': Dir('#').abspath,
             'build_type': env['build'].build_type,
             'build_dir': env['build'].build_dir}
    return """
#include "ChasteBuildRoot.hpp" 

const char* ChasteBuildRootDir() 
{ 
    return "%(chaste_root)s/"; 
}

std::string ChasteComponentBuildDir(const std::string& rComponent)
{
    return std::string(ChasteBuildRootDir()) + rComponent + "/build/%(build_dir)s/";
}

std::string ChasteBuildDirName()
{
    return "%(build_dir)s";
}

std::string ChasteBuildType()
{
    return "%(build_type)s";
}

""" % subst

def _GenerateCppFromValue(env, target, source):
    """An Action to generate a source file from a value node.

    Use like:
    env.Command('global/src/Version.cpp', [Value(GetVersionCpp(templateFilePath, env))], GenerateCppFromValue)
    or:
    env.Command('global/src/ChasteBuildRoot.cpp', [Value(GetChasteBuildRootCpp())], GenerateCppFromValue)
    """
    out = open(target[0].path, "w")
    out.write(source[0].get_contents())
    out.close()
GenerateCppFromValue = SCons.Action.Action(_GenerateCppFromValue, "Generating $TARGET from build information.")

def RecordBuildInfo(env, build_type, static_libs, use_chaste_libs):
    """Record key build information for the provenance system."""
    # TODO: Add library versions used, etc?
    build_info = build_type
    if use_chaste_libs:
        libtype = ['shared', 'static'][static_libs]
        build_info += ', ' + libtype + ' libraries'
    else:
        build_info += ', no Chaste libraries'
    env['CHASTE_BUILD_INFO'] = build_info

def CreateXsdBuilder(build, buildenv):
    """'Builder' for running xsd to generate parser code from an XML schema."""
    # Check if  'xsd' is really CodeSynthesis xsd...
    if not SCons.Script.GetOption('clean'):
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
        # And to assist transitioning to the new builder...
        for path in glob.glob('heart/src/io/ChasteParameters*.?pp'):
            try:
                os.remove(path)
            except OSError:
                pass

    def RunXsd(target, source, env):
        """Action for running XSD."""
        schema_file = str(source[0])
        output_dir = os.path.dirname(target[0].abspath)
        command = ' '.join([build.tools['xsd'], 'cxx-tree',
                            '--generate-serialization',
                            '--output-dir', output_dir,
                            '--hxx-suffix', '.hpp', '--cxx-suffix', '.cpp',
                            '--prologue-file', 'heart/src/io/XsdPrologue.txt',
                            '--epilogue-file', 'heart/src/io/XsdEpilogue.txt',
                            '--namespace-regex', "'X.* $Xchaste::parametersX'",
                            '--namespace-regex', "'X.* https://chaste.comlab.ox.ac.uk/nss/parameters/(.+)Xchaste::parameters::v$1X'",
                            schema_file])
        os.system(command)
    XsdAction = buildenv.Action(RunXsd, "Running xsd on $SOURCE")
    def XsdEmitter(target, source, env):
        hpp = os.path.splitext(str(target[0]))[0] + '.hpp'
        return (target + [hpp], source)
    # Add XSD as a source of .cpp files
    c_file, cxx_file = SCons.Tool.createCFileBuilders(buildenv)
    cxx_file.add_action('.xsd', XsdAction)
    cxx_file.add_emitter('.xsd', XsdEmitter)

def CreatePyCmlBuilder(build, buildenv):
    """Create a builder for running PyCml to generate C++ source code from CellML.
    
    First step just runs with default args.  Later step will be to customize the
    process with a separate options file (e.g. <cellml_file>.opts).
    """
    def IsDynamicSource(source):
        parts = source[0].srcnode().path.split(os.path.sep)
        return (parts[1] == 'dynamic' or
                (parts[1] == 'build' and parts[3] == 'dynamic'))
    def RunPyCml(target, source, env):
        script = os.path.join(Dir('#').abspath, 'python', 'ConvertCellModel.py')
        args = ['--normal', '--opt', '--cvode', '-A',
                '--output-dir', os.path.dirname(target[0].abspath)]
        if IsDynamicSource(source):
            args.remove('--cvode')
            args.append('-y')
        # TODO: Strip CVODE if not being used
        # TODO: Backward Euler if .out file available
        command = [script] + args + [str(source[0])]
        retcode = subprocess.call(command)
        assert retcode == 0, "PyCml execution failed"
    PyCmlAction = buildenv.Action(RunPyCml)
    def PyCmlEmitter(target, source, env):
        base, ext = os.path.splitext(str(target[0]))
        target.extend([base + '.hpp',
                       base + 'Opt.cpp',
                       base + 'Opt.hpp'])
        if not IsDynamicSource(source):
            target.extend([base + 'Cvode.cpp',
                           base + 'Cvode.hpp',
                           base + 'CvodeOpt.cpp',
                           base + 'CvodeOpt.hpp'])
        return (target, source)
    # Add PyCml as a source of .cpp files
    c_file, cxx_file = SCons.Tool.createCFileBuilders(buildenv)
    cxx_file.add_action('.cellml', PyCmlAction)
    cxx_file.add_emitter('.cellml', PyCmlEmitter)


def DoProjectSConscript(projectName, chasteLibsUsed, otherVars):
    """Main logic for a project's SConscript file.
    
    The aim of this method is that a project's SConscript file should be able to be as
    simple as:
        import os
        Import("*")
        project_name = os.path.basename(os.path.dirname(os.path.dirname(os.getcwd())))
        chaste_libs_used = comp_deps['core']
        result = SConsTools.DoProjectSConscript(project_name, chaste_libs_used, globals())
        Return("result")
    """
    # Commonly used variables
    env = otherVars['env']
    use_chaste_libs = otherVars['use_chaste_libs']
    # Note that because we are using SCons' variant dir functionality, and the build
    # dir is created before the SConscript files are executed, that the working dir
    # will be set to <project>/build/<something>.
    curdir = os.getcwd()
    # Look for .cpp files within the project's src folder
    os.chdir('../..') # This is so .o files are built in <project>/build/<something>/
    files, extra_cpppath = FindSourceFiles(env, 'src', includeRoot=True)
    # Look for source files that tests depend on under <project>/test/.
    testsource, test_cpppath = FindSourceFiles(env, 'test', ['data'])
    extra_cpppath.extend(test_cpppath)
    del test_cpppath
    # Move back to the build dir
    os.chdir(curdir)

    # Look for files containing a test suite
    # A list of test suites to run will be found in a test/<name>TestPack.txt
    # file, one per line.
    # Alternatively, a single test suite may have been specified on the command line.
    testfiles = FindTestsToRun(otherVars['build'], otherVars['BUILD_TARGETS'],
                               otherVars['single_test_suite'],
                               otherVars['single_test_suite_dir'],
                               otherVars['all_tests'],
                               project=projectName)

    # Add extra source and test folders to CPPPATH only for this project
    if extra_cpppath:
        newenv = env.Clone()
        newenv.Prepend(CPPPATH=extra_cpppath)
        # Make sure both envs reference the same dict *object*,
        # so updates in one env are reflected in all.
        newenv['CHASTE_OBJECTS'] = env['CHASTE_OBJECTS']
        env = newenv

    # Libraries to link against (TODO: only add project libs if sources exist)
    all_libs = ['test'+projectName, projectName] + chasteLibsUsed + otherVars['other_libs']

    if use_chaste_libs:
        # Build the library for this project
        project_lib = env.StaticLibrary(projectName, files)
        
        # Build the test library for this project
        test_lib = env.StaticLibrary('test'+projectName, testsource)
    else:
        # Build the object files for this project
        project_lib = test_lib = None
        for source_file in files + testsource:
            objs = env.StaticObject(source_file)
            key = os.path.join('projects', projectName, source_file)
            #print projectName, "source", key
            env['CHASTE_OBJECTS'][key] = objs

    # Make test output depend on shared libraries, so if implementation changes
    # then tests are re-run.
    lib_deps = [project_lib, test_lib] # only this project's libraries
    #lib_deps.extend(map(lambda lib: '#lib/lib%s.so' % lib, chasteLibsUsed)) # all Chaste libs used

    # Collect a list of test log files to use as dependencies for the test
    # summary generation
    test_log_files = []

    # Build and run tests of this project
    if not use_chaste_libs:
        env['TestBuilder'] = \
            lambda target, source: env.Program(target, source,
                        LIBS=otherVars['other_libs'],
                        LIBPATH=otherVars['other_libpaths'])
    for testfile in testfiles:
        prefix = os.path.splitext(testfile)[0]
        #print projectName, 'test', prefix
        test_hpp = os.path.join('test', testfile)
        runner_cpp = env.Test(prefix+'Runner.cpp', test_hpp)
        runner_exe = env.File(prefix+'Runner').abspath
        if use_chaste_libs:
            runner_exe = env.Program(runner_exe, runner_cpp,
                                     LIBS = all_libs,
                                     LIBPATH = ['#/lib', '.'] + otherVars['other_libpaths'])
        else:
            runner_obj = env.StaticObject(runner_cpp)
            runner_dummy = runner_exe+'.dummy'
            runner_exe = env.File(ExeName(env, runner_exe))
            env.BuildTest(runner_dummy, runner_obj, RUNNER_EXE=runner_exe)
            env.AlwaysBuild(runner_dummy)
            env.Depends(runner_exe, runner_dummy)
        # Make sure we build the test unless the user says otherwise
        env.Default(runner_exe)
        if not otherVars['compile_only']:
            log_file = env.File(prefix+'.log')
            if use_chaste_libs:
                env.Depends(log_file, lib_deps)
            else:
                env.Depends(log_file, runner_dummy)
            test_log_files.append(log_file)
            env.RunTest(log_file, runner_exe)
            if otherVars['force_test_runs']:
                env.AlwaysBuild(log_file)
    
    return test_log_files

def DoComponentSConscript(component, otherVars):
    """Main logic for a Chaste component's SConscript file.
    
    The aim of this method is that a component's SConscript file should be able to be as
    simple as:
        import os
        Import("*")
        toplevel_dir = os.path.basename(os.path.dirname(os.path.dirname(os.getcwd())))
        result = SConsTools.DoComponentSConscript(toplevel_dir, globals())
        Return("result")
    """
    # Commonly used variables
    env = otherVars['env']
    use_chaste_libs = otherVars['use_chaste_libs']
    # Note that this is executed from within the <component>/build/<something>/ folder
    curdir = os.getcwd()
    # Look for source files within the <component>/src folder
    os.chdir('../..') # This is so .o files are built in <component>/build/<something>/
    files, _ = FindSourceFiles(env, 'src', ignoreDirs=['broken'])
    # Look for source files that tests depend on under test/.
    # We also need to add any subfolders to the CPPPATH, so they are searched
    # for #includes.
    testsource, test_cpppath = FindSourceFiles(env, 'test',
                                               ignoreDirs=['data'],
                                               includeRoot=True)
    # Find any source files that should get made into dynamically loadable modules.
    dyn_source, _ = FindSourceFiles(env, 'dynamic')
    # Move back to the buid dir
    os.chdir(curdir)

    # Look for files containing a test suite
    # A list of test suites to run will be found in a test/<name>TestPack.txt
    # file, one per line.
    # Alternatively, a single test suite may have been specified on the command line.
    testfiles = FindTestsToRun(otherVars['build'], otherVars['BUILD_TARGETS'],
                               otherVars['single_test_suite'],
                               otherVars['single_test_suite_dir'],
                               otherVars['all_tests'],
                               component=component)
    
    # Add test folders to CPPPATH only for this component
    if test_cpppath:
        newenv = env.Clone()
        newenv.Prepend(CPPPATH=test_cpppath)
        # Make sure both envs reference the same dict *object*,
        # so updates in one env are reflected in all.
        newenv['CHASTE_OBJECTS'] = env['CHASTE_OBJECTS']
        env = newenv
    
    # Build any dynamically loadable modules
    dyn_libs = []
    for s in dyn_source:
        # Note: if building direct from CellML, there will be more than 1 target
        dyn_objs = otherVars['dynenv'].SharedObject(source=s)
        for o in dyn_objs:
            so_lib = otherVars['dynenv'].OriginalSharedLibrary(source=o)
            so_dir = os.path.abspath(os.path.join(curdir, '..', '..', os.path.dirname(s)))
            dyn_libs.append(otherVars['dynenv'].Install(so_dir, so_lib))
    
    # Build and install the library for this component
    if use_chaste_libs:
        if otherVars['static_libs']:
            lib = env.StaticLibrary(component, files)
            lib = env.Install('#lib', lib)
            libpath = '#lib'
        else:
            if files:
                lib = env.SharedLibrary(component, files)
            else:
                lib = None
            libpath = '#linklib'
        # Build the test library for this component
        env.StaticLibrary('test'+component, testsource)
    else:
        # Don't build libraries - tests will link against object files directly
        lib = None
        for source_file in files + testsource:
            objs = env.StaticObject(source_file)
            key = os.path.join(component, str(source_file))
            #print component, "source", key
            env['CHASTE_OBJECTS'][key] = objs
    
    # Determine libraries to link against.
    # Note that order does matter!
    if lib:
        chaste_libs = [component] + otherVars['comp_deps'][component]
    else:
        chaste_libs = otherVars['comp_deps'][component]
    all_libs = ['test'+component] + chaste_libs + otherVars['other_libs']
    
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
                        LIBS=otherVars['other_libs'],
                        LIBPATH=otherVars['other_libpaths'])
    
    for testfile in testfiles:
        prefix = os.path.splitext(testfile)[0]
        #print component, 'test', prefix
        test_hpp = os.path.join('test', testfile)
        runner_cpp = env.Test(prefix+'Runner.cpp', test_hpp)
        runner_exe = env.File(prefix+'Runner').abspath
        if use_chaste_libs:
            runner_exe = env.Program(runner_exe, runner_cpp,
                                     LIBS = all_libs,
                                     LIBPATH = [libpath, '.'] + otherVars['other_libpaths'])
        else:
            runner_obj = env.StaticObject(runner_cpp)
            runner_dummy = runner_exe+'.dummy'
            runner_exe = env.File(ExeName(env, runner_exe))
            env.BuildTest(runner_dummy, runner_obj, RUNNER_EXE=runner_exe)
            env.AlwaysBuild(runner_dummy)
            env.Depends(runner_exe, runner_dummy)
        # Make sure we build the test unless the user says otherwise
        env.Default(runner_exe)
        if not otherVars['compile_only']:
            log_file = env.File(prefix+'.log')
            if use_chaste_libs:
                env.Depends(log_file, lib_deps)
            else:
                env.Depends(log_file, runner_dummy)
            if dyn_libs:
                # All tests should depend on dynamically loadable modules, just in case
                env.Depends(log_file, dyn_libs)
            test_log_files.append(log_file)
            env.RunTest(log_file, runner_exe)
            if otherVars['force_test_runs']:
                env.AlwaysBuild(log_file)
    
    return (test_log_files, lib)
