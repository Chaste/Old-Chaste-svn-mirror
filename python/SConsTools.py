
"""Copyright (C) University of Oxford, 2008

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
import glob
import SCons.Defaults

# Compatability with Python 2.3
try:
    set = set
except NameError:
    import sets
    set = sets.Set

def IsTemplateCpp(filepath): 
    """Check if the .cpp file defines a templated class. 
     
    We assume this is the case if one of the first 2 lines starts '#ifndef '.
    
    TODO: Has this been silently broken by the copyright blocks?
    """ 
    fp = open(filepath) 
    line = fp.next() 
    template = line.startswith('#ifndef ') 
    if not template: 
        line = fp.next() 
        template = line.startswith('#ifndef ') 
    fp.close() 
    return template

def FindSourceFiles(rootDir, ignoreDirs=[], dirsOnly=False, includeRoot=False):
    """Look for source files under rootDir.
    
    Returns 2 lists: the first of source (.cpp) files, and the second
    of the directories in which they may be found.
    
    Optionally:
     * specify ignoreDirs to not search within particular folder names
     * set dirsOnly to True to only find source directories.  In this case
       only a single list is returned
     * set includeRoot to True to include the rootDir in the returned folder list
    """
    source_files = []
    source_dirs = []
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
                filepath = os.path.join(dirpath, filename)
                if filename[-4:] == '.cpp' and not IsTemplateCpp(filepath):
                    source_files.append(filepath)
    if dirsOnly:
        return source_dirs
    else:
        return source_files, source_dirs

def BuildTest(target, source, env):
    """A builder for test executables.

    Takes as a single source the object file compiled for the test,
    and uses its implicit dependencies to work out which other object
    files must be linked against.  For each header file, if it is a
    Chaste header and has a corresponding 'real' source file (as
    determined by IsTemplateCpp), then we should link with the
    associated object file.  This analysis is recursive - we analyse
    each object file in the same way.

    It requires a few attributes of the environment:
     * env['CHASTE_COMPONENTS']
       A list of the Chaste components.
     * env['CHASTE_OBJECTS']
       A dictionary mapping cpp file paths (relative to the Chaste root)
       to object file nodes.
     * env['RUNNER_EXE']
       The test runner executable path.
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
                component = hdr.split(os.path.sep)[0]
                if component in env['CHASTE_COMPONENTS']:
                    #print "Header", hdr
                    # Does it have a source file?
                    base, ext = os.path.splitext(hdr)
                    cpp_file = base + '.cpp'
                    if os.path.exists(cpp_file) and \
                           not IsTemplateCpp(cpp_file):
                        # Find the object file and analyse it
                        obj = env['CHASTE_OBJECTS'][cpp_file]
                        objects.append(obj)
                        process(obj)

    for o in source:
        process(o)
    # Build the test itself
    runner = env['RUNNER_EXE']
    print "Building", runner, "from", pns(source+objects)
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

    Returns an iterable of header file leaf names.
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
        if singleTestSuiteDir == component:
            testfiles.add(singleTestSuite)
            # Remove any old test output file to force a re-run
            try:
                os.remove(singleTestSuite[:-4] + '.log')
            except OSError:
                pass
    else:
        # Are we building this component/project?
        test_this_comp = False
        root_dir = SCons.Defaults.DefaultEnvironment().Dir('#').abspath
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
