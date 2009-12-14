
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

"""Chaste Build System

This module is designed to be imported by both the build scripts and the
web interface to test results.  Given a name representing a build type
(a valid value of the build argument to scons) it determines what compile
tools & flags to use, and also how to interpret the status string of a test
suite.
"""

import os

class BuildType(object):
    """
    Base class for all objects representing a build type.
    Also gives the default build options.
    """
  
    def __init__(self, buildType):
        """
        Do any setup.
        Here we set member variables for each method to use.
        """
        self.build_type = buildType
        self._compiler_type = 'gcc'
        self._cc_flags = ['-Wall', '-Werror']
        self._link_flags = []
        self._include_flag = ['-isystem']
        self._test_packs = ['Continuous']
        self._revision = ''
        self.build_dir = 'default'
        self._num_processes = 1
        self._hostConfigSettings = {}
        self.using_dealii = False
        self.dealii_debugging = False
        self.is_optimised = False
        self.is_profile = False
        self.is_production = False
        # Where test output will go
        import socket
        machine_fqdn = socket.getfqdn()
        self.output_dir = os.path.join(self.GetTestReportDir(),
                                       machine_fqdn+'.'+buildType)
        # Where tools such as mpirun can be found;
        # by default assume they're on the PATH.
        # The SConstruct file can then override these as appropriate.
        self.tools = {'mpirun': 'mpirun', 'mpicxx': 'mpicxx',
                      'ar': 'ar', 'cxx': 'cxx',
                      'valgrind': 'valgrind',
                      'xsd': 'xsd',
                      'gprof': 'gprof', 'pprof': 'pprof',
                      'rm': 'rm', 'cat': 'cat'}

    def SetNumProcesses(self, np):
        """Set the number of parallel processes to run."""
        assert np > 0, 'Cannot run fewer than 1 process!'
        self._num_processes = np

    def CompilerType(self):
        """
        Return the type of compiler tools to use.
        Currently recognised strings are 'gcc' and 'intel'.
        """
        return self._compiler_type
    
    def CcFlags(self):
        """
        Return the CC flags to use, as a string.
        Note that this does not cover include paths or library search paths.
        """
        return ' '.join(self._cc_flags)
    
    def LinkFlags(self):
        """
        Return the linker flags to use, as a string.
        Note that this does not cover library search paths or what to link with.
        """
        return ' '.join(self._link_flags)

    def IncludeFlag(self):
        """
        Return the flags to use for include files.
        """
        return ' '.join(self._include_flag)

    
    def TestPacks(self):
        """
        Return a list of the test packs to run as part of this build.
        """
        return self._test_packs

    def AddTestPacks(self, *packs):
        """
        Adds each string argument to the list of test packs to be run.
        """
        for pack in packs:
            if not pack in self._test_packs:
                self._test_packs.append(pack)
    
    def ClearTestPacks(self):
        "Empty the list of test packs to be run."
        self._test_packs = []

    def StatusColour(self, status):
        """
        Return a colour string indicating whether the given status string
        represents a 'successful' test suite under this build type.
        """
        # By default, 'OK' is ok and anything else isn't.
        if status == 'OK':
            return 'green'
        elif status == 'MPI':
            return 'orange'
        else:
            return 'red'
        
    def DisplayStatus(self, status):
        """
        Return a (more) human readable version of the given status string.
        """
        if status == 'OK':
            return 'All tests passed'
        elif status == 'Unknown':
            return 'Test output unrecognised'
        elif status == 'MPI':
            return 'MPI semaphore error'
        else:
            return status.replace('_', '/') + ' tests failed (RED)'

    def EncodeStatus(self, exitCode, logFile):
        """Encode the output from a test program as a status string.
        
        Parses the output looking for a line
        'Failed (\d+) of (\d+) tests?'; if one is found then the
        testsuite failed and the status string is '\1_\2'.
        Otherwise if the output contains as many 'OK!' lines as
        the number of processes running then the test suite is 
        deemed to have passed.
        If neither type of line is found (e.g. due to premature termination)
        then the status is 'Unknown'.
        """
        status = 'Unknown'
        
        import re
        failed_tests = re.compile('Failed (\d+) of (\d+) tests?')
        ok, ok_count = re.compile('OK!'), 0
        infrastructure_ok = re.compile('Infrastructure test passed ok.')
        mpi_error = 'semget failed for setnum = '

        first_line = True
        for line in logFile:
            if first_line:
                first_line = False
                if line.find(mpi_error) != -1:
                    status = 'MPI'
                    break
            m = failed_tests.match(line)
            if m:
                status = '%d_%d' % (int(m.group(1)), int(m.group(2)))
                break
            m = ok.match(line)
            if m:
                ok_count += 1
            m = infrastructure_ok.match(line)
            if m:
                ok_count = self._num_processes
                break
        
        if ok_count > 0 and status == 'Unknown':
            # All tests passed on all processes
            status = 'OK'
        return status

    def SetRevision(self, revision):
        """
        Set the subversion revision number of the code that is being built.
        revision will be '' if we don't know or don't care.
        """
        self._revision = revision
    
    def GetTestReportDir(self):
        """
        Return the base directory in which to store the output from all
        the tests. Files with names that include status info will be
        saved in a subdirectory named 'machine.buildtype'.

        Note: various places assume this includes a trailing slash.
        Note2: the builder script also has this path hardcoded.
        """
        return 'testoutput/'
    
    def GetTestRunnerCommand(self, exefile, exeflags=''):
        """Return the command to be used to run a test suite.
        
        exefile is the filename of the test executable.
        exeflags are any flags to be passed to the executable.
        
        The default behaviour is just to run the given exectuable.
        If self._num_processes > 1, then mpirun is used to run the
        executable in parallel.
        """
        cmd = exefile + ' ' + exeflags
        if self._num_processes > 1:
            cmd = ' '.join([self.tools['mpirun'], '-np',
                            str(self._num_processes), cmd])
        return cmd
        
    def ResultsFileName(self, dir, testsuite, status, runtime):
        """
        Return the path to a results file.
        dir is the directory in which files should be put.
        testsuite is the name of the test suite giving these results.
        status is an encoded status string summarising the results of the test.
        runtime is the time taken for the test to complete, in seconds (as a floating point no.).
        """
        leafname = testsuite + '.' + status
        if runtime >= 0:
            leafname = leafname + '.' + str(int(runtime))
        pathname = os.path.join(dir, leafname)
        return pathname
    
    def GetInfoFromResultsFileName(self, leafname):
        """
        Extract the metadata held within the name of a results file.
        This returns a dictionary, with keys 'testsuite', 'status' and 'runtime'.
        testsuite is the name of the test suite.
        status is the encoded status string.
        runtime is the run time for the test suite in seconds.
        """
        # Components are separated by '.'
        i2 = leafname.rfind('.')
        i1 = leafname.rfind('.', 0, i2)
        if i1 == -1:
            # No runtime info available
            runtime = -1
            i1, i2 = i2, len(leafname)
        else:
            runtime = int(leafname[i2+1:])
        return {'testsuite': leafname[:i1],
                'status': leafname[i1+1:i2],
                'runtime': runtime}

    def UseDealii(self, use_dealii):
        """Set whether this build should link against Deal.II.
        
        Several things need to change if we do:
         * The build_dir - alter this directly
         * The libraries linked against, and search paths for them
                - set a flag that SConstruct can check
         * The default test pack - use DealiiContinuous
        """
        self.using_dealii = use_dealii
        self.build_dir = 'dealii_' + self.build_dir
        if 'Continuous' in self._test_packs:
            self._test_packs[self._test_packs.index('Continuous')] = 'DealiiContinuous'

    def GetDealiiLibraries(self, dealii_basepath):
        """Return a list of Deal.II libraries to link against.
        
        This method is provided so that optimised builds can use
        the optimised libraries.
        """
        metis_libs = ['metis']
        dealii_libs = ['deal_II_1d', 'deal_II_2d', 'deal_II_3d', 'lac', 'base']
        dealii_petsc = dealii_basepath + 'lib/libpetsc'
        if self.dealii_debugging:
            dealii_libs = map(lambda s: s + '.g', dealii_libs)
            dealii_petsc = dealii_petsc + '.g'
        dealii_petsc = dealii_petsc + '.so'
        return dealii_libs + metis_libs

    def GetSignature(self):
        """Return the signature of this build object for SCons.

        This determines when tests need to be re-run because a different build
        is being used.
        """
        return '*'.join([self.build_type, self.build_dir])

    def SetHostConfig(self, configString):
        """Parse hostconfig settings from a build type option string.
        
        This method extracts prefered versions of libraries from a string
        with format "libraryName1=version1,libraryName2=version2", where
        version numbers are given as "1-2-3".  Currently recognised library
        names are 'petsc', 'boost', 'hdf5' and 'xsd', but support from the
        machine-specific hostconfig file is needed too.  Prefered versions
        can be retrieved using GetPreferedVersions.
        """
        items = configString.split(',')
        config = {'petsc': '2.3',
                  'boost': '1.34',
                  'hdf5': '1.6',
                  'xsd': '3.2'}
        for item in items:
            key, val = item.split('=')
            config[key] = val.replace('-', '.')
        self._hostConfigSettings = config
        return

    def GetPreferedVersions(self):
        """Get the prefered versions of libraries parsed with SetHostConfig.
        
        Returns a dictionary mapping library name to version string, e.g.
        {'petsc': '2.3', 'boost': '1.34'}
        """
        return self._hostConfigSettings

Gcc = BuildType

class GccDebug(Gcc):
    """
    gcc compiler with debug enabled.
    """
    def __init__(self, *args, **kwargs):
        Gcc.__init__(self, *args, **kwargs)
        self._cc_flags.append('-g')
        self.build_dir = 'debug'


class Coverage(GccDebug):
    """
    gcc compiler with options to allow for coverage testing.
    """
    def __init__(self, *args, **kwargs):
        GccDebug.__init__(self, *args, **kwargs)
        self._cc_flags.extend(['-fprofile-arcs', '-ftest-coverage'])
        self._link_flags.extend(['-fprofile-arcs', '-ftest-coverage'])
        self.build_dir = 'coverage'
        self._num_processes = 2
        #self._test_packs.append('Failing')
        #self.UseDealii(True)

    def UseDealii(self, use_dealii):
        """Set whether this build should link against Deal.II.

        Extends the base method so both continuous test packs are run.
        """
        super(Coverage, self).UseDealii(use_dealii)
        self._test_packs.append('Continuous')

    def GetTestRunnerCommand(self, exefile, exeflags=''):
        "Run test on 1 processor then on 2 processors"
        return exefile + ' ' + exeflags + '; ' + self.tools['mpirun'] + \
                ' -np ' + str(self._num_processes) + ' ' + exefile + ' ' + exeflags

    def DisplayStatus(self, status):
        """
        Return a (more) human readable version of the given status string.
        """
        if status == 'OK':
            s = 'All lines covered'
        elif status == 'Unknown':
            s = 'Output unrecognised'
        else:
            if status.startswith('ignore_'):
                s = 'Unterminated COVERAGE_IGNORE block. '
                status = status[7:]
            else:
                s = ''
            if status.startswith('warn_'):
                s = s + status[5:].replace('_', '/') + " lines 'spuriously' uncovered"
            else:
                s = s + status.replace('_', '/') + ' lines marked uncovered (RED)'
        return s

    def StatusColour(self, status):
        """
        Return a colour string indicating whether the given status string
        represents a 'successful' test suite under this build type.
        """
        # 'OK' is green, warnings are orange, otherwise red
        if status == 'OK':
            return 'green'
        elif status.startswith('warn_'):
            return 'orange'
        else:
            return 'red'

class DoxygenCoverage(GccDebug):
    """Check for documentation coverage/problems."""
    def DisplayStatus(self, status):
        """
        Return a (more) human readable version of the given status string.
        """
        if status == 'OK':
            s = 'No Doxygen problems found'
        elif status == 'Unknown':
            s = 'Output unrecognised'
        else:
            s = ''
            if status.startswith('warn_'):
                s = s + status[5:].split('_')[0] + " Doxygen warnings"
            else:
                s = s + status.split('_')[0] + ' Doxygen errors (RED)'
        return s
    
    def GetTestRunnerCommand(self, exefile, exeflags=''):
        "We don't actually run any tests in this build..."
        return ''

class CovTool(Coverage):
    """Coverage testing using the covtool software."""
    def __init__(self, *args, **kwargs):
        # Don't call Coverage.__init__ as it adds flags we don't want
        GccDebug.__init__(self, *args, **kwargs)
        self.build_dir = 'covtool'
        self._num_processes = 2
        self._cc_flags.append('-v')
        # Note: the following can be overridden by hostconfig
        self.tools['cov++'] = 'cov++'
        self.tools['covannotate'] = 'covannotate.exe'
        self.tools['covmerge'] = 'covmerge.exe'

    def UseCovTool(self, other_includepaths, other_libs):
        """Switch to using the covtool toolchain."""
#        covcmd = self.tools['cov++'] + ' -EXT .cpp .c++ -DIAG -VER' + \
#            ''.join(map(lambda p: ' -skip '+p,
#                        other_includepaths + ['/usr'])) # Don't instrument other people's code
#        self.tools['mpicxx'] += ' -CC="%s"' % covcmd
        self.tools['mpicxx'] = ' '.join([self.tools['cov++'], '-CMD',
                                         self.tools['mpicxx'], self.tools['mpicxx'],
                                         '-EXT .cpp .c++'])
        self.tools['mpicxx'] += ' -DIAG -VER' # for debugging
        # Don't instrument other people's code
        self.tools['mpicxx'] += ''.join(map(lambda p: ' -skip '+p,
                                            other_includepaths + ['/usr']))
        # Must be last
        other_libs.append('covtoolhelper')


class Profile(GccDebug):
    """
    gcc compiler with profiling enabled (and optimisation).
    Uses -O2 rather than -O3 since we don't want inlining when profiling.
    """
    def __init__(self, *args, **kwargs):
        GccDebug.__init__(self, *args, **kwargs)
        self._cc_flags.extend(['-O2', '-pg']) 
        self._link_flags.append('-pg')
        self._test_packs = ['Profile']
        self.build_dir = 'profile'
        self.is_profile = True
    
    def GetTestRunnerCommand(self, exefile, exeflags=''):
        "Run test with a profiler and rename gmon.out"
        return exefile + ' ' + exeflags + ' ; ' + self.tools['gprof'] + ' ' + exefile

    def SetNumProcesses(self, np):
        """Can't run profiling in parallel (yet)."""
        raise ValueError("The profiling builds cannot be run in parallel.")

class LineProfile(Profile):
    def __init__(self, *args, **kwargs):
        Profile.__init__(self, *args, **kwargs)
        self.build_dir = 'line_profile'

    def GetTestRunnerCommand(self, exefile, exeflags=''):
        "Run test with a profiler and rename gmon.out"
        return exefile + ' ' + exeflags + ' ; ' + self.tools['gprof'] + ' -l ' + exefile

class GoogleProfile(GccDebug):
    """
    gcc compiler with profiling enabled (and optimisation).
    """
    def __init__(self, *args, **kwargs):
        GccDebug.__init__(self, *args, **kwargs)
        self._cc_flags.append('-O3')
        self._link_flags.append('-lprofiler')
        self._test_packs = ['Profile']
        self.build_dir = 'google_profile'
        self.is_profile = True
 
    def ParseGraphFilename(self, filename):
        "Remove the string 'Runner.gif' from the end of a filename, thus returning test_suite name"
        return filename[:-10]

    def GetTestRunnerCommand(self, exefile, exeflags=''):
        "Run test with a profiler and rename gmon.out"
        base = os.path.basename(exefile)
        profileFile = '/tmp/'+base+'.prof'
        return "export HOME='.'; export CPUPROFILE=" + profileFile + '; ' \
            + exefile + ' ' + exeflags    + ' ; ' \
            + self.tools['pprof'] + ' -gif ' + exefile + ' ' + profileFile \
            + ' > ' + self.output_dir+'/'+base+'.gif ; ' \
            + self.tools['rm'] + ' ' + profileFile
    
    def SetNumProcesses(self, np):
        """Can't run profiling in parallel (yet)."""
        raise ValueError("The profiling builds cannot be run in parallel.")

    def StatusColour(self, status):
        """
        Return a colour string indicating whether the given status string
        represents a 'successful' test suite under this build type.
        """
        prof = False
        if status[-5:] == '_prof':
            prof = True
            status = status[:-5]
        base_col = super(GoogleProfile, self).StatusColour(status)
        if prof and base_col == 'green':
            return 'orange'
        else:
            return base_col
        
    def DisplayStatus(self, status):
        """
        Return a (more) human readable version of the given status string.
        """
        ret = ''
        if status[-5:] == '_prof':
            ret = 'Profiler failed.  (RED)'
            status = status[:-5]
        return ret + super(GoogleProfile, self).DisplayStatus(status)

    def EncodeStatus(self, exitCode, logFile):
        """
        Encode the output from a test program as a status string.
        If the exit code is zero then all tests passed, and the status
        is 'OK'. Otherwise the output must be parsed looking for a line
        'Failed (\d+) of (\d+) tests?' and the status string is '\1_\2'.
        Return the encoded status.
        """
        status = super(GoogleProfile, self).EncodeStatus(exitCode, logFile)
        if exitCode:
            status = status + '_prof'
        return status


class Parallel(GccDebug):
    """
    Run using mpi run for tests which run in a parallel environment
    """
    def __init__(self, *args, **kwargs):
        GccDebug.__init__(self, *args, **kwargs)
        self._test_packs = ['Parallel']
        self._num_processes = 2
    
    def GetTestRunnerCommand(self, exefile, exeflags=''):
        "Run test with a two processor environment"
        return self.tools['mpirun'] + ' -np ' + str(self._num_processes) \
            + ' ' + exefile + ' ' + exeflags

class Parallel3(Parallel):
    """
    Run using mpi run for tests which run in a parallel environment
    """
    def __init__(self, *args, **kwargs):
        Parallel.__init__(self, *args, **kwargs)
        self._num_processes = 3

class Parallel4(Parallel):
    """
    Run using mpi run for tests which run in a parallel environment
    """
    def __init__(self, *args, **kwargs):
        Parallel.__init__(self, *args, **kwargs)
        self._num_processes = 4

class Parallel10(Parallel):
    """
    Run using mpi run for tests which run in a parallel environment
    """
    def __init__(self, *args, **kwargs):
        Parallel.__init__(self, *args, **kwargs)
        self._num_processes = 10
        

class MemoryTesting(GccDebug):
    """
    Compile using gcc with debugging turned on, and run tests under valgrind.
    """
    _petsc_flags = "-malloc_debug -malloc_dump -memory_info"
    _valgrind_flags = "--tool=memcheck --log-file=%s --track-fds=yes --leak-check=yes --num-callers=50 --suppressions=chaste.supp"
#    _valgrind_flags +=" --gen-suppressions=all"
    def __init__(self, *args, **kwargs):
        GccDebug.__init__(self, *args, **kwargs)
        #self._cc_flags.append('-DPETSC_MEMORY_TRACING')
        #self.build_dir += '_mem'

    def GetTestRunnerCommand(self, exefile, exeflags=''):
        "Run all tests using valgrind to check for memory leaks."
        test_suite = os.path.basename(exefile)
        log_prefix = self.GetTestReportDir() + test_suite
        cmd = ' '.join([self.tools['valgrind'], self._valgrind_flags % log_prefix,
                                        exefile, exeflags, self._petsc_flags,
                                        ';', self.tools['cat'], log_prefix + '*',
                                        ';', self.tools['rm'], log_prefix + '*'])
        return cmd
    
    def SetNumProcesses(self, np):
        """Can't run profiling in parallel (yet)."""
        raise ValueError("Use ParallelMemoryTesting to run memory tests in parallel.")
    
    def StatusColour(self, status):
        """
        Return a colour string indicating whether the given status string
        represents a 'successful' test suite under this build type.
        """
        if status == 'OK':
            return 'green'
        elif status == 'Warn':
            return 'orange'
        else:
            return 'red'
        
    def DisplayStatus(self, status):
        "Return a (more) human readable version of the given status string."
        if status == 'OK':
            return 'No leaks found'
        elif status == 'Unknown':
            return 'Test output unrecognised'
        elif status == 'Warn':
            return 'Possible leak found'
        else:
            return 'Memory leaks found (RED)'

    def EncodeStatus(self, exitCode, logFile, outputLines=None):
        """
        Encode the output from a test program as a status string.
        The output from valgrind needs to be parsed to check for a leak summary.
        If one is found the status is 'Leaky', otherwise 'OK'.
        Return the encoded status.
        """
        status = 'Unknown'
        
        # Regexps to check for
        import re
        invalid = re.compile('==\d+== Invalid ')
        glibc = re.compile('__libc_freeres')
        leaks = re.compile('==\d+== LEAK SUMMARY:')
        lost = re.compile('==\d+==\s+(definitely|indirectly|possibly) lost: ([0-9,]+) bytes in ([0-9,]+) blocks.')
        petsc = re.compile('\[0]Total space allocated (\d+) bytes')
        uninit = re.compile('==\d+== (Conditional jump or move depends on uninitialised value\(s\)|Use of uninitialised value)')
        open_files = re.compile('==(\d+)== Open (?:file descriptor|AF_UNIX socket) (?![012])(\d+): (?!(?:/home/bob/eclipse/lockfile|/dev/urandom))(.*)')
        
        if outputLines is None:
            outputLines = logFile.readlines()
        for lineno in range(len(outputLines)):
            m = petsc.match(outputLines[lineno])
            if m and int(m.group(1)) > 0:
                # PETSc Vec or Mat allocated and not destroyed
                status = 'Leaky'
                break
                
            m = uninit.match(outputLines[lineno])
            if m:
                # Uninitialised values problem
                status = 'Uninit'
                break
        
            m = invalid.match(outputLines[lineno])
            if m:
                # Invalid read/write/free()/etc. found. This is bad, unless it's glibc's fault.
                match = glibc.search(outputLines[lineno+3])
                if not match:
                    status = 'Leaky'
                    break
                    
            m = leaks.match(outputLines[lineno])
            if m:
                # Check we have really lost some memory
                # (i.e. ignore 'still reachable' memory)
                status = 'OK'
                lineno += 1
                match = lost.match(outputLines[lineno])
                while match:
                    blocks = match.group(3).replace(',', '')
                    if int(blocks) > 0:
                        # Hack for chaste-bob
                        bytes = match.group(2).replace(',', '')
                        if match.group(1) == 'indirectly' and \
                           int(bytes) == 240 and \
                           int(blocks) == 10:
                            status = 'Warn'
                        #HDF5 Test giving unexpected indirect loss
                        elif match.group(1) == 'indirectly' and \
                           int(bytes) == 256 and \
                           int(blocks) == 14:
                            status = 'Warn'   
                        else:
                            status = 'Leaky'
                        break
                    lineno += 1
                    match = lost.match(outputLines[lineno])
                break
                
            m = open_files.match(outputLines[lineno])
            if m:
                # There's a file open that shouldn't be.
                # Descriptors 0, 1 and 2 are ok, as are names /dev/urandom
                # and /home/bob/eclipse/lockfile, and the log files.
                # All these OK files are inherited from the parent process.
                if not outputLines[lineno+1].strip().endswith("<inherited from parent>"):
                    status = 'Openfile'
                    break
        else:
            # No leak summary found
            status = 'OK'
        return status


class ParallelMemoryTesting(MemoryTesting, Parallel):
    """
    """
    def __init__(self, *args, **kwargs):
        Parallel.__init__(self, *args, **kwargs)

    def GetTestRunnerCommand(self, exefile, exeflags=''):
        "Run test within a two processor environment"
        cmd = self.tools['mpirun'] + ' -np ' + str(self._num_processes) + ' ' + \
                MemoryTesting.GetTestRunnerCommand(self, exefile, exeflags)
        return cmd

    def EncodeStatus(self, exitCode, logFile):
        """
        Encode the output from a test program as a status string.
        The output is sorted by process ID so that checking context in the valgrind
        output (from a single process) works as expected.
        The output from valgrind needs to be parsed to check for a leak summary.
        If one is found the status is 'Leaky', otherwise 'OK'.
        Return the encoded status.
        """
        # First stably sort the output by process id
        import re
        pid = re.compile('==(\d+)==')
        def cmp(l1, l2):
            m1, m2 = pid.match(l1), pid.match(l2)
            if m1:
                pid1 = int(m1.group(1))
            else:
                pid1 = 0
            if m2:
                pid2 = int(m2.group(1))
            else:
                pid2 = 0
            if pid1 == pid2: return 0
            elif pid1 < pid2: return -1
            else: return 1

        output_lines = logFile.readlines()
        output_lines.sort(cmp)
        
        # Now use the parsing from the superclass
        return MemoryTesting.EncodeStatus(self, exitCode, logFile, outputLines=output_lines)


class GccOpt(Gcc):
    """
    gcc compiler with some optimisations enabled.
    """
    def __init__(self, *args, **kwargs):
        Gcc.__init__(self, *args, **kwargs)
        self._cc_flags = ['-O3']
        self.build_dir = 'optimised'
        self.is_optimised = True

class GccOptP4(GccOpt):
    """
    gcc compiler with optimisations for Pentium 4.
    """
    def __init__(self, *args, **kwargs):
        GccOpt.__init__(self, *args, **kwargs)
        self._cc_flags.extend(['-march=pentium4', '-mmmx', '-msse', '-msse2', '-mfpmath=sse'])
        self.build_dir = 'optimised_P4'


class Intel(BuildType):
    """Intel compiler tools."""
    def __init__(self, *args, **kwargs):
        BuildType.__init__(self, *args, **kwargs)
        self._compiler_type = 'intel'
        # Turn off some warnings, and report warnings as errors
        #self._cc_flags = ['-wr470', '-wr186', '-pc64'] #Emulates a 64-bit (not 80-bit) FPU
        self._cc_flags = ['-Werror']
        self._link_flags = ['-static-libcxa']
        self.build_dir = 'intel'
        # Intel compiler uses optimisation by default
        self.is_optimised = True

    def SetReporting(self, vec=1):
        """Set the reporting level.
        
        vec controls the vectoriser report, and is the number to put after
            -vec_report. Default is 1 to indicate vectorised loops; use 3 to
            find out why loops aren't vectorised.
        """
        # Remove any current reporting first
        for i, flag in enumerate(self._cc_flags):
            if flag.startswith('-vec_report'):
                del self._cc_flags[i]
                break
        self._cc_flags.append('-vec_report' + str(vec))

class GccPower(Gcc):
  "GNU compiler on IBM POWER architecture."
  def __init__(self, *args, **kwargs):
    Gcc.__init__(self, *args, **kwargs)
    self.tools['mpicxx'] = 'mpCC'
    self.tools['mpirun'] = ''
    self._cc_flags.append('-compiler g++')
    self._cc_flags.append('-m64')
    self._link_flags.append('-compiler g++')
    self._link_flags.append('-m64')
    self.build_dir = 'gccpower'

class GccPowerDebug(GccPower):
  "GNU compiler on IBM POWER with debugging."
  def __init__(self, *args, **kwargs):
    GccPower.__init__(self, *args, **kwargs)
    self._cc_flags.append('-g')
    self._link_flags.append('-g')
    self.build_dir = 'gccpowerdebug'

class CrayGcc(BuildType):
  "Cray XT platform."
  def __init__(self, *args, **kwargs):
    Gcc.__init__(self, *args, **kwargs)
    self.tools['mpicxx'] = 'CC'
    self._cc_flags.append('-DMPICH_IGNORE_CXX_SEEK')
    self._cc_flags.append('-g')
    self.build_dir = 'craygcc'

class Vacpp(BuildType):
  "IBM Visual Age C++ compiler"
  def __init__(self, *args, **kwargs):
    BuildType.__init__(self, *args, **kwargs)
    self.tools['mpicxx'] = 'mpCC'
    self.tools['mpirun'] = 'poe'
    self._cc_flags = ['-q64 -qhalt=w']
    self._link_flags = ['-q64']
    self._include_flag = ['-I']
    self.build_dir = 'vacpp'

class VacppOpt(Vacpp):
  "Optmized IBM build"
  def __init__(self, *args, **kwargs):
    Vacpp.__init__(self, *args, **kwargs)
    self._cc_flags.append('-qarch=auto -qstrict -qhot -O3')
    self.build_dir = 'vacppopt'

class Fle(BuildType):
  "Intel compiler tools on FLE cluster."
  def __init__(self, *args, **kwargs):
    BuildType.__init__(self, *args, **kwargs)
    self._compiler_type = 'intel'
    # Turn off some warnings
    self._cc_flags = ['-i-dynamic', '-wr470', '-wr186']
    self._cc_flags.extend(['-O3', '-xW'])
    self._cc_flags.extend(['-DNDEBUG'])
    self._link_flags = ['-static-libgcc']
    self.build_dir = 'fle'
    # Intel compiler uses optimisation by default
    self.is_optimised = True

  def GetTestRunnerCommand(self, exefile, exeflags=''):
    "Run test with a single processor environment"
    return self.tools['mpirun'] + ' -machinefile /home/southern/.mpihosts' + ' -np ' \
             + str(self._num_processes) + ' ' + exefile + ' ' + exeflags

class FleDebug(Fle):
    "Intel compilers with debugging enabled on FLE cluster."
    def __init__(self, *args, **kwargs):
        Fle.__init__(self, *args, **kwargs)
        self._cc_flags.append('-g')
        self.build_dir = 'fle_debug'

class FleProfile(Fle):
  "Intel compilers with no optimisation on FLE cluster."
  def __init__(self, *args, **kwargs):
    Fle.__init__(self, *args, **kwargs)
    self._cc_flags.extend(['-p', '-g'])
    self._link_flags.extend(['-p', '-g'])
    self.build_dir = 'fle_profile'

class FleItcProfile(Fle):
  "Intel compilers with no optimisation on FLE cluster."
  def __init__(self, *args, **kwargs):
    Fle.__init__(self, *args, **kwargs)
    self._cc_flags.extend(['-DITC'])
    self._link_flags.extend(['-DITC', '-L/opt/intel/ict/3.0/itac/7.0/itac/lib_mpich', '-lVT', '-ldwarf', '-lelf', '-lnsl', '-lm', '-ldl', '-lpthread'])
    self.build_dir = 'fle_itc_profile'

class FleNonopt(Fle):
  "Intel compilers with no optimisation on FLE cluster."
  def __init__(self, *args, **kwargs):
    Fle.__init__(self, *args, **kwargs)
    self._cc_flags = ['-i-dynamic', '-wr470', '-wr186', '-O0']
    self.build_dir = 'fle_nonopt'
    self.is_optimised = False
				    
class FleMemoryTesting(FleDebug):
    """
    Compile using intel compilers with debugging turned on, and run tests under valgrind on FLE cluster.
    """
    _petsc_flags = "-malloc_debug -malloc_dump -memory_info"
    _valgrind_flags = "--tool=memcheck --log-file=%s --track-fds=yes --leak-check=yes --num-callers=50 --suppressions=chaste.supp --suppressions=fle.supp"

    def __init__(self, *args, **kwargs):
        FleDebug.__init__(self, *args, **kwargs)
        #self._cc_flags.append('-DPETSC_MEMORY_TRACING')
        #self.build_dir += '_mem'

    def GetTestRunnerCommand(self, exefile, exeflags=''):
        "Run all tests using valgrind to check for memory leaks."
        test_suite = os.path.basename(exefile)
        log_prefix = self.GetTestReportDir() + test_suite
        cmd = ' '.join([self.tools['valgrind'], self._valgrind_flags % log_prefix,
                                        exefile, exeflags, self._petsc_flags,
                                        ';', self.tools['cat'], log_prefix + '*',
                                        ';', self.tools['rm'], log_prefix + '*'])
        return cmd
    
    def SetNumProcesses(self, np):
        """Can't run profiling in parallel (yet)."""
        raise ValueError("Use ParallelMemoryTesting to run memory tests in parallel.")
    
    def StatusColour(self, status):
        """
        Return a colour string indicating whether the given status string
        represents a 'successful' test suite under this build type.
        """
        if status == 'OK':
            return 'green'
        elif status == 'Warn':
            return 'orange'
        else:
            return 'red'
        
    def DisplayStatus(self, status):
        "Return a (more) human readable version of the given status string."
        if status == 'OK':
            return 'No leaks found'
        elif status == 'Unknown':
            return 'Test output unrecognised'
        elif status == 'Warn':
            return 'Possible leak found'
        else:
            return 'Memory leaks found (RED)'

    def EncodeStatus(self, exitCode, logFile, outputLines=None):
        """
        Encode the output from a test program as a status string.
        The output from valgrind needs to be parsed to check for a leak summary.
        If one is found the status is 'Leaky', otherwise 'OK'.
        Return the encoded status.
        """
        status = 'Unknown'
        
        # Regexps to check for
        import re
        invalid = re.compile('==\d+== Invalid ')
        glibc = re.compile('__libc_freeres')
        leaks = re.compile('==\d+== LEAK SUMMARY:')
        lost = re.compile('==\d+==\s+(definitely|indirectly|possibly) lost: ([0-9,]+) bytes in ([0-9,]+) blocks.')
        petsc = re.compile('\[0]Total space allocated (\d+) bytes')
        uninit = re.compile('==\d+== (Conditional jump or move depends on uninitialised value\(s\)|Use of uninitialised value)')
        open_files = re.compile('==(\d+)== Open (?:file descriptor|AF_UNIX socket) (?![012])(\d+): (?!(?:/home/bob/eclipse/lockfile|/dev/urandom))(.*)')
        
        if outputLines is None:
            outputLines = logFile.readlines()
        for lineno in range(len(outputLines)):
            m = petsc.match(outputLines[lineno])
            if m and int(m.group(1)) > 0:
                # PETSc Vec or Mat allocated and not destroyed
                status = 'Leaky'
                break
                
            m = uninit.match(outputLines[lineno])
            if m:
                # Uninitialised values problem
                status = 'Uninit'
                break
        
            m = invalid.match(outputLines[lineno])
            if m:
                # Invalid read/write/free()/etc. found. This is bad, unless it's glibc's fault.
                match = glibc.search(outputLines[lineno+3])
                if not match:
                    status = 'Leaky'
                    break
                    
            m = leaks.match(outputLines[lineno])
            if m:
                # Check we have really lost some memory
                # (i.e. ignore 'still reachable' memory)
                status = 'OK'
                lineno += 1
                match = lost.match(outputLines[lineno])
                while match:
                    blocks = match.group(3).replace(',', '')
                    if int(blocks) > 0:
                        # Hack for chaste-bob
                        bytes = match.group(2).replace(',', '')
                        if match.group(1) == 'indirectly' and \
                           int(bytes) == 240 and \
                           int(blocks) == 10:
                            status = 'Warn'
                        else:
                            status = 'Leaky'
                        break
                    lineno += 1
                    match = lost.match(outputLines[lineno])
                break
                
            m = open_files.match(outputLines[lineno])
            if m:
                # There's a file open that shouldn't be.
                # Descriptors 0, 1 and 2 are ok, as are names /dev/urandom
                # and /home/bob/eclipse/lockfile, and the log files.
                # All these OK files are inherited from the parent process.
                if not outputLines[lineno+1].strip().endswith("<inherited from parent>"):
                    status = 'Openfile'
                    break
        else:
            # No leak summary found
            status = 'OK'
        return status

class IntelNonopt(Intel):
    """Intel compilers with no optimisation."""
    def __init__(self, *args, **kwargs):
        Intel.__init__(self, *args, **kwargs)
        self._cc_flags.extend(['-O0', '-xK'])
        self.build_dir = 'intel_nonopt'
        self.is_optimised = False

class IntelP3(Intel):
    """Intel compilers optimised for Pentium 3."""
    def __init__(self, *args, **kwargs):
        Intel.__init__(self, *args, **kwargs)
        self._cc_flags.extend(['-xK', '-O3', '-ip', '-ipo0', '-ipo_obj'])
        self._link_flags.append('-ipo')
        self.build_dir = 'intel_p3'

class IntelP4(Intel):
    """Intel compilers optimised for Pentium 4."""
    def __init__(self, *args, **kwargs):
        Intel.__init__(self, *args, **kwargs)
        # -xN = P4 + SSE2 while -xW = P4 (suitable for AMD CPUs)
        self._cc_flags.extend(['-xW', '-O3', '-ip', '-ipo1'])
        self._link_flags.extend(['-ipo1'])
        self.build_dir = 'intel_p4'

class IntelProduction(IntelP4):
    """Intel Production version optimised for Pentium 4."""
    def __init__(self, *args, **kwargs):
        IntelP4.__init__(self, *args, **kwargs)
        self.build_dir = 'intel_production'
        self._cc_flags.append('-DNDEBUG')
        self.is_production = True

class Vtune(IntelProduction):
    """Production build with debug symbols for vtune analyser."""
    def __init__(self, *args, **kwargs):
        super(Vtune, self).__init__(*args, **kwargs)
        self.build_dir = 'intel_vtune'
        self._cc_flags.append('-g')

class IntelProductionParallel4(IntelProduction):
    """Intel production build, run tests in parallel on 4 nodes"""
    def __init__(self, *args, **kwargs):
        IntelProduction.__init__(self, *args, **kwargs)
        self._test_packs = ['Parallel']
        self._num_processes = 4
    
    def GetTestRunnerCommand(self, exefile, exeflags=''):
        "Run test with a 4 processor environment"
        return self.tools['mpirun'] + ' -np ' + str(self._num_processes) \
            + ' ' + exefile + ' ' + exeflags

class StyleCheck(GccDebug):
    """Check the code against Effective C++ style guidelines."""
    def __init__(self, *args, **kwargs):
        GccDebug.__init__(self, *args, **kwargs)
        self._cc_flags = ['-Weffc++']
        self.build_dir = 'style_check'
        self._test_packs.extend(['Failing', 'Profile', 'Nightly'])
        
    def GetTestRunnerCommand(self, exefile, exeflags=''):
        """This build shouldn't be used to run tests."""
        return ""



# Define mappings between arguments on the command line and BuildType objects.
def GetBuildType(buildType):
    """
    Given a string representing a build type, create and return an instance of
    the appropriate BuildType subclass.
    Components of the string are separated by '_'. The first component is the
    basic BuildType, and further components can customise that.
    """
    parts = buildType.split('_')
    classname = parts[0]
    extras = parts[1:]
    
    if classname == '' or classname == 'default':
        # Default build type
        classname = 'GccDebug'
    elif classname == 'failing':
        # Check failing tests
        classname = 'GccDebug'
        extras = ['onlytests', 'Failing'] + extras
    exec "obj = %s('%s')" % (classname, buildType)
    
    for extra in extras:
        if extra == 'report':
            if isinstance(obj, Intel):
                obj.SetReporting(vec=3)
        elif extra == 'onlytests':
            obj.ClearTestPacks()
        elif extra == 'traceksp':
            obj._cc_flags.append('-DTRACE_KSP')
        elif extra == 'ndebug':
            obj._cc_flags.append('-DNDEBUG')
            obj.build_dir += '_ndebug'
        elif extra == 'fpe':
            obj._cc_flags.append('-DTEST_FOR_FPE')
            obj.build_dir += '_fpe'
        elif extra == 'dealii':
            obj.UseDealii(True)
        elif extra == 'debug':
            obj.dealii_debugging = True
        elif extra == 'warn':
            try:
                obj._cc_flags.remove('-Werror')
                obj.build_dir += '_warn'
            except ValueError:
                pass
        elif extra.startswith('hostconfig'):
            obj.SetHostConfig(extra[11:])
            obj.build_dir += '_' + extra
        else:
            try:
                np = int(extra)
                # If it's an integer, assume it sets the number of
                # parallel processes to run
                obj.SetNumProcesses(np)
            except ValueError:
                # Assume it's a test pack
                obj.AddTestPacks(extra)
    
    return obj
