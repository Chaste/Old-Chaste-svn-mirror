# Chaste Build System Scripts

# This module is designed to be imported by both the build scripts and the
# web interface to test results. Given a name representing a build type
# (a valid value of the build argument to scons) it determines what compile
# tools & flags to use, and also how to interpret the status string of a test
# suite.

import os

class BuildType(object):
  """
  Base class for all objects representing a build type.
  Also gives the default build options.
  """
  
  def __init__(self):
    """
    Do any setup.
    Here we set member variables for each method to use.
    """
    self._compiler_type = 'gcc'
    self._cc_flags = ''
    self._link_flags = ''
    self._test_packs = ['Continuous']
    self._revision = ''
    self.build_dir = 'default'
    self._num_processes = 1
    self.using_dealii = False
    self.is_optimised = False
  
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
    return self._cc_flags
  
  def LinkFlags(self):
    """
    Return the linker flags to use, as a string.
    Note that this does not cover library search paths or what to link with.
    """
    return self._link_flags
  
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
    else:
      return status.replace('_', '/') + ' tests failed'

  def EncodeStatus(self, exitCode, logFile):
    """
    Encode the output from a test program as a status string.
    Parses the output looking for a line
    'Failed (\d+) of (\d+) tests?'; if one is found then the
    testsuite failed and the status string is '\1_\2'.
    Otherwise if the output contains as many 'OK!' lines as
    the number of processes running then the test suite is 
    deemed to have passed.
    If neither type of line is found (e.g. due to premature termination)
    then the status is 'Unknown'.
    Return the encoded status.
    """
    status = 'Unknown'
    
    import re
    failed_tests = re.compile('Failed (\d+) of (\d+) tests?')
    ok, ok_count = re.compile('OK!'), 0
    infrastructure_ok = re.compile('Infrastructure test passed ok.')
    
    for line in logFile:
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
    """
    Return the command to be used to run a test suite.
    exefile is the filename of the test executable.
    exeflags are any flags to be passed to the executable.
    The default is just to run the given exectuable.
    """
    return exefile + ' ' + exeflags
    
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
    """
    self.using_dealii = use_dealii
    self.build_dir = 'dealii_' + self.build_dir

  def GetDealiiLibraries(self, dealii_basepath):
    """Return a list of Deal.II libraries to link against.
    
    This method is provided so that optimised builds can use
    the optimised libraries.
    """
    metis_libs = ['metis']
    dealii_libs = ['deal_II_1d', 'deal_II_2d', 'deal_II_3d', 'lac', 'base']
    dealii_petsc = dealii_basepath + 'lib/libpetsc'
    if not self.is_optimised:
      dealii_libs = map(lambda s: s + '.g', dealii_libs)
      dealii_petsc = dealii_petsc + '.g'
    dealii_petsc = dealii_petsc + '.so'
    return dealii_libs + metis_libs


Gcc = BuildType

class GccDebug(Gcc):
  """
  gcc compiler with debug enabled.
  """
  def __init__(self):
    Gcc.__init__(self)
    self._cc_flags += ' -g'
    self.build_dir = 'debug'
    
class Coverage(GccDebug):
  """
  gcc compiler with options to allow for coverage testing.
  """
  def __init__(self):
    GccDebug.__init__(self)
    self._cc_flags += ' -fprofile-arcs -ftest-coverage'
    self._link_flags += ' -fprofile-arcs -ftest-coverage'
    self.build_dir = 'coverage'
    self._num_processes = 2
    self._test_packs.append('Failing')
    #self.UseDealii(True)

  def GetTestRunnerCommand(self, exefile, exeflags=''):
    "Run test on 1 processor then on 2 processors"
    return exefile + ' ' + exeflags + '; mpirun -np ' + str(self._num_processes) + ' ' + exefile + ' ' + exeflags

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
        s = s + status.replace('_', '/') + ' lines uncovered'
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

class Profile(GccDebug):
  """
  gcc compiler with profiling enabled (and optimisation).
  """
  def __init__(self):
    GccDebug.__init__(self)
    self._cc_flags += ' -O3 -pg'
    self._link_flags += ' -pg'
    self._test_packs = ['Profile']
    self.build_dir = 'profile'
  
  def GetTestRunnerCommand(self, exefile, exeflags=''):
    "Run test with a profiler and rename gmon.out"
    return exefile + ' ' + exeflags + ' ; gprof ' + exefile


class GoogleProfile(GccDebug):
  """
  gcc compiler with profiling enabled (and optimisation).
  """
  def __init__(self):
    GccDebug.__init__(self)
    self._cc_flags += ' -O3'
    self._link_flags += ' -lprofiler'
    self._test_packs = ['Profile']
    self.build_dir = 'google_profile'
    import socket
    machine = socket.getfqdn()
    self.output_dir = os.path.join(self.GetTestReportDir(), machine+'.'+self.__class__.__name__)
 
  def ParseGraphFilename(self, filename):
    "Remove the string 'Runner.gif' from the end of a filename, thus returning test_suite name"
    return filename[:-10]

  def GetTestRunnerCommand(self, exefile, exeflags=''):
    "Run test with a profiler and rename gmon.out"
    base = os.path.basename(exefile)
    profileFile = '/tmp/'+base+'.prof'
    return "export HOME='.' ;export CPUPROFILE=" + profileFile + '; ' \
      + exefile + ' ' + exeflags  + ' ; pprof -gif ' \
      + exefile + ' ' + profileFile + ' > ' +self.output_dir+'/'+base+'.gif ; rm ' \
      + profileFile

  def StatusColour(self, status):
    """
    Return a colour string indicating whether the given status string
    represents a 'successful' test suite under this build type.
    """
    prof = False
    if status[-5:] == '_prof':
      prof = True
      status = status[:-5]
    # By default, 'OK' is ok and anything else isn't.
    if status == 'OK':
      if prof:
        return 'orange'
      else:
        return 'green'
    else:
      return 'red'
    
  def DisplayStatus(self, status):
    """
    Return a (more) human readable version of the given status string.
    """
    ret = ''
    if status[-5:] == '_prof':
      ret = 'Profiler failed. '
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
  def __init__(self):
    GccDebug.__init__(self)
    self._test_packs = ['Parallel']
    self._num_processes = 2
  
  def GetTestRunnerCommand(self, exefile, exeflags=''):
    "Run test with a two processor environment"
    return 'mpirun -np ' + str(self._num_processes) + ' ' + exefile + ' ' + exeflags


class Parallel10(Parallel):
  """
  Run using mpi run for tests which run in a parallel environment
  """
  def __init__(self):
    Parallel.__init__(self)
    self._num_processes = 10
    

class MemoryTesting(GccDebug):
  """
  Compile using gcc with debugging turned on, and run tests under valgrind.
  """
  _petsc_flags = "-trmalloc -trdebug -trdump"

  def __init__(self):
    GccDebug.__init__(self)
    #self._cc_flags = self._cc_flags + ' -DPETSC_MEMORY_TRACING'
    self._valgrind_flags = "--tool=memcheck --log-fd=1 --track-fds=yes --leak-check=yes"

  def GetTestRunnerCommand(self, exefile, exeflags=''):
    "Run all tests using valgrind to check for memory leaks."
    return 'valgrind  ' + self._valgrind_flags + ' ' + \
      exefile + ' ' + exeflags + ' ' + self._petsc_flags

  def DisplayStatus(self, status):
    "Return a (more) human readable version of the given status string."
    if status == 'OK':
      return 'No leaks found'
    elif status == 'Unknown':
      return 'Test output unrecognised'
    else:
      return 'Memory leaks found'

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
    lost = re.compile('==\d+==\s+(definitely|indirectly|possibly) lost: (\d+) bytes in (\d+) blocks.')
    petsc = re.compile('\[0]Total space allocated (\d+) bytes')
    uninit = re.compile('==\d+== (Conditional jump or move depends on uninitialised value\(s\)|Use of uninitialised value)')
    open_files = re.compile('==\d+== Open (?:file descriptor|AF_UNIX socket) (?![012])(\d+): (?!(?:/home/bob/eclipse/lockfile|/dev/urandom))(.*)')
    
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
          if int(match.group(3)) > 0:
            status = 'Leaky'
            break
          lineno += 1
          match = lost.match(outputLines[lineno])
        break
        
      m = open_files.match(outputLines[lineno])
      if m:
        # There's a file open that shouldn't be.
        # Descriptors 0, 1 and 2 are ok, as are names /dev/urandom
        # and /home/bob/eclipse/lockfile, and the log file.
        file_name = m.group(2)
        if not file_name.endswith(logFile.name):
          status = 'Openfile'
          break
    else:
      # No leak summary found
      status = 'OK'
    return status


class ParallelMemoryTesting(MemoryTesting, Parallel):
  """
  """
  def __init__(self):
    Parallel.__init__(self)

  def GetTestRunnerCommand(self, exefile, exeflags=''):
    "Run test within a two processor environment"
    return 'mpirun -np ' + str(self._num_processes) + ' -dbg=valgrind ' + \
      exefile + ' ' + exeflags + ' ' + self._petsc_flags

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
  def __init__(self):
    Gcc.__init__(self)
    self._cc_flags = '-O3'
    self.build_dir = 'optimised'
    self.is_optimised = True

class GccOptP4(GccOpt):
  """
  gcc compiler with optimisations for Pentium 4.
  """
  def __init__(self):
    GccOpt.__init__(self)
    self._cc_flags = self._cc_flags+' -march=pentium4 -mmmx -msse -msse2 -mfpmath=sse'
    self.build_dir = 'optimised_P4'
    
class Intel(BuildType):
  "Intel compiler tools."
  def __init__(self):
    BuildType.__init__(self)
    self._compiler_type = 'intel'
    # Turn off some warnings
    self._cc_flags = '-wr470 -wr186'
    self._link_flags = '-static-libcxa'
    self.build_dir = 'intel'
    # Intel compiler uses optimisation by default
    self.is_optimised = True

  def SetReporting(self, vec=1):
    """
    Set the reporting level.
    vec controls the vectoriser report, and is the number to put after
      -vec_report. Default is 1 to indicate vectorised loops; use 3 to
      find out why loops aren't vectorised.
    """
    # Remove any current reporting
    i = self._cc_flags.find('-vec_report')
    if i > -1:
      self._cc_flags = self._cc_flags[:i] + self._cc_flags[i+13:]
    self._cc_flags = self._cc_flags + ' -vec_report' + vec

class IntelNonopt(Intel):
  "Intel compilers with no optimisation."
  def __init__(self):
    Intel.__init__(self)
    self._cc_flags = self._cc_flags + ' -O0 -xK'
    self.build_dir = 'intel_nonopt'
    self.is_optimised = False

class IntelP3(Intel):
  "Intel compilers optimised for Pentium 3."
  def __init__(self):
    Intel.__init__(self)
    self._cc_flags = self._cc_flags + ' -xK -O3 -ip -ipo0 -ipo_obj'
    self._link_flags = self._link_flags + ' -ipo'
    self.build_dir = 'intel_p3'

class IntelP4(Intel):
  "Intel compilers optimised for Pentium 4."
  def __init__(self):
    Intel.__init__(self)
    self._cc_flags = self._cc_flags + ' -xN -O3 -ip -ipo0 -ipo_obj -static'
    self._link_flags = self._link_flags + ' -ipo -lsvml -L/opt/intel_cc_80/lib -static'
    self.build_dir = 'intel_p4'





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
  exec "obj = %s()" % classname
  
  for extra in extras:
    if extra == 'report':
      if issubclass(obj, Intel):
        obj.SetReporting(vec=3)
    elif extra == 'onlytests':
      obj.ClearTestPacks()
    elif extra == 'ndebug':
      obj._cc_flags += ' -DNDEBUG'
      obj.build_dir += '_ndebug'
    elif extra == 'dealii':
      obj.UseDealii(True)
    else:
      # Assume it's a test pack
      obj.AddTestPacks(extra)
  
  return obj
