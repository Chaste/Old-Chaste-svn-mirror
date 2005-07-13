# Chaste Build System Scripts

# This module is designed to be imported by both the build scripts and the
# web interface to test results. Given a name representing a build type
# (a valid value of the build argument to scons) it determines what compile
# tools & flags to use, and also how to interpret the status string of a test
# suite.

class BuildType:
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

  def IsGoodStatus(self, status):
    """
    Check the given status string to see if it represents a 'successful'
    test suite under this build type. Return True if so.
    """
    # By default, 'OK' is ok and anything else isn't.
    return status == 'OK'
    
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

  def EncodeStatus(self, exitCode, outputLines):
    """
    Encode the output from a test program as a status string.
    If the exit code is zero then all tests passed, and the status
    is 'OK'. Otherwise the output must be parsed looking for a line
    'Failed (\d+) of (\d+) tests?' and the status string is '\1_\2'.
    Return the encoded status.
    """
    status = 'Unknown'
    if exitCode:
      # At least one test failed. Find out how many by parsing the output.
      import re
      failed_tests = re.compile('Failed (\d+) of (\d+) tests?')
      for line in outputLines:
        m = failed_tests.match(line)
        if m:
          status = '%d_%d' % (int(m.group(1)), int(m.group(2)))
          break
    else:
      # All tests passed
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
    """
    return 'testoutput/'
  
  def GetTestRunnerCommand(self, exefile):
    """
    Return the command to be used to run a test suite.
    exefile is the filename of the test executable.
    The default is just to run the given exectuable.
    """
    return exefile


class GccDebug(BuildType):
  """
  gcc compiler with debug enabled.
  """
  def __init__(self):
    BuildType.__init__(self)
    self._cc_flags = '-g'

class MemoryTesting(BuildType):
  """
  Compile using gcc with debugging turned on, and run tests under valgrind.
  """
  def GetTestRunnerCommand(self, exefile):
    "Run all tests using valgrind to check for memory leaks."
    return 'valgrind --tool=memcheck --log-fd=1 --track-fds=yes --leak-check=full ' + exefile

  def DisplayStatus(self, status):
    "Return a (more) human readable version of the given status string."
    if status == 'OK':
      return 'No leaks found'
    elif status == 'Unknown':
      return 'Test output unrecognised'
    else:
      return 'Memory leaks found'

  def EncodeStatus(self, exitCode, outputLines):
    """
    Encode the output from a test program as a status string.
    The output from valgrind needs to be parsed to check for a leak summary.
    If one is found the status is 'Leaky', otherwise 'OK'.
    Return the encoded status.
    """
    status = 'Unknown'
    import re
    leaks = re.compile('==\d+== LEAK SUMMARY:')
    for line in outputLines:
      m = leaks.match(line)
      if m:
        status = 'Leaky'
        break
    else:
      # No leak summary found
      status = 'OK'
    return status

class GccOpt(BuildType):
  """
  gcc compiler with some optimisations enabled.
  """
  def __init__(self):
    BuildType.__init__(self)
    self._cc_flags = '-O3'

class GccOptP4(GccOpt):
  """
  gcc compiler with optimisations for Pentium 4.
  """
  def __init__(self):
    GccOpt.__init__(self)
    self._cc_flags = self._cc_flags+' -march=pentium4 -mmx -msse -msse2 -mfpmath=sse'
    
class GccProfiled(BuildType):
  """
  gcc compiler with profiling.
  """
  def _init__(self):
    BuildType.__init__(self)
    self._cc_flags = '-pg'
    self._link_flags = '-pg'

class Intel(BuildType):
  "Intel compiler tools."
  def __init__(self):
    BuildType.__init__(self)
    self._compiler_type = 'intel'
    # Turn off some warnings
    self._cc_flags = '-wr470 -wr186'
    self._link_flags = '-static-libcxa'

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

class IntelP3(Intel):
  "Intel compilers optimised for Pentium 3."
  def __init__(self):
    Intel.__init__(self)
    self._cc_flags = self._cc_flags + ' -xK -O3 -ip -ipo0 -ipo_obj'
    self._link_flags = self._link_flags + ' -ipo'

class IntelP4(Intel):
  "Intel compilers optimised for Pentium 4."
  def __init__(self):
    Intel.__init__(self)
    self._cc_flags = self._cc_flags + ' -xN -O3 -ip -ipo0 -ipo_obj -static'
    self._link_flags = self._link_flags + ' -ipo -lsvml -L/opt/intel_cc_80/lib -static'





# Define mappings between arguments on the command line and BuildType objects.
def GetBuildType(buildType):
  """
  Given a string representing a build type, create and return an instance of
  the appropriate BuildType subclass.
  Components of the string are separated by ':'. The first component is the
  basic BuildType, and further components can customise that.
  """
  parts = buildType.split(':')
  classname = parts[0]
  extras = parts[1:]
  
  if classname == '' or classname == 'default':
    classname = 'BuildType'
  exec "obj = %s()" % classname
  
  for extra in extras:
    if extra == 'report':
      if issubclass(obj, Intel):
        obj.SetReporting(vec=3)
    elif extra == 'notests':
      obj.ClearTestPacks()
    else:
      # Assume it's a test pack
      obj.AddTestPacks(extra)
  
  return obj
