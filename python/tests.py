# Chaste tests display script.

# This program is designed to run from mod_python, and displays results
# from tests run on Chaste in a user-friendly manner.

# Requests to https://comlab2.lsi.ox.ac.uk/tests.py will be directed to
# the index method.
# Requests to https://comlab2.lsi.ox.ac.uk/tests.py/methodname will be
# directed to the appropriate method.

#####################################################################
##                        Configuration.                           ##
#####################################################################

class _conf:
  # This is the directory where test results are stored.
  _tests_dir = '/var/www/chaste_test_data'

  # Base URL for source browsing
  _source_browser_url = 'https://comlab2.lsi.ox.ac.uk/cgi-bin/trac.cgi/browser/'

  # Base URL for this script
  _our_url = '/tests.py'

  # Repository URL or path, username, and password
  #_svn_repos = 'https://comlab2.lsi.ox.ac.uk/svn/chaste/trunk'
  _svn_repos = 'file:///var/svn/chaste/trunk'
  _svn_user  = 'dev'
  _svn_pass  = 'deve10p3r'


import os, time, pysvn

_m = None

#####################################################################
##                          Webpages                               ##
#####################################################################

def index(req):
  """The main test display page.
  
  This displays a summary of the most recent tests.
  """
  _importCodeModule()
  return _m.index(req)


def testsuite(req, type, revision, machine, buildType, testsuite, status):
  """
  Display the results for the given testsuite, by passing the file back
  to the user.
  """
  _importCodeModule()
  return _m.testsuite(req, type, revision, machine, buildType, testsuite, status)

def recent(req, type=''):
  """Display brief summaries of recent builds of the given type.
  
  Returns a string representing part of a webpage.
  """
  _importCodeModule()
  return _m.recent(req, type)

def summary(req, type, revision, machine, buildType):
  """Display a summary of a build.
  
  Returns a string representing part of a webpage.
  """
  _importCodeModule()
  return _m.summary(req, type, revision, machine, buildType)

def buildType(req, buildType, revision=None):
  """
  Display information on the compiler settings, etc. used to build a set
  of tests.
  buildType is the user-friendly name describing these settings, such as
  can be passed to scons build=buildType.
  revision is the code revision of the set of tests, in case the 
  definition of buildType has changed since.
  """
  _importCodeModule()
  return _m.buildType(req, buildType, revision)


#####################################################################
##                    Helper functions.                            ##
#####################################################################

def _importModuleFromSvn(module_name, module_filepath,
                         revision=pysvn.Revision(pysvn.opt_revision_kind.head)):
  """
  Use pysvn and imp to import the requested revision of the given
  module from the repository.
  module_name is the name to give the module.
  module_filepath is the path to the module file within the trunk
  directory of the repository.
  By default import the latest version.
  Return the module object.
  """
  filepath = _conf._svn_repos + module_filepath
  client = _svnClient()  
  module_text = client.cat(filepath, revision)
  return _importCode(module_text, module_name)

def _svnClient():
  "Return a pysvn.Client object for communicating with the svn repository."
  client = pysvn.Client()
  def get_login(realm, username, may_save):
    return True, _conf._svn_user, _conf._svn_pass, True
  client.callback_get_login = get_login
  # We trust our server even though it has a dodgy certificate
  # Might want to make this a bit more secure?
  def ssl_server_trust_prompt(trust_dict):
    return True, trust_dict['failures'], True
  client.callback_ssl_server_trust_prompt = ssl_server_trust_prompt
  return client

def _importCode(code, name, add_to_sys_modules=0):
  """
  Import dynamically generated code as a module. code is the
  object containing the code (a string, a file handle or an
  actual compiled code object, same types as accepted by an
  exec statement). The name is the name to give to the module,
  and the final argument says wheter to add it to sys.modules
  or not. If it is added, a subsequent import statement using
  name will return this module. If it is not added to sys.modules
  import will try to load it in the normal fashion.
  Code from the Python Cookbook.
  
  import foo
  
  is equivalent to
  
  foofile = open("/path/to/foo.py")
  foo = importCode(foofile,"foo",1)
  
  Returns a newly generated module.
  """
  import sys, imp
  
  module = imp.new_module(name)
  
  exec code in module.__dict__
  if add_to_sys_modules:
    sys.modules[name] = module
    
  return module



def _importCodeModule():
  "Import our code from the repository."
  global _m
  if _m is None:
    _m = _importModuleFromSvn('DisplayTests', '/python/DisplayTests.py')

    # Transfer the configuration data
    for item in dir(_conf):
      if item[1] != '_':
        exec '_m.%s = _conf.%s' % (item, item)
