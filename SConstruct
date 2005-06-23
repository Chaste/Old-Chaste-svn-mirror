# Controlling scons build script for Chaste

import sys
sys.path.append('python')
import BuildTypes

build_type = ARGUMENTS.get('build', 'default')
build = BuildTypes.GetBuildType(build_type)
build.SetRevision(ARGUMENTS.get('revision', ''))
Export('build', 'build_type')

# Specify test_summary=1 to scons to generate a summary html page
test_summary = ARGUMENTS.get('test_summary', 0)
Export('test_summary')

# Specify system_name=finarfin to scons to change default paths
system_name = ARGUMENTS.get('system_name', '')

## PETSc library paths
if system_name == 'finarfin':
  # Finarfin (Debian sarge)
  petsc_base = '/home/jonc/work/dtc/courses/IB/petsc-2.2.1/'
  petsc_inc = '-I'+petsc_base+'include '
  petsc_bmake = '-I'+petsc_base+'bmake/linux-gnu '
  petsc_mpi = '-I'+petsc_base+'include/mpiuni '
  petsc_incs = petsc_inc+petsc_bmake+petsc_mpi
  
  petsc_libpath = petsc_base+'lib/libg_c++/linux-gnu/'
else:
  # DTC (default)
  petsc_base = '../../petsc-2.2.1-with-mpi/'
  petsc_inc = '-I'+petsc_base+'include '
  petsc_bmake = '-I'+petsc_base+'bmake/linux-gnu '
  petsc_mpi = '-I'+petsc_base+'include/mpiuni '
  petsc_incs = petsc_inc+petsc_bmake+petsc_mpi
  
  petsc_libpath = '#'+petsc_base+'lib/libg_c++/linux-gnu/'

Export("petsc_base", "petsc_inc", "petsc_bmake", "petsc_mpi", "petsc_incs", "petsc_libpath")


## C++ build tools & MPI runner
if system_name == 'finarfin':
  mpirun = '/usr/bin/mpirun'
  if build.CompilerType() == 'intel':
    # Use intel compiler
    mpicxx = '/usr/bin/mpicxx -CC=icpc'
    cxx    = '/opt/intel_cc_80/bin/icpc'
    ar     = '/opt/intel_cc_80/bin/xiar'
  else:
    # Use gcc
    mpicxx = '/usr/bin/mpicxx'
    cxx    = '/usr/bin/g++'
    ar     = '/usr/bin/ar'
else:
  mpicxx = '../../mpi/bin/mpicxx'
  mpirun = '../../mpi/bin/mpirun'
  cxx = '/usr/bin/g++'
  ar = '/usr/bin/ar'

Export("mpicxx", "mpirun", "cxx", "ar")


## Any extra CCFLAGS and LINKFLAGS
extra_flags = build.CcFlags()
link_flags  = build.LinkFlags()

Export("extra_flags", "link_flags")

# Search path for #includes
import glob, os
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

SConscript('maths/SConscript', build_dir='maths/build', duplicate=0)
SConscript('mesh/SConscript', build_dir='mesh/build', duplicate=0)
SConscript('global/SConscript', build_dir='global/build', duplicate=0)
SConscript('io/SConscript', build_dir='io/build', duplicate=0)
SConscript('ode/SConscript', build_dir='ode/build', duplicate=0)
SConscript('pde/SConscript', build_dir='pde/build', duplicate=0)
SConscript('coupled/SConscript', build_dir='coupled/build', duplicate=0)

#SConscript('global/SConscript', build_dir='build')
#SConscript('io/SConscript', build_dir='build')
#SConscript('ode/SConscript', build_dir='build')
#SConscript('pde/SConscript', build_dir='build')
#SConscript('coupled/SConscript', build_dir='build')
