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
build_type = ARGUMENTS.get('build', '')
use_intel = (build_type[:3] == 'icc')
if system_name == 'finarfin':
  mpirun = '/usr/bin/mpirun'
  if use_intel:
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
extra_flags, link_flags = '', ''

if build_type == 'gcc_p4':
  extra_flags, link_flags = ' -O3 -march=pentium4 -mmmx -msse -msse2 -mfpmath=sse ', ''
elif build_type == 'gcc_opt':
  extra_flags, link_flags = ' -O3 ', ''
elif build_type == 'gcc_prof':
  extra_flags, link_flags = ' -pg ', ' -pg ' # gcc profiling
elif use_intel:
  common_flags = ' -static-libcxa -wr470 -wr186 ' # (Turn off some warnings)
  if build_type == 'icc_O0':
    extra_flags = common_flags + '-O0 -xK '
    link_flags = common_flags
  elif build_type == 'icc_p3':
    extra_flags = common_flags + '-xK -O3 -ip -ipo0 -ipo_obj '
    link_flags = common_flags + '-ipo '
  elif build_type == 'icc_p4':
    extra_flags = common_flags + '-xN -O3 -ip -ipo0 -ipo_obj -static '
    link_flags = common_flags + '-ipo -lsvml -L/opt/intel_cc_80/lib -static '
  else:
    extra_flags = common_flags
    link_flags = common_flags
  if ARGUMENTS.get('report', 0):
    extra_flags = extra_flags + '-vec_report3 '


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

SConscript('global/SConscript', build_dir='datawriters/build', duplicate=0)
SConscript('io/SConscript', build_dir='io/build', duplicate=0)
SConscript('ode/SConscript', build_dir='odes/build', duplicate=0)
SConscript('pde/SConscript', build_dir='pdes/build', duplicate=0)
SConscript('coupled/SConscript', build_dir='heart/build', duplicate=0)
