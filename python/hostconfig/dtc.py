# Configuration for core Chaste machines

import os

petsc_2_3_path = '/usr/apps/chaste-plugins/petsc-2.3.2-p4/'
petsc_build_name = 'linux-gnu'

#other_includepaths = ['../../../xsd-2.3.1-i686-linux-gnu/libxsd', '/usr/apps/chaste-plugins/hdf5/include']
other_includepaths = ['/usr/apps/chaste-plugins/hdf5/include']
other_libpaths = [os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu/'),
                  '/usr/apps/chaste-plugins/hdf5/lib']
blas_lapack = ['f2clapack', 'f2cblas']
other_libraries = ['boost_serialization', 'xerces-c', 'z', 'hdf5']

tools = {'mpicxx': '/usr/apps/chaste-plugins/mpi/bin/mpicxx',
         'mpirun': '/usr/apps/chaste-plugins/mpi/bin/mpirun'}

#These are unused - possibly set erroneously
petsc_2_2_path = '/usr/petsc-2.2.1/'
petsc_build_name_profile = 'linux-gnu-profile'
petsc_build_name_optimized = 'linux-gnu-opt'
petsc_build_name_production = 'linux-intel-opt-mkl'
dealii_path = '../../../deal.II/'
metis_path = '../../../metis-4.0/'
intel_path = '/opt/intel/cc/9.1.039'
blas_lapack_production = ['mkl_lapack', 'mkl', 'svml']
#tools = {'texttest': '/home/chaste/texttest-3.10/source/bin/texttest.py'}
