# Configuration for Joe's machines

import os

petsc_2_2_path = None
petsc_2_3_path = '/home/jmpf/petsc-2.3.1-p16/'
petsc_build_name = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
petsc_build_name_production = 'linux-intel-opt-mkl'
dealii_path = None
metis_path = None
intel_path = '/opt/intel/cc/9.1.039/'

other_includepaths = ['/home/jmpf/xsd-2.3.1-i686-linux-gnu/libxsd']
other_libpaths = [os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu/'), '/opt/intel/mkl/9.1.023/lib/32']
blas_lapack = ['f2clapack', 'f2cblas']
blas_lapack_production = ['mkl_lapack', 'mkl']
other_libraries = ['boost_serialization', 'xerces-c']

tools = {'mpirun': '/home/jmpf/mpi/bin/mpirun',
         'mpicxx': '/home/jmpf/mpi/bin/mpicxx'}
