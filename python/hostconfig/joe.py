# Configuration for Joe's machines

import os

petsc_2_2_path = None
petsc_2_3_path = '/home/jmpf/petsc-2.3.1-p16/'
petsc_build_name='linux-gnu'
dealii_path = None
metis_path = None
intel_path = '/opt/intel/cc/9.1.039'

other_includepaths = ['/home/jmpf/xsd-2.3.1-i686-linux-gnu/libxsd']
other_libpaths = [os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu/')]
blas_lapack = ['f2clapack', 'f2cblas']
other_libraries = ['boost_serialization', 'xerces-c']

tools = {'mpirun': '/home/jmpf/mpi/bin/mpirun',
         'mpicxx': '/home/jmpf/mpi/bin/mpicxx'}
