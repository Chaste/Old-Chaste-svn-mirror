# Configuration for core Chaste machines

import os

petsc_2_2_path = '../../../petsc-2.2.1/'
petsc_2_3_path = '../../../petsc-2.3.2-p4/'
petsc_build_name = 'linux-gnu'
petsc_build_name_profile = 'linux-gnu-profile'
petsc_build_name_optimized = 'linux-intel-opt-mkl'
dealii_path = '../../../deal.II/'
metis_path = '../../../metis-4.0/'
intel_path = '/opt/intel/cc/9.1.039'

other_includepaths = ['../../../xsd-2.3.1-i686-linux-gnu/libxsd']
other_libpaths = ['/opt/intel/mkl/9.1.023/lib/32', '/opt/intel/cc/9.1.039/lib/']
blas_lapack = ['mkl_lapack', 'mkl']
other_libraries = ['boost_serialization', 'xerces-c']


tools = {}
