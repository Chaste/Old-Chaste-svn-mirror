# Configuration for core Chaste machines

import os

petsc_2_2_path = '../../../petsc-2.2.1/'
petsc_2_3_path = '../../../petsc-2.3.2-p4/'
petsc_build_name = 'linux-gnu'
petsc_build_name_profile = 'linux-gnu-profile'
petsc_build_name_optimized = 'linux-gnu-opt'
petsc_build_name_production = 'linux-intel-opt-mkl'
dealii_path = '../../../deal.II-5.2.0/'
metis_path = '../../../metis-4.0/'
intel_path = '/opt/intel/cce/9.1.039'

other_includepaths = ['../../../xsd-2.3.1-i686-linux-gnu/libxsd', '../../../hdf5/include']
other_libpaths = [os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu/'),  '../../../lib', '/opt/intel/mkl/9.1.023/lib/em64t', '../../../hdf5/lib']
blas_lapack = ['f2clapack', 'f2cblas']
blas_lapack_production = ['mkl_lapack', 'mkl', 'svml']
other_libraries = ['boost_serialization', 'xerces-c', 'z', 'hdf5']

tools = {}
