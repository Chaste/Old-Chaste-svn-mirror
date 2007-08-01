# Configuration for Zuse

petsc_2_2_path = '/home/zuse/system/software/petsc-2.2.1/'
petsc_2_3_path = None
dealii_path = None
metis_path = None
intel_path = '/opt/intel/cc/9.1.039/lib'

other_includepaths = ['/home/zuse/system/joe/xerces-c-suse_80_AMD_64-gcc_32/include',
                      '/home/zuse/system/joe/xsd-2.3.1-i686-linux-gnu/libxsd',
                      '/home/zuse/system/joe/boost_v1_34/include/boost-1_34']
other_libpaths = ['/opt/intel/mkl/8.0/lib/em64t/',
                  '/home/zuse/system/joe/xerces-c-suse_80_AMD_64-gcc_32/lib',
                  '/home/zuse/system/joe/boost_v1_34/lib']
blas_lapack = []
other_libraries = ['boost_serialization-gcc32', 'xerces-c']

tools = {'mpicxx': '/home/zuse/system/software/mpich-gcc/bin/mpicxx',
         'mpirun': '/home/zuse/system/software/mpich-gcc/bin/mpirun'}
