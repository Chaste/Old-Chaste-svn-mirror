# Configuration for Oxford's Maths Institute

petsc_2_2_path = None
petsc_2_3_path = '/usr/lib/petsc/'
petsc_build_name = 'linux-gnu-c-opt'
petsc_build_name_optimized = 'linux-gnu-c-opt'
dealii_path = None
metis_path = None
intel_path = '/scratch/chaste/intel/cc/9.1.039/'

other_includepaths = ['/scratch/chaste/xsd-2.3.1-i686-linux-gnu/libxsd', '/scratch/chaste/hdf5/include']
other_libpaths = ['/scratch/chaste/hdf5/lib']
blas_lapack = ['lapack', 'blas-3']
other_libraries = ['boost_serialization', 'xerces-c', 'z', 'hdf5'] 

ccflags = ''

tools = {}
