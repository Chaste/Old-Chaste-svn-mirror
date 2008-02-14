# Configuration for core Chaste machines

petsc_2_2_path = '/opt/petsc-2.2.1-with-mpi/'
petsc_2_3_path = None
petsc_build_name='linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
dealii_path = None
metis_path = None
intel_path = None

other_includepaths = ['/opt/boost/include/boost-1_33_1','/opt/hdf5/include']
other_libpaths = ['/opt/boost/lib/','/opt/hdf5/lib']
blas_lapack = ['lapack', 'blas']
other_libraries = ['boost_serialization-gcc', 'xerces-c', 'z', 'hdf5']

tools = {'mpirun': '/opt/mpi/bin/mpirun',
         'mpicxx': '/opt/mpi/bin/mpicxx'}
