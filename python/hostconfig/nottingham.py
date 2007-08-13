# Configuration for core Chaste machines

petsc_2_2_path = '/opt/petsc-2.2.1-with-mpi/'
petsc_2_3_path = None
petsc_build_name='linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
dealii_path = None
metis_path = None
intel_path = None

other_includepaths = ['/opt/boost/include/boost-1_33_1']
other_libpaths = ['/opt/boost/lib/']
blas_lapack = ['lapack', 'blas']
other_libraries = ['boost_serialization-gcc', 'xerces-c']

tools = {'mpirun': '/opt/mpi/bin/mpirun',
         'mpicxx': '/opt/mpi/bin/mpicxx'}
