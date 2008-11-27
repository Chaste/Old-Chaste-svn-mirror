# Configuration for cluster3 at FLE

import os

petsc_2_3_path = '/opt/petsc/'
petsc_build_name = 'linux-x86_64'
petsc_build_name_optimized = 'linux-x86_64'
dealii_path = None
metis_path = None
intel_path = '/opt/intel/fce/9.1.040/'
icpc = 'icpc'

other_includepaths = ['/home/southern/boost/include/boost-1_34',
                      '/opt/hypre/include',
                      '/share/apps/opt/opt/include',
                      '/home/southern/xsd-2.3.1-i686-linux-gnu/libxsd',
                      '/home/southern/xerces-c-src_2_7_0/include',
		      '/opt/infinipath-mpi/include' ]
other_libpaths = ['/opt/intel/mkl/9.0/lib/em64t',
                  '/usr/X11R6/lib64',
                  '/opt/hypre/lib',
                  '/home/southern/boost/lib',
                  '/home/southern/lib',
		  '/opt/infinipath-mpi/lib64' ]
blas_lapack = ['mkl_lapack64', 'mkl_em64t']
other_libraries = ['mkl_em64t', 'mkl', 'X11', 'HYPRE', 'boost_serialization-gcc34-1_34',
                     'xerces-c', 'svml', 'imf', 'irc']

tools = {'mpirun': '/opt/infinipath-mpi/bin/mpirun',
         'mpicxx': '/opt/infinipath-mpi/bin/mpicxx -CC=icpc'}
