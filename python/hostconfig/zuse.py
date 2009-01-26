# Configuration for Zuse

"""Copyright (C) University of Oxford, 2005-2009

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.
"""

petsc_2_2_path = '/home/zuse/system/software/petsc-2.2.1/'
petsc_2_3_path = None
petsc_build_name = 'linux-mpich-gnu-mkl'
petsc_build_name_optimized = 'linux-mpich-gnu-mkl'
dealii_path = None
metis_path = None
intel_path = '/opt/intel/cc/9.1.039/lib'
icpc='icpc'

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
