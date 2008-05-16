# Configuration for core Chaste machines

"""Copyright (C) University of Oxford, 2008

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

petsc_2_2_path = '/opt/petsc-2.2.1-with-mpi/'
petsc_2_3_path = None
petsc_build_name='linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
dealii_path = None
metis_path = None
intel_path = None
icpc='icpc'

other_includepaths = ['/opt/boost/include/boost-1_33_1','/opt/hdf5/include']
other_libpaths = ['/opt/boost/lib/','/opt/hdf5/lib']
blas_lapack = ['lapack', 'blas']
other_libraries = ['boost_serialization-gcc', 'xerces-c', 'z', 'hdf5']

tools = {'mpirun': '/opt/mpi/bin/mpirun',
         'mpicxx': '/opt/mpi/bin/mpicxx'}
