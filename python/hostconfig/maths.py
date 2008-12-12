# Configuration for Oxford's Maths Institute

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

petsc_2_2_path = None
petsc_2_3_path = '/usr/lib/petsc/'
petsc_build_name = 'linux-gnu-c-opt'
petsc_build_name_optimized = 'linux-gnu-c-opt'
dealii_path = None
metis_path = None
intel_path = '/scratch/chaste/intel/cc/9.1.039/'
icpc='icpc'

other_includepaths = ['/scratch/chaste/xsd-2.3.1-i686-linux-gnu/libxsd', '/scratch/chaste/hdf5/include']
other_libpaths = ['/scratch/chaste/hdf5/lib']
blas_lapack = ['lapack', 'blas-3']
other_libraries = ['boost_serialization', 'xerces-c', 'z', 'hdf5'] 

use_cvode = True
if use_cvode:
    other_includepaths.append('/scratch/chaste/cvode/include')
    other_libpaths.append('/scratch/chaste/cvode/lib')
    other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])
ccflags = ''

tools = {}
