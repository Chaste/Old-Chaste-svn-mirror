# Configuration for Ozzy's Comlab machine

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
petsc_2_3_path = '/home/ozzy/petsc-2.3.2-p10/'
petsc_build_name = 'linux-gnu'
petsc_build_name_profile = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
dealii_path = None
metis_path = None
intel_path = None #'/opt/intel/cce/10.0.025/'
#icpc = 'icpc -gcc-version=413 -I /usr/include/c++/4.1.3/x86_64-linux-gnu/ -I/usr/include/c++/4.1.3/'

ldflags = ' /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf '

other_includepaths = [] 
other_libpaths = ['/usr/lib', '/home/ozzy/lib']
blas_lapack = []
blas_lapack_production = ['mkl_lapack', 'mkl', 'svml']
other_libraries = ['boost_serialization', 'xerces-c', 'hdf5']

use_cvode = False
if use_cvode:
    other_includepaths.append('/home/ozzy/cvode/include')
    other_libpaths.append('home/ozzy/cvode/lib')
    other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])

tools = {'texttest': '/home/chaste/texttest-3.10/source/bin/texttest.py',
         'mpirun': '/home/ozzy/mpi/bin/mpirun',
         'mpicxx': '/home/ozzy/mpi/bin/mpicxx'}
    


