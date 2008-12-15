# Configuration for Gary's machine in physiology

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

import os

petsc_2_2_path = None
petsc_2_3_path = '/home/pc-gary/wolf1768/petsc-2.3.3-p15'
petsc_build_name='linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
petsc_build_name_production = 'linux-intel-opt-mkl'
dealii_path = None
metis_path = None
intel_path = None
icpc = 'icpc -gcc-version=431 -I /usr/include/c++/4.3/x86_64-suse-linux/ -I/usr/include/c++/4.3/'

other_includepaths = ['/home/pc-gary/wolf1768/include',
                      '/home/pc-gary/wolf1768/xsd-2.3.1-i686-linux-gnu/libxsd']

other_libpaths = [os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu/'),
                  '/home/pc-gary/wolf1768/lib']

blas_lapack = ['f2clapack', 'f2cblas']  # For some reason they have to be in this order!

blas_lapack_production = ['mkl_lapack', 'mkl', 'svml']

other_libraries = ['boost_serialization', 'xerces-c', 'z', 'hdf5']

use_cvode = True
if use_cvode:
    other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])

tools = {'mpirun': '/home/pc-gary/wolf1768/bin/mpirun',
         'mpicxx': '/home/pc-gary/wolf1768/bin/mpicxx'}



