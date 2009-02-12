# Configuration for Joe's machines

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

import os

petsc_2_2_path = None
petsc_2_3_path = '/home/jmpf/petsc-2.3.1-p16/'
#petsc_2_3_path = '/home/jmpf/petsc-2.3.2-p4/'
petsc_build_name = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
petsc_build_name_production = 'linux-intel-opt-mkl'
dealii_path = None
metis_path = None
intel_path = '' #opt/intel/Compiler/11.0/074/bin/ia32/iccvars_ia32.sh 
icpc = 'icpc'


other_includepaths = ['/home/jmpf/xsd-2.3.1-i686-linux-gnu/libxsd', '/home/jmpf/hdf5/include']
other_libpaths = [os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu/'),  
                    '/home/jmpf/hdf5/lib', '/home/jmpf/xerces/lib'] #/opt/intel/Compiler/11.0/074/mkl/tools/environment/mklvars32.sh
blas_lapack = ['f2clapack', 'f2cblas']
blas_lapack_production = ['mkl_lapack', 'mkl', 'svml']
other_libraries = ['boost_serialization', 'xerces-c', 'hdf5', 'z']

tools = {'mpirun': '/home/jmpf/mpi/bin/mpirun',
         'mpicxx': '/home/jmpf/mpi/bin/mpicxx'}


use_vtk = False
if use_vtk:
    other_libraries.extend([ 'vtkIO'])
    other_includepaths.extend(['/usr/include/vtk-5.0'])
    