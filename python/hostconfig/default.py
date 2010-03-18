# Configuration

"""Copyright (C) University of Oxford, 2005-2010

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
import sys

############################################################
# TO CONFIGURE YOUR MACHINE: if you have followed the manual
# installation instructions, edit the definition of
# chaste_libs_path below to point to where you installed the
# dependencies.
############################################################

#EDIT HERE
#For a simple installation all paths will be below this directory
chaste_libs_path = '/home/scratch/chaste-libs/'
#EDIT HERE

if not os.path.exists(chaste_libs_path) or not os.path.isdir(chaste_libs_path):
    print >>sys.stderr, "Chaste dependencies folder", chaste_libs_path, \
        "not found; please edit python/hostconfig/default.py"
    sys.exit(1)

petsc_2_2_path = ''
petsc_2_3_path = chaste_libs_path+'petsc-2.3.3-p15/'
petsc_build_name = 'linux-gnu'
petsc_build_name_profile = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
dealii_path = ''
metis_path = chaste_libs_path+'/metis-4.0/'
parmetis_path = chaste_libs_path+'/ParMetis-3.1/'
intel_path = '/opt/intel/cc/9.1.039/lib'
icpc = 'icpc'

other_includepaths = [chaste_libs_path+'hdf5/include',
                      chaste_libs_path+'xerces/include',
                      chaste_libs_path+'boost/include/boost-1_34_1',
                      chaste_libs_path+'xsd-3.2.0-i686-linux-gnu/libxsd']

other_libpaths = [chaste_libs_path+'lib',
                  chaste_libs_path+'boost/lib', 
                  chaste_libs_path+'xerces/lib',
                  chaste_libs_path+'hdf5/lib',
                  os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu')]


blas_lapack = ['f2clapack', 'f2cblas']
other_libraries = ['boost_serialization', 'xerces-c', 'hdf5', 'z']
# Note that boost serialization sometimes has a different name:
# other_libraries = ['boost_serialization-gcc41', 'xerces-c', 'hdf5', 'z', 'metis']


# use_vtk set to false initially. Change to True if VTK development libraries are 
# available.
use_vtk = False
# use_parmetis set to false initially. Change to True if parMETIS is available.
use_parmetis = False


tools = {'mpirun': chaste_libs_path+'mpi/bin/mpirun',
         'mpicxx': chaste_libs_path+'mpi/bin/mpicxx',
         'xsd': chaste_libs_path+'xsd-3.2.0-i686-linux-gnu/bin/xsd'}

#Extra libraries for VTK output
if use_vtk:
    other_includepaths.append(chaste_libs_path+'/Vtk5/include/vtk-5.2')
    other_libpaths.append(chaste_libs_path+'/Vtk5/lib/vtk-5.2')
    other_libraries.extend(['vtkFiltering','vtkIO',  'vtkCommon', 'vtksys', 'vtkzlib'])

# We need to make sure that we pick the right library since ParMETIS 3.1 requires a symbol not defined 
# in METIS 4.0.
if use_parmetis:
    other_includepaths.append(parmetis_path)
    other_libpaths.append(parmetis_path)
    other_libraries.extend(['parmetis', 'metis'])
else:
    other_includepaths.append(os.path.join(metis_path,'Lib'))
    other_libpaths.append(metis_path)
    other_libraries.extend(['metis'])

