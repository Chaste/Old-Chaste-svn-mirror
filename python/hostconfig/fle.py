# Configuration for cluster3 at FLE

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

petsc_2_3_path = '/opt/petsc/'
petsc_build_name = 'linux-x86_64'
petsc_build_name_optimized = 'linux-x86_64'
dealii_path = None
metis_path = None
intel_path = '/opt/intel/fce/10.1.013/'
icpc = 'icpc'

other_includepaths = ['/home/southern/boost/include/boost-1_34',
                      '/opt/hypre/include',
                      '/share/apps/opt/opt/include',
                      '/home/southern/xsd-2.3.1-i686-linux-gnu/libxsd',
                      '/home/southern/xerces-c-src_2_7_0/include',
		      '/home/southern/hdf5/include',
		      '/home/southern/deal.II/deal.II/include',
		      '/home/southern/deal.II/base/include', 
		      '/opt/infinipath-mpi/include' ]
other_libpaths = ['/opt/intel/mkl/9.0/lib/em64t',
                  '/usr/X11R6/lib64',
                  '/opt/hypre/lib',
		  '/home/southern/lib',
                  '/home/southern/boost/lib',
		  '/home/southern/hdf5/lib',
		  '/opt/infinipath-mpi/lib64' ]
blas_lapack = ['mkl_lapack64', 'mkl_em64t']
other_libraries = ['mkl_em64t', 'mkl', 'X11', 'HYPRE', 'boost_serialization-gcc34-1_34',
                     'xerces-c', 'svml', 'imf', 'irc', 'z', 'hdf5' ]

tools = {}
#tools = {'mpirun': '/opt/infinipath-mpi/bin/mpirun',
#	 'mpicxx': '/opt/infinipath-mpi/bin/mpicxx'}
