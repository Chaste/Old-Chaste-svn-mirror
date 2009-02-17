# Configuration for cluster3 at FLE

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

petsc_2_3_path = '/share/apps/petsc-2.3.3-p15-intel-10.1.013-openmpi-intel/'
petsc_build_name = 'linux'
petsc_build_name_optimized = 'linux'
dealii_path = None
metis_path = None
intel_path = '/opt/intel/fce/10.1.013/'
icpc = 'icpc'

other_includepaths = ['/home/southern/boost/include/boost-1_33_1',
                      '/share/apps/opt/opt/include',
                      '/home/southern/xsd-2.3.1-i686-linux-gnu/libxsd',
                      '/home/southern/xerces-c-src_2_7_0/include',
		      '/home/southern/hdf5-openmpi/include',
		      '/home/southern/deal.II/deal.II/include',
		      '/home/southern/deal.II/base/include',
		      '/usr/mpi/intel/openmpi-2.2/include' ]

other_libpaths = ['/opt/intel/mkl/9.0/lib/em64t',
                  '/usr/X11R6/lib64',
		  '/home/southern/lib',
                  '/home/southern/boost/lib',
		  '/home/southern/hdf5-openmpi/lib',
		  '/usr/mpi/intel/openmpi-2.2/lib64',
		  '/share/apps/ParMetis-3.1-intel-10.1.013-openmpi-intel/lib', 
		  os.path.join(petsc_2_3_path, 'externalpackages/hypre-2.0.0/linux/lib'),
		  os.path.join(petsc_2_3_path, 'externalpackages/Prometheus-1.8.6/linux/lib') ]
		  
blas_lapack = ['mkl_lapack64', 'mkl_em64t']
other_libraries = ['mkl_em64t', 'mkl', 'X11', 'HYPRE', 'boost_serialization-gcc-1_33_1',
                     'xerces-c', 'svml', 'imf', 'irc', 'z', 'hdf5', 'promfei', 'prometheus', 
		     'parmetis', 'metis' ]

tools = {}
#tools = {'mpirun': '/opt/infinipath-mpi/bin/mpirun',
#	 'mpicxx': '/opt/infinipath-mpi/bin/mpicxx'}
