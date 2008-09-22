
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

"""Host-specific configuration.

This module contains the logic to pick up the right configuration file
for the machine it is run on, based on the machine's hostname.

Each configuration file should be a Python program providing certain variables
specifying where to find libraries and tools.  These variables are:
 * petsc_2_2_path - path to PETSc 2.2 install (None if not present)
 * petsc_2_3_path - path to PETSc 2.3 install (None if not present)
 * dealii_path    - path to Deal.II install (None if not present)
 * metis_path     - path to METIS install (None if not present)
 * intel_path     - path to Intel compiler installation

 * other_includepaths - list of paths containing other header files
 * other_libpaths     - list of paths containing other libraries, including metis, xsd, and boost
 * other_libraries    - list of other libraries to link against.  This *must* include
                        lapack, blas, and boost_serialization, as their names vary across systems.

 * ccflags - any extra compiler flags needed, as a string

 * tools - a dictionary mapping executable names to their absolute paths, for tools
           not found on $PATH.

Any non-absolute paths will be considered relative to the root of the Chaste install.
"""

import os
import socket
import sys

machine_fqdn = socket.getfqdn()

if machine_fqdn in ["userpc30.comlab.ox.ac.uk", "userpc33.comlab.ox.ac.uk"]:
    import joe as conf
elif machine_fqdn in [ "userpc60.comlab.ox.ac.uk", "userpc58.comlab.ox.ac.uk" ]:
    import chaste32 as conf
elif machine_fqdn in ["userpc44.comlab.ox.ac.uk", "userpc59.comlab.ox.ac.uk", "userpc36.comlab.ox.ac.uk" ]:
    import chaste64 as conf
elif machine_fqdn == "chaste-bob.comlab.ox.ac.uk":
    import chastebob as conf
elif machine_fqdn == "userpc36.comlab.ox.ac.uk":
    import migb as conf
elif machine_fqdn == "zuse.osc.ox.ac.uk":
    import zuse as conf
elif machine_fqdn.endswith(".comlab.ox.ac.uk"):
    import comlab as conf
elif machine_fqdn.endswith(".dtc.ox.ac.uk"):
    import dtc as conf
elif machine_fqdn.lower().startswith('finarfin'):
    import finarfin as conf
elif machine_fqdn.endswith(".maths.nottingham.ac.uk"):
    import nottingham as conf
elif machine_fqdn.endswith(".maths.ox.ac.uk"):
    import maths as conf
elif machine_fqdn.startswith('alex-laptop'):
    import alexf as conf
else:
    import default as conf
# For debugging
#for name in dir(conf):
#    if name[0] != '_':
#        print name, '=', getattr(conf, name)

# This is a bit ugly at present: SConstruct calls configure() to fill
# these global variables in, then reads them directly.
# Note that the order in which things are added to these lists often matters!
libpaths = []
incpaths = []
libraries = []

def do_petsc(version, optimised, profile=False, production=False, includes_only=False):
    """Determine PETSc include and library paths.

    The locations vary depending on the version of PETSc, and possibly
    whether optimised libraries are to be used.

    The version can be given as 2_2 or 2_3 to choose PETSc minor version.
    If a host doesn't support 2.3, we attempt to use 2.2 instead.  A
    ValueError is raised if 2.2 isn't present when asked for.

    Set optimised to True to use optimised builds of the libraries rather
    than debug builds.
    Set profile to True to use profile builds of PETSc.
    """
    if version == '2_3' and conf.petsc_2_3_path is None:
        # Use 2.2 instead
        version = '2_2'
    if version == '2_2' and conf.petsc_2_2_path is None:
        # Raise a friendly error
        raise ValueError('PETSc 2.2 required, but no path given in the host config.')
    if version == '2_2':
        petsc_base = os.path.abspath(conf.petsc_2_2_path)
        # Gracefully fall back to optimised/non-opt if the requested one isn't there
        if optimised:
            dirs = ['libO_c++', 'libg_c++']
        else:
            dirs = ['libg_c++', 'libO_c++']
        for d in dirs:
            libpath = os.path.join(petsc_base, 'lib', d, conf.petsc_build_name)
            if os.path.exists(libpath): break
        else:
            raise ValueError('No PETSc 2.2 libraries found.')
        incpaths.append(os.path.join(petsc_base, 'bmake', conf.petsc_build_name))
    else:
        petsc_base = os.path.abspath(conf.petsc_2_3_path)
        if production:
            build_name = conf.petsc_build_name_production
        elif profile:
            optimised = False
            build_name = conf.petsc_build_name_profile
        elif optimised:
            build_name = conf.petsc_build_name_optimized
        else:
            build_name = conf.petsc_build_name
        libpath = os.path.join(petsc_base, 'lib', build_name)
        incpaths.append(os.path.join(petsc_base, 'bmake', build_name))
    incpaths.append(os.path.join(petsc_base, 'include'))
    if not includes_only:
        libpaths.append(libpath)
        libraries.extend(['petscts', 'petscsnes', 'petscksp', 'petscdm', 
                          'petscmat', 'petscvec', 'petsc'])

def do_metis():
    """Add METIS include and library paths."""
    if conf.metis_path is None:
        raise ValueError('METIS required, but no path given in the host config.')
    libpath = os.path.abspath(conf.metis_path)
    incpath = os.path.join(libpath, 'Lib') # Yes, METIS is odd!
    libpaths.append(libpath)
    incpaths.append(incpath)
    libraries.append('metis')

def do_dealii(build):
    """Add Deal.II include & library paths, and libraries.

    Deal.II uses different library *names* to distinguish optimised versions.
    """
    if conf.dealii_path is None:
        raise ValueError('Deal.II required, but no path given in the host config.')
    base = os.path.abspath(conf.dealii_path)
    # Check Deal.II version
    version = open(os.path.join(base, 'Version')).read()
    if not version.startswith('6.'):
        # Older Deal.II requires PETSc 2.2
        do_petsc('2_2', build.is_optimised)
    else:
        # Just pick up the header files
        do_petsc('2_3', build.is_optimised, includes_only=True)
    # Add Deal.II libraries
    libpaths.append(os.path.join(base, 'lib'))
    relative_incpaths = ['base/include', 'lac/include', 'deal.II/include']
    incpaths.extend(map(lambda relpath: os.path.join(base, relpath),
                        relative_incpaths))
    libs = ['deal_II_1d', 'deal_II_2d', 'deal_II_3d', 'lac', 'base']
    if version.startswith('6.'):
        libs.append('petscall')
    if build.dealii_debugging:
        libs = map(lambda s: s + '.g', libs)
    libraries.extend(libs)

def configure(build):
    """Given a build object (BuildTypes.BuildType instance), configure the build."""
    if build.using_dealii:
        do_dealii(build)
        do_metis()
        libraries.extend(['blas', 'lapack']) # Use versions provided with Deal.II
    else:
        do_petsc('2_3', build.is_optimised, build.is_profile, build.is_production)
        if build.is_production:
            libraries.extend(conf.blas_lapack_production)
        else:
            libraries.extend(conf.blas_lapack)
    if build.CompilerType() == 'intel':
        intel_path = os.path.abspath(conf.intel_path)
        libpaths.append(os.path.join(intel_path, 'lib'))
    incpaths.extend(conf.other_includepaths)
    libpaths.extend(map(os.path.abspath, conf.other_libpaths))
    libraries.extend(conf.other_libraries)

    build.tools.update(conf.tools)

    if build.CompilerType() == 'intel':
        # Switch to use Intel toolchain
	build.tools['mpicxx'] += ' -CC="'+conf.icpc+'"'
        build.tools['cxx'] = os.path.join(intel_path, 'bin', 'icpc')
        build.tools['ar'] = os.path.join(intel_path, 'bin', 'xiar')

def ccflags():
    try:
        return conf.ccflags
    except:
        return ''
