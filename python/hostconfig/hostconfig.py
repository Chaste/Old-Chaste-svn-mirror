
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

"""Host-specific configuration.

This module contains the logic to set various build parameters based on a
host configuration file.  Configuration files for developer machines may be
stored in a 'machines' package.  If a machine doesn't have a specific file,
default.py is used.

Each configuration file should be a Python program providing certain variables
specifying where to find libraries and tools.  In each case, if the library or
tool is not present, the variable should be set to None.  These variables are:
 * petsc_2_2_path - path to PETSc 2.2 install
 * petsc_2_3_path - path to PETSc 2.3 install
 * petsc_3_0_path - path to PETSc 3.0 install
 * dealii_path    - path to Deal.II install
 * metis_path     - path to METIS install
 * intel_path     - path to Intel compiler installation

These flags indicate whether to use certain optional external libraries.
Eventually we will use SCons' Configure functionality to make this unnecessary.
If the variable is not present, it will default to False.  If the variable is
set to True, then the appropriate paths (see below) should be specified, too.
 * use_cvode  - whether to use CVODE
 * use_vtk - whether to use VTK development libraries
 * use_adaptivity - whether to use adaptivity library
 
 * other_includepaths - list of paths containing other header files
 * other_libpaths     - list of paths containing other libraries, including
                        metis, xsd, and boost
 * other_libraries    - list of other libraries to link against
 * blas_lapack        - the names of the blas and lapack libraries on this system
                        (a 2 element list)

 * ccflags - any extra compiler flags needed, as a string.

 * tools - a dictionary mapping executable names to their absolute paths, for tools
           not found on $PATH.
 * icpc - a special case tool: allows you to change the Intel compiler command, e.g.
          to add flags specific to 64-bit machines.

Any non-absolute paths will be considered relative to the root of the Chaste install.
"""

import os
import sys

# Do we have any machine-specific config?
try:
    import machines
    conf = machines.config_module()
except ImportError:
    # How about distro-specific config?
    try:
        fp = open('/etc/issue')
        distro = fp.read().split()[0].lower()
        fp.close()
        conf = __import__(distro)
    except (ImportError, IOError):
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

    The version can be given as 2.2, 2.3 or 3.0 to choose PETSc version.
    If a host doesn't support 3.0 we attempt to use 2.3 or 2.2 respectively.
    A ValueError is raised if 2.2 isn't present when asked for.

    Set optimised to True to use optimised builds of the libraries rather
    than debug builds.
    Set profile to True to use profile builds of PETSc.
    """
    if  os.environ.get('XTPE_COMPILE_TARGET', ''):
        return

    conf.petsc_2_2_path = getattr(conf, 'petsc_2_2_path', None)
    conf.petsc_2_3_path = getattr(conf, 'petsc_2_3_path', None)
    conf.petsc_3_0_path = getattr(conf, 'petsc_3_0_path', None)
    requested_version = version
    if version == '3.0' and (conf.petsc_3_0_path is None or 
                             not os.path.isdir(conf.petsc_3_0_path)):
        # Use 2.3 instead
        version = '2.3'
    if version == '2.3' and (conf.petsc_2_3_path is None or 
                             not os.path.isdir(conf.petsc_2_3_path)):
        # Use 2.2 instead
        version = '2.2'
    if version == '2.2' and (conf.petsc_2_2_path is None or 
                             not os.path.isdir(conf.petsc_2_2_path)):
        # Raise a friendly error
        raise ValueError('PETSc %s requested, but no path for this or an earlier version given in the host config.' % requested_version)
    if version == '2.2':
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
    elif version == '2.3':
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
    else: #version == '3.0'
        petsc_base = os.path.abspath(conf.petsc_3_0_path)
        if production:
            build_name = conf.petsc_build_name_production
        elif profile:
            optimised = False
            build_name = conf.petsc_build_name_profile
        elif optimised:
            build_name = conf.petsc_build_name_optimized
        else:
            build_name = conf.petsc_build_name
        libpath = os.path.join(petsc_base, build_name, 'lib')
        incpaths.append(os.path.join(petsc_base, build_name, 'include'))
    incpaths.append(os.path.join(petsc_base, 'include'))
    if not includes_only:
        libpaths.append(libpath)
        libraries.extend(['petscts', 'petscsnes', 'petscksp', 'petscdm', 
                          'petscmat', 'petscvec', 'petsc'])
        if sys.platform == 'cygwin':
            libraries.extend(['gdi32', 'user32', 'advapi32', 'kernel32', 'dl'])

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
        do_petsc('2.2', build.is_optimised)
    else:
        # Just pick up the header files
        do_petsc('2.3', build.is_optimised, includes_only=True)
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

def optional_library_defines():
    """
    Work out what optional libraries have been asked for,
    and return the appropriate #define flags, as a list.
    """
    possible_flags = {'cvode': 'CHASTE_CVODE', 'vtk': 'CHASTE_VTK', 'adaptivity': 'CHASTE_ADAPTIVITY'}
    actual_flags = []
    for libname, symbol in possible_flags.iteritems():
        if getattr(conf, 'use_' + libname, False):
            actual_flags.append('-D' + symbol)
    return actual_flags

def configure(build):
    """Given a build object (BuildTypes.BuildType instance), configure the build."""
    if build.using_dealii:
        do_dealii(build)
        #do_metis()
        libraries.extend(conf.other_libraries) # Some of "other_libraries" may depend on BLAS/LAPACK, make sure they are included before them.
        libraries.extend(['blas', 'lapack']) # Use versions provided with Deal.II
    else:
        prefs = build.GetPreferedVersions()
        if prefs:
            if hasattr(conf, 'SetPreferedVersions') and callable(conf.SetPreferedVersions):
                conf.SetPreferedVersions(prefs)
            petsc_version = prefs['petsc'][:3]
        else:
            petsc_version = '3.0'
        do_petsc(petsc_version, build.is_optimised, build.is_profile, build.is_production) # PETSc links against some objects defined in "other_libraries"
        libraries.extend(conf.other_libraries) # Some of "other_libraries" may depend on BLAS/LAPACK, make sure they are included before them.
        if build.is_production:
            libraries.extend(conf.blas_lapack_production)
        else:
            libraries.extend(conf.blas_lapack)
    if build.CompilerType() == 'intel':
        intel_path = os.path.abspath(conf.intel_path)
        libpaths.append(os.path.join(intel_path, 'lib'))
    incpaths.extend(conf.other_includepaths)
    libpaths.extend(map(os.path.abspath, conf.other_libpaths))    

    build.tools.update(conf.tools)

    if build.CompilerType() == 'intel':
        # Switch to use Intel toolchain
        if hasattr(conf, 'icpc'):
            build.tools['mpicxx'] += ' -CC="'+conf.icpc+'"'
        build.tools['cxx'] = os.path.join(intel_path, 'bin', 'icpc')
        build.tools['ar'] = os.path.join(intel_path, 'bin', 'xiar')

    if hasattr(conf, 'ModifyBuild') and callable(conf.ModifyBuild):
        conf.ModifyBuild(build)

def ccflags():
    opt_lib_flags = optional_library_defines()
    conf_flags = getattr(conf, 'ccflags', '')
    return conf_flags + ' ' + ' '.join(opt_lib_flags)

def ldflags():
    return getattr(conf, 'ldflags', '')

