#!/usr/bin/env python

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

"""
A little helper script to facilitate calling PyCml for common Chaste usage scenarios.
"""

import os
import sys

options = ['--conf=config.xml',
           '--use-chaste-stimulus',
           '--convert-interfaces',
           '--Wu',
           '--row-lookup-method']
# Other possible options:
#   --assume-valid --no-member-vars
#   --lt-index-uses-floor

# Use external PyCml if requested
if 'PYCML_DIR' in os.environ and os.path.isdir(os.environ['PYCML_DIR']):
    print 'Using external PyCml from PYCML_DIR =', os.environ['PYCML_DIR']
    pycml_dir = os.environ['PYCML_DIR']
else:
    pycml_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pycml')

# Poor man's argument parsing
our_args = ['--cvode', '--normal', '--opt']
found_any_arg = False
for arg in our_args:
    if arg in sys.argv:
        sys.argv.remove(arg)
        exec("use_%s = True" % arg[2:])
        found_any_arg = True
    else:
        exec("use_%s = False" % arg[2:])
if not found_any_arg:
    use_cvode = False
    use_normal = True
    use_opt = True

# What options should be passed on to PyCml?
try:
    end = sys.argv.index('--')
    args = sys.argv[1:end]
    options.extend(sys.argv[end+1:])
except ValueError:
    args = sys.argv[1:]

options.extend([arg for arg in args if arg[0] == '-'])

models = []
for model in [arg for arg in args if arg[0] != '-']:
    models.append(os.path.abspath(model))

# The main workhorse function
def convert(model):
    model_dir = os.path.dirname(model)
    model_base = os.path.basename(model)
    model_base = os.path.splitext(model_base)[0]
    class_name = "Cell" + model_base + "FromCellML"

    command_base = './translate.py %(opts)s -c %(classname)s %(model)s -o %(outfile)s'

    if use_normal:
        # Basic class
        cmd = command_base % {'opts': ' '.join(options),
                              'classname': class_name,
                              'model': model,
                              'outfile': os.path.join(model_dir, model_base + '.hpp')}
        print cmd
        os.system(cmd)

        if use_opt:
            # With optimisation
            cmd = command_base % {'opts': ' '.join(['-p -l'] + options),
                                  'classname': class_name + 'Opt',
                                  'model': model,
                                  'outfile': os.path.join(model_dir, model_base + 'Opt.hpp')}
            print cmd
            os.system(cmd)
    
    if use_cvode:
        cmd = command_base % {'opts': ' '.join(['-t CVODE'] + options),
                              'classname': class_name + 'Cvode',
                              'model': model,
                              'outfile': os.path.join(model_dir, model_base + 'Cvode.hpp')}
        print cmd
        os.system(cmd)

        if use_opt:
            # With optimisation
            cmd = command_base % {'opts': ' '.join(['-p -l', '-t CVODE'] + options),
                                  'classname': class_name + 'CvodeOpt',
                                  'model': model,
                                  'outfile': os.path.join(model_dir, model_base + 'CvodeOpt.hpp')}
            print cmd
            os.system(cmd)


os.chdir(pycml_dir)

for model in models:
    convert(model)
