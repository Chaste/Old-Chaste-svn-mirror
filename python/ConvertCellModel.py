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
import subprocess
import sys

options = ['--conf=config.xml',
           '--use-chaste-stimulus',
           '--convert-interfaces',
           '--Wu',
           '--row-lookup-method']
# Other possible options:
#   --assume-valid --no-member-vars
#   --lt-index-uses-floor

validation_options = ['-u', '--Wu']

# Use external PyCml if requested
if 'PYCML_DIR' in os.environ and os.path.isdir(os.environ['PYCML_DIR']):
    print 'Using external PyCml from PYCML_DIR =', os.environ['PYCML_DIR']
    pycml_dir = os.environ['PYCML_DIR']
else:
    pycml_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pycml')

# Poor man's argument parsing
our_args = ['--cvode', '--normal', '--opt', '--backward-euler']
def arg2varname(arg):
    return arg[2:].replace('-', '_')
number_of_args = 0
for arg in our_args:
    if arg in sys.argv:
        sys.argv.remove(arg)
        exec("use_%s = True" % arg2varname(arg))
        number_of_args += 1
    else:
        exec("use_%s = False" % arg2varname(arg))
if number_of_args == 0:
    use_cvode = False
    use_normal = True
    use_opt = True
    use_backward_euler = False
    number_of_args = 2
if number_of_args > 1:
    options.append('--assume-valid')

if '--output-dir' in sys.argv:
    i = sys.argv.index('--output-dir')
    output_dir = sys.argv[i+1]
    sys.argv[i:i+2] = []
else:
    output_dir = None

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
    
def do_cmd(cmd):
    """Print and execute a command."""
    print cmd
    rc = subprocess.call(cmd, shell=True)
    if rc:
        sys.exit(rc)

def convert(model, output_dir):
    """The main workhorse function."""
    model_dir = os.path.dirname(model)
    model_base = os.path.basename(model)
    model_base = os.path.splitext(model_base)[0]
    class_name = "Cell" + model_base + "FromCellML"
    if output_dir is None:
        output_dir = model_dir

    if number_of_args > 1:
        # Run validation separately
        cmd = './validator.py ' + ' '.join(validation_options) + ' ' + model
        do_cmd(cmd)

    command_base = './translate.py %(opts)s -c %(classname)s %(model)s -o %(outfile)s'

    if use_normal:
        # Basic class
        cmd = command_base % {'opts': ' '.join(options),
                              'classname': class_name,
                              'model': model,
                              'outfile': os.path.join(output_dir, model_base + '.cpp')}
        do_cmd(cmd)

    if use_opt and (use_normal or not use_cvode):
        # With optimisation
        cmd = command_base % {'opts': ' '.join(['-p -l'] + options),
                              'classname': class_name + 'Opt',
                              'model': model,
                              'outfile': os.path.join(output_dir, model_base + 'Opt.cpp')}
        do_cmd(cmd)
    
    if use_cvode:
        cmd = command_base % {'opts': ' '.join(['-t CVODE'] + options),
                              'classname': class_name + 'Cvode',
                              'model': model,
                              'outfile': os.path.join(output_dir, model_base + 'Cvode.cpp')}
        do_cmd(cmd)

        if use_opt:
            # With optimisation
            cmd = command_base % {'opts': ' '.join(['-p -l', '-t CVODE'] + options),
                                  'classname': class_name + 'CvodeOpt',
                                  'model': model,
                                  'outfile': os.path.join(output_dir, model_base + 'CvodeOpt.cpp')}
            do_cmd(cmd)
    
    if use_backward_euler:
        maple_output = os.path.splitext(model)[0] + '.out'
        be_opts = ['-j', maple_output, '-p', '-l']
        cmd = command_base % {'opts': ' '.join(be_opts + options),
                              'classname': class_name + 'BackwardEuler',
                              'model': model,
                              'outfile': os.path.join(output_dir, model_base + 'BackwardEuler.cpp')}
        do_cmd(cmd)


os.chdir(pycml_dir)

for model in models:
    convert(model, output_dir)
