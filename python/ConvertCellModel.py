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

import operator
import optparse
import os
import subprocess
import sys


# Use external PyCml if requested
if 'PYCML_DIR' in os.environ and os.path.isdir(os.environ['PYCML_DIR']):
    print 'Using external PyCml from PYCML_DIR =', os.environ['PYCML_DIR']
    pycml_dir = os.environ['PYCML_DIR']
else:
    pycml_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pycml')

# Options that we will supply to PyCml anyway
essential_options = ['--conf=config.xml',
                     '--use-chaste-stimulus',
                     '--convert-interfaces']
validation_options = ['-u', '--Wu']
# Options supplied if the user doesn't give a config file
default_options = []

# Parse command-line options
class OptionParser(optparse.OptionParser):
    def _process_args(self, largs, rargs, values):
        """
        Override to catch BadOption errors and pass these options on to PyCml.
        """
        while rargs:
            arg = rargs[0]
            if arg == "--":
                del rargs[0]
                return
            elif arg[0:2] == "--":
                try:
                    self._process_long_opt(rargs, values)
                except optparse.BadOptionError:
                    largs.append(arg)
            elif arg[0] == "-":
                largs.append(rargs.pop(0))
            elif self.allow_interspersed_args:
                largs.append(arg)
                del rargs[0]
            else:
                return                  # stop now, leave this arg in rargs
usage = 'usage: %prog [options] <file1.cellml> ...'
parser = OptionParser(usage=usage)
parser.add_option('--opt', action='store_true', default=False,
                  help="apply default optimisations to all generated code")
parser.add_option('--normal', action='store_true', default=False,
                  help="generate standard cell model code, suitable for"
                  " use with Chaste's ODE solvers")
parser.add_option('--cvode', action='store_true', default=False,
                  help="generate a subclass of AbstractCvodeCell,"
                  " i.e. that can be simulated with CVODE")
parser.add_option('--backward-euler', action='store_true', default=False,
                  help="generate a version of the cell model that can be"
                  " solved using a backward Euler method.  Requires the"
                  " presence of a .out file accompanying the CellML.")
parser.add_option('--output-dir', action='store',
                  help="directory to place output files in")
parser.add_option('--config-file',
                  action='store',
                  help="pathname of configuration file.  This can be used to"
                  " override options supplied on the command line, and will"
                  " also be passed on to PyCml itself, as will any options"
                  " not listed above.")
options, args = parser.parse_args()

# Read further arguments from config file?
if options.config_file:
    essential_options.append('--conf=' + options.config_file)
    # Parse the config file and extract any options
    sys.path[0:0] = [pycml_dir]
    import pycml
    rules = [pycml.bt.ws_strip_element_rule(u'*')]
    config_doc = pycml.amara_parse(options.config_file, rules=rules)
    if hasattr(config_doc.pycml_config, 'command_line_args'):
        config_args = map(str, config_doc.pycml_config.command_line_args.arg)
        options, extra_args = parser.parse_args(config_args, options)
        args.extend(extra_args)
    del config_doc

# If no options supplied, default to --normal --opt
option_names = ['opt', 'normal', 'cvode', 'backward_euler']
number_of_options = len(filter(None, operator.attrgetter(*option_names)(options)))
if number_of_options == 0:
    options.opt = True
    options.normal = True
    number_of_options = 2

# Whether to do a separate validation run
if number_of_options > 1:
    essential_options.append('--assume-valid')

# What options should be passed on to PyCml?
pycml_options = filter(lambda a: not a.endswith('.cellml'), args)
if not pycml_options:
    pycml_options = default_options
pycml_options.extend(essential_options)

# Models to process
models = []
for model in filter(lambda a: a.endswith('.cellml'), args):
    if os.path.exists(model):
        models.append(os.path.abspath(model))
    else:
        print >>sys.stderr, "Skipping", model, "because it does not exist"

def do_cmd(cmd):
    """Print and execute a command."""
    print ' '.join(cmd)
    rc = subprocess.call(cmd)
    if rc:
        sys.exit(rc)

def add_out_opts(options, output_dir, classname, filename):
    return options + ['-c', classname, '-o', os.path.join(output_dir, filename)]

def convert(model, output_dir):
    """The main workhorse function."""
    model_dir = os.path.dirname(model)
    model_base = os.path.basename(model)
    model_base = os.path.splitext(model_base)[0]
    class_name = "Cell" + model_base + "FromCellML"
    if not output_dir:
        output_dir = model_dir

    if number_of_options > 1:
        # Run validation separately
        cmd = ['./validator.py'] + validation_options + [model]
        do_cmd(cmd)

    command_base = ['./translate.py', model] + pycml_options

    if options.normal:
        # Basic class
        cmd = add_out_opts(command_base, output_dir, class_name, model_base + '.cpp')
        do_cmd(cmd)

    if options.opt and (options.normal or number_of_options == 1):
        # Normal with optimisation
        cmd = add_out_opts(command_base + ['-p', '-l'], output_dir,
                           class_name + 'Opt', model_base + 'Opt.cpp')
        do_cmd(cmd)
    
    if options.cvode:
        # For use with CVODE
        cmd = add_out_opts(command_base + ['-t', 'CVODE'], output_dir,
                           class_name + 'Cvode', model_base + 'Cvode.cpp')
        do_cmd(cmd)

        if options.opt:
            # With optimisation
            cmd = add_out_opts(command_base + ['-p', '-l', '-t', 'CVODE'],
                               output_dir,
                               class_name + 'CvodeOpt',
                               model_base + 'CvodeOpt.cpp')
            do_cmd(cmd)
    
    if options.backward_euler:
        maple_output = os.path.splitext(model)[0] + '.out'
        be_opts = ['-j', maple_output, '-p', '-l']
        cmd = add_out_opts(command_base + ['-j', maple_output, '-p', '-l'],
                           output_dir,
                           class_name + 'BackwardEuler',
                           model_base + 'BackwardEuler.cpp')
        do_cmd(cmd)


# TODO: This is bad for scons -j
os.chdir(pycml_dir)

for model in models:
    convert(model, options.output_dir)
