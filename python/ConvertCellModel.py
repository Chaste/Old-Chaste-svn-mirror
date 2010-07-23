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

import optparse
import os
import subprocess
import sys
import tempfile


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
    def __init__(self, *args, **kwargs):
        if 'our_short_options' in kwargs:
            self.__our_short_options = kwargs['our_short_options']
            del kwargs['our_short_options']
        else:
            self.__our_short_options = []
        self.__our_short_options.append('-h')
        optparse.OptionParser.__init__(self, *args, **kwargs)
        
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
                if arg in self.__our_short_options:
                    self._process_short_opts(rargs, values)
                else:
                    largs.append(rargs.pop(0))
            elif self.allow_interspersed_args:
                largs.append(arg)
                del rargs[0]
            else:
                return                  # stop now, leave this arg in rargs
usage = 'usage: %prog [options] <file1.cellml> ...'
parser = OptionParser(usage=usage, our_short_options=['-y'])
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
parser.add_option('--show-outputs', action='store_true', default=False,
                  help="don't actually run PyCml, just show what files would"
                  " be generated, one per line")
parser.add_option('--config-file',
                  action='store',
                  help="pathname of configuration file.  This can be used to"
                  " override options supplied on the command line, and will"
                  " also be passed on to PyCml itself, as will any options"
                  " not listed above.")
parser.add_option('-y', '--dll', '--dynamically-loadable',
                  dest='dynamically_loadable',
                  action='store_true', default=False,
                  help="add code to allow the model to be compiled to a "
                  "shared library and dynamically loaded.  If this option is "
                  "given, only one type of output will be generated (so you "
                  "can't combine, e.g. --cvode --normal).")
options, args = parser.parse_args()

option_names = ['opt', 'normal', 'cvode', 'backward_euler']
def arg2name(arg):
    return str(arg)[2:].replace('-', '_')

# Read further arguments from config file?
if options.config_file:
    # Parse the config file and extract any options
    import xml.dom.minidom
    def getText(nodelist):
        rc = []
        for node in nodelist:
            if node.nodeType == node.TEXT_NODE:
                rc.append(node.data)
        return ''.join(rc)
    config_doc = xml.dom.minidom.parse(options.config_file)
    config_modified = False
    cl_args_elt = config_doc.getElementsByTagName('command_line_args')
    if cl_args_elt:
        cl_args_elt = cl_args_elt[0]
        arg_elts = cl_args_elt.getElementsByTagName('arg')
        config_args = []
        # Strip from the file any arguments only understood by this script
        for arg_elt in arg_elts:
            arg = getText(arg_elt.childNodes)
            config_args.append(arg)
            if arg2name(arg) in option_names:
                cl_args_elt.removeChild(arg_elt)
                config_modified = True
        if config_modified:
            # If the config file supplied such arguments, then pretend there weren't
            # any on the command line, since we can't turn options off.
            for option in option_names:
                setattr(options, option, False)
            # An empty command_line_args isn't allowed
            for node in cl_args_elt.childNodes:
                if node.nodeType == node.ELEMENT_NODE:
                    break
            else:
                for node in config_doc.childNodes:
                    if node.nodeType == node.ELEMENT_NODE:
                        node.removeChild(cl_args_elt)
        # Parse additional options
        options, extra_args = parser.parse_args(config_args, options)
        args.extend(extra_args)
    if config_modified:
        # Write a new config file
        fp, config_path = tempfile.mkstemp(suffix='.xml', text=True)
        config_file = os.fdopen(fp, 'w')
        config_doc.writexml(config_file)
        config_file.close()
        essential_options.append('--conf=' + config_path)
    else:
        # Just pass on the one we were given
        essential_options.append('--conf=' + options.config_file)
    del config_doc


# If no options supplied, default to --normal --opt
number_of_options = len(filter(None, [getattr(options, opt_name) for opt_name in option_names]))
if number_of_options == 0:
    options.normal = True
    if not options.dynamically_loadable:
        options.opt = True
        number_of_options = 2
    else:
        number_of_options = 1

# Check for .so creation
if options.dynamically_loadable:
    if number_of_options > 1:
        parser.error("Only one output type may be specified if creating a dynamic library")
    essential_options.append('-y')

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

def tidy_up():
    """Clean up temporary file, if created."""
    if options.config_file and config_modified:
        os.remove(config_path)

def do_cmd(cmd, outputs):
    """Print and execute a command.
    
    If the command fails, remove any generated outputs and exit.
    """
    if options.show_outputs:
        for output in outputs:
            print output
    else:
        print ' '.join(cmd)
        rc = subprocess.call(cmd)
        if rc:
            for output in outputs:
                try:
                    os.remove(output)
                except OSError:
                    pass
            tidy_up()
            sys.exit(rc)

def add_out_opts(base_options, output_dir, classname, file_base, file_extra=''):
    """Add options specifying output path and class name.
    
    Returns extended options list and list of output file paths.
    """
    if options.dynamically_loadable:
        filename = file_base + '.cpp'
    else:
        filename = file_base + file_extra + '.cpp'
    cpp_path = os.path.join(output_dir, filename)
    full_options = base_options + ['-c', classname, '-o', cpp_path]
    outputs = [cpp_path, cpp_path[:-3] + 'hpp']
    return (full_options, outputs)

def convert(model, output_dir):
    """The main workhorse function."""
    model_dir = os.path.dirname(model)
    model_base = os.path.basename(model)
    model_base = os.path.splitext(model_base)[0]
    class_name = "Cell" + model_base.replace('-', '_') + "FromCellML"
    if not output_dir:
        output_dir = model_dir

    if number_of_options > 1:
        # Run validation separately
        cmd = ['./validator.py'] + validation_options + [model]
        do_cmd(cmd, [])

    command_base = ['./translate.py', model] + pycml_options

    if options.normal:
        # Basic class
        cmd, outputs = add_out_opts(command_base, output_dir, class_name, model_base)
        do_cmd(cmd, outputs)

    if options.opt and (options.normal or number_of_options == 1):
        # Normal with optimisation
        cmd, outputs = add_out_opts(command_base + ['-p', '-l'], output_dir,
                                    class_name + 'Opt', model_base, 'Opt')
        do_cmd(cmd, outputs)
    
    if options.cvode:
        # For use with CVODE
        cmd, outputs = add_out_opts(command_base + ['-t', 'CVODE'], output_dir,
                                    class_name + 'Cvode', model_base, 'Cvode')
        do_cmd(cmd, outputs)

        if options.opt:
            # With optimisation
            cmd, outputs = add_out_opts(command_base + ['-p', '-l', '-t', 'CVODE'],
                                        output_dir,
                                        class_name + 'CvodeOpt',
                                        model_base, 'CvodeOpt')
            do_cmd(cmd, outputs)
    
    if options.backward_euler:
        maple_output = os.path.splitext(model)[0] + '.out'
        be_opts = ['-j', maple_output, '-p', '-l']
        cmd, outputs = add_out_opts(command_base + ['-j', maple_output, '-p', '-l'],
                                    output_dir,
                                    class_name + 'BackwardEuler',
                                    model_base, 'BackwardEuler')
        do_cmd(cmd, outputs)


# TODO #1493: This is bad for scons -j
os.chdir(pycml_dir)

for model in models:
    convert(model, options.output_dir)

tidy_up()
