#!/usr/bin/env python

import os
import sys

options = ['--conf=config.xml', '--use-chaste-stimulus', '--convert-interfaces', '--Wu']
# Other possible options:
#   --assume-valid --no-member-vars
#   --row-lookup-method
#   --lt-index-uses-floor

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

def convert(model):
    model_dir = os.path.dirname(model)
    model_base = os.path.basename(model)
    model_base = os.path.splitext(model_base)[0]
    class_name = "Cell" + model_base + "FromCellML"

    command_base = './translate.py %(opts)s -c %(classname)s %(model)s -o %(outfile)s'

    # Basic class
    cmd = command_base % {'opts': ' '.join(options),
                          'classname': class_name,
                          'model': model,
                          'outfile': os.path.join(model_dir, model_base + '.hpp')}
    print cmd
    os.system(cmd)

    # With optimisation
    cmd = command_base % {'opts': ' '.join(['-p -l'] + options),
                          'classname': class_name + 'Opt',
                          'model': model,
                          'outfile': os.path.join(model_dir, model_base + 'Opt.hpp')}
    print cmd
    os.system(cmd)

os.chdir(os.environ['PYCML_DIR'])
for model in models:
    convert(model)
