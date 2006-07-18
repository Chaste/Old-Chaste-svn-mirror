#!/usr/bin/env python

# Beautify source files using astyle

exts = ['.cpp', '.hpp']
dir_ignores = ['build', 'cxxtest', 'testoutput', 'doc', 'anim']
chaste_dir = '.'
astyle_options_file = 'astylerc'

import os

os.chdir(chaste_dir)

# AStyle command line template
cmd = "astyle --options=%s %%s" % astyle_options_file

for root, dirs, files in os.walk(chaste_dir):
    # Check for ignored dirs
    for dir in dir_ignores:
        if dir in dirs:
            dirs.remove(dir)
    # Check for source files
    for file in files:
        name, ext = os.path.splitext(file)
        if ext in exts:
            # Run astyle
            os.system(cmd % os.path.join(root, file))
            
