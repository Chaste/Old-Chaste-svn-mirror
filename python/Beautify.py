#!/usr/bin/env python

"""Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>."""

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
            
