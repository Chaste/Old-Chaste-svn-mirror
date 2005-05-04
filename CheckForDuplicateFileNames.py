#!/usr/bin/env python

# Check through all the Chaste source directories, noting any source files
# that have the same name.

exts = ['.cpp', '.hpp']
dir_ignores = ['build']
chaste_dir = '.'

import os

# Dictionary mapping file names to locations
source_files = {}

for root, dirs, files in os.walk(chaste_dir):
    # Check for ignored dirs
    for dir in dir_ignores:
        if dir in dirs:
            dirs.remove(dir)
    # Check for source files
    for file in files:
        name, ext = os.path.splitext(file)
        if ext in exts:
            if source_files.has_key(file):
                # We've already found a file with this name
                source_files[file].append(os.path.join(root, file))
            else:
                # This is the first occurence of this name
                source_files[file] = [os.path.join(root, file)]

# Now check dictionary for duplicates
for file in source_files:
    if len(source_files[file]) > 1:
        print "Duplicate occurrences of",file,":"
        for loc in source_files[file]:
            print "  ",loc
