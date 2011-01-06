#!/usr/bin/env python


"""Copyright (C) University of Oxford, 2005-2011

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


# This script checks all the test folders to see if there are any
# Test*.hpp files that aren't listed in a test pack.
# It expects to be run from the trunk of the Chaste distribution.

import glob
import os
import re

# Compatibility with Python 2.3
try:
    set = set
except NameError:
    import sets
    set = sets.Set

# Pre-2.6 compatibility (on posix)
try:
    relpath = os.path.relpath
except AttributeError:
    def relpath(path, start=os.path.curdir):
        """Return a relative version of a path"""
    
        if not path:
            raise ValueError("no path specified")
        
        start_list = os.path.abspath(start).split(os.path.sep)
        path_list = os.path.abspath(path).split(os.path.sep)
        
        # Work out how much of the filepath is shared by start and path.
        i = len(os.path.commonprefix([start_list, path_list]))
    
        rel_list = [os.path.pardir] * (len(start_list)-i) + path_list[i:]
        return os.path.join(*rel_list)


chaste_dir = '.'

suite_res = {'.hpp': re.compile(r'class\s+(\w+)\s*:\s*public\s+((::)?\s*CxxTest\s*::\s*)?\w*TestSuite\s+$'),
             '.py': re.compile(r'class\s+\w+\(unittest\.TestCase\):\s+$')}

def IsTestFile(test_dir, test_file_path):
    """Does the given file define a test suite?"""
    is_test = False
    test_file = os.path.basename(test_file_path)
    test_ext = os.path.splitext(test_file)[1]
    if test_file[:4] == 'Test' and test_ext in suite_res.keys():
        fp = open(os.path.join(test_dir, test_file_path))
        for line in fp:
            m = suite_res[test_ext].match(line)
            if m:
                is_test = True
                break
        fp.close()
    print test_dir, test_file, test_ext, is_test
    return is_test

test_packs  = set()  # Names of test packs found
orphans     = set()  # Names of any orphaned test files
found_tests = set()  # Names of tests found in test packs

# First get a list of all tests in all test packs
test_pack_files = glob.glob('*/test/*TestPack.txt')
for test_pack_file in test_pack_files:
    # Add to list of test packs?
    test_pack = os.path.basename(test_pack_file)[:-12]
    test_packs.add(test_pack)
    # Add all tests in this file
    fp = file(test_pack_file)
    for line in fp:
        line = line.strip()
        if line:
            found_tests.add(line)
    fp.close()


# Now check for orphaned tests in each top-level dir
test_dirs = glob.glob('*/test/')

local_found_tests = {} # Names of tests found in test packs in each folder

for test_dir in test_dirs:
    local_found_tests[test_dir] = set()
    test_pack_files = glob.glob(test_dir + '*TestPack.txt')
    for test_pack_file in test_pack_files:
        # Update list of tests that should be in this folder
        fp = file(test_pack_file)
        for line in fp:
            line = line.strip()
            if line:
                local_found_tests[test_dir].add(line)
        fp.close()
    # Check for orphans in this folder
    for dirpath, dirnames, filenames in os.walk(test_dir):
        for dirname in dirnames[:]:
            if dirname in ['.svn', 'data']:
                dirnames.remove(dirname)
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            filepath = relpath(filepath, test_dir)
            if IsTestFile(test_dir, filepath):
                if not filepath in local_found_tests[test_dir]:
                    orphans.add(os.path.join(test_dir, filepath))
                else:
                    local_found_tests[test_dir].remove(filepath)

# Output the names of test packs found
if test_packs:
    print "Test packs found:"
    for test_pack in test_packs:
        print " ", test_pack
    print

# Compute a list of tests listed in test packs without .hpp files
not_found = []
for test_dir in local_found_tests.keys():
    for test_file in local_found_tests[test_dir]:
        not_found.append(test_dir + test_file)

# Display results
if orphans or not_found:
    if orphans:
        print "Orphaned tests found:"
        for orphan in orphans:
            print " ", orphan
        print
    if not_found:
        print "Tests that don't exist:"
        for test in not_found:
            print " ", test
        print
    print "The next line is for the benefit of the test summary scripts."
    n_orphans, n_found = len(orphans), len(found_tests)
    print "Failed",n_orphans,"of",n_orphans+n_found,"tests"

    # Return a non-zero exit code if problems were found
    import sys
    sys.exit(n_orphans + len(not_found))
else:
    print "Infrastructure test passed ok."
  
