#!/usr/bin/env python

# This script checks all the test folders to see if there are any
# Test*.hpp files that aren't listed in a test pack.
# It expects to be run from the trunk of the Chaste distribution.

import os, glob, re

chaste_dir = '.'

suite_re = re.compile( r'\bclass\s+(\w+)\s*:\s*public\s+((::)?\s*CxxTest\s*::\s*)?\w*TestSuite\b' )
def IsTestFile(test_dir, test_file):
  """Does the given file define a test suite?"""
  is_test = False
  if test_file[-4:] == '.hpp' and test_file[:4] == 'Test':
    fp = open(os.path.join(test_dir, test_file))
    for line in fp:
      m = suite_re.search( line )
      if m:
        is_test = True
        break
    fp.close()
  return is_test

test_packs  = []  # Names of test packs found
orphans     = []  # Names of any orphaned test files
found_tests = []  # Names of tests found in test packs

# First get a list of all tests in all test packs
test_pack_files = glob.glob('*/test/*TestPack.txt')
for test_pack_file in test_pack_files:
  # Add to list of test packs?
  test_pack = os.path.basename(test_pack_file)[:-12]
  if not test_pack in test_packs:
    test_packs.append(test_pack)
  # Add all tests in this file
  fp = file(test_pack_file)
  for line in fp:
    line = line.strip()
    if line and not line in found_tests:
      found_tests.append(line)
  fp.close()


# Now check for orphaned tests in each top-level dir
test_dirs = glob.glob('*/test/')

local_found_tests = {} # Names of tests found in test packs in each folder

for test_dir in test_dirs:
  local_found_tests[test_dir] = []
  test_pack_files = glob.glob(test_dir + '*TestPack.txt')
  for test_pack_file in test_pack_files:
    # Update list of tests that should be in this folder
    fp = file(test_pack_file)
    for line in fp:
      line = line.strip()
      if line and not line in local_found_tests[test_dir]:
        local_found_tests[test_dir].append(line)
    fp.close()
  # Check for orphans in this folder
  for filename in os.listdir(test_dir):
    if IsTestFile(test_dir, filename):
      if not filename in local_found_tests[test_dir]:
        orphans.append(os.path.join(test_dir, filename))
      else:
        local_found_tests[test_dir].remove(filename)

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
  
