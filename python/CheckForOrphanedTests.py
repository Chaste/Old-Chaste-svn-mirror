#!/usr/bin/env python

# This script checks all the test folders to see if there are any
# Test*.hpp files that aren't listed in a test pack.
# It expects to be run from the trunk of the Chaste distribution.

chaste_dir = '.'

def IsTestFile(test_dir, test_file):
  "Does the given file define a test suite?"
  is_test = False
  if test_file[-4:] == '.hpp' and test_file[:4] == 'Test':
    fp = open(os.path.join(test_dir, test_file))
    for line in fp:
      if line.find('CxxTest::TestSuite') > -1:
        is_test = True
        break
    fp.close()
  return is_test


import os, glob

test_packs  = []  # Names of test packs found
orphans     = []  # Names of any orphaned test files
found_tests = []  # Names of tests found in test packs

# First get a list of all tests in test packs
test_pack_files = glob.glob('*/test/*TestPack.txt')

for test_pack_file in test_pack_files:
  # Add to list of test packs?
  test_pack = os.path.basename(test_pack_file)[:-12]
  if not test_pack in test_packs:
    test_packs.append(test_pack)
  # Add all tests in this file
  fp = file(test_pack_file)
  for line in fp:
    if not line in found_tests:
      found_tests.append(line.strip())
  fp.close()

# Now check for orphaned tests
test_dirs = glob.glob('*/test/')

for test_dir in test_dirs:
  for file in os.listdir(test_dir):
    if IsTestFile(test_dir, file) and not file in found_tests:
      orphans.append(os.path.join(test_dir, file))

# Output the names of test packs found
if test_packs:
  print "Test packs found:"
  for test_pack in test_packs:
    print " ", test_pack
  print

# Output any orphaned tests found
if orphans:
  print "Orphaned tests found:"
  for orphan in orphans:
    print " ", orphan
  print
  print "The next line is for the benefit of the test summary scripts."
  n_orphans, n_found = len(orphans), len(found_tests)
  print "Failed",n_orphans,"of",n_orphans+n_found,"tests"

  # Return a non-zero exit code if orphans were found
  import sys
  sys.exit(n_orphans)
