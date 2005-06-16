#!/usr/bin/env python

# Strip out non-interesting output from the build process, and only
# display errors, with minimal context.

# The first argument is expected to be the log filename, which
# can be a URL (e.g. http://comlab2.lsi.ox.ac.uk/out, which is
# the default). Use '-' to read from stdin.

import sys, re, urllib

if len(sys.argv) > 1:
  logfile = sys.argv[1]
else:
  logfile = "http://comlab2.lsi.ox.ac.uk/out"

print "Analysing build logs from",logfile,"..."

# Context variables
line_compile, lineno_compile = "", -1
line_testsuite, lineno_testsuite = "", -1
line_test, lineno_test = "", -1
line_doxygen, lineno_doxygen = "", -1

# Regular expressions
res_compile = []
res_compile.append(re.compile("^/usr/bin/g\+\+"))
res_compile.append(re.compile("^cxxtest/cxxtestgen.pl"))
res_compile.append(re.compile("^ar"))
res_compile.append(re.compile("^ranlib"))
res_compile.append(re.compile("^[a-zA-Z/\.]*mpicxx "))
res_testsuite = []
res_testsuite.append(re.compile("^ \*\*\*\*\* [\w\.]+ \*\*\*\*\*"))
res_test = []
res_test.append(re.compile("^Entering .*"))
res_doxygen = []
res_doxygen.append(re.compile("^[A-Z](\w+ )+\w+\.\.\.$"))
res_doxygen.append(re.compile("^Preprocessing "))
res_doxygen.append(re.compile("^Parsing file "))
res_doxygen.append(re.compile("^Generating "))
#res_doxygen.append(re.compile("^Notice: Output directory "))
res_doxygen.append(re.compile("^Reading and parsing tag files"))
res_doxygen.append(re.compile("^Searching for files to exclude"))
res_doxygen.append(re.compile("^Read \d+ bytes"))
res_doxygen.append(re.compile("^Freeing entry tree"))
res_doxygen.append(re.compile("^Determining which enums are documented"))
res_doxygen.append(re.compile("^Adding members to member groups\."))
res_other = []
res_other.append(re.compile("^Passed"))
res_other.append(re.compile("^Running \d+ test"))
res_other.append(re.compile("^OK!"))
res_other.append(re.compile("^Failed "))
res_other.append(re.compile("^Success rate:"))
res_other.append(re.compile("^$"))
res_other.append(re.compile("^[a-zA-Z/\.]*testrunner "))
res_other.append(re.compile("^[a-zA-Z/\.]*mpirun "))
res_other.append(re.compile("^[a-zA-Z/\.]*paralleltestrunner "))

if logfile == '-':
  fp = sys.stdin
else:
  fp = urllib.urlopen(logfile)

line, lineno = fp.readline(), 0
done_context = False
while line:
  ignore = False

  # Check for ignored lines
  for regexp in res_compile:
  # Check for a compilation command
    if regexp.match(line):
      ignore, done_context = True, False
      line_compile = line
      lineno_compile = lineno
      break
  if not ignore:
    for regexp in res_testsuite:
      # Check for the start of a test suite
      if regexp.match(line):
        ignore, done_context = True, False
        line_testsuite = line
        lineno_testsuite = lineno
        break
  if not ignore:
    for regexp in res_test:
      # Check for the start of a test
      if regexp.match(line):
        ignore, done_context = True, False
        line_test = line
        lineno_test = lineno
        break
  if not ignore:
    for regexp in res_doxygen:
      # Doxygen output
      if regexp.match(line):
        ignore, done_context = True, False
        line_doxygen = line
        lineno_doxygen = lineno
        break
  if not ignore:
    for regexp in res_other:
      # Any other ignored lines
      if regexp.match(line):
        ignore, done_context = True, False
        break

  if not ignore:
    if not done_context:
      l = [lineno_compile, lineno_testsuite, lineno_test, lineno_doxygen]
      i = l.index(max(l))
      if i == 0 and line_compile:
        print line_compile,
      elif i == 3 and line_doxygen:
        pass
      else:
        if line_testsuite: print line_testsuite,
        if line_test: print line_test,
      done_context = True
    print line,
  line = fp.readline()
  lineno = lineno + 1

fp.close()
