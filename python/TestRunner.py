#!/usr/bin/env python

# Script to run a test executable and record the status
# (whether all tests passed or how many failed) as part of
# the output filename.

# The test executable should be run and the output captured,
# and stored to a file.

# This script expects 3 or 4 arguments.
# The first is the executable to run.
# The second is the name of the output .log file, in which the output
#  of the program will be stored.
# The 3rd is the short name for the type of build performed.   
# The 4th is the directory in which to store copies of the log files
#  with altered names encoding the status, as determined by the
#  BuildTypes module. If not provided then this isn't done.
#  Files will be stored in a subfolder machine.build_type.
#  The .log file name, without extension, will be prepended to the status.

import os, sys

def help():
  print "Usage:",sys.argv[0],"<test exe> <.log file> <build type> [output dir]"

# Check for valid arguments.
if len(sys.argv) < 4:
  print "Syntax error: insufficient arguments."
  help()
  sys.exit(1)

exefile, logfile, build_type = sys.argv[1:4]
if len(sys.argv) > 4:
  outputdir = sys.argv[4]
  if not os.path.isdir(outputdir):
    print "Output directory",outputdir,"does not exist."
    help()
    sys.exit(1)
  import socket, BuildTypes
  build = BuildTypes.GetBuildType(build_type)
else:
  outputdir = None
  build = None

# Run the test program and record output & exit code
test_fp = os.popen(exefile, 'r')
test_output = []
for line in test_fp:
  test_output.append(line)
  print line,
exit_code = test_fp.close()

# Write output to the log file
log_fp = file(logfile, 'w')
log_fp.writelines(test_output)
log_fp.close()

if outputdir:
  print "Storing..."
  # Get the test status and copy log file
  machine  = socket.getfqdn()
  test_dir = os.path.join(outputdir, machine+'.'+build_type)
  if not os.path.exists(test_dir):
    os.mkdir(test_dir)
  if not os.path.isdir(test_dir):
    print test_dir,"is not a directory; unable to copy output."
    sys.exit(1)
  test_name = os.path.splitext(logfile)[0]
  status    = build.EncodeStatus(exit_code, test_output)
  print test_name, status
  copy_to   = os.path.join(test_dir, test_name+'.'+status)
  print copy_to
  fp = file(copy_to, 'w')
  fp.writelines(test_output)
  fp.close()
