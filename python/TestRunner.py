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


# Script to run a test executable and record the status
# (whether all tests passed or how many failed) as part of
# the output filename.

# The test executable should be run and the output captured,
# and stored to a file.

# This script expects 3 - 5 arguments.
# The first is the executable to run.
# The second is the name of the output .log file, in which the output
#  of the program will be stored.
# The 3rd is the short name for the type of build performed.
# The 4th, if present, is any extra flags to pass to the executable.

# The directory in which to store copies of the log files
# with altered names encoding the status, as determined by the
# BuildTypes module, is determined from the build object.
# It will be created if necessary.
# The .log file basename, without extension, will be prepended to the status.

import os, sys, time

def help():
    print "Usage:",sys.argv[0],\
        "<test exe> <.log file> <build type> [run time flags] [--no-stdout]"

def run_test(exefile, logfile, build, run_time_flags='', echo=True):
    """Actually run the given test."""
    # Find out how we're supposed to run tests under this build type
    if exefile.startswith("python/CheckFor"):
        command = exefile + ' 2>&1'
    elif exefile.startswith('python/test/Test'):
        command = './python/TestPythonCode.py ' + exefile + ' 2>&1'
    else:
        command = build.GetTestRunnerCommand(exefile, '2>&1 ' + run_time_flags)

    # Run the test program and record output & exit code
    log_fp = file(logfile, 'w')
    if not log_fp:
        raise IOError("Unable to open log file")
    start_time = time.time()
    test_fp = os.popen(command, 'r')
    for line in test_fp:
        log_fp.write(line)
        if echo:
            print line,
    exit_code = test_fp.close()
    end_time = time.time()
    log_fp.close()

    #print "Time",end_time,start_time

    if True:
        import glob
        # Get the test status and copy log file
        test_dir = build.output_dir
        if not os.path.exists(test_dir):
            os.makedirs(test_dir)
        if not os.path.isdir(test_dir):
            print test_dir, "is not a directory; unable to copy output."
            sys.exit(1)
        test_name = os.path.splitext(os.path.basename(logfile))[0]
        log_fp = open(logfile, 'r')
        status = build.EncodeStatus(exit_code, log_fp)
        log_fp.close()
        #print test_name, status
        # Remove any old copies of results from this test
        oldfiles = glob.glob(os.path.join(test_dir, test_name+'.*'))
        for oldfile in oldfiles:
            os.remove(oldfile)
        # Copy the new results
        copy_to = build.ResultsFileName(dir=test_dir, testsuite=test_name, status=status,
                                        runtime=end_time-start_time)
        #print copy_to
        os.system("/bin/cp -f " + logfile + " " + copy_to)
        return copy_to


if __name__ == '__main__':
    # Check for valid arguments.
    if '--no-stdout' in sys.argv:
        echo = False
        sys.argv.remove('--no-stdout')
    else:
        echo = True

    if len(sys.argv) < 4:
        print "Syntax error: insufficient arguments."
        help()
        sys.exit(1)

    exefile, logfile, build_type = sys.argv[1:4]

    # Get any extra command line args for the test
    if len(sys.argv) > 4:
        run_time_flags = sys.argv[4]
    else:
        run_time_flags = ''

    # Get a build object for this build type
    import BuildTypes
    build = BuildTypes.GetBuildType(build_type)

    run_test(exefile, logfile, build, run_time_flags, echo)

# Builder function for running via SCons
def get_build_function(build, run_time_flags=''):
    """Return a function that can be used as a Builder by SCons."""
    
    def build_function(target, source, env):
        # Set up the environment from env['ENV']
        os.environ.update(env['ENV'])
        # Run the test
        log = run_test(str(source[0]), str(target[0]), build, run_time_flags)
        # Note the extra dependency of the copied log file
        #env.SideEffect(log, target)
        env.Depends(os.path.join(os.path.dirname(log), 'index.html'), log)
        return None

    return build_function
