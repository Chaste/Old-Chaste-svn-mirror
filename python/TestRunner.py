#!/usr/bin/env python

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
    if exefile in ["python/CheckForOrphanedTests.py",
                   "python/CheckForDuplicateFileNames.py"]:
        command = exefile
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
            os.mkdir(test_dir)
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
        run_test(str(source[0]), str(target[0]), build, run_time_flags)
        return None

    return build_function
