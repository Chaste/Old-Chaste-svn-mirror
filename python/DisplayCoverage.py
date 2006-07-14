#!/usr/bin/env python

# Script to run gcov on source files after a Coverage build has been done,
# and summarise the results.

# The script takes arguments:
#  <output_dir>  The directory in which to generate summary files and
#                an index page.
#  <build_type>  The build type used; defaults to Coverage.

import sys
import glob
import os
import itertools
import BuildTypes

# Arguments to gcov
#  -l  Create long file names for included source files.
#  -p  Preserve complete path information in the names of generated .gcov files.
gcov_flags = ' -lp '

# Get output dir and build type object
if len(sys.argv) < 2:
    print "Syntax error."
    print "Usage:",sys.argv[0],"<test output dir> [<build type>]"
    sys.exit(1)
output_dir = sys.argv[1]
if len(sys.argv) > 2:
    build_type = sys.argv[2]
else:
    build_type = 'Coverage'
build = BuildTypes.GetBuildType(build_type)

# Remove any old output files/test results from output_dir
for filename in os.listdir(output_dir):
    os.remove(os.path.join(output_dir, filename))

# Find .gcda files to determine which source files to run gcov on
# First, find appropriate build directories
build_dirs = glob.glob('*/build/' + build.build_dir)
# Now find .gcda files within there
gcda_files = []
for build_dir in build_dirs:
    for dirpath, dirnames, filenames in os.walk(build_dir):
       for filename in filenames:
           if filename[-5:] == '.gcda':
               gcda_files.append({'dir': dirpath, 'file': filename})

# Run gcov on all the .cpp files which have .gcda files.
for gcda_file in gcda_files:
    # For added interest, the source file to process is in different locations
    # depending on whether it is a test or not.
    if gcda_file['file'][:4] == 'Test' or \
           gcda_file['dir'][-5:] == '/test':
        #gcda_file['dir'].find('/test/') != -1:
        # .cpp file is in the same folder
        os.system('gcov -o ' + gcda_file['dir'] + gcov_flags +
                  os.path.join(gcda_file['dir'], gcda_file['file'][:-4] + 'cpp'))
    else:
        # .cpp file is contained within the Chaste source tree
        # gcda_file['dir'] should look something like mesh/build/coverage/src/reader
        # We then want to look in mesh/src/reader
        try:
            start, end = gcda_file['dir'].split('src')
        except:
            print gcda_file
            raise
        toplevel, junk = start.split('build')
        # Get rid of slashes (or system equivalent)
        toplevel = os.path.dirname(toplevel)
        if end: end = end.split(os.path.sep, 1)[1]
        # Run gcov
        os.system('gcov -o ' + gcda_file['dir'] + gcov_flags +
                  os.path.join(toplevel, 'src', end, gcda_file['file'][:-4] + 'cpp'))

# Now find all our source files
src_dirs = glob.glob('*/src')
src_files = []
for src_dir in src_dirs:
    for dirpath, dirnames, filenames in os.walk(src_dir):
       for filename in filenames:
           if filename[-4:] in ['.cpp', '.hpp']:
               src_files.append({'dir': dirpath, 'file': filename})

for src_file in src_files:
    # Mangle the name like gcov does
    mangled_dir = src_file['dir'].replace(os.path.sep, '#')
    # Find .gcov files relating to this source file
    gcov_files = glob.glob('*' + mangled_dir + '#' + src_file['file'] + '.gcov')
    # Open all the files, and an output file
    gcov_fps = [open(gcov_file) for gcov_file in gcov_files]
    out_file_name = os.path.join(output_dir, mangled_dir + '#' + src_file['file'])
    out_file_name = out_file_name.replace('#', '-')
    out_file = open(out_file_name, 'w')
    # Now go through them line by line in lock-step,
    # aggregating line execution counts
    covered_line_count, missed_line_count, warn, ignore = 0, 0, True, False
    for lines in itertools.izip(*gcov_fps):
        aggregated_count = 0
        for line in lines:
            count, line_no, src_line = line.split(':', 2)
            count, line_no = count.strip(), line_no.strip()
            if src_line.strip() == '#define COVERAGE_IGNORE':
                ignore = True
            elif src_line.strip() == '#undef COVERAGE_IGNORE':
                ignore = False
            if count == '-' or line_no == 0:
                out_file.write(line)
                break
            if count != '#####':
                aggregated_count += int(count)
        else:
            if aggregated_count == 0:
                src_line_stripped = src_line.strip()
                if not (ignore or src_line_stripped == '}' or
                        (src_line_stripped.startswith('return') and src_line_stripped[6] in [';', ' '])):
                    warn = False
                aggregated_count = '#####'
                missed_line_count += 1
            else:
                covered_line_count += 1
            out_file.write("%9s:%5s:%s" % (aggregated_count, line_no, src_line))
    # Output a summary
    if not gcov_files:
        # No gcov files found for this source file.
        # This may not be an error, if the source file in question is an .hpp file with
        # an associated .cpp file containing all the code for the class.
        ##print src_file
        if src_file['file'][-4:] == '.hpp' and \
            os.path.exists(os.path.join(src_file['dir'], src_file['file'][:-3]+'cpp')):
            status = '' # So output file will be deleted
        else:
            out_file.write("This source file wasn't used at all!\n\nFailed 1 of 1 test\n")
            status = "1_1"
    elif missed_line_count == 0:
        out_file.write('\nOK!\n\n')
        status = 'OK'
    else:
        counts = (missed_line_count, missed_line_count+covered_line_count)
        out_file.write('\nFailed %d of %d tests\n\n' % counts)
        status = "%d_%d" % counts
        if warn:
            status = 'warn_' + status
    # Close all files
    [fp.close() for fp in gcov_fps]
    out_file.close()
    # Alter file name to indicate summary
    if status:
        os.rename(out_file_name, out_file_name + '.' + status + '.0')
    else:
        os.remove(out_file_name)

# Now remove .gcov files from the Chaste root directory
for filename in os.listdir('.'):
    if filename[-5:] == '.gcov':
        os.remove(filename)

# And generate a summary page
os.system('python python/DisplayTests.py '+output_dir+' '+build_type)
