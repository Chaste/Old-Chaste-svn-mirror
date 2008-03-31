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
src_dirs.remove('apps/src')
src_dirs.remove('dealii/src')
src_files = []
for src_dir in src_dirs:
    for dirpath, dirnames, filenames in os.walk(src_dir):
       for filename in filenames:
           if filename[-4:] in ['.cpp', '.hpp']:
               src_files.append({'dir': dirpath, 'file': filename})

def coverage_ignore(src_file):
    """Whether to ignore the fact that a source file is not used.
    
    If a file contains only typedefs, for example, this is not an error. 
    For .hpp files we check this by looking for the presence of either
    'template' or 'class' at the start of a line.  If neither are found,
    we assume the file contains no real code.
    
    This will only work if header files don't contain non-template function
    definitions, which should be the case if we're being good programmers.
    Unfortunately the boost serialization tweaking file "TemplatedExport.hpp"
    has some templated definitions which are not code, for this reason we only
    scrape the file for "template" or "class" definitions that are not surrounded
    by COVERAGE_IGNORE.
    """
    ignore = False
    if src_file['file'] == 'triangle.cpp':
        # We don't try to cover other people's code
        ignore = True
    if src_file['file'][-4:] == '.hpp':
        ignore = True
        fp = open(os.path.join(src_file['dir'], src_file['file']))
        code = True
        for line in fp:
            if line.find('#define COVERAGE_IGNORE') != -1:
                code = False
            elif line.find('#undef COVERAGE_IGNORE') != -1:
                code = True
	    if code and (line.startswith('template') or line.startswith('class ')):
                ignore = False
                break
        fp.close()
    return ignore

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
        maybe_not_code, really_uncovered = False, False
        for line in lines:
            count, line_no, src_line = line.split(':', 2)
            count, line_no = count.strip(), line_no.strip()
            if src_line.find('#define COVERAGE_IGNORE') != -1:
                ignore = True
            elif src_line.find('#undef COVERAGE_IGNORE') != -1:
                ignore = False
            if line_no == 0:
                # This is a gcov header line; what it is doesn't matter
                out_file.write(line)
                break
            if count == '-':
                # This line "isn't code".  This may be because it's blank, a comment, or
                # similar.  Or it may be because it's within a templated method that hasn't
                # been instantiated in this particular execution, but it might be in another.
                maybe_not_code = True
            elif count == '#####':
                # The line was really uncovered here, so it must be code
                really_uncovered = True
            else:
                aggregated_count += int(count)
        else:
            if aggregated_count == 0:
                if maybe_not_code and not really_uncovered:
                    # This really wasn't a code line (or the template is *never* instantiated).
                    # Would be nice to differentiate these cases, but doing so is decidedly
                    # non-trivial.
                    aggregated_count = '-'
                else:
                    src_line_stripped = src_line.strip()
                    if not (ignore or src_line_stripped in ['{', '}', 'NEVER_REACHED;'] or
                            (src_line_stripped.startswith('return') and
                             src_line_stripped[6] in [';', ' ']) or
                            (src_line_stripped.startswith('catch ') and #Line is catch(...)
                             src_line_stripped[-1] == ')')
                             or (len(src_line_stripped)>0 and src_line_stripped[-1] == ')') ): #Method definition (possibly). Currently overlaps with previous 'catch' ignore
                        warn = False
                        aggregated_count = '#####'
                        #print 'Full details of coverage: ', src_line_stripped,'\t',src_file,'\t',aggregated_count,'\t', line_no,'\t', src_line
                    else:
                        aggregated_count = 'ignored'
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
    elif not ignore and missed_line_count == 0:
        out_file.write('\nOK!\n\n')
        status = 'OK'
    else:
        counts = (missed_line_count, missed_line_count+covered_line_count)
        out_file.write('\nFailed %d of %d tests\n\n' % counts)
        status = "%d_%d" % counts
        if warn:
            status = 'warn_' + status
        if ignore:
            status = 'ignore_' + status
    if coverage_ignore(src_file):
        # All special case ignorable files (not just ones with partial coverage)
        status = ''
    
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
