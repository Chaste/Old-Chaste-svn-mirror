#!/usr/bin/env python

# Copyright (C) Oxford University 2008
# 
# This file is part of Chaste.
# 
# CHASTE is free software: you can redistribute it and/or modify it
# under the terms of the Lesser GNU General Public License as
# published by the Free Software Foundation, either version 2.1 of the
# License, or (at your option) any later version.
#
# Chaste is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with Chaste.  If not, see <http://www.gnu.org/licenses/>.

"""Script to convert a tutorial test file into wiki code

Notes: 
1. The script starts converting everything after the first #define line
2. The script stops converting after an #endif. If this isn't the last
   line the script needs changing.
3. All C-style block comments '/*' to '*/' are converted to wiki text.
4. All other lines, including C++ style comments '//', are kept as code lines
5. In C-block comments (ie wiki text), all whitespace is removed. Hence 
   indentations for bulleted lists, subsections etc won't work
6. To print an empty line in the wiki page, write EMPTYLINE in the block
   comment, with nothing else
7. Lines inside a block comment which start with a '*', i.e.
     /* my comment is
      * two lines long */
   are ok, the initial '*' is removed.
"""


import sys
import glob
import os
import itertools

if len(sys.argv)!=3:
    print "Syntax error."
    print "Usage:",sys.argv[0],"<test_file> <output_file>"
    sys.exit(1)

test_file = sys.argv[1]
out_file_name = sys.argv[2]


if test_file[-4:] != '.hpp' or os.path.basename(test_file)[:4] != 'Test':
    print "Syntax error:",test_file,"does not appear to be a test file"
    sys.exit(1)

out_file = open(out_file_name, 'w')

out_file.write('This tutorial is automatically generated from the file '+test_file+'.\n')
out_file.write('Note that the code is given in full at the bottom of the page.\n\n\n')

# 'state machine' state variables
Parsing = False
ST_NONE, ST_TEXT, ST_CODE = 0, 1, 2 # Status codes
Status = ST_NONE

code_store = [] # a store of every code-line, to print at the end


in_file = open(test_file)
# TODO: Use a generator so this loop only processes interesting lines
for line in in_file:
    # remove trailing whitespace, including '\n'
    line = line.rstrip()
    # Remove all whitespace and save to a new string.
    # We don't remove the initial whitespace as it will
    # be needed for code lines
    stripped_line = line.strip()

    # We stop processing input after an #endif - if the test contains
    # an #endif before the final one the script needs changing.
    if stripped_line.startswith('#endif'):
	if Status is ST_CODE:
            # close code block
	    out_file.write('}}}\n')
        Parsing = False

    # if in Parsing mode
    if Parsing:
	if Status is ST_TEXT:
	    # we are still in a comment line, so strip it
	    line = line.strip()
	
	# check if the line is a new text line
        if stripped_line.startswith('/*'):
	    # assert not already text
	    assert Status is not ST_TEXT, 'Nested comment' 
	    # remove all whitespace and the '/*'
            line = line.strip()[2:]
	    line = line.strip()
	    # if the last line was code, print }}}
            if Status is ST_CODE:
     	        out_file.write('}}}\n')
	    # set the status as text
            Status = ST_TEXT
	elif Status is ST_TEXT and stripped_line.startswith('*'):
            # we are in a comment, so get rid of whitespace and the initial '*'
            line = line.strip()[1:]
	    line = line.strip()
	elif Status is ST_NONE and len(line.strip()) > 0:
	    # Line has content and isn't a comment => it's code
	    out_file.write('{{{\n#!c\n')
            Status = ST_CODE
	
	# check if comment ends
	if stripped_line.endswith('*/'):
	    # comment ended, verify was text
	    assert Status is ST_TEXT, 'Unexpected end comment'
	    # get rid of whitespace and '*/'
            line = line.strip()[:-2]
	    # set status as unknown 
            Status = ST_NONE

	# if the line is a comment justing saying 'EMPTYLINE', we print a blank line
	if line.strip()=='EMPTYLINE':
	    out_file.write('\n') 
	# else we print the line if is it non-empty
	elif len(line.strip()) > 0:
            out_file.write(line+'\n')
	
	# if the line is a code line we store it,
        # unless there would be two consecutive empty lines
        if Status is ST_CODE:
            if len(line.strip()) > 0 or len(code_store[-1].strip()) > 0:
		code_store.append(line)
	
    # we start processing lines AFTER the first #define..
    if stripped_line.startswith('#define'):
        Parsing = True


in_file.close()

# write out all the code
out_file.write('\n\n= Code =\nThe full code is given below\n{{{\n#!c\n')
for line in code_store:
    out_file.write(line+'\n')
out_file.write('}}}\n\n')

out_file.close()

