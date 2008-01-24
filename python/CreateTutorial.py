#!/usr/bin/env python

# Script to convert a tutorial test file into wiki code
#
# Notes: 
# 1. The script starts converting everything after the first #define line
# 2. The script stops converting after an #endif. If this isn't the last
#    line the script needs changing.
# 3. All C-style block comments '/*' to '*/' are converted to wiki text.
# 4. All other lines, including C++ style comments '//', are kept as code lines
# 5. In C-block comments (ie wiki text), all whitespace is removed. Hence 
#    indentations for bulleted lists, subsections etc won't work
# 6. To print an empty line in the wiki page, write EMPTYLINE in the block
#    comment, with nothing else
# 7. Lines inside a block comment which start with a '*', i.e.
#      /* my comment is
#       * two lines long */
#    are ok, the initial '*' is removed.



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

# todo: get rid of path in file name then do this
#if test_file[:4]!='Test':
#    print "Syntax error:",test_file,"does not appear to be a test file"
#    sys.exit(1)

if test_file[-4:]!='.hpp':
    print "Syntax error:",test_file,"does not appear to be a test file"
    sys.exit(1)

out_file = open(out_file_name, 'w')

out_file.write('This tutorial is automatically generated from the file '+test_file+'.\n')
out_file.write('Note that the code is given in full at the bottom of the page.\n\n\n')

Parsing=False
Parse_from_next_line=False
Status=0 #0 for none, 1 for Text, 2 for Code.
code_store=[''] # a store of every code-line, to print at the end


in_file = open(test_file)
for line in in_file:
    # remove trailing whitespace, including '\n'
    line = line.rstrip()
    # Remove all whitespace and save to a new string.
    # We don't remove the initial whitespace as it will
    # be needed for code lines
    stripped_line = line.strip()
	
    # we start printing line AFTER the first #define..
    if stripped_line.startswith('#define'):
        Parse_from_next_line=True

    # ..and stop after an #endif - if the test contains
    # an #endif before the final one the script needs changing
    if stripped_line.startswith('#endif'):
	if Status==2:
            # close code block
	    out_file.write('}}}\n')
        Parsing=False

    # if in Parsing mode
    if Parsing:
	if Status==1:
	    # we are still in a comment line, so strip it
	    line = line.strip()
	
	# check if the line is a new text line
        if stripped_line.startswith('/*'):
	    # assert not already text
	    assert Status!=1, 'StatusError1' 
	    # remove all whitespace and the '/*'
	    line = line.strip()
	    line = line[2:]           
	    line = line.strip()
	    # if line line was code, print }}}
            if Status==2:             
     	        out_file.write('}}}\n') 
	    # set the status as text
            Status=1
	elif stripped_line.startswith('*') and Status==1:
            # we are in a comment, so get rid of whitespace and the initial '*'
	    line = line.strip()
	    line = line[1:]
	    line = line.strip()
	elif Status==0 and len(line.strip())!=0:
	    # status=none, and newline isn't a comment=>code
	    out_file.write('{{{\n#!c\n')
            Status=2        
	
	# check if comment ends
	if stripped_line.endswith('*/'):
	    # comment ended, verify was text
	    assert Status==1, 'StatusError2'
	    # get rid of whitespace and '*/'
	    line = line.strip()
	    line = line[:-2]
	    # set status as unknown 
            Status=0                  

	# if the line is a comment justing saying 'EMPTYLINE', we print a blank line
	if line.strip()=='EMPTYLINE':
	    out_file.write('\n') 
	# else we print the line if is it non-empty
	elif len(line.strip())!=0:
            out_file.write(line+'\n')
	
	# if the line is a code line we store it, unless there would be two consecutive
        # empty lines
	if Status==2:
            if not len(line.strip())==0 or not len(code_store[-1].strip())==0:
		code_store.append(line)

    if Parse_from_next_line:
	Parsing=True

in_file.close()

# write out all the code
out_file.write('\n\n= Code =\nThe full code is given below\n{{{\n#!c\n')
for line in code_store:
    out_file.write(line+'\n')
out_file.write('}}}\n\n')

out_file.close()

