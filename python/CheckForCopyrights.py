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



# Check, apply or modify the copyright notice
import os, sys
exts = ['.cpp', '.hpp', '.py', '.java']
dir_ignores = ['build', 'cxxtest', 'testoutput', 'doc', 'projects']
exclusions = ['python/pycml/enum.py', 'python/pycml/pyparsing.py', 'python/pycml/schematron.py']

apply_update =  '-update' in sys.argv
apply_new = '-new' in sys.argv

chaste_dir = '.'
if '-dir' in sys.argv: 
    i = sys.argv.index('-dir')
    chaste_dir = os.path.realpath(sys.argv[i+1])

#This next variable is unused -- it shows the previous deprecated notice
#as a mechanism for rolling the notice forward
previous_depricated_notice="""Copyright (C) University of Oxford, 2005-2009

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

depricated_notice="""Copyright (C) University of Oxford, 2005-2009

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

current_notice="""Copyright (C) University of Oxford, 2005-2010

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

#py_depricated_notice=''
#for line in depricated_notice.splitlines():#
#	py_depricated_notice+=''.join(['# ',line,'\n'])
#py_current_notice=''
#for line in current_notice.splitlines():
#	py_current_notice+=''.join(['# ',line,'\n'])

py_current_notice='\"\"\"'+current_notice+'\"\"\"\n'
cpp_current_notice='/*\n\n'+current_notice+'\n*/'
cpp_depricated_notice='/*\n\n'+depricated_notice+'\n*/'


pycml_notice=" Processed by pycml - CellML Tools in Python"
xsd2_notice="// Copyright (C) 2005-2007 Code Synthesis Tools CC"
xsd3_notice="// Copyright (C) 2005-2008 Code Synthesis Tools CC"
triangle_notice="""/*  Copyright 1993, 1995, 1997, 1998, 2002, 2005                             */
/*  Jonathan Richard Shewchuk                                                */"""

def CheckForCopyrightNotice(findStr, fileIn):
    """Test if the (possibly multi-line) string findStr is contained anywhere in fileIn."""
    fileIn.seek(0)
    file_text = fileIn.read()
    return (file_text.find(findStr) >= 0)
    
def ReplaceStringInFile(findStr, repStr, filePath):
   """Replaces all findStr by repStr in file filePath"""
   tempName = filePath+'~'
   input = open(filePath)
   output = open(tempName, 'w')

   s = input.read()
   output.write(s.replace(findStr, repStr))
   output.close()
   input.close()
   os.rename(tempName, filePath)
   print 'Notice: replaced deprecated copyright notice in', filePath

def HeadAppendStringInFile(appendString, filePath):
   """Adds appendStr to the top of file filePath"""
   tempName = filePath+'~'
   input = open(filePath)
   output = open(tempName, 'w')

   s = input.read()
   output.write(appendString)
   output.write(s)
   output.close()
   input.close()
   os.rename(tempName, filePath)
   print 'Notice: applied copyright notice in ', filePath
    

   

def InspectFile(fileName):
    file_in = open(fileName)
    if fileName[-21:] == 'CheckForCopyrights.py':
    	#Can't really check this one, since it knows all the licences
    	return True
    valid_notice = False
    if (CheckForCopyrightNotice(cpp_current_notice, file_in) or
        CheckForCopyrightNotice(py_current_notice, file_in)):
        #print 'Found current notice in '+file_name
        valid_notice=True
    if (CheckForCopyrightNotice(pycml_notice, file_in) or
        CheckForCopyrightNotice(xsd2_notice, file_in) or
        CheckForCopyrightNotice(xsd3_notice, file_in) or
        CheckForCopyrightNotice(triangle_notice, file_in)):
        #print 'Found 3rd party notice in '+file_name
        if valid_notice:
            print "Multiple notices on", file_name
            return False
        else:
            return True
    if valid_notice:
        return True
    if CheckForCopyrightNotice(cpp_depricated_notice, file_in):
        print 'Found deprecated copyright notice for', fileName
        if apply_update:
           	ReplaceStringInFile(cpp_depricated_notice, cpp_current_notice, fileName)
        	return True
        else:
            print 'Fix this by doing: python python/CheckForCopyrights.py -update'
            return False
    
    print 'Found no copyright notice for', fileName
    if apply_new:
        if fileName[-3:] == '.py':
          	print 'Not implemented for .py files'
          	return False
        else:
	        HeadAppendStringInFile(cpp_current_notice+"\n\n", fileName)
        return True
    else:
        print 'Fix this by doing: python python/CheckForCopyrights.py -new'
        return False

num_no_copyrights = 0
num_copyrights = 0
chaste_dir_len = len(os.path.join(chaste_dir, ''))
for root, dirs, files in os.walk(chaste_dir):
    relative_root = root[chaste_dir_len:]
    # Check for ignored dirs
    for dir in dir_ignores:
        if dir in dirs:
            dirs.remove(dir)
    # Check for source files
    for file in files:
        relative_path = os.path.join(relative_root, file)
        name, ext = os.path.splitext(file)
        if ((ext in exts or file=='SConscript' or file=='SConstruct') and
            relative_path not in exclusions):
            file_name = os.path.join(root, file)
            if InspectFile(file_name) == False:
                num_no_copyrights += 1
            else:
                num_copyrights += 1

# Let the test summary script know
if chaste_dir == ".":
    dir = os.getcwd()
else:
    dir = chaste_dir
    
print "Copyright test run over ",dir," (",num_no_copyrights+num_copyrights,") files"
if num_no_copyrights > 0:
    print
    print "The next line is for the benefit of the test summary scripts."
    print "Failed",num_no_copyrights,"of",num_no_copyrights+num_copyrights,"tests"

    # Return a non-zero exit code if orphans were found
    import sys
    sys.exit(num_no_copyrights)
else:
    print "Infrastructure test passed ok."

